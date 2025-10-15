#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"
#include "../Hydro/kernel.h"


#ifdef GM_STARDENSITY
#define CM_PER_KPC                ((double)3.085678e21)
#define CM_PER_PC                 ((double)3.085678e18)
#define SOLAR_MASS_IN_CODE_UNITS  (SOLAR_MASS / All.UnitMass_in_g)

#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#define LOCK_PARTNODEDRIFT
#define UNLOCK_PARTNODEDRIFT


struct kernel_density
{
  double dx, dy, dz;
  double r;
  double dvx, dvy, dvz; 
  double wk, dwk;
  double hinv, hinv3, hinv4;
  double mj_wk, mj_dwk_r;
};

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct densdata_in_stars
{
  MyLongDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *DensDataInStars, *DensDataGetStars;

static struct densdata_out_stars
{
  MyLongDouble Rho;
  MyLongDouble Ngb;
#ifdef EVALPOTENTIAL
  MyLongDouble SmoothPot;
#endif
}
 *DensDataResultStars, *DensDataOutStars;

void out2particle_density_stars(struct densdata_out_stars *, int , int );
void particle2in_density_stars(struct densdata_in_stars *, int );


void particle2in_density_stars(struct densdata_in_stars *in, int i)
{
  int k;

  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

#ifdef GM_STARDENSITY_FIXAPERTURE
  P[i].StarHsml = CM_PER_KPC*GM_STARDENSITY_FIXAPERTURE/All.UnitLength_in_cm;
  if(All.ComovingIntegrationOn)
	   P[i].StarHsml *= All.Time/All.HubbleParam;
#endif
  in->Hsml = P[i].StarHsml;

}

void out2particle_density_stars(struct densdata_out_stars *out, int i, int mode)
{
  int k, j;

  ASSIGN_ADD(P[i].StarNumNgb, out->Ngb, mode);

  if(P[i].Type == 4)
    {
      ASSIGN_ADD(P[i].StarDensity, out->Rho, mode);
#ifdef EVALPOTENTIAL
      ASSIGN_ADD(P[i].SmoothPot, out->SmoothPot, mode);
#endif
  

    }

}

/*! \file density_stars.c
 *  \brief SPH stellar density computation and smoothing length determination
 *
 *  Identical to "density.c", but used to compute the stellar density around star particle  
 *  and the smoothing length
 *  Number of neighbours and deviations are the same of the gas particle 
 *
 *  If EVALPOTENTIAL is ON,  it also evaluated the gravitational potential
 *  SPH-smoothed on STAR PARTICLES ONLY 
 */


/*! 
This function is a clone of density(). Evaluates the SPH smoothed stellar density for each star particles.
Also finds the smoothing lenght.
 */
void density_stars(void)
{
  MyFloat *Left, *Right;
  int i, j, k, ii, ndone, ndone_flag, npleft, iter = 0;
  int ngrp, recvTask, place;
  long long ntot;
  double fac;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  double desnumngb, desnumngbdev;
  int save_NextParticle;
  long long n_exported = 0;
  int redo_particle;

  int numstars=0, totnumstars=0;


  CPU_Step[CPU_STARDENSMISC] += measure_time();

  int NTaskTimesNumPart;



  /*  check that enough stars do exist for the neighbour search 
      (at high redshift this is not true */
  /* I ask for 10 times the number of neighbours */
  for(i=0; i<NumPart; i++)
    if(P[i].Type==4)
      {
      numstars++;
      if (P[i].StarHsml>=0.)
	P[i].StarHsml=All.SofteningTable[4];
      }
  MPI_Allreduce(&numstars, &totnumstars, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  if (totnumstars < 10 * All.DesNumNgb)
    {
      if (ThisTask==0)
	printf(" density_stars: not enough star particles for neighbour search (we have %d stars)\n",totnumstars); fflush(stdout);
      return;
    }
  else
    if (ThisTask==0)
      printf(" ...%d stars available\n",totnumstars); fflush(stdout);

  
  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));


  
#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
  int il;
#pragma omp parallel for private(i)
  for(il = 0; il < NActivePart; il++)
    {
      i = ActiveParticleList[il];
#else /* KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif
      if(density_isactive_stars(i))
	{
	  Left[i] = Right[i] = 0;
	}			// if
    }				// for

  /* allocate buffers to arrange communication */

#ifdef KD_BUFFER_MANAGEMENT
  size_t MyBufferSize = KD_BUFFER_MANAGEMENT * FreeBytes / (1024.0 * 1024.0);
#else
  size_t MyBufferSize = All.BufferSize;
#endif

  All.BunchSize =
    (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					   sizeof(struct densdata_in_stars) + sizeof(struct densdata_out_stars) +
					   sizemax(sizeof(struct densdata_in_stars), sizeof(struct densdata_out_stars))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

#ifdef KD_BUFFER_MANAGEMENT
  if(ThisTask == 0)
    printf("Density_stars: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
	   FreeBytes / (1024.0 * 1024.0));
#endif

  t0 = second();

  desnumngb = All.DesNumNgb*GM_STARDENSITY_NEIGH;
  desnumngbdev = All.MaxNumNgbDeviation;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {

      NextParticle = FirstActiveParticle;	/* begin with this index */

      do
	{
	  BufferFullFlag = 0;
	  Nexport = 0;
	  save_NextParticle = NextParticle;

	  tstart = second();

#ifdef _OPENMP
#pragma omp parallel
#endif
	  {
#ifdef _OPENMP
	    int mainthreadid = omp_get_thread_num();
#else
	    int mainthreadid = 0;
#endif
	    density_evaluate_primary_stars(&mainthreadid);	/* do local particles and prepare export list */
	  }


	  tend = second();
	  timecomp1 += timediff(tstart, tend);

	  if(BufferFullFlag)
	    {
	      int last_nextparticle = NextParticle;

	      NextParticle = save_NextParticle;

	      while(NextParticle >= 0)
		{
		  if(NextParticle == last_nextparticle)
		    break;

		  if(ProcessedFlag[NextParticle] != 1)
		    break;

		  ProcessedFlag[NextParticle] = 2;

		  NextParticle = NextActiveParticle[NextParticle];
		}

	      if(NextParticle == save_NextParticle)
		{
		  /* in this case, the buffer is too small to process even a single particle */
		  printf("Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n", ThisTask, P[NextParticle].Type,
			 P[NextParticle].Pos[0], P[NextParticle].Pos[1], P[NextParticle].Pos[2],
			 P[NextParticle].Mass);
		  if(P[NextParticle].Type == 0)
		    printf("   rho=%g hsml=%g\n", P[NextParticle].StarDensity, P[NextParticle].StarHsml);

		  endrun(12998);
		}


	      int new_export = 0;

	      for(j = 0, k = 0; j < Nexport; j++)
		if(ProcessedFlag[DataIndexTable[j].Index] != 2)
		  {
		    if(k < j + 1)
		      k = j + 1;

		    for(; k < Nexport; k++)
		      if(ProcessedFlag[DataIndexTable[k].Index] == 2)
			{
			  int old_index = DataIndexTable[j].Index;

			  DataIndexTable[j] = DataIndexTable[k];
			  DataNodeList[j] = DataNodeList[k];
			  DataIndexTable[j].IndexGet = j;
			  new_export++;

			  DataIndexTable[k].Index = old_index;
			  k++;
			  break;
			}
		  }
		else
		  new_export++;

	      Nexport = new_export;

	    }


	  n_exported += Nexport;

	  for(j = 0; j < NTask; j++)
	    Send_count[j] = 0;
	  for(j = 0; j < Nexport; j++)
	    Send_count[DataIndexTable[j].Task]++;

	  MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

	  tstart = second();

	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);

	  tend = second();
	  timewait1 += timediff(tstart, tend);

	  for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	    {
	      Nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

	  DensDataGetStars = (struct densdata_in_stars *) mymalloc("DensDataGetStars", Nimport * sizeof(struct densdata_in_stars));
	  DensDataInStars = (struct densdata_in_stars *) mymalloc("DensDataInStars", Nexport * sizeof(struct densdata_in_stars));

	  /* prepare particle data for export */
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      particle2in_density_stars(&DensDataInStars[j], place);

	      memcpy(DensDataInStars[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	    }
	  /* exchange particle data */
	  tstart = second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&DensDataInStars[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in_stars), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGetStars[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in_stars), MPI_BYTE,
				   recvTask, TAG_DENS_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(DensDataInStars);
	  DensDataResultStars =
	    (struct densdata_out_stars *) mymalloc("DensDataResultStars", Nimport * sizeof(struct densdata_out_stars));
	  DensDataOutStars =
	    (struct densdata_out_stars *) mymalloc("DensDataOutStars", Nexport * sizeof(struct densdata_out_stars));

	  report_memory_usage(&HighMark_sphdensity, "STAR_DENSITY");

	  /* now do the particles that were sent to us */

	  tstart = second();

	  NextJ = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
	  {
#ifdef _OPENMP
	    int mainthreadid = omp_get_thread_num();
#else
	    int mainthreadid = 0;
#endif
	    density_evaluate_secondary_stars(&mainthreadid);
	  }


	  tend = second();
	  timecomp2 += timediff(tstart, tend);

	  if(NextParticle < 0)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  tstart = second();
	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
	  tend = second();
	  timewait2 += timediff(tstart, tend);


	  /* get the result */
	  tstart = second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&DensDataResultStars[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_out_stars),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DensDataOutStars[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_out_stars),
				   MPI_BYTE, recvTask, TAG_DENS_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }
	  tend = second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the result to the local particles */
	  tstart = second();
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
	      out2particle_density_stars(&DensDataOutStars[j], place, 1);
	    }
	  tend = second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(DensDataOutStars);
	  myfree(DensDataResultStars);
	  myfree(DensDataGetStars);
	}
      while(ndone < NTask);


      /* do final operations on results */
      tstart = second();

      npleft = 0;

#ifdef KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP
      int il;
#pragma omp parallel for reduction(+:npleft) private(i,k,desnumngb,desnumngbdev,redo_particle,ii,fac) if(NActivePart>80)
      for(il = 0; il < NActivePart; il++)
	{
	  i = ActiveParticleList[il];
#else /* KD_ACTIVE_PARTICLE_LIST_FOR_OPENMP */
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
#endif
	  if(density_isactive_stars(i))
	    {

	      P[i].SmoothPot /= P[i].StarDensity; /* smoothed gravitational potential */
 
	      
	      /* now check whether we had enough neighbours */
	      desnumngb = All.DesNumNgb*GM_STARDENSITY_NEIGH;
	      desnumngbdev = All.MaxNumNgbDeviation;


	      redo_particle = 0;

#if !defined(GM_STARDENSITY_FIXAPERTURE) /* if this is defined, the stellar Hsml is fixed*/
	      if(P[i].StarNumNgb < (desnumngb - desnumngbdev) ||
		 P[i].StarNumNgb > (desnumngb + desnumngbdev)   )
		redo_particle = 1;
#endif

#ifdef MAXHSML
	      if(P[i].StarNumNgb < (desnumngb - desnumngbdev) && P[i].StarHsml >= All.MaxHsml)
		{
		  P[i].StarHsml = All.MaxHsml;
		  redo_particle = 0;
		}
#endif

	      if(redo_particle)
		{
		  /* need to redo this particle */
		  npleft++;
		  /* added variable my_continue because the openmp multithreading does not allow "continue" statement here */
		  short my_continue = 0;

		  if(!my_continue)
		    if(Left[i] > 0 && Right[i] > 0)
		      if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
			{
			  /* this one should be ok */
			  npleft--;
			  P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
			  my_continue = 1;
			}
		  if(!my_continue)
		    if(P[i].StarNumNgb < (desnumngb - desnumngbdev))
		      Left[i] = DMAX(P[i].StarHsml, Left[i]);
		    else
		      {
			if(Right[i] != 0)
			  {
			    if(P[i].StarHsml < Right[i])
			      Right[i] = P[i].StarHsml;
			  }
			else
			  Right[i] = P[i].StarHsml;
		      }

		  if(!my_continue)
		    if(iter >= MAXITER - 10)
		      {
			printf
			  ("i=%d task=%d ID=%llu StarHsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			   i, ThisTask, (unsigned long long) P[i].ID, P[i].StarHsml, Left[i], Right[i],
			   (float) P[i].NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1],
			   P[i].Pos[2]);
			fflush(stdout);
		      }

		  if(!my_continue)
		    if(Right[i] > 0 && Left[i] > 0)
		      P[i].StarHsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		    else
		      {
			if(Right[i] == 0 && Left[i] == 0)
			  {
			    char buf[1000];
			    printf("Right[i] == 0 && Left[i] == 0 P[i].StarHsml=%g\n", P[i].StarHsml);
			    exit(-11732);
			  }

			if(Right[i] == 0 && Left[i] > 0)
			  {
			    P[i].StarHsml *= 1.26;
			  }

			if(Right[i] > 0 && Left[i] == 0)
			  {
			      P[i].StarHsml /= 1.26;
			  }
		      }

		}		// redo particle
	      else
		P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
	    }			// if
	}			// for
      tend = second();
      timecomp1 += timediff(tstart, tend);

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
	{
	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("star density ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
   }
  while(ntot > 0);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Right);
  myfree(Left);
  myfree(Ngblist);


  /* mark as active again */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].TimeBin < 0)
	P[i].TimeBin = -P[i].TimeBin - 1;
    }



  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_STARCOMPUTE] += timecomp;
  CPU_Step[CPU_STARWAIT] += timewait;
  CPU_Step[CPU_STARCOMM] += timecomm;
  CPU_Step[CPU_STARDENSMISC] += timeall - (timecomp + timewait + timecomm);
}


/*! 
 *  Identical to density_evaluate, but for stellar densities around stars
 *  
 */
int density_evaluate_stars(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		     int *ngblist)
{
  int k, j, n;
  int startnode, numngb_inbox, listindex = 0;
  double r2, h2, u, mass_j;

  struct kernel_density kernel;
  struct densdata_in_stars local;
  struct densdata_out_stars out;


  memset(&out, 0, sizeof(struct densdata_out_stars));


  if(mode == 0)
    particle2in_density_stars(&local, target);
  else
    local = DensDataGetStars[target];

  h2 = local.Hsml * local.Hsml;

  kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGetStars[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  /* XXXX controlla che sta roba qui sotto trovi tutto e non solo il  gas */
	  numngb_inbox =
	    ngb_treefind_variable_threads_stars(local.Pos, local.Hsml, target, &startnode, mode, exportflag,
					  exportnodecount, exportindex, ngblist);

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = ngblist[n];

	      kernel.dx = local.Pos[0] - P[j].Pos[0];
	      kernel.dy = local.Pos[1] - P[j].Pos[1];
	      kernel.dz = local.Pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      kernel.dx = NEAREST_X(kernel.dx);
	      kernel.dy = NEAREST_Y(kernel.dy);
	      kernel.dz = NEAREST_Z(kernel.dz);
#endif
	      r2 = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;

	      if(r2 < h2)
		{
		  kernel.r = sqrt(r2);
		  u = kernel.r * kernel.hinv;

		  kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);

		  mass_j = P[j].Mass;

		  kernel.mj_wk = FLT(mass_j * kernel.wk);
		  out.Rho += kernel.mj_wk;
		  out.Ngb += FLT(NORM_COEFF * kernel.wk / kernel.hinv3);    /* NORM_COEFF= 4.0/3 * PI = 4.188790204786 */

#ifdef EVALPOTENTIAL
		  out.SmoothPot += mass_j * kernel.wk * P[j].Potential;
#endif
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DensDataGetStars[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    out2particle_density_stars(&out, target, 0);
  else
    DensDataResultStars[target] = out;

  return 0;
}


void *density_evaluate_primary_stars(void *p)
{
  int thread_id = *(int *) p;
  int i, j;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	if(BufferFullFlag != 0 || NextParticle < 0)
	  {
	    exitFlag = 1;
	  }
	else
	  {
	    i = NextParticle;
	    ProcessedFlag[i] = 0;
	    NextParticle = NextActiveParticle[NextParticle];
	  }
      }
      UNLOCK_NEXPORT;
      if(exitFlag)
	break;

      if(density_isactive_stars(i))
	{
	  if(density_evaluate_stars(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *density_evaluate_secondary_stars(void *p)
{
  int thread_id = *(int *) p;

  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;


  while(1)
    {
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	j = NextJ;
	NextJ++;
      }
      UNLOCK_NEXPORT;

      if(j >= Nimport)
	break;

      density_evaluate_stars(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}


int density_isactive_stars(int n)
{
  if(P[n].TimeBin < 0)
    return 0;

  if(P[n].Type == 4)
    return 1;

  return 0;
}

int ngb_treefind_variable_threads_stars(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				  int mode, int *exportflag, int *exportnodecount, int *exportindex,
				  int *ngblist)
{
  int numngb, no, nexp, p, task;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;

#ifndef DO_NOT_BRACH_IF
  MyLongDouble dx, dy, dz, dist;
#ifdef PERIODIC
  MyLongDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  /* GM: note: this is the only difference with the normal  function
	     I prefer to clone it to make of this file a standalone module */

	  if(P[p].Type != 4)
	    continue;

	  if(P[p].Ti_current != ti_Current)
	    {
	      LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      drift_particle(p, ti_Current);
	      UNLOCK_PARTNODEDRIFT;
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif
	  ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      exportflag[task] = target;
		      exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(exportnodecount[task] == NODELISTLENGTH)
		    {
		      int exitFlag = 0;
		      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
		      {
			if(Nexport >= bunchSize)
			  {
			    /* out of buffer space. Need to discard work for this particle and interrupt */
			    BufferFullFlag = 1;
			    exitFlag = 1;
			  }
			else
			  {
			    nexp = Nexport;
			    Nexport++;
			  }
		      }
		      UNLOCK_NEXPORT;
		      if(exitFlag)
			return -1;

		      exportnodecount[task] = 0;
		      exportindex[task] = nexp;
		      DataIndexTable[nexp].Task = task;
		      DataIndexTable[nexp].Index = target;
		      DataIndexTable[nexp].IndexGet = nexp;
		    }

		  DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;

		}

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    {
	      LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      force_drift_node(no, ti_Current);
	      UNLOCK_PARTNODEDRIFT;
	    }

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}

/* This function does the BH repositioning, on the maximum of the stellar density,
   with the required options */
void blackhole_stardens_reposition(int n)
{
  double pbh[3], stardens_in_msun_per_pc3;
#ifdef GM_USE_STARDENS_AND_POTMIN
  double dd1, dd2;
#endif
#ifdef GM_BHMOVE_WT13
  double dxvers, dyvers, dzvers, vmod, rr, shift1, shift2, dt_drift;
#endif

  if(BPP(n).BH_StarDens>0 && P[n].Mass>0.) /* second condition for the (unlikely) case when BH has no valid star particle neighbour */
    {                                        /* in that case, it does not move */

      /* converting stardensity in Msun/pc^3 */
      stardens_in_msun_per_pc3 = BPP(n).BH_StarDens / pow( All.UnitLength_in_cm, 3);  // codemass/cm3
      stardens_in_msun_per_pc3 *=   CM_PER_PC * CM_PER_PC *  CM_PER_PC   / SOLAR_MASS_IN_CODE_UNITS;     // Msun/pc^3
      if (All.ComovingIntegrationOn) //to physical
	stardens_in_msun_per_pc3 *= All.HubbleParam*All.HubbleParam / pow(All.Time,3);


      if(stardens_in_msun_per_pc3 > GM_STARDENSITY_THRESH) /* the BH only is repositioned if max star density is above threshold, i.e. it does not move if it is in an underdense environment */
	{
#ifdef GM_BHMOVE_WT13 // shifts BHs according to formula (21) in  Wurster & Thacker 2013
	  rr = (P[n].Pos[0] - BPP(n).BH_StarDensPos[0]) * (P[n].Pos[0] - BPP(n).BH_StarDensPos[0]) +
	    (P[n].Pos[1] - BPP(n).BH_StarDensPos[1]) * (P[n].Pos[1] - BPP(n).BH_StarDensPos[1]) +
	    (P[n].Pos[2] - BPP(n).BH_StarDensPos[2]) * (P[n].Pos[2] - BPP(n).BH_StarDensPos[2]);
	  dxvers= (BPP(n).BH_StarDensPos[0] - P[n].Pos[0])/rr;
	  dyvers= (BPP(n).BH_StarDensPos[1] - P[n].Pos[1])/rr;
	  dzvers= (BPP(n).BH_StarDensPos[2] - P[n].Pos[2])/rr;

	  vmod= sqrt(P[n].Vel[0]*P[n].Vel[0] +  P[n].Vel[1]*P[n].Vel[1] +  P[n].Vel[2]*P[n].Vel[2]);
	  if(All.ComovingIntegrationOn)
	    dt_drift = get_drift_factor(P[n].Ti_current, All.Ti_Current);
	  else
	    dt_drift = (All.Ti_Current - P[n].Ti_current) * All.Timebase_interval;

	  shift1= 0.01*All.SofteningTable[5];
	  shift2 = 0.03*vmod*dt_drift;
	  if(shift1<=shift2)
	    {		  
	      pbh[0] = P[n].Pos[0] + shift1*dxvers;
	      pbh[1] = P[n].Pos[1] + shift1*dyvers;
	      pbh[2] = P[n].Pos[2] + shift1*dzvers;
	    }
	  else
	    {
	      pbh[0] = P[n].Pos[0] + shift2*dxvers;
	      pbh[1] = P[n].Pos[1] + shift2*dyvers;
	      pbh[2] = P[n].Pos[2] + shift2*dzvers;
	    }
	  
	}
#else
	  for( int k=0; k<3; k++)
	    pbh[k]=BPP(n).BH_StarDensPos[k];
        }
#endif
	      
	  for(int k = 0; k < 3; k++)                   
	    {
	      P[n].Pos[k] = pbh[k];
#ifdef GM_USE_VELS_ON_REPOSITIONING
	      P[n].Vel[k] = BPP(n).BH_StarDensVel[k];
#endif
	    }

	} //ends the density threshold if
  
  
}

#endif /* closes GM_STARDENSITY */


  /********************** Config Options **************************************************
 #--------------------------------------- BH merger fix switches
GM_USE_ABSVAL_IN_VREL      #BH mergers only happen if their relative vel is smaller than
                           #All.BlackHoleMergeCsndFrac set in the param file. In this case
                           #SOUND speed is NOT used

 #--------------------------------------- Stellar density switches
GM_STARDENSITY                  #master switch for the evaluation of stellar density around stars
#GM_STARDENSITY_FIXAPERTURE=15   #if set, a PHYSICAL radius is used for the stellar density
                                #otherwise, the a fixed number of stars is used.
                                #in output you have the number of neighbours (1st case) or the stellar 
                                #hsml (2nd). You always have RHOS (rho star) in output
GM_STARDENSITY_NEIGH=0.5          #number of star neighbour for the above case, WITH RESCPECT TO THE GAS NEIGHBOURS
                                # maybe <1. e.g. if 200 Ngb for gas and this is 0.5, will require 100 stars

GM_REPOSITION_ON_STARDENSITY_MAX #master switch for continuosly pinning BHs on the top of the 
                                 #star particle having the maximum stellar density
                                 #note that this requires All.BlackHoleRepositioningOn set to 1 in parameter file
GM_STARDENSITY_THRESH=0.004      #density threshold for repositioning BHs (msun/pc^3)
                                 #Bovi 2017: in solar neighbourhood is 0.04 Msun/pc^3. This  is one tenth of solar
GM_VREL_REPOSITION=60            #relative velocity threshold for repositioning BHs (km/s)
#GM_BHMOVE_WT13                   #alternative repositioning as in Wurster & Thacker 2013       
GM_USE_VELS_ON_REPOSITIONING     #the BH velocity is set to the star particle velocity upon which the BH moves   
**********************************End of Config options*********************************/
