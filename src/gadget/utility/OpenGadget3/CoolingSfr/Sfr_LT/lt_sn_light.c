/*
* @file
* This file is part of the developer version of GADGET3 and contains
* the license conditions for its usage.
*
* @author GADGET-development team, led by Volker Springel and Klaus Dolag.
*
* @section LICENSE
* Copyright (c) 2016, Volker Springel, Klaus Dolag, and all contributing authors
* (see change logs). All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Received source code may be modified and used as convenient.
*
* 2. Redistributions of source code or in binary form is only possible with
*    explicit agreement of the copyright holders.
*
* 3. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*
* 4. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* 5. Neither the name of the copyright holder nor the names of its
*    contributors may be used to endorse or promote products derived from this
*    software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#if defined(LT_EJECTA_IN_HOTPHASE) && defined(GM_MUPPI)
#error LT_EJECTA_IN_HOTPHASE AND GM_MUPPI cannot be active simultaneously!
#endif



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"
#include "../../Hydro/kernel.h"

#include "lt_error_codes.h"

#ifdef LT_STELLAREVOLUTION

#include "lt_sn.h"

#if defined(GL_CR_DUST) && (!defined(GL_DWEK) && !defined(GL_GIVENGRAIN))
#error with GL_CR_DUST either GL_DWEK and GL_GIVENGRAIN have to be defined
#endif

#if defined(GL_DWEK) && defined(GL_GIVENGRAIN)
#error only one of GL_DWEK and GL_GIVENGRAIN can be defined
#endif

#if defined(GL_CR_ALMOSTGIVENGRAIN) && !defined(GL_GIVENGRAIN)
#error nonsense to use GL_CR_ALMOSTGIVENGRAIN without GL_GIVENGRAIN
#endif



static int ndone;

long double local_produced_metals[LT_NMetP];

//!
/*!

  \brief This routine evolves all active stars and star-forming gas particles
  \param mode instructs evolve_SN whether to evolve all the particles or the active only

  \return the number of evolved particles

  This routine evolves all active stars and star-forming gas
  particles in order to advance them on the chemical timeline and
  to spread the released metals and energy released on their neighbourhood

*/

unsigned int evolve_SN(void)
{
  int starsnum, gasnum;
  int i, j, k, i_list;		/*!< general counters */
  int ndone_flag;		/*!< to account for the bunch size */
  int nimport, nexport;
  int recvTask, ngrp;		/*!< used in the communication structure */
  double Metals[LT_NMet], energy, egystep, LMMass;	/*!< stores the result from stellarevolution() */
  double tstart, tend;		/*!< used to account for used cpu-time */

#ifdef GL_CR_DUST
  double dDustL[LT_NMet];	// Timestep increments of dust, large grains
#endif /* GL_CR_DUST */

  /* : ---------------------------------------------- */
  /* : preliminary operations                         */
  /* : ---------------------------------------------- */


  long double metals_before[LT_NMetP];
  get_metals_checksum(0, metals_before);
  for(int j = 0; j < LT_NMetP; j++)
    local_produced_metals[j] = 0;


  /* if this option is set, the release of metals is
     halted below
     All.Below_this_redshift_stop_cooling */
#ifdef LT_STOP_MET_BELOW_Z
  if((1.0 / All.Time - 1.0) < All.Below_this_redshift_stop_cooling)
    return 0;
#endif

  /* find the total number of active stars and sf gas
     particles */
  tstart = second();
  count_evolving_stars(&starsnum, &gasnum);
  tend = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    printf("EXTRA TIMER: SN count_evolving = %g\n", timediff(tstart, tend));
#endif

  sumup_large_ints(1, &starsnum, &tot_starsnum);
  sumup_large_ints(1, &N_stars, &tot_allstarsnum);
  sumup_large_ints(1, &gasnum, &tot_gasnum);


/*   MPI_Allreduce(&starsnum, &tot_starsnum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD); */
/*   MPI_Allreduce(&gasnum, &tot_gasnum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD); */


  if(tot_starsnum + tot_gasnum == 0)	/* if no particles are active, exit */
    {
      myfree(NextChemActiveParticle);
      return 0;
    }


  if(ThisTask == 0)
    printf("calculating stellar evolution for %llu stars (over %llu) and %llu gas partcles\n", tot_starsnum,
	   tot_allstarsnum, tot_gasnum);
  fflush(stdout);

  for(i = 0; i < INUM; i++)
    infos[i] = 0;

  /* : ----------------------------------------------- */
  /* : allocate memory to exchange data                */
  /* : ----------------------------------------------- */

  /* allocate the space for the list of neighbours */
  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
  /* find available room in memory                 */
  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) /
	   (sizeof(struct data_index) + sizeof(struct data_nodelist) + 2 * sizeof(struct metaldata_index)));

  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  MetalDataIn =
    (struct metaldata_index *) mymalloc("MetalDataIn", All.BunchSize * sizeof(struct metaldata_index));

  /* now, all the evolving stars have the right spreading lenght and
   *  have stored the total weight of neighbours.
   *
   * the nex step is to communicate the stellar evolution details and
   *  the total weight in order to spread over the neighbours
   */

  i_list = starsnum + gasnum - 1;

  /* MAIN (outer) cycle                 */
  /* .................................. */

  //  i = FirstChemActiveParticle;

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Recv_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0; i_list >= 0; i_list--)
	{
	  i = NextChemActiveParticle[i_list];
	  if(P[i].Type & EVOLVE)
	    {
	      P[i].Type &= ~EVOLVE;
	      /* here we account for the stellar evolution and do
	         operations for local particles
	       */
#ifdef GL_CR_DUST
	      for(k = 0; k < LT_NMet; k++)
		dDustL[k] = 0.;	//needed to overcome non-assignment when under minimum metallicity for forming dust

	      /* XXXX debug outpt
	         for(k=0; k<LT_NMet; k++)
	         printf(" DUST L before computation, TAsk %d, part %d, metal %d, value %e\n",
	         ThisTask, i, k, SphP[i].DustL[k]); fflush(stdout);
	       */
	      if((LMMass =
		  perform_stellarevolution_operations(i, &nexport, &Metals[0], &dDustL[0], &energy,
						      &egystep)) < 0)
#else
	      if((LMMass =
		  perform_stellarevolution_operations(i, &nexport, &Metals[0], &energy, &egystep)) < 0)
#endif
		{
		  P[i].Type |= EVOLVE;
		  break;
		}
	      else
		{
		  if(P[i].Type == 4)
		    MetP[P[i].pt.MetID].weight = 0;
		}
	    }
	}
      tend = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	printf("EXTRA TIMER: SN evolve local = %g\n", timediff(tstart, tend));
#endif

      tstart = second();
      MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

      STD_SORT(MetalDataIn, nexport, sizeof(struct metaldata_index), metaldata_index_compare);

      tend = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	printf("EXTRA TIMER: SN sort = %g\n", timediff(tstart, tend));
#endif

      tstart = second();
#ifndef IMPORT_ALLtoALLv_FROM_G4
      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);
#else
#ifdef IMPORT_ALLtoALLv_FROM_G4_WHEN_TIMEBIN_ACTIVE_AFTER
      int non_active_timebin = 0;
      for(j = TIMEBINS - 1; j >= 0; j--)
	{
	  if(TimeBinActive[j])
	    {
	      break;
	    }
	  else
	    {
	      non_active_timebin++;
	    }
	}
      if(non_active_timebin > IMPORT_ALLtoALLv_FROM_G4_WHEN_TIMEBIN_ACTIVE_AFTER)
	{
#endif
	  MPI_Win_create(Recv_count, NTask * sizeof(MPI_INT), sizeof(MPI_INT), MPI_INFO_NULL,
			 MYMPI_COMM_WORLD, &win);
	  MPI_Win_fence(0, win);
	  for(i = 0; i < NTask; ++i)
	    {
	      if(Send_count[i] > 0)
		MPI_Put(&Send_count[i], 1, MPI_INT, i, ThisTask, 1, MPI_INT, win);
	    }
	  MPI_Win_fence(0, win);
	  MPI_Win_free(&win);
#ifdef IMPORT_ALLtoALLv_FROM_G4_WHEN_TIMEBIN_ACTIVE_AFTER
	}
      else
	{
	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);
	}
#endif //ACTIVE_AFTER
#endif
      tend = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	printf("EXTRA TIMER: SN wait1 = %g\n", timediff(tstart, tend));
#endif
      //timewait1 += timediff(tstart, tend);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      /* preparing particle data for export has been done during the local evaluation */

      /* exchange particle data */

      MetalDataGet =
	(struct metaldata_index *) mymalloc("MetalDataGet", nimport * sizeof(struct metaldata_index));

      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  //      sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&MetalDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct metaldata_index), MPI_BYTE,
			       recvTask, TAG_SMOOTH_A,
			       &MetalDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct metaldata_index), MPI_BYTE,
			       recvTask, TAG_SMOOTH_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	printf("EXTRA TIMER: SN comm1 = %g\n", timediff(tstart, tend));
#endif
      //timecommsumm1 += timediff(tstart, tend);

      /* now do the particles that were sent to us */

      tstart = second();
      for(j = 0; j < nimport; j++)
#ifdef GL_CR_DUST
	spread_evaluate(j, BUFFER, MetalDataGet[j].Metals, MetalDataGet[j].Dust,
			MetalDataGet[j].numsnIa, MetalDataGet[j].numsnII, MetalDataGet[j].LMMass,
			MetalDataGet[j].energy, MetalDataGet[j].egystep, &ngrp, &ngrp);
#else
	spread_evaluate(j, BUFFER, MetalDataGet[j].Metals, MetalDataGet[j].LMMass, MetalDataGet[j].energy,
			MetalDataGet[j].egystep, &ngrp, &ngrp);
#endif

      tend = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	printf("EXTRA TIMER: SN evolve external = %g\n", timediff(tstart, tend));
#endif
      /*timecomp2 += timediff(tstart, tend); */

      if(i_list < 0)		/* flag the actual active particle */
	ndone_flag = 1;
      else
	ndone_flag = 0;

      tstart = second();
      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
      tend = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	printf("EXTRA TIMER: SN wait2 = %g\n", timediff(tstart, tend));
#endif
      /* timewait2 += timediff(tstart, tend); */

      myfree(MetalDataGet);

    }
  while(ndone < NTask);

  myfree(MetalDataIn);
  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  tstart = second();
#ifdef _OPENMP
#pragma omp parallel for
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0 && SphP[i].MassRes > 0)
#ifdef GM_MUPPI
      {
	if(P[i].Type == 0)
	  if(is_chemically_active(i))
	    SphP[i].mstar = 0.0;

	if(P[i].Type == 0 && SphP[i].MassRes > 0)
#else
      if(P[i].Type == 0 && SphP[i].MassRes > 0)
#endif
	{
#ifdef GM_MUPPI
	  if(SphP[i].MultiPhase == 0)	/* otherwise, this is accounted for in MUPPI_main */
	    {
#endif
#if GADGET_HYDRO == HYDRO_MFM
	      mfm_add_mass(i,SphP[i].MassRes);
#endif
	      P[i].Mass += SphP[i].MassRes;
	      SphP[i].MassRes = 0;
#ifdef GM_MUPPI
	    }
#endif

	}
#ifdef GM_MUPPI
}
#endif


tend = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
if(ThisTask == 0)
  printf("EXTRA TIMER: SN final = %g\n", timediff(tstart, tend));
#endif

myfree(NextChemActiveParticle);

search_for_metalspread = 0;

long double metals_after[LT_NMetP];
long double tot_metals[LT_NMetP];
get_metals_checksum(0, metals_after);
get_metals_checksum(1, tot_metals);

long double all_produced_metals[NTask * LT_NMetP];
MPI_Gather(local_produced_metals, LT_NMetP * sizeof(long double), MPI_BYTE, all_produced_metals,
	   LT_NMetP * sizeof(long double), MPI_BYTE, 0, MYMPI_COMM_WORLD);


if(ThisTask == 0)
  {
    int failures_step = 0;
    int failures_global = 0;
    for(int j = 0; j < LT_NMetP; j++)
      {
	for(int i = 1; i < NTask; i++)
	  all_produced_metals[j] += all_produced_metals[i * LT_NMetP + j];
	All.MetalsCheckSum[j] += all_produced_metals[j];

	/* if( j == 0 ) */
	/*   continue; */

	metals_after[j] -= metals_before[j];

	long double diff_step =
	  fabsl(all_produced_metals[j] >
		0 ? (metals_after[j] -
		     all_produced_metals[j]) / all_produced_metals[j] : (long double) (metals_after[j] > 0));

	long double diff_global =
	  fabsl(All.MetalsCheckSum[j] >
		0 ? (tot_metals[j] -
		     All.MetalsCheckSum[j]) / All.MetalsCheckSum[j] : (long double) (metals_after[j] > 0));

	failures_step += (diff_step > 0.0001);
	failures_global += (diff_global > 0.0001);
      }

    if(failures_step)
      {
	printf("=== CHECKSUM step %d time %g : ", All.NumCurrentTiStep, All.Time);
	for(int j = 0; j < LT_NMetP; j++)
	  printf("[ %Lg %Lg %Lg ] ", metals_after[j], all_produced_metals[j],
		 (all_produced_metals[j] >
		  0 ? (metals_after[j] - all_produced_metals[j]) / all_produced_metals[j] : 0.0));
	printf("\n");
      }
    if(failures_global)
      {
	printf("=== CHECKSUM global %d time %g : ", All.NumCurrentTiStep, All.Time);
	for(int j = 0; j < LT_NMetP; j++)
	  printf("[ %Lg %Lg %Lg ] ", tot_metals[j], All.MetalsCheckSum[j],
		 (All.MetalsCheckSum[j] >
		  0 ? (tot_metals[j] - All.MetalsCheckSum[j]) / All.MetalsCheckSum[j] : 0.0));
	printf("\n");
      }

  }

MPI_Barrier(MYMPI_COMM_WORLD);



return tot_starsnum;
}

/*
   =======================================================
     END of   M a i n   R o u t i n e
     .........................................
   =======================================================
*/


/* :::.............................................................. */
/* ***************************************************************** */


/* NOTE: not used now */
int iterate(int n, int ntrue, MyFloat grad, double *deltax)
{
  MyFloat old, Lold, L, L3, rho;

  int i = 0;

  Lold = *deltax;		/* at the begin *deltax contains the guess for L */
  rho = n / (*deltax * *deltax * *deltax);
  L = *deltax * pow(ntrue / n, 0.33333333);
  L3 = L * L * L;
  *deltax = L - *deltax;
  old = *deltax / 2;

  while(fabs(*deltax - old) / (*deltax) > 0.001)
    {
      i++;
      old = *deltax;
      L3 = L * L * L;
      *deltax = (ntrue - rho * L3) / (grad * L3);
      L = Lold + *deltax;
    }
  return i;
}


/* / ---------------------------------------------- \
 * |                                                |
 * | SECTION IIa                                    |
 * |                                                |
 * | ........................                       |
 * |                                                |
 * | Neighbours search                              |
 * | Stellar Evolution Driver                       |
 * | Spreading                                      |
 * |                                                |
 * |    > stellarevolution                          |
 * |    > spread                                    |
 * |                                                |
 * \ ---------------------------------------------- / */



/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */


int metaldata_index_compare(const void *a, const void *b)
{
  if(((struct metaldata_index *) a)->Task < (((struct metaldata_index *) b)->Task))
    return -1;

  if(((struct metaldata_index *) a)->Task > (((struct metaldata_index *) b)->Task))
    return +1;

  if(((struct metaldata_index *) a)->Index < (((struct metaldata_index *) b)->Index))
    return -1;

  if(((struct metaldata_index *) a)->Index > (((struct metaldata_index *) b)->Index))
    return +1;

  if(((struct metaldata_index *) a)->IndexGet < (((struct metaldata_index *) b)->IndexGet))
    return -1;

  if(((struct metaldata_index *) a)->IndexGet > (((struct metaldata_index *) b)->IndexGet))
    return +1;

  return 0;
}


/*
   * ..................... *
   :                       :
   :   Evolving Stars      :
   * ..................... *
*/

#ifdef GL_CR_DUST
double perform_stellarevolution_operations(int i, int *nexport, double *Metals, double *dDustL,
					   double *energy, double *egystep)
#else
double perform_stellarevolution_operations(int i, int *nexport, double *Metals, double *energy,
					   double *egystep)
#endif
{
  int j;
  int Yset, YZbin;
  double Z;
  double tstart, tend;
  double LMMass;
  float metals[LT_NMet];
  int ret;

#ifdef GL_CR_DUST
  float dust[LT_NMet];
  double numsnIa, numsnII;
#endif /* GL_CR_DUST */

#ifdef LT_TRACK_CONTRIBUTES
  NULL_CONTRIB(&contrib);
  NULL_EXPCONTRIB(&IIcontrib);
  NULL_EXPCONTRIB(&Iacontrib);
  NULL_EXPCONTRIB(&AGBcontrib);
#endif

  for(j = 0; j < LT_NMet; j++)
    Metals[j] = 0;

  *energy = 0;
  *egystep = 0;
  LMMass = 0;

  int mySFi, myIMFi;
  IMF_Type *IMFp;
  SF_Type *SFp;

  get_SF_index(i, &mySFi, &myIMFi);
  SFp = (SF_Type *) & SFs[mySFi];
  IMFp = (IMF_Type *) & IMFs[myIMFi];

  tstart = second();
#ifndef LT_LOCAL_IRA
  if(P[i].Type & 4)
#endif
    {
      //printf("calling stellarevolution for %llu\n", (unsigned long long)P[i].ID); /* if this is an active star, calculate its evolution and update fields */
#ifdef GL_CR_DUST
      LMMass = stellarevolution(i, &Metals[0], energy, egystep, &dDustL[0], &numsnIa, &numsnII);
#else
      LMMass = stellarevolution(i, &Metals[0], energy, egystep);
#endif
      for(int j = 0; j < LT_NMetP; j++)
	local_produced_metals[j] += (j != Hyd ? Metals[j] : 0);

    }

#ifndef LT_LOCAL_IRA
  else if(SFp->nonZeroIRA > 0)
    {				/* if this is a star-forming gas particle accountfor the IRA part, if any */

      Yset = IMFp->YSet;
      Z = get_metallicity(i, -1);
      for(YZbin = IIZbins_dim[Yset] - 1; Z < IIZbins[Yset][YZbin] && YZbin > 0; YZbin--)
	;
      for(j = 0; j < LT_NMet; j++)
	{
	  Metals[j] = SphP[i].mstar * SnII_ShortLiv_Yields[Yset][j][YZbin];
	}

#ifndef LT_SNegy_IN_HOTPHASE
      *energy = 0.;		/* already used in effective model */
      *egystep = 0.;
#endif
#if defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
      LMMass = 0;
#endif

#ifdef LT_TRACK_CONTRIBUTES
      for(j = 0; j < LT_NMetP; j++)
	{
	  IIcontrib[j] = 1;
	  Iacontrib[j] = AGBcontrib[j] = 0;
	}
#endif
    }
#endif

#ifdef LT_TRACK_CONTRIBUTES
  pack_contrib(&contrib, myIMFi, IIcontrib, Iacontrib, AGBcontrib);
#endif

  tend = second();
  infos[SN_Calc] += timediff(tstart, tend);

  tstart = second();

  for(j = 0; j < LT_NMet; j++)
    {
#ifdef GL_CR_DUST
      metals[j] = (float) Metals[j];
      dust[j] = (float) dDustL[j];
#else
      metals[j] = (float) Metals[j];
#endif
    }
#ifdef GL_CR_DUST
  ret =
    spread_evaluate(i, LOCAL, &metals[0], &dust[0], numsnIa, numsnII, LMMass, *energy, *egystep, nexport,
		    Send_count);
#else
  ret = spread_evaluate(i, LOCAL, &metals[0], LMMass, *energy, *egystep, nexport, Send_count);
#endif

  tend = second();
  infos[SN_Spread] += timediff(tstart, tend);

  if(ret >= 0)
    return LMMass;
  else
    return -1;
}


/* This routine return the stellar evolution products for a given star
 * it also upgrade the mass of the star and the last chemical times to the
 * current values.
 */
#ifdef GL_CR_DUST
double stellarevolution(int i, double *metals, double *energy, double *egystep, double *dDustL,
			double *NumsnIa, double *NumsnII)
#else
double stellarevolution(int i, double *metals, double *energy, double *egystep)
#endif
{
  int I, Yset;

/* #if !(defined (UM_CHEMISTRY) && defined (UM_METAL_COOLING)) */
  double mymetals[LT_NMet], sum_mymetals[LT_NMet];
/* #else */
/*   double mymetals[LT_NMet+1], sum_mymetals[LT_NMet+1]; */
/* #endif */
  double myenergy;
  double mass;
  double starlifetime, mylifetime, delta_evolve_time, prec_evolve_time;
  double inf_mass, sup_mass;
  double numsn;
  int chem_step;
  double NextChemTime, LMmass = 0;
  int j, ti_min;

#if defined(GL_CR_DUST)
  double SnIImet[LT_NMet], SnIamet[LT_NMet], AGBmet[LT_NMet];
  double wk;
#endif


  I = P[i].pt.MetID;

  int mySFi, myIMFi;

  /* find the IMF associated to this particle */
  get_SF_index(i, &mySFi, &myIMFi);

  /* find the yield set associated to this IMF */
  Yset = IMFs[myIMFi].YSet;

  for(j = 0; j < LT_NMet; j++)
    mymetals[j] = sum_mymetals[j] = 0;
  *energy = 0.;
#if defined(GL_CR_DUST)
  for(j = 0; j < LT_NMet; j++)
    SnIImet[j] = SnIamet[j] = AGBmet[j] = 0.;
#endif
  *egystep = 0.;

  myenergy = 0.;
  /* calculate the lifetime of the star */
  if(MPP(i).StellarAge < (float) All.Time)
    starlifetime = get_age(MPP(i).StellarAge) - All.Time_Age;
  else
    /* due to float-double conversion possible round off */
    starlifetime = 0;

  if(!All.ComovingIntegrationOn)
    starlifetime *= -1.0;


  if(UseSnII && (MetP[I].LastChemTime < All.mean_lifetime) && (Yset < All.II_Nset_ofYields))
    {

      mylifetime = starlifetime;

      if((mylifetime >= SFs[mySFi].ShortLiv_TimeTh) &&
	 (mylifetime - MetP[I].LastChemTime >= All.MinChemTimeStep))
	{
	  prec_evolve_time = MetP[I].LastChemTime;
	  if(mylifetime > All.mean_lifetime)
	    mylifetime = All.mean_lifetime;

	  inf_mass = dying_mass(mylifetime);
	  sup_mass = dying_mass(prec_evolve_time);

	  mass = get_SnII_product(i, Yset, &mymetals[0], &myenergy, inf_mass, sup_mass, &numsn);

	  if(mass >= 0)
	    {
	      P[i].Mass -= mass;

	      for(j = 0; j < LT_NMet; j++)
		sum_mymetals[j] += mymetals[j];

	      if(myenergy > 0)
		{
		  *egystep = (*energy * *egystep + myenergy * (mylifetime - prec_evolve_time)) /
		    (*energy + myenergy);
		  *energy += myenergy;
		}
	      numsn = numsn * UnitMassFact / All.HubbleParam;
#if defined(GL_CR_DUST)
	      for(j = 0; j < LT_NMet; j++)
		SnIImet[j] = mymetals[j];
	      *NumsnII = numsn;
#endif //GL_CR_DUST

	    }

	  MetP[I].LastChemTime = mylifetime;
	}

#ifdef LT_TRACK_CONTRIBUTES
      for(j = 0; j < LT_NMetP; j++)
	IIcontrib[j] = mymetals[j];
#endif

    }



  if(UseSnIa || UseAGB)
    {
      if(starlifetime > All.mean_lifetime)
	{
	  prec_evolve_time = MetP[I].LastChemTime;
	  delta_evolve_time = starlifetime - prec_evolve_time;

	  if(delta_evolve_time >= All.MinChemTimeStep)
	    {
	      inf_mass = dying_mass(prec_evolve_time + delta_evolve_time);
	      sup_mass = dying_mass(prec_evolve_time);

	      if(UseSnIa && (Yset < All.Ia_Nset_ofYields))
		{
		  mass =
		    get_SnIa_product(i, Yset, &mymetals[0], &myenergy, prec_evolve_time, delta_evolve_time);

		  if(mass >= 0)
		    {
		      P[i].Mass -= mass;

		      if(myenergy > 0)
			{
			  *egystep = (*energy * *egystep + myenergy * delta_evolve_time) /
			    (*energy + myenergy);
			  *energy += myenergy;
			}

		      for(j = 0; j < LT_NMet; j++)
			{
			  sum_mymetals[j] += mymetals[j];
			}

		    }

#if defined(GL_CR_DUST)
		  for(j = 0; j < LT_NMet; j++)
		    SnIamet[j] = mymetals[j];
		  *NumsnIa = myenergy / All.SnIaEgy;	// GLG 081016
		  if(*NumsnIa < 1e-30)
		    *NumsnIa = 0.0;	// GLG 081016
#endif //GL_CR_DUST


#ifdef LT_TRACK_CONTRIBUTES
		  for(j = 0; j < LT_NMetP; j++)
		    Iacontrib[j] = (float) mymetals[j];
#endif


		}

	      if(UseAGB && (Yset < All.AGB_Nset_ofYields))
		{
		  mass = get_AGB_product(i, Yset, &mymetals[0], inf_mass, sup_mass, &numsn);

		  if(mass >= 0)
		    {
		      P[i].Mass -= mass;


		      for(j = 0; j < LT_NMet; j++)
			{
			  sum_mymetals[j] += mymetals[j];
			  /* keep track of ejecta that are not injected explosively */
			  LMmass += mymetals[j];
			}
#if defined(GL_CR_DUST)
		      for(j = 0; j < LT_NMet; j++)
			AGBmet[j] = mymetals[j];
#endif


		    }

#ifdef LT_TRACK_CONTRIBUTES
		  for(j = 0; j < LT_NMetP; j++)
		    AGBcontrib[j] = (float) mymetals[j];
#endif

		  MetP[I].LastChemTime = starlifetime;
		}
	    }

	}
    }


#ifdef GL_CR_DUST
  /* coefficients in lt_sn.h */

// GL 020217 begin : calculate condensation efficiencies for C and O following Dwek 98
#ifdef GL_DWEK			// GL 030217

// first AGB channel
//  GLG 050517 in this Dwek case for simplicity dcondeffAGB_var is actually the mass going to dust in AGB channels,
// not simply the condensation efficiency
  dcondeffAGB_var[ind_O] = 0.0;
  if(AGBmet[ind_C] / AGBmet[ind_O] <= 0.75)	// Dwek eq 23 silicate grains are produced
    //GLG 050517 0.75 not 1 because quantities are masses, not numbers
    {
      for(j = 0; j < LT_NMet; j++)
	{
	  dcondeffAGB_var[j] = dcondeffAGB[j] * AGBmet[j];	//GLG 050517
	}
      dcondeffAGB_var[ind_C] = 0.0;
      dcondeffAGB_var[ind_O] = 0.0;
      for(j = 0; j < n_sil; j++)
	{
	  dcondeffAGB_var[ind_O] += dcondeffAGB[ind_sil[j]] * AGBmet[ind_sil[j]] / mu_sil[j];
	}
      dcondeffAGB_var[ind_O] *= mu_O;
      dcondeffAGB_var[ind_C] = 0.0;
    }
  else				// dwek eq 22 carbon grains are produced
    {
      for(j = 0; j < LT_NMet; j++)
	dcondeffAGB_var[j] = 0.0;
      dcondeffAGB_var[ind_C] = dcondeffAGB[ind_C] * (AGBmet[ind_C] - 0.75 * AGBmet[ind_O]);	// 0.75 is the ratio mu_C/mu_O
    }

  // then SNae Dwek eq 24 and
  for(j = 0; j < LT_NMet; j++)	//GLG 060517
    {
      dcondeffIa_var[j] = dcondeffIa[j] * SnIamet[j];
      dcondeffII_var[j] = dcondeffII[j] * SnIImet[j];
    }


  dcondeffIa_var[ind_O] = 0.0;
  dcondeffII_var[ind_O] = 0.0;
  for(j = 0; j < n_sil; j++)
    {
      //GLG 060517 dcondeffI?_var[ind_sil[j]] already multiplied abobe *SnI?met[ind_sil[j]]
      dcondeffIa_var[ind_O] += dcondeffIa_var[ind_sil[j]] / mu_sil[j];
      dcondeffII_var[ind_O] += dcondeffII_var[ind_sil[j]] / mu_sil[j];
    }
  dcondeffIa_var[ind_O] *= mu_O;
  dcondeffII_var[ind_O] *= mu_O;


// GL 020217 end



// GLG 20111q7 fix (porfin!) the problem that this scheme gives problems with SNIa, producing mostly Fe
  if(dcondeffIa_var[ind_O] > SnIamet[ind_O])
    {
      float cor;
      cor = SnIamet[ind_O] / dcondeffIa_var[ind_O];
      dcondeffIa_var[ind_O] = SnIamet[ind_O];
      for(j = 0; j < n_sil; j++)
	dcondeffIa_var[ind_sil[j]] *= cor;
    }


  for(j = 0; j < LT_NMet; j++)
    {
// GL 020217 following uses dcondeffAGB_var instead of dcondeffAGB
      dDustL[j] = dcondeffAGB_var[j] + dcondeffIa_var[j] + dcondeffII_var[j];
      // GLG 060517  dcondeff???_var are alredy masses
      sum_mymetals[j] -= dDustL[j];	//GLG 030217 prima calcolato di nuovo...suppongo totalmente inutile
      /* XXXX WARNING: non abbiamo sottratto individualmente le polveri dai vari canali, il tracking
         in questa versione NON funzionera' correttamente (o cmq terra' conto della produzione globale comprese le polveri) */
      if(sum_mymetals[j] < 0.)
	{
	  printf
	    ("   DUST WARNING: Task %d metal %d is negative (%e), contrib met %f %f %f contrib dust %e %e %e\n",
	     ThisTask, j, sum_mymetals[j], dcondeffAGB_var[j], dcondeffIa_var[j], dcondeffII_var[j],
	     AGBmet[j], SnIamet[j], SnIImet[j]);
	  fflush(stdout);
	  sum_mymetals[j] = 0.0;

	}
    }


#endif /* GL_DWEK */


// GL 030217 begin


#ifdef GL_GIVENGRAIN		// GL 030217

  // find the number of grain that can be formed in each channel, set by the minimum
  // between the elements of the number of assumed chemical compound "molecule" that can be formed


  nummol_AGB = effsil_AGB * (AGBmet[ind_O] / mu_O - AGBmet[ind_C] / mu_C) / n_ato_O;
  if(nummol_AGB < 0)
    nummol_AGB = 0;
  nummol_II = effsil_II * SnIImet[ind_O] / mu_O / n_ato_O;
  nummol_Ia = effsil_Ia * SnIamet[ind_O] / mu_O / n_ato_O;

  for(j = 0; j < n_sil; j++)
    {


#ifdef GL_CR_ALMOSTGIVENGRAIN	// GL 29072019

      if(ind_sil[j] == ind_Mg || ind_sil[j] == ind_Fe)
	{

	  wk =
	    effsil_AGB * (AGBmet[ind_sil[ind_Mg]] / mu_sil[ind_Mg] +
			  AGBmet[ind_sil[ind_Fe]] / mu_sil[ind_Fe]) / (n_ato_sil[ind_Mg] + n_ato_sil[ind_Fe]);
	  if(wk < nummol_AGB)
	    {
	      nummol_AGB = wk;
	    }

	  wk =
	    effsil_II * (SnIImet[ind_sil[ind_Mg]] / mu_sil[ind_Mg] +
			 SnIImet[ind_sil[ind_Fe]] / mu_sil[ind_Fe]) / (n_ato_sil[ind_Mg] + n_ato_sil[ind_Fe]);
	  if(wk < nummol_II)
	    {
	      nummol_II = wk;
	    }

	  wk =
	    effsil_Ia * (SnIamet[ind_sil[ind_Mg]] / mu_sil[ind_Mg] +
			 SnIamet[ind_sil[ind_Fe]] / mu_sil[ind_Fe]) / (n_ato_sil[ind_Mg] + n_ato_sil[ind_Fe]);
	  if(wk < nummol_Ia)
	    {
	      nummol_Ia = wk;
	    }

	}
      else
	{
	  wk = effsil_AGB * AGBmet[ind_sil[j]] / mu_sil[j] / n_ato_sil[j];
	  if(wk < nummol_AGB)
	    {
	      nummol_AGB = wk;
	    }

	  wk = effsil_II * SnIImet[ind_sil[j]] / mu_sil[j] / n_ato_sil[j];
	  if(wk < nummol_II)
	    {
	      nummol_II = wk;
	    }

	  wk = effsil_Ia * SnIamet[ind_sil[j]] / mu_sil[j] / n_ato_sil[j];
	  if(wk < nummol_Ia)
	    {
	      nummol_Ia = wk;
	    }
	}
#else /* GL_CR_ALMOSTGIVENGRAIN */
      wk = effsil_AGB * AGBmet[ind_sil[j]] / mu_sil[j] / n_ato_sil[j];
      if(wk < nummol_AGB)
	{
	  nummol_AGB = wk;
	}

      wk = effsil_II * SnIImet[ind_sil[j]] / mu_sil[j] / n_ato_sil[j];
      if(wk < nummol_II)
	{
	  nummol_II = wk;
	}

      wk = effsil_Ia * SnIamet[ind_sil[j]] / mu_sil[j] / n_ato_sil[j];
      if(wk < nummol_Ia)
	{
	  nummol_Ia = wk;
	}

#endif /* GL_CR_ALMOSTGIVENGRAIN */


    }

// GL 030217 note-catarsis: in fortran 90 or later the previous search would be one line each X:
// nummol_X=MINVAL(AGBmet(ind_sil)/mu(:)/n_ato_sil(:))
// but we are smart people and then we use c for numerical programmimg:-)


// compute dcondeffAGB_var only for C. In this GIVENGRAIN case only C is treated with dcondeff*
// first AGB channel
  dcondeffAGB_var[ind_O] = 0.0;
  if(AGBmet[ind_C] / AGBmet[ind_O] <= 0.75)	// Dwek eq 23 silicate grains are produced
    //GLG 050517 0.75 (the ratio mu_C/mu_O) and not 1 because quantities are masses, not numbers
    {
      dcondeffAGB_var[ind_C] = 0.0;
    }
  else				// dwek eq 22 carbon grains are produced
    {
      dcondeffAGB_var[ind_C] = dcondeffAGB[ind_C];
    }


  for(j = 0; j < LT_NMet; j++)
    dDustL[j] = 0.;

// GL 020217 following uses dcondeffAGB_var instead of dcondeffAGB
// In this GIVENGRAIN case only C is treated with dcondeff*0
  dDustL[ind_C] =
    dcondeffAGB_var[ind_C] * (AGBmet[ind_C] - 0.75 * AGBmet[ind_O]) + dcondeffIa[ind_C] * SnIamet[ind_C] +
    dcondeffII[ind_C] * SnIImet[ind_C];
// 0.75 is the ratio mu_C/mu_O
  sum_mymetals[ind_C] -= dDustL[ind_C];
/* XXXX WARNING: non abbiamo sottratto individualmente le polveri dai vari canali, il tracking
in questa versione NON funzionera' correttamente (o cmq terra' conto della produzione globale comprese le polveri) */

  for(j = 0; j < n_sil; j++)
    {
      dDustL[ind_sil[j]] = (nummol_AGB + nummol_Ia + nummol_II) * mu_sil[j] * n_ato_sil[j];
      sum_mymetals[ind_sil[j]] -= dDustL[ind_sil[j]];
    }
  dDustL[ind_O] = (nummol_AGB + nummol_Ia + nummol_II) * mu_O * n_ato_O;	// O not included in previous loop
  sum_mymetals[ind_O] -= dDustL[ind_O];

//XXXXX questo check forse ancora dentro un ciclo completo?
  for(j = 0; j < LT_NMet; j++)
    {
      if(sum_mymetals[j] < 0.)
	{
	  printf("   DUST WARNING: Task %d metal %d is negative (%e), coeff %f %f %f contrib %e %e %e\n",
		 ThisTask, j, sum_mymetals[j], dcondeffAGB_var[j], dcondeffIa[j], dcondeffII[j],
		 AGBmet[j], SnIamet[j], SnIImet[j]);
	  fflush(stdout);
	  sum_mymetals[j] = 0.0;

	}

    }


#endif /* GL_GIVENGRAIN */


// GL 030217 end
#endif /* GL_CR_DUST */


  NextChemTime = get_NextChemTime(starlifetime, mySFi, 0x0);
  j = get_chemstep_bin(All.Time, All.Time_Age - NextChemTime, &chem_step, i);

/*   double Time_Age_new, Time_new; */
/*   Time_new = All.TimeBegin * exp( (All.Ti_Current + (1 << j)) * All.Timebase_interval); */
/*   Time_Age_new = get_age(Time_new); */
/*   printf(">>> %llu @ %d @ %g has chemtimestep=%g, j=%d => new a: %u %u %g %g dt: %g frac: %g\n",  */
/* 	 (unsigned long long int)P[i].ID, All.NumCurrentTiStep, starlifetime,All.Time_Age - NextChemTime, j,  */
/* 	 All.Ti_Current, 1<<j, All.Time, Time_new, */
/* 	 All.Time_Age - Time_Age_new,  */
/* 	 (NextChemTime - Time_Age_new)/(All.Time_Age - NextChemTime)); */
  /*       if(TimeBinActive[j] == 0) */
  /*         { */
  /*           while(TimeBinActive[j] == 0 && j > 0) */
  /*             j--; */
  /*           chem_step = j ? (1 << j) : 0; */
  /*         } */
  if(All.Ti_Current >= TIMEBASE)
    chem_step = j = 0;

  if((TIMEBASE - All.Ti_Current) < chem_step)
    {
      chem_step = TIMEBASE - All.Ti_Current;
      ti_min = TIMEBASE;
      while(ti_min > chem_step)
	ti_min >>= 1;
      chem_step = ti_min;
      j = get_timestep_bin(chem_step);
    }

  if(j != MetP[I].ChemTimeBin)
    {
      TimeBinCountStars[MetP[I].ChemTimeBin]--;
      TimeBinCountStars[j]++;
      MetP[I].ChemTimeBin = j;
    }


#ifdef LT_TRACK_CONTRIBUTES
  for(j = 0; j < LT_NMetP; j++)
    if(sum_mymetals[j] > 0)
      {
	IIcontrib[j] = (float) ((double) IIcontrib[j] / sum_mymetals[j]);
	AGBcontrib[j] = (float) ((double) AGBcontrib[j] / sum_mymetals[j]);
	Iacontrib[j] = (float) ((double) Iacontrib[j] / sum_mymetals[j]);
      }
    else
      {
	IIcontrib[j] = 0;
	AGBcontrib[j] = 0;
	Iacontrib[j] = 0;
      }
#endif

  if(P[i].Mass <= 0)
    {
      /*.. shouldn't occour (oh, really??) */
      printf("  @@@@@@@@@ %i %i %llu %i %g %g %g\n", ThisTask, i, (unsigned long long) P[i].ID,
	     All.NumCurrentTiStep, P[i].Mass, MPP(i).StellarAge, starlifetime);
      fflush(stdout);
    }

  for(j = 0; j < LT_NMet; j++)
    {
      metals[j] = (float) sum_mymetals[j];
    }

  return (float) LMmass;
}



/*
   * ..................... *
   :                       :
   :   Spreading           :
   * ..................... *
*/



#ifdef GL_CR_DUST
int spread_evaluate(int target, int mode, float *metals, float *dust, double numsnIa, double numsnII,
		    float LMmass, double energy, double egystep, int *nexport, int *nsend_local)
#else
int spread_evaluate(int target, int mode, float *metals, float LMmass, double energy, double egystep,
		    int *nexport, int *nsend_local)
#endif
{
  int I, startnode, listindex = 0, Type;
  int numngb_inbox;
  MyLongDouble Pos[3];
  MyFloat L;
  double weight;
  double L2;
  double add_mass;
  double egyfrac;
  double linv, linv3, linv4;

#ifdef LT_EJECTA_IN_HOTPHASE
  double a3inv;
  double dt, current_egy, current_egyhot, hotmass;

#ifdef LT_SNegy_IN_HOTPHASE
  double x_hotejecta, sn_spec_egy;
#endif
#if defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
  double v;
#endif
#endif


#ifdef LT_TRACK_CONTRIBUTES
  float contrib_metals[LT_NMetP];
#endif

#if defined(LT_SEv_INFO) && defined(LT_EJECTA_IN_HOTPHASE)
  double tot_spreadegy = 0, spreadegy_ratio, tot_agbspreadegy = 0;
  double agb_frac;
  double egy_ratio = 0, x_ratio = 0;
#endif

#if defined(LT_ZAGE) || defined(GM_MUPPI)
  double metal_add_mass;
#endif
#ifdef LT_ZAGE_LLV
  double llvmetal_add_mass;
#endif

  int mySFi;
  SF_Type *SFp;

/* #if defined (UM_CHEMISTRY) && defined (UM_METAL_COOLING) */
/*   float myFillEl_mu; */
/* #endif */

#ifdef LT_EJECTA_IN_HOTPHASE
  a3inv = 1.0 / (All.Time * All.Time * All.Time);
#endif

  startnode = All.MaxPart;

  add_mass = 0;
#ifdef LT_ZAGE
  metal_add_mass = 0;
#endif
#ifdef LT_ZAGE_LLV
  llvmetal_add_mass = metals[Iron];
#endif

  for(int k = 0; k < LT_NMet; k++)
    {
      add_mass += metals[k];
#ifdef LT_ZAGE
      if(k != Hyd && k != Hel)
	metal_add_mass += metals[k];
#endif
    }
  if(add_mass == 0)
    /* may happen if a star is evolved when too less
     * time has elapsed since the last evolution;
     * altghough not dangerous, this calls for a more
     * clever timing! */
    return (0);

  if(mode)
    {
      Type = MetalDataGet[target].Type;
      for(int j = 0; j < 3; j++)
	Pos[j] = MetalDataGet[target].Pos[j];
      L = MetalDataGet[target].L;
      weight = MetalDataGet[target].weight;
#ifdef LT_TRACK_CONTRIBUTES
      contrib = MetalDataGet[target].contrib;
#endif
      mySFi = MetalDataGet[target].SFi;
      SFp = (SF_Type *) & SFs[mySFi];
    }
  else
    {
      Type = P[target].Type & 4;
      L = P[target].Hsml;
      if(Type)
	{
	  I = P[target].pt.MetID;
	  weight = MetP[I].weight;
	}
      else
	weight = SphP[target].Density;

      for(int j = 0; j < 3; j++)
	Pos[j] = P[target].Pos[j];

      mySFi = (int) (SFp - SFs);
    }

#ifdef LT_EJECTA_IN_HOTPHASE
  agb_frac = LMmass / add_mass;
#endif

#ifdef LT_SNegy_IN_HOTPHASE
  sn_spec_egy = energy / add_mass;
  energy = 0;
  egystep = 0.;
#endif

  kernel_hinv(L, &linv, &linv3, &linv4);
  L2 = L * L;


  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = MetalDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifndef KD_FRICTION
	  numngb_inbox = ngb_treefind_variable(&Pos[0], L, target, &startnode, mode, nexport, nsend_local);
#else
	  numngb_inbox = ngb_treefind_variable(&Pos[0], L, target, &startnode, mode, nexport, nsend_local, 0);
#endif

	  if(numngb_inbox < 0)
	    return -1;

	  if(mode == LOCAL)
	    {
	      int n = *nexport - 1;
	      while((n >= 0) && (DataIndexTable[n].Index == target))
		{
		  memcpy(MetalDataIn[n].Pos, Pos, sizeof(MyLongDouble) * 3);
		  memcpy(MetalDataIn[n].Metals, metals, sizeof(float) * LT_NMet);
#ifdef GL_CR_DUST
		  memcpy(MetalDataIn[n].Dust, dust, sizeof(float) * LT_NMet);
		  MetalDataIn[n].numsnIa = numsnIa;
		  MetalDataIn[n].numsnII = numsnII;
#endif

		  MetalDataIn[n].L = L;
		  MetalDataIn[n].weight = weight;
		  MetalDataIn[n].energy = energy;
		  MetalDataIn[n].egystep = egystep;
		  MetalDataIn[n].SFi = mySFi;
#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
		  MetalDataIn[n].LMMass = LMmass;
#endif

		  memcpy(MetalDataIn[n].NodeList,
			 DataNodeList[DataIndexTable[n].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
		  MetalDataIn[n].Type = Type;
		  MetalDataIn[n].Task = DataIndexTable[n].Task;
		  MetalDataIn[n].Index = target;
		  MetalDataIn[n].IndexGet = n;
#ifdef LT_TRACK_CONTRIBUTES
		  MetalDataIn[n].contrib = contrib;	/* contrib is globa �l in this scope; it has been set in perform_stellarevolution_operations() */
#endif
		  n--;
		}
	    }

	  for(int n = 0; n < numngb_inbox; n++)
	    {
	      int j = Ngblist[n];

#ifdef BLACK_HOLES
	      if(P[j].Mass == 0)
		continue;
#endif

#if defined(GM_MUPPI) && defined(WINDS)
	      if(SphP[j].DelayTime0 > 0.0)
		continue;
#endif

	      double dx = NEAREST_X(Pos[0] - P[j].Pos[0]);
	      double dy = NEAREST_Y(Pos[1] - P[j].Pos[1]);
	      double dz = NEAREST_Z(Pos[2] - P[j].Pos[2]);
	      double dist = dx * dx + dy * dy + dz * dz;

	      if(dist < L2)
		{
		  double myweight = 1.0;
		  double u;
		  {
		    double dwk;
		    switch (Type)
		      {

		      case 0:
			dist = sqrt(dist);
			u = dist * linv;
			kernel_main(u, linv3, linv4, &myweight, &dwk, 0);
			break;

		      case 4:
#if defined( LT_USE_KERNEL_IN_WEIGHT ) || defined( LT_USE_DENSITY_IN_WEIGHT)
			dist = sqrt(dist);
			u = dist * linv;
			kernel_main(u, linv3, linv4, &myweight, &dwk, 0);
#if defined( LT_USE_DENSITY_IN_WEIGHT)
			myweight /= SphP[j].Density;
#endif
#endif
			break;

		      default:
			// obvisouly that should not happen
			break;
		      }
		  }


#if !defined(GM_MUPPI)

		  myweight *= P[j].Mass;

#else // GM_MUPPI defined
		  // GM: check this absolutely not sure! this is not supposed to be gas!

		  myweight *= ((SphP[j].MultiPhase > 0) ? SphP[j].Mgas0 : P[j].Mass);

#if defined(LT_USE_DENSITY_IN_WEIGHT)
		  myweight *= P[j].Mass;
		  myweight = (Type ? myweight / SphP[j].Density : myweight);
#endif

#endif /* ends GM_MUPPI */


		  double wfac = myweight / weight;




		  if(wfac > (1 + 1e-3))	/* > ============================================== < */
		    /* >  relative weight                               < */
		    {
		      printf
			("Something strange arose when calculating weight in metal distribution fpr part type %d [ %lf %lf %lf ] %d:\n"
			 "weightsum: %.8g <  wfac : %.8g\n  *  spreading length : %.8g\n"
			 "Task : %i, neighbour : %i, neigh id: %llu, neigh type: %i, neigh dist: %g, "
			 "particle : %i, part ID: %llu, mode: %i, numngb: %i\n", Type, Pos[0], Pos[1], Pos[2],
			 All.NumCurrentTiStep, weight, wfac, L, ThisTask, j, (unsigned long long) P[j].ID,
			 P[j].Type, u, target, (unsigned long long) ((mode == LOCAL) ? P[target].ID : 0),
			 mode, numngb_inbox);

		      fflush(stdout);
		      wfac = 1.0;
		      /*
		         endrun(101010);
		       */
		    }

		  else if(wfac < 0)	/* > ============================================== < */
		    /* >  should not happen                             < */
		    {
		      if(wfac > -5e-3)
			{
			  printf
			    ("warning: particle has got a negative weight factor! possibly due to a round-off error.\n"
			     "         [%d][%d][%d] %g %g %g %g\n", ThisTask, mode, j, u, myweight, weight,
			     wfac);
			  wfac = 0;
			}
		      else
			{
			  printf
			    ("warning: particle has got a negative weight factor! too large for being a round-off error.\n"
			     "         [%d][%d][%d] %g %g %g %g\n", ThisTask, mode, j, u, myweight, weight,
			     wfac);
			  endrun(919193);
			}
		    }

#ifdef LT_TRACK_CONTRIBUTES	/* > ============================================== < */
		  /* >  update tracks                                 < */
		  for(int k = 0; k < LT_NMetP; k++)
		    contrib_metals[k] = metals[k] * wfac;

		  update_contrib(&SphP[j].contrib, &SphP[j].Metals[0], &contrib, &contrib_metals[0]);
#endif


		  for(int k = 0; k < LT_NMetP; k++)	/* > ============================================== < */
		    /* >  check consistency                             < */
		    {
		      /* Hydrogen is the last element, we don't store it in the Metals array */
		      if((SphP[j].Metals[k] += metals[k] * wfac) < 0)
			{
			  printf(" \n ooops... it's a shame! %i %i %i %i %i %i %g %g %g \n",
				 ThisTask, mode, j, target, k, Type, metals[k], wfac, weight);
			  endrun(333333);
			}
#ifdef GL_CR_DUST
		      SphP[j].DustL[k] += dust[k] * wfac;
#endif

		    }

#ifdef GL_CR_DUST
		  SphP[j].numSnIa += numsnIa * wfac;
		  SphP[j].numSnII += numsnII * wfac;
#endif

		  double zage_term;

#ifdef LT_ZAGE			/* > ============================================== < */
		  /* >  update ZAGE                                   < */
		  zage_term = (cosmic_time - All.Time_Age) * metal_add_mass * wfac / P[j].Mass;
#ifdef LT_LOGZAGE
		  zage_term = log10(zage_term);
#endif
		  SphP[j].ZAge += zage_term;
		  SphP[j].ZAgeW += metal_add_mass * wfac / P[j].Mass;
#endif
#ifdef LT_ZAGE_LLV		/* > ============================================== < */
		  /* >  update ZAGE for Fe                            < */
		  zage_term = (cosmic_time - All.Time_Age) * llvmetal_add_mass * wfac / P[j].Mass;
#ifdef LT_LOGZAGE
		  zage_term = log10(zage_term);
#endif
		  SphP[j].ZAge_llv += zage_term;
		  SphP[j].ZAgeW_llv += llvmetal_add_mass * wfac / P[j].Mass;

/*                   SphP[j].ZAge += (cosmic_time - All.Time_Age) * get_metallicity(j, -1); */
/*                   SphP[j].ZAgeW += get_metallicity(j, -1); */

#endif

		  /* > ============================================== < */
		  /* >  FEEDBACK                                      < */
		  /* ===============================================================
		   *
		   * a feedback form
		   * ejecta from sn (not from agb!) are put in the hot phase of the
		   * gas; thie means:
		   *  (1) the current intrinsic energy is update supposing that
		   *      the added mass has the same erg/g than the hot phase
		   *  (2) the entropy or the entropy change rate is update accordingly
		   *
		   * also, you can choose to put the ejecta into the hot phase with
		   * some their own specific energy, either using the specific energy
		   * of the supernovae (1^51 erg / ejecta mass for each SN) or a
		   * specific energy that you specify in paramfile. These two options
		   * correspond to switching on either LT_SNegy_IN_HOTPHASE or
		   * LT_HOT_EJECTA.
		   *
		   * =============================================================== */

#ifdef LT_EJECTA_IN_HOTPHASE
		  /* calculate the current specific energy */
		  if(P[j].Ti_endstep == All.Ti_Current)
		    {
		      /* if this is an active particle, calculate from the end-step quantities */
		      dt = (P[j].Ti_endstep - P[j].Ti_begstep) * All.Timebase_interval;
		      current_egy = DMAX(All.MinEgySpec, (SphP[j].Entropy + SphP[j].DtEntropy * dt) /
					 GAMMA_MINUS1 * pow(SphP[j].Density * a3inv, GAMMA_MINUS1));
		    }
		  else
		    {
		      dt = 0;
		      current_egy =
			SphP[j].Entropy / GAMMA_MINUS1 * pow(SphP[j].Density * a3inv, GAMMA_MINUS1);
		    }

		  /* calculate the current specific energy of the hot phase */
		  current_egyhot = (current_egy - All.EgySpecCold * SphP[j].x) / (1 - SphP[j].x);
		  hotmass = P[j].Mass * (1 - SphP[j].x);

		  /* update the mass and the cold fraction */
		  if(SphP[j].x > 0)
		    {
		      x_ratio = SphP[j].x;
#ifndef LT_LOCAL_IRA
		      if(dist == 0 && !Type && mode == 0)
			/* the particle itself; should be sufficient dist == 0 */
			{
			  SphP[j].x = (SphP[j].x * P[j].Mass - add_mass) /
			    (P[j].Mass - (1 - wfac) * add_mass);
			  SphP[j].MassRes -= (1 - wfac) * add_mass;
			}
		      else
#endif
			{
			  SphP[j].x *= P[j].Mass / (P[j].Mass + add_mass * wfac);
			  SphP[j].MassRes += add_mass * wfac;
			}

		      x_ratio /= SphP[j].x;
		    }
		  else
		    {
		      x_ratio = 0;
#ifndef LT_LOCAL_IRA
		      if(dist == 0 && !Type && mode == 0)
			/* the particle itself; should be sufficient dist == 0 */
			SphP[j].MassRes -= (1 - wfac) * add_mass;
		      else
#endif
			SphP[j].MassRes += add_mass * wfac;
		    }

#if defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
		  v = (add_mass - LMmass) * wfac / (hotmass + add_mass * wfac);
#ifdef LT_HOT_EJECTA
		  current_egyhot = current_egyhot * (1 - v) + All.EgySpecEjecta * v;
#endif
#ifdef LT_SNegy_IN_HOTPHASE
		  current_egyhot = current_egyhot * (1 - v) + sn_spec_egy * v;
#endif
#endif
		  egy_ratio = current_egy;
		  current_egy = current_egyhot * (1 - SphP[j].x) + All.EgySpecCold * SphP[j].x;
		  egy_ratio /= current_egy;


		  if(dt > 0 && SphP[j].DtEntropy != 0)
		    {
		      SphP[j].DtEntropy =
			(current_egy * GAMMA_MINUS1 / pow(SphP[j].Density * a3inv, GAMMA_MINUS1) -
			 SphP[j].Entropy) / dt;
		      if(SphP[j].DtEntropy < -0.5 * SphP[j].Entropy / dt)
			SphP[j].DtEntropy = -0.5 * SphP[j].Entropy / dt;
		    }
		  else
		    SphP[j].Entropy = current_egy * GAMMA_MINUS1 / pow(SphP[j].Density * a3inv, GAMMA_MINUS1);
#if GADGET_HYDRO == HYDRO_MFM
		  SphP[j].InternalEnergy = eos->GetUFromEntropy(SphP[j].Density, SphP[i].Entropy, All.a3inv);
#endif

		  /* ===============================================================
		   *
		   * end of feedback
		   * =============================================================== */

#else /*  LT_EJECTA_IN_HOTPHASE */
		  /* > ============================================== < */
		  /* >  NORMAL FEEDBACK                               < */

		  {
		    double TotalEgy = SphP[j].EgyRes + energy * wfac;
#if defined(GM_MUPPI)
		    SphP[j].EgyStep =
		      (TotalEgy >
		       0 ? (energy * wfac * egystep +
			    SphP[j].EgyRes * SphP[j].EgyStep) / TotalEgy : SphP[j].EgyStep);
#endif

		    SphP[j].EgyRes = TotalEgy;
		  }

		  /* here we update the mass of the receveing particle if LT_EJECTA_IN_HOTPHASE
		     has not been used */

#ifndef LT_LOCAL_IRA
		  if(dist == 0 && !Type && mode == 0)
		    /* the particle itself; should be sufficient dist == 0 */
		    {

#ifdef GM_MUPPI
		      SphP[j].M_sf -= add_mass;
		      SphP[j].MassRes += add_mass * wfac;
#else
		      SphP[j].MassRes -= (1 - wfac) * add_mass;
#endif

		    }
		  else
#endif
		    SphP[j].MassRes += add_mass * wfac;
#endif /* > ============================================== < */
		  /* >  end of feedback                               < */


		  /* > ============================================== < */
		  /* >  collect INFOS                                 < */

		}
	    }
	}			/* closes inner while */
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = MetalDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }				/* closes outer while */


  return numngb_inbox;
}



/* / ---------------------------------------------- \
 * |                                                |
 * |   here we include lt_sn_calc.c                 |
 * |   it contains all the routines that actually   |
 * |   calculate the produced elements              |
 * |                                                |
 * |   don't include it in the OBJS section of the  |
 * |   Makefile                                     |
 * |                                                |
 * | SECTION III                                    |
 * |                                                |
 * | ........................                       |
 * |                                                |
 * | General Routines                               |
 * |                                                |
 * |    > get_SnIa_product                          |
 * |    > nRSnIa                                    |
 * |    > get_AGB_product                           |
 * |    > get_SnII_product                          |
 * |    > nRSnII                                    |
 * |    > mRSnII                                    |
 * |    > zmRSnII                                   |
 * |    > ztRSnII                                   |
 * |    > ejectaSnII                                |
 * |                                                |
 * \ ---------------------------------------------- / */


#include "lt_sn_calc.c"

/* / ---------------------------------------------- \
 * |                                                |
 * | SECTION IV                                     |
 * |                                                |
 * | ........................                       |
 * |                                                |
 * | General Routines                               |
 * |                                                |
 * |    > initialize_star_lifetimes                 |
 * |    > dm_dt                                     |
 * |    > lifetime                                  |
 * |    > dying_mass                                |
 * |    > sec_dist                                  |
 * |    > get_metallicity_stat                      |
 * |    > write_metallicity_stat                    |
 * |    > get_metals_sumcheck                       |
 * |                                                |
 * \ ---------------------------------------------- / */

void initialize_star_lifetimes(void)
     /*- initialize stellar lifetimes for given range of star masses -*/
{
  int i;

  All.mean_lifetime = lifetime(All.Mup);
  All.sup_lifetime = lifetime(All.MBms);
  All.inf_lifetime = All.sup_lifetime;
  if(ThisTask == 0)
    {
      printf("\nstellar mean lifetime:   mean (%5.4g Msun) = %g Gyrs\n\n", All.Mup, All.mean_lifetime);
      fflush(stdout);
    }

  for(i = 0; i < IMFs_dim; i++)
    {
      IMFs[i].inf_lifetime = lifetime(IMFs[i].MU);
      if(IMFs[i].Mm > All.MBms)
	IMFs[i].sup_lifetime = lifetime(IMFs[i].Mm);
      else
	IMFs[i].sup_lifetime = All.sup_lifetime;

      if(All.inf_lifetime > IMFs[i].inf_lifetime)
	All.inf_lifetime = IMFs[i].inf_lifetime;

      if(ThisTask == 0)
	{
	  printf("\n"
		 "   IMF %3d      inf (%5.4g Msun) = %g Gyrs\n"
		 "                sup (%5.4g Msun) = %g Gyrs\n",
		 i, IMFs[i].MU, IMFs[i].inf_lifetime, IMFs[i].Mm, IMFs[i].sup_lifetime);
	  fflush(stdout);
	}
    }

  if(ThisTask == 0)
    {
      printf("\nstellar mean lifetime:   mean (%5.4g Msun) = %g Gyrs\n\n", All.Mup, All.mean_lifetime);
      fflush(stdout);
    }
}


double INLINE_FUNC dm_dt(double m, double t)
{
  /* t is in Gyr */
  /* the last factor 1/agefact normalize in 1/yr: otherwise
     the results would be 1/Gyr */
  /*
     if(t > 0.0302233)
     return -0.37037 * m / (t - 0.012);
     else
     return -0.54054 * m / (t - 0.003);
   */

#ifdef LT_PM_LIFETIMES
  /* padovani & matteucci 1993 */
  if(t > 0.039765318659064693)
    return -m / t * (1.338 - 0.1116 * (9 + log10(t)));
  else
    return -0.545045045045 * m / (t - 0.003);
#endif

#ifdef LT_MM_LIFETIMES
  /* maeder & meynet 1989 */
  if(m <= 1.3)
    return -m / t / 0.6545;
  if(m > 1.3 && m <= 3)
    return -m / t / 3.7;
  if(m > 3 && m <= 7)
    return -m / t / 2.51;
  if(m > 7 && m <= 15)
    return -m / t / 1.78;
  if(m > 15 && m <= 53.054)
    return -m / t / 0.86;
  if(m > 53.054)
    return -0.54054054054 * m / (t - 0.003);
#endif

}

double INLINE_FUNC lifetime(double mass)
{
  /* calculates lifetime for a given mass  */
  /*                                       */
  /* padovani & matteucci (1993) approach  */
  /* move to gibson (1997) one             */
  /*                                       */
  /* mass is intended in solar units, life */
  /* time in Gyr                           */

  /* if(mass > 100) */
  /*   return 0; */
  /* else */

  if(mass < 0.6)
    return 160;			/* should be INF */


  /* padovani & matteucci 1993 */
#ifdef LT_PM_LIFETIMES
  if(mass <= 6.6)
    return pow(10, ((1.338 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) / 0.1116) - 9);
  else
    return 1.2 * pow(mass, -1.85) + 0.003;
#endif

  /*
     if(mass <= 8)
     return 5 * pow(mass, -2.7) + 0.012;
     else
     return 1.2 * pow(mass, -1.85) + 0.003;
   */

#ifdef LT_MM_LIFETIMES
  /* maeder & meynet 1989 */
  if(mass <= 1.3)
    return pow(10, -0.6545 * log10(mass) + 1);
  if(mass > 1.3 && mass <= 3)
    return pow(10, -3.7 * log10(mass) + 1.35);
  if(mass > 3 && mass <= 7)
    return pow(10, -2.51 * log10(mass) + 0.77);
  if(mass > 7 && mass <= 15)
    return pow(10, -1.78 * log10(mass) + 0.17);
  if(mass > 15 && mass <= 53.054)
    return pow(10, -0.86 * log10(mass) - 0.94);
  if(mass > 53.054)
    return 1.2 * pow(mass, -1.85) + 0.003;
#endif
}

double INLINE_FUNC dying_mass(double time)
{
  /* calculates mass dying at some time */
  /*                                    */
  /* time is time_in_Gyr                */

  if((time < All.inf_lifetime) || (time > All.sup_lifetime))
    return 0;

  /*
     if(time > 0.0302233)
     return pow((time - 0.012) / 5, -0.37037);
     else
     return pow((time - 0.003) / 1.2, -0.54054);
   */

#ifdef LT_PM_LIFETIMES
  /* padovani & matteucci 1993 */
  if(time > 0.039765318659064693)
    return pow(10, 7.764 - (1.79 - pow(1.338 - 0.1116 * (9 + log10(time)), 2)) / 0.2232);
  else
    return pow((time - 0.003) / 1.2, -1.0 / 1.85);
#endif

#ifdef LT_MM_LIFETIMES
  /* maeder & meynet 1989 */
  if(time >= 8.4221714076)
    return pow(10, (1 - log10(time)) / 0.6545);
  if(time < 8.4221714076 && time >= 0.38428316376)
    return pow(10, (1.35 - log10(time)) / 3.7);
  if(time < 0.38428316376 && time >= 0.044545508363)
    return pow(10, (0.77 - log10(time)) / 2.51);
  if(time < 0.044545508363 && time >= 0.01192772338)
    return pow(10, (0.17 - log10(time)) / 1.78);
  if(time < 0.01192772338 && time >= 0.0037734864318)
    return pow(10, -(0.94 + log10(time)) / 0.86);
  if(time < 0.0037734864318)
    return pow((time - 0.003) / 1.2, -0.54054);
#endif
}


double INLINE_FUNC sec_dist(double gamma, double mu)
{
  /* calculates secondary distribution function */
  /* for Sn type Ia                             */
  /* as far, we take gamma=2, so this function  */
  /* isn't used                                 */

  return pow(2, 1 + gamma) * (1 + gamma) * pow(mu, gamma);
}



/* .........................................
   * ------------------------------------- *
   |                                       |
   |   Chemical Stepping                   |
   * ------------------------------------- *
*/


double INLINE_FUNC get_NextChemTime(double lifetime, int mySFi, int *bin)
     /*!<
        returns the next look-back Time (in Gyr) for the
        chemical evolution. *bin will contain the ordinal
        number of the array's timestep in which the
        star currently lives.

        NOTE:: mySFi MUST be set to the appropriate value
      */
{
  int i;
  double diff;

  for(i = 1; i < Nsteps[mySFi]; i++)
    if(lifetime < SNtimesteps[mySFi][0][i + 1])
      break;

  if(bin != 0x0)
    *bin = i;

  if(i == Nsteps[mySFi])
    {
      if(bin != 0x0)
	*bin = Nsteps[mySFi] - 1;
      return FOREVER;
    }


  if((diff =
      (SNtimesteps[mySFi][0][i + 1] - lifetime) / (SNtimesteps[mySFi][0][i + 1] -
						   SNtimesteps[mySFi][0][i])) <= 0.4)
    /* if less than the 40% of the current chem timestep is left, use all the next timestep */
    diff = All.Time_Age - (SNtimesteps[mySFi][0][i + 2] - lifetime);
  else
    diff = All.Time_Age - (SNtimesteps[mySFi][0][i + 1] - lifetime);

  if(diff < 0)
    /* would mean at negative redshift, then set the last time right before z=0 */
    diff = 1e-4;

  return diff;
}


double INLINE_FUNC get_da_dota(double y, void *param)
     /* gives da/(da/dt) */
{
  //#ifdef LT_USE_SYSTEM_HUBBLE_FUNCTION
  return 1 / (y * sqrt(All.Omega0 * y * y * y + All.OmegaLambda +
		       (1 - All.Omega0 - All.OmegaLambda) * y * y));
  //#else
  //return 1.0 / (y * hubble_function(y) / All.Hubble);
  //#endif

}

double INLINE_FUNC get_age(double a)
     /* gives the look-back time corresponding to the expansion factor a */
{
  static double sec_in_Gyr = (86400.0 * 365.0 * 1e9);
  double result, error;

  if(All.ComovingIntegrationOn)
    {

      F.function = &get_da_dota;
      F.params = NULL;

      if((gsl_status =
	  gsl_integration_qag(&F, 1, 1 / a, 1e-4, 1e-5, gsl_ws_dim, qag_INT_KEY, w, &result, &error)))
	{
	  printf(">>>>> [%3d] qag integration error %d in get_age\n", ThisTask, gsl_status);
	  endrun(LT_ERROR_INTEGRATION_ERROR);
	}
      return 1.0 / (HUBBLE * All.HubbleParam) * result / sec_in_Gyr;
    }
  else
    {
      return a * All.UnitTime_in_s / sec_in_Gyr;
    }
}


/* :::........................................................................ */
/* *************************************************************************** */

#ifdef DOUBLEPRECISION
double INLINE_FUNC myfloor(double v)
{
  return floor(v);
}
#else
float INLINE_FUNC myfloor(double v)
{
  return floorf((float) v);
}
#endif



/* / ---------------------------------------------- \
 * |                                                |
 * | SECTION II                                     |
 * |                                                |
 * | ........................                       |
 * |                                                |
 * | Searching for Stars and Main Cycle             |
 * |                                                |
 * |    > count_evolving_stars                      |
 * |    > evolve_SN                                 |
 * |                                                |
 * \ ---------------------------------------------- / */


/* ...........................................................
   * ............................... *
   :                                 :
   :   Searching for Evolving Stars  :
   * ............................... *
*/


int is_chemically_active(int i)
{
  int ret_value = 0;

  switch (P[i].Type)
    {
    case 0:
#if !defined(LT_LOCAL_IRA)
#ifdef GM_MUPPI
      ret_value = (SphP[i].MultiPhase > 0 && SphP[i].mstar > 0 && SphP[i].DelayTime0 <= 0.0);
#else
      ret_value = (SphP[i].mstar > 0);
#endif
#endif
      break;
    case 4:
      ret_value = (TimeBinActive[MetP[P[i].pt.MetID].ChemTimeBin] && MPP(i).StellarAge < (MyFloat) All.Time);
      /* if( (MetP[P[i].pt.MetID].NextChemTime >= All.Time_Age) ) */
      /*  || ((All.Time_Age - MetP[I].NextChemTime) / (lifetime - MetP[I].LastChemTime) < 0.05) ) */
      /* condition 1 :
         the look-back time NextChemTime is larger than the present look-backtime
         condition 2 :
         the look-back time NextChemTime is lower than the present look-backtime, but
         the difference between the twos is less than 0.05 of the requested chemical timestep
       */
      break;
    }

  return ret_value;
}


void count_evolving_stars(int *num_of_stars, int *num_of_gas)
{
  int starsnum, gasnum, idx;

  NextChemActiveParticle = (int *) mymalloc("NextChemActiveParticle", NumPart * sizeof(int));

  starsnum = gasnum = idx = 0;

#pragma omp parallel
  {
    int flag, idx_l;
    long int i;
#pragma omp for reduction(+:gasnum,starsnum)
    for(i = 0; i < NumPart; i++)
      {
	flag = is_chemically_active(i);

	if(flag)
	  {
	    // capture old value of idx and update idx
#pragma omp atomic capture
	    idx_l = idx++;

	    NextChemActiveParticle[idx_l] = i;

	    if(P[i].Type == 0)
	      gasnum++;
	    else
	      starsnum++;
	    P[i].Type |= EVOLVE;

	  }
      }
  }
  *num_of_gas = gasnum;
  *num_of_stars = starsnum;

  return;
}


/* :::........................................................................ */
/* *************************************************************************** */


/* / ---------------------------------------------- \
 * |                                                |
 * | SECTION  V                                     |
 * |                                                |
 * | ........................                       |
 * |                                                |
 * | Initialization Routines                        |
 * |                                                |
 * |    > calculate_effective_yields                |
 * |    > calculate_FactorSN                        |
 * |    > init_SN                                   |
 * |    > build_SN_Stepping                         |
 * |    > get_Egy_and_Beta                          |
 * |    > calculate_ShortLiving_related             |
 * |    > setup_SF_related                          |
 * |    > TestStellarEvolution                      |
 * |                                                |
 * \ ---------------------------------------------- / */



#include "lt_sn_init.c"





void fsolver_error_handler(const char *reason, const char *file, int line, int err)
{
  if(err == GSL_EINVAL)
    {
      gsl_status = err;
      return;
    }
  return;
}

#include "lt_test_suite.c"

#endif

/*
 *
 * SCRATCH.NOTES AREA
 *

 *
 * END of SCRATCH.NOTES AREA
 *
 */

/*
 * here below and example of how to convert a time interval in code discretized time steps
 */

/* int INLINE_FUNC convert_time_to_timesteps(double start_a, double start_time, double delta_time) */
/* { */
/*   double delta_a; */

/*   /\* get delta_expansion_factor when moving from start_time to start_time + delta_time   * */
/*    * start_a is the expansion factor corresponding to start_time                         *\/ */
/*   delta_a = gsl_spline_eval(spline, cosmic_time - start_time + delta_time * 1.01, accel) - start_a; */

/*   /\* converts in code steps *\/ */
/*   return (int) (log(delta_a / start_a + 1) / All.Timebase_interval); */
/* } */
