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
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"
#include "../../Hydro/kernel.h"

#ifdef GM_MUPPI

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>
 

extern int NextParticle;
extern int Nexport, Nimport;
extern int BufferFullFlag;
extern int NextJ;
extern int TimerFlag;


/*! \file muppi_communications.c
*  \brief Redistribution of thermal and kinetic energy for MUPPI
*
*  This file contains the communication of energy among (source) MUPPI
*  particles, and (target) SPH neighbours.
*  Note, it's done with a SCATTER scheme.
*/


struct muppidata_in
{
  MyLongDouble Pos[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyAtLeastDouble E_out;
  MyAtLeastDouble GradDens[3];
  MyAtLeastDouble Eout_norm;
  MyIDType Nnp;

#ifdef MV_GM_STELLAR_KIN_FB2
  MyAtLeastDouble Ekin_norm_fountain;
#endif
#ifdef PARTICLE_DEBUG
  MyIDType ID;			/*!< particle identifier */
#endif

  int NodeList[NODELISTLENGTH];
}
 *MuppiDataIn, *MuppiDataGet;


struct muppidata_out
{
  MyAtLeastDouble Eout_norm;

#ifdef MV_GM_STELLAR_KIN_FB2
  MyAtLeastDouble Ekin_norm_fountain;
#endif

}
 *MuppiDataResult, *MuppiDataOut;


/***************************************/
/*   thermal energy exchange functions */
/***************************************/
void thermalenergy(void)
{
  int ngrp, ndone, ndone_flag;
  int recvTask, place;
  MyAtLeastDouble timecomp, timewait, timecomm;
  MyAtLeastDouble timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  MyAtLeastDouble tstart, tend, t0, t1;

  int save_NextParticle;

  long long n_exported = 0;

  /* allocate buffers to arrange communication */

  int NTaskTimesNumPart;
  int mainthreadid = 0;

  NTaskTimesNumPart = NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct muppidata_in) +
					     sizeof(struct muppidata_out) +
					     sizemax(sizeof(struct muppidata_in),
						     sizeof(struct muppidata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  t0 = second();

  NextParticle = FirstActiveParticle;	/* begin with this index */

  do
    {
      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;

      for(int j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();

      thermalenergy_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */

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
	      endrun(12998);
	    }


	  int new_export = 0;

	  
	  for(int j = 0, k = 0; j < Nexport; j++)
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

      for(int j = 0; j < NTask; j++)
	Send_count[j] = 0;
      for(int j = 0; j < Nexport; j++)
	Send_count[DataIndexTable[j].Task]++;

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
#endif

      tstart = second();
      
      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);

      tend = second();
      timewait1 += timediff(tstart, tend);
      {
	int j;
      for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}
      }
      MuppiDataGet = (struct muppidata_in *) mymalloc("MuppiDataGet", Nimport * sizeof(struct muppidata_in));
      MuppiDataIn = (struct muppidata_in *) mymalloc("MuppiDataIn", Nexport * sizeof(struct muppidata_in));

      /* prepare particle data for export */
      for(int j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(int k = 0; k < 3; k++)
	    {
	      MuppiDataIn[j].Pos[k] = P[place].Pos[k];
	      MuppiDataIn[j].GradDens[k] = SphP[place].GradDens[k];
	    }
	  MuppiDataIn[j].Hsml = P[place].Hsml;
	  MuppiDataIn[j].E_out = SphP[place].E_out;
	  MuppiDataIn[j].Eout_norm = SphP[place].Eout_norm;
	  MuppiDataIn[j].Nnp = SphP[place].Nnp;

#ifdef PARTICLE_DEBUG
	  MuppiDataIn[j].ID = P[place].ID;
#endif



#ifdef MV_GM_STELLAR_KIN_FB2
	  MuppiDataIn[j].Ekin_norm_fountain = SphP[place].Ekin_norm_fountain;
#endif

	  memcpy(MuppiDataIn[j].NodeList,
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
		  MPI_Sendrecv(&MuppiDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct muppidata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &MuppiDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct muppidata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);


      myfree(MuppiDataIn);
      MuppiDataResult =
	(struct muppidata_out *) mymalloc("MuppiDataResult", Nimport * sizeof(struct muppidata_out));
      MuppiDataOut =
	(struct muppidata_out *) mymalloc("MuppiDataOut", Nexport * sizeof(struct muppidata_out));

      /* now do the particles that were sent to us */

      tstart = second();

      NextJ = 0;

	  {
	    int mainthreadid = 0;
	    thermalenergy_evaluate_secondary(&mainthreadid);
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


      myfree(MuppiDataOut);
      myfree(MuppiDataResult);
      myfree(MuppiDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


#ifdef GM_MUPPI_DEBUG
  /* energy statistics */
  double TotEgyOut, TotEgyInThermal, TotEgyInKinetic;
  double TotEgyOutCum, TotEgyInThermalCum, TotEgyInKineticCum;
  All.Egy_outCum += All.Egy_out;
  All.Egy_in_thermalCum += All.Egy_in_thermal;
  All.Egy_in_kineticCum += All.Egy_in_kinetic;


  MPI_Allreduce(&All.Egy_out, &TotEgyOut, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&All.Egy_in_thermal, &TotEgyInThermal, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&All.Egy_in_kinetic, &TotEgyInKinetic, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);

  MPI_Allreduce(&All.Egy_outCum, &TotEgyOutCum, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&All.Egy_in_thermalCum, &TotEgyInThermalCum, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&All.Egy_in_kineticCum, &TotEgyInKineticCum, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);


  if(ThisTask==0) {
    FILE *f;
    
    f=fopen("FBenergy.txt","a");
    printf(" ...FB ENERGY STATISTICS:\n");
    printf("    Total SN energy in timestep, cumulative: %e\n",TotEgyOutCum);
    printf("....RECEIVED THERMAL energy: %e  (expected %e); fraction %f;\n", 
	   TotEgyInThermalCum, TotEgyOutCum*All.FracEgyOut, TotEgyInThermalCum/(TotEgyOutCum*All.FracEgyOut)); 
    printf("....RECEIVED KINETIC energy: %e  (expected %e); fraction %f;\n", 
	   TotEgyInKineticCum, TotEgyOutCum*All.FracEgyKin, TotEgyInKineticCum/(TotEgyOutCum*All.FracEgyKin)); 
    printf(".........................\n");
    printf("    Total SN energy in timestep, differential: %e\n",TotEgyOut);
    printf("....RECEIVED THERMAL energy: %e  (expected %e); fraction %f;\n", 
	   TotEgyInThermal, TotEgyOut*All.FracEgyOut, TotEgyInThermal/(TotEgyOut*All.FracEgyOut)); 
    printf("....RECEIVED KINETIC energy: %e  (expected %e); fraction %f;\n", 
	   TotEgyInKinetic, TotEgyOut*All.FracEgyKin, TotEgyInKinetic/(TotEgyOut*All.FracEgyKin)); 
    printf(".........................\n");
    fflush(stdout);
    fprintf(f,"%f %e %e   %e %e %e   %e %e %e   %e %e %e   %e %e %e\n", All.Time, TotEgyOut, TotEgyOutCum,
	    TotEgyInThermal, TotEgyOut*All.FracEgyOut, TotEgyInThermal/(TotEgyOut*All.FracEgyOut),
	    TotEgyInKinetic, TotEgyOut*All.FracEgyKin, TotEgyInKinetic/(TotEgyOut*All.FracEgyKin),
	    TotEgyInThermalCum, TotEgyOutCum*All.FracEgyOut, TotEgyInThermalCum/(TotEgyOutCum*All.FracEgyOut),
	    TotEgyInKineticCum, TotEgyOut*All.FracEgyKin, TotEgyInKineticCum/(TotEgyOutCum*All.FracEgyKin) );
    fclose(f);
  }
  All.Egy_out = All.Egy_in_thermal=All.Egy_in_kinetic = 0.0;
#endif

  int i;
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	SphP[i].wait = 0;
      }
  
  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1); /* total */

  timecomp = timecomp1 + timecomp2; /* computation */
  timewait = timewait1 + timewait2; /* waiting other procs */
  timecomm = timecommsumm1 + timecommsumm2; /* communications */

  CPU_Step[CPU_MUPPINET] += timeall; /* all toghether */
}



/*! This function is cloned from hydra.c
*  evaluates the normalization for MUPPI energy exchanges
*/
int thermalenergy_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
			   int *ngblist)
{
  int startnode, numngb, listindex = 0;
  int j, n, nangle;
  MyLongDouble *pos;
  MyAtLeastDouble *grad;
  MyFloat h;
  MyAtLeastDouble mass_j;

  MyAtLeastDouble dx, dy, dz;
  MyAtLeastDouble h2, hinv, hinv3, hinv4;
  MyAtLeastDouble wk, r, r2, u, dwk;

  MyAtLeastDouble r_a, r2_a;  
  MyAtLeastDouble u_a, wk_a, dwk_a;			
  MyAtLeastDouble pto_rx, pto_ry, pto_rz, dir_norm;   
  MyAtLeastDouble grad_tot, theta, phi, theta_1, theta_2, phi_1, phi_2;
  MyAtLeastDouble alpha, beta, alphap, alpham, betap, betam, eout_norm, e_out, E_tot_rec, fbenergy;


#ifdef PARTICLE_DEBUG
  MyIDType ID;			/*!< particle identifier */
#endif

  MyAtLeastDouble energy_kick;
#ifdef  MV_GM_STELLAR_KIN_FB2
  MyAtLeastDouble ekin_norm_fountain;
#endif

  MyIDType Nnp;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = P[target].Hsml;
      grad = SphP[target].GradDens;
      eout_norm = SphP[target].Eout_norm;
      e_out = SphP[target].E_out;
#ifdef PARTICLE_DEBUG
      ID=P[target].ID;			/*!< particle identifier */
#endif


#ifdef MV_GM_STELLAR_KIN_FB2
      ekin_norm_fountain = SphP[target].Ekin_norm_fountain;
#endif

      Nnp = SphP[target].Nnp;

    } 
  else
    {
      pos = MuppiDataGet[target].Pos;
      h = MuppiDataGet[target].Hsml;
      grad = MuppiDataGet[target].GradDens;
      eout_norm = MuppiDataGet[target].Eout_norm;
      e_out = MuppiDataGet[target].E_out;
#ifdef PARTICLE_DEBUG
      ID=MuppiDataGet[target].ID;			/*!< particle identifier */
#endif

#ifdef MV_GM_STELLAR_KIN_FB2
      ekin_norm_fountain = MuppiDataGet[target].Ekin_norm_fountain;
#endif

      Nnp = MuppiDataGet[target].Nnp;
    }

  h2 = h * h;
  hinv = 1.0 / h; 
  hinv3 = hinv * hinv * hinv;
  hinv4 = hinv3 * hinv;

  grad_tot = sqrt(grad[0]* grad[0] + grad[1]* grad[1] + 
		  grad[2]* grad[2]);

  theta = atan2(-grad[1],-grad[0]);
  phi =   acos(-grad[2]/grad_tot);

  theta_1 =  theta - ConeAngle_FB/2;
  theta_2 =  theta + ConeAngle_FB/2;

  phi_1 = phi - ConeAngle_FB/2;
  phi_2 = phi + ConeAngle_FB/2;
  
  nangle=0;


  /* initialize variables before SPH loop is started */
  h2 = h * h;
  /* Now start the actual SPH computation for this particle */
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = MuppiDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{

	  numngb =
	    ngb_treefind_variable_threads(pos, h, target, &startnode, mode, exportflag, exportnodecount,
					  exportindex, ngblist);
	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = ngblist[n];

	      pto_rx = -grad[2]*(P[j].Pos[1]-pos[1]) + grad[1]*(P[j].Pos[2]-pos[2]);
	      pto_ry = -grad[0]*(P[j].Pos[2]-pos[2]) + grad[2]*(P[j].Pos[0]-pos[0]);    
	      pto_rz = -grad[1]*(P[j].Pos[0]-pos[0]) + grad[0]*(P[j].Pos[1]-pos[1]);
	      dir_norm = grad[0]*grad[0] +grad[1]*grad[1] + grad[2]*grad[2];
 
	      r2= (pto_rx*pto_rx + pto_ry*pto_ry + pto_rz*pto_rz)/dir_norm;
	      r = sqrt(r2);

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      dx = NEAREST_X(dx);
	      dy = NEAREST_Y(dy);
	      dz = NEAREST_Z(dz);
#endif
	      dx = -dx;
	      dy = -dy;
	      dz = -dz;

	      r2_a = dx * dx + dy * dy + dz * dz;
	      r_a = sqrt(r2_a);

	      alpha = atan2(dy, dx);
	      beta= acos(dz/r_a);	 

	      alphap = alpha + 2.*M_PI;
	      alpham = alpha - 2.*M_PI;
	      betap = beta + M_PI;
	      betam = beta - M_PI;

	      if(SphP[j].MultiPhase>0) 
		{
		  mass_j = SphP[j].M_h + SphP[j].M_c;
		}
	      else 
		{
		  mass_j = P[j].Mass;
		  if(SphP[j].NMF>0) mass_j -= SphP[j].M_sf; //WRONG?
		}


	      /* here it distributes thermal energy */
	      /* if there is at least one particle within the cone, then distribute thermal energy */
	      if (eout_norm>0.0)
		{

		  if(r2 < h2 && r2_a < h2 && r_a > 0.0  &&
		     ( (theta_1 < alpha && alpha < theta_2) || (theta_1 < alphap && alphap < theta_2) ||  (theta_1 < alpham && alpham < theta_2) )&&
		     ( (phi_1 < beta && beta < phi_2) || (phi_1 < betap && betap < phi_2) || (phi_1 < betam && betam < phi_2) ) )
		      {

			nangle++;

			u = r * hinv;
			
			kernel_hinv(h, &hinv, &hinv3, &hinv4);
			kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);

			fbenergy = mass_j/ SphP[j].Density / eout_norm * e_out * wk * All.FracEgyOut;
			  
			
			SphP[j].E_rec += fbenergy;
			E_tot_rec +=  fbenergy;
			
#ifdef  GM_MUPPI_DEBUG
			All.Egy_in_thermal += fbenergy;
#endif

		      }
		}

	      /* if there are no particles in the cone, then give all THERMAL energy to the 
		 nearest particle to the axis (with ID Nnp).
		 Kinetic energy to fountain is not distributed here. */

	      else if (P[j].ID == Nnp)
		{

		  fbenergy = e_out * All.FracEgyOut;

		  SphP[j].E_rec += fbenergy;
		  E_tot_rec +=  fbenergy;
#ifdef  GM_MUPPI_DEBUG
		  All.Egy_in_thermal += fbenergy;
#endif
		}

#if defined(MV_GM_STELLAR_KIN_FB2)

	      kernel_hinv(h, &hinv, &hinv3, &hinv4);
	      u_a = r_a * hinv;          				/* milena */
	      kernel_main(u_a, hinv3, hinv4, &wk_a, &dwk_a, -1);       
	      if((r_a < h) && (r_a > 0.0000001))    /* milena */
		{

		  if (ekin_norm_fountain>0.0 && SphP[j].DelayTime>0.0 && SphP[j].wait!=1)
		    {
		      MyAtLeastDouble energy_kick;

		      energy_kick = mass_j/ SphP[j].Density / ekin_norm_fountain * e_out * wk_a * All.FracEgyKin;
		      
		      SphP[j].E_kin +=  energy_kick; 

#ifdef  GM_MUPPI_DEBUG
		      All.Egy_in_kinetic += energy_kick;
#endif

		      E_tot_rec += energy_kick;
		      
		    }
		}
#endif       /* MV_GM_STELLAR_KIN_FB2    milena */

	      

	    }    /* end of ngb loop, milena*/
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = MuppiDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  return 0;
}



void *thermalenergy_evaluate_primary(void *p)
{
  MyAtLeastDouble Eout;
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
      if(exitFlag)
	break;

      /* HERE WE HAVE THE CONTROL ON WHAT PARTICLE TO RUN */
      if(P[i].Type==0)
	Eout = SphP[i].E_out;
      else
	Eout = -1.0;
      if(P[i].Type == 0 &&  Eout>1.e-20)  
	{
	  if(thermalenergy_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}


void *thermalenergy_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;

  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;


  while(1)
    {
      {
	j = NextJ;
	NextJ++;
      }

      if(j >= Nimport)
	  break;

      if(MuppiDataGet[j].E_out <= 1.e-20) continue;
      thermalenergy_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);

    }

  return NULL;

}

#endif         /* closes GM_MUPPI */
