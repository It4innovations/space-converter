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
#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"
#include "../Hydro/kernel.h"
#include "blackhole.h"
#include "../System/communication.h"

/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

#ifdef BLACK_HOLES

void *blackhole_evaluate_swallow_secondary(void *p);
void *blackhole_evaluate_swallow_primary(void *p);

/***********************************************************************************/
/************************** Data Structures for loop (accreation) ******************/
/***********************************************************************************/

struct blackholeswallowdata_in
{
  MyLongDouble Pos[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat BH_Mass;
  MyIDType ID;
  int NodeList[NODELISTLENGTH];
};

struct blackholeswallowdata_out
{
  MyAtLeastDouble BH_accreted_Mass;
  MyAtLeastDouble BH_accreted_BH_Mass;
  MyAtLeastDouble AccretedMomentum[3];
  int BH_CountProgs;
  MyFloat Accreted_Age;
};

void particle2in_blackholeswallow(struct blackholeswallowdata_in *in, int i)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  in->Hsml = P[i].Hsml;
  in->Mass = P[i].Mass;
  in->BH_Mass = BPP(i).BH_Mass;
  in->ID = P[i].ID;
}

void out2particle_blackholeswallow(struct blackholeswallowdata_out *out, int i, int mode)
{
  ASSIGN_ADD(BPP(i).BH_accreted_Mass, out->BH_accreted_Mass, 1);
  ASSIGN_ADD(BPP(i).BH_accreted_BHMass, out->BH_accreted_BH_Mass, 1);

  for(int k = 0; k < 3; k++)
    ASSIGN_ADD(BPP(i).BH_accreted_momentum[k], out->AccretedMomentum[k], 1);

  ASSIGN_ADD(BPP(i).BH_CountProgs, out->BH_CountProgs, 1);

  if(BPP(i).StellarAge > out->Accreted_Age)
    BPP(i).StellarAge = out->Accreted_Age;
}


static struct blackholeswallowdata_in *BlackholeSwallowDataIn, *BlackholeSwallowDataGet;
static struct blackholeswallowdata_out *BlackholeSwallowDataResult, *BlackholeSwallowDataOut;

static int N_gas_swallowed, N_BH_swallowed;
static int N_gas_slices_swallowed;

void blackhole_swallow_loop(void)
{
  int ndone_flag, ndone;
  int Ntot_gas_swallowed, Ntot_BH_swallowed;
  int save_NextParticle;
  MPI_Status status;
  int Ntot_gas_slices_swallowed;

  N_gas_swallowed = N_BH_swallowed = 0;
  N_gas_slices_swallowed = 0;

  int NTaskTimesNumPart;

  NTaskTimesNumPart = ((int) (maxThreads)) * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  size_t MyBufferSize = KD_BUFFER_MANAGEMENT * FreeBytes / (1024.0 * 1024.0);

  All.BunchSize =
    (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					   sizeof(struct blackholeswallowdata_in) +
					   sizeof(struct blackholeswallowdata_out) +
					   sizemax(sizeof(struct blackholeswallowdata_in),
						   sizeof(struct blackholeswallowdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  VERBOSE(0, "BH: Start swallowing of gas particles and black holes");

  NextParticle = FirstActiveParticle;	/* first particle for this task */

  do
    {

      /* do local particles and prepare export list */
      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int thread_id = omp_get_thread_num();
#else
	int thread_id = 0;
#endif
	blackhole_evaluate_swallow_primary(&thread_id);	/* do local particles and prepare export list */
      }

      prepare_sendrec_offset_send(save_NextParticle, 1);	//1 means panic on bufer too small

      BlackholeSwallowDataGet =
	(struct blackholeswallowdata_in *) mymalloc("BlackholeDataGet",
						    Nimport * sizeof(struct blackholeswallowdata_in));
      BlackholeSwallowDataIn =
	(struct blackholeswallowdata_in *) mymalloc("BlackholeDataIn",
						    Nexport * sizeof(struct blackholeswallowdata_in));

      for(int j = 0; j < Nexport; j++)
	{
	  int place = DataIndexTable[j].Index;
	  particle2in_blackholeswallow(&BlackholeSwallowDataIn[j], place);
	  memcpy(BlackholeSwallowDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  int recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&BlackholeSwallowDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholeswallowdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeSwallowDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholeswallowdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MYMPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(BlackholeSwallowDataIn);
      BlackholeSwallowDataResult =
	(struct blackholeswallowdata_out *) mymalloc("BlackholeDataResult",
						     Nimport * sizeof(struct blackholeswallowdata_out));
      BlackholeSwallowDataOut =
	(struct blackholeswallowdata_out *) mymalloc("BlackholeDataOut",
						     Nexport * sizeof(struct blackholeswallowdata_out));

      /* now do the particles that were sent to us */
      NextJ = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int thread_id = omp_get_thread_num();
#else
	int thread_id = 0;
#endif
	blackhole_evaluate_swallow_secondary(&thread_id);	/* do local particles and prepare export list */
      }

      if(NextParticle < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

      /* get the result */
      for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  int recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&BlackholeSwallowDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholeswallowdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeSwallowDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholeswallowdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MYMPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the particles */
      for(int j = 0; j < Nexport; j++)
	{
	  int place = DataIndexTable[j].Index;
	  out2particle_blackholeswallow(&BlackholeSwallowDataOut[j], place, 1);
	}

      myfree(BlackholeSwallowDataOut);
      myfree(BlackholeSwallowDataResult);
      myfree(BlackholeSwallowDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&N_gas_slices_swallowed, &Ntot_gas_slices_swallowed, 1, MPI_INT, MPI_SUM, 0, MYMPI_COMM_WORLD);

  VERBOSE(0,
	  "BH: Accretion done, %d gas slices swallowed, %d last-slices swallowed (gas particles disappeared), %d BH particles swallowed",
	  Ntot_gas_slices_swallowed, Ntot_gas_swallowed, Ntot_BH_swallowed);

}



int blackhole_evaluate_swallow(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
			       int *ngblist)
{
  struct blackholeswallowdata_in local;
  struct blackholeswallowdata_out out;
  memset(&out, 0, sizeof(struct blackholeswallowdata_out));

  if(mode == 0)
    particle2in_blackholeswallow(&local, target);
  else
    local = BlackholeSwallowDataGet[target];

  out.Accreted_Age = 1;

  int startnode, listindex = 0;
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeSwallowDataGet[target].NodeList[listindex];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  int numngb =
	    ngb_treefind_blackhole_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag,
					   exportnodecount,
					   exportindex, ngblist);

	  PANIC_IF(numngb > NumPart / DENSITY_LESS_NGB_FACTOR,
		   "black hole: Number of neighbours found (%d) exceeds allowed number (%d) derived from NumPart/DENSITY_LESS_NGB_FACTOR, to continue, increase DENSITY_LESS_NGB_FACTOR (currently %f, max possible val would be 1) !!",
		   numngb, (int) (NumPart / DENSITY_LESS_NGB_FACTOR), DENSITY_LESS_NGB_FACTOR);


	  if(numngb < 0)
	    return -1;

	  for(int n = 0; n < numngb; n++)
	    {
	      int j = ngblist[n];

	      if(P[j].Type == 5 && P[j].Mass > 0)	/* we have a black hole merger */
		{
		  if(BPP(j).SwallowID == local.ID)
		    {
#ifdef _OPENMP
#pragma omp critical(_bh5)
#endif
		      {
			if(All.BlackHoleDetails >= 2)
			  fprintf(FdBlackHolesDetails,
				  "ThisTask=%d, time=%g: id=%llu swallows %llu (%g %g)\n",
				  ThisTask, All.Time, (unsigned long long) local.ID,
				  (unsigned long long) P[j].ID, local.BH_Mass, BPP(j).BH_Mass);

			out.BH_accreted_Mass += FLT(BPP(j).BH_Mass);
			out.BH_accreted_BH_Mass += FLT(BPP(j).BH_Mass);
			if(All.BlackHoleIgnoreMomentum == 0)
			  for(int k = 0; k < 3; k++)
			    out.AccretedMomentum[k] += FLT(BPP(j).BH_Mass * P[j].Vel[k]);

			out.BH_CountProgs += BPP(j).BH_CountProgs;

			int bin = P[j].TimeBin;
			TimeBin_BH_mass[bin] -= BPP(j).BH_Mass;
			TimeBin_BH_dynamicalmass[bin] -= P[j].Mass;
			TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
			if(BPP(j).BH_Mass > 0)
			  TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;

			P[j].Mass = 0;
			BPP(j).BH_Mass = 0;
			BPP(j).BH_Mdot = 0;

			out.Accreted_Age = BPP(j).StellarAge;

			N_BH_swallowed++;

		      }		// critical
		    }		// end of BH merger
		}

	      if(P[j].Type == 0)
		{
		  if(SPP[j].SwallowID == local.ID)
		    {
#ifdef _OPENMP
#pragma omp critical(_bh0)
#endif
		      {
			double myaccreted_mass = P[j].Mass;
			int number_of_slices_generated = 0;
			if(All.bits > 0 && All.BlackHoleAccretionSlicesOn == 1)
			  {
			    number_of_slices_generated = (P[j].ID >> (sizeof(MyIDType) * 8 - All.bits));
			    myaccreted_mass /= (GENERATIONS - number_of_slices_generated);
			  }
			if(number_of_slices_generated == (GENERATIONS - 1))
			  {
			    myaccreted_mass = P[j].Mass;
			    P[j].Mass = 0;
#ifdef LT_STELLAREVOLUTION
			    for(int k = 0; k < LT_NMetP; k++)
			      SphP[j].Metals[k] = 0;
#endif
			    N_gas_swallowed++;
			  }
			else
			  {
#ifdef LT_STELLAREVOLUTION
			    MyAtLeastDouble dec_factor = myaccreted_mass / P[j].Mass;
			    for(int k = 0; k < LT_NMetP; k++)
			      SphP[j].Metals[k] *= dec_factor;
#endif
#if GADGET_HYDRO == HYDRO_MFM
			    mfm_add_mass(j,-myaccreted_mass);
#endif
			    P[j].Mass -= myaccreted_mass;
			    P[j].ID += ((MyIDType) 1 << (sizeof(MyIDType) * 8 - All.bits));
			  }

			out.BH_accreted_Mass += myaccreted_mass;
			if(All.BlackHoleIgnoreMomentum == 0)
			  for(int k = 0; k < 3; k++)
			    out.AccretedMomentum[k] += FLT(myaccreted_mass * P[j].Vel[k]);

			N_gas_slices_swallowed++;
			SPP[j].SwallowID = 0;

		      }		// critical
		    }
		}
	    }			// for
	}			// while
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeSwallowDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }				// while

  /* Now collect the result at the right place */
  if(mode == 0)
    out2particle_blackholeswallow(&out, target, 0);
  else
    BlackholeSwallowDataResult[target] = out;

  return 0;
}

void *blackhole_evaluate_swallow_primary(void *p)
{
  int thread_id = *((int *) p);
  int i, j;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;

  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));

  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
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
#ifdef _OPENMP
	    do
	      {
		i = NextParticle;
		NextParticle = NextActiveParticle[i];
		if(NextParticle < 0)
		  break;
		ProcessedFlag[i] = 1;
	      }
	    while(P[i].Type != 5);
	    ProcessedFlag[i] = 0;
#else
	    i = NextParticle;
	    ProcessedFlag[NextParticle] = 0;
	    NextParticle = NextActiveParticle[i];
#endif
	  }
      }
      if(exitFlag)
	break;

      if(P[i].Type == 5)
	{
	  if(BPP(i).SwallowID == 0)
	    if(blackhole_evaluate_swallow(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	      {
		break;
	      }
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}

void *blackhole_evaluate_swallow_secondary(void *p)
{
  int thread_id = *((int *) p);
  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));

  while(1)
    {
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	j = NextJ;
	NextJ++;
      }

      if(j >= Nimport)
	break;

      blackhole_evaluate_swallow(j, 1, &dummy, &dummy, &dummy, ngblist);
    }
  return NULL;
}


#endif
