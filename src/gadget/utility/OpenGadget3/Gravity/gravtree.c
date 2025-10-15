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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/sem.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"

double Ewaldcount, Costtotal;
long long N_nodesinlist;
int Ewald_iter;			/* global in file scope, for simplicity */

#if !(defined ACC_INCLUDE_GRAVITY) && defined(ACC_GRAVITY)	//we reinclude this same file if ACC is active

#include "../OpenACC/gpuallvars_branch.h"
#include "../OpenACC/gravtree.h"

#endif


/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */


#ifdef STATICBRANDT
inline double OmegaR(int i, int mode)
{
  double r, m;
  r =
    sqrt((P[i].Pos[0] - LONG_X / 2.) * (P[i].Pos[0] - LONG_X / 2.) +
	 (P[i].Pos[1] - LONG_Y / 2.) * (P[i].Pos[1] - LONG_Y / 2.));
  if(r > 0)
    {
      switch (mode)
	{
	case -1:
	  m = 1;
	  break;
	case 0:
	  if(r < 7)
	    m = 1.0;
	  else
	    m = abs(sin((r - 7) * M_PI / (2 * All.Alfa2H)) + 1);
	  if(r > 8)
	    m = 1.0;
	  break;
	case 1:
	  m = (BRANDT_OmegaBr) / sqrt(1 + (r / BRANDT_R0) * (r / BRANDT_R0));
	  break;
	case 2:
	  if(r < BRANDT_R0)
	    m = BRANDT_OmegaBr;
	  else
	    m = BRANDT_OmegaBr * BRANDT_R0 / r;
	  break;
	default:
	  VERBOSE_ALL(0, "Wrong Immplementation of Omega\n");
	  exit(80085);
	  break;
	};
    }
  else
    m = 0;
  return m;
};
#endif

void sum_top_level_node_costfactors(void);



/*! This function computes the gravitational forces for all active particles.
 *  If needed, a new tree is constructed, otherwise the dynamically updated
 *  tree is used.  Particles are only exported to other processors when really
 *  needed, thereby allowing a good use of the communication buffer.
 */
void gravity_tree(void)
{
  long long n_exported = 0;
  int i, j, maxnumnodes, iter = 0;
  double t0, t1;
  double timeall = 0, timetree1 = 0, timetree2 = 0;
  double timetree, timewait, timecomm;
  double timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
  double sum_costtotal, ewaldtot;
  double maxt, sumt, maxt1, sumt1, maxt2, sumt2, sumcommall, sumwaitall;
  double plb, plb_max;

#ifdef FIXEDTIMEINFIRSTPHASE
  int counter;
  double min_time_first_phase, min_time_first_phase_glob;
#endif
#ifndef NOGRAVITY
  int k, ewald_max, diff, save_NextParticle;
  int ndone, ndone_flag, ngrp;
  int place;
  int recvTask;
  double tstart, tend;
  MPI_Status status;

#endif

#ifndef GRAVITY_CENTROID
  CPU_Step[CPU_MISC] += measure_time();

  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();

  /* construct tree if needed */
  if(TreeReconstructFlag)
    {
      VERBOSE(2, "Tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

      on_end_of_timestep();	/* used by OpenACC, empty otherwise */



      CPU_Step[CPU_MISC] += measure_time();

      force_treebuild(NumPart, NULL);

      CPU_Step[CPU_TREEBUILD] += measure_time();

      TreeReconstructFlag = 0;

      on_beginning_of_timestep();


      VERBOSE(2, "Tree construction done.\n");
    }
#endif


#ifdef ACC_GRAVITY
  int is_primary_loop_on_accelerator = 0;
  int should_drift = 1;
  acc_all.go_gpu = NActivePart > ACC_ACTIVEPART_CPU;

  openacc_gravity_drift(is_primary_loop_on_accelerator, should_drift);

#endif //ACC


#ifndef NOGRAVITY


#ifdef ADAPTGRAVSOFT
  ags_density();
  ags_force_update_hmax();
#endif

  /* allocate buffers to arrange communication */
  VERBOSE(2, "Begin tree force.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  size_t MyBufferSize = 0.5 * KD_BUFFER_MANAGEMENT * FreeBytes / (1024.0 * 1024.0);



  All.BunchSize =
    (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					   sizeof(struct gravdata_in) + sizeof(struct gravdata_out) +
					   sizemax(sizeof(struct gravdata_in), sizeof(struct gravdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));

  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeListIn", All.BunchSize * sizeof(struct data_nodelist));
#ifdef AR_NO_NODELIST_IN_DATASTRUCT
  struct data_nodelist *DataNodeListOut =
    (struct data_nodelist *) mymalloc("DataNodeListOut", All.BunchSize * sizeof(struct data_nodelist));
#endif

  VERBOSE(2, "Gravity: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
	  FreeBytes / (1024.0 * 1024.0));

  VERBOSE(2, "All.BunchSize=%d\n", All.BunchSize);

  Ewaldcount = 0;
  Costtotal = 0;
  N_nodesinlist = 0;


  CPU_Step[CPU_TREEMISC] += measure_time();
  t0 = second();

#if defined(PERIODIC) && !defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC)
  ewald_max = 1;
#else
  ewald_max = 0;
#endif

#ifdef SCF_HYBRID
  int scf_counter, max_scf_counter = 1;

  if(SCF_HYBRID == 2)
    max_scf_counter = 0;
  /*
     calculates the following forces (depending on SCF_HYBRID value)
     STAR<->STAR, STAR->DM (scf_counter=0)
     DM<->DM (scf_counter=1)
   */

  for(scf_counter = 0; scf_counter <= max_scf_counter; scf_counter++)
    {
      /* set DM mass to zero and set gravsum to zero */
      if(scf_counter == 0)
	{
	  for(i = 0; i < NumPart; i++)
	    {
	      if(P[i].Type == 1)	/* DM particle */
		P[i].Mass = 0.0;

	      for(j = 0; j < 3; j++)
		P[i].GravAccelSum[j] = 0.0;
	    }
	}
      /* set stellar mass to zero */
      if(scf_counter == 1)
	{
	  for(i = 0; i < NumPart; i++)
	    {
	      if(P[i].Type == 2)	/* stellar particle */
		P[i].Mass = 0.0;
	    }
	}

      /* particle masses changed, so reconstruct tree */
      VERBOSE(2, "SCF Tree construction %d\n", scf_counter);
      force_treebuild(NumPart, NULL);
      VERBOSE(2, "done.\n");
#endif

#ifdef KD_TREEDOMAIN_UPDATE_FREQENCY
      if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency / All.Time * All.TotNumPart)
#else
      if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)
#endif
	{
	  /* we have a fresh tree and would like to measure gravity cost */

	  /* find the closest level */
	  for(i = 1, TakeLevel = 0, diff = abs(All.LevelToTimeBin[0] - All.HighestActiveTimeBin);
	      i < GRAVCOSTLEVELS; i++)
	    {
	      if(diff > abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin))
		{
		  TakeLevel = i;
		  diff = abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin);
		}
	    }

	  if(diff != 0)		/* we have not found a matching slot */
	    {

	      if(All.HighestOccupiedTimeBin - All.HighestActiveTimeBin < GRAVCOSTLEVELS)	/* we should have space */
		{
		  /* clear levels that are out of range */
		  for(i = 0; i < GRAVCOSTLEVELS; i++)
		    {
		      if(All.LevelToTimeBin[i] > All.HighestOccupiedTimeBin)
			All.LevelToTimeBin[i] = 0;
		      if(All.LevelToTimeBin[i] < All.HighestOccupiedTimeBin - (GRAVCOSTLEVELS - 1))
			All.LevelToTimeBin[i] = 0;
		    }
		}

	      for(i = 0, TakeLevel = -1; i < GRAVCOSTLEVELS; i++)
		{
		  if(All.LevelToTimeBin[i] == 0)
		    {
		      All.LevelToTimeBin[i] = All.HighestActiveTimeBin;
		      TakeLevel = i;
		      break;
		    }
		}

	      if(TakeLevel < 0 && All.HighestOccupiedTimeBin - All.HighestActiveTimeBin < GRAVCOSTLEVELS)	/* we should have space */
		PANIC("TakeLevel < 0, even though we should have a slot");
	    }
	}
      else
	{
	  /* in this case we do not measure gravity cost. Check whether this time-level
	     has previously mean measured. If yes, then delete it so to make sure that it is not out of time */

	  for(i = 0; i < GRAVCOSTLEVELS; i++)
	    if(All.LevelToTimeBin[i] == All.HighestActiveTimeBin)
	      All.LevelToTimeBin[i] = 0;

	  TakeLevel = -1;
	}




#ifdef ACC_GRAVITY
      struct gravdata_out *LitePOut;
      double *PGravCosts;
      double *NodeGravCosts;


      int NodesSize = MaxNodes;	//Numnodestree;
      openacc_gravity_init(LitePOut, PGravCosts, NodeGravCosts, NodesSize);
#endif //#ifdef ACC_GRAVITY

      if(TakeLevel >= 0)
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i = 0; i < NumPart; i++)
	  P[i].GravCost[TakeLevel] = 0;


      for(Ewald_iter = 0; Ewald_iter <= ewald_max; Ewald_iter++)
	{
	  /* this pragma is ignored if there is no openacc support. But we need that extra block {} here in gravtree */
	  PUSH_T("ewald");
	  NextParticle = FirstActiveParticle;	/* beginn with this index */


	  PRAGMA_GRAVITY_DATA_REGION()	//vuoto se non c'e'  ACC, ma se no,crea un blocco data di openacc
	  {

#ifdef ACC_GRAVITY

	    //note: costTotal_GPU_I  is defined in OpenACC/*.h and zeroed in the gravity init function
	    openacc_gravity_primary_loop(LitePOut, PGravCosts, NodeGravCosts, NodesSize);

#endif

	    do
	      {
		iter++;
		BufferFullFlag = 0;
		Nexport = 0;
		save_NextParticle = NextParticle;

		tstart = second();


		PUSH_T("primary");
/* do local particles and prepare export list */
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
		  int mainthreadid = ThisThread;
#ifdef ACC_GRAVITY
		  //gravity_primary_loop(&mainthreadid);
		  gravity_primary_loop_vect(&mainthreadid, is_primary_loop_on_accelerator, should_drift);
		  //                      gravity_primary_loop_vect(&mainthreadid, 0, should_drift);
		  NextParticle = -1;
#else
		  gravity_primary_loop(&mainthreadid);
#endif
		}
		POP_T();


		tend = second();
		timetree1 += timediff(tstart, tend);


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
		    for(j = 0; j < NTask; j++)
		      Recv_count[j] = 0;
		    MPI_Win_create(Recv_count, NTask * sizeof(int), sizeof(int), MPI_INFO_NULL,
				   MYMPI_COMM_WORLD, &win);
		    MPI_Win_fence(0, win);
		    for(j = 0; j < NTask; ++j)
		      {
			if(Send_count[j] > 0)
			  MPI_Put(&Send_count[j], 1, MPI_INT, j, ThisTask, 1, MPI_INT, win);
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

		GravDataGet =
		  (struct gravdata_in *) mymalloc("GravDataGet", Nimport * sizeof(struct gravdata_in));
		GravDataIn =
		  (struct gravdata_in *) mymalloc("GravDataIn", Nexport * sizeof(struct gravdata_in));

		/* prepare particle data for export */

		for(j = 0; j < Nexport; j++)
		  {
		    place = DataIndexTable[j].Index;

#ifdef GRAVITY_CENTROID
		    if(P[place].Type == 0)
		      {
			for(k = 0; k < 3; k++)
			  GravDataIn[j].Pos[k] = SphP[place].Center[k];
		      }
		    else
		      {
			for(k = 0; k < 3; k++)
			  GravDataIn[j].Pos[k] = P[place].Pos[k];
		      }
#else
		    for(k = 0; k < 3; k++)
		      GravDataIn[j].Pos[k] = P[place].Pos[k];
#endif
		    GravDataIn[j].Type = P[place].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		    if(P[place].Type == 0)
		      GravDataIn[j].Soft = P[place].Hsml;
#endif
		    GravDataIn[j].OldAcc = P[place].OldAcc;

#ifdef ADAPTGRAVSOFT
		    GravDataIn[j].AGS_zeta = P[place].AGS_zeta;
		    GravDataIn[j].AGS_omega = P[place].AGS_omega;
		    GravDataIn[j].AGS_Hsml = P[place].AGS_Hsml;
		    GravDataIn[j].Mass = P[place].Mass;
#endif

		    memcpy(GravDataIn[j].NodeList,
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
			    MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]],
					 Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
					 recvTask, TAG_GRAV_A,
					 &GravDataGet[Recv_offset[recvTask]],
					 Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
					 recvTask, TAG_GRAV_A, MYMPI_COMM_WORLD, &status);

			  }
		      }
		  }
		tend = second();
		timecommsumm1 += timediff(tstart, tend);


		myfree(GravDataIn);
		GravDataResult =
		  (struct gravdata_out *) mymalloc("GravDataResult", Nimport * sizeof(struct gravdata_out));
		GravDataOut =
		  (struct gravdata_out *) mymalloc("GravDataOut", Nexport * sizeof(struct gravdata_out));

		report_memory_usage(&HighMark_gravtree, "GRAVTREE");

		PUSH_T("secondary");

		/* now do the particles that were sent to us */
		tstart = second();
		NextJ = 0;

#ifdef ACC_GRAVITY
		openacc_gravity_secondary_loop(should_drift, NodesSize);

#else //#ifdef ACC_GRAVITY
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
		  int mainthreadid = ThisThread;
		  gravity_secondary_loop(&mainthreadid);
		}
#endif //#ifdef ACC_GRAVITY
		POP_T();
		tend = second();
		timetree2 += timediff(tstart, tend);

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
			    MPI_Sendrecv(&GravDataResult[Recv_offset[recvTask]],
					 Recv_count[recvTask] * sizeof(struct gravdata_out),
					 MPI_BYTE, recvTask, TAG_GRAV_B,
					 &GravDataOut[Send_offset[recvTask]],
					 Send_count[recvTask] * sizeof(struct gravdata_out),
					 MPI_BYTE, recvTask, TAG_GRAV_B, MYMPI_COMM_WORLD, &status);
			  }
		      }

		  }
		tend = second();
		timecommsumm2 += timediff(tstart, tend);

		/* add the results to the local particles */
		tstart = second();
		for(j = 0; j < Nexport; j++)
		  {
		    place = DataIndexTable[j].Index;

		    for(k = 0; k < 3; k++)
		      P[place].GravAccel[k] += GravDataOut[j].Acc[k];

#ifdef DISTORTIONTENSORPS
		    int i1, i2;
		    for(i1 = 0; i1 < 3; i1++)
		      for(i2 = 0; i2 < 3; i2++)
			P[place].tidal_tensorps[i1][i2] += GravDataOut[j].tidal_tensorps[i1][i2];
#endif

#ifdef EVALPOTENTIAL
		    P[place].Potential += GravDataOut[j].Potential;
#endif

#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
		    P[place].AGS_corr += GravDataOut[j].AGS_corr;
#endif
		  }
		tend = second();
		timetree1 += timediff(tstart, tend);

		myfree(GravDataOut);
		myfree(GravDataResult);
		myfree(GravDataGet);
	      }
	    while(ndone < NTask);
	    POP_T();
#ifdef ACC_GRAVITY
#pragma acc wait(1)
#endif
	  }			//exit GPU region
#ifdef ACC_GRAVITY

	  if(acc_all.go_gpu)
	    {
	      //truct gravdata_out *LitePOut, double *PGravCosts, double*NodeGravCosts, double costTotal, int NodesSize)
	      openacc_gravity_join_results(LitePOut, PGravCosts, NodeGravCosts, NodesSize);
	    }
#endif
	  POP_T();
	}			/* Ewald_iter */

#ifdef SCF_HYBRID
      /* restore particle masses */
      for(i = 0; i < NumPart; i++)
	P[i].Mass = P[i].MassBackup;


      /* add up accelerations from tree to AccelSum */
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
	  /* ignore STAR<-DM contribution */
	  if(scf_counter == 1 && P[i].Type == 2)
	    {
	      continue;
	    }
	  else
	    {
	      for(j = 0; j < 3; j++)
		P[i].GravAccelSum[j] += P[i].GravAccel[j];
	    }
	}
    }				/* scf_counter */

  /* set acceleration to summed up accelerations */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = P[i].GravAccelSum[j];
    }
#endif



#ifdef ACC_GRAVITY

  if(acc_all.go_gpu)
    {

      myfree(LitePOut);
      myfree(PGravCosts);
      myfree(NodeGravCosts);
    }
#endif

  myfree(DataNodeList);
  myfree(DataIndexTable);


  /* assign node cost to particles */
  if(TakeLevel >= 0)
    {
      sum_top_level_node_costfactors();

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(i = 0; i < NumPart; i++)
	{
#ifdef NEUTRINOS
	  if((P[i].Type != 2) || (All.Time > All.Time_tree_on_nu))
#endif
	    {
	      int no = Father[i];

	      while(no >= 0)
		{
		  if(Nodes[no].u.d.mass > 0)
		    P[i].GravCost[TakeLevel] += Nodes[no].GravCost * P[i].Mass / Nodes[no].u.d.mass;

		  no = Nodes[no].u.d.father;
		}
	    }
	}
    }

  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      double fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
	}
    }
#endif
#endif

  /* to prevent that we overwrite OldAcc in the first evaluation for 2lpt ICs */
  if(!(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0))
    {
#ifdef _OPENMP
      int il;
#pragma omp parallel for private(i)
      for(il = 0; il < NActivePart; il++)
	{
	  i = ActiveParticleList[il];
#else /* _OPENMP */
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
#endif
	  double ax, ay, az;
#ifdef PMGRID
	  ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
	  ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
	  az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
	  ax = P[i].GravAccel[0];
	  ay = P[i].GravAccel[1];
	  az = P[i].GravAccel[2];
#endif
	  P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
	}			// for
    }				// if

  if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
    {
      if(!(All.Ti_Current == 0 && RestartFlag == 0))
	if(All.TypeOfOpeningCriterion == 1)
	  All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */
    }
  else
    {
      if(All.TypeOfOpeningCriterion == 1)
	All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */
    }

  /*  muliply by G */
#ifdef _OPENMP
  int il;
#pragma omp parallel for private(i,j,k)	// if(NActivePart>40)
  for(il = 0; il < NActivePart; il++)
    {
      i = ActiveParticleList[il];
#else /* _OPENMP */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif
      for(j = 0; j < 3; j++)
	{
	  P[i].GravAccel[j] *= All.G;
	}


#ifdef DISTORTIONTENSORPS
      /*
         Diaganol terms of tidal tensor need correction, because tree is running over
         all particles -> also over target particle -> extra term -> correct it
       */
      /* 3D -> full forces */
      P[i].tidal_tensorps[0][0] +=
	P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] *
		     All.ForceSoftening[P[i].Type]) * 10.666666666667;

      P[i].tidal_tensorps[1][1] +=
	P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] *
		     All.ForceSoftening[P[i].Type]) * 10.666666666667;

      P[i].tidal_tensorps[2][2] +=
	P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] *
		     All.ForceSoftening[P[i].Type]) * 10.666666666667;

      if(All.ComovingIntegrationOn)
	{
	  P[i].tidal_tensorps[0][0] -= All.TidalCorrection / All.G;
	  P[i].tidal_tensorps[1][1] -= All.TidalCorrection / All.G;
	  P[i].tidal_tensorps[2][2] -= All.TidalCorrection / All.G;
	}
      /*now muliply by All.G */
      int j1, j2;
      for(j1 = 0; j1 < 3; j1++)
	for(j2 = 0; j2 < 3; j2++)
	  P[i].tidal_tensorps[j1][j2] *= All.G;
#endif /* DISTORTIONTENSORPS */

#ifdef EVALPOTENTIAL
      /* remove self-potential */
      P[i].Potential += P[i].Mass / All.SofteningTable[P[i].Type];

      if(All.ComovingIntegrationOn)
	if(All.PeriodicBoundariesOn)
	  P[i].Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
	    pow(All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G), 1.0 / 3);

      P[i].Potential *= All.G;

#ifdef PMGRID
      P[i].Potential += P[i].PM_Potential;	/* add in long-range potential */
#endif

      if(All.ComovingIntegrationOn)
	{
#ifndef PERIODIC
	  double fac, r2;

	  fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

	  for(k = 0, r2 = 0; k < 3; k++)
	    r2 += P[i].Pos[k] * P[i].Pos[k];

	  P[i].Potential += fac * r2;
#endif
	}
      else
	{
	  double fac, r2;

	  fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;

	  if(fac != 0)
	    {
	      for(k = 0, r2 = 0; k < 3; k++)
		r2 += P[i].Pos[k] * P[i].Pos[k];

	      P[i].Potential += fac * r2;
	    }
	}
#endif
    }				// for

  /* Finally, the following factor allows a computation of a cosmological simulation
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      double fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	for(j = 0; j < 3; j++)
	  P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif


  VERBOSE(2, "tree is done.\n");

#else /* gravity is switched off */
  t0 = second();

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    for(j = 0; j < 3; j++)
      P[i].GravAccel[j] = 0;

#ifdef DISTORTIONTENSORPS
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      P[i].tidal_tensorps[0][0] = 0.0;
      P[i].tidal_tensorps[0][1] = 0.0;
      P[i].tidal_tensorps[0][2] = 0.0;
      P[i].tidal_tensorps[1][0] = 0.0;
      P[i].tidal_tensorps[1][1] = 0.0;
      P[i].tidal_tensorps[1][2] = 0.0;
      P[i].tidal_tensorps[2][0] = 0.0;
      P[i].tidal_tensorps[2][1] = 0.0;
      P[i].tidal_tensorps[2][2] = 0.0;
    }
#endif
#endif /* end of NOGRAVITY */

#ifdef NOGRAVITY
  int k;
#endif

#ifdef CONSTANT_GRAVITY_Y
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef SPH_BND_PARTICLES
      if(P[i].ID != 0)
#endif
	P[i].GravAccel[1] += CONSTANT_GRAVITY_Y;
    }
#endif

#ifdef SCFPOTENTIAL
  MyDouble xs, ys, zs;
  MyDouble pots, axs, ays, azs;

  VERBOSE(2, "Starting SCF calculation...\n");

  /* reset the expansion coefficients to zero */
  SCF_reset();
#ifdef SCF_HYBRID
  /*
     calculate SCF coefficients for local DM particles.
     sum them up from all processors, so every processor
     sees the same expansion coefficients
   */
  SCF_calc_from_particles();

  /* sum up local coefficients */
  MPI_Allreduce(sinsum, sinsum_all, (SCF_NMAX + 1) * (SCF_LMAX + 1) * (SCF_LMAX + 1), MPI_DOUBLE, MPI_SUM,
		MYMPI_COMM_WORLD);
  MPI_Allreduce(cossum, cossum_all, (SCF_NMAX + 1) * (SCF_LMAX + 1) * (SCF_LMAX + 1), MPI_DOUBLE, MPI_SUM,
		MYMPI_COMM_WORLD);

  /* update local coefficients to global coefficients -> every processor has now complete SCF expansion */
  SCF_collect_update();
  VERBOSE(2, "calculated and collected coefficients.\n");

#else
  long old_seed, global_seed_min, global_seed_max;

  /*
     resample coefficients for expansion
     make sure that every processors sees the SAME potential,
     i.e. has the same seed to generate coefficients
   */
  old_seed = scf_seed;
  SCF_calc_from_random(&scf_seed);
  /* check that all cpus have the same random seed (min max must be the same) */
  MPI_Allreduce(&scf_seed, &global_seed_max, 1, MPI_LONG, MPI_MAX, MYMPI_COMM_WORLD);
  MPI_Allreduce(&scf_seed, &global_seed_min, 1, MPI_LONG, MPI_MIN, MYMPI_COMM_WORLD);
  VERBOSE(2, "sampled coefficients with old/new seed = %ld/%ld         min/max=%ld/%ld\n", old_seed, scf_seed,
	  global_seed_min, global_seed_max);

#endif


  /* get accelerations for all active particles based on current expansion */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      /* convert to unit sphere */
      to_unit(P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], &xs, &ys, &zs);
      /* OR: not */
      //xs = P[i].Pos[0]; ys = P[i].Pos[1]; zs = P[i].Pos[2];

      /* evaluate potential and acceleration */
      SCF_evaluate(xs, ys, zs, &pots, &axs, &ays, &azs);

      /* scale to system size and add to acceleration */
#ifdef SCF_HYBRID
      /*
         add missing STAR<-DM force from SCF (was excluded in tree above)
       */
      if(P[i].Type == 2 || SCF_HYBRID == 2)
	{
#endif
	  /* scale */
	  P[i].GravAccel[0] += All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * axs;
	  P[i].GravAccel[1] += All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * ays;
	  P[i].GravAccel[2] += All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * azs;
	  /* OR: not */
	  //P[i].GravAccel[0] += All.G * axs;
	  //P[i].GravAccel[1] += All.G * ays;
	  //P[i].GravAccel[2] += All.G * azs;

#ifdef DEBUG
	  if(P[i].ID == 150000)
	    {
	      VERBOSE_ALL(2, "SCF-ACCEL (scf)   %d  (%g|%g|%g)\n", All.NumCurrentTiStep,
			  All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * axs,
			  All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * ays,
			  All.G * SCF_HQ_MASS / (SCF_HQ_A * SCF_HQ_A) * azs);
	      /* analyic potential of zeroth order of expansion */
	      sphere_acc(xs, ys, zs, &axs, &ays, &azs);
	      VERBOSE_ALL(2, "SCF-ACCEL (exact) %d  (%g|%g|%g)\n", All.NumCurrentTiStep, All.G * axs,
			  All.G * ays, All.G * azs);
	    }
#endif

#ifdef SCF_HYBRID
	}
#endif
    }

  VERBOSE(2, "done.\n");

#endif


#ifdef STATICISO
  {
    double r, m;
    double dx, dy, dz;

    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
      {
	dx = P[i].Pos[0];
	dy = P[i].Pos[1];
	dz = P[i].Pos[2];

	r = sqrt(dx * dx + dy * dy + dz * dz);

	if(r > ISO_R200)
	  m = ISO_M200;
	else
	  m = ISO_M200 * r / ISO_R200;

#ifdef ISO_FRACTION
	m *= ISO_FRACTION;
#endif
	if(r > 0)
	  {
	    P[i].GravAccel[0] += -All.G * m * dx / r / (r * r + ISO_Eps * ISO_Eps);
	    P[i].GravAccel[1] += -All.G * m * dy / r / (r * r + ISO_Eps * ISO_Eps);
	    P[i].GravAccel[2] += -All.G * m * dz / r / (r * r + ISO_Eps * ISO_Eps);
	  }
      }
  }
#endif

#ifdef DISKPOT
  gravity_tree_subfunc_diskpot();
#endif



#ifdef GROWING_DISK_POTENTIAL
  {
    double mdisk, dx, dy, dz, r, z, aR, az;

    growing_disk_init();

    mdisk = get_disk_mass(All.Time);

    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
      {
	dx = P[i].Pos[0];
	dy = P[i].Pos[1];
	dz = P[i].Pos[2];

	r = sqrt(dx * dx + dy * dy);
	z = fabs(dz);

	get_disk_forces(r, z, &aR, &az);

	aR *= mdisk;
	az *= mdisk;

	if(r > 0)
	  {
	    P[i].GravAccel[0] += -dx / r * aR;
	    P[i].GravAccel[1] += -dy / r * aR;
	    P[i].GravAccel[2] += -dz / z * az;
	  }
      }
  }
#endif


#ifdef STATICNFW
  double r, m;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      double dx, dy, dz;

#ifdef NFW_BOXCENTERED
      dx = P[i].Pos[0] - boxHalf_X;
      dy = P[i].Pos[1] - boxHalf_Y;
      dz = P[i].Pos[2] - boxHalf_Z;
#else
      dx = P[i].Pos[0];
      dy = P[i].Pos[1];
      dz = P[i].Pos[2];
#endif

      r = sqrt(dx * dx + dy * dy + dz * dz);

      m = enclosed_mass(r);
#ifdef NFW_DARKFRACTION
      m *= NFW_DARKFRACTION;
#endif
      if(r > 0)
	{
	  P[i].GravAccel[0] += -All.G * m * dx / (r * r * r);
	  P[i].GravAccel[1] += -All.G * m * dy / (r * r * r);
	  P[i].GravAccel[2] += -All.G * m * dz / (r * r * r);

#ifdef DISTORTIONTENSORPS
	  double R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
	  double Rs = R200 / NFW_C;
	  double K = All.G * NFW_M200 / (Rs * (log(1 + NFW_C) - NFW_C / (1 + NFW_C)));
	  double r_red = r / Rs;
	  double x, y, z;

	  x = P[i].Pos[0];
	  y = P[i].Pos[1];
	  z = P[i].Pos[2];

	  P[i].tidal_tensorps[0][0] +=
	    -(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (1 / r - x * x / (r * r * r)) -
	      K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
		   2.0 * Rs * log(1 + r_red) / (r * r * r)) * x * x / (r * r));
	  P[i].tidal_tensorps[0][1] +=
	    -(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (0 - x * y / (r * r * r)) -
	      K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
		   2.0 * Rs * log(1 + r_red) / (r * r * r)) * x * y / (r * r));
	  P[i].tidal_tensorps[0][2] +=
	    -(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (0 - x * z / (r * r * r)) -
	      K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
		   2.0 * Rs * log(1 + r_red) / (r * r * r)) * x * z / (r * r));
	  P[i].tidal_tensorps[1][1] +=
	    -(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (1 / r - y * y / (r * r * r)) -
	      K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
		   2.0 * Rs * log(1 + r_red) / (r * r * r)) * y * y / (r * r));
	  P[i].tidal_tensorps[1][2] +=
	    -(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (0 - y * z / (r * r * r)) -
	      K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
		   2.0 * Rs * log(1 + r_red) / (r * r * r)) * y * z / (r * r));
	  P[i].tidal_tensorps[2][2] +=
	    -(-K * (1.0 / (r * (1 + r_red)) - log(1 + r_red) / (r * r_red)) * (1 / r - z * z / (r * r * r)) -
	      K * (-2.0 / (r * r * (1 + r_red)) - 1.0 / (r * (1 + r_red) * (1 + r_red) * Rs) +
		   2.0 * Rs * log(1 + r_red) / (r * r * r)) * z * z / (r * r));

	  P[i].tidal_tensorps[1][0] += P[i].tidal_tensorps[0][1];
	  P[i].tidal_tensorps[2][0] += P[i].tidal_tensorps[0][2];
	  P[i].tidal_tensorps[2][1] += P[i].tidal_tensorps[1][2];
#endif

	}
    }
#endif

#ifdef STATICPLUMMER
  int l;
  double r;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);

      for(l = 0; l < 3; l++)
	P[i].GravAccel[l] += -P[i].Pos[l] / pow(r * r + 1, 1.5);

#ifdef DISTORTIONTENSORPS
      double x, y, z, r2, f, f2;

      x = P[i].Pos[0];
      y = P[i].Pos[1];
      z = P[i].Pos[2];

      r2 = r * r;;
      f = pow(r2 + 1, 1.5);
      f2 = pow(r2 + 1, 2.5);


      P[i].tidal_tensorps[0][0] += -1.0 / f + 3.0 * x * x / f2;
      P[i].tidal_tensorps[0][1] += -0.0 / f + 3.0 * x * y / f2;
      P[i].tidal_tensorps[0][2] += -0.0 / f + 3.0 * x * z / f2;
      P[i].tidal_tensorps[1][1] += -1.0 / f + 3.0 * y * y / f2;
      P[i].tidal_tensorps[1][2] += -0.0 / f + 3.0 * y * z / f2;
      P[i].tidal_tensorps[2][2] += -1.0 / f + 3.0 * z * z / f2;
      P[i].tidal_tensorps[1][0] += P[i].tidal_tensorps[0][1];
      P[i].tidal_tensorps[2][0] += P[i].tidal_tensorps[0][2];
      P[i].tidal_tensorps[2][1] += P[i].tidal_tensorps[1][2];
#endif
    }
#endif



#ifdef STATICHQ
  double r, m, a;
  int l;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);

      a = pow(All.G * HQ_M200 / (100 * All.Hubble * All.Hubble), 1.0 / 3) / HQ_C *
	sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));

      m = HQ_M200 * pow(r / (r + a), 2);
#ifdef HQ_DARKFRACTION
      m *= HQ_DARKFRACTION;
#endif
      if(r > 0)
	{
	  for(l = 0; l < 3; l++)
	    P[i].GravAccel[l] += -All.G * m * P[i].Pos[l] / (r * r * r);

#ifdef DISTORTIONTENSORPS
	  double x, y, z, r2, r3, f, f2, f3;

	  x = P[i].Pos[0];
	  y = P[i].Pos[1];
	  z = P[i].Pos[2];

	  r2 = r * r;
	  r3 = r * r2;
	  f = r + a;
	  f2 = f * f;
	  f3 = f2 * f;


	  P[i].tidal_tensorps[0][0] +=
	    All.G * (2.0 * HQ_M200 / (r2 * f3) * x * x + HQ_M200 / (r3 * f2) * x * x - HQ_M200 / (r * f2));
	  P[i].tidal_tensorps[0][1] +=
	    All.G * (2.0 * HQ_M200 / (r2 * f3) * x * y + HQ_M200 / (r3 * f2) * x * y);
	  P[i].tidal_tensorps[0][2] +=
	    All.G * (2.0 * HQ_M200 / (r2 * f3) * x * z + HQ_M200 / (r3 * f2) * x * z);
	  P[i].tidal_tensorps[1][1] +=
	    All.G * (2.0 * HQ_M200 / (r2 * f3) * y * y + HQ_M200 / (r3 * f2) * y * y - HQ_M200 / (r * f2));
	  P[i].tidal_tensorps[1][2] +=
	    All.G * (2.0 * HQ_M200 / (r2 * f3) * y * z + HQ_M200 / (r3 * f2) * y * z);
	  P[i].tidal_tensorps[2][2] +=
	    All.G * (2.0 * HQ_M200 / (r2 * f3) * z * z + HQ_M200 / (r3 * f2) * z * z - HQ_M200 / (r * f2));
	  P[i].tidal_tensorps[1][0] += P[i].tidal_tensorps[0][1];
	  P[i].tidal_tensorps[2][0] += P[i].tidal_tensorps[0][2];
	  P[i].tidal_tensorps[2][1] += P[i].tidal_tensorps[1][2];
#endif
	}
    }
#endif

#ifdef STATICBRANDT

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {

/* note there is no acceleration in z */
      P[i].GravAccel[2] = 0;
      P[i].GravAccel[0] -= OmegaR(i, BRANDT_MODE) * OmegaR(i, BRANDT_MODE) * (P[i].Pos[0] - LONG_X / 2.);
      P[i].GravAccel[1] -= OmegaR(i, BRANDT_MODE) * OmegaR(i, BRANDT_MODE) * (P[i].Pos[1] - LONG_Y / 2.);
    }

#endif

#ifdef KEPLER_DISK
  MyFloat x00 = 4.0, y00 = 4.0;	/* 2D center of orbit: the is hard-coded for the relevant test problem */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      // add the constant field
      MyFloat r = pow(pow(P[i].Pos[1] - y00, 2.) + pow(P[i].Pos[0] - x00, 2.), 0.5);
      if(r <= 0.35)
	{
	  P[i].GravAccel[0] +=
	    -(P[i].Pos[0] - x00) * pow(r / 0.35,
				       2) / pow(pow(P[i].Pos[1] - y00, 2.) + pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[0] +=
	    +(P[i].Pos[0] - x00) * (0.35 - r) / 0.35 / pow(pow(P[i].Pos[1] - y00, 2.) +
							   pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[1] +=
	    -(P[i].Pos[1] - y00) * pow(r / 0.35,
				       2) / pow(pow(P[i].Pos[1] - y00, 2.) + pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[1] +=
	    +(P[i].Pos[1] - y00) * (0.35 - r) / 0.35 / pow(pow(P[i].Pos[1] - y00, 2.) +
							   pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[2] += 0;
	}
      else if((r > 0.35) && (r < 2.1))
	{
	  P[i].GravAccel[0] +=
	    -(P[i].Pos[0] - x00) / pow(pow(P[i].Pos[1] - y00, 2.) + pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[1] +=
	    -(P[i].Pos[1] - y00) / pow(pow(P[i].Pos[1] - y00, 2.) + pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[2] += 0;
	}
      else if(r >= 2.1)
	{
	  P[i].GravAccel[0] +=
	    -(P[i].Pos[0] - x00) * (1 + (r - 2.1) / 0.1) / pow(pow(P[i].Pos[1] - y00, 2.) +
							       pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[1] +=
	    -(P[i].Pos[1] - y00) * (1 + (r - 2.1) / 0.1) / pow(pow(P[i].Pos[1] - y00, 2.) +
							       pow(P[i].Pos[0] - x00, 2.), 1.5);
	  P[i].GravAccel[2] += 0;
	}
    }
#endif


  /* Now the force computation is finished */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  /*  gather some diagnostic information */

  timetree = timetree1 + timetree2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;


  MPI_Reduce(&timetree, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&timetree, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&timetree1, &sumt1, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&timetree1, &maxt1, 1, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&timetree2, &sumt2, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&timetree2, &maxt2, 1, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&timewait, &sumwaitall, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&timecomm, &sumcommall, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&Costtotal, &sum_costtotal, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&Ewaldcount, &ewaldtot, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);

  sumup_longs(1, &n_exported, &n_exported);
  sumup_longs(1, &N_nodesinlist, &N_nodesinlist);


  All.TotNumOfForces += GlobNumForceUpdate;

  plb = (NumPart / ((double) All.TotNumPart)) * NTask;
  MPI_Reduce(&plb, &plb_max, 1, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&Numnodestree, &maxnumnodes, 1, MPI_INT, MPI_MAX, 0, MYMPI_COMM_WORLD);

  CPU_Step[CPU_TREEMISC] += timeall - (timetree + timewait + timecomm);
  CPU_Step[CPU_TREEWALK1] += timetree1;
  CPU_Step[CPU_TREEWALK2] += timetree2;
  CPU_Step[CPU_TREESEND] += timecommsumm1;
  CPU_Step[CPU_TREERECV] += timecommsumm2;
  CPU_Step[CPU_TREEWAIT1] += timewait1;
  CPU_Step[CPU_TREEWAIT2] += timewait2;



#ifdef FIXEDTIMEINFIRSTPHASE
  MPI_Reduce(&min_time_first_phase, &min_time_first_phase_glob, 1, MPI_DOUBLE, MPI_MIN, 0, MYMPI_COMM_WORLD);
  VERBOSE(2, "FIXEDTIMEINFIRSTPHASE=%g  min_time_first_phase_glob=%g\n",
	  FIXEDTIMEINFIRSTPHASE, min_time_first_phase_glob);
#endif

  if(ThisTask == 0)
    {
      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g (%g) iter= %d\n",
	      (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      n_exported / ((double) GlobNumForceUpdate), N_nodesinlist / ((double) n_exported + 1.0e-10),
	      iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fprintf(FdTimings, "work-load balance: %g (%g %g) rel1to2=%g   max=%g avg=%g\n",
	      maxt / (1.0e-6 + sumt / NTask), maxt1 / (1.0e-6 + sumt1 / NTask),
	      maxt2 / (1.0e-6 + sumt2 / NTask), sumt1 / (1.0e-6 + sumt1 + sumt2), maxt, sumt / NTask);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart + NTopnodes));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", GlobNumForceUpdate / (sumt + 1.0e-20),
	      GlobNumForceUpdate / (1.0e-6 + maxt * NTask),
	      ((double) (sum_costtotal)) / (1.0e-20 + GlobNumForceUpdate),
	      ((double) ewaldtot) / (1.0e-20 + GlobNumForceUpdate));
      fprintf(FdTimings, "\n");

      fflush(FdTimings);
    }

  CPU_Step[CPU_TREEMISC] += measure_time();

  double costtotal_new = 0, sum_costtotal_new;
  if(TakeLevel >= 0)
    {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:costtotal_new)
#endif
      for(i = 0; i < NumPart; i++)
	costtotal_new += P[i].GravCost[TakeLevel];
      MPI_Reduce(&costtotal_new, &sum_costtotal_new, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
      if(ThisTask == 0)
	if(sum_costtotal)
	  VERBOSE_ALL(1, "relative error in the total number of tree-gravity interactions = %g\n",
		      (sum_costtotal - sum_costtotal_new) / sum_costtotal);
      /* can be non-zero if THREAD_SAFE_COSTS is not used (and due to round-off errors). */
    }

}


void *gravity_primary_loop(void *p)
{
  int i, j, ret;
  int thread_id = *(int *) p;

  int *exportflag, *exportnodecount, *exportindex;

  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

#ifdef FIXEDTIMEINFIRSTPHASE
  int counter = 0;
  double tstart;

  if(thread_id == 0)
    {
      tstart = second();
    }
#endif

  while(1)
    {
      int exitFlag = 0;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
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

      if(exitFlag)
	break;

#if !defined(PMGRID)
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
      if(Ewald_iter)
	{
#ifdef ACC_GRAVITY
	  ret = force_treeevaluate_ewald_correction(i, 0, exportflag, exportnodecount, exportindex, 1);
#else
#ifdef ACC
	  ret = force_treeevaluate_ewald_correction(i, 0, exportflag, exportnodecount, exportindex, 1);
#else
	  ret = force_treeevaluate_ewald_correction(i, 0, exportflag, exportnodecount, exportindex);
#endif
#endif
	  if(ret >= 0)
	    {

#ifdef _OPENMP
#pragma omp atomic
#endif
	      Ewaldcount += ret;	/* note: ewaldcount may be slightly incorrect for multiple threads if buffer gets filled up */
	    }
	  else
	    break;		/* export buffer has filled up */
	}
      else
#endif
	{
#ifdef ACC_GRAVITY
	  ret = force_treeevaluate(i, 0, exportflag, exportnodecount, exportindex, 1);
#else
#ifdef ACC
	  ret = force_treeevaluate(i, 0, exportflag, exportnodecount, exportindex, 1);
#else
	  ret = force_treeevaluate(i, 0, exportflag, exportnodecount, exportindex);
#endif
#endif
	  if(ret < 0)
	    break;		/* export buffer has filled up */

#ifdef _OPENMP
#pragma omp atomic
#endif
	  Costtotal += ret;
	}
#else

#ifdef NEUTRINOS
      if((P[i].Type != 2) || (All.Time > All.Time_tree_on_nu))
#endif
	{
#ifdef ACC_PROTO
	  ret = force_treeevaluate_shortrange(i, 0, exportflag, exportnodecount, exportindex, 1);	//1= drift!
#else
	  ret = force_treeevaluate_shortrange(i, 0, exportflag, exportnodecount, exportindex);
#endif
	  if(ret < 0)
	    break;		/* export buffer has filled up */

#ifdef _OPENMP
#pragma omp atomic
#endif
	  Costtotal += ret;

	}

#endif

      ProcessedFlag[i] = 1;	/* particle successfully finished */

#ifdef FIXEDTIMEINFIRSTPHASE
      if(thread_id == 0)
	{
	  counter++;
	  if((counter & 255) == 0)
	    {
	      if(timediff(tstart, second()) > FIXEDTIMEINFIRSTPHASE)
		{
		  TimerFlag = 1;
		  break;
		}
	    }
	}
      else
	{
	  if(TimerFlag)
	    break;
	}
#endif
    }

  return NULL;
}


void *gravity_secondary_loop(void *p)
{
  int j, nodesinlist, dummy, ret;

  while(1)
    {
#ifdef _OPENMP
#pragma omp atomic capture
#endif
      j = NextJ++;

      if(j >= Nimport)
	break;

#if !defined(PMGRID)
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
      if(Ewald_iter)
	{
#ifdef ACC_GRAVITY
	  int cost = force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy, &dummy, 1);

#else
#ifdef ACC
	  int cost = force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy, &dummy, 1);
#else
	  int cost = force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy, &dummy);
#endif
#endif
#ifdef _OPENMP
#pragma omp atomic
#endif
	  Ewaldcount += cost;
	}
      else
#endif
	{
#ifdef ACC_GRAVITY
	  ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy, 1);

#else
#ifdef ACC
	  ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy, 1);
#else
	  ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy);
#endif
#endif
#ifdef _OPENMP
#pragma omp critical(__cost__)
#endif
	  {
	    N_nodesinlist += nodesinlist;
	    Costtotal += ret;
	  }
	}
#else
#ifdef ACC_PROTO
      ret = force_treeevaluate_shortrange(j, 1, &nodesinlist, &dummy, &dummy, 1);
#else
      ret = force_treeevaluate_shortrange(j, 1, &nodesinlist, &dummy, &dummy);
#endif
#ifdef _OPENMP
#pragma omp critical(__cost__)
#endif
      {
	N_nodesinlist += nodesinlist;
	Costtotal += ret;
      }
#endif
    }

  return NULL;
}



void sum_top_level_node_costfactors(void)
{
  int i;

  double *costlist = (double *) mymalloc("costlist", NTopnodes * sizeof(double));
  double *costlist_all = (double *) mymalloc("costlist_all", NTopnodes * sizeof(double));

  for(i = 0; i < NTopnodes; i++)
    costlist[i] = Nodes[All.MaxPart + i].GravCost;

  MPI_Allreduce(costlist, costlist_all, NTopnodes, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);

  for(i = 0; i < NTopnodes; i++)
    Nodes[All.MaxPart + i].GravCost = costlist_all[i];

  myfree(costlist_all);
  myfree(costlist);
}





/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
	All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
	All.SofteningTable[0] = All.SofteningGas;

      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
	All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
	All.SofteningTable[1] = All.SofteningHalo;

      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
	All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
	All.SofteningTable[2] = All.SofteningDisk;

      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
	All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
	All.SofteningTable[3] = All.SofteningBulge;

      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
	All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
	All.SofteningTable[4] = All.SofteningStars;

      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
	All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
	All.SofteningTable[5] = All.SofteningBndry;
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];

}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
    */
int data_index_compare(const void *a, const void *b)
{
  if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task))
    return -1;

  if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task))
    return +1;

  if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index))
    return -1;

  if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index))
    return +1;

  if(((struct data_index *) a)->IndexGet < (((struct data_index *) b)->IndexGet))
    return -1;

  if(((struct data_index *) a)->IndexGet > (((struct data_index *) b)->IndexGet))
    return +1;

  return 0;
}

static void msort_dataindex_with_tmp(struct data_index *b, size_t n, struct data_index *t)
{
  struct data_index *tmp;
  struct data_index *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_dataindex_with_tmp(b1, n1, t);
  msort_dataindex_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(b1->Task < b2->Task || (b1->Task == b2->Task && b1->Index <= b2->Index))
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(struct data_index));

  memcpy(b, t, (n - n2) * sizeof(struct data_index));
}

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
  const size_t size = n * s;

  struct data_index *tmp = (struct data_index *) mymalloc("struct data_index *tmp", size);

  msort_dataindex_with_tmp((struct data_index *) b, n, tmp);

  myfree(tmp);
}



#ifdef ACC_GRAVITY


//gcc wants reduction variables to be public if omp region is opened outside function
 /*inline */ void *gravity_secondary_loop_vect(void *p, int should_drift)
{
  int nodesinlist, dummy, ret;


  double costTotal_CPU_II = 0.;
#pragma omp for reduction(+:costTotal_CPU_II) reduction(+:Ewaldcount) reduction(+:N_nodesinlist)
  for(int j = 0; j < Nimport; j++)
    {
#if !defined(PMGRID)
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
      if(Ewald_iter)
	{
	  int cost = force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy, &dummy, should_drift);
	  Ewaldcount += cost;
	}
      else
#endif
	{
	  ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy, should_drift);
	  N_nodesinlist += nodesinlist;
	  costTotal_CPU_II += ret;
	}
#else
      ret = force_treeevaluate_shortrange(j, 1, &nodesinlist, &dummy, &dummy, should_drift);
      N_nodesinlist += nodesinlist;
      costTotal_CPU_II += ret;
#endif

    }				//end for

  //  VERBOSE(-1,"costTotal_CPU_II=%g;  Costtotal=%g (before update)",costTotal_CPU_II, Costtotal);
#pragma omp atomic
  Costtotal += costTotal_CPU_II;
  //  VERBOSE(-1,"costTotal_CPU_II=%g;  Costtotal=%g (after update)",costTotal_CPU_II, Costtotal);

  return NULL;
}				//end /*inline*/



/*inline*/ void *gravity_primary_loop_vect(void *p, int is_primary_on_gpu, int should_drift)
{
  int ret;
  int thread_id = *(int *) p;
  int *exportflag, *exportnodecount, *exportindex;

  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;


  double costTotal_CPU_I = 0.;
  int chunk_counter = 0;
  for(int j = 0; j < NTask; j++)
    exportflag[j] = -1;
  int mode = 0;
#ifdef ACC_GRAVITY
  if(is_primary_on_gpu)
    mode = 3;
#endif
#pragma omp for reduction(+:costTotal_CPU_I)
  for(int chunk_counter = 0; chunk_counter < NActivePart; chunk_counter++)
    {
      if(BufferFullFlag != 0)
	{
	  //wait. GCC doesn't want break in parallel regions
#ifdef ACC_PGI
	  break;
#endif
	}
      else
	{
	  int target_id = chunk_counter;
	  int target = ActiveParticleList[target_id];
	  int ret = 0;
	  if(ProcessedFlag[target] == 1)
	    continue;
	  if(ProcessedFlag[target] == 2)
	    continue;
	  ProcessedFlag[target] = 0;
#if !defined(PMGRID)
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
	  if(Ewald_iter)
	    {
	      ret =
		force_treeevaluate_ewald_correction(target, mode, exportflag, exportnodecount, exportindex,
						    should_drift);
	    }
	  else
#endif
	    {
	      ret = force_treeevaluate(target, mode, exportflag, exportnodecount, exportindex, should_drift);


	    }
#else
#ifdef NEUTRINOS
	  if((P[target].Type != 2) || (All.Time > All.Time_tree_on_nu))
#endif

	    ret =
	      force_treeevaluate_shortrange(target, mode, exportflag, exportnodecount, exportindex,
					    should_drift);
#endif


	  if(ret < 0)
	    {
	      BufferFullFlag = 1;
#ifdef ACC_PGI
	      break;
#endif
	    }
	  else
	    {
	      ProcessedFlag[target] = 1;


	      if(mode == 0)
		costTotal_CPU_I += ret;
	    }
	}

    }
  NextParticle = -1;
  if(mode == 0)
    {
      // only if complete I phase is in GPU then we add costtoal, otherwise it is accounted by GPU
#pragma acc atomic
      Costtotal += costTotal_CPU_I;
    }
  if(Nexport > All.BunchSize)
    Nexport = All.BunchSize;
  return NULL;
}



inline void openacc_gravity_join_results(struct gravdata_out *LitePOut, double *PGravCosts,
					 double *NodeGravCosts, int NodesSize)
{
#pragma omp parallel for
  for(int counter = 0; counter < NActivePart; counter++)
    {
      int i1 = 0;
      int i2 = 0;
      int i = ActiveParticleList[counter];
#ifdef NEUTRINOS
      if(!((P[i].Type != 2) || (All.Time > All.Time_tree_on_nu)))
	{
	  continue;
	}
#endif
      P[i].GravAccel[0] += LitePOut[counter].Acc[0];
      P[i].GravAccel[1] += LitePOut[counter].Acc[1];
      P[i].GravAccel[2] += LitePOut[counter].Acc[2];
#ifdef EVALPOTENTIAL
      P[i].Potential += LitePOut[counter].Potential;
#endif
#ifdef DISTORTIONTENSORPS
      for(i1 = 0; i1 < 3; i1++)
	for(i2 = 0; i2 < 3; i2++)
	  P[i].tidal_tensorps[i1][i2] += LitePOut[counter].tidal_tensorps[i1][i2];
#endif
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
      P[i].AGS_corr += LitePOut[counter].AGS_corr;
#endif
    }
  double PGravCostsTotal = 0.;
#pragma omp parallel for reduction(+:PGravCostsTotal)
  for(int counter = 0; counter < NumPart; counter++)
    {
      P[counter].GravCost[TakeLevel] += PGravCosts[counter];
      PGravCostsTotal += PGravCosts[counter];
    }
  Costtotal += costTotal_GPU_I;	//finally add I phase (probably async) gpu grav cost
#pragma omp parallel for
  for(int counter = 0; counter < NodesSize; counter++)
    Nodes_base[counter].GravCost += NodeGravCosts[counter];
}

#endif
