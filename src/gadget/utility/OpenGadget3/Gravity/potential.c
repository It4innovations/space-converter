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


#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"


/*! \file potential.c
 *  \brief Computation of the gravitational potential of particles
 */

#if !defined(EVALPOTENTIAL) && (defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL))

/*! This function computes the gravitational potential for ALL the particles.
 *  First, the (short-range) tree potential is computed, and then, if needed,
 *  the long range PM potential is added.
 */
void compute_potential(void)
{
  int i;

#ifndef NOGRAVITY
  int j, k, ret, recvTask;
  int ndone, ndone_flag, dummy;
  int ngrp, place, nexport, nimport;
  double fac;
  MPI_Status status;
  double r2;

  if(All.ComovingIntegrationOn)
    set_softenings();

  VERBOSE(2 "Start computation of potential for all particles...\n");


  CPU_Step[CPU_MISC] += measure_time();


  if(TreeReconstructFlag)
    {
      VERBOSE(2, "Tree construction.\n");

      CPU_Step[CPU_MISC] += measure_time();

#if defined(SFR) || defined(BLACK_HOLES)
      rearrange_particle_sequence();
#endif

      force_treebuild(NumPart, NULL);

      CPU_Step[CPU_TREEBUILD] += measure_time();

      TreeReconstructFlag = 0;

      VERBOSE(2, "Tree construction done.\n");
    }


  /* allocate buffers to arrange communication */
  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct gravdata_in) + sizeof(struct potdata_out) +
					     sizemax(sizeof(struct gravdata_in),
						     sizeof(struct potdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

#ifndef SUBFIND_RESHUFFLE_AND_POTENTIAL
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_current != All.Ti_Current)
      drift_particle(i, All.Ti_Current);
#endif
  i = 0;			/* beginn with this index */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      for(nexport = 0; i < NumPart; i++)
	{
#ifndef PMGRID
	  ret = force_treeevaluate_potential(i, 0, &nexport, Send_count);
#else
	  ret = force_treeevaluate_potential_shortrange(i, 0, &nexport, Send_count);
#endif
	  if(ret < 0)
	    break;		/* export buffer has filled up */
	}

      MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      GravDataGet = (struct gravdata_in *) mymalloc("GravDataGet", nimport * sizeof(struct gravdata_in));
      GravDataIn = (struct gravdata_in *) mymalloc("GravDataIn", nexport * sizeof(struct gravdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    GravDataIn[j].Pos[k] = P[place].Pos[k];

	  GravDataIn[j].Type = P[place].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(P[place].Type == 0)
	    GravDataIn[j].Soft = SphP[place].Hsml;
#endif
#ifdef ADAPTGRAVSOFT
	  GravDataIn[j].AGS_Hsml = P[place].AGS_Hsml;
#endif
	  GravDataIn[j].OldAcc = P[place].OldAcc;

	  for(k = 0; k < NODELISTLENGTH; k++)
	    GravDataIn[j].NodeList[k] = DataNodeList[DataIndexTable[j].IndexGet].NodeList[k];
	}


      /* exchange particle data */

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
			       recvTask, TAG_POTENTIAL_A,
			       &GravDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
			       recvTask, TAG_POTENTIAL_A, MYMPI_COMM_WORLD, &status);
		}
	    }
	}

      myfree(GravDataIn);
      PotDataResult = (struct potdata_out *) mymalloc("PotDataResult", nimport * sizeof(struct potdata_out));
      PotDataOut = (struct potdata_out *) mymalloc("PotDataOut", nexport * sizeof(struct potdata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	{
#ifndef PMGRID
	  force_treeevaluate_potential(j, 1, &dummy, &dummy);
#else
	  force_treeevaluate_potential_shortrange(j, 1, &dummy, &dummy);
#endif
	}

      if(i >= NumPart)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&PotDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct potdata_out),
			       MPI_BYTE, recvTask, TAG_POTENTIAL_B,
			       &PotDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct potdata_out),
			       MPI_BYTE, recvTask, TAG_POTENTIAL_B, MYMPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the results to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  //P[place].p.dPotential += PotDataOut[j].Potential;
	  P[place].Potential += PotDataOut[j].Potential;
	}

      myfree(PotDataOut);
      myfree(PotDataResult);
      myfree(GravDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);

#ifndef ADAPTGRAVSOFT
  /* add correction to exclude self-potential */

  for(i = 0; i < NumPart; i++)
    {
      /* remove self-potential */
      P[i].Potential += P[i].Mass / All.SofteningTable[P[i].Type];

      if(All.ComovingIntegrationOn)
	if(All.PeriodicBoundariesOn)
	  P[i].Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
	    pow(All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G), 1.0 / 3);
    }
#endif

  /* multiply with the gravitational constant */

  for(i = 0; i < NumPart; i++)
    P[i].Potential *= All.G;


#ifdef PMGRID

#ifdef PERIODIC
  pmpotential_periodic();
#ifdef PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(1);	/* try again */
    }
  if(i == 1)
    endrun(88686);
#endif
#else
  i = pmpotential_nonperiodic(0);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(0);	/* try again */
    }
  if(i == 1)
    endrun(88687);
#ifdef PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();

      i = pmpotential_nonperiodic(1);
    }
  if(i != 0)
    endrun(88688);
#endif
#endif

#endif



  if(All.ComovingIntegrationOn)
    {
#ifndef PERIODIC
      fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	{
	  for(k = 0, r2 = 0; k < 3; k++)
	    r2 += P[i].Pos[k] * P[i].Pos[k];

	  P[i].Potential += fac * r2;
	}
#endif
    }
  else
    {
      fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;
      if(fac != 0)
	{
	  for(i = 0; i < NumPart; i++)
	    {
	      for(k = 0, r2 = 0; k < 3; k++)
		r2 += P[i].Pos[k] * P[i].Pos[k];

	      P[i].Potential += fac * r2;
	    }
	}
    }


  VERBOSE(2, "potential done.\n");



#else
  for(i = 0; i < NumPart; i++)
    P[i].Potential = 0;
#endif

  CPU_Step[CPU_POTENTIAL] += measure_time();
}


#endif
