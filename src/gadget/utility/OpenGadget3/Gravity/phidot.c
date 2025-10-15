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
/*
 * phidot.c
 *
 *  Created on: Aug 20, 2013
 *      Author: margarita
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>


#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"

#ifdef PHIDOT

void phidot(void)
{
  int i, save_NextParticle;
  int j, k, ret, recvTask;
  int ndone, ndone_flag, dummy;
  int ngrp, place, nexport, nimport;
  double fac;
  double tstart, tend;
  MPI_Status status;
  double r2;

  if(All.ComovingIntegrationOn)
    set_softenings();

  VERBOSE(3, "Start computation of phi dot for all gas particles...\n");


  tstart = second();

  /* allocate buffers to arrange communication */
  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct gravdata_in) + sizeof(struct phidotdata_out) +
					     sizemax(sizeof(struct gravdata_in),
						     sizeof(struct phidotdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  NextParticle = 0;		/* beginn with this index */

  do
    {
      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;


#ifdef _OPENMP
#pragma omp parallel
#endif
      {
	int mainthreadid = ThisThread;
	phidot_primary_loop(&mainthreadid);	/* do local particles and prepare export list */
      }



      if(BufferFullFlag)
	{
	  int last_nextparticle = NextParticle;

	  NextParticle = save_NextParticle;

	  while(NextParticle < N_gas)
	    {
	      if(NextParticle == last_nextparticle)
		break;

	      if(ProcessedFlag[NextParticle] != 1)
		break;

	      ProcessedFlag[NextParticle] = 2;

	      NextParticle++;
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

      for(j = 0; j < NTask; j++)
	Send_count[j] = 0;
      for(j = 0; j < Nexport; j++)
	Send_count[DataIndexTable[j].Task]++;


      MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);


      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);


      for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      GravDataGet = (struct gravdata_in *) mymalloc("GravDataGet", Nimport * sizeof(struct gravdata_in));
      GravDataIn = (struct gravdata_in *) mymalloc("GravDataIn", Nexport * sizeof(struct gravdata_in));

      /* prepare particle data for export */

      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    GravDataIn[j].Pos[k] = P[place].Pos[k];

	  GravDataIn[j].Type = P[place].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(P[place].Type == 0)
	    GravDataIn[j].Soft = P[place].Hsml;
#endif

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


      myfree(GravDataIn);
      PhiDotDataResult =
	(struct phidotdata_out *) mymalloc("PhiDotDataResult", Nimport * sizeof(struct phidotdata_out));
      PhiDotDataOut =
	(struct phidotdata_out *) mymalloc("PhiDotDataOut", Nexport * sizeof(struct phidotdata_out));


      /* now do the particles that were sent to us */

      NextJ = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
	int mainthreadid = ThisThread;
	phidot_secondary_loop(&mainthreadid);
      }


      if(NextParticle >= N_gas)
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
		  MPI_Sendrecv(&PhiDotDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct phidotdata_out),
			       MPI_BYTE, recvTask, TAG_GRAV_B,
			       &PhiDotDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct phidotdata_out),
			       MPI_BYTE, recvTask, TAG_GRAV_B, MYMPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the results to the local particles */
      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  SphP[place].PhiDot += PhiDotDataOut[j].PhiDot;
	}

      myfree(PhiDotDataOut);
      myfree(PhiDotDataResult);
      myfree(GravDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);

  /* multiply with the gravitational constant */

  for(i = 0; i < N_gas; i++)
    SphP[i].PhiDot *= All.G;

  tend = second();

  VERBOSE(4, "phi dot calculation done, took %g seconds\n", timediff(tstart, tend));



}

void *phidot_primary_loop(void *p)
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


  while(1)
    {
      int exitFlag = 0;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      if(BufferFullFlag != 0 || NextParticle >= N_gas)
	{
	  exitFlag = 1;
	}
      else
	{
	  i = NextParticle;
	  ProcessedFlag[i] = 0;
	  NextParticle++;
	}

      if(exitFlag)
	break;

      if(P[i].Type == 0)
	{
	  ret = phidot_treeevaluate(i, 0, exportflag, exportnodecount, exportindex);
	  if(ret < 0)
	    break;		/* export buffer has filled up */
	}
      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;
}


void *phidot_secondary_loop(void *p)
{
  int j, nodesinlist, dummy;

  while(1)
    {
#ifdef _OPENMP
#pragma omp atomic capture
#endif
      j = NextJ++;

      if(j >= Nimport)
	break;

      phidot_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy);

    }

  return NULL;
}

#endif
