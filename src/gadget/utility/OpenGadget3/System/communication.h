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

#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"


/*  memcpy(DensDataIn[j].NodeList,
	 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
*/

static inline void copy_data_nodelist(){
  for(int j = 0; j < Nexport; j++){

    memcpy(&DataNodeListIn[j],
	   &DataNodeList[DataIndexTable[j].IndexGet], sizeof(struct data_nodelist));


  }
}

template <typename data_in> 
static inline   void   templated_exchange_data_in_get(data_in *DataIn,data_in *DataGet, int exchange_datanodelist){
  int status, recvTask;
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)    {
    recvTask = ThisTask ^ ngrp;
    if(recvTask < NTask){
      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)	    {
	/* get the particles */
	MPI_Sendrecv(&DataIn[Send_offset[recvTask]],
		     Send_count[recvTask] * sizeof(data_in), MPI_BYTE,
		     recvTask, TAG_HYDRO_A,
		     &DataGet[Recv_offset[recvTask]],
		     Recv_count[recvTask] * sizeof(data_in), MPI_BYTE,
		     recvTask, TAG_HYDRO_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if(exchange_datanodelist){


	   MPI_Sendrecv(&DataNodeListIn[Send_offset[recvTask]],
		     Send_count[recvTask] * sizeof(struct data_nodelist), MPI_BYTE,  
		     recvTask, TAG_GRAV_A,
		     &DataNodeListGet[Recv_offset[recvTask]],
		     Recv_count[recvTask] * sizeof(struct data_nodelist), MPI_BYTE,
		     recvTask, TAG_GRAV_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	

      }
    }
  }

}


static inline int prepare_sendrec_offset_send(const int save_nextparticle, const int panic_on_error)
{
  long long n_exported = 0;
  if(BufferFullFlag)
    {
      int last_nextparticle = NextParticle;
      
      NextParticle = save_nextparticle;

      while(NextParticle >= 0)
	{
	  if(NextParticle == last_nextparticle)
	    break;
	  
	  if(ProcessedFlag[NextParticle] != 1)
	    break;

	  ProcessedFlag[NextParticle] = 2;

	  NextParticle = NextActiveParticle[NextParticle];
	}

      if(NextParticle == save_nextparticle)
	{
	  /* in this case, the buffer is too small to process even a single particle */
	  printf("Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n", ThisTask, P[NextParticle].Type,
		 P[NextParticle].Pos[0], P[NextParticle].Pos[1], P[NextParticle].Pos[2],
		 P[NextParticle].Mass);
	  if(P[NextParticle].Type == 0)
	    printf("   rho=%g hsml=%g NgB= %g/%d NextParticle=%d/%d/%d\n", SphP[NextParticle].Density, P[NextParticle].Hsml,P[NextParticle].NumNgb,P[NextParticle].TrueNGB,NextParticle,N_gas,NumPart);
	  if(panic_on_error)
	    PANIC("  the buffer is too small to process even a single particle ");
	  else 
	    return 0;
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
  
  MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

#ifndef IMPORT_ALLtoALLv_FROM_G4
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);
#else
#ifdef IMPORT_ALLtoALLv_FROM_G4_WHEN_TIMEBIN_ACTIVE_AFTER
  int non_active_timebin = 0;
  for(int j = TIMEBINS - 1; j >= 0; j--)
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
      VERBOSE(-1,"Sendreceive performed with MPI windows");
#endif
      for(int j = 0; j < NTask; j++)
	Recv_count[j] = 0;
      MPI_Win_create(Recv_count, NTask * sizeof(int), sizeof(int), MPI_INFO_NULL,
		     MYMPI_COMM_WORLD, &win);
      MPI_Win_fence(0, win);
      for(int i = 0; i < NTask; ++i)
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

  Nimport = 0;
  Recv_offset[0] = 0;
  Send_offset[0] = 0;
  for(int j = 0; j < NTask; j++)
    {
      Nimport += Recv_count[j];
      
      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }
  return 1;
}

static inline int prepare_sendrec_offset_send_conduction(const int save_nextparticle, const int panic_on_error)
{
  long long n_exported = 0;
  if(BufferFullFlag)
    {
      int last_nextparticle = NextParticle;
      
      NextParticle = save_nextparticle;

      while(NextParticle < N_gas)
	{
	  if(NextParticle == last_nextparticle)
	    break;
	  
	  if(ProcessedFlag[NextParticle] != 1)
	    break;

	  ProcessedFlag[NextParticle] = 2;

	  NextParticle++;
	}

      if(NextParticle == save_nextparticle)
	{
	  /* in this case, the buffer is too small to process even a single particle */
	  printf("Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n", ThisTask, P[NextParticle].Type,
		 P[NextParticle].Pos[0], P[NextParticle].Pos[1], P[NextParticle].Pos[2],
		 P[NextParticle].Mass);
	  if(P[NextParticle].Type == 0)
	    printf("   rho=%g hsml=%g NextParticle=%d/%d/%d\n", SphP[NextParticle].Density, P[NextParticle].Hsml,NextParticle,N_gas,NumPart);
	  if(panic_on_error)
            PANIC("  the buffer is too small to process even a single particle ");
          else
            return 0;


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
  
  MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

#ifndef IMPORT_ALLtoALLv_FROM_G4
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);
#else
#ifdef IMPORT_ALLtoALLv_FROM_G4_WHEN_TIMEBIN_ACTIVE_AFTER
  int non_active_timebin = 0;
  for(int j = TIMEBINS - 1; j >= 0; j--)
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
      for(int j = 0; j < NTask; j++)
	Recv_count[j] = 0;
      MPI_Win_create(Recv_count, NTask * sizeof(int), sizeof(int), MPI_INFO_NULL,
		     MYMPI_COMM_WORLD, &win);
      MPI_Win_fence(0, win);
      for(int i = 0; i < NTask; ++i)
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

  Nimport = 0;
  Recv_offset[0] = 0;
  Send_offset[0] = 0;
  for(int j = 0; j < NTask; j++)
    {
      Nimport += Recv_count[j];
      
      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }
  return 1;
}


#endif
