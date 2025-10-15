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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"



static struct aux_data
{
  MyIDType ID;
  MyIDType GrNr;
  int OriginTask;
  int OriginIndex;
  int FinalTask;
}
 *Aux;


static int compare_Aux_GrNr_ID(const void *a, const void *b)
{
  if(((struct aux_data *) a)->GrNr < (((struct aux_data *) b)->GrNr))
    return -1;

  if(((struct aux_data *) a)->GrNr > (((struct aux_data *) b)->GrNr))
    return +1;

  if(((struct aux_data *) a)->ID < (((struct aux_data *) b)->ID))
    return -1;

  if(((struct aux_data *) a)->ID > (((struct aux_data *) b)->ID))
    return +1;

  return 0;
}

static int compare_Aux_OriginTask_OriginIndex(const void *a, const void *b)
{
  if(((struct aux_data *) a)->OriginTask < (((struct aux_data *) b)->OriginTask))
    return -1;

  if(((struct aux_data *) a)->OriginTask > (((struct aux_data *) b)->OriginTask))
    return +1;

  if(((struct aux_data *) a)->OriginIndex < (((struct aux_data *) b)->OriginIndex))
    return -1;

  if(((struct aux_data *) a)->OriginIndex > (((struct aux_data *) b)->OriginIndex))
    return +1;

  return 0;
}



/*  The following function should do the same thing as a call of

       parallel_sort(P, NumPart, sizeof(struct particle_data), io_compare_P_GrNr_ID)

    but reducing the peak memory consumption through use of an auxiliary array.

*/

#if defined(ORDER_SNAPSHOTS_BY_ID) || defined(SUBFIND_READ_FOF) || defined(SUBFIND_RESHUFFLE_CATALOGUE)

void parallel_sort_special_P_GrNr_ID(void)
{
  int i, j, Nimport, ngrp, sendTask, recvTask;

  Aux = mymalloc("Aux", NumPart * sizeof(struct aux_data));

  for(i = 0; i < NumPart; i++)
    {
      Aux[i].GrNr = P[i].GrNr;
      Aux[i].ID = P[i].ID;
      Aux[i].OriginTask = ThisTask;
      Aux[i].OriginIndex = i;
    }

  parallel_sort(Aux, NumPart, sizeof(struct aux_data), compare_Aux_GrNr_ID);

  for(i = 0; i < NumPart; i++)
    Aux[i].FinalTask = ThisTask;

  parallel_sort(Aux, NumPart, sizeof(struct aux_data), compare_Aux_OriginTask_OriginIndex);


  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < NumPart; i++)
    Send_count[Aux[i].FinalTask]++;


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

  if(Nimport != NumPart)
    PANIC("Nimport != NumPart");


  struct particle_data *pbuf;

  pbuf = mymalloc("pbuf", NumPart * sizeof(struct particle_data));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < NumPart; i++)
    pbuf[Send_offset[Aux[i].FinalTask] + Send_count[Aux[i].FinalTask]++] = P[i];


  memcpy(&P[Recv_offset[ThisTask]], &pbuf[Send_offset[ThisTask]],
	 Send_count[ThisTask] * sizeof(struct particle_data));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the particles */
	      MPI_Sendrecv(&pbuf[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct particle_data), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &P[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct particle_data), MPI_BYTE,
			   recvTask, TAG_DENS_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(pbuf);


  qsort(P, NumPart, sizeof(struct particle_data), io_compare_P_GrNr_ID);


  myfree(Aux);
}

#endif
