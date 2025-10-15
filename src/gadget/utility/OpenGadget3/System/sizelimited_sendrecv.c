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

#ifdef MPISENDRECV_SIZELIMIT


#undef MPI_Sendrecv


int MPI_Sizelimited_Sendrecv(void *sendbuf, size_t sendcount, MPI_Datatype sendtype,
			     int dest, int sendtag, void *recvbuf, size_t recvcount,
			     MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
			     MPI_Status * status)
{
  int iter = 0, size_sendtype, size_recvtype, send_now, recv_now;
  int count_limit;


  if(dest != source)
    endrun(3);

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  if(dest == ThisTask)
    {
      memcpy(recvbuf, sendbuf, recvcount * size_recvtype);
      return 0;
    }

  count_limit = (int) ((((long long) MPISENDRECV_SIZELIMIT) * 1024 * 1024) / size_sendtype);

  while(sendcount > 0 || recvcount > 0)
    {
      if(sendcount > count_limit)
	{
	  send_now = count_limit;
	  if(iter == 0)
	    {
	      printf("imposing size limit on MPI_Sendrecv() on task=%d (send of size=%d)\n",
		     ThisTask, sendcount * size_sendtype);
	      fflush(stdout);
	    }
	  iter++;
	}
      else
	send_now = sendcount;

      if(recvcount > count_limit)
	recv_now = count_limit;
      else
	recv_now = recvcount;

      MPI_Sendrecv(sendbuf, send_now, sendtype, dest, sendtag,
		   recvbuf, recv_now, recvtype, source, recvtag, comm, status);

      sendcount -= send_now;
      recvcount -= recv_now;

      sendbuf += send_now * size_sendtype;
      recvbuf += recv_now * size_recvtype;
    }

  return 0;
}

#endif
