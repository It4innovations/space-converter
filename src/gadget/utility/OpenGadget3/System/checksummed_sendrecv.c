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

#ifdef MPISENDRECV_CHECKSUM

#undef MPI_Sendrecv


int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		       int dest, int sendtag, void *recvbufreal, int recvcount,
		       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status)
{
  int checksumtag = 1000, errtag = 2000;
  int i, iter = 0, err_flag, err_flag_imported, size_sendtype, size_recvtype, Local_ThisTask, Local_NTask;
  long long sendCheckSum, recvCheckSum, importedCheckSum;
  unsigned char *p, *buf, *recvbuf;

  if(dest != source)
    endrun(3);

  MPI_Comm_rank(comm, &Local_ThisTask);
  MPI_Comm_size(comm, &Local_NTask);

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  if(dest == Local_ThisTask)
    {
      memcpy(recvbufreal, sendbuf, recvcount * size_recvtype);
      return 0;
    }


  if(!(buf = (unsigned char*)mymalloc("buf", recvcount * size_recvtype + 1024)))
    endrun(6);

  for(i = 0, p = buf; i < recvcount * size_recvtype + 1024; i++)
    *p++ = 255;

  recvbuf = buf + 512;

  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
	       recvbuf, recvcount, recvtype, source, recvtag, comm, status);

  for(i = 0, p = buf; i < 512; i++, p++)
    {
      if(*p != 255)
	{
	  printf
	    ("MPI-ERROR: Task=%d/%s: Recv occured before recv buffer. message-size=%d from %d, i=%d c=%d\n",
	     Local_ThisTask, getenv("HOST"), recvcount, dest, i, *p);
	  fflush(stdout);
	  endrun(6);
	}
    }

  for(i = 0, p = recvbuf + recvcount * size_recvtype; i < 512; i++, p++)
    {
      if(*p != 255)
	{
	  printf
	    ("MPI-ERROR: Task=%d/%s: Recv occured after recv buffer. message-size=%d from %d, i=%d c=%d\n",
	     Local_ThisTask, getenv("HOST"), recvcount, dest, i, *p);
	  fflush(stdout);
	  endrun(6);
	}
    }


  for(i = 0, p = (unsigned char*)sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
    sendCheckSum += *p;

  importedCheckSum = 0;

  if(dest > Local_ThisTask)
    {
      if(sendcount > 0)
	MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
      if(recvcount > 0)
	MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, comm, status);
    }
  else
    {
      if(recvcount > 0)
	MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, comm, status);
      if(sendcount > 0)
	MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
    }

  checksumtag++;

  for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount * size_recvtype; i++, p++)
    recvCheckSum += *p;


  err_flag = err_flag_imported = 0;

  if(recvCheckSum != importedCheckSum)
    {
      printf
	("MPI-ERROR: Receive error on task=%d/%s from task=%d, message size=%d, sendcount=%d checksums= %d %d  %d %d. Try to fix it...\n",
	 Local_ThisTask, getenv("HOST"), source, recvcount, sendcount, (int) (recvCheckSum >> 32),
	 (int) recvCheckSum, (int) (importedCheckSum >> 32), (int) importedCheckSum);
      fflush(stdout);

      err_flag = 1;
    }

  if(dest > Local_ThisTask)
    {
      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
    }
  else
    {
      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
    }
  errtag++;

  if(err_flag > 0 || err_flag_imported > 0)
    {
      printf("Task=%d is on %s, wants to send %d and has checksum=%d %d of send data\n",
	     Local_ThisTask, getenv("HOST"), sendcount, (int) (sendCheckSum >> 32), (int) sendCheckSum);
      fflush(stdout);

      do
	{
	  sendtag++;
	  recvtag++;

	  for(i = 0, p = (unsigned char*)recvbuf; i < recvcount * size_recvtype; i++, p++)
	    *p = 0;

	  if((iter & 1) == 0)
	    {
	      if(dest > Local_ThisTask)
		{
		  if(sendcount > 0)
		    MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, comm);
		  if(recvcount > 0)
		    MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, comm, status);
		}
	      else
		{
		  if(recvcount > 0)
		    MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, comm, status);
		  if(sendcount > 0)
		    MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, comm);
		}
	    }
	  else
	    {
	      if(iter > 5)
		{
		  printf("we're trying to send each byte now on task=%d (iter=%d)\n", Local_ThisTask, iter);
		  fflush(stdout);
		  if(dest > Local_ThisTask)
		    {
		      for(i = 0, p = (unsigned char*)sendbuf; i < sendcount * size_sendtype; i++, p++)
			MPI_Ssend(p, 1, MPI_BYTE, dest, i, comm);
		      for(i = 0, p = (unsigned char*)recvbuf; i < recvcount * size_recvtype; i++, p++)
			MPI_Recv(p, 1, MPI_BYTE, dest, i, comm, status);
		    }
		  else
		    {
		      for(i = 0, p = (unsigned char*)recvbuf; i < recvcount * size_recvtype; i++, p++)
			MPI_Recv(p, 1, MPI_BYTE, dest, i, comm, status);
		      for(i = 0, p = (unsigned char*)sendbuf; i < sendcount * size_sendtype; i++, p++)
			MPI_Ssend(p, 1, MPI_BYTE, dest, i, comm);
		    }
		}
	      else
		{
		  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
			       recvbuf, recvcount, recvtype, source, recvtag, comm, status);
		}
	    }

	  importedCheckSum = 0;

	  for(i = 0, p = (unsigned char*)sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
	    sendCheckSum += *p;

	  printf("Task=%d gas send_checksum=%d %d\n", Local_ThisTask, (int) (sendCheckSum >> 32),
		 (int) sendCheckSum);
	  fflush(stdout);

	  if(dest > Local_ThisTask)
	    {
	      if(sendcount > 0)
		MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
	      if(recvcount > 0)
		MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag,
			 comm, status);
	    }
	  else
	    {
	      if(recvcount > 0)
		MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag,
			 comm, status);
	      if(sendcount > 0)
		MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, comm);
	    }

	  for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount; i++, p++)
	    recvCheckSum += *p;

	  err_flag = err_flag_imported = 0;

	  if(recvCheckSum != importedCheckSum)
	    {
	      printf
		("MPI-ERROR: Again (iter=%d) a receive error on task=%d/%s from task=%d, message size=%d, checksums= %d %d  %d %d. Try to fix it...\n",
		 iter, Local_ThisTask, getenv("HOST"), source, recvcount, (int) (recvCheckSum >> 32),
		 (int) recvCheckSum, (int) (importedCheckSum >> 32), (int) importedCheckSum);
	      fflush(stdout);
	      err_flag = 1;
	    }

	  if(dest > Local_ThisTask)
	    {
	      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
	      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
	    }
	  else
	    {
	      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, comm, status);
	      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, comm);
	    }

	  if(err_flag == 0 && err_flag_imported == 0)
	    break;

	  errtag++;
	  checksumtag++;
	  iter++;
	}
      while(iter < 10);

      if(iter >= 10)
	{
	  char buf[1000];
	  int length;
	  FILE *fd;

	  sprintf(buf, "send_data_%d.dat", Local_ThisTask);
	  fd = fopen(buf, "w");
	  length = sendcount * size_sendtype;
	  fwrite(&length, 1, sizeof(int), fd);
	  fwrite(sendbuf, sendcount, size_sendtype, fd);
	  fclose(fd);

	  sprintf(buf, "recv_data_%d.dat", Local_ThisTask);
	  fd = fopen(buf, "w");
	  length = recvcount * size_recvtype;
	  fwrite(&length, 1, sizeof(int), fd);
	  fwrite(recvbuf, recvcount, size_recvtype, fd);
	  fclose(fd);

	  printf("MPI-ERROR: Even 10 trials proved to be insufficient on task=%d/%s. Stopping\n",
		 Local_ThisTask, getenv("HOST"));
	  fflush(stdout);
	  endrun(10);
	}
    }

  memcpy(recvbufreal, recvbuf, recvcount * size_recvtype);

  myfree(buf);

  return 0;
}

#endif
