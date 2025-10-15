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
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "utilities.h"
#include "proto.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

void myprintf(const char *format, ...)
{
#ifdef PARALLEL
  printf("Task %03d: ", ThisTask);
#endif
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

char *util_fgets(char *str, int num, FILE * stream, char *file, int line)
{
  char *ret = fgets(str, num, stream);
  if(ret == NULL)
    {
      printf("error: fgets in file %s at line %d\n", file, line);
      endrun(200);
    }

  return ret;
}

size_t util_fread(void *ptr, size_t size, size_t count, FILE * stream, char *file, int line)
{
  size_t result = fread(ptr, size, count, stream);
  if(result != count)
    {
      printf("error: fread in file %s at line %d\n", file, line);
      endrun(201);
    }

  return result;
}



//ideally we can use this MPI_alltoall instead of copying the below bunch of #ifdef every time.
void MPI_Alltoall_hybrid(int Send_count, int one, int MPI_INT, int Recv_count, one_2, int MPI_INT, int MYMPI_COMM_WORLD){

#ifdef IMPORT_ALLtoALLv_FROM_G4
  return MPI_Alltoall(Send_count, one, MPI_INT, Recv_count, one_2, MPI_INT, MYMPI_COMM_WORLD);
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
      MPI_Win_create(Recv_count, NTask * sizeof(MPI_INT), sizeof(MPI_INT), MPI_INFO_NULL,
		     MYMPI_COMM_WORLD, &win);
      MPI_Win_fence(0, win);
      for(j = 0; j < NTask; ++j)
	{
	  if(Send_count[j] > 0)
	    MPI_Put(&Send_count[j], 1, MPI_INT, j, ThisTask, 1, MPI_INT, win);
	}
      MPI_Win_fence(0, win);
      MPI_Win_free(&win);
      return 0;
#ifdef IMPORT_ALLtoALLv_FROM_G4_WHEN_TIMEBIN_ACTIVE_AFTER
    }
  else
    {
      return MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);
    }
#endif //ACTIVE_AFTER
#endif
  
}

