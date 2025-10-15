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
#include <unistd.h>

#ifdef IMPOSE_PINNING

#define __USE_GNU

#include <sched.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"

#ifndef SOCKETS
#define SOCKETS   4		/* Setting for our 4 x AMD Magny-Cours nodes */
#define MAX_CORES 48		/* Setting for our 4 x AMD Magny-Cours nodes */
#endif


 /* SOCKETS   = 12 for 2 x Intel 6-core nodes with hypethreading  */
 /* MAX_CORES = 24 for 2 x Intel 6-core nodes with hypethreading  */



static cpu_set_t cpuset;


void get_core_set(void)
{
  CPU_ZERO(&cpuset);
  sched_getaffinity(getpid(), sizeof(cpuset), &cpuset);
}


void pin_to_core_set(void)
{
  int core, task, num_threads;
  static cpu_set_t cpuset_new;

  char *p = getenv("OMP_NUM_THREADS");
  if(p)
    num_threads = atoi(p);
  else
    num_threads = 1;



  CPU_ZERO(&cpuset_new);

  int corestart = 0;

  for(task = 0, core = corestart; task < NTask * num_threads; task++)
    {
      while(!CPU_ISSET(core, &cpuset))
	{
	  core += SOCKETS;
	  if(core >= MAX_CORES)
	    {
	      corestart++;
	      if(corestart >= SOCKETS)
		corestart = 0;

	      core = corestart;
	    }
	}

      if((task / num_threads) == ThisTask)
	CPU_SET(core, &cpuset_new);

      core += SOCKETS;
      if(core >= MAX_CORES)
	{
	  corestart++;
	  if(corestart >= SOCKETS)
	    corestart = 0;

	  core = corestart;
	}
    }

  sched_setaffinity(getpid(), sizeof(cpuset_new), &cpuset_new);
}


void report_pinning(void)
{
  cpu_set_t cpuset;
  int i;
  char buf[MAX_CORES + 1];

  CPU_ZERO(&cpuset);
  sched_getaffinity(getpid(), sizeof(cpuset), &cpuset);

  for(i = 0; i < MAX_CORES; i++)
    if(CPU_ISSET(i, &cpuset))
      buf[i] = '1';
    else
      buf[i] = '-';
  buf[MAX_CORES] = 0;

  for(i = 0; i < NTask; i++)
    {
      if(ThisTask == i)
	printf("Task=%02d: %s\n", ThisTask, buf);
      fflush(stdout);
      MPI_Barrier(MYMPI_COMM_WORLD);
    }
}

#endif
