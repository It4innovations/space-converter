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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#include "allvars.h"
#include "proto.h"


/*  This function aborts the simulations. If a single processors
 *  wants an immediate termination,  the function needs to be
 *  called with ierr>0. A bunch of MPI-error messages will also
 *  appear in this case.
 *  For ierr=0, MPI is gracefully cleaned up, but this requires
 *  that all processors call endrun().
 */
#ifdef ACC_RUN
#include "../OpenACC/gpuallvars_branch.h"
#endif

#ifndef LT_STELLAREVOLUTION
void endrun(int ierr)
{
#ifdef _DEBUG
    assert(0);
#endif

#ifdef ACC_RUN
  if(acc_all.entered_gpu)
    update_openacc_structures_exit();
#endif
  if(ierr)
    {
      printf("task %d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
      fflush(stdout);
      MPI_Abort(MYMPI_COMM_WORLD, ierr);
      exit(0);
    }

  MPI_Finalize();
  exit(0);
}
#else

void EndRun(int ierr, const char *func, const char *file, const int line)
{
#ifdef _DEBUG
    assert(0);
#endif

#ifdef ACC_RUN
  if(acc_all.entered_gpu)
    update_openacc_structures_exit();
#endif

  if(ierr)
    {
      printf("task %d: endrun called with an error level of %d from func %s in file %s at line %d\n\n\n",
	     ThisTask, ierr, func, file, line);
      fflush(stdout);
      MPI_Abort(MYMPI_COMM_WORLD, ierr);
      exit(0);
    }

  MPI_Finalize();
  exit(0);
}



#endif

void EndRunString(const int do_output, const char *func, const char *file, const int line, const char *f_,
		  ...)
{
#ifdef ACC_RUN
  if(acc_all.entered_gpu)
    update_openacc_structures_exit();
#endif

  if(do_output)
    {
      va_list args;

      va_start(args, f_);
      printf("\n\n====== Gadget Panic. ThisTask: %d File: %s Function: %s Line:%d. Last words: ", ThisTask,
	     file, func, line);
      vprintf(f_, args);
      va_end(args);
      printf("\n\n");
      fflush(stdout);
#ifdef PANIC_PERPETUAL
      printf("\n\n====== Looping perpetually in panic\n");
      while(1)
	;
#endif
    }
  MPI_Abort(MYMPI_COMM_WORLD, 1);
  MPI_Finalize();
  exit(1);			//exit 1 is more right

}
