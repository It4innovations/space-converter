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
#include <time.h>

#ifndef _WIN32
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#endif

#ifndef GADGET3_IO_LIB
#include <gsl/gsl_rng.h>
#endif

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"



#ifdef DEBUG
#include <fenv.h>
void enable_core_dumps_and_fpu_exceptions(void)
{
  struct rlimit rlim;
  extern int feenableexcept(int __excepts);

  /* enable floating point exceptions */

  /*
     feenableexcept(FE_DIVBYZERO | FE_INVALID);
   */

  /* Note: FPU exceptions appear not to work properly
   * when the Intel C-Compiler for Linux is used
   */

  /* set core-dump size to infinity */
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  /* MPICH catches the signales SIGSEGV, SIGBUS, and SIGFPE....
   * The following statements reset things to the default handlers,
   * which will generate a core file.
   */
  /*
     signal(SIGSEGV, catch_fatal);
     signal(SIGBUS, catch_fatal);
     signal(SIGFPE, catch_fatal);
     signal(SIGINT, catch_fatal);
   */

  signal(SIGSEGV, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGFPE, SIG_DFL);
  signal(SIGINT, SIG_DFL);

  /* Establish a handler for SIGABRT signals. */
  signal(SIGABRT, catch_abort);
}


void catch_abort(int sig)
{
  MPI_Finalize();
  exit(0);
}

void catch_fatal(int sig)
{
  terminate_processes();
  MPI_Finalize();

  signal(sig, SIG_DFL);
  raise(sig);
}


void terminate_processes(void)
{
  pid_t my_pid;
  char buf[500], hostname[500], *cp;
  char commandbuf[500];
  FILE *fd;
  int i, pid;

  sprintf(buf, "%s%s", All.OutputDir, "PIDs.txt");

  my_pid = getpid();

  if((fd = fopen(buf, "r")))
    {
      for(i = 0; i < NTask; i++)
	{
	  int ret;

	  ret = fscanf(fd, "%s %d", hostname, &pid);

	  cp = hostname;
	  while(*cp)
	    {
	      if(*cp == '.')
		*cp = 0;
	      else
		cp++;
	    }

	  if(my_pid != pid)
	    {
	      sprintf(commandbuf, "ssh %s kill -ABRT %d", hostname, pid);
	      printf("--> %s\n", commandbuf);
	      fflush(stdout);
#ifndef NOCALLSOFSYSTEM
	      ret = system(commandbuf);
#endif
	    }
	}

      fclose(fd);
    }
}

void write_pid_file(void)
{
  pid_t my_pid;
  char mode[8], buf[500];
  FILE *fd;
  int i;

  my_pid = getpid();

  sprintf(buf, "%s%s", All.OutputDir, "PIDs.txt");

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  for(i = 0; i < NTask; i++)
    {
      if(ThisTask == i)
	{
	  if(ThisTask == 0)
	    sprintf(mode, "w");
	  else
	    sprintf(mode, "a");

	  if((fd = fopen(buf, mode)))
	    {
	      fprintf(fd, "%s %d\n", getenv("HOST"), (int) my_pid);
	      fclose(fd);
	    }
	}

      MPI_Barrier(MYMPI_COMM_WORLD);
    }
}
#endif

#ifdef PAUSE_RUN_TO_ATTACH_DEBUGGER
void pause_run_to_attach_debugger()
{
  int continue_run = 0, my_pid;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  my_pid = getpid();
#if ( PAUSE_RUN_TO_ATTACH_DEBUGGER < 0 )
  if( ! getenv("GADGET_PAUSE_RUN_TO_ATTACH_DEBUGGER") ) {
    if(ThisTask == 0)
      printf("GADGET_PAUSE_RUN_TO_ATTACH_DEBUGGER environment variable unset: will not wait for debugger to attach.\n", my_pid, hostname);
    return;
  }
#endif
  printf("PID %d on %s ready for debugger to attach\n", my_pid, hostname);
  fflush(stdout);

  int *all_pids;
  if(ThisTask == 0)
    all_pids = (int *) malloc(sizeof(int) * NTask);

  MPI_Gather(&my_pid, 1, MPI_INT, all_pids, 1, MPI_INT, 0, MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      FILE *fd = fopen("pid_list_for_debugger.txt", "w");
      int j;
      for(j = 0; j < NTask; ++j)
	fprintf(fd, "%i \n", all_pids[j]);
      fclose(fd);
      free(all_pids);
      printf("PID file written\n");
      fflush(stdout);
    }

  while(continue_run == 0)	/* continue_run needsto be set to 1 ("set var continue_run = 1") by debugger to continue */
    sleep(1);
}
#endif

/*
double get_random_number(unsigned int id)
{
  return RndTable[(id % RNDTABLE)];
}
*/

double get_random_number(MyIDType id)
{
  return RndTable[(int) (id % RNDTABLE)];
}

#ifndef GADGET3_IO_LIB
void set_random_numbers(void)
{
  int i;

  for(i = 0; i < RNDTABLE; i++)
    RndTable[i] = gsl_rng_uniform(random_generator);
}
#endif

/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
#ifdef WALLCLOCK
  return MPI_Wtime();
#else
  return ((double) clock()) / CLOCKS_PER_SEC;
#endif

  /* note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}

double measure_time(void)	/* strategy: call this at end of functions to account for time in this function, and before another (nontrivial) function is called */
{
  double t, dt;

  t = second();
  dt = t - WallclockTime;
  WallclockTime = t;

  return dt;
}

double report_time(void)	/* strategy: call this to measure sub-times of functions */
{
  double t, dt;

  t = second();
  dt = t - WallclockTime;

  return dt;
}


/* returns the time difference between two measurements
 * obtained with second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)			/* overflow has occured (for systems with 32bit tick counter) */
    {
#ifdef WALLCLOCK
      dt = 0;
#else
      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}




#ifdef X86FIX

#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#define _FPU_EXTENDED 0x0300
#define _FPU_DOUBLE   0x0200

void x86_fix(void)
{
  unsigned short dummy, new_cw;
  unsigned short *old_cw;

  old_cw = &dummy;

  _FPU_GETCW(*old_cw);
  new_cw = (*old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
}

#endif


void minimum_large_ints(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
		MYMPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = src[j];

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      if(res[j] > numlist[i * n + j])
	res[j] = numlist[i * n + j];

  myfree(numlist);
}

void sumup_large_ints(int n, int *src, long long *res)
{
  int i, j, *numlist;

  numlist = (int *) mymalloc("numlist", NTask * n * sizeof(int));
  MPI_Allgather(src, n, MPI_INT, numlist, n, MPI_INT, MYMPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}

void sumup_longs(int n, long long *src, long long *res)
{
  int i, j;
  long long *numlist;

  numlist = (long long *) mymalloc("numlist", NTask * n * sizeof(long long));
  MPI_Allgather(src, n * sizeof(long long), MPI_BYTE, numlist, n * sizeof(long long), MPI_BYTE,
		MYMPI_COMM_WORLD);

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(i = 0; i < NTask; i++)
    for(j = 0; j < n; j++)
      res[j] += numlist[i * n + j];

  myfree(numlist);
}

size_t sizemax(size_t a, size_t b)
{
  if(a < b)
    return b;
  else
    return a;
}

#ifndef _WIN32
void report_VmRSS(void)
{
  pid_t my_pid;
  FILE *fd;
  char buf[1024];

  my_pid = getpid();

  sprintf(buf, "/proc/%d/status", my_pid);

  if((fd = fopen(buf, "r")))
    {
      while(1)
	{
	  if(fgets(buf, 500, fd) != buf)
	    break;

	  if(strncmp(buf, "VmRSS", 5) == 0)
	    {
	      printf("ThisTask=%d: %s", ThisTask, buf);
	    }
	  if(strncmp(buf, "VmSize", 6) == 0)
	    {
	      printf("ThisTask=%d: %s", ThisTask, buf);
	    }
	}
      fclose(fd);
    }
}
#endif

#ifndef GADGET3_IO_LIB
long long report_comittable_memory(long long *MemTotal,
				   long long *Committed_AS, long long *SwapTotal, long long *SwapFree)
{
  FILE *fd;
  char buf[1024];

  if((fd = fopen("/proc/meminfo", "r")))
    {
      while(1)
	{
	  if(fgets(buf, 500, fd) != buf)
	    break;

	  if(bcmp(buf, "MemTotal", 8) == 0)
	    {
	      *MemTotal = atoll(buf + 10);
	    }
	  if(strncmp(buf, "Committed_AS", 12) == 0)
	    {
	      *Committed_AS = atoll(buf + 14);
	    }
	  if(strncmp(buf, "SwapTotal", 9) == 0)
	    {
	      *SwapTotal = atoll(buf + 11);
	    }
	  if(strncmp(buf, "SwapFree", 8) == 0)
	    {
	      *SwapFree = atoll(buf + 10);
	    }
	}
      fclose(fd);
    }

  return (*MemTotal - *Committed_AS);
}
#endif

int FullestTask, EmptiestTask;

#ifndef GADGET3_IO_LIB
void mpi_report_comittable_memory(void)
{
  long long *sizelist, maxsize[6], minsize[6];
  double avgsize[6];
  int i, imem, mintask[6], maxtask[6];
  long long Mem[6];
  char label[512];

  Mem[0] = report_comittable_memory(&Mem[1], &Mem[2], &Mem[3], &Mem[4]);
  Mem[5] = Mem[1] - Mem[0];

  for(imem = 0; imem < 6; imem++)
    {
      sizelist = (long long *) malloc(NTask * sizeof(long long));
      MPI_Allgather(&Mem[imem], sizeof(long long), MPI_BYTE, sizelist, sizeof(long long), MPI_BYTE,
		    MPI_COMM_WORLD);

      for(i = 1, mintask[imem] = 0, maxtask[imem] = 0, maxsize[imem] = minsize[imem] =
	  sizelist[0], avgsize[imem] = sizelist[0]; i < NTask; i++)
	{
	  if(sizelist[i] > maxsize[imem])
	    {
	      maxsize[imem] = sizelist[i];
	      maxtask[imem] = i;
	    }
	  if(sizelist[i] < minsize[imem])
	    {
	      minsize[imem] = sizelist[i];
	      mintask[imem] = i;
	    }
	  avgsize[imem] += sizelist[i];
	}

      free(sizelist);
    }

  if(ThisTask == 0)
    {
      printf
	("___________________________________________________________________________________________________________________\n");
      for(imem = 0; imem < 6; imem++)
	{
	  switch (imem)
	    {
	    case 0:
	      sprintf(label, "AvailMem");
	      break;
	    case 1:
	      sprintf(label, "Total Mem");
	      break;
	    case 2:
	      sprintf(label, "Committed_AS");
	      break;
	    case 3:
	      sprintf(label, "SwapTotal");
	      break;
	    case 4:
	      sprintf(label, "SwapFree");
	      break;
	    case 5:
	      sprintf(label, "AllocMem");
	      break;
	    }
	  printf
	    ("%s:\t Largest = %10.2f Mb (on task=%d), Smallest = %10.2f Mb (on task=%d), Average = %10.2f Mb\n",
	     label, maxsize[imem] / (1024.0), maxtask[imem], minsize[imem] / (1024.0), mintask[imem],
	     avgsize[imem] / (1024.0 * NTask));
	}
      printf
	("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }


  FullestTask = maxtask[2];
  EmptiestTask = mintask[2];


  if(ThisTask == maxtask[2])
    {
      char name[MPI_MAX_PROCESSOR_NAME];
      int len;
      MPI_Get_processor_name(name, &len);
      printf("Task=%d has the maximum commited memory and is host: %s\n", ThisTask, name);
      printf
	("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }

  fflush(stdout);
}
#endif

void report_task_list(void)
{

  if(ThisTask == FullestTask)
    {
      char name[MPI_MAX_PROCESSOR_NAME];
      int len;
      MPI_Get_processor_name(name, &len);

      char command[10000];
      sprintf(command, "top -b -n 1 -o RES > processes_largest_memuse_%d_%s.txt", ThisTask, name);
      system(command);
      sprintf(command, "free > free_largest_memuse_%d_%s.txt", ThisTask, name);
      system(command);
      sprintf(command, " ps -eo pid,pmem,rss,vsize,comm --sort vsize > ps_largest_memuse_%d_%s.txt", ThisTask,
	      name);
      system(command);
    }



  if(ThisTask == EmptiestTask)
    {
      char name[MPI_MAX_PROCESSOR_NAME];
      int len;
      MPI_Get_processor_name(name, &len);

      char command[10000];
      sprintf(command, "top -b -n 1 -o RES > processes_smallest_memuse_%d_%s.txt", ThisTask, name);
      system(command);
      sprintf(command, "free > free_smallest_memuse_%d_%s.txt", ThisTask, name);
      system(command);
      sprintf(command, " ps -eo pid,pmem,rss,vsize,comm --sort vsize > ps_smallest_memuse_%d_%s.txt",
	      ThisTask, name);
      system(command);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  fflush(stdout);
}



void write_ps_files(void)
{
  long long pid;
  FILE *fd, *fd2, *ofd;
  char buf[1024], buf2[1024];

  char command[10000], filename[200], memfilename[200], outfile[200];

  if(ThisTask == FullestTask || ThisTask == EmptiestTask)
    {
      sprintf(filename, "ps_file_%i.txt", ThisTask);
      sprintf(command, "ps -eo pid > %s", filename);
      system(command);

      if(ThisTask == FullestTask)
	sprintf(outfile, "memory_largest_%i.txt", ThisTask);

      if(ThisTask == EmptiestTask)
	sprintf(outfile, "memory_smallest_%i.txt", ThisTask);

      if(!(ofd = fopen(outfile, "w")))
	endrun(4532);


      if((fd = fopen(filename, "r")))
	{
	  fgets(buf, 500, fd);

	  while(1)
	    {
	      if(fgets(buf, 500, fd) != buf)
		break;
	      pid = atoll(buf);


	      sprintf(memfilename, "/proc/%lld/status", pid);

	      if((fd2 = fopen(memfilename, "r")))
		{
		  while(1)
		    {
		      if(fgets(buf2, 500, fd2) != buf2)
			break;

		      if(strncmp(buf2, "VmRSS", 5) == 0)
			{
			  fprintf(ofd, "PID %lld: %s \n", pid, buf2);
			}
		      if(strncmp(buf2, "VmSize", 6) == 0)
			{
			  fprintf(ofd, "PID %lld: %s \n", pid, buf2);
			}
		    }
		  fclose(fd2);
		}

	    }
	  fclose(fd);
	}
      fclose(ofd);
    }
}
