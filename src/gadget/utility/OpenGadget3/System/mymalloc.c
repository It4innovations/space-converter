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

#ifndef GADGET3_IO_LIB
#include <gsl/gsl_math.h>
#endif

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"


#define MAXBLOCKS 500
#define MAXCHARS  16

static size_t TotBytes;
static void *Base;

static unsigned long Nblocks;

static void **Table;
static size_t *BlockSize;
static char *MovableFlag;
static void ***BasePointers;

static char *VarName;
static char *FunctionName;
static char *FileName;
static int *LineNumber;


void mymalloc_init()
{
  size_t n;

  BlockSize = (size_t *) malloc(MAXBLOCKS * sizeof(size_t));
  Table = (void **) malloc(MAXBLOCKS * sizeof(void *));
  MovableFlag = (char *) malloc(MAXBLOCKS * sizeof(char));
  BasePointers = (void ***) malloc(MAXBLOCKS * sizeof(void **));
  VarName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FunctionName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FileName = (char *) malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  LineNumber = (int *) malloc(MAXBLOCKS * sizeof(int));

  memset(VarName, 0, MAXBLOCKS * MAXCHARS);
  memset(FunctionName, 0, MAXBLOCKS * MAXCHARS);
  memset(FileName, 0, MAXBLOCKS * MAXCHARS);

#ifdef KD_FFTW_MEMORY_BALANCE
  n = All.MaxMemSize * ((size_t) 1024 * 1024);
  printf("Task %d: dofftw=%d\n", ThisTask, dofftw);
#else
  n = All.MaxMemSize * ((size_t) 1024 * 1024);
#endif

#ifndef KD_ITTEREATE_MEMORYUSAGE

  if(!(Base = malloc(n)))
    {
      printf("Failed to allocate memory for `Base' (%d Mbytes).\n", All.MaxMemSize);
      endrun(122);
    }

#else

  n -= MAXBLOCKS * (sizeof(size_t) + sizeof(void *) + sizeof(char) + sizeof(void **) +
		    3 * MAXCHARS * sizeof(char) + sizeof(int));

  while(!(Base = malloc(n)) && n > All.MaxMemSize * ((size_t) 1024 * 1024 / 2))
    {
      printf("Task %d: Failed to allocate memory for `Base' (%d Mbytes).\n", ThisTask,
	     (int) (n / ((size_t) 1024 * 1024)));
      n = n - ((size_t) 50 * 1024 * 1024);
    }

  if(n <= All.MaxMemSize * ((size_t) 1024 * 1024 / 2))
    {
      printf("Task %d: Gave up to allocate memory for `Base' (%d Mbytes).\n", ThisTask,
	     (int) (n / ((size_t) 1024 * 1024)));
      endrun(122);
    }
#endif

  TotBytes = FreeBytes = n;

  AllocatedBytes = 0;
  Nblocks = 0;
  HighMarkBytes = 0;
}

void report_detailed_memory_usage_of_largest_task(size_t * OldHighMarkBytes, const char *label,
						  const char *func, const char *file, int line)
{
  size_t *sizelist, maxsize, minsize;
  double avgsize;
  int i, task;

  sizelist = (size_t *) mymalloc("sizelist", NTask * sizeof(size_t));
  MPI_Allgather(&AllocatedBytes, sizeof(size_t), MPI_BYTE, sizelist, sizeof(size_t), MPI_BYTE,
		MYMPI_COMM_WORLD);

  for(i = 1, task = 0, maxsize = minsize = sizelist[0], avgsize = sizelist[0]; i < NTask; i++)
    {
      if(sizelist[i] > maxsize)
	{
	  maxsize = sizelist[i];
	  task = i;
	}
      if(sizelist[i] < minsize)
	{
	  minsize = sizelist[i];
	}
      avgsize += sizelist[i];
    }

  myfree(sizelist);

  if(maxsize > 1.1 * (*OldHighMarkBytes))
    {
      *OldHighMarkBytes = maxsize;

      avgsize /= NTask;

      if(ThisTask == task)
	{
	  printf
	    ("\nAt '%s', %s()/%s/%d: Largest Allocation = %g Mbyte (on task=%d), Smallest = %g Mbyte, Average = %g Mbyte\n\n",
	     label, func, file, line, maxsize / (1024.0 * 1024.0), task, minsize / (1024.0 * 1024.0),
	     avgsize / (1024.0 * 1024.0));
	  dump_memory_table();
	}
      fflush(stdout);
    }
  MPI_Barrier(MYMPI_COMM_WORLD);

}




void dump_memory_table(void)
{
  unsigned int i;
  size_t totBlocksize = 0;

  printf("------------------------ Allocated Memory Blocks----------------------------------------\n");
  printf("Task   Nr F          Variable      MBytes   Cumulative         Function/File/Linenumber\n");
  printf("----------------------------------------------------------------------------------------\n");
  for(i = 0; i < Nblocks; i++)
    {
      totBlocksize += BlockSize[i];

      printf("%4d %4d %d  %16s  %10.4f   %10.4f  %s()/%s/%d\n",
	     ThisTask, i, MovableFlag[i], VarName + i * MAXCHARS, BlockSize[i] / (1024.0 * 1024.0),
	     totBlocksize / (1024.0 * 1024.0), FunctionName + i * MAXCHARS,
	     FileName + i * MAXCHARS, LineNumber[i]);
    }
  printf("----------------------------------------------------------------------------------------\n");
}

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line)
{
  //VERBOSE(-1,"%s %d %s %s %d", varname, n, func, file, line);
  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks >= MAXBLOCKS)
    {
      printf("Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n", ThisTask,
	     func, file, line, MAXBLOCKS);
      endrun(813);
    }

  if(n > FreeBytes)
    {
      dump_memory_table();
      printf
	("\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB).\n",
	 ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));
      endrun(812);
    }
  Table[Nblocks] = (char *) Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;
  MovableFlag[Nblocks] = 0;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}


void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file,
				int line)
{
  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks >= MAXBLOCKS)
    {
      printf("Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n", ThisTask,
	     func, file, line, MAXBLOCKS);
      endrun(816);
    }

  if(n > FreeBytes)
    {
      dump_memory_table();
      printf
	("\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB).\n",
	 ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));
      endrun(817);
    }
  Table[Nblocks] = (char *) Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks] = n;
  MovableFlag[Nblocks] = 1;
  BasePointers[Nblocks] = (void **) ptr;

  Nblocks += 1;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}



void myfree_fullinfo(void *p, const char *func, const char *file, int line)
{
  //VERBOSE(-1, "%s %s %d", func, file , line);
  if(Nblocks == 0)
    endrun(76878);

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();
      printf("Task=%d: Wrong call of myfree() at %s()/%s/line %d: not the last allocated block!\n", ThisTask,
	     func, file, line);
      fflush(stdout);
      endrun(814);
    }

  Nblocks -= 1;
  AllocatedBytes -= BlockSize[Nblocks];
  FreeBytes += BlockSize[Nblocks];
}



void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line)
{
  unsigned int i;

  if(Nblocks == 0)
    endrun(768728);

  /* first, let's find the block */
  unsigned int nr;

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();
      printf
	("Task=%d: Wrong call of myfree_movable() from %s()/%s/line %d - this block has not been allocated!\n",
	 ThisTask, func, file, line);
      fflush(stdout);
      endrun(8152);
    }

  if(nr < Nblocks - 1)		/* the block is not the last allocated block */
    {
      /* check that all subsequent blocks are actually movable */
      for(i = nr + 1; i < Nblocks; i++)
	if(MovableFlag[i] == 0)
	  {
	    dump_memory_table();
	    printf
	      ("Task=%d: Wrong call of myfree_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable allocated blocks\n",
	       ThisTask, func, file, line, nr);
	    fflush(stdout);
	    endrun(81252);
	  }
    }


  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  size_t offset = -BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += (size_t) BlockSize[i];

  if(nr < Nblocks - 1)
    memmove((char *) Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] = (char *) Table[i] + offset;
      *BasePointers[i] = (char *) *BasePointers[i] + offset;
    }

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i - 1] = Table[i];
      BasePointers[i - 1] = BasePointers[i];
      BlockSize[i - 1] = BlockSize[i];
      MovableFlag[i - 1] = MovableFlag[i];

      strncpy(VarName + (i - 1) * MAXCHARS, VarName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FunctionName + (i - 1) * MAXCHARS, FunctionName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FileName + (i - 1) * MAXCHARS, FileName + i * MAXCHARS, MAXCHARS - 1);
      LineNumber[i - 1] = LineNumber[i];
    }

  Nblocks -= 1;
}






void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks == 0)
    endrun(76879);

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();
      printf("Task=%d: Wrong call of myrealloc() at %s()/%s/line %d - not the last allocated block!\n",
	     ThisTask, func, file, line);
      fflush(stdout);
      endrun(815);
    }

  AllocatedBytes -= BlockSize[Nblocks - 1];
  FreeBytes += BlockSize[Nblocks - 1];

  if(n > FreeBytes)
    {
      dump_memory_table();
      printf
	("Task=%d: Not enough memory in myremalloc(n=%g MB) at %s()/%s/line %d. previous=%g FreeBytes=%g MB\n",
	 ThisTask, n / (1024.0 * 1024.0), func, file, line, BlockSize[Nblocks - 1] / (1024.0 * 1024.0),
	 FreeBytes / (1024.0 * 1024.0));
      endrun(812);
    }
  Table[Nblocks - 1] = (char *) Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  AllocatedBytes += n;
  BlockSize[Nblocks - 1] = n;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[Nblocks - 1];
}

void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  unsigned int i;

  if((n % 64) > 0)
    n = (n / 64 + 1) * 64;

  if(n < 64)
    n = 64;

  if(Nblocks == 0)
    endrun(768799);

  /* first, let's find the block */
  unsigned int nr;

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();
      printf
	("Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line %d - this block has not been allocated!\n",
	 ThisTask, func, file, line);
      fflush(stdout);
      endrun(8151);
    }

  if(nr < Nblocks - 1)		/* the block is not the last allocated block */
    {
      /* check that all subsequent blocks are actually movable */
      for(i = nr + 1; i < Nblocks; i++)
	if(MovableFlag[i] == 0)
	  {
	    dump_memory_table();
	    printf
	      ("Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable allocated blocks\n",
	       ThisTask, func, file, line, nr);
	    fflush(stdout);
	    endrun(8152);
	  }
    }


  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  if(n > FreeBytes)
    {
      dump_memory_table();
      printf
	("Task=%d: at %s()/%s/line %d: Not enough memory in myremalloc_movable(n=%g MB). previous=%g FreeBytes=%g MB\n",
	 ThisTask, func, file, line, n / (1024.0 * 1024.0), BlockSize[nr] / (1024.0 * 1024.0),
	 FreeBytes / (1024.0 * 1024.0));
      endrun(812);
    }

  size_t offset = n - BlockSize[nr];
  size_t length = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove((char *) Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] = (char *) Table[i] + offset;

      *BasePointers[i] = (char *) *BasePointers[i] + offset;
    }

  FreeBytes -= n;
  AllocatedBytes += n;
  BlockSize[nr] = n;

  if(AllocatedBytes > HighMarkBytes)
    HighMarkBytes = AllocatedBytes;

  return Table[nr];
}
