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

#include "allvars.h"
#include "proto.h"

#include <gsl/gsl_heapsort.h>


#ifdef THRUST_SORT
#include "../OpenACC/thrust_sort.h"
#endif

static struct peano_hilbert_data
{
  peanokey key;
  int index;
}
 *__restrict__ mp;

static int *__restrict__ Id;

void peano_hilbert_order(void)
{
  int i;
#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
  double t0_kd, t1_kd;
  t0_kd = second();
#endif


  if(ThisTask == 0)
    printf("begin Peano-Hilbert order...\n");

  if(N_gas)
    {
      mp = (struct peano_hilbert_data *) mymalloc("mp", sizeof(struct peano_hilbert_data) * N_gas);
      Id = (int *) mymalloc("Id", sizeof(int) * N_gas);

      int N_gas_4 = (N_gas / 4) * 4;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int ii = 0; ii < N_gas_4; ii++)
	{
	  mp[ii].index = ii;
	  mp[ii].key = Key[ii];
	}
      for(int ii = N_gas_4; ii < N_gas; ii++)
	{
	  mp[ii].index = ii;
	  mp[ii].key = Key[ii];
	}

#ifdef MYSORT
#ifdef _OPENMP
      if(N_gas > OMP_SORT_THRESH)
	serial_sort_omp(mp, N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
      else
	mysort_peano(mp, N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
#else
      mysort_peano(mp, N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
#endif
#else
      STD_SORT(mp, N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
#endif

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:    sorting gas keys took %g sec\n", timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int ii = 0; ii < N_gas_4; ii++)
	Id[mp[ii].index] = ii;
      for(int ii = N_gas_4; ii < N_gas; ii++)
	Id[mp[ii].index] = ii;

#ifdef _OPENMP
      kd_reorder_gas2();
#else
      reorder_gas();
#endif

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:    reodering gas particles took %g sec\n", timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

      myfree(Id);
      myfree(mp);
    }


  if(NumPart - N_gas > 0)
    {
      mp =
	(struct peano_hilbert_data *) mymalloc("mp", sizeof(struct peano_hilbert_data) * (NumPart - N_gas));
      mp -= (N_gas);

      Id = (int *) mymalloc("Id", sizeof(int) * (NumPart - N_gas));
      Id -= (N_gas);

      int NumPart_4 = N_gas + ((NumPart - N_gas) / 4) * 4;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int ii = N_gas; ii < NumPart_4; ii += 4)
	{
	  mp[ii].index = ii;
	  mp[ii + 1].index = ii + 1;
	  mp[ii + 2].index = ii + 2;
	  mp[ii + 3].index = ii + 3;
	  mp[ii].key = Key[ii];
	  mp[ii + 1].key = Key[ii + 1];
	  mp[ii + 2].key = Key[ii + 2];
	  mp[ii + 3].key = Key[ii + 3];
	}
      for(int ii = NumPart_4; ii < NumPart; ii++)
	{
	  mp[ii].index = ii;
	  mp[ii].key = Key[ii];
	}


#ifdef THRUST_SORT

      int *indexes = (int *) mymalloc("indexes", NumPart * sizeof(int));
      unsigned long long *keys =
	(unsigned long long *) mymalloc("keys", NumPart * sizeof(unsigned long long));

#pragma omp parallel for
      for(int i = 0; i < NumPart; i++)
	{
	  indexes[i] = i;
	  keys[i] = Key[i];
	}

      thrust_device_keys_sort(keys, indexes, NumPart);

#pragma omp parallel for
      for(int i = 0; i < NumPart; i++)
	{
	  int j = indexes[i];
	  mp[i].index = j;
	  mp[i].key = Key[j];
	}




      myfree(keys);
      myfree(indexes);
#else //else THRUST_SORT
#ifdef MYSORT
#ifdef _OPENMP
      if(NumPart - N_gas > OMP_SORT_THRESH)
	serial_sort_omp(mp + N_gas, NumPart - N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
      else
	mysort_peano(mp + N_gas, NumPart - N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
#else
      mysort_peano(mp + N_gas, NumPart - N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
#endif
#else
      STD_SORT(mp + N_gas, NumPart - N_gas, sizeof(struct peano_hilbert_data), peano_compare_key);
#endif
#endif //end THRUST_SORT

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:    sorting non gas particles keys took %g sec\n", timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int ii = N_gas; ii < NumPart_4; ii += 4)
	{
	  Id[mp[ii].index] = ii;
	  Id[mp[ii + 1].index] = ii + 1;
	  Id[mp[ii + 2].index] = ii + 2;
	  Id[mp[ii + 3].index] = ii + 3;
	}
      for(int ii = NumPart_4; ii < NumPart; ii++)
	Id[mp[ii].index] = ii;

#ifdef _OPENMP
      kd_reorder_particles2();
#else
      reorder_particles();
#endif

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:    reodering non gas particles took %g sec\n", timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

      Id += N_gas;
      myfree(Id);
      mp += N_gas;
      myfree(mp);
    }

  if(ThisTask == 0)
    printf("Peano-Hilbert done.\n");
}


int peano_compare_key(const void *a, const void *b)
{
  if(((struct peano_hilbert_data *) a)->key < (((struct peano_hilbert_data *) b)->key))
    return -1;

  if(((struct peano_hilbert_data *) a)->key > (((struct peano_hilbert_data *) b)->key))
    return +1;

  return 0;
}

void reorder_gas(void)
{
  int i;
  struct particle_data Psave, Psource;
  struct sph_particle_data SphPsave, SphPsource;
  int idsource, idsave, dest;
#ifdef WRITE_KEY_FILES
  peanokey Ksource, Ksave;
#endif

  for(i = 0; i < N_gas; i++)
    {
      if(Id[i] != i)
	{
#ifdef WRITE_KEY_FILES
	  Ksource = Key[i];
#endif
	  Psource = P[i];
	  SphPsource = SphP[i];

	  idsource = Id[i];
	  dest = Id[i];

	  do
	    {
#ifdef WRITE_KEY_FILES
	      Ksave = Key[dest];
#endif
	      Psave = P[dest];
	      SphPsave = SphP[dest];
	      idsave = Id[dest];

#ifdef WRITE_KEY_FILES
	      Key[dest] = Ksource;
#endif
	      P[dest] = Psource;
	      SphP[dest] = SphPsource;
	      Id[dest] = idsource;

	      if(dest == i)
		break;

	      Psource = Psave;
	      SphPsource = SphPsave;
	      idsource = idsave;
#ifdef WRITE_KEY_FILES
	      Ksource = Ksave;
#endif

	      dest = idsource;
	    }
	  while(1);
	}
    }
}


void reorder_particles(void)
{
  int i;
  struct particle_data Psave, Psource;
  int idsource, idsave, dest;
#ifdef WRITE_KEY_FILES
  peanokey Ksource, Ksave;
#endif

  for(i = N_gas; i < NumPart; i++)
    {
      if(Id[i] != i)
	{
#ifdef WRITE_KEY_FILES
	  Ksource = Key[i];
#endif
	  Psource = P[i];
	  idsource = Id[i];

	  dest = Id[i];

	  do
	    {
#ifdef WRITE_KEY_FILES
	      Ksave = Key[dest];
#endif
	      Psave = P[dest];
	      idsave = Id[dest];

#ifdef WRITE_KEY_FILES
	      Key[dest] = Ksource;
#endif
	      P[dest] = Psource;
	      Id[dest] = idsource;
#ifdef LT_STELLAREVOLUTION
	      if(P[dest].Type == 4)
		MetP[P[dest].pt.MetID].PID = dest;
#endif
#if defined(BLACK_HOLES)
	      if(P[dest].Type == 5)
		BHP[P[dest].pt.BHID].PID = dest;
#endif

#ifdef AXION_DM
	      if(P[dest].Type == AXION_TYPE)
		APP(dest).PID = dest;
#endif

	      if(dest == i)
		break;

	      Psource = Psave;
	      idsource = idsave;
#ifdef WRITE_KEY_FILES
	      Ksource = Ksave;
#endif

	      dest = idsource;
	    }
	  while(1);
	}
    }
}

#ifdef _OPENMP
void kd_reorder_gas2(void)
{
  struct particle_data *__restrict__ Psave;
  struct sph_particle_data *__restrict__ SphPsave;
#ifdef WRITE_KEY_FILES
  peanokey *__restrict__ Ksave;
#endif
#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
  double t1_kd, t0_kd = second();
#endif

  if(FreeBytes < N_gas * (sizeof(struct particle_data) + sizeof(struct sph_particle_data) + sizeof(peanokey)))
    {
      reorder_gas();
    }
  else
    {
      Psave = (struct particle_data *) mymalloc("Psave", N_gas * sizeof(struct particle_data));
      SphPsave = (struct sph_particle_data *) mymalloc("Psave", N_gas * sizeof(struct sph_particle_data));
#ifdef WRITE_KEY_FILES
      Ksave = (peanokey *) mymalloc("Ksave", N_gas * sizeof(peanokey));
#endif

      if(ThisTask == 0)
	printf("Peano: using %g MB of remaining %g MB for non inline sorting of sph particles and keys\n",
	       N_gas * (sizeof(struct particle_data) + sizeof(struct sph_particle_data) +
			sizeof(peanokey)) / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));

      int N_gas_4 = (N_gas / 4) * 4;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i0 = 0; i0 < N_gas_4; i0 += 4)
	{
	  Psave[Id[i0]] = P[i0];
	  Psave[Id[i0 + 1]] = P[i0 + 1];
	  Psave[Id[i0 + 2]] = P[i0 + 2];
	  Psave[Id[i0 + 3]] = P[i0 + 3];

	  SphPsave[Id[i0]] = SphP[i0];
	  SphPsave[Id[i0 + 1]] = SphP[i0 + 1];
	  SphPsave[Id[i0 + 2]] = SphP[i0 + 2];
	  SphPsave[Id[i0 + 3]] = SphP[i0 + 3];

#ifdef WRITE_KEY_FILES
	  Ksave[Id[i0]] = Key[i0];
	  Ksave[Id[i0 + 1]] = Key[i0 + 1];
	  Ksave[Id[i0 + 2]] = Key[i0 + 2];
	  Ksave[Id[i0 + 3]] = Key[i0 + 3];
#endif
	}
      for(int i0 = N_gas_4; i0 < N_gas; i0++)
	{
	  Psave[Id[i0]] = P[i0];
	  SphPsave[Id[i0]] = SphP[i0];
#ifdef WRITE_KEY_FILES
	  Ksave[Id[i0]] = Key[i0];
#endif
	}

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:       reodering gas particles (stage 1) took %g sec\n",
	       timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

#ifdef _OPENMP
#pragma omp parallel for
      for(int i0 = 0; i0 < N_gas_4; i0 += 4)
	{
	  P[i0] = Psave[i0];
	  P[i0 + 1] = Psave[i0 + 1];
	  P[i0 + 2] = Psave[i0 + 2];
	  P[i0 + 3] = Psave[i0 + 3];

	  SphP[i0] = SphPsave[i0];
	  SphP[i0 + 1] = SphPsave[i0 + 1];
	  SphP[i0 + 2] = SphPsave[i0 + 2];
	  SphP[i0 + 3] = SphPsave[i0 + 3];

#ifdef WRITE_KEY_FILES
	  Key[i0] = Ksave[i0];
	  Key[i0 + 1] = Ksave[i0 + 1];
	  Key[i0 + 2] = Ksave[i0 + 2];
	  Key[i0 + 3] = Ksave[i0 + 3];
#endif
	}
      for(int i0 = N_gas_4; i0 < N_gas; i0++)
	{
	  P[i0] = Psave[i0];
	  SphP[i0] = SphPsave[i0];
#ifdef WRITE_KEY_FILES
	  Key[i0] = Ksave[i0];
#endif
	}

#else
      memcpy(P, Psave, N_gas * sizeof(struct particle_data));
      memcpy(SphP, SphPsave, N_gas * sizeof(struct sph_particle_data));
#ifdef WRITE_KEY_FILES
      memcpy(Key, Ksave, N_gas * sizeof(peanokey));
#endif
#endif // _OPENMP

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:       reodering gas particles (stage 2) took %g sec\n",
	       timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

#ifdef WRITE_KEY_FILES
      myfree(Ksave);
#endif
      myfree(SphPsave);
      myfree(Psave);
    }
}


void kd_reorder_gas2_alternative(void)
{
  struct particle_data *Psave;
  struct sph_particle_data *SphPsave;
  int *FlagFree;
  int *FlagCopied;
#ifdef WRITE_KEY_FILES
  peanokey *Ksave;
#endif
#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
  double t1_kd, t0_kd = second();
#endif

  if(FreeBytes <
     N_gas * (sizeof(struct particle_data) + sizeof(struct sph_particle_data) + sizeof(peanokey) +
	      2 * sizeof(int)))
    {
      reorder_gas();
    }
  else
    {
      Psave = (struct particle_data *) mymalloc("Psave", N_gas * sizeof(struct particle_data));
      SphPsave = (struct sph_particle_data *) mymalloc("Psave", N_gas * sizeof(struct sph_particle_data));
      FlagFree = (int *) mymalloc("FlagFree", N_gas * sizeof(int));
      FlagCopied = (int *) mymalloc("FlagCopied", N_gas * sizeof(int));
#ifdef WRITE_KEY_FILES
      Ksave = (peanokey *) mymalloc("Ksave", (NumPart - N_gas) * sizeof(peanokey));
#endif

      if(ThisTask == 0)
	printf("Peano: using %g MB of remaining %g MB for non inline sorting of sph particles and keys\n",
	       N_gas * (sizeof(struct particle_data) + sizeof(struct sph_particle_data) + sizeof(peanokey) +
			2 * sizeof(int)) / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));

      memset(FlagFree, 0, N_gas * sizeof(int));
      memset(FlagCopied, 0, N_gas * sizeof(int));

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:       reodering gas particles (zeoring) took %g sec\n",
	       timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i0 = 0; i0 < N_gas; i0++)
	{
	  if(i0 != Id[i0])
	    {
	      if(FlagFree[Id[i0]] == 1)
		{
		  if(FlagCopied[i0] == 1)
		    {
		      P[Id[i0]] = Psave[i0];
		      SphP[Id[i0]] = SphPsave[i0];
#ifdef WRITE_KEY_FILES
		      Key[Id[i0]] = Ksave[i0];
#endif
		    }
		  else
		    {
		      P[Id[i0]] = P[i0];
		      SphP[Id[i0]] = SphP[i0];
		      Key[Id[i0]] = Key[i0];
		    }
		  FlagFree[i0] == 1;
		}
	      else
		{
		  Psave[Id[i0]] = P[Id[i0]];
		  SphPsave[Id[i0]] = SphP[Id[i0]];
#ifdef WRITE_KEY_FILES
		  Ksave[Id[i0]] = Key[Id[i0]];
#endif
		  FlagCopied[Id[i0]] = 1;

		  if(FlagCopied[i0] == 1)
		    {
		      P[Id[i0]] = Psave[i0];
		      SphP[Id[i0]] = SphPsave[i0];
#ifdef WRITE_KEY_FILES
		      Key[Id[i0]] = Ksave[i0];
#endif
		    }
		  else
		    {
		      P[Id[i0]] = P[i0];
		      SphP[Id[i0]] = SphP[i0];
		      Key[Id[i0]] = Key[i0];
		    }
		  FlagFree[i0] == 1;
		}
	    }
	}

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:       reodering gas particles (stage 2) took %g sec\n",
	       timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

#ifdef WRITE_KEY_FILES
      myfree(Ksave);
#endif
      myfree(FlagCopied);
      myfree(FlagFree);

      myfree(SphPsave);
      myfree(Psave);
    }
}


void kd_reorder_particles2(void)
{
  struct particle_data *__restrict__ Psave;
#ifdef WRITE_KEY_FILES
  peanokey *__restrict__ Ksave;
#endif
#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
  double t1_kd, t0_kd = second();
#endif

  if(FreeBytes < (NumPart - N_gas) * (sizeof(struct particle_data) + sizeof(peanokey)))
    {
      reorder_particles();
    }
  else
    {
      Psave = (struct particle_data *) mymalloc("Psave", (NumPart - N_gas) * sizeof(struct particle_data));
#ifdef WRITE_KEY_FILES
      Ksave = (peanokey *) mymalloc("Ksave", (NumPart - N_gas) * sizeof(peanokey));
#endif

      if(ThisTask == 0)
	printf("Peano: using %g MB of remaining %g MB for non inline sorting of particles and keys\n",
	       (NumPart - N_gas) * (sizeof(struct particle_data) + sizeof(peanokey)) / (1024.0 * 1024.0),
	       FreeBytes / (1024.0 * 1024.0));

#pragma omp parallel for
      for(int i1 = N_gas; i1 < NumPart; i1++)
	{
	  int i0 = Id[i1] - N_gas;
	  Psave[i0] = P[i1];
#ifdef WRITE_KEY_FILES
	  Ksave[i0] = Key[i1];
#endif

#ifdef LT_STELLAREVOLUTION
	  if(Psave[i0].Type == 4)
	    MetP[Psave[i0].pt.MetID].PID = Id[i1];
#endif
#if defined(BLACK_HOLES)
	  if(Psave[i0].Type == 5)
	    BHP[Psave[i0].pt.BHID].PID = Id[i1];
#endif

#ifdef AXION_DM
	  if(Psave[i0].Type == AXION_TYPE)
	    APP(Id[i1]).PID = Id[i1];
#endif
	}

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:       reodering particles (stage 1) took %g sec\n",
	       timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif

#ifdef _OPENMP
#pragma omp parallel for
      for(int i0 = 0; i0 < NumPart - N_gas; i0++)
	{
	  P[i0 + N_gas] = Psave[i0];
#ifdef WRITE_KEY_FILES
	  Key[i0 + N_gas] = Ksave[i0];
#endif
	}
#else
      memcpy(&P[N_gas], Psave, (NumPart - N_gas) * sizeof(struct particle_data));
#ifdef WRITE_KEY_FILES
      memcpy(&Key[N_gas], Ksave, (NumPart - N_gas) * sizeof(peanokey));
#endif
#endif

#ifdef KD_EXTRA_TIMER_OUTPUT_PEANO
      t1_kd = second();
      if(ThisTask == 0)
	printf("EXTRA TIMER PEANO:       reodering particles (stage 2) took %g sec\n",
	       timediff(t0_kd, t1_kd));
      t0_kd = second();
#endif


#ifdef WRITE_KEY_FILES
      myfree(Ksave);
#endif
      myfree(Psave);

    }
}
#endif


/*  The following rewrite of the original function
 *  peano_hilbert_key_old() has been written by MARTIN REINECKE. 
 *  It is about a factor 2.3 - 2.5 faster than Volker's old routine!
 */
const unsigned char rottable3[48][8] = {
  {36, 28, 25, 27, 10, 10, 25, 27},
  {29, 11, 24, 24, 37, 11, 26, 26},
  {8, 8, 25, 27, 30, 38, 25, 27},
  {9, 39, 24, 24, 9, 31, 26, 26},
  {40, 24, 44, 32, 40, 6, 44, 6},
  {25, 7, 33, 7, 41, 41, 45, 45},
  {4, 42, 4, 46, 26, 42, 34, 46},
  {43, 43, 47, 47, 5, 27, 5, 35},
  {33, 35, 36, 28, 33, 35, 2, 2},
  {32, 32, 29, 3, 34, 34, 37, 3},
  {33, 35, 0, 0, 33, 35, 30, 38},
  {32, 32, 1, 39, 34, 34, 1, 31},
  {24, 42, 32, 46, 14, 42, 14, 46},
  {43, 43, 47, 47, 25, 15, 33, 15},
  {40, 12, 44, 12, 40, 26, 44, 34},
  {13, 27, 13, 35, 41, 41, 45, 45},
  {28, 41, 28, 22, 38, 43, 38, 22},
  {42, 40, 23, 23, 29, 39, 29, 39},
  {41, 36, 20, 36, 43, 30, 20, 30},
  {37, 31, 37, 31, 42, 40, 21, 21},
  {28, 18, 28, 45, 38, 18, 38, 47},
  {19, 19, 46, 44, 29, 39, 29, 39},
  {16, 36, 45, 36, 16, 30, 47, 30},
  {37, 31, 37, 31, 17, 17, 46, 44},
  {12, 4, 1, 3, 34, 34, 1, 3},
  {5, 35, 0, 0, 13, 35, 2, 2},
  {32, 32, 1, 3, 6, 14, 1, 3},
  {33, 15, 0, 0, 33, 7, 2, 2},
  {16, 0, 20, 8, 16, 30, 20, 30},
  {1, 31, 9, 31, 17, 17, 21, 21},
  {28, 18, 28, 22, 2, 18, 10, 22},
  {19, 19, 23, 23, 29, 3, 29, 11},
  {9, 11, 12, 4, 9, 11, 26, 26},
  {8, 8, 5, 27, 10, 10, 13, 27},
  {9, 11, 24, 24, 9, 11, 6, 14},
  {8, 8, 25, 15, 10, 10, 25, 7},
  {0, 18, 8, 22, 38, 18, 38, 22},
  {19, 19, 23, 23, 1, 39, 9, 39},
  {16, 36, 20, 36, 16, 2, 20, 10},
  {37, 3, 37, 11, 17, 17, 21, 21},
  {4, 17, 4, 46, 14, 19, 14, 46},
  {18, 16, 47, 47, 5, 15, 5, 15},
  {17, 12, 44, 12, 19, 6, 44, 6},
  {13, 7, 13, 7, 18, 16, 45, 45},
  {4, 42, 4, 21, 14, 42, 14, 23},
  {43, 43, 22, 20, 5, 15, 5, 15},
  {40, 12, 21, 12, 40, 6, 23, 6},
  {13, 7, 13, 7, 41, 41, 22, 20}
};

const unsigned char subpix3[48][8] = {
  {0, 7, 1, 6, 3, 4, 2, 5},
  {7, 4, 6, 5, 0, 3, 1, 2},
  {4, 3, 5, 2, 7, 0, 6, 1},
  {3, 0, 2, 1, 4, 7, 5, 6},
  {1, 0, 6, 7, 2, 3, 5, 4},
  {0, 3, 7, 4, 1, 2, 6, 5},
  {3, 2, 4, 5, 0, 1, 7, 6},
  {2, 1, 5, 6, 3, 0, 4, 7},
  {6, 1, 7, 0, 5, 2, 4, 3},
  {1, 2, 0, 3, 6, 5, 7, 4},
  {2, 5, 3, 4, 1, 6, 0, 7},
  {5, 6, 4, 7, 2, 1, 3, 0},
  {7, 6, 0, 1, 4, 5, 3, 2},
  {6, 5, 1, 2, 7, 4, 0, 3},
  {5, 4, 2, 3, 6, 7, 1, 0},
  {4, 7, 3, 0, 5, 6, 2, 1},
  {6, 7, 5, 4, 1, 0, 2, 3},
  {7, 0, 4, 3, 6, 1, 5, 2},
  {0, 1, 3, 2, 7, 6, 4, 5},
  {1, 6, 2, 5, 0, 7, 3, 4},
  {2, 3, 1, 0, 5, 4, 6, 7},
  {3, 4, 0, 7, 2, 5, 1, 6},
  {4, 5, 7, 6, 3, 2, 0, 1},
  {5, 2, 6, 1, 4, 3, 7, 0},
  {7, 0, 6, 1, 4, 3, 5, 2},
  {0, 3, 1, 2, 7, 4, 6, 5},
  {3, 4, 2, 5, 0, 7, 1, 6},
  {4, 7, 5, 6, 3, 0, 2, 1},
  {6, 7, 1, 0, 5, 4, 2, 3},
  {7, 4, 0, 3, 6, 5, 1, 2},
  {4, 5, 3, 2, 7, 6, 0, 1},
  {5, 6, 2, 1, 4, 7, 3, 0},
  {1, 6, 0, 7, 2, 5, 3, 4},
  {6, 5, 7, 4, 1, 2, 0, 3},
  {5, 2, 4, 3, 6, 1, 7, 0},
  {2, 1, 3, 0, 5, 6, 4, 7},
  {0, 1, 7, 6, 3, 2, 4, 5},
  {1, 2, 6, 5, 0, 3, 7, 4},
  {2, 3, 5, 4, 1, 0, 6, 7},
  {3, 0, 4, 7, 2, 1, 5, 6},
  {1, 0, 2, 3, 6, 7, 5, 4},
  {0, 7, 3, 4, 1, 6, 2, 5},
  {7, 6, 4, 5, 0, 1, 3, 2},
  {6, 1, 5, 2, 7, 0, 4, 3},
  {5, 4, 6, 7, 2, 3, 1, 0},
  {4, 3, 7, 0, 5, 2, 6, 1},
  {3, 2, 0, 1, 4, 5, 7, 6},
  {2, 5, 1, 6, 3, 4, 0, 7}
};

/*! This function computes a Peano-Hilbert key for an integer triplet (x,y,z),
  *  with x,y,z in the range between 0 and 2^bits-1.
  */
peanokey peano_hilbert_key(int x, int y, int z, int bits)
{
  int mask;
  unsigned char rotation = 0;
  peanokey key = 0;

  for(mask = 1 << (bits - 1); mask > 0; mask >>= 1)
    {
      unsigned char pix = ((x & mask) ? 4 : 0) | ((y & mask) ? 2 : 0) | ((z & mask) ? 1 : 0);

      key <<= 3;
      key |= subpix3[rotation][pix];
      rotation = rottable3[rotation][pix];
    }

  return key;
}



peanokey morton_key(int x, int y, int z, int bits)
{
  int mask;
  peanokey morton = 0;

  for(mask = 1 << (bits - 1); mask > 0; mask >>= 1)
    {
      morton <<= 3;
      morton += ((z & mask) ? 4 : 0) + ((y & mask) ? 2 : 0) + ((x & mask) ? 1 : 0);
    }

  return morton;
}


peanokey peano_and_morton_key(int x, int y, int z, int bits, peanokey *morton_key)
{
  int mask;
  unsigned char rotation = 0;
  peanokey key = 0;
  peanokey morton = 0;


  for(mask = 1 << (bits - 1); mask > 0; mask >>= 1)
    {
      unsigned char pix = ((x & mask) ? 4 : 0) | ((y & mask) ? 2 : 0) | ((z & mask) ? 1 : 0);

      key <<= 3;
      key |= subpix3[rotation][pix];
      rotation = rottable3[rotation][pix];

      morton <<= 3;
      morton += ((z & mask) ? 4 : 0) + ((y & mask) ? 2 : 0) + ((x & mask) ? 1 : 0);
    }

  *morton_key = morton;

  return key;
}




static int quadrants[24][2][2][2] = {
  /* rotx=0, roty=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rotx=1, roty=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rotx=2, roty=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rotx=3, roty=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rotx=4, roty=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rotx=5, roty=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
  11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static int rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static int roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };


peanokey peano_hilbert_key_old(int x, int y, int z, int bits)
{
  int i, quad, bitx, bity, bitz;
  int mask, rotation, rotx, roty, sense;
  peanokey key;


  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;


  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}
    }

  return key;
}


peanokey peano_and_morton_key_old(int x, int y, int z, int bits, peanokey *morton_key)
{
  int i, quad, bitx, bity, bitz;
  int mask, rotation, rotx, roty, sense;
  peanokey key, morton;


  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;
  morton = 0;

  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}

      morton <<= 3;
      morton += (bitz << 2) + (bity << 1) + bitx;
    }

  *morton_key = morton;

  return key;
}


peanokey morton_key_old(int x, int y, int z, int bits)
{
  int i, bitx, bity, bitz;
  int mask;
  peanokey morton;


  mask = 1 << (bits - 1);
  morton = 0;

  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      morton <<= 3;
      morton += (bitz << 2) + (bity << 1) + bitx;
    }

  return morton;
}


static void msort_peano_with_tmp(struct peano_hilbert_data *b, size_t n, struct peano_hilbert_data *t)
{
  struct peano_hilbert_data *tmp;
  struct peano_hilbert_data *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_peano_with_tmp(b1, n1, t);
  msort_peano_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(b1->key <= b2->key)
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(struct peano_hilbert_data));
  memcpy(b, t, (n - n2) * sizeof(struct peano_hilbert_data));
}

void mysort_peano(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
  const size_t size = n * s;

  struct peano_hilbert_data *tmp =
    (struct peano_hilbert_data *) mymalloc("struct peano_hilbert_data *tmp", size);

  msort_peano_with_tmp((struct peano_hilbert_data *) b, n, tmp);

  myfree(tmp);
}
