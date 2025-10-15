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

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"
#ifdef ACC
#include "../OpenACC/gpuallvars_branch.h"
#endif






/* This routine allocates memory for
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{

#ifdef KD_EXTRA_TIMER_OUTPUT_TREEBUILD
  double t1_kd,t0_kd_merk,t0_kd=second();
#endif

  size_t bytes;

  double bytes_tot = 0;

  int NTaskTimesThreads;

  NTaskTimesThreads = maxThreads * NTask;

  Exportflag = (int *) mymalloc("Exportflag", NTaskTimesThreads * sizeof(int));
  Exportindex = (int *) mymalloc("Exportindex", NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *) mymalloc("Exportnodecount", NTaskTimesThreads * sizeof(int));

  Send_count = (int *) mymalloc("Send_count", sizeof(int) * NTask);
  Send_offset = (int *) mymalloc("Send_offset", sizeof(int) * NTask);
  Recv_count = (int *) mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc("Recv_offset", sizeof(int) * NTask);
#ifdef IMPORT_ALLtoALLv_FROM_G4
  //  printf("Task %d: WinCreate: %d %d\n",ThisTask,NTask*sizeof(MPI_INT), sizeof(MPI_INT));
  //  MPI_Win_create(Recv_count, NTask*sizeof(MPI_INT), sizeof(MPI_INT), MPI_INFO_NULL, MYMPI_COMM_WORLD, &win);
#endif

  if(ThisTask == 0)
    printf("Max number of particles per task: %llu\n", (unsigned long long)All.MaxPart);

  ProcessedFlag = (unsigned char *) mymalloc("ProcessedFlag", bytes = All.MaxPart * sizeof(unsigned char));
  bytes_tot += bytes;

  NextActiveParticle = (int *) mymalloc("NextActiveParticle", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

#ifdef _OPENMP
  ActiveParticleList = (int *) mymalloc("ActiveParticleList", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;
#endif

  NextInTimeBin = (int *) mymalloc("NextInTimeBin", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  PrevInTimeBin = (int *) mymalloc("PrevInTimeBin", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;


  if(All.MaxPart > 0)
    {
      if(!(P = (struct particle_data *) mymalloc("P", bytes = All.MaxPart * sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}

#ifdef GADGET3_IO_LIB
      memset(P, -1, bytes);
#endif

      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for particle storage.\n\n", bytes_tot / (1024.0 * 1024.0));
    }



  if(All.MaxPartSph > 0)
    {
      bytes_tot = 0;

      if(!
	 (SphP =
	  (struct sph_particle_data *) mymalloc("SphP", bytes =
						All.MaxPartSph * sizeof(struct sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
#ifdef GADGET3_IO_LIB
      memset(SphP, -1, bytes);
#endif

      bytes_tot += bytes;

#ifdef AR_SOA_LITE //with dictionary this will became automatic
      AR_ALLOC(SphPSoa.wakeup, short int, All.MaxPartSph, bytes_tot);
#endif
      
      
      
      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of SPH data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }

#ifdef AXION_DM

  if(All.MaxPartAx > 0)
    {
      bytes_tot = 0;

      if(!
	 (AxP =
	  (struct ax_particle_data *) mymalloc("AxP", bytes =
						All.MaxPartAx * sizeof(struct ax_particle_data))))
	{
	  printf("failed to allocate memory for `AxP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of AX data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }

#endif

#ifdef LT_STELLAREVOLUTION
  if(All.MaxPartMet > 0)
    {
      bytes_tot = 0;

      if(!
	 (MetP =
	  (struct met_particle_data *) mymalloc("MetP", bytes =
						All.MaxPartMet * sizeof(struct met_particle_data))))
	{
	  printf("failed to allocate memory for `MetP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of MET data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }
#endif

#if defined(BLACK_HOLES)
  if(All.MaxPartBH > 0)
    {
      bytes_tot = 0;

      if(!
	 (BHP =
	  (struct bh_particle_data *) mymalloc("BHP", bytes =
					       All.MaxPartBH * sizeof(struct bh_particle_data))))
	{
	  printf("failed to allocate memory for `BHP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of BH data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }
#endif


#ifdef WRITE_KEY_FILES
  if(All.MaxPart > 0)
    {
      bytes_tot = 0;
      KeyIndex = (peanokey *) mymalloc("KeyIndex", bytes = All.MaxPart * sizeof(peanokey));
      bytes_tot += bytes;
      NPartPerKey = (int *) mymalloc("KeyNpart", bytes = All.MaxPart * sizeof(int));
      bytes_tot += bytes;
      PartKeyOffset = (int *) mymalloc("KeyOffset", bytes = All.MaxPart * sizeof(int));
      bytes_tot += bytes;
      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of KEY data.\n\n", bytes_tot / (1024.0 * 1024.0));
    }
#endif


  // TODO : (David) Find convenient location to destroy objects (i.e. deallocate memory)
#if GADGET_HYDRO == HYDRO_MFM
  eos = new AdiabaticEos<NUMDIMS>(MyFloat(GAMMA));
  fluxSolver = new FluxSolver<NUMDIMS>(eos);
//#if MFM_SLOPE_LIMITER == NULL_LIMITER
//  slopeLimiter = new NullSlopeLimiter<NUMDIMS>();
//#else
//  slopeLimiter = new ZeroSlopeLimiter<NUMDIMS>();
//#endif
#endif

  
#ifdef AR_TREEBUILD_PARALLEL
  double tree_alloc_factor = DMAX(All.TreeAllocFactor, 1.);
  int _MaxNodes = (int) (All.MaxPart * tree_alloc_factor)+ MaxTopNodes;
  if(!( NodesLockList = (omp_lock_t *) mymalloc("NodesLockList", bytes = (_MaxNodes+ 1)* sizeof(omp_lock_t )))){
    printf("failed to allocate memory for %d tree-nodes locks (%g MB).\n", _MaxNodes, bytes / (1024.0 * 1024.0));
    endrun(3);
  }
  bytes_tot += bytes;

  for(int i=0; i<_MaxNodes +1; i++){ //note loop gets worse if parallelised
    omp_init_lock(&NodesLockList[i]);
  }

  VERBOSE(-1,"Allocated %d locks: my TreeAllocFactor=%g; All.MaxPart=%d; NTopnodes=%d",_MaxNodes +1, tree_alloc_factor,All.MaxPart,NTopnodes);
  NumberOfAllocatedNodeLocks = _MaxNodes;
#ifdef KD_EXTRA_TIMER_OUTPUT_TREEBUILD
  t1_kd=second();
  if(ThisTask ==0)
    printf("EXTRA TIMER TREEBUILD:   initialising openmp locks took %g sec\n",timediff(t0_kd,t1_kd));
  t0_kd=second();
#endif
  
  
#endif
  
  
  
}
