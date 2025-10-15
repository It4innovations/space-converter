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
#include <errno.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/file.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"
#include "../CodeBase/domain.h"

static FILE *fd;

static void in(int *x, int modus);
static void byten(void *x, size_t n, int modus);

int old_MaxPart = 0, new_MaxPart;


/* This function reads or writes the restart files.
 * Each processor writes its own restart file, with the
 * I/O being done in parallel. To avoid congestion of the disks
 * you can tell the program to restrict the number of files
 * that are simultaneously written to NumFilesWrittenInParallel.
 *
 * If modus>0  the restart()-routine reads,
 * if modus==0 it writes a restart file.
 */
void restart(int modus)
{
  char buf[200], buf_bak[200], buf_mv[500];
  double save_PartAllocFactor;

  int nprocgroup, masterTask, groupTask;
  struct global_data_all_processes all_task0;
  int nmulti = MULTIPLEDOMAINS;

#ifdef LT_STELLAREVOLUTION
  int save_NumFilesPerSnapshot, save_NumFilesWrittenInParallel;
  double safe_SFfactor;
  double save_LLv_Step_Prec, save_SnII_Step_Prec;
  int buffer;
  double save_SofteningGasMaxPhys,
    save_SofteningHaloMaxPhys,
    save_SofteningDiskMaxPhys,
    save_SofteningBulgeMaxPhys, save_SofteningStarsMaxPhys, save_SofteningBndryMaxPhys;
  double save_MinChemTimeStep, save_MinChemSpreadL, save_MaxChemSpreadL;
#endif

#if defined(BLACK_HOLES)
  int bhbuffer;
#endif

  double t0_rst, t1_rst;

  t0_rst = second();

  if(ThisTask == 0 && modus == 0)
    {
      sprintf(buf, "%s/restartfiles", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MYMPI_COMM_WORLD);

  sprintf(buf, "%s/restartfiles/%s.%d", All.OutputDir, All.RestartFile, ThisTask);
  sprintf(buf_bak, "%s/restartfiles/%s.%d.bak", All.OutputDir, All.RestartFile, ThisTask);
  sprintf(buf_mv, "mv %s %s", buf, buf_bak);

  if((NTask < All.NumFilesWrittenInParallel))
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(2131);
    }

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    {
      nprocgroup++;
    }

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))
	{
	  if(!modus)
	    {
#ifndef NOCALLSOFSYSTEM
	      int ret;

	      ret = system(buf_mv);	/* move old restart files to .bak files */
#else
#ifdef CRENAME
	      if(rename(buf, buf_bak) != 0)
		{
		  printf("WARNING: Could not back up restart file %s\n%s\n", buf, strerror(errno));
		}
#endif
#endif
	    }
	}
    }

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	{
	  if(modus)
	    {
	      if(!(fd = fopen(buf, "r")))
		{
		  if(!(fd = fopen(buf_bak, "r")))
		    {
		      printf("Restart file '%s' nor '%s' found.\n", buf, buf_bak);
		      endrun(7870);
		    }
		}
	    }
	  else
	    {
	      if(!(fd = fopen(buf, "w")))
		{
		  printf("Restart file '%s' cannot be opened.\n", buf);
		  endrun(7878);
		}
	    }


	  save_PartAllocFactor = All.PartAllocFactor;

#ifdef LT_STELLAREVOLUTION
	  save_NumFilesPerSnapshot = All.NumFilesPerSnapshot;
	  save_NumFilesWrittenInParallel = All.NumFilesWrittenInParallel;

	  save_SofteningGasMaxPhys = All.SofteningGasMaxPhys;
	  save_SofteningHaloMaxPhys = All.SofteningHaloMaxPhys;
	  save_SofteningDiskMaxPhys = All.SofteningDiskMaxPhys;
	  save_SofteningBulgeMaxPhys = All.SofteningBulgeMaxPhys;
	  save_SofteningStarsMaxPhys = All.SofteningStarsMaxPhys;
	  save_SofteningBndryMaxPhys = All.SofteningBndryMaxPhys;

	  save_SnII_Step_Prec = All.SnII_Step_Prec;
	  save_LLv_Step_Prec = All.LLv_Step_Prec;

	  save_MinChemTimeStep = All.MinChemTimeStep;
	  save_MinChemSpreadL = All.MinChemSpreadL;
	  save_MaxChemSpreadL = All.MaxChemSpreadL;

	  safe_SFfactor = All.SFfactor;

	  if(!modus)
	    *(float *) &buffer = (float) All.Time;
	  in(&buffer, modus);
	  if(!modus)
	    buffer = sizeof(struct global_data_all_processes);
	  in(&buffer, modus);
	  if(modus && buffer != sizeof(struct global_data_all_processes))
	    {
	      printf
		("in file <%s> :: sizes of the current All structure and of the stored All structure are different! (%d vs %d bytes)\n",
		 buf, (int) buffer, (int) sizeof(struct global_data_all_processes));
	      endrun(23);
	    }
#endif

	  /* common data  */
	  byten(&All, sizeof(struct global_data_all_processes), modus);

#ifdef LT_STELLAREVOLUTION
	  All.NumFilesPerSnapshot = save_NumFilesPerSnapshot;
	  All.NumFilesWrittenInParallel = save_NumFilesWrittenInParallel;

	  All.SnII_Step_Prec = save_SnII_Step_Prec;
	  All.LLv_Step_Prec = save_LLv_Step_Prec;
	  All.MinChemTimeStep = save_MinChemTimeStep;

	  All.SofteningGasMaxPhys = save_SofteningGasMaxPhys;
	  All.SofteningHaloMaxPhys = save_SofteningHaloMaxPhys;
	  All.SofteningDiskMaxPhys = save_SofteningDiskMaxPhys;
	  All.SofteningBulgeMaxPhys = save_SofteningBulgeMaxPhys;
	  All.SofteningStarsMaxPhys = save_SofteningStarsMaxPhys;
	  All.SofteningBndryMaxPhys = save_SofteningBndryMaxPhys;

	  All.MinChemSpreadL = save_MinChemSpreadL;
	  All.MaxChemSpreadL = save_MaxChemSpreadL;
#endif

	  if(ThisTask == 0 && modus > 0)
	    all_task0 = All;

	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MYMPI_COMM_WORLD);
	    }



	  if(modus)		/* read */
	    {
	      if(All.PartAllocFactor != save_PartAllocFactor)
		{
		  old_MaxPart = All.MaxPart;	/* old MaxPart */

		  if(ThisTask == 0)
		    printf("PartAllocFactor changed: %f/%f , adapting bounds ...\n",
			   All.PartAllocFactor, save_PartAllocFactor);

		  All.PartAllocFactor = save_PartAllocFactor;
		  All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));
		  All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));
#ifdef INHOMOG_GASDISTR_HINT
		  All.MaxPartSph = All.MaxPart;
#endif

		  new_MaxPart = All.MaxPart;

#ifdef LT_STELLAREVOLUTION
		  if(ThisTask == 0)
		    printf("All.TotNumPart=%llu, All.TotN_gas=%llu, All.TotN_stars=%llu \n",
			   (unsigned long long) All.TotNumPart, (unsigned long long) All.TotN_gas,
			   (unsigned long long) All.TotN_stars);
		  if(All.TotN_stars == 0)
		    All.MaxPartMet =
		      All.PartAllocFactor * (All.TotN_gas / NTask) * safe_SFfactor * All.Generations;
		  else
		    All.MaxPartMet =
		      All.PartAllocFactor * (All.TotN_stars / NTask +
					     (All.TotN_gas / NTask) * safe_SFfactor * All.Generations);
		  if(ThisTask == 0)
		    printf("All.MaxPart=%llu, All.MaxPartSph=%llu, All.MaxPartMet=%llu \n",
			   (unsigned long long) All.MaxPart, (unsigned long long) All.MaxPartSph,
			   (unsigned long long) All.MaxPartMet);
		  All.SFfactor = safe_SFfactor;
#endif
#if defined(BLACK_HOLES)
		  if(All.TotBHs == 0)
		    All.MaxPartBH = All.PartAllocFactor * (All.TotN_gas / NTask) * All.BHfactor;
		  else
		    All.MaxPartBH = All.PartAllocFactor * (All.TotBHs / NTask +
							   (All.TotN_gas / NTask) * All.BHfactor);
#endif

		  save_PartAllocFactor = -1;
		}

#ifdef LT_STELLAREVOLUTION
	      if(ThisTask == 0)
		printf("SFfactor changed: %f/%f , adapting bounds ...\n", All.SFfactor, safe_SFfactor);

	      if(All.SFfactor != safe_SFfactor)
		{
		  All.MaxPartMet = (All.MaxPartMet / All.SFfactor) * safe_SFfactor;
		  All.SFfactor = safe_SFfactor;
		}
#endif

	      if(all_task0.Time != All.Time)
		{
		  printf("The restart file on task=%d is not consistent with the one on task=0\n", ThisTask);
		  fflush(stdout);
		  endrun(16);
		}

	      allocate_memory();
	    }

	  in(&NumPart, modus);
#ifdef LT_STELLAREVOLUTION
	  if(!modus)
	    buffer = sizeof(struct particle_data);
	  in(&buffer, modus);
	  if(modus && buffer != sizeof(struct particle_data))
	    {
	      printf
		("in file <%s> :: sizes of the current particle data and of the stored particle data are different! (%d vs %d bytes)\n",
		 buf, buffer, (int) sizeof(struct particle_data));
	      endrun(23);
	    }
#endif
	  if(NumPart > All.MaxPart)
	    {
	      printf
		("it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		 NumPart / (((double) All.TotNumPart) / NTask));
	      printf("fatal error\n");
	      endrun(22);
	    }

	  if(modus)		/* read */
	    {
	      if(old_MaxPart)
		All.MaxPart = old_MaxPart;	/* such that tree is still valid */
	    }

	  /* Particle data  */
	  byten(&P[0], (size_t) NumPart * sizeof(struct particle_data), modus);

	  in(&N_gas, modus);
#ifdef LT_STELLAREVOLUTION
	  if(!modus)
	    buffer = sizeof(struct sph_particle_data);
	  in(&buffer, modus);
	  if(modus && buffer != sizeof(struct sph_particle_data))
	    {
	      printf
		("in file <%s> :: sizes of the current sph particle data and of the stored sph particle data are different! (%d vs %d bytes)\n",
		 buf, buffer, (int) sizeof(struct sph_particle_data));
	      endrun(23);
	    }
#endif
	  if(N_gas > 0)
	    {
	      if(N_gas > All.MaxPartSph)
		{
		  printf
		    ("SPH: it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		     N_gas / (((double) All.TotN_gas) / NTask));
		  printf("fatal error\n");
		  endrun(222);
		}
	      /* Sph-Particle data  */
	      byten(&SphP[0], (size_t) N_gas * sizeof(struct sph_particle_data), modus);
	      
	    }

#if defined(AR_SOA_LITE)  
#ifndef AR_SOA_LITE_FROM_RESTART_AOS
	  //unfortunately we must clone the byten for each property
	  byten(SphPSoa.wakeup, (size_t) N_gas * sizeof(short int), modus);
#endif
#endif
#ifdef AR_SOA_LITE_FROM_RESTART_AOS //if you want to restart with soa-lite active
	      for(int i=0;i<N_gas; i++) SPHP(i,wakeup) = SphP[i].wakeup;
#endif

	  /* write state of random number generator */
	  byten(gsl_rng_state(random_generator), gsl_rng_size(random_generator), modus);
	  byten(&SelRnd, sizeof(SelRnd), modus);

	  /* write flags for active timebins */
	  byten(TimeBinActive, TIMEBINS * sizeof(int), modus);

	  /* now store relevant data for tree */
#ifdef SFR
	  in(&Stars_converted, modus);
#endif

#ifdef LT_STELLAREVOLUTION
	  in(&N_stars, modus);
	  if(!modus)
	    buffer = sizeof(struct met_particle_data);
	  in(&buffer, modus);
	  if(modus && buffer != sizeof(struct met_particle_data))
	    {
	      printf
		("in file <%s> :: sizes of the current (star) met particle data and of the stored (star) met particle data are different! (%d vs %d bytes)\n",
		 buf, buffer, (int) sizeof(struct met_particle_data));
	      endrun(23);
	    }


	  if(N_stars > 0)
	    {
	      if(N_stars > All.MaxPartMet)
		{
		  printf
		    ("MET: it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		     N_stars / (((double) All.TotN_stars) / NTask));
		  printf("fatal error\n");
		  endrun(2222);
		}
	      /* Sph-Particle data  */
	      byten(&MetP[0], (size_t) N_stars * sizeof(struct met_particle_data), modus);
	    }

	  byten(&All.MetalsCheckSum, sizeof(long double)*LT_NMetP, modus);
#endif

#if defined(BLACK_HOLES)
	  in(&N_BHs, modus);
	  if(!modus)
	    bhbuffer = sizeof(struct bh_particle_data);
	  in(&bhbuffer, modus);
	  if(modus && bhbuffer != sizeof(struct bh_particle_data))
	    {
	      printf
		("in file <%s> :: sizes of the current bh particle data and of the stored bh particle data are different! (%d vs %d bytes)\n",
		 buf, bhbuffer, (int) sizeof(struct bh_particle_data));
	      endrun(23);
	    }

	  if(N_BHs > 0)
	    {
	      if(N_BHs > All.MaxPartBH)
		{
		  printf
		    ("BH: it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		     N_BHs / (((double) All.TotBHs) / NTask));
		  printf("fatal error\n");
		  endrun(2222);
		}
	      /* Sph-Particle data  */
	      byten(&BHP[0], (size_t) N_BHs * sizeof(struct bh_particle_data), modus);
	    }
#endif

	  /* now store relevant data for tree */

	  in(&nmulti, modus);
	  if(modus != 0 && nmulti != MULTIPLEDOMAINS)
	    {
	      if(ThisTask == 0)
		printf
		  ("Looks like you changed MULTIPLEDOMAINS from %d to %d.\nWe will need to discard tree stored in restart files and construct a new one.\n",
		   nmulti, (int) MULTIPLEDOMAINS);

	      /* In this case we must do a new domain decomposition! */
	    }
	  else
	    {
#ifndef NOTREEINRESTARTFILE
	      in(&NTopleaves, modus);
	      in(&NTopnodes, modus);

	      if(modus)		/* read */
		{
		  domain_allocate();
		  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
		}

	      in(&Numnodestree, modus);

	      if(Numnodestree > MaxNodes)
		{
		  printf
		    ("Tree storage: it seems you have reduced(!) 'PartAllocFactor' below the value needed to load the restart file (task=%d). "
		     "Numnodestree=%d  MaxNodes=%d\n", ThisTask, Numnodestree, MaxNodes);
		  endrun(221);
		}

	      byten(Nodes_base, Numnodestree * sizeof(struct NODE), modus);
	      byten(Extnodes_base, Numnodestree * sizeof(struct extNODE), modus);

	      byten(Father, NumPart * sizeof(int), modus);

	      byten(Nextnode, NumPart * sizeof(int), modus);
	      byten(Nextnode + All.MaxPart, NTopnodes * sizeof(int), modus);

	      byten(DomainStartList, NTask * MULTIPLEDOMAINS * sizeof(int), modus);
	      byten(DomainEndList, NTask * MULTIPLEDOMAINS * sizeof(int), modus);
	      byten(TopNodes, NTopnodes * sizeof(struct topnode_data), modus);
	      byten(DomainTask, NTopnodes * sizeof(int), modus);
	      byten(DomainNodeIndex, NTopleaves * sizeof(int), modus);

	      byten(DomainCorner, 3 * sizeof(double), modus);
	      byten(DomainCenter, 3 * sizeof(double), modus);
	      byten(&DomainLen, sizeof(double), modus);
	      byten(&DomainFac, sizeof(double), modus);
#endif
	    }

	  fclose(fd);
	}
      else			/* wait inside the group */
	{
	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MYMPI_COMM_WORLD);
	    }
	}

      MPI_Barrier(MYMPI_COMM_WORLD);
    }

  if(modus != 0)
#ifndef NOTREEINRESTARTFILE
    if(nmulti != MULTIPLEDOMAINS)	/* in this case we must force a domain decomposition */
#endif
      {
	if(ThisTask == 0)
	  printf("Doing extra domain decomposition\n");

	domain_Decomposition(0, 0);
	if(TreeReconstructFlag)
	  {
	    if(ThisTask == 0)
	      printf("Tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

	    force_treebuild(NumPart, NULL);

	    CPU_Step[CPU_TREEBUILD] += measure_time();

	    TreeReconstructFlag = 0;

	    if(ThisTask == 0)
	      printf("Tree construction done.\n");
	  }

      }





  t1_rst = second();

  if(ThisTask == 0)
    printf("restatfile in modus %d took %g sec\n", modus, timediff(t0_rst, t1_rst));

}



/* reads/writes n bytes
 */
void byten(void *x, size_t n, int modus)
{
  if(modus)
    my_fread(x, n, 1, fd);
  else
    my_fwrite(x, n, 1, fd);
}


/* reads/writes one int
 */
void in(int *x, int modus)
{
  if(modus)
    my_fread(x, 1, sizeof(int), fd);
  else
    my_fwrite(x, 1, sizeof(int), fd);
}
