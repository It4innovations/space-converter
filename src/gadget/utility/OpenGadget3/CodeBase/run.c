
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
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"
#include "greentree.h"


#ifdef ACC_RUN
#include "../OpenACC/gpuallvars_branch.h"
#endif

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */

void domain_decomposition_intensity_decision();
void domain_decomposition_intensity_execute();


void run(void)
{
  PUSH_T("run initial part");
  double t0, t1;
#ifdef FPE_TRAPPING
  feenableexcept(FE_INVALID);
#endif

  CPU_Step[CPU_MISC] += measure_time();

  if(RestartFlag != 1)		/* need to compute forces at initial synchronization time, unless we restarted from restart files */
    {

      output_log_messages();

      domain_Decomposition(0, 0);
#ifdef RECOMPOSE_DOMAIN_FACTOR
      LastDomainUpdate = 0;
#endif

      set_non_standard_physics_for_current_time();

      compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */

      compute_densities();	/* densities for synchronous particles */

#ifdef GM_STARDENSITY		/*...and stellar densities for star particles */
      if(ThisTask == 0)
	{
	  printf("Start stellar density computation...\n");
	  fflush(stdout);
	}
      density_stars();
#endif


#if GADGET_HYDRO == HYDRO_SPH || GADGET_HYDRO == HYDRO_PESPH
      compute_hydro_accelerations();	/* hydro-accels for synchronous particles */
#elif GADGET_HYDRO == HYDRO_MFM
      compute_gradients();	// Compute gradients once smoothing lengths are known

      compute_limiters();	// Compute limiters once all gradients are known

      find_timesteps();		// find-timesteps for flux calculation

      compute_fluxes();		// Compute Godunov fluxes
#endif

      calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */

#if (defined(fSIDM) || defined(rSIDM)) && defined(mSIDM_TIMESTEP)
      msidm_init_omega();
#endif

    }
  POP_T();

/*
									 main timestep iteration loop
*/
  while(1)
    {
      PUSH_T("run while loop");
#ifdef GDE_BIGFLOAT
      xErrNo = 0;
#endif

#ifdef KD_COOLING_TEST
      savepositions(++All.SnapshotFileCount);	/* this will be overwritten if All.TimeMax is increased and the run is continued */
      break;
#endif

      compute_statistics();	/* regular statistics outputs (like total energy) */

      create_snapshot_if_desired();

      write_cpu_log();		/* output some CPU usage log-info (accounts for everything needed up to the current sync-point) */

#ifdef WRITE_ANOTHER_SNAPSHOT
#ifdef OUTPUT_LIGHTCONES
      if(ThisTask == 0)
	{
	  printf("WRITE_ANOTHER_SNAPSHOT\n");
	}
      lightcone_output(All.SnapshotFileCount++);
#else
      savepositions(All.SnapshotFileCount++);	/* this will be overwritten if All.TimeMax is increased and the run is continued */
#endif
      break;
#endif

      if(All.Ti_Current >= TIMEBASE)	/* check whether we reached the final time */
	{
	  if(ThisTask == 0)
	    {
	      printf("\nFinal time=%g reached. Simulation ends.\n", All.TimeMax);
	    }

#ifndef DISABLE_OUTPUT_AT_FINAL_TIME
	  restart(0);		/* write a restart file to allow continuation of the run for a larger value of TimeMax */

	  if(All.Ti_lastoutput != All.Ti_Current)	/* make a snapshot at the final time in case none has produced at this time */
	    {
#ifdef OUTPUT_LIGHTCONES
	      if(ThisTask == 0)
		{
		  printf("All.Ti_Current >= TIMEBASE\n");
		}
	      lightcone_output(All.SnapshotFileCount++);
#else
	      savepositions(All.SnapshotFileCount++);	/* this will be overwritten if All.TimeMax is increased and the run is continued */
#endif
	    }
#endif
	  break;
	}

      find_timesteps();		/* find-timesteps */


      do_first_halfstep_kick();	/* half-step kick at beginning of timestep for synchronous particles */

      find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
					 * If needed, this function will also write an output file
					 * at the desired time. */

      domain_decomposition_intensity_decision();

      output_log_messages();	/* write some info to log-files */

      set_non_standard_physics_for_current_time();	/* update auxiliary physics for current time */

      domain_decomposition_intensity_execute();

      compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */

#ifdef AR_TREE_STATS
      ar_tree_stats();		//provides information on number of particles per tree-level
#endif
      compute_densities();	/* densities for synchronous particles */
      /* GM WARNING : Supernovae evaluation is now here! */

#ifdef GM_STARDENSITY		/*...and stellar densities for star particles */
      if(ThisTask == 0)
	{
	  printf("Start stellar density computation...\n");
	  fflush(stdout);
	}
      density_stars();
#endif


#if GADGET_HYDRO == HYDRO_SPH || GADGET_HYDRO == HYDRO_PESPH
      compute_hydro_accelerations();	/* hydro-accels for synchronous particles */
#elif GADGET_HYDRO == HYDRO_MFM
      compute_gradients();	// Compute gradients once smoothing lenghts are known

      compute_limiters();	// Compute limiters once all gradients are known

      compute_fluxes();		// Compute Godunov fluxes
#endif

      do_second_halfstep_kick();	/* this does the half-step kick at the end of the timestep */

      calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */


#if defined(fSIDM) || defined(rSIDM)
      do_msidm_full_kick();	/* computes new velocities due to self-interactions */
#endif

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	{
	  printf("EXTRA TIMER: non standard physics took %g sec\n", timediff(t0, t1));
	}
#endif

      int stopflag = check_stop_condition();
      if(stopflag > 0)
	{
	  return;
	}

#ifdef GDE_BIGFLOAT
      if(xErrNo != 0)
	{
	  PANIC("GDE_BIGFLOAT error\n");
	}
#endif

      POP_T();
    }

}

#ifdef DOMAIN_DECOMPOSITION_ON_TB

void domain_decomposition_intensity_decision()
{


  long long tot_count[TIMEBINS];
  long long tot_cumulative[TIMEBINS];
  int i;
  sumup_large_ints(TIMEBINS, TimeBinCount, tot_count);

  for(i = 1, tot_cumulative[0] = tot_count[0]; i < TIMEBINS; i++)
    tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];

  AllDerived.DoDomainDecompositionTBFull = 0;
  AllDerived.DoDomainDecompositionTBDiffuse = 0;
  AllDerived.DoDomainDecompositionTBFix = 0;

  for(int i = TIMEBINS - 1; i >= 0; i--)
    {

      if(tot_count[i] > 0)
	{
	  if(tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart)
	    {

	      AllDerived.DoDomainDecompositionTBFull = i;
	      AllDerived.DoDomainDecompositionTBDiffuse = AllDerived.DoDomainDecompositionTBFull - 1;
	      AllDerived.DoDomainDecompositionTBFix = AllDerived.DoDomainDecompositionTBDiffuse - 1;

	    }
	}
    }

}

void domain_decomposition_intensity_execute()
{

  if(TimeBinActive[AllDerived.DoDomainDecompositionTBFull])
    {
      VERBOSE(0, "Domain Decomposition Intensity: Full");
      domain_Decomposition(0, 0);	/* do domain decomposition if step is big enough, and set new list of active particles  */
      return;
    }

  force_update_tree();		/* update tree dynamically with kicks of last step so that it can be reused */

  if(TimeBinActive[AllDerived.DoDomainDecompositionTBDiffuse])
    {
      VERBOSE(0, "Domain Decomposition Intensity: Diffuse");
#ifdef RECOMPOSE_DOMAIN_DIFFUSE
      domain_diffuse();
      domain_recomposition();
#else
      VERBOSE(0, "No RECOMPOSE_DOMAIN_DIFFUSE flag\n");
#endif
    }
  else if(TimeBinActive[AllDerived.DoDomainDecompositionTBFix])
    {
      VERBOSE(0, "Domain Decomposition Intensity: Fix");
#ifdef RECOMPOSE_DOMAIN
      domain_recomposition();
#else
      VERBOSE(0, "No RECOMPOSE_DOMAIN flag\n");
#endif

    }

  make_list_of_active_particles();	/* now we can set the new chain list of active particles */

}

#else //DOMAIN_DECOMPOSITION_ON_TB


void domain_decomposition_intensity_decision()
{

  // this function is empty if you don't have DOMAIN_DECOMPOSITION_ON_TB set up.
}

void domain_decomposition_intensity_execute()
{				//good old domain decomposition
#ifdef KD_TREEDOMAIN_UPDATE_FREQENCY
  if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency / All.Time * All.TotNumPart)
#else
  if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)	/* check whether we have a big step */
#endif
    {
      domain_Decomposition(0, 0);	/* do domain decomposition if step is big enough, and set new list of active particles  */
#ifdef RECOMPOSE_DOMAIN_FACTOR
      LastDomainUpdate = 0;
#endif
    }
  else
    {
      force_update_tree();	/* update tree dynamically with kicks of last step so that it can be reused */

#ifdef RECOMPOSE_DOMAIN
#ifdef RECOMPOSE_DOMAIN_FACTOR
#ifdef KD_TREEDOMAIN_UPDATE_FREQENCY
      int ThisDomainUpdate = floor(GlobNumForceUpdate /
				   (All.TreeDomainUpdateFrequency / All.Time * All.TotNumPart /
				    RECOMPOSE_DOMAIN_FACTOR));
#else
      int ThisDomainUpdate = floor(GlobNumForceUpdate /
				   (All.TreeDomainUpdateFrequency * All.TotNumPart /
				    RECOMPOSE_DOMAIN_FACTOR));
#endif
      if(ThisDomainUpdate != LastDomainUpdate)
	{
	  LastDomainUpdate = ThisDomainUpdate;
#endif

	  if(ThisTask == 0)
	    printf(" [AR] Recomposing the domain...\n");
	  double t0 = second();
#ifdef RECOMPOSE_DOMAIN_DIFFUSE
	  domain_diffuse();
#endif
	  domain_recomposition();
	  double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
	  if(ThisTask == 0)
	    {
	      printf("EXTRA TIMER: recomposing domain took %g sec\n", timediff(t0, t1));
	    }
#endif

#ifdef RECOMPOSE_DOMAIN_FACTOR
	}
#endif
#endif //RECOMPOSE_DOMAIN

      make_list_of_active_particles();	/* now we can set the new chain list of active particles */

    }

}
#endif //DOMAIN_DECOMPOSITION_ON_TB

int check_stop_condition(void)
{
  double t0 = second();
  /* Check whether we need to interrupt the run */
  int stopflag = 0;
  if(ThisTask == 0)
    {
      FILE *fd;
      char stopfname[1000];
      sprintf(stopfname, "%sstop", All.OutputDir);
      if((fd = fopen(stopfname, "r")))	/* Is the stop-file present? If yes, interrupt the run. */
	{
	  fclose(fd);
	  stopflag = 1;
	  unlink(stopfname);
	}

      if(CPUThisRun > 0.85 * All.TimeLimitCPU)	/* are we running out of CPU-time ? If yes, interrupt run. */
	{
	  printf("reaching time-limit. stopping.\n");
	  stopflag = 2;
	}
    }
  /*for testing purposes it is quite useful to interrupt the run after a well deifned a */
#ifdef AR_STOP_AT_A
  if(All.Time >= AR_STOP_AT_A)
    stopflag = 1;
#endif
  MPI_Bcast(&stopflag, 1, MPI_INT, 0, MYMPI_COMM_WORLD);

  if(stopflag)
    {
#ifndef DISABLE_OUTPUT_AT_FINAL_TIME
      restart(0);		/* write restart file */
#endif
      MPI_Barrier(MYMPI_COMM_WORLD);

      if(stopflag == 2 && ThisTask == 0)
	{
	  FILE *fd;
	  char contfname[1000];
	  sprintf(contfname, "%scont", All.OutputDir);
	  if((fd = fopen(contfname, "w")))
	    fclose(fd);
	  if(All.ResubmitOn)
	    execute_resubmit_command();
	}
      POP_T();
      return stopflag;
    }

  if(ThisTask == 0)
    {
      /* is it time to write one of the regularly space restart-files? */
      if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
	{
	  All.TimeLastRestartFile = CPUThisRun;
	  stopflag = 3;
	}
      else
	{
	  stopflag = 0;
	}
    }

  MPI_Bcast(&stopflag, 1, MPI_INT, 0, MYMPI_COMM_WORLD);

  if(stopflag == 3)
    {
#ifndef DISABLE_OUTPUT_AT_FINAL_TIME
      restart(0);		/* write an occasional restart file */
#endif
      stopflag = 0;
      All.TimeLastRestartFile += report_time();
    }

  set_random_numbers();		/* draw a new list of random numbers */

  report_memory_usage(&HighMark_run, "RUN");

  double t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: checking stop flag took %g sec\n", timediff(t0, t1));
    }
#endif

  return stopflag;
}



void set_non_standard_physics_for_current_time(void)
{
  double t0 = second();

#if defined(RADIATIVE_RATES) || defined(RADIATION)
  init_rad(All.Time);
#endif

#ifdef COOLING
  IonizeParams();		/* set UV background for the current time */
#endif

  double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: set non standard physics took %g sec\n", timediff(t0, t1));
    }
#endif
}



void calculate_non_standard_physics(void)
{
  PUSH_T("non standard physics");
  double tstart = second();
  double t0, t1;

#ifdef CONDUCTION
  if(All.Conduction_Ti_endstep == All.Ti_Current)
    {
      t0 = second();
      conduction();
      t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	{
	  printf("EXTRA TIMER: conduction took %g sec\n", timediff(t0, t1));
	}
#endif

    }
#endif


#ifdef BLACK_HOLES
  /***** black hole accretion and feedback *****/
  t0 = second();
  blackhole();
  t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    printf("EXTRA TIMER: BH took %g sec\n", timediff(t0, t1));
#endif
#endif

#ifdef LMB_SPECTRAL_CRs
  t0 = second();
  evolve_spectral_crs();
  t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    printf("EXTRA TIMER: CRs took %g sec\n", timediff(t0, t1));
#endif

#ifdef LMB_SPECTRAL_CRs_DIFFUSION
  if(All.CRDiffusion_Ti_endstep == All.Ti_Current)
    {
      t0 = second();
      spectral_cr_diffusion();
      t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
      if(ThisTask == 0)
	{
	  printf("EXTRA TIMER: CR diffusion took %g sec\n", timediff(t0, t1));
	}
#endif
    }
#endif // LMB_SPECTRAL_CRs_DIFFUSION

#endif



#ifdef ACC_RUN
  on_end_of_timestep();

#endif


#if defined(BLACK_HOLES) || defined(VARIABLE_WINDS)
#ifdef FOF
  /* this will find new black hole seed halos and/or assign host halo masses for the variable wind model */
  if(All.Time >= All.TimeNextOnTheFlyFoF)
    {
      if(All.BlackHoleSeedStarMassFraction > 0)
	fof_fof(-2);
      else
	fof_fof(-1);

      if(All.ComovingIntegrationOn)
	All.TimeNextOnTheFlyFoF *= All.TimeBetOnTheFlyFoF;
      else
	All.TimeNextOnTheFlyFoF += All.TimeBetOnTheFlyFoF;
    }
#endif
#endif

#ifdef COOLING	/**** radiative cooling and star formation *****/

#if GADGET_SFR == SFR_NONE
  t0 = second();
  cl_CoolingOnly();
  t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    printf("EXTRA TIMER: cooling only took %g sec\n", timediff(t0, t1));
#endif

#endif

#ifdef GM_MUPPI
  if(ThisTask == 0)
    {
      printf("Start MUPPI normalization..\n");
      fflush(stdout);
    }
  muppinorm();
  if(ThisTask == 0)
    {
      printf("..MUPPI normalization done.\n");
      fflush(stdout);
    }
#endif


#ifdef SFR
  t0 = second();
  cooling_and_starformation();
  t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    printf("EXTRA TIMER: cooling & sfr took %g sec\n", timediff(t0, t1));
#endif

#if defined(GM_MUPPI)
  if(ThisTask == 0)
    {
      printf("Start MUPPI energy distribution..\n");
      fflush(stdout);
    }
  thermalenergy();
  if(ThisTask == 0)
    {
      printf("..MUPPI energy redistribution done.\n");
      fflush(stdout);
    }
#endif

#else
  cl_CoolingOnly();
#endif

  CPU_Step[CPU_COOLINGSFR] += measure_time();
#endif /*ends COOLING */


#ifdef MOL_CLOUDS
  if(ThisTask == 0)
    {
      printf("Starting MOL_CLOUDS...\n");
      fflush(stdout);
    }

  CPU_Step[CPU_MISC] += measure_time();

  do_mol_clouds();

  if(ThisTask == 0)
    {
      printf("done.\n");
      fflush(stdout);
    }
#endif

#ifdef SCF_HYBRID
  SCF_do_center_of_mass_correction(0.75, 10.0 * SCF_HQ_A, 0.01, 1000);
#endif

#ifdef AA_RELAX_ATMOSPHERE
  ++relaxCounter;

  if(relaxCounter > AA_RELAX_ATMOSPHERE)
    {
      if(ThisTask == 0)
	{
	  printf("Resetting internal energy to atmospheric profile\n");
	}

      const double u0 = 2442715.9;
      const double y0 = 4000.0;
      int i;
      double a3inv;

      if(All.ComovingIntegrationOn)
	{
	  a3inv = 1 / (All.Time * All.Time * All.Time);
	}
      else
	{
	  a3inv = 1.0;
	}

      for(i = 0; i < N_gas; i++)
	{
	  SphP[i].Entropy =
	    GAMMA_MINUS1 / pow(SphP[i].Density * a3inv, GAMMA_MINUS1) * u0 * (1 - P[i].Pos[1] / y0);
	  SphP[i].EntropyPred = SphP[i].Entropy;

#ifdef AA_RELAX_ATMOSPHERE_V
	  P[i].Vel[0] = P[i].Vel[1] = P[i].Vel[2] = 0.0;
	  SphP[i].VelPred[0] = SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0.0;
#endif
	}


      relaxCounter = 0;
    }
#endif

  double tend = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: non standard physics took %g sec\n", timediff(tstart, tend));
    }
#endif

  POP_T();
}


void compute_statistics(void)
{
  double t0 = second();

  if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
    {
#ifdef COMPUTE_POTENTIAL_ENERGY
      compute_potential();
#endif
      energy_statistics();	/* compute and output energy statistics */

#ifdef SCFPOTENTIAL
      SCF_write(0);
#endif

      All.TimeLastStatistics += All.TimeBetStatistics;
    }

  double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: statistics took %g sec\n", timediff(t0, t1));
    }
#endif
}


void execute_resubmit_command(void)
{
  char buf[1000];
  sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
  system(buf);
#endif
}


void create_snapshot_if_desired(void)
{
#ifdef SYNCRONIZ_OUTPUT
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
      {
#ifdef OUTPUTPOTENTIAL
#if !defined(EVALPOTENTIAL) || (defined(EVALPOTENTIAL) && defined(RECOMPUTE_POTENTIAL_ON_OUTPUT))
	domain_Decomposition(0, 0);

	compute_potential();
#endif
#endif

#ifdef OUTPUT_LIGHTCONES
	if(ThisTask == 0)
	  printf("create_snapshot_if_desired\n");
	lightcone_output(All.SnapshotFileCount++);
#else
	savepositions(All.SnapshotFileCount++);
#endif
	All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
      }
#endif
}


/*! This function finds the next synchronization point of the system
 * (i.e. the earliest point of time any of the particles needs a force
 * computation), and drifts the system to this point of time.  If the
 * system dirfts over the desired time of a snapshot file, the
 * function will drift to this moment, generate an output, and then
 * resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  double t0 = second();
  PUSH_T("drift");

  int n, i, prev;
  integertime dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;
  int highest_active_bin, highest_occupied_bin;
  double timeold;

  CPU_Step[CPU_MISC] += measure_time();

  timeold = All.Time;

  All.NumCurrentTiStep++;	/* we are now moving to the next sync point */

  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE, highest_occupied_bin = 0; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      highest_occupied_bin = n;
	      dt_bin = (((integertime) 1) << n);	// = 2^(n+1)
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }
	  else
	    {
	      dt_bin = 0;
	      ti_next_for_bin = All.Ti_Current;
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &ti_next_kick, &ti_next_kick_global);
#else
  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MYMPI_COMM_WORLD);
#endif


#ifndef SYNCRONIZ_OUTPUT
  if(ThisTask == 0)
    printf("ti_next_kick_global = %i, All.Ti_nextoutput = %i\n", ti_next_kick_global, All.Ti_nextoutput);
  while(ti_next_kick_global >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
      All.Ti_Current = All.Ti_nextoutput;

      if(All.ComovingIntegrationOn)
	All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      else
	All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

      set_cosmo_factors_for_current_time();

#ifdef TIMEDEPGRAV
      All.G = All.Gini * dGfak(All.Time);
#endif

      move_particles(All.Ti_nextoutput);

      CPU_Step[CPU_DRIFT] += measure_time();

#ifdef OUTPUTPOTENTIAL
#if !defined(EVALPOTENTIAL) || (defined(EVALPOTENTIAL) && defined(RECOMPUTE_POTENTIAL_ON_OUTPUT))
      domain_Decomposition(0, 0);

      compute_potential();
#endif
#endif

#ifdef LT_SEvDbg
      get_metals_sumcheck(9);
#endif

#ifdef OUTPUT_LIGHTCONES

      if(ThisTask == 0)
	printf("find_next_sync_point_and_drift\n");

      lightcone_output(All.SnapshotFileCount++);
#else
      savepositions(All.SnapshotFileCount++);	/* write snapshot file */
#endif

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }
#endif
  CPU_Step[CPU_MISC] += measure_time();



  All.Previous_Ti_Current = All.Ti_Current;
  All.Ti_Current = ti_next_kick_global;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  set_cosmo_factors_for_current_time();

#ifdef TIMEDEPGRAV
  All.G = All.Gini * dGfak(All.Time);
#endif

  All.TimeStep = All.Time - timeold;


#ifdef WINDTUNNEL
  All.WindCurrentX += All.WindVel * All.TimeStep;
#endif

#ifdef LT_STELLAREVOLUTION
  All.Time_Age = get_age(All.Time);
#endif


  /* mark the bins that will be active */
  for(n = 1, TimeBinActive[0] = 1, NumForceUpdate = TimeBinCount[0], highest_active_bin = 0; n < TIMEBINS;
      n++)
    {
      dt_bin = (((integertime) 1) << n);
      if((ti_next_kick_global % dt_bin) == 0)
	{
	  TimeBinActive[n] = 1;
	  NumForceUpdate += TimeBinCount[n];
	  if(TimeBinCount[n])
	    highest_active_bin = n;
	}
      else
	TimeBinActive[n] = 0;
    }

  sumup_large_ints(1, &NumForceUpdate, &GlobNumForceUpdate);
  MPI_Allreduce(&highest_active_bin, &All.HighestActiveTimeBin, 1, MPI_INT, MPI_MAX, MYMPI_COMM_WORLD);
  MPI_Allreduce(&highest_occupied_bin, &All.HighestOccupiedTimeBin, 1, MPI_INT, MPI_MAX, MYMPI_COMM_WORLD);

  if(GlobNumForceUpdate == All.TotNumPart)
    {
      Flag_FullStep = 1;
      if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
	PANIC("Something is wrong with the time bins.\n");
    }
  else
    Flag_FullStep = 0;


  /* move the new set of active/synchronized particles */
  /* Note: We do not yet call make_list_of_active_particles(), since we
   * may still need to old list in the dynamic tree update
   */
#ifndef KD_DRIFT_ALL_PARTICLES
  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	    drift_particle(i, All.Ti_Current);
	}
    }
#else
#pragma omp parallel for
  for(int i = 0; i < NumPart; i++)
    drift_particle(i, All.Ti_Current);
#pragma omp parallel for
  for(int no = All.MaxPart; no < All.MaxPart + Numnodestree; no++)
    force_drift_node(no, All.Ti_Current);
#endif

  CPU_Step[CPU_DRIFT] += measure_time();

  POP_T();
  double t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: find next and drift took %g sec\n", timediff(t0, t1));
    }
#endif

}


void make_list_of_active_particles(void)
{
  double t0 = second();
  PUSH_T("make list active part");

  int i, n, prev;
  /* make a link list with the particles in the active time bins */
  FirstActiveParticle = -1;

  CPU_Step[CPU_MISC] += measure_time();

  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	    {
#ifdef BLACK_HOLES
	      if(P[i].Mass == 0)
		continue;
#endif
	      if(prev == -1)
		FirstActiveParticle = i;

	      if(prev >= 0)
		NextActiveParticle[prev] = i;

	      prev = i;
	    }
	}
    }

  if(prev >= 0)
    NextActiveParticle[prev] = -1;

  CPU_Step[CPU_TREEMISC] += measure_time();

#ifdef _OPENMP
  build_active_particle_list();
#endif

  POP_T();
  double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: make list of active particles took %g sec\n", timediff(t0, t1));
    }
#endif

}

#ifdef _OPENMP
void build_active_particle_list()
{
  NActivePart = 0;
  for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      ActiveParticleList[NActivePart] = i;
      NActivePart++;
    }
}
#endif

/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
integertime find_next_outputtime(integertime ti_curr)
{
  int i, iter = 0;
  integertime ti, ti_next;
  double next, time;

  DumpFlag = 1;
  ti_next = -1;


  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn)
		ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
	      else
		ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

#ifdef PROCESS_TIMES_OF_OUTPUTLIST
	      /* first, determine maximum output interval based on All.MaxSizeTimestep */
	      integertime timax = (integertime) (All.MaxSizeTimestep / All.Timebase_interval);

	      /* make it a power 2 subdivision */
	      integertime ti_min = TIMEBASE;
	      while(ti_min > timax)
		ti_min >>= 1;
	      timax = ti_min;

	      double multiplier = ti / ((double) timax);

	      /* now round this to the nearest multiple of timax */
	      ti = ((integertime) (multiplier + 0.5)) * timax;
#endif

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount)
			All.SnapshotFileCount = i;
		    }

		  if(ti_next > ti)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount)
			All.SnapshotFileCount = i;
		    }
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      time = All.TimeOfFirstSnapshot;

      iter = 0;

      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(110);
	    }
	}
      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn)
	    ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
	  else
	    ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }


  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g  (DumpFlag=%d)\n\n", next, DumpFlag);

    }

  return ti_next;
}




/*! This routine writes for every synchronisation point in the timeline information to two log-files:
 * In FdInfo, we just list the timesteps that have been done, while in
 * FdTimebins we inform about the distribution of particles over the timebins, and which timebins are active on this step.
 * code is stored.
 */
void output_log_messages(void)
{
  double t0 = second();

  double z;
  int i, j;
  long long tot, tot_sph;
  long long tot_count[TIMEBINS];
  long long tot_count_sph[TIMEBINS];
  long long tot_cumulative[TIMEBINS];
  int weight, corr_weight;
  double sum, avg_CPU_TimeBin[TIMEBINS], frac_CPU_TimeBin[TIMEBINS];

#ifdef GM_MUPPI
  double Mstar, Xstar, Ystar, Zstar;
  double totMstar, totXstar, totYstar, totZstar;
  int Nstar, totNstar;
#endif

  sumup_large_ints(TIMEBINS, TimeBinCount, tot_count);
  sumup_large_ints(TIMEBINS, TimeBinCountSph, tot_count_sph);

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nSync-Point %d, Time: %g, Redshift: %g, Nf = %d%09d, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z,
		  (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  printf("\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
		 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z, All.TimeStep,
		  log(All.Time) - log(All.Time - All.TimeStep));
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nSync-Point %d, Time: %g, Nf = %d%09d, Systemstep: %g\n", All.NumCurrentTiStep,
		  All.Time, (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep);
	  printf("\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
	  fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
		  All.TimeStep);
	  fflush(FdInfo);
	}

      for(i = 1, tot_cumulative[0] = tot_count[0]; i < TIMEBINS; i++)
	tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];


      for(i = 0; i < TIMEBINS; i++)
	{
	  for(j = 0, sum = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++)
	    sum += All.CPU_TimeBinMeasurements[i][j];
	  if(All.CPU_TimeBinCountMeasurements[i])
	    avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];
	  else
	    avg_CPU_TimeBin[i] = 0;
	}

      for(i = All.HighestOccupiedTimeBin, weight = 1, sum = 0; i >= 0 && tot_count[i] > 0; i--, weight *= 2)
	{
	  if(weight > 1)
	    corr_weight = weight / 2;
	  else
	    corr_weight = weight;

	  frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
	  sum += frac_CPU_TimeBin[i];
	}

      for(i = All.HighestOccupiedTimeBin; i >= 0 && tot_count[i] > 0; i--)
	{
	  if(sum)
	    frac_CPU_TimeBin[i] /= sum;
	}


      printf
	("Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
      fprintf(FdTimebin,
	      "Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
      for(i = TIMEBINS - 1, tot = tot_sph = 0; i >= 0; i--)
	if(tot_count_sph[i] > 0 || tot_count[i] > 0)
	  {
	    printf(" %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%",
		   TimeBinActive[i] ? 'X' : ' ',
		   i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
		   i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
		   (i == All.HighestActiveTimeBin) ? '<' : ' ',
		   (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
		   avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);
	    fprintf(FdTimebin,
		    " %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%",
		    TimeBinActive[i] ? 'X' : ' ', i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
		    i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
		    (i == All.HighestActiveTimeBin) ? '<' : ' ',
		    (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
		    avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);


#ifdef DOMAIN_DECOMPOSITION_ON_TB
	    printf(" %c%c%c", AllDerived.DoDomainDecompositionTBFull == i ? 'D' : ' ',
		   AllDerived.DoDomainDecompositionTBDiffuse == i ? 'd' : ' ',
		   AllDerived.DoDomainDecompositionTBFix == i ? 'f' : ' ');
	    fprintf(FdTimebin, " %c%c%c", AllDerived.DoDomainDecompositionTBFull == i ? 'D' : ' ',
		    AllDerived.DoDomainDecompositionTBDiffuse == i ? 'd' : ' ',
		    AllDerived.DoDomainDecompositionTBFix == i ? 'f' : ' ');
#endif


	    printf("\n");
	    fprintf(FdTimebin, "\n");
	    if(TimeBinActive[i])
	      {
		tot += tot_count[i];
		tot_sph += tot_count_sph[i];
	      }
	  }
      printf("               ------------------------\n");
      fprintf(FdTimebin, "               ------------------------\n");

#ifdef PMGRID
      if(All.PM_Ti_endstep == All.Ti_Current)
	{
	  printf("PM-Step. Total: %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
	  fprintf(FdTimebin, "PM-Step. Total: %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
	}
      else
#endif
	{
	  printf("Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
	  fprintf(FdTimebin, "Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);

	}
      fprintf(FdTimebin, "\n");
      fflush(FdTimebin);

    }


#ifdef GM_MUPPI
  /* GM: we must decide if to leave this, or drop it */
  Mstar = Xstar = Ystar = Zstar = 0.0;
  totXstar = totYstar = totZstar = totMstar = 0.0;
  Nstar = totNstar = 0;
  All.Nstar = 0;
  for(i = 1; i < NumPart + 1; i++)
    if(P[i].Type == 4)
      {
	Xstar += P[i].Pos[0] * P[i].Mass;
	Ystar += P[i].Pos[1] * P[i].Mass;
	Zstar += P[i].Pos[2] * P[i].Mass;
	Mstar += P[i].Mass;
	Nstar++;
      }
  MPI_Allreduce(&Xstar, &totXstar, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Ystar, &totYstar, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Zstar, &totZstar, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Mstar, &totMstar, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Nstar, &totNstar, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  if(totNstar > 0)
    {
      All.StarCM[0] = totXstar / totMstar;
      All.StarCM[1] = totYstar / totMstar;
      All.StarCM[2] = totZstar / totMstar;
      All.Mstar = totMstar;
      All.Nstar = totNstar;
    }
  else
    {
      All.StarCM[0] = 0.0;
      All.StarCM[1] = 0.0;
      All.StarCM[2] = 0.0;
      All.Mstar = 0.0;
      All.Nstar = 0;
    }

  if(ThisTask == 0)
    {
      printf(" \n ===Star summary:\n");
      printf(" %d star particles for a mass of %g Msol\n", All.Nstar,
	     All.Mstar * All.UnitMass_in_g / SOLAR_MASS);
      printf(" CM position: %f %f %f\n\n", All.StarCM[0], All.StarCM[1], All.StarCM[2]);
    }
#endif


  output_extra_log_messages();

  double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: output log took %g sec\n", timediff(t0, t1));
    }
#endif

}




void write_cpu_log(void)
{
  double tstart = second();
  double max_CPU_Step[CPU_PARTS], avg_CPU_Step[CPU_PARTS], t0, t1, tsum;
  int i;

  CPU_Step[CPU_MISC] += measure_time();

  for(i = 1, CPU_Step[0] = 0; i < CPU_PARTS4SUM; i++)
    CPU_Step[0] += CPU_Step[i];

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	avg_CPU_Step[i] /= NTask;

      put_symbol(0.0, 1.0, '#');

      for(i = 1, tsum = 0.0; i < CPU_PARTS4SUM; i++)
	{
	  if(max_CPU_Step[i] > 0)
	    {
	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_Symbol[i]);
	      tsum += t1 - t0;

	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_SymbolImbalance[i]);
	      tsum += t1 - t0;
	    }
	}

      put_symbol(tsum / max_CPU_Step[0], 1.0, '-');

      fprintf(FdBalance, "Step=%7d  sec=%10.3f  Nf=%2d%09d  %s\n", All.NumCurrentTiStep, max_CPU_Step[0],
	      (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), CPU_String);
      fflush(FdBalance);

#ifdef KD_MONITOR_PERFORMANCE
      if(GlobNumForceUpdate < 10000)
	{
	  if(max_CPU_Step[0] > KD_MONITOR_PERFORMANCE)
	    {
	      count_performance++;
	      if(count_performance > 10)
		{
		  printf("KD PERFORMANCE WARNING: %d: %d %f \n", count_performance, (int) GlobNumForceUpdate,
			 max_CPU_Step[0]);
		  FILE *this_fd;
		  this_fd = fopen("stop", "w");
		  fclose(this_fd);
		}
	    }
	  else
	    count_performance = 0;
	}
#endif

      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
	{
	  All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
	  memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0],
		  &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1],
		  (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
	}

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements
							    [All.HighestActiveTimeBin]++] = max_CPU_Step[0];
    }

  CPUThisRun += CPU_Step[0];

  for(i = 0; i < CPU_PARTS; i++)
    CPU_Step[i] = 0;

  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	All.CPU_Sum[i] += avg_CPU_Step[i];

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);
      fprintf(FdCPU,
	      "total                      %10.2f  %5.1f%%\n"
	      "treegrav                   %10.2f  %5.1f%%\n"
	      "   treebuild               %10.2f  %5.1f%%\n"
	      "   treeupdate              %10.2f  %5.1f%%\n"
	      "   treewalk                %10.2f  %5.1f%%\n"
	      "   treecomm                %10.2f  %5.1f%%\n"
	      "   treeimbal               %10.2f  %5.1f%%\n"
	      "pmgrav                     %10.2f  %5.1f%%\n"
	      "sph                        %10.2f  %5.1f%%\n"
	      "   density                 %10.2f  %5.1f%%\n"
	      "   denscomm                %10.2f  %5.1f%%\n"
	      "   densimbal               %10.2f  %5.1f%%\n"
	      "   hydrofrc                %10.2f  %5.1f%%\n"
	      "   hydcomm                 %10.2f  %5.1f%%\n"
	      "   hydmisc                 %10.2f  %5.1f%%\n"
	      "   hydnetwork              %10.2f  %5.1f%%\n"
	      "   hydimbal                %10.2f  %5.1f%%\n" "   hmaxupdate              %10.2f  %5.1f%%\n"
#if GADGET_HYDRO == HYDRO_MFM
	      "mfm                        %10.2f  %5.1f%%\n"
	      "   gradients               %10.2f  %5.1f%%\n"
	      "      mfmcompute           %10.2f  %5.1f%%\n"
	      "      mfmcomm              %10.2f  %5.1f%%\n"
	      "      mfmmisc              %10.2f  %5.1f%%\n"
	      "      mfmimbal             %10.2f  %5.1f%%\n"
	      "   limiter                 %10.2f  %5.1f%%\n"
	      "      mfmcompute           %10.2f  %5.1f%%\n"
	      "      mfmcomm              %10.2f  %5.1f%%\n"
	      "      mfmmisc              %10.2f  %5.1f%%\n"
	      "      mfmimbal             %10.2f  %5.1f%%\n"
	      "   fluxes                  %10.2f  %5.1f%%\n"
	      "      mfmcompute           %10.2f  %5.1f%%\n"
	      "      mfmcomm              %10.2f  %5.1f%%\n"
	      "      mfmmisc              %10.2f  %5.1f%%\n" "      mfmimbal             %10.2f  %5.1f%%\n"
#endif
	      "domain                     %10.2f  %5.1f%%\n"
	      "potential                  %10.2f  %5.1f%%\n"
	      "predict                    %10.2f  %5.1f%%\n"
	      "kicks                      %10.2f  %5.1f%%\n"
	      "i/o                        %10.2f  %5.1f%%\n"
	      "peano                      %10.2f  %5.1f%%\n" "sfrcool                    %10.2f  %5.1f%%\n"
#ifdef GM_MUPPI
	      "    MUPPI     %10.2f  %5.1f%%\n" "    MP comm   %10.2f  %5.1f%%\n"
#endif
	      "blackholes                 %10.2f  %5.1f%%\n" "fof/subfind                %10.2f  %5.1f%%\n"
#ifdef SUBTIMERS
	      "  fof                      %10.2f  %5.1f%%\n"
	      "  subfind                  %10.2f  %5.1f%%\n"
	      "    treebuild_species      %10.2f  %5.1f%%\n"
	      "    treebuild              %10.2f  %5.1f%%\n"
	      "    smoothlength           %10.2f  %5.1f%%\n"
	      "    density                %10.2f  %5.1f%%\n"
	      "    DM_density             %10.2f  %5.1f%%\n"
	      "    save_density           %10.2f  %5.1f%%\n"
	      "    exchange               %10.2f  %5.1f%%\n"
	      "    proc_collective_halos  %10.2f  %5.1f%%\n"
	      "    sort_local             %10.2f  %5.1f%%\n"
	      "    proc_local_groups      %10.2f  %5.1f%%\n"
	      "    unsort_local           %10.2f  %5.1f%%\n"
	      "    exchange_return        %10.2f  %5.1f%%\n"
	      "    domain                 %10.2f  %5.1f%%\n"
	      "    overdensity_masses     %10.2f  %5.1f%%\n"
	      "    contamination          %10.2f  %5.1f%%\n" "    save_final             %10.2f  %5.1f%%\n"
#endif
	      "smoothing                  %10.2f  %5.1f%%\n"
#ifdef CONDUCTION
	      "conduction                 %10.2f  %5.1f%%\n"
#endif
#ifdef ADAPTGRAVSOFT
	      "ags_dens                   %10.2f  %5.1f%%\n"
	      "ags_d_comm                 %10.2f  %5.1f%%\n"
	      "ags_d_imbal                %10.2f  %5.1f%%\n"
	      "ags_d_misc                 %10.2f  %5.1f%%\n" "ags_hmaxupd                %10.2f  %5.1f%%\n"
#endif
#if defined(fSIDM) || defined(rSIDM)
	      "mSIDM                      %10.2f  %5.1f%%\n"
	      "   comp                    %10.2f  %5.1f%%\n"
	      "      local                %10.2f  %5.1f%%\n"
	      "      parallel             %10.2f  %5.1f%%\n"
	      "      frequent scatter     %10.2f  %5.1f%%\n"
	      "         drag force        %10.2f  %5.1f%%\n"
	      "         random            %10.2f  %5.1f%%\n"
	      "      rare scatter         %10.2f  %5.1f%%\n"
	      "   prep                    %10.2f  %5.1f%%\n"
	      "   comm                    %10.2f  %5.1f%%\n"
	      "   wait                    %10.2f  %5.1f%%\n"
	      "   netw                    %10.2f  %5.1f%%\n" "   misc                    %10.2f  %5.1f%%\n"
#endif
#ifdef GM_STARDENSITY
	      "misc                       %10.2f  %5.1f%%\n"
	      "star density               %10.2f  %5.1f%%\n"
	      "    star density compute   %10.2f  %5.1f%%\n"
	      "    star density wait      %10.2f  %5.1f%%\n" "    star density comm      %10.2f  %5.1f%%\n",
#else
	      "misc                       %10.2f  %5.1f%%\n",
#endif
	      All.CPU_Sum[CPU_ALL], 100.0,
	      All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]
	      + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
	      + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]
	      + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE]
	      + All.CPU_Sum[CPU_TREEMISC],
	      (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]
	       + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
	       + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]
	       + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE]
	       + All.CPU_Sum[CPU_TREEMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEBUILD],
	      (All.CPU_Sum[CPU_TREEBUILD]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEUPDATE],
	      (All.CPU_Sum[CPU_TREEUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2],
	      (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV],
	      (All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2],
	      (All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MESH],
	      (All.CPU_Sum[CPU_MESH]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT]
	      + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
	      + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
	      + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC] + All.CPU_Sum[CPU_HYDNETWORK],
	      (All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT]
	       + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
	       + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
	       + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC] +
	       All.CPU_Sum[CPU_HYDNETWORK]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DENSCOMPUTE],
	      (All.CPU_Sum[CPU_DENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DENSCOMM],
	      (All.CPU_Sum[CPU_DENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DENSWAIT],
	      (All.CPU_Sum[CPU_DENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDCOMPUTE],
	      (All.CPU_Sum[CPU_HYDCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDCOMM],
	      (All.CPU_Sum[CPU_HYDCOMM]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDMISC],
	      (All.CPU_Sum[CPU_HYDMISC]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDNETWORK],
	      (All.CPU_Sum[CPU_HYDNETWORK]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_HYDWAIT],
	      (All.CPU_Sum[CPU_HYDWAIT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_TREEHMAXUPDATE],
	      (All.CPU_Sum[CPU_TREEHMAXUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
#if GADGET_HYDRO == HYDRO_MFM
	      All.CPU_Sum[CPU_MFMGRADCOMPUTE] + All.CPU_Sum[CPU_MFMGRADCOMM]
	      + All.CPU_Sum[CPU_MFMGRADMISC] + All.CPU_Sum[CPU_MFMGRADWAIT]
	      + All.CPU_Sum[CPU_MFMLIMCOMPUTE] + All.CPU_Sum[CPU_MFMLIMCOMM]
	      + All.CPU_Sum[CPU_MFMLIMMISC] + All.CPU_Sum[CPU_MFMLIMWAIT]
	      + All.CPU_Sum[CPU_MFMFLUXCOMPUTE] + All.CPU_Sum[CPU_MFMFLUXCOMM]
	      + All.CPU_Sum[CPU_MFMFLUXMISC] + All.CPU_Sum[CPU_MFMFLUXWAIT],
	      (All.CPU_Sum[CPU_MFMGRADCOMPUTE] + All.CPU_Sum[CPU_MFMGRADCOMM]
	       + All.CPU_Sum[CPU_MFMGRADMISC] + All.CPU_Sum[CPU_MFMGRADWAIT]
	       + All.CPU_Sum[CPU_MFMLIMCOMPUTE] + All.CPU_Sum[CPU_MFMLIMCOMM]
	       + All.CPU_Sum[CPU_MFMLIMMISC] + All.CPU_Sum[CPU_MFMLIMWAIT]
	       + All.CPU_Sum[CPU_MFMFLUXCOMPUTE] + All.CPU_Sum[CPU_MFMFLUXCOMM]
	       + All.CPU_Sum[CPU_MFMFLUXMISC] + All.CPU_Sum[CPU_MFMFLUXWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMGRADCOMPUTE] + All.CPU_Sum[CPU_MFMGRADCOMM]
	      + All.CPU_Sum[CPU_MFMGRADMISC] + All.CPU_Sum[CPU_MFMGRADWAIT],
	      (All.CPU_Sum[CPU_MFMGRADCOMPUTE] + All.CPU_Sum[CPU_MFMGRADCOMM]
	       + All.CPU_Sum[CPU_MFMGRADMISC] + All.CPU_Sum[CPU_MFMGRADWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMGRADCOMPUTE],
	      (All.CPU_Sum[CPU_MFMGRADCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMGRADCOMM],
	      (All.CPU_Sum[CPU_MFMGRADCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMGRADMISC],
	      (All.CPU_Sum[CPU_MFMGRADMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMGRADWAIT],
	      (All.CPU_Sum[CPU_MFMGRADWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMLIMCOMPUTE] + All.CPU_Sum[CPU_MFMLIMCOMM]
	      + All.CPU_Sum[CPU_MFMLIMMISC] + All.CPU_Sum[CPU_MFMLIMWAIT],
	      (All.CPU_Sum[CPU_MFMLIMCOMPUTE] + All.CPU_Sum[CPU_MFMLIMCOMM]
	       + All.CPU_Sum[CPU_MFMLIMMISC] + All.CPU_Sum[CPU_MFMLIMWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMLIMCOMPUTE],
	      (All.CPU_Sum[CPU_MFMLIMCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMLIMCOMM],
	      (All.CPU_Sum[CPU_MFMLIMCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMLIMMISC],
	      (All.CPU_Sum[CPU_MFMLIMMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMLIMWAIT],
	      (All.CPU_Sum[CPU_MFMLIMWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMFLUXCOMPUTE] + All.CPU_Sum[CPU_MFMFLUXCOMM]
	      + All.CPU_Sum[CPU_MFMFLUXMISC] + All.CPU_Sum[CPU_MFMFLUXWAIT],
	      (All.CPU_Sum[CPU_MFMFLUXCOMPUTE] + All.CPU_Sum[CPU_MFMFLUXCOMM]
	       + All.CPU_Sum[CPU_MFMFLUXMISC] + All.CPU_Sum[CPU_MFMFLUXWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMFLUXCOMPUTE],
	      (All.CPU_Sum[CPU_MFMFLUXCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMFLUXCOMM],
	      (All.CPU_Sum[CPU_MFMFLUXCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMFLUXMISC],
	      (All.CPU_Sum[CPU_MFMFLUXMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MFMFLUXWAIT], (All.CPU_Sum[CPU_MFMFLUXWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
	      All.CPU_Sum[CPU_DOMAIN],
	      (All.CPU_Sum[CPU_DOMAIN]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_POTENTIAL],
	      (All.CPU_Sum[CPU_POTENTIAL]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_DRIFT],
	      (All.CPU_Sum[CPU_DRIFT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_TIMELINE],
	      (All.CPU_Sum[CPU_TIMELINE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_SNAPSHOT],
	      (All.CPU_Sum[CPU_SNAPSHOT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_PEANO],
	      (All.CPU_Sum[CPU_PEANO]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_COOLINGSFR],
	      (All.CPU_Sum[CPU_COOLINGSFR]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef GM_MUPPI
	      All.CPU_Sum[CPU_MUPPI],
	      (All.CPU_Sum[CPU_MUPPI]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MUPPINET], (All.CPU_Sum[CPU_MUPPINET]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
	      All.CPU_Sum[CPU_BLACKHOLES],
	      (All.CPU_Sum[CPU_BLACKHOLES]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_FOF],
	      (All.CPU_Sum[CPU_FOF]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef SUBTIMERS
	      All.CPU_Sum[CPU_FOF_TOT], (All.CPU_Sum[CPU_FOF_TOT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_TOT], (All.CPU_Sum[CPU_SUBFIND_TOT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_TREEBUILD_SPECIES],
	      (All.CPU_Sum[CPU_SUBFIND_TREEBUILD_SPECIES]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_TREEBUILD],
	      (All.CPU_Sum[CPU_SUBFIND_TREEBUILD]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_SMOOTHINGLENGTH],
	      (All.CPU_Sum[CPU_SUBFIND_SMOOTHINGLENGTH]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_DENSITY],
	      (All.CPU_Sum[CPU_SUBFIND_DENSITY]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_DMDENSITY],
	      (All.CPU_Sum[CPU_SUBFIND_DMDENSITY]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_SAVE_DENSITY],
	      (All.CPU_Sum[CPU_SUBFIND_SAVE_DENSITY]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_EXCHANGE],
	      (All.CPU_Sum[CPU_SUBFIND_EXCHANGE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_COLLHALOS],
	      (All.CPU_Sum[CPU_SUBFIND_COLLHALOS]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_SORTLOCAL],
	      (All.CPU_Sum[CPU_SUBFIND_SORTLOCAL]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_LOCALGROUPS],
	      (All.CPU_Sum[CPU_SUBFIND_LOCALGROUPS]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_UNSORTLOCAL],
	      (All.CPU_Sum[CPU_SUBFIND_UNSORTLOCAL]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_EXCHANGERETURN],
	      (All.CPU_Sum[CPU_SUBFIND_EXCHANGERETURN]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_DOMAIN], (All.CPU_Sum[CPU_SUBFIND_DOMAIN]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_MASSES], (All.CPU_Sum[CPU_SUBFIND_MASSES]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_CONT], (All.CPU_Sum[CPU_SUBFIND_CONT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SUBFIND_SAVEFINAL],
	      (All.CPU_Sum[CPU_SUBFIND_SAVEFINAL]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
	      All.CPU_Sum[CPU_SMTHCOMPUTE] + All.CPU_Sum[CPU_SMTHWAIT] + All.CPU_Sum[CPU_SMTHCOMM] +
	      All.CPU_Sum[CPU_SMTHMISC],
	      (All.CPU_Sum[CPU_SMTHCOMPUTE] + All.CPU_Sum[CPU_SMTHWAIT] + All.CPU_Sum[CPU_SMTHCOMM] +
	       All.CPU_Sum[CPU_SMTHMISC]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef CONDUCTION
	      All.CPU_Sum[CPU_CONDUCTION], (All.CPU_Sum[CPU_CONDUCTION]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#ifdef ADAPTGRAVSOFT
	      All.CPU_Sum[CPU_AGSDENSCOMPUTE], (All.CPU_Sum[CPU_AGSDENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_AGSDENSCOMM], (All.CPU_Sum[CPU_AGSDENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_AGSDENSWAIT], (All.CPU_Sum[CPU_AGSDENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_AGSDENSMISC], (All.CPU_Sum[CPU_AGSDENSMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_AGSTREEHMAXUPD], (All.CPU_Sum[CPU_AGSTREEHMAXUPD]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
#if defined(fSIDM) || defined(rSIDM)
	      All.CPU_Sum[CPU_MSIDMCOMPUTE] + All.CPU_Sum[CPU_MSIDMWAIT] + All.CPU_Sum[CPU_MSIDMCOMM]
	      + All.CPU_Sum[CPU_MSIDMNETWORK] + All.CPU_Sum[CPU_MSIDMMISC],
	      (All.CPU_Sum[CPU_MSIDMCOMPUTE] + All.CPU_Sum[CPU_MSIDMWAIT] + All.CPU_Sum[CPU_MSIDMCOMM]
	       + All.CPU_Sum[CPU_MSIDMNETWORK] + All.CPU_Sum[CPU_MSIDMMISC]) / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMCOMPUTE],
	      All.CPU_Sum[CPU_MSIDMCOMPUTE] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMLOCAL],
	      All.CPU_Sum[CPU_MSIDMLOCAL] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMPARALLEL],
	      All.CPU_Sum[CPU_MSIDMPARALLEL] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_FSIDMSCATTER],
	      All.CPU_Sum[CPU_FSIDMSCATTER] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_FSIDMSCATTER_DF],
	      All.CPU_Sum[CPU_FSIDMSCATTER_DF] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_FSIDMSCATTER_RC],
	      All.CPU_Sum[CPU_FSIDMSCATTER_RC] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_RSIDMSCATTER],
	      All.CPU_Sum[CPU_RSIDMSCATTER] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMPREPP],
	      All.CPU_Sum[CPU_MSIDMPREPP] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMCOMM],
	      All.CPU_Sum[CPU_MSIDMCOMM] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMWAIT],
	      All.CPU_Sum[CPU_MSIDMWAIT] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMNETWORK],
	      All.CPU_Sum[CPU_MSIDMNETWORK] / All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_MSIDMMISC], All.CPU_Sum[CPU_MSIDMMISC] / All.CPU_Sum[CPU_ALL] * 100,
#endif
#ifdef GM_STARDENSITY
	      All.CPU_Sum[CPU_STARDENSMISC] + All.CPU_Sum[CPU_STARCOMPUTE] + All.CPU_Sum[CPU_STARWAIT] +
	      All.CPU_Sum[CPU_STARCOMM],
	      (All.CPU_Sum[CPU_STARDENSMISC] + All.CPU_Sum[CPU_STARCOMPUTE] + All.CPU_Sum[CPU_STARWAIT] +
	       All.CPU_Sum[CPU_STARCOMM]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_STARCOMPUTE],
	      (All.CPU_Sum[CPU_STARCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_STARWAIT],
	      (All.CPU_Sum[CPU_STARWAIT]) / All.CPU_Sum[CPU_ALL] * 100, All.CPU_Sum[CPU_STARCOMM],
	      (All.CPU_Sum[CPU_STARCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
	      All.CPU_Sum[CPU_MISC], (All.CPU_Sum[CPU_MISC]) / All.CPU_Sum[CPU_ALL] * 100);
      fprintf(FdCPU, "\n");
      fflush(FdCPU);
    }

  double tend = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  if(ThisTask == 0)
    {
      printf("EXTRA TIMER: cpu logs took %g sec\n", timediff(tstart, tend));
    }
#endif
}


void put_symbol(double t0, double t1, char c)
{
  int i, j;

  i = (int) (t0 * CPU_STRING_LEN + 0.5);
  j = (int) (t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j < 0)
    j = 0;
  if(i >= CPU_STRING_LEN)
    i = CPU_STRING_LEN;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    CPU_String[i++] = c;

  CPU_String[CPU_STRING_LEN] = 0;
}


/*! This routine first calls a computation of various global
 * quantities of the particle distribution, and then writes some
 * statistics about the energies in the various particle components to
 * the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

#if defined(MAGNETIC_STATISTICS)
      fprintf(FdEnergy, " %g", SysState.EnergyMag);
#if defined(TRACEDIVB)
      fprintf(FdEnergy, " %g", SysState.DivBerr);
#endif
#endif

      fprintf(FdEnergy, " \n");
      fflush(FdEnergy);
    }
}




void output_extra_log_messages(void)
{

  if(ThisTask == 0)
    {
#ifdef LT_STELLAREVOLUTION_NOT_YET_FIXED
      long long tot, tot_sph, tot_stars;
      long long tot_count[TIMEBINS], tot_count_sph[TIMEBINS], tot_count_stars[TIMEBINS];
      int i;

      sumup_large_ints(TIMEBINS, TimeBinCountStars, tot_count_stars);

      printf("Occupied timebins: non-sph         sph       stars    dt\n");
      for(i = TIMEBINS - 1, tot = tot_sph = tot_stars = 0; i >= 0; i--)
	if(tot_count_sph[i] > 0 || tot_count[i] > 0 || tot_count_stars[i] > 0)
	  {
	    printf(" %c  bin=%2d     %2d%09d %2d%09d %2d%09d   %6g\n",
		   TimeBinActive[i] ? 'X' : ' ',
		   i,
		   (int) ((tot_count[i] - tot_count_sph[i]) / 1000000000),
		   (int) ((tot_count[i] - tot_count_sph[i]) % 1000000000),
		   (int) (tot_count_sph[i] / 1000000000),
		   (int) (tot_count_sph[i] % 1000000000),
		   (int) (tot_count_stars[i] / 1000000000),
		   (int) (tot_count_stars[i] % 1000000000),
		   i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0);
	    if(TimeBinActive[i])
	      {
		tot += tot_count[i];
		tot_sph += tot_count_sph[i];
		tot_stars += tot_count_stars[i];
	      }
	  }
#endif


#ifdef CHEMISTRY
      printf("Abundances elec: %g, HM: %g, H2I: %g, H2II: %g\n",
	     SphP[1].elec, SphP[1].HM, SphP[1].H2I, SphP[1].H2II);
#endif

#ifdef XXLINFO
      if(Flag_FullStep == 1)
	{
	  fprintf(FdXXL, "%d %g ", All.NumCurrentTiStep, All.Time);
#ifdef MAGNETIC
	  fprintf(FdXXL, "%e ", MeanB);
#ifdef TRACEDIVB
	  fprintf(FdXXL, "%e ", MaxDivB);
#endif
#endif
#ifdef TIME_DEP_ART_VISC
	  fprintf(FdXXL, "%f ", MeanAlpha);
#endif
	  fprintf(FdXXL, "\n");
	  fflush(FdXXL);
	}
#endif

#ifdef DARKENERGY
      if(All.ComovingIntegrationOn == 1)
	{
	  double hubble_a;

	  hubble_a = hubble_function(All.Time);
	  fprintf(FdDE, "%d %g %e ", All.NumCurrentTiStep, All.Time, hubble_a);
#ifndef TIMEDEPDE
	  fprintf(FdDE, "%e ", All.DarkEnergyParam);
#else
	  fprintf(FdDE, "%e %e ", get_wa(All.Time), DarkEnergy_a(All.Time));
#endif
#ifdef TIMEDEPGRAV
	  fprintf(FdDE, "%e %e", dHfak(All.Time), dGfak(All.Time));
#endif
	  fprintf(FdDE, "\n");
	  fflush(FdDE);
	}
#endif
    }
}


void check_node_info(const char *func, const char *file, int linenr)
{
  int maxpart = All.MaxPart;

  MPI_Barrier(MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Checking node data (function %s in file %s at line %d) ...\n", func, file, linenr);

  for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      int no = Father[i];
      PANIC_IF(no < All.MaxPart, "error at check 1: i=%d, no=%d < All.MaxPart=%d\n", i, no, All.MaxPart);
      PANIC_IF(no >= maxpart + MaxNodes + 1, "error: at check 2: i=%d, no=%d >= maxpart=%d\n", i, no,
	       maxpart + MaxNodes + 1);

      while(no >= 0)
	{
	  if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* top-level tree-node reached */
	    break;

	  PANIC_IF(no == Nodes[no].u.d.father,
		   "ettot at check 3: Fixed point  no == Nodes[no].u.d.father  detected, %d/%d\n", no,
		   Nodes[no].u.d.father);
	  no = Nodes[no].u.d.father;
	  PANIC_IF(no < All.MaxPart, "error: at check 4: i=%d, no=%d < All.MaxPart=%d\n", i, no, All.MaxPart);
	  PANIC_IF(no >= maxpart + MaxNodes + 1, "error: at check 5: i=%d, no=%d >= maxpart=%d\n", i, no,
		   maxpart + MaxNodes + 1);
	}
    }
  if(ThisTask == 0)
    printf("Node check done ...\n");

  MPI_Barrier(MYMPI_COMM_WORLD);

}


void check_particles_info(const char *func, const char *file, int linenr)
{
  int i, k, vok = 0, pok = 0, vsph = 0, mok = 0, vstar = 0;
  double vv;

  MPI_Barrier(MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Checking particle data (function %s in file %s at line %d) ...\n", func, file, linenr);

#ifdef _OPENMP
#pragma omp parallel for private(k,vv) reduction(+:vok,pok,vsph,vstar)
#endif
  for(i = 0; i < NumPart; i++)
    {

      for(k = 0; k < 3; k++)
	{
	  if(P[i].Vel[k] > -1e8 && P[i].Vel[k] < 1e8)
	    {
	      vv = sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
	      if(vv > 15000)
		{
		  printf
		    ("task=%d: WARNING: Large velocity for particle %d ID %llu v[%d]=%g, renormalizing it !!\n",
		     ThisTask, i, (unsigned long long) P[i].ID, k, vv);
		  fflush(stdout);
		  P[i].Vel[0] = P[i].Vel[0] / vv * 10000;
		  P[i].Vel[1] = P[i].Vel[1] / vv * 10000;
		  P[i].Vel[2] = P[i].Vel[2] / vv * 10000;
		}
	      vok++;
	    }
	  else
	    {
	      printf
		("task=%d:  strange value in velocity in for particle %d ID %llu , type=%d, mass=%g, v[%d]=%g\n",
		 ThisTask, i, (unsigned long long) P[i].ID, P[i].Type, (float) P[i].Mass, k,
		 (float) P[i].Vel[k]);
	      if(P[i].Type == 0)
		printf("        vred[%d]=%g  a_hydro[%d]=%g\n", k, (float) SphP[i].VelPred[k], k,
		       (double) SphP[i].HydroAccel[k]);
	      fflush(stdout);
	      endrun(712401);
	    }

	  if(P[i].Pos[k] > -10000 && P[i].Pos[k] < All.BoxSize + 10000)
	    pok++;
	  else
	    {
	      printf("task=%d:  strange value in position in for particle %d ID %llu x[%d]=%g\n", ThisTask, i,
		     (unsigned long long) P[i].ID, k, (double) P[i].Pos[k]);
	      fflush(stdout);
	      endrun(712402);
	    }
	}

      double massDMpart;

      if(All.MassTable[1] > 0)
	massDMpart = All.MassTable[1];
      else
#ifdef BLACK_HOLES
	massDMpart = All.massDMpart;
#else
      if(All.MassTable[0] > 0)
	massDMpart = All.MassTable[0] * 10;
      else
	massDMpart = P[0].Mass * 10;
#endif

      if((P[i].Mass > massDMpart / 5000 && P[i].Mass < massDMpart * 10) || P[i].Type == 2 || P[i].Type == 3
	 || P[i].Type == 5)
	mok++;
      else
	{
	  printf("task=%d:  strange value in Mass in for particle %d ID %llu Type %d mass=%g\n", ThisTask, i,
		 (unsigned long long) P[i].ID, P[i].Type, (float) P[i].Mass);
	  fflush(stdout);
	  if(P[i].Mass > 0)
	    endrun(712403);
	}

      if(P[i].Type == 0)
	{
	  if(SphP[i].Density >= 0 && SphP[i].Density < 1e10)
	    vsph++;
	  else
	    printf("task=%d: Gas Particle id=%llu,m=%e strange Density value: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass, (double) SphP[i].Density);

	  if((SphP[i].EntropyPred > -1e20 && SphP[i].EntropyPred < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange predicted Entropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass, (double) SphP[i].EntropyPred);

	  if((SphP[i].Entropy > -1e20 && SphP[i].Entropy < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange Entropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass, (double) SphP[i].Entropy);

	  if((SphP[i].DtEntropy > -1e20 && SphP[i].DtEntropy < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange DtEntropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass, (double) SphP[i].DtEntropy);

#if GADGET_HYDRO == HYDRO_SPH
	  if((SphP[i].Pressure > -1e20 && SphP[i].Pressure < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu,m=%e strange Pressure value in hydro: %e,%e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass, (double) SphP[i].EntropyPred,
		   (float) SphP[i].Pressure);
#elif GADGET_HYDRO == HYDRO_PESPH
	  if((SphP[i].Pressure > -1e20 && SphP[i].Pressure < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu,m=%e strange Pressure value in hydro: %e,%e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass, (double) SphP[i].EntropyPred,
		   (float) SphP[i].Pressure);
#endif
	  if(SphP[i].HydroAccel[0] > -1e20 && SphP[i].HydroAccel[0] < 1e20 &&
	     SphP[i].HydroAccel[1] > -1e20 && SphP[i].HydroAccel[1] < 1e20 &&
	     SphP[i].HydroAccel[2] > -1e20 && SphP[i].HydroAccel[2] < 1e20)
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu,m=%e strange acceleration value in hydro: %e,%e,%e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass,
		   (double) SphP[i].HydroAccel[0], (double) SphP[i].HydroAccel[1],
		   (double) SphP[i].HydroAccel[2]);
	}

      if(P[i].Type > 5 || P[i].Type < 0)
	{
	  printf("task=%d:  P[i=%d].Type=%d\n", ThisTask, i, P[i].Type);
	  endrun(712411);
	}

#ifdef LT_STELLAREVOLUTION
      if(P[i].Type == 4)
	{
	  if(MetP[P[i].pt.MetID].PID != i)
	    {
	      printf("task=%d:  error in cross-indexes for star-particle %d ID %llu\n", ThisTask, i,
		     (unsigned long long) P[i].ID);
	      fflush(stdout);
	      endrun(712412);
	    }

	  if(MetP[P[i].pt.MetID].weight >= 0 && MetP[P[i].pt.MetID].weight < 1e3)
	    vstar++;
	  else
	    printf("task=%d: Star Particle id=%llu,m=%e strange weight for star: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, (float) P[i].Mass,
		   (float) MetP[P[i].pt.MetID].weight);
	}
#endif

#if defined(BLACK_HOLES)
      if(P[i].Type == 5)
	if(BHP[P[i].pt.BHID].PID != i)
	  {
	    printf("task=%d:  error in cross-indexes for bh-particle %d ID %llu\n", ThisTask, i,
		   (unsigned long long) P[i].ID);
	    fflush(stdout);
	    endrun(712413);
	  }
#endif
    }

  if(ThisTask == 0)
    printf
      ("Positions, Velocities, Masses and  SphVelocities fine for (%d,%d,%d,%d) of %d/%d/%d cases on task 0 (and %d stars)...\n\n",
       pok, vok, mok, vsph / 6, NumPart * 3, NumPart, N_gas, vstar);

  MPI_Barrier(MYMPI_COMM_WORLD);
}


int first_timestep = 0;


void on_beginning_of_timestep()
{
  //  if(ThisTask==0)  dump_memory_table();
#ifdef ACC_RUN
  VERBOSE(-1, "Alloc WalkLieP %d", NumPart);
  WalkLiteP = (struct gravity_walkdata *) mymalloc("WalkLiteP", NumPart * sizeof(struct gravity_walkdata));
#ifdef LT_STELLAREVOLUTION
  Acc_SFs_PhysDensThresh_table =
    (double *) mymalloc("Acc_SFs_PhysDensThresh_table", SFs_dim * ZBins * sizeof(double));
  for(int i = 0; i < SFs_dim; i++)
    {
      for(int j = 0; j < ZBins; j++)
	{
	  Acc_SFs_PhysDensThresh_table[i * ZBins + j] = SFs[i].PhysDensThresh[j];
	}
    }
#endif

#endif
}

void on_end_of_timestep()
{
#ifdef ACC_RUN
  if(acc_all.inside_data_region)
    {
      update_openacc_structures_exit();
    }


#ifdef LT_STELLAREVOLUTION
  myfree(Acc_SFs_PhysDensThresh_table);
#endif
  VERBOSE(-1, "Free WalkLieP %d", NumPart);
  myfree(WalkLiteP);

#endif
}
