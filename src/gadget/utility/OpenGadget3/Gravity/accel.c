
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



/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */

#ifdef DEBUG_PARTICLE_GRAVITY_ACCELERATION
void debug_particle(double Costtotal)
{


  /* for instance: add in Config.sh: DEBUG_PARTICLE_GRAVITY_ACCELERATION=717 and the particle 717 will be tracked */
  {

    int PARTICELLONA, mycounter;
    for(PARTICELLONA = 0; PARTICELLONA < NumPart; PARTICELLONA++)
      if(P[PARTICELLONA].ID == DEBUG_PARTICLE_GRAVITY_ACCELERATION)
	{
	  VERBOSE_ALL
	    (2,
	     "\n\t Debug:true, Time:%f, Task:%d, ID=%d, Acc:[%lg,%lg,%lg], CostTotal:%e, P.GravCost:%e, NActivePart:%d \n",
	     All.Time, ThisTask, DEBUG_PARTICLE_GRAVITY_ACCELERATION, P[PARTICELLONA].g.dGravAccel[0],
	     P[PARTICELLONA].g.dGravAccel[1], P[PARTICELLONA].g.dGravAccel[2], Costtotal,
	     P[PARTICELLONA].GravCost[TakeLevel], NActivePart);
	}

  }

}
#endif

void compute_grav_accelerations(void)
{
  double t0 = second();
  PUSH_T("compute grav acceleration");

  CPU_Step[CPU_MISC] += measure_time();
  VERBOSE(2, "Start gravity force computation...\n");






#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      long_range_force();
      CPU_Step[CPU_MESH] += measure_time();
    }
#endif


#ifdef ACC_RUN
  on_beginning_of_timestep();
#endif

#ifndef ONLY_PM

  gravity_tree();		/* computes gravity accel. on CPU */
#ifdef DEBUG_PARTICLE_GRAVITY_ACCELERATION
  debug_particle(0);
#endif

  /* For the first timestep, we redo it
   * to allow usage of relative opening
   * criterion for consistent accuracy.
   */



  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    {
      gravity_tree();		/* computes gravity accel. */
#ifdef DEBUG_PARTICLE_GRAVITY_ACCELERATION
      debug_particle(0);
#endif

    }
#endif


  if(All.Ti_Current == 0 && RestartFlag == 0 && header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
    second_order_ics();		/* produces the actual ICs from the special second order IC file */


#ifdef FORCETEST
  gravity_forcetest();
#endif

  VERBOSE(2, "gravity force computation done.\n");


  POP_T();
  double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  VERBOSE(5, "EXTRA TIMER: grav accel took %g sec\n", timediff(t0, t1));

#endif

}


void compute_densities(void)
{
  PUSH_T("density");

#if (defined(fSIDM) || defined(rSIDM)) && !defined(ADAPTGRAVSOFT)
  VERBOSE(2, "Start density computation for PartType == 1\n");
  pt1_density();		/* compute hsml for PartType==1, used for mSIDM */
  pt1_force_update_hmax();
#endif

  if(All.TotN_gas > 0)
    {
      VERBOSE(2, "Start density computation...\n");

      density();		/* computes density, and pressure */

#if (defined(SMOOTH_DIVB)  || defined(SMOOTH_ROTB) || defined(BSMOOTH) || defined(LT_USE_DENSITY_IN_WEIGHT) || defined(LT_SMOOTH_Z) || defined(LT_SMOOTH_XCLD) || defined(LT_TRACK_WINDS) || defined(BSMOOTH_TIME)) || defined(VSMOOTH)
      {
	double t0 = second();
	smoothed_values();
	double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
	VERBOSE(0, "EXTRA TIMER: density smoothing took %g sec\n", timediff(t0, t1));
#endif
      }
#endif

#ifdef LT_STELLAREVOLUTION
      {
	VERBOSE(2, "Start supernovae computation...\n");


	double t0 = second();
	evolve_SN();
	double t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
	VERBOSE(0, "EXTRA TIMER: SN evolve took %g sec\n", timediff(t0, t1));
#endif
      }
#endif

#ifdef GL_CR_DUST
      {
	VERBOSE(2, "Start dust computation...\n");

	double t0 = second();
	evolve_dust();
	double t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
	VERBOSE(0, "EXTRA TIMER: dust evolve took %g sec\n", timediff(t0, t1));
#endif
      }
#endif

#if defined(SNIA_HEATING)
      snIa_heating();
#endif

    }

#ifdef AXION_DM
  if(All.TotN_axion > 0)
    {
      VERBOSE(2, "Start axion density computation: first cycle...\n");


      ax_density(0);

      VERBOSE(2, "Start axion density computation: second cycle...\n");


      // MN: MPI_BARRIER

      ax_density(1);

    }
#endif

  POP_T();
}


#if GADGET_HYDRO == HYDRO_SPH || GADGET_HYDRO == HYDRO_PESPH
void compute_hydro_accelerations(void)
{
  PUSH_T("hydro accel");

  if(All.TotN_gas > 0)
    {
      VERBOSE(2, "Start hydro-force computation...\n");


#ifdef ACC
      on_end_of_timestep();
#endif
      force_update_hmax();	/* update smoothing lengths in tree */

#ifdef ACC
      on_beginning_of_timestep();

#endif

      hydro_force();		/* adds hydrodynamical accelerations  and computes du/dt  */

#ifdef WINDTUNNEL
      set_injection_accel();
#endif
      VERBOSE(2, "hydro force computation done.\n");

    }

#ifdef AXION_DM
  if(All.TotN_axion > 0)
    {
      VERBOSE(2, "Start axion hydro-force computation...\n");


      ax_force_update_hmax();	/* update smoothing lengths in tree */
      ax_hydro_force();		/* adds hydrodynamical accelerations  and computes du/dt  */

      VERBOSE(2, "hydro axion force computation done.\n");

    }
#endif

  POP_T();
}

#elif GADGET_HYDRO == HYDRO_MFM
void compute_gradients(void)
{
  PUSH_T("gradients");

  double t0, t1;
  if(All.TotN_gas > 0)
    {
      VERBOSE(2, "Start gradient computation...\n");
      
      force_update_hmax();	/* update smoothing lengths in tree */


      gradients();		/* computes gradients and matrixes */

      VERBOSE(2, "gradient computation done.\n");
    }

  POP_T();
}


void compute_fluxes(void)
{
  PUSH_T("fluxes");

  double t0, t1;
  if(All.TotN_gas > 0)
    {
      VERBOSE(2, "Start flux computation...\n");


      mfm_fluxes();		// computes godunov fluxes

      VERBOSE(2, "flux computation done.\n");
    }

  POP_T();
}


void compute_limiters(void)
{
  PUSH_T("slope_limiters");

  double t0, t1;
  if(All.TotN_gas > 0)
    {
      VERBOSE(2, "Start slope limiter computation...\n");


      mfm_limiters();		// computes slope limiter terms
      
      VERBOSE(2, "slope limiter computation done.\n");
    }

  POP_T();
}
#endif
