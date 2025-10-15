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

/*! \file timestep.c
 *  \brief routines for 'kicking' particles in
 *  momentum space and assigning new timesteps
 */

static double dt_displacement = 0;
static double modgrav_timestep = 0;


#ifdef LT_STELLAREVOLUTION
static double time_convert_factor;
#endif // LT_STELLAREVOLUTION

#ifdef GADGET3_IO_LIB
double INLINE_FUNC hubble_function(double a)
{
	double hubble_a;

#ifdef EXTERNALHUBBLE
	hubble_a = hubble_function_external(a);
#else
	hubble_a = All.Omega0 / (a * a * a) + OMEGAK / (a * a)
#ifdef INCLUDE_RADIATION
		/*Note OMEGAR is defined to be 0 if INCLUDE_RADIATION is not on */
		+ OMEGAR / (a * a * a * a)
#ifdef KSPACE_NEUTRINOS_2
		/*Add massive neutrinos, possibly slightly relativistic, to the evolution */
		+ OmegaNu(a) - All.OmegaNu / (a * a * a)
#endif
#endif
#ifdef DARKENERGY
		+ DarkEnergy_a(a);
#else
		+ All.OmegaLambda;
#endif
	hubble_a = All.Hubble * sqrt(hubble_a);
#endif
#ifdef TIMEDEPGRAV
	hubble_a *= dHfak(a);
#endif
	return (hubble_a);
}
#endif

void set_cosmo_factors_for_current_time(void)
{

  if(All.ComovingIntegrationOn)
    {
      All.cf_atime = All.Time;
      All.cf_a2inv = 1 / (All.Time * All.Time);
      All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
      All.cf_afac1 = pow(All.Time, 3 * GAMMA_MINUS1);
      All.cf_afac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      All.cf_afac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      All.cf_hubble_a = hubble_function(All.Time);
    }
  else
    {
      All.cf_atime = 1;
      All.cf_a2inv = 1;
      All.cf_a3inv = 1;
      All.cf_afac1 = 1;
      All.cf_afac2 = 1;
      All.cf_afac3 = 1;
      All.cf_hubble_a = 1;
    }
}


/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
#ifndef GADGET3_IO_LIB

void find_timesteps(void)
{
  double tstart = second();
  CPU_Step[CPU_MISC] += measure_time();

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin || dt_displacement == 0)
    {
      find_dt_displacement_constraint(All.cf_hubble_a * All.cf_atime * All.cf_atime);
    }

#ifdef LT_STELLAREVOLUTION
  time_convert_factor = (1000 * SEC_PER_MEGAYEAR) /	/* chemical time is in Gyr */
    All.UnitTime_in_s *		/* convert in code units   */
    All.HubbleParam;
#endif // LT_STELLAREVOLUTION

#ifdef MAKEGLASS
  do_glass_making_step();
#endif // MAKEGLASS

#if defined(FORCE_EQUAL_TIMESTEPS) || defined(KD_REDUCE_TIMESTEP_BH_AND_STARS)
  integertime ti_min_glob;
  {
    integertime ti_min = TIMEBASE;
    integertime ti_step;
    for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
      {
	if(P[i].Type == 0 && P[i].Mass != 0)
	  {
	    double aphys;
	    ti_step = get_timestep(i, &aphys, 0);

	    if(ti_step < ti_min)
	      {
		ti_min = ti_step;
	      }
	  }
      }

    if(ti_min > (dt_displacement / All.Timebase_interval))
      {
	ti_min = (dt_displacement / All.Timebase_interval);
      }

#ifdef KD_REDUCE_TIMESTEP_BH_AND_STARS
    ti_step = ti_min;		/* we only need the minimum timestep, no need to bring it to power of 2 */
#else // KD_REDUCE_TIMESTEP_BH_AND_STARS
    ti_step = TIMEBASE;
    while(ti_step > ti_min)
      {
	ti_step >>= 1;
      }
#endif //KD_REDUCE_TIMESTEP_BH_AND_STARS

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
    minimum_large_ints(1, &ti_step, &ti_min_glob);
#else // ENLARGE_DYNAMIC_RANGE_IN_TIME
    MPI_Allreduce(&ti_step, &ti_min_glob, 1, MPI_INT, MPI_MIN, MYMPI_COMM_WORLD);
#endif // ENLARGE_DYNAMIC_RANGE_IN_TIME
  }
#endif // if defined(FORCE_EQUAL_TIMESTEPS) || defined(KD_REDUCE_TIMESTEP_BH_AND_STARS)

#ifdef RELAXOBJECT
  determine_relaxfac();
#endif // RELAXOBJECT

  int same = 0, changed = 0;
  /* Now assign new timesteps  */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:same,changed)	// if(NActivePart>500)
  for(int il = 0; il < NActivePart; il++)
    {
      int bin, binold, prev, next;
      integertime ti_step;
      int i = ActiveParticleList[il];
#else /* _OPENMP */
  for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      int bin, binold, prev, next;
      integertime ti_step;
#endif // _OPENMP
#ifdef FORCE_EQUAL_TIMESTEPS
      ti_step = ti_min_glob;
#else // FORCE_EQUAL_TIMESTEPS
      double aphys;
      ti_step = get_timestep(i, &aphys, 0);
#endif // FORCE_EQUAL_TIMESTEPS

#ifdef KD_REDUCE_TIMESTEP_BH_AND_STARS
      if(P[i].Type == 4 || P[i].Type == 5)
	{
	  if(ti_step > ti_min_glob)
	    {
	      ti_step = (int) sqrt((double) ti_step * (double) ti_min_glob);	/* geometric mean of global and local timestep */
	    }
	}
#endif // KD_REDUCE_TIMESTEP_BH_AND_STARS

      /* make it a power 2 subdivision */
      integertime ti_min = TIMEBASE;
      for(ti_min; ti_min > ti_step;)
	{
	  ti_min >>= 1;
	}

      /* old version
         ti_min = TIMEBASE;
         while(ti_min > ti_step)
         ti_min >>= 1;
       */

      ti_step = ti_min;

      bin = get_timestep_bin(ti_step);
      binold = P[i].TimeBin;

      if(bin > binold)		/* timestep wants to increase */
	{
	  while(TimeBinActive[bin] == 0 && bin > binold)
	    {			/* make sure the new step is synchronized */
	      bin--;
	    }

	  ti_step = bin ? (((integertime) 1) << bin) : 0;
	}

      if(All.Ti_Current >= TIMEBASE)	/* we here finish the last timestep. */
	{
	  ti_step = 0;
	  bin = 0;
	}

      if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	{
	  PANIC("we are beyond the end of the timeline");	/* should not happen */

	  for(integertime ti_min = TIMEBASE; ti_min > TIMEBASE - All.Ti_Current;)
	    {
	      ti_min >>= 1;
	      ti_step = ti_min;
	    }

	  /* old version
	     ti_step = TIMEBASE - All.Ti_Current;
	     ti_step = ti_min;
	     while(ti_min > ti_step)
	     ti_min >>= 1;
	   */

	}
      if(bin != binold)
	{
	  changed++;
	  if(P[i].Type == 0)
	    {
#pragma omp critical(_fix_timebins_sph_)
	      {
#ifdef SFR
		TimeBinSfr[binold] -= SphP[i].Sfr;
		TimeBinSfr[bin] += SphP[i].Sfr;
#endif // SFR
	      }
	    }

#ifdef BLACK_HOLES
	  if(P[i].Type == 5)
	    {
#pragma omp critical(_fix_timebins_bhs_)
	      {
		TimeBin_BH_mass[binold] -= BPP(i).BH_Mass;
		TimeBin_BH_dynamicalmass[binold] -= P[i].Mass;
		TimeBin_BH_Mdot[binold] -= BPP(i).BH_Mdot;
		if(BPP(i).BH_Mass > 0)
		  TimeBin_BH_Medd[binold] -= BPP(i).BH_Mdot / BPP(i).BH_Mass;
		TimeBin_BH_mass[bin] += BPP(i).BH_Mass;
		TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
		TimeBin_BH_Mdot[bin] += BPP(i).BH_Mdot;
		if(BPP(i).BH_Mass > 0)
		  TimeBin_BH_Medd[bin] += BPP(i).BH_Mdot / BPP(i).BH_Mass;
	      }
	    }
#endif // BLACK_HOLES

#pragma omp critical(_fix_timebins_)
	  {
	    TimeBinCount[binold]--;
	    prev = PrevInTimeBin[i];
	    next = NextInTimeBin[i];

	    if(FirstInTimeBin[binold] == i)
	      {
		FirstInTimeBin[binold] = next;
	      }
	    if(LastInTimeBin[binold] == i)
	      {
		LastInTimeBin[binold] = prev;
	      }
	    if(prev >= 0)
	      {
		NextInTimeBin[prev] = next;
	      }
	    if(next >= 0)
	      {
		PrevInTimeBin[next] = prev;
	      }

	    if(TimeBinCount[bin] > 0)
	      {
		PrevInTimeBin[i] = LastInTimeBin[bin];
		NextInTimeBin[LastInTimeBin[bin]] = i;
		NextInTimeBin[i] = -1;
		LastInTimeBin[bin] = i;
	      }
	    else
	      {
		FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
		PrevInTimeBin[i] = NextInTimeBin[i] = -1;
	      }
	    TimeBinCount[bin]++;
	    if(P[i].Type == 0)
	      {
		TimeBinCountSph[binold]--;
		TimeBinCountSph[bin]++;
	      }
	  }
	  P[i].TimeBin = bin;
	}
      else
	{
	  same++;
	}

      integertime ti_step_old;
      if(WAKEUP > 0)
	{
	  ti_step_old = P[i].dt_step;
	  P[i].dt_step = ti_step;
	}
      else
	{
	  ti_step_old = binold ? (((integertime) 1) << binold) : 0;
	}
      P[i].Ti_begstep += ti_step_old;

    }


#ifdef CONDUCTION
  if(All.Conduction_Ti_endstep == All.Ti_Current)
    {
      integertime ti_step = TIMEBASE;
      while(ti_step > (All.MaxSizeConductionStep / All.Timebase_interval))
	{
	  ti_step >>= 1;
	}
      while(ti_step > (All.MaxSizeTimestep / All.Timebase_interval))
	{
	  ti_step >>= 1;
	}

      if(ti_step > (All.Conduction_Ti_endstep - All.Conduction_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.Conduction_Ti_endstep) % ti_step) > 0)
	    {
	      ti_step = All.Conduction_Ti_endstep - All.Conduction_Ti_begstep;	/* leave at old step */
	    }
	}

      if(All.Ti_Current == TIMEBASE)
	{			/* we here finish the last timestep. */
	  ti_step = 0;
	}

      All.Conduction_Ti_begstep = All.Conduction_Ti_endstep;
      All.Conduction_Ti_endstep = All.Conduction_Ti_begstep + ti_step;
    }
#endif // CONDUCTION

#ifdef LMB_SPECTRAL_CRs_DIFFUSION
  if(All.CRDiffusion_Ti_endstep == All.Ti_Current)
    {
      integertime ti_step = TIMEBASE;
      while(ti_step > (All.MaxSizeCRDiffusionStep / All.Timebase_interval))
	{
	  ti_step >>= 1;
	}
      while(ti_step > (All.MaxSizeTimestep / All.Timebase_interval))
	{
	  ti_step >>= 1;
	}

      if(ti_step > (All.CRDiffusion_Ti_endstep - All.CRDiffusion_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.CRDiffusion_Ti_endstep) % ti_step) > 0)
	    ti_step = All.CRDiffusion_Ti_endstep - All.CRDiffusion_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      All.CRDiffusion_Ti_begstep = All.CRDiffusion_Ti_endstep;
      All.CRDiffusion_Ti_endstep = All.CRDiffusion_Ti_begstep + ti_step;
    }
#endif // LMB_SPECTRAL_CRs_DIFFUSION

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      integertime ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	{
	  ti_step >>= 1;
	}

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  integertime bin = get_timestep_bin(ti_step);
	  integertime binold = get_timestep_bin(All.PM_Ti_endstep - All.PM_Ti_begstep);

	  while(TimeBinActive[bin] == 0 && bin > binold)
	    {			/* make sure the new step is synchronized */
	      bin--;
	    }

	  ti_step = bin ? (((integertime) 1) << bin) : 0;
	}

      if(All.Ti_Current == TIMEBASE)
	{			/* we here finish the last timestep. */
	  ti_step = 0;
	}

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;
    }
#endif // PMGRID

  if(WAKEUP > 0)
    {
      double t0 = second();
      process_wake_ups();
      double t1 = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
      VERBOSE(5, "EXTRA TIMER: wakeup took %g sec\n", timediff(t0, t1));
#endif // KD_EXTRA_TIMER_OUTPUT
    }

#ifdef KD_EXTRA_TIMER_OUTPUT
  VERBOSE(5, "Find Next Timestep: same=%d changed=%d\n", same, changed);
#endif // KD_EXTRA_TIMER_OUTPUT

  CPU_Step[CPU_TIMELINE] += measure_time();

  double tend = second();
#ifdef KD_EXTRA_TIMER_OUTPUT
  VERBOSE(5, "EXTRA TIMER: finding timesteps took %g sec\n", timediff(tstart, tend));
#endif // KD_EXTRA_TIMER_OUTPUT

}

#endif



/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
integertime get_timestep(int p,	/*!< particle index */
			 double *aphys,	/*!< acceleration (physical units) */
			 int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
					   aphys */ )
{

#if (defined(fSIDM) || defined(rSIDM)) && defined(mSIDM_TIMESTEP)
  double dt_msidm;
  if(P[p].Type == 1)
    {
      dt_msidm = msidm_scatter.calc_dt(P[p]);
#ifndef mSIDM_VDEP0
      P[p].omega_frequent = 0.0;
#ifdef rSIDM
      P[p].omega_rare = 0.0;
#endif // rSIDM
#endif // mSIDM_VDEP0
    }
#endif // if (defined(fSIDM) || defined(rSIDM)) && defined(mSIDM_TIMESTEP)

  double ac;
  if(flag == 0)
    {
      double ax, ay, az;
      ax = All.cf_a2inv * P[p].GravAccel[0];
      ay = All.cf_a2inv * P[p].GravAccel[1];
      az = All.cf_a2inv * P[p].GravAccel[2];
#ifdef PMGRID
      ax += All.cf_a2inv * P[p].GravPM[0];
      ay += All.cf_a2inv * P[p].GravPM[1];
      az += All.cf_a2inv * P[p].GravPM[2];
#endif // PMGRID

#if defined(AB_TURB) || defined(TURB_DRIVING)
      if(P[p].Type == 0)
	{
	  ax += SphP[p].TurbAccel[0];
	  ay += SphP[p].TurbAccel[1];
	  az += SphP[p].TurbAccel[2];
	}
#endif // if defined(AB_TURB) || defined(TURB_DRIVING)

      if(P[p].Type == 0)
	{
	  ax += All.cf_afac2 * SphP[p].HydroAccel[0];
	  ay += All.cf_afac2 * SphP[p].HydroAccel[1];
	  az += All.cf_afac2 * SphP[p].HydroAccel[2];
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    {
      ac = *aphys;
    }

  if(ac == 0)
    {
      ac = 1.0e-30;
    }

  double dt = 0;
  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;

	  dt /= All.cf_hubble_a;	/* convert dloga to physical timestep  */
#ifdef ADAPTGRAVSOFT
	  ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * P[p].AGS_Hsml / 2.8 / (dt * dt);
#else // ADAPTGRAVSOFT
	  ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / (dt * dt);
#endif // ADAPTGRAVSOFT
	  *aphys = ac;
	  return flag;
	}
#ifdef ADAPTGRAVSOFT
      dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * P[p].AGS_Hsml / 2.8 / ac);
#else // ADAPTGRAVSOFT
      dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac);
#endif // ADAPTGRAVSOFT

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
      if(P[p].Type == 0)
	dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * P[p].Hsml / 2.8 / ac);
#else // ADAPTIVE_GRAVSOFT_FORGAS_HSML
      if(P[p].Type == 0)
	dt =
	  sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] *
	       pow(P[p].Mass / All.ReferenceGasMass, 1.0 / 3) / ac);
#endif // ADAPTIVE_GRAVSOFT_FORGAS_HSML
#endif // ADAPTIVE_GRAVSOFT_FORGAS
      break;

    default:
      fprintf(stderr, "\n !!!2@@@!!! \n");
      endrun(888);
      fprintf(stderr, "\n !!!2@@@!!! \n");
      break;
    }

#ifdef KD_BH_TMESTEP_FAC
  if(P[p].Type == 5)
    {
      dt /= KD_BH_TMESTEP_FAC;
    }
#endif // KD_BH_TMESTEP_FAC

  double csnd = 0, dt_divv = 0, dt_courant = 0;
  if(IS_SPH_ALIVE(P[p]))
    {
#if GADGET_HYDRO == HYDRO_SPH
#ifndef LMB_SPECTRAL_CRs
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].Density);
#else // LMB_SPECTRAL_CRs
      csnd = sqrt((GAMMA * (SphP[p].Pressure - (SphP[p].CRpPressure + SphP[p].CRePressure))
		   + All.CR_Gamma * (SphP[p].CRpPressure + SphP[p].CRePressure)) / SphP[p].Density);
#endif // LMB_SPECTRAL_CRs

#elif GADGET_HYDRO == HYDRO_PESPH
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].Density);
#elif GADGET_HYDRO == HYDRO_MFM
      csnd = eos->SoundSpeed(SphP[p].Density, SphP[p].InternalEnergyPred, All.cf_a3inv);
#endif //  GADGET_HYDRO == HYDRO_SPH

      if(All.ComovingIntegrationOn)
	{
	  dt_courant = 2 * All.CourantFac * All.Time * P[p].Hsml / (All.cf_afac3 * SphP[p].MaxSignalVel);
	}
      else
	{
	  dt_courant = 2 * All.CourantFac * P[p].Hsml / SphP[p].MaxSignalVel;
	}
      //#ifdef AR_MAX_DELAYTIME_AFFECTS_SIGNALVEL
      //if(SphP[p].DelayTime<=AR_MAX_DELAYTIME_AFFECTS_SIGNALVEL) //if winds are decoupled why their timesteps should be influenced
      //#endif
      if(dt_courant < dt)
	{
	  dt = dt_courant;
	}

#ifdef DEDNER_TIMESTEP
      double dt_ded;
      if(SphP[p].divB != 0)
	{
	  dt_ded = sqrt(SphP[p].Density / (MU0_1 * SphP[p].divB * SphP[p].divB));
	  if(P[p].ID != 0)
	    {
	      if(P[p].TimeBin > 10 && TimeBinCountSph[P[p].TimeBin] > All.TotN_gas / 100.)
		{
		  dt = dt_ded < dt ? dt_ded : dt;
		}
	    }
	}
#endif // DEDNER_TIMESTEP

      /* make sure that the velocity divergence does not imply a too large change of density or smoothing length in the step */
      if(SphP[p].DivVel != 0)
	{
	  dt_divv = 1.5 / fabs(All.cf_a2inv * SphP[p].DivVel);
#ifdef DO_NOT_PATCH_DIVV_TIMESTEP
	  if(dt_divv < dt)
	    {
	      dt = dt_divv;
	    }
#else // DO_NOT_PATCH_DIVV_TIMESTEP
	  if(dt_divv < dt / 2)
	    {
	      VERBOSE_ALL(2, "WARNING divv timestep: Part-ID=%llu  dt=%g dtc=%g dtv=%g  hsml=%g\n",
			  (unsigned long long) P[p].ID, dt, dt_courant, dt_divv, (float) P[p].Hsml);
	    }
	  if(dt_divv < dt && dt_divv > dt / 4)
	    {
	      dt = dt_divv;
	    }
#endif // DO_NOT_PATCH_DIVV_TIMESTEP
	}


#ifdef LMB_SPECTRAL_CRs_TIMESTEP
      /*  Function in CosmicRays/bp_cr_timestep.h  */
      //MyAtLeastDouble dt_cr = compute_bp_cr_timestep_from_cooling_times(p);

      if(SphP[p].DivVel != 0)
	{
	  // The factor 1.4e6 comes from tests
	  MyAtLeastDouble dt_cr = fabs(SphP[p].DivVel) / (All.cf_a2inv * 1.4e6);
	  if(dt_cr < dt)
	    {
	      dt = dt_cr;
	    }
	}

#endif // LMB_SPECTRAL_CRs_TIMESTEP

#ifdef MYFALSE
      dt_viscous = All.CourantFac * SphP[p].MaxViscStep / All.cf_hubble_a;	/* to convert dloga to physical dt */

      if(dt_viscous < dt)
	{
	  dt = dt_viscous;
	}
#endif // MYFALSE

#ifdef NAVIERSTOKES_TIMESTEP
      double dt_NS;
      if(fabs(SphP[p].ViscEntropyChange))
	{
	  dt_NS = NAVIERSTOKES_TIMESTEP * SphP[p].Entropy / SphP[p].ViscEntropyChange / All.cf_hubble_a;

	  if(dt_NS < dt)
	    {
	      dt = dt_NS;
	    }
	}
#endif // NAVIERSTOKES_TIMESTEP
    }

#ifdef BLACK_HOLES
  double dt_accr;
  if(P[p].Type == 5)
    {
      if(BPP(p).BH_Mdot > 0 && BPP(p).BH_Mass > 0)
	{
	  dt_accr = 0.25 * BPP(p).BH_Mass / BPP(p).BH_Mdot;
	  if(dt_accr < dt)
	    {
	      dt = dt_accr;
	    }
	}

      double dt_ngbs =
	(BPP(p).BH_TimeBinGasNeighbor ? (1 << BPP(p).BH_TimeBinGasNeighbor) : 0) * All.Timebase_interval /
	All.cf_hubble_a;

      if(dt > dt_ngbs && dt_ngbs > 0)
	{
	  dt = 1.01 * dt_ngbs;
	}
    }
#endif // BLACK_HOLES

#ifdef BH_BUBBLES
  if(P[p].Type == 5)
    {
      if(BPP(p).BH_Mdot > 0 && BPP(p).BH_Mass > 0)
	{
#ifdef UNIFIED_FEEDBACK
	  double meddington;
	  meddington = (4 * M_PI * GRAVITY * C_LIGHT * PROTONMASS /
			(All.BlackHoleRadiativeEfficiency * C_LIGHT * C_LIGHT * THOMPSON)) * BPP(p).BH_Mass *
	    All.UnitTime_in_s;
	  if(BPP(p).BH_Mdot < All.RadioThreshold * meddington)
#endif // UNIFIED_FEEDBACK
	    dt_accr = (All.BlackHoleRadioTriggeringFactor - 1) * BPP(p).BH_Mass / BPP(p).BH_Mdot;
	  if(dt_accr < dt)
	    {
	      dt = dt_accr;
	    }
	}
    }
#endif // BH_BUBBLES

#ifdef NONEQUILIBRIUM
  double dt_cool, dt_elec;
  /* another criterion given by the local cooling time */
  if(P[p].Type == 0)
    {
      dt_cool = fabs(SphP[p].t_cool);	/* still in yrs */
      dt_cool *= SEC_PER_YEAR;	/* in seconds */
      dt_cool /= All.UnitTime_in_s;
      dt_cool *= All.HubbleParam;	/* internal units */

      dt_cool = All.Epsilon * dt_cool;

#ifndef UM_CONTINUE
      if(dt_cool > 0 && dt_cool < dt)
	{
	  dt = dt_cool;
	}
#else // UM_CONTINUE
#ifdef WINDS
      if(dt_cool > 0 && dt_cool < dt && SphP[p].DelayTime < 0)
	{
	  dt = dt_cool;
	}
#else // WINDS
      if(dt_cool > 0 && dt_cool < dt)
	{
	  dt = dt_cool;
	}
#endif // WINDS
#endif // UM_CONTINUE

      /* yet another criterion given by the electron number density change */
      dt_elec = fabs(SphP[p].t_elec);	/* still in yrs */
      dt_elec *= SEC_PER_YEAR;	/* in seconds */
      dt_elec /= All.UnitTime_in_s;
      dt_elec *= All.HubbleParam;	/* internal units */
      dt_elec = All.Epsilon * dt_elec;

#ifndef UM_CONTINUE
      if(dt_elec > 0 && dt_elec < dt)
	{
	  dt = dt_elec;
	}
#else // UM_CONTINUE
      if(dt_elec > 0 && dt_elec < dt && SphP[p].DelayTime < 0)
	{
	  dt = dt_elec;
	}
#endif // UM_CONTINUE
    }
#endif // NONEQUILIBRIUM

#if defined(LT_STELLAREVOLUTION) && defined(LT_CONSTRAIN_DYNAMICAL_TIMESTEP)
  double dt_chem;
  if(P[i].Type == 4)
    {
      dt_chem = MetP[P[i].MetID].NextChemTime - All.Time_Age;	/* this is the next chemical time in Gyr */
      dt_chem *= time_convert_factor;	/* convert in code physical units        */

      if(dt_chem > 0 && dt_chem < dt)
	{
	  dt = dt_chem;
	}
    }
#endif // if defined(LT_STELLAREVOLUTION) && defined(LT_CONSTRAIN_DYNAMICAL_TIMESTEP)

#if (defined(fSIDM) || defined(rSIDM)) && defined(mSIDM_TIMESTEP)
  if(P[p].Type == 1)
    {
      dt = std::min(dt, dt_msidm);
    }
#endif // if (defined(fSIDM) || defined(rSIDM)) && defined(mSIDM_TIMESTEP)

  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     All.cf_hubble_a=1.
   */
  dt *= All.cf_hubble_a;

#ifdef ONLY_PM
  dt = All.MaxSizeTimestep;
#endif // ONLY_PM

  if(dt >= All.MaxSizeTimestep)
    {
      dt = All.MaxSizeTimestep;
    }

  if(dt >= dt_displacement)
    {
      dt = dt_displacement;
    }

#ifdef CONDUCTION
  if(P[p].Type == 0)
    {
      if(dt >= All.MaxSizeConductionStep)
	{
	  dt = All.MaxSizeConductionStep;
	}
    }
#endif // CONDUCTION

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      VERBOSE_ALL(2, "warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
	  VERBOSE_ALL
	    (5, "Part-ID=%llu  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
	     (unsigned long long) P[p].ID, dt, dt_courant * All.cf_hubble_a, ac, (double) P[p].Pos[0],
	     (double) P[p].Pos[1], (double) P[p].Pos[2], (float) P[p].Hsml, csnd,
	     sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac) *
	     All.cf_hubble_a, All.SofteningTable[P[p].Type]);
#ifdef NAVIERSTOKES_TIMESTEP
	  VERBOSE_ALL
	    (5, "NAVIERSTOKES:  dt_NS=%g  A=%g  rho=%g  dotAvisc=%g  dtold=%g, meanpath=%g \n",
	     dt_NS * All.cf_hubble_a, SphP[p].Entropy, SphP[p].Density,
	     SphP[p].ViscEntropyChange, (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval,
	     All.IonMeanFreePath *
	     pow((SphP[p].Entropy * pow(SphP[p].Density * All.cf_a3inv, GAMMA_MINUS1) / GAMMA_MINUS1),
		 2.0) / SphP[p].Density);
	  VERBOSE_ALL(5, "    Stressd=(%g|%g|%g) \n", SphP[p].u.s.StressDiag[0], SphP[p].u.s.StressDiag[1],
		      SphP[p].u.s.StressDiag[2]);
	  VERBOSE_ALL(5, "    Stressoffd=(%g|%g|%g) \n", SphP[p].u.s.StressOffDiag[0],
		      SphP[p].u.s.StressOffDiag[1], SphP[p].u.s.StressOffDiag[2]);
#endif // NAVIERSTOKES_TIMESTEP
	}
      else
	{
	  VERBOSE_ALL(5, "Part-ID=%llu  dt=%g ac=%g xyz=(%g|%g|%g)\n", (unsigned long long) P[p].ID, dt, ac,
		      (double) P[p].Pos[0], (double) P[p].Pos[1], (double) P[p].Pos[2]);
	}
      fflush(stdout);
      fprintf(stderr, "\n @ fflush \n");
      endrun(888);
#endif // NOSTOP_WHEN_BELOW_MINTIMESTEP
      dt = All.MinSizeTimestep;
    }

  integertime ti_step = (integertime) (dt / All.Timebase_interval);

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
  if(ti_step == 0)
    {
      VERBOSE(0, "\nError: A timestep of size zero was assigned on the integer timeline!\n"
	      "We better stop.\n"
	      "Task=%d Part-ID=%llu dt=%g dt_elec=%g dt_cool=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g)\n\n",
	      ThisTask, (unsigned long long) P[p].ID, dt, (float) SphP[p].t_elec, (float) SphP[p].t_cool,
	      All.Timebase_interval, (int) ti_step, ac, (double) P[p].Pos[0], (double) P[p].Pos[1],
	      (double) P[p].Pos[2]);
      endrun(818);
    }

#ifdef GM_STARDENSITY
  if(P[p].Type == 4)
    {
      VERBOSE_ALL(5, "    Stellar density, hsml: %g %g\n", P[p].StarDensity, P[p].StarHsml);
    }
#endif // GM_STARDENSITY

#endif // if defined(CHEMISTRY) || defined(UM_CHEMISTRY)

  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      VERBOSE_ALL(0, "\nError: A timestep of size zero was assigned on the integer timeline!\n"
		  "We better stop.\n"
		  "Task=%d Part-ID=%llu Type=%d mass=%g dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
		  ThisTask, (unsigned long long) P[p].ID, P[p].Type, (float) P[p].Mass, dt, dt_courant,
		  dt_divv, dt_displacement, All.Timebase_interval, ti_step, ac, (double) P[p].Pos[0],
		  (double) P[p].Pos[1], (double) P[p].Pos[2], (double) P[p].GravAccel[0],
		  (double) P[p].GravAccel[1], (double) P[p].GravAccel[2]);
#ifdef PMGRID
      VERBOSE_ALL(0, "pm_force=(%g|%g|%g)\n", (double) P[p].GravPM[0], (double) P[p].GravPM[1],
		  (double) P[p].GravPM[2]);
#endif // PMGRID

      if(P[p].Type == 0)
	{
	  VERBOSE_ALL(5, "hydro-frc=(%g|%g|%g) dens=%g hsml=%g  dTdE=%g\n", (double) SphP[p].HydroAccel[0],
		      (double) SphP[p].HydroAccel[1], (double) SphP[p].HydroAccel[2],
		      (double) SphP[p].Density, (float) P[p].Hsml, (float) SphP[p].DtEntropy);
	}

#ifdef LMB_SPECTRAL_CRs
      if(!ThisTask)
	if(P[p].Type == 0)
	  report_cr_properties(p);

#endif // LMB_SPECTRAL_CRs

#ifdef GM_STARDENSITY
      if(P[p].Type == 4)
	{
	  VERBOSE_ALL(5, "    Stellar density, hsml: %g %g\n", P[p].StarDensity, P[p].StarHsml);
	}
#endif // GM_STARDENSITY

#ifdef BLACK_HOLES
      if(P[p].Type == 5)
	{
	  VERBOSE_ALL
	    (5,
	     "BH: mass = %g, md = %g, prog=%d, dens=%g entr=%g v=(%g,%g,%g) macc=%g mBHacc=%e p=(%g,%g,%g)\n",
	     BPP(p).BH_Mass, BPP(p).BH_Mdot, BPP(p).BH_CountProgs, BPP(p).BH_Density, BPP(p).BH_Entropy,
	     BPP(p).BH_SurroundingGasVel[0], BPP(p).BH_SurroundingGasVel[1], BPP(p).BH_SurroundingGasVel[2],
	     BPP(p).BH_accreted_Mass, BPP(p).BH_accreted_BHMass, BPP(p).BH_accreted_momentum[0],
	     BPP(p).BH_accreted_momentum[1], BPP(p).BH_accreted_momentum[2]);
#ifdef KD_FRICTION_DYNAMIC
	  VERBOSE_ALL(5, "BH: v=(%g,%g,%g) dens=%g sigma=%g, bmax=%g\n",
		      BPP(p).BH_SurroundingVel[0], BPP(p).BH_SurroundingVel[1], BPP(p).BH_SurroundingVel[2],
		      BPP(p).BH_SurroundingDensity, BPP(p).BH_sigma, BPP(p).BH_bmax);
#endif // KD_FRICTION_DYNAMIC
	}
#endif // BLACK_HOLES

#ifdef DIVBFORCE3
      if(P[p].Type == 0)
	{
	  VERBOSE_ALL(5, "mhd-frc=(%g|%g|%g) mhd-corr=(%g|%g|%g)\n",
		      (float) SphP[p].magacc[0], (float) SphP[p].magacc[1], (float) SphP[p].magacc[2],
		      (float) SphP[p].magcorr[0], (float) SphP[p].magcorr[1], (float) SphP[p].magcorr[2]);
	}
#endif // DIVBFORCE3
      endrun(818);
    }

  return ti_step;
}

/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      int count[6];
      long long count_sum[6];
      double v[6], v_sum[6], mim[6], min_mass[6];
      for(int type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(int i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(P[i].Mass > 0)
	    {
	      if(mim[P[i].Type] > P[i].Mass)
		mim[P[i].Type] = P[i].Mass;
	    }
	  count[P[i].Type]++;
	}
      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MYMPI_COMM_WORLD);
      sumup_large_ints(6, count, count_sum);

#ifdef SFR
      /* add star and gas particles together to treat them on equal footing, using the original gas particle
         spacing. */
      v_sum[0] += v_sum[4];
      count_sum[0] += count_sum[4];
      v_sum[4] = v_sum[0];
      count_sum[4] = count_sum[0];
#ifdef BLACK_HOLES
      v_sum[0] += v_sum[5];
      count_sum[0] += count_sum[5];
      v_sum[5] = v_sum[0];
      count_sum[5] = count_sum[0];
      min_mass[5] = min_mass[0];
#endif // BLACK_HOLES
#endif // SFR

      for(int type = 0; type < 6; type++)
	{
	  double dmean;
	  if(count_sum[type] > 0)
	    {
	      if(type == 0 || (type == 4 && All.StarformationOn))
		{
		  dmean =
		    pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
			1.0 / 3);
		}
	      else
		{
		  dmean =
		    pow(min_mass[type] /
			((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
			1.0 / 3);
		}

#ifdef BLACK_HOLES
	      if(type == 5)
		{
		  dmean =
		    pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
			1.0 / 3);
		}
#endif // BLACK_HOLES
	      double dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

	      double asmth = 0;
#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		{
		  asmth = All.Asmth[1];
		}
#endif // PLACEHIGHRESREGION
	      if(asmth < dmean)
		{
		  dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
		}
#endif // PMGRID

	      VERBOSE(5, "type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		      type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);

#ifdef NEUTRINOS
	      if(type != 2)	/* don't constrain the step to the neutrinos */
#endif // NEUTRINOS
		if(dt < dt_displacement)
		  {
		    dt_displacement = dt;
		  }
	    }
	}

      VERBOSE(5, "displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);

    }
}

int get_timestep_bin(integertime ti_step)
{
  if(ti_step == 0)
    {
      return 0;
    }

  if(ti_step == 1)
    {
      PANIC("time-step of integer size 1 not allowed\n");
    }

  int bin = -1;
  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}

#ifdef MAKEGLASS
void do_glass_making_step(void)
{
  int dispmax = 0, disp2sum = 0;
  for(int i = 0; i < NumPart; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  P[i].GravAccel[j] *= -1;
#ifdef PMGRID
	  P[i].GravPM[j] *= -1;
	  P[i].GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
#endif // PMGRID
	}

      disp = sqrt(P[i].GravAccel[0] * P[i].GravAccel[0] + P[i].GravAccel[1] * P[i].GravAccel[1]
		  + P[i].GravAccel[2] * P[i].GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	{
	  dispmax = disp;
	}
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MYMPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    {
      fac = dmean / globmax;
    }
  else
    {
      fac = 1.0;
    }

  VERBOSE(5, "\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n", dmean, globmax,
	  sqrt(globdisp2sum / All.TotNumPart));

  for(int i = 0; i < NumPart; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].GravAccel[j] = 0;
	}
    }
}
#endif // MAKEGLASS

#ifdef RELAXOBJECT
void determine_relaxfac(void)
{
  if(All.Time < 0.2 * All.TimeMax)
    {
      All.RelaxFac = 1. / All.RelaxBaseFac;
    }
  else
    {
      if(All.Time > 0.8 * All.TimeMax)
	{
	  All.RelaxFac = 0.;
	}
      else
	{
	  All.RelaxFac =
	    1. / (All.RelaxBaseFac * pow(10., (All.Time - 0.2 * All.TimeMax) / (0.6 * All.TimeMax) * 3.));
	}
    }
}
#endif // RELAXOBJECT

#ifndef GADGET3_IO_LIB
void process_wake_ups(void)
{
  /* find the next kick time */
  int ti_next_kick = TIMEBASE;
  for(int i = 0; i < TIMEBINS; i++)
    {
      if(TimeBinCount[i])
	{
	  int ti_next_for_bin = All.Ti_Current;
	  if(i > 0)
	    {
	      int dt_bin = (((integertime) 1) << i);
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  // The whole function here needs to be addapted to use integertime, so the missing ti_next_kick_global here is to remind for this !!
  minimum_large_ints(1, &ti_next_kick, &ti_next_kick_global);
#else // ENLARGE_DYNAMIC_RANGE_IN_TIME
  int ti_next_kick_global;
  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MYMPI_COMM_WORLD);
#endif /* / ENLARGE_DYNAMIC_RANGE_IN_TIME */

  VERBOSE(5, "predicting next timestep: %g\n",
	  (ti_next_kick_global - All.Ti_Current) * All.Timebase_interval);

  int n = 0;
  //  VERBOSE(5,"ThisTask%d wakup\n",ThisTask);
#ifdef LT_STELLAREVOLUTION
  double sfr_local = 0;
#endif
#pragma omp parallel for	// private()  //reduction(+:n)
  for(int i = 0; i < N_gas; i++)
    {
      if(P[i].Type != 0)
	{
	  continue;
	}

#ifdef WINDS			// wind particles do not wakeup, yet they should stay in the hydroloop to keep their MaxWindVelocity updated for when they'll wake up
      if(SphP[i].DelayTime > 0.)
	{
	  continue;
	}
#endif // WINDS

#ifdef LT_STELLAREVOLUTION
#pragma omp atomic
      sfr_local += SphP[i].Sfr;
#endif // LT_STELLAREVOLUTION

      if(SPHP(i, wakeup) == TIMEBINS)
	{
	  continue;
	}

      int binold = P[i].TimeBin;
      if(TimeBinActive[binold])
	{
	  continue;
	}

      int bin = SPHP(i, wakeup) < binold ? SPHP(i, wakeup) : binold;

      if(bin + 1 < binold)
	{

	  VERBOSE_ALL(2, "Task %d (WAKEUP DEBUG): old=%d new=%d\n", ThisTask, binold, bin);

#pragma omp critical(_wakeup_)
	  {
	    TimeBinCount[binold]--;
	    TimeBinCountSph[binold]--;

#ifdef LT_STELLAREVOLUTION
	    TimeBinSfr[binold] -= SphP[i].Sfr;
	    TimeBinSfr[bin] += SphP[i].Sfr;
	    //GM: this fix the double-count for the Sfr of the woken-up particles
#endif // LT_STELLAREVOLUTION

	    int prev = PrevInTimeBin[i];
	    int next = NextInTimeBin[i];

	    if(FirstInTimeBin[binold] == i)
	      {
		FirstInTimeBin[binold] = next;
	      }
	    if(LastInTimeBin[binold] == i)
	      {
		LastInTimeBin[binold] = prev;
	      }
	    if(prev >= 0)
	      {
		NextInTimeBin[prev] = next;
	      }
	    if(next >= 0)
	      {
		PrevInTimeBin[next] = prev;
	      }
	    if(TimeBinCount[bin] > 0)
	      {
		PrevInTimeBin[i] = LastInTimeBin[bin];
		NextInTimeBin[LastInTimeBin[bin]] = i;
		NextInTimeBin[i] = -1;
		LastInTimeBin[bin] = i;
	      }
	    else
	      {
		FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
		PrevInTimeBin[i] = NextInTimeBin[i] = -1;
	      }
	    TimeBinCount[bin]++;
	    TimeBinCountSph[bin]++;

	    P[i].TimeBin = bin;

	    if(TimeBinActive[bin])
	      {
		NumForceUpdate++;
	      }

	    /* correct quantities predicted for a longer timestep */
	    int ti_step_old = P[i].dt_step;

	    int dt_step = ti_next_kick_global - P[i].Ti_current;
	    // P[i].dt_step = dt_step * 2; 
	    P[i].dt_step = (((integertime) 1) << bin);

	    /*
	       dt_entr = (-ti_step_old / 2 + dt_step / 2) * All.Timebase_interval;
	     */
	    int time0 = P[i].Ti_begstep;
	    int time1_old = P[i].Ti_begstep + ti_step_old;
	    int time1_new = P[i].Ti_begstep + dt_step;

	    /* This part has still to be adapted ...
	       #ifdef PMGRID
	       if(All.ComovingIntegrationOn)
	       dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
	       get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
	       else
	       dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
	       #endif
	     */

	    double dt_entr, dt_gravkick, dt_magkick, dt_hydrokick;

	    if(All.ComovingIntegrationOn)
	      {
		dt_entr = (-(time1_old - time0) / 2 + (time1_new - time0) / 2) * All.Timebase_interval;
		dt_gravkick =
		  -get_gravkick_factor(time0, time1_old) / 2 + get_gravkick_factor(time0, time1_new) / 2;
		dt_hydrokick =
		  -get_hydrokick_factor(time0, time1_old) / 2 + get_hydrokick_factor(time0, time1_new) / 2;
#ifdef MAGNETIC
		dt_magkick = -get_magkick_factor(time0, time1_old) / 2
		  + get_magkick_factor(time0, time1_new) / 2;
#endif // MAGNETIC
	      }
	    else
	      {
		dt_entr = dt_gravkick = dt_magkick = dt_hydrokick = (-(time1_old - time0) / 2
								     + (time1_new -
									time0) / 2) * All.Timebase_interval;
	      }

	    /* This may now work in comoving runs */
	    /* WARNING: this velocity correction is inconsistent,
	     * as the position of the particle was calculated with a "wrong" velocity before  */
#ifdef AA_BND_PARTICLES
	    if(P[i].ID < AA_BND_PARTICLES)
	      {
#endif // AA_BND_PARTICLES
		for(int k = 0; k < 3; k++)
		  {
		    P[i].Vel[k] += P[i].GravAccel[k] * dt_gravkick;
		  }

		for(int k = 0; k < 3; k++)
		  {
		    P[i].Vel[k] += SphP[i].HydroAccel[k] * dt_hydrokick;
		  }
#ifdef AA_BND_PARTICLES
	      }
#endif // AA_BND_PARTICLES

#ifdef MAGNETIC
	    for(int k = 0; k < 3; k++)
	      {
		SphP[i].B[k] += SphP[i].DtB[k] * dt_magkick;
	      }
#endif // MAGNETIC

	    SphP[i].Entropy += SphP[i].DtEntropy * dt_entr;
	    // drift_particle(i, All.Ti_Current);
	    P[i].Ti_begstep = P[i].Ti_current;

#if GADGET_HYDRO == HYDRO_SPH
	    SphP[i].Pressure = get_pressure(i);
#elif GADGET_HYDRO == HYDRO_PESPH
	    SphP[i].Pressure = get_pressure(i);
#elif GADGET_HYDRO == HYDRO_MFM
	    SphP[i].Pressure = get_pressure(i);
#endif // if GADGET_HYDRO == HYDRO_SPH
	  }
#pragma omp atomic
	  n++;
	}
    }

  long long ntot = 0;
  sumup_large_ints(1, &n, &ntot);

#ifdef LT_STELLAREVOLUTION
  double sfr_total;
  MPI_Allreduce(&sfr_local, &sfr_total, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);

  VERBOSE(5, "%d%09d particles woken up, total sfr = %g.\n", (int) (ntot / 1000000000),
	  (int) (ntot % 1000000000), sfr_total);

#else // LT_STELLAREVOLUTION
  VERBOSE(5, "%d%09d particles woken up.\n", (int) (ntot / 1000000000), (int) (ntot % 1000000000));
#endif // LT_STELLAREVOLUTION
}
#endif