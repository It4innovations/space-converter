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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3, atime;

#ifdef MAGNETIC
  double a2_fac, gauss2gadget = 1.0;
#endif


#ifdef LMB_SPECTRAL_CRs
  double Pinit_integral;
  double Pth0;
#endif

#ifdef CR_INITPRESSURE
  double cr_pressure, q_phys, C_phys[NUMCRPOP];
#endif
#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
  int ifunc;
  double min_t_cool, max_t_cool;
  double min_t_elec, max_t_elec;
  double a_start, a_end;
#endif
#ifdef BLACK_HOLES
  int count_holes = 0;
#endif

#if defined(GM_MUPPI) && defined(MV_GM_AGNMUPPI_OUTPUT)
  All.TotE_AGNM_tot = All.TotE_AGNM_cool = All.TotE_AGNM_sfr_c = All.TotE_AGNM_sfr_h = 0.0;
#endif


#ifdef START_WITH_EXTRA_NGBDEV
  double MaxNumNgbDeviationMerk;
#endif

#ifdef DISTORTIONTENSORPS
  int i1, i2;
#endif

#ifdef LT_STELLAREVOLUTION
  double this_age, NextChemTime;
  int IMFi, SFi;
  unsigned int I;
  int chem_step, ti_min;
#endif

  double t0_ics, t1_ics;

#ifdef LT_ADD_GAL_TO_SUB
#ifdef OBSERVER_FRAME
#ifdef INTERP_OBSERVER_FRAME
  flag_readCB07obs_redshifts = 0;
#endif
#endif
#endif

  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();

  if(RestartFlag == -1)
    {
#ifdef BLACK_HOLES
      if(ThisTask == 0)
	{
	  printf("BH: Evaluating Settings of the parameters\n");
	  report_internal_model_for_blackhole();
	}
#endif

#ifdef LMB_SPECTRAL_CRs
      if(ThisTask == 0)
	{
	  report_internal_model_for_spectral_crs();
	}
#endif

#ifdef LMB_CR_OUTPUT_SYNCHROTRON
      if(ThisTask == 0)
	{
	  report_synchrotron_kernel();
	}
#endif

      MPI_Barrier(MYMPI_COMM_WORLD);

      PANIC("Ending OpenGadget after explaining all physics modules");
    }


#if defined(SFR) || defined(BLACK_HOLES)
  for(All.bits = 0; GENERATIONS > (1 << All.bits); All.bits++);
#endif

  if(RestartFlag == 3 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");
      endrun(0);
    }

  if(RestartFlag == 4 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if snapshot should be converted\n");
      endrun(0);
    }

  if(RestartFlag == 5 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf
	  ("Need to give the snapshot number if power spectrum and two-point correlation function should be calculated\n");
      endrun(0);
    }

  if(RestartFlag == 6 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf
	  ("Need to give the snapshot number if velocity power spectrum for the gas cells should be calculated\n");
      endrun(0);
    }

  t0_ics = second();
  switch (All.ICFormat)
    {
    case 1:
    case 2:
    case 3:
    case 4:
      if(RestartFlag >= 2 && RestartSnapNum >= 0)
	{
	  char fname[1000];

	  if(All.NumFilesPerSnapshot > 1)
	    sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase,
		    RestartSnapNum);
	  else
	    sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
	  read_ic(fname);

	}
      else
	{
	  if(All.ICFormat == 4)
	    read_ic_cluster(All.InitCondFile);
	  else
	    read_ic(All.InitCondFile);
	}
      break;

    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }
  t1_ics = second();

  if(ThisTask == 0)
    printf("reading ICs took %g sec\n", timediff(t0_ics, t1_ics));

  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();



#ifdef WINDTUNNEL
  All.WindCurrentX = -0.5 * All.WindDmean;
#endif


#ifdef MOL_CLOUDS
  if(ThisTask == 0)
    {
      printf("Init MOL_CLOUDS...\n");
      fflush(stdout);
    }
  init_mol_clouds();
  if(ThisTask == 0)
    {
      printf("done.\n");
      fflush(stdout);
    }
#endif


#ifdef SCFPOTENTIAL
  if(ThisTask == 0)
    {
      printf("Init SCF...\n");
      fflush(stdout);
    }
  SCF_init();
  if(ThisTask == 0)
    {
      printf("Initial random seed = %ld\n", scf_seed);
      printf("done.\n");
      fflush(stdout);
    }
#endif


#ifdef COOLING
#ifndef LT_METAL_COOLING_WAL
  IonizeParams();
#endif
#endif

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
  InitChem();
#endif

#ifdef LMB_SPECTRAL_CRs
  init_spectral_crs_bounds();
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = All.Time * All.Time * All.Time;
      atime = All.Time;
#ifdef MAGNETIC
#ifndef MU0_UNITY
      gauss2gadget *=
	sqrt(All.UnitTime_in_s * All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g /
	     (All.HubbleParam * All.HubbleParam));
#endif
      a2_fac = (All.Time * All.Time);
#endif
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = 1;
      atime = 1;
#ifdef MAGNETIC
#ifndef MU0_UNITY
      gauss2gadget *= sqrt(All.UnitTime_in_s * All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g);
#endif
      a2_fac = 1;
#endif
    }


  set_softenings();


#ifdef SIDM
  sidm_Init_Particles();
#endif

#if defined(OUTPUT_LONGRANGE_POTENTIAL)
  All.PotentialFileCount = 0;
#endif

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    {
      if(RestartSnapNum < 0)
	{
	  char *underscore = strrchr(All.InitCondFile, '_');
	  if(!underscore)
	    {
	      char buf[1000];
	      sprintf(buf, "Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n",
		      All.InitCondFile);
	      PANIC(buf);
	    }
	  else
	    {
	      All.SnapshotFileCount = atoi(underscore + 1) + 1;
#if defined(OUTPUT_LONGRANGE_POTENTIAL)
	      All.PotentialFileCount = atoi(underscore + 1) + 1;
#endif
	    }
	}
      else
	{
	  All.SnapshotFileCount = RestartSnapNum + 1;
	  All.SnapshotFileCount--;
#if defined(OUTPUT_LONGRANGE_POTENTIAL)
	  All.PotentialFileCount = RestartSnapNum + 1;
#endif
	}
    }

  All.TotNumOfForces = 0;
#if defined(MAGNETIC) && defined(BSMOOTH)
#ifdef SETMAINTIMESTEPCOUNT
  All.MainTimestepCounts = All.MainTimestepCountIni;
#else
  All.MainTimestepCounts = 0;
#endif
#endif

  All.TopNodeAllocFactor = 0.008;
  All.TreeAllocFactor = 0.35;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#if defined(BLACK_HOLES) || defined(VARIABLE_WINDS)
  All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif

  for(i = 0; i < GRAVCOSTLEVELS; i++)
    All.LevelToTimeBin[i] = 0;

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < GRAVCOSTLEVELS; j++)
      P[i].GravCost[j] = 0;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }


  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].GravAccel[j] = 0;
	}

      /* DISTORTION PARTICLE SETUP */
#ifdef DISTORTIONTENSORPS
      /*init tidal tensor for first output (not used for calculation) */
      for(i1 = 0; i1 < 3; i1++)
	for(i2 = 0; i2 < 3; i2++)
	  P[i].tidal_tensorps[i1][i2] = 0.0;

      /* find caustics by sign analysis of configuration space distortion */
      P[i].last_determinant = 1.0;

#ifdef OUTPUT_LAST_CAUSTIC
      /* all entries zero -> no caustic yet */
      P[i].lc_Time = 0.0;
      P[i].lc_Pos[0] = 0.0;
      P[i].lc_Pos[1] = 0.0;
      P[i].lc_Pos[2] = 0.0;
      P[i].lc_Vel[0] = 0.0;
      P[i].lc_Vel[1] = 0.0;
      P[i].lc_Vel[2] = 0.0;
      P[i].lc_rho_normed_cutoff = 0.0;

      P[i].lc_Dir_x[0] = 0.0;
      P[i].lc_Dir_x[1] = 0.0;
      P[i].lc_Dir_x[2] = 0.0;
      P[i].lc_Dir_y[0] = 0.0;
      P[i].lc_Dir_y[1] = 0.0;
      P[i].lc_Dir_y[2] = 0.0;
      P[i].lc_Dir_z[0] = 0.0;
      P[i].lc_Dir_z[1] = 0.0;
      P[i].lc_Dir_z[2] = 0.0;

      P[i].lc_smear_x = 0.0;
      P[i].lc_smear_y = 0.0;
      P[i].lc_smear_z = 0.0;
#endif


#ifdef PMGRID
      /* long range tidal field init */
      P[i].tidal_tensorpsPM[0][0] = 0;
      P[i].tidal_tensorpsPM[0][1] = 0;
      P[i].tidal_tensorpsPM[0][2] = 0;
      P[i].tidal_tensorpsPM[1][0] = 0;
      P[i].tidal_tensorpsPM[1][1] = 0;
      P[i].tidal_tensorpsPM[1][2] = 0;
      P[i].tidal_tensorpsPM[2][0] = 0;
      P[i].tidal_tensorpsPM[2][1] = 0;
      P[i].tidal_tensorpsPM[2][2] = 0;
#endif

      for(i1 = 0; i1 < 6; i1++)
	for(i2 = 0; i2 < 6; i2++)
	  {
	    if((i1 == i2))
	      P[i].distortion_tensorps[i1][i2] = 1.0;
	    else
	      P[i].distortion_tensorps[i1][i2] = 0.0;
	  }

      /* for cosmological simulations we do init here, not read from ICs */
      if(All.ComovingIntegrationOn)
	{
#ifndef GDE_READIC
	  /* no caustic passages in the beginning */
	  P[i].caustic_counter = 0.0;
#ifndef GDE_LEAN
	  /* Lagrange time of particle */
	  P[i].a0 = All.TimeBegin;
	  /* approximation: perfect Hubble Flow -> peculiar sheet orientation is exactly zero */
	  for(i1 = 0; i1 < 3; i1++)
	    for(i2 = 0; i2 < 3; i2++)
	      GDE_VMATRIX(i, i1, i2) = 0.0;
	  /* approximation: initial sream density equals background density */
	  P[i].init_density = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#else
	  All.GDEInitStreamDensity = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#endif
#endif
	}

#ifndef GDE_LEAN
      /* annihilation stuff */
      P[i].s_1_last = 1.0;
      P[i].s_2_last = 1.0;
      P[i].s_3_last = 1.0;
      P[i].second_deriv_last = 0.0;
      P[i].rho_normed_cutoff_last = 1.0;

      P[i].s_1_current = 1.0;
      P[i].s_2_current = 1.0;
      P[i].s_3_current = 1.0;
      P[i].second_deriv_current = 0.0;
      P[i].rho_normed_cutoff_current = 1.0;

      P[i].annihilation = 0.0;
      P[i].analytic_caustics = 0.0;
      P[i].analytic_annihilation = 0.0;
#endif

      if(All.ComovingIntegrationOn)
	P[i].stream_density = GDE_INITDENSITY(i) / (All.TimeBegin * All.TimeBegin * All.TimeBegin);
      else
	P[i].stream_density = GDE_INITDENSITY(i);

#endif /* DISTORTIONTENSORPS */

#ifdef KEEP_DM_HSML_AS_GUESS
      if(RestartFlag != 1)
	P[i].DM_Hsml = -1;
#endif

#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_begstep = 0;
      P[i].Ti_current = 0;
      P[i].TimeBin = 0;

      if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
	P[i].OldAcc = 0;	/* Do not zero in 2lpt case as masses are stored here */

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
      P[i].Potential = 0;
#endif
#ifdef STELLARAGE
      if(RestartFlag == 0)
	{
	  if(P[i].Type == 4)
	    MPP(i).StellarAge = 0;

#if defined(BLACK_HOLES)
	  if(P[i].Type == 5)
	    BPP(i).StellarAge = 0;
#endif
	}
#endif

#if defined(KD_MAIN_HALO) && defined(SUBFIND)
      if(RestartFlag == 0)
	P[i].SubLevel = 0;
#endif

#ifdef METALS
      if(RestartFlag == 0)
	P[i].Metallicity = 0;

#endif

#ifdef LT_STELLAREVOLUTION
      if(RestartFlag != 5)
	{
	  All.Time_Age = get_age(All.Time);

	  if(P[i].Type == 4)
	    {
	      this_age = get_age(MPP(i).StellarAge) - All.Time_Age;
	      SFi = get_SF_index(i, &SFi, &IMFi);

	      NextChemTime = get_NextChemTime(this_age, SFi, 0x0);
	      j = get_chemstep_bin(All.Time, All.Time_Age - NextChemTime, &chem_step, i);

	      I = P[i].pt.MetID;

	      if(All.Ti_Current >= TIMEBASE)
		chem_step = j = 0;

	      if((TIMEBASE - All.Ti_Current) < chem_step)
		{
		  chem_step = TIMEBASE - All.Ti_Current;
		  ti_min = TIMEBASE;
		  while(ti_min > chem_step)
		    ti_min >>= 1;
		  chem_step = ti_min;
		  j = get_timestep_bin(chem_step);
		}

	      if(j != MetP[I].ChemTimeBin)
		{
		  TimeBinCountStars[MetP[I].ChemTimeBin]--;
		  TimeBinCountStars[j]++;
		  MetP[I].ChemTimeBin = j;
		}
	      MetP[I].LastChemTime = this_age;

	      if(this_age < SFs[SFi].ShortLiv_TimeTh)
		MetP[I].LastChemTime = SFs[SFi].ShortLiv_TimeTh;
	    }
	}
#endif


#ifdef BLACK_HOLES
      if(P[i].Type == 5)
	{
	  count_holes++;

	  if(RestartFlag == 0)
	    BPP(i).BH_Mass = All.SeedBlackHoleMass;
	}
#endif
    }

#ifdef BLACK_HOLES
  MPI_Allreduce(&count_holes, &All.TotBHs, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
#endif

  for(i = 0; i < TIMEBINS; i++)
    TimeBinActive[i] = 1;

  reconstruct_timebins();

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef CONDUCTION
  All.Conduction_Ti_endstep = All.Conduction_Ti_begstep = 0;
#endif

#ifdef LMB_SPECTRAL_CRs_DIFFUSION
  All.CRDiffusion_Ti_endstep = All.CRDiffusion_Ti_begstep = 0;
#endif

#ifdef  GM_MUPPI_DEBUG
  All.Egy_out = All.Egy_in_thermal = All.Egy_in_kinetic = 0.0;
  All.Egy_outCum = All.Egy_in_thermalCum = All.Egy_in_kineticCum = 0.0;
#endif


  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      SphP[i].EntropyPred = SphP[i].Entropy;

      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = 0;
	}

      SphP[i].DtEntropy = 0;


#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      SphP[i].Gamma = GAMMA;	/* set universal value */
      SphP[i].t_cool = 0;
      SphP[i].t_elec = 0;
#endif



      if(RestartFlag == 0)
	{
#ifndef READ_HSML
	  P[i].Hsml = 0;
#endif
	  SphP[i].Density = -1;
#ifdef COOLING
#ifndef UM_CHEMISTRY
	  SphP[i].elec = 1.0;
#endif
#endif
	  SphP[i].DivVel = 0;
	}
#ifdef WINDS

      if(RestartFlag == 0)
	SphP[i].DelayTime = 0;
#ifdef VARIABLE_WINDS
      SphP[i].HostHaloMass = 0;
#endif
#endif
#ifdef SFR
      SphP[i].Sfr = 0;
#endif
#ifdef STATICBRANDT_INIT
      P[i].Vel[0] = 0.0;
      P[i].Vel[1] = 0.0;
      P[i].Vel[2] = 0.0;

      P[i].Vel[0] = -OmegaR(i, BRANDT_MODE) * (P[i].Pos[1] - LONG_Y / 2.);
      P[i].Vel[1] = OmegaR(i, BRANDT_MODE) * (P[i].Pos[0] - LONG_X / 2.);

      if(P[i].Type == 0)
	for(j = 0; j < 3; j++)
	  SphP[i].VelPred[j] = P[i].Vel[j];
#endif
#if defined (UM_CHEMISTRY) && defined (UM_CHEMISTRY_INISET)
      SphP[i].elec = All.Startelec;

      SphP[i].HI = All.StartHI;
      SphP[i].HII = All.StartHII;
      SphP[i].HM = All.StartHM;

      SphP[i].HeI = All.StartHeI;
      SphP[i].HeII = All.StartHeII;
      SphP[i].HeIII = All.StartHeIII;

      SphP[i].H2I = All.StartH2I;
      SphP[i].H2II = All.StartH2II;

      SphP[i].HD = All.StartHD;
      SphP[i].DI = All.StartDI;
      SphP[i].DII = All.StartDII;

      SphP[i].HeHII = All.StartHeHII;
#endif
#ifdef MAGNETIC
#if defined BINISET && !defined(MAGNETICZERO)
      if(RestartFlag == 0)
	{			/* Set only when starting from ICs */
	  SphP[i].BPred[0] = All.BiniX;
	  SphP[i].BPred[1] = All.BiniY;
	  SphP[i].BPred[2] = All.BiniZ;
	}
#endif /*BINISET*/
	for(j = 0; j < 3; j++)
	{
	  SphP[i].BPred[j] *= a2_fac * gauss2gadget;
	  SphP[i].B[j] = SphP[i].BPred[j];
	}
#ifdef TIME_DEP_MAGN_DISP
#ifdef HIGH_MAGN_DISP_START
      SphP[i].Balpha = All.ArtMagDispConst;
#else
      SphP[i].Balpha = All.ArtMagDispMin;
#endif
      SphP[i].DtBalpha = 0.0;
#endif
#ifdef AB_ART_DISP
      int k;
      for(j = 0; j < 3; j++)
	for(k = 0; k < 3; k++)
	  SphP[i].ArtMagDispMatrix[j][k] = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
      SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#endif
#ifdef TIME_DEP_ART_COND
      SphP[i].Calpha = TIME_DEP_ART_COND * All.ArtCondMin;
#endif
#ifdef TIME_DEP_ART_VISC
#ifdef HIGH_ART_VISC_START
      if(HIGH_ART_VISC_START > 0)
	SphP[i].alpha = HIGH_ART_VISC_START;
      else
	SphP[i].alpha = All.ArtBulkViscConst;
#else
      SphP[i].alpha = All.AlphaMin;
#endif
      SphP[i].Dtalpha = 0.0;
#ifdef AB_ART_VISC
      SphP[i].visc_R = 0.0;
      SphP[i].visc_OldDvel = 0.0;
#endif
#endif

      SphP[i].Injected_BH_Energy = 0;

#ifdef GM_MUPPI
      if(RestartFlag == 0)
	{
	  SphP[i].M_sf = 0.0;
	  SphP[i].M_h = 0.0;
	  SphP[i].E_h = 0.0;
	  SphP[i].Egy_tot_0 = 0.0;
	  SphP[i].out_part = 0;
	  SphP[i].MultiPhase = 0;
	  SphP[i].tdyn = 0.0;
	  SphP[i].clock = 0.0;
	  SphP[i].E_out = 0.0;
	  SphP[i].E_rec = 0.0;
	  SphP[i].Eout_norm = 0.0;
	  SphP[i].GradDens[0] = 0.0;
	  SphP[i].GradDens[1] = 0.0;
	  SphP[i].GradDens[2] = 0.0;

	  SphP[i].E_kin = 0.0;
	  SphP[i].xkin = 0.0;
	  SphP[i].ykin = 0.0;
	  SphP[i].zkin = 0.0;

	  SphP[i].NoArtCond = 0;

#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
	  SphP[i].CouplingBHEnHot = 0.;
	  SphP[i].CouplingBHEnCold = 0.;
#endif


	}
#endif
    }


#ifdef TWODIMS
  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[2] = 0;
      P[i].Vel[2] = 0;

      P[i].GravAccel[2] = 0;

      if(P[i].Type == 0)
	{
	  SphP[i].VelPred[2] = 0;
	  SphP[i].HydroAccel[2] = 0;
	}
    }
#endif

#ifdef ONEDIM
  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[1] = P[i].Pos[2] = 0;
      P[i].Vel[1] = P[i].Vel[2] = 0;

      P[i].GravAccel[1] = P[i].GravAccel[2] = 0;

      if(P[i].Type == 0)
	{
	  SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0;
	  SphP[i].HydroAccel[1] = SphP[i].HydroAccel[2] = 0;
	}
    }
#endif

#ifdef ASSIGN_NEW_IDS
  assign_unique_ids();
#endif

#ifndef NOTEST_FOR_IDUNIQUENESS
  test_id_uniqueness();
#endif

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  TreeReconstructFlag = 1;


#ifdef SHIFT_BY_HALF_BOX
  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      P[i].Pos[j] += 0.5 * All.BoxSize;
#endif


  domain_Decomposition(0, 0);	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;

#ifdef START_WITH_EXTRA_NGBDEV
  MaxNumNgbDeviationMerk = All.MaxNumNgbDeviation;
  All.MaxNumNgbDeviation = All.MaxNumNgbDeviationStart;
#endif

#ifdef ACC_RUN
  on_beginning_of_timestep();
#endif


#if GADGET_HYDRO == HYDRO_MFM
  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      if(header.flag_entropy_instead_u == 0)
	{
	  // TODO : Also add restart option later
	  SphP[i].dQ = MyFloat(0.0);
	  SphP[i].dQdt = MyFloat(0.0);
	  SphP[i].InternalEnergy = SphP[i].Entropy;
	  SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
	  SphP[i].DtInternalEnergy = (MyFloat) 0.0;
	}
    }
#endif



  if(RestartFlag != 3 && RestartFlag != 5)
    setup_smoothinglengths();


#ifdef ACC_RUN
  on_end_of_timestep();
#endif

#ifdef ADAPTGRAVSOFT
  if(RestartFlag != 3 && RestartFlag != 5)
    ags_setup_smoothinglengths();
#endif

#ifdef START_WITH_EXTRA_NGBDEV
  All.MaxNumNgbDeviation = MaxNumNgbDeviationMerk;
#endif


  /* at this point, the entropy variable actually contains the
   * internal energy, read in from the initial conditions file.
   * Once the density has been computed, we can convert to entropy.
   */
  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      if(header.flag_entropy_instead_u == 0)
	{

#if !defined(ISOTHERM_EQS)

	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy !\n");

	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
	  SphP[i].EntropyPred = SphP[i].Entropy;
#endif
	}


      SphP[i].DtEntropy = 0;

      SphP[i].DivVel = 0;

#ifdef LMB_SPECTRAL_CRs
#ifndef LMB_SPECTRAL_CRs_FROM_ICS
      Pth0 = SphP[i].Entropy * pow(SphP[i].Density / a3, GAMMA);
      init_spectral_crs_pressure(i, Pth0);
      init_spectral_crs_arrays(i);
#else // LMB_SPECTRAL_CRs_FROM_ICS
      init_spectral_crs_from_file(i);
#endif // LMB_SPECTRAL_CRs_FROM_ICS
#endif

    }

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)

  if(ThisTask == 0)
    {
      printf("Initial abundances (for P[1].ID=%llu):\n", (unsigned long long) P[1].ID);

      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);
      printf("HM=%g, H2I=%g, H2II=%g, elec=%g\n", SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec);

#if defined (UM_CHEMISTRY) && defined(UM_HD_COOLING)
      printf("HD=%g,  DI=%g, DII=%g\n", SphP[1].HD, SphP[1].DI, SphP[1].DII);
      printf("HeHII=%g\n", SphP[1].HeHII);
#endif

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g,\ndensity=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].Density, SphP[N_gas - 1].Entropy);



#if defined(UM_METAL_COOLING) && defined (LT_METAL_COOLING)
      printf("::  mass of C=%g, O=%g, Si=%g, Fe=%g, N=%g\n",
	     SphP[N_gas - 1].Metals[Carbon], SphP[N_gas - 1].Metals[Oxygen], SphP[N_gas - 1].Metals[Silicon],
	     SphP[N_gas - 1].Metals[Iron], SphP[N_gas - 1].Metals[Nitrogen]);

#endif
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;


  for(i = 0; i < N_gas; i++)	/* Init Chemistry: */
    {
#ifdef CHEMISTRY
      a_start = All.Time;
      a_end = All.Time + 0.001;	/*  0.001 as an arbitrary value */
      ifunc = compute_abundances(0, i, a_start, a_end);
#endif


#ifdef UM_CHEMISTRY		/* Init um_Chemistry/Z: */
#if defined (UM_METAL_COOLING) && defined (LT_METAL_COOLING)
      um_ZsPoint = SphP[i].Metals;
      um_mass = P[i].Mass;
#endif
      double um_u;

      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */
      ifunc = compute_abundances(0, i, a_start, a_end, &um_u);

#endif

      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);
    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif


  if(RestartFlag == 3)
    {
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
      if(ThisTask == 0)
	printf("SUBFIND_RESHUFFLE_AND_POTENTIAL: Calculating potential energy before reshuffling...\n");
#ifdef PMGRID
      long_range_init_regionsize();
#endif
      compute_potential();
      if(ThisTask == 0)
	printf("potential energy done.\n");

#endif

#ifdef ADAPTGRAVSOFT
#ifdef PMGRID
      /* If PMGRID was on during the run, then at least a Maximum softening must have also existed:
       * we then need to take that into account during the computation of the potential energy of the particles */
      long_range_init_regionsize();
      All.AGS_MaxSoft[0] = AGS_MAXSOFT * All.Asmth[0];
      All.AGS_MinSoft[0] = AGS_MINSOFT * All.Asmth[0];
#ifdef PLACEHIGHRESREGION
      All.AGS_MaxSoft[1] = AGS_MAXSOFT * All.Asmth[1];
      All.AGS_MinSoft[1] = AGS_MINSOFT * All.Asmth[1];

#endif
#endif
      if(ThisTask == 0)
	printf("*ADAPTGRAVSOFT* Computation of softening lengths... \n");
      ags_setup_smoothinglengths();
      if(ThisTask == 0)
	printf("*ADAPTGRAVSOFT* Computation of softening lengths done. \n");
#endif

#ifdef FOF
      fof_fof(RestartSnapNum);
#endif

#ifdef LMB_SPECTRAL_CRs_POST_PROCESS

      double t0 = second();
      // CR injection
      inject_crs_postprocess();

      double t1 = second();

#ifdef LMB_CR_TIMER
      if(ThisTask == 0)
	{
	  printf("EXTRA TIMER: CR injection took %g sec\n", timediff(t0, t1));
	}
#endif

      t0 = second();
      // spectral evolution
      evolve_spectral_crs();
      t1 = second();

#ifdef LMB_CR_TIMER
      if(ThisTask == 0)
	{
	  printf("EXTRA TIMER: CR evolution took %g sec\n", timediff(t0, t1));
	}
#endif

      // write the snapshot
      savepositions(RestartSnapNum);

#endif

      endrun(0);
    }

  if(RestartFlag == 5)
    {
      /* calculating powerspec and twopoint function */
#ifdef PMGRID
      long_range_init_regionsize();
#ifdef PERIODIC
      int n, n_type[6];
      long long ntot_type_all[6];
      /* determine global and local particle numbers */
      for(n = 0; n < 6; n++)
	n_type[n] = 0;
      for(n = 0; n < NumPart; n++)
	n_type[P[n].Type]++;
      sumup_large_ints(6, n_type, ntot_type_all);

      calculate_power_spectra(RestartSnapNum, ntot_type_all);
#endif
#endif
      force_treebuild(NumPart, NULL);
      endrun(0);
    }


  if(RestartFlag == 6)
    {
      endrun(0);
    }


#ifdef CHEMISTRY
  if(ThisTask == 0)
    {
      printf("Initial abundances: \n");
      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);

      printf("HM=%g, H2I=%g, H2II=%g, elec=%g, %d\n",
	     SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec, P[1].ID);

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g, density=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].Density, SphP[N_gas - 1].Entropy);
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;

  for(i = 0; i < N_gas; i++)
    {
      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */

      ifunc = compute_abundances(0, i, a_start, a_end);


      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);

    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif

  if(RestartFlag == 4)
    {
      All.Time = All.TimeBegin = header.time;
      sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
      if(ThisTask == 0)
	printf("Start writing file %s\n", All.SnapshotFileBase);
      printf("RestartSnapNum %d\n", RestartSnapNum);

      All.TopNodeAllocFactor = 0.008;

      savepositions(RestartSnapNum);
      endrun(0);
    }
}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));
#ifdef TIMEDEPGRAV
  omega *= All.Gini / All.G;
#endif

#if defined KSPACE_NEUTRINOS || defined KSPACE_NEUTRINOS_2
  omega += All.OmegaNu;
#endif

  if(fabs(omega - All.Omega0) > 1.0e-2)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
#ifndef LMB_SPECTRAL_CRs_POST_PROCESS
	  printf("\nI better stop.\n");
#else
	  printf
	    ("\nBut since we're only running the CR component in post-processing we don't have to care!\nYay!\n");
#endif
	  fflush(stdout);
	}
#ifndef LMB_SPECTRAL_CRs_POST_PROCESS
      endrun(1);
#endif
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

#ifndef READ_HSML
#ifndef TWODIMS
#ifndef ONEDIM
	  P[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  P[i].Hsml = All.DesNumNgb * (P[i].Mass / Nodes[no].u.d.mass) * Nodes[no].len;
#endif
#else
	  P[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	  if(All.SofteningTable[0] != 0 && P[i].Hsml > 200.0 * All.SofteningTable[0])
	    P[i].Hsml = All.SofteningTable[0];
#endif
	}

    }

#if (defined(fSIDM) || defined(rSIDM)) && !defined(ADAPTGRAVSOFT)
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; ++i)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

	  if(P[i].Type == 1)
	    P[i].Hsml =
	      pow(3.0 / (4 * M_PI) * All.DesNumNgb_pt1 * P[i].Mass / Nodes[no].u.d.mass,
		  1.0 / 3) * Nodes[no].len;
	}
    }
  pt1_density();
#endif

#ifdef BLACK_HOLES
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 5)
	  P[i].Hsml = All.SofteningTable[5];
    }
#endif

#if  defined(LT_STELLAREVOLUTION)
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 4)
	  P[i].Hsml = All.SofteningTable[4];
    }
#endif

  density();

}


void assign_unique_ids(void)
{
  int i, *numpartlist;
  MyIDType idfirst;

  numpartlist = (int *) mymalloc("numpartlist", NTask * sizeof(int));

  MPI_Allgather(&NumPart, 1, MPI_INT, numpartlist, 1, MPI_INT, MYMPI_COMM_WORLD);

  idfirst = 1;

  for(i = 0; i < ThisTask; i++)
    idfirst += numpartlist[i];

  for(i = 0; i < NumPart; i++)
    {
      P[i].ID = idfirst;
      idfirst++;
    }

  myfree(numpartlist);
}

#ifdef ADAPTGRAVSOFT

void ags_setup_smoothinglengths(void)
{
  int i, no, p;

#ifdef  AGS_OUTPUTGRAVSOFT
  /* If the gravitational softening lengths are stored in the outputfiles, then we'd rather read those than
   * recompute them from scratch when the code is restarted from a snapshot.*/
  //well, actually...
  if(RestartFlag == 0 || RestartFlag == 2)
    //  if(RestartFlag == 0)
    {
      for(i = 0; i < NumPart; i++)
	{
	  no = Father[i];

	  while(10 * All.AGS_DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

	  P[i].AGS_Hsml =
	    pow(3.0 / (4 * M_PI) * All.AGS_DesNumNgb * P[i].Mass / Nodes[no].u.d.mass,
		1.0 / 3) * Nodes[no].len;

	}
    }
  ags_density();

#else

  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	{
	  no = Father[i];

	  while(10 * All.AGS_DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

	  P[i].AGS_Hsml =
	    pow(3.0 / (4 * M_PI) * All.AGS_DesNumNgb * P[i].Mass / Nodes[no].u.d.mass,
		1.0 / 3) * Nodes[no].len;

	}
    }
  ags_density();

#endif

}

#endif

void test_id_uniqueness(void)
{
  int i;
  double t0, t1;
  MyIDType *ids, *ids_first;

  if(ThisTask == 0)
    {
      printf("Testing ID uniqueness...\n");
      fflush(stdout);
    }

  if(NumPart == 0)
    {
      printf("need at least one particle per cpu\n");
      endrun(8);
    }

  t0 = second();

#ifndef SPH_BND_PARTICLES
  ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
  ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));

  for(i = 0; i < NumPart; i++)
    ids[i] = P[i].ID;

#ifdef ALTERNATIVE_PSORT
  init_sort_ID(ids, NumPart);
#else
  parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);
#endif

  for(i = 1; i < NumPart; i++)
    if(ids[i] == ids[i - 1])
      {
#ifdef LONGIDS
	printf("non-unique ID=%d%09d found on task=%d (i=%d NumPart=%d)\n",
	       (int) (ids[i] / 1000000000), (int) (ids[i] % 1000000000), ThisTask, i, NumPart);

#else
	printf("non-unique ID=%d found on task=%d   (i=%d NumPart=%d)\n", (int) ids[i], ThisTask, i, NumPart);
#endif
	endrun(12);
      }

  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MYMPI_COMM_WORLD);

  if(ThisTask < NTask - 1)
    if(ids[NumPart - 1] == ids_first[ThisTask + 1])
      {
	printf("non-unique ID=%d found on task=%d\n", (int) ids[NumPart - 1], ThisTask);
	endrun(13);
      }

  myfree(ids_first);
  myfree(ids);
#endif

  t1 = second();

  if(ThisTask == 0)
    {
      printf("success.  took=%g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
}

int compare_IDs(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}
