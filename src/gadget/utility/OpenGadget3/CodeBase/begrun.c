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
#include <sys/types.h>
#include <sys/stat.h>

#ifndef _WIN32
#include <unistd.h>
#endif

#ifndef GADGET3_IO_LIB
#include <gsl/gsl_rng.h>
#endif

#include <ctype.h>

#ifndef _WIN32
#include <libgen.h>		/* libraries to easily build path of included files */
#endif
#include <string.h>


#include "allvars.h"
#include "proto.h"
#include "../Blackholes/blackhole_begrun.h"


#ifdef GDE_BIGFLOAT
#if (GDE_BIGFLOAT==1)
#include <iostream>
#endif
#endif

#ifdef OPENACC
#include <openacc.h>
#endif

#define MAX_PATH_LENGTH 300	//maximum path length of a included parameter file


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;


#if defined(OPENACC)		//||defined(ACC)
  // basic GPU selection - assumes (number GPUs / node >= MPI ranks / node )
  const int num_devices = acc_get_num_devices(acc_device_nvidia);
  if(num_devices > 0)
    {
      int mydevice = ThisTask % num_devices;	// 1 device per MPI rank, less MPI_ranks then GPUs / node
      acc_set_device_num(mydevice, acc_device_nvidia);
      acc_init(acc_device_nvidia);
      printf("OpenACC init: Mpi Rank %d/%d using GPU %d/%d\n", ThisTask, NTask, mydevice, num_devices);
#pragma omp parallel
      {
	acc_set_device_num(mydevice, acc_device_nvidia);
	acc_init(acc_device_nvidia);
      }
    }
  else
    {
      printf("OpenACC init:  Mpi Rank %d/%d No GPU found.\n", ThisTask, NTask);
    }
  if(ThisTask == 0)
    printf("\nUsing %d GPU(s) per MPI rank via OpenACC\n", acc_get_num_devices(acc_device_nvidia));
#endif


#ifdef _OPENMP
  int tid;
#endif
  if(ThisTask == 0)
    {
      /*    printf("\nThis is P-Gadget, version `%s', svn-revision `%s'.\n", GADGETVERSION, svn_version()); */
      printf("\nThis is Open3Gadget, version %s.\n", GADGETVERSION);
      printf("\nRunning on %d MPI tasks.\n", NTask);
#ifdef _OPENMP
#pragma omp parallel private(tid)
      {
#pragma omp master
	{
	  printf("\nUsing %d OpenMP threads\n", omp_get_num_threads());
	  printf("\nUsing %d omp threads for fftw\n", omp_get_num_threads());
	}
	/*
	   tid = omp_get_thread_num();
	   printf("Hello from thread = %d\n", tid);
	 */
      }
#endif

      printf("\nCode was compiled with settings:\n\n");

#ifndef GADGET3_IO_LIB
      output_compile_time_options();
#endif

#ifndef DOUBLEPRECISION
      printf("\nCode was compile with no extra precission settings\n");
#else
      printf("\nCode was compiled with DOUBLEPRECISSION set to %d\n", DOUBLEPRECISION);
#endif
      printf("   Size MyFloat is                  %d  [bytes]\n", (int) sizeof(MyFloat));
      printf("   Size MyDouble is                 %d  [bytes]\n", (int) sizeof(MyDouble));
      printf("   Size MyIDType is                 %d  [bytes]\n", (int) sizeof(MyIDType));
      printf("   Size MyLongDouble is             %d  [bytes]\n", (int) sizeof(MyLongDouble));
      printf("   Size MyAtLeastDouble is          %d  [bytes]\n\n", (int) sizeof(MyAtLeastDouble));

      printf("   Size of All structure            %d  [bytes]\n", (int) sizeof(struct global_data_all_processes));

      printf("   Size of particle structure       %d  [bytes]\n", (int) sizeof(struct particle_data));
      printf("   Size of sph particle structure   %d  [bytes]\n", (int) sizeof(struct sph_particle_data));
#ifdef LT_STELLAREVOLUTION
      printf("   Size of metal particle structure %d  [bytes]\n\n", (int) sizeof(struct met_particle_data));
#endif
#ifdef BLACK_HOLES
      printf("   Size of BH particle structure    %d  [bytes]\n\n", (int) sizeof(struct bh_particle_data));
#endif
      printf("\n");

#ifdef GDE_BIGFLOAT
#if (GDE_BIGFLOAT==1)
      MyBig tmp;
      tmp.SetMax();
      std::cout << "GDE_BIGFLOAT: Largest floating point number for GDE: " << tmp << std::endl;

#endif
#endif
    }

#if defined(X86FIX) && defined(SOFTDOUBLEDOUBLE)
  x86_fix();			/* disable 80bit treatment of internal FPU registers in favour of proper IEEE 64bit double precision arithmetic */
#endif

  read_parameter_file(ParameterFile, NULL, NULL, NULL, 0);	/* ... read in parameters for this run, parameter 1 implies to set default values. */

  mymalloc_init();

#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

  set_units();

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif

#ifdef OUTPUT_LIGHTCONES
  lc_step_data_present = 0;
  read_step_file();
#endif

#ifdef GROWING_DISK_POTENTIAL
  growing_disk_init();
#endif


#ifdef COOLING
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
#ifndef LT_STELLAREVOLUTION

#ifndef GADGET3_IO_LIB
  InitCool();
#endif

#endif
#endif

#ifdef UM_CHEMISTRY
  printf("Initialize UVB rates to be interpolated, Task %d...\n", ThisTask);
  gJH0uvb = gJHe0uvb = gJHepuvb = 0;
  epsH0uvb = epsHe0uvb = epsHepuvb = 0;

  HIshielding = 1;
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
#ifndef GADGET3_IO_LIB
  printf("Initialize chemistry..\n");
  InitChem();
#endif
#endif

#if defined(SFR) && !defined(LT_STELLAREVOLUTION) && !defined(GM_MUPPI) && !defined(EB_SFR_MAGNETIC)
#ifndef GADGET3_IO_LIB
  init_clouds();
#endif
#endif

#ifdef PERIODIC
#ifndef GADGET3_IO_LIB
  ewald_init();
#endif
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
  inverse_boxSize = 1. / boxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
  inverse_boxSize_X = 1. / boxSize_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
  inverse_boxSize_Y = 1. / boxSize_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
  inverse_boxSize_Z = 1. / boxSize_Z;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
  All.ViscSource = All.ViscSource0 / log((GAMMA + 1) / (GAMMA - 1));
#ifdef AB_ART_VISC
  All.DecayTime = 1 / All.DecayLength;
#else
  All.DecayTime = 1 / All.DecayLength * sqrt((GAMMA - 1) / (2 * GAMMA));
#endif
#endif

#ifndef GADGET3_IO_LIB
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

  set_random_numbers();
#endif

#ifdef PMGRID
#ifndef ADAPTGRAVSOFT
  if(RestartFlag != 3 && RestartFlag != 4)
#endif
    long_range_init();
#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
  long_range_init();
#endif
#endif

#ifdef SUBFIND
  GrNr =-1;
#endif

#ifdef SIDM
  sidm_Init_CrossSection();
#endif

#if defined(fSIDM) || defined(rSIDM)
  msidm_init();
#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(ThisTask == 0)
    printf("Openinglog files ...\n");
  open_outputfiles();

  if(RestartFlag == -1 || RestartFlag == 0 || RestartFlag == 2 || RestartFlag == 3 || RestartFlag == 4
     || RestartFlag == 5 || RestartFlag == 6)
    {
#ifndef GADGET3_IO_LIB
      init();			/* ... read in initial model */
#endif
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

#ifndef GADGET3_IO_LIB
      restart(RestartFlag);	/* ... read restart file. Note: This also resets
				   all variables in the struct `All'.
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.
				 */

      set_random_numbers();
#endif

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MinGasHsmlFractional = all.MinGasHsmlFractional;
      All.MinGasTemp = all.MinGasTemp;

#ifdef MAXHSML
      All.MaxHsml = all.MaxHsml;
#endif
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesPerSnapshot = all.NumFilesPerSnapshot;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);
      memcpy(All.OutputListFlag, all.OutputListFlag, sizeof(char) * All.OutputListLength);

#if defined(LT_STELLAREVOLUTION)
      memcpy(All.MetalsCheckSum, all.MetalsCheckSum, sizeof(long double) * LT_NMetP);
#endif

#if defined(GM_MUPPI) && defined(LT_STELLAREVOLUTION)
      memcpy(All.beta_sf, all.beta_sf, sizeof(double) * 10);
      memcpy(All.E_SN_51, all.E_SN_51, sizeof(double) * 10);
      memcpy(All.f_rest, all.f_rest, sizeof(double) * 10);

      All.MuppiOn = all.MuppiOn;
      All.MuppiDebugOn = all.MuppiDebugOn;
      All.CountParticlesInConeOn = all.CountParticlesInConeOn;
      All.StellarKinFB2On = all.StellarKinFB2On;
      All.StellarKinFB2OutputOn = all.StellarKinFB2OutputOn;
#endif


#ifdef TIME_DEP_ART_VISC
      All.ViscSource = all.ViscSource;
      All.ViscSource0 = all.ViscSource0;
      All.DecayTime = all.DecayTime;
      All.DecayLength = all.DecayLength;
      All.AlphaMin = all.AlphaMin;
#endif

#ifdef ARTIFICIAL_CONDUCTIVITY
      All.ArtCondConstant = all.ArtCondConstant;
#ifdef TIME_DEP_ART_COND
      All.ArtCondMin = all.ArtCondMin;
#endif
#endif

#if defined(MAGNETIC_DISSIPATION)
      All.ArtMagDispConst = all.ArtMagDispConst;
#ifdef TIME_DEP_MAGN_DISP
      All.ArtMagDispMin = all.ArtMagDispMin;
      All.ArtMagDispSource = all.ArtMagDispSource;
      All.ArtMagDispTime = all.ArtMagDispTime;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
      All.DivBcleanParabolicSigma = all.DivBcleanParabolicSigma;
      All.DivBcleanHyperbolicSigma = all.DivBcleanHyperbolicSigma;
      All.DivBcleanQ = all.DivBcleanQ;
#endif

#ifdef BLACK_HOLES
      change_param_values_blackhole(&All, &all, &FdParamChangeLog);
#endif

#ifdef LT_STELLAREVOLUTION
      All.MaxChemSpreadL = all.MaxChemSpreadL;
#ifdef LT_ADD_GAL_TO_SUB
      strcpy(All.BC_SED_Path, all.BC_SED_Path);
#endif
#endif

#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

#ifdef ADAPTGRAVSOFT
      /* Allow the tolerance over the number of neighbours to vary during the run:
       * If it was initially set to a very strict value, convergence in ngb-iteration may at some point fail */
      All.AGS_MaxNumNgbDeviation = all.AGS_MaxNumNgbDeviation;
#endif

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.TimebinFile, all.TimebinFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);


#ifdef RELAXOBJECT
      All.RelaxBaseFac = all.RelaxBaseFac;
#endif


      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);
    }

#ifdef LT_STELLAREVOLUTION
  if(All.TestSuite && (ThisTask == 0))
    {
      printf(">>>>>>>>>>>>>>>>>>>>>\n" "  testing SE\n" ">>>>>>>>>>>>>>>>>>>>>\n");

#ifndef GADGET3_IO_LIB
      TestStellarEvolution();
#endif
    }

  if (RestartFlag != 1) {
	  // if All.MetalsCheckSum was not read from restart file
	  // we initialize it here
#ifndef GADGET3_IO_LIB
	  get_metals_checksum(1, All.MetalsCheckSum);
#endif
  }
#endif

  char contfname[1000];
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  //  open_outputfiles();

#ifdef LMB_SPECTRAL_CRs
  // write headers for output files
  if(All.CR_Details > 0)
    {
      write_cr_file_headers();
    }
#endif

#ifdef PMGRID
  long_range_init_regionsize();

#ifdef ADAPTGRAVSOFT
  /* Setting the maximum\minimum allowed gravitational softening length */
  All.AGS_MaxSoft[0] = AGS_MAXSOFT * All.Asmth[0];
  All.AGS_MinSoft[0] = AGS_MINSOFT * All.Asmth[0];
  if(ThisTask == 0)
    printf("*ADAPTGRAVSOFT* MinSoft, MaxSoft= %g %g\n", All.AGS_MinSoft[0], All.AGS_MaxSoft[0]);
#ifdef PLACEHIGHRESREGION
  All.AGS_MaxSoft[1] = AGS_MAXSOFT * All.Asmth[1];
  All.AGS_MinSoft[1] = AGS_MINSOFT * All.Asmth[1];
  if(ThisTask == 0)
    {
      printf("*ADAPTGRAVSOFT* for Highresolution Region: MinSoft, MaxSoft= %g %g\n", All.AGS_MinSoft[1],
	     All.AGS_MaxSoft[1]);
      printf("ASmth[0] = %.2f ASmth[1] = %.2f RCut[0] = %.2f RCut[1] = %.2f\n", All.Asmth[0], All.Asmth[1],
	     All.Rcut[0], All.Rcut[1]);
    }
#endif
#endif

#endif

#ifndef GADGET3_IO_LIB
  reconstruct_timebins();
#endif


#ifdef TWODIMS
  int i;

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


#ifndef GADGET3_IO_LIB
  if(All.ComovingIntegrationOn)
    init_drift_table();

  if(RestartFlag == 2)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
  else if(RestartFlag == 1)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);

#endif

#if defined(OUTPUT_LONGRANGE_POTENTIAL)
  All.Ti_nextpotdump = All.Ti_nextoutput;
  All.dumped_after_snap = 0;	/* need to do first potential dump */
#endif

  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT
  double coulomb_log;
#endif
#endif
#ifdef STATICNFW
  double Mtot;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;
#ifdef TIMEDEPGRAV
  All.Gini = All.G;
  All.G = All.Gini * dGfak(All.TimeBegin);
#endif

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

#ifdef DISTORTIONTENSORPS
  /* 5.609589206e23 is the factor to convert from g to GeV/c^2, the rest comes from All.UnitDensity_in_cgs */
  All.UnitDensity_in_Gev_per_cm3 = 5.609589206e23 / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g;
#endif
  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;
#ifdef KSPACE_NEUTRINOS_2
  /*Set OmegaNu from the neutrino mass. */
  All.OmegaNu = OmegaNu(1);
#endif

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
#ifdef DISTORTIONTENSORPS
      printf("Annihilation radiation units:\n");
      printf("UnitDensity_in_Gev_per_cm3 = %g\n", All.UnitDensity_in_Gev_per_cm3);
#endif
#ifdef INCLUDE_RADIATION
#ifdef KSPACE_NEUTRINOS_2
      printf("Omega_R = %g Omega_Nu = %g\n", OMEGAR, All.OmegaNu);
#else
      printf("Omega_R = %g\n", OMEGAR);
#endif
#endif

      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#ifdef LT_STELLAREVOLUTION
  All.MinChemTimeStep /= 86400 * 365 * 1e9;
  /*
     Add here Spreading Length in the case of fixed distance
   */
  All.SnIaEgy /= All.UnitEnergy_in_cgs;
  All.SnIIEgy /= All.UnitEnergy_in_cgs;

#ifdef LT_HOT_EJECTA
  /* All.EgySpecEjecta is supposed to be given in paramfile as Km/sec        */
  /* then we have to transform in erg / g through the 1e10 conversion factor */
  All.EgySpecEjecta = pow(All.EgySpecEjecta, 2) * 1e10 * All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

  Hyd = LT_NMet - 1;
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();

#ifndef GADGET3_IO_LIB
  init_SN();
#endif

#endif


#if defined(SFR) && ( (!defined(LT_STELLAREVOLUTION)) ||  defined(GM_MUPPI) )

#ifndef GADGET3_IO_LIB
  set_units_sfr();
#endif

#endif


#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)

#ifdef NAVIERSTOKES
  /* Braginskii-Spitzer shear viscosity parametrization */
  /* mu = 0.406 * m_p^0.5 * (k_b* T)^(5/2) / e^4 / logLambda  [g/cm/s] */
  /* eta = frac * mu */

  meanweight = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* assuming full ionization */

#ifdef NAVIERSTOKES_CONSTANT
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 0.406 * pow(PROTONMASS, 0.5) * pow(BOLTZMANN * All.ShearViscosityTemperature, 5. / 2.) / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;	/*in cgs units */

  if(ThisTask == 0)
    printf("Constant shear viscosity in cgs units: eta = %g\n", All.NavierStokes_ShearViscosity);

  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g / All.HubbleParam;	/* in internal code units */

  if(ThisTask == 0)
    printf("Constant shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);

#else
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 0.406 * pow(PROTONMASS, 0.5) * pow((meanweight * PROTONMASS * GAMMA_MINUS1), 5. / 2.) / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;	/*in cgs units */
  /*T = mu*m_p*(gamma-1)/k_b * E * UnitEnergy/UnitMass */

  All.NavierStokes_ShearViscosity *= pow((All.UnitEnergy_in_cgs / All.UnitMass_in_g), 5. / 2.);	/* now energy can be multiplied later in the internal code units */
  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g / All.HubbleParam;	/* in internal code units */

  if(ThisTask == 0)
    printf("Variable shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);

#endif

#ifdef NAVIERSTOKES_BULK
  if(ThisTask == 0)
    printf("Costant bulk viscosity in internal code units: zeta = %g\n", All.NavierStokes_BulkViscosity);
#endif

#ifdef NAVIERSTOKES_VISCOSITY_SATURATION
  /* calculate ion mean free path assuming complete ionization:                                                                                                                                                                                
     ion mean free path for hydrogen is similar to that of helium,                                                                                                                                                                             
     thus we calculate only for hydrogen */
  /* l_i = 3^(3/2)*(k*T)^2 / (4*\pi^(1/2)*ni*(Z*e)^4*lnL) */

  All.IonMeanFreePath = pow(3.0, 1.5) / (4.0 * sqrt(M_PI) * pow(ELECTRONCHARGE, 4) * LOG_LAMBDA);

  All.IonMeanFreePath *= pow(meanweight * PROTONMASS * GAMMA_MINUS1, 2) * pow((All.UnitEnergy_in_cgs / All.UnitMass_in_g), 2);	/*kT -> u */

  All.IonMeanFreePath /= (HYDROGEN_MASSFRAC / PROTONMASS) *
    (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
  /* n_H = rho * Hfr / mp *//* now is cgs units *///changed / to * in front of the unitdensity                                                                                                                                                 

  All.IonMeanFreePath *= All.HubbleParam / All.UnitLength_in_cm;
  /* in internal code units */
#endif

#endif



#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT

  meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  /* assuming full ionization */

  coulomb_log = 37.8;
  /* according to Sarazin's book */

  All.ConductionCoeff *=
    (1.84e-5 / coulomb_log * pow(meanweight / k_B * GAMMA_MINUS1, 2.5) * erg / (s * deg * cm));
  /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003
   * ( ApJ 582:162-169, Eq. (5) )
   */

  /* Note: Because we replace \nabla(T) in the conduction equation with
   * \nable(u), our conduction coefficient is not the usual kappa, but
   * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with
   * another factor of (meanweight / k_B * GAMMA_MINUS1).
   */
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  All.ConductionCoeff /= All.HubbleParam;

#ifdef CONDUCTION_SATURATION
  All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
    / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)
    / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam)
    * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;
#endif

#endif /* CONDUCTION_CONSTANT */

#ifdef CONDUCTION_INCLUDEMAGNETIC_PERPENDICULAR
  // constants from omega_g * tau
  All.ConductionPerpendicularConstants =
    pow(meanweight * GAMMA_MINUS1, 1.5) * pow(ELECTRONMASS * All.HubbleParam / All.UnitMass_in_g, 0.5);
  All.ConductionPerpendicularConstants /= 16.0 * M_PI * M_PI * pow(ELECTRONCHARGE * All.HubbleParam / pow(All.UnitEnergy_in_cgs * All.UnitLength_in_cm, 0.5), 3.0) * 2.9979e10 / All.UnitVelocity_in_cm_per_s;	//c
#endif

#endif /* CONDUCTION */

#ifdef LMB_SPECTRAL_CRs

#ifndef LMB_DPP
  All.CR_Dpp *= All.UnitTime_in_s;
#endif

#ifdef LMB_SPECTRAL_CRs_DIFFUSION
  All.CR_Kappa_10k /= (All.UnitVelocity_in_cm_per_s * All.UnitLength_in_cm);
#endif // LMB_SPECTRAL_CRs_SPACIAL_DIFFUSION

#endif // LMB_SPECTRAL_CRs

#ifdef STATICNFW
  R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs = R200 / NFW_C;
  Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200 = 10 * All.Hubble * R200;
  if(ThisTask == 0)
    printf("V200= %g\n", V200);

  fac = 1.0;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);

  /* fac = M200 / Mtot */
  fac = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
#endif
}

#ifdef STATICNFW
/*! auxiliary function for static NFW potential
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  /* use unsoftened NFW if NFW_Eps=0 */
  if(NFW_Eps > 0.0)
    if(R > Rs * NFW_C)
      R = Rs * NFW_C;

  if(NFW_Eps > 0.0)
    {
      return fac * 4 * M_PI * RhoCrit * Dc *
	(-
	 (Rs * Rs * Rs *
	  (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) +
	   NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1)) + (Rs * Rs * Rs * (Rs -
													NFW_Eps
													* Rs -
													(2 *
													 NFW_Eps
													 -
													 1) *
													(R +
													 Rs) *
													log(R
													    +
													    Rs)
													+
													NFW_Eps
													*
													NFW_Eps
													* (R +
													   Rs)
													*
													log(R
													    +
													    NFW_Eps
													    *
													    Rs)))
	 / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
    }
  else				/* analytic NFW */
    {
      return fac * 4 * M_PI * RhoCrit * Dc *
	(-(Rs * Rs * Rs * (1 + log(Rs))) + Rs * Rs * Rs * (Rs + (R + Rs) * log(R + Rs)) / (R + Rs));
    }
}
#endif



/*!  This function opens various log-files that report on the status and
 *   performance of the simulation. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  if(ThisTask == 0)
    mkdir(All.OutputDir, 02755);
  MPI_Barrier(MYMPI_COMM_WORLD);

#ifdef BLACK_HOLES
  /* Note: This is done by everyone */
  if(ThisTask == 0)
    {
      sprintf(buf, "%sblackhole_details", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MYMPI_COMM_WORLD);

  if(All.BlackHoleDetails >= 1)
    {
      sprintf(buf, "%sblackhole_details/blackhole_details_%d.txt", All.OutputDir, ThisTask);
      if(!(FdBlackHolesDetails = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
    }
#endif

#ifdef MV_GM_AGNMUPPI_OUTPUT
  sprintf(buf, "%s%s", All.OutputDir, "EgyAGN_global.txt");
  if(!(FdEgyAGNtot = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdEgyAGNtot, "# AGN energy output: global information\n");
      fprintf(FdEgyAGNtot,
	      "# 1.time       2.TotE_AGNM_tot       3.TotE_AGNM_cool       4.TotE_AGNM_sfr_c       5.TotE_AGNM_sfr_h        6.TotE_AGNM_sfr_E      7.TotE_AGNM_sfr_E_aft        8.TotE_AGNM_sfr_BHEcoldD       9.Tot_n_AGNM           10.Tot_n_AGNcool\n");
      fflush(FdEgyAGNtot);
    }
#endif


#ifdef LMB_SPECTRAL_CRs
/* Note: This is done by everyone */
  if(ThisTask == 0)
    {
      sprintf(buf, "%sCRs_details", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MYMPI_COMM_WORLD);

  if(All.CR_Details >= 1)
    {
      sprintf(buf, "%sCRs_details/CRs_details_adiabatic_%d.txt", All.OutputDir, ThisTask);
      if(!(FdCRsDetailsAdiabatic = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      sprintf(buf, "%sCRs_details/CRs_details_radiative_%d.txt", All.OutputDir, ThisTask);
      if(!(FdCRsDetailsRadiative = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}

      sprintf(buf, "%sCRs_details/CRs_details_Dpp_%d.txt", All.OutputDir, ThisTask);
      if(!(FdCRsDetailsDpp = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}

#ifdef JD_SHOCK
      sprintf(buf, "%sCRs_details/CRs_details_shock_injection_%d.txt", All.OutputDir, ThisTask);
      if(!(FdCRsDetailsShockInjection = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
#endif

#if (defined(LMB_CR_SFR_INJECTION) || defined(LMB_SPECTRAL_CRs_SN_SEEDING))
      sprintf(buf, "%sCRs_details/CRs_details_sfr_injection_%d.txt", All.OutputDir, ThisTask);
      if(!(FdCRsDetailsSfrInjection = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
#endif
    }

#endif


#ifdef LT_TRACK_WINDS
  sprintf(buf, "%s%s%d%s", All.OutputDir, "trckwinds.", ThisTask, ".txt");
  if(!(FdTrackW = fopen(buf, mode)))
    {
      printf("[%d] error in opening file '%s'\n", ThisTask, buf);
      endrun(1);
    }
#endif

#ifdef MV_GM_STELLAR_KIN_FB2_OUTPUT
  sprintf(buf, "%s%s_%03d.txt", All.OutputDir, "winds", ThisTask);
  if(!(FdWind = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdWind, "# Statistics of wind gas particles\n");
      fprintf(FdWind,
	      "#    time  x y z vx vy vz mass rho pres temp id vkick/All.Time SphP[i].E_kin SphP[i].xkin SphP[i].ykin SphP[i].zkin DelayTime\n");
      fflush(FdWind);
    }
#endif

#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION_OUTPUT
  if(ThisTask == 0)
    {
      sprintf(buf, "%s%s_%03d.txt", All.OutputDir, "SF_Molecular_fraction", ThisTask);
      if(!(FdMolFrac = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      else
	{
	  fprintf(FdMolFrac, "# SF molecular fractions - output\n");
	  fprintf(FdMolFrac,
		  "# 1.case   2.time    3.particleID    4.H2_frac_Krumholz    5.H2_frac_Blitz&R    6.metallicity    7.number_dens hot(1)/cold(2)    8.temperature hot(1)/cold(2)\n");
	  fflush(FdMolFrac);
	}
    }
#endif

#ifdef MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z_OUTPUT
  sprintf(buf, "%s%s_%03d.txt", All.OutputDir, "Hypernovae", ThisTask);
  if(!(FdHypernovae = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdHypernovae,
	      "# Info about MP star-forming particles with low Z injecting higher E_out from early SNe\n");
      fprintf(FdHypernovae,
	      "# 1.time    2.ID     3.x     4.y     5.z      6.rho     7.temp     8.mass     9.metallicity     10.E_out\n");
      fflush(FdHypernovae);
    }
#endif

#ifdef MV_GM_AGNMUPPI_OUTPUT
  sprintf(buf, "%s%s_%03d.txt", All.OutputDir, "EgyAGN", ThisTask);
  if(!(FdEgyAGN = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdEgyAGN, "# AGN energy output\n");
      fprintf(FdEgyAGN,
	      "# 1.time    2.x    3.y    4.z    5.vx    6.vy    7.vz    8.mass    9.rho    10.pres    11.temp    12.id    13.flag    14.Injected_BH_Energy    15.MultiPhase    16.DelayTime      17.CouplingBHEnHot    18.CouplingBHEnCold\n");
      fflush(FdEgyAGN);
    }
#endif

#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
  sprintf(buf, "%s%s_%03d.txt", All.OutputDir, "CovFacAGN", ThisTask);
  if(!(FdCovFacAGN = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdCovFacAGN, "# AGN covering factors-coupling factors output for MP particles\n");
      fprintf(FdCovFacAGN,
	      "# 1.time    2.ID    3.filling_fac_hot    4.hsml_Muppi      5.Injected_BH_Energy      6.SphP[i].CouplingBHEnHot     7.SphP[i].CouplingBHEnCold    8.SphP[i].M_h    9.SphP[i].M_c    10.SphP[i].T_h    11.n_h    12.n_c    13.SphP[i].Sfr     14.hsml \n");
      fflush(FdCovFacAGN);
    }
#endif

#if defined (MV_GM_AGNMUPPI) && !defined(MV_AGNMUPPI_COOLING_OFF_OUTPUT)
  sprintf(buf, "%s%s_%03d.txt", All.OutputDir, "AGNmuppi", ThisTask);
  if(!(FdAGNmuppi = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdAGNmuppi, "# MV_GM_AGNMUPPI_OUTPUT - Properties of MP particles when AGNmuppi is on\n");
      fprintf(FdAGNmuppi,
	      "# 1.time    2.timestep   3.Egy_to_Hot_phase_from_AGNfb    4.Egy_to_Cold_phase_from_AGNfb    5.M_c_in_timestep    6.M_h_in_timestep   7.T_h_in_timestep    8.M_c_fin_timestep    9.M_h_fin_timestep   10.T_h_fin_timestep   11.InitialBHEcold   12.y[4]    13.tcool    14.Sph_density     15.Pressione_sph_preAGN_(cgs/K_B)      16.Pressione_postAGN\n");
      fflush(FdAGNmuppi);
    }
#endif

#if defined(MV_AGNMUPPI_COOLING_OFF) && defined(MV_AGNMUPPI_COOLING_OFF_OUTPUT)
  sprintf(buf, "%s%s_%03d.txt", All.OutputDir, "CoolOff", ThisTask);
  if(!(FdCoolOffAGNmuppi = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdCoolOffAGNmuppi,
	      "# MV_AGNMUPPI_COOLING_OFF - Properties of MP particles with cooling (eventually) switched off\n");
      fprintf(FdCoolOffAGNmuppi,
	      "# 1.time    2.timestep    3.timeCoolOff   4.Egy_to_Hot_phase_from_AGNfb    5.Egy_to_Cold_phase_from_AGNfb    6.M_c_in_timestep    7.M_h_in_timestep    8.T_h_in_timestep    9.M_c_fin_timestep    10.M_h_fin_timestep   11.T_h_fin_timestep   12.InitialBHEcold   13.y[4]    14.tcool    15.Sph_density     16.Pressione_sph_preAGN_(cgs/K_B)      17.Pressione_postAGN\n");
      fflush(FdCoolOffAGNmuppi);
    }
#endif



  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimebinFile);
  if(!(FdTimebin = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "paramchanges.txt");
  if(!(FdParamChangeLog = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  fprintf(FdBalance, "\n");
  fprintf(FdBalance, "Treewalk1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK1],
	  CPU_SymbolImbalance[CPU_TREEWALK1]);
  fprintf(FdBalance, "Treewalk2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK2],
	  CPU_SymbolImbalance[CPU_TREEWALK2]);
  fprintf(FdBalance, "Treewait1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT1],
	  CPU_SymbolImbalance[CPU_TREEWAIT1]);
  fprintf(FdBalance, "Treewait2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT2],
	  CPU_SymbolImbalance[CPU_TREEWAIT2]);
  fprintf(FdBalance, "Treesend       = '%c' / '%c'\n", CPU_Symbol[CPU_TREESEND],
	  CPU_SymbolImbalance[CPU_TREESEND]);
  fprintf(FdBalance, "Treerecv       = '%c' / '%c'\n", CPU_Symbol[CPU_TREERECV],
	  CPU_SymbolImbalance[CPU_TREERECV]);
  fprintf(FdBalance, "Treebuild      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEBUILD],
	  CPU_SymbolImbalance[CPU_TREEBUILD]);
  fprintf(FdBalance, "Treeupdate     = '%c' / '%c'\n", CPU_Symbol[CPU_TREEUPDATE],
	  CPU_SymbolImbalance[CPU_TREEUPDATE]);
  fprintf(FdBalance, "Treehmaxupdate = '%c' / '%c'\n", CPU_Symbol[CPU_TREEHMAXUPDATE],
	  CPU_SymbolImbalance[CPU_TREEHMAXUPDATE]);
  fprintf(FdBalance, "Treemisc =       '%c' / '%c'\n", CPU_Symbol[CPU_TREEMISC],
	  CPU_SymbolImbalance[CPU_TREEMISC]);
  fprintf(FdBalance, "Domain decomp  = '%c' / '%c'\n", CPU_Symbol[CPU_DOMAIN],
	  CPU_SymbolImbalance[CPU_DOMAIN]);
  fprintf(FdBalance, "Density compute= '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMPUTE],
	  CPU_SymbolImbalance[CPU_DENSCOMPUTE]);
  fprintf(FdBalance, "Density imbal  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSWAIT],
	  CPU_SymbolImbalance[CPU_DENSWAIT]);
  fprintf(FdBalance, "Density commu  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMM],
	  CPU_SymbolImbalance[CPU_DENSCOMM]);
  fprintf(FdBalance, "Density misc   = '%c' / '%c'\n", CPU_Symbol[CPU_DENSMISC],
	  CPU_SymbolImbalance[CPU_DENSMISC]);
  fprintf(FdBalance, "Hydro compute  = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMPUTE],
	  CPU_SymbolImbalance[CPU_HYDCOMPUTE]);
  fprintf(FdBalance, "Hydro imbalance= '%c' / '%c'\n", CPU_Symbol[CPU_HYDWAIT],
	  CPU_SymbolImbalance[CPU_HYDWAIT]);
  fprintf(FdBalance, "Hydro comm     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMM],
	  CPU_SymbolImbalance[CPU_HYDCOMM]);
  fprintf(FdBalance, "Hydro misc     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDMISC],
	  CPU_SymbolImbalance[CPU_HYDMISC]);
#if GADGET_HYDRO == HYDRO_MFM
  fprintf(FdBalance, "MFM grad comp  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMGRADCOMPUTE],
	  CPU_SymbolImbalance[CPU_MFMGRADCOMPUTE]);
  fprintf(FdBalance, "MFM grad imbal = '%c' / '%c'\n", CPU_Symbol[CPU_MFMGRADWAIT],
	  CPU_SymbolImbalance[CPU_MFMGRADWAIT]);
  fprintf(FdBalance, "MFM grad comm  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMGRADCOMM],
	  CPU_SymbolImbalance[CPU_MFMGRADCOMM]);
  fprintf(FdBalance, "MFM grad misc  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMGRADMISC],
	  CPU_SymbolImbalance[CPU_MFMGRADMISC]);
  fprintf(FdBalance, "MFM lim comp  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMLIMCOMPUTE],
	  CPU_SymbolImbalance[CPU_MFMLIMCOMPUTE]);
  fprintf(FdBalance, "MFM lim imbal = '%c' / '%c'\n", CPU_Symbol[CPU_MFMLIMWAIT],
	  CPU_SymbolImbalance[CPU_MFMLIMWAIT]);
  fprintf(FdBalance, "MFM lim comm  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMLIMCOMM],
	  CPU_SymbolImbalance[CPU_MFMLIMCOMM]);
  fprintf(FdBalance, "MFM lim misc  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMLIMMISC],
	  CPU_SymbolImbalance[CPU_MFMLIMMISC]);
  fprintf(FdBalance, "MFM flux comp  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMFLUXCOMPUTE],
	  CPU_SymbolImbalance[CPU_MFMFLUXCOMPUTE]);
  fprintf(FdBalance, "MFM flux imbal = '%c' / '%c'\n", CPU_Symbol[CPU_MFMFLUXWAIT],
	  CPU_SymbolImbalance[CPU_MFMFLUXWAIT]);
  fprintf(FdBalance, "MFM flux comm  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMFLUXCOMM],
	  CPU_SymbolImbalance[CPU_MFMFLUXCOMM]);
  fprintf(FdBalance, "MFM flux misc  = '%c' / '%c'\n", CPU_Symbol[CPU_MFMFLUXMISC],
	  CPU_SymbolImbalance[CPU_MFMFLUXMISC]);
#endif
  fprintf(FdBalance, "Drifts         = '%c' / '%c'\n", CPU_Symbol[CPU_DRIFT], CPU_SymbolImbalance[CPU_DRIFT]);
  fprintf(FdBalance, "Blackhole      = '%c' / '%c'\n", CPU_Symbol[CPU_BLACKHOLES],
	  CPU_SymbolImbalance[CPU_BLACKHOLES]);
  fprintf(FdBalance, "Kicks          = '%c' / '%c'\n", CPU_Symbol[CPU_TIMELINE],
	  CPU_SymbolImbalance[CPU_TIMELINE]);
  fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL],
	  CPU_SymbolImbalance[CPU_POTENTIAL]);
  fprintf(FdBalance, "PM             = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
  fprintf(FdBalance, "Peano-Hilbert  = '%c' / '%c'\n", CPU_Symbol[CPU_PEANO], CPU_SymbolImbalance[CPU_PEANO]);
  fprintf(FdBalance, "Cooling & SFR  = '%c' / '%c'\n", CPU_Symbol[CPU_COOLINGSFR],
	  CPU_SymbolImbalance[CPU_COOLINGSFR]);
  fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT],
	  CPU_SymbolImbalance[CPU_SNAPSHOT]);
  fprintf(FdBalance, "FoF            = '%c' / '%c'\n", CPU_Symbol[CPU_FOF], CPU_SymbolImbalance[CPU_FOF]);
  fprintf(FdBalance, "Miscellaneous  = '%c' / '%c'\n", CPU_Symbol[CPU_MISC], CPU_SymbolImbalance[CPU_MISC]);
  fprintf(FdBalance, "\n");

#ifdef SCFPOTENTIAL
  sprintf(buf, "%s%s", All.OutputDir, "scf_coeff.txt");
  if(!(FdSCF = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#ifdef GM_MUPPI
  sprintf(buf, "%s%s", All.OutputDir, "exit.txt");
  if(!(FdExit = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      fprintf(FdExit, "# Statistics of exits from Multi-Phase\n");
      fprintf(FdExit,
	      "#    time  density    clock   freeze    spawn  forceSP     wind    --    depl      gsl   sf run      mnc      tms   --    hot     cold    alminef ineff\n");
    }
#ifdef GM_COUNT_PARTICLES_IN_CONE
  sprintf(buf, "%s%s", All.OutputDir, "particles_in_cone.txt");
  if(!(FdCone = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif
#endif


#endif


#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LMB_SPECTRAL_CRs
  sprintf(buf, "%s%s", All.OutputDir, "CRs.txt");
  if(!(FdCRs = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#endif

#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif

#ifdef XXLINFO
  sprintf(buf, "%s%s", All.OutputDir, "xxl.txt");
  if(!(FdXXL = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdXXL, "nstep time ");
#ifdef MAGNETIC
	  fprintf(FdXXL, "<|B|> ");
#ifdef TRACEDIVB
	  fprintf(FdXXL, "max(divB) ");
#endif
#endif
#ifdef TIME_DEP_ART_VISC
	  fprintf(FdXXL, "<alpha> ");
#endif
	  fprintf(FdXXL, "\n");
	  fflush(FdXXL);
	}
    }
#endif

#ifdef DARKENERGY
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdDE, "nstep time H(a) ");
#ifndef TIMEDEPDE
	  fprintf(FdDE, "w0 Omega_L ");
#else
	  fprintf(FdDE, "w(a) Omega_L ");
#endif
#ifdef TIMEDEPGRAV
	  fprintf(FdDE, "dH dG ");
#endif
	  fprintf(FdDE, "\n");
	  fflush(FdDE);
	}
    }
#endif

#ifdef LT_EXTEGY_INFO
  sprintf(buf, "%s%s", All.OutputDir, "extegy.txt");
  if(!(FdExtEgy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_SEv_INFO_DETAILS_onSPREAD
  sprintf(buf, "%s%s", All.OutputDir, "spinfo.txt");
  if(!(FdSPinfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_STELLAREVOLUTION

#ifdef LT_SEv_INFO

#ifdef LT_SEvDbg
  sprintf(buf, "%s%s", All.OutputDir, "met_sumcheck.txt");
  if(!(FdMetSumCheck = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

  sprintf(buf, "%s%s", All.OutputDir, "metals.txt");
  if(!(FdMetals = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "sn.txt");
  if(!(FdSn = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "sn_lost.txt");
  if(!(FdSnLost = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#ifdef WINDS
  sprintf(buf, "%s%s", All.OutputDir, "winds.txt");
  if(!(FdWinds = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#endif
#endif


}




/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
#ifdef BLACK_HOLES
  if(All.BlackHoleDetails >= 1)
    fclose(FdBlackHolesDetails);	/* needs to be done by everyone */
#endif

#ifdef LMB_SPECTRAL_CRs
  if(All.CR_Details >= 1)
    {
      fclose(FdCRsDetailsAdiabatic);	/* needs to be done by everyone */
      fclose(FdCRsDetailsRadiative);	/* needs to be done by everyone */
      fclose(FdCRsDetailsDpp);	/* needs to be done by everyone */
      fclose(FdCRsDetailsShockInjection);	/* needs to be done by everyone */
#ifdef LMB_CR_SFR_INJECTION
      fclose(FdCRsDetailsSfrInjection);	/* needs to be done by everyone */
#endif
    }

#endif

#ifdef LT_TRACK_WINDS
  fclose(FdTrackW);
#endif

#ifdef MV_GM_STELLAR_KIN_FB2_OUTPUT
  fclose(FdWind);
#endif

#ifdef MV_GM_AGNMUPPI
  fclose(FdAGNmuppi);
#ifdef MV_GM_AGNMUPPI_OUTPUT
  fclose(FdEgyAGN);
#endif
#endif
#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
  fclose(FdCovFacAGN);
#endif
#if defined(MV_AGNMUPPI_COOLING_OFF) && defined(MV_AGNMUPPI_COOLING_OFF_OUTPUT)
  fclose(FdCoolOffAGNmuppi);
#endif
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION_OUTPUT
  fclose(FdMolFrac);
#endif
#ifdef MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z_OUTPUT
  fclose(FdHypernovae);
#endif


  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdTimebin);
  fclose(FdBalance);
  fclose(FdParamChangeLog);

#ifdef SCFPOTENTIAL
  fclose(FdSCF);
#endif

#ifdef SFR
  fclose(FdSfr);
#ifdef GM_MUPPI
  fclose(FdExit);
#ifdef GM_COUNT_PARTICLES_IN_CONE
  fclose(FdCone);
#endif
#endif
#endif


#ifdef BLACK_HOLES
  fclose(FdBlackHoles);
#endif

#ifdef LMB_SPECTRAL_CRs
  fclose(FdCRs);
#endif

#ifdef XXLINFO
  fclose(FdXXL);
#endif

#ifdef DARKENERGY
  fclose(FdDE);
#endif



#ifdef LT_SEv_INFO_DETAILS_onSPREAD
  fclose(FdSPinfo);
#endif

#ifdef LT_EXTEGY_INFO
  fclose(FdExtEgy);
#endif

#ifdef LT_STELLAREVOLUTION

#ifdef LT_SEv_INFO
#ifdef LT_SEvDbg
  fclose(FdMetSumCheck);
#endif
  fclose(FdSn);
  fclose(FdSnLost);
  fclose(FdMetals);
#ifdef WINDS
  fclose(FdWinds);
#endif
#endif
#endif


}





/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
#define REAL 1
#define STRING 2
#define INT 3

void read_parameter_file(char *fname, char *tag[], void **addr, int *id, int nt)
{




  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j;



  int included = tag != NULL;



  if(!included)
    {
      tag = (char **) malloc(MAXTAGS * sizeof(char *));
      addr = (void **) malloc(MAXTAGS * sizeof(void *));
      id = (int *) malloc(MAXTAGS * sizeof(int));
      for(int i = 0; i < MAXTAGS; i++)
	{
	  tag[i] = (char *) malloc(MAXTAGLEN * sizeof(char));
	}
    }


  int pnum, errorFlag = 0;

  if(!included)
    {
      All.StarformationOn = 0;	/* defaults */
    }



  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      if(!included)
	{

	  nt = 0;

	  strcpy(tag[nt], "InitCondFile");
	  addr[nt] = All.InitCondFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "OutputDir");
	  addr[nt] = All.OutputDir;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "SnapshotFileBase");
	  addr[nt] = All.SnapshotFileBase;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "EnergyFile");
	  addr[nt] = All.EnergyFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "CpuFile");
	  addr[nt] = All.CpuFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "InfoFile");
	  addr[nt] = All.InfoFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "TimingsFile");
	  addr[nt] = All.TimingsFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "TimebinFile");
	  addr[nt] = All.TimebinFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "RestartFile");
	  addr[nt] = All.RestartFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "ResubmitCommand");
	  addr[nt] = All.ResubmitCommand;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "OutputListFilename");
	  addr[nt] = All.OutputListFilename;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "OutputListOn");
	  addr[nt] = &All.OutputListOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "Omega0");
	  addr[nt] = &All.Omega0;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "OmegaBaryon");
	  addr[nt] = &All.OmegaBaryon;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "OmegaLambda");
	  addr[nt] = &All.OmegaLambda;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "HubbleParam");
	  addr[nt] = &All.HubbleParam;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BoxSize");
	  addr[nt] = &All.BoxSize;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "PeriodicBoundariesOn");
	  addr[nt] = &All.PeriodicBoundariesOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "MaxMemSize");
	  addr[nt] = &All.MaxMemSize;
	  id[nt++] = INT;

	  strcpy(tag[nt], "TimeOfFirstSnapshot");
	  addr[nt] = &All.TimeOfFirstSnapshot;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CpuTimeBetRestartFile");
	  addr[nt] = &All.CpuTimeBetRestartFile;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "TimeBetStatistics");
	  addr[nt] = &All.TimeBetStatistics;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "TimeBegin");
	  addr[nt] = &All.TimeBegin;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "TimeMax");
	  addr[nt] = &All.TimeMax;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "TimeBetSnapshot");
	  addr[nt] = &All.TimeBetSnapshot;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
	  addr[nt] = &All.UnitVelocity_in_cm_per_s;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "UnitLength_in_cm");
	  addr[nt] = &All.UnitLength_in_cm;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "UnitMass_in_g");
	  addr[nt] = &All.UnitMass_in_g;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "TreeDomainUpdateFrequency");
	  addr[nt] = &All.TreeDomainUpdateFrequency;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "ErrTolIntAccuracy");
	  addr[nt] = &All.ErrTolIntAccuracy;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "ErrTolTheta");
	  addr[nt] = &All.ErrTolTheta;
	  id[nt++] = REAL;

#ifdef SIDM
	  strcpy(tag[nt], "CrossSectionPerMass_in_cgs");
	  addr[nt] = &All.CrossSectionPerMass_in_cgs;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "PeakSigma");
	  addr[nt] = &All.PeakSigma;
	  id[nt++] = REAL;
#ifdef SIDM_INELASTIC
	  strcpy(tag[nt], "Splitting");
	  addr[nt] = &All.PeakSigma;
	  id[nt++] = REAL;
#endif
#endif

#ifdef fSIDM
	  strcpy(tag[nt], "FSIDM_CrossSectionPerMass_in_cgs");
	  addr[nt] = &All.FSIDM_CrossSectionPerMass_in_cgs;
	  id[nt++] = REAL;
#ifdef mSIDM_TIMESTEP
	  strcpy(tag[nt], "FSIDM_tfact");
	  addr[nt] = &All.FSIDM_tfact;
	  id[nt++] = REAL;
#endif
#endif
#ifdef rSIDM
	  strcpy(tag[nt], "RSIDM_CrossSectionPerMass_in_cgs");
	  addr[nt] = &All.RSIDM_CrossSectionPerMass_in_cgs;
	  id[nt++] = REAL;
#ifdef mSIDM_TIMESTEP
	  strcpy(tag[nt], "RSIDM_tfact");
	  addr[nt] = &All.RSIDM_tfact;
	  id[nt++] = REAL;
#endif
#ifdef rSIDM_DISSIPATION
	  strcpy(tag[nt], "RSIDM_f_diss");
	  addr[nt] = &All.RSIDM_f_diss;
	  id[nt++] = REAL;
#endif
#endif
#ifdef mSIDM_TIMESTEP
#if !defined(mSIDM_VDEP0)
	  strcpy(tag[nt], "MSIDM_OmegaInit");
	  addr[nt] = &All.MSIDM_OmegaInit;
	  id[nt++] = REAL;
#endif
#endif
#if defined(fSIDM) || defined(rSIDM)
	  strcpy(tag[nt], "MSIDM_SEED");
	  addr[nt] = &All.MSIDM_SEED;
	  id[nt++] = INT;

	  strcpy(tag[nt], "MSIDM_MaxIntOrder");
	  addr[nt] = &All.MSIDM_MaxIntOrder;
	  id[nt++] = INT;

	  strcpy(tag[nt], "MSIDM_IntOrder");
	  addr[nt] = &All.MSIDM_IntOrder;
	  id[nt++] = INT;

#ifdef mSIDM_VDEP0
	  strcpy(tag[nt], "MSIDM_vdep_w");
	  addr[nt] = &All.MSIDM_vdep_w;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MSIDM_vdep_alpha");
	  addr[nt] = &All.MSIDM_vdep_alpha;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MSIDM_vdep_beta");
	  addr[nt] = &All.MSIDM_vdep_beta;
	  id[nt++] = REAL;
#endif

#ifndef ADAPTGRAVSOFT
	  strcpy(tag[nt], "DesNumNgb_pt1");
	  addr[nt] = &All.DesNumNgb_pt1;
	  id[nt++] = INT;

	  strcpy(tag[nt], "MaxNumNgbDeviation_pt1");
	  addr[nt] = &All.MaxNumNgbDeviation_pt1;
	  id[nt++] = REAL;
#endif
#endif

#ifdef WINDTUNNEL
	  strcpy(tag[nt], "WindGrid");
	  addr[nt] = &All.WindGrid;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindVel");
	  addr[nt] = &All.WindVel;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindDmean");
	  addr[nt] = &All.WindDmean;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindEntr");
	  addr[nt] = &All.WindEntr;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindDens");
	  addr[nt] = &All.WindDens;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindPmass");
	  addr[nt] = &All.WindPmass;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindSizeOfInjectionRegion");
	  addr[nt] = &All.WindSizeOfInjectionRegion;
	  id[nt++] = REAL;

#endif

#ifdef SUBFIND
	  strcpy(tag[nt], "ErrTolThetaSubfind");
	  addr[nt] = &All.ErrTolThetaSubfind;
	  id[nt++] = REAL;
#endif

	  strcpy(tag[nt], "ErrTolForceAcc");
	  addr[nt] = &All.ErrTolForceAcc;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MinGasHsmlFractional");
	  addr[nt] = &All.MinGasHsmlFractional;
	  id[nt++] = REAL;

#if defined(ISOTHERM_EQS)
	  strcpy(tag[nt], "IsoSoundSpeed");
	  addr[nt] = &All.IsoSoundSpeed;
	  id[nt++] = REAL;
#endif


#ifdef MAXHSML
	  strcpy(tag[nt], "MaxHsml");
	  addr[nt] = &All.MaxHsml;
	  id[nt++] = REAL;
#endif

	  strcpy(tag[nt], "MaxSizeTimestep");
	  addr[nt] = &All.MaxSizeTimestep;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MinSizeTimestep");
	  addr[nt] = &All.MinSizeTimestep;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MaxRMSDisplacementFac");
	  addr[nt] = &All.MaxRMSDisplacementFac;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "ArtBulkViscConst");
	  addr[nt] = &All.ArtBulkViscConst;
	  id[nt++] = REAL;

#ifdef ARTIFICIAL_CONDUCTIVITY
	  strcpy(tag[nt], "ArtCondConstant");
	  addr[nt] = &All.ArtCondConstant;
	  id[nt++] = REAL;
#ifdef TIME_DEP_ART_COND
	  strcpy(tag[nt], "ArtCondMin");
	  addr[nt] = &All.ArtCondMin;
	  id[nt++] = REAL;
#endif
#endif

	  strcpy(tag[nt], "CourantFac");
	  addr[nt] = &All.CourantFac;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "DesNumNgb");
	  addr[nt] = &All.DesNumNgb;
	  id[nt++] = INT;

#ifdef NEUTRINOS
	  strcpy(tag[nt], "Time_tree_on_nu");
	  addr[nt] = &All.Time_tree_on_nu;
	  id[nt++] = REAL;
	  All.tree_nu_on = 0;
#endif

#ifdef KSPACE_NEUTRINOS
	  strcpy(tag[nt], "KspaceNeutrinoSeed");
	  addr[nt] = &All.KspaceNeutrinoSeed;
	  id[nt++] = INT;

	  strcpy(tag[nt], "Nsample");
	  addr[nt] = &All.Nsample;
	  id[nt++] = INT;

	  strcpy(tag[nt], "SphereMode");
	  addr[nt] = &All.SphereMode;
	  id[nt++] = INT;

	  strcpy(tag[nt], "KspaceDirWithTransferfunctions");
	  addr[nt] = All.KspaceDirWithTransferfunctions;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "KspaceBaseNameTransferfunctions");
	  addr[nt] = All.KspaceBaseNameTransferfunctions;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "PrimordialIndex");
	  addr[nt] = &All.PrimordialIndex;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "Sigma8");
	  addr[nt] = &All.Sigma8;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "OmegaNu");
	  addr[nt] = &All.OmegaNu;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
	  addr[nt] = &All.InputSpectrum_UnitLength_in_cm;
	  id[nt++] = REAL;

#endif

#if defined KSPACE_NEUTRINOS_2
	  strcpy(tag[nt], "KspaceTransferFunction");
	  addr[nt] = All.KspaceTransferFunction;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "TimeTransfer");
	  addr[nt] = &All.TimeTransfer;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "OmegaBaryonCAMB");
	  addr[nt] = &All.OmegaBaryonCAMB;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
	  addr[nt] = &All.InputSpectrum_UnitLength_in_cm;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MNue");
	  addr[nt] = &(All.MNu[0]);
	  id[nt++] = REAL;
	  strcpy(tag[nt], "MNum");
	  addr[nt] = &(All.MNu[1]);
	  id[nt++] = REAL;
	  strcpy(tag[nt], "MNut");
	  addr[nt] = &(All.MNu[2]);
	  id[nt++] = REAL;
#endif



#ifdef SUBFIND
	  strcpy(tag[nt], "DesLinkNgb");
	  addr[nt] = &All.DesLinkNgb;
	  id[nt++] = INT;
#endif

	  strcpy(tag[nt], "MaxNumNgbDeviation");
	  addr[nt] = &All.MaxNumNgbDeviation;
	  id[nt++] = REAL;

#ifdef START_WITH_EXTRA_NGBDEV
	  strcpy(tag[nt], "MaxNumNgbDeviationStart");
	  addr[nt] = &All.MaxNumNgbDeviationStart;
	  id[nt++] = REAL;
#endif

	  strcpy(tag[nt], "ComovingIntegrationOn");
	  addr[nt] = &All.ComovingIntegrationOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "ICFormat");
	  addr[nt] = &All.ICFormat;
	  id[nt++] = INT;

	  strcpy(tag[nt], "SnapFormat");
	  addr[nt] = &All.SnapFormat;
	  id[nt++] = INT;

	  strcpy(tag[nt], "NumFilesPerSnapshot");
	  addr[nt] = &All.NumFilesPerSnapshot;
	  id[nt++] = INT;

	  strcpy(tag[nt], "NumFilesWrittenInParallel");
	  addr[nt] = &All.NumFilesWrittenInParallel;
	  id[nt++] = INT;

	  strcpy(tag[nt], "ResubmitOn");
	  addr[nt] = &All.ResubmitOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "CoolingOn");
	  addr[nt] = &All.CoolingOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "StarformationOn");
	  addr[nt] = &All.StarformationOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "TypeOfTimestepCriterion");
	  addr[nt] = &All.TypeOfTimestepCriterion;
	  id[nt++] = INT;

	  strcpy(tag[nt], "TypeOfOpeningCriterion");
	  addr[nt] = &All.TypeOfOpeningCriterion;
	  id[nt++] = INT;

	  strcpy(tag[nt], "TimeLimitCPU");
	  addr[nt] = &All.TimeLimitCPU;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningHalo");
	  addr[nt] = &All.SofteningHalo;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningDisk");
	  addr[nt] = &All.SofteningDisk;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningBulge");
	  addr[nt] = &All.SofteningBulge;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningGas");
	  addr[nt] = &All.SofteningGas;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningStars");
	  addr[nt] = &All.SofteningStars;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningBndry");
	  addr[nt] = &All.SofteningBndry;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningHaloMaxPhys");
	  addr[nt] = &All.SofteningHaloMaxPhys;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningDiskMaxPhys");
	  addr[nt] = &All.SofteningDiskMaxPhys;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningBulgeMaxPhys");
	  addr[nt] = &All.SofteningBulgeMaxPhys;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningGasMaxPhys");
	  addr[nt] = &All.SofteningGasMaxPhys;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningStarsMaxPhys");
	  addr[nt] = &All.SofteningStarsMaxPhys;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SofteningBndryMaxPhys");
	  addr[nt] = &All.SofteningBndryMaxPhys;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BufferSize");
	  addr[nt] = &All.BufferSize;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "PartAllocFactor");
	  addr[nt] = &All.PartAllocFactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "GravityConstantInternal");
	  addr[nt] = &All.GravityConstantInternal;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "InitGasTemp");
	  addr[nt] = &All.InitGasTemp;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MinGasTemp");
	  addr[nt] = &All.MinGasTemp;
	  id[nt++] = REAL;

#ifdef DISTORTIONTENSORPS
	  strcpy(tag[nt], "TidalCorrection");
	  addr[nt] = &All.TidalCorrection;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "DM_velocity_dispersion");
	  addr[nt] = &All.DM_velocity_dispersion;
	  id[nt++] = REAL;
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
	  strcpy(tag[nt], "ReferenceGasMass");
	  addr[nt] = &All.ReferenceGasMass;
	  id[nt++] = REAL;
#endif


#ifdef NAVIERSTOKES
	  strcpy(tag[nt], "FractionSpitzerViscosity");
	  addr[nt] = &All.FractionSpitzerViscosity;
	  id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_CONSTANT
	  strcpy(tag[nt], "ShearViscosityTemperature");
	  addr[nt] = &All.ShearViscosityTemperature;
	  id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_BULK
	  strcpy(tag[nt], "NavierStokes_BulkViscosity");
	  addr[nt] = &All.NavierStokes_BulkViscosity;
	  id[nt++] = REAL;
#endif


#ifdef CHEMISTRY
	  strcpy(tag[nt], "Epsilon");
	  addr[nt] = &All.Epsilon;
	  id[nt++] = REAL;
#endif


#ifdef CONDUCTION
	  strcpy(tag[nt], "ConductionEfficiency");
	  addr[nt] = &All.ConductionCoeff;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MaxSizeConductionStep");
	  addr[nt] = &All.MaxSizeConductionStep;
	  id[nt++] = REAL;
#endif

#ifdef LMB_SPECTRAL_CRs_DIFFUSION
	  strcpy(tag[nt], "CR_DiffusionKappa10k");
	  addr[nt] = &All.CR_Kappa_10k;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CR_AlphaKappa");
	  addr[nt] = &All.CR_AlphaKappa;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MaxSizeCRDiffusionStep");
	  addr[nt] = &All.MaxSizeCRDiffusionStep;
	  id[nt++] = REAL;
#endif

	  strcpy(tag[nt], "LevelOfStrickness");
	  addr[nt] = &All.LevelOfStrickness;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHolesOn");
	  addr[nt] = &All.BlackHolesOn;
	  id[nt++] = INT;

#ifdef BLACK_HOLES

	  strcpy(tag[nt], "BlackHoleThermalFeedbackOn");
	  addr[nt] = &All.BlackHoleThermalFeedbackOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleKineticFeedbackOn");
	  addr[nt] = &All.BlackHoleKineticFeedbackOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "massDMpart");
	  addr[nt] = &All.massDMpart;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "TimeBetOnTheFlyFoF");
	  addr[nt] = &All.TimeBetOnTheFlyFoF;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleFoFMinMassForNewSeed");
	  addr[nt] = &All.MinFoFMassForNewSeed;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleSeedMass");
	  addr[nt] = &All.SeedBlackHoleMass;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleSeedStarMassFraction");
	  addr[nt] = &All.BlackHoleSeedStarMassFraction;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleSeedDMFraction");
	  addr[nt] = &All.BlackHoleSeedDMFraction;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleSeedGasFraction");
	  addr[nt] = &All.BlackHoleSeedGasFraction;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleFormationFactor");
	  addr[nt] = &All.BHfactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleNgbFactor");
	  addr[nt] = &All.BlackHoleNgbFactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleMaxAccretionRadius");
	  addr[nt] = &All.BlackHoleMaxAccretionRadius;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleSwallowGasOn");
	  addr[nt] = &All.BlackHoleSwallowGasOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleIgnoreMomentum");
	  addr[nt] = &All.BlackHoleIgnoreMomentum;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleAccretionSlicesOn");
	  addr[nt] = &All.BlackHoleAccretionSlicesOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleDetails");
	  addr[nt] = &All.BlackHoleDetails;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleRepositioningOn");
	  addr[nt] = &All.BlackHoleRepositioningOn;
	  id[nt++] = INT;


	  strcpy(tag[nt], "BlackHoleHotColdAccretionOn");
	  addr[nt] = &All.BlackHoleHotColdAccretionOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleAccretionFactor");
	  addr[nt] = &All.BlackHoleAccretionFactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleColdAccretionFactor");
	  addr[nt] = &All.BlackHoleColdAccretionFactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleVelocityInBondiOn");
	  addr[nt] = &All.BlackHoleVelocityInBondiOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleVariableAccretionFactorOn");
	  addr[nt] = &All.BlackHoleVariableAccretionFactorOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleVariableAccretionSlope");
	  addr[nt] = &All.BlackHoleVariableAccretionSlope;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleVariableEfficiencyOn");
	  addr[nt] = &All.BlackHoleVariableEfficiencyOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleRadiativeEfficiency");
	  addr[nt] = &All.BlackHoleRadiativeEfficiency;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadiativeEfficiencySlope");
	  addr[nt] = &All.BlackHoleRadiativeEfficiencySlope;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadiativeEfficiencyNorm");
	  addr[nt] = &All.BlackHoleRadiativeEfficiencyNorm;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadiativeEfficiencySlopeMdot");
	  addr[nt] = &All.BlackHoleRadiativeEfficiencySlopeMdot;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadiativeEfficiencySlopeMbh");
	  addr[nt] = &All.BlackHoleRadiativeEfficiencySlopeMbh;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHolePhysMdotUpperLimit");
	  addr[nt] = &All.BlackHolePhysMdotUpperLimit;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadiativeEfficiencyMax");
	  addr[nt] = &All.BlackHoleRadiativeEfficiencyMax;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadiativeEfficiencyMin");
	  addr[nt] = &All.BlackHoleRadiativeEfficiencyMin;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleEddingtonFactor");
	  addr[nt] = &All.BlackHoleEddingtonFactor;
	  id[nt++] = REAL;


	  strcpy(tag[nt], "BlackHoleOutflowModelOn");
	  addr[nt] = &All.BlackHoleOutflowModelOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleFeedbackFactor");
	  addr[nt] = &All.BlackHoleFeedbackFactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadioModeFeedbackBoost");
	  addr[nt] = &All.BlackHoleRadioModeFeedbackBoost;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleRadioModeTreshold");
	  addr[nt] = &All.BlackHoleRadioModeTreshold;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleLimitFeedbackOn");
	  addr[nt] = &All.BlackHoleLimitFeedbackOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleLimitMaxTemp");
	  addr[nt] = &All.BlackHoleLimitMaxTemp;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleColdGasTemperatureTresh");
	  addr[nt] = &All.BlackHoleColdGasTemperatureTresh;
	  id[nt++] = REAL;


	  strcpy(tag[nt], "BlackHoleFrictionForceOn");
	  addr[nt] = &All.BlackHoleFrictionForceOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleFrictionForceDynaicsOn");
	  addr[nt] = &All.BlackHoleFrictionForceDynaicsOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "BlackHoleFrictionForceMaxCorrect");
	  addr[nt] = &All.BlackHoleFrictionForceMaxCorrect;
	  id[nt++] = REAL;


	  strcpy(tag[nt], "BlackHoleMergeCsndFrac");
	  addr[nt] = &All.BlackHoleMergeCsndFrac;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleMergeDistFrac");
	  addr[nt] = &All.BlackHoleMergeDistFrac;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BlackHoleMergeBindingFrac");
	  addr[nt] = &All.BlackHoleMergeBindingFrac;
	  id[nt++] = REAL;



#endif /* BLACK_HOLES */

#if defined (UM_CHEMISTRY) && defined (UM_CHEMISTRY_INISET)
	  /* read the composition from the parameter file */
	  strcpy(tag[nt], "START_elec");
	  addr[nt] = &All.Startelec;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HI");
	  addr[nt] = &All.StartHI;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HII");
	  addr[nt] = &All.StartHII;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HM");
	  addr[nt] = &All.StartHM;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HeI");
	  addr[nt] = &All.StartHeI;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HeII");
	  addr[nt] = &All.StartHeII;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HeIII");
	  addr[nt] = &All.StartHeIII;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_H2I");
	  addr[nt] = &All.StartH2I;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_H2II");
	  addr[nt] = &All.StartH2II;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HD");
	  addr[nt] = &All.StartHD;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_DI");
	  addr[nt] = &All.StartDI;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_DII");
	  addr[nt] = &All.StartDII;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "START_HeHII");
	  addr[nt] = &All.StartHeHII;
	  id[nt++] = REAL;
#endif

#ifdef SFR
	  strcpy(tag[nt], "CritOverDensity");
	  addr[nt] = &All.CritOverDensity;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CritPhysDensity");
	  addr[nt] = &All.CritPhysDensity;
	  id[nt++] = REAL;

#ifndef LT_STELLAREVOLUTION
	  strcpy(tag[nt], "FactorSN");
	  addr[nt] = &All.FactorSN;
	  id[nt++] = REAL;
	  strcpy(tag[nt], "FactorEVP");
	  addr[nt] = &All.FactorEVP;
	  id[nt++] = REAL;
#endif

	  strcpy(tag[nt], "TempSupernova");
	  addr[nt] = &All.TempSupernova;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "TempClouds");
	  addr[nt] = &All.TempClouds;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MaxSfrTimescale");
	  addr[nt] = &All.MaxSfrTimescale;
	  id[nt++] = REAL;

#ifdef GM_MUPPI
	  strcpy(tag[nt], "FracEgyIn");
	  addr[nt] = &All.FracEgyIn;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "FracEgyOut");
	  addr[nt] = &All.FracEgyOut;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "FracEgyKin");
	  addr[nt] = &All.FracEgyKin;
	  id[nt++] = REAL;

#if defined (MV_GM_AGNMUPPI) && !defined(MV_GM_COVER_FACT_MCLS)
	  strcpy(tag[nt], "CouplingBHEnHot");
	  addr[nt] = &All.CouplingBHEnHot;
	  id[nt++] = REAL;
	  strcpy(tag[nt], "CouplingBHEnCold");
	  addr[nt] = &All.CouplingBHEnCold;
	  id[nt++] = REAL;
#endif


	  /* Here we have our new flags for MUPPI */
	  /* I initialize them to "On" since if we are here we
	     want to use Muppi */
	  All.MuppiOn = 1;
	  All.MuppiDebugOn = 0;
	  All.CountParticlesInConeOn = 1;	/* GM: this will need to default at 0 */
	  All.StellarKinFB2On = 1;
	  All.StellarKinFB2OutputOn = 1;	/* GM: this will need to default at 0 */

	  strcpy(tag[nt], "MuppiOn");
	  addr[nt] = &All.MuppiOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "MuppiDebugOn");
	  addr[nt] = &All.MuppiDebugOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "CountParticlesInConeOn");
	  addr[nt] = &All.CountParticlesInConeOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "StellarKinFB2");
	  addr[nt] = &All.StellarKinFB2On;
	  id[nt++] = INT;

	  strcpy(tag[nt], "StellarKinFB2Output");
	  addr[nt] = &All.StellarKinFB2OutputOn;
	  id[nt++] = INT;
#endif


#ifdef WINDS
	  strcpy(tag[nt], "WindEfficiency");
	  addr[nt] = &All.WindEfficiency;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindEnergyFraction");
	  addr[nt] = &All.WindEnergyFraction;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindFreeTravelMaxTimeFactor");
	  addr[nt] = &All.WindFreeTravelMaxTimeFactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "WindFreeTravelDensFac");
	  addr[nt] = &All.WindFreeTravelDensFac;
	  id[nt++] = REAL;

#ifdef VARIABLE_WINDS
	  strcpy(tag[nt], "VariableWindVelFactor");
	  addr[nt] = &All.VariableWindVelFactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "VariableWindSpecMomentum");
	  addr[nt] = &All.VariableWindSpecMomentum;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "HaloConcentrationNorm");
	  addr[nt] = &All.HaloConcentrationNorm;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "HaloConcentrationSlope");
	  addr[nt] = &All.HaloConcentrationSlope;
	  id[nt++] = REAL;
#endif
#endif /* End Winds */

#endif

#if defined(SNIA_HEATING)
	  strcpy(tag[nt], "SnIaHeatingRate");
	  addr[nt] = &All.SnIaHeatingRate;
	  id[nt++] = REAL;
#endif



#ifdef DARKENERGY
#ifndef TIMEDEPDE
	  strcpy(tag[nt], "DarkEnergyParam");
	  addr[nt] = &All.DarkEnergyParam;
	  id[nt++] = REAL;
#endif
#endif

#ifdef RESCALEVINI
	  strcpy(tag[nt], "VelIniScale");
	  addr[nt] = &All.VelIniScale;
	  id[nt++] = REAL;
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
	  strcpy(tag[nt], "DarkEnergyFile");
	  addr[nt] = All.DarkEnergyFile;
	  id[nt++] = STRING;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
	  strcpy(tag[nt], "ViscositySourceScaling");
	  addr[nt] = &All.ViscSource0;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "ViscosityDecayLength");
	  addr[nt] = &All.DecayLength;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "ViscosityAlphaMin");
	  addr[nt] = &All.AlphaMin;
	  id[nt++] = REAL;
#endif

#if defined(MAGNETIC_DISSIPATION)
	  strcpy(tag[nt], "ArtificialMagneticDissipationConstant");
	  addr[nt] = &All.ArtMagDispConst;
	  id[nt++] = REAL;

#ifdef TIME_DEP_MAGN_DISP
	  strcpy(tag[nt], "ArtificialMagneticDissipationMin");
	  addr[nt] = &All.ArtMagDispMin;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "ArtificialMagneticDissipationSource");
	  addr[nt] = &All.ArtMagDispSource;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "ArtificialMagneticDissipationDecaytime");
	  addr[nt] = &All.ArtMagDispTime;
	  id[nt++] = REAL;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
	  strcpy(tag[nt], "DivBcleaningParabolicSigma");
	  addr[nt] = &All.DivBcleanParabolicSigma;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "DivBcleaningHyperbolicSigma");
	  addr[nt] = &All.DivBcleanHyperbolicSigma;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "DivBcleaningQ");
	  addr[nt] = &All.DivBcleanQ;
	  id[nt++] = REAL;

#endif



#ifdef VSMOOTH
	  strcpy(tag[nt], "VSmoothScale");
	  addr[nt] = &All.VSmoothScale;
	  id[nt++] = REAL;
#endif

#ifdef MAGNETIC
#ifdef BINISET
	  strcpy(tag[nt], "BiniX");
	  addr[nt] = &All.BiniX;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BiniY");
	  addr[nt] = &All.BiniY;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BiniZ");
	  addr[nt] = &All.BiniZ;
	  id[nt++] = REAL;
#ifdef MAGNETICZERO
	  All.BiniX = 0.;
	  All.BiniY = 0.;
	  All.BiniZ = 0.;
#endif
#endif

#if defined(BSMOOTH)
	  strcpy(tag[nt], "BSmoothInt");
	  addr[nt] = &All.BSmoothInt;
	  id[nt++] = INT;
	  strcpy(tag[nt], "BSmoothFrac");
	  addr[nt] = &All.BSmoothFrac;
	  id[nt++] = REAL;

#ifdef SETMAINTIMESTEPCOUNT
	  strcpy(tag[nt], "MainTimestepCount");
	  addr[nt] = &All.MainTimestepCountIni;
	  id[nt++] = INT;
#endif
#endif
#ifdef BSMOOTH_TIME
	  All.BSmoothInt = 1;
	  All.BSmoothFrac = 1;
#endif

#ifdef MAGNETIC_SN_SEEDING
	  strcpy(tag[nt], "SnSeedRadius");
	  addr[nt] = &All.SnSeedRadius;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SnSeedBubble");
	  addr[nt] = &All.SnSeedBubble;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SnSeedField");
	  addr[nt] = &All.SnSeedField;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SnSeedSoftening");
	  addr[nt] = &All.SnSeedSoftening;
	  id[nt++] = REAL;
#endif

#ifdef MAGNETIC_DIFFUSION
	  strcpy(tag[nt], "MagneticEta");
	  addr[nt] = &All.MagneticEta;
	  id[nt++] = REAL;
#endif
#ifdef MAGNETIC_DIFFUSION_LIMIT
	  strcpy(tag[nt], "MagneticDiffSpeed");
	  addr[nt] = &All.MagneticDiffSpeed;
	  id[nt++] = REAL;
#endif
#endif /* MAGNETIC */



#ifdef RELAXOBJECT
	  strcpy(tag[nt], "RelaxBaseFac");
	  addr[nt] = &All.RelaxBaseFac;
	  id[nt++] = REAL;
#endif

	  strcpy(tag[nt], "SpectralCRsOn");
	  addr[nt] = &All.SpectralCRsOn;
	  id[nt++] = INT;

#ifdef LMB_SPECTRAL_CRs

	  strcpy(tag[nt], "CR_ProtonsOn");
	  addr[nt] = &All.CR_ProtonsOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "CR_ElectronsOn");
	  addr[nt] = &All.CR_ElectronsOn;
	  id[nt++] = INT;

	  strcpy(tag[nt], "CR_Details");
	  addr[nt] = &All.CR_Details;
	  id[nt++] = INT;

	  strcpy(tag[nt], "CR_TrackParticle");
	  addr[nt] = &All.CR_TrackParticle;
	  id[nt++] = INT;

	  strcpy(tag[nt], "CR_AdiabaticHighUpdateRate");
	  addr[nt] = &All.CR_AdiabaticHighUpdateRate;
	  id[nt++] = INT;

	  strcpy(tag[nt], "CR_SpektrumHighUpdateRate");
	  addr[nt] = &All.CR_SpektrumHighUpdateRate;
	  id[nt++] = INT;

#ifdef JD_SHOCK
	  strcpy(tag[nt], "CR_DSA_Model");
	  addr[nt] = &All.CR_DSA_Model;
	  id[nt++] = INT;

#ifdef LMB_CRs_CALCULATE_P_INJ
	  strcpy(tag[nt], "CR_Chi_pinj");
	  addr[nt] = &All.CR_Chi_pinj;
	  id[nt++] = REAL;
#endif // LMB_CRs_CALCULATE_P_INJ

#endif // JD_SHOCK

#ifdef LMB_CR_PROTONS
	  strcpy(tag[nt], "CRp_MinMomentum");
	  addr[nt] = &All.CRp_pmin;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CRp_MaxMomentum");
	  addr[nt] = &All.CRp_pmax;
	  id[nt++] = REAL;
#endif // LMB_CR_PROTONS
#ifdef LMB_CR_ELECTRONS
	  strcpy(tag[nt], "CRe_MinMomentum");
	  addr[nt] = &All.CRe_pmin;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CRe_MaxMomentum");
	  addr[nt] = &All.CRe_pmax;
	  id[nt++] = REAL;
#endif // LMB_CR_ELECTRONS

	  strcpy(tag[nt], "CR_MaxSlope");
	  addr[nt] = &All.CR_SlopeBig;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CR_MomentumDiffusionCoefficient");
	  addr[nt] = &All.CR_Dpp;
	  id[nt++] = REAL;

	  //      strcpy(tag[nt], "CR_EpsilonB");
	  // addr[nt] = &All.CR_EpsilonB;
	  // id[nt++] = REAL;

	  strcpy(tag[nt], "CR_SlopeSoftening");
	  addr[nt] = &All.CR_DSlope;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CR_ElectronInjectionRatio");
	  addr[nt] = &All.CRe_Ratio;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "CR_InjectionStartTime");
	  addr[nt] = &All.CR_T0Inj;
	  id[nt++] = REAL;

#ifdef LMB_SPECTRAL_CRs_SEED
#ifdef LMB_CR_PROTONS
	  strcpy(tag[nt], "CRp_InitSlope");
	  addr[nt] = &All.CRp_InitSlope;
	  id[nt++] = REAL;
#endif // LMB_CR_PROTONS
#ifdef LMB_CR_ELECTRONS
	  strcpy(tag[nt], "CRe_InitSlope");
	  addr[nt] = &All.CRe_InitSlope;
	  id[nt++] = REAL;
#endif // LMB_CR_ELECTRONS
#endif
#ifdef LMB_SPECTRAL_CRs_ARTIFICIAL_CONDUCTIVITY
	  strcpy(tag[nt], "CR_ArtCondConstant");
	  addr[nt] = &All.CR_ArtCondConstant;
	  id[nt++] = REAL;
#endif
#endif

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
	  strcpy(tag[nt], "Epsilon");
	  addr[nt] = &All.Epsilon;
	  id[nt++] = REAL;
#endif

#ifdef LT_STELLAREVOLUTION

#ifdef LT_STOP_COOL_BELOW_Z
	  strcpy(tag[nt], "CoolStop_redshift");
	  addr[nt] = &All.Below_this_redshift_stop_cooling;
	  id[nt++] = REAL;
#endif

#ifdef LT_TRACK_CONTRIBUTES
	  strcpy(tag[nt], "PackPowerBase");
	  addr[nt] = &PowerBase;
	  id[nt++] = INT;
#endif

#ifdef LT_POPIII
	  strcpy(tag[nt], "PopIII_IMF_idx");
	  addr[nt] = &All.PopIII_IMF_idx;
	  id[nt++] = INT;

	  strcpy(tag[nt], "PopIII_Zlimit");	/* Z threshold for PIII formation */
	  addr[nt] = &All.PopIII_Zlimit;	/* in units of solar abundance */
	  id[nt++] = REAL;
#endif

#ifdef LT_STARBURSTS
	  strcpy(tag[nt], "StarBurstCondition");
	  addr[nt] = &All.StarBurstCondition;
	  id[nt++] = INT;

	  strcpy(tag[nt], "SB_Density_Thresh");
	  addr[nt] = &All.SB_Density_Thresh;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SB_DEntropy_Thresh");
	  addr[nt] = &All.SB_DEntropy_Thresh;
	  id[nt++] = REAL;
#endif

#ifdef LT_SMOOTH_Z
#if defined(LT_SMOOTH_SIZE) && !defined(LT_SMOOTH_NGB)
	  strcpy(tag[nt], "SmoothRegionSize");
	  addr[nt] = &All.SmoothRegionSize;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SmoothRegionSizeMax");
	  addr[nt] = &All.SmoothRegionSizeMax;
	  id[nt++] = REAL;
#endif
#if defined(LT_SMOOTH_NGB) && !defined(LT_SMOOTH_SIZE)
	  strcpy(tag[nt], "DesNumNgnSmooth");
	  addr[nt] = &All.DesNumNgbSmooth;
	  id[nt++] = INT;
#endif
#endif

	  strcpy(tag[nt], "TestSuite");
	  addr[nt] = &All.TestSuite;
	  id[nt++] = INT;

	  /* the minimum number of neighbours used to spread
	   * metals from Stars.
	   */
	  strcpy(tag[nt], "InfNeighNum");
	  addr[nt] = &All.NeighInfNum;
	  id[nt++] = INT;

	  /* the desired number of neighbours used to spread
	   * metals from Stars.
	   */
	  strcpy(tag[nt], "DesNumNgbSN");
	  addr[nt] = &All.DesNumNgbSN;
	  id[nt++] = INT;

	  /* the allowed deviation on number of neighbours used
	   * to spread metals from Stars.
	   */
	  strcpy(tag[nt], "SpreadNumNgbDev");
	  addr[nt] = &All.SpreadNumNgbDev;
	  id[nt++] = INT;

	  /* note: this is used just to calculate All.MaxPartMet in
	   * read_file(); then, it should give a reasonable "average"
	   * of the generations of each specified IMFs.
	   *
	   * WARNING: this will be overridden by read_imfs() !!!
	   */
	  strcpy(tag[nt], "Generations");
	  addr[nt] = &All.Generations;
	  id[nt++] = INT;

	  /* the baryon fraction that you expect to end in stars
	   */
	  strcpy(tag[nt], "SFfactor");
	  addr[nt] = &All.SFfactor;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "IMFsFileName");
	  addr[nt] = All.IMFfilename;
	  id[nt++] = STRING;


	  strcpy(tag[nt], "SFsFileName");
	  addr[nt] = All.SFfilename;
	  id[nt++] = STRING;

	  /* the minimum mass for CC supernovae (usually 8 Msun) */
	  strcpy(tag[nt], "SnII_InfMass");
	  addr[nt] = &All.Mup;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MinBinSystemMass");
	  addr[nt] = &All.MBm;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MaxBinSystemMass");
	  addr[nt] = &All.MBM;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "BinarySystemFrac");
	  addr[nt] = &All.BinFrac;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MinStarMass_inBinSystem");
	  addr[nt] = &All.MBms;
	  id[nt++] = REAL;

#ifdef LT_AVOID_ENRICH_SFGAS
	  /* this value is set in Msun/yr ;
	   * Gas Particles having a sfr larger that this
	   * value will not receive metals from stars.
	   */
	  strcpy(tag[nt], "SFTh_for_enrich");
	  addr[nt] = &All.Enrich_SFGas_Th;
	  id[nt++] = REAL;
#endif

#ifdef LT_MOD_EFFM
	  strcpy(tag[nt], "ModSEffCrit");
	  addr[nt] = &All.Mod_SEff;
	  id[nt++] = INT;
#endif

#ifdef LT_HOT_EJECTA
	  strcpy(tag[nt], "EgySpecEjecta");
	  addr[nt] = &All.EgySpecEjecta;
	  id[nt++] = REAL;
#endif

	  strcpy(tag[nt], "SnIaDataFileName");
	  addr[nt] = All.SnIaDataFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "SnIaDataSetNum");
	  addr[nt] = &All.Ia_Nset_ofYields;
	  id[nt++] = INT;

	  strcpy(tag[nt], "SnIaEnergy");
	  addr[nt] = &All.SnIaEgy;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "LongLiving_Step_Prec");
	  addr[nt] = &All.LLv_Step_Prec;
	  id[nt++] = REAL;

	  /*
	     strcpy(tag[nt], "SnIa_Remnant");
	     addr[nt] = &All.SnIaRemn;
	     id[nt++] = REAL;
	   */

	  strcpy(tag[nt], "SnIIdataFileName");
	  addr[nt] = All.SnIIDataFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "SnIIDataSetNum");
	  addr[nt] = &All.II_Nset_ofYields;
	  id[nt++] = INT;

	  strcpy(tag[nt], "SnIIEnergy");
	  addr[nt] = &All.SnIIEgy;
	  id[nt++] = REAL;

	  /*
	     strcpy(tag[nt], "ChemTimeStepII");
	     addr[nt] = &All.ChemTimeStepII;
	     id[nt++] = REAL;

	     strcpy(tag[nt], "LongCTStepII");
	     addr[nt] = &All.LongChemTimeStepII;
	     id[nt++] = REAL;
	   */

	  strcpy(tag[nt], "metIRAThMass");
	  addr[nt] = &All.metIRA_ThMass;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "egyIRAThMass");
	  addr[nt] = &All.egyIRA_ThMass;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "SnII_Step_Prec");
	  addr[nt] = &All.SnII_Step_Prec;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "AGBdataFileName");
	  addr[nt] = All.AGBDataFile;
	  id[nt++] = STRING;

	  strcpy(tag[nt], "AGBDataSetNum");
	  addr[nt] = &All.AGB_Nset_ofYields;
	  id[nt++] = INT;

	  strcpy(tag[nt], "MinChemTimeStep");
	  addr[nt] = &All.MinChemTimeStep;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "MinSpreadLength");
	  addr[nt] = &All.MinChemSpreadL;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "Z_toset_SF_DensTh");
	  addr[nt] = &All.referenceZ_toset_SF_DensTh;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "Zdependent_SFTh");
	  addr[nt] = &All.SFTh_Zdep;
	  id[nt++] = INT;

	  strcpy(tag[nt], "MaxSpreadLength");
	  addr[nt] = &All.MaxChemSpreadL;
	  id[nt++] = REAL;

#endif

#ifdef LT_METAL_COOLING_WAL
	  strcpy(tag[nt], "WalCoolTables_Path");
	  addr[nt] = All.WalCool_CoolTables_path;
	  id[nt++] = STRING;
#endif

#ifdef LT_ADD_GAL_TO_SUB
	  strcpy(tag[nt], "BC_SED_Path");
	  addr[nt] = All.BC_SED_Path;
	  id[nt++] = STRING;
#endif

#ifdef GENERATE_GAS_IN_ICS
#ifdef GENERATE_GAS_TG
	  strcpy(tag[nt], "GenGasRefFac");
	  addr[nt] = &All.GenGasRefFac;
	  id[nt++] = INT;
#endif
#endif

#ifdef ADAPTGRAVSOFT
	  strcpy(tag[nt], "AGS_DesNumNgb");
	  addr[nt] = &All.AGS_DesNumNgb;
	  id[nt++] = INT;

	  strcpy(tag[nt], "AGS_MaxNumNgbDeviation");
	  addr[nt] = &All.AGS_MaxNumNgbDeviation;
	  id[nt++] = REAL;
#endif

#ifdef SPHERICAL_BOUNDARY
	  strcpy(tag[nt], "SB_type");
	  addr[nt] = &All.spherical_boundary.type;
	  id[nt++] = INT;

	  strcpy(tag[nt], "SB_radius");
	  addr[nt] = &All.spherical_boundary.radius;
	  id[nt++] = REAL;
#endif

#if GADGET_HYDRO == HYDRO_MFM
	  // default: energy-entropy switch turned off
	  All.EkinSwitchFraction = (double) 0.0;
	  All.EpotSwitchFraction = (double) 0.0;

	  strcpy(tag[nt], "EkinSwitchFraction");
	  addr[nt] = &All.EkinSwitchFraction;
	  id[nt++] = REAL;

	  strcpy(tag[nt], "EpotSwitchFraction");
	  addr[nt] = &All.EpotSwitchFraction;
	  id[nt++] = REAL;
#endif
	}			// if include_tags == NULL
      /*
         if we are including a paremter file , we need the parent tag array because
         the code set the tag[i]==0 in order to check if the tag was assigned or not.
       */


      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      printf("Obtaining parameters from file '%s':\n", fname);
	      while(!feof(fd))
		{

		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)

		    if(strcmp("Include", buf1) == 0)
		      {
			char include_file_name[MAX_PATH_LENGTH];
#ifdef _WIN32
			strcpy(include_file_name, buf2);
#else
			sprintf(include_file_name, "%s/%s", dirname(strdup(fname)), buf2);
#endif

			PANIC_IF(strlen(include_file_name) >= MAX_PATH_LENGTH,
				 "Inclusion path exeeded limit of %d chars: %s\n", MAX_PATH_LENGTH,
				 include_file_name);

			read_parameter_file(include_file_name, tag, addr, id, nt);
			j = -2;	//j=-2 signal that we didn't fund any parameter because we are including stuff
			break;
		      }
		    else if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  fprintf(stdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  fprintf(stdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  fprintf(stdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else if(j == -2)
		    {
		      //we just included a parameter file.
		    }
		  else
		    {

		      fprintf(stdout,
			      "begun: WARNING: Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 2;
		    }
		}
	      fclose(fd);
	      fclose(fdout);
	      printf("\n");

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
	      int ret;

	      ret = system(buf3);
#endif
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}


      if(included)
	{			//we are an included parameter file, our job ends here!
	  return;
	}

#ifdef LMB_SPECTRAL_CRs
      printf("LMB_SPECTRAL_CRs model switched on\n");
#ifdef LMB_CR_PROTONS
      printf("Proton bins = %i\n", LMB_CR_PROTONS);
#else
      if(All.CR_ProtonsOn == 1)
	{
	  printf
	    ("Error. Protons are switched on in parameter file, but code was compiled without LMB_CR_PROTONS.\n");
	  errorFlag = 1;
	}
#endif // LMB_CR_PROTONS
#ifdef LMB_CR_ELECTRONS
      printf("Electron bins = %i\n", LMB_CR_ELECTRONS);
#else
      if(All.CR_ElectronsOn == 1)
	{
	  printf
	    ("Error. Electrons are switched on in parameter file, but code was compiled without LMB_CR_ELECTRONS.\n");
	  errorFlag = 1;
	}
#endif // LMB_CR_ELECTRONS
#ifdef JD_SHOCK
      if(All.CR_DSA_Model == 0)
	printf("DSA model by Kang et al. 2007: ApJ, 669, 729\n");
      if(All.CR_DSA_Model == 1)
	printf("DSA model by Kang&Ryu 2013: ApJ, 764, 95\n");
      if(All.CR_DSA_Model == 2)
	printf("DSA model by Ryu et al. 2019: https://arxiv.org/pdf/1905.04476v2.pdf\n");
      if(All.CR_DSA_Model == 3)
	printf("DSA model by Caprioli&Spitkovsky 2014: ApJ, 783, 91\n");
      if(All.CR_DSA_Model == 4)
	printf("Constant efficiency like in Pfrommer+2016: MNRAS, 465, 4500-4529\n");
#endif // JD_SHOCK

#endif // LMB_SPECTRAL_CRs

      for(i = 0; i < nt; i++)
	{
#ifdef LT_STELLAREVOLUTION
	  if(*tag[i] &&
	     (strcmp(tag[i], "metIRAThMass") != 0 && strcmp(tag[i], "egyIRAThMass") != 0) &&
	     strcmp(tag[i], "WindEfficiency") != 0 && strcmp(tag[i], "WindEnergyFraction") != 0 &&
	     strcmp(tag[i], "LocalSpreadFactor") != 0 && strcmp(tag[i], "TestSuite"))
#else
	  if(*tag[i])
#endif
#if GADGET_HYDRO == HYDRO_MFM
	    if(strcmp(tag[i], "EkinSwitchFraction") != 0 && strcmp(tag[i], "EpotSwitchFraction") != 0)
#endif
	      {
		printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
		errorFlag = 1;
	      }
	}
#ifndef GADGET3_IO_LIB
      if(All.OutputListOn)
	{
	  int res_outputlist = read_outputlist(All.OutputListFilename);
	  PANIC_IF(res_outputlist != 0, "Unable to process '%s'", All.OutputListFilename);
	}
      else
#endif
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MYMPI_COMM_WORLD);


  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  PANIC_IF(errorFlag == 1
	   || (errorFlag == 2
	       && All.LevelOfStrickness == 0),
	   " Problems in parameter file detected, set LevelOfStrickness >=1 to continue in case of warnings");



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched on.\n");
	  printf("You must set `CoolingOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.CoolingOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched off.\n");
	  printf("You must set `CoolingOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#ifdef SFR

#ifndef MOREPARAMS
  if(ThisTask == 0)
    {
      printf("Code was compiled with SFR, but not with MOREPARAMS.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif

  if(All.StarformationOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched on.\n");
	  printf("You must set `StarformationOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("You try to use the code with star formation enabled,\n");
	  printf("but you did not switch on cooling.\nThis mode is not supported.\n");
	}
      endrun(0);
    }
#else
  if(All.StarformationOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched off.\n");
	  printf("You must set `StarformationOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif



#ifdef METALS
#ifndef SFR
  if(ThisTask == 0)
    {
      printf("Code was compiled with METALS, but not with SFR.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef MOREPARAMS
#ifdef TIME_DEP_ART_VISC
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIME_DEP_ART_VISC, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DARKENERGY, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef TIMEDEPDE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif


#ifndef MAGNETIC
#ifdef TRACEDIVB
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TRACEDIVB, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DBOUTPUT
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DBOUTPUT, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MAGFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGFORCE, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef BSMOOTH
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with BSMOOTH, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MU0_UNITY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MU0_UNITY, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef MAGNETIC_DISSIPATION
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGNETIC_DISSIPATION, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef TIME_DEP_MAGN_DISP
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIME_DEP_MAGN_DISP, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef DIVBCLEANING_DEDNER
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBCLEANING_DEDNER, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#endif

#ifndef MAGFORCE
#if defined(DIVBFORCE) || defined(DIVBFORCE3)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBFORCE, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#if defined(DIVBFORCE) && defined(DIVBFORCE3)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBFORCE and with DIVBFORCE3.\n");
      fprintf(stdout, "This will lead to no correction at all, better stop.\n");
    }
  endrun(0);
#endif

#if defined(LT_STELLAREVOLUTION)

#if !defined(LT_PM_LIFETIMES) && !defined(LT_MM_LIFETIMES)
  if(ThisTask == 0)
    printf("LT_STELLAREVOLUTION requires either LT_PM_LIFETIMES or LT_MM_LIFETIMES\n");
  endrun(0);
#endif

#if (defined(LT_WIND_VELOCITY) || defined(LT_HYDROWINDS) || defined(LT_DECOUPLE_POSTWINDS_FROM_SF)) && !defined(WINDS)
  if(ThisTask == 0)
    printf("LT_WIND_VELOCITY, LT_HYDROWINDS and LT_DECOUPLE_POSTWINDS_FROM_SF require WIND\n");
  endrun(0);
#endif


#endif

#if defined(LT_METAL_COOLING)

#if defined(LT_METAL_COOLING_on_SMOOTH_Z) && !defined(LT_SMOOTH_Z)
  if(ThisTask == 0)
    printf("LT_METAL_COOLING_on_SMOOTH_Z requires LT_SMOOTH_Z\n ");
  endrun(0);
#endif

#if defined(LT_METAL_COOLING_on_SMOOTH_Z) && defined(LT_METAL_COOLING_WAL) && !defined(LT_SMOOTH_ALLMETALS)
  if(ThisTask == 0)
    printf("LT_METAL_COOLING_on_SMOOTH_Z requires LT_SMOOTH_ALLMETALS if LT_METAL_COOLING_WAL is defined\n ");
  endrun(0);
#endif


#if defined(CONDUCTION_INCLUDEMAGNETIC) && (!defined(CONDUCTION) || !defined(MAGNETIC))
  if(ThisTask == 0)
    printf("CONDUCTION_INCLUDEMAGNETIC requires CONDUCTION and MAGNETIC\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_INCLUDEMAGNETIC_EASY) && !defined(CONDUCTION_INCLUDEMAGNETIC)
  if(ThisTask == 0)
    printf("CONDUCTION_INCLUDEMAGNETIC_EASY requires CONDUCTION_INCLUDEMAGNETIC\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_INCLUDEMAGNETIC_PERPENDICULAR) && !defined(CONDUCTION_INCLUDEMAGNETIC)
  if(ThisTask == 0)
    printf("CONDUCTION_INCLUDEMAGNETIC_PERPENDICULAR requires CONDUCTION_INCLUDEMAGNETIC\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_INCLUDEMAGNETIC_ISOTROPISED) && !defined(CONDUCTION_INCLUDEMAGNETIC)
  if(ThisTask == 0)
    printf("CONDUCTION_INCLUDEMAGNETIC_ISOTROPISED requires CONDUCTION_INCLUDEMAGNETIC\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_BICGSTAB) && !defined(CONDUCTION_INCLUDEMAGNETIC)
  if(ThisTask == 0)
    printf("CONDUCTION_BICGSTAB requires CONDUCTION_INCLUDEMAGNETIC\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_INCLUDEMAGNETIC_CHECK_ANGLE) && !defined(CONDUCTION_INCLUDEMAGNETIC)
  if(ThisTask == 0)
    printf("CONDUCTION_INCLUDEMAGNETIC_CHECK_ANGLE requires CONDUCTION_INCLUDEMAGNETIC\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_INCLUDEMAGNETIC_PERPENDICULAR) && defined(CONDUCTION_INCLUDEMAGNETIC_EASY)
  if(ThisTask == 0)
    printf
      ("CONDUCTION_INCLUDEMAGNETIC_PERPENDICULAR can not be used with CONDUCTION_INCLUDEMAGNETIC_EASY\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_INCLUDEMAGNETIC_ISOTROPISED) && defined(CONDUCTION_INCLUDEMAGNETIC_EASY)
  if(ThisTask == 0)
    printf("CONDUCTION_INCLUDEMAGNETIC_ISOTROPISED can not be used with CONDUCTION_INCLUDEMAGNETIC_EASY\n ");
  endrun(0);
#endif

#if defined(CONDUCTION_INCLUDEMAGNETIC_CHECK_ANGLE) && defined(CONDUCTION_INCLUDEMAGNETIC_EASY)
  if(ThisTask == 0)
    printf("CONDUCTION_INCLUDEMAGNETIC_CHECK_ANGLE can not be used with CONDUCTION_INCLUDEMAGNETIC_EASY\n ");
  endrun(0);
#endif

#endif

/* ---------------------------------------------------- */
/* --------- Place here obsolete switches !!!! -------- */
/* ---------------------------------------------------- */

/*
#ifdef SUBFIND_ID_OVERBOOK_FACTOR
#error "SUBFIND_ID_OVERBOOK_FACTOR switch obsolete, remove it from Config.sh !"
#endif
*/

#ifdef SUBFIND_GROUP_SUBCOMM
#error "SUBFIND_GROUP_SUBCOMM switch obsolete, remove it from Config.sh !"
#endif

#ifdef KD_CHOOSE_PSUBFIND_LIMIT
#error "KD_CHOOSE_PSUBFIND_LIMIT switch obsolete, remove it from Config.sh !"
#endif

#ifdef SAVE_MASS_TAB
#error "SAVE_MASS_TAB switch obsolete, remove it from Config.sh !"
#endif

#ifdef FFTW3
#error "FTTW3 switch obsolete, remove it from Config.sh !"
#endif

#ifdef PM_FROM_L3
#error "PM_FROM_L3 switch obsolete, remove it from Config.sh !"
#endif

#ifdef PM_FROM_L3_COMBINED
#error "PM_FROM_L3_COMBINED switch obsolete, remove it from Config.sh !"
#endif

#ifdef FFTW3_OMP
#error "FFTW3_OMP switch obsolete, remove it from Config.sh !"
#endif

#ifdef NOTYPEPREFIX_FFTW
#error "NOTYPEPREFIX_FFTW switch obsolete, remove it from Config.sh !"
#endif

#ifdef OPENMP
#error "OPENMP switch obsolete, remove it from Config.sh !"
#endif

#ifdef KD_HMAX_ESTIMATE
#error "KD_HMAX_ESTIMATE switch obsolete, remove it from Config.sh !"
#endif

#ifdef AR_GREEN_TREE_HYDRA
#error "AR_GREEN_TREE_HYDRA switch obsolete, remove it from Config.sh !"
#endif

#ifdef AR_GREEN_TREE_DENSITY
#error "AR_GREEN_TREE_DENSITY switch obsolete, remove it from Config.sh !"
#endif

#ifdef AR_GREEN_TREE_CONDUCTION
#error "AR_GREEN_TREE_CONDUCTION switch obsolete, remove it from Config.sh !"
#endif

#ifdef AR_NODE_REORDERED
#error "AR_NODE_REORDERED switch obsolete, remove it from Config.sh !"
#endif

#ifdef KD_NEW_ORDER_IN_PARTICLE_STRUCTURE
#error "KD_NEW_ORDER_IN_PARTICLE_STRUCTURE switch obsolete, remove it from Config.sh !"
#endif


  if(ThisTask == 0)
    {
      printf("BegRun: Check BlackHoleSetting ...\n");
      check_and_report_settings_for_blackhole();
#if defined(fSIDM) || defined(rSIDM)
      check_and_report_settings_for_msidm();
#endif
    }

  free(addr);
  free(id);
  for(int i = 0; i < MAXTAGS; i++)
    {
      free(tag[i]);
    }
  free(tag);


}

#undef REAL
#undef STRING
#undef INT
#undef MAXTAGS

/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;
  int count, flag;
  char buf[512];

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;

  while(1)
    {
      if(fgets(buf, 500, fd) != buf)
	break;

      count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);

      if(count == 1)
	flag = 1;

      if(count == 1 || count == 2)
	{
	  if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
	    {
	      if(ThisTask == 0)
		printf("\ntoo many entries in output-list. You should increase MAXLEN_OUTPUTLIST=%d.\n",
		       (int) MAXLEN_OUTPUTLIST);
	      endrun(13);
	    }

	  All.OutputListFlag[All.OutputListLength] = flag;
	  All.OutputListLength++;
	}
    }

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = (long long) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);
  else
    ti_end = (long long) ((TimeMax_new - All.TimeBegin) / All.Timebase_interval);

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif
#ifdef CONDUCTION
      All.Conduction_Ti_begstep /= 2;
      All.Conduction_Ti_endstep /= 2;
#endif
#ifdef LMB_SPECTRAL_CRs_DIFFUSION
      All.CRDiffusion_Ti_begstep /= 2;
      All.CRDiffusion_Ti_endstep /= 2;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_current /= 2;

	  if(P[i].TimeBin > 0)
	    {
	      P[i].TimeBin--;
	      if(P[i].TimeBin <= 0)
		{
		  printf("Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
		  endrun(8765);
		}
	    }
	}

      All.Ti_nextlineofsight /= 2;
    }

  All.TimeMax = TimeMax_new;
}
