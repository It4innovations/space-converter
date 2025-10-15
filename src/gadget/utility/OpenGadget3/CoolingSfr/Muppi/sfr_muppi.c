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
#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"
#include "../../Gravity/forcetree.h"
#include "muppi_proto.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/* GM WARNING: all the variables defined as <<double>> in this file,
are still used as <<double>>, and not as <<MyAtLeastDouble>>, somewhere else.
Note that the structure All contains <<double>>s. The same holds for function 
calls (formal parameters declared or casted <<double>> are required by
the function definitions) */

#if defined(MV_GM_AGNMUPPI) && !(defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK))
#error MV_GM_AGNMUPPI requires BH_THERMALFEEDBACK or BH_KINETICFEEDBACK
#endif


/* VERY IMPORTANT: COMPILE THIS CODE IN DOUBLE PRECISION OR IT WILL NOT WORK */
#if defined(GM_MUPPI) && !defined(DOUBLEPRECISION)	/* XXXX ctrl for OpenGadget */
#error GM_MUPPI *REQUIRES* DOUBLEPRECISION
#endif

#ifdef COOLING
#ifdef SFR
#ifdef GM_MUPPI

#define checkline(x) check_line(x, __FUNCTION__, __FILE__, __LINE__)
void check_line(const char *text, const char *func, const char *file, int line);


/* if incompatible directives are chosen, the compilation stops with an error */
//#if defined(SFR_METALS) || defined(MHM) || defined(BG_SFR) || defined(COSMIC_RAYS) || defined(QUICK_LYALPHA) || defined(BLACK_HOLES) || defined(SFR_FEEDBACK) // XXXX ctrl for OpenGadget
//#error PLEASE CHECK PRECOMPILER DIRECTIVES WHEN MUPPI IS USED
//#endif

/* ALL muppi parameters inside the following header. */
#include "muppi.h"

/* GLOBAL functions variables for this module */
/* note. Pay attention to these, when rewriting for OpenMP */


MyAtLeastDouble ESN_tot;
MyAtLeastDouble tempthresh = 0.0;
MyAtLeastDouble E_SN_51, beta_sf, f_rest;
MyAtLeastDouble dt, dtime, a3inv, Time;	/* XXXX put these in auto? */
int stars_spawned, stars_converted, bits;

int nmf;
int nexit_dens, nexit_clock, nexit_mnc, nexit_freeze, nexit_depl, nexit_spawn, nexit_gsl, nexit_sfraw,
  nexit_tms, nexit_force, nexit_wind;

MyAtLeastDouble Sum_sm, Sum_mass_stars;

size_t neqz;			/* this is 4 without AGNMUPPI and 5 with it; initialized in MUPPI_main */

/* AGN-related counters */
#ifdef MV_GM_AGNMUPPI
#if defined (MV_GM_AGNMUPPI_OUTPUT)
double E_AGNM_tot, E_AGNM_cool, E_AGNM_sfr_h, E_AGNM_sfr, E_AGNM_sfr_c;
double E_AGNM_sfr_dE, E_AGNM_sfr_dE_aft, E_AGNM_sfr_BHEcoldD;
int n_AGNM, n_AGNcool;
#endif
#endif


/* note: here I use <<double>> because I want to use MPI_DOUBLE reductions */
#ifdef COUNT_PARTICLES_IN_CONE
int histog[21], tothistog[21], totnactive;
double Energy_total, Energy_used, Energy_total_1, Energy_used_1;
double S_Energy_total, S_Energy_used, S_Energy_total_1, S_Energy_used_1;
#endif

#ifdef LT_STELLAREVOLUTION
int SFi;
#ifdef LT_METAL_COOLING
double Zcool;
#endif
#ifdef LT_METAL_COOLING_WAL	/* XXXX warning, these are currently defined <<double>>
				   in lt_sft */
static double *mMetallicities, mDZ, mRedshift;
#if defined(GL_DUST_COOLING)
static double DL, DS;		// particle's dust to gas ratio for small and big grains
#endif // GL_DUST_COOLING
#endif
#endif


/* Table of exit codes from multi-phase

   MultiPhase =   0   low density exit         -  normal behaviour
   MultiPhase =  -1   clock exit               -  normal behaviour
   MultiPhase =  -2   force single phase exit  -  special behaviour
   MultiPhase =  -3   freezing exit            -  special behaviour
   MultiPhase =  -4   depletion of cold phase  -  traumatic case
   MultiPhase =  -5   GSL error                -  traumatic case
   MultiPhase =  -6   mass is not conserver    -  traumatic case  REMOVED
   MultiPhase =  -7   wind exit                -  special behaviour
   MultiPhase =  -8   too many steps           -  traumatic case
   MultiPhase = -10   star formation runaway   -  traumatic case

   MultiPhase -= 100  a hot particle in the case of phase decoupling

*/


/*
 * This routine performs standard cooling for single-phase particles 
 * and implements the MUlti-Phase Particle Integrator (MUPPI) model
 */

void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  MyAtLeastDouble ascale = 1, hubble_a = 0;
  MyAtLeastDouble time_hubble_a;
  int mfin, mfout;

#ifndef  LT_METAL_COOLING_WAL
  MyAtLeastDouble ne_in;
#endif

#ifdef LT_STELLAREVOLUTION
  int IMFi, GENERATIONS_BUNCH;
#endif


#ifdef MV_GM_AGNMUPPI
#if defined (MV_GM_AGNMUPPI_OUTPUT)
  E_AGNM_tot = E_AGNM_cool = E_AGNM_sfr_h = E_AGNM_sfr = E_AGNM_sfr_c = 0.0;
  E_AGNM_sfr_dE = E_AGNM_sfr_dE_aft = E_AGNM_sfr_BHEcoldD = 0.0;
  n_AGNM = n_AGNcool = 0;
#endif
#endif


  /* initializing global counters */
  stars_spawned = stars_converted = 0;
  Sum_sm = Sum_mass_stars = 0;

  nmf = 0;
  nexit_dens = nexit_clock = nexit_freeze = nexit_depl = nexit_spawn = nexit_gsl = nexit_sfraw = nexit_mnc =
    nexit_tms = nexit_force = nexit_wind = 0;


#ifdef COUNT_PARTICLES_IN_CONE
  for(int i = 0; i < 21; i++)
    histog[i] = 0;
  Energy_total = Energy_used = Energy_total_1 = Energy_used_1 = S_Energy_total = S_Energy_used =
    S_Energy_total_1 = S_Energy_used_1 = 0.0;
#endif

  for(int bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinActive[bin])
      TimeBinSfr[bin] = 0;

  if(All.ComovingIntegrationOn)
    {

      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY		/* XXXX check with OpenGadget */
	+ DarkEnergy_a(All.Time);
#else
	+ All.OmegaLambda;
#endif /* DARKENERGY */

      hubble_a = All.Hubble * sqrt(hubble_a);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;

    }
  else
    a3inv = ascale = time_hubble_a = 1;


  /*  compute the number of bits necessary to store number of stars generated */
#ifndef LT_STELLAREVOLUTION
  for(bits = 0; GENERATIONS > (1 << bits); bits++);
#else
  bits = All.StarBits;
#endif


#ifdef LT_METAL_COOLING_WAL	/* initializing cooling at this redshift */
  mRedshift = 1.0 / ascale - 1;
  get_cool_redshift(mRedshift, &mDZ);
  mMetallicities = (double *) mymalloc("Metallicities", sizeof(double) * get_cool_n_el());
  WalCool_tables_load(mRedshift);
#endif

  mfin = mfout = 0;

  /* setting SPH temperatures for all gas particles (to be used later) */
  for(int i = 0; i < N_gas; i++)
    {
      double sum_old_new_energies, temp = 0.0;

      sum_old_new_energies = SphP[i].EntropyPred / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1) * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;	/* this is independent of h */

#ifndef LT_METAL_COOLING_WAL
      rho = SphP[i].Density * a3inv * All.UnitDensity_in_cgs;	/* we transform this in h true */
      rho *= All.HubbleParam * All.HubbleParam;
      ne_in = SphP[i].elec;
      temp = convert_u_to_tempSTD(sum_old_new_energies, rho, &ne_in);	/* convert_u_to_temp requires h-true cgs units */
#else
      temp = convert_u_to_tempMET(sum_old_new_energies, mRedshift, mDZ, mMetallicities);
#endif
      SphP[i].Temperature = temp;	/* XXXX warning duplicated Temp/Temperature */

    }


  for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      /* begin of loop over gas particles */
      if(P[i].Type == 0)
	{
	  int flag;
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;

	  if(All.ComovingIntegrationOn)
	    {
	      dtime = All.Time * dt / time_hubble_a;
	      Time = All.Time / time_hubble_a;
	    }
	  else
	    {
	      dtime = dt;
	      Time = All.Time;
	    }

	  /* initializing IMF related quantities  (gas particle could have different IMFs */
#ifdef LT_STELLAREVOLUTION
	  get_SF_index(i, &SFi, &IMFi);

	  E_SN_51 = (MyAtLeastDouble) All.E_SN_51[SFi];
	  beta_sf = (MyAtLeastDouble) All.beta_sf[SFi];
	  f_rest = (MyAtLeastDouble) 0.0;	/* in this case, restoration is done by evolve_sn */
	  GENERATIONS_BUNCH = All.Generations / SFs[SFi].Generations;
	  SphP[i].mstar = 0.0;
#else
	  E_SN_51 = myE_SN_51;
	  beta_sf = mybeta_sf;
	  f_rest = myf_rest;
#endif

#ifdef LT_METAL_COOLING_WAL	/* initializing cooling of THIS gas particle */
	  set_metallicities(i, mMetallicities, (double) a3inv);	/* note: sets the physical density too! */
#if defined(GL_DUST_COOLING)	// compute DSG for cooling and/or scaling SF; also Z for scaling SF
	  if(P[i].Mass > 0)
	    {
	      DL = get_metalmass(SphP[i].DustL) / P[i].Mass;
	      DS = get_metalmass(SphP[i].DustS) / P[i].Mass;
	    }
	  else
	    {
	      DL = 0;
	      DS = 0;
	    }
#endif //GL_DUST_COOLING
#elif defined(LT_METAL_COOLING)
	  Zcool = get_metallicity_solarunits(get_metallicity(i, Iron));
#endif


	  /* check whether conditions for star formation are fulfilled.
	   *  
	   * flag=1  normal cooling
	   * flag=0  star formation
	   */
	  flag = check_sf_or_cooling(i);	/* default is flag=1, normal cooling */


#ifdef MV_GM_STELLAR_KIN_FB2
	  SphP[i].DelayTime0 = 0.0;
#endif

#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
	  if(SphP[i].Injected_BH_Energy > 0. && SphP[i].MultiPhase > 0)
	    {
	      MyAtLeastDouble Vol_tot, T_h_m, fill_h_m, Rho_h_m, n_h_m, Rho_c_m, n_c_m;
	      Vol_tot = (P[i].Mass - SphP[i].M_sf) / (SphP[i].Density * a3inv) / pow(All.HubbleParam, 3);
	      compute_ism_properties(SphP[i].M_h / All.HubbleParam,
				     SphP[i].E_h / All.HubbleParam,
				     SphP[i].M_c / All.HubbleParam,
				     Vol_tot, &T_h_m, &fill_h_m, &Rho_h_m, &Rho_c_m, &n_h_m, &n_c_m);
	      if(fill_h_m < 1.0)
		{
		  SphP[i].CouplingBHEnCold = ((1.0 - fill_h_m) * (pow(3. * Vol_tot / 4. / M_PI, 1. / 3.) / MV_GM_COVER_FACT_MCLS));	/* covering factor = filling factor * (hsml della p.lla MP / Giant Molecular Cloud Linear Size) , MV_GM_COVER_FACT_MCLS=parametro del modello, milena  MV+GM change particle volume for covering factor */
		}
	      else
		{
		  SphP[i].CouplingBHEnCold = 0.0;
		}
	      SphP[i].CouplingBHEnHot = 1.0 - SphP[i].CouplingBHEnCold;

	      if(SphP[i].CouplingBHEnCold > 1.0)
		{
		  SphP[i].CouplingBHEnCold = 1.0;
		  SphP[i].CouplingBHEnHot = 1.0 - SphP[i].CouplingBHEnCold;
		}
	      fprintf(FdCovFacAGN,
		      "%g    %d   %e  %e     %e        %e    %e       %e  %e  %e  %e  %e  %e    %e\n",
		      All.Time, P[i].ID, fill_h_m, pow(3. * Vol_tot / 4. / M_PI, 1. / 3.),
		      SphP[i].Injected_BH_Energy, SphP[i].CouplingBHEnHot, SphP[i].CouplingBHEnCold,
		      SphP[i].M_h, SphP[i].M_c, T_h_m, n_h_m, n_c_m, SphP[i].Sfr, P[i].Hsml);
	      fflush(FdCovFacAGN);	//milena


	    }
#endif


	  /* protection against start-up dt=0 */
	  if(dt <= 0.0 || !P[i].TimeBin)
	    flag = -1;


	  if(flag == 1)		/* normal implicit isochoric cooling */
	    cooling(i, dtime, flag);

	  if(flag == 0)		/* active star formation */
	    starformation(i, GENERATIONS_BUNCH, &mfin, &mfout, dtime, flag);

	}

    }				/* end of main loop over active particles */

  final_statistics(time_hubble_a, mfin, mfout);	/* producing sfr data for sfr file, particles in cones, etc */

#ifdef LT_METAL_COOLING_WAL
  myfree(mMetallicities);
#endif
}


/* This function performs a check over gas particles, to estabilish if it must go into starformation, 
   or only do isocoric cooling */
static inline int check_sf_or_cooling(int i)
{
  int flag = 1;


  if(SphP[i].MultiPhase > 0)
    flag = 0;

#if defined(MV_GM_STELLAR_KIN_FB2)
  /* particles in fountain do not go multi-phase */
  else if(SphP[i].DelayTime > 0.0)
    {
      flag = 1;
#ifdef LT_STELLAREVOLUTION
      SphP[i].DelayTime0 = SphP[i].DelayTime;
#endif
      SphP[i].DelayTime -= dtime / All.HubbleParam;	/* delay time is in true h */

      /* density threshold for fountain */
      if(SphP[i].Density * a3inv * All.HubbleParam * All.HubbleParam <
	 FOUNTAIN_DENSITY_THR * All.PhysDensThresh)
	SphP[i].DelayTime = 0;
    }
#endif
  else if(SphP[i].MultiPhase == 0)
    {
      if(SphP[i].Density * a3inv * All.HubbleParam * All.HubbleParam >= All.PhysDensThresh
	 && SphP[i].Temperature <= TEMPSFTHRES)
	{
	  flag = 0;
	  SphP[i].NoArtCond = 0;
	}
      if(SphP[i].NoArtCond == 1)
	{
	  if((SphP[i].Temperature <= 1.2 * SphP[i].AveT) && (SphP[i].Temperature >= 0.8 * SphP[i].AveT))
	    SphP[i].NoArtCond = 0;
	}

    }

  if(All.ComovingIntegrationOn && SphP[i].Density < All.OverDensThresh)	/* overdensity threshold */
    flag = 1;

#ifdef MV_GM_AGNEGY_NO_MP
  if(SphP[i].Injected_BH_Energy > 0.)
    {

      SphP[i].M_sf = 0.0;
      SphP[i].MultiPhase = 0;
      SphP[i].M_h = 0.0;
      SphP[i].M_c = 0.0;
      SphP[i].E_h = 0.0;
      SphP[i].E_out = 0.0;
#ifdef LT_STELLAREVOLUTION
      SphP[i].mstar = 0.0;
#endif
      flag = 1;			//particle is never multiphase if it takes BH energy
    }
#endif


  return flag;

}



/* This function performs isocoric cooling */
/* Called vy cooling_and_starformation() */
static inline void cooling(int i, MyAtLeastDouble dtime, int flag)
{
  MyAtLeastDouble ne, MyEgy;
  double unew, temp, sum_old_new_energies, rho;

#ifndef LT_METAL_COOLING_WAL
  double ne_in;
#endif


  SphP[i].Sfr = 0;

  ne = SphP[i].elec;		/* electron abundance (gives ionization state and mean molecular weight) */


  unew = DMAX(All.MinEgySpec,	/* min temperature allowed expressed as energy per unit mass */
	      (SphP[i].Entropy) /	/* rate of change of entropy modified at line 281 */
	      GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1));	/* new thermal energy */

#if !defined(LT_METAL_COOLING) && !defined(LT_METAL_COOLING_WAL)
  unew = DoCoolingSTD(unew, SphP[i].Density * a3inv, dtime, &ne);
#else
#if defined(LT_METAL_COOLING)
  unew = DoCoolingZ(unew, (double) SphP[i].Density * a3inv, dtime, &ne, Zcool, &temp);
#endif
#if defined(LT_METAL_COOLING_WAL)
#ifdef GL_DUST_COOLING
  unew =
    DoCoolingMET(unew, (double) SphP[i].Density * a3inv, mMetallicities, DL, DS, mRedshift, mDZ, dtime,
		 &temp);
#else
  unew = DoCoolingMET(unew, (double) SphP[i].Density * a3inv, mMetallicities, mRedshift, mDZ, dtime, &temp);
#endif // GL_DUST_COOLING
#endif
#endif

#ifdef LT_STELLAREVOLUTION
  /* Here is where energy from SNIa is given */
  MyEgy = (dtime > SphP[i].EgyStep ? SphP[i].EgyRes : SphP[i].EgyRes * dtime / SphP[i].EgyStep);
  unew += MyEgy / P[i].Mass;
  SphP[i].EgyRes -= MyEgy;
#endif

#if defined(MV_GM_AGNMUPPI) || defined(MV_GM_AGNEGY_NO_MP)
  if(SphP[i].Injected_BH_Energy > 0.)
    {
      if(P[i].Mass == 0)
	SphP[i].Injected_BH_Energy = 0;
      else
	{
	  unew += SphP[i].Injected_BH_Energy / P[i].Mass;

#if defined (MV_GM_AGNMUPPI_OUTPUT)
	  E_AGNM_cool += SphP[i].Injected_BH_Energy;	//milena
	  E_AGNM_tot += SphP[i].Injected_BH_Energy;	//milena
	  n_AGNcool++;


	  // KLAUS: is this print inside threads zone OK?
	  fprintf(FdEgyAGN, "%g   %f %f %f   %f %f %f   %e   %e %e %e   %d   %e %e %e   %e    %e    %e\n",
		  All.Time, P[i].Pos[0], P[i].Pos[0], P[i].Pos[2],
		  P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass,
		  SphP[i].Density, SphP[i].Pressure, SphP[i].Temperature, P[i].ID,
		  flag, SphP[i].Injected_BH_Energy, SphP[i].MultiPhase, SphP[i].DelayTime, 1.0, 0.0);
	  fflush(FdEgyAGN);	//milena
#endif
	}

#ifdef LT_MAX_TEMP_FEEDBACK
      temp = u_to_temp_fac * unew;
      if(temp > LT_MAX_TEMP_FEEDBACK)
	unew = LT_MAX_TEMP_FEEDBACK / u_to_temp_fac;
#endif
    }
#endif


  SphP[i].elec = ne;

  SphP[i].Entropy = unew / pow((double) SphP[i].Density * a3inv, GAMMA_MINUS1) * GAMMA_MINUS1;


  SphP[i].EntropyPred = SphP[i].Entropy;
  SphP[i].Pressure = get_pressure(i);

  sum_old_new_energies = (SphP[i].Entropy) / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1) * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;	/* this is independent of h */
  rho = SphP[i].Density * a3inv * All.UnitDensity_in_cgs;	/* we transform this in h true */
  rho *= All.HubbleParam * All.HubbleParam;

#ifndef LT_METAL_COOLING_WAL
  ne_in = SphP[i].elec;
  temp = convert_u_to_tempSTD(sum_old_new_energies, rho, &ne_in);	/* convert_u_to_temp requires h-true cgs units */
#else
  temp = convert_u_to_tempMET(sum_old_new_energies, mRedshift, mDZ, mMetallicities);
#endif
  SphP[i].Temperature = temp;	/* WARNING dupicate Temp/Temperature */


}





/* This function performs star formation with MUPPI */
/* Called by cooling_and_starformation() */
static inline void starformation(int i, int GENERATIONS_BUNCH, int *mfin, int *mfout, MyAtLeastDouble dtime,
				 int flag)
{

  MyAtLeastDouble sm;
  MyAtLeastDouble tstart, tend;



  sm = SphP[i].M_sf;		/* stellar mass before calling MUPPI */

  /****************** multiphase particle integration ******************/
  tstart = second();

  MUPPI_main(i, SphP[i].Density * a3inv, mfin, mfout);


  tend = second();
  CPU_Step[CPU_MUPPI] += timediff(tstart, tend);
  /* this should exactly measure the integration time */

  /* various counters */
  if(SphP[i].MultiPhase > 0)
    nmf++;

  /* exit status */
  /* previously:
     if(SphP[i].MultiPhase==  0) nexit_dens++;
     if(SphP[i].MultiPhase== -1) nexit_clock++;
     if(SphP[i].MultiPhase== -2) nexit_force++;
     if(SphP[i].MultiPhase== -3) nexit_freeze++;
     if(SphP[i].MultiPhase== -4) nexit_depl++;
     if(SphP[i].MultiPhase== -5) nexit_gsl++;
     if(SphP[i].MultiPhase== -6) nexit_mnc++;
     if(SphP[i].MultiPhase== -7) nexit_wind++;
     if(SphP[i].MultiPhase== -8) nexit_tms++;
     if(SphP[i].MultiPhase==-10) nexit_sfraw++;
   */

  if(SphP[i].MultiPhase == 0)
    nexit_dens++;


  switch (-SphP[i].MultiPhase)
    {
    case 1:
      nexit_clock++;
      break;
    case 2:
      nexit_clock++;
      break;
    case 3:
      nexit_force++;
      break;
    case 4:
      nexit_freeze++;
      break;
    case 5:
      nexit_gsl++;
      break;
    case 6:
      nexit_mnc++;
      break;
    case 7:
      nexit_wind++;
      break;
    case 8:
      nexit_tms++;
      break;
    case 10:
      nexit_sfraw++;
      break;
    }



  if(SphP[i].MultiPhase <= 0)
    SphP[i].NoArtCond = 1;


#ifdef GM_COUNT_PARTICLES_IN_CONE
  Energy_total += SphP[i].E_out * All.FracEgyKin;
  Energy_used += SphP[i].Ekin_total;
  if(SphP[i].nnorm > 0)
    {
      Energy_total_1 += SphP[i].E_out * All.FracEgyKin;
      Energy_used_1 += SphP[i].Ekin_total;
    }
#endif

#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)

#if  defined(MV_GM_AGNMUPPI_OUTPUT)

  E_AGNM_sfr_h += SphP[i].Injected_BH_Energy * SphP[i].CouplingBHEnHot;	//milena
  E_AGNM_sfr_c += SphP[i].Injected_BH_Energy * SphP[i].CouplingBHEnCold;	//milena
  E_AGNM_tot += SphP[i].Injected_BH_Energy * SphP[i].CouplingBHEnHot +
    fprintf(FdEgyAGN, "%g   %f %f %f   %f %f %f   %e   %e %e %e   %d   %e %e %e   %e     %e   %e\n",
	    All.Time, P[i].Pos[0], P[i].Pos[0], P[i].Pos[2],
	    P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass,
	    SphP[i].Density, SphP[i].Pressure, SphP[i].Temperature, P[i].ID,
	    flag, SphP[i].Injected_BH_Energy, SphP[i].MultiPhase,
	    SphP[i].DelayTime, SphP[i].CouplingBHEnHot, SphP[i].CouplingBHEnCold);
  fflush(FdEgyAGN);		//milena

#endif

  SphP[i].Injected_BH_Energy = 0.;

#endif



#if defined(MV_GM_STELLAR_KIN_FB2)
  /* selects particles that must be kicked off as a fountain/wind */
  if(SphP[i].MultiPhase == 0 || SphP[i].MultiPhase == -1 || SphP[i].MultiPhase == -7
     || SphP[i].MultiPhase == -2)
    if(get_random_number(P[i].ID) < FOUNTAIN_PROBABILITY)
      {
	if(DELAY_IN_TDYN > 0.0)
	  SphP[i].DelayTime = DELAY_IN_TDYN * SphP[i].tdyn;
	else
	  SphP[i].DelayTime = -DELAY_IN_TDYN * YEAR_IN_CODE_UNITS - CLOCK * SphP[i].tdyn;
	SphP[i].wait = 1;
#ifdef LT_STELLAREVOLUTION
	SphP[i].mstar = 0.0;
#endif
      }
  /* as long as SphP[i].DelayTime>0.0 the particle will not enter MP and will receive kinetic energy */
#endif


  /*  STATISTICAL STARFORMATION */
  make_a_star(i, sm, GENERATIONS_BUNCH);

  /* resets MP particle if MP phase is over */
  if(SphP[i].MultiPhase <= 0)
    {
      SphP[i].MultiPhase = 0;
      SphP[i].M_sf = 0.0;
      SphP[i].E_out = 0.0;
      SphP[i].tdyn = 0.0;

      SphP[i].E_kin = 0.0;
      SphP[i].xkin = 0.0;
      SphP[i].ykin = 0.0;
      SphP[i].zkin = 0.0;

#ifdef LT_STELLAREVOLUTION
      SphP[i].mstar = 0.0;
#endif
    }


#ifdef GM_COUNT_PARTICLES_IN_CONE
  /* istogramma delle particelle a cui si da energia */
  if(SphP[i].nnorm < 0)
    histog[20]++;		/* should not happen */
  else if(SphP[i].nnorm < 10)
    histog[SphP[i].nnorm]++;
  else if(SphP[i].nnorm < 100)
    histog[SphP[i].nnorm / 10 + 9]++;
  else
    histog[19]++;
#endif

#ifdef GM_MUPPI_DEBUG
  All.Egy_out += SphP[i].E_out;
#endif


}



/* This function actually evaluate the probability of making a star, and produces it when needed.
   Called by starformation()*/
static inline void make_a_star(int i, MyAtLeastDouble sm, int GENERATIONS_BUNCH)
{

  MyAtLeastDouble p, xp, prob, mass_of_star;
  double NextChemTime;
  int index_of_star, bin, number_of_stars_generated, chem_step;



  /* star spawning is based on the acquired stellar mass AFTER restoration in the IRA */
  /* SphP[i].M_sf set by MUPPI */
  p = (SphP[i].M_sf - sm) / P[i].Mass;	/* probabilistic treatment for spawning new stars, pag 11 SH */

  Sum_sm += P[i].Mass * (1 - exp(-p));
  if(p < 0.0)
    p = 0.0;

  /* the higher the mass of formed stars the higher the probability to spawn a new star */

  /* the upper bits of the gas particle ID store how many stars this gas
     particle gas already generated */
  if(bits == 0)
    number_of_stars_generated = 0;
  else
    number_of_stars_generated = (P[i].ID >> (32 - bits));	/* here it counts the number of already generated  stars,
								 * the number of stars is incremented later */

#ifndef LT_STELLAREVOLUTION
  mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);
  /* if GENERATION=2, one gas particle forms star particles of equal mass */
  /* Volker's definition:
     SphP[i].Sfr = (1 - All.FactorSN) * cloudmass / tsfr *                      eq. 1 pag 3   SH multi
     (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
   */

  /* compute SFR only if some star has formed - this is to avoid to form stars after errors */
#else
  mass_of_star = P[i].Mass / (All.Generations - number_of_stars_generated) * GENERATIONS_BUNCH;
#endif
  if(SphP[i].M_sf > sm)
    {

#ifdef LT_STELLAREVOLUTION
      SphP[i].mstar = (SphP[i].M_sf - sm);	/* this tells to the chemical part 
						   to account for metals continuosly produced
						   by the MP particles */
      /* we give LT stellar mass BEFORE restoration */
#endif

      /* setting particle's SFR */
      SphP[i].Sfr = (SphP[i].M_sf - sm) / dtime / SOLAR_MASS_IN_CODE_UNITS * YEAR_IN_CODE_UNITS;
#ifndef LT_STELLAREVOLUTION
      SphP[i].Sfr /= (1. - f_rest);	/* to have sfr before rest if stellarevolution not used */
#endif

    }
  /* M_sf in g, dtime physical internal time units */
  else
    {
      SphP[i].Sfr = 0.0;
#ifdef LT_STELLAREVOLUTION
      SphP[i].mstar = 0.0;
#endif
    }
  TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;

  prob = P[i].Mass / mass_of_star * (1 - exp(-p));	/*  eq. 39 pag 11 ... spawning a new star */
  xp = get_random_number(P[i].ID + 1);

  if(xp < prob)			/* ok, make a star */
    {

#ifdef LT_STELLAREVOLUTION
      NextChemTime = get_NextChemTime(0, SFi, 0x0);
      if((number_of_stars_generated + GENERATIONS_BUNCH) >= All.Generations)
#else
      if(number_of_stars_generated == (GENERATIONS - 1))
#endif
	index_of_star = convert_gas_into_star_particle(i);
      else
	index_of_star = produce_star_particle(i, mass_of_star, GENERATIONS_BUNCH);

#ifdef LT_STELLAREVOLUTION
      bin = get_chemstep_bin(All.Time, All.Time_Age - NextChemTime, &chem_step, index_of_star);
      MetP[P[index_of_star].pt.MetID].ChemTimeBin = bin;

      TimeBinCountStars[bin]++;
    }
#endif
//  printf("Task %d make_a_star completed i=%d N_stars=%d stars_Spawned=%d index_of_stars=%d\n", ThisTask, i, N_stars, stars_spawned, index_of_star); fflush(stdout);

}


/* This function converts a gas particle in a star particle. 
   Called by make_a_star() if probability check is passed and the gas particle be converted. */
static inline int convert_gas_into_star_particle(int i)
{
  int index_of_star;

  /* here we convert the gas particle itself into a star */
  Stars_converted++;
  stars_converted++;

  Sum_mass_stars += P[i].Mass;	/* total mass in stars */
  P[i].Type = 4;
  TimeBinCountSph[P[i].TimeBin]--;
  TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;

#ifdef LT_STELLAREVOLUTION
  if(P[i].Hsml > All.MaxChemSpreadL)
    P[i].Hsml = All.MaxChemSpreadL;

  if(N_stars + 1 >= All.MaxPartMet)
    {
      printf
	("On Task=%d with NumPart=%d we try to convert %d particles. Sorry, no space left...(All.MaxPartMet=%d)\n",
	 ThisTask, NumPart, stars_converted, All.MaxPartMet);
      fflush(stdout);
      endrun(8889);
    }

  P[i].pt.MetID = N_stars;
  MetP[N_stars].PID = i;
  MetP[N_stars].iMass = P[i].Mass;	/* mass of stars BEFORE rest */

  for(int j = 0; j < LT_NMetP; j++)
    MetP[N_stars].Metals[j] = SphP[i].Metals[j];
#ifdef GL_CR_DUST
  for(int j = 0; j < LT_NMetP; j++)
    MetP[N_stars].Metals[j] += SphP[i].DustL[j] + SphP[i].DustS[j];
#endif

  MetP[N_stars].LastChemTime = SNtimesteps[SFi][0][1];
  P[i].Mass *= (1 - SFs[SFi].metFactorSN);
#endif


#ifdef STELLARAGE
#ifdef LT_STELLAREVOLUTION	/* function takes into account comoving/physical setup */
  MPP(i).StellarAge = (MyFloat) All.Time;
  if(MPP(i).StellarAge < 0.0)
    MPP(i).StellarAge = 0.0;	/* fix negative stardate in isolation */

  printf(" Task %d Star particle %d Age %g  All.Time %g\n", ThisTask, i, MPP(i).StellarAge, All.Time);
  fflush(stdout);

#else
  MPP(i).StellarAge = All.Time;
#endif
#endif


  index_of_star = i;
  N_stars++;

  return index_of_star;
}


/* This function spawns a new particle. Called by make_a_star() if probability check is passed and
   the gas particle must not be converted. */
static inline int produce_star_particle(int i, MyAtLeastDouble mass_of_star, int GENERATIONS_BUNCH)
{
  int index_of_star;
  MyAtLeastDouble mass_excess, factor;

  /* here we spawn a new star particle */
  if(NumPart + stars_spawned >= All.MaxPart)	/* the gas particle still exists, with less energy and mass */
    /* in the hot and cold phase */
    {
      printf
	("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
	 ThisTask, NumPart, stars_spawned, All.MaxPart);
      fflush(stdout);
      endrun(8888);
    }

  P[NumPart + stars_spawned] = P[i];	/*  stars_spawned incremented at the end of the loop */
  P[NumPart + stars_spawned].Type = 4;


  NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
  FirstActiveParticle = NumPart + stars_spawned;
  NumForceUpdate++;

  TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;

  PrevInTimeBin[NumPart + stars_spawned] = i;
  NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
  if(NextInTimeBin[i] >= 0)
    PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
  NextInTimeBin[i] = NumPart + stars_spawned;
  if(LastInTimeBin[P[i].TimeBin] == i)
    LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;
  P[NumPart + stars_spawned].Mass = mass_of_star;

#ifndef LT_STELLAREVOLUTION
  P[i].ID += ((unsigned) 1 << (32 - bits));	/* !!! add one to the upper bits (those for storing num of *) of the part ID */
  number_of_stars_generated = (P[i].ID >> (32 - bits));
#else
  factor = P[NumPart + stars_spawned].Mass / P[i].Mass;
  P[i].ID += (GENERATIONS_BUNCH << (32 - bits));
#endif

  Sum_mass_stars += P[NumPart + stars_spawned].Mass;
  P[i].Mass -= P[NumPart + stars_spawned].Mass;


  /* Removes the mass of the spawned star from the stellar mass component */
  mass_excess = mass_of_star - SphP[i].M_sf;
  SphP[i].M_sf -= mass_of_star;

  /* if the stellar mass is not enough, take away mass from the cold component */
  /* (without stocastic model it should never happen) */
  if(SphP[i].M_sf < SOLAR_MASS_IN_CODE_UNITS_h100 || SphP[i].MultiPhase == -10)	/* in solar masses!! */
    {
      SphP[i].M_sf = SOLAR_MASS_IN_CODE_UNITS_h100;
      mass_excess += SOLAR_MASS_IN_CODE_UNITS_h100;	/* :-) */
      SphP[i].M_c -= mass_excess;

      /* if the cold mass is depopulated as well, the particle can't stay multiphase */
      if(SphP[i].M_c < 0.0 || SphP[i].MultiPhase == -10)
	{
	  SphP[i].M_sf = 0.0;
	  SphP[i].MultiPhase = 0;
	  SphP[i].M_h = 0.0;
	  SphP[i].M_c = 0.0;
	  SphP[i].E_h = 0.0;
	  SphP[i].E_out = 0.0;

	  SphP[i].E_kin = 0.0;
	  SphP[i].xkin = 0.0;
	  SphP[i].ykin = 0.0;
	  SphP[i].zkin = 0.0;
#ifdef LT_STELLAREVOLUTION
	  SphP[i].mstar = 0.0;
#endif

	  nexit_spawn++;
	}
    }


#ifdef LT_STELLAREVOLUTION
  if(N_stars + 1 >= All.MaxPartMet)
    {
      printf
	("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPartMet=%d)\n",
	 ThisTask, NumPart, stars_spawned, All.MaxPartMet);
      fflush(stdout);
      endrun(8889);
    }

  if(P[NumPart + stars_spawned].Hsml > All.MaxChemSpreadL)
    P[NumPart + stars_spawned].Hsml = All.MaxChemSpreadL;

  MetP[N_stars].iMass = P[NumPart + stars_spawned].Mass;

  P[NumPart + stars_spawned].Mass *= (1 - SFs[SFi].metFactorSN);
  P[NumPart + stars_spawned].pt.MetID = (unsigned int) N_stars;
  MetP[N_stars].PID = NumPart + stars_spawned;
  MetP[N_stars].LastChemTime = SNtimesteps[SFi][0][1];


  for(int j = 0; j < LT_NMetP; j++)
    {
      MetP[N_stars].Metals[j] = SphP[i].Metals[j] * factor;
      SphP[i].Metals[j] *= (1 - factor);
#ifdef GL_CR_DUST
      MetP[N_stars].Metals[j] += factor * (SphP[i].DustL[j] + SphP[i].DustS[j]);
      SphP[i].DustL[j] -= factor * SphP[i].DustL[j];
      SphP[i].DustS[j] -= factor * SphP[i].DustS[j];
#endif
    }
#endif




#ifdef STELLARAGE
  MPP(NumPart + stars_spawned).StellarAge = (MyFloat) All.Time;
  if(MPP(NumPart + stars_spawned).StellarAge < 0.0)
    MPP(NumPart + stars_spawned).StellarAge = 0.0;	/* fix negative stardate in isolation */
  printf("   star %d stellarage %g %g\n", NumPart + stars_spawned, MPP(NumPart + stars_spawned).StellarAge,
	 MPP(i).StellarAge);
  fflush(stdout);
#endif


  index_of_star = NumPart + stars_spawned;
  N_stars++;
  force_add_star_to_tree(i, NumPart + stars_spawned);
  stars_spawned++;


  return index_of_star;
}


void final_statistics(MyAtLeastDouble time_hubble_a, int mfin, int mfout)
{

  int i, bin;
  int tot_spawned, tot_converted;
  int texit_dens, texit_clock, texit_mnc, texit_freeze, texit_depl, texit_spawn, texit_gsl, texit_sfraw,
    texit_tms, texit_force, texit_wind, totnmf, totmfin, totmfout;
  double sfrrate, rate, totsfrrate, rate_in_msunperyear, total_sm, total_sum_mass_stars;


  texit_dens = texit_clock = texit_freeze = texit_depl = texit_spawn = texit_gsl = texit_sfraw = texit_mnc =
    texit_tms = texit_force = texit_wind = 0;
  totnmf = totmfin = totmfout = 0;


#ifdef COUNT_PARTICLES_IN_CONE
  MPI_Allreduce(&histog[0], &tothistog[0], 21, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Energy_total, &S_Energy_total, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Energy_used, &S_Energy_used, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Energy_total_1, &S_Energy_total_1, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&Energy_used_1, &S_Energy_used_1, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
#endif

  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  MPI_Allreduce(&nmf, &totnmf, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&mfin, &totmfin, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&mfout, &totmfout, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  MPI_Allreduce(&nexit_dens, &texit_dens, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_clock, &texit_clock, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_freeze, &texit_freeze, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_depl, &texit_depl, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_spawn, &texit_spawn, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_gsl, &texit_gsl, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_sfraw, &texit_sfraw, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_mnc, &texit_mnc, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_tms, &texit_tms, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_force, &texit_force, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&nexit_wind, &texit_wind, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);


#ifdef MV_GM_AGNMUPPI_OUTPUT
  {
    double TotE_AGNM_tot, TotE_AGNM_cool, TotE_AGNM_sfr_h, TotE_AGNM_sfr_c;
    double TotE_AGNM_sfr_dE, TotE_AGNM_sfr_dE_aft, TotE_AGNM_sfr_BHEcoldD;
    int Tot_n_AGNM, Tot_n_AGNcool;
    MPI_Allreduce(&E_AGNM_tot, &TotE_AGNM_tot, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&E_AGNM_cool, &TotE_AGNM_cool, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&E_AGNM_sfr_h, &TotE_AGNM_sfr_h, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&E_AGNM_sfr_c, &TotE_AGNM_sfr_c, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&E_AGNM_sfr_dE, &TotE_AGNM_sfr_dE, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&E_AGNM_sfr_dE_aft, &TotE_AGNM_sfr_dE_aft, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&E_AGNM_sfr_BHEcoldD, &TotE_AGNM_sfr_BHEcoldD, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&n_AGNM, &Tot_n_AGNM, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(&n_AGNcool, &Tot_n_AGNcool, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);


    All.TotE_AGNM_tot += TotE_AGNM_tot;
    All.TotE_AGNM_cool += TotE_AGNM_cool;
    All.TotE_AGNM_sfr_c += TotE_AGNM_sfr_c;
    All.TotE_AGNM_sfr_h += TotE_AGNM_sfr_h;
    if(ThisTask == 0)
      {
	fprintf(FdEgyAGNtot, "%g         %e         %e          %e     %e             %e     %e             %e         %d        %d\n", All.Time, TotE_AGNM_tot, TotE_AGNM_cool, TotE_AGNM_sfr_c, TotE_AGNM_sfr_h, TotE_AGNM_sfr_dE, TotE_AGNM_sfr_dE_aft,	//sommati su tutte le particelle MP che han preso energia da AGN
		TotE_AGNM_sfr_BHEcoldD, Tot_n_AGNM, Tot_n_AGNcool);
	fflush(FdEgyAGNtot);

      }
  }
#endif



  if(ThisTask == 0)
    {
      printf(" ----> MultiPhase particles: %d ", totnmf);
      printf("\n\n  --->MF: %d, MF-->:, %d\n", totmfin, totmfout);
      printf
	("\n\n Exits from MP( Time %f ), density: %d, clock: %d, freeze: %d, spawining: %d, force SP: %d, wind:     %d, depleted: %d, gsl failure: %d, SF runaway: %d, mass: %d, long int: %d\n",
	 All.Time, texit_dens, texit_clock, texit_freeze, texit_spawn, texit_force, texit_wind, texit_depl,
	 texit_gsl, texit_sfraw, texit_mnc, texit_tms);
      fflush(stdout);

      fprintf(FdExit, " %8f %8d %8d %8d %8d %8d %8d     %8d %8d %8d %8d %8d\n", All.Time,
	      texit_dens, texit_clock, texit_freeze, texit_spawn, texit_force, texit_wind, texit_depl,
	      texit_gsl, texit_sfraw, texit_mnc, texit_tms);
      fflush(FdExit);

#ifdef GM_COUNT_PARTICLES_IN_CONE
      totnactive = 0;

      printf(" Histogram of particles in cone: %f ", All.Time);
      for(i = 0; i < 20; i++)
	{

	  totnactive += tothistog[i];
	  printf(" (N=%d): %d ", (i < 10 ? i : (i - 9) * 10), tothistog[i]);
	}
      printf("  NEGATIVI: %d", tothistog[20]);
      printf("  TOTALE: %d\n", totnactive);
      fflush(stdout);

      fprintf(FdCone, " %8f ", All.Time);

      fprintf(FdCone, " %15g %15g %15g %15g",
	      S_Energy_total, S_Energy_used, S_Energy_total_1, S_Energy_used_1);

      for(i = 0; i < 21; i++)
	fprintf(FdCone, " %8d ", tothistog[i]);
      fprintf(FdCone, "   %8d\n", totnactive);

      fflush(FdCone);

#endif

    }

  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("\n----> spawned %d stars, converted %d gas particles into stars\n\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}

      All.TotNumPart += tot_spawned;	/* updating number of particles of different types */
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;
#ifdef LT_STELLAREVOLUTION
      All.TotN_stars += tot_spawned + tot_converted;
#endif

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */
      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Reduce(&Sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&Sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    {

      /* computation of sfr */
      if(All.TimeStep > 0)
	rate = total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;
      /* convert to solar masses per yr */
      rate_in_msunperyear = rate / SOLAR_MASS_IN_CODE_UNITS * YEAR_IN_CODE_UNITS;
#ifndef LT_STELLAREVOLUTION
      /* converto into SFR BEFORE restoration */
      rate_in_msunperyear /= (1. - f_rest);
#endif

      fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear,
	      total_sum_mass_stars);

      fflush(FdSfr);
    }


}


MyAtLeastDouble get_starformation_rate(int i)
{
  /* returns the particle's SFR */

  return (MyAtLeastDouble) SphP[i].Sfr;

}



void rearrange_particle_sequence(void)
{
  int i, j, flag = 0, flag_sum;

  struct particle_data psave;

#ifdef BLACK_HOLES
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
#endif

  if(Stars_converted)
    {
      N_gas -= Stars_converted;
      Stars_converted = 0;

      for(i = 0; i < N_gas; i++)
	if(P[i].Type != 0)
	  {
	    for(j = N_gas; j < NumPart; j++)
	      if(P[j].Type == 0)
		break;

	    if(j >= NumPart)
	      endrun(181170);
#ifdef LT_STELLAREVOLUTION
	    if(P[i].Type == 4)
	      MetP[P[i].pt.MetID].PID = j;
#endif

#if defined(BLACK_HOLES)
	    if(P[i].Type == 5)
	      BHP[P[i].pt.BHID].PID = j;
#endif

	    psave = P[i];
	    P[i] = P[j];
	    SphP[i] = SphP[j];
	    P[j] = psave;
	  }
      flag = 1;
    }

#ifdef BLACK_HOLES
  count_elim = 0;
  count_gaselim = 0;

  for(i = 0; i < NumPart; i++)
    if(P[i].Mass == 0)
      {
	TimeBinCount[P[i].TimeBin]--;

	if(TimeBinActive[P[i].TimeBin])
	  NumForceUpdate--;

	if(P[i].Type == 0)
	  {
	    TimeBinCountSph[P[i].TimeBin]--;

	    P[i] = P[N_gas - 1];
	    SphP[i] = SphP[N_gas - 1];

	    P[N_gas - 1] = P[NumPart - 1];
	    if(P[N_gas - 1].Type == 4)
	      MetP[P[N_gas - 1].pt.MetID].PID = N_gas - 1;

	    if(P[N_gas - 1].Type == 5)
	      BHP[P[N_gas - 1].pt.BHID].PID = N_gas - 1;

	    N_gas--;

	    count_gaselim++;
	  }
	else if(P[i].Type == 4)
	  {
	    j = P[i].pt.MetID;
	    MetP[j] = MetP[N_stars - 1];
	    P[MetP[j].PID].pt.MetID = j;
	    N_stars--;

	    P[i] = P[NumPart - 1];

	    if(P[i].Type == 4)
	      MetP[P[i].pt.MetID].PID = i;

	    if(P[i].Type == 5)
	      BHP[P[i].pt.BHID].PID = i;
	  }
	else if(P[i].Type == 5)
	  {
	    j = P[i].pt.BHID;
	    BHP[j] = BHP[N_BHs - 1];
	    P[BHP[j].PID].pt.BHID = j;
	    N_BHs--;

	    P[i] = P[NumPart - 1];

	    if(P[i].Type == 4)
	      MetP[P[i].pt.MetID].PID = i;

	    if(P[i].Type == 5)
	      BHP[P[i].pt.BHID].PID = i;
	  }
	else
	  {
	    P[i] = P[NumPart - 1];

	    if(P[i].Type == 4)
	      MetP[P[i].pt.MetID].PID = i;

	    if(P[i].Type == 5)
	      BHP[P[i].pt.BHID].PID = i;
	  }

	NumPart--;
	i--;

	count_elim++;
      }

  MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  if(count_elim)
    flag = 1;

  if(ThisTask == 0)
    {
      printf("Blackholes: Eliminated %d gas particles and merged away %d black holes.\n",
	     tot_gaselim, tot_elim - tot_gaselim);
      fflush(stdout);
    }

  All.TotNumPart -= tot_elim;
  All.TotN_gas -= tot_gaselim;
  All.TotBHs -= tot_elim - tot_gaselim;
#endif

  MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  if(flag_sum)
    reconstruct_timebins();
}


void set_units_sfr(void)
{

#ifndef LT_STELLAREVOLUTION

  MyAtLeastDouble meanweight;

  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;	/* This is in h true */

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
  All.UnitSurfDens_in_Msunpc2 =
    1.0 / (SOLAR_MASS / All.UnitMass_in_g / pow(CM_PER_PC / All.UnitLength_in_cm, 2));
#endif


  if(ThisTask == 0)
    {
      printf("OverDensThresh= %g\nPhysDensThresh= %g (internal units)\n", All.OverDensThresh,
	     All.PhysDensThresh);
      printf("\n MUPPI parameters:\n");
#ifdef LT_STELLAREVOLUTION
      printf("f_star=%f, f_evap=%f,f_rest=%f, beta_sn=%g, E_SN_51=%g\n",
	     f_star, f_evap, All.f_rest[0], All.beta_sf[0], All.E_SN_51[0]);
#else
      printf("f_star=%f, f_evap=%f,f_rest=%f, beta_sn=%g, E_SN_51=%g\n",
	     f_star, f_evap, myf_rest, mybeta_sf, myE_SN_51);
#endif

      printf(" Frac Egy in, out, kin: %f %f %f\n", All.FracEgyIn, All.FracEgyOut, All.FracEgyKin);
      fflush(stdout);
    }
}



/* This funcion is the integrator of our MUPPI PDEs system */
void MUPPI_main(int i, MyAtLeastDouble Rho, int *mfin, int *mfout)
{
  /* This routine is a driver for MUltiphase Particle Integration 
   * input: i (particle ID)
   *        Rho (density in code units)
   * *mfin, *mfout: number of particles entering/exiting MultiPhase 
   */

  /* NB: this routine works in h true */

  if(P[i].ID == 1000921906 && All.Time > 0.709132)
    printf("Eccomi in Muppi\n");

  //  checkline("Entering MUPPI_main");

  MyAtLeastDouble dblMuppiVars[nMuppiDbls];
  //  MyAtLeastDouble *dblMuppiVars;
  //  dblMuppiVars = (MyAtLeastDouble*)mymalloc("MUPPI DBLS",sizeof(MyAtLeastDouble) * nMuppiDbls);

  //  int minternal, mpart_ID, mnstep, mmfin, mmfout, mncall;
  int intMuppiVars[nMuppiInts];
  //  int *intMuppiVars;
  //  intMuppiVars = (int*)mymalloc("MUPPI INTS",sizeof(int) * nMuppiInts);

  /* all variables of MUPPI_main and daughter functions defined in muppi.h and allocated here */
  /* WARNING: these areas MUST keep their names in functions, because the code uses alias
     defined in muppi.h on these names */

  //  if(dblMuppiVars==NULL) printf(" NULL DOUBLES!\n"); fflush(stdout);
  //  if(intMuppiVars==NULL) printf(" NULL INTS!\n"); fflush(stdout);


#ifdef MV_GM_AGNMUPPI
  neqz = 5;
#else
  neqz = 4;
#endif


  double y[neqz];		/* NOTE THAT THIS MUST BE <<double>> */
  struct myparams param;
  int derivatives(double, const double[], double *, void *);
  int jac(double, const double[], double *, double[], void *);
  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(T, neqz);
  gsl_odeiv2_control *c = gsl_odeiv2_control_standard_new(1.0e-4, 1.0e-8, 1.0, 1.0);
  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(neqz);
  gsl_odeiv2_system sst = { derivatives, jac, neqz, (void *) &param };


  /* NOTE: we are using gsl integration functions, and they want <<double>> args */

  initialize_muppi_vars(i, &Rho, dblMuppiVars, intMuppiVars);

  if(SphP[i].MultiPhase == 0)	/* first MUPPI step */
    start_MUPPI(i, dblMuppiVars, intMuppiVars);
  else				/* subsequent MUPPI steps */
    {
      mT_h = SE_TO_T * SphP[i].E_h / SphP[i].M_h;
    }

  /* exit condition from multiphase: density threshold */
  if(check_exit_muppi(i, Rho, dblMuppiVars, intMuppiVars) == 1)
    {
      gsl_odeiv2_evolve_free(e);
      gsl_odeiv2_control_free(c);
      gsl_odeiv2_step_free(s);

      //      CheckMuppiMemory(i, Rho, mfin, mfout, dblMuppiVars, intMuppiVars);

      //      myfree(intMuppiVars);
      //      myfree(dblMuppiVars);
      return;
    }


  /* initializes the integration */
  /* ALSO: sets all ISM properties. Currently not all used, here for future possible uses */
  muppi_initialize_integration(i, Rho, dblMuppiVars, intMuppiVars, y, &param);


  /*
     {
     char bfr[1024];
     FILE *fdbg;

     if(All.Time>0.0507727)
     {
     sprintf(bfr,"dbg/dbg%02d.txt",ThisTask);
     fdbg=fopen(bfr,"a");
     fprintf(fdbg, "%f %d %u %e %e %e %d %d %e %e %e %e %e %d %e %e %e %e %e %e %e %d %e %e\n",
     All.Time, i, P[i].ID, SphP[i].Egy_tot_0, 
     SphP[i].tdyn,SphP[i].clock, SphP[i].MultiPhase, SphP[i].NMF,
     SphP[i].Ekin_rec[0], SphP[i].Ekin_rec[1], SphP[i].Ekin_rec[2],
     SphP[i].Ekin_norm, SphP[i].Ekin_norm_fountain, SphP[i].wait, SphP[i].DelayTime0,
     SphP[i].Mgas0, SphP[i].GradDens[0], SphP[i].GradDens[1], SphP[i].GradDens[2],
     SphP[i].E_out, SphP[i].E_rec, SphP[i].NoArtCond, SphP[i].Temp, SphP[i].EgyStep); fflush(fdbg);
     fclose(fdbg);
     }
     } 
   */

  /* MAIN INTEGRATION LOOP */
  /* calls the integration routine */
  //  checkline("beginning integration");
  while(mx1 < mx2)
    {
      mhh_old = mhh;

      /* gsl integrator */
      /* WARNING: here we must be sure that arguments are <<double>> */
      double dmx1, dmx2, dmhh;
      dmx1 = (double) mx1;
      dmx2 = (double) mx2;
      dmhh = (double) mhh;
      int status = gsl_odeiv2_evolve_apply(e, c, s, &sst, &dmx1, dmx2, &dmhh, y);
      mnstep += 1;
      mx1 = (MyAtLeastDouble) dmx1;
      mhh = (MyAtLeastDouble) dmhh;

      /* check exits */
      if(muppi_integration_check_exits(i, status, dblMuppiVars, intMuppiVars, y) < 0)
	break;

    }				/* end of WHILE cycle on integration timesteps */


  /* if an exit condition arose in integration, exec it */
  if(SphP[i].MultiPhase < 0)
    {
      muppi_exec_exits(i, dblMuppiVars, intMuppiVars, y);


      gsl_odeiv2_evolve_free(e);
      gsl_odeiv2_control_free(c);
      gsl_odeiv2_step_free(s);

      *mfin += mmfin;
      *mfout += mmfout;

      //      myfree(intMuppiVars);
      //      myfree(dblMuppiVars);
      return;
    }

  /* othewise make final operations, conversions and the like */
  *mfin += mmfin;
  *mfout += mmfout;
  muppi_final_step_operations(i, Rho, dblMuppiVars, intMuppiVars, y);

  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);
  fflush(stdout);

#ifdef GM_MUPPI_DEBUG
  check_muppi_nan(i);
#endif

  //  myfree(intMuppiVars);
  //  myfree(dblMuppiVars);
  //  checkline("MUPPI_main ENDS"); fflush(stdout);

}


int jac(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  printf("This integration method should not call this function.\n");
  printf("Better stop here.\n");
  endrun(767676);
  return GSL_SUCCESS;
}


/* this function is the core of the model */
/* if normalization is active, ALL energies and times here
   are normalized, only in OUTPUT and when they enters
   calculation of other quantities (...Th...) they are converted
   back to CGS */
int derivatives(double x, const double y[], double *dydx, void *par)
{
  int j, ichk;
  MyAtLeastDouble E_h, M_h, fill_h, T_h;
  MyAtLeastDouble M_s, M_c, Rho_h, n_c, Rho_c, n_h;
  MyAtLeastDouble dMcool, dEcool, dEsn, dEhydro, dM;
  MyAtLeastDouble dMevap, dMrest, dMsf;
  MyAtLeastDouble tdyn;
  MyAtLeastDouble tcool, Fcoll;

  MyAtLeastDouble FBVol_tot;
  int FBinternal, FBpart_ID, FBfiltered, FBncall;

  MyAtLeastDouble *dblMuppiVars;
  int *intMuppiVars;

  struct myparams *param;

  MyAtLeastDouble dumRho = 0.0;
  int dummfin = 0, dummfout = 0;

#ifdef MV_GM_AGNMUPPI
  double BHEcold_used, dBHEcold_used;
  double BHEcoldD_surplus, dMBH;	//change in cold mass due to the BH energy
  double E_c_AGN_used, dBHEcold;
#endif


#ifdef MV_UPDATE_TDYN_FOR_STARFORMATION
  MyAtLeastDouble upd_tdyn;
#endif

#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
  MyAtLeastDouble sdensity, metallicity, metal_mass, rho, Fcoll_BeR;
#endif


  param = (struct myparams *) par;

  /* local parameters */
  param->ncall += 1;
  FBncall = param->ncall;
  dEhydro = param->dE;
  dM = param->dM;
  FBVol_tot = param->VolTot;
  FBinternal = param->internal;
  FBpart_ID = param->part_ID;
  FBfiltered = param->filtered;
  dblMuppiVars = param->dblsmuppi;
  intMuppiVars = param->intsmuppi;

#ifdef MV_GM_AGNMUPPI
  dBHEcold = mBHEcoldD;		//mile6
  dBHEcold_used = 0.0;		//qui accumuliamo l'energia USATA 
#endif


  /* local names of variables */
  M_s = (MyAtLeastDouble) y[0];
  M_c = (MyAtLeastDouble) y[1];
  M_h = (MyAtLeastDouble) y[2];
  E_h = (MyAtLeastDouble) y[3];
#ifdef MV_GM_AGNMUPPI
  E_c_AGN_used = y[4];		//mile6
#endif

#ifdef MV_GM_AGNMUPPI
  ichk =
    derivatives_initial_check(x, E_h, M_h, M_c, dEhydro, FBVol_tot, M_s, E_c_AGN_used, dBHEcold, FBpart_ID,
			      FBinternal, &FBfiltered, dydx, dblMuppiVars);
#else
  ichk =
    derivatives_initial_check(x, E_h, M_h, M_c, dEhydro, FBVol_tot, M_s, FBpart_ID, FBinternal, &FBfiltered,
			      dydx, dblMuppiVars);
#endif


  if(ichk == -1)
    return GSL_FAILURE;
  else if(ichk == -2)
    {
      param->filtered = FBfiltered;
      return GSL_SUCCESS;
    }

  /* ISM */
  compute_ism_properties(M_h, E_h, M_c, FBVol_tot, &T_h, &fill_h, &Rho_h, &Rho_c, &n_h, &n_c);
  param->fill_c = 1.0 - fill_h;



  /* this it the behaviour of the system in case of maximal starburst */
  /* SFR is directly connected to the hot phase */
  /* cold phase is evacuated */
  if(SphP[FBinternal].t_startMP < 0.0)
    {
      derivatives_maximal_starburst(dM, n_h, T_h, E_h, M_h, Rho_h, M_c, FBpart_ID, FBinternal, y, dydx,
				    dblMuppiVars);
      return GSL_SUCCESS;
    }

  /* cooling time */
  tcool = cooling_time(T_h, E_h / M_h, Rho_h, n_h);

  /* exit if cooling time is unreasonable */
  if(tcool < 0.0 || n_h <= 0.0)
    {
      for(j = 0; j < 4; j++)
	dydx[j] = 0.;

#ifdef GM_MUPPI_WARNINGS
      printf
	("WARNING: task %d  part %d got wrong cooling time: Th=%g Eh=%g Mh=%g Rhoh=%g nh=%g tcool=%g call:%d\n",
	 ThisTask, FBpart_ID, T_h, E_h, M_h, Rho_h, n_h, tcool, FBncall);
#endif

      return GSL_SUCCESS;
    }


  /* star formation rate */
  if(M_c < SOLAR_MASS_IN_CODE_UNITS)
    dMsf = 0.0;			/* no SF if there is no cold phase (less than 1 Msun) */
  else
    {

      /* dynamical time */
#ifdef MV_UPDATE_TDYN_FOR_STARFORMATION
      compute_ism_properties(y[2], y[3], y[1], mVol_tot, &T_h, &fill_h, &Rho_h, &Rho_c, &n_h, &n_c);
      upd_tdyn = 5.15e7 / sqrt(mu_c * n_c) * YEAR_IN_CODE_UNITS;
#endif

      if(SphP[FBinternal].tdyn == 0.0)
	tdyn = 5.15e7 / sqrt(mu_c * n_c) * YEAR_IN_CODE_UNITS;
      else
	tdyn = SphP[FBinternal].tdyn;

      /* molecular fraction */
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
      rho = mRho;
      /* surface density in Sobolev approximation */
      sdensity = rho * rho / sqrt(SphP[FBinternal].GradDens[0] * SphP[FBinternal].GradDens[0] +
				  SphP[FBinternal].GradDens[1] * SphP[FBinternal].GradDens[1] +
				  SphP[FBinternal].GradDens[2] * SphP[FBinternal].GradDens[2]);
      sdensity /= pow(All.HubbleParam, 3) / pow(All.Time, 4);
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION_CLUMPINESS
      double eps_dave = 0.125 / All.HubbleParam / 1000;	// comoving used in Dave, converted to Mpc
      double mgas_dave = 2.85e5 / All.HubbleParam;	// gas particle mass, in solar masses
      double rho_max_dave_com = mgas_dave / pow(eps_dave, 3);	// maximum comoving density from Dave
      double mgas_aq5 = 4.11e5;	// mass resolution for AqC5, in solar masses
      double rho_max_aq5_com = mgas_aq5 / pow((All.SofteningTable[0] / All.HubbleParam), 3);	// density for AqC5, to compute a scaling factor
      double cfak_to_use = rho_max_dave_com / rho_max_aq5_com * MV_KRUMHOLZ_MOLECULAR_FRACTION_CLUMPINESS;
      if(cfak_to_use < 1.0)
	cfak_to_use = 1;
      sdensity *= cfak_to_use;
#endif
      metal_mass = 0.0;
      for(j = 1; j < LT_NMetP; j++)
	metal_mass = metal_mass + SphP[FBinternal].Metals[j];
      metallicity = metal_mass / P[FBinternal].Mass;

      Fcoll = molecular_fraction_Krum(metallicity, sdensity);

      Fcoll_BeR = molecular_fraction(n_c * T_c);
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION_OUTPUT
      /* fprintf(FdMolFrac,"# 1.case   2.time    3.particleID    4.H2_frac_Krumholz    5.H2_frac_Blitz&R */
      /*               6.metallicity    7.number_dens hot(1)/cold(2)    8.temperature hot(1)/cold(2)\n"); */
      /* to compare different predictions on a particle basis, for particles with different properties */
      /* case 2 (1st column) means no FORCE SINGLE PHASE  */
      if(ThisTask == 0)
	{
	  fprintf(FdMolFrac, "%d    %g    %d    %g    %g     %g     %g    %g \n",
		  2, All.Time, mpart_ID, Fcoll, Fcoll_BeR, metallicity, n_c, T_c);
	  fflush(FdMolFrac);
	}
#endif
#else
      Fcoll = molecular_fraction(n_c * T_c);
#endif

#ifdef MV_UPDATE_TDYN_FOR_STARFORMATION
      dMsf = Fcoll * M_c / upd_tdyn * f_star;
#else
      dMsf = Fcoll * M_c / tdyn * f_star;
#endif

      /* this is for output purposes only */
      SphP[FBinternal].Fcoll = Fcoll;

    }

  /* restoration and evaporation */
  dMrest = f_rest * dMsf;
  dMevap = f_evap * dMsf;

#ifdef MV_GM_AGNMUPPI
  compute_ism_properties(M_h, E_h, M_c, mVol_tot, &T_h, &fill_h, &Rho_h, &Rho_c, &n_h, &n_c);	//mile1

  dMBH = dBHEcold * GAMMA_MINUS1 * mu_c * PROTONMASS / BOLTZMANN / (T_h - T_c);	//mile1: prima non c'era -T_c
  // dBHEcold is in code units (*All.UnitVelocity_in_cm_per_s ^2), PROTONMASS/BOLTZMASS/T_h in cgs, the conversion 
  //   is done once only in the assignment of FB.BHEcoldD
  // dMBH is the amount of cold mass brought to the temperature of the hot phase (in dtime)

  // dBHEcold_used is the BH energy required to bring to T_h a mass dMDB*dt_rungekutta of cold gas
  // We limit it to M_c/dt_rungekutta (FB.hh) if cold mass to be heated is more than what exists
  if(dMBH * mhh > M_c)
    dMBH = M_c / mhh;
  dBHEcold_used =
    dMBH / GAMMA_MINUS1 / mu_c / PROTONMASS * BOLTZMANN * (T_h -
							   T_c) / All.UnitVelocity_in_cm_per_s /
    All.UnitVelocity_in_cm_per_s;
  // [dMBH / GAMMA_MINUS1 / mu_c / PROTONMASS * BOLTZMANN * (T_h-T_c)] e' in cgs
  // dBHEcold_used e' ora in code units
#endif


  /* cooling rate */
  if(fabs(tcool) < SMALL)	/* tcool=0 means that the particle is heated by the UV background */
    {
      dMcool = 0.0;
      dEcool = 0.0;
    }
  else
    {
      dMcool = M_h / tcool;
      dEcool = E_h / tcool;
    }

  /* heating rate */
  dEsn = E_SN_51 * All.FracEgyIn * dMsf / beta_sf;

#if defined(MV_GM_AGNMUPPI) && defined(MV_AGNMUPPI_COOLING_OFF)
  if((mt_startCoolOff != -99.) && (mt_startCoolOff < MV_AGNMUPPI_COOLING_OFF * YEAR_IN_CODE_UNITS))	//mile3
    {
      //cosi' spengo il cooling nello stesso momento in cui la particella riceve energia
      //cioe' all'entrata dello stage MP, NON al timestep successivo
      dydx[0] = dMsf - dMrest;
      dydx[1] = -dMsf - dMevap - dMBH;
      dydx[2] = dMevap + dMrest + dM + dMBH;
      dydx[3] = dEsn + dEhydro + dBHEcold_used;	//mile6
      dydx[4] = dBHEcold_used;
    }
  else
    {
      dydx[0] = dMsf - dMrest;
      dydx[1] = dMcool - dMsf - dMevap - dMBH;
      dydx[2] = dMevap + dMrest - dMcool + dM + dMBH;
      dydx[3] = dEsn - dEcool + dEhydro + dBHEcold_used;	//mile6
      dydx[4] = dBHEcold_used;
    }

#else //means: if NOT( defined(MV_GM_AGNMUPPI) && defined(MV_AGNMUPPI_COOLING_OFF) )

#ifdef MV_GM_AGNMUPPI		//if only MV_GM_AGNMUPPI defined
  dydx[0] = dMsf - dMrest;
  dydx[1] = dMcool - dMsf - dMevap - dMBH;
  dydx[2] = dMevap + dMrest - dMcool + dM + dMBH;
  dydx[3] = dEsn - dEcool + dEhydro + dBHEcold_used;
  dydx[4] = dBHEcold_used;
#else //original system of equations
  dydx[0] = dMsf - dMrest;
  dydx[1] = dMcool - dMsf - dMevap;
  dydx[2] = dMevap + dMrest - dMcool + dM;
  dydx[3] = dEsn - dEcool + dEhydro;
#endif

#endif //if defined(MV_GM_AGNMUPPI) && defined(MV_AGNMUPPI_COOLING_OFF)


  /* system of equations ORIG 
     dydx[0] = (double) (dMsf - dMrest); 
     dydx[1] = (double) (dMcool - dMsf - dMevap);
     dydx[2] = (double) (dMevap + dMrest - dMcool + dM); 
     dydx[3] = (double) (dEsn - dEcool + dEhydro);
   */

  /* XXXX forcing zeroes */
  if(fabs(dydx[0]) < SMALL)
    dydx[0] = 0.0;
  if(fabs(dydx[1]) < SMALL)
    dydx[1] = 0.0;
  if(fabs(dydx[2]) < SMALL)
    dydx[2] = 0.0;
  if(fabs(dydx[3]) < SMALL)
    dydx[3] = 0.0;
#ifdef MV_GM_AGNMUPPI
  if(fabs(dydx[4]) < SMALL)
    dydx[4] = 0.0;
#endif

  return GSL_SUCCESS;
}


/* This function is called from derivatives
   Checks that nothing strange happens in input */

#ifdef MV_GM_AGNMUPPI
static inline int derivatives_initial_check(double x, MyAtLeastDouble E_h, MyAtLeastDouble M_h,
					    MyAtLeastDouble M_c, MyAtLeastDouble dEhydro,
					    MyAtLeastDouble FBVol_tot, MyAtLeastDouble M_s,
					    MyAtLeastDouble E_c_AGN_used, MyAtLeastDouble dBHEcold,
					    int FBpart_ID, int FBinternal, int *FBfiltered, double *dydx,
					    MyAtLeastDouble * dblMuppiVars)
#else
static inline int derivatives_initial_check(double x, MyAtLeastDouble E_h, MyAtLeastDouble M_h,
					    MyAtLeastDouble M_c, MyAtLeastDouble dEhydro,
					    MyAtLeastDouble FBVol_tot, MyAtLeastDouble M_s, int FBpart_ID,
					    int FBinternal, int *FBfiltered, double *dydx,
					    MyAtLeastDouble * dblMuppiVars)
#endif
{
  int j;
  MyAtLeastDouble Mass_tol;


  /* first filter: exclude NANs */
#ifdef MV_GM_AGNMUPPI
  if(isnan(E_h) || isinf(E_h) || isnan(M_h) || isinf(M_h) || isnan(M_c) || isinf(M_c) || isnan(dEhydro) || isinf(dEhydro) || isnan(mVol_tot) || isinf(mVol_tot) || isnan(M_s) || isinf(M_s) || isnan(E_c_AGN_used) || isinf(E_c_AGN_used) || isnan(dBHEcold) || isinf(dBHEcold))	//mile6
#else
  if(isnan(E_h) || isinf(E_h) || isnan(M_h) || isinf(M_h) || isnan(M_c) || isinf(M_c)
     || isnan(dEhydro) || isinf(dEhydro) || isnan(mVol_tot) || isinf(mVol_tot) || isnan(M_s) || isinf(M_s))
#endif
    {
      for(j = 0; j < neqz; j++)
	dydx[j] = 0.;

#ifdef GM_MUPPI_WARNINGS
      printf("WARNING: task %d got a NAN for particle %ld (originally %ld)\n", ThisTask,
	     (long int) P[FBinternal].ID, (long int) FBpart_ID);
      printf("         time: %f  %f\n", All.Time, x);
      printf("         checked variables: E_h=%g, M_h=%g, M_c=%g, M_s=%g, dEhydro=%g, FB.Vol_tot=%g\n",
	     E_h, M_h, M_c, M_s, dEhydro, FBVol_tot);
      printf("         state of particle: MultiPhase=%d, tdyn=%g\n",
	     SphP[FBinternal].MultiPhase, SphP[FBinternal].tdyn);
      printf("                            t_startMP=%g\n", SphP[FBinternal].t_startMP);
#if defined(MV_GM_STELLAR_KIN_FB2_OUTPUT)
      printf("                            DelayTime=%g\n", SphP[FBinternal].DelayTime);
#endif
#endif

      return -1;
    }

  /* second check on masses and energy; null derivatives if masses or energy are negative or too large */
  Mass_tol = 1.5 * P[FBinternal].Mass / All.HubbleParam;	/* P[].Mass is in h100 */
#ifdef MV_GM_AGNMUPPI
  if(E_h < SMALLNEG || E_c_AGN_used < SMALLNEG	//mile6 
     || M_c <= SMALLNEG || M_c > Mass_tol || M_h <= SMALLNEG || M_h > Mass_tol || M_s < SMALLNEG
     || M_s > Mass_tol)
#else
  if(E_h < SMALLNEG || M_c <= SMALLNEG || M_c > Mass_tol || M_h <= SMALLNEG || M_h > Mass_tol
     || M_s < SMALLNEG || M_s > Mass_tol)
#endif
    {
      *FBfiltered += 1;

      for(j = 0; j < neqz; j++)
	dydx[j] = 0.;

      return -2;
    }

  return 0;
}


/* This function is called from derivatives
   Estabilish the behaviour in case of maximal starburst */
static inline void derivatives_maximal_starburst(int dM, MyAtLeastDouble n_h, MyAtLeastDouble T_h,
						 MyAtLeastDouble E_h, MyAtLeastDouble M_h,
						 MyAtLeastDouble Rho_h, MyAtLeastDouble M_c, int FBpart_ID,
						 int FBinternal, const double *y, double *dydx,
						 MyAtLeastDouble * dblMuppiVars)
{

  int j;
  MyAtLeastDouble Fcoll, tdyn, dMevap, dMsf, dMrest;
#ifdef  MV_UPDATE_TDYN_FOR_STARFORMATION
  MyAtLeastDouble fill_h, Rho_c, upd_tdyn, n_c;
#endif
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
  MyAtLeastDouble sdensity, metallicity, metal_mass, rho, Fcoll_BeR;
#endif

  //  checkline("");

#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
  rho = mRho;
  /* surface density in Sobolev approximation */
  sdensity = rho * rho / sqrt(SphP[FBinternal].GradDens[0] * SphP[FBinternal].GradDens[0] +
			      SphP[FBinternal].GradDens[1] * SphP[FBinternal].GradDens[1] +
			      SphP[FBinternal].GradDens[2] * SphP[FBinternal].GradDens[2]);
  sdensity /= pow(All.HubbleParam, 3) / pow(All.Time, 4);
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION_CLUMPINESS
  double eps_dave = 0.125 / All.HubbleParam / 1000;	// comoving used in Dave, converted to Mpc
  double mgas_dave = 2.85e5 / All.HubbleParam;	// gas particle mass, in solar masses
  double rho_max_dave_com = mgas_dave / pow(eps_dave, 3);	// maximum comoving density from Dave
  double mgas_aq5 = 4.11e5;	// mass resolution for AqC5, in solar masses
  double rho_max_aq5_com = mgas_aq5 / pow((All.SofteningTable[0] / All.HubbleParam), 3);	// density for AqC5, to compute a scaling factor
  double cfak_to_use = rho_max_dave_com / rho_max_aq5_com * MV_KRUMHOLZ_MOLECULAR_FRACTION_CLUMPINESS;
  if(cfak_to_use < 1.0)
    cfak_to_use = 1;
  sdensity *= cfak_to_use;
#endif
  metal_mass = 0.0;
  for(j = 1; j < LT_NMetP; j++)
    metal_mass = metal_mass + SphP[FBinternal].Metals[j];
  metallicity = metal_mass / P[FBinternal].Mass;

  Fcoll = molecular_fraction_Krum(metallicity, sdensity);

  Fcoll_BeR = molecular_fraction(n_h * T_h);
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION_OUTPUT
  /* fprintf(FdMolFrac,"# 1.case   2.time    3.particleID    4.H2_frac_Krumholz    5.H2_frac_Blitz&R */
  /*                   6.metallicity    7.number_dens hot(1)/cold(2)    8.temperature hot(1)/cold(2)\n"); */
  /* to compare different predictions on a particle basis, for particles with different properties */
  /* case 1 (1st column) means FORCE SINGLE PHASE */
  if(ThisTask == 0)
    {
      fprintf(FdMolFrac, "%d    %g    %d    %g    %g     %g     %g    %g \n",
	      1, All.Time, FBpart_ID, Fcoll, Fcoll_BeR, metallicity, n_h, T_h);
      fflush(FdMolFrac);
    }
#endif
#else //KRUMHOLZ not defined
  Fcoll = molecular_fraction(n_h * T_h);
#endif //ends MV_KRUMHOLZ_MOLECULAR_FRACTION


  Fcoll = molecular_fraction(n_h * T_h);

  if(isnan(Fcoll) || isinf(Fcoll) || Fcoll > 1.0 || Fcoll < SMALL)
    {
      for(j = 0; j < 4; j++)
	dydx[j] = 0.;

#ifdef GM_MUPPI_WARNINGS
      printf("WARNING: task %d  part %d got wrong Fcoll (%g) in FSP: Th=%g Eh=%g Mh=%g Rhoh=%g nh=%g\n",
	     ThisTask, FBpart_ID, Fcoll, T_h, E_h, M_h, Rho_h, n_h);
#endif

      return;

    }


#ifdef MV_UPDATE_TDYN_FOR_STARFORMATION
  compute_ism_properties(y[2], y[3], y[1], mVol_tot, &T_h, &fill_h, &Rho_h, &Rho_c, &n_h, &n_c);
  upd_tdyn = 5.15e7 / sqrt(mu_c * n_h * T_h / T_c) * YEAR_IN_CODE_UNITS;
#endif

  tdyn = SphP[FBinternal].tdyn;

  dMevap = M_c / tdyn;		/* the cold phase is evacuated */

#ifdef MV_UPDATE_TDYN_FOR_STARFORMATION
  dMsf = Fcoll * f_star * M_h / upd_tdyn;	/* star formation is computed on the hot phase */
#else
  dMsf = Fcoll * f_star * M_h / tdyn;	/* star formation is computed on the hot phase */
#endif

  dMrest = f_rest * dMsf;	/* restoration is as usual */


  /* system of equations */
  /* no cooling term, thermal energy is lost proportional to 
     loss in hot phase mass */
  dydx[0] = (double) (dMsf - dMrest);
  dydx[1] = (double) (-dMevap);
  dydx[2] = (double) (dMevap - dMsf + dMrest + dM);
  dydx[3] = (double) (dydx[2] / M_h * E_h);

  /* forcing zeroes */
  if(fabs(dydx[0]) < SMALL)
    dydx[0] = 0.0;
  if(fabs(dydx[1]) < SMALL)
    dydx[1] = 0.0;
  if(fabs(dydx[2]) < SMALL)
    dydx[2] = 0.0;
  if(fabs(dydx[3]) < SMALL)
    dydx[3] = 0.0;

  return;

}





#ifdef GM_MUPPI_DEBUG
/* DEBUG function! */
/* Checks whether one of the MUPPI variables is NAN or INF */
void check_muppi_nan(int i)
{
  unsigned int bits, id;

  //    checkline("");
  for(bits = 0; GENERATIONS > (1 << bits); bits++);
  id = P[i].ID << bits;
  id = id >> bits;

  if(isnan(SphP[i].M_c) || isnan(SphP[i].M_h) || isnan(SphP[i].E_h) ||
     isnan(SphP[i].Egy_tot_0) || isnan(SphP[i].tdyn) ||
     isnan(SphP[i].M_sf) || isnan(SphP[i].E_out) ||
     isnan(SphP[i].Eout_norm) || isnan(SphP[i].E_rec) ||
     isnan(SphP[i].Density) || isnan(SphP[i].Pressure) ||
     isnan(SphP[i].Entropy) || isnan(SphP[i].GradDens[0]) ||
     isnan(SphP[i].GradDens[1]) || isnan(SphP[i].GradDens[2]))
    {
      printf("\n\n ************ NAN FOUND ***************\n");
      printf("   Task: %d, Particle: %d, Time: %f\n", ThisTask, id, All.Time);
      printf(" MUPPI: NMF=%d, MultiPhase=%d\n", SphP[i].NMF, SphP[i].MultiPhase);
      printf(" M_c=%g M_h=%g M_sf=%g E_h=%g\n", SphP[i].M_c, SphP[i].M_h, SphP[i].M_sf, SphP[i].E_h);
      printf(" tdyn=%g\n", SphP[i].tdyn);
      printf(" Egy_tot_0=%g \n", SphP[i].Egy_tot_0);
      printf("Eout=%g Eout_norm=%g E_rec=%g\n", SphP[i].E_out, SphP[i].Eout_norm, SphP[i].E_rec);
      printf("\nTask: %d part_i: %d ID:%u\n", ThisTask, i, id);
      printf(" Mass=%g Pos= %f, %f, %f Vel= %f, %f,%f\n", P[i].Mass,
	     P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2]);
      printf("Pres=%f Dens=%g Entropy=%g dEntropy=%g\n", SphP[i].Pressure,
	     SphP[i].Density, SphP[i].Entropy, SphP[i].Entropy);
      printf("Grad rho= %g %g %g\n", SphP[i].GradDens[0], SphP[i].GradDens[1], SphP[i].GradDens[2]);
    }

  if(isinf(SphP[i].M_c) || isinf(SphP[i].M_h) || isinf(SphP[i].E_h) ||
     isinf(SphP[i].Egy_tot_0) || isinf(SphP[i].tdyn) ||
     isinf(SphP[i].M_sf) || isinf(SphP[i].E_out) ||
     isinf(SphP[i].Eout_norm) || isinf(SphP[i].E_rec) ||
     isinf(SphP[i].Density) || isinf(SphP[i].Pressure) ||
     isinf(SphP[i].Entropy) || isinf(SphP[i].GradDens[0]) ||
     isinf(SphP[i].GradDens[1]) || isinf(SphP[i].GradDens[2]) || SphP[i].Density <= 0.0)
    {
      printf("\n\n ************ INF FOUND ***************\n");
      printf("   Task: %d, Particle: %d, Time: %f\n", ThisTask, id, All.Time);
      printf(" MUPPI: NMF=%d, MultiPhase=%d\n", SphP[i].NMF, SphP[i].MultiPhase);
      printf(" M_c=%g M_h=%g M_sf=%g E_h=%g\n", SphP[i].M_c, SphP[i].M_h, SphP[i].M_sf, SphP[i].E_h);
      printf(" tdyn=%g\n", SphP[i].tdyn);
      printf(" Egy_tot_0=%g \n", SphP[i].Egy_tot_0);
      printf("Eout=%g Eout_norm=%g Erec=%g\n", SphP[i].E_out, SphP[i].Eout_norm, SphP[i].E_rec);
      printf("\nTask: %d part_i: %d ID:%d\n", ThisTask, i, id);
      printf(" Mass=%g Pos= %f, %f, %f Vel= %f, %f,%f\n", P[i].Mass,
	     P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2]);
      printf("Pres=%f Dens=%g Entropy=%g dEntropy=%g\n", SphP[i].Pressure,
	     SphP[i].Density, SphP[i].Entropy, SphP[i].Entropy);
      printf("Grad rho= %g %g %g\n", SphP[i].GradDens[0], SphP[i].GradDens[1], SphP[i].GradDens[2]);
    }
}
#endif


/* This function initializes values needed by MUPPI.
   Called from MUPPI_main */
static inline void initialize_muppi_vars(int i, MyAtLeastDouble * Rho, MyAtLeastDouble * dblMuppiVars,
					 int *intMuppiVars)
{

  int bits, id;
#ifdef LT_STELLAREVOLUTION
  MyAtLeastDouble MyEgy;
#endif


  /* zeroing MUPPI RAM */
  //  checkline("");

  mMass_tot = 0.0;
  mEgy_tot = 0.0;
  mM_h_0 = 0.0;
  mM_c_0 = 0.0;
  mx1 = 0.0;
  mx2 = 0.0;
  mhh = 0.0;
  mT_h = 0.0;
  mfill_h = 0.0;
  mRho_h = 0.0;
  mn_h = 0.0;
  mRho_c = 0.0;
  mn_c = 0.0;
  mM_sf0 = 0.0;
  mdE = 0.0;
  mdM = 0.0;
  mVol_tot = 0.0;
  mfill_c = 0.0;
  mTime = 0.0;
  mdtime = 0.0;
  mhh_old = 0.0;
  mTotMassWtRest = 0.0;

  mnstep = 0;
  mncall = 0;
  mmfin = 0;
  mmfout = 0;



  /* particle ID */
  for(bits = 0; GENERATIONS > (1 << bits); bits++);
  id = P[i].ID << bits;
  id = id >> bits;
  mpart_ID = id;
  minternal = i;



  /* Mass and energy */
  mMass_tot = P[i].Mass;
  mEgy_tot = (SphP[i].Entropy) / GAMMA_MINUS1 * pow(*Rho, GAMMA_MINUS1);	/* specific internal energy */

#ifdef LT_STELLAREVOLUTION
  /* Tot Mass including restored mass from outside */
  /* used if a non-clock exit from MP is called */
  mTotMassWtRest = SphP[i].M_h + SphP[i].M_c + SphP[i].M_sf + SphP[i].MassRes;
#endif

  /* Converting relevant quantities into h true */
  *Rho *= All.HubbleParam * All.HubbleParam;
  mdtime = dtime / All.HubbleParam;
  if(dtime <= 0.0)
    printf("AAAAAAHHHHH! dtime=%e\n", dtime);
  fflush(stdout);



  mMass_tot /= All.HubbleParam;
  mTime = All.Time / All.HubbleParam;

  SphP[i].M_h /= All.HubbleParam;
  SphP[i].M_c /= All.HubbleParam;
  SphP[i].M_sf /= All.HubbleParam;
  SphP[i].E_h /= All.HubbleParam;

  /* storing information in the FB structure */
  mVol_tot = (mMass_tot - SphP[i].M_sf) / *Rho;

#ifdef LT_STELLAREVOLUTION	/* account for energy from long-living stars */
  MyEgy = (mdtime > SphP[i].EgyStep ? SphP[i].EgyRes : SphP[i].EgyRes * mdtime / SphP[i].EgyStep);
  mEgy_tot += MyEgy / (mMass_tot - SphP[i].M_sf);
  SphP[i].EgyRes -= MyEgy;
  mdM = SphP[i].MassRes / All.HubbleParam / mdtime;
  SphP[i].MassRes = 0.0;	//this is the mass rest evaluated at the PREVIOUS time step, we set now it to zero
  //note that it INCLUDES mass res from SnIa that should be evaluated as the energy above
#else
  mdM = 0.0;
#endif

#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
  mRho = *Rho;
#endif

  /* These vars are in CODE unit, not converted */
#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
  mBHEhotOrig = (MyAtLeastDouble) SphP[i].Injected_BH_Energy * SphP[i].CouplingBHEnHot;
  mBHEcoldOrig = (MyAtLeastDouble) SphP[i].Injected_BH_Energy * SphP[i].CouplingBHEnCold;
#endif

  return;
}


/* This function initializes MUPPI at the first time step.
   Called from MUPPI_main */
static inline void start_MUPPI(int i, MyAtLeastDouble * dblMuppiVars, int *intMuppiVars)
{

  //  checkline("");

  /* -------------- Particle just entered Multi-Phase -------------- */
  /* in this case MUPPI is only initialized */
  SphP[i].NMF++;		/* number of times the particle entered the multiphase stage */
  SphP[i].M_sf = FracS * mMass_tot;	/* mass in stars */
  SphP[i].M_h = FracH * mMass_tot;	/* mass in hot phase */
  SphP[i].M_c = FracC * mMass_tot;	/* mass in cold phase */
#ifdef LT_STELLAREVOLUTION	//need to repeat to avoid zeroes at the first muppi step
  mTotMassWtRest = SphP[i].M_h + SphP[i].M_c + SphP[i].M_sf + SphP[i].MassRes;
  mTotMassWtRest *= All.HubbleParam;
#endif

  SphP[i].tdyn = 0.0;
  SphP[i].clock = 0.0;

#if defined(MV_GM_STELLAR_KIN_FB2)
  SphP[i].DelayTime = 0.0;
#endif

  SphP[i].t_startMP = 0.0;


  /* hot phase thermal energy; cold phase thermal energy is subtracted (very small!) */
  mT_h =
    SE_TO_T * mEgy_tot * (SphP[i].M_h + SphP[i].M_c) / SphP[i].M_h -
    SphP[i].M_c / SphP[i].M_h * rap_mol_weight * T_c;

  SphP[i].E_h = T_TO_SE * mT_h * SphP[i].M_h;

  mmfin = mmfin + 1;		// printf("    startMUPPI ....%d\n", mmfin); fflush(stdout);                  /* number of particles that get multiphase */
#ifdef GM_MUPPI_DEBUG
  check_muppi_nan(i);
#endif
}

/* This function checks if the gas particle must exit MUPPI.
   Called from MUPPI_main */
static inline int check_exit_muppi(int i, MyAtLeastDouble Rho, MyAtLeastDouble * dblMuppiVars,
				   int *intMuppiVars)
{

  if(Rho < All.PhysDensThresh / 1.5)	/* now Rho and All.PhysDensThresh are both in h true */
    {
      /* SphP[i].M_sf and SphP[i].E_out are not zeroed because they are used later
       * for star formation and energy distribution */

      SphP[i].MultiPhase = 0;



      SphP[i].M_h = 0.0;
      SphP[i].M_c = 0.0;
      SphP[i].E_h = 0.0;
      SphP[i].clock = 0.0;
      mmfout += 1;
      /* Transforms non-nulled quantities in h100 */
      SphP[i].M_sf *= All.HubbleParam;
      mTime *= All.HubbleParam;
      SphP[i].E_out *= All.HubbleParam;

#ifdef GM_MUPPI_DEBUG
      check_muppi_nan(i);
#endif
      return 1;
    }
  return 0;
}

/* This function is called from MUPPI_main during the integration, to check it it must be stopped.*/
static inline int muppi_integration_check_exits(int i, int status, MyAtLeastDouble * dblMuppiVars,
						int *intMuppiVars, double *y)
{
  MyAtLeastDouble tt;

  //  checkline("");
  /* checks that the integration has succeded, otherwise break */
  if(status != GSL_SUCCESS)
    {
      SphP[i].MultiPhase = -5;	/* MultiPhase=-5 is the errorcode for GSL error */
      return -1;
    }

  mT_h = SE_TO_T * y[3] / y[2];

  /* particle is freezing: break */
  if(mT_h < T_c || y[3] < 0 || y[2] < 0)
    {
      SphP[i].MultiPhase = -3;	/* MultiPhase=-3  is the errorcode for freezing */
      return -1;
    }

  /* particle is depleted of cold phase: break */
  if(SphP[i].t_startMP >= -0.1 && y[1] < SMALL && y[0] > SMALL)	/* particles where single phase is forced are excluded */
    {
      SphP[i].MultiPhase = -4;	/* MultiPhase=-4  is the errorcode for depletion of cold phase */
      return -1;
    }


  /* SF RUNAWAY! break */
  if((y[0] > mMass_tot / GENERATIONS ||
      y[2] + y[1] + y[0] < 1.2 * mMass_tot / GENERATIONS) && (y[2] + y[1]) / y[0] < 0.01)
    {
      SphP[i].MultiPhase = -10;
      return -1;
    }

  /*  old EXIT_MP_FOR_WIND_PARTICLES, now default */
  /* Particle in wind? break */
  /* Particles in wind (fountain) have accumulated significant stellar mass */
  /*    but have a high (>0.1) hot gas fraction */
  if(y[2] > 0.1 * y[1] && y[0] > 0.03 * y[1])
    {
      SphP[i].MultiPhase = -7;
      return -1;
    }


  /* too many steps? break */
  if(mnstep > 1000)
    {
      SphP[i].MultiPhase = -8;
      return -1;
    }

  if(SphP[i].tdyn == 0.0)
    if(y[1] > NC_COLDMASS_TH * y[2])
      {
	compute_ism_properties(y[2], y[3], y[1], mVol_tot, &mT_h, &mfill_h, &mRho_h, &mRho_c, &mn_h, &mn_c);

	SphP[i].tdyn = 5.15e7 / sqrt(mu_c * mn_c) * YEAR_IN_CODE_UNITS;
	SphP[i].clock = CLOCK * SphP[i].tdyn;
      }

  if(SphP[i].t_startMP >= -0.1)
    {
      SphP[i].t_startMP += mdtime;

      compute_ism_properties(y[2], y[3], y[1], mVol_tot, &mT_h, &mfill_h, &mRho_h, &mRho_c, &mn_h, &mn_c);
      tt = 5.15e7 / sqrt(mu_c * mn_h * mT_h / T_c) * YEAR_IN_CODE_UNITS;

      /* if multi-phase is not set after a dynamical time */
      /* and if pressure is higher than PRESFMOL */
      /* then trigger a maximal starburst */
      if(SphP[i].tdyn == 0.0 && mT_h < 12000. && mT_h * mn_h > PRESFMOL && SphP[i].t_startMP > tt)
	{

	  SphP[i].tdyn = tt;
	  SphP[i].clock = CLOCK * SphP[i].tdyn;

	  SphP[i].t_startMP = -1.0;

	}
    }
  return 0;

}


/* This function is called from MUPPI_main, to complete the exit from MP
   if such a condition arised during integration */
static inline void muppi_exec_exits(int i, MyAtLeastDouble * dblMuppiVars, int *intMuppiVars, double *y)
{
  //  checkline("");
  /* exit from multi-phase */
  /* SphP[i].M_sf and SphP[i].E_out are zeroed because SFR and energy redistribution are switched off */
  if(SphP[i].MultiPhase != -10)	/* if the particle exits for star formation runaway then do not ignore these stars */
    SphP[i].M_sf = 0.0;

  SphP[i].M_h = 0.0;
  SphP[i].M_c = 0.0;
  SphP[i].E_h = 0.0;
  SphP[i].E_out = 0.0;

  SphP[i].E_kin = 0.0;
  SphP[i].xkin = 0.0;
  SphP[i].ykin = 0.0;
  SphP[i].zkin = 0.0;

#ifdef LT_STELLAREVOLUTION
  P[i].Mass = mTotMassWtRest;	/* to avoid mass non-conservations */
  SphP[i].mstar = 0.0;
#endif
  SphP[i].clock = 0.0;
  //  mmfout += 1; /* number of parts that exit MF */

#ifdef GM_MUPPI_WARNINGS
  /* warnings for remarkable traumatic exits (not for freezing) */
  if(SphP[i].MultiPhase == -4)
    printf
      ("WARNING: Task %d particle %d t=[%f %f] - DEPLETED OF COLD PHASE, fill_c = %g y=[%g %g %g %g]...\n",
       ThisTask, mpart_ID, mx1, mx2, mfill_c, y[0], y[1], y[2], y[3]);
  if(SphP[i].MultiPhase == -5)
    printf("WARNING: Task %d particle %d t=[%f %f] - GSL ERROR...\n", ThisTask, mpart_ID, mx1, mx2);
  if(SphP[i].MultiPhase == -6)
    printf("WARNING: Task %d particle %d t=[%f %f] - MASS NOT CONSERVED y=[%g %g %g %g], Mgas=%g...\n",
	   ThisTask, mpart_ID, mx1, mx2, y[0], y[1], y[2], y[3], mMass_tot);
  if(SphP[i].MultiPhase == -8)
    printf
      ("WARNING: Task %d particle %d t=[%f %f] - INTEGRATION TOO LONG y=[%g %g %g %g], initial hh=%g  hh=%g...\n",
       ThisTask, mpart_ID, mx1, mx2, y[0], y[1], y[2], y[3], mhh_old, mhh);
  if(SphP[i].MultiPhase == -10)
    printf("WARNING: Task %d particle %d t=[%f %f] - STARFORMATION RUNAWAY y=[%g %g %g %g]...\n", ThisTask,
	   mpart_ID, mx1, mx2, y[0], y[1], y[2], y[3]);
  fflush(stdout);
#endif


#ifdef GM_MUPPI_DEBUG
  check_muppi_nan(i);
#endif

  /* Transforms non-nulled quantities in h100 */
  SphP[i].M_sf *= All.HubbleParam;
  mTime *= All.HubbleParam;

  mmfout = mmfout + 1;
  return;
}

/* This function exec the final operations on particle when the integration is successfully finished.
   Called from MUPPI_main */
static inline void muppi_final_step_operations(int i, MyAtLeastDouble Rho, MyAtLeastDouble * dblMuppiVars,
					       int *intMuppiVars, double *y)
{
  MyAtLeastDouble Rho_h, Rho_c, T_h, fill_h, T_med, mu_med, n_h, n_c;


  //  checkline("");
  /* assigns variables again */
  SphP[i].M_sf = (MyAtLeastDouble) y[0];
  SphP[i].M_c = (MyAtLeastDouble) y[1];
  SphP[i].M_h = (MyAtLeastDouble) y[2];
  SphP[i].E_h = (MyAtLeastDouble) y[3];


#ifdef MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z
  int j;
  MyAtLeastDouble metallicity_esn, metal_mass_esn;
  metal_mass_esn = 0.0;
  for(j = 1; j < LT_NMetP; j++)
    metal_mass_esn = metal_mass_esn + SphP[i].Metals[j];
  metallicity_esn = metal_mass_esn / P[i].Mass;
  metallicity_esn /= SOLAR_METALLICITY;
  if(metallicity_esn < MIN_METALLICITY)
    {
      SphP[i].E_out =
	MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z * E_SN_51 * (SphP[i].M_sf - mM_sf0) / (1.0 - f_rest) / beta_sf;
      /* XXXX
         printf("Task %d E_out %e M_sf %e M_sf0 %e\n", ThisTask, SphP[i].E_out, SphP[i].M_sf, mM_sf0); fflush(stdout);
       */
#ifdef MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z_OUTPUT
      fprintf(FdHypernovae,
	      "%f       %d        %f  %f  %f         %e        %e         %g         %g        %g      \n",
	      All.Time, mpart_ID, P[i].Pos[0], P[i].Pos[0], P[i].Pos[2], SphP[i].Density, SphP[i].Temperature,
	      P[i].Mass, metallicity_esn, SphP[i].E_out);
      fflush(FdHypernovae);
#endif
    }
  else
    SphP[i].E_out = E_SN_51 * (SphP[i].M_sf - mM_sf0) / (1.0 - f_rest) / beta_sf;
#else
  SphP[i].E_out = E_SN_51 * (SphP[i].M_sf - mM_sf0) / (1.0 - f_rest) / beta_sf;
#endif


#ifdef MV_GM_AGNMUPPI		//GM-chk qui diamo il surplus di energia, tutto assieme
  if(mInitialBHEcold - y[4] > 1.e-6)
    SphP[i].E_h += mInitialBHEcold - y[4];
#endif


  compute_ism_properties(SphP[i].M_h, SphP[i].E_h, SphP[i].M_c, mVol_tot, &T_h, &fill_h, &Rho_h, &Rho_c, &n_h,
			 &n_c);


  /* average (mass-weighted) molecular weight and temperature  */
  mu_med = (mu_c * SphP[i].M_c + mu_h * SphP[i].M_h) / (SphP[i].M_c + SphP[i].M_h);
  T_med = mu_med / (SphP[i].M_h + SphP[i].M_c) * ((SphP[i].M_h * T_h) / mu_h + (SphP[i].M_c * T_c) / mu_c);

  mEgy_tot = BOLTZMANN * T_med / GAMMA_MINUS1 / mu_med / PROTONMASS;
  mEgy_tot *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;	/* non dipende da h */

#ifdef MV_GM_AGNMUPPI

#if defined (MV_GM_AGNMUPPI) && !defined(MV_AGNMUPPI_COOLING_OFF_OUTPUT)
  if((mBHEhotOrig > 0.) || (mBHEcoldOrig > 0.))
    {
      // remember, mtcool is the initial one  
      fprintf(FdAGNmuppi, "%g %e    %e %e      %e %e %e      %e %e %e      %e %e      %e %e       %e %e\n",
	      All.Time, mdtime,
	      mBHEhotOrig, mBHEcoldOrig,
	      mM_c_0, mM_h_0, mT_h_0,
	      SphP[i].M_c, SphP[i].M_h, T_h,
	      mInitialBHEcold, y[4],
	      mtcool, SphP[i].Density,
	      SphP[i].Pressure * All.UnitPressure_in_cgs / BOLTZMANN * pow(All.HubbleParam, 3), n_h * T_h);
      fflush(FdAGNmuppi);	//milena   
    }
#endif

#ifdef MV_AGNMUPPI_COOLING_OFF
  mt_startCoolOff += mdtime;	//mile3: questo e' il tempo trascorso dall'ultima volta in cui la fase cold della
#endif //particella MP considerata ha ricevuto energia dall'AGN#endif

#ifdef MV_AGNMUPPI_COOLING_OFF_OUTPUT
  if((mBHEhotOrig > 0.) || (mBHEcoldOrig > 0.))
    {
      // remember, mtcool is the initial one
      fprintf(FdCoolOffAGNmuppi,
	      "%g   %e %e     %e %e      %e %e %e      %e %e %e      %e %e      %e %e       %e %e\n",
	      All.Time, mdtime, mt_startCoolOff, mBHEhotOrig, mBHEcoldOrig, mM_c_0, mM_h_0, mT_h_0,
	      SphP[i].M_c, SphP[i].M_h, T_h, mInitialBHEcold, y[4], mtcool, SphP[i].Density,
	      SphP[i].Pressure * All.UnitPressure_in_cgs / BOLTZMANN * pow(All.HubbleParam, 3), n_h * T_h);
      fflush(FdCoolOffAGNmuppi);	//milena
    }
#endif


#endif //ifdef MV_GM_AGNMUPPI


  /* This is the (non-specific) energy after MUPPI - 
     every change in next timestep will be due to hydrodynamics
     (including changes in entropy and volume) and will be taken into account in the dE term */
  SphP[i].Egy_tot_0 = mEgy_tot * (SphP[i].M_h + SphP[i].M_c);
  SphP[i].Entropy = 0.0;

  /* Specific entropy is given in h100 */
  SphP[i].Entropy = mEgy_tot * GAMMA_MINUS1 /
    pow(Rho / All.HubbleParam / All.HubbleParam * (SphP[i].M_h + SphP[i].M_c) / (mM_h_0 + mM_c_0),
	GAMMA_MINUS1);


  SphP[i].EntropyPred = SphP[i].Entropy;
  SphP[i].Pressure = get_pressure(i);



#ifdef LT_STELLAREVOLUTION
  /* in this case we must reset Mass to account for 
     changes (increase due to FB.dM, and the decrease
     of M_sf evaluated by stellar evolution */
  P[i].Mass = (SphP[i].M_c + SphP[i].M_h + SphP[i].M_sf) * All.HubbleParam;
  if(P[i].Mass <= 0.0)
    {
      printf(" WARNING, negative mass in MUPPI!\n");
      fflush(stdout);
    }
#endif

  /* exit MP after CLOCK * tdyn */
  if(SphP[i].clock > 0.0)
    {
      SphP[i].clock -= mdtime;
      if(SphP[i].clock < 0.0)
	{
	  /* exit from multi-phase */
	  /* SphP[i].M_sf and SphP[i].E_out are not zeroed to have SFR and distribution of energy */
	  SphP[i].MultiPhase = -1;	/* MultiPhase=-1 implies regular exit by clock */

	  if(SphP[i].t_startMP < 0.0)
	    SphP[i].MultiPhase = -2;	/* MultiPhase=-2 implies regular exit during forced single phase regime */

	  SphP[i].M_h = 0.0;
	  SphP[i].M_c = 0.0;
	  SphP[i].E_h = 0.0;
	  SphP[i].clock = 0.0;
	  mmfout = mmfout + 1;	/* number of particles that go out of MP */
	}
    }


  /* converting to h100 */
  SphP[i].M_sf *= All.HubbleParam;
  SphP[i].M_c *= All.HubbleParam;
  SphP[i].M_h *= All.HubbleParam;
  SphP[i].E_h *= All.HubbleParam;
  SphP[i].E_out *= All.HubbleParam;
  mTime *= All.HubbleParam;
  //  { int mff=0, moo=0;
  //    CheckMuppiMemory(i, Rho, &mff, &moo, intMuppiVars, dblMuppiVars);
  //  }
}

/* This function initializes MUPPI integration
   Called from MUPPI_main */
static inline void muppi_initialize_integration(int i, MyAtLeastDouble Rho, MyAtLeastDouble * dblMuppiVars,
						int *intMuppiVars, double *y, struct myparams *param)
{
  MyAtLeastDouble tstar_guess, tcool;

  //  checkline("");
  /* number of time steps in multiphase */
  SphP[i].MultiPhase++;

  mM_sf0 = SphP[i].M_sf;
  mM_h_0 = SphP[i].M_h;
  mM_c_0 = SphP[i].M_c;

  /* first timestep, advance Th to S at endstep and no dS */
  if(SphP[i].MultiPhase == 1)
    mdE = 0.0;
  else
    {				/* next timesteps, dS advances S to endstep */
      mdE = (mEgy_tot * (mMass_tot - SphP[i].M_sf) - SphP[i].Egy_tot_0) / mdtime;
    }

#ifdef MV_GM_AGNMUPPI

#ifdef MV_GM_AGNMUPPI_OUTPUT
  if(mBHEcoldOrig > 0. || mBHEhotOrig > 0.)
    {
      E_AGNM_sfr_dE += mdE * mdtime;	//milena, valore prima dell'incremento
      n_AGNM++;
    }
#endif


  mBHEhot = mBHEhotOrig * All.HubbleParam;
  mInitialBHEcold = mBHEcoldOrig * All.HubbleParam;	//GM-chk ricordiamo l'en cold prima dei cambi di UdM
  mdE += mBHEhot / mdtime;	//BH energy coupled with the hot phase
  mBHEcoldD = mInitialBHEcold / mdtime * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s;
  /* the unit conversion is explained in derivatives, search dMBH */
  // perche' alla riga 2356 dBHEcold=FB.BHEcoldD viene diviso per un'energia che e' in cgs
#ifdef MV_AGNMUPPI_COOLING_OFF
  if(mBHEcoldOrig > 0.)
    mt_startCoolOff = 0.;	//mile3: definisco questo tempo NON negativo solo per particelle MP che ricevono
  //AGN feedback energy alla componente cold
  else
    mt_startCoolOff = -99.;
#endif

#ifdef MV_GM_AGNMUPPI_OUTPUT
  if(mBHEcoldOrig > 0. || mBHEhotOrig > 0.)
    {
      E_AGNM_sfr_dE_aft += mdE * mdtime;	//milena, valore dopo l'incremento
      E_AGNM_sfr_BHEcoldD += mBHEcoldOrig;	//milena, valore prima dell'incremento
    }
#endif


#endif //MV_GM_AGNMUPPI


  /* particle ISM */
  compute_ism_properties(SphP[i].M_h, SphP[i].E_h, SphP[i].M_c, mVol_tot, &mT_h, &mfill_h, &mRho_h, &mRho_c,
			 &mn_h, &mn_c);
  mfill_c = 1.0 - mfill_h;


  /* computation of cooling time */
  tcool = cooling_time(mT_h, SphP[i].E_h / SphP[i].M_h, mRho_h, mn_h);
  //  checkline("back to muppi_initialize_integration");

  /* integration interval */
  mx1 = 0.;
  mx2 = mdtime;

  if(SphP[i].tdyn > 0.0)
    tstar_guess = SphP[i].tdyn / f_star;
  else
    tstar_guess = 5.15e7 / sqrt(Rho * All.UnitDensity_in_cgs / PROTONMASS) * YEAR_IN_CODE_UNITS / 0.1;	/* all in h true */

  /* first integration time-step */
  if(tcool > SMALL)
    {
      mhh = tcool / SCALE_HH;
      if(mhh > tstar_guess / SCALE_HH)
	mhh = tstar_guess / SCALE_HH;
    }
  else
    mhh = tstar_guess / SCALE_HH;

  if(mhh > 0.25 * mdtime)
    mhh = 0.25 * mdtime;
  if(mhh < 0.001 * mdtime)
    mhh = 0.001 * mdtime;

  /* initial conditions for integration */

  y[0] = (double) SphP[i].M_sf;
  y[1] = (double) SphP[i].M_c;
  y[2] = (double) SphP[i].M_h;
  y[3] = (double) SphP[i].E_h;
#ifdef MV_GM_AGNMUPPI
  y[4] = 0.0;
#endif


  /* derivative function parameters */
  param->dE = mdE;
  param->dM = mdM;
  param->VolTot = mVol_tot;
  param->internal = minternal;
  param->part_ID = mpart_ID;
  param->filtered = 0;
  param->fill_c = mfill_c;
  param->ncall = 0;
  param->dblsmuppi = dblMuppiVars;
  param->intsmuppi = intMuppiVars;

}


/* compute the properties of multi-phase ISM in pressure equilibrium */
/* the quantities M_h, M_c, E_h and Vol_tot are assumed to be given in code units */
/* Rho_h and Rho_c are returned in code units, T_h, n_h and n_c in cgs */
/* works in h true */
//static inline void compute_ism_properties(MyAtLeastDouble M_h, MyAtLeastDouble E_h, MyAtLeastDouble M_c, MyAtLeastDouble Vol_tot,
void compute_ism_properties(MyAtLeastDouble M_h, MyAtLeastDouble E_h, MyAtLeastDouble M_c,
			    MyAtLeastDouble Vol_tot, MyAtLeastDouble * T_h, MyAtLeastDouble * fill_h,
			    MyAtLeastDouble * Rho_h, MyAtLeastDouble * Rho_c, MyAtLeastDouble * n_h,
			    MyAtLeastDouble * n_c)
{
  MyAtLeastDouble Vol_h, Vol_c, Frac_h, Frac_c;


  //  checkline("");
  *T_h = SE_TO_T * E_h / M_h;
  Frac_h = M_h / (M_c + M_h);
  Frac_c = 1.0 - Frac_h;
  *fill_h = 1.0 / (1.0 + Frac_c / Frac_h * mu_h / mu_c * T_c / *T_h);
  Vol_h = Vol_tot * *fill_h;
  Vol_c = Vol_tot * (1.0 - *fill_h);
  *Rho_h = M_h / Vol_h;
  *n_h = *Rho_h * All.UnitDensity_in_cgs / mu_h / PROTONMASS;
  if(*fill_h < 1.0)
    {
      *Rho_c = M_c / Vol_c;
      *n_c = *Rho_c * All.UnitDensity_in_cgs / mu_c / PROTONMASS;
    }
  else
    {
      *Rho_c = 0.0;
      *n_c = 0.0;
    }
  return;
}

/* compute cooling time */
/* sE_h, Rho_h and cooling_time are in code units, T_h and n_h in cgs */
/* works in h true */
static inline MyAtLeastDouble cooling_time(MyAtLeastDouble T_h, MyAtLeastDouble sE_h, MyAtLeastDouble Rho_h,
					   MyAtLeastDouble n_h)
{
  double ne;
  MyAtLeastDouble tcool, Rho_cgs, sE_cgs;
  MyAtLeastDouble Temperature;

  //  checkline("");
  Rho_cgs = Rho_h * All.UnitDensity_in_cgs;
  sE_cgs = sE_h * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
  ne = (double) n_h;
  if(T_h > TMAXCOOL)
    Temperature = TMAXCOOL;
  else
    Temperature = T_h;

#if !defined(LT_METAL_COOLING) && !defined(LT_METAL_COOLING_WAL)
  tcool = get_cooling_time_fromT((double) Temperature, (double) sE_cgs, (double) Rho_cgs, &ne);
#else
#if defined(LT_METAL_COOLING)
  tcool = get_cooling_time_fromT((double) Temperature, (double) sE_cgs, (double) Rho_cgs, &ne, Zcool);
#endif

#if defined(LT_METAL_COOLING_WAL)
#ifdef GL_DUST_COOLING
  tcool =
    get_cooling_time_fromT((double) Temperature, (double) sE_cgs, (double) Rho_cgs, &ne, mRedshift, mDZ, DL,
			   DS, mMetallicities);
#else
  tcool =
    get_cooling_time_fromT((double) Temperature, (double) sE_cgs, (double) Rho_cgs, &ne, mRedshift, mDZ,
			   mMetallicities);
#endif // GL_DUST_COOLING
#endif //LT_METAL_COOLING_WAL

#endif // !defined(LT_METAL_COOLING) && !defined(LT_METAL_COOLING_WAL)

  tcool /= All.UnitTime_in_s;
  //  checkline("Exiting");
  return tcool;
}


/* molecular fraction, see Krumholz, McKee, Tumlinson (2009, 699) */

#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
static inline MyAtLeastDouble molecular_fraction_Krum(MyAtLeastDouble metallicity, MyAtLeastDouble sdens)
{
  MyAtLeastDouble Fcoll_Krum, ss, chi, delta;

  /* metallicity is needed in solar units */
  metallicity /= SOLAR_METALLICITY;
  if(metallicity < MIN_METALLICITY)
#ifndef MV_EARLY_FB_HIGH_MOL_FRACTION
    metallicity = MIN_METALLICITY;
#else
    metallicity = metallicity;
#endif

  /* surface density is needed in Msun pc^-2 */
  sdens *= All.UnitSurfDens_in_Msunpc2;

  chi = 0.77 * (1 + 3.1 * pow(metallicity, 0.365));
  ss = log(1 + 0.6 * chi) / (0.04 * sdens * metallicity);
  delta = 0.0712 * pow(0.1 / ss + 0.675, -2.8);

  if(metallicity < MIN_METALLICITY)	/* this relation should be never satisfied if MV_EARLY_FB_HIGH_MOL_FRACTION is off */
#ifndef MV_EARLY_FB_HIGH_MOL_FRACTION
    Fcoll_Krum = 1 - pow(1 + pow(3. / 4. * ss / (1 + delta), -5.), -1. / 5);
#else
    Fcoll_Krum = MV_EARLY_FB_HIGH_MOL_FRACTION;
#endif
  else
    Fcoll_Krum = 1 - pow(1 + pow(3. / 4. * ss / (1 + delta), -5.), -1. / 5);

  return Fcoll_Krum;
}
#endif



/* molecular fraction, see Blitz & Rosolowski (2006) */
/* pressure must be physical (no h100) */
static inline MyAtLeastDouble molecular_fraction(MyAtLeastDouble pressure)
{
  MyAtLeastDouble Fcoll;

  if(pressure <= 0.0)
    return 0.0;
  Fcoll = 1. / (1. + (PRESFMOL / pressure));

  return Fcoll;
}


/* This function is only here for debugging purposes, please don't use it unless
   you know what you are doing */
void checkmuppimemory(int i, MyAtLeastDouble Rho, int *mfin, int *mfout, MyAtLeastDouble * dblMuppiVars,
		      int *intMuppiVars, const char *func, const char *file, int line)
{

  FILE *f;
  char bfr[1024];

  /*
     if((unsigned)P[i].ID == 1000045260) 
     {
     printf("Particle reached! %d\n", i); fflush(stdout);
     }
   */

  sprintf(bfr, "particles/part%d.dat", (unsigned int) P[i].ID);
  f = fopen(bfr, "a");

  fprintf(f, "%f  %e %e %e %e  %e  %e %e %e %e %d %d\n", All.Time, SphP[i].E_h, SphP[i].M_h, SphP[i].M_c,
	  SphP[i].M_sf, Rho, mdE, mdM, mMass_tot, mEgy_tot, SphP[i].NMF, SphP[i].MultiPhase);
  fflush(f);

  fclose(f);

  sprintf(bfr, "particles/egypart%d.dat", (unsigned int) P[i].ID);
  f = fopen(bfr, "a");

  fprintf(f, "%f  %d   %e %e   %e %e %e   %e\n", All.Time, SphP[i].out_part, SphP[i].Eout_norm,
	  SphP[i].Ekin_norm, SphP[i].Ekin_rec[0], SphP[i].Ekin_rec[1], SphP[i].Ekin_rec[2],
	  SphP[i].Ekin_norm_fountain);
  fflush(stdout);


  fclose(f);


  return;
}

/* debugging funcion for print debugging */
inline void check_line(const char *text, const char *func, const char *file, int line)
{
  //  printf("Task=%d: %s reached: %s()/%s/line %d\n", ThisTask, text, func, file, line); fflush(stdout);

}
#endif /* closes GM_MUPPI */


#endif /* closes SFR */
#endif /* closes COOLING */
