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
#include <ctype.h>
#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"
#include "../../Gravity/forcetree.h"

#if defined(LT_SEv_INFO) || defined(LT_ZAGE)
#ifdef LT_STELLAREVOLUTION
#ifdef LT_SEv_INFO
double mZ[LT_NMet], tot_mZ[LT_NMet];

#define STAT_INUM  6
#define ERatio     0
#define EFRatio    1
#define XRatio     2
#define NoREFratio 3
#define OEgy       4
#define OEgy_I     5

#define REScount   6
#define NoREScount 7

double StatSum[STAT_INUM + 2], tot_StatSum[STAT_INUM + 2];

double StatMin[STAT_INUM], tot_StatMin[STAT_INUM];

double StatMax[STAT_INUM], tot_StatMax[STAT_INUM];

#ifdef LT_EXTEGY_INFO
static int EEInfo_grain;
#endif

#endif

#ifdef LT_TRACK_CONTRIBUTES
thread__ Contrib contrib;
thread__ float contrib_metals[LT_NMetP];
thread__ float IIcontrib[LT_NMetP], Iacontrib[LT_NMetP], AGBcontrib[LT_NMetP];

#define NULL_CONTRIB(c)  {bzero((void*)c, sizeof(Contrib));}
#define NULL_EXPCONTRIB(c)  {bzero((void*)c, sizeof(float)*LT_NMetP);}
#endif

#ifdef LT_METAL_COOLING_WAL
double Metallicities[LT_NMet], DZ, Redshift;
double REDSHIFT_DEFAULT = 0.0;
#endif
/*
 * This routine does cooling and star formation for
 * the effective multi-phase model.
 */


#ifndef GM_MUPPI
void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i, j, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  double Zcool;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, mass_of_star;
  double Sum_sm, Total_sm, sm, rate, Sum_mass_stars, Total_sum_mass_stars;
  double p, prob;
  double cloudmass;
  double factorEVP;
  double tsfr, trelax;
  double egyhot, egyeff, egycurrent, tcool, x, y, rate_in_msunperyear;
  double sfrrate, totsfrrate, dmax1, dmax2;
  double NonMetalMass;
  double myFactorEVP, myPhysDensThresh, myFactorSN, myEgySpecSN, myMaxSfrTimescale;
  double SNEgy;
  double mstar, factor;
  MyIDType GENERATIONS_BUNCH;
  int Yset, YZbin;
  double Temperature;

  double rhomean;

#if defined(LT_EXTEGY_INFO) || defined(LT_SEv_INFO)
  double tstart, tend;
#endif
  double exp_minus_p;
  double NextChemTime;
  int index_of_star, chem_step;

#ifdef WINDS
  double v;

  double norm, dir[3];

#ifdef ISOTROPICWINDS
  double theta, phi;
#endif

#ifdef LT_HOT_WINDS
  FLOAT u_hotwinds;
#endif

#ifdef LT_TRACK_WINDS
  double ne_guess, tw_temp;

  char trackwinds_flag;
#endif

#ifdef LT_SEv_INFO
  int windn, tot_windn;

  double windv_min, windv_max, windv, tot_windv_min, tot_windv_max, tot_windv;
#endif
#endif


#ifdef LT_EXTEGY_INFO
  double r;
#endif

#if defined(LT_TEMP_THRESH_FOR_MULTIPHASE) || defined(LT_TRACK_WINDS)
  double ne_guess;
#endif


#if defined(QUICK_LYALPHA) || defined(BH_THERMALFEEDBACK) || defined (BH_KINETICFEEDBACK) || defined(MODIFIED_SFR)
  double temp;
#endif
#if defined(QUICK_LYALPHA) || defined(BH_THERMALFEEDBACK) || defined (BH_KINETICFEEDBACK) || defined(MODIFIED_SFR) || defined(LT_STELLAREOLUTION)
  double u_to_temp_fac;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

#ifdef MODIFIED_SFR
  double SFRTempThresh;

  SFRTempThresh = 5.0e5 / u_to_temp_fac;
#endif

#ifdef FLTROUNDOFFREDUCTION
#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      SphP[i].i.Injected_BH_Energy = FLT(SphP[i].i.dInjected_BH_Energy);
#endif
#endif

  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinActive[bin])
      TimeBinSfr[bin] = 0;

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
      rhomean = 3 * hubble_a * hubble_a / 8 / M_PI / All.G;	/* mean density at time a in unit's code */
    }
  else
    {
      a3inv = ascale = time_hubble_a = 1;
      rhomean = 1;
    }

/* XXXX test & check GM */
#ifdef LT_METAL_COOLING_WAL
  if(All.ComovingIntegrationOn)
    Redshift = 1.0 / ascale - 1;
  else
    Redshift = REDSHIFT_DEFAULT;

  get_cool_redshift(Redshift, &DZ);
  if(LT_NMet < get_cool_n_el())
    endrun(1009);
  WalCool_tables_load(Redshift);
#endif

  for(i = 0; i < SFs_dim; i++)
    {
      sum_sm[i] = sum_mass_stars[i] = 0;
      total_sm[i] = total_sum_mass_stars[i] = 0;
      Sum_sm = Sum_mass_stars = Total_sm = Total_sum_mass_stars = 0;
      sfrrates[i] = totsfrrates[i] = 0;
    }
  for(i = 0; i < NBINS_SFR_by_DENSITY; i++)
    sfrrates_bydensity[i] = totsfrrates_bydensity[i] = 0;
  for(i = 0; i < NBINS_SFR_by_Z; i++)
    sfrrates_byZ[i] = totsfrrates_byZ[i] = 0;
  stars_spawned = stars_converted = 0;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      int mySFi, myIMFi;

#ifdef BLACK_HOLES
      if(P[i].Type == 0 && P[i].Mass > 0)
#else
      if(P[i].Type == 0)
#endif
	{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

#ifndef LT_LOCAL_IRA
	  SphP[i].mstar = 0;
#endif

#ifdef LT_METAL_COOLING_WAL
	  //      char buffer[1000];

	  set_metallicities(i, &Metallicities[0], a3inv);	/* note: sets the physical density too! */
	  Zcool = get_metallicity_solarunits(get_metallicity(i, Iron));

/* 	  if(All.NumCurrentTiStep > 0) */
/* 	    printf("\t cooling: metallicities : %d ** %d\n", ThisTask, All.NumCurrentTiStep); fflush(stdout); */

/* 	  if(P[i].ID == DBG_ID) */
/* 	    { */
/* 	      sprintf(buffer, " metallicities [%d][%u] %d: ", ThisTask, P[i].ID, get_cool_n_el()); */
/* 	      for(j = 0; j < get_cool_n_el(); j++) */
/* 		sprintf(&buffer[strlen(buffer)], " (%d) %e ", j, Metallicities[j]); fflush(stdout); */
/* 	      printf("%s\n", buffer);fflush(stdout); */
/* 	    } */
#else
#ifndef LT_METAL_COOLING_on_SMOOTH_Z
	  Zcool = get_metallicity_solarunits(get_metallicity(i, Iron));
#else
	  double metalmass;
	  metalmass = get_metalmass(SphP[i].Metals);
	  if(metalmass > 0)
	    {
#ifdef LT_ZSMOOTH_ALLMETALS
	      Zcool = get_metallicity_solarunits(SphP[i].Zsmooth[Iron]);
#else
	      Zcool = get_metallicity_solarunits(SphP[i].Zsmooth * SphP[i].Metals[Iron] / metalmass);
#endif
	    }
	  else
	    {
	      Zcool = NO_METAL;
	    }
#endif
#endif

/* #if defined (UM_METAL_COOLING) */
/* 	  um_ZsPoint = SphP[i].Metals; */
/* 	  um_mass = P[i].Mass; */
/* /\*         um_FillEl_mu = SphP[i].FillEl_mu; *\/ */
/* #endif */

	  get_SF_index(i, &mySFi, &myIMFi);

	  GENERATIONS_BUNCH = All.Generations / SFs[mySFi].Generations;

	  Yset = IMFs[myIMFi].YSet;
	  myFactorSN = SFs[mySFi].FactorSN;
	  myEgySpecSN = SFs[mySFi].EgySpecSN;
	  myMaxSfrTimescale = SFs[mySFi].MaxSfrTimescale;
	  for(YZbin = IIZbins_dim[Yset] - 1; Zcool < IIZbins[Yset][YZbin] && YZbin > 0; YZbin--)
	    ;

	  if(SFs[mySFi].SFTh_Zdep)
	    {
	      /*
	         if the effective model is allowed to depend on metallicity,
	         thresholds and evaporation factors differ for different Z.
	       */
	      getindex(&CoolZvalue[0], 0, ZBins - 1, &Zcool, &flag);

	      if(flag == 0 || flag == ZBins - 1)
		{
		  myFactorEVP = SFs[mySFi].FEVP[flag];
		  myPhysDensThresh = SFs[mySFi].PhysDensThresh[flag];
		}
	      else
		{
		  p = (Zcool - CoolZvalue[flag + 1]) / (CoolZvalue[flag + 1] - CoolZvalue[flag]);
		  /* interpolate */
		  myFactorEVP = SFs[mySFi].FEVP[flag] * (1 - p) + SFs[mySFi].FEVP[flag + 1] * p;
		  myPhysDensThresh =
		    SFs[mySFi].PhysDensThresh[flag] * (1 - p) + SFs[mySFi].PhysDensThresh[flag + 1] * p;
		}
	    }
	  else
	    {
	      myFactorEVP = SFs[mySFi].FEVP[0];
	      myPhysDensThresh = SFs[mySFi].PhysDensThresh[0];
	    }

	  /* check whether conditions for star formation are fulfilled.
	   *
	   * f=1  normal cooling
	   * f=0  star formation
	   */

	  /* default is normal cooling */
	  flag = 1;

	  double fac_entr_to_u = pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
	  double uold = DMAX(All.MinEgySpec, SphP[i].Entropy * fac_entr_to_u);

#ifndef MODIFIED_SFR
	  if(SphP[i].d.Density * a3inv >= myPhysDensThresh)
	    {
	      flag = 0;
#ifdef GM_BH_EVAPORATE_COLD_CLOUDS
	      factorEVP = pow(SphP[i].d.Density * a3inv / myPhysDensThresh, -0.8) * myFactorEVP;
	      egyhot = myEgySpecSN / (1 + factorEVP) + All.EgySpecCold;
	      tsfr = sqrt(myPhysDensThresh / (SphP[i].d.Density * a3inv)) * myMaxSfrTimescale;
	      tcool =
		GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, Redshift, DZ, &Metallicities[0],
			       &Temperature);
	      y = tsfr / tcool * egyhot / (myFactorSN * myEgySpecSN - (1 - myFactorSN) * All.EgySpecCold);
	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      double uc_new;
	      /* All BH energy goes to heat cold gas, can be changed... */
	      uc_new = All.EgySpecCold * x + SphP[i].i.Injected_BH_Energy / (x * P[i].Mass);

	      /* if specific energy of COLD gas is larger than average specific energy of particle,
	         all cold gas has been evaporated and particle is not MP anymore */
	      if(uc_new > uold)
		flag = 1;

	      // GM: I set a temp thresh to avoid that heated gas re-enters mp the timestep after being heated
#if !defined(LT_METAL_COOLING_WAL)
	      Temperature = convert_u_to_tempSTD(uold * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
						 SphP[i].d.Density * All.UnitDensity_in_cgs *
						 All.HubbleParam * All.HubbleParam * a3inv, &ne_guess);
#else
	      Temperature = convert_u_to_tempMET(uold * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
						 Redshift, DZ, &Metallicities[0]);
#endif

	      if(Temperature > 250000.0)	//value taken from phase diagrams of one test cluster
		flag = 1;
#endif


	    }




#else
	  if((SphP[i].d.Density * a3inv >= myPhysDensThresh) && (uold < SFRTempThresh))
	    flag = 0;
#endif

	  x = 0;

#ifdef LT_TEMP_THRESH_FOR_MULTIPHASE
	  ignore_failure_in_convert_u = 1;
#if !defined(LT_METAL_COOLING_WAL)
	  Temperature = convert_u_to_tempSTD(uold * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
					     SphP[i].d.Density * All.UnitDensity_in_cgs * All.HubbleParam *
					     All.HubbleParam * a3inv, &ne_guess);
#else
	  Temperature = convert_u_to_tempMET(uold * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
					     Redshift, DZ, &Metallicities[0]);
#endif
	  ignore_failure_in_convert_u = 0;

	  if(Temperature > All.MultiPhaseTempThresh)
	    flag = 1;
#endif


	  if(All.ComovingIntegrationOn)
	    if(SphP[i].d.Density < All.OverDensThresh)
	      flag = 1;

#ifdef BLACK_HOLES
	  if(P[i].Mass == 0)
	    flag = 1;
#endif

	  if(P[i].Mass > 0)
	    SNEgy = SphP[i].EgyRes / P[i].Mass;
	  else
	    SNEgy = 0;

#ifdef WINDS

#ifdef LT_TRACK_WINDS
	  trackwinds_flag = '\0';
#endif

	  if(SphP[i].DelayTime > 0)
	    {
	      flag = 1;		/* only normal cooling for particles in the wind */
	      SphP[i].DelayTime -= dtime;
	      if((SphP[i].DelayTime < 0) ||
		 ((SphP[i].DelayTime > 0) &&
		  (SphP[i].d.Density * a3inv <
		   0.5 * All.WindFreeTravelDensFac * SFs[mySFi].PhysDensThresh[0])))
		{
		  SphP[i].DelayTime = 0;
#ifdef LT_TRACK_WINDS
		  trackwinds_flag = 'R';
#endif
		}

#ifdef LT_DECOUPLE_POSTWINDS_FROM_SF
	      if(SphP[i].DelayTime == 0)
		{
#ifndef LT_WIND_VELOCITY
		  v = sqrt(2 * SFs[mySFi].WindEnergyFraction / SFs[mySFi].WindEfficiency *
			   (SFs[mySFi].totFactorSN / (1 - SFs[mySFi].totFactorSN) * myEgySpecSN + SNEgy));
#else
		  v = LT_WIND_VELOCITY;
#endif
		  SphP[i].DelayTimeSF = All.WindFreeTravelDensFac / v / 2;
		}
#endif
	    }

#ifdef UM_CONTINUE
	  if(SphP[i].DelayTime > 0)
	    continue;		/* neglect cooling in winds */
#endif

#ifdef LT_DECOUPLE_POSTWINDS_FROM_SF
	  if(SphP[i].DelayTimeSF > 0)
	    {
	      flag = 1;
	      SphP[i].DelayTimeSF -= dtime;
	      if((SphP[i].DelayTimeSF < 0) ||
		 ((SphP[i].DelayTimeSF > 0) &&
		  (SphP[i].d.Density * a3inv <
		   0.5 * All.WindFreeTravelDensFac * SFs[mySFi].PhysDensThresh[0])))
		{
		  SphP[i].DelayTimeSF = 0;
#ifdef LT_TRACK_WINDS
		  trackwinds_flag = 'S';
#endif
		}
	    }
#endif

#ifdef LT_TRACK_WINDS
	  if(trackwinds_flag != '\0')
	    {
	      ne_guess = 1.0;
	      ignore_failure_in_convert_u = 1;
#if !defined(LT_METAL_COOLING_WAL)
	      tw_temp = convert_u_to_tempSTD(uold * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
					     SphP[i].d.Density * All.UnitDensity_in_cgs * All.HubbleParam *
					     All.HubbleParam * a3inv, &ne_guess);
#else
	      tw_temp = convert_u_to_tempMET(uold * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
					     Redshift, DZ, &Metallicities[0]);
#endif

	      fprintf(FdTrackW,
		      "[%c] %8.6e id: %10u %2d pos: %8.6e %8.6e %8.6e rho: %8.6e %8.6e Z: %8.6e %8.6e"
		      " h: %8.6e %8.6e %8.6e\n", trackwinds_flag, All.Time, (P[i].ID << All.StarBits),
		      (P[i].ID >> (sizeof(MyIDType) * 8 - All.StarBits)), P[i].Pos[0], P[i].Pos[1],
		      P[i].Pos[2],
		      SphP[i].d.Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * a3inv,
		      tw_temp, get_metallicity(i, -1), x, PPP[i].Hsml, SphP[i].AvgHsml, P[i].p.Potential);
	      ignore_failure_in_convert_u = 0;
	    }
#endif

#endif


#ifdef QUICK_LYALPHA
	  temp = u_to_temp_fac * SphP[i].Entropy /
	    GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

	  if(SphP[i].d.Density > All.OverDensThresh && temp < 1.0e5)
	    flag = 0;
	  else
	    flag = 1;
#endif

/* 	  if(All.NumCurrentTiStep> 0) */
/* 	    printf("\t cooling: flag is %d : %d ** %d\n", flag, ThisTask, All.NumCurrentTiStep); fflush(stdout); */
#if !defined(NOISMPRESSURE) && !defined(QUICK_LYALPHA)
	  if(flag == 1)		/* normal implicit isochoric cooling */
#endif
	    {
	      SphP[i].Sfr = 0;

	      ne = SphP[i].elec;

#ifdef RT_COOLING_PHOTOHEATING
	      unew = radtransfer_cooling_photoheating(i, dtime);	/* !!! here unew = dtemp/molecular_weight actually !!! */

	      if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
		{
		  /* note: the adiabatic rate has been already added in ! */

		  if(dt > 0)
		    {
		      SphP[i].e.DtEntropy = unew * (BOLTZMANN / All.UnitEnergy_in_cgs) /
			pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / (PROTONMASS / All.UnitMass_in_g) / dt;

		      if(SphP[i].e.DtEntropy < -0.5 * SphP[i].Entropy / dt)
			SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt;
		    }
		}
#else

#ifdef LT_METAL_COOLING_WAL
//            WalCool_set_PID(P[i].ID);
#ifndef UM_CHEMISTRY
	      unew =
		DoCoolingMET(uold, SphP[i].d.Density * a3inv, &Metallicities[0], Redshift, DZ, dtime,
			     &Temperature);
#else
	      unew =
		Um_DoCooling(uold, SphP[i].d.Density * a3inv, dtime, &ne, &Metallicities[0], Redshift, DZ, i,
			     flag);
#endif
#else
#ifdef LT_METAL_COOLING
#ifdef UM_CHEMISTRY
	      unew = Um_DoCooling(uold, SphP[i].d.Density * a3inv, dtime, &ne, Zcool, i, flag);
#else
	      unew = DoCoolingZ(uold, SphP[i].d.Density * a3inv, dtime, &ne, Zcool, &Temperature);
#endif
#else
	      unew = DoCoolingSTD(uold, SphP[i].d.Density * a3inv, dtime, &ne);
#endif /* closes LT_METAL_COOLING */
#endif /* closes LT_METAL_COOLING_WAL */
#endif /* closes RT_COOLING_PHOTOHEATING */

	      unew += SNEgy;

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
	      if(SphP[i].i.Injected_BH_Energy)
		{
		  if(P[i].Mass == 0)
		    SphP[i].i.Injected_BH_Energy = 0;
		  else
		    unew += SphP[i].i.Injected_BH_Energy / P[i].Mass;
#ifdef LT_MAX_TEMP_FEEDBACK
		  temp = u_to_temp_fac * unew;
		  if(temp > LT_MAX_TEMP_FEEDBACK)
		    unew = LT_MAX_TEMP_FEEDBACK / u_to_temp_fac;
#endif

#ifdef FLTROUNDOFFREDUCTION
		  SphP[i].i.dInjected_BH_Energy = 0;
#else
		  SphP[i].i.Injected_BH_Energy = 0;
#endif
		}
#endif


	      SphP[i].elec = ne;
	      SphP[i].Entropy = unew / fac_entr_to_u;

	      SphP[i].EntropyPred = SphP[i].Entropy;
#if GADGET_HYDRO == HYDRO_MFM
	      SphP[i].InternalEnergy = unew;
	      SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
	      mfm_synchronize_entropy();
#endif

	      SphP[i].Pressure = get_pressure(i);
	    }

	  if(flag == 0)		/* active star formation */
	    {
#if !defined(QUICK_LYALPHA)
	      tsfr = sqrt(myPhysDensThresh / (SphP[i].d.Density * a3inv)) * myMaxSfrTimescale;
	      factorEVP = pow(SphP[i].d.Density * a3inv / myPhysDensThresh, -0.8) * myFactorEVP;

	      egyhot = myEgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = SphP[i].elec;

#ifndef LT_METAL_COOLING_WAL
#ifndef LT_METAL_COOLING
#ifndef RT_COOLING_PHOTOHEATING
	      tcool = GetCoolingTimeSTD(egyhot, SphP[i].d.Density * a3inv, &ne);
#else
	      tcool = rt_GetCoolingTime(i, egyhot, SphP[i].d.Density * a3inv);
#endif
#else
#ifndef UM_CHEMISTRY
	      /*Z = get_metallicity(i); already done */
	      if((tcool = GetCoolingTimeZ(egyhot, SphP[i].d.Density * a3inv, &ne, Zcool, &Temperature)) == 0)
		tcool = 1e13;
#else
	      if((tcool = Um_GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne, Zcool, i)) == 0)
		tcool = 1e13;
#endif
#endif
#else
	      tcool =
		GetCoolingTimeMET(egyhot, SphP[i].d.Density * a3inv, Redshift, DZ, &Metallicities[0],
				  &Temperature);
#endif

	      SphP[i].elec = ne;

	      y = tsfr / tcool * egyhot / (myFactorSN * myEgySpecSN - (1 - myFactorSN) * All.EgySpecCold);

	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      cloudmass = x * P[i].Mass;

	      if(tsfr < dtime)
		tsfr = dtime;

#ifndef LT_SMOOTH_XCLD
	      sm = dtime / tsfr * cloudmass;
#else
	      if(x > 0)
		sm = dtime / tsfr * SphP[i].XCLDsmooth * P[i].Mass;
	      else
		sm = 0;
#endif

	      p = sm / P[i].Mass;

	      sum_sm[mySFi] += P[i].Mass * (1 - exp(-p));

	      SphP[i].Sfr = (1 - myFactorSN) * cloudmass / tsfr *
		(All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

	      TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
	      sfrrates[mySFi] += SphP[i].Sfr;

	      if(SphP[i].Sfr > 0)
		{
		  int byrho_idx, byZ_idx;
		  byrho_idx = (int) (log10(SphP[i].d.Density * a3inv / rhomean) - SFR_by_DENSITY_LOGMIN);
		  if(byrho_idx < 0)
		    byrho_idx = 0;
		  if(byrho_idx > NBINS_SFR_by_DENSITY)
		    byrho_idx = NBINS_SFR_by_DENSITY;
		  sfrrates_bydensity[byrho_idx] += SphP[i].Sfr;

		  byZ_idx = (int) (Zcool - (NO_METAL));
		  if(byZ_idx < 0)
		    byZ_idx = 0;
		  if(byZ_idx > NBINS_SFR_by_Z)
		    byZ_idx = NBINS_SFR_by_Z;
		  sfrrates_byZ[byZ_idx] += SphP[i].Sfr;
		}

	      exp_minus_p = exp(-p);
#ifndef LT_STOCHASTIC_IRA
	      mstar = sm * exp_minus_p;
#else
	      mstar = P[i].Mass * (1 - exp_minus_p);
#endif


	      if(dt > 0)
		{
		  if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
		    {
		      trelax = tsfr * (1 - x) / x / (myFactorSN * (1 + factorEVP));
		      egycurrent =
			SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;

#ifndef LT_SNegy_IN_HOTPHASE
		      egycurrent += SNEgy;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
		      if(SphP[i].i.Injected_BH_Energy > 0)
			{
			  egycurrent += SphP[i].i.Injected_BH_Energy / P[i].Mass;
#ifdef LT_MAX_TEMP_FEEDBACK
			  temp = u_to_temp_fac * egycurrent;
			  if(temp > LT_MAX_TEMP_FEEDBACK)
			    egycurrent = LT_MAX_TEMP_FEEDBACK / u_to_temp_fac;
#endif

			  if(egycurrent > egyeff)
			    {

#if !defined(LT_METAL_COOLING) && !defined(LT_METAL_COOLING_WAL) && !defined(RT_COOLING_PHOTOHEATING)	/* normal cooling */
			      if((tcool = GetCoolingTimeSTD(egyhot, SphP[i].d.Density * a3inv, &ne)) == 0)
				tcool = 1e13;
#else

#if defined(LT_METAL_COOLING)	/* metal cooling S & D */
			      ne = 1.0;
#ifndef UM_CHEMISTRY		/* low-T metal & molecular cooling */
			      if((tcool =
				  GetCoolingTimeZ(egycurrent, SphP[i].d.Density * a3inv, &ne, Zcool,
						  &Temperature)) == 0)
				tcool = 1e13;
#else
			      ne = SphP[i].elec;
			      if((tcool =
				  Um_GetCoolingTime(egycurrent, SphP[i].d.Density * a3inv, &ne, Zcool,
						    i)) == 0)
				tcool = 1e13;
#endif
#endif

#if defined(LT_METAL_COOLING_WAL)	/* WAL metal cooling */
			      tcool =
				GetCoolingTimeMET(egycurrent, SphP[i].d.Density * a3inv, Redshift, DZ,
						  &Metallicities[0], &Temperature);
#endif

#if defined(RT_COOLING_PHOTOHEATING)	/* RT cooling */
			      tcool = rt_GetCoolingTime(i, egycurrent, SphP[i].d.Density * a3inv);
#endif

#endif
			      if(tcool < trelax && tcool > 0)
				trelax = tcool;
			    }

			  SphP[i].i.Injected_BH_Energy = 0;
			}
#endif



#if !defined(NOISMPRESSURE)
		      SphP[i].Entropy =
			(egyeff +
			 (egycurrent -
			  egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 /
			pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

		      SphP[i].e.DtEntropy = 0;
#if GADGET_HYDRO == HYDRO_MFM
		      SphP[i].InternalEnergy = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));
#endif
#endif

		      SphP[i].EntropyPred = SphP[i].Entropy;
#if GADGET_HYDRO == HYDRO_MFM
		      SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
		      mfm_synchronize_entropy();
#endif
		      SphP[i].Pressure = get_pressure(i);

#ifdef LT_EXTEGY_INFO
		      if(EEInfo_grain > SEvInfo_GRAIN)
			{
			  tstart = second();
			  if(SphP[i].EgyRes > 0)
			    {
			      StatSum[REScount]++;

			      r = SNEgy / SFs[mySFi].TOT_erg_per_g;
			      if(r < StatMin[OEgy] && r > 0)
				StatMin[OEgy] = r;
			      if(r > StatMax[OEgy])
				StatMax[OEgy] = r;
			      StatSum[OEgy] += r;

			      r *= P[i].Mass / mstar;
			      if(r < StatMin[OEgy_I] && r > 0)
				StatMin[OEgy_I] = r;
			      if(r > StatMax[OEgy_I])
				StatMax[OEgy_I] = r;
			      StatSum[OEgy_I] += r;

			    }

			  tend = second();
			  CPU_ee_info = timediff(tstart, tend);
			}
#endif
		    }
		}

	      if(NonProcOn_II[Yset])
		NonMetalMass = P[i].Mass - get_metalmass(SphP[i].Metals);

#if defined(LT_LOCAL_IRA)
	      if(UseSnII)
		{
		  tot_mass = 0;
#if defined (LT_ZAGE) || defined(LT_POPIII)
		  metal_add_mass = 0;
#endif
		  if(SFs[mySFi].nonZeroIRA)
		    {

		      for(j = 0; j < LT_NMetP; j++)
			{
			  /* add the newly formed elements */

			  mz = mstar * SnII_ShortLiv_Yields[Yset][j][YZbin];

			  /* subtract the fraction locked in stars at birth
			   *   m_star * (All.MassFrac_inIRA - myFactorSN) * SphP[i].Metals[j] / P[i].Mass
			   */
			  /*
			     mz -= (1 - exp_minus_p) * SphP[i].Metals[j] * (All.MassFrac_inIRA - myFactorSN);
			   */

			  if(NonProcOn_II[Yset])
			    mz += sm * exp_minus_p * SphP[i].Metals[j] / NonMetalMass *
			      (SnII_ShortLiv_Yields[Yset][Hyd][YZbin] +
			       SnII_ShortLiv_Yields[Yset][Hel][YZbin]);

			  /*
			   *           MassFrac_inIRA = mass fraction of stars in IRA Sn
			   *           FactorSN = restored fraction by IRA Sn
			   *           MassFra_inIRA - FactorSN = locked fraction
			   */

			  tot_mass += mz;
#if defined (LT_ZAGE) || defined(LT_POPIII)
			  if(j != Hel)
			    metal_add_mass += mz;
#endif

#ifdef LT_ZAGE_LLV
			  if(j == Iron)
			    llvmetal_add_mass = mz;
#endif

			  if(SphP[i].Metals[j] + mz < 0)
			    {
			      printf(" \n\n !!! GULP !!! \n %i %i %i :: %g %g %g %g\n", ThisTask, P[i].ID, j,
				     SphP[i].Metals[j], mz, sm, p);
			      endrun(101);
			    }
			  else
			    {
#ifndef LT_TRACK_CONTRIBUTES
			      SphP[i].Metals[j] += mz;
#else
			      contrib_metals[j] = mz;
#endif
#ifdef LT_SEv_INFO
			      mZ[j] += mz;
#endif
			    }

			}

		      if(tot_mass > P[i].Mass)
			{
			  printf(" \n\n !!! GULP 2 !!! \n %i %i %i :: %g %g %g %g\n", ThisTask, P[i].ID, j,
				 tot_mass, P[i].Mass, sm, p);
			  endrun(102);
			}
		      tot_mass += mstar * SnII_ShortLiv_Yields[Yset][LT_NMetP][YZbin];

#ifdef LT_EJECTA_IN_HOTPHASE
		      SphP[i].XColdCloud = (P[i].Mass * SphP[i].XColdCloud - mstar) / P[i].Mass;
		      x = tot_mass / P[i].Mass;
		      egyeff = (egyhot * (1 - SphP[i].XColdCloud + x) + All.EgySpecCold * SphP[i].XColdCloud);
#if !defined(NOISMPRESSURE)
		      SphP[i].Entropy =
			(egyeff +
			 (egycurrent -
			  egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 /
			pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
#if GADGET_HYDRO == HYDRO_MFM
		      SphP[i].InternalEnergy = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));
#endif
#endif
#endif

#if defined(LT_TRACK_CONTRIBUTES)
		      NULL_CONTRIB(&contrib);
		      NULL_EXPCONTRIB(&IIcontrib);
		      NULL_EXPCONTRIB(&Iacontrib);
		      NULL_EXPCONTRIB(&AGBcontrib);

		      for(j = 0; j < LT_NMetP; j++)
			{
			  IIcontrib[j] = 1;
			  Iacontrib[j] = AGBcontrib[j] = 0;
			}

		      pack_contrib(&contrib, mySFi, IIcontrib, Iacontrib, AGBcontrib);

		      update_contrib(&SphP[i].contrib, &SphP[i].Metals[0], &contrib, &contrib_metals[0]);
		      for(j = 0; j < LT_NMetP; j++)
			SphP[i].Metals[j] += contrib_metals[j];
#endif /*  <<<----   actually add metals   */

		      double zage_term;

#if defined(LT_ZAGE)
		      zage_term = (cosmic_time - All.Time_Age) * metal_add_mass / P[i].Mass;
#ifdef LT_LOGZAGE
		      zage_term = log10(zage_term);
#endif
		      SphP[i].ZAge += zage_term;
		      SphP[i].ZAgeW += metal_add_mass / P[i].Mass;
#endif
#if defined(LT_ZAGE_LLV)
		      zage_term = (cosmic_time - All.Time_Age) * llvmetal_add_mass / P[i].Mass;
#ifdef LT_LOGZAGE
		      zage_term = log10(zage_term);
#endif
		      SphP[i].ZAge_llv += zage_term;
		      SphP[i].ZAgeW_llv += llvmetal_add_mass / P[i].Mass;
		      /*            SphP[i].ZAge += (cosmic_time - All.Time_Age) * get_metallicity(j, -1); */
		      /*            SphP[i].ZAgeW += get_metallicity(j, -1); */
#endif
		    }
		}
#else /* closes ifdef(LT_LOCAL_IRA) */
	      if(SFs[mySFi].nonZeroIRA)
		SphP[i].mstar = mstar;
#endif

#endif /* belongs to ifndef(QUICK_LYALPHA) */

	      /* the upper bits of the gas particle ID store how man stars this gas
	         particle gas already generated                                     */
	      number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - All.StarBits));

	      mass_of_star = P[i].Mass / (All.Generations - number_of_stars_generated) * GENERATIONS_BUNCH;


#ifndef QUICK_LYALPHA
	      prob = P[i].Mass / mass_of_star * (1 - exp(-p));
#else
	      prob = 2.0;	/* this will always cause a star creation event */
#endif /* ends to QUICK_LYALPHA */

	      if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
		{
		  NextChemTime = get_NextChemTime(0, mySFi, 0x0);	/* !! note: this is look-back time */

		  if((number_of_stars_generated + GENERATIONS_BUNCH) >= All.Generations)
		    {
		      /* here we turn the gas particle itself into a star */
		      Stars_converted++;
		      stars_converted++;

		      sum_mass_stars[mySFi] += P[i].Mass;

		      P[i].Type = 4;
		      TimeBinCountSph[P[i].TimeBin]--;
		      TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;

		      if(PPP[i].Hsml > All.MaxChemSpreadL)
			PPP[i].Hsml = All.MaxChemSpreadL;

#ifdef LT_SEvDbg
		      if(ThisTask == 0 && FirstID == 0)
			if(get_random_number(P[NumPart + stars_spawned].ID + 3) > 0.5)
			  FirstID = P[i].ID;
#endif
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

		      MetP[N_stars].iMass = P[i].Mass;
#ifdef STELLARAGE
		      MetP[N_stars].StellarAge = All.Time;
#endif
		      for(j = 0; j < LT_NMetP; j++)
			MetP[N_stars].Metals[j] = SphP[i].Metals[j];

#ifdef GL_CR_DUST
		      for(j = 0; j < LT_NMetP; j++)
			MetP[N_stars].Metals[j] += SphP[i].DustL[j] + SphP[i].DustS[j];
#endif

#ifdef LT_SMOOTHZ_IN_IMF_SWITCH
		      MetP[N_stars].Zsmooth = SphP[i].Zsmooth;
#endif

#ifdef LT_TRACK_CONTRIBUTES
		      MetP[N_stars].contrib = SphP[i].contrib;
#endif

#ifdef LT_ZAGE
#ifndef LT_LOGZAGE
		      MetP[N_stars].ZAge = SphP[i].ZAge / SphP[i].ZAgeW;
#else
		      MetP[N_stars].ZAge = exp(10, SphP[i].ZAge) / SphP[i].ZAgeW;
#endif
#endif

#ifdef LT_ZAGE_LLV
#ifndef LT_LOGZAGE
		      MetP[N_stars].ZAge_llv = SphP[i].ZAge_llv / SphP[i].ZAgeW_llv;
#else
		      MetP[N_stars].ZAge_llv = exp(10, SphP[i].ZAge_llv) / SphP[i].ZAgeW_llv;
#endif
#endif

		      MetP[N_stars].LastChemTime = SNtimesteps[mySFi][0][1];

		      index_of_star = i;

		      P[i].Mass *= (1 - SFs[mySFi].metFactorSN);

		      N_stars++;
		    }
		  else
		    {
		      /* here we spawn a new star particle */

		      if(NumPart + stars_spawned >= All.MaxPart)
			{
			  printf
			    ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			     ThisTask, NumPart, stars_spawned, All.MaxPart);
			  fflush(stdout);
			  endrun(8888);
			}

		      if(N_stars + 1 >= All.MaxPartMet)
			{
			  printf
			    ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPartMet=%d)\n",
			     ThisTask, NumPart, stars_spawned, All.MaxPartMet);
			  fflush(stdout);
			  endrun(8889);
			}

		      P[NumPart + stars_spawned] = P[i];
		      P[NumPart + stars_spawned].Type = 4;

		      if(P[NumPart + stars_spawned].Hsml > All.MaxChemSpreadL)
			P[NumPart + stars_spawned].Hsml = All.MaxChemSpreadL;

		      P[i].ID += (GENERATIONS_BUNCH << (sizeof(MyIDType) * 8 - All.StarBits));

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
		      MetP[N_stars].iMass = P[NumPart + stars_spawned].Mass;

		      factor = P[NumPart + stars_spawned].Mass / P[i].Mass;

		      P[i].Mass -= P[NumPart + stars_spawned].Mass;
		      sum_mass_stars[mySFi] += P[NumPart + stars_spawned].Mass;
		      P[NumPart + stars_spawned].Mass *= (1 - SFs[mySFi].metFactorSN);

#ifdef LT_SEvDbg
		      if(ThisTask == 0 && FirstID == 0)
			if(get_random_number(P[NumPart + stars_spawned].ID + 3) > 0.5)
			  FirstID = P[NumPart + stars_spawned].ID;
#endif
		      P[NumPart + stars_spawned].pt.MetID = N_stars;
		      MetP[N_stars].PID = NumPart + stars_spawned;
#ifdef LT_SMOOTHZ_IN_IMF_SWITCH
		      MetP[N_stars].Zsmooth = SphP[i].Zsmooth;
#endif
#ifdef LT_ZAGE
		      MetP[N_stars].ZAge = SphP[i].ZAge / SphP[i].ZAgeW;
#endif
#ifdef LT_ZAGE_LLV
		      MetP[N_stars].ZAge_llv = SphP[i].ZAge_llv / SphP[i].ZAgeW_llv;
#endif
		      MetP[N_stars].LastChemTime = SNtimesteps[mySFi][0][1];

		      index_of_star = NumPart + stars_spawned;

		      for(j = 0; j < LT_NMetP; j++)
			{
#ifdef LT_LOCAL_IRA
			  SphP[i].Metals[j] = DMAX(SphP[i].Metals[j] - mstar * SnII_ShortLiv_Yields[Yset][j][YZbin], 0);	/* to avoid round-off errors */
#endif
			  MetP[N_stars].Metals[j] = SphP[i].Metals[j] * factor;
			  SphP[i].Metals[j] *= (1 - factor);
#ifdef LT_LOCAL_IRA
			  SphP[i].Metals[j] += mstar * SnII_ShortLiv_Yields[Yset][j][YZbin];
#endif
			}
#ifdef STELLARAGE
		      MetP[N_stars].StellarAge = All.Time;
#endif

#ifdef LT_TRACK_CONTRIBUTES
		      MetP[N_stars].contrib = SphP[i].contrib;
#endif
		      N_stars++;

		      force_add_star_to_tree(i, NumPart + stars_spawned);

		      stars_spawned++;
		    }

		  bin = get_chemstep_bin(All.Time, All.Time_Age - NextChemTime, &chem_step, index_of_star);
		  /*                 if(TimeBinActive[bin] == 0) */
		  /*                   { */
		  /*                     while(TimeBinActive[bin] == 0 && bin > 0) */
		  /*                       bin--; */
		  /*                     chem_step = bin ? (1 << bin) : 0; */
		  /*                   }                 */
		  MetP[P[index_of_star].pt.MetID].ChemTimeBin = bin;

		  TimeBinCountStars[bin]++;
		}

#ifdef WINDS
	      /* Here comes the wind model */

	      if(P[i].Type == 0 &&	/* to protect using a particle that has been turned into a star */
		 SFs[mySFi].WindEfficiency > 0)
		{
		  p = SFs[mySFi].WindEfficiency * sm / P[i].Mass;

		  prob = 1 - exp(-p);

		  if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
		    {


#ifndef LT_WIND_VELOCITY
		      v =
			sqrt(2 * SFs[mySFi].WindEnergyFraction / SFs[mySFi].WindEfficiency *
			     (SFs[mySFi].totFactorSN / (1 - SFs[mySFi].totFactorSN) * myEgySpecSN + SNEgy));

#ifdef LT_SEv_INFO

		      if(v < windv_min)
			windv_min = v;
		      if(v > windv_max)
			windv_max = v;
		      windv += v;
		      windn++;

#endif
#else
		      v = LT_WIND_VELOCITY;
#ifdef LT_SEv_INFO
		      SFs[mySFi].WindEnergyFraction = v * v * SFs[mySFi].WindEfficiency /
			(2 * (SFs[mySFi].totFactorSN / (1 - SFs[mySFi].totFactorSN) * myEgySpecSN + SNEgy));

		      if(SFs[mySFi].WindEnergyFraction < windv_min)
			windv_min = SFs[mySFi].WindEnergyFraction;
		      if(SFs[mySFi].WindEnergyFraction > windv_max)
			windv_max = SFs[mySFi].WindEnergyFraction;
		      windv += SFs[mySFi].WindEnergyFraction;
		      windn++;
#endif
#endif

#ifdef LT_TRACK_WINDS
		      ne_guess = 1.0;
#ifndef LT_METAL_COOLING_WAL
		      tw_temp = convert_u_to_tempSTD(DMAX(All.MinEgySpec,
							  (SphP[i].Entropy +
							   SphP[i].e.DtEntropy * dt) / GAMMA_MINUS1 *
							  pow(SphP[i].d.Density * a3inv,
							      GAMMA_MINUS1)) * All.UnitPressure_in_cgs /
						     All.UnitDensity_in_cgs,
						     SphP[i].d.Density * All.UnitDensity_in_cgs *
						     All.HubbleParam * All.HubbleParam * a3inv, &ne_guess);
#else
		      tw_temp = convert_u_to_tempMET(DMAX(All.MinEgySpec,
							  (SphP[i].Entropy +
							   SphP[i].e.DtEntropy * dt) / GAMMA_MINUS1 *
							  pow(SphP[i].d.Density * a3inv,
							      GAMMA_MINUS1)) * All.UnitPressure_in_cgs /
						     All.UnitDensity_in_cgs, Redshift, DZ, &Metallicities[0]);
#endif

		      fprintf(FdTrackW,
			      "[G] %8.6e id: %10u %2d pos: %8.6e %8.6e %8.6e rho: %8.6e %8.6e Z: %8.6e %8.6e"
			      " h: %8.6e %8.6e %8.6e\n", All.Time, (P[i].ID << All.StarBits),
			      (P[i].ID >> (32 - All.StarBits)), P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
			      SphP[i].d.Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam *
			      a3inv, tw_temp, get_metallicity(i, -1), x, PPP[i].Hsml, SphP[i].AvgHsml,
			      P[i].p.Potential);
#endif

#ifdef ISOTROPICWINDS
		      theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
		      phi = 2 * M_PI * get_random_number(P[i].ID + 4);

		      dir[0] = sin(theta) * cos(phi);
		      dir[1] = sin(theta) * sin(phi);
		      dir[2] = cos(theta);
#else
		      dir[0] = P[i].g.GravAccel[1] * P[i].Vel[2] - P[i].g.GravAccel[2] * P[i].Vel[1];
		      dir[1] = P[i].g.GravAccel[2] * P[i].Vel[0] - P[i].g.GravAccel[0] * P[i].Vel[2];
		      dir[2] = P[i].g.GravAccel[0] * P[i].Vel[1] - P[i].g.GravAccel[1] * P[i].Vel[0];
#endif

		      for(j = 0, norm = 0; j < 3; j++)
			norm += dir[j] * dir[j];

		      norm = sqrt(norm);
		      if(get_random_number(P[i].ID + 5) < 0.5)
			norm = -norm;

		      if(norm != 0)
			{
			  for(j = 0; j < 3; j++)
			    dir[j] /= norm;

			  for(j = 0; j < 3; j++)
			    {
			      P[i].Vel[j] += v * ascale * dir[j];
			      SphP[i].VelPred[j] += v * ascale * dir[j];
			    }

#ifndef LT_HYDROWINDS
#ifdef UM_WIND_DELAYTIME
			  SphP[i].DelayTime = DMAX(2. * PPP[i].Hsml, All.WindFreeTravelLength) / v;
#else
			  SphP[i].DelayTime = All.WindFreeTravelMaxTimeFactor / hubble_function(All.Time);
#endif
#endif
#ifdef LT_HOT_WINDS
			  u_hotwinds = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.16e7;
			  u_hotwinds *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
			  u_hotwinds /= 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));

			  SphP[i].Entropy =
			    GAMMA_MINUS1 * u_hotwinds / pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);
#if GADGET_HYDRO == HYDRO_MFM
			  SphP[i].InternalEnergy = u_hotwinds;	// consider the modified density here ???
#endif
#endif
			}
		    }
		}
#endif
	    }

	  SphP[i].EgyRes = 0;

	  SphP[i].XColdCloud = (float) x;
	  SphP[i].Temperature = (float) Temperature;

	}

    }				/* end of main loop over active particles */

#ifdef LT_SEvDbg
  if(checkFirstID == 0)
    {
      MPI_Bcast(&FirstID, 1, MPI_UNSIGNED, 0, MYMPI_COMM_WORLD);
      if(FirstID > 0)
	checkFirstID = 1;
    }
#endif

#ifdef LT_TRACK_WINDS
  fflush(FdTrackW);
#endif


  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("SFR: spawned %d stars, converted %d gas particles into stars\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}


      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;
      All.TotN_stars += tot_spawned + tot_converted;

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */

      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

#ifdef _OPENMP
  if(stars_spawned > 0)
    build_active_particle_list();
#endif

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(sfrrates, totsfrrates, SFs_dim, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(sfrrates_bydensity, totsfrrates_bydensity, NBINS_SFR_by_DENSITY, MPI_DOUBLE, MPI_SUM,
		MYMPI_COMM_WORLD);
  MPI_Allreduce(sfrrates_byZ, totsfrrates_byZ, NBINS_SFR_by_DENSITY, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);

  for(i = 0; i < SFs_dim; i++)
    {
      Sum_sm += sum_sm[i];
      Sum_mass_stars += sum_mass_stars[i];
    }
#ifdef LT_SEv_INFO
  MPI_Allreduce(&Sum_sm, &Total_sm, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
#else
  MPI_Reduce(&Sum_sm, &Total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
#endif
  MPI_Reduce(&Sum_mass_stars, &Total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&sum_sm[0], &total_sm[0], SFs_dim, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars[0], &total_sum_mass_stars[0], SFs_dim, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	{
	  rate = Total_sm / (All.TimeStep / time_hubble_a);

	  /* convert to solar masses per yr */

	  rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

	  fprintf(FdSfr, "%10.8e %10.8e %10.8e %10.8e %10.8e ", All.Time, Total_sm, totsfrrate,
		  rate_in_msunperyear, Total_sum_mass_stars);
	  for(i = 0; i < SFs_dim; i++)
	    fprintf(FdSfr, "%10.8e %10.8e ", totsfrrates[i], total_sm[i] / (All.TimeStep / time_hubble_a) *
		    (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR));

	  for(i = 0; i < SFs_dim; i++)
	    fprintf(FdSfr, "%10.8e ", total_sum_mass_stars[i]);

	  fprintf(FdSfrDens, "%10.8e ", All.Time);
	  for(i = 0; i < NBINS_SFR_by_DENSITY; i++)
	    fprintf(FdSfrDens, "%10.8e ", totsfrrates_bydensity[i]);

	  fprintf(FdSfrZ, "%10.8e ", All.Time);
	  for(i = 0; i < NBINS_SFR_by_Z; i++)
	    fprintf(FdSfrZ, "%10.8e ", totsfrrates_byZ[i]);
	}
      else
	{
	  fprintf(FdSfr, "%10.8e %10.8e 0 0 %10.8e ", All.Time, Total_sm, Total_sum_mass_stars);
	  for(i = 0; i < SFs_dim; i++)
	    fprintf(FdSfr, "0 ");
	  for(i = 0; i < SFs_dim; i++)
	    fprintf(FdSfr, "%10.8e ", total_sum_mass_stars[i]);

	  fprintf(FdSfrDens, "%10.8e ", All.Time);
	  for(i = 0; i < NBINS_SFR_by_DENSITY; i++)
	    fprintf(FdSfrDens, "0 ");
	  for(i = 0; i < NBINS_SFR_by_Z; i++)
	    fprintf(FdSfrZ, "0 ");
	}
      fprintf(FdSfr, "\n");
      fprintf(FdSfrDens, "\n");
      fprintf(FdSfrZ, "\n");

      fflush(FdSfr);
    }

#ifdef LT_EXTEGY_INFO

  /* collect some infos about external energy in the reservoir SphP[].EgyRes.
     could be used for any egy source (AGN etc).
   */
  if(EEInfo_grain > SEvInfo_GRAIN)
    {
      tstart = second();

      MPI_Reduce(&StatSum[OEgy], &tot_StatSum[OEgy], 3, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
      MPI_Reduce(&StatMin[OEgy], &tot_StatMin[OEgy], 3, MPI_DOUBLE, MPI_MIN, 0, MYMPI_COMM_WORLD);
      MPI_Reduce(&StatMax[OEgy], &tot_StatMax[OEgy], 3, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);

      if(ThisTask == 0)
	{
	  fprintf(FdExtEgy, "%8.6lg %8.6lg %9d ", All.Time, tot_StatSum[REScount], N_stars);

	  if(tot_StatSum[REScount])
	    {
	      if(tot_StatSum[REScount] > 2)
		{
		  tot_StatSum[OEgy] = (tot_StatSum[OEgy] -=
				       tot_StatMin[OEgy] + tot_StatMax[OEgy]) / (tot_StatSum[REScount] - 2);
		  tot_StatSum[OEgy_I] = (tot_StatSum[OEgy_I] -=
					 tot_StatMin[OEgy_I] + tot_StatMax[OEgy]) / (tot_StatSum[REScount] -
										     2);
		}
	      else
		{
		  tot_StatSum[OEgy] /= tot_StatSum[REScount];
		  tot_StatSum[OEgy_I] /= tot_StatSum[REScount];
		}

	      fprintf(FdExtEgy, "%8.6lg %8.6lg %8.6lg " "%8.6lg %8.6lg %8.6lg ",
		      /* (3-5) {min,max,mean} ratio between egy of reservoir and egy due to IRA SnII */
		      tot_StatMin[OEgy], tot_StatMax[OEgy], tot_StatSum[OEgy],
		      tot_StatMin[OEgy_I], tot_StatMax[OEgy_I], tot_StatSum[OEgy_I]);

	    }

	  fprintf(FdExtEgy, "\n");
	  fflush(FdExtEgy);
	}
      tend = second();
      CPU_ee_info += timediff(tstart, tend);

      MPI_Reduce(&CPU_ee_info, &sumCPU_ee_info, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
      if(ThisTask == 0)
	All.CPU_EE_info += sumCPU_ee_info / NTask;

      EEInfo_grain = 0;
    }
#endif /* closes LT_EXTEGY_INFO */

#ifdef LT_SEv_INFO

  /* collect some informations */

  tstart = second();
#ifdef LT_LOCAL_IRA
  if(UseSnII)
    {
      /* write metals produced by SnII in IRA */

      if(Total_sm > 0)
	{
	  MPI_Reduce(&mZ[0], &tot_mZ[0], LT_NMet, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
	  if(ThisTask == 0)
	    {
	      for(j = 0; j < SFs_dim; j++)
		{
		  fprintf(FdSn, "IRA %1d %.8lg 0 0 %g", j, All.Time,
			  total_sm[j] * All.UnitMass_in_g / SOLAR_MASS * 1);
		  for(i = 0; i < LT_NMet; i++)
		    fprintf(FdSn, "%.8lg ", tot_mZ[i] * All.UnitMass_in_g / SOLAR_MASS);
		  fprintf(FdSn, "\n");
		}
	      fflush(FdSn);
	    }
	}

    }
#endif

#ifdef WINDS
  /* write the minimum, maximum and mean wind velocity (to check the external egy contrib) */
  MPI_Allreduce(&windn, &tot_windn, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  if(tot_windn > 0)
    {
      if(windv_min == 1e6)
	windv_min = 0;
      MPI_Reduce(&windv, &tot_windv, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
      MPI_Reduce(&windv_min, &tot_windv_min, 1, MPI_DOUBLE, MPI_MIN, 0, MYMPI_COMM_WORLD);
      MPI_Reduce(&windv_max, &tot_windv_max, 1, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);
      if(ThisTask == 0)
	{
	  if(tot_windn > 1)
	    tot_windv = (tot_windv - tot_windv_max) / (tot_windn - 1);

	  fprintf(FdWinds, "%g %i %g %g %g\n", All.Time, tot_windn, tot_windv_min, tot_windv_max, tot_windv);
	  fflush(FdWinds);
	}
    }
#endif

  tend = second();
/*   CPU_sn_info = timediff(tstart, tend); */

/*   MPI_Reduce(&CPU_sn_info, &sumCPU_sn_info, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD); */
/*   if(ThisTask == 0) */
/*     All.CPU_SN_info += sumCPU_sn_info / NTask; */

#endif /* closes LT_SEv_INFO  */

  //MPI_Reduce(&CPU_eff_iter, &sumCPU_eff_iter, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  //if(ThisTask == 0)
  //All.CPU_Eff_Iter += sumCPU_eff_iter / NTask;

}

double get_starformation_rate(int i, float *Temperature, float *xclouds)
{
  double rateOfSF;
  double a3inv, dt;
  int mySFi, myIMFi, flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;
  double Z;
  double myFactorEVP, myPhysDensThresh, myEgySpecSN, myFactorSN, myMaxSfrTimescale;
  double temp;

#ifdef LT_METAL_COOLING_on_SMOOTH_Z
  double metalmass;
#endif

#ifdef UM_CHEMISTRY
  double u;
#endif

#ifdef LT_METAL_COOLING_WAL
  double Metallicities[LT_NMet], DZ, Redshift;
#endif

#ifdef WINDS
  if(SphP[i].DelayTime > 0)
    return 0;
#endif

  *xclouds = 0;

  mySFi = get_SF_index(i, &mySFi, &myIMFi);
  myFactorSN = SFs[mySFi].FactorSN;
  myEgySpecSN = SFs[mySFi].EgySpecSN;
  myMaxSfrTimescale = SFs[mySFi].MaxSfrTimescale;


  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;

#ifdef LT_METAL_COOLING_WAL
  if(LT_NMet < get_cool_n_el())
    endrun(2009);
  set_metallicities(i, &Metallicities[0], a3inv);
  if(All.ComovingIntegrationOn)
    {
      Redshift = 1.0 / All.Time - 1;
      get_cool_redshift(Redshift, &DZ);
    }
  else
    {
      Redshift = 0;
      DZ = 0;
    }
#else
#ifndef LT_METAL_COOLING_on_SMOOTH_Z
  Z = get_metallicity_solarunits(get_metallicity(i, Iron));
#else
  double metalmass = get_metalmass(SphP[i].Metals);
  if(metalmass > 0)
    {
#ifdef LT_ZSMOOTH_ALLMETALS
      Z = get_metallicity_solarunits(SphP[i].Zsmooth[Iron]);
#else
      Z = get_metallicity_solarunits(SphP[i].Zsmooth * SphP[i].Metals[Iron] / metalmass);
#endif
    }
  else
    Z = NO_METAL;
#endif
#endif

#ifdef LT_METAL_COOLING_WAL
  if(SFs[mySFi].SFTh_Zdep)
    {
#if defined (UM_METAL_COOLING)
      um_ZsPoint = SphP[i].Metals;
      um_mass = P[i].Mass;
#endif

      getindex(&CoolZvalue[0], 0, ZBins - 1, &Z, &flag);
      if(flag == 0 || flag == ZBins - 1)
	{
	  myFactorEVP = SFs[mySFi].FEVP[flag];
	  myPhysDensThresh = SFs[mySFi].PhysDensThresh[flag];
	}
      else
	{
	  x = (Z - CoolZvalue[flag + 1]) / (CoolZvalue[flag + 1] - CoolZvalue[flag]);
	  myFactorEVP = SFs[mySFi].FEVP[flag] * (1 - x) + SFs[mySFi].FEVP[flag + 1] * x;
	  myPhysDensThresh =
	    SFs[mySFi].PhysDensThresh[flag] * (1 - x) + SFs[mySFi].PhysDensThresh[flag + 1] * x;
	}
    }
  else
#endif
    {
      myFactorEVP = SFs[mySFi].FEVP[0];
      myPhysDensThresh = SFs[mySFi].PhysDensThresh[0];

#ifdef UM_MYTH
      myPhysDensThresh = UM_MYTH;	//1.e20,  0.480065; //<-new th
#endif
    }

  flag = 1;			/* default is normal cooling */

  if(SphP[i].d.Density * a3inv >= myPhysDensThresh)
    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].d.Density < All.OverDensThresh)
      flag = 1;

#ifdef LT_TEMP_THRESH_FOR_MULTIPHASE
  ignore_failure_in_convert_u = 1;
#else
  if(flag == 1)
    {
#endif

      ne = SphP[i].elec;
#ifdef UM_CHEMISTRY
      Um_Compute_MeanMolecularWeight(i);
      u = DMAX(All.MinEgySpec, SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1));
      *Temperature =
	GAMMA_MINUS1 / BOLTZMANN * PROTONMASS * SphP[i].Um_MeanMolecularWeight * u * All.UnitPressure_in_cgs /
	All.UnitDensity_in_cgs;
      if(*Temperature >= T_SUP_INTERPOL_LIMIT)
#endif

#ifndef LT_METAL_COOLING_WAL
	*Temperature = convert_u_to_tempSTD(DMAX(All.MinEgySpec,
						 (SphP[i].Entropy +
						  SphP[i].e.DtEntropy * dt) / GAMMA_MINUS1 *
						 pow(SphP[i].d.Density * a3inv,
						     GAMMA_MINUS1)) * All.UnitEnergy_in_cgs /
					    All.UnitMass_in_g,
					    SphP[i].d.Density * All.UnitDensity_in_cgs * All.HubbleParam *
					    All.HubbleParam * a3inv, &ne);
#else
	*Temperature = convert_u_to_tempMET(DMAX(All.MinEgySpec,
						 (SphP[i].Entropy +
						  SphP[i].e.DtEntropy * dt) / GAMMA_MINUS1 *
						 pow(SphP[i].d.Density * a3inv,
						     GAMMA_MINUS1)) * All.UnitEnergy_in_cgs /
					    All.UnitMass_in_g, Redshift, DZ, &Metallicities[0]);
#endif

#ifdef LT_TEMP_THRESH_FOR_MULTIPHASE
      ignore_failure_in_convert_u = 0;
      if(*Temperature > All.MultphaseTempThresh || flag == 1)
	return 0;
#else
      return 0;
    }
#endif

  tsfr = sqrt(myPhysDensThresh / (SphP[i].d.Density * a3inv)) * myMaxSfrTimescale;

  factorEVP = pow(SphP[i].d.Density * a3inv / myPhysDensThresh, -0.8) * myFactorEVP;

  egyhot = myEgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = SphP[i].elec;

#ifndef RT_COOLING_PHOTOHEATING
#ifndef LT_METAL_COOLING_WAL
#ifdef LT_METAL_COOLING
#ifndef UM_CHEMISTRY
  if((tcool = GetCoolingTimeZ(egyhot, SphP[i].d.Density * a3inv, &ne, Z, &temp)) == 0)
    tcool = 1e13;
#else
  if((tcool = Um_GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne, Z, i)) == 0)
    tcool = 1e13;
#endif
#else
  tcool = GetCoolingTimeSTD(egyhot, SphP[i].d.Density * a3inv, &ne);
#endif
#else /* LT_METAL_COOLING_WAL */
  tcool = GetCoolingTimeMET(egyhot, SphP[i].d.Density * a3inv, Redshift, DZ, &Metallicities[0], &temp);
#endif
#else /* RT_COOLING_PHOTOHEATING */
  tcool = rt_GetCoolingTime(i, egyhot, SphP[i].d.Density * a3inv);
#endif

  y = tsfr / tcool * egyhot / (myFactorSN * myEgySpecSN - (1 - myFactorSN) * All.EgySpecCold);

  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  cloudmass = x * P[i].Mass;

  *xclouds = x;

  rateOfSF = (1 - myFactorSN) * cloudmass / tsfr;

  /* convert to solar masses per yr */

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

#ifdef SUBFIND

double get_starformation_rate_subfind(int i, float *Temperature, float *xclouds)
{
  double rateOfSF;
  double a3inv, dt;
  int mySFi, myIMFi, flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;
  double Z;
  double myFactorEVP, myPhysDensThresh, myEgySpecSN, myFactorSN, myMaxSfrTimescale;
  double temp;

#ifdef LT_METAL_COOLING_on_SMOOTH_Z
  double metalmass;
#endif

#ifdef UM_CHEMISTRY
  double u;
#endif

#ifdef LT_METAL_COOLING_WAL
  double Metallicities[LT_NMet], DZ, Redshift;
#endif

#ifdef WINDS
  if(SphP[P[i].origindex].DelayTime > 0)
    return 0;
#endif

  *xclouds = 0;

  mySFi = get_SF_index_subfind(i, &mySFi, &myIMFi);
  myFactorSN = SFs[mySFi].FactorSN;
  myEgySpecSN = SFs[mySFi].EgySpecSN;
  myMaxSfrTimescale = SFs[mySFi].MaxSfrTimescale;


  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;

#ifdef LT_METAL_COOLING_WAL
  if(LT_NMet < get_cool_n_el())
    endrun(3009);
  set_metallicities_subfind(i, &Metallicities[0], a3inv);
  if(All.ComovingIntegrationOn)
    {
      Redshift = 1.0 / All.Time - 1;
      get_cool_redshift(Redshift, &DZ);
    }
  else
    {
      Redshift = 0;
      DZ = 0;
    }
#else
#ifndef LT_METAL_COOLING_on_SMOOTH_Z
  Z = get_metallicity_solarunits(get_metallicity_subfind(i));
#else
  double metalmass = get_metalmass(SphP[P[i].origindex].Metals);
  if(metalmass > 0)
    Z =
      get_metallicity_solarunits(SphP[P[i].origindex].Zsmooth * SphP[P[i].origindex].Metals[Iron] /
				 metalmass);
  else
    Z = NO_METAL;
#endif
#endif

#ifdef LT_METAL_COOLING_WAL
  if(SFs[mySFi].SFTh_Zdep)
    {
#if defined (UM_METAL_COOLING)
      um_ZsPoint = SphP[P[i].origindex].Metals;
      um_mass = P[i].Mass;
#endif

      getindex(&CoolZvalue[0], 0, ZBins - 1, &Z, &flag);
      if(flag == 0 || flag == ZBins - 1)
	{
	  myFactorEVP = SFs[mySFi].FEVP[flag];
	  myPhysDensThresh = SFs[mySFi].PhysDensThresh[flag];
	}
      else
	{
	  x = (Z - CoolZvalue[flag + 1]) / (CoolZvalue[flag + 1] - CoolZvalue[flag]);
	  myFactorEVP = SFs[mySFi].FEVP[flag] * (1 - x) + SFs[mySFi].FEVP[flag + 1] * x;
	  myPhysDensThresh =
	    SFs[mySFi].PhysDensThresh[flag] * (1 - x) + SFs[mySFi].PhysDensThresh[flag + 1] * x;
	}
    }
  else
#endif
    {
      myFactorEVP = SFs[mySFi].FEVP[0];
      myPhysDensThresh = SFs[mySFi].PhysDensThresh[0];

#ifdef UM_MYTH
      myPhysDensThresh = UM_MYTH;	//1.e20,  0.480065; //<-new th
#endif
    }

  flag = 1;			/* default is normal cooling */

  if(SphP[P[i].origindex].d.Density * a3inv >= myPhysDensThresh)
    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[P[i].origindex].d.Density < All.OverDensThresh)
      flag = 1;

#ifdef LT_TEMP_THRESH_FOR_MULTIPHASE
  ignore_failure_in_convert_u = 1;
#else
  if(flag == 1)
    {
#endif

      ne = SphP[P[i].origindex].elec;
#ifdef UM_CHEMISTRY
      Um_Compute_MeanMolecularWeight(i);
      u =
	DMAX(All.MinEgySpec,
	     SphP[P[i].origindex].Entropy / GAMMA_MINUS1 * pow(SphP[P[i].origindex].d.Density * a3inv,
							       GAMMA_MINUS1));
      *Temperature =
	GAMMA_MINUS1 / BOLTZMANN * PROTONMASS * SphP[P[i].origindex].Um_MeanMolecularWeight * u *
	All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
      if(*Temperature >= T_SUP_INTERPOL_LIMIT)
#endif

#ifndef LT_METAL_COOLING_WAL
	*Temperature = convert_u_to_tempSTD(DMAX(All.MinEgySpec,
						 (SphP[P[i].origindex].Entropy +
						  SphP[P[i].origindex].e.DtEntropy * dt) / GAMMA_MINUS1 *
						 pow(SphP[P[i].origindex].d.Density * a3inv,
						     GAMMA_MINUS1)) * All.UnitEnergy_in_cgs /
					    All.UnitMass_in_g,
					    SphP[P[i].origindex].d.Density * All.UnitDensity_in_cgs *
					    All.HubbleParam * All.HubbleParam * a3inv, &ne);
#else
	*Temperature = convert_u_to_tempMET(DMAX(All.MinEgySpec,
						 (SphP[P[i].origindex].Entropy +
						  SphP[P[i].origindex].e.DtEntropy * dt) / GAMMA_MINUS1 *
						 pow(SphP[P[i].origindex].d.Density * a3inv,
						     GAMMA_MINUS1)) * All.UnitEnergy_in_cgs /
					    All.UnitMass_in_g, Redshift, DZ, &Metallicities[0]);
#endif

#ifdef LT_TEMP_THRESH_FOR_MULTIPHASE
      ignore_failure_in_convert_u = 0;
      if(*Temperature > All.MultphaseTempThresh || flag == 1)
	return 0;
#else
      return 0;
    }
#endif

  tsfr = sqrt(myPhysDensThresh / (SphP[P[i].origindex].d.Density * a3inv)) * myMaxSfrTimescale;

  factorEVP = pow(SphP[P[i].origindex].d.Density * a3inv / myPhysDensThresh, -0.8) * myFactorEVP;

  egyhot = myEgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = SphP[P[i].origindex].elec;

#ifndef RT_COOLING_PHOTOHEATING
#ifndef LT_METAL_COOLING_WAL
#ifdef LT_METAL_COOLING
#ifndef UM_CHEMISTRY
  if((tcool = GetCoolingTimeZ(egyhot, SphP[P[i].origindex].d.Density * a3inv, &ne, Z, &temp)) == 0)
    tcool = 1e13;
#else
  if((tcool = Um_GetCoolingTime(egyhot, SphP[P[i].origindex].d.Density * a3inv, &ne, Z, i)) == 0)
    tcool = 1e13;
#endif
#else
  tcool = GetCoolingTimeSTD(egyhot, SphP[P[i].origindex].d.Density * a3inv, &ne);
#endif
#else /* LT_METAL_COOLING_WAL */
  tcool =
    GetCoolingTimeMET(egyhot, SphP[P[i].origindex].d.Density * a3inv, Redshift, DZ, &Metallicities[0], &temp);
#endif
#else /* RT_COOLING_PHOTOHEATING */
  tcool = rt_GetCoolingTime(i, egyhot, SphP[P[i].origindex].d.Density * a3inv);
#endif

  y = tsfr / tcool * egyhot / (myFactorSN * myEgySpecSN - (1 - myFactorSN) * All.EgySpecCold);

  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  cloudmass = x * P[i].Mass;

  *xclouds = x;

  rateOfSF = (1 - myFactorSN) * cloudmass / tsfr;

  /* convert to solar masses per yr */

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

#endif

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

	    if(P[i].Type == 4)
	      MetP[P[i].pt.MetID].PID = j;

#if defined(BLACK_HOLES) && defined(DETACH_BLACK_HOLES)
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

#if defined(DETACH_BLACK_HOLES)
	    if(P[N_gas - 1].Type == 5)
	      BHP[P[N_gas - 1].pt.BHID].PID = N_gas - 1;
#endif

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
#if defined(DETACH_BLACK_HOLES)
	    if(P[i].Type == 5)
	      BHP[P[i].pt.BHID].PID = i;
#endif
	  }
#if defined(DETACH_BLACK_HOLES)
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
#endif
	else
	  {
	    P[i] = P[NumPart - 1];

	    if(P[i].Type == 4)
	      MetP[P[i].pt.MetID].PID = i;
#if defined(DETACH_BLACK_HOLES)
	    if(P[i].Type == 5)
	      BHP[P[i].pt.BHID].PID = i;
#endif
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
#endif // closes GM_MUPPI


int write_eff_model(int num, int mySFi)
{
  char string[200];

  FILE *file;

  int i;

  if(ThisTask == 0)
    {
      if(num < 0)
	sprintf(string, "%s.%02d", "eff_model.dat", mySFi);
      else
	sprintf(string, "eff_model.dat.%02d.%03d", mySFi, num);


      /* write data on a file to save time in the next run */
      if((file = fopen(string, "w")) == NULL)
	{
	  printf("it's impossible to write in <%s> \n", string);
	  fclose(file);
	  return -1;
	}
      else
	{
	  fprintf(file, "%d\n", ZBins);
	  for(i = 0; i < ZBins; i++)
	    fprintf(file, "%lg %lg %lg\n", CoolZvalue[i], SFs[mySFi].PhysDensThresh[i], SFs[mySFi].FEVP[i]);
	  fclose(file);
	}
    }
  return 0;
}

int read_eff_model(int num, int mySFi)
{
  int i, j, k, Zbins;

  char string[200];

  FILE *file;

  k = 0;

  if(num < 0)
    sprintf(string, "eff_model.dat.%02d", mySFi);
  else
    sprintf(string, "eff_model.dat.%02d.%03d", mySFi, num);

  if((file = fopen(string, "r")) != 0x0)
    {
      if(fscanf(file, "%d", &Zbins) < 1)
	return k;
      if(Zbins != ZBins)
	{
	  printf("# of bins in %s is different than that found in cooling tables!\n", string);
	  endrun(20002);
	}
      for(i = j = 0; i < Zbins; i++)
	{
	  j += fscanf(file, "%*g %lg %lg", &SFs[mySFi].PhysDensThresh[i], &SFs[mySFi].FEVP[i]);
	  printf("[Fe/H]: %-4.3lg - RhoTh: %-9.7lg  - fEVP: %-9.7lg\n",
		 CoolZvalue[i], SFs[mySFi].PhysDensThresh[i], SFs[mySFi].FEVP[i]);
	}
      if(j != Zbins * 2)
	k = 0;
      else
	k = 1;
    }
  return k;
}

void init_clouds_cm(int mode, double *PhDTh, double *fEVP, double EgySN, double FSN, double SFt,
		    int Zbins, double *ZArray)
{
  int i;

  int *offset, *counts;

  if(mode > 0)
    {
      if(ThisTask == 0)
	printf("\n\ninitialize effective model.. \n"
	       "it is needed to recalculate metallicity dependence for effective model.. it will take a while\n");
      fflush(stdout);
      for(i = 0; i < Zbins; i++)
	PhDTh[i] = 0;
    }

  /* distribute work over processor */
  offset = (int *) mymalloc("offset", sizeof(int) * NTask);
  counts = (int *) mymalloc("counts", sizeof(int) * NTask);

  offset[0] = 0;
  counts[0] = Zbins / NTask + (0 < (int) fmod(Zbins, NTask));
  for(i = 1; i < NTask; i++)
    {
      if((counts[i] = Zbins / NTask + (i < (int) fmod(Zbins, NTask))))
	offset[i] = offset[i - 1] + counts[i - 1];
      else
	offset[i] = 0;
    }

  /* ThisTask will calculate params for its range of Z */
  for(i = 0; i < counts[ThisTask]; i++)
    init_clouds(mode, EgySN, FSN, SFt, ZArray[i + offset[ThisTask]], &PhDTh[i + offset[ThisTask]],
		&fEVP[i + offset[ThisTask]]);

  MPI_Barrier(MYMPI_COMM_WORLD);
  if(mode > 0)
    {
      /* collect informations */
      MPI_Allgatherv((void *) &PhDTh[offset[ThisTask]], counts[ThisTask], MPI_DOUBLE,
		     (void *) &PhDTh[0], counts, offset, MPI_DOUBLE, MYMPI_COMM_WORLD);

      MPI_Allgatherv((void *) &fEVP[offset[ThisTask]], counts[ThisTask], MPI_DOUBLE,
		     (void *) &fEVP[0], counts, offset, MPI_DOUBLE, MYMPI_COMM_WORLD);

      MPI_Barrier(MYMPI_COMM_WORLD);
    }

  myfree(counts);
  myfree(offset);

  return;
}



void init_clouds(int mode, double egySN, double fSN, double SFt, double Z, double *PhDTh, double *fEVP)
{
#ifdef LT_METAL_COOLING
  int Zbin;
#endif
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight, temp;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

#ifdef UM_MET_IN_LT_COOLING
  float dummy_Metals[LT_NMet];

  memset(dummy_Metals, 0, sizeof(float) * LT_NMet);
  um_ZsPoint = dummy_Metals;
#endif


  /*
   * calculating effective model parameters you have different choices
   * in the case you have settled on LT_METAL_COOLING:
   *   (a) make A0 parameter (onset of thermal instability) dependent on
   *       metal cooling
   *   (b) leave A0 parameter independent of metal cooling
   * the following code trigger this choice as stated by compiler directives
   */

#if defined(LT_METAL_COOLING)
  for(Zbin = ZBins - 1; Z < CoolZvalue[Zbin] && Zbin > 0; Zbin--)
    ;
#endif

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  A0 = *fEVP;

#ifndef LT_METAL_COOLING_WAL
  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;		/* to be guaranteed to get z=0 rate */
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }
#endif

  if(mode)
    {

      egyhot = egySN / A0;

      if(All.ComovingIntegrationOn)
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);


#ifndef LT_METAL_COOLING_WAL

      ne = 1.0;
      SetZeroIonization();

#ifdef LT_METAL_COOLING
/* #ifndef UM_CHEMISTRY */
/*       if((tcool = GetCoolingTimeZ(egyhot, dens, &ne, Z)) == 0) */
/* 	tcool = 1e13; */
/* #else */
/*       /\* what about Um_GetCoolingTime() here ???? I don't know!!! *\/ */
/* #endif */
      if((tcool = GetCoolingTimeZ(egyhot, dens, &ne, Z, &temp)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTimeSTD(egyhot, dens, &ne);
#endif
#else
      if((tcool = GetCoolingTimeMET(egyhot, dens, Redshift, DZ, &Metallicities[0], &temp)) == 0)
	tcool = 1e13;
#endif

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      *PhDTh = x / pow(1 - x, 2) * (fSN * egySN - (1 - fSN) * All.EgySpecCold) / (SFt * coolrate);

      if(mode == 2)
	return;
    }

  dens = *PhDTh * 10;

  do
    {


      tsfr = sqrt(*PhDTh / (dens)) * SFt;
      factorEVP = pow(dens / *PhDTh, -0.8) * *fEVP;
      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

      ne = 0.5;

#ifndef LT_METAL_COOLING_WAL
#ifdef LT_METAL_COOLING
/* #ifndef UM_CHEMISTRY */
/*       if((tcool = GetCoolingTimeZ(egyhot, dens, &ne, Z)) == 0) */
/* 	tcool = 1e13; */
/* #else */
/*       /\* what about Um_GetCoolingTime() here ???? I don't know!!! *\/ */
/* #endif */
      if((tcool = GetCoolingTimeZ(egyhot, dens, &ne, Z, &temp)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTimeSTD(egyhot, dens, &ne);
#endif
#else
      tcool = GetCoolingTimeMET(egyhot, dens, Redshift, DZ, &Metallicities[0], &temp);
#endif


      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
      if(y < 1e-4)
	x = y;
      else
	x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      peff = GAMMA_MINUS1 * dens * egyeff;
      fac = 1 / (log(dens * 1.025) - log(dens));
      dens *= 1.025;
      neff = -log(peff) * fac;

      tsfr = sqrt(*PhDTh / (dens)) * SFt;
      factorEVP = pow(dens / *PhDTh, -0.8) * *fEVP;
      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

      ne = 0.5;

#ifndef LT_METAL_COOLING_WAL
#ifdef LT_METAL_COOLING
/* #ifndef UM_CHEMISTRY */
/*       if((tcool = GetCoolingTimeZ(egyhot, dens, &ne, Z)) == 0) */
/* 	tcool = 1e13; */
/* #else */
/*       /\* what about Um_GetCoolingTime() here ???? I don't know!!! *\/ */
/* #endif */
      if((tcool = GetCoolingTimeZ(egyhot, dens, &ne, Z, &temp)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTimeSTD(egyhot, dens, &ne);
#endif
#else
      tcool = GetCoolingTimeMET(egyhot, dens, Redshift, DZ, &Metallicities[0], &temp);
#endif

      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
      if(y < 1e-4)
	x = y;
      else
	x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      peff = GAMMA_MINUS1 * dens * egyeff;

      neff += log(peff) * fac;
    }
  while(neff > 4.0 / 3);

  thresholdStarburst = dens;
  integrate_sfr(*PhDTh, *fEVP, egySN, fSN, SFt, Z);

  sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

  printf("\n[%i] .: [Fe] = %.3g :. \n"
	 "\nA0= %g  \n"
	 "Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n"
	 "EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n"
	 "tcool=%g dens=%g egyhot=%g\n"
	 "Run-away sets in for dens=%g\n"
	 "Dynamic range for quiescent star formation= %g\n\n"
	 "Isotherm sheet central density: %g   z0=%g\n", ThisTask,
#ifdef LT_METAL_COOLING
	 CoolZvalue[Zbin],
#else
	 (double) 0,
#endif
	 A0,
	 *PhDTh,
	 *PhDTh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs), x, tcool, dens,
	 egyhot, thresholdStarburst, thresholdStarburst / *PhDTh,
	 M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
	 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
  fflush(stdout);

#ifndef LT_METAL_COOLING_WAL
  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }
#endif

  return;
}

void integrate_sfr(double PhDTh, double fEVP, double egySN, double fSN, double SFt, double Z)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, y, P, P2, x2, y2, tsfr2, factorEVP2, egyhot2, tcool2, drho, dq;
  double meanweight, u4, z, tsfr, tcool, egyhot, factorEVP, egyeff, egyeff2, temp;
  char buff[30];
  FILE *fd;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#ifndef LT_METAL_COOLING_WAL
  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;		/* to be guaranteed to get z=0 rate */
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }
#endif

  if(ThisTask == 0)
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = PhDTh; rho <= 1000 * PhDTh; rho *= 1.1)
    {
      tsfr = sqrt(PhDTh / rho) * SFt;

      factorEVP = pow(rho / PhDTh, -0.8) * fEVP;

      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

      ne = 1.0;

#ifndef LT_METAL_COOLING_WAL
#ifdef LT_METAL_COOLING
/* #ifndef UM_CHEMISTRY */
/*       if((tcool = GetCoolingTimeZ(egyhot, rho, &ne, Z)) == 0) */
/* 	tcool = 1e13; */
/* #else */
/*       /\* what about Um_GetCoolingTime() here ???? I don't know!!! *\/ */
/* #endif */
      if((tcool = GetCoolingTimeZ(egyhot, rho, &ne, Z, &temp)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTimeSTD(egyhot, rho, &ne);
#endif
#else
      tcool = GetCoolingTimeMET(egyhot, rho, Redshift, DZ, &Metallicities[0], &temp);
#endif
      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
      if(y < 1e-4)
	x = y;
      else
	x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      P = GAMMA_MINUS1 * rho * egyeff;

      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g\n", rho, P);
	}
    }

  if(ThisTask == 0)
    fclose(fd);


  sprintf(buff, "sfrrate.%03d.txt", sfrrate_filenum);
  if(ThisTask == 0)
    fd = fopen(buff, "w");
  else
    fd = 0;

  for(rho0 = PhDTh; rho0 <= 10000 * PhDTh; rho0 *= 1.02)
    {
      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
	{
	  if(rho > PhDTh)
	    {
	      tsfr = sqrt(PhDTh / rho) * SFt;

	      factorEVP = pow(rho / PhDTh, -0.8) * fEVP;

	      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

	      ne = 1.0;

#ifndef LT_METAL_COOLING_WAL
#ifdef LT_METAL_COOLING
/* #ifndef UM_CHEMISTRY */
/*               if((tcool = GetCoolingTimeZ(egyhot, rho, &ne, Z)) == 0) */
/*                 tcool = 1e13; */
/* #else */
/*       /\* what about Um_GetCoolingTime() here ???? I don't know!!! *\/ */
/* #endif */
	      if((tcool = GetCoolingTimeZ(egyhot, rho, &ne, Z, &temp)) == 0)
		tcool = 1e13;
#else
	      tcool = GetCoolingTimeSTD(egyhot, rho, &ne);
#endif
#else
	      tcool = GetCoolingTimeMET(egyhot, rho, Redshift, DZ, &Metallicities[0], &temp);
#endif

	      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
	      if(y < 1e-4)
		x = y;
	      else
		x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      P = P1 = GAMMA_MINUS1 * rho * egyeff;

	      rho2 = 1.1 * rho;
	      tsfr2 = sqrt(PhDTh / rho2) * SFt;
	      factorEVP2 = pow(rho2 / PhDTh, -0.8) * fEVP;
	      egyhot2 = egySN / (1 + factorEVP2) + All.EgySpecCold;

#ifndef LT_METAL_COOLING_WAL
#ifdef LT_METAL_COOLING
/* #ifndef UM_CHEMISTRY */
/* 	      if((tcool2 = GetCoolingTimeZ(egyhot2, rho2, &ne, Z)) == 0) */
/* 		tcool2 = 1e13; */
/* #else */
/*       /\* what about Um_GetCoolingTime() here ???? I don't know!!! *\/ */
/* #endif */
	      if((tcool2 = GetCoolingTimeZ(egyhot2, rho2, &ne, Z, &temp)) == 0)
		tcool2 = 1e13;
#else
	      tcool2 = GetCoolingTimeSTD(egyhot2, rho2, &ne);
#endif
#else
	      tcool2 = GetCoolingTimeMET(egyhot2, rho2, Redshift, DZ, &Metallicities[0], &temp);
#endif
	      y2 = tsfr2 / tcool2 * egyhot2 / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
	      if(y2 < 1e-4)
		x2 = y2;
	      else
		x2 = 1 + 1 / (2 * y2) - sqrt(1 / y2 + 1 / (4 * y2 * y2));
	      egyeff2 = egyhot2 * (1 - x2) + All.EgySpecCold * x2;
	      P2 = GAMMA_MINUS1 * rho2 * egyeff2;

	      gam = log(P2 / P1) / log(rho2 / rho);
	    }
	  else
	    {
	      tsfr = 0;

	      P = GAMMA_MINUS1 * rho * u4;
	      gam = 1.0;


	      sigma_u4 += rho * dz;
	    }



	  drho = q;
	  dq = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

	  sigma += rho * dz;
	  if(tsfr > 0)
	    {
	      sigmasfr += (1 - SFs[sfrrate_filenum].FactorSN) * rho * x / tsfr * dz;
	    }

	  rho += drho * dz;
	  q += dq * dz;
	}


      sigma *= 2;		/* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;


      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
	}
    }


#ifndef LT_METAL_COOLING_WAL
  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }
#endif

  if(ThisTask == 0)
    fclose(fd);
}

#endif


#ifdef UM_CHECK
void Um_cooling_check(void)
{
  int i, ifunc;
  double a_start, a_end, um_u;

  for(i = 0; i < N_gas; i++)
    {
      a_start = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
      a_end = All.TimeBegin * exp(P[i].Ti_endstep * All.Timebase_interval);

      if(ThisTask == 0 && i == 1)
	{
	  printf("--- in  Um_cooling_check: start = %g, end = %g\n", a_start, a_end);
	  printf
	    ("--- From  Um_cooling_check: Step %d, Time: %g, Systemstep: %g, Dloga=log(Time-Systemstep): %g\n",
	     All.NumCurrentTiStep, All.Time, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	}

      ifunc = compute_abundances(1, i, a_start, a_end, &um_u);
    }
}
#endif

#endif
