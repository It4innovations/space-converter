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

#if !defined(LT_STELLAREVOLUTION) && !defined(GM_MUPPI)
#ifdef COOLING

/*
 * This routine does cooling and star formation for
 * the effective multi-phase model.
 */

#ifdef SFR

double cold_function(double p, double delta)
{
  float g;
  g = 1 + sqrt(p) / 2 / delta - sqrt(sqrt(p) / delta + p / (4 * delta * delta)) - p;
  return g;
}


double find_root_x(double delta)
{
  double a, b, fa, fb, xx, yy, root;

  a = 0.;
  b = 1.;
  int mycount = 0;
  fa = cold_function(a, delta);
  fb = cold_function(b, delta);
  if(fa * fb >= 0)		/* No root, should never happen ! */
    {
      if(fabs(fa) < 1e-3)
	{
	  printf("Task=%d: WARNING, root at beginning of range: a=%.5f, b=%.5f, fa=%.5f, fb=%.5f \n",
		 ThisTask, a, b, fa, fb);
	  b = a;
	}
      else if(fabs(fb) < 1e-3)
	{
	  printf("Task=%d: WARNING, root at end of range: a=%.5f, b=%.5f, fa=%.5f, fb=%.5f \n", ThisTask, a,
		 b, fa, fb);
	  a = b;
	}
      else
	{
	  printf("Task=%d: ERROR, no root possible: a=%.5f, b=%.5f, fa=%.5f, fb=%.5f \n", ThisTask, a, b, fa,
		 fb);
	  //            root = 0.15;
	  fflush(stdout);
	  endrun(8762019);
	}
    }

  while(fabs(a - b) > 1e-5 && mycount <= 1000)
    {
      xx = (a + b) / 2.;
      yy = cold_function(xx, delta);
      if((yy >= 0 && fa >= 0) || (yy <= 0 && fa <= 0))
	{
	  a = xx;
	  fa = yy;
	}
      else
	{
	  b = xx;
	  fb = yy;
	}
      mycount++;
    }
  root = (a + b) / 2.;
  return root;
}

double calculate_dens_threshold(int i)
{
  double ne, tcool, coolrate, egyhot, meanweight, u4, x, A0, PhysDensThresh_cgs, PhysDensThresh,
    PhysDensThresh_old, P4, f_mol4;
  double e_unit, delta, Po, f_mol, t_dyn_int, P_b, rho_gas, Pressure, coolrate_cgs, dens, lambda;
  double f_star = 0.1;

  A0 = All.FactorEVP;
  egyhot = All.EgySpecSN / A0;
  ne = SphP[i].elec;

  dens = 1.0e5 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  tcool = GetCoolingTimeSTD(egyhot, dens, &ne);
  lambda = egyhot / tcool / dens;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;	// cgs
  P4 = u4 * dens * GAMMA_MINUS1 * All.UnitDensity_in_cgs;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
  Po = 35000 * BOLTZMANN;
  f_mol4 = 1 / (1 + Po / P4);

  x = (egyhot - u4) / (egyhot - All.EgySpecCold);
  t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / 32 / All.G / dens);


  PhysDensThresh =
    x / pow(1 - x,
	    2) * f_star * f_mol4 * (All.FactorSN * All.EgySpecSN - (1 -
								    All.FactorSN) * All.EgySpecCold) /
    t_dyn_int / lambda;

  return PhysDensThresh;
}


void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  unsigned int bits;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, mass_of_star;
  double sum_sm, total_sm, sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p, prob;
  double cloudmass;
  double factorEVP;
  double tsfr, trelax;
  double egyhot, egyeff, egycurrent, tcool, x, y, rate_in_msunperyear;
  double sfrrate, totsfrrate;
  double e_unit, delta, root, Po, f_mol, t_dyn, t_dyn_int, P_b, rho_gas, Pressure;

#ifdef WINDS
  int j;
  double v;
  double norm, dir[3];

#ifdef ISOTROPICWINDS
  double theta, phi;
#endif

#ifdef VARIABLE_WINDS
  double z, r200c, v_esc, c_halo, wind_energy, wind_momentum, wind_mass;
  double rhocrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  if(All.ComovingIntegrationOn)
    z = 1 / All.Time - 1;
  else
    z = 0;

  rhocrit *= (All.Omega0 * pow(1 + z, 3) + (1 - All.Omega0 - All.OmegaLambda) * pow(1 + z, 2) + All.OmegaLambda);	/* physical critical density at redshift z */
#endif
#endif

#ifdef METALS
  double w;
#endif

  double temp, u_to_temp_fac;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

#ifdef MODIFIED_SFR

  double SFRTempThresh;

  SFRTempThresh = 5.0e5 / u_to_temp_fac;

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
    }
  else
    a3inv = ascale = hubble_a = time_hubble_a = 1;



  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;

  for(bits = 0; GENERATIONS > (1 << bits); bits++);

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

	  /* check whether conditions for star formation are fulfilled.
	   *  
	   * f=1  normal cooling
	   * f=0  star formation
	   */
	  flag = 1;		/* default is normal cooling */




#if !defined(MODIFIED_SFR) && !defined(EB_SFR_MAGNETIC)
	  if(SphP[i].Density * a3inv >= All.PhysDensThresh)
	    flag = 0;
#elif defined(MODIFIED_SFR) && !defined(EB_SFR_MAGNETIC)
	  if((SphP[i].Density * a3inv >= All.PhysDensThresh)
	     && (SphP[i].Entropy * pow(SphP[i].Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1 < SFRTempThresh))
	    flag = 0;
#else
	  SphP[i].DensityThreshold = calculate_dens_threshold(i);

	  if(SphP[i].Density * a3inv >= SphP[i].DensityThreshold)
	    flag = 0;
#endif

	  if(All.ComovingIntegrationOn)
	    if(SphP[i].Density < All.OverDensThresh)
	      flag = 1;

#ifdef BLACK_HOLES
	  if(P[i].Mass == 0)
	    flag = 1;
#endif

#ifdef WINDS
	  if(SphP[i].DelayTime > 0)
	    flag = 1;		/* only normal cooling for particles in the wind */

	  if(SphP[i].DelayTime > 0)
	    SphP[i].DelayTime -= dtime;

	  if(SphP[i].DelayTime > 0)
	    {
#ifdef QUICK_LYALPHA
	      if(SphP[i].Density < All.WindFreeTravelDensFac * All.OverDensThresh)
		SphP[i].DelayTime = 0;
#else
	      if(SphP[i].Density * a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh)
		SphP[i].DelayTime = 0;
#endif
	    }
	  if(SphP[i].DelayTime < 0)
	    SphP[i].DelayTime = 0;

#endif


#ifdef QUICK_LYALPHA
	  temp = u_to_temp_fac * (SphP[i].Entropy) /
	    GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);

	  if(SphP[i].Density > All.OverDensThresh && temp < 1.0e5)
	    flag = 0;
	  else
	    flag = 1;
#ifdef WINDS
	  if(SphP[i].DelayTime > 0)
	    flag = 1;		/* only normal cooling for particles in the wind */
#endif
#endif

#ifdef MAGNETIC
	  x = 0.;
#endif

#ifdef MAGNETIC_SN_SEEDING
	  SphP[i].MagSeed[0] = SphP[i].MagSeed[1] = SphP[i].MagSeed[2] = 0.0;
#endif

#ifdef LMB_SPECTRAL_CRs_SN_SEEDING
	  SphP[i].CRsSNe = 0.0;
#endif

#ifdef EB_SFR_MAGNETIC
	  e_unit = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);
	  Po = 35000 * BOLTZMANN;
#ifdef MAGNETIC
	  P_b =
	    sqrt(SphP[i].BPred[0] * SphP[i].BPred[0] + SphP[i].BPred[1] * SphP[i].BPred[1] +
		 SphP[i].BPred[2] * SphP[i].BPred[2]) / (8 * M_PI);
#else
	  P_b = 0;
#endif
	  Pressure =
	    SphP[i].Pressure * e_unit * All.HubbleParam * All.HubbleParam / pow(All.UnitLength_in_cm,
										3) + P_b;
	  f_mol = 1 / (1 + Po / Pressure);

#endif


#if !defined(NOISMPRESSURE) && !defined(QUICK_LYALPHA)
	  if(flag == 1)		/* normal implicit isochoric cooling */
#endif
	    {
	      SphP[i].Sfr = 0;


	      ne = SphP[i].elec;	/* electron abundance (gives ionization state and mean molecular weight) */

	      double fac_entr_to_u = pow(SphP[i].Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
	      double uold = DMAX(All.MinEgySpec, SphP[i].Entropy * fac_entr_to_u);

#ifndef RT_COOLING_PHOTOHEATING
	      unew = DoCoolingSTD(uold, SphP[i].Density * a3inv, dtime, &ne);
#else
	      unew = uold + dt * fac_entr_to_u * (rt_DoHeating(i, dt) + rt_DoCooling(i, dt));
#endif


#ifdef UM_CHEMISTRY
	      unew = Um_DoCooling(uold, SphP[i].Density * a3inv, dtime, &ne, i, flag);
#endif



#ifdef BLACK_HOLES
	      if(All.BlackHoleThermalFeedbackOn == 1 || All.BlackHoleKineticFeedbackOn == 1)
		if(SphP[i].Injected_BH_Energy)
		  {
		    if(P[i].Mass == 0)
		      SphP[i].Injected_BH_Energy = 0;
		    else
		      unew += SphP[i].Injected_BH_Energy / P[i].Mass;

		    temp = u_to_temp_fac * unew;

		    if(temp > 5.0e9)
		      unew = 5.0e9 / u_to_temp_fac;

		    SphP[i].Injected_BH_Energy = 0;
		  }
#endif

#if GADGET_HYDRO == HYDRO_MFM
	      SphP[i].InternalEnergy = unew;
#endif
	      SphP[i].Entropy = unew / fac_entr_to_u;
	      SphP[i].elec = ne;

#if GADGET_HYDRO == HYDRO_MFM
	      SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
	      mfm_synchronize_entropy(i);
#endif
	      SphP[i].EntropyPred = SphP[i].Entropy;
	      SphP[i].Pressure = get_pressure(i);
	    }


	  if(flag == 0)		/* active star formation */
	    {
#if !defined(QUICK_LYALPHA)
#ifndef EB_SFR_MAGNETIC
	      factorEVP = pow(SphP[i].Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;
#else
	      factorEVP = pow(SphP[i].Density * a3inv / SphP[i].DensityThreshold, -0.8) * All.FactorEVP;
#endif
	      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = SphP[i].elec;
#ifndef UM_CHEMISTRY
	      tcool = GetCoolingTimeSTD(egyhot, SphP[i].Density * a3inv, &ne);
#else
	      tcool = Um_GetCoolingTime(egyhot, SphP[i].Density * a3inv, &ne, i);
#endif

	      SphP[i].elec = ne;

#ifndef EB_SFR_MAGNETIC
	      tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * a3inv)) * All.MaxSfrTimescale;
	      y =
		tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
#else
	      double f_star = 0.1;
	      delta =
		sqrt(3 * M_PI / (32 * All.G * SphP[i].Density * a3inv)) * (egyhot / tcool) / (f_star * f_mol *
											      (All.FactorSN *
											       All.EgySpecSN -
											       (1 -
												All.
												FactorSN) *
											       All.
											       EgySpecCold));
	      x = find_root_x(delta);
	      t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / (32 * All.G * SphP[i].Density * a3inv));
#endif

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

#ifdef MAGNETIC_SFR_EGY
	      egyeff +=
		(SphP[i].BPred[0] * SphP[i].BPred[0] + SphP[i].BPred[1] * SphP[i].BPred[1] +
		 SphP[i].BPred[2] * SphP[i].BPred[2]) / (2. * MU0) * x;
#endif

	      cloudmass = x * P[i].Mass;

	      if(tsfr < dtime)
		tsfr = dtime;
#ifndef EB_SFR_MAGNETIC
	      sm = (1 - All.FactorSN) * dtime / tsfr * cloudmass;	/* amount of stars expect to form */
#else
	      sm = (1 - All.FactorSN) * f_star * dtime / t_dyn_int * cloudmass * f_mol;
#endif
	      p = sm / P[i].Mass;

	      sum_sm += P[i].Mass * (1 - exp(-p));


	      if(dt > 0)
		{
		  if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
		    {
#ifndef EB_SFR_MAGNETIC
		      trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
#else
		      trelax = t_dyn_int * (1 - x) / x / f_mol / f_star / (All.FactorSN * (1 + factorEVP));
#endif
		      egycurrent =
			SphP[i].Entropy * pow(SphP[i].Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;

#ifdef BLACK_HOLES
		      if(All.BlackHoleThermalFeedbackOn == 1 || All.BlackHoleKineticFeedbackOn == 1)
			if(SphP[i].Injected_BH_Energy > 0)
			  {
			    egycurrent += SphP[i].Injected_BH_Energy / P[i].Mass;

			    temp = u_to_temp_fac * egycurrent;

			    if(temp > 5.0e9)
			      egycurrent = 5.0e9 / u_to_temp_fac;

			    if(egycurrent > egyeff)
			      {
				tcool = GetCoolingTimeSTD(egycurrent, SphP[i].Density * a3inv, &ne);

				if(tcool < trelax && tcool > 0)
				  trelax = tcool;
			      }

			    SphP[i].Injected_BH_Energy = 0;
			  }
#endif

#ifdef MAGNETIC_SN_SEEDING
		      double gauss2gadget =
			All.Time * All.Time * sqrt(All.UnitTime_in_s * All.UnitTime_in_s *
						   All.UnitLength_in_cm / All.UnitMass_in_g);
		      if(All.ComovingIntegrationOn)
			gauss2gadget /= All.HubbleParam;

		      double mag_seed = All.SnSeedField * pow((All.SnSeedRadius / All.SnSeedBubble), 2)
			* pow((All.SnSeedBubble / P[i].Hsml), 3);
		      double mag_norm = pow(P[i].Hsml, 3) * sqrt(0.5 * pow(All.SnSeedSoftening,
									   3) * (1.0 +
										 pow(All.SnSeedSoftening,
										     3)));

		      double n_acc[3];
		      n_acc[0] = P[i].GravAccel[0] + SphP[i].HydroAccel[0];
		      n_acc[1] = P[i].GravAccel[1] + SphP[i].HydroAccel[1];
		      n_acc[2] = P[i].GravAccel[2] + SphP[i].HydroAccel[2];

		      double n_abs = sqrt(n_acc[0] * n_acc[0] + n_acc[1] * n_acc[1] + n_acc[2] * n_acc[2]);
#ifndef EB_SFR_MAGNETIC
		      double sn_eff = sqrt(0.008 * cloudmass * All.UnitMass_in_g / SOLAR_MASS * dtime / tsfr);
#else
		      double sn_eff =
			sqrt(0.008 * f_mol * f_star * cloudmass * All.UnitMass_in_g / SOLAR_MASS * dtime /
			     t_dyn_int);
#endif

		      if(sn_eff < 1.0)
			{
			  if(get_random_number(P[i].ID + 6) <= sn_eff * sn_eff)
			    sn_eff = 1.0;
			  else
			    sn_eff = 0.0;
			}

		      double magseed = mag_norm * sn_eff * mag_seed * gauss2gadget / dtime;

		      SphP[i].MagSeed[0] = magseed * n_acc[0] / n_abs;
		      SphP[i].MagSeed[1] = magseed * n_acc[1] / n_abs;
		      SphP[i].MagSeed[2] = magseed * n_acc[2] / n_abs;

#ifdef LMB_SPECTRAL_CRs_SN_SEEDING

		      /*
		         #JPO_is_never_wrong      #U$$$ requested this
		         Inject 10% of SN energy into CRs.
		         Needs to be converted from erg to UnitEnergy per Mass.
		       */
#ifndef MAGNETIC_SN_SEEDING
		      double sn_eff;
#endif
#ifndef EB_SFR_MAGNETIC
		      sn_eff = sqrt(0.008 * cloudmass * All.UnitMass_in_g / SOLAR_MASS * dtime / tsfr);
#else
		      sn_eff =
			sqrt(0.008 * f_mol * f_star * cloudmass * All.UnitMass_in_g / SOLAR_MASS * dtime /
			     t_dyn_int);
#endif
		      if(sn_eff < 1.0)
			{
			  if(get_random_number(P[i].ID + 6) <= sn_eff * sn_eff)
			    sn_eff = 1.0;
			  else
			    sn_eff = 0.0;
			}

		      // energy to be injected
		      double SN_Energy = 1.e50 / All.UnitEnergy_in_cgs / P[i].Mass;	// All.EgySpecSN;
		      // compute SN seed energy and adiabatically expand over the kernel
		      double cr_seed =
			All.FactorSN * SN_Energy * pow((All.SnSeedRadius / All.SnSeedBubble),
						       2) * pow((All.SnSeedBubble / P[i].Hsml), 3);

		      /// save to particle -> to be injected
			  SphP[i].CRsSNe = cr_seed * sn_eff;
#endif // LMB_SPECTRAL_CRs_SN_SEEDING

#endif // MAGNETIC_SN_SEEDING

#if !defined(NOISMPRESSURE)

#ifdef ALTERNATIVE_SFR
		      double dens = SphP[i].Density;
		      if(egycurrent > egyeff)
			{
			  ne = SphP[i].elec;
			  unew = DoCoolingSTD(egycurrent, dens * a3inv, dtime, &ne);

			  if(unew < egyeff)
			    {
#ifdef SLOW_RELAX_TO_EOS
			      tcool = GetCoolingTimeSTD(egycurrent, dens * a3inv, &ne);
			      if(tcool < trelax && tcool > 0)
				trelax = tcool;

			      unew = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));
#else
			      unew = egyeff;
#endif
			    }
			}
		      else
			unew = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));

#if GADGET_HYDRO == HYDRO_MFM
		      SphP[i].InternalEnergy = unew;
#endif
		      SphP[i].Entropy = unew * GAMMA_MINUS1 / pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
#else
#if GADGET_HYDRO == HYDRO_MFM
		      SphP[i].InternalEnergy = egyeff + (egycurrent - egyeff) * exp(-dtime / trelax);
#endif
		      SphP[i].Entropy =
			(egyeff +
			 (egycurrent -
			  egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 /
			pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
#endif
#if GADGET_HYDRO == HYDRO_MFM
		      SphP[i].dQdt[NUMDIMS+1] =
			eos->GetDtQFromDtEntropy((MyFloat) 0.0, SphP[i].Entropy, SphP[i].InternalEnergy,
						 P[i].Mass, SphP[i].DivVel, All.cf_hubble_a);
		      SphP[i].DtInternalEnergy =
			eos->GetDtUFromDtEntropy((MyFloat) 0.0, SphP[i].Entropy, SphP[i].InternalEnergy,
						 SphP[i].DivVel, All.cf_hubble_a, All.ComovingIntegrationOn);
#endif
		      SphP[i].DtEntropy = 0;
#endif

#if GADGET_HYDRO == HYDRO_MFM
		      SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
		      mfm_synchronize_entropy(i);
#endif
		      SphP[i].EntropyPred = SphP[i].Entropy;
		      SphP[i].Pressure = get_pressure(i);
		    }
		}



	      /* the upper bits of the gas particle ID store how man stars this gas
	         particle gas already generated */

	      if(bits == 0)
		number_of_stars_generated = 0;
	      else
		number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));

	      mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);


#ifndef EB_SFR_MAGNETIC
	      SphP[i].Sfr = (1 - All.FactorSN) * cloudmass / tsfr *
		(All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
#else
	      SphP[i].Sfr = (1 - All.FactorSN) * f_star * f_mol * cloudmass / t_dyn_int *
		(All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

#endif

	      TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
#ifdef METALS
	      w = get_random_number(P[i].ID);
	      P[i].Metallicity += w * METAL_YIELD * (1 - exp(-p));
#endif

	      prob = P[i].Mass / mass_of_star * (1 - exp(-p));
#else /* belongs to ifndef(QUICK_LYALPHA) */

	      prob = 2.0;	/* this will always cause a star creation event */
	      if(bits == 0)
		number_of_stars_generated = 0;
	      else
		number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));

	      mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);

	      SphP[i].Sfr = 0;

#ifdef QUICK_LYALPHA_FEEDBACK
	      egyhot = All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) * mass_of_star;
	      sm = mass_of_star;
#endif

#endif /* ends to QUICK_LYALPHA */

#ifdef DO_NOT_CREATE_STAR_PARTICLES
	      prob = -1;
#endif

	      if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
		{
		  if(number_of_stars_generated == (GENERATIONS - 1))
		    {
		      /* here we turn the gas particle itself into a star */
		      Stars_converted++;
		      stars_converted++;

		      sum_mass_stars += P[i].Mass;

		      P[i].Type = 4;
		      TimeBinCountSph[P[i].TimeBin]--;
		      TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;

#ifdef STELLARAGE
		      MPP(i).StellarAge = All.Time;
#endif
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

		      P[NumPart + stars_spawned] = P[i];
		      P[NumPart + stars_spawned].Type = 4;
#ifdef SNIA_HEATING
		      PPP[NumPart + stars_spawned].Hsml = All.SofteningTable[0];
#endif

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

		      P[i].ID += ((MyIDType) 1 << (sizeof(MyIDType) * 8 - bits));

		      P[NumPart + stars_spawned].Mass = mass_of_star;
#if GADGET_HYDRO == HYDRO_MFM
		      mfm_add_mass(i, -P[NumPart + stars_spawned].Mass);
#endif
		      P[i].Mass -= P[NumPart + stars_spawned].Mass;
		      sum_mass_stars += P[NumPart + stars_spawned].Mass;
#ifdef STELLARAGE
		      MPP(NumPart + stars_spawned).StellarAge = All.Time;
#endif
		      force_add_star_to_tree(i, NumPart + stars_spawned);

		      stars_spawned++;
		    }
		}

#ifdef METALS
	      if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		P[i].Metallicity += (1 - w) * METAL_YIELD * (1 - exp(-p));
#endif



#ifdef WINDS
	      /* Here comes the wind model */

	      if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		{
#ifdef VARIABLE_WINDS
		  if(SphP[i].HostHaloMass > 0 && sm > 0)
		    {
		      r200c = pow(SphP[i].HostHaloMass / (4 * M_PI / 3.0 * 200 * rhocrit), 1.0 / 3.0);	/* physical r_200,crit value, assuming FoF mass = M_200,crit */
		      v_esc = sqrt(All.G * SphP[i].HostHaloMass / r200c);	/* physical circular velocity at r_200,crit */
		      c_halo =
			All.HaloConcentrationNorm * pow(SphP[i].HostHaloMass, All.HaloConcentrationSlope);
		      v_esc *= sqrt(2 * c_halo / (log(1 + c_halo) - c_halo / (1 + c_halo)));	/* physical escape velocity of halo */
		      v = All.VariableWindVelFactor * v_esc;	/* physical wind velocity */

		      wind_momentum = sm * All.VariableWindSpecMomentum;
		      wind_energy =
			sm * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN);

		      wind_mass =
			(wind_energy +
			 sqrt(wind_energy * wind_energy + v * v * wind_momentum * wind_momentum)) / (v * v);
		      /* wind mass for this particle, assuming the wind is first given the energy wind_energy and then the momentum wind_momentum */
		      p = wind_mass / P[i].Mass;

		      /*printf("Mfof %f r200_c %f c_halo %f v_esc %f v %f sm %e wind_mass %e ml %f \n", SphP[i].HostHaloMass, r200c, c_halo, v_esc, v, sm, wind_mass, wind_mass/sm); */
		    }
		  else
		    {
		      v = 0;
		      p = 0;
		    }
#else
		  p = All.WindEfficiency * sm / P[i].Mass;
#endif

		  prob = 1 - exp(-p);

		  if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
		    {
#ifndef VARIABLE_WINDS
		      v =
			sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			     All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency);
#endif
#ifdef ISOTROPICWINDS
		      theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
		      phi = 2 * M_PI * get_random_number(P[i].ID + 4);

		      dir[0] = sin(theta) * cos(phi);
		      dir[1] = sin(theta) * sin(phi);
		      dir[2] = cos(theta);
#else
		      dir[0] = P[i].GravAccel[1] * P[i].Vel[2] - P[i].GravAccel[2] * P[i].Vel[1];
		      dir[1] = P[i].GravAccel[2] * P[i].Vel[0] - P[i].GravAccel[0] * P[i].Vel[2];
		      dir[2] = P[i].GravAccel[0] * P[i].Vel[1] - P[i].GravAccel[1] * P[i].Vel[0];
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

			  SphP[i].DelayTime = All.WindFreeTravelMaxTimeFactor / hubble_a;
			}
#if defined(MAGNETIC)
		      x = 0.;
#endif
		    }
		}
#endif /* WINDS */
	    }
#if defined(MAGNETIC)
	  SphP[i].XColdCloud = x;
#endif
	}			/* End of If Type = 0 */
    }				/* end of main loop over active particles */


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

  MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	rate = total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;

      /* convert to solar masses per yr */

      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear,
	      total_sum_mass_stars);
      fflush(FdSfr);
    }
}

double get_starformation_rate(int i)
{
  double rateOfSF;
  double a3inv;
  int flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;
  double e_unit, delta, a, b, fa, fb, root, xx, yy, Po, f_mol, t_dyn, t_dyn_int, P_b, rho_gas, Pressure;

#ifdef WINDS
  if(SphP[i].DelayTime > 0)
    return 0;
#endif

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  flag = 1;			/* default is normal cooling */


#ifndef EB_SFR_MAGNETIC
  if(SphP[i].Density * a3inv >= All.PhysDensThresh)
    flag = 0;
#else
  if(SphP[i].Density * a3inv >= SphP[i].DensityThreshold)
    flag = 0;
#endif


  if(All.ComovingIntegrationOn)
    if(SphP[i].Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

#ifndef EB_SFR_MAGNETIC
  factorEVP = pow(SphP[i].Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;
#else
  factorEVP = pow(SphP[i].Density * a3inv / SphP[i].DensityThreshold, -0.8) * All.FactorEVP;
#endif


  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = SphP[i].elec;
#ifndef UM_CHEMISTRY
  tcool = GetCoolingTimeSTD(egyhot, SphP[i].Density * a3inv, &ne);
#else
  tcool = Um_GetCoolingTime(egyhot, SphP[i].Density * a3inv, &ne, i);
#endif

#ifndef EB_SFR_MAGNETIC
  tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * a3inv)) * All.MaxSfrTimescale;
  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
#else
  e_unit = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);
  Po = 35000 * BOLTZMANN;
  double f_star = 0.1;
#ifdef MAGNETIC
  P_b =
    sqrt(SphP[i].BPred[0] * SphP[i].BPred[0] + SphP[i].BPred[1] * SphP[i].BPred[1] +
	 SphP[i].BPred[2] * SphP[i].BPred[2]) / (8 * M_PI);
#else
  P_b = 0;
#endif
  Pressure =
    SphP[i].Pressure * e_unit * All.HubbleParam * All.HubbleParam / pow(All.UnitLength_in_cm, 3) + P_b;
  f_mol = 1 / (1 + Po / Pressure);
  delta =
    sqrt(3 * M_PI / (32 * All.G * SphP[i].Density * a3inv)) * (egyhot / tcool) / (f_star * f_mol *
										  (All.FactorSN *
										   All.EgySpecSN - (1 -
												    All.
												    FactorSN)
										   * All.EgySpecCold));
  x = find_root_x(delta);
  t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / (32 * All.G * SphP[i].Density * a3inv));
#endif

  cloudmass = x * P[i].Mass;

#ifndef EB_SFR_MAGNETIC
  rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;
#else
  rateOfSF = (1 - All.FactorSN) * f_star * cloudmass * f_mol / t_dyn_int;
#endif


  /* convert to solar masses per yr */

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}


#endif /* SFR */

#if defined(SFR) || defined(BLACK_HOLES) || defined(WINDTUNNEL)
void rearrange_particle_sequence(void)
{
  int i, j, flag = 0, flag_sum;
  struct particle_data psave;

#if defined(BLACK_HOLES) || defined(WINDTUNNEL)
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
#endif

#ifdef SFR
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
#endif

#if defined(BLACK_HOLES) || defined(WINDTUNNEL)
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

	    if(P[N_gas - 1].Type == 5)
	      BHP[P[N_gas - 1].pt.BHID].PID = N_gas - 1;

	    N_gas--;

	    count_gaselim++;
	  }
	else if(P[i].Type == 5)
	  {
	    j = P[i].pt.BHID;
	    BHP[j] = BHP[N_BHs - 1];
	    P[BHP[j].PID].pt.BHID = j;
	    N_BHs--;

	    P[i] = P[NumPart - 1];

	    if(P[i].Type == 5)
	      BHP[P[i].pt.BHID].PID = i;
	  }
	else
	  {
	    P[i] = P[NumPart - 1];
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
#if defined(BLACK_HOLES)
  All.TotBHs -= tot_elim - tot_gaselim;
#endif
#endif

  MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  if(flag_sum)
    reconstruct_timebins();
}
#endif /* closing of SFR-conditional */



#if defined(SFR)
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;
  double e_unit, delta, a, b, fa, fb, root, xx, yy, Po, f_mol4, t_dyn, t_dyn_int, P_b, rho_gas, P4, f_mol_eff,
    Pb, Po_int;

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
      e_unit = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

#ifndef EB_SFR_MAGNETIC
      dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#else
      dens = 1.0e5 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#endif

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
      tcool = GetCoolingTimeSTD(egyhot, dens, &ne);

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);
#ifndef EB_SFR_MAGNETIC
      All.PhysDensThresh =
	x / pow(1 - x,
		2) * (All.FactorSN * All.EgySpecSN - (1 -
						      All.FactorSN) * All.EgySpecCold) /
	(All.MaxSfrTimescale * coolrate);
#else
      double f_star = 0.1;
      Po = 35000 * BOLTZMANN;
      Po_int = Po / e_unit * pow(All.UnitLength_in_cm, 3) / pow(All.HubbleParam, 2);	// parameter Po in internal units
      P4 = u4 * dens * GAMMA_MINUS1 * e_unit * All.UnitDensity_in_cgs * pow(All.HubbleParam, 2);
      f_mol4 = 1 / (1 + Po / P4);
      t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / 32 / All.G / dens);
      All.PhysDensThresh =
	x / pow(1 - x,
		2) * f_star * f_mol4 * (All.FactorSN * All.EgySpecSN - (1 -
									All.FactorSN) * All.EgySpecCold) /
	t_dyn_int / coolrate;

#endif

      if(ThisTask == 0)
	{
	  printf("\nA0= %g  \n", A0);
	  printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
	  printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
	  printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);
	}

      dens = All.PhysDensThresh * 10;


      do
	{
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
	  tcool = GetCoolingTimeSTD(egyhot, dens, &ne);
#ifndef EB_SFR_MAGNETIC
	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
#else
	  delta =
	    sqrt(3 * M_PI / (32 * All.G * dens)) * (egyhot / tcool) / (f_star * f_mol4 *
								       (All.FactorSN * All.EgySpecSN -
									(1 -
									 All.FactorSN) * All.EgySpecCold));
	  x = find_root_x(delta);
	  t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / (32 * All.G * dens));
#endif

	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  fac = 1 / (log(dens * 1.025) - log(dens));
	  dens *= 1.025;

	  neff = -log(peff) * fac;

	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
	  f_mol_eff = 1 / (1 + Po_int / peff);

	  ne = 0.5;
	  tcool = GetCoolingTimeSTD(egyhot, dens, &ne);
#ifndef EB_SFR_MAGNETIC
	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
#else
	  delta =
	    sqrt(3 * M_PI / (32 * All.G * dens)) * (egyhot / tcool) / (f_star * f_mol_eff *
								       (All.FactorSN * All.EgySpecSN -
									(1 -
									 All.FactorSN) * All.EgySpecCold));
	  x = find_root_x(delta);
	  t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / (32 * All.G * dens));
#endif

	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  neff += log(peff) * fac;
	}
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;

      if(ThisTask == 0)
	{
	  printf("Run-away sets in for dens=%g\n", thresholdStarburst);
	  printf("Dynamic range for quiescent star formation= %g\n", thresholdStarburst / All.PhysDensThresh);
	  fflush(stdout);
	}

      integrate_sfr();

      if(ThisTask == 0)
	{
	  sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

	  printf("Isotherm sheet central density: %g   z0=%g\n",
		 M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
		 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
	  fflush(stdout);

	}

      if(All.ComovingIntegrationOn)
	{
	  All.Time = All.TimeBegin;
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

#ifdef WINDS
      if(All.WindEfficiency > 0)
	if(ThisTask == 0)
	  printf("Windspeed: %g\n",
		 sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) /
		      All.WindEfficiency));
#endif
    }
}


void integrate_sfr(void)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, y, P, P2, x2, y2, tsfr2, factorEVP2, egyhot2, tcool2, drho, dq;
  double meanweight, u4, z, tsfr, tcool, egyhot, factorEVP, egyeff, egyeff2;
  double delta, delta2, a, b, fa, fb, root, xx, yy, Po, f_mol, t_dyn, t_dyn2, t_dyn_int, t_dyn_int2, P_b,
    rho_gas, Pressure, f_mol_eff, Po_int, e_unit, P4, f_mol4, dens;
  FILE *fd;

#ifndef EB_SFR_MAGNETIC
  dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#else
  dens = 1.0e5 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#endif


  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
  Po = 35000 * BOLTZMANN;
  e_unit = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);
  Po_int = Po / e_unit * pow(All.UnitLength_in_cm, 3) / pow(All.HubbleParam, 2);
  P4 = u4 * dens * GAMMA_MINUS1 * e_unit * All.UnitDensity_in_cgs * pow(All.HubbleParam, 2);
  f_mol4 = 1 / (1 + Po / P4);


  double f_star = 0.1;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;		/* to be guaranteed to get z=0 rate */
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(ThisTask == 0)
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = All.PhysDensThresh; rho <= 1000 * All.PhysDensThresh; rho *= 1.1)
    {
      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

      ne = 1.0;
      tcool = GetCoolingTimeSTD(egyhot, rho, &ne);
#ifndef EB_SFR_MAGNETIC
      y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
#else
      delta =
	sqrt(3 * M_PI / (32 * All.G * rho)) * (egyhot / tcool) / (f_star * f_mol4 *
								  (All.FactorSN * All.EgySpecSN -
								   (1 - All.FactorSN) * All.EgySpecCold));
      x = find_root_x(delta);
      t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / (32 * All.G * rho));
#endif

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      P = GAMMA_MINUS1 * rho * egyeff;
      f_mol_eff = 1 / (1 + Po_int / P);

      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g\n", rho, P);
	}
    }

  if(ThisTask == 0)
    fclose(fd);


  if(ThisTask == 0)
    fd = fopen("sfrrate.txt", "w");
  else
    fd = 0;

  for(rho0 = All.PhysDensThresh; rho0 <= 10000 * All.PhysDensThresh; rho0 *= 1.02)
    {
      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
	{
	  if(rho > All.PhysDensThresh)
	    {
	      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

	      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

	      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = 1.0;
	      tcool = GetCoolingTimeSTD(egyhot, rho, &ne);
#ifndef EB_SFR_MAGNETIC
	      y =
		tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
#else
	      delta =
		sqrt(3 * M_PI / (32 * All.G * rho)) * (egyhot / tcool) / (f_star * f_mol_eff *
									  (All.FactorSN * All.EgySpecSN -
									   (1 -
									    All.FactorSN) * All.EgySpecCold));
	      x = find_root_x(delta);
	      t_dyn_int = 1 / sqrt(x) * sqrt(3 * M_PI / (32 * All.G * rho));
#endif
	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      P = P1 = GAMMA_MINUS1 * rho * egyeff;
	      f_mol_eff = 1 / (1 + Po_int / P);

	      rho2 = 1.1 * rho;
	      tsfr2 = sqrt(All.PhysDensThresh / rho2) * All.MaxSfrTimescale;
	      factorEVP2 = pow(rho2 / All.PhysDensThresh, -0.8) * All.FactorEVP;
	      egyhot2 = All.EgySpecSN / (1 + factorEVP2) + All.EgySpecCold;
	      tcool2 = GetCoolingTimeSTD(egyhot2, rho2, &ne);

#ifndef EB_SFR_MAGNETIC
	      y2 =
		tsfr2 / tcool2 * egyhot2 / (All.FactorSN * All.EgySpecSN -
					    (1 - All.FactorSN) * All.EgySpecCold);
	      x2 = 1 + 1 / (2 * y2) - sqrt(1 / y2 + 1 / (4 * y2 * y2));
#else
	      delta2 =
		sqrt(3 * M_PI / (32 * All.G * rho2)) * (egyhot2 / tcool2) / (f_star * f_mol_eff *
									     (All.FactorSN * All.EgySpecSN -
									      (1 -
									       All.FactorSN) *
									      All.EgySpecCold));
	      x2 = find_root_x(delta);
	      t_dyn_int2 = 1 / sqrt(x2) * sqrt(3 * M_PI / (32 * All.G * rho2));
#endif
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
#ifndef EB_SFR_MAGNETIC
	      sigmasfr += (1 - All.FactorSN) * rho * x / tsfr * dz;
#else
	      sigmasfr += (1 - All.FactorSN) * rho * f_star * f_mol * x / t_dyn_int * dz;
#endif
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


  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(ThisTask == 0)
    fclose(fd);
}




#endif /* closing of SFR-conditional */


#if defined(SFR)
void set_units_sfr(void)
{
  double meanweight;

  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

}


#endif /* closes SFR */

#endif /* closes COOLING */

#endif /* closes LT_STELLAREVOLUTION */
