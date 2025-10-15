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
#include "../Hydro/kernel.h"
#include "blackhole.h"
#include "../System/communication.h"

/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

#ifdef BLACK_HOLES

MyAtLeastDouble hubble_a, ascale;

extern void blackhole_feedback_loop();
extern void blackhole_swallow_loop();

void blackhole(void)
{

  VERBOSE(0, "BH: Beginning black-hole calculation");

  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    hubble_a = ascale = 1;

/***********************************************************************************/
/************************** Local computation of values ****************************/
/***********************************************************************************/

  /* Let's first compute the Mdot values */

  for(int n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {
	MyAtLeastDouble dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / hubble_a;
	MyAtLeastDouble mdot = 0;
	MyAtLeastDouble meddington = 0;

	if(All.BlackHoleDetails >= 1)
	  fprintf(FdBlackHolesDetails, "BHGROWTH: idBH=%llu (%2.7f %2.7f %2.7f) time=%g mbh=%g (",
		  (unsigned long long) P[n].ID, (double) P[n].Pos[0], (double) P[n].Pos[1],
		  (double) P[n].Pos[2], All.Time, (double) BPP(n).BH_Mass);

	if(All.BlackHoleHotColdAccretionOn == 0)
	  Get_Bondi_Parameter(n, BPP(n).BH_Density, BPP(n).BH_SurroundingGasVel, BPP(n).BH_Entropy,
			      All.BlackHoleAccretionFactor, dt, &mdot);
	else
	  {
	    if(All.BlackHoleDetails >= 1)
	      fprintf(FdBlackHolesDetails, " rho=%g (cold: ", (double) BPP(n).BH_Density);
	    Get_Bondi_Parameter(n, BPP(n).BH_ColdDensity, BPP(n).BH_SurroundingColdGasVel,
				BPP(n).BH_ColdEntropy, All.BlackHoleColdAccretionFactor, dt, &mdot);

	    MyAtLeastDouble hotmdot;
	    if(All.BlackHoleDetails >= 1)
	      fprintf(FdBlackHolesDetails, " hot: ");
	    Get_Bondi_Parameter(n, BPP(n).BH_HotDensity, BPP(n).BH_SurroundingHotGasVel, BPP(n).BH_HotEntropy,
				All.BlackHoleAccretionFactor, dt, &hotmdot);
	    mdot += hotmdot;
	  }

       #if defined(KD_DYNFRIC) && !defined(AD_DYNFRIC)
        if(All.BlackHoleFrictionForceOn == 1)
          ApplyBlackHoleFriction(n,dt);
       #endif 

	if(All.BlackHoleVariableEfficiencyOn == 2)
	  {
	    meddington = 1.;	//to avoid redefining the definition of BlackHoleTotalFeedbackEfficiency
	    MyAtLeastDouble mdotphys = mdot * All.UnitMass_in_g / All.UnitTime_in_s / SOLAR_MASS * SEC_PER_YEAR;	//Convert in Msun/yr
	    MyAtLeastDouble mbh = BPP(n).BH_Mass * All.UnitMass_in_g / SOLAR_MASS;	//Convert in Msun
	    if(mdot >= All.BlackHolePhysMdotUpperLimit)	//Upper limit for Mdot (in Msun/yr) due to observational constraints.
	      {
		mdot =
		  All.BlackHolePhysMdotUpperLimit / (All.UnitMass_in_g / All.UnitTime_in_s / SOLAR_MASS *
						     SEC_PER_YEAR);
	      }
	  }
	else
	  {
	    meddington = Mdot_Eddington(BPP(n).BH_Mass);
	    if(All.BlackHoleEddingtonFactor > 0 && mdot > All.BlackHoleEddingtonFactor * meddington)
	      mdot = All.BlackHoleEddingtonFactor * meddington;
	  }

	MyAtLeastDouble epsTot = BlackHoleTotalEfficiency(mdot / meddington, BPP(n).BH_Mass);

	BPP(n).BH_Mass += (1. - epsTot) * mdot * dt;

	BPP(n).BH_TotalFeedbackEfficiency = BlackHoleTotalFeedbackEfficiency(mdot / meddington, BPP(n).BH_Mass);	// Here the quantities are passed in internal units.

	BPP(n).BH_Mdot = mdot;

	if(All.BlackHoleDetails >= 1)
	  fprintf(FdBlackHolesDetails, ") MdotEdd=%g EpsTotFeed=%g EpsTot=%g dt=%g mdot=%g Mnow=%g\n",
		  (double) meddington, (double) BPP(n).BH_TotalFeedbackEfficiency, (double) epsTot,
		  (double) dt, (double) BPP(n).BH_Mdot, (double) BPP(n).BH_Mass);

      }

/***********************************************************************************/
/************************** FIRST loop (accreation) ********************************/
/***********************************************************************************/

  blackhole_feedback_loop();

/***********************************************************************************/
/************************** SECOND loop (swallow) **********************************/
/***********************************************************************************/

  blackhole_swallow_loop();

/***********************************************************************************/
/************************** Final computation of values ****************************/
/***********************************************************************************/

  if(All.BlackHoleRepositioningOn)
    for(int n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
      if(P[n].Type == 5)
#ifdef GM_REPOSITION_ON_STARDENSITY_MAX
	blackhole_stardens_reposition(n);
#else
	if(BPP(n).BH_MinPot < 0.5 * BHPOTVALUEINIT)
	  for(int k = 0; k < 3; k++)
	    P[n].Pos[k] = BPP(n).BH_MinPotPos[k];
#endif

  for(int n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  TimeBin_BH_mass[n] = 0;
	  TimeBin_BH_dynamicalmass[n] = 0;
	  TimeBin_BH_Mdot[n] = 0;
	  TimeBin_BH_Medd[n] = 0;
	}
    }

  for(int n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {
	if(BPP(n).BH_accreted_Mass > 0 && P[n].Mass > 0)
	  {
	    if(All.BlackHoleIgnoreMomentum == 0)
	      for(int k = 0; k < 3; k++)
		P[n].Vel[k] = (P[n].Vel[k] * P[n].Mass + BPP(n).BH_accreted_momentum[k]) /
		  (P[n].Mass + BPP(n).BH_accreted_Mass);

	    P[n].Mass += BPP(n).BH_accreted_Mass;
	    BPP(n).BH_Mass += BPP(n).BH_accreted_BHMass;
	  }
	BPP(n).BH_accreted_Mass = 0;
	if(All.BlackHoleIgnoreMomentum == 1)
	  BPP(n).BH_accreted_BHMass = 0;

	int bin = P[n].TimeBin;
	TimeBin_BH_mass[bin] += BPP(n).BH_Mass;
	TimeBin_BH_dynamicalmass[bin] += P[n].Mass;
	TimeBin_BH_Mdot[bin] += BPP(n).BH_Mdot;
	if(BPP(n).BH_Mass > 0)
	  TimeBin_BH_Medd[bin] += BPP(n).BH_Mdot / BPP(n).BH_Mass;
      }

  double mass_holes = 0, mass_real = 0, mdot = 0, mdotedd = 0;
  for(int bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      {
	mass_holes += TimeBin_BH_mass[bin];
	mass_real += TimeBin_BH_dynamicalmass[bin];
	mdot += TimeBin_BH_Mdot[bin];
	mdotedd += TimeBin_BH_Medd[bin];
      }

  double total_mass_holes, total_mass_real, total_mdot, total_mdotedd;
  MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&mdotedd, &total_mdotedd, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      /* convert to solar masses per yr */
      double mdot_in_msun_per_year =
	total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      /* Note that for mass dependent radiative efficiency the statistics about eddington
         ratios is screwed up. Here we just take the efficiency for a typical, 1e8 Msol BH */
      total_mdotedd *= 1.0 / Mdot_Eddington(total_mdot);

      fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
	      All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year,
	      total_mass_real, total_mdotedd);
      fflush(FdBlackHoles);
    }

  if(All.BlackHoleDetails >= 1)
    fflush(FdBlackHolesDetails);

  CPU_Step[CPU_BLACKHOLES] += measure_time();
}

#endif
