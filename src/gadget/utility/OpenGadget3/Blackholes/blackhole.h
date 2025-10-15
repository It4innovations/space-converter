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

#ifndef BLACK_HOLES_INCLUDE
#define BLACK_HOLES_INCLUDE

#ifdef BLACK_HOLES

extern MyAtLeastDouble hubble_a, ascale;

#define BHPOTVALUEINIT 1.0e30

/***********************************************************************************/
/************************** Helper functions for model parameters ******************/
/***********************************************************************************/


static inline MyAtLeastDouble BlackHoleMaximumRadiativeEfficiency(MyFloat bh_mass)
{
  if(All.BlackHoleVariableEfficiencyOn == 2)
    {
      return (MyAtLeastDouble) All.BlackHoleRadiativeEfficiency;
    }
  else if(All.BlackHoleVariableEfficiencyOn == 1)
    {
      MyAtLeastDouble efficiency = All.BlackHoleRadiativeEfficiencyNorm * pow(bh_mass * 100 / All.HubbleParam,
									      All.
									      BlackHoleRadiativeEfficiencySlope);
      if(efficiency > All.BlackHoleRadiativeEfficiencyMax)	// theoretical maximum value (spin=1)
	efficiency = All.BlackHoleRadiativeEfficiencyMax;
      if(efficiency < All.BlackHoleRadiativeEfficiencyMin)
	efficiency = All.BlackHoleRadiativeEfficiencyMin;
      return efficiency;
    }
  else
    return (MyAtLeastDouble) All.BlackHoleRadiativeEfficiency;
}


static inline MyAtLeastDouble BlackHoleOutflowEfficiency(MyFloat f_Edd, MyFloat bh_mass)
{
  if(All.BlackHoleOutflowModelOn == 1)
    {
      if(f_Edd < 0.05)
	{
	  MyAtLeastDouble beta = 0.5;	// slope
	  MyAtLeastDouble gamma = pow(10., -4.0) * pow(0.05, (-2.8431 - beta));
	  return 0.1 - (gamma * 0.1 * pow(f_Edd, beta));
	}
      else
	return (MyAtLeastDouble) pow(10.0, -4.0) * 0.1 * pow(f_Edd, -2.8431);
//The slope 2.8431 was estimated numerically: -5-log(0.05)/(log(0.05))
    }
  else
    return (MyAtLeastDouble) 0;
}


static inline MyAtLeastDouble BlackHoleRadiativeEfficiency(MyFloat f_Edd, MyFloat bh_mass)
{
  if(All.BlackHoleVariableEfficiencyOn == 1)
    {
      if(f_Edd < 0.05)
	{
	  MyAtLeastDouble beta = 0.5;	// slope
	  MyAtLeastDouble gamma = pow(10., -4.0) * pow(0.05, (-2.8431 - beta));
	  return gamma * BlackHoleMaximumRadiativeEfficiency(bh_mass) * pow(f_Edd, beta);
	}
      else
	{
	  return BlackHoleMaximumRadiativeEfficiency(bh_mass) * (1.0 -
								 (pow(10.0, -4.0) * pow(f_Edd, -2.8431)));
	}
    }
  else if(All.BlackHoleVariableEfficiencyOn == 2)
    {
      MyAtLeastDouble mdotphys = f_Edd * All.UnitMass_in_g / All.UnitTime_in_s / SOLAR_MASS * SEC_PER_YEAR;	//Convert in Msun/yr (the quantities are passed in internal units; f_Edd is actually mdot - mdotEdd set to 1 before passed))
      MyAtLeastDouble mbh = bh_mass * All.UnitMass_in_g / SOLAR_MASS;	//Convert in Msun/yr
      MyAtLeastDouble er = All.BlackHoleRadiativeEfficiencyNorm * pow(mdotphys, All.BlackHoleRadiativeEfficiencySlopeMdot) * pow(mbh, All.BlackHoleRadiativeEfficiencySlopeMbh);	//A*mdot^alpha*M^beta;
      if(er >= 0.42)
	{
	  return 0.42;		// theoretical maximum value (spin=0.998)
	}
      else
	{
	  return er;
	}
    }
  else
    return (MyAtLeastDouble) All.BlackHoleRadiativeEfficiency;
}


/* Following Davis&Laor 2010, astro-ph 1012.3213 and Chelouche 2013, astro-ph 1305:6507 */
/* e.g. eps_r or eps_0 + eps_r in Steinborn et al. 2015, eq. (12) */
static inline MyAtLeastDouble BlackHoleTotalEfficiency(MyFloat f_Edd, MyFloat bh_mass)
{
  return BlackHoleOutflowEfficiency(f_Edd, bh_mass) + BlackHoleRadiativeEfficiency(f_Edd, bh_mass);
  // std model: BlackHoleOutflowEfficiency = 0 and BlackHoleRadiativeEfficiency = eps_r  -> returns eps_r
}


/* eps_r * eps_f (where eps_f can be boosted in radio mode) or eps_0 + eps_r * eps_f in Steinborn et al. 2015, eq. (9) */
static inline MyAtLeastDouble BlackHoleTotalFeedbackEfficiency(MyFloat f_Edd, MyFloat bh_mass)
{
  MyAtLeastDouble feedback_factor = All.BlackHoleFeedbackFactor;
  if(All.BlackHoleOutflowModelOn == 0)
    if(All.BlackHoleRadioModeFeedbackBoost > 1)
      if(f_Edd < All.BlackHoleRadioModeTreshold)
	feedback_factor *= All.BlackHoleRadioModeFeedbackBoost;

  return BlackHoleOutflowEfficiency(f_Edd, bh_mass) + BlackHoleRadiativeEfficiency(f_Edd,
										   bh_mass) * feedback_factor;
  // std model: BlackHoleOutflowEfficiency = 0 and BlackHoleRadiativeEfficiency = eps_r  ->  returns eps_r * eps_f
}


static inline MyAtLeastDouble Mdot_Eddington(MyFloat bh_mass)
{
  MyAtLeastDouble const_coeff =
    (4 * M_PI * GRAVITY * C_LIGHT * PROTONMASS / (C_LIGHT * C_LIGHT * THOMPSON)) * bh_mass *
    All.UnitTime_in_s / All.HubbleParam;
  return const_coeff / BlackHoleMaximumRadiativeEfficiency(bh_mass);
}


static inline void Get_Bondi_Parameter(const int n, const MyFloat Rho, const MyFloat GasVel[3],
				       const MyFloat Entropy, const double AccretionFactor,
				       MyAtLeastDouble dt, MyAtLeastDouble *mdot)
{
  MyAtLeastDouble rho = Rho, bhvel;
  if(All.BlackHoleVelocityInBondiOn)
    bhvel =
      sqrt(pow(P[n].Vel[0] - GasVel[0], 2) + pow(P[n].Vel[1] - GasVel[1], 2) +
	   pow(P[n].Vel[2] - GasVel[2], 2));
  else
    bhvel = 0;

  if(All.ComovingIntegrationOn)
    {
      bhvel /= All.Time;
      rho /= pow(All.Time, 3);
    }
  MyAtLeastDouble soundspeed = sqrt(GAMMA * Entropy * pow(rho, GAMMA_MINUS1));

  MyAtLeastDouble accretionfactor = AccretionFactor;
  if(All.BlackHoleVariableAccretionFactorOn == 1)
    {
#ifdef SFR
      double rho_threshold;
#ifdef EB_SFR_MAGNETIC
      rho_threshold = SphP[n].DensityThreshold;
#else
      rho_threshold = All.PhysDensThresh;
#endif

      if(rho > rho_threshold)
	accretionfactor = pow(rho / rho_threshold, All.BlackHoleVariableAccretionSlope);
      else
	accretionfactor = 1.0;
#else
      accretionfactor = 1.0;
#endif
    }

  MyAtLeastDouble norm = pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);
  if(norm > 0)
    *mdot = 4. * M_PI * accretionfactor * All.G * All.G * pow(BPP(n).BH_Mass, 2) * (rho) / norm;
  else
    *mdot = 0;

  if(All.BlackHoleDetails >= 1)
    fprintf(FdBlackHolesDetails, " mdot=%g rho=%g csnd=%g bhvel=%g ", (double) *mdot, (double) rho,
	    (double) soundspeed, (double) bhvel);

#ifdef GM_MV_ANGMOM		/* Schaye+ MNRAS 446, 521–554 (2015) (eagle) */
  /* VphiMod as in Rosas-Guevara+ MNRAS 454, 1038–1057 (2015) */
  MyAtLeastDouble AngMomLimiter = (1 / GM_MV_ANGMOM) * pow(soundspeed / BPP(n).VphiMod, 3);
  if(AngMomLimiter < 1.0)
    {
      *mdot *= AngMomLimiter;
      if(All.BlackHoleDetails >= 1)
	fprintf(FdBlackHolesDetails, " AML=%g ", (double) AngMomLimiter);
    }
#endif

#if defined(LB_PRESSURE_DEPENDENT_ACCRETION) && defined(LT_NMet)
  if((All.BlackHoleLimitFeedbackOn & 2) == 2)
    {
      MyAtLeastDouble M_feedback = BPP(n).BH_SurroundingGasMass;
      MyAtLeastDouble Lambda_input = All.BlackHoleFeedbackFactor * BlackHoleRadiativeEfficiency(*mdot / Mdot_Eddington(BPP(n).BH_Mass), BPP(n).BH_Mass) * C_LIGHT * C_LIGHT * 4. * M_PI * All.BlackHoleAccretionFactor * GRAVITY * GRAVITY * BPP(n).BH_Mass * BPP(n).BH_Mass * PROTONMASS * PROTONMASS * All.UnitMass_in_g / (pow(soundspeed * All.UnitVelocity_in_cm_per_s, 3.0) * M_feedback);	//richtige Formel!
      MyAtLeastDouble Redshift, DZ;

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

      MyAtLeastDouble Temperature = BPP(n).BH_SurroundingTemperature / rho;
      MyAtLeastDouble Metallicity[LT_NMet];
      Metallicity[0] = BPP(n).BH_SurroundingMetallicity / rho;
      int q;
      for(q = 1; q < LT_NMet; q++)
	Metallicity[q] = 0.0;
      MyAtLeastDouble u_eq =
	GetUFromLambda(Lambda_input, rho, &Metallicity[0], Redshift, DZ, dt, &Temperature, 1e-12, 1e-6);
      MyAtLeastDouble PhysDensThresh = BPP(n).BH_SfrThreshold / rho;
      MyAtLeastDouble P_ref = (GAMMA_MINUS1) * PhysDensThresh * u_eq * All.UnitPressure_in_cgs;
      MyAtLeastDouble P_external = BPP(n).BH_SurroundingGasPressure / rho;
      MyAtLeastDouble factor = pow((P_external / P_ref), 2.0);
      if(P_external < P_ref)
	{
	  printf
	    ("WARNING: Plimiter: Lambda=%g, u_eq=%g, GAMMA_MINUS1=%g, PysDensThresh=%g, Metallicity=%g, Temperature=%g\n",
	     (double) Lambda_input, (double) u_eq, GAMMA_MINUS1, (double) PhysDensThresh, Metallicity[0],
	     (double) Temperature);
	  printf("P_ref=%g, P_ext=%g, factor=%g\n", (double) P_ref, (double) P_external, (double) factor);
	  *mdot *= pow((P_external / P_ref), 2.0);
	  if(All.BlackHoleDetails >= 1)
	    fprintf(FdBlackHolesDetails, " PextLim=%g ", (double) pow((P_external / P_ref), 2.0));
	}
    }
#endif

}


static inline void ApplyBlackHoleFriction(const int n, const MyAtLeastDouble dt)
{
  if(BPP(n).BH_Mass > 0 && dt > 0)
    if(P[n].Hsml < 0.5 * All.BlackHoleMergeDistFrac * All.ForceSoftening[5])
      {
	/* averaged value for colomb logarithm and integral over the distribution function */
	/* fac_friction = log(lambda) * [erf(x) - 2*x*exp(-x^2)/sqrt(pi)]                  */
	/*       lambda = b_max * v^2 / G / (M+m)                                          */
	/*        b_max = Size of system (e.g. Rvir)                                       */
	/*            v = Relative velocity of BH with respect to the environment          */
	/*            M = Mass of BH                                                       */
	/*            m = individual mass elements composing the large system (e.g. m<<M)  */
	/*            x = v/sqrt(2)/sigma                                                  */
	/*        sigma = width of the max. distr. of the host system                      */
	/*                (e.g. sigma = v_disp / 3                                         */
	MyAtLeastDouble relvel = 0, fac_friction = 10;
	for(int k = 0; k < 3; k++)
	  relvel += pow(P[n].Vel[k] - BPP(n).BH_SurroundingVel[k], 2);
	if(All.BlackHoleFrictionForceDynaicsOn == 1)
	  {
	    MyAtLeastDouble a_erf = 8 * (M_PI - 3) / (3 * M_PI * (4. - M_PI));
	    MyAtLeastDouble x = sqrt(relvel) / sqrt(2) / BPP(n).BH_sigma;
	    /* First term is aproximation of the error function */
	    fac_friction =
	      x / fabs(x) * sqrt(1 - exp(-x * x * (4 / M_PI + a_erf * x * x) / (1 + a_erf * x * x))) -
	      2 * x / sqrt(M_PI) * exp(-x * x);
	    MyAtLeastDouble lambda = BPP(n).BH_bmax * relvel / All.G / BPP(n).BH_Mass;
	    /*
	       printf("Task %d: x=%e, log(lambda)=%e, facerf=%e m=%e, sigma=%e\n",
	       ThisTask,x,log(lambda),fac_friction,P[n].BH_Mass,P[n].BH_sigma);
	     */
	    fac_friction *= log(lambda);
	    if(fac_friction < 0)
	      fac_friction = 0;
	  }
	else
	  fac_friction = 10;

	fac_friction *=
	  4 * M_PI * All.G * All.G * BPP(n).BH_SurroundingDensity * BPP(n).BH_Mass / relvel / sqrt(relvel);
	MyAtLeastDouble accgrv = 0, accfrc = 0;
	for(int k = 0; k < 3; k++)
	  {
	    accgrv += pow(P[n].GravAccel[k], 2);
	    accfrc += fac_friction * (pow(P[n].Vel[k] - BPP(n).BH_SurroundingVel[k], 2));
	  }
	accgrv = sqrt(accgrv);
	accfrc = sqrt(accfrc);
	if(accgrv > 0)
	  {
	    if(accfrc / accgrv > All.BlackHoleFrictionForceMaxCorrect)	/* Restrict friction force to be relatively mild */
	      {
		if(All.BlackHoleDetails >= 2)
		  fprintf(FdBlackHolesDetails,
			  "FRICTION: Large friction !! id = %llu fac = %e, vrel=%e, acc=(%e,%e,%e), adcc=(%e,%e,%e)\n",
			  (unsigned long long) P[n].ID, (double) fac_friction, (double) sqrt(relvel),
			  (double) P[n].GravAccel[0], (double) P[n].GravAccel[1], (double) P[n].GravAccel[2],
			  (float) fac_friction * (P[n].Vel[0] - BPP(n).BH_SurroundingVel[0]),
			  (float) fac_friction * (P[n].Vel[1] - BPP(n).BH_SurroundingVel[1]),
			  (float) fac_friction * (P[n].Vel[2] - BPP(n).BH_SurroundingVel[2]));
		fac_friction *= (accgrv / accfrc * All.BlackHoleFrictionForceMaxCorrect);
	      }
	    for(int k = 0; k < 3; k++)
	      P[n].GravAccel[k] -= fac_friction * (P[n].Vel[k] - BPP(n).BH_SurroundingVel[k]);
	  }
      }
}


static inline int ngb_treefind_blackhole_threads(MyLongDouble searchcenter[3], MyFloat hsml, int target,
						 int *startnode, int mode, int *exportflag,
						 int *exportnodecount, int *exportindex, int *ngblist)
{
  struct NODE *current;
#ifdef PERIODIC
  MyLongDouble xtmp;
#endif
  int numngb = 0;
  int no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  int p = no;
	  no = Nextnode[no];

	  if(All.BlackHoleRepositioningOn == 0 && All.BlackHoleFrictionForceOn == 0)
	    if(P[p].Type != 0 && P[p].Type != 5)
	      continue;

	  MyAtLeastDouble dist = hsml;
	  MyAtLeastDouble dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  MyAtLeastDouble dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  MyAtLeastDouble dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      int task;
	      if(exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		{
		  exportflag[task] = target;
		  exportnodecount[task] = NODELISTLENGTH;
		}
	      if(exportnodecount[task] == NODELISTLENGTH)
		{
		  int exitFlag = 0, nexp;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
		  {
		    if(Nexport >= All.BunchSize)
		      {
			BufferFullFlag = 1;
			exitFlag = 1;
		      }
		    else
		      {
			nexp = Nexport;
			Nexport++;
		      }
		  }

		  if(exitFlag)
		    return -1;

		  exportnodecount[task] = 0;
		  exportindex[task] = nexp;
		  DataIndexTable[nexp].Task = task;
		  DataIndexTable[nexp].Index = target;
		  DataIndexTable[nexp].IndexGet = nexp;
		}

	      DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
		DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

	      if(exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  MyAtLeastDouble dist = hsml + 0.5 * current->len;;
	  MyAtLeastDouble dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  MyAtLeastDouble dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  MyAtLeastDouble dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}

/*static inline int ngb_treefind_blackhole(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			   int *nexport, int *nsend_local)
{
  struct NODE *current;
#ifdef PERIODIC
  MyLongDouble xtmp;
#endif
  int nexport_save = *nexport;
  int numngb = 0;
  int no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	// single particle
	{
	  int p = no;
	  no = Nextnode[no];

	  if(All.BlackHoleRepositioningOn == 0)
	    if(P[p].Type != 0 && P[p].Type != 5)
	      continue;

	  MyAtLeastDouble dist = hsml;
	  MyAtLeastDouble dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  MyAtLeastDouble dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  MyAtLeastDouble dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	// pseudo particle
	    {
	      if(mode == 1)
		endrun(12312);

	      int task;
	      if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= All.BunchSize)
		    {
		      *nexport = nexport_save;
		      for(task = 0; task < NTask; task++)
			nsend_local[task] = 0;
		      for(no = 0; no < nexport_save; no++)
			nsend_local[DataIndexTable[no].Task]++;
		      return -1;
		    }
		  Exportnodecount[task] = 0;
		  Exportindex[task] = *nexport;
		  DataIndexTable[*nexport].Task = task;
		  DataIndexTable[*nexport].Index = target;
		  DataIndexTable[*nexport].IndexGet = *nexport;
		  *nexport = *nexport + 1;
		  nsend_local[task]++;
		}

	      DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	// we reached a top-level node again, which means that we are done with the branch
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  no = current->u.d.sibling;	// in case the node can be discarded

	  MyAtLeastDouble dist = hsml + 0.5 * current->len;;
	  MyAtLeastDouble dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  MyAtLeastDouble dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  MyAtLeastDouble dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  // now test against the minimal sphere enclosing everything
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	// ok, we need to open the node
	}
    }

  *startnode = -1;
  return numngb;
}*/



#endif // BLACK_HOLES

#endif // BLACK_HOLES_INCLUDE
