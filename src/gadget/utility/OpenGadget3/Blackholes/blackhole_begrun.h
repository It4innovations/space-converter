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

#include "blackhole.h"

void check_and_report_settings_for_blackhole()
{
#ifndef BLACK_HOLES
  PANIC_IF(All.BlackHolesOn != 0,
	   "Blackholes switched on in param file but BLACK_HOLES was not set at compile time");
#else
  // Check if On/Off settings are correct
  if(All.BlackHolesOn == 0)
    {
      VERBOSE(0, "BH: Blackhole model is available, but switched off\n");
    }
  else
    {
      if(All.BlackHolesOn == 1)
	{
	  printf("BH: Blackhole model is activated (refer to Springel & DiMattheo 2005)\n");
	}
      else
	{
	  PANIC("BH: Blackhole model is available, therefore BlackHolesOn has to be 0/1 !");
	}
      // Check if other needed modules are acailable 
#ifndef FOF
      VERBOSE(0,
	      "BH: WARNING: FOF was not set at compile time, therefore no new black holes will be seeded !\n");
#endif
#ifndef SUBFIND
      PANIC_IF(All.BlackHoleFrictionForceDynaicsOn == 1,
	       "BH: FrictionForceDynaics is switched on but needs SUBFIND to be set at compiletime !");
#endif
#ifndef EVALPOTENTIAL
      PANIC_IF(All.BlackHoleRepositioningOn == 1,
	       "BH: Repositioning is switched on but needs EVALPOTENTIAL to be set at compiletime !");
      PANIC_IF(All.BlackHoleSeedStarMassFraction > 0,
	       "BH: Seeding on Stars implies to seed on minimum binding energy and  needs EVALPOTENTIAL to be set at compiletime !");
#endif


      // Check validity of On/Off switches
      PANIC_IF(All.BlackHoleSwallowGasOn < 0
	       || All.BlackHoleSwallowGasOn > 2, "BH: BlackHoleSwallowGasOn has to be 0/1/2 !");
      PANIC_IF(All.BlackHoleDetails < 0
	       || All.BlackHoleDetails > 3, "BH: BlackHoleDetails has to be 0,1,2 or 3 !");
      PANIC_IF(All.BlackHoleRepositioningOn < 0
	       || All.BlackHoleRepositioningOn > 1, "BH: BlackHoleRepositioningOn has to be 0/1 !");
      PANIC_IF(All.BlackHoleHotColdAccretionOn < 0
	       || All.BlackHoleHotColdAccretionOn > 1, "BH: BlackHoleHotColdAccretionOn has to be 0/1 !");
      PANIC_IF(All.BlackHoleVelocityInBondiOn < 0
	       || All.BlackHoleVelocityInBondiOn > 1, "BH: BlackHoleVelocityInBondiOn has to be 0/1 !");
      PANIC_IF(All.BlackHoleVariableAccretionFactorOn < 0
	       || All.BlackHoleVariableAccretionFactorOn > 1,
	       "BH: BlackHoleVariableAccretionFactorOn has to be 0/1 !");
      PANIC_IF(All.BlackHoleVariableEfficiencyOn < 0
	       || All.BlackHoleVariableEfficiencyOn > 2,
	       "BH: BlackHoleVariableEfficiencyOn has to be 0,1 or 2 !");
      PANIC_IF(All.BlackHoleOutflowModelOn < 0
	       || All.BlackHoleOutflowModelOn > 1, "BH: BlackHoleOutflowModelOn has to be 0/1 !");
      PANIC_IF(All.BlackHoleLimitFeedbackOn < 0
	       || All.BlackHoleLimitFeedbackOn > 3, "BH: BlackHoleLimitFeedbackOn has to be 0/1/2/3 !");
      PANIC_IF(All.BlackHoleFrictionForceOn < 0
	       || All.BlackHoleFrictionForceOn > 1, "BH: BlackHoleFrictionForceOn has to be 0/1 !");
      PANIC_IF(All.BlackHoleFrictionForceDynaicsOn < 0
	       || All.BlackHoleFrictionForceDynaicsOn > 1,
	       "BH: BlackHoleFrictionForceDynaicsOn has to be 0/1 !");
#ifdef LB_PRESSURE_DEPENDENT_ACCRETION
      PANIC_IF((All.BlackHoleLimitFeedbackOn & 2) != 2,
	       "BH: LB_PRESSURE_DEPENDENT_ACCRETION is defined, but BlackHoleLimitFeedbackOn is 0 or 1 !");
#else
      PANIC_IF((All.BlackHoleLimitFeedbackOn & 2) == 2,
	       "BH: BlackHoleLimitFeedbackOn implies pressure dependent limiter for accretion, but LB_PRESSURE_DEPENDENT_ACCRETION is not defined!");
#endif

      // Check combinations:
      if(All.BlackHoleOutflowModelOn == 1 && All.BlackHoleRadioModeFeedbackBoost > 1)
	if(All.LevelOfStrickness == 0)
	  PANIC("BH: RadioModeFeedbackBoost and OutflowModel can not be used together !");
	else
	  VERBOSE(0, "BH: WARNING: RadioModeFeedbackBoost will be ignored as OutflowModel is switched on !");

      if(All.BlackHoleOutflowModelOn == 1 && All.BlackHoleRadioModeTreshold != 0.05)
	if(All.LevelOfStrickness == 0)
	  PANIC
	    ("BH: BlackHoleRadioModeTreshold is different from 0,05 used in OutflowModel. Will be inored but can lead to confusion !");
	else
	  VERBOSE(0,
		  "BH: WARNING: BlackHoleRadioModeTreshold will be ignored as OutflowModel is using 0.05 in any case !");

      PANIC_IF(All.BlackHoleFrictionForceDynaicsOn == 1
	       && All.BlackHoleFrictionForceOn == 0,
	       "   FrictionForceDynaics needs FrictionForce to be turned on !");

      // Give credits: 
      VERBOSE_CONDITION(All.BlackHoleFrictionForceOn == 1, 0,
			"BH: BlackHoleFrictionForce is on (refer to Hirschmann+ 2014)");
      VERBOSE_CONDITION(All.BlackHoleRepositioningOn == 0, 0,
			"BH: BlackHoleRepositioning is off (refer to Hirschmann+ 2014)");
      VERBOSE_CONDITION(All.BlackHoleVariableEfficiencyOn == 1, 0,
			"BH: BlackHoleVariableEfficiency is on (refer to Steinborn+ 2015)");
      VERBOSE_CONDITION(All.BlackHoleVariableEfficiencyOn == 2, 0,
			"BH: BlackHoleVariableEfficiency is on (refer to Sala+ 2022)");
      VERBOSE_CONDITION(All.BlackHoleSeedStarMassFraction > 0, 0,
			"BH: Blackholes are seeded onto stellar body (refer to Hirschmann+ 2014)");
      VERBOSE_CONDITION(All.BlackHoleHotColdAccretionOn == 1, 0,
			"BH: BlackHoleHotColdAccretion is on (refer to Steinborn+ 2015)");
      VERBOSE_CONDITION(All.BlackHoleOutflowModelOn == 1, 0,
			"BH: BlackHoleOutflowModel is on (refer to Steinborn+ 2015)");
      VERBOSE_CONDITION((All.BlackHoleLimitFeedbackOn & 2) == 2, 0,
			"BH: BlackHoleLimitFeedback implies pressure dependent accretion limiter (refer to Vogelsberger+ 2014)");
      VERBOSE_CONDITION(All.BlackHoleRadioModeFeedbackBoost > 1, 0,
			"BH: BlackHoleRadioModeFeedbackBoost is on (refer to Fabjan+ 2010)");
      VERBOSE_CONDITION(All.BlackHoleVariableAccretionFactorOn == 1, 0,
			"BH: BlackHoleVariableAccretionFactor is on (refer to Booth&Schaye 2013)");

      // Warning for expert modules:
      if(All.BlackHoleFrictionForceOn == 1)
	if(All.LevelOfStrickness == 0)
	  PANIC("BH: FrictionForce is an expert module and only avalable for LevelOfStrickness >= 1 !");
	else
	  VERBOSE(0,
		  "BH: WARNING: FrictionForce is an expert module, check results carefully and use on own risc !\n");

      if(All.BlackHoleRepositioningOn == 0)
	if(All.LevelOfStrickness == 0)
	  PANIC
	    ("BH: Not using Repositioning is an expert module and only avalable for LevelOfStrickness >= 1 !");
	else
	  VERBOSE(0,
		  "BH: WARNING: Not using Repositioning is an expert module, check results carefully and use on own risc !\n");

      if(All.BlackHoleVariableAccretionFactorOn == 1)
	if(All.LevelOfStrickness == 0)
	  PANIC
	    ("BH: Not using VariableAccretionFactor is an expert module and only avalable for LevelOfStrickness >= 1 !");
	else
	  VERBOSE(0,
		  "BH: WARNING: using VariableAccretionFactor is an expert module, check results carefully and use on own risc !\n");

      if(All.BlackHoleThermalFeedbackOn == 0 && All.BlackHoleKineticFeedbackOn == 0)
	if(All.LevelOfStrickness == 0)
	  PANIC
	    ("BH: You are running with black holes but without any feedback, please switch to LevelOfStrickness >= 1 !");
	else
	  VERBOSE(0,
		  "BH: WARNING: You are running with black holes but without any feedback, use on own risc !\n");



      // General warnings for expert modules:
      VERBOSE_CONDITION(All.BlackHoleDetails >= 1, 0,
			"BH: WARNING: Details are written by every CPU which can lead to file system problems for very large number of MPI ranks !\n");
      VERBOSE_CONDITION(All.BlackHoleVelocityInBondiOn == 0, 0,
			"BH: WARNING: Not using VelocityInBondi is not recomendet !\n");
      VERBOSE_CONDITION(All.BlackHoleVelocityInBondiOn == 0, 0,
			"BH: WARNING: Not using VelocityInBondi is not recomendet !\n");
      VERBOSE_CONDITION(All.BlackHoleIgnoreMomentum == 0 && All.BlackHoleRepositioningOn == 0, 0,
			"BH: WARNING: accounting for Momentum and not Repositioning is not recomendet !");
    }
#endif
}


void change_param_values_blackhole(struct global_data_all_processes *All,
				   struct global_data_all_processes *all, FILE **logfile)
{

  if(ThisTask == 0)
    printf("BH: Checking paramtere changes ...\n");

  // Very dangerous changes
  if(All->BlackHolesOn != all->BlackHolesOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHolesOn changed from %d to %d \n",
		All->Time, All->BlackHolesOn, all->BlackHolesOn);
      PANIC_IF(All->LevelOfStrickness <= 2,
	       " changing BlackHolesOn is not at all recomended, needs LevelOfStrickness >= 3\n");
      All->BlackHolesOn = all->BlackHolesOn;
    }

#ifdef BLACK_HOLES

  // Changing parfameters of models
  if(All->BlackHoleMaxAccretionRadius != all->BlackHoleMaxAccretionRadius)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleMaxAccretionRadius changed from %f to %f \n",
		All->Time, All->BlackHoleMaxAccretionRadius, all->BlackHoleMaxAccretionRadius);
      All->BlackHoleMaxAccretionRadius = all->BlackHoleMaxAccretionRadius;
    }

  if(All->BlackHoleRadiativeEfficiency != all->BlackHoleRadiativeEfficiency)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadiativeEfficiency changed from %f to %f \n",
		All->Time, All->BlackHoleRadiativeEfficiency, all->BlackHoleRadiativeEfficiency);
      All->BlackHoleRadiativeEfficiency = all->BlackHoleRadiativeEfficiency;
    }

  if(All->BlackHoleRadiativeEfficiencySlope != all->BlackHoleRadiativeEfficiencySlope)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadiativeEfficiencySlope changed from %f to %f \n",
		All->Time, All->BlackHoleRadiativeEfficiencySlope, all->BlackHoleRadiativeEfficiencySlope);
      All->BlackHoleRadiativeEfficiencySlope = all->BlackHoleRadiativeEfficiencySlope;
    }

  if(All->BlackHoleRadiativeEfficiencyNorm != all->BlackHoleRadiativeEfficiencyNorm)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadiativeEfficiencyNorm changed from %f to %f \n",
		All->Time, All->BlackHoleRadiativeEfficiencyNorm, all->BlackHoleRadiativeEfficiencyNorm);
      All->BlackHoleRadiativeEfficiencyNorm = all->BlackHoleRadiativeEfficiencyNorm;
    }

  if(All->BlackHoleRadiativeEfficiencySlopeMdot != all->BlackHoleRadiativeEfficiencySlopeMdot)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadiativeEfficiencyNorm changed from %f to %f \n",
		All->Time, All->BlackHoleRadiativeEfficiencySlopeMdot,
		all->BlackHoleRadiativeEfficiencySlopeMdot);
      All->BlackHoleRadiativeEfficiencyNorm = all->BlackHoleRadiativeEfficiencyNorm;
    }

  if(All->BlackHoleRadiativeEfficiencyNorm != all->BlackHoleRadiativeEfficiencyNorm)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadiativeEfficiencyNorm changed from %f to %f \n",
		All->Time, All->BlackHoleRadiativeEfficiencyNorm, all->BlackHoleRadiativeEfficiencyNorm);
      All->BlackHoleRadiativeEfficiencyNorm = all->BlackHoleRadiativeEfficiencyNorm;
    }

  if(All->BlackHoleRadiativeEfficiencyMax != all->BlackHoleRadiativeEfficiencyMax)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadiativeEfficiencyMax changed from %f to %f \n",
		All->Time, All->BlackHoleRadiativeEfficiencyMax, all->BlackHoleRadiativeEfficiencyMax);
      All->BlackHoleRadiativeEfficiencyMax = all->BlackHoleRadiativeEfficiencyMax;
    }

  if(All->BlackHoleRadiativeEfficiencyMin != all->BlackHoleRadiativeEfficiencyMin)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadiativeEfficiencyMin changed from %f to %f \n",
		All->Time, All->BlackHoleRadiativeEfficiencyMin, all->BlackHoleRadiativeEfficiencyMin);
      All->BlackHoleRadiativeEfficiencyMin = all->BlackHoleRadiativeEfficiencyMin;
    }

  if(All->BlackHoleMergeCsndFrac != all->BlackHoleMergeCsndFrac)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleMergeCsndFrac changed from %f to %f \n",
		All->Time, All->BlackHoleMergeCsndFrac, all->BlackHoleMergeCsndFrac);
      All->BlackHoleMergeCsndFrac = all->BlackHoleMergeCsndFrac;
    }

  if(All->BlackHoleMergeDistFrac != all->BlackHoleMergeDistFrac)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleMergeDistFrac changed from %f to %f \n",
		All->Time, All->BlackHoleMergeDistFrac, all->BlackHoleMergeDistFrac);
      All->BlackHoleMergeDistFrac = all->BlackHoleMergeDistFrac;
    }

  if(All->BlackHoleMergeBindingFrac != all->BlackHoleMergeBindingFrac)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleMergeBindingFrac changed from %f to %f \n",
		All->Time, All->BlackHoleMergeBindingFrac, all->BlackHoleMergeBindingFrac);
      All->BlackHoleMergeBindingFrac = all->BlackHoleMergeBindingFrac;
    }

  // Changing between models, potentially dangerous
  if(All->BlackHoleFrictionForceOn != all->BlackHoleFrictionForceOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleFrictionForceOn changed from %d to %d \n",
		All->Time, All->BlackHoleFrictionForceOn, all->BlackHoleFrictionForceOn);
      PANIC_IF(All->LevelOfStrickness <= 1,
	       " changing BlackHoleFrictionForceOn is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleFrictionForceOn = all->BlackHoleFrictionForceOn;
    }

  if(All->BlackHoleFrictionForceOn != all->BlackHoleFrictionForceOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleFrictionForceOn changed from %d to %d \n",
		All->Time, All->BlackHoleFrictionForceOn, all->BlackHoleFrictionForceOn);
      PANIC_IF(All->LevelOfStrickness <= 1,
	       " changing BlackHoleFrictionForceOn is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleFrictionForceOn = all->BlackHoleFrictionForceOn;
    }

  if(All->BlackHoleSeedStarMassFraction != all->BlackHoleSeedStarMassFraction)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleSeedStarMassFractio changed from %f to %f \n",
		All->Time, All->BlackHoleSeedStarMassFraction, all->BlackHoleSeedStarMassFraction);
      if(All->BlackHoleSeedStarMassFraction == 0 || all->BlackHoleSeedStarMassFraction == 0)
	PANIC_IF(All->LevelOfStrickness <= 1,
		 " changing seedin mechanism is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleSeedStarMassFraction = all->BlackHoleSeedStarMassFraction;
    }

  if(All->BlackHoleSeedDMFraction != all->BlackHoleSeedDMFraction)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleSeedDMFraction changed from %f to %f \n",
		All->Time, All->BlackHoleSeedDMFraction, all->BlackHoleSeedDMFraction);
      if(All->BlackHoleSeedStarMassFraction > 0)
	PANIC_IF(All->LevelOfStrickness <= 1,
		 " changing parameters of seeding mechanism is not recomended, needs LevelOfStrickness >= 2\n");
      if(All->BlackHoleSeedStarMassFraction == 0)
	PANIC_IF(All->LevelOfStrickness <= 1,
		 " changing parameters of seeding mechanism while switched off is meaningless, needs LevelOfStrickness >= 2\n");
      All->BlackHoleSeedDMFraction = all->BlackHoleSeedDMFraction;
    }

  if(All->BlackHoleSeedGasFraction != all->BlackHoleSeedGasFraction)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleSeedGasFraction changed from %f to %f \n",
		All->Time, All->BlackHoleSeedGasFraction, all->BlackHoleSeedGasFraction);
      if(All->BlackHoleSeedStarMassFraction > 0)
	PANIC_IF(All->LevelOfStrickness <= 1,
		 " changing parameters of seeding mechanism is not recomended, needs LevelOfStrickness >= 2\n");
      if(All->BlackHoleSeedStarMassFraction == 0)
	PANIC_IF(All->LevelOfStrickness <= 1,
		 " changing parameters of seeding mechanism while switched off is meaningless, needs LevelOfStrickness >= 2\n");
      All->BlackHoleSeedGasFraction = all->BlackHoleSeedGasFraction;
    }

  if(All->BlackHoleRepositioningOn != all->BlackHoleRepositioningOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleFrictionForceOn changed from %d to %d \n",
		All->Time, All->BlackHoleRepositioningOn, all->BlackHoleRepositioningOn);
      PANIC_IF(All->LevelOfStrickness <= 1,
	       " changing BlackHoleRepositioningOn is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleRepositioningOn = all->BlackHoleRepositioningOn;
    }

  if(All->BlackHoleHotColdAccretionOn != all->BlackHoleHotColdAccretionOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleHotColdAccretionOn changed from %d to %d \n",
		All->Time, All->BlackHoleHotColdAccretionOn, all->BlackHoleHotColdAccretionOn);
      PANIC_IF(All->LevelOfStrickness <= 1,
	       " changing BlackHoleHotColdAccretionOn is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleHotColdAccretionOn = all->BlackHoleHotColdAccretionOn;
    }

  if(All->BlackHoleOutflowModelOn != all->BlackHoleOutflowModelOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleOutflowModelOn changed from %d to %d \n",
		All->Time, All->BlackHoleOutflowModelOn, all->BlackHoleOutflowModelOn);
      PANIC_IF(All->LevelOfStrickness <= 1,
	       " changing BlackHoleOutflowModelOn is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleOutflowModelOn = all->BlackHoleOutflowModelOn;
    }

  if(All->BlackHoleLimitFeedbackOn != all->BlackHoleLimitFeedbackOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleLimitFeedbackOn changed from %d to %d \n",
		All->Time, All->BlackHoleLimitFeedbackOn, all->BlackHoleLimitFeedbackOn);
      PANIC_IF(All->LevelOfStrickness <= 0,
	       " changing BlackHoleLimitFeedbackOn can is not recomended, needs LevelOfStrickness >= 1\n");
      All->BlackHoleLimitFeedbackOn = all->BlackHoleLimitFeedbackOn;
    }

  if(All->BlackHoleRadioModeFeedbackBoost != all->BlackHoleRadioModeFeedbackBoost)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleRadioModeFeedbackBoost changed from %f to %f \n",
		All->Time, All->BlackHoleRadioModeFeedbackBoost, all->BlackHoleRadioModeFeedbackBoost);
      PANIC_IF(All->LevelOfStrickness <= 1,
	       " changing BlackHoleRadioModeFeedbackBoost is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleRadioModeFeedbackBoost = all->BlackHoleRadioModeFeedbackBoost;
    }

  if(All->BlackHoleVariableAccretionFactorOn != all->BlackHoleVariableAccretionFactorOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleVariableAccretionFactorOn changed from %d to %d \n",
		All->Time, All->BlackHoleVariableAccretionFactorOn, all->BlackHoleVariableAccretionFactorOn);
      PANIC_IF(All->LevelOfStrickness <= 1,
	       " changing BlackHoleVariableAccretionFactorOn is not recomended, needs LevelOfStrickness >= 2\n");
      All->BlackHoleVariableAccretionFactorOn = all->BlackHoleVariableAccretionFactorOn;
    }

  if(All->BlackHoleThermalFeedbackOn != all->BlackHoleThermalFeedbackOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleThermalFeedbackOn changed from %d to %d \n",
		All->Time, All->BlackHoleThermalFeedbackOn, all->BlackHoleThermalFeedbackOn);
      PANIC_IF(All->LevelOfStrickness <= 0,
	       " changing BlackHoleThermalFeedback is not recomended, needs LevelOfStrickness >= 1\n");
      All->BlackHoleThermalFeedbackOn = all->BlackHoleThermalFeedbackOn;
    }

  if(All->BlackHoleKineticFeedbackOn != all->BlackHoleKineticFeedbackOn)
    {
      if(ThisTask == 0)
	fprintf(*logfile, "Time %f: BlackHoleKineticFeedbackOn changed from %d to %d \n",
		All->Time, All->BlackHoleKineticFeedbackOn, all->BlackHoleKineticFeedbackOn);
      PANIC_IF(All->LevelOfStrickness <= 0,
	       " changing BlackHoleKineticFeedback is not recomended, needs LevelOfStrickness >= 1\n");
      All->BlackHoleKineticFeedbackOn = all->BlackHoleKineticFeedbackOn;
    }

#endif
  if(ThisTask == 0)
    fflush(*logfile);
}


void report_internal_model_for_blackhole()
{
#ifdef BLACK_HOLES
  FILE *fmodel;

  int N_BIN_BH_MASS = 20;
  MyAtLeastDouble MIN_BH_MASS = 5.0;
  MyAtLeastDouble MAX_BH_MASS = 9.0;
  int N_BIN_FEDD = 20;
  MyAtLeastDouble MIN_FEDD = -8.0;
  MyAtLeastDouble MAX_FEDD = 1.3;
  if(All.BlackHoleVariableEfficiencyOn == 2)
    {
      MIN_BH_MASS = 5.0;
      MAX_BH_MASS = 10.0;
      MIN_FEDD = -3.0;		//In this model the variable is reinterpreted as mdot, in Log(Msun/yr)
      MAX_FEDD = 2;
    }
  MyFloat frac;
  MyFloat mbh;

  if(!(fmodel = fopen("BH_MODEL_Mbh_Medd_MaxRadEff.txt", "w")))
    PANIC("Could not open <BH_MODEL_Medd_Mbh.txt> for writing !");
  fprintf(fmodel, "# Mbh/1e10Msun Medd MaxRadEff\n");
  for(int i = 0; i < N_BIN_BH_MASS; i++)
    {
      frac = (double) i / (N_BIN_BH_MASS - 1);
      mbh = pow(10.0, MIN_BH_MASS + frac * (MAX_BH_MASS - MIN_BH_MASS)) / 1e10;
      MyAtLeastDouble medd = Mdot_Eddington(mbh);
      MyAtLeastDouble MaxRadEff = BlackHoleMaximumRadiativeEfficiency(mbh);
      fprintf(fmodel, "%1.2e %1.2e %1.2e\n", mbh, medd, MaxRadEff);
    }
  fclose(fmodel);

  if(!(fmodel = fopen("BH_MODEL_RadEff.txt", "w")))
    PANIC("Could not open <BH_MODEL_RadEff.txt> for writing !");
  (All.BlackHoleVariableEfficiencyOn == 2) ? fprintf(fmodel,
						     "# Mbh(1e10Msun) mdot(Msun/yr) -> \n") : fprintf(fmodel,
												      "# Mbh(1e10Msun) f_edd -> \n");

  fprintf(fmodel, "         ");
  for(int j = 0; j < N_BIN_FEDD; j++)
    {
      frac = (double) j / (N_BIN_FEDD - 1);
      MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
      fprintf(fmodel, "%1.2e ", fedd);
    }
  fprintf(fmodel, "\n");
  for(int i = 0; i < N_BIN_BH_MASS; i++)
    {
      frac = (double) i / (N_BIN_BH_MASS - 1);
      mbh = pow(10.0, MIN_BH_MASS + frac * (MAX_BH_MASS - MIN_BH_MASS)) / 1e10;
      fprintf(fmodel, "%1.2e ", mbh);
      for(int j = 0; j < N_BIN_FEDD; j++)
	{
	  frac = (double) j / (N_BIN_FEDD - 1);
	  MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
	  if(All.BlackHoleVariableEfficiencyOn == 2)
	    fedd /= (All.UnitMass_in_g / All.UnitTime_in_s / SOLAR_MASS * SEC_PER_YEAR);	//convert back Msun/yr in internal units to pass to the function.

	  MyAtLeastDouble RadEff = BlackHoleRadiativeEfficiency(fedd, mbh);
	  fprintf(fmodel, "%1.2e ", RadEff);
	}
      fprintf(fmodel, "\n");
    }
  fclose(fmodel);

  if(!(fmodel = fopen("BH_MODEL_OutflowEff.txt", "w")))
    PANIC("Could not open <BH_MODEL_OutflowEff.txt> for writing !");
  (All.BlackHoleVariableEfficiencyOn == 2) ? fprintf(fmodel,
						     "# Mbh(1e10Msun) mdot(Msun/yr) -> \n") : fprintf(fmodel,
												      "# Mbh(1e10Msun) f_edd -> \n");
  fprintf(fmodel, "         ");
  for(int j = 0; j < N_BIN_FEDD; j++)
    {
      frac = (double) j / (N_BIN_FEDD - 1);
      MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
      fprintf(fmodel, "%1.2e ", fedd);
    }
  fprintf(fmodel, "\n");
  for(int i = 0; i < N_BIN_BH_MASS; i++)
    {
      frac = (double) i / (N_BIN_BH_MASS - 1);
      mbh = pow(10.0, MIN_BH_MASS + frac * (MAX_BH_MASS - MIN_BH_MASS)) / 1e10;
      fprintf(fmodel, "%1.2e ", mbh);
      for(int j = 0; j < N_BIN_FEDD; j++)
	{
	  frac = (double) j / (N_BIN_FEDD - 1);
	  MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
	  if(All.BlackHoleVariableEfficiencyOn == 2)
	    fedd /= (All.UnitMass_in_g / All.UnitTime_in_s / SOLAR_MASS * SEC_PER_YEAR);	//convert back Msun/yr in internal units to pass to the function.
	  MyAtLeastDouble RadEff = BlackHoleOutflowEfficiency(fedd, mbh);
	  fprintf(fmodel, "%1.2e ", RadEff);
	}
      fprintf(fmodel, "\n");
    }
  fclose(fmodel);

  if(!(fmodel = fopen("BH_MODEL_TotEff.txt", "w")))
    PANIC("Could not open <BH_MODEL_TotEff.txt> for writing !");
  (All.BlackHoleVariableEfficiencyOn == 2) ? fprintf(fmodel,
						     "# Mbh(1e10Msun) mdot(Msun/yr) -> \n") : fprintf(fmodel,
												      "# Mbh(1e10Msun) f_edd -> \n");
  fprintf(fmodel, "         ");
  for(int j = 0; j < N_BIN_FEDD; j++)
    {
      frac = (double) j / (N_BIN_FEDD - 1);
      MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
      fprintf(fmodel, "%1.2e ", fedd);
    }
  fprintf(fmodel, "\n");
  for(int i = 0; i < N_BIN_BH_MASS; i++)
    {
      frac = (double) i / (N_BIN_BH_MASS - 1);
      mbh = pow(10.0, MIN_BH_MASS + frac * (MAX_BH_MASS - MIN_BH_MASS)) / 1e10;
      fprintf(fmodel, "%1.2e ", mbh);
      for(int j = 0; j < N_BIN_FEDD; j++)
	{
	  frac = (double) j / (N_BIN_FEDD - 1);
	  MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
	  if(All.BlackHoleVariableEfficiencyOn == 2)
	    fedd /= (All.UnitMass_in_g / All.UnitTime_in_s / SOLAR_MASS * SEC_PER_YEAR);	//convert back Msun/yr in internal units to pass to the function.
	  MyAtLeastDouble RadEff = BlackHoleTotalEfficiency(fedd, mbh);
	  fprintf(fmodel, "%1.2e ", RadEff);
	}
      fprintf(fmodel, "\n");
    }
  fclose(fmodel);

  if(!(fmodel = fopen("BH_MODEL_TotFeedEff.txt", "w")))
    PANIC("Could not open <BH_MODEL_TotFeedEff.txt> for writing !");
  (All.BlackHoleVariableEfficiencyOn == 2) ? fprintf(fmodel, "# Mbh mdot -> \n") : fprintf(fmodel,
											   "# Mbh f_edd -> \n");
  fprintf(fmodel, "         ");
  for(int j = 0; j < N_BIN_FEDD; j++)
    {
      frac = (double) j / (N_BIN_FEDD - 1);
      MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
      fprintf(fmodel, "%1.2e ", fedd);
    }
  fprintf(fmodel, "\n");
  for(int i = 0; i < N_BIN_BH_MASS; i++)
    {
      frac = (double) i / (N_BIN_BH_MASS - 1);
      mbh = pow(10.0, MIN_BH_MASS + frac * (MAX_BH_MASS - MIN_BH_MASS)) / 1e10;
      fprintf(fmodel, "%1.2e ", mbh);
      for(int j = 0; j < N_BIN_FEDD; j++)
	{
	  frac = (double) j / (N_BIN_FEDD - 1);
	  MyFloat fedd = pow(10.0, MIN_FEDD + frac * (MAX_FEDD - MIN_FEDD));
	  if(All.BlackHoleVariableEfficiencyOn == 2)
	    fedd /= (All.UnitMass_in_g / All.UnitTime_in_s / SOLAR_MASS * SEC_PER_YEAR);	//convert back Msun/yr in internal units to pass to the function.

	  MyAtLeastDouble RadEff = BlackHoleTotalFeedbackEfficiency(fedd, mbh);
	  fprintf(fmodel, "%1.2e ", RadEff);
	}
      fprintf(fmodel, "\n");
    }
  fclose(fmodel);



#endif
}
