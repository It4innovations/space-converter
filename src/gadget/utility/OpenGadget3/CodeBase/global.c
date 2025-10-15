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

#include "allvars.h"
#include "proto.h"


/* This routine computes various global properties of the particle
 * distribution and stores the result in the struct `SysState'.
 * Currently, not all the information that's computed here is 
 * actually used (e.g. momentum is not really used anywhere),
 * just the energies are written to a log-file every once in a while.
 */
void compute_global_quantities_of_system(void)
{
  int dt_step;
  struct state_of_system sys;
  double a1, a2, a3;
  double entr = 0, egyspec, vel[3];
  double dt_entr, dt_gravkick, dt_hydrokick;

#if defined(MAGNETIC_STATISTICS)
  double dt_mag, mag[3];
  mag[0] = mag[1] = mag[2] = 0.0;
  sys.EnergyMag = 0.0;
#if defined(TRACEDIVB)
  sys.DivBerr = 0;
#endif
#endif
#if defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB)
  double volweight = 0, massweight = 0, Volweight = 0, Massweight = 0;
#endif

  if(All.ComovingIntegrationOn)
    {
      a1 = All.Time;
      a2 = All.Time * All.Time;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      a1 = a2 = a3 = 1;
    }


  for(int n = 0; n < 6; n++)
    {
      sys.MassComp[n] = sys.EnergyKinComp[n] = sys.EnergyPotComp[n] = sys.EnergyIntComp[n] = 0;

      for(int j = 0; j < 4; j++)
	sys.CenterOfMassComp[n][j] = sys.MomentumComp[n][j] = sys.AngMomentumComp[n][j] = 0;
    }

  for(int i = 0; i < NumPart; i++)
    {
      sys.MassComp[P[i].Type] += P[i].Mass;

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
      sys.EnergyPotComp[P[i].Type] += 0.5 * P[i].Mass * P[i].Potential / a1;
#endif

      if(WAKEUP > 0)
	{
	  dt_step = P[i].dt_step;
	}
      else
	{
	  dt_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
	}

      if(All.ComovingIntegrationOn)
	{
	  dt_entr = (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;
	  dt_gravkick = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
	    get_gravkick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
	  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
	    get_hydrokick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
#if defined(MAGNETIC_STATISTICS)
	  dt_mag = get_magkick_factor(P[i].Ti_begstep, All.Ti_Current) -
	    get_magkick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
#endif
	}
      else
	dt_entr = dt_gravkick = dt_hydrokick =
#if defined(MAGNETIC_STATISTICS)
	  dt_mag =
#endif
	  (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

      for(int j = 0; j < 3; j++)
	{
	  vel[j] = P[i].Vel[j] + P[i].GravAccel[j] * dt_gravkick;
	  if(P[i].Type == 0)
	    vel[j] += SphP[i].HydroAccel[j] * dt_hydrokick;

	}
      if(P[i].Type == 0)
	{
	  entr = SphP[i].Entropy + SphP[i].DtEntropy * dt_entr;
#if defined(MAGNETIC_STATISTICS)
	  for(int j = 0; j < 3; j++)
	    mag[j] = SphP[i].b2.BPred[j] + SphP[i].DtB[j] * dt_mag;
#endif
	}
#ifdef PMGRID
      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
	  get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkick = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;

      for(int j = 0; j < 3; j++)
	{
	  vel[j] += P[i].GravPM[j] * dt_gravkick;
	}
#endif

      sys.EnergyKinComp[P[i].Type] +=
	0.5 * P[i].Mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / a2;

      if(P[i].Type == 0)
	{
#if defined(ISOTHERM_EQS)
	  egyspec = entr;
#else
	  egyspec = entr / (GAMMA_MINUS1) * pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif

#if defined(MAGNETIC_STATISTICS)
	  sys.EnergyMag +=
	    P[i].Mass / SphP[i].Density * MU0_1 / 2.0 * (mag[0] * mag[0] + mag[1] * mag[1] + mag[2] * mag[2]);
#endif
#if defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB)
	  if(P[i].ID != 0
	     && (SphP[i].b2.BPred[0] * SphP[i].b2.BPred[0] + SphP[i].b2.BPred[1] * SphP[i].b2.BPred[1] +
		 SphP[i].b2.BPred[2] * SphP[i].b2.BPred[2]) != 0)
	    sys.DivBerr +=
	      fabs(SphP[i].divB) * P[i].Hsml * P[i].Mass / (SphP[i].b2.BPred[0] * SphP[i].b2.BPred[0] +
							    SphP[i].b2.BPred[1] * SphP[i].b2.BPred[1] +
							    SphP[i].b2.BPred[2] * SphP[i].b2.BPred[2] +
							    1E-2) / SphP[i].Density;
#endif
#if defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB)
	  volweight += P[i].Mass / SphP[i].Density;
	  massweight += P[i].Mass;
#endif
	  sys.EnergyIntComp[0] += P[i].Mass * egyspec;
	}

      for(int j = 0; j < 3; j++)
	{
	  sys.MomentumComp[P[i].Type][j] += P[i].Mass * vel[j];
	  sys.CenterOfMassComp[P[i].Type][j] += P[i].Mass * P[i].Pos[j];
	}

      sys.AngMomentumComp[P[i].Type][0] += P[i].Mass * (P[i].Pos[1] * vel[2] - P[i].Pos[2] * vel[1]);
      sys.AngMomentumComp[P[i].Type][1] += P[i].Mass * (P[i].Pos[2] * vel[0] - P[i].Pos[0] * vel[2]);
      sys.AngMomentumComp[P[i].Type][2] += P[i].Mass * (P[i].Pos[0] * vel[1] - P[i].Pos[1] * vel[0]);
    }


  /* some the stuff over all processors */
  MPI_Reduce(&sys.MassComp[0], &SysState.MassComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyPotComp[0], &SysState.EnergyPotComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyIntComp[0], &SysState.EnergyIntComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyKinComp[0], &SysState.EnergyKinComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&sys.MomentumComp[0][0], &SysState.MomentumComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0,
	     MYMPI_COMM_WORLD);
  MPI_Reduce(&sys.AngMomentumComp[0][0], &SysState.AngMomentumComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0,
	     MYMPI_COMM_WORLD);
  MPI_Reduce(&sys.CenterOfMassComp[0][0], &SysState.CenterOfMassComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0,
	     MYMPI_COMM_WORLD);

#if defined(MAGNETIC_STATISTICS)
  MPI_Reduce(&sys.EnergyMag, &SysState.EnergyMag, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
#endif
#if defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB)
  MPI_Reduce(&sys.DivBerr, &SysState.DivBerr, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
#endif
#if defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB)
  MPI_Reduce(&volweight, &Volweight, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(&massweight, &Massweight, 1, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {

#if (defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB))
      SysState.DivBerr /= Volweight;
#endif

      for(int i = 0; i < 6; i++)
	SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] +
	  SysState.EnergyPotComp[i] + SysState.EnergyIntComp[i];

      SysState.Mass = SysState.EnergyKin = SysState.EnergyPot = SysState.EnergyInt = SysState.EnergyTot = 0;

      for(int j = 0; j < 3; j++)
	SysState.Momentum[j] = SysState.AngMomentum[j] = SysState.CenterOfMass[j] = 0;

      for(int i = 0; i < 6; i++)
	{
	  SysState.Mass += SysState.MassComp[i];
	  SysState.EnergyKin += SysState.EnergyKinComp[i];
	  SysState.EnergyPot += SysState.EnergyPotComp[i];
	  SysState.EnergyInt += SysState.EnergyIntComp[i];
	  SysState.EnergyTot += SysState.EnergyTotComp[i];

	  for(int j = 0; j < 3; j++)
	    {
	      SysState.Momentum[j] += SysState.MomentumComp[i][j];
	      SysState.AngMomentum[j] += SysState.AngMomentumComp[i][j];
	      SysState.CenterOfMass[j] += SysState.CenterOfMassComp[i][j];
	    }
	}

      for(int i = 0; i < 6; i++)
	for(int j = 0; j < 3; j++)
	  if(SysState.MassComp[i] > 0)
	    SysState.CenterOfMassComp[i][j] /= SysState.MassComp[i];

      for(int j = 0; j < 3; j++)
	if(SysState.Mass > 0)
	  SysState.CenterOfMass[j] /= SysState.Mass;


      for(int i = 0; i < 6; i++)
	{
	  SysState.CenterOfMassComp[i][3] = SysState.MomentumComp[i][3] = SysState.AngMomentumComp[i][3] = 0;
	  for(int j = 0; j < 3; j++)
	    {
	      SysState.CenterOfMassComp[i][3] +=
		SysState.CenterOfMassComp[i][j] * SysState.CenterOfMassComp[i][j];
	      SysState.MomentumComp[i][3] += SysState.MomentumComp[i][j] * SysState.MomentumComp[i][j];
	      SysState.AngMomentumComp[i][3] +=
		SysState.AngMomentumComp[i][j] * SysState.AngMomentumComp[i][j];
	    }
	  SysState.CenterOfMassComp[i][3] = sqrt(SysState.CenterOfMassComp[i][3]);
	  SysState.MomentumComp[i][3] = sqrt(SysState.MomentumComp[i][3]);
	  SysState.AngMomentumComp[i][3] = sqrt(SysState.AngMomentumComp[i][3]);
	}

      SysState.CenterOfMass[3] = SysState.Momentum[3] = SysState.AngMomentum[3] = 0;

      for(int j = 0; j < 3; j++)
	{
	  SysState.CenterOfMass[3] += SysState.CenterOfMass[j] * SysState.CenterOfMass[j];
	  SysState.Momentum[3] += SysState.Momentum[j] * SysState.Momentum[j];
	  SysState.AngMomentum[3] += SysState.AngMomentum[j] * SysState.AngMomentum[j];
	}

      SysState.CenterOfMass[3] = sqrt(SysState.CenterOfMass[3]);
      SysState.Momentum[3] = sqrt(SysState.Momentum[3]);
      SysState.AngMomentum[3] = sqrt(SysState.AngMomentum[3]);
    }

  /* give everyone the result, maybe the want to do something with it */
  MPI_Bcast(&SysState, sizeof(struct state_of_system), MPI_BYTE, 0, MYMPI_COMM_WORLD);
}
