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
#include "../System/communication.h"
#include "Matrix.hpp"


struct kernel_flux
{
  MyAtLeastDouble dx, dy, dz;
  MyAtLeastDouble r, vsig, sound_i, sound_j;
  MyAtLeastDouble dvx, dvy, dvz, vdotr2;
  MyAtLeastDouble wk_i, wk_j, dwk_i, dwk_j;
  MyAtLeastDouble h_i, h_j, dwk_ij, rho_ij_inv;
};


/*! Structure for communication during the flux computation. Holds data that is sent to other processors.
 */
static struct fluxdata_in
{
  MyIDType ID;
  MyLongDouble Pos[3];
  MyFloat GravAccel[3];
#ifdef PMGRID
  MyFloat GravPM[3];
#endif
  integertime Ti_begstep;
  MyFloat Hsml;
  MyFloat SoundSpeed;
  MyFloat Volume;
  WFluidVector<NUMDIMS> Wprim;
  MyFloat AlphaSlope[NUMDIMS+2];
  Vector<NUMDIMS> gradW[NUMDIMS+2];
  SquareMatrix<NUMDIMS> Bgrad;
  short int TimeBin;
  int NodeList[NODELISTLENGTH];
  MyFloat ConditionNumber;
}
*FluxDataIn, *FluxDataGet;


/*! Structure for communication during the flux computation. Returns partial results from node.
 */
static struct fluxdata_out
{
  QFluidVector<NUMDIMS> dQ;
  QFluidVector<NUMDIMS> dQdt;
}
 *FluxDataResult, *FluxDataOut;


/*! Create reduced particle data structure for communicating to other nodes for gradient computation */
void particle2in_flux(struct fluxdata_in *in, int i)
{
  in->ID = P[i].ID;
  for(int k = 0; k < 3; k++) {
    in->Pos[k] = P[i].Pos[k];
    in->GravAccel[k] = P[i].GravAccel[k];
#ifdef PMGRID
    in->GravPM[k] = P[i].GravPM[k];
#endif
  }
  in->Hsml = P[i].Hsml;
  in->TimeBin = P[i].TimeBin;
  in->Ti_begstep = P[i].Ti_begstep;
  const MyFloat u = SphP[i].InternalEnergyPred;
  in->SoundSpeed = eos->SoundSpeed(P[i].Mass*SphP[i].NumDens, u, All.cf_a3inv);
  in->Volume = MyFloat(1.0) / SphP[i].NumDens;
  for(int k = 0; k < NUMDIMS+2; k++) {
    in->Wprim[k] = SphP[i].Wprim[k];
    in->AlphaSlope[k] = SphP[i].AlphaSlope[k];
    in->gradW[k] = SphP[i].gradW[k];
  }
  in->Bgrad = SphP[i].Bgrad;
  in->ConditionNumber = SphP[i].ConditionNumber;
}


/*! Create data structure with partial sums/results to send back to particle host node */
void out2particle_flux(struct fluxdata_out *out, int i, int mode)
{
  if(P[i].Type == 0)
    {
      for (int k = 0; k < NUMDIMS+2; k++) {
        ASSIGN_ADD(SphP[i].dQ[k], out->dQ[k], mode);
        ASSIGN_ADD(SphP[i].dQdt[k], out->dQdt[k], mode);
      }
    }
}


// Returns 1 if particle is active (i.e. beginning of timestep).  Otherwise returns 0.
int flux_isactive(int n)
{
  if(P[n].TimeBin < 0)
    return 0;

  if(P[n].Type == 0)
    {
#ifdef WINDS
      if(SphP[n].DelayTime > 0)
	return 0;
#endif
      return 1;
    }
  
  return 0;
}

/*! combined function to calculate all relevant hydro changes based on MFM Godunov fluxes.*/
void mfm_fluxes_to_code(int i, integertime tstart_step, integertime tend_step)
{
  double dt_hydrokick, dt_entr_phys, dt_entr, dt_drift;

  if(All.ComovingIntegrationOn)
    {
      dt_hydrokick = get_hydrokick_factor(tstart_step, tend_step);
      dt_entr = (tend_step - tstart_step) * All.Timebase_interval;
      dt_entr_phys = (tend_step - tstart_step) * All.Timebase_interval / All.cf_hubble_a;
      dt_drift = get_drift_factor(tstart_step, tend_step);
    }
  else
    {
      dt_hydrokick = dt_entr = dt_entr_phys = dt_drift = (tend_step - tstart_step) * All.Timebase_interval;
    }

  // Update HydroAccel used in code
  for (int k = 0; k < NUMDIMS; k++) SphP[i].HydroAccel[k] = SphP[i].dQdt[k] * All.cf_atime /P[i].Mass * dt_entr_phys/dt_hydrokick;


  // Update DtEntropy / DtInternalEnergy used in code
  SphP[i].DtInternalEnergy = SphP[i].dQdt[NUMDIMS+1]/P[i].Mass * dt_entr_phys/dt_entr; // dQdt is internal energy in physical units
  if(All.ComovingIntegrationOn)
    {
      SphP[i].DtInternalEnergy -= 3*GAMMA_MINUS1 * SphP[i].InternalEnergyPred * All.cf_hubble_a * dt_entr_phys/dt_entr; // adiabatic hubble flux correction
    }
  SphP[i].DtEntropy = eos->GetEntropyFromU(SphP[i].Density, SphP[i].dQdt[NUMDIMS+1]*dt_entr_phys/P[i].Mass, All.cf_a3inv) / dt_entr;
  SphP[i].DtEntropy += GAMMA_MINUS1* eos->GetEntropyFromU(SphP[i].Density, SphP[i].InternalEnergyPred, All.cf_a3inv) * SphP[i].DivVel * dt_drift / dt_entr; // reduce by adiabatic energy change (only rough approximation!)
  
  //Energy-Entropy switch from gizmo.
  MyFloat grav_acc[3];
  MyFloat e_thermal,e_kinetic,e_potential;
  e_potential=0;
  for(int k = 0; k < 3; k++)
    {
      grav_acc[k] = All.cf_a2inv * P[i].GravAccel[k];
#ifdef PMGRID
      grav_acc[k] += All.cf_a2inv * P[i].GravPM[k];
#endif
      e_potential += grav_acc[k]*grav_acc[k];
    }
  e_potential = P[i].Mass * sqrt(e_potential) * ((MyFloat)0.5/*is KERNEL_CORE_SIZE*/*P[i].Hsml*All.cf_atime); // = M*|a_grav|*h (physical)
  e_kinetic = 0.5 * P[i].Mass * All.cf_a2inv * SphP[i].MaxVelSquareNgb;
  MyFloat dEnt_Gravity = 0; //only for MFV
  MyFloat InternalEnergy = SphP[i].InternalEnergyPred;
  MyFloat dEnt = InternalEnergy + SphP[i].DtInternalEnergy*dt_entr + dEnt_Gravity;
  e_thermal = DMAX(0.5*InternalEnergy, dEnt) * P[i].Mass;
  // this is the actual switch. in many cases, only one is sufficient (KE or potential).
  if(All.EpotSwitchFraction*e_potential + All.EkinSwitchFraction*e_kinetic > e_thermal) // recommended values on the order of 1e-2, use 0 for both to ignore the switch.
    {
      SphP[i].DtInternalEnergy = (SphP[i].Pressure/SphP[i].Density) * SphP[i].DivVel*All.cf_a2inv * dt_entr_phys/dt_entr;
      if(All.ComovingIntegrationOn)
	{
	  SphP[i].DtInternalEnergy -= 3*GAMMA_MINUS1 * SphP[i].InternalEnergyPred * All.cf_hubble_a * dt_entr_phys/dt_entr;
	}
      SphP[i].DtEntropy = (MyFloat) 0.0; // assume purely adiabatic changes
    }

  if(dt_entr_phys == 0)
    {
      for(int k = 0; k < NUMDIMS; k++) SphP[i].HydroAccel[k] = (MyFloat)0.0;
      SphP[i].DtInternalEnergy = (MyFloat)0.0;
      SphP[i].DtEntropy = (MyFloat)0.0;
    }
}
/*! synchronize entropy with internal energy, which is used for the evolution in mfm */
void mfm_synchronize_entropy(int i)
{
  SphP[i].Entropy = eos->GetEntropyFromU(SphP[i].Density, SphP[i].InternalEnergy, All.cf_a3inv);
  SphP[i].EntropyPred = eos->GetEntropyFromU(SphP[i].Density, SphP[i].InternalEnergyPred, All.cf_a3inv);
}

/*! scale dQdt according to change in mass */
void mfm_update_mass(int i, MyFloat old_mass, MyFloat new_mass)
{
  MyFloat factor = new_mass / old_mass;
  // components with information about HydroAccel
  for(int j = 0; j < NUMDIMS; j++) {
    SphP[i].dQdt[j] *= factor;
  }
  // component with information about InternalEnergy / Entropy
  SphP[i].dQdt[NUMDIMS+1] *= factor;
}
/*! scale dQdt according to change in mass */
void mfm_add_mass(int i, MyFloat mass_change)
{
  if ((P[i].Type == 0) && (P[i].Mass > 0) && (P[i].Mass + mass_change > 0)) {
    MyFloat factor = (P[i].Mass + mass_change) / P[i].Mass;
    // components with information about HydroAccel
    for(int j = 0; j < NUMDIMS; j++) {
      SphP[i].dQdt[j] *= factor;
    }
    // component with information about InternalEnergy / Entropy
    SphP[i].dQdt[NUMDIMS+1] *= factor;
  }
}

/*! \file mfm_fluxes.cpp
 *  \brief Compute Godunov fluxes for all active particles on this step
 */


// This function computes the gradients for all active particles
void mfm_fluxes(void)
{
  int ndone, ndone_flag, dummy;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  long long n_exported = 0;

  // Display information message that this step is executed on Task 0
  VERBOSE(3, "Calculating MFM fluxes\n");
  
  // TODO : Check with Klaus about DENSITY_LESS_NGB_FACTOR
  int NTaskTimesNumPart = ((int) (maxThreads)) * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  size_t MyBufferSize = KD_BUFFER_MANAGEMENT * FreeBytes / (1024.0 * 1024.0);
  
  All.BunchSize =
    (int) ((MyBufferSize  * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct fluxdata_in) + sizeof(struct fluxdata_out) +
                                             sizemax(sizeof(struct fluxdata_in), sizeof(struct fluxdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  VERBOSE(2, "MFM fluxes: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
	  FreeBytes / (1024.0 * 1024.0));
  
  CPU_Step[CPU_MFMFLUXMISC] += measure_time();
  t0 = second();


  // begin with this index
  NextParticle = FirstActiveParticle;

  int group_particles = 1;
  do
    {
      BufferFullFlag = 0;
      Nexport = 0;
      int save_NextParticle = NextParticle;

      tstart = second();

#ifdef _OPENMP
#pragma omp parallel
#endif
      flux_evaluate_primary(&ThisThread);  // do local particles and prepare export list
      
      tend = second();
      timecomp1 += timediff(tstart, tend);

      /* if we are not grouping particles and SendRecv fails, then we panic as there is no memory to process even a single particles */
      int panic_on_error = !group_particles;
      int is_send_receive_ok = prepare_sendrec_offset_send(save_NextParticle, panic_on_error); //1 means panic on bufer too small
      if(!is_send_receive_ok) {
	group_particles = 0;
	continue;
      }

      FluxDataGet =
        (struct fluxdata_in *) mymalloc("FluxDataGet", Nimport * sizeof(struct fluxdata_in));
      FluxDataIn =
        (struct fluxdata_in *) mymalloc("FluxDataIn", Nexport * sizeof(struct fluxdata_in));

      /* prepare particle data for export */
      for(int j = 0; j < Nexport; j++)
        {
          int place = DataIndexTable[j].Index;

          particle2in_flux(&FluxDataIn[j], place);

          memcpy(FluxDataIn[j].NodeList,
                 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

      // exchange particle data
      tstart = second();
      for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          int sendTask = ThisTask;
          int recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  // get the particles */
                  MPI_Sendrecv(&FluxDataIn[Send_offset[recvTask]],
                   Send_count[recvTask] * sizeof(struct fluxdata_in), MPI_BYTE,
                   recvTask, TAG_FLUX_A,
                   &FluxDataGet[Recv_offset[recvTask]],
                   Recv_count[recvTask] * sizeof(struct fluxdata_in), MPI_BYTE,
                   recvTask, TAG_FLUX_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      tend = second();
      timecommsumm1 += timediff(tstart, tend);

      myfree(FluxDataIn);
      FluxDataResult =
       (struct fluxdata_out *) mymalloc("FluxDataResult", Nimport * sizeof(struct fluxdata_out));
      FluxDataOut =
        (struct fluxdata_out *) mymalloc("FluxDataOut", Nexport * sizeof(struct fluxdata_out));


      // now do the particles that were sent to us
      tstart = second();
      NextJ = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
      flux_evaluate_secondary(&ThisThread);
      
      tend = second();
      timecomp2 += timediff(tstart, tend);

      if(NextParticle < 0)
        ndone_flag = 1;
      else
        ndone_flag = 0;
      
      tstart = second();
      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
      tend = second();
      timewait2 += timediff(tstart, tend);

      // get the result
      tstart = second();
      for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          int sendTask = ThisTask;
          int recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* send the results */
                  MPI_Sendrecv(&FluxDataResult[Recv_offset[recvTask]],
                               Recv_count[recvTask] * sizeof(struct fluxdata_out),
                               MPI_BYTE, recvTask, TAG_FLUX_B,
                               &FluxDataOut[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct fluxdata_out),
                               MPI_BYTE, recvTask, TAG_FLUX_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
      fflush(stdout);
      
      tend = second();
      timecommsumm2 += timediff(tstart, tend);


      // add the result to the local particles
      tstart = second();
      for(int j = 0; j < Nexport; j++)
        {
          int place = DataIndexTable[j].Index;
          out2particle_flux(&FluxDataOut[j], place, 1);
        }
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(FluxDataOut);
      myfree(FluxDataResult);
      myfree(FluxDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  //QFluidVector<NUMDIMS> dQtot;

  // do final operations on results
  tstart = second();
#ifdef _OPENMP
  int il;
#pragma omp parallel for
  for(il = 0; il < NActivePart; il ++)
    {
      int i = ActiveParticleList[il];
#else
  for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif
      if (flux_isactive(i))
        {
	    double dt_entr_phys;
	    integertime ti_step, tstart_step, tend_step;
	    ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
	    tstart_step = P[i].Ti_begstep /* + ti_step / 2*/;   /* midpoint of step */
	    tend_step = P[i].Ti_begstep + ti_step; /* end of step */
	    
	    if(All.ComovingIntegrationOn)
	      {
		dt_entr_phys = (tend_step - tstart_step) * All.Timebase_interval / All.cf_hubble_a;
	      }
	    else
	      {
		dt_entr_phys = (tend_step - tstart_step) * All.Timebase_interval;
	      }
          
	  for (int k = 0; k < NUMDIMS; k++)
		{
		  SphP[i].dQdt[NUMDIMS+1] -= (SphP[i].Wprim[k]/All.cf_atime + (MyFloat)0.5 * SphP[i].dQ[k]/P[i].Mass) * SphP[i].dQdt[k]; // from total energy flux to internal energy flux (physical units) // without the dQ term, it produces too high internal energies for the shock.
		  SphP[i].dQ[NUMDIMS+1] -= (SphP[i].Wprim[k]/All.cf_atime + (MyFloat)0.5 * SphP[i].dQ[k]/P[i].Mass) * SphP[i].dQ[k]; // from total energy flux to internal energy flux (physical units) // without the dQ term, it produces too high internal energies for the shock.
		}
	  if(dt_entr_phys > 0)
	  SphP[i].dQdt = SphP[i].dQ * ((MyFloat)1.0 / dt_entr_phys); // this should be replaced at some point by a more direct calculation without going via the timesteps.
	  else
	    SphP[i].dQdt = (MyFloat)0.0;
#ifdef WINDS
	  if (SphP[i].DelayTime > 0) /* decouple wind particles */
	    {
	      SphP[i].dQdt = (MyFloat)0.0;
	      SphP[i].dQdt[NUMDIMS+1] = eos->GetDtQFromDtEntropy((MyFloat)0.0,SphP[i].Entropy,SphP[i].InternalEnergy, P[i].Mass, SphP[i].DivVel,All.cf_hubble_a);
	    }
#endif
	  mfm_fluxes_to_code(i, tstart_step, tend_step);
        }

    }

  tend = second();
  timecomp1 += timediff(tstart, tend);

  // collect some timing information
  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_MFMFLUXCOMPUTE] += timecomp;
  CPU_Step[CPU_MFMFLUXWAIT] += timewait;
  CPU_Step[CPU_MFMFLUXCOMM] += timecomm;
  CPU_Step[CPU_MFMFLUXMISC] += timeall - (timecomp + timewait + timecomm);

}



// This function represents the core of the flux computation. The
// target particle may either be local, or reside in the communication buffer.
int flux_evaluate(int target, int mode, int *exportflag, int *exportnodecount,
                      int *exportindex, int *ngblist)
{
  int n, listindex = 0;
  int startnode;
  MyAtLeastDouble h, h2, hinv, hinv3, hinv4;
  MyAtLeastDouble hinvi, hinvj;

  struct kernel_flux kernel;
  struct fluxdata_in local;
  struct fluxdata_out out;
  memset(&out, 0, sizeof(struct fluxdata_out));

  if(mode == 0) {
    particle2in_flux(&local, target);
  }
  else {
    local = FluxDataGet[target];
  }

  double dt_gravkick, dt_hydrokick, dt_entr_phys, dt_entr_inv;
  integertime ti_step, tstart_step, tend_step;
  ti_step = local.TimeBin ? (((integertime) 1) << local.TimeBin) : 0;
  tstart_step = local.Ti_begstep /*+ ti_step / 2*/;   /* midpoint of step */
  tend_step = local.Ti_begstep + ti_step; /* end of step */
  
  if(All.ComovingIntegrationOn)
    {
      dt_gravkick = get_gravkick_factor(tstart_step, tend_step);
      dt_hydrokick = get_hydrokick_factor(tstart_step, tend_step);
      dt_entr_phys = (tend_step - tstart_step) * All.Timebase_interval / All.cf_hubble_a; // used for fluxes, that are all in physical units.
    }
  else
    {
      dt_gravkick = dt_hydrokick = dt_entr_phys = (tend_step - tstart_step) * All.Timebase_interval;
    }
  dt_entr_inv = MyFloat(1.0)/dt_entr_phys;
  h2 = local.Hsml * local.Hsml;
  kernel.h_i = local.Hsml;
  //kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

  if(mode == 0)
    {
      startnode = All.MaxPart;    // root node
    }
  else
    {
      startnode = FluxDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;    // open it
    }


  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          int numngb_inbox = ngb_treefind_pairs_threads(&(local.Pos[0]), local.Hsml, target,
                                                       &startnode, mode, exportflag,
                                                       exportnodecount, exportindex, ngblist);

          if(numngb_inbox < 0)
            return -1;

          for(int n = 0; n < numngb_inbox; n++)
            {
              int j = ngblist[n];

#ifdef MFM_BIN_FLUXES
              // Skip neighbours on shorter steps
              if (P[j].TimeBin > local.TimeBin) continue; // does not work for dQdt yet
#endif
	      
	      // avoid self-interaction
	      //if (local.ID == P[j].ID) continue;
              
#if defined(BLACK_HOLES)
              if(P[j].Mass == 0)
                continue;
#endif
#ifdef WINDS
	      if (SphP[j].DelayTime > 0) /* decouple wind particles */
		continue;
#endif

              Vector<NUMDIMS> dr;
              for (int k = 0; k < NUMDIMS; k++) dr[k] = P[j].Pos[k] - local.Pos[k];
#ifdef PERIODIC
              dr[0] = NEAREST_X(dr[0]);
#if NUMDIMS > 1
              dr[1] = NEAREST_Y(dr[1]);
#endif
#if NUMDIMS == 3
              dr[2] = NEAREST_Z(dr[2]);
#endif
#endif
              MyAtLeastDouble r2 = dr.GetMagnitudeSquared();
              kernel.h_j = P[j].Hsml;

              // Try to avoid self-interaction terms
              if (sqrt(r2)/kernel.h_i < 1.0e-8 || sqrt(r2)/kernel.h_j < 1.0e-8) continue;

              if(r2 < kernel.h_i * kernel.h_i || r2 < kernel.h_j * kernel.h_j)
                {
                  kernel.r = sqrt(r2);

                  if(kernel.r < kernel.h_i)
                    {
                      kernel_hinv(kernel.h_i, &hinvi, &hinv3, &hinv4);
                      MyAtLeastDouble u = kernel.r * hinvi;
                      kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, 0);
                    }
                  else
                    {
                      hinvi = 0;
                      kernel.dwk_i = 0;
                      kernel.wk_i = 0;
                    }

                  if(kernel.r < kernel.h_j)
                    {
                      kernel_hinv(kernel.h_j, &hinvj, &hinv3, &hinv4);
                      MyAtLeastDouble u = kernel.r * hinvj;
                      kernel_main(u, hinv3, hinv4, &kernel.wk_j, &kernel.dwk_j, 0);
                    }
                  else
                    {
                      hinvj = 0;
                      kernel.dwk_j = 0;
                      kernel.wk_j = 0;
                    }

                  // Skip interaction if both kernel terms are zero
                  //if (kernel.wk_i < 1.0e-10 && kernel.wk_j < 1.0e-10) continue;

                  // Compute psi-tilda values (for integral gradient terms)
                  MyFloat vol_j = MyFloat(1.0) / SphP[j].NumDens;
                  Vector<NUMDIMS> psitildai;
                  Vector<NUMDIMS> psitildaj;
		  
                  for (int k = 0; k < NUMDIMS; k++) {
                    for (int kk = 0; kk < NUMDIMS; kk++) {
                      psitildai[k] -= (SphP[j].Bgrad(k,kk)*dr[kk]*kernel.wk_j*hinvj*vol_j);
                      psitildaj[k] += (local.Bgrad(k,kk)*dr[kk]*kernel.wk_i*hinvi*local.Volume);
                    }
                  }

                  // Compute the face pseudo-area, Aij, and corresponding unit vector, Aunit

                  Vector<NUMDIMS> Aij = psitildaj*local.Volume - psitildai*vol_j;
                  MyFloat Amag = Aij.GetMagnitude();
                  Vector<NUMDIMS> Aunit = Aij.GetUnitVector();
		  		  
                  // Compute position and velocity of the working face (i.e. where the Riemann problem is solved)
                  //Vector<NUMDIMS> rface;
                  Vector<NUMDIMS> vface;
		  //s_star_ij = 0; /*first order*/
		  MyFloat s_star_ij = 0.5 * kernel.r * (P[j].Hsml - local.Hsml) / (local.Hsml + P[j].Hsml); /*second order*/
		  MyFloat s_i = s_star_ij - 0.5 * kernel.r;
		  MyFloat s_j = s_star_ij + 0.5 * kernel.r;
                  Vector<NUMDIMS> drFacei,drFacej;
                  for (int k = 0; k < NUMDIMS; k++) {
                    //rface[k] = MyFloat(0.5)*(2.0*local.Pos[k] + dr[k]);
		    //drFacej[k] = rface[k] - local.Pos[k] - dr[k];
		    //drFacei[k] = rface[k] - local.Pos[k];
		    drFacej[k] = -dr[k]/kernel.r * s_j;
		    drFacei[k] = -dr[k]/kernel.r * s_i;
 		    vface[k] = (s_j*local.Wprim[k] - s_i*SphP[j].Wprim[k]) / kernel.r; // allows for face to be off-center (second-order)
                  }

                  if( (local.ConditionNumber > 1e12) || (SphP[j].ConditionNumber > 1e12 /*gizmo value*/) || (Vector<NUMDIMS>::DotProduct(Aunit,drFacei) < 0.0) )
		    {
		      MyFloat wt_i, wt_j;
		      Aunit = drFacei.GetUnitVector();
#ifdef COOLING
		      if(abs(local.Volume - vol_j)/std::min(local.Volume, vol_j) > 1.25)
			{
			  wt_i = wt_j = (local.Volume*P[j].Hsml + vol_j*local.Hsml) / (local.Hsml + P[j].Hsml);
			}
#else
		      if(abs(local.Volume - vol_j)/std::min(local.Volume, vol_j) > 1.5)
			{
			  wt_i = wt_j = (MyFloat)2.0 * local.Volume * vol_j / (local.Volume + vol_j);
			}
#endif
		      else
			{
			  wt_i = local.Volume;
			  wt_j = vol_j;
			}
		      Amag = -(wt_i*local.Volume*kernel.dwk_i + wt_j*vol_j*kernel.dwk_j) / kernel.r;
		    }

		  // Skip interaction if both kernel terms are zero
                  if (Amag == 0) continue;

                  WFluidVector<NUMDIMS> dW;
                  WFluidVector<NUMDIMS> Wdot;
                  Vector<NUMDIMS> gradWi[NUMDIMS+2],gradWj[NUMDIMS+2];

                  // Compute primitive fluid vectors for the left-hand state of the Riemann problem
                  WFluidVector<NUMDIMS> Wi(local.Wprim);
                  for (int k = 0; k < NUMDIMS + 2; k++)
                  {
		    gradWi[k] = local.gradW[k];
                    gradWi[k] *= local.AlphaSlope[k];
                    dW[k] = Vector<NUMDIMS>::DotProduct(gradWi[k], drFacei);
		    Wi[k] += dW[k];
                  }

                  // Compute primitive fluid vectors for the right-hand state of the Riemann problem
		  const MyFloat u = SphP[j].InternalEnergyPred;
                  const MyFloat SoundSpeed_j = eos->SoundSpeed(P[j].Mass*SphP[j].NumDens, u , All.cf_a3inv);
                  WFluidVector<NUMDIMS> Wj(SphP[j].Wprim);
		  for (int k = 0; k < NUMDIMS + 2; k++) 
                  {
		    gradWj[k] = SphP[j].gradW[k];
                    gradWj[k] *= SphP[j].AlphaSlope[k];
                    dW[k] = Vector<NUMDIMS>::DotProduct(gradWj[k], drFacej);
                    Wj[k] += dW[k];
                  }


#if (defined(SLOPE_LIMITER_GIZMO_FANCY) || defined(SLOPE_LIMITER_GIZMO_PAIRWISE))
		  MyFloat phimin, phimax, phimed,phimax_eff,phimin_eff,fac,phimed_max,phimed_min;
		  MyFloat fac_minmax = 0.5;
		  MyFloat fac_meddev = 0.375;
#if SLOPE_LIMITER_TOLERANCE == 2
		  fac_minmax = 0.75;
		  fac_meddev = 0.4;
#elif SLOPE_LIMITER_TOLERANCE == 0		  
		  fac_minmax = 0.0;
		  fac_meddev = 0.0;
#endif
#if defined(SLOPE_LIMITER_GIZMO_PAIRWISE)
		  fac_minmax = 0.5;
		  fac_meddev = 0.25;
#endif

		  WFluidVector<NUMDIMS> Wi_orig(local.Wprim);
		  WFluidVector<NUMDIMS> Wj_orig(SphP[j].Wprim);
		  for (int k = 0; k < NUMDIMS; k++)
		    {
		      Wi[k] -= vface[k];
		      Wj[k] -= vface[k];
		      Wi_orig[k] -= vface[k];
		      Wj_orig[k] -= vface[k];
		    }
		  
		  for (int k = 0; k < NUMDIMS + 2; k++)
		    {
		      phimed = 0.5 * (Wi_orig[k] + Wj_orig[k]);
		      if(Wi_orig[k] < Wj_orig[k])
			{phimax = Wj_orig[k]; phimin = Wi_orig[k]; }
		      else
			{phimax = Wi_orig[k]; phimin = Wj_orig[k]; }
		      fac = fac_minmax * (phimax - phimin);
		      phimax_eff = phimax + fac;
		      phimin_eff = phimin - fac;
		      if (phimax < 0) {if(phimax_eff>0) phimax_eff=phimax*phimax/(phimax-(phimax_eff-phimax));}
		      if (phimin > 0) {if(phimin_eff<0) phimin_eff=phimin*phimin/(phimin+(phimin-phimin_eff));}
		      fac = fac_meddev * (phimax - phimin);
		      phimed_max = phimed + fac;
		      phimed_min = phimed - fac;
		      if(phimed_max > phimax_eff) phimed_max = phimax_eff;
		      if(phimed_min < phimin_eff) phimed_min = phimin_eff;
		      if(local.Wprim[k] < SphP[j].Wprim[k])
			{
			  if(Wj[k] < phimin_eff) Wj[k] = phimin_eff;
			  if(Wj[k] > phimed_max) Wj[k] = phimed_max;
			  if(Wi[k] > phimax_eff) Wi[k] = phimax_eff;
			  if(Wi[k] < phimed_min) Wi[k] = phimed_min;
			} else {
			  if(Wj[k] > phimax_eff) Wj[k] = phimax_eff;
			  if(Wj[k] < phimed_min) Wj[k] = phimed_min;
			  if(Wi[k] < phimin_eff) Wi[k] = phimin_eff;
			  if(Wi[k] > phimed_max) Wi[k] = phimed_max;
		      }

		      /*if(Vector<NUMDIMS>::DotProduct(local.gradW[k], drFacei) != 0) for (int kk = 0; kk < NUMDIMS; kk++) gradWi[k][kk] = local.gradW[k][kk]* (Wi[k]-local.Wprim[k]) / Vector<NUMDIMS>::DotProduct(local.gradW[k], drFacei);
			if(Vector<NUMDIMS>::DotProduct(SphP[j].gradW[k], drFacej) != 0) for (int kk = 0; kk < NUMDIMS; kk++) gradWj[k][kk] = SphP[j].gradW[k][kk]* (Wj[k]-SphP[j].Wprim[k]) / Vector<NUMDIMS>::DotProduct(SphP[j].gradW[k], drFacej);*/
		    }
		  for (int k = 0; k < NUMDIMS; k++)
		    {
		      Wi[k] += vface[k];
		      Wj[k] += vface[k];
		      //Wi_orig[k] += vface[k];
		      //Wj_orig[k] += vface[k];
		    }
#endif

#ifndef NOMUSCL
		  // Calculate the predicted left-hand state at the half-timestep
                  Wdot = WFluidVector<NUMDIMS>::ComputePrimitiveTimeDerivative(local.SoundSpeed*local.SoundSpeed, Wi, gradWi, vface);
                  for (int k = 0; k < NUMDIMS+2; k++) Wi[k] += MyFloat(0.5)*dt_entr_phys*Wdot[k];

                  // Calculate the predicted right-hand state at the half-timestep
                  Wdot = WFluidVector<NUMDIMS>::ComputePrimitiveTimeDerivative(SoundSpeed_j*SoundSpeed_j, Wj, gradWj, vface);
		  for (int k = 0; k < NUMDIMS+2; k++) Wj[k] += MyFloat(0.5)*dt_entr_phys*Wdot[k]; //need to check this!  
                  for (int k = 0; k < NUMDIMS; k++) {
		    Wi[k] += MyFloat(0.5)*dt_gravkick*local.GravAccel[k];
		    Wj[k] += MyFloat(0.5)*dt_gravkick*P[j].GravAccel[k]; // compare to the Muscl implementation in Gandalf.
		  }
#ifdef PMGRID
		  for (int k = 0; k < NUMDIMS; k++) {
		    Wi[k] += MyFloat(0.5)*dt_gravkick*local.GravPM[k];
		    Wj[k] += MyFloat(0.5)*dt_gravkick*P[j].GravPM[k];
		  }
#endif
		  s_star_ij = 0.5 * kernel.r * (P[j].Hsml - local.Hsml) / (local.Hsml + P[j].Hsml); //second order
                  s_i = s_star_ij - 0.5 * kernel.r;
                  s_j = s_star_ij + 0.5 * kernel.r;
                  for (int k = 0; k < NUMDIMS; k++) {
		    s_i += 0.5 * dt_entr_phys * local.Wprim[k]*dr[k] / kernel.r;
		    s_j += 0.5 * dt_entr_phys * SphP[j].Wprim[k]*dr[k] / kernel.r;
                  }
                  for (int k = 0; k < NUMDIMS; k++) {
                    //vface[k] = MyFloat(0.5)*(Wi[k] + Wj[k]);
                    vface[k] = (s_j*Wi[k] - s_i*Wj[k]) / kernel.r; // allows for face to be off-center (second-order)
                  }
#endif //NOMUSCL
		  for (int k = NUMDIMS; k < NUMDIMS+2; k++) if (Wi[k] <= (MyFloat)1.0e-30 or !std::isfinite(Wi[k])) Wi[k]=(MyFloat)1.0e-30;
		  for (int k = NUMDIMS; k < NUMDIMS+2; k++) if (Wj[k] <= (MyFloat)1.0e-30 or !std::isfinite(Wj[k])) Wj[k]=(MyFloat)1.0e-30;

		  if(All.ComovingIntegrationOn)
		    {
		      // Convert to physical units
		      for (int k = 0; k < NUMDIMS; k++)
			{
			  vface[k] /= All.cf_atime; // convert to physical units
			  Wi[k] /= All.cf_atime;
			  Wj[k] /= All.cf_atime;
			}
		      int irho = NUMDIMS;
		      int ipress = NUMDIMS+1;
		      Wi[irho] *= All.cf_a3inv; //density
		      Wj[irho] *= All.cf_a3inv;
		      Wi[ipress] *= All.cf_a3inv / All.cf_afac1; //pressure
		      Wj[ipress] *= All.cf_a3inv / All.cf_afac1;
		      
		      Amag *= All.cf_atime*All.cf_atime;
		    }
		  
		  // Calculate Godunov flux using the selected Riemann solver
		  FluxVector<NUMDIMS> flux = fluxSolver->ComputeFluxes(true, Aunit, vface, Wi, Wj);

		  QFluidVector<NUMDIMS> Qflux = flux*Amag*dt_entr_phys;
		  
		  out.dQ -= Qflux;
                  out.dQdt -= flux*Amag;

#ifdef MFM_BIN_FLUXES
                  // If neighbour is in longer bin, then record opposite flux here
                  if (P[j].TimeBin < local.TimeBin)
		    {
     		    #ifdef _OPENMP
		    #pragma omp critical(_fluxes_)
		    #endif
		    {
                    SphP[j].dQ += Qflux;

                    // If neighbour is also at the beginning of its step, then record dQdt
                    if (flux_isactive(j))
                    {
                      SphP[j].dQdt += flux*Amag; // this does not work yet. First just use dQ, but solve later.
                    }
		    } // only in combination with timebin criterion in the beginning! Have to double check this!
		    }
#endif
                }
	      
	    }
        }
      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = FluxDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Nodes[startnode].u.d.nextnode;    // open it
            }
        }

    }

  
  if(mode == 0)
    {
      out2particle_flux(&out, target, 0);
      //out2particle_flux(&out, target, 1); // this casues some numerical problem in a parallelized version ???
    }
  else
    {
      FluxDataResult[target] = out;
    }

  return 0;
}


void *flux_evaluate_primary(void *p)
{
  int thread_id = *(int *) p;
  int i;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;

#ifdef DENSITY_LESS_NGB_FACTOR
  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
#else
  ngblist = Ngblist + thread_id * NumPart;
#endif
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  // Note: exportflag is local to each thread
  for(int j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
        if(BufferFullFlag != 0 || NextParticle < 0)
          {
            exitFlag = 1;
          }
        else
          {
#ifdef _OPENMP
            do
              {
                i = NextParticle;
                NextParticle = NextActiveParticle[NextParticle];
                if(NextParticle < 0)
                  break;
                ProcessedFlag[i] = 1;
              }
            while(!flux_isactive(i));
              ProcessedFlag[i] = 0;
#else
            i = NextParticle;
            ProcessedFlag[i] = 0;
            NextParticle = NextActiveParticle[NextParticle];
#endif
          }
      }
      if(exitFlag)
	break;

      if(flux_isactive(i))
        {
          if(flux_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
            break;		/* export buffer has filled up */
        }

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *flux_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;

  int j, dummy, *ngblist;

#ifdef DENSITY_LESS_NGB_FACTOR
  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
#else
  ngblist = Ngblist + thread_id * NumPart;
#endif


  while(1)
    {
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
        j = NextJ;
        NextJ++;
      }

      if (j >= Nimport) break;

      flux_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}
