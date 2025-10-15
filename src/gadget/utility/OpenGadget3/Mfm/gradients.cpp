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


/* Structure for properties only needed within this function.
 */
static struct grad_store_out{
  Vector<NUMDIMS> gradWSph[NUMDIMS+2]; /*!< SPH computation of primitive gradients */
  MyFloat neighbour_mass;
} * GradStoreOut;


/*! Structure for communication during the gradient computation. Holds data that is sent to other processors.
 */
static struct graddata_in
{
  MyLongDouble Pos[3];
  //MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Volume;
  MyFloat NumDens;
  MyFloat SoundSpeed;
  WFluidVector<NUMDIMS> Wprim;
  int NodeList[NODELISTLENGTH];
}
*GradDataIn, *GradDataGet;


/*! Structure for communication during the gradient computation. Returns partial results from node.
 */
static struct graddata_out
{
  SquareMatrix<NUMDIMS> Ematrix;
  Vector<NUMDIMS> gradW[NUMDIMS+2];
  Vector<NUMDIMS> gradWSph[NUMDIMS+2];
  MyFloat MaxSignalVel;
#if MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009 || MFM_SLOPE_LIMITER == LIMITER_SCALAR || MFM_SLOPE_LIMITER == LIMITER_GIZMO
  WFluidVector<NUMDIMS> Wmin;
  WFluidVector<NUMDIMS> Wmax;
#endif
#if MFM_SLOPE_LIMITER == LIMITER_SCALAR || MFM_SLOPE_LIMITER == LIMITER_GIZMO
  MyFloat DistNgbSqdMax;
#endif
  Vector<NUMDIMS> Vel_smooth;
  MyFloat neighbour_mass;
  MyFloat MaxVelSquareNgb;
}
 *GradDataResult, *GradDataOut;

/*! Create reduced particle data structure for communicating to other nodes for gradient computation */
void particle2in_gradient(struct graddata_in *in, int i)
{
  for(int k = 0; k<3; k++) {
    in->Pos[k] = P[i].Pos[k];
  }
  in->Hsml = P[i].Hsml;
  in->NumDens = SphP[i].NumDens;
  in->Volume  = MyFloat(1.0) / SphP[i].NumDens;
  const MyFloat u = SphP[i].InternalEnergyPred;
  in->SoundSpeed = eos->SoundSpeed(P[i].Mass*SphP[i].NumDens, u, All.cf_a3inv);
  for(int k = 0; k < NUMDIMS+2; k++) {
    in->Wprim[k] = SphP[i].Wprim[k];
  }
  assert(std::isnormal(in->Hsml));
  assert(std::isnormal(in->NumDens));
}


/*! Create data structure with partial sums/results to send back to particle host node */
void out2particle_gradient(struct graddata_out *out, int i, int mode)
{
  if(P[i].Type == 0)
    {
      for (int j = 0; j < NUMDIMS; j++) {
        for (int k = 0; k < NUMDIMS; k++) {
          ASSIGN_ADD(SphP[i].Bgrad(j,k), out->Ematrix(j,k), mode);
        }
      }
      for (int j = 0; j < NUMDIMS+2; j++) {
        for (int k = 0; k < NUMDIMS; k++) {
          ASSIGN_ADD(SphP[i].gradW[j][k], out->gradW[j][k], mode);
          ASSIGN_ADD(GradStoreOut[i].gradWSph[j][k], out->gradWSph[j][k], mode);
        }
      }
      ASSIGN_MAX(SphP[i].MaxSignalVel, out->MaxSignalVel, mode);
#if MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009 || MFM_SLOPE_LIMITER == LIMITER_SCALAR || MFM_SLOPE_LIMITER == LIMITER_GIZMO
      for (int k = 0; k < NUMDIMS + 2; k++) {
        ASSIGN_MIN(SphP[i].Wmin[k], out->Wmin[k], mode);
        ASSIGN_MAX(SphP[i].Wmax[k], out->Wmax[k], mode);
      }
#endif
#if MFM_SCALAR_LIMITER == LIMITER_SCALAR || MFM_SLOPE_LIMITER == LIMITER_GIZMO
      ASSIGN_MAX(SphP[i].DistNgbSqdMax, out->DistNgbSqdMax, mode);
#endif
      ASSIGN_MAX(SphP[i].MaxVelSquareNgb, out->MaxVelSquareNgb, mode);
      for (int k = 0; k < NUMDIMS; k++) {
	ASSIGN_ADD(SphP[i].Vel_smooth[k], out->Vel_smooth[k], mode);
      }
      ASSIGN_ADD(GradStoreOut[i].neighbour_mass, out->neighbour_mass, mode);
    }
}


// Returns 1 if particle is active (i.e. beginning of timestep).  Otherwise returns 0.
int gradient_isactive(int n)
{
  if(P[n].TimeBin < 0)
    return 0;


  if(P[n].Type == 0)
    {
#ifdef WINDS
      if (SphP[n].DelayTime > 0) /* decouple wind particles */
        return 0;
#endif
    return 1;
    }

  return 0;
}


/*! \file gradients.cpp
 *  \brief Gradient computation for extrapolating state vectors to Riemann surface for 2nd order
 */


// This function computes the gradients for all active particles
void gradients(void)
{
  int // i, j, k,
    ndone, ndone_flag, dummy;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  long long n_exported = 0;

  // Display information message that this step is executed on Task 0
  VERBOSE(3, "Updating gradients for primitive fluid vector (Wprim)\n");
  
  // Reset wakeups
if (WAKEUP >0){
  #ifdef _OPENMP
    #pragma omp parallel for
  #endif
  for(int i = 0; i < N_gas; i++){
    SphP[i].wakeup = TIMEBINS;
  }
}

  GradStoreOut = (grad_store_out *) mymalloc("GradStoreOut", NumPart * sizeof(grad_store_out));

  // TODO : Check with Klaus about DENSITY_LESS_NGB_FACTOR
  int NTaskTimesNumPart = ((int) (maxThreads)) * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  size_t MyBufferSize = KD_BUFFER_MANAGEMENT * FreeBytes / (1024.0 * 1024.0);
  
  All.BunchSize =
    (int) ((MyBufferSize  * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct graddata_in) + sizeof(struct graddata_out) +
                                             sizemax(sizeof(struct graddata_in), sizeof(struct graddata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  VERBOSE(2, "Gradients: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
	  FreeBytes / (1024.0 * 1024.0));
  
  CPU_Step[CPU_MFMGRADMISC] += measure_time();
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
        gradient_evaluate_primary(&ThisThread);  // do local particles and prepare export list

      tend = second();
      timecomp1 += timediff(tstart, tend);

      /* if we are not grouping particles and SendRecv fails, then we panic as there is no memory to process even a single particles */
      int panic_on_error = !group_particles;
      int is_send_receive_ok = prepare_sendrec_offset_send(save_NextParticle, panic_on_error); // 1 means panic on buffer too small
      if(!is_send_receive_ok) {
	group_particles = 0;
	continue;
      }

      GradDataGet =
        (struct graddata_in *) mymalloc("GradDataGet", Nimport * sizeof(struct graddata_in));
      GradDataIn =
        (struct graddata_in *) mymalloc("GradDataIn", Nexport * sizeof(struct graddata_in));

      /* prepare particle data for export */
      for(int j = 0; j < Nexport; j++)
        {
          int place = DataIndexTable[j].Index;

          particle2in_gradient(&GradDataIn[j], place);

          memcpy(GradDataIn[j].NodeList,
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
                  MPI_Sendrecv(&GradDataIn[Send_offset[recvTask]],
                   Send_count[recvTask] * sizeof(struct graddata_in), MPI_BYTE,
                   recvTask, TAG_GRAD_A,
                   &GradDataGet[Recv_offset[recvTask]],
                   Recv_count[recvTask] * sizeof(struct graddata_in), MPI_BYTE,
                   recvTask, TAG_GRAD_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      tend = second();
      timecommsumm1 += timediff(tstart, tend);

      myfree(GradDataIn);
      GradDataResult =
       (struct graddata_out *) mymalloc("GradDataResult", Nimport * sizeof(struct graddata_out));
      GradDataOut =
        (struct graddata_out *) mymalloc("GradDataOut", Nexport * sizeof(struct graddata_out));


      // now do the particles that were sent to us
      tstart = second();
      NextJ = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
      gradient_evaluate_secondary(&ThisThread);

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
                  MPI_Sendrecv(&GradDataResult[Recv_offset[recvTask]],
                               Recv_count[recvTask] * sizeof(struct graddata_out),
                               MPI_BYTE, recvTask, TAG_GRAD_B,
                               &GradDataOut[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct graddata_out),
                               MPI_BYTE, recvTask, TAG_GRAD_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      tend = second();
      timecommsumm2 += timediff(tstart, tend);

      // add the result to the local particles
      tstart = second();
      for(int j = 0; j < Nexport; j++)
        {
          int place = DataIndexTable[j].Index;
          out2particle_gradient(&GradDataOut[j], place, 1);
        }
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(GradDataOut);
      myfree(GradDataResult);
      myfree(GradDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

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

#ifdef WINDS
      if (SphP[i].DelayTime > 0)
	{
	  MyAtLeastDouble windspeed, hsml_c;
#if !defined(LT_WIND_VELOCITY) && !defined(LT_STELLAREVOLUTION)
	  windspeed = sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			   All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency) * All.Time;
#else
#ifdef LT_WIND_VELOCITY
	  windspeed = LT_WIND_VELOCITY * All.Time;
#else /* !GM_MUPPI */
	  int IMFi, SFi;
	  get_SF_index(i, &SFi, &IMFi);
	  windspeed = sqrt(2 * SFs[SFi].WindEnergyFraction * SFs[SFi].totFactorSN *
			   SFs[SFi].EgySpecSN / (1 - SFs[SFi].totFactorSN) / SFs[SFi].WindEfficiency) *
	    All.Time;
#endif /* ends WIND VELOCITY */
#endif /* ends !WIND_VELOCITY && !LT_STELLAREVOLUTION */
	  MyAtLeastDouble fac_mu;
	  if(All.ComovingIntegrationOn)
	    fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
	  else
	    fac_mu = 1.0;
	  windspeed *= fac_mu;
#ifdef QUICK_LYALPHA
	  hsml_c = pow(All.WindFreeTravelDensFac * All.OverDensThresh / SphP[i].Density, (1. / 3.));
#else
	  MyAtLeastDouble myPhysDensThresh;
#ifdef LT_STELLAREVOLUTION
	  myPhysDensThresh = get_PhysDensThresh(i);
#elif defined(EB_SFR_MAGNETIC)
	  myPhysDensThresh = SphP[i].DensityThreshold;
#else
	  myPhysDensThresh = All.PhysDensThresh;
#endif
	  hsml_c = pow(All.WindFreeTravelDensFac * myPhysDensThresh /
		       (SphP[i].Density * All.cf_a3inv), (1. / 3.));
#endif
	  SphP[i].MaxSignalVel = hsml_c * DMAX((2 * windspeed), SphP[i].MaxSignalVel);
	}
#endif
      
      if (gradient_isactive(i))
        {
          SquareMatrix<NUMDIMS> Egrad(SphP[i].Bgrad);
          MyFloat sqdConditionNumber;
	  MyFloat det_tmp = SphP[i].Bgrad.GetDeterminant();
	  if (std::isnormal(det_tmp) and fabs(det_tmp) > 1e-30)
	    {
	      SphP[i].Bgrad.InvertMatrix();
	  
	      // Check the accuracy of the integral gradients (using the square of the condition number)
	      MyFloat modE = MyFloat(0.0);
	      MyFloat modB = MyFloat(0.0);
	      for (int k=0; k<NUMDIMS; k++) {
		for (int kk=0; kk<NUMDIMS; kk++) {
		  modE += Egrad(k,kk)*Egrad(k,kk);
		  modB += SphP[i].Bgrad(k,kk)*SphP[i].Bgrad(k,kk);
		}
	      }
	      sqdConditionNumber = modE*modB / (MyFloat) (NUMDIMS*NUMDIMS);
	    }
	  else
	    {
	      // ensure the matrix is invertable
	      MyFloat correction = 0;
	      for (int k = 0; k<NUMDIMS; k++) correction += SphP[i].Bgrad(k,k);
	      correction = ((correction>=0)-(correction<0)) /*as sign, but returning 1 for 0*/ * std::max(1.05 * abs(correction) / NUMDIMS / 1.0e4, 1e-15);

	      while(true)
		{
		  for (int k = 0; k<NUMDIMS; k++) SphP[i].Bgrad(k,k) += correction;
		  if (std::isnormal(SphP[i].Bgrad.GetDeterminant()) and fabs(SphP[i].Bgrad.GetDeterminant()) >= 1.0e-30) break;
		  correction *= 1.2;
		}
	      SphP[i].Bgrad.InvertMatrix();

	      sqdConditionNumber = 1e30; // some high number
	    }
	  SphP[i].ConditionNumber = sqrt(sqdConditionNumber);
	  
          // If integral gradients are ok, then continue to use to compute gradients
          if (sqdConditionNumber < 1.0e4) {
            for (int var=0; var<NUMDIMS+2; var++) {
              SphP[i].gradW[var] = SphP[i].Bgrad*SphP[i].gradW[var];
            }
          }
          // Otherwise, use less accurate but more stable SPH-like gradients
          else {
            for (int var=0; var<NUMDIMS+2; var++) { 
              SphP[i].gradW[var] = GradStoreOut[i].gradWSph[var];
            }
          }
	  for (int k = 0; k < NUMDIMS+2; k++) { // avoid nan values
	    for (int kk = 0; kk < NUMDIMS; kk++) {
	      if (!std::isnormal(SphP[i].gradW[k][kk])) SphP[i].gradW[k][kk] = (MyFloat)0.0;
	    }
	  }

	  SphP[i].Vel_smooth /= GradStoreOut[i].neighbour_mass;

          // Update SphP.HydroAccel (used in timestepping)
#if NUMDIMS == 1
          Vector<NUMDIMS> WVel(SphP[i].Wprim[0]);
#elif NUMDIMS == 2
          Vector<NUMDIMS> WVel(SphP[i].Wprim[0], SphP[i].Wprim[1]);
#elif NUMDIMS == 3
          Vector<NUMDIMS> WVel(SphP[i].Wprim[0], SphP[i].Wprim[1], SphP[i].Wprim[2]);
#endif
	  const MyFloat u = SphP[i].InternalEnergyPred;
          //const MyFloat u = eos->GetUFromEntropy(SphP[i].Density, SphP[i].Entropy, 1);
          WFluidVector<NUMDIMS> Wcurrent(WVel, SphP[i].Wprim[NUMDIMS], SphP[i].Wprim[NUMDIMS+1]);
          WFluidVector<NUMDIMS> Wdot = WFluidVector<NUMDIMS>::ComputePrimitiveTimeDerivative
            (pow(eos->SoundSpeed(P[i].Mass*SphP[i].NumDens, u, All.cf_a3inv), 2), Wcurrent, SphP[i].gradW, SphP[i].Vel_smooth);

          for (int k=0; k<NUMDIMS; k++) {
            SphP[i].HydroAccel[k] = Wdot[k];
          }
	  SphP[i].dQ = (MyFloat)0.0; // zero out before flux calculation
	  SphP[i].dQdt = (MyFloat)0.0;
	}
    }

  tend = second();
  timecomp1 += timediff(tstart, tend);

  myfree(GradStoreOut);
  
  // collect some timing information
  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_MFMGRADCOMPUTE] += timecomp;
  CPU_Step[CPU_MFMGRADWAIT] += timewait;
  CPU_Step[CPU_MFMGRADCOMM] += timecomm;
  CPU_Step[CPU_MFMGRADMISC] += timeall - (timecomp + timewait + timecomm);

}



// This function represents the core of the gradient computation. The
// target particle may either be local, or reside in the communication buffer.
int gradient_evaluate(int target, int mode, int *exportflag, int *exportnodecount,
                      int *exportindex, int *ngblist)
{
  int listindex = 0;
  int startnode;
  double h, h2, hinv, hinv3, hinv4;

  struct kernel_density kernel;
  struct graddata_in local;
  struct graddata_out out;
  memset(&out, 0, sizeof(struct graddata_out));

  if(mode == 0) {
    particle2in_gradient(&local, target);
  }
  else {
    local = GradDataGet[target];
  }

  h2 = local.Hsml * local.Hsml;
  kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

  if(mode == 0)
    {
      startnode = All.MaxPart;    // root node
    }
  else
    {
      startnode = GradDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;    // open it
    }


  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          int numngb_inbox = ngb_treefind_variable_threads(&(local.Pos[0]), local.Hsml, target,
                                                       &startnode, mode, exportflag,
                                                       exportnodecount, exportindex, ngblist);

          if(numngb_inbox < 0)
            return -1;

          for(int n = 0; n < numngb_inbox; n++)
            {
              int j = ngblist[n];

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
              double r2 = dr.GetMagnitudeSquared();

              if(r2 < h2)
                {
                  kernel.r = sqrt(r2);
                  double u = kernel.r * kernel.hinv;

                  kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                  double wk = kernel.wk*kernel.hinv*local.Volume;
                  double dwk = kernel.dwk*local.Volume;

                  for (int k = 0; k < NUMDIMS; k++) {
                    for (int kk = 0; kk < NUMDIMS; kk++) {
                      out.Ematrix(k,kk) += dr[k]*dr[kk]*wk;
                    }
                  }

                  for (int var = 0; var < NUMDIMS + 2; var++) {
                    for (int k = 0; k < NUMDIMS; k++) {
                      out.gradW[var][k] += dr[k]*(SphP[j].Wprim[var] - local.Wprim[var])*wk;
                    }
                  }

                  // Special slope limiter terms
#if MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009 || MFM_SLOPE_LIMITER == LIMITER_SCALAR || MFM_SLOPE_LIMITER == LIMITER_GIZMO
                  for (int k = 0; k < NUMDIMS + 2; k++) {
                    out.Wmin[k] = std::min(out.Wmin[k], SphP[j].Wprim[k]);
                    out.Wmax[k] = std::max(out.Wmax[k], SphP[j].Wprim[k]);
                  }
#endif
#if MFM_SLOPE_LIMITER == LIMITER_SCALAR || MFM_SLOPE_LIMITER == LIMITER_GIZMO
		  out.DistNgbSqdMax = std::max(out.DistNgbSqdMax, MyFloat(r2));
#endif
		  for (int k = 0; k < NUMDIMS; k++) {
		    out.Vel_smooth[k] += P[j].Mass*SphP[j].Wprim[k];
		  }
		  out.neighbour_mass+=P[j].Mass;
		  
                  if(kernel.r > 0)
                    {
                      const double invr = 1.0 / kernel.r;
                      for (int var = 0; var < NUMDIMS + 2; var++) {
                        for (int k = 0; k < NUMDIMS; k++) {
                          out.gradWSph[var][k] -= dr[k] * (SphP[j].Wprim[var] - local.Wprim[var]) * dwk * invr;
                        }
                      }

                      Vector<NUMDIMS> dv;
                      for (int k = 0; k < NUMDIMS; k++) dv[k] = SphP[j].Wprim[k] - local.Wprim[k];
                      const MyFloat dvdr = Vector<NUMDIMS>::DotProduct(dv, dr);

		      out.MaxVelSquareNgb = std::max(out.MaxVelSquareNgb, dv.GetMagnitudeSquared());
		      
                      // Need to compute signal speed properly here
		      const MyFloat u = SphP[j].InternalEnergyPred;
		      //const MyFloat u = eos->GetUFromEntropy(SphP[j].NumDens*P[j].Mass, SphP[j].Entropy, All.cf_a3inv);
                      double vsigij = local.SoundSpeed + eos->SoundSpeed(P[j].Mass*SphP[j].NumDens, u, All.cf_a3inv) - std::min(MyAtLeastDouble(0.0), dvdr / kernel.r);
                      out.MaxSignalVel = std::max(double(out.MaxSignalVel), vsigij);

                      if (WAKEUP >0){
                        if(vsigij > WAKEUP * SphP[j].MaxSignalVel){
                          if(SphP[j].wakeup > P[target].TimeBin)

                            #if defined(_OPENMP) && !defined(ACC_INCLUDE_HYDRA)
                              #pragma omp critical(__wakeup__)
                            #endif
                          {
                            if(SphP[j].wakeup > P[target].TimeBin)
                        	  SphP[j].wakeup = P[target].TimeBin;
                        	}
                        }
                      }

                    }
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = GradDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }

    }


  if(mode == 0)
    {
      out2particle_gradient(&out, target, 0);
    }
  else
    {
      GradDataResult[target] = out;
    }

  return 0;
}


void *gradient_evaluate_primary(void *p)
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
            while(!gradient_isactive(i));
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

      if(gradient_isactive(i))
        {
          if(gradient_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
            break;		/* export buffer has filled up */
        }

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *gradient_evaluate_secondary(void *p)
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

      gradient_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}
