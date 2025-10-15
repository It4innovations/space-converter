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


/*! Structure for communication during the slope limiter computation. Holds data that is sent to other processors.
 */
static struct slopedata_in
{
  MyIDType ID;
  MyFloat Hsml;
  MyLongDouble Pos[3];
  WFluidVector<NUMDIMS> Wprim;
  Vector<NUMDIMS> gradW[NUMDIMS+2];
  int NodeList[NODELISTLENGTH];
#if MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009
  WFluidVector<NUMDIMS> dWmin;
  WFluidVector<NUMDIMS> dWmax;
#endif
}
*SlopeDataIn, *SlopeDataGet;


/*! Structure for communication during the slope limiter computation. Returns partial results from node.
 */
static struct slopedata_out
{
  MyFloat AlphaSlope[NUMDIMS+2];
}
 *SlopeDataResult, *SlopeDataOut;


/*! Create reduced particle data structure for communicating to other nodes for gradient computation */
void particle2in_slope(struct slopedata_in *in, int i)
{
  in->ID = P[i].ID;
  in->Hsml = P[i].Hsml;
  for(int k = 0; k<3; k++) {
    in->Pos[k] = P[i].Pos[k];
  }
  for(int k = 0; k < NUMDIMS+2; k++) {
    in->Wprim[k] = SphP[i].Wprim[k];
    in->gradW[k] = SphP[i].gradW[k];
#if MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009
    in->dWmax[k] = SphP[i].Wmax[k] - SphP[i].Wprim[k];
    in->dWmin[k] = SphP[i].Wmin[k] - SphP[i].Wprim[k];
#endif
  }
}


/*! Create data structure with partial sums/results to send back to particle host node */
void out2particle_slope(struct slopedata_out *out, int i, int mode)
{
  if(P[i].Type == 0) {
#if MFM_SLOPE_LIMITER == LIMITER_TVD_SCALAR || MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009
    for (int k = 0; k < NUMDIMS + 2; k++) {
      ASSIGN_MIN(SphP[i].AlphaSlope[k], out->AlphaSlope[k], mode);
    }
#endif
  }
}


// Returns 1 if particle is active (i.e. beginning of timestep).  Otherwise returns 0.
int slope_isactive(int n)
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


/*! \file mfm_limiter.cpp
 *  \brief Compute slope limiter terms for all active particles on this step
 */


// This function computes the gradients for all active particles
void mfm_limiters(void)
{
  int ndone, ndone_flag, dummy;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  long long n_exported = 0;

  // Display information message that this step is executed on Task 0
  VERBOSE(3, "Calculating slope limiter for primitive fluid vector (Wprim)\n");
  
#if MFM_SLOPE_LIMITER == LIMITER_TVD_SCALAR || MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009 // skip the additional neighbor loop for all limiters, where it is not necessary.
  int NTaskTimesNumPart = ((int) (maxThreads)) * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  size_t MyBufferSize = KD_BUFFER_MANAGEMENT * FreeBytes / (1024.0 * 1024.0);
  
  All.BunchSize =
    (int) ((MyBufferSize  * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct slopedata_in) + sizeof(struct slopedata_out) +
                                             sizemax(sizeof(struct slopedata_in), sizeof(struct slopedata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  VERBOSE(2, "MFM limiters: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
	  FreeBytes / (1024.0 * 1024.0));
#endif // MFM_SLOPE_LIMITER == LIMITER_TVD_SCALAR || MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009
  
  CPU_Step[CPU_MFMLIMMISC] += measure_time();
  t0 = second();

#if MFM_SLOPE_LIMITER == LIMITER_TVD_SCALAR || MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009 // skip the additional neighbor loop for all limiters, where it is not necessary.
  // begin with this index
  NextParticle = FirstActiveParticle;

  do
    {
      BufferFullFlag = 0;
      Nexport = 0;
      int save_NextParticle = NextParticle;

      tstart = second();

#ifdef _OPENMP
#pragma omp parallel
#endif
      slope_evaluate_primary(&ThisThread);  // do local particles and prepare export list

      tend = second();
      timecomp1 += timediff(tstart, tend);

      /* if we are not grouping particles and SendRecv fails, then we panic as there is no memory to process even a single particles */
      int panic_on_error = !group_particles;
      int is_send_receive_ok = prepare_sendrec_offset_send(save_NextParticle, panic_on_error); //1 means panic on bufer too small
      if(!is_send_receive_ok) {
	group_particles = 0;
	continue;
      }
      
      SlopeDataGet =
        (struct slopedata_in *) mymalloc("SlopeDataGet", Nimport * sizeof(struct slopedata_in));
      SlopeDataIn =
        (struct slopedata_in *) mymalloc("SlopeDataIn", Nexport * sizeof(struct slopedata_in));

      /* prepare particle data for export */
      for(int j = 0; j < Nexport; j++)
        {
          int place = DataIndexTable[j].Index;

          particle2in_slope(&SlopeDataIn[j], place);

          memcpy(SlopeDataIn[j].NodeList,
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
                  MPI_Sendrecv(&SlopeDataIn[Send_offset[recvTask]],
                   Send_count[recvTask] * sizeof(struct slopedata_in), MPI_BYTE,
                   recvTask, TAG_SLOPE_A,
                   &SlopeDataGet[Recv_offset[recvTask]],
                   Recv_count[recvTask] * sizeof(struct slopedata_in), MPI_BYTE,
                   recvTask, TAG_SLOPE_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      tend = second();
      timecommsumm1 += timediff(tstart, tend);

      myfree(SlopeDataIn);
      SlopeDataResult =
       (struct slopedata_out *) mymalloc("SlopeDataResult", Nimport * sizeof(struct slopedata_out));
      SlopeDataOut =
        (struct slopedata_out *) mymalloc("SlopeDataOut", Nexport * sizeof(struct slopedata_out));


      // now do the particles that were sent to us
      tstart = second();
      NextJ = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
      slope_evaluate_secondary(&ThisThread);

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
                  MPI_Sendrecv(&SlopeDataResult[Recv_offset[recvTask]],
                               Recv_count[recvTask] * sizeof(struct slopedata_out),
                               MPI_BYTE, recvTask, TAG_SLOPE_B,
                               &SlopeDataOut[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct slopedata_out),
                               MPI_BYTE, recvTask, TAG_SLOPE_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
          out2particle_slope(&SlopeDataOut[j], place, 1);
        }
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(SlopeDataOut);
      myfree(SlopeDataResult);
      myfree(SlopeDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
#endif // MFM_SLOPE_LIMITER == LIMITER_TVD_SCALAR || MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009

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
      if (slope_isactive(i))
        {
#if MFM_SLOPE_LIMITER == LIMITER_TVD_SCALAR

#elif MFM_SLOPE_LIMITER == LIMITER_SCALAR
          // TODO : Factor of 1.0 refers to extend of M4 kernel: replace with whatever is needed
          const MyFloat drMax = std::max(sqrt(SphP[i].DistNgbSqdMax), MyFloat(1.0)*P[i].Hsml);
          for (int k = 0; k < NUMDIMS + 2; k++) {
            const MyFloat gradW = SphP[i].gradW[k].GetMagnitude();
            const MyFloat dW = drMax*gradW;
            const MyFloat dWmax = SphP[i].Wmax[k] - SphP[i].Wprim[k];
            const MyFloat dWmin= SphP[i].Wprim[k] - SphP[i].Wmin[k];
            SphP[i].AlphaSlope[k] = MyFloat(1.0);
            if (dW != 0.0) {
              SphP[i].AlphaSlope[k] = std::max(MyFloat(0.0), std::min(MyFloat(1.0), std::min(dWmax/dW, dWmin/dW)));
            }
          }
#elif MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009

#elif MFM_SLOPE_LIMITER == LIMITER_NULL
	  for (int k = 0; k < NUMDIMS + 2; k++) {
            SphP[i].AlphaSlope[k] = MyFloat(1.0);
          }
#elif MFM_SLOPE_LIMITER == LIMITER_ZERO_SLOPES
          for (int k = 0; k < NUMDIMS + 2; k++) {
            SphP[i].AlphaSlope[k] = MyFloat(0.0);
          }
#elif MFM_SLOPE_LIMITER == LIMITER_GIZMO
	  MyFloat a_limiter = 0.25;
	  //if(SphP[i].ConditionNumber > 100) a_limiter=std::min(0.5, 0.25 + 0.25*(SphP[i].ConditionNumber-100)/100);
	  MyFloat shoot_tol = 0.0;
#if defined(NOGRAVITY) && !defined(MAGNETIC)
	  shoot_tol = 0.0/*0.1*/;
#endif	  
#if SLOPE_LIMITER_TOLERANCE == 2
	  a_limiter *= 2.0;
	  shoot_tol = 0.125;
#elif SLOPE_LIMITER_TOLERANCE == 0
	  a_limiter *= 0.5;
	  shoot_tol = 0.0;
#endif
	  MyFloat h0 = MyFloat(1.0) / (a_limiter * P[i].Hsml);
	  MyFloat drMax = std::max(sqrt(SphP[i].DistNgbSqdMax), MyFloat(1.0)*P[i].Hsml);
	  for (int k = 0; k < NUMDIMS + 2; k++) {
	    MyFloat d_abs = SphP[i].gradW[k].GetMagnitude();
	    SphP[i].AlphaSlope[k] = MyFloat(1.0);
            if (d_abs > 0) {
	      MyFloat this_shoot_tol = (k == NUMDIMS) ? 0.0 : /*0.0*/shoot_tol;
	      MyFloat cfac = h0 / d_abs;
	      const MyFloat dWmax = SphP[i].Wmax[k] - SphP[i].Wprim[k];
	      const MyFloat dWmin = SphP[i].Wprim[k] - SphP[i].Wmin[k];
	      MyFloat dWabsmax = std::max(dWmax,dWmin);
	      MyFloat dWabsmin = std::min(dWmax,dWmin);
	      MyFloat f_corr_over = std::min(dWabsmin + this_shoot_tol*dWabsmax,dWabsmax);
	      cfac *= f_corr_over;
	      if (k >= NUMDIMS) { //positivity-preserving for pressure and density
	        MyFloat fmin = std::min(SphP[i].Wprim[k], std::max(MyFloat(0), std::max(MyFloat(1e-20)*SphP[i].Wprim[k], std::min(MyFloat(0.5)*(SphP[i].Wprim[k]-dWmin),SphP[i].Wprim[k]-f_corr_over))));
	        cfac = std::min(cfac, (SphP[i].Wprim[k]-fmin)/drMax/d_abs );
	      }
	      SphP[i].AlphaSlope[k] = std::min(MyFloat(1.0), cfac);
	    } 
	  }
#endif
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

  CPU_Step[CPU_MFMLIMCOMPUTE] += timecomp;
  CPU_Step[CPU_MFMLIMWAIT] += timewait;
  CPU_Step[CPU_MFMLIMCOMM] += timecomm;
  CPU_Step[CPU_MFMLIMMISC] += timeall - (timecomp + timewait + timecomm);

}
  


// This function represents the core of the gradeint computation. The
// target particle may either be local, or reside in the communication buffer.
int slope_evaluate(int target, int mode, int *exportflag, int *exportnodecount,
                   int *exportindex, int *ngblist)
{
  int listindex = 0;
  int startnode;

  struct slopedata_in local;
  struct slopedata_out out;
  memset(&out, 0, sizeof(struct slopedata_out));

  for (int k = 0; k < NUMDIMS + 2; k++) out.AlphaSlope[k] = 1.0;

  if(mode == 0) {
    particle2in_slope(&local, target);
  }
  else {
    local = SlopeDataGet[target];
  }

  if(mode == 0)
    {
      startnode = All.MaxPart;    // root node
    }
  else
    {
      startnode = SlopeDataGet[target].NodeList[0];
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

              // Skip self-interaction
              //if (local.ID == SphP[j].ID) continue;


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

#if MFM_SLOPE_LIMITER == LIMITER_TVD_SCALAR
              for (int k = 0; k < NUMDIMS + 2; k++) {
                MyFloat dW = Vector<NUMDIMS>::DotProduct(local.gradW[k], dr);
                MyFloat dWcell = SphP[j].Wprim[k] - local.Wprim[k];
                if (dW != 0.0) {
                  out.AlphaSlope[k] = std::min(out.AlphaSlope[k], std::max((MyFloat)0.0, std::min((MyFloat)1.0, dWcell / dW)));
                }
              }
#elif MFM_SLOPE_LIMITER == LIMITER_SPRINGEL_2009
              for (int k = 0; k < NUMDIMS + 2; k++) {
                MyFloat dW = Vector<NUMDIMS>::DotProduct(local.gradW[k], dr);
                if (dW > 0.0) {
                  out.AlphaSlope[k] = std::min(out.AlphaSlope[k], local.dWmax[k] / dW);
                }
                else if (dW < 0.0) {
                  out.AlphaSlope[k] = std::min(out.AlphaSlope[k], local.dWmin[k] / dW);
                }
              }
#elif MFM_SLOPE_LIMITER == LIMITER_NULL
	      
#elif MFM_SLOPE_LIMITER == LIMITER_ZERO_SLOPES
              for (int k = 0; k < NUMDIMS + 2; k++) {
                out.AlphaSlope[k] = 0.0;
              }
#endif
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = SlopeDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Nodes[startnode].u.d.nextnode;    // open it
            }
          }

    }


  if(mode == 0)
    {
      out2particle_slope(&out, target, 0);
      //out2particle_slope(&out, target, 1);
    }
  else
    {
      SlopeDataResult[target] = out;
    }

  return 0;
}


void *slope_evaluate_primary(void *p)
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
            while(!slope_isactive(i));
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

      if(slope_isactive(i))
        {
          if(slope_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
            break;		/* export buffer has filled up */
        }

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *slope_evaluate_secondary(void *p)
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

      slope_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}
