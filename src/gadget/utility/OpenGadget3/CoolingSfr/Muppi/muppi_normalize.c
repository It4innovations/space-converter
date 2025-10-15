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


#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"
#include "../../Hydro/kernel.h"
#include "../../System/tags.h"

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>


/*! \file muppi_normalize.c
*  \brief Computation normalizations for MUPPI
*
*/


struct kernel_muppinorm
{
  MyAtLeastDouble dx, dy, dz, r;
  MyAtLeastDouble wk_i, wk_j, dwk_i, dwk_j;
  MyAtLeastDouble h_i, h_j;
};

struct muppinormdata_in
{
  MyLongDouble Pos[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyLongDouble Density;

#if defined(PARTICLE_DEBUG) 
  MyIDType ID;			/*!< particle identifier */
#endif

  MyAtLeastDouble E_out;
  MyAtLeastDouble GradDens[3];
  MyAtLeastDouble Eout_norm;
  MyAtLeastDouble NnpD;
  MyIDType Nnp;
#ifdef MV_GM_STELLAR_KIN_FB2
  MyAtLeastDouble Ekin_norm_fountain;
#endif
#ifdef GM_COUNT_PARTICLES_IN_CONE
  MyAtLeastDouble Ekin_total;
#endif


#ifndef DONOTUSENODELIST
  int NodeList[NODELISTLENGTH];
#endif
} *MuppiNormDataIn, *MuppiNormDataGet;


struct muppinormdata_out
{
  MyAtLeastDouble Eout_norm;
#ifdef MV_GM_STELLAR_KIN_FB2
  MyAtLeastDouble Ekin_norm_fountain;
#ifdef GM_COUNT_PARTICLES_IN_CONE
  MyAtLeastDouble Ekin_total;
  int nnorm; 
#endif
#endif

  MyAtLeastDouble NnpD;
  MyIDType Nnp;
} *MuppiNormDataResult, *MuppiNormDataOut;

static inline void particle2in_muppinorm(struct muppinormdata_in *in, int i);
static inline void out2particle_muppinorm(struct muppinormdata_out *out, int i, int mode);



static inline void particle2in_muppinorm(struct muppinormdata_in *in, int i)
{
  int k;

  for(k = 0; k < 3; k++)
    {
      in->Pos[k] = P[i].Pos[k];
      in->GradDens[k] = SphP[i].GradDens[k];
    }
  in->Hsml = P[i].Hsml;

  in->Mass = P[i].Mass - SphP[i].M_sf;
  in->NnpD = SphP[i].NnpD;
#ifdef GM_COUNT_PARTICLES_IN_CONE
  in->E_out = SphP[i].E_out;
#endif

  in->Density = SphP[i].Density;


#if defined(PARTICLE_DEBUG) 
  in->ID = P[i].ID;
#endif

}

static inline void out2particle_muppinorm(struct muppinormdata_out *out, int i, int mode)
{

  ASSIGN_ADD(SphP[i].Eout_norm, out->Eout_norm, mode);


#ifdef MV_GM_STELLAR_KIN_FB2
  ASSIGN_ADD(SphP[i].Ekin_norm_fountain, out->Ekin_norm_fountain, mode);

#ifdef GM_COUNT_PARTICLES_IN_CONE
  ASSIGN_ADD(SphP[i].Ekin_total, out->Ekin_total, mode);
#endif

#endif


#ifdef GM_COUNT_PARTICLES_IN_CONE
  ASSIGN_ADD(SphP[i].nnorm, out->nnorm, mode);
#endif
  if(out->Nnp != -1 && out->NnpD > SphP[i].NnpD)
    {
      SphP[i].Nnp = out->Nnp;
      SphP[i].NnpD = out->NnpD;
    }

}



void muppinorm(void)
{
  int i, j, k, ngrp, ndone, ndone_flag;
  int recvTask, place;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0, timenetwork = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;

  int save_NextParticle;

  long long n_exported = 0;


  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      {
#if defined(GM_COUNT_PARTICLES_IN_CONE)
	SphP[i].nnorm = 0;
#endif

	SphP[i].Nnp = -1;
	SphP[i].NnpD = -100;
      }

  /* allocate buffers to arrange communication */


  int NTaskTimesNumPart;

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct muppinormdata_in) +
					     sizeof(struct muppinormdata_out) +
					     sizemax(sizeof(struct muppinormdata_in),
						     sizeof(struct muppinormdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  CPU_Step[CPU_HYDMISC] += measure_time();
  t0 = second();

  NextParticle = FirstActiveParticle;	/* beginn with this index */

  do
    {

      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();

      {
	int mainthreadid = 0;
	muppinorm_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
      }

      tend = second();
      timecomp1 += timediff(tstart, tend);

      if(BufferFullFlag)
	{
	  int last_nextparticle = NextParticle;

	  NextParticle = save_NextParticle;

	  while(NextParticle >= 0)
	    {
	      if(NextParticle == last_nextparticle)
		break;

	      if(ProcessedFlag[NextParticle] != 1)
		break;

	      ProcessedFlag[NextParticle] = 2;

	      NextParticle = NextActiveParticle[NextParticle];
	    }

	  if(NextParticle == save_NextParticle)
	    {
	      /* in this case, the buffer is too small to process even a single particle */
	      endrun(12998);
	    }


	  int new_export = 0;

	  for(j = 0, k = 0; j < Nexport; j++)
	    if(ProcessedFlag[DataIndexTable[j].Index] != 2)
	      {
		if(k < j + 1)
		  k = j + 1;

		for(; k < Nexport; k++)
		  if(ProcessedFlag[DataIndexTable[k].Index] == 2)
		    {
		      int old_index = DataIndexTable[j].Index;

		      DataIndexTable[j] = DataIndexTable[k];
		      DataNodeList[j] = DataNodeList[k];
		      DataIndexTable[j].IndexGet = j;
		      new_export++;

		      DataIndexTable[k].Index = old_index;
		      k++;
		      break;
		    }
	      }
	    else
	      new_export++;

	  Nexport = new_export;

	}


      n_exported += Nexport;

      for(j = 0; j < NTask; j++)
	Send_count[j] = 0;
      for(j = 0; j < Nexport; j++)
	Send_count[DataIndexTable[j].Task]++;

      MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

      tstart = second();

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);

      tend = second();
      timewait1 += timediff(tstart, tend);

      for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      MuppiNormDataGet = (struct muppinormdata_in *) mymalloc("MuppiNormDataGet", Nimport * sizeof(struct muppinormdata_in));
      MuppiNormDataIn = (struct muppinormdata_in *) mymalloc("MuppiNormDataIn", Nexport * sizeof(struct muppinormdata_in));

      /* prepare particle data for export */

      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  particle2in_muppinorm(&MuppiNormDataIn[j], place);
#ifndef DONOTUSENODELIST
	  memcpy(MuppiNormDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif

	}

      /* exchange particle data */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&MuppiNormDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct muppinormdata_in), MPI_BYTE,
			       recvTask, TAG_MPPIN_A,
			       &MuppiNormDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct muppinormdata_in), MPI_BYTE,
			       recvTask, TAG_MPPIN_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);


      myfree(MuppiNormDataIn);
      MuppiNormDataResult =
	(struct muppinormdata_out *) mymalloc("MuppiNormDataResult", Nimport * sizeof(struct muppinormdata_out));
      MuppiNormDataOut =
	(struct muppinormdata_out *) mymalloc("MuppiNormDataOut", Nexport * sizeof(struct muppinormdata_out));


      report_memory_usage(&HighMark_mppin, "MUPPI_NORMALIZATION");

      /* now do the particles that were sent to us */

      tstart = second();

      NextJ = 0;

      {
	int mainthreadid = 0;
	muppinorm_evaluate_secondary(&mainthreadid);
      }

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


      /* get the result */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&MuppiNormDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct muppinormdata_out),
			       MPI_BYTE, recvTask, TAG_MPPIN_B,
			       &MuppiNormDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct muppinormdata_out),
			       MPI_BYTE, recvTask, TAG_MPPIN_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm2 += timediff(tstart, tend);



      /* add the result to the local particles */
      tstart = second();
      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  out2particle_muppinorm(&MuppiNormDataOut[j], place, 1);
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(MuppiNormDataOut);
      myfree(MuppiNormDataResult);
      myfree(MuppiNormDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


  /* do final operations on results */

#ifdef GM_COUNT_PARTICLES_IN_CONE
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	if(SphP[i].Ekin_norm_fountain>0.0)
	  SphP[i].Ekin_total /= SphP[i].Ekin_norm_fountain;
	else
	  SphP[i].Ekin_total = 0.0;
      }
#endif



  /* collect some timing information */
  
  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_HYDCOMPUTE] += timecomp;
  CPU_Step[CPU_HYDWAIT] += timewait;
  CPU_Step[CPU_HYDCOMM] += timecomm;
  CPU_Step[CPU_HYDNETWORK] += timenetwork;
  CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm + timenetwork);
}




/*! This function is the 'core' of the SPH force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*/
int muppinorm_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		   int *ngblist)
{
  int startnode, numngb, listindex = 0;
  int j, k, n, l;

  /* avoid unused variable compiler warnings */
  (void) l;
  (void) k;

  double hinv, hinv3, hinv4, r2, u;

  struct kernel_muppinorm kernel;
  struct muppinormdata_in local;
  struct muppinormdata_out out;
  memset(&out, 0, sizeof(struct muppinormdata_out)); /* I am strongly critical on this */

  if(mode == 0)
    particle2in_muppinorm(&local, target);
  else
    local = MuppiNormDataGet[target];


  kernel.h_i = local.Hsml;

  MyAtLeastDouble angle;
  MyAtLeastDouble rn, r2n, wk_n, dwk_n, un, h_i2;
  MyAtLeastDouble pto_rx, pto_ry, pto_rz, dir_norm, r_a, r2_a;
  MyAtLeastDouble grad_tot, theta, phi, theta_1, theta_2, phi_1, phi_2;
  MyAtLeastDouble alpha, beta, alphap, alpham, betap, betam;


  local.Eout_norm = 0.0;

#ifdef MV_GM_STELLAR_KIN_FB2
  local.Ekin_norm_fountain = 0.0;
#endif

  
#ifdef GM_COUNT_PARTICLES_IN_CONE
  out.nnorm = 0;
#endif
  grad_tot = sqrt(local.GradDens[0] * local.GradDens[0] +
		  local.GradDens[1] * local.GradDens[1] + local.GradDens[2] * local.GradDens[2]);

  theta = atan2(-local.GradDens[1], -local.GradDens[0]);
  phi = acos(-local.GradDens[2] / grad_tot);

  theta_1 = theta - ConeAngle_FB / 2;
  theta_2 = theta + ConeAngle_FB / 2;

  phi_1 = phi - ConeAngle_FB / 2;
  phi_2 = phi + ConeAngle_FB / 2;


  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
#ifndef DONOTUSENODELIST
      startnode = MuppiNormDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
#else
      startnode = All.MaxPart;	/* root node */
#endif
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb =
	    ngb_treefind_pairs_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
				       exportnodecount, exportindex, ngblist);

	  if(numngb < 0)
	    return -1;


	  for(n = 0; n < numngb; n++)
	    {
	      j = ngblist[n];



#ifdef NOWINDTIMESTEPPING
#ifdef WINDS
	      if(P[j].Type == 0)
		if(SphP[j].DelayTime > 0)	/* ignore the wind particles */
		  continue;
#endif
#endif

	      kernel.dx = local.Pos[0] - P[j].Pos[0];
	      kernel.dy = local.Pos[1] - P[j].Pos[1];
	      kernel.dz = local.Pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      kernel.dx = NEAREST_X(kernel.dx);
	      kernel.dy = NEAREST_Y(kernel.dy);
	      kernel.dz = NEAREST_Z(kernel.dz);
#endif

	      r2 = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;
	      kernel.h_j = P[j].Hsml;

	      /* here I evaluate the normalization for MP energy redistribution.
	         I do this here for all MP particles to avoid a communication
	         round later. */

	      pto_rx =
		-local.GradDens[2] * (P[j].Pos[1] - local.Pos[1]) + local.GradDens[1] * (P[j].Pos[2] -
											 local.Pos[2]);
	      pto_ry =
		-local.GradDens[0] * (P[j].Pos[2] - local.Pos[2]) + local.GradDens[2] * (P[j].Pos[0] -
											 local.Pos[0]);
	      pto_rz =
		-local.GradDens[1] * (P[j].Pos[0] - local.Pos[0]) + local.GradDens[0] * (P[j].Pos[1] -
											 local.Pos[1]);
	      dir_norm =
		local.GradDens[0] * local.GradDens[0] + local.GradDens[1] * local.GradDens[1] +
		local.GradDens[2] * local.GradDens[2];

	      r2n = (pto_rx * pto_rx + pto_ry * pto_ry + pto_rz * pto_rz) / dir_norm;
	      rn = sqrt(r2n);

	      r2_a = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;
	      r_a = sqrt(r2_a);

	      alpha = atan2(-kernel.dy, -kernel.dx);
	      beta = acos(-kernel.dz / r_a);
	      alphap = alpha + 2. * M_PI;
	      alpham = alpha - 2. * M_PI;
	      betap = beta + M_PI;
	      betam = beta - M_PI;

	      
	      if(r2 < kernel.h_i * kernel.h_i || r2 < kernel.h_j * kernel.h_j)
		{   
		  kernel.r = sqrt(r2);
		  if(kernel.r > 0)
		    {

		      if(kernel.r < kernel.h_i)
			{
			  kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
			  u = kernel.r * hinv;
			  kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, 1);
			}
		      else
			{
			  kernel.dwk_i = 0;
			  kernel.wk_i = 0;
			}

		      if(kernel.r < kernel.h_j)
			{
			  kernel_hinv(kernel.h_j, &hinv, &hinv3, &hinv4);
			  u = kernel.r * hinv;
			  kernel_main(u, hinv3, hinv4, &kernel.wk_j, &kernel.dwk_j, 1);
			}
		      else
			{
			  kernel.dwk_j = 0;
			  kernel.wk_j = 0;
			}

		      if(rn < kernel.h_i)
			{
			  kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
			  un = rn * hinv;
			  kernel_main(un, hinv3, hinv4, &wk_n, &dwk_n, -1);
			}
		      else
			wk_n = 0;

		      angle = (local.GradDens[0] * kernel.dx +
			       local.GradDens[1] * kernel.dy +
			       local.GradDens[2] * kernel.dz) / grad_tot / sqrt(r2);

		      h_i2 = kernel.h_i * kernel.h_i;
		      if(r2 < h_i2 && angle > local.NnpD)    
			{
			  local.NnpD = angle;
			  local.Nnp = P[j].ID;
			  out.NnpD = angle;
			  out.Nnp = P[j].ID;
			}


		      if(r2n < h_i2 && r2 < h_i2 && rn > 0.00001 &&
			 ((theta_1 < alpha && alpha < theta_2) || (theta_1 < alphap && alphap < theta_2)
			  || (theta_1 < alpham && alpham < theta_2)) && ((phi_1 < beta && beta < phi_2)
			  || (phi_1 < betap && betap < phi_2)
			  || (phi_1 < betam && betam < phi_2)))
			{
	
			  //note, added r2 check because here _tree_pairs is used
			  //note2, Eout_norm evaluated for ALL SPH PARTS

			  out.Eout_norm += (P[j].Mass - SphP[j].M_sf) * wk_n / SphP[j].Density;
			} 

#if defined(MV_GM_STELLAR_KIN_FB2)

		      if((kernel.r < kernel.h_i) && (kernel.r > 0.0000001))  
			{
			  if(SphP[j].DelayTime > 0.0)
			    {
			      kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4); 
			      u = kernel.r * hinv;

			      kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, -1);
			      out.Ekin_norm_fountain += (P[j].Mass - SphP[j].M_sf) * kernel.wk_i / SphP[j].Density;         /* milena, normalizzazione */
			    }
			  
		
#ifdef GM_COUNT_PARTICLES_IN_CONE
		      if(SphP[j].DelayTime > 0.0)
			out.Ekin_total +=  (P[j].Mass - SphP[j].M_sf)/ SphP[j].Density * local.E_out * kernel.wk_i * All.FracEgyKin;
#endif
#endif

#ifdef GM_COUNT_PARTICLES_IN_CONE
#if defined(COUNT_FOUNTAIN_PARTICLES_IN_CONE) && defined(MV_GM_STELLAR_KIN_FB2)
		      if(SphP[j].DelayTime > 0.0)
#endif
			out.nnorm++;
#endif
			}  


		    } //if(kernel.r>0..
		} //if(r2<...  
	    } //for(n=0, n<numnbg
	} //while(startnode..

#ifndef DONOTUSENODELIST
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = MuppiNormDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
#endif
    }


  /* Now collect the result at the right place */
  if(mode == 0)
    out2particle_muppinorm(&out, target, 0);
  else
    MuppiNormDataResult[target] = out;

  return 0;
}

void *muppinorm_evaluate_primary(void *p)
{
  int thread_id = *(int *) p;

  int i, j;

  int *exportflag, *exportnodecount, *exportindex, *ngblist;


  ngblist = Ngblist + thread_id * NumPart;
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
      {
	if(BufferFullFlag != 0 || NextParticle < 0)
	  {
	    exitFlag = 1;
	  }
	else
	  {
	    i = NextParticle;
	    ProcessedFlag[i] = 0;
	    NextParticle = NextActiveParticle[NextParticle];
	  }
      }
      if(exitFlag)
	break;

#ifdef BLACK_HOLES
      if(P[i].Type == 0 && P[i].Mass > 0)
#else
      if(P[i].Type == 0)
#endif
	{
	  if(muppinorm_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *muppinorm_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;

  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;


  while(1)
    {
      {
	j = NextJ;
	NextJ++;
      }

      if(j >= Nimport)
	break;

      muppinorm_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}

