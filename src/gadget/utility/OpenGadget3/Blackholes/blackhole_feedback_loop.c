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

void *blackhole_evaluate_primary(void *p);
void *blackhole_evaluate_secondary(void *p);

/***********************************************************************************/
/************************** Data Structures for loop (accreation) ******************/
/***********************************************************************************/

struct blackholedata_in
{
  MyLongDouble Pos[3];
  MyFloat BH_Density;
  MyFloat BH_Mdot;
  MyFloat Dt;
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat BH_Mass;
  MyFloat Vel[3];
  MyFloat Csnd;
  MyIDType ID;
  int NodeList[NODELISTLENGTH];
  MyFloat BH_TotalFeedbackEfficiency;
  MyFloat Potential;
#ifdef KD_USE_IPOT_FOR_BH_MERGER
  MyFloat IPotential;
#endif
};

struct blackholedata_out
{
  MyLongDouble BH_MinPotPos[3];
  MyFloat BH_MinPot;
  #if defined AD_DYNFRIC
     MyFloat BHDF_DynFric[3];
     int  BHDF_ngb;
  #endif
#if defined(GM_REPOSITION_ON_STARDENSITY_MAX)
  MyLongDouble BH_StarDensPos[3];
  MyLongDouble BH_StarDens;
#ifdef GM_USE_VELS_ON_REPOSITIONING
  MyLongDouble BH_StarDensVel[3];
#endif
#endif
};

void particle2in_blackhole(struct blackholedata_in *in, int i)
{
  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = P[i].Pos[k];
      in->Vel[k] = P[i].Vel[k];
    }
  in->BH_Density = BPP(i).BH_Density;
  in->BH_Mdot = BPP(i).BH_Mdot;
  in->Dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a;
  in->Hsml = P[i].Hsml;
  in->Mass = P[i].Mass;
  in->BH_Mass = BPP(i).BH_Mass;
  in->Csnd = sqrt(GAMMA * BPP(i).BH_Entropy *
		  pow(BPP(i).BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
  in->ID = P[i].ID;
  in->BH_TotalFeedbackEfficiency = BPP(i).BH_TotalFeedbackEfficiency;

  in->Potential = P[i].Potential;
#ifdef KD_USE_IPOT_FOR_BH_MERGER
  in->IPotential = P[i].IPotential;
#endif
}

void out2particle_blackhole(struct blackholedata_out *out, int i, int mode)
{
  int k;
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	BPP(i).BH_MinPotPos[k] = out->BH_MinPotPos[k];
      BPP(i).BH_MinPot = out->BH_MinPot;
     #ifdef AD_DYNFRIC
	for(k=0; k<3; k++) {
	  ASSIGN_ADD(BPP(i).BHDF_DynFric[k], out->BHDF_DynFric[k],mode); }
	  ASSIGN_ADD(BPP(i).BHDF_ngb, out->BHDF_ngb, mode);
     #endif
#ifdef GM_REPOSITION_ON_STARDENSITY_MAX
      for(k = 0; k < 3; k++)
	{
	  BPP(i).BH_StarDensPos[k] = out->BH_StarDensPos[k];
#ifdef GM_USE_VELS_ON_REPOSITIONING
	  BPP(i).BH_StarDensVel[k] = out->BH_StarDensVel[k];
#endif
	}
      BPP(i).BH_StarDens = out->BH_StarDens;
#endif
    }
  else
    {
      if(BPP(i).BH_MinPot > out->BH_MinPot)
	{
	  for(k = 0; k < 3; k++)
	    BPP(i).BH_MinPotPos[k] = out->BH_MinPotPos[k];
	  BPP(i).BH_MinPot = out->BH_MinPot;
	}
#ifdef GM_REPOSITION_ON_STARDENSITY_MAX
      if(BPP(i).BH_StarDens < out->BH_StarDens)
	{
	  BPP(i).BH_StarDens = out->BH_StarDens;
	  for(k = 0; k < 3; k++)
	    {
	      BPP(i).BH_StarDensPos[k] = out->BH_StarDensPos[k];
#ifdef GM_USE_VELS_ON_REPOSITIONING
	      BPP(i).BH_StarDensVel[k] = out->BH_StarDensVel[k];

#endif

	    }
	}
#endif
     #ifdef AD_DYNFRIC
      for(k=0; k<3; k++) {
        ASSIGN_ADD(BPP(i).BHDF_DynFric[k], out->BHDF_DynFric[k],mode); }
      ASSIGN_ADD(BPP(i).BHDF_ngb, out->BHDF_ngb, mode);
     #endif


    }
}

static struct blackholedata_in *BlackholeDataIn, *BlackholeDataGet;
static struct blackholedata_out *BlackholeDataResult, *BlackholeDataOut;

void blackhole_feedback_loop(void)
{
  int ndone_flag, ndone;
  int save_NextParticle;
  MPI_Status status;

  VERBOSE(0, "BH: Compute feedback and prepare swallowing of gas particles and black holes");

  /* allocate buffers to arrange communication */

  int NTaskTimesNumPart;

  NTaskTimesNumPart = ((int) (maxThreads)) * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  size_t MyBufferSize = KD_BUFFER_MANAGEMENT * FreeBytes / (1024.0 * 1024.0);

  All.BunchSize =
    (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					   sizeof(struct blackholedata_in) +
					   sizeof(struct blackholedata_out) +
					   sizemax(sizeof(struct blackholedata_in),
						   sizeof(struct blackholedata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  VERBOSE(1, "BH: using %g MB for buffering, %g MB remaining for other buffers\n", MyBufferSize / 1.0,
	  FreeBytes / (1024.0 * 1024.0));

  /** Let's first spread the feedback energy, and determine which particles may be swalled by whom */

  NextParticle = FirstActiveParticle;	/* first particle for this task */

  do
    {

      /* do local particles and prepare export list */
      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;


#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int thread_id = omp_get_thread_num();
#else
	int thread_id = 0;
#endif
	blackhole_evaluate_primary(&thread_id);	/* do local particles and prepare export list */
      }

      prepare_sendrec_offset_send(save_NextParticle, 1);	//1 means panic on bufer too small

      BlackholeDataGet =
	(struct blackholedata_in *) mymalloc("BlackholeDataGet", Nimport * sizeof(struct blackholedata_in));
      BlackholeDataIn =
	(struct blackholedata_in *) mymalloc("BlackholeDataIn", Nexport * sizeof(struct blackholedata_in));

      for(int j = 0; j < Nexport; j++)
	{
	  int place = DataIndexTable[j].Index;
	  particle2in_blackhole(&BlackholeDataIn[j], place);
	  memcpy(BlackholeDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}


      /* exchange particle data */
      for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  int recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MYMPI_COMM_WORLD, &status);
		}
	    }
	}

      myfree(BlackholeDataIn);

      BlackholeDataResult =
	(struct blackholedata_out *) mymalloc("BlackholeDataResult",
					      Nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut =
	(struct blackholedata_out *) mymalloc("BlackholeDataOut", Nexport * sizeof(struct blackholedata_out));


      /* now do the particles that were sent to us */
      NextJ = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int thread_id = omp_get_thread_num();
#else
	int thread_id = 0;
#endif
	blackhole_evaluate_secondary(&thread_id);	/* do local particles and prepare export list */
      }

      if(NextParticle < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

      /* get the result */
      for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  int recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MYMPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the particles */
      for(int j = 0; j < Nexport; j++)
	{
	  int place = DataIndexTable[j].Index;
	  out2particle_blackhole(&BlackholeDataOut[j], place, 1);
	}

      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

#ifdef AD_DYNFRIC
#ifdef _OPENMP
  int il;
#pragma omp parallel for if(NActivePart>80)
  for(il = 0; il < NActivePart; il++)
    {
      int i = ActiveParticleList[il];
#else 
  for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#endif  // _OPENMP
      if(P[i].Type == 5 && P[i].Mass >0 && BPP(i).BHDF_ngb>0 )
	{
	  
	  for(int k = 0; k < 3 ; k++)
	    {
	      BPP(i).BHDF_DynFric[k] = 3/2 * All.G * All.G /
		(All.ForceSoftening[5]*All.ForceSoftening[5]*All.ForceSoftening[5]) *
		BPP(i).BHDF_DynFric[k];
	      P[i].GravAccel[k] += BPP(i).BHDF_DynFric[k];
	    }

	  fprintf(FdBlackHolesDetails,
		  "DYNFRIC:At time : %g  on ID: %llu Particles found: %d Smoothing Length : "
		  "%g Dynamical friction : (%g %g %g) GravAccel (%g %g %g)  \n",
		  All.Time, (unsigned long int)P[i].ID, BPP(i).BHDF_ngb,
		  P[i].Hsml,BPP(i).BHDF_DynFric[0], BPP(i).BHDF_DynFric[1], BPP(i).BHDF_DynFric[2],
		  P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2]);
	}
    }
#endif // AD_DYNFRIC
}

int blackhole_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		       int *ngblist)
{
  struct blackholedata_in local;
  struct blackholedata_out out;
  memset(&out, 0, sizeof(struct blackholedata_out));
#if defined(GM_REPOSITION_ON_STARDENSITY_MAX)
  MyLongDouble stardenspos[3] = { 0., 0., 0. }, stardens = -1.0;
#ifdef GM_USE_VELS_ON_REPOSITIONING
  MyLongDouble stardensvel[3] = { 0., 0., 0. };
#endif
#endif


  out.BH_MinPot = BHPOTVALUEINIT;

  if(mode == 0)
    particle2in_blackhole(&local, target);
  else
    local = BlackholeDataGet[target];

  /* initialize variables before SPH loop is started */
  MyAtLeastDouble h_i2 = local.Hsml * local.Hsml;
  

  /* Now start the actual SPH computation for this particle */
  int startnode, listindex = 0;
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[listindex];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  int numngb =
	    ngb_treefind_blackhole_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag,
					   exportnodecount,
					   exportindex, ngblist);

	  PANIC_IF(numngb > NumPart / DENSITY_LESS_NGB_FACTOR,
		   "black_hole: Number of neighbours found (%d) exceeds allowed number (%d) derived from NumPart/DENSITY_LESS_NGB_FACTOR, to continue, increase DENSITY_LESS_NGB_FACTOR (currently %f, max possible val would be 1) !!",
		   numngb, (int) (NumPart / DENSITY_LESS_NGB_FACTOR), DENSITY_LESS_NGB_FACTOR);

	  if(numngb < 0)
	    return -1;

	  for(int n = 0; n < numngb; n++)
	    {
	      int j = ngblist[n];
#if defined(WINDS) && defined(KD_EXCLUDE_WIND_IN_BH_FEEDBACK_LOOP)
	      int Valid = 1;
	      if(P[j].Mass <= 0)
		Valid = 0;
	      if(P[j].Type == 0)
		if(SphP[j].DelayTime > 0)
		  Valid = 0;
	      if(Valid == 1)
#else
	      if(P[j].Mass > 0)
#endif
		{
		  if(local.Mass > 0)
		    {
		      MyAtLeastDouble dx = NEAREST_X(local.Pos[0] - P[j].Pos[0]);
		      MyAtLeastDouble dy = NEAREST_Y(local.Pos[1] - P[j].Pos[1]);
		      MyAtLeastDouble dz = NEAREST_Z(local.Pos[2] - P[j].Pos[2]);

		      MyAtLeastDouble r2 = dx * dx + dy * dy + dz * dz;
                     #ifdef AD_DYNFRIC
		      MyAtLeastDouble h_soft = All.ForceSoftening[5]*All.ForceSoftening[5];
                      if(All.BlackHoleFrictionForceOn > 0)

                        if(r2 < h_soft)
                          {
                            if(P[j].Type == 4 || P[j].Type == 1)
                              {

                                out.BHDF_ngb ++ ;
                                MyAtLeastDouble vrel_2 = 0;
                                MyAtLeastDouble vrel[3]={0};
                                for (int k=0; k<3; k++)
                                  { vrel_2 += (P[j].Vel[k]-local.Vel[k])*(P[j].Vel[k]-local.Vel[k]);
                                    vrel[k]=P[j].Vel[k]-local.Vel[k];
                                  }
                                MyAtLeastDouble Lambda = All.ForceSoftening[5]*vrel_2/(All.G*(P[j].Mass+local.Mass));
                                for (int k=0; k<3 ; k++)
                                  {
                                    out.BHDF_DynFric[k]+=log(1+Lambda*Lambda)*P[j].Mass*(local.Mass+P[j].Mass)*vrel[k]/vrel_2/sqrt(vrel_2);
                                  }
                              }
                          }
                     #endif


		      if(r2 < h_i2)
			{
			  MyAtLeastDouble vrel;
			  if(All.BlackHoleRepositioningOn == 1)	/* if this option is switched on, we may also encounter dark matter particles or stars */
			    if(P[j].Potential < out.BH_MinPot)
			      if(P[j].Type == 0 || P[j].Type == 1 || P[j].Type == 4 || P[j].Type == 5)
				{
				  vrel = 0;
				  for(int k = 0; k < 3; k++)
				    vrel += (P[j].Vel[k] - local.Vel[k]) * (P[j].Vel[k] - local.Vel[k]);	/* compute relative velocities */
				  vrel = sqrt(vrel) / ascale;
				  if(vrel <= 0.25 * local.Csnd)
				    {
				      out.BH_MinPot = P[j].Potential;
				      for(int k = 0; k < 3; k++)
					out.BH_MinPotPos[k] = P[j].Pos[k];
				    }
				}


#ifdef GM_REPOSITION_ON_STARDENSITY_MAX
			  if(All.BlackHoleRepositioningOn == 1)	/* if this option is switched on, stars are used here */
			    if(P[j].Type == 4)	//here we only consider stars
			      if(P[j].StarDensity > out.BH_StarDens)
				{
				  {
				    /* compute relative velocities */

				    for(int k = 0, vrel = 0; k < 3; k++)
				      vrel += (P[j].Vel[k] - local.Vel[k]) * (P[j].Vel[k] - local.Vel[k]);

				    vrel = sqrt(vrel) / ascale;

				    if(vrel <= GM_VREL_REPOSITION)
				      {
					out.BH_StarDens = P[j].StarDensity;
					for(int k = 0; k < 3; k++)
					  {
					    out.BH_StarDensPos[k] = P[j].Pos[k];
#ifdef GM_USE_VELS_ON_REPOSITIONING
					    out.BH_StarDensVel[k] = P[j].Vel[k];
#endif

					  }
					if(All.BlackHoleDetails >= 2)
					  fprintf(FdBlackHolesDetails,
						  "ThisTask=%d, time=%g: reposition on stellar densiy max BH %llu vrel=%g dens=%g \n",
						  ThisTask, All.Time, (unsigned long long) local.ID, vrel,
						  out.BH_StarDens);
				      }
				  }

				}
#endif /* ends GM_REPOSITION_ON_STARDENSITY_MAX */


			  if(P[j].Type == 5)	/* we have a black hole merger */
			    if(local.ID != P[j].ID)
			      {
				vrel = 0;
				for(int k = 0; k < 3; k++)
				  vrel += (P[j].Vel[k] - local.Vel[k]) * (P[j].Vel[k] - local.Vel[k]);	/* compute relative velocity of BHs */
				vrel = sqrt(vrel) / ascale;
				MyAtLeastDouble rrel = sqrt(r2);
				MyAtLeastDouble brel =
				  fabs(P[j].Potential - local.Potential) / ascale + vrel * vrel;
				// BlackHoles should meger if fullfill soundspeed, distance and binding criteria (if set)
				if(
#ifdef GM_USE_ABSVAL_IN_VREL
				    ((vrel < All.BlackHoleMergeCsndFrac / ascale)
				     || (All.BlackHoleMergeCsndFrac <= 0))
#else
				    ((vrel < All.BlackHoleMergeCsndFrac * local.Csnd)
				     || (All.BlackHoleMergeCsndFrac <= 0))
#endif
				    && ((rrel < All.BlackHoleMergeDistFrac * All.SofteningTable[5])
					|| (All.BlackHoleMergeDistFrac <= 0))
				    && ((brel < All.BlackHoleMergeBindingFrac * local.Csnd * local.Csnd)
					|| All.BlackHoleMergeBindingFrac <= 0))
				  {
#ifdef _OPENMP
#pragma omp critical(_bh5)
#endif
				    {	// Second condition captures the right selection in rare cases of multiple (>2) swallow posibilities
				      MyAtLeastDouble P_j_Pot = P[j].Potential, l_Pot = local.Potential;
#ifdef KD_USE_IPOT_FOR_BH_MERGER
				      P_j_Pot = P[j].IPotential;
				      l_Pot = local.IPotential;
#endif
#ifdef KD_COMBINE_POT_AND_MASS_FOR_MERGER
				      if(P_j_Pot > 0)
					P_j_Pot /= log(BPP(j).BH_Mass);
				      else
					P_j_Pot *= log(BPP(j).BH_Mass);
				      if(l_Pot > 0)
					l_Pot /= log(local.BH_Mass);
				      else
					l_Pot *= log(local.BH_Mass);
#endif
				      if(P_j_Pot > l_Pot && BPP(j).SwallowPot < l_Pot)
					{
					  BPP(j).SwallowID = local.ID;
					  BPP(j).SwallowPot = l_Pot;
					}
				    }
				  }
				else
				  {
				    if(All.BlackHoleDetails >= 2)
				      fprintf(FdBlackHolesDetails,
					      "SWALLOW: ThisTask=%d, time=%g: id=%llu would like to swallow %llu, but vrel=%g csnd=%g rrel=%g soft=%g brel=%g\n",
					      ThisTask, All.Time, (unsigned long long) local.ID,
					      (unsigned long long) P[j].ID, vrel, local.Csnd, rrel,
					      All.SofteningTable[5], brel);
				  }

			      }	// local.ID != P[j].ID

			  if(P[j].Type == 0)	/* here we have a gas particle */
			    {
			      MyAtLeastDouble r = sqrt(r2);
			      MyAtLeastDouble hinv, hinv3, hinv4;
			      kernel_hinv(local.Hsml, &hinv, &hinv3, &hinv4);
			      MyAtLeastDouble u = r * hinv;
			      MyAtLeastDouble wk, dwk;
			      kernel_main(u, hinv3, hinv4, &wk, &dwk, 0);

#ifdef _OPENMP
#pragma omp critical(_bh0)
#endif
			      {
				int number_of_slices_generated = 0;
				if(All.BlackHoleSwallowGasOn >= 1)	/* compute accretion probability */
				  {
				    MyAtLeastDouble acc_mass = P[j].Mass;
				    if(All.bits > 0 && All.BlackHoleAccretionSlicesOn == 1)
				      {
					number_of_slices_generated =
					  (P[j].ID >> (sizeof(MyIDType) * 8 - All.bits));
					acc_mass = acc_mass / (GENERATIONS - number_of_slices_generated);
				      }

				    MyAtLeastDouble sm = local.BH_Mdot * local.Dt;
				    double p = -1;
				    if(local.BH_Mass > local.Mass)
				      p = P[j].Mass / acc_mass * (1 - exp(-sm / P[j].Mass));

				    double w = get_random_number(P[j].ID);	/* compute random number, uniform in [0,1] */

				    if(All.BlackHoleSwallowGasOn == 2)
				      {
#ifdef LT_STELLAREVOLUTION
					MyFloat xclouds = SphP[j].XColdCloud;
					if(xclouds > 0.01)	//Here we want to swallow only the star forming particles
#else
					// FIX IT: need to compute xclouds and Temperature here !!
					// get_starformation_rate(j, &Temperature, &xclouds);
					PANIC
					  ("BlackHoleSwallowGasOn == 2, need xcloud to decide which particles to swallow, but LT_STELLAREVOLUTION is off. xclouds computation not yet implemented in this case.");
#endif
					if(w < p)
					  if(SPP[j].SwallowID < local.ID)
					    SPP[j].SwallowID = local.ID;
				      }
				    else if(All.BlackHoleSwallowGasOn == 1)	//In this case we swallow the particles regardless their state
				      {
					if(w < p)
					  if(SPP[j].SwallowID < local.ID)
					    SPP[j].SwallowID = local.ID;
				      }
				  }

				if(P[j].Mass > 0)
				  {
				    MyAtLeastDouble energy =
				      local.BH_TotalFeedbackEfficiency * local.BH_Mdot *
				      local.Dt * pow(C_LIGHT / All.UnitVelocity_in_cm_per_s, 2);

				    if(All.BlackHoleLimitMaxTemp > 0)	/* Limit maximum temperature */
				      {
					MyAtLeastDouble u_to_temp_fac =
					  (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN *
					  GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
					MyAtLeastDouble temp = u_to_temp_fac * energy / P[j].Mass;
					if(temp > All.BlackHoleLimitMaxTemp)
					  {
					    VERBOSE(1,
						    "BH: WARNING: task=%d: id=%llu, m=%e large Injected_BH_Energy in blackhole: %e+%e would be temp=%e\n",
						    ThisTask, (unsigned long long) P[j].ID, (float) P[j].Mass,
						    (float) SphP[j].Injected_BH_Energy, (double) energy,
						    (double) temp);
					    energy = All.BlackHoleLimitMaxTemp / u_to_temp_fac * P[j].Mass;
					  }
				      }

				    if(local.BH_Density > 0)
				      {
					MyAtLeastDouble meddington = Mdot_Eddington(local.BH_Mass);
					int was_heated = 1;
					if(All.BlackHoleLimitFeedbackOn == 0 || (All.BlackHoleLimitFeedbackOn & 2) == 2 || (local.BH_Mdot / meddington) < 0.05)	/* no restriction in heating */
					  SphP[j].Injected_BH_Energy +=
					    FLT(energy * P[j].Mass * wk / local.BH_Density);
					else
					  {
					    MyFloat Temperature, xclouds;
#ifdef LT_STELLAREVOLUTION
					    Temperature = SphP[j].Temperature;
					    xclouds = SphP[j].XColdCloud;
#else
					    // FIX IT: need to compute Temperature here !!
					    // get_starformation_rate(j, &Temperature, &xclouds);
					    PANIC
					      ("Trying to limit feedback but LT_STELLAREVOLUTION is off. Temperature computation not yet implemented in this case.");
#endif
					    if((All.BlackHoleLimitFeedbackOn & 1) == 1
					       && Temperature > All.BlackHoleColdGasTemperatureTresh
					       && xclouds < 0.01)
					      SphP[j].Injected_BH_Energy +=
						FLT(energy * P[j].Mass * wk / local.BH_Density);
					    else
					      was_heated = 0;
					  }

					if(SPP[j].SwallowID == local.ID && was_heated == 1)
					  SphP[j].Injected_BH_Energy -=
					    FLT(energy * P[j].Mass /
						(GENERATIONS -
						 number_of_slices_generated) * wk / local.BH_Density);

				      }	// if(local.BH_Density > 0)
				    if(All.BlackHoleDetails >= 3)
				      fprintf(FdBlackHolesDetails,
					      "ENERGY: idBH=%llu time=%g: mbh=%g mdot=%g dt=%g: idGas=%llu gets %g (of %g available)\n",
					      (unsigned long long) local.ID, All.Time, local.BH_Mass,
					      local.BH_Mdot, local.Dt, (unsigned long long) P[j].ID,
					      SphP[j].Injected_BH_Energy, energy);
				  }	// if(P[j].Mass > 0)
			      }	// critical
			    }	// if(P[j].Type == 0)
			}	// if(r2 < h_i2)
		    }		// if(mass > 0)
		}		// (P[j].Mass > 0)
	    }			// for
	}			// while

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }				// while



  /* Now collect the result at the right place */
  if(mode == 0)
    out2particle_blackhole(&out, target, 0);
  else
    BlackholeDataResult[target] = out;

  return 0;
}

void *blackhole_evaluate_primary(void *p)
{
  int thread_id = *((int *) p);
  int i;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;

  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));

  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
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
		NextParticle = NextActiveParticle[i];
		if(NextParticle < 0)
		  break;
		ProcessedFlag[i] = 1;
	      }
	    while(P[i].Type != 5);
	    ProcessedFlag[i] = 0;
#else
	    i = NextParticle;
	    ProcessedFlag[NextParticle] = 0;
	    NextParticle = NextActiveParticle[i];
#endif
	  }
      }
      if(exitFlag)
	break;

      if(P[i].Type == 5)
	{
	  if(blackhole_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    {
	      break;		// buffer was too small to add this BH, we exit without flagging it as processed.
	    }
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */
    }
  return NULL;
}

void *blackhole_evaluate_secondary(void *p)
{
  int thread_id = *((int *) p);
  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));

  while(1)
    {
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	j = NextJ;
	NextJ++;
      }

      if(j >= Nimport)
	break;

      blackhole_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);

    }

  return NULL;

}



#endif
