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
#ifndef __GREENTREE_H__
#define __GREENTREE_H__

#include <stdlib.h>
#include <stdio.h>

#ifdef AR_TREE_STATS
#include <map>
#endif

#ifndef KD_EXTEND_SEARCH
#define KD_EXTEND_SEARCH 1.00001
#endif

#define PTWO(x) x*x


/*

is_active functions for hydro and conduction criteria. Maybe they should sttay in their respective modules (as density_isactive).

*/

static inline int dummy_isactive(int n)
{
  return 1;			// every particle is a good particle
}

static inline int hydro_isactive(int n)
{
#ifdef BLACK_HOLES
  return (P[n].Type == 0 && P[n].Mass > 0);
#else
  return (P[n].Type == 0);
#endif
}

#if defined(fSIDM) || defined(rSIDM)
static inline int msidm_isactive(int n)
{
  return (P[n].Type == 1);
};
#endif


static inline int conduction_isactive(int n)
{
  return P[n].Type == 0 && P[n].Mass > 0;
}

static inline int cr_diffusion_isactive(int n)
{
  return P[n].Type == 0 && P[n].Mass > 0;
}





/*
  criterias for breaking the grouping of particles, e.g. KD_FRICTION will check for different types.
  I cannot hard code in the following routines because for instance density DataGet has Type while hydro DataGet has not. 
 */


//will be either particle_data or DataGet depending on the I or II phase
static inline int group_criteria_distance(MyFloat dist_from_center_2, MyFloat max_allowed_dist_2, void *P,
					  int next_particle, int first_particle)
{
  return dist_from_center_2 <= max_allowed_dist_2;
}


#ifdef KD_FRICTION
template < typename group_particle_data >	//will be either particle_data or DataGet depending on the I or II phase
static inline int group_criteria_w_ptype(MyFloat dist_from_center_2, MyFloat max_allowed_dist_2,
					 group_particle_data P, int next_particle, int first_particle)
{
  if(dist_from_center_2 > max_allowed_dist_2)
    return 0;
  return P[first_particle].Type == P[next_particle].Type;
}
#endif


/*
  criteras for filtering out particles in tree-walks
*/
static inline int treewalk_filter_gas_only(particle_data *P, int p, void *DataGet, int next_particle)
{
  return P[p].Type > 0;
}


#ifdef KD_FRICTION
template < typename group_particle_data >
  static inline int treewalk_filter_friction_criteria(particle_data *P, int p, group_particle_data P_In,
						      int next_particle)
{
  int targettype = P_In[next_particle].Type;

#if !defined(KD_FRICTION) && !defined(AXION_DM)
  return P[p].Type > 0;

#else //else of !defined(KD_FRICTION) && !defined(AXION_DM)
#ifndef AXION_DM
  if(P[p].Type > 0 && targettype != 5)
    return 1;
#else
  if(!((1 << P[p].Type) & (targettype)))
    return 1;
#endif
#endif //end else  !defined(KD_FRICTION) && !defined(AXION_DM)
  return 0;
}
#endif


/*
   For debugging purposes it is useful to know which module is executing one of the following functions.
   This struct is  a parameter of most of the functions below.
*/
enum Module
{
  Module_Gravity,
  Module_Density,
  Module_Hydro,
  Module_Conduction,
  Module_mSIDM,
  Module_Density_pt1,
  Module_CRDiffusion
};


/*
  struct to identify neighbour finding criteria
*/
enum NeighbourSearchCriteria
{
  NeighbourSearchCriteria_Density,	/* tree walk only takes particles within a given searching radius h */
  NeighbourSearchCriteria_Hydro,	/* tree walk takes particles whose smothing length intercept a given searching radius h */
#if defined(fSIDM) || defined(rSIDM)
  NeighbourSearchCriteria_mSIDM,
#ifndef ADAPTGRAVSOFT
  NeighbourSearchCriteria_Density_pt1
#endif
#endif
};

/*
  how to choose NextParticle? NextPArticleLoop_ActiveParticleList uses Active PArticle List,
  while NextPArticleLoop_ZeroToNgas counts from 0 to N_gas (as done in conduction.c)
*/
enum NextParticleLoop
{
  NextParticleLoop_ActiveParticleList,
  NextParticleLoop_ZeroToNgas,
};

/*
  Types of Neighbours SEarch
  Ideally, in the future, you can switch to 'good old Gadget2 tree walk' for debugging purposes
 */
enum TreeWalkType
{
  TreeWalkType_GreenTreeWalk,
  TreeWalkType_Gadget3TreeWalk
};

enum IterationPhase
{
  IterationPhase_IPhaseComplete,	/* neighbour search + fills of export buffer */
  IterationPhase_IIPhase,	/* second phase */
  IterationPhase_IPhaseExportBuffer,	/* only export buffer (used when doing neighbour search on GPU) */
  IterationPhase_IPhaseComputation,	/* only neighbour search (possibly used when mergin OpenACC code to this routines) */
};


//ultra unified neighbour finding routine.
template < typename TreewakFilterCriteria, typename data_in > static inline int ar_green_ngb_openmp(MyLongDouble searchcenter[3], MyFloat phsml, MyFloat pradius, int *startnode, enum IterationPhase mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int *targets, int n_targets, int should_drift, enum NeighbourSearchCriteria tree_walk_search_type, int max_ngbs, enum Module module, TreewakFilterCriteria treewalk_filter_criteria, data_in DataGet)	/* max_ngbs==0 => all neighbours */
{

#ifdef DO_NOT_BRACH_IF
#error "tree walk never tested with DO_NOT_BRACH_IF - but it may work btw"
#endif

  int no, p, numngb, task, nexp;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  MyFloat hsml = phsml;		//* KD_EXTEND_SEARCH;
  MyFloat radius = pradius;	// * KD_EXTEND_SEARCH;
  MyFloat hsmlpradius = hsml + radius;
  integertime ti_Current = All.Ti_Current;
#ifndef DO_NOT_BRACH_IF
  MyLongDouble dist, dx, dy, dz;
#ifdef PERIODIC
  MyLongDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  numngb = 0;

  no = *startnode;

  while(no >= 0)
    {

      if(max_ngbs > 0 && numngb >= max_ngbs)
	{			/* GPU processes neighbours in chunks of (typically) 32 */
	  *startnode = no;
#ifndef DO_NOT_BRACH_IF
	  return numngb;
#else
	  return ngb_filter_pairs(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
	}

      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

#if defined(fSIDM) || defined(rSIDM)
#ifndef ADAPTGRAVSOFT
	  if(tree_walk_search_type == NeighbourSearchCriteria_mSIDM ||
	     tree_walk_search_type == NeighbourSearchCriteria_Density_pt1)
#else
	  if(tree_walk_search_type == NeighbourSearchCriteria_mSIDM)
#endif
	    {
	      if(P[p].Type != 1)
		continue;
	    }
	  else
#endif
	    {
	      if(treewalk_filter_criteria(P, p, DataGet, targets[0]))
		continue;
	    }

#ifdef BLACK_HOLES
	  if(P[p].Mass == 0)
	    continue;
#endif

	  if(should_drift)
	    if(P[p].Ti_current != ti_Current)
	      {
#pragma omp critical(_partnodedrift_)
		drift_particle(p, ti_Current);
	      }

	  if(mode == IterationPhase_IPhaseExportBuffer)
	    continue;

#ifndef DO_NOT_BRACH_IF

	  if(tree_walk_search_type == NeighbourSearchCriteria_Hydro)
	    dist = DMAX(P[p].Hsml, hsml) + radius;
	  else if(tree_walk_search_type == NeighbourSearchCriteria_Density)
	    dist = hsmlpradius;
#if defined(fSIDM) || defined(rSIDM)
#ifndef ADAPTGRAVSOFT
	  // this is only an approximate criterion
	  else if(tree_walk_search_type == NeighbourSearchCriteria_mSIDM)
	    dist = P[p].Hsml + hsml + radius;
	  else if(tree_walk_search_type == NeighbourSearchCriteria_Density_pt1)
	    dist = hsmlpradius;
#else
	  // this is only an approximate criterion
	  else if(tree_walk_search_type == NeighbourSearchCriteria_mSIDM)
	    dist = P[p].AGS_Hsml + hsml + radius;
#endif
#endif
	  else
	    PANIC("TreeWalkType got wrong");

	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif
#ifdef AR_GRAVCOST_IN_DENSITY
	  P[p].GravCost[TakeLevel] += 0.2;
#endif



	  ngblist[numngb++] = p;	/* Note: unlike in previous versions of the code, the buffer
					   can hold up to all particles */

	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      int i_target = 0;

	      PANIC_IF(mode == IterationPhase_IIPhase,
		       "mode =%d task=%d no=%d maxpart=%d maxnodes=%d +=%d startnode=%d\n", mode, ThisTask,
		       no, All.MaxPart, maxNodes, maxPart + maxNodes, *startnode);


	      if(targets[i_target] >= 0)	/* if no target is given, export will not occur */
		{


		  if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != targets[i_target])
		    {
		      exportflag[task] = targets[i_target];
		      exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(exportnodecount[task] == NODELISTLENGTH)
		    {
		      int exitFlag = 0;


#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
		      {
			if(Nexport + n_targets - 1 >= bunchSize)
			  {
			    // out of buffer space. Need to discard work for this particle and interrupt
			    BufferFullFlag = 1;
			    exitFlag = 1;
			  }
			else
			  {
			    nexp = Nexport;
			    Nexport += n_targets;
			  }
		      }

		      if(exitFlag)
			{
			  return -1;
			}
		      /*
		         #pragma omp atomic capture
		         nexp = Nexport++;
		       */
		      /*                      if(nexp>=bunchSize)
		         {
		         BufferFullFlag=1;
		         Nexport = bunchSize;
		         return -1;
		         }
		       */
		      if(exitFlag)
			return -1;
		      exportnodecount[task] = 0;
		      exportindex[task] = nexp;	//per questo task, indice della coppia nexp->particella
		      for(i_target = 0; i_target < n_targets; i_target++)
			{

			  DataIndexTable[nexp + i_target].Task = task;
			  DataIndexTable[nexp + i_target].Index = targets[i_target];
			  DataIndexTable[nexp + i_target].IndexGet = nexp + i_target;
			}
		    }

		  for(i_target = 0; i_target < n_targets; i_target++)
		    DataNodeList[exportindex[task] + i_target].NodeList[exportnodecount[task]] =
		      DomainNodeIndex[no - (maxPart + maxNodes)];
		  exportnodecount[task]++;
		  if(exportnodecount[task] < NODELISTLENGTH)
		    for(i_target = 0; i_target < n_targets; i_target++)
		      DataNodeList[exportindex[task] + i_target].NodeList[exportnodecount[task]] = -1;


		}

	      no = Nextnode[no - maxNodes];
	      continue;

	    }

	  current = &Nodes[no];

	  if(mode == IterationPhase_IIPhase)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
#error "NOT TESTED WITH DO NOT BRANCH IF"
		  if(tree_walk_search_type == NeighbourSearchCriteria.NeighbourSearchCriteria_Hydro)
		    return ngb_filter_pairs(numngb, ngblist, &vcenter, &box, &hbox, hsmlpradius);
		  if(tree_walk_search_type == NeighbourSearchCriteria.NeighbourSearchCriteria_Density)
		    return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsmlpradius);

#endif
		}
	    }


	  if(should_drift)
	    if(current->Ti_current != ti_Current)
	      {
#pragma omp critical
		force_drift_node(no, ti_Current);
	      }




	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }


	  /* on I phase hybrid in CPU, we stay only in the toplevel nodes to fill the export buffer */
	  if(mode == IterationPhase_IPhaseExportBuffer && !(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL)))
	    {
	      no = Nodes[no].u.d.sibling;
	      continue;
	    }

#ifndef DO_NOT_BRACH_IF


	  if(tree_walk_search_type == NeighbourSearchCriteria_Hydro)
	    dist = DMAX(Extnodes[no].hmax, hsml) + 0.5 * current->len + radius;
	  if(tree_walk_search_type == NeighbourSearchCriteria_Density)
	    dist = hsmlpradius + 0.5 * current->len;
#if defined(fSIDM) || defined(rSIDM)
	  // this is only an approximate criterion
	  if(tree_walk_search_type == NeighbourSearchCriteria_mSIDM)
	    dist = current->maxsoft + hsml + 0.5 * current->len + radius;
#ifndef ADAPTGRAVSOFT
	  if(tree_walk_search_type == NeighbourSearchCriteria_Density_pt1)
	    dist = hsmlpradius + 0.5 * current->len;
#endif
#endif
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, DMAX(Extnodes[no].hmax, hsml) + radius);
#endif
	}
    }


  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_pairs(numngb, ngblist, &vcenter, &box, &hbox, hsmlpradius);
#endif
}

/*
// GreenTree grouping on the first phase
//It takes an 'is_active' criteria, a NextParticleLoop criteria (to know how to 'increment' NextParticle'),
// a list of 'group_ids' (output) where to store the `group_len` grouepd particles.
max_r2 and &max_h are output variables storing the maximum distnace square between particles and the maximum h.
group_particles is a flag: 0 do not group (good old Gadget2 tree walk), 1 do the green tree
module stores the current module, for debugging purposes.
*/
template < typename IsActiveFunction, typename GroupCriteria >
  static inline int ar_pick_next_particle_to_group_in_i_phase(IsActiveFunction is_active,
							      NextParticleLoop next_particle_loop,
							      int *group_ids, int &group_len,
							      MyFloat & max_r2, MyFloat & max_h,
							      int group_particles, Module module,
							      int &processed_particles,
							      GroupCriteria group_criteria)
{
  int first_group_particle = -1;
  MyLongDouble xtmp, dx, dy, dz, dist2, dist, cluster_dist_2 = -1.;

  group_len = 0;
  if(BufferFullFlag != 0 || NextParticle < 0)
    return -1;

  while(1)
    {				/* we search for the first particle of the group, but it must be active! */
      if((next_particle_loop == NextParticleLoop_ActiveParticleList && NextParticle < 0)
	 || (next_particle_loop == NextParticleLoop_ZeroToNgas && NextParticle == N_gas))
	{
	  /*we reached the last particle without even starting to group them.. gadget iteration is done! */
	  return -1;
	}
      if(!is_active(NextParticle))
	{			/* skip the particle */
	  ProcessedFlag[NextParticle] = 1;
	  if(next_particle_loop == NextParticleLoop_ActiveParticleList)
	    NextParticle = NextActiveParticle[NextParticle];
	  if(next_particle_loop == NextParticleLoop_ZeroToNgas)
	    NextParticle++;
	  processed_particles++;
	  continue;
	}

      max_r2 = 0.;
      first_group_particle = NextParticle;

#if (defined(fSIDM) || defined(rSIDM)) && defined(ADAPTGRAVSOFT)
      max_h = P[first_group_particle].AGS_Hsml;
#else
      max_h = P[first_group_particle].Hsml;
#endif
      group_ids[group_len] = first_group_particle;
      group_len = 1;

#if (defined(fSIDM) || defined(rSIDM)) && defined(ADAPTGRAVSOFT)
      cluster_dist_2 = P[first_group_particle].AGS_Hsml * AR_GREEN_TREE_HSML_FRACTION;
#else
      cluster_dist_2 = P[first_group_particle].Hsml * AR_GREEN_TREE_HSML_FRACTION;
#endif
      cluster_dist_2 *= cluster_dist_2;
      ProcessedFlag[first_group_particle] = 0;
      if(next_particle_loop == NextParticleLoop_ActiveParticleList)
	NextParticle = NextActiveParticle[NextParticle];
      if(next_particle_loop == NextParticleLoop_ZeroToNgas)
	NextParticle++;

      break;
    }
  //  return 1;
  if(!group_particles)
    return 1;

#ifdef AR_ONE_PARTICLE_PER_GROUP_I
  return 1;
#endif

  /* I prefear to do black holes alone - may not be necessary */
  if(P[first_group_particle].Type == 5)
    return 1;

  while(1)
    {
      if((next_particle_loop == NextParticleLoop_ActiveParticleList && NextParticle < 0)
	 || (next_particle_loop == NextParticleLoop_ZeroToNgas && NextParticle == N_gas))
	{

	  break;
	}
      if(!is_active(NextParticle))
	{
	  ProcessedFlag[NextParticle] = 1;
	  if(next_particle_loop == NextParticleLoop_ActiveParticleList)
	    NextParticle = NextActiveParticle[NextParticle];
	  if(next_particle_loop == NextParticleLoop_ZeroToNgas)
	    NextParticle++;

	  continue;
	}



#ifdef AR_GREEN_TREE_MAX_GROUP_SIZE
      if(group_len >= AR_GREEN_TREE_MAX_GROUP_SIZE)
	break;
#endif

      /* I prefear to do black holes alone - may not be necessary */
      if(P[NextParticle].Type == 5)
	break;

      dx = NGB_PERIODIC_LONG_X(P[first_group_particle].Pos[0] - P[NextParticle].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(P[first_group_particle].Pos[1] - P[NextParticle].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(P[first_group_particle].Pos[2] - P[NextParticle].Pos[2]);
      dist2 = dx * dx + dy * dy + dz * dz;

      if(group_criteria(dist2, cluster_dist_2, P, NextParticle, first_group_particle))
	{			/* found a new particle to group */


	  if(dist2 > max_r2)
	    max_r2 = dist2;	//update max group size
#if (defined(fSIDM) || defined(rSIDM)) && defined(ADAPTGRAVSOFT)
	  if(P[NextParticle].AGS_Hsml > max_h)
	    max_h = P[NextParticle].AGS_Hsml;
#else
	  if(P[NextParticle].Hsml > max_h)
	    max_h = P[NextParticle].Hsml;
#endif
	  group_ids[group_len] = NextParticle;	//add particle to group
	  group_len++;
	  ProcessedFlag[NextParticle] = 0;
	  if(next_particle_loop == NextParticleLoop_ActiveParticleList)
	    NextParticle = NextActiveParticle[NextParticle];
	  if(next_particle_loop == NextParticleLoop_ZeroToNgas)
	    NextParticle++;


	}
      else
	{
	  break;		/* ...next! */
	}

    }				/* while ends */
  return 1;			/*ok! */
}

/*
group guest particles.
DataGet is a templated DensityDataGet/HydroDataGet etc..
first_group_particle is an output array with the list grouepd particles
group_len is the output size of the group
max_r2 and   max_h output variables with respectively the maximum distance square or the smoothing lenght
module is for debugging purpose
*/
template < typename data_in, typename GroupCriteria >
  static inline int ar_pick_next_particle_to_group_in_ii_phase(data_in DataGet, int &first_group_particle,
							       int &group_len, MyFloat & max_r2,
							       MyFloat & max_h, Module module,
							       GroupCriteria group_criteria)
{
  first_group_particle = -1;
  MyLongDouble xtmp, dx, dy, dz, dist2, dist, cluster_dist_2 = -1.;
  /*
     if(NextJ<Nimport){
     first_group_particle=NextJ;
     group_len=1;
     NextJ++;
     max_r2=0.;
     max_h = DataGet[first_group_particle].Hsml;
     return 1;
     }else{
     return -1;
     }
   */
  while(1)
    {
      int n_nodes;
      int same_nodes;
      if(NextJ >= Nimport)
	{			/* we reached the last particle: end import */
	  if(first_group_particle == -1)
	    return -1;
	  break;
	}

      if(first_group_particle == -1)
	{


	  first_group_particle = NextJ;
	  group_len++;
	  NextJ++;
	  n_nodes = 0;



	  while(DataNodeListGet[first_group_particle].NodeList[n_nodes] >= 0)
	    {			/* II phase group also particles with same NodeList */
	      n_nodes++;
	      if(n_nodes == NODELISTLENGTH)
		{
#if defined(fSIDM) || defined(rSIDM)
		  if(module == Module_mSIDM)
		    {
		      std::
			cout << "mSIDM Warning: pseudo particle, NODELISTLENGTH is too small (" <<
			NODELISTLENGTH << ")" << std::endl;
#ifndef mSIDM_IGNOREWARNING
		      PANIC_ANY("This can be resolved by increasing NODELISTLENGTH_mSIDM in Config.sh.");
#endif
		    }
#endif
		  break;
		}

	    }

	  max_h = DataGet[first_group_particle].Hsml;
	  max_r2 = 0.;
	  cluster_dist_2 = DataGet[first_group_particle].Hsml * AR_GREEN_TREE_HSML_FRACTION;
	  cluster_dist_2 *= cluster_dist_2;

#if defined(LT_STELLAREVOLUTION) || defined(KD_FRICTION_DYNAMIC)
	  if(DataGet[first_group_particle].Type == 5)
	    break;		/* black holes are don alone */
#endif
	  continue;
	}

      data_in D = &DataGet[first_group_particle];
#if defined(LT_STELLAREVOLUTION) || defined(KD_FRICTION_DYNAMIC)
      //    if(D->Type==5)break; /* black hole will be processed in the next group - alone*/
#endif
      same_nodes = 1;

      int i_node = 0;

      while(DataNodeListGet[NextJ].NodeList[i_node] >= 0)
	{			/* check if the current particle share the same node list of all particles in the group  */
	  if(DataNodeListGet[first_group_particle].NodeList[i_node] !=
	     DataNodeListGet[NextJ].NodeList[i_node])
	    {
	      same_nodes = 0;
	      break;
	    }
	  i_node++;

	  if(i_node > n_nodes)
	    {
	      same_nodes = 0;
	      break;
	    }
	  if(i_node == NODELISTLENGTH)
	    break;
	}
      if(i_node != n_nodes)
	{
	  same_nodes = 0;
	}

      if(same_nodes != 1)
	break;


      data_in N = &DataGet[NextJ];

#ifdef AR_ONE_PARTICLE_PER_GROUP_II
      break;
#endif
      if(same_nodes != 1)
	break;
      dx = NGB_PERIODIC_LONG_X(D->Pos[0] - N->Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(D->Pos[1] - N->Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(D->Pos[2] - N->Pos[2]);
      dist2 = dx * dx + dy * dy + dz * dz;

      if(group_criteria(dist2, cluster_dist_2, DataGet, NextJ, first_group_particle) && same_nodes == 1)
	{
	  if(dist2 > max_r2)
	    max_r2 = dist2;
	  if(N->Hsml > max_h)
	    max_h = N->Hsml;
	  group_len++;
	  NextJ++;

	}
      else
	break;

    }

  return 1;			/*ok! */
}


/*
Take the list of grouped particles and a possible list of nodes (in case of the secondary phase), and do:
- for each node :
- - perform both a single tree walk for the whole group
- - for each particle in the group:
- - -  do evaluate().
AR.
 */



template < typename IsActiveFunction, typename EvaluateFunction, typename DataIn, typename TreeWalkFilter >
  static inline int ar_green_tree_inner_loop(DataIn P_In, IsActiveFunction is_active,
					     EvaluateFunction evaluate,
					     NeighbourSearchCriteria tree_walk_search_type, int *targets,
					     int n_targets, int thread_id, MyFloat max_r, MyFloat max_h,
					     IterationPhase mode, int should_drift, int *nodes_list,
					     Module module, void *data, TreeWalkFilter tree_walk_filter)
{

  int *exportflag, *exportnodecount, *exportindex, *ngblist;

#ifdef DENSITY_LESS_NGB_FACTOR
  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
#else
  ngblist = Ngblist + thread_id * NumPart;
#endif

  if(mode != IterationPhase_IIPhase)
    {
      exportflag = Exportflag + thread_id * NTask;
      exportnodecount = Exportnodecount + thread_id * NTask;
      exportindex = Exportindex + thread_id * NTask;
    }
  else
    {
      exportflag = exportnodecount = exportindex = NULL;
    }
  int is_loop_ok = 1;



  int i_node = 0;		/*on the secondary phase this will increase with the varioujs NODELISTs, and will signal the evaluate function that it has to to DensDataOut instead of assigning. */



  do
    {				/* loop over the list of nodes, that is not empty in the II phase */

      int extended_numngb_inbox = -2;
      int evaluate_in_node = All.MaxPart;
      int target0 = targets[0];	/* the tree walk is performed around the first particle, ofc over a length of max_h + max_r */

      if(i_node >= NODELISTLENGTH)	/* we finished the nodes to evaluate */
	return 1;


      if(nodes_list != NULL)
	{			/*if we got a node list, we use it */
	  evaluate_in_node = nodes_list[i_node];
	}

      if(evaluate_in_node < 0)
	return 1;		/* we finished the nodes to evaluate */

      if(evaluate_in_node > All.MaxPart + MaxNodes)
	{
	  for(int i = 0; i < NODELISTLENGTH; i++)
	    {
	      printf("%d ", nodes_list[i]);

	    }
	  printf("\n");
	}
      PANIC_IF(evaluate_in_node > All.MaxPart + MaxNodes,
	       "Bad node no=%d > %d+%d = %d; i_node=%d<NODELISTLENGTH=%d, nodes=%p",
	       evaluate_in_node, All.MaxPart, MaxNodes, All.MaxPart + MaxNodes, i_node, NODELISTLENGTH,
	       nodes_list);
      PANIC_IF(Nodes[evaluate_in_node].u.d.nextnode >= All.MaxPart + MaxNodes,
	       "Very bad nextnode = %d > %d; mode=%d", Nodes[evaluate_in_node].u.d.nextnode,
	       All.MaxPart + MaxNodes, mode);

      if(nodes_list != NULL)
	{			/* enter the node if mode==1 */
	  evaluate_in_node = Nodes[evaluate_in_node].u.d.nextnode;
	}

      int i_startnode = evaluate_in_node;


      MyLongDouble *pos = P_In[target0].Pos;


      extended_numngb_inbox = ar_green_ngb_openmp(pos, max_h, max_r, &i_startnode, mode,
						  exportflag, exportnodecount, exportindex,
						  ngblist, targets, n_targets, should_drift,
						  tree_walk_search_type, 0, module, tree_walk_filter, P_In);
      PANIC_IF(mode == IterationPhase_IPhaseExportBuffer
	       && extended_numngb_inbox > 0,
	       "mode == IterationPhase_IPhaseExportBuffer && extended_numngb_inbox = %d>0",
	       extended_numngb_inbox);

      /*
         MyFloat x = max_h + max_r;
         extended_numngb_inbox = ar_green_ngb_treefind_variable_threads(pos, &x, &i_startnode, (mode==IterationPhase::IterationPhase_IPhaseComplete)?0:1,
         exportflag, exportnodecount, exportindex,
         ngblist, targets, n_targets, should_drift);
       */
      PANIC_IF(extended_numngb_inbox < 0
	       && mode == IterationPhase_IIPhase, "Error in ngb search; maxh=%f maxr=%f", max_h, max_r);

      if(extended_numngb_inbox < 0)
	return -1;



#ifdef DENSITY_LESS_NGB_FACTOR
      PANIC_IF(extended_numngb_inbox > ((int) (NumPart / DENSITY_LESS_NGB_FACTOR)),
	       "mode=%d,  numngb = %d max=numpart/density_less_ngb_factor=%d \n", mode, extended_numngb_inbox,
	       ((int) (NumPart / DENSITY_LESS_NGB_FACTOR)));
#endif
      PANIC_IF(extended_numngb_inbox > ((int) (NumPart)), " numngb=%d max=numpart=%d \n",
	       extended_numngb_inbox, ((int) (NumPart)));

      int i_target = 0;



#ifdef GREEN_TREE_OUTER_LOOP_NUM_THREADS
      int num_threads_inner_loop = maxThreads / GREEN_TREE_OUTER_LOOP_NUM_THREADS;
#pragma omp parallel for num_threads(num_threads_inner_loop)
#error ""
#endif
      for(i_target = 0; i_target < n_targets; i_target++)
	{			/* loop over the grouped particles */
	  //    VERBOSE(-1,"1. evaluate itarget=%d",i_target);
	  if(!is_loop_ok)
	    continue;
	  int target = targets[i_target];
	  PANIC_IF(mode != IterationPhase_IIPhase
		   && !is_active(target), "Processing a non active particle. mode=%d, target=%d is_active=%d",
		   mode, target, is_active(target));
	  //      VERBOSE(-1,"2. evaluate %d i_startnode=%d maxpart=%d",target, i_startnode, All.MaxPart);


	  if(mode != IterationPhase_IPhaseExportBuffer &&
	     evaluate(target, mode, ngblist, extended_numngb_inbox, &i_startnode, evaluate_in_node, i_node,
		      data) < 0)
	    {
	      is_loop_ok = 0;
	      PANIC_IF(mode == IterationPhase_IIPhase, "tragic error in evaluate");

	    }
	  else
	    {
	      if(mode == IterationPhase_IPhaseComplete || mode == IterationPhase_IPhaseExportBuffer)

		ProcessedFlag[target] = 1;	/* particle successfully finished */

	    }

	}
      if(is_loop_ok != 1)
	return -1;		/* buffer is full */
      if(nodes_list != NULL)
	{
	  i_node++;
	}
    }
  while((nodes_list != NULL && nodes_list[i_node] >= 0 && i_node < NODELISTLENGTH));
  return is_loop_ok ? 1 : -1;

}

/*

enum NextParticleLoop{
  NextParticleLoop_ActiveParticleList,
  NextParticleLoop_ZeroToNgas,
};

*/

template < typename IsActive, typename Evaluate, typename GroupCriteria, typename TreeWalkFilter >
  static inline void *ar_green_evaluate_primary(IsActive is_active, Evaluate evaluate, Module module,
						TreeWalkType search_type, NextParticleLoop next_particle_loop,
						NeighbourSearchCriteria tree_walk_type,
						GroupCriteria group_criteria, TreeWalkFilter treewalk_filter,
						void *data, int is_primary_on_gpu, int should_drift,
						int group_particles)
{
  int thread_id = *(int *) data;

  int i, j;
  IterationPhase mode = IterationPhase_IPhaseComplete;
  if(is_primary_on_gpu)
    mode = IterationPhase_IPhaseExportBuffer;



  int *recyclingbuffer, *exportflag;

#ifdef DENSITY_LESS_NGB_FACTOR
  recyclingbuffer = RecyclingBuffer + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
#else
  recyclingbuffer = RecyclingBuffer + thread_id * NumPart;
#endif

  if(search_type == TreeWalkType_Gadget3TreeWalk)
    group_particles = 0;
  int *group_ids;
  MyLongDouble xtmp;
  int k;
  group_ids = recyclingbuffer;

  exportflag = Exportflag + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  int max_group_len = -1;
  int skipped_active = 0;
  int processed_particles = 0;
  int exit_flag = 0;

  /* loop until either we did all particls or buffer is full */

  int grouping_iter = -1;
  PANIC_IF(thread_id == -123, "GREPME in ar_green_evaluate_primary");


  while(1)
    {
      grouping_iter++;		/* for debugging purpose */
      int break_em_all = 0;
      int group_len = 0;
      int first_group_particle = -1;
      MyFloat max_r2, max_h, max_r1, max_h_p_r;

      if((next_particle_loop == NextParticleLoop_ActiveParticleList && NextParticle < 0)
	 || (next_particle_loop == NextParticleLoop_ZeroToNgas && NextParticle == N_gas))
	{
	  return NULL;
	}


      //note that processed_particles is passsed as reference, see signature of ar_pick_next_particle_to_group_in_i_phase.
#pragma omp critical
      exit_flag =
	ar_pick_next_particle_to_group_in_i_phase(is_active, next_particle_loop, recyclingbuffer, group_len,
						  max_r2, max_h, group_particles, module, processed_particles,
						  group_criteria);
      /*note: ar_pick_next_particle_to_group get paramas by reference */



      if(exit_flag == -1)
	return NULL;



      if(max_group_len < 0 || max_group_len < group_len)
	max_group_len = group_len;
      max_r1 = sqrt(max_r2);
      max_h_p_r = max_r1 + max_h;



      if(group_len > 0)
	{
	  int was_green_outerloop_ok =
	    ar_green_tree_inner_loop(P, is_active, evaluate, tree_walk_type, group_ids, group_len, thread_id,
				     max_r1, max_h, mode, should_drift, NULL, module, data, treewalk_filter);
	  if(was_green_outerloop_ok != 1)
	    return NULL;	/* ar_first_phase_inner_loop returns -1 on error, as tradition of Gadget */

	}


      continue;
    }


  return NULL;

}

template < typename IsActive, typename Evaluate >
  static inline void *ar_green_evaluate_primary(IsActive is_active, Evaluate evaluate, Module module,
						TreeWalkType search_type, NextParticleLoop next_particle_loop,
						NeighbourSearchCriteria tree_walk_type, void *data,
						int is_primary_on_gpu, int should_drift, int group_particles
						//,GroupCriteria group_criteria, TreeWalkFilter tree_walk_filter
  )
{


  return ar_green_evaluate_primary(is_active, evaluate, module, search_type, next_particle_loop,
				   tree_walk_type, group_criteria_distance, treewalk_filter_gas_only, data,
				   is_primary_on_gpu, should_drift, group_particles);
}

template < typename data_in, typename Evaluate, typename GroupCriteria, typename TreeWalkFilter >
  static inline void *ar_green_evaluate_secondary(data_in DataGet, Evaluate evaluate, Module module,
						  TreeWalkType search_type,
						  NeighbourSearchCriteria ngb_criteria,
						  GroupCriteria group_criteria,
						  TreeWalkFilter treewalk_filter, void *data,
						  int should_drift)
{
  int thread_id = *(int *) data;

  int j, dummy, *ngblist, *recyclingbuffer, processed_particles = 0;

  //auto group_criteria = boh1;
  //auto treewalk_filter = boh2;

#ifdef DENSITY_LESS_NGB_FACTOR
  ngblist = Ngblist + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
  recyclingbuffer = RecyclingBuffer + thread_id * ((int) (NumPart / DENSITY_LESS_NGB_FACTOR));
#else
  ngblist = Ngblist + thread_id * NumPart;
  recyclingbuffer = RecyclingBuffer + thread_id * NumPart;
#endif

  //  printf("Numport=%d NextJ=%d\n",Nimport, NextJ);
  int *group_ids;
  group_ids = recyclingbuffer;
  int exitFlag = 0;
  int i_node, max_group_len = -1;
  enum IterationPhase mode = IterationPhase_IIPhase;
  int exit_flag;
  while(1)
    {
      int n_nodes = 0;

      int same_nodes = 1;
      int group_len = 0;

      MyFloat max_h;
      MyFloat max_r2;
      MyFloat max_r1;
      MyFloat max_r_p_h;
      int first_group_particle;

      /*note: ar_pick_next_particle_to_group get paramas by reference */


#pragma omp critical(_nexport_)
      exit_flag =
	ar_pick_next_particle_to_group_in_ii_phase(DataGet, first_group_particle, group_len, max_r2, max_h,
						   module, group_criteria);


      if(exit_flag != 1)
	return NULL;

      if(first_group_particle == -1 || first_group_particle >= Nimport)
	return NULL;		/* w finished */
      max_r1 = sqrt(max_r2);
      max_r_p_h = max_r1 + max_h;

      if(max_group_len < 0 || max_group_len < group_len)
	max_group_len = group_len;	/* we do statistics */

      for(int k = 0; k < group_len; k++)
	{
	  group_ids[k] = first_group_particle;
	  first_group_particle++;
	}

      int *nodes = DataNodeListGet[group_ids[0]].NodeList;
      int was_green_outerloop_ok =
	ar_green_tree_inner_loop(DataGet, dummy_isactive, evaluate, ngb_criteria, group_ids, group_len,
				 thread_id, max_r1, max_h, IterationPhase_IIPhase, should_drift, nodes,
				 module, data, treewalk_filter);
      PANIC_IF(was_green_outerloop_ok != 1, "II phase failed miserably");
      if(exit_flag == -1)
	return NULL;
    }

  return NULL;

}


template < typename data_in, typename Evaluate >
  static inline void *ar_green_evaluate_secondary(data_in DataGet, Evaluate evaluate, Module module,
						  TreeWalkType search_type,
						  NeighbourSearchCriteria ngb_criteria, void *data,
						  int should_drift)
{
  return ar_green_evaluate_secondary(DataGet, evaluate, module, search_type, ngb_criteria,
				     group_criteria_distance, treewalk_filter_gas_only, data, should_drift);
}





/*
  we do statistics of the tree here below
*/

#ifdef AR_TREE_STATS

inline static int ar_is_tree_index_root(int n)
{
  return n == All.MaxPart;
}

inline static int ar_get_tree_index_levels(int n)
{
  int levels = 0;
  while(!ar_is_tree_index_root(n))
    {
      n = Nodes[n].u.d.father;
      levels++;
    }
  return levels;
}


inline static void ar_tree_stats(void)
{


  std::map < int, int >mappa;
  int nlevels;

  int max_level = 0;

  for(int i = 0; i < NumPart; i++)
    {

      int level = ar_get_tree_index_levels(Father[i]);
      if(mappa.find(level) == mappa.end())
	{
	  mappa[level] = 0;
	}

      mappa[level] += 1;	//mappa[level] + 1;
      if(level > max_level)
	max_level = level;
    }

  for(int i = 0; i <= max_level; i++)
    {
      if(mappa.find(i) != mappa.end())
	{
	  VERBOSE(0, "Level %d has %d particles", i, mappa[i]);
	}
    }

}

#endif

#endif /* __GREENTREE_H_ */
