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
#include <time.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"
#include "../System/vector.h"
#include "../System/communication.h"

/*! \file ngb.c
 *  \brief neighbour search by means of the tree
 *
 *  This file contains routines for neighbour finding.  We use the
 *  gravity-tree and a range-searching technique to find neighbours.
 */





#ifdef DO_NOT_BRACH_IF

#ifdef __xlC__
#pragma alloca
#define ALLOC_STACK(n) alloca(n)
#elif defined(__GNUC__)
#define ALLOC_STACK(n) alloca(n)
#elif defined(__INTEL_COMPILER)
#define ALLOC_STACK(n) _alloca(n)
#else
#define ALLOC_STACK(n) alloca(n)
#endif


int ngb_filter_pairs(long long numngb, int list[], t_vector * center, t_vector * box, t_vector * hbox,
		     MyFloat hsml) __attribute__((noinline));
int ngb_filter_variables(long long numngb, int list[], t_vector * center, t_vector * box, t_vector * hbox,
			 MyFloat hsml) __attribute__((noinline));

int ngb_filter_pairs(long long numngb, int list[], t_vector *center, t_vector *box, t_vector *hbox,
		     MyFloat hsml)
{
#ifdef POWER6
  long long numngb_old = numngb;
  long long *comp;
  long long no;
  if(!(comp = ALLOC_STACK(numngb_old * sizeof(long long))))
    {
      VERBOSE_ALL(0,
		  "Failed to allocate additional memory for `comp' (%d Mbytes), switch off 'DO_NOT_BRACH_IF'.\n",
		  numngb_old * sizeof(long long));
      endrun(124);
    }
  // mymalloc is nod thread save !!
  // comp = (long long *) mymalloc("NgbFilter", numngb_old * sizeof(long long));
#else
  int numngb_old = numngb;
  int *comp;
  int no;
  if(!(comp = ALLOC_STACK(numngb_old * sizeof(long long))))
    {
      VERBOSE_ALL(0,
		  "Failed to allocate additional memory for `comp' (%d Mbytes), switch off 'DO_NOT_BRACH_IF'.\n",
		  numngb_old * sizeof(long long));
      endrun(124);
    }
  // mymalloc is nod thread save !!
  // comp = (int *) mymalloc("NgbFilter", numngb_old * sizeof(int));
#endif

  // first compute all the distances

  numngb = 0;
  for(no = 0; no < numngb_old; no++)
    {
      int p = list[no];
      MyLongDouble dx, dy, dz, d2;
      MyLongDouble dist = DMAX(P[p].Hsml, hsml);

      dx = NGB_PERIODIC_LONG(P[p].Pos[0] - center->d[0], box->d[0], hbox->d[0]);
      dy = NGB_PERIODIC_LONG(P[p].Pos[1] - center->d[1], box->d[1], hbox->d[1]);
      dz = NGB_PERIODIC_LONG(P[p].Pos[2] - center->d[2], box->d[2], hbox->d[2]);
      d2 = dx * dx + dy * dy + dz * dz;
#ifdef POWER6
      comp[no] = (long long) __fsel(d2 - dist * dist, 0.0, 1.0);
#else
      comp[no] = (d2 < dist * dist);
#endif
    }

  // then filter out the distant particles
  if(numngb_old > 0)
    for(no = 0; no < numngb_old; no++)
      {
	if(comp[no])
	  list[numngb++] = list[no];
      }

  //  myfree(comp);

#ifdef NGB_DEBUG
  VERBOSE_ALL(0, "ngb_treefind_pairs: numngb before/after filter: %d / %d\n", numngb_old, numngb);
#endif
#
  return (int) numngb;
}

int ngb_filter_variables(long long numngb, int list[], t_vector *center, t_vector *box, t_vector *hbox,
			 MyFloat dist)
{
#ifdef POWER6
  long long numngb_old = numngb;
  long long *comp;
  long long no;
  if(!(comp = ALLOC_STACK(numngb_old * sizeof(long long))))
    {
      VERBOSE_ALL(0,
		  "Failed to allocate additional memory for `comp' (%d Mbytes), switch off 'DO_NOT_BRACH_IF'.\n",
		  numngb_old * sizeof(long long));
      endrun(124);
    }
  // mymalloc is nod thread save !!
  //  comp = (long long *) mymalloc("NgbFilter", numngb_old * sizeof(long long));
#else
  int numngb_old = numngb;
  int *comp;
  int no;
  if(!(comp = ALLOC_STACK(numngb_old * sizeof(long long))))
    {
      VERBOSE_ALL(0,
		  "Failed to allocate additional memory for `comp' (%d Mbytes), switch off 'DO_NOT_BRACH_IF'.\n",
		  numngb_old * sizeof(long long));
      endrun(124);
    }
  // mymalloc is nod thread save !!
  //  comp = (int *) mymalloc("NgbFilter", numngb_old * sizeof(int));
#endif

  // first compute all the distances

  numngb = 0;
  for(no = 0; no < numngb_old; no++)
    {
      int p = list[no];
      MyLongDouble dx, dy, dz, d2;

      dx = NGB_PERIODIC_LONG(P[p].Pos[0] - center->d[0], box->d[0], hbox->d[0]);
      dy = NGB_PERIODIC_LONG(P[p].Pos[1] - center->d[1], box->d[1], hbox->d[1]);
      dz = NGB_PERIODIC_LONG(P[p].Pos[2] - center->d[2], box->d[2], hbox->d[2]);
      d2 = dx * dx + dy * dy + dz * dz;
#ifdef POWER6
      comp[no] = (long long) __fsel(d2 - dist * dist, 0.0, 1.0);
#else
      comp[no] = (d2 < dist * dist);
#endif
    }

  // then filter out the distant particles
  if(numngb_old > 0)
    for(no = 0; no < numngb_old; no++)
      {
	if(comp[no])
	  list[numngb++] = list[no];
      }

  //  myfree(comp);

#ifdef NGB_DEBUG
  VERBOSE_ALL(0, "ngb_treefind_vars: numngb before/after filter: %d / %d\n", numngb_old, numngb);
#endif
#
  return numngb;
}


#endif // DO_NOT_BRACH_IF


/*! This routine finds all neighbours `j' that can interact with the
 *  particle `i' in the communication buffer.
 *
 *  Note that an interaction can take place if
 *  \f$ r_{ij} < h_i \f$  OR if  \f$ r_{ij} < h_j \f$.
 *
 *  In the range-search this is taken into account, i.e. it is guaranteed that
 *  all particles are found that fulfil this condition, including the (more
 *  difficult) second part of it. For this purpose, each node knows the
 *  maximum h occuring among the particles it represents.
 */
int ngb_treefind_pairs(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
		       int mode, int *nexport, int *nsend_local)
{
  int no, p, numngb, task, nexport_save;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
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

  nexport_save = *nexport;

  numngb = 0;

  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

#ifdef BLACK_HOLES
	  if(P[p].Mass == 0)
	    continue;
#endif

	  if(P[p].Ti_current != ti_Current)
	    drift_particle(p, ti_Current);

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(P[p].Hsml, hsml);

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
	  Ngblist[numngb++] = p;	/* Note: unlike in previous versions of the code, the buffer
					   can hold up to all particles */
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(23131);

	      if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= bunchSize)
		    {
		      *nexport = nexport_save;
		      if(nexport_save == 0)
			endrun(13003);	/* in this case, the buffer is too small to process even a single particle */
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
		DomainNodeIndex[no - (maxPart + maxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_pairs(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    force_drift_node(no, ti_Current);

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(Extnodes[no].hmax, hsml) + 0.5 * current->len;

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
	  no = ngb_check_node(current, &vcenter, &box, &hbox, DMAX(Extnodes[no].hmax, hsml));
#endif
	}
    }


  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_pairs(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}

int ngb_treefind_pairs_threads(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
			       int mode, int *exportflag, int *exportnodecount, int *exportindex,
			       int *ngblist)
{
  int no, p, numngb, task, nexp;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
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
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;
#ifdef BLACK_HOLES
	  if(P[p].Mass == 0)
	    continue;
#endif

	  if(P[p].Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      drift_particle(p, ti_Current);
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(P[p].Hsml, hsml);

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

	  ngblist[numngb++] = p;	/* Note: unlike in previous versions of the code, the buffer
					   can hold up to all particles */
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
#ifdef DONOTUSENODELIST
	      if(mode == 1)
		{
		  no = Nextnode[no - maxNodes];
		  continue;
		}
#endif
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      exportflag[task] = target;
		      exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(exportnodecount[task] == NODELISTLENGTH)
		    {
#pragma omp atomic capture
		      nexp = Nexport++;
		      if(nexp >= bunchSize)
			{
			  BufferFullFlag = 1;
			  Nexport = bunchSize;
			  return -1;
			}
		      exportnodecount[task] = 0;
		      exportindex[task] = nexp;
		      DataIndexTable[nexp].Task = task;
		      DataIndexTable[nexp].Index = target;
		      DataIndexTable[nexp].IndexGet = nexp;
		    }

#ifndef DONOTUSENODELIST
		  DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
#endif
		}

	      no = Nextnode[no - maxNodes];
	      continue;

	    }

	  current = &Nodes[no];

#ifndef DONOTUSENODELIST
	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_pairs(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }
#endif

	  if(current->Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
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

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(Extnodes[no].hmax, hsml) + 0.5 * current->len;

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
	  no = ngb_check_node(current, &vcenter, &box, &hbox, DMAX(Extnodes[no].hmax, hsml));
#endif
	}
    }


  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_pairs(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}





/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
#if !defined(KD_FRICTION)
int ngb_treefind_variable(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			  int *nexport, int *nsend_local)
#else
int ngb_treefind_variable(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			  int *nexport, int *nsend_local, int targettype)
#endif
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
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

  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

#if !defined(KD_FRICTION)
	  if(P[p].Type > 0)
	    continue;
#else
#ifdef KD_FRACTION
	  if(P[p].Type > 0 && targettype != 5)
	    continue;
#else
	  if(P[p].Type != targettype)
	    continue;
#endif
#endif

#ifdef BLACK_HOLES
	  if(P[p].Mass == 0)
	    continue;
#endif

	  if(P[p].Ti_current != ti_Current)
	    drift_particle(p, ti_Current);

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;

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

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= bunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13004);	/* in this case, the buffer is too small to process even a single particle */
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
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    force_drift_node(no, ti_Current);

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
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
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }


  //printf("%d\n", numngb);

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}

int ngb_treefind_fof_nearest_openmp(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				    int mode, int *exportflag, int *exportnodecount, int *exportindex,
				    int *ngblist, int MyFOF_PRIMARY_LINK_TYPES)
{
  int numngb, no, nexp, p, task;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;

#ifndef DO_NOT_BRACH_IF
  MyLongDouble dx, dy, dz, dist;
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

      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(!((1 << P[p].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
	    continue;

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
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
	  ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      exportflag[task] = target;
		      exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(exportnodecount[task] == NODELISTLENGTH)
		    {
#pragma omp atomic capture
		      nexp = Nexport++;
		      if(nexp >= bunchSize)
			{
			  BufferFullFlag = 1;
			  Nexport = bunchSize;
			  return -1;
			}
		      exportnodecount[task] = 0;
		      exportindex[task] = nexp;
		      DataIndexTable[nexp].Task = task;
		      DataIndexTable[nexp].Index = target;
		      DataIndexTable[nexp].IndexGet = nexp;
		    }

		  DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;

		}

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
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

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
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
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}


/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */

int ngb_treefind_variable_threads_noexport(MyLongDouble searchcenter[3], MyFloat hsml, int target,
					   int *startnode, int mode, int *ngblist, int debug)
{
  int numngb, no, nexp, p, task;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;

#ifndef DO_NOT_BRACH_IF
  MyLongDouble dx, dy, dz, dist;
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
      if(debug)
	VERBOSE_ALL(0, "ThisTask%d no=%d\n", ThisTask, no);

      if(no < maxPart)		/* single particle */
	{
	  if(debug)
	    VERBOSE_ALL(0, "ThisTask%d no<maxPart = %d\n", ThisTask, maxPart);
	  p = no;
	  no = Nextnode[no];


	  if(P[p].Type > 0)
	    continue;


#ifdef BLACK_HOLES
	  if(P[p].Mass == 0)
	    continue;
#endif

	  if(P[p].Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      drift_particle(p, ti_Current);
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
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
	  ngblist[numngb++] = p;
	  if(debug)
	    VERBOSE(0, "ThisTask%d dx = %f; dy = %f; dz = %f; [%f==%f] dist=%f; ID=%llu; numngb=%d\n",
		    ThisTask, dx, dy, dz, dx * dx + dy * dy + dz * dz, dist * dist, dist,
		    (unsigned long long) P[p].ID, numngb);
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(debug)
		VERBOSE_ALL(0, "ThisTask%d no >= maxPart + maxNode = %d + %d = %d\n", ThisTask, maxPart,
			    maxNodes, maxPart + maxNodes);
	      no = Nextnode[no - maxNodes];
	      continue;
	    }
	  if(debug)
	    VERBOSE_ALL(0, "ThisTask%d else\n", ThisTask);

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    {

#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      force_drift_node(no, ti_Current);
	    }

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  if(debug)
		    VERBOSE_ALL(0, "ThisTask%d no -> & continue%d\n", ThisTask, no);
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
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
	  if(debug)
	    VERBOSE_ALL
	      (0,
	       "ThisTask%d dx = %f; dy = %f; dz = %f; [%f==%f] dist=%f  = %f + .5*%f + %f * %f) ; numngb=%d; no-> %d\n",
	       ThisTask, dx, dy, dz, dx * dx + dy * dy + dz * dz, dist * dist, dist, hsml, current->len,
	       FACT1, current->len, numngb, no);
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}



/*

  This function will run a special neighbour search on a given position and hsml and ensure that witin r<hsml it finds a subset of the neighboursp rovided 
  in ngblist.
  
  Call it anywere you have doubts about your neighbours.

  AR.
*/
int test_ngb_inconsitency(int *ngblist, int numngb, int *check_ngblist, MyLongDouble *local_pos,
			  MyFloat local_hsml, int startnode, int mode, int target)
{

  int _startnode = startnode;
  int check_numngb_inbox = 0;
  MyAtLeastDouble h2 = local_hsml * local_hsml;

  //local_pos[0]+=0.01; //trigger a problem
  //local_hsml*=1.2;

  check_numngb_inbox =
    ngb_treefind_variable_threads_noexport(local_pos, local_hsml, target, &_startnode, mode, check_ngblist,
					   0);
  //VERBOSE_OUTPUT("%d vs %d",numngb, check_numngb_inbox);
  PANIC_IF((numngb < 0) && (check_numngb_inbox > 0), "(numngb_inbox<0)&&(check_numngb_inbox>0)");

  PANIC_IF(((check_numngb_inbox) > NumPart / DENSITY_LESS_NGB_FACTOR),
	   "density: Number of check neighbours found (%d) exceeds allowed number (%d) derived "
	   "from NumPart/DENSITY_LESS_NGB_FACTOR, to continue, decreas DENSITY_LESS_NGB_FACTOR "
	   "(currently %f, min possible val would be 1) !!\n", check_numngb_inbox,
	   (int) (NumPart / DENSITY_LESS_NGB_FACTOR), (double) DENSITY_LESS_NGB_FACTOR);

  PANIC_IF(check_numngb_inbox > numngb,
	   "P-Gadget NGB search found %d ngbs, while greentree found %d. Mode %d from"
	   "(POS=[%f, %f, %f], HSML=%f) on Task%d", check_numngb_inbox, numngb, mode, local_pos[0],
	   local_pos[1], local_pos[2], local_hsml, ThisTask);

  for(int _i = 0; _i < check_numngb_inbox; _i++)
    {
      int found = 0;
      int _ngb = check_ngblist[_i];
      MyLongDouble xtmp;
      MyAtLeastDouble dist, dx, dy, dz;
      dx = NGB_PERIODIC_LONG_X(local_pos[0] - P[_ngb].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(local_pos[1] - P[_ngb].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(local_pos[2] - P[_ngb].Pos[2]);

      if(dx * dx + dy * dy + dz * dz > h2)	//non-interesting ngb
	continue;
      for(int _j = 0; _j < numngb; _j++)
	{
	  if(_ngb == ngblist[_j])
	    {
	      found = 1;
	      break;
	    }
	}
      if(found == 0)
	{

	  int _startnode = startnode;
	  check_numngb_inbox =
	    ngb_treefind_variable_threads_noexport(local_pos, local_hsml, target, &_startnode, mode,
						   check_ngblist, 1);



	  VERBOSE_ALL(2, "P[%d] (ID=%llu; POS=[%f, %f, %f]) was found as ngb with P-Gadget3 "
		      "NGB search but it was not found with GREEN_NTREE from on phase %d from"
		      "(POS=[%f, %f, %f], HSML=%f) on Task%d; ngblist=%p; check_numngb_inbox=%d, check_ngblist=%p;  numngb=%d;",
		      _ngb, (unsigned long long) P[_ngb].ID, P[_ngb].Pos[0], P[_ngb].Pos[1], P[_ngb].Pos[2],
		      mode,
		      local_pos[0], local_pos[1], local_pos[2], local_hsml,
		      ThisTask, ngblist, check_numngb_inbox, check_ngblist, numngb);
	  return 0;
	}
    }
  return 1;
}



/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
#if !defined(KD_FRICTION)
int ngb_treefind_variable_threads(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				  int mode, int *exportflag, int *exportnodecount, int *exportindex,
				  int *ngblist)
#else
int ngb_treefind_variable_threads(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				  int mode, int *exportflag, int *exportnodecount, int *exportindex,
				  int *ngblist, int targettype)
#endif
{
  int numngb, no, nexp, p, task;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;

#ifndef DO_NOT_BRACH_IF
  MyLongDouble dx, dy, dz, dist;
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

      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

#if !defined(KD_FRICTION)
	  if(P[p].Type > 0)
	    continue;
#else
#ifdef KD_FRACTION
	  if(P[p].Type > 0 && targettype != 5)
	    continue;
#else
	  if(P[p].Type != targettype)
	    continue;
#endif
#endif

#ifdef BLACK_HOLES
	  if(P[p].Mass == 0)
	    continue;
#endif

	  if(P[p].Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      drift_particle(p, ti_Current);
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
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
	  ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      exportflag[task] = target;
		      exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(exportnodecount[task] == NODELISTLENGTH)
		    {
#pragma omp atomic capture
		      nexp = Nexport++;
		      if(nexp >= bunchSize)
			{
			  BufferFullFlag = 1;
			  Nexport = bunchSize;
			  return -1;
			}
		      exportnodecount[task] = 0;
		      exportindex[task] = nexp;
		      DataIndexTable[nexp].Task = task;
		      DataIndexTable[nexp].Index = target;
		      DataIndexTable[nexp].IndexGet = nexp;
		    }

		  DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;

		}

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
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

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
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
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}






/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
#define n_last_nodes 20

int ngb_treefind_variable_threads_check(MyLongDouble searchcenter[3], MyFloat hsml, int target,
					int *startnode, int mode)
{
  int numngb, no, nexp, p, task;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;

#ifndef DO_NOT_BRACH_IF
  MyLongDouble dx, dy, dz, dist;
#ifdef PERIODIC
  MyLongDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif
  int last_nodes[n_last_nodes];
  int i_last_node = 0;
  for(int i = 0; i < n_last_nodes; i++)
    {
      last_nodes[i] = -1;
    }
  numngb = 0;
  no = *startnode;
  while(no >= 0)
    {

      for(int i = 0; i < n_last_nodes; i++)
	{
	  if(last_nodes[i] != -1 && last_nodes[i] == no)
	    {
	      for(int j = 0; j < n_last_nodes; j++)
		VERBOSE_ALL(5, " %d ", last_nodes[j]);
	      PANIC_ANY
		("Tree is broken. no=%d and i=%d i_last_node=%d target=%d search_center=[%f %f %f] hsml=%f P[target].ID=%llu Type=%d numngb=0",
		 no, i, i_last_node, target, searchcenter[0], searchcenter[1], searchcenter[2], hsml,
		 (unsigned long long) P[target].ID, P[target].Type);
	    }
	}
      last_nodes[i_last_node] = no;
      i_last_node = (i_last_node + 1) % n_last_nodes;


      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

#if !defined(KD_FRICTION)
	  if(P[p].Type > 0)
	    continue;
#else
	  /*
	     #ifdef KD_FRACTION
	     if(P[p].Type > 0 && targettype != 5)
	     continue;
	     #else
	     if(P[p].Type != targettype)
	     continue;
	     #endif
	   */
#endif

#ifdef BLACK_HOLES
	  if(P[p].Mass == 0)
	    continue;
#endif

	  if(P[p].Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      drift_particle(p, ti_Current);
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
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
	  numngb++;
	  //      ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);


	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    {
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
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

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
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
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}





/*! Allocates memory for the neighbour list buffer.
 */
void ngb_init(void)
{

}


/*! This function constructs the neighbour tree. To this end, we actually need
 *  to construct the gravitational tree, because we use it now for the
 *  neighbour search.
 */
void ngb_treebuild(void)
{
  VERBOSE(3, "Begin Ngb-tree construction.\n");

  CPU_Step[CPU_MISC] += measure_time();

  force_treebuild(NumPart, NULL);

  CPU_Step[CPU_TREEBUILD] += measure_time();

  VERBOSE(4, "Ngb-Tree contruction finished \n");
}



int ngb_treefind_fof_primary(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local, int MyFOF_PRIMARY_LINK_TYPES)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  MyLongDouble dx, dy, dz, dist, r2;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;

#ifndef DO_NOT_BRACH_IF
#ifdef PERIODIC
  MyLongDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(!((1 << P[p].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
	    continue;

	  if(mode == 0)
	    continue;

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
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
	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(mode == 0)
		{
		  if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= bunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13005);	/* in this case, the buffer is too small to process even a single particle */
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
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      if(mode == -1)
		{
		  *nexport = 1;
		}

	      no = Nextnode[no - maxNodes];
	      continue;

	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(mode == 0)
	    {
	      if(!(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL)))	/* we have a node with only local particles, can skip branch */
		{
		  no = current->u.d.sibling;
		  continue;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
#ifndef DO_NOT_BRACH_IF
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
#else
	  dx = NGB_PERIODIC_LONG(current->center[0] - searchcenter[0], box.d[0], hbox.d[0]);
#endif
	  if(dx > dist)
	    continue;
#ifndef DO_NOT_BRACH_IF
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
#else
	  dy = NGB_PERIODIC_LONG(current->center[1] - searchcenter[1], box.d[1], hbox.d[1]);
#endif
	  if(dy > dist)
	    continue;
#ifndef DO_NOT_BRACH_IF
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
#else
	  dz = NGB_PERIODIC_LONG(current->center[2] - searchcenter[2], box.d[2], hbox.d[2]);
#endif
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
	    continue;

	  if((current->u.d.bitflags & ((1 << BITFLAG_TOPLEVEL) + (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))) == 0)	/* only use fully local nodes */
	    {
	      /* test whether the node is contained within the sphere */
	      dist = hsml - FACT2 * current->len;
	      if(dist > 0)
		if(r2 < dist * dist)
		  {
		    if(current->u.d.bitflags & (1 << BITFLAG_INSIDE_LINKINGLENGTH))	/* already flagged */
		      {
			/* sufficient to return only one particle inside this cell */

			p = current->u.d.nextnode;
			while(p >= 0)
			  {
			    if(p < maxPart)
			      {
				if(((1 << P[p].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
				  {
#ifndef DO_NOT_BRACH_IF
				    dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
				    dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
				    dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
				    if(dx * dx + dy * dy + dz * dz > hsml * hsml)
				      break;
#endif
				    Ngblist[numngb++] = p;
				    break;
				  }
				p = Nextnode[p];
			      }
			    else if(p >= maxPart + maxNodes)
			      p = Nextnode[p - maxNodes];
			    else
			      p = Nodes[p].u.d.nextnode;
			  }
			continue;
		      }
		    else
		      {
			/* flag it now */
			current->u.d.bitflags |= (1 << BITFLAG_INSIDE_LINKINGLENGTH);
		      }
		  }
	    }

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}



/* find all particles of type FOF_PRIMARY_LINK_TYPES in smoothing length in order to find nearest one */
int ngb_treefind_fof_nearest(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local, int MyFOF_PRIMARY_LINK_TYPES)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
#ifndef DO_NOT_BRACH_IF
  MyLongDouble dx, dy, dz, dist;
#ifdef PERIODIC
  MyLongDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

#define FACT2 0.86602540

  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(!((1 << P[p].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
	    continue;

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
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

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(123192);

	      if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= bunchSize)
		    {
		      *nexport = nexport_save;
		      if(nexport_save == 0)
			endrun(13005);	/* in this case, the buffer is too small to process even a single particle */
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
		DomainNodeIndex[no - (maxPart + maxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
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
	  if((dx * dx + dy * dy + dz * dz) > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}











#ifdef ADAPTGRAVSOFT
#include"AdaptGravSoft/ads_ngb.c"
#endif
