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


void force_update_tree(void)
{

  double t0 = second();
  //  PUSH_T("force update tree");

  int i, j;

  VERBOSE(2, "kicks will prepare for dynamic update of tree\n");

  CPU_Step[CPU_MISC] += measure_time();


  GlobFlag++;
  DomainNumChanged = 0;
  DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));


  /* note: the current list of active particles still refers to that
   * synchronized at the previous time.
   */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      force_kick_node(i, P[i].dp);	/* kick the parent nodes with this momentum
					   difference, also updated maximum velocity, softening and soundspeed, if needed */
      for(j = 0; j < 3; j++)
	P[i].dp[j] = 0;
    }

  force_finish_kick_nodes();
  myfree(DomainList);

  CPU_Step[CPU_TREEUPDATE] += measure_time();

  VERBOSE(2, "Tree has been updated dynamically.\n");


  //  POP_T();
  double t1 = second();

#ifdef KD_EXTRA_TIMER_OUTPUT
  VERBOSE(2, "EXTRA TIMER: force update took %g sec\n", timediff(t0, t1));
#endif
}


void force_kick_node(int i, MyFloat *dp)
{
  int j, no;
  MyFloat v, vmax;

#ifdef NEUTRINOS
  if((P[i].Type == 2) && (All.Time <= All.Time_tree_on_nu))
    {
      VERBOSE(0, "Situation 1\n");
      return;
    }
#endif

  for(j = 0, vmax = 0; j < 3; j++)
    if((v = fabs(P[i].Vel[j])) > vmax)
      vmax = v;

  no = Father[i];

  PANIC_IF(no < All.MaxPart, "error: at i=%d, no=%d < All.MaxPart=%d\n", no, i, All.MaxPart);

  while(no >= 0)
    {
      force_drift_node(no, All.Ti_Current);

      for(j = 0; j < 3; j++)
	{
	  Extnodes[no].dp[j] += dp[j];
	}

      if(Extnodes[no].vmax < vmax)
	Extnodes[no].vmax = vmax;

      Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);

      Extnodes[no].Ti_lastkicked = All.Ti_Current;

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* top-level tree-node reached */
	{
	  if(Extnodes[no].Flag != GlobFlag)
	    {
	      Extnodes[no].Flag = GlobFlag;
	      DomainList[DomainNumChanged++] = no;
	    }
	  break;
	}

      PANIC_IF(no == Nodes[no].u.d.father,
	       "Task=%d Fixed point  no == Nodes[no].u.d.father  detected, no=%d\n", ThisTask, no);
      no = Nodes[no].u.d.father;
    }
}






void force_finish_kick_nodes(void)
{
  int i, j, no, ta, totDomainNumChanged;
  int *domainList_all;
  int *counts, *counts_dp, *offset_list, *offset_dp, *offset_vmax;
  MyLongDouble *domainDp_loc, *domainDp_all;

  MyFloat *domainVmax_loc, *domainVmax_all;

  /* share the momentum-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  counts_dp = (int *) mymalloc("counts_dp", sizeof(int) * NTask);
  offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
  offset_dp = (int *) mymalloc("offset_dp", sizeof(int) * NTask);
  offset_vmax = (int *) mymalloc("offset_vmax", sizeof(int) * NTask);

  domainDp_loc = (MyLongDouble *) mymalloc("domainDp_loc", DomainNumChanged * 3 * sizeof(MyLongDouble));

  domainVmax_loc = (MyFloat *) mymalloc("domainVmax_loc", DomainNumChanged * sizeof(MyFloat));

  for(i = 0; i < DomainNumChanged; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  domainDp_loc[i * 3 + j] = Extnodes[DomainList[i]].dp[j];
	}
      domainVmax_loc[i] = Extnodes[DomainList[i]].vmax;
    }

  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MYMPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_dp[0] = 0, offset_vmax[0] = 0; ta < NTask;
      ta++)
    {
      totDomainNumChanged += counts[ta];
      if(ta > 0)
	{
	  offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
	  offset_dp[ta] = offset_dp[ta - 1] + counts[ta - 1] * 3 * sizeof(MyLongDouble);
	  offset_vmax[ta] = offset_vmax[ta - 1] + counts[ta - 1] * sizeof(MyFloat);
	}
    }

  VERBOSE(2, "I exchange kick momenta for %d top-level nodes out of %d\n", totDomainNumChanged, NTopleaves);


  domainDp_all = (MyLongDouble *) mymalloc("domainDp_all", totDomainNumChanged * 3 * sizeof(MyLongDouble));

  domainVmax_all = (MyFloat *) mymalloc("domainVmax_all", totDomainNumChanged * sizeof(MyFloat));

  domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

  MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
		 domainList_all, counts, offset_list, MPI_INT, MYMPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    {
      counts_dp[ta] = counts[ta] * 3 * sizeof(MyLongDouble);
      counts[ta] *= sizeof(MyFloat);
    }


  MPI_Allgatherv(domainDp_loc, DomainNumChanged * 3 * sizeof(MyLongDouble), MPI_BYTE,
		 domainDp_all, counts_dp, offset_dp, MPI_BYTE, MYMPI_COMM_WORLD);

  MPI_Allgatherv(domainVmax_loc, DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
		 domainVmax_all, counts, offset_vmax, MPI_BYTE, MYMPI_COMM_WORLD);


  /* construct momentum kicks in top-level tree */
  for(i = 0; i < totDomainNumChanged; i++)
    {
      no = domainList_all[i];

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))	/* to avoid that the local one is kicked twice */
	no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  force_drift_node(no, All.Ti_Current);

	  for(j = 0; j < 3; j++)
	    {
	      Extnodes[no].dp[j] += domainDp_all[3 * i + j];
	    }

	  if(Extnodes[no].vmax < domainVmax_all[i])
	    Extnodes[no].vmax = domainVmax_all[i];

	  Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);
	  Extnodes[no].Ti_lastkicked = All.Ti_Current;

	  no = Nodes[no].u.d.father;
	}
    }

  myfree(domainList_all);
  myfree(domainVmax_all);
  myfree(domainDp_all);
  myfree(domainVmax_loc);
  myfree(domainDp_loc);
  myfree(offset_vmax);
  myfree(offset_dp);
  myfree(offset_list);
  myfree(counts_dp);
  myfree(counts);
}



void force_drift_node(int no, int time1)
{
  int j, time0;
  double dt_drift, dt_drift_hmax, fac;
#ifdef AR_XMAS_TREE		//xmas tree has holes
  if(Nodes[no].allocated != 1)
    return;
#endif
  if(time1 == Nodes[no].Ti_current)
    return;

  time0 = Extnodes[no].Ti_lastkicked;

  if(Nodes[no].u.d.bitflags & (1 << BITFLAG_NODEHASBEENKICKED))
    {
      if(Extnodes[no].Ti_lastkicked != Nodes[no].Ti_current)
	{
	  VERBOSE(0, "Task=%d Extnodes[no].Ti_lastkicked=%d  Nodes[no].Ti_current=%d\n",
		  ThisTask, (int) Extnodes[no].Ti_lastkicked, (int) Nodes[no].Ti_current);
	  PANIC("inconsistency in drift node");
	}

      if(Nodes[no].u.d.mass)
	fac = 1 / Nodes[no].u.d.mass;
      else
	fac = 0;

      for(j = 0; j < 3; j++)
	{
	  Extnodes[no].vs[j] += fac * FLT(Extnodes[no].dp[j]);
	  Extnodes[no].dp[j] = 0;
	}
      Nodes[no].u.d.bitflags &= (~(1 << BITFLAG_NODEHASBEENKICKED));
    }

  if(All.ComovingIntegrationOn)
    {
      dt_drift_hmax = get_drift_factor(Nodes[no].Ti_current, time1);
      dt_drift = dt_drift_hmax;
    }
  else
    {

      dt_drift_hmax = (time1 - Nodes[no].Ti_current) * All.Timebase_interval;
      dt_drift = dt_drift_hmax;
    }

  for(j = 0; j < 3; j++)
    Nodes[no].u.d.s[j] += Extnodes[no].vs[j] * dt_drift;
  Nodes[no].len += 2 * Extnodes[no].vmax * dt_drift;

  Nodes[no].Ti_current = time1;
}





/*! This function updates the hmax-values in tree nodes that hold SPH
 *  particles. These values are needed to find all neighbors in the
 *  hydro-force computation.  Since the Hsml-values are potentially changed
 *  in the SPH-denity computation, force_update_hmax() should be carried
 *  out just before the hydrodynamical SPH forces are computed, i.e. after
 *  density().
 */
void force_update_hmax(void)
{
  int i, no, ta, totDomainNumChanged;
  int *domainList_all;
  int *counts, *offset_list, *offset_hmax;
  MyFloat *domainHmax_loc, *domainHmax_all;
  int OffsetSIZE = 1;
  integertime ti_step, time0, time1;
  double dt_drift_hmax;
  MyFloat hmax_est;

  GlobFlag++;

  DomainNumChanged = 0;
  DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
#ifdef BLACK_HOLES
    if(P[i].Type == 0 && P[i].Mass > 0)
#else
    if(P[i].Type == 0)
#endif
      {
	no = Father[i];

	while(no >= 0)
	  {
	    force_drift_node(no, All.Ti_Current);
#ifdef KD_FIX_HMAXUPDATE
	    time0 = P[i].Ti_current;
	    time1 = All.Ti_Current;
	    if(time0 != time1)
	      {
		dt_drift_hmax = get_drift_factor(time0, time1);
		hmax_est = P[i].Hsml * exp(0.333333333333 * SphP[i].DivVel * dt_drift_hmax);
	      }
	    else
	      hmax_est = P[i].Hsml;
#else
	    ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
	    time0 = All.Ti_Current;
	    time1 = P[i].Ti_begstep + ti_step;
	    if(All.ComovingIntegrationOn)
	      dt_drift_hmax = get_drift_factor(time0, time1);
	    else
	      dt_drift_hmax = (time1 - time0) * All.Timebase_interval;
	    hmax_est = P[i].Hsml * exp(0.333333333333 * SphP[i].DivVel * dt_drift_hmax);
#endif
	    if(P[i].Hsml > Extnodes[no].hmax || hmax_est > Extnodes[no].hmax)
	      {
		if(P[i].Hsml > Extnodes[no].hmax)
		  Extnodes[no].hmax = P[i].Hsml;

		if(hmax_est > Extnodes[no].hmax)
		  Extnodes[no].hmax = hmax_est;
		if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node */
		  {
		    if(Extnodes[no].Flag != GlobFlag)
		      {
			Extnodes[no].Flag = GlobFlag;
			DomainList[DomainNumChanged++] = no;
		      }
		    break;
		  }
	      }
	    else
	      break;

	    no = Nodes[no].u.d.father;
	  }
      }

  /* share the hmax-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
  offset_hmax = (int *) mymalloc("offset_hmax", sizeof(int) * NTask);

  domainHmax_loc = (MyFloat *) mymalloc("domainHmax_loc", DomainNumChanged * OffsetSIZE * sizeof(MyFloat));

  for(i = 0; i < DomainNumChanged; i++)
    {
      domainHmax_loc[OffsetSIZE * i] = Extnodes[DomainList[i]].hmax;
    }


  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MYMPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_hmax[0] = 0; ta < NTask; ta++)
    {
      totDomainNumChanged += counts[ta];
      if(ta > 0)
	{
	  offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
	  offset_hmax[ta] = offset_hmax[ta - 1] + counts[ta - 1] * OffsetSIZE * sizeof(MyFloat);
	}
    }

  VERBOSE(2, "Hmax exchange: %d topleaves out of %d\n", totDomainNumChanged, NTopleaves);

  domainHmax_all = (MyFloat *) mymalloc("domainHmax_all", totDomainNumChanged * OffsetSIZE * sizeof(MyFloat));
  domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

  MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
		 domainList_all, counts, offset_list, MPI_INT, MYMPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    counts[ta] *= OffsetSIZE * sizeof(MyFloat);

  MPI_Allgatherv(domainHmax_loc, OffsetSIZE * DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
		 domainHmax_all, counts, offset_hmax, MPI_BYTE, MYMPI_COMM_WORLD);


  for(i = 0; i < totDomainNumChanged; i++)
    {
      no = domainList_all[i];

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))	/* to avoid that the hmax is updated twice */
	no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  force_drift_node(no, All.Ti_Current);


	  if(domainHmax_all[OffsetSIZE * i] > Extnodes[no].hmax)
	    Extnodes[no].hmax = domainHmax_all[OffsetSIZE * i];
	  else
	    break;

	  no = Nodes[no].u.d.father;
	}
    }


  myfree(domainList_all);
  myfree(domainHmax_all);
  myfree(domainHmax_loc);
  myfree(offset_hmax);
  myfree(offset_list);
  myfree(counts);
  myfree(DomainList);

  CPU_Step[CPU_TREEHMAXUPDATE] += measure_time();
}
