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


#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"



/*  This function computes the gravitational forces for all active particles.
 *  A new tree is constructed, if the number of force computations since
 *  it's last construction exceeds some fraction of the total
 *  particle number, otherwise tree nodes are dynamically updated if needed.
 */

#ifdef FORCETEST

void gravity_forcetest(void)
{
  int iter = 0, nthis, task, nloc, ntot;
  double tstart, tend, timetree = 0;
  int i, j, ndone, ndone_flag, ngrp, place;

#ifndef NOGRAVITY
  int k, nexport, nimport;
  int sendTask, recvTask;
  double fac1;
  MPI_Status status;
#endif
  double costtotal, *costtreelist;
  double maxt, sumt, *timetreelist;
  double fac;
  char buf[200];

#ifdef PMGRID
  if(All.PM_Ti_endstep != All.Ti_Current)
    return;
#endif

  if(All.ComovingIntegrationOn)
    set_softenings();		/* set new softening lengths */

  for(i = FirstActiveParticle, nloc = 0; i >= 0; i = NextActiveParticle[i])
    {
      if(get_random_number(P[i].ID) < FORCETEST)
	{
	  P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as selected */
	  nloc++;
	}
    }

  MPI_Allreduce(&nloc, &ntot, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct gravdata_in) + sizeof(struct gravdata_out) +
					     sizemax(sizeof(struct gravdata_in),
						     sizeof(struct gravdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  costtotal = 0;


  i = FirstActiveParticle;	/* beginn with this index */

  do
    {
      iter++;

      for(j = 0; j < NTask; j++)
	Send_count[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].TimeBin < 0)
	  {
	    if(nexport + NTask - 1 <= All.BunchSize)
	      {
		for(task = 0; task < NTask; task++)
		  if(task != ThisTask)
		    {
		      DataIndexTable[nexport].Task = task;
		      DataIndexTable[nexport].Index = i;
		      DataIndexTable[nexport].IndexGet = nexport;
		      nexport++;
		      Send_count[task]++;
		    }

		costtotal += force_treeevaluate_direct(i, 0);
	      }
	    else
	      break;
	  }
      tend = second();
      timetree += timediff(tstart, tend);

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif
      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      GravDataGet = (struct gravdata_in *) mymalloc("GravDataGet", nimport * sizeof(struct gravdata_in));
      GravDataIn = (struct gravdata_in *) mymalloc("GravDataIn", nexport * sizeof(struct gravdata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    GravDataIn[j].Pos[k] = P[place].Pos[k];

	  GravDataIn[j].Type = P[place].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(P[place].Type == 0)
	    GravDataIn[j].Soft = SphP[place].Hsml;
#endif
	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
			       recvTask, TAG_GRAV_A,
			       &GravDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
			       recvTask, TAG_GRAV_A, MYMPI_COMM_WORLD, &status);
		}
	    }
	}

      myfree(GravDataIn);
      GravDataResult =
	(struct gravdata_out *) mymalloc("GravDataResult", nimport * sizeof(struct gravdata_out));
      GravDataOut = (struct gravdata_out *) mymalloc("GravDataOut", nexport * sizeof(struct gravdata_out));

      for(j = 0; j < nimport; j++)
	costtotal += force_treeevaluate_direct(j, 1);


      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

      /* get the result */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&GravDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct gravdata_out),
			       MPI_BYTE, recvTask, TAG_GRAV_B,
			       &GravDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct gravdata_out),
			       MPI_BYTE, recvTask, TAG_GRAV_B, MYMPI_COMM_WORLD, &status);
		}
	    }
	}

      /* add the results to the local particles */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    P[place].GravAccelDirect[k] += GravDataOut[j].Acc[k];
	}

      myfree(GravDataOut);
      myfree(GravDataResult);
      myfree(GravDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);


  /* now add things for comoving integration */

  if(All.ComovingIntegrationOn)
    {
#ifndef PERIODIC
      fac1 = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(P[i].TimeBin < 0)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccelDirect[j] += fac1 * P[i].Pos[j];
#endif
    }


  /*  muliply by G */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].TimeBin < 0)
      for(j = 0; j < 3; j++)
	P[i].GravAccelDirect[j] *= All.G;



  /* Finally, the following factor allows a computation of cosmological simulation
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      fac1 = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(P[i].TimeBin < 0)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccelDirect[j] += fac1 * P[i].Pos[j];
    }

  /* now output the forces to a file */

  for(nthis = 0; nthis < NTask; nthis++)
    {
      if(nthis == ThisTask)
	{
	  sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
	  if(!(FdForceTest = fopen(buf, "a")))
	    {
	      VERBOSE_ALL(0, "error in opening file '%s'\n", buf);
	      endrun(17);
	    }
	  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	    if(P[i].TimeBin < 0)
	      {
#ifndef PMGRID
		fprintf(FdForceTest, "%d %g %g %g %g %g %g %g %g %g %g %g\n",
			P[i].Type, All.Time, All.Time - TimeOfLastTreeConstruction,
			P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
			P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2],
			P[i].g.GravAccel[0], P[i].g.GravAccel[1], P[i].g.GravAccel[2]);
#else
		fprintf(FdForceTest, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
			P[i].Type, All.Time, All.Time - TimeOfLastTreeConstruction,
			P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
			P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2],
			P[i].g.GravAccel[0], P[i].g.GravAccel[1], P[i].g.GravAccel[2],
			P[i].GravPM[0] + P[i].g.GravAccel[0],
			P[i].GravPM[1] + P[i].g.GravAccel[1], P[i].GravPM[2] + P[i].g.GravAccel[2]);
#endif
	      }
	  fclose(FdForceTest);
	}
      MPI_Barrier(MYMPI_COMM_WORLD);
    }

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].TimeBin < 0)
      P[i].TimeBin = -P[i].TimeBin - 1;

  /* Now the force computation is finished */



  timetreelist = mymalloc("timetreelist", sizeof(double) * NTask);
  costtreelist = mymalloc("costtreelist", sizeof(double) * NTask);

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MYMPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, costtotal = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];
	}

      fprintf(FdTimings, "DIRECT Nf= %d    part/sec=%g | %g  ia/part=%g \n", ntot, ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);
    }

  myfree(costtreelist);
  myfree(timetreelist);
}

#endif
