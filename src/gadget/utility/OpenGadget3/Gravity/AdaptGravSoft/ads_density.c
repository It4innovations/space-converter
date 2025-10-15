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

void ags_density(void)
{
  MyFloat *Left, *Right;
  int i, j, k, ndone, ndone_flag, npleft, dt_step, iter = 0;
  int ngrp, sendTask, recvTask, place;
  long long ntot;
  double dmax1, dmax2, fac;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double dt_entr, tstart, tend, t0, t1;
  double desnumngb;
  int save_NextParticle;
  long long n_exported = 0;
  double maxsoft, minsoft;

  CPU_Step[CPU_AGSDENSMISC] += measure_time();

  int NTaskTimesNumPart;

  NTaskTimesNumPart = NumPart;
/*****************************************check******************************************/
#ifdef NUM_THREADS
  NTaskTimesNumPart = NUM_THREADS * NumPart;
#endif

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

#ifndef AGS_UPDATEALLPARTICLES
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
#else
  for(i = 0; i < NumPart; i++)
#endif
    {
      if(ags_density_isactive(i))
	{
	  Left[i] = Right[i] = 0;

	}
    }

/* allocate buffers to arrange communication */

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct AGS_densdata_in) +
					     sizeof(struct AGS_densdata_out) +
					     sizemax(sizeof(struct AGS_densdata_in),
						     sizeof(struct AGS_densdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  t0 = second();

  desnumngb = All.AGS_DesNumNgb;

/* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {

#ifndef AGS_UPDATEALLPARTICLES
      NextParticle = FirstActiveParticle;	/* begin with this index */
#else
      NextParticle = 0;
#endif

      do
	{
	  BufferFullFlag = 0;
	  Nexport = 0;
	  save_NextParticle = NextParticle;
	  tstart = second();

#ifdef NUM_THREADS
	  pthread_t mythreads[NUM_THREADS - 1];

	  int threadid[NUM_THREADS - 1];

	  pthread_attr_t attr;

	  pthread_attr_init(&attr);
	  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	  pthread_mutex_init(&mutex_nexport, NULL);
	  pthread_mutex_init(&mutex_partnodedrift, NULL);

	  TimerFlag = 0;

	  for(j = 0; j < NUM_THREADS - 1; j++)
	    {
	      threadid[j] = j + 1;
	      pthread_create(&mythreads[j], &attr, ags_density_evaluate_primary, &threadid[j]);
	    }
#endif
	  int mainthreadid = 0;

	  ags_density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */

#ifdef NUM_THREADS
	  for(j = 0; j < NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);
#endif

	  tend = second();
	  timecomp1 += timediff(tstart, tend);

	  if(BufferFullFlag)
	    {
	      int last_nextparticle = NextParticle;

	      NextParticle = save_NextParticle;

#ifndef AGS_UPDATEALLPARTICLES
	      while(NextParticle >= 0)
		{
		  if(NextParticle == last_nextparticle)
		    break;

		  if(ProcessedFlag[NextParticle] != 1)
		    break;

		  ProcessedFlag[NextParticle] = 2;

		  NextParticle = NextActiveParticle[NextParticle];
		}
#else
	      while(NextParticle < NumPart)
		{
		  if(NextParticle == last_nextparticle)
		    break;

		  if(ProcessedFlag[NextParticle] != 1)
		    break;

		  ProcessedFlag[NextParticle] = 2;

		  NextParticle = NextParticle + 1;
		}
#endif

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

#ifdef MYSORT
	  mysort_dataindex(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
#else
	  qsort(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
#endif

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

	  AGS_DensDataGet =
	    (struct AGS_densdata_in *) mymalloc("AGS_DensDataGet", Nimport * sizeof(struct AGS_densdata_in));
	  AGS_DensDataIn =
	    (struct AGS_densdata_in *) mymalloc("AGS_DensDataIn", Nexport * sizeof(struct AGS_densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      AGS_DensDataIn[j].Pos[0] = P[place].Pos[0];
	      AGS_DensDataIn[j].Pos[1] = P[place].Pos[1];
	      AGS_DensDataIn[j].Pos[2] = P[place].Pos[2];
	      AGS_DensDataIn[j].AGS_Hsml = P[place].AGS_Hsml;

	      memcpy(AGS_DensDataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));


	    }
	  /* exchange particle data */
	  tstart = second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&AGS_DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct AGS_densdata_in), MPI_BYTE,
				   recvTask, TAG_AGS_DENS_A,
				   &AGS_DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct AGS_densdata_in), MPI_BYTE,
				   recvTask, TAG_AGS_DENS_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(AGS_DensDataIn);
	  AGS_DensDataResult =
	    (struct AGS_densdata_out *) mymalloc("AGS_DensDataResult",
						 Nimport * sizeof(struct AGS_densdata_out));
	  AGS_DensDataOut =
	    (struct AGS_densdata_out *) mymalloc("AGS_DensDataOut",
						 Nexport * sizeof(struct AGS_densdata_out));
/*****************************************check************************************************/
	  report_memory_usage(&HighMark_agsdensity, "AGS_DENSITY");

	  /* now do the particles that were sent to us */

	  tstart = second();

	  NextJ = 0;

#ifdef NUM_THREADS
	  for(j = 0; j < NUM_THREADS - 1; j++)
	    pthread_create(&mythreads[j], &attr, ags_density_evaluate_secondary, &threadid[j]);
#endif
	  ags_density_evaluate_secondary(&mainthreadid);

#ifdef NUM_THREADS
	  for(j = 0; j < NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);

	  pthread_mutex_destroy(&mutex_partnodedrift);
	  pthread_mutex_destroy(&mutex_nexport);
	  pthread_attr_destroy(&attr);
#endif

	  tend = second();
	  timecomp2 += timediff(tstart, tend);

#ifndef AGS_UPDATEALLPARTICLES
	  if(NextParticle < 0)
#else
	  if(NextParticle == NumPart)
#endif
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
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&AGS_DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct AGS_densdata_out),
				   MPI_BYTE, recvTask, TAG_AGS_DENS_B,
				   &AGS_DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct AGS_densdata_out),
				   MPI_BYTE, recvTask, TAG_AGS_DENS_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
#ifdef AGS_OUTPUTNGBS
	      P[place].AGS_REALNumNgb += AGS_DensDataOut[j].AGS_REALNgb;
#endif
	      P[place].AGS_NumNgb += AGS_DensDataOut[j].AGS_Ngb;
	      P[place].AGS_Density += AGS_DensDataOut[j].AGS_Rho;
	      P[place].AGS_zeta += AGS_DensDataOut[j].AGS_zeta;
	      P[place].AGS_omega += AGS_DensDataOut[j].AGS_omega;

	    }
	  tend = second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(AGS_DensDataOut);
	  myfree(AGS_DensDataResult);
	  myfree(AGS_DensDataGet);
	}
      while(ndone < NTask);


      /* do final operations on results */
      tstart = second();
#ifndef AGS_UPDATEALLPARTICLES
      for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
#else
      for(i = 0, npleft = 0; i < NumPart; i++)
#endif
	{
	  if(ags_density_isactive(i))
	    {

	      if(P[i].AGS_Density > 0)
		{

		  P[i].AGS_zeta *= P[i].AGS_Hsml / (NUMDIMS * P[i].AGS_Density);
		  P[i].AGS_zeta *= -1;

		  P[i].AGS_omega *= P[i].AGS_Hsml / (NUMDIMS * P[i].AGS_Density);

		  if(P[i].AGS_omega > -0.9)
		    P[i].AGS_omega = 1 / (P[i].AGS_omega + 1);
		  else
		    P[i].AGS_omega = 1;

		}



#ifndef WAKEUP
	      dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
#else
	      dt_step = P[i].dt_step;
#endif
	      dt_entr = (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;



	      /* now check whether we had enough neighbours */

	      desnumngb = All.AGS_DesNumNgb;
#ifdef PMGRID
#ifdef PLACEHIGHRESREGION
	      if(pmforce_is_particle_high_res(P[i].Type, P[i].Pos))
		{
		  minsoft = All.AGS_MinSoft[1];
		  maxsoft = All.AGS_MaxSoft[1];
		}
	      else
		{
		  minsoft = All.AGS_MinSoft[0];
		  maxsoft = All.AGS_MaxSoft[0];
		}
#else
	      minsoft = All.AGS_MinSoft[0];
	      maxsoft = All.AGS_MaxSoft[0];
#endif
	      if((P[i].AGS_NumNgb < (desnumngb - All.AGS_MaxNumNgbDeviation)
		  && (P[i].AGS_Hsml < maxsoft || maxsoft == 0)) ||
		 (P[i].AGS_NumNgb > (desnumngb + All.AGS_MaxNumNgbDeviation) && P[i].AGS_Hsml > minsoft))
#else
	      if(P[i].AGS_NumNgb < (desnumngb - All.AGS_MaxNumNgbDeviation) ||
		 (P[i].AGS_NumNgb > (desnumngb + All.AGS_MaxNumNgbDeviation)))
#endif

		{
		  /* need to redo this particle */
		  npleft++;

		  if(Left[i] > 0 && Right[i] > 0)
		    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
		      {

			/* this one should be ok */
			npleft--;
			P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
			continue;
		      }

		  if(P[i].AGS_NumNgb < (desnumngb - All.AGS_MaxNumNgbDeviation))
		    Left[i] = DMAX(P[i].AGS_Hsml, Left[i]);
		  else
		    {
		      if(Right[i] != 0)
			{
			  if(P[i].AGS_Hsml < Right[i])
			    Right[i] = P[i].AGS_Hsml;
			}
		      else
			Right[i] = P[i].AGS_Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      VERBOSE
			(3,"i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, P[i].AGS_Hsml, Left[i], Right[i],
			 (float) P[i].AGS_NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);

		    }

		  if(Right[i] > 0 && Left[i] > 0)
		    P[i].AGS_Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		  else
		    {
		      if(Right[i] == 0 && Left[i] == 0)
			endrun(8188);	/* can't occur */

		      if(Right[i] == 0 && Left[i] > 0)
			{
			  if(fabs(P[i].AGS_NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (P[i].AGS_NumNgb -
					 desnumngb) / (NUMDIMS * P[i].AGS_NumNgb) * P[i].AGS_omega;


			      if(fac < 1.26)
				P[i].AGS_Hsml *= fac;
			      else
				P[i].AGS_Hsml *= 1.26;
			    }
			  else
			    P[i].AGS_Hsml *= 1.26;
			}

		      if(Right[i] > 0 && Left[i] == 0)
			{
			  if(fabs(P[i].AGS_NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (P[i].AGS_NumNgb -
					 desnumngb) / (NUMDIMS * P[i].AGS_NumNgb) * P[i].AGS_omega;


			      if(fac > 1 / 1.26)
				P[i].AGS_Hsml *= fac;
			      else
				P[i].AGS_Hsml /= 1.26;
			    }
			  else
			    P[i].AGS_Hsml /= 1.26;
			}
		    }


		}
	      else
		{
#ifdef PMGRID
		  if(P[i].AGS_NumNgb < (desnumngb - All.AGS_MaxNumNgbDeviation) ||
		     P[i].AGS_NumNgb > (desnumngb + All.AGS_MaxNumNgbDeviation))
		    {
		      P[i].AGS_zeta = 0;
		      P[i].AGS_omega = 1.;
		    }

#endif
		  P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
		}

	    }

	}


      tend = second();
      timecomp1 += timediff(tstart, tend);

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
	{
	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      VERBOSE(2,"ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	    }

	  if(iter > MAXITER)
	    {
	      VERBOSE(0,"failed to converge in neighbour iteration in ags_density()\n");
	      endrun(1155);
	    }
	}
    }
  while(ntot > 0);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Right);
  myfree(Left);
  myfree(Ngblist);


/* mark as active again */
#ifndef AGS_UPDATEALLPARTICLES
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
#else
  for(i = 0; i < NumPart; i++)
#endif
    {
      if(P[i].TimeBin < 0)
	P[i].TimeBin = -P[i].TimeBin - 1;

    }


/* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp;
  CPU_Step[CPU_AGSDENSWAIT] += timewait;
  CPU_Step[CPU_AGSDENSCOMM] += timecomm;
  CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);

}

int ags_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
			 int *ngblist)
{
  int j, n;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h, h2, hinv, hinv3, hinv4;
  double u3, u4, u5;
  MyLongDouble rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyLongDouble weighted_numngb;
#ifdef AGS_OUTPUTNGBS
  int real_numngb = 0;
#endif
  MyDouble *pos;
  double dPhidh, dWdh;
  MyFloat zeta, omega;

  rho = weighted_numngb = 0;

  dPhidh = dWdh = 0;
  zeta = omega = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = P[target].AGS_Hsml;

    }
  else
    {
      pos = AGS_DensDataGet[target].Pos;
      h = AGS_DensDataGet[target].AGS_Hsml;

    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
#ifndef  ONEDIM
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv;
#endif
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;



  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = AGS_DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{

	  numngb_inbox =
	    ags_ngb_treefind_variable_threads(pos, h, target, &startnode, mode, exportflag, exportnodecount,
					      exportindex, ngblist);


	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {

	      j = ngblist[n];

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      dx = NEAREST_X(dx);
	      dy = NEAREST_Y(dy);
	      dz = NEAREST_Z(dz);
#endif
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(r2 < h2)
		{
		  numngb++;

		  r = sqrt(r2);

		  u = r * hinv;
		  u3 = u * u * u;
		  u4 = u3 * u;
		  u5 = u4 * u;


		  if(u < 0.5)
		    {
		      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		      dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);

		      /* compute quantities for zeta and omega: CUBIC kernel is assumed */

//
		      dPhidh = hinv * hinv * (-16. * u * u + 48. * u4 - 38.4 * u5 + 2.8);

//
//                      dPhidh = hinv * hinv *  (-16. * u * u + 48. * pow(u,4.) - 38.4 * pow(u,5.) + 2.8);
//                      dWdh  = hinv4 * KERNEL_COEFF_1 * (-3. + 30. * u * u - 36. * u * u * u);
//
		    }
		  else
		    {
		      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);

		      /* compute quantities for zeta and omega: CUBIC kernel is assumed */

		      dPhidh = hinv * hinv * (-32. * u * u + 64. * u3 - 48. * u4 + 12.8 * u5 + 3.2);

//
//                    dPhidh = hinv * hinv * (-32. * u * u + 64. * u * u * u - 48. * pow(u,4.) + 12.8 * pow(u,5.) + 3.2);
//                    dWdh  =  hinv4 * KERNEL_COEFF_1 * (-6. + 24. * u - 30. * u * u + 12. * u * u * u);
//
		    }

		  mass_j = P[j].Mass;


		  rho += FLT(wk);

		  zeta += FLT(mass_j * dPhidh);

//                                omega += FLT(dWdh);
		  omega += FLT(-1 * (NUMDIMS * hinv * wk + u * dwk));

		  weighted_numngb += FLT(NORM_COEFF * wk / hinv3);	/* 4.0/3 * PI = 4.188790204786 */

#ifdef AGS_OUTPUTNGBS
		  real_numngb += 1;
#endif
		}
	    }
	}


      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = AGS_DensDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    {
#ifdef AGS_OUTPUTNGBS
      P[target].AGS_REALNumNgb = real_numngb;
#endif
      P[target].AGS_NumNgb = weighted_numngb;
      P[target].AGS_Density = rho;
      P[target].AGS_zeta = zeta;
      P[target].AGS_omega = omega;
    }
  else
    {

      AGS_DensDataResult[target].AGS_Rho = rho;
#ifdef AGS_OUTPUTNGBS
      AGS_DensDataResult[target].AGS_REALNgb = real_numngb;
#endif
      AGS_DensDataResult[target].AGS_Ngb = weighted_numngb;
      AGS_DensDataResult[target].AGS_zeta = zeta;
      AGS_DensDataResult[target].AGS_omega = omega;

    }

  return 0;

}


void *ags_density_evaluate_primary(void *p)
{
  int thread_id = *(int *) p;
  int i, j;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;
  double minsoft, maxsoft;

  ngblist = Ngblist + thread_id * NumPart;
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {

      LOCK_NEXPORT;
#ifndef AGS_UPDATEALLPARTICLES
      if(BufferFullFlag != 0 || NextParticle < 0)
#else
      if(BufferFullFlag != 0 || NextParticle == NumPart)
#endif
	{
	  UNLOCK_NEXPORT;
	  break;
	}

      i = NextParticle;
      ProcessedFlag[i] = 0;

#ifndef AGS_UPDATEALLPARTICLES
      NextParticle = NextActiveParticle[NextParticle];
#else
      NextParticle += 1;
#endif

      UNLOCK_NEXPORT;

      if(ags_density_isactive(i))
	{

#ifdef PMGRID
	  /* All.AGS_MaxSoft=0 when this routine is called for the first time, namely
	   * in ags_setup_smoothinglengths(). This is because the value for Asmth
	   * is set for the first time AFTER the init() call in begrun().
	   * We drop the requirement for the softening not to exceed the desired maximum value the first
	   * time ags_density() is called, since nothing interesting will take place before ags_density()
	   * is called again within run().*/

#ifdef PLACEHIGHRESREGION
	  if(pmforce_is_particle_high_res(P[i].Type, P[i].Pos))
	    {
	      minsoft = All.AGS_MinSoft[1];
	      maxsoft = All.AGS_MaxSoft[1];
	    }
	  else
	    {
	      minsoft = All.AGS_MinSoft[0];
	      maxsoft = All.AGS_MaxSoft[0];
	    }
#else
	  minsoft = All.AGS_MinSoft[0];
	  maxsoft = All.AGS_MaxSoft[0];
#endif
	  if((P[i].AGS_Hsml > maxsoft) && maxsoft != 0)
	    P[i].AGS_Hsml = maxsoft;

	  if(P[i].AGS_Hsml < minsoft)
	    P[i].AGS_Hsml = minsoft;
#endif

	  if(ags_density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}
      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}

void *ags_density_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;

  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;


  while(1)
    {
      LOCK_NEXPORT;
      j = NextJ;
      NextJ++;
      UNLOCK_NEXPORT;

      if(j >= Nimport)
	break;

      ags_density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}


int ags_density_isactive(int n)
{
  if(P[n].TimeBin < 0)
    return 0;
  else
    return 1;
}
