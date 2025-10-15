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

#if !defined(GADGET3_IO_LIB) || defined(GADGET3_IO_LIB_MPI)
#	include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#ifdef GADGET3_IO_LIB
#	include "gadget3_io_lib.h"
#endif

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"

#ifdef GM_MUPPI
#include "../CoolingSfr/Muppi/muppi.h"
#endif

/* This function reads initial conditions that are in the default file format
 * of Gadget, i.e. snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with restartflag==2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameterfile.  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasTemp>0 is given, the gas temperature will be initialzed to this
 * value assuming a mean colecular weight either corresponding to complete
 * neutrality, or full ionization.
 */

#ifdef AUTO_SWAP_ENDIAN_READIC
int swap_file = 8;
#endif

#ifdef LT_STELLAREVOLUTION
double time_age, this_age;
int Zs_present = 0, ti_step, IMFi;

#ifdef LT_ZAGE
int ZAge_present = 0;
#endif
#ifdef LT_ZAGE_LLV
int ZAge_llv_present = 0;
#endif
#endif

int count = 0;

void read_ic(char *fname)
{
  int i, j, k, num_files, rest_files, ngroups, gr, filenr, masterTask, lastTask, groupMaster;
  double u_init, molecular_weight;
  char buf[500];

  CPU_Step[CPU_MISC] += measure_time();

#ifdef RESCALEVINI
  if(ThisTask == 0 && RestartFlag == 0)
    {
      fprintf(stdout, "\nRescaling v_ini !\n\n");
      fflush(stdout);
    }
#endif

  NumPart = 0;
  N_gas = 0;
  All.TotNumPart = 0;
  N_stars = 0;
  N_BHs = 0;

  num_files = find_files(fname);

  rest_files = num_files;

  if(ThisTask == 0)
    printf("num_files = %i\n", num_files);

  while(rest_files > NTask)
    {
      sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
      if(All.ICFormat == 3)
	sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));

      ngroups = NTask / All.NumFilesWrittenInParallel;
      if((NTask % All.NumFilesWrittenInParallel))
	ngroups++;
      groupMaster = (ThisTask / ngroups) * ngroups;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if(ThisTask == (groupMaster + gr))	/* ok, it's this processor's turn */

#ifdef GADGET3_IO_LIB
	    read_file(buf, ThisTask, ThisTask, num_files);
#else
		read_file(buf, ThisTask, ThisTask);
#endif
	  MPI_Barrier(MYMPI_COMM_WORLD);
	}

      rest_files -= NTask;
    }


  if(rest_files > 0)
    {
      distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(num_files > 1)
	{
	  sprintf(buf, "%s.%d", fname, filenr);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.%d.hdf5", fname, filenr);
	}
      else
	{
	  sprintf(buf, "%s", fname);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.hdf5", fname);
	}

      ngroups = rest_files / All.NumFilesWrittenInParallel;
      if((rest_files % All.NumFilesWrittenInParallel))
	ngroups++;

      for(gr = 0; gr < ngroups; gr++)
	{

	  if((filenr / num_files) == gr)	/* ok, it's this processor's turn */
#ifdef GADGET3_IO_LIB
	    read_file(buf, masterTask, lastTask, num_files);
#else
	    read_file(buf, masterTask, lastTask);
#endif		

	  MPI_Barrier(MYMPI_COMM_WORLD);
	}
    }

  myfree(CommBuffer);

  // Re arrange particles (e.g. bring gas to the front)
  // First step: Find the last gas partilce
  int index_gas_to_be = 0;
  while(index_gas_to_be < N_gas)
    {
      if(P[index_gas_to_be].Type != 0)
	break;
      index_gas_to_be++;
    }
  // Second step: Find the next gas particle which has to be moved
  int index_gas_to_move = index_gas_to_be;
  while(index_gas_to_move < NumPart)
    {
      if(P[index_gas_to_move].Type == 0)
	break;
      index_gas_to_move++;
    }

  struct particle_data psave;

  // Third step: Move all gas particles tio fromt
  while(index_gas_to_be < N_gas && index_gas_to_move < NumPart)
    {
      // repair pointers for extra data for stars and BHs
#ifdef LT_STELLAREVOLUTION
      if(P[index_gas_to_be].Type == 4)
	MetP[P[index_gas_to_be].pt.MetID].PID = index_gas_to_move;
#endif
#if defined(BLACK_HOLES)
      if(P[index_gas_to_be].Type == 5)
	BHP[P[index_gas_to_be].pt.BHID].PID = index_gas_to_move;
#endif
      // swap particles
      psave = P[index_gas_to_be];
      P[index_gas_to_be] = P[index_gas_to_move];
      P[index_gas_to_move] = psave;

      // find next free spot for gas particles
      while(index_gas_to_be < N_gas)
	{
	  if(P[index_gas_to_be].Type != 0)
	    break;
	  index_gas_to_be++;
	}
      // find next gas particle to be move
      while(index_gas_to_move < NumPart)
	{
	  if(P[index_gas_to_move].Type == 0)
	    break;
	  index_gas_to_move++;
	}
    }

  // here we do a sanity check if there unexpected gas particles left
  while(index_gas_to_move < NumPart)
    {
      if(P[index_gas_to_move].Type == 0)
	{
	  printf("Task %d: SOMETHING WRONG ! Seems to be gas at %d, N_gas=%d, N_tot=%d, type=%d count=%d\n",
		 ThisTask, index_gas_to_move, N_gas, NumPart, P[index_gas_to_move].Type, count);
	  for(i = 0; i < N_gas; i++)
	    if(P[i].Type != 0)
	      {
		printf("Task %d: EVEN MORE WRONG ! Seems to be a non gas at %d, N_gas=%d type=%d \n",
		       ThisTask, i, N_gas, P[i].Type);
		endrun(3894573895);
	      }
	  endrun(8357837);
	}
      index_gas_to_move++;
    }

  if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
    {
      /* this makes sure that masses are initialized in the case that the mass-block
         is empty for this particle type */
      for(i = 0; i < NumPart; i++)
	{
	  if(All.MassTable[P[i].Type] != 0)
	    P[i].Mass = All.MassTable[P[i].Type];
	}
    }


#ifdef GENERATE_GAS_IN_ICS
  int count;

  if(RestartFlag == 0)
    {
      header.flag_entropy_instead_u = 0;

      check_particles();

      for(i = 0, count = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 1)
	    count++;

	  if(P[i].Type == 0)
	    {
	      printf
		("Task %d: ICs already contain gas particles but flag GENERATE_GAS_IN_ICS is set. We better stop.\n",
		 ThisTask);
	      endrun(81234091);
	    }
	}

      if(count != NumPart)
	{
	  printf("Task=%d not only type 1 partcles %d/%d !!\n", ThisTask, count, NumPart);
	}

      if(NumPart + count > All.TotNumPart)
	{
	  printf("Task=%d ends up getting more particles (%d) than allowed (%ld)\n",
		 ThisTask, NumPart + count, All.TotNumPart);
	  endrun(222);
	}

      if(N_gas + count > All.MaxPartSph)
	{
	  printf("Task=%d ends up getting more SPH particles (%d) than allowed (%d)\n",
		 ThisTask, N_gas, All.MaxPartSph);
	  endrun(111);
	}

      double fac = All.OmegaBaryon / All.Omega0;
      double rho = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      for(i = 0; i < count; i++)
	{
	  P[i + count] = P[i];
	  P[i + count].ID += 300000000000;
	  P[i].Type = 0;

	  double d = pow(P[i].Mass / rho, 1.0 / 3);
	  double a = 0.5 * All.OmegaBaryon / All.Omega0 * d;
	  double b = 0.5 * (All.Omega0 - All.OmegaBaryon) / All.Omega0 * d;

	  P[i].Mass = All.MassTable[1] * fac;
	  P[i + count].Mass = All.MassTable[1] * (1 - fac);
	  P[i + count].Pos[0] += a;
	  P[i + count].Pos[1] += a;
	  P[i + count].Pos[2] += a;
	  P[i].Pos[0] -= b;
	  P[i].Pos[1] -= b;
	  P[i].Pos[2] -= b;

	}

      NumPart += count;
      N_gas += count;

      All.MassTable[0] = 0;
      All.MassTable[1] *= (1 - fac);

#ifdef PEDANTIC_CHECK
      if(ThisTask == 0)
	{
	  printf(" Checking Particles after adding gas %d/%d/%d\n", NumPart, N_gas, NumPart - N_gas);
	  int count_gas = 0, count_dm = 0;
	  for(i = 0; i < NumPart; i++)
	    {
	      if(P[i].Type == 0 && i >= N_gas)
		printf("  Gas paricle at %d >= %d !!\n", i, N_gas);
	      if(P[i].Type == 0 && P[i].Mass != P[0].Mass)
		printf("  Gas mass wrong (%g/%g) !!\n", P[i].Mass, P[0].Mass);
	      if(P[i].Type == 1 && P[i].Mass != P[N_gas].Mass)
		printf("  DM mass wrong (%g/%g) !!\n", P[i].Mass, P[N_gas].Mass);
	      if(P[i].Type == 0)
		count_gas++;
	      if(P[i].Type == 1)
		count_dm++;
	    }
	  printf(" Checking Particles after adding gas %d/%d/%d done ...\n", NumPart, count_gas, count_dm);
	}
#endif
    }
#endif



#if defined(BLACK_HOLES)
  if(RestartFlag == 0)
    {
      All.MassTable[5] = 0;
    }
#endif

#ifdef SFR
  if(RestartFlag == 0)
    {
      if(All.MassTable[4] == 0 && All.MassTable[0] > 0)
	{
	  All.MassTable[0] = 0;
	  All.MassTable[4] = 0;
	}
    }
#endif

#ifdef LT_STELLAREVOLUTION

#ifndef GADGET3_IO_LIB
  All.Time_Age = get_age(All.Time);
#endif

  /* note: initialization of timings for stellar evolution occurs in init.c */

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  if(!Zs_present)
	    {
	      //              memset(SphP[i].Metals, 0, sizeof(float) * LT_NMetP);
	      SphP[i].Metals[Hel] = P[i].Mass * (1 - HYDROGEN_MASSFRAC);
	    }

#ifdef GL_CR_DUST
	  if(RestartFlag == 0)
	    {
	      memset(SphP[i].DustL, 0, sizeof(float) * LT_NMetP);
	      memset(SphP[i].DustS, 0, sizeof(float) * LT_NMetP);
	    }

	  SphP[i].numSnIa=0.; /* XXXX may be a problem when restarting from snapshots */
	  SphP[i].numSnII=0.;
#endif


	}
      else if(P[i].Type == 4)
	{
	  //MetP[N_star_idx].PID = i;
	  //P[i].pt.MetID = N_star_idx;

	  if(!Zs_present)
	    {
	      memset(MPP(i).Metals, 0, sizeof(float) * LT_NMetP);
	      MPP(i).Metals[Hel] = P[i].Mass * (1 - HYDROGEN_MASSFRAC);
	    }
	}
    }

  if(!Zs_present)
    /* metal array is not present; we initialize Helium to the cosmic fraction */
    {
      if(ThisTask == 0)
	printf("Helium Fraction initialized to the Cosmic Value %8.5g\n", 1 - HYDROGEN_MASSFRAC);
    }

#ifdef PEDANTIC_CHECK
  if(ThisTask == 0)
    {
      printf(" Checking Particles after setting metals %d/%d/%d\n", NumPart, N_gas, NumPart - N_gas);
      int count_gas = 0, count_dm = 0;
      for(i = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 0 && i >= N_gas)
	    printf("  Gas paricle at %d >= %d !!\n", i, N_gas);
	  if(P[i].Type == 0 && P[i].Mass != P[0].Mass)
	    printf("  Gas mass wrong (%g/%g) !!\n", P[i].Mass, P[0].Mass);
	  if(P[i].Type == 1 && P[i].Mass != P[N_gas].Mass)
	    printf("  Gas mass wrong (%g/%g) !!\n", P[i].Mass, P[N_gas].Mass);
	  if(P[i].Type == 0)
	    count_gas++;
	  if(P[i].Type == 1)
	    count_dm++;
	}
      printf(" Checking Particles after setting metals %d/%d/%d done ...\n", NumPart, count_gas, count_dm);
    }
#endif
#endif

  u_init = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * DMAX(All.InitGasTemp, 1e-20);
  u_init *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;	/* unit conversion */

  if(All.InitGasTemp > 1.0e4)	/* assuming FULL ionization */
    molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  else				/* assuming NEUTRAL GAS */
    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);

  u_init /= molecular_weight;

  All.InitGasU = u_init;

  if(RestartFlag == 0)
    for(i = 0; i < N_gas; i++)
      if(SphP[i].Entropy == 0)
	SphP[i].Entropy = All.InitGasU;


  //  for(i = 0; i < N_gas; i++)
  //    SphP[i].Entropy = DMAX(All.MinEgySpec, SphP[i].Entropy);

#ifdef PEDANTIC_CHECK
  if(ThisTask == 0)
    {
      printf(" Checking Particles after setting gas entropy %d/%d/%d\n", NumPart, N_gas, NumPart - N_gas);
      int count_gas = 0, count_dm = 0;
      for(i = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 0 && i >= N_gas)
	    printf("  Gas paricle at %d >= %d !!\n", i, N_gas);
	  if(P[i].Type == 0 && P[i].Mass != P[0].Mass)
	    printf("  Gas mass wrong (%g/%g) !!\n", P[i].Mass, P[0].Mass);
	  if(P[i].Type == 1 && P[i].Mass != P[N_gas].Mass)
	    printf("  Gas mass wrong (%g/%g) !!\n", P[i].Mass, P[N_gas].Mass);
	  if(P[i].Type == 0)
	    count_gas++;
	  if(P[i].Type == 1)
	    count_dm++;
	}
      printf(" Checking Particles after setting gas entropy %d/%d/%d done ...\n", NumPart, count_gas,
	     count_dm);
    }
#endif


  MPI_Barrier(MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("reading done.\n");
      fflush(stdout);
    }

  if(ThisTask == 0)
    {
      printf("Total number of particles :  %llu/%d\n", All.TotNumPart, NumPart);
      printf("Task %d: Total number of particles : %d\n", ThisTask, NumPart);
      printf("Task %d: Total number of Gas part. : %d\n", ThisTask, N_gas);
      printf("Task %d: Total number of Star part.: %d\n", ThisTask, N_stars);
      printf("Task %d: Total number of BH part.  : %d\n\n", ThisTask, N_BHs);

      fflush(stdout);
    }

  CPU_Step[CPU_SNAPSHOT] += measure_time();

#ifdef PEDANTIC_CHECK
  if(ThisTask == 0)
    {
      printf(" Checking Particles after reading %d/%d/%d\n", NumPart, N_gas, NumPart - N_gas);
      int count_gas = 0, count_dm = 0;
      for(i = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 0 && i >= N_gas)
	    printf("  Gas paricle at %d >= %d !!\n", i, N_gas);
	  if(P[i].Type == 0 && P[i].Mass != P[0].Mass)
	    printf("  Gas mass wrong (%g/%g) !!\n", P[i].Mass, P[0].Mass);
	  if(P[i].Type == 1 && P[i].Mass != P[N_gas].Mass)
	    printf("  Gas mass wrong (%g/%g) !!\n", P[i].Mass, P[N_gas].Mass);
	  if(P[i].Type == 0)
	    count_gas++;
	  if(P[i].Type == 1)
	    count_dm++;
	}
      printf(" Checking Particles after reading %d/%d/%d done ...\n", NumPart, count_gas, count_dm);
    }
#endif

}


/*! This function reads out the buffer that was filled with particle data.
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type, int offset_gas, int offset_stars,
		       int offset_bhs)
{
  int n, k, i, j;
  MyInputFloat *fp;
  MyIDType *ip;
  float *fp_single;
#ifdef GM_MUPPI
  int *ip_single;
#endif

#ifdef AUTO_SWAP_ENDIAN_READIC
  int vt, vpb;
  char *cp;
#endif

  fp = (MyInputFloat *) CommBuffer;
  fp_single = (float *) CommBuffer;
  ip = (MyIDType *) CommBuffer;
#ifdef GM_MUPPI
  ip_single = (int *) CommBuffer;
#endif

#ifdef AUTO_SWAP_ENDIAN_READIC
  if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP)
    {
      cp = (char *) CommBuffer;
      vt = get_datatype_in_block(blocknr);
      vpb = get_values_per_blockelement(blocknr);
      if(vt == 2)
	swap_Nbyte(cp, pc * vpb, 8);
      else
	{
#ifdef INPUT_IN_DOUBLEPRECISION
	  if(vt == 1)
	    swap_Nbyte(cp, pc * vpb, 8);
	  else
#endif
	    swap_Nbyte(cp, pc * vpb, 4);
	}
    }
#endif


  switch (blocknr)
    {
// Note that we use the reading of the postion array to also initialize the type
// as well as setting the pointers for the extra data structures for stars and BHs
    case IO_POS:		/* positions */
      for(n = 0; n < pc; n++)
	{
	  for(k = 0; k < 3; k++)
	    P[offset + n].Pos[k] = *fp++;
	}

#if defined(DS_SHIFT_BOX) && defined(PERIODIC)
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  P[offset + n].Pos[k] += All.BoxSize / 2;
#endif

      for(n = 0; n < pc; n++)
	{
	  P[offset + n].Type = type;	/* initialize type here as well */
	  if(type == 0)
	    count++;

#ifdef LT_STELLAREVOLUTION
	  if(P[offset + n].Type == 4)	/* initialize pointer to extra informations for stars */
	    {
	      MetP[offset_stars + n].PID = offset + n;
	      P[offset + n].pt.MetID = offset_stars + n;
	    }
#endif
#if defined(BLACK_HOLES)
	  if(P[offset + n].Type == 5)	/* initialize pointer to extra informations for BHs */
	    {
	      BHP[offset_bhs + n].PID = offset + n;
	      P[offset + n].pt.BHID = offset_bhs + n;
	    }
#endif
	}
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
#ifdef RESCALEVINI
	  /* scaling v to use same IC's for different cosmologies */
	  if(RestartFlag == 0)
	    P[offset + n].Vel[k] = (*fp++) * All.VelIniScale;
	  else
	    P[offset + n].Vel[k] = *fp++;
#else
	  P[offset + n].Vel[k] = *fp++;
#endif
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; n++)
	{
	  P[offset + n].ID = *ip++;
	}
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; n++)
	P[offset + n].Mass = *fp++;
      break;


    case IO_SHEET_ORIENTATION:	/* initial particle sheet orientation */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
      for(n = 0; n < pc; n++)
	{
#ifndef GDE_LEAN
	  P[offset + n].V_matrix[0][0] = *fp++;
	  P[offset + n].V_matrix[0][1] = *fp++;
	  P[offset + n].V_matrix[0][2] = *fp++;
	  P[offset + n].V_matrix[1][0] = *fp++;
	  P[offset + n].V_matrix[1][1] = *fp++;
	  P[offset + n].V_matrix[1][2] = *fp++;
	  P[offset + n].V_matrix[2][0] = *fp++;
	  P[offset + n].V_matrix[2][1] = *fp++;
	  P[offset + n].V_matrix[2][2] = *fp++;
#else
	  *fp += 8;
#endif
	}
#endif
      break;

    case IO_INIT_DENSITY:	/* initial stream density */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
      for(n = 0; n < pc; n++)
	GDE_INITDENSITY(offset + n) = *fp++;
      break;
#endif

    case IO_CAUSTIC_COUNTER:	/* initial caustic counter */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
      for(n = 0; n < pc; n++)
	P[offset + n].caustic_counter = *fp++;
      break;
#endif

    case IO_SECONDORDERMASS:
      for(n = 0; n < pc; n++)
	{
	  P[offset + n].OldAcc = P[offset + n].Mass;	/* use this to temporarily store the masses in the 2plt IC case */
	  P[offset + n].Mass = *fp++;
	}
      break;

    case IO_U:			/* temperature */
      for(n = 0; n < pc; n++)
	{
#ifdef JD_RELAX_CLUSTERS
	  SphP[offset_gas + n].U = *fp;
#endif
	  SphP[offset_gas + n].Entropy = *fp++;
	}
      break;

    case IO_RHO:		/* density */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Density = *fp++;
      break;

    case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(CHEMISTRY) || defined(UM_CHEMISTRY)
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].elec = *fp++;
#endif
      break;

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)	//now that's weird. But I am afraid of removing this #if, I am not sure where it ends. AR.
    case IO_NH:		/* neutral hydrogen abundance */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HI = *fp++;
      break;

    case IO_HII:		/* ionized hydrogen abundance */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HII = *fp++;
      break;

    case IO_HeI:		/* neutral Helium */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HeI = *fp++;
      break;

    case IO_HeII:		/* ionized Heluum */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HeII = *fp++;
      break;

    case IO_HeIII:		/* double ionised Helium */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HeIII = *fp++;
      break;
#endif

    case IO_H2I:		/* H2 molecule */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].H2I = *fp++;
      break;

    case IO_H2II:		/* ionised H2 molecule */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].H2II = *fp++;
      break;

    case IO_HM:		/* H minus */
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HM = *fp++;
      break;

    case IO_HeHII:		/* HeH+ */
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HeHII = *fp++;
#endif
      break;

    case IO_HD:		/* HD */
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].HD = *fp++;
#endif
      break;

    case IO_DI:		/* D */
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].DI = *fp++;
#endif
      break;

    case IO_DII:		/* D plus */
#if defined (UM_CHEMISTRY)
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].DII = *fp++;
#endif
      break;

#else
    case IO_NH:		/* neutral hydrogen abundance */
    case IO_HII:		/* ionized hydrogen abundance */
    case IO_HeI:		/* neutral Helium */
    case IO_HeII:		/* ionized Heluum */
    case IO_HeIII:		/* double ionised Helium */
    case IO_H2I:		/* H2 molecule */
    case IO_H2II:		/* ionised H2 molecule */
    case IO_HM:		/* H minus */
    case IO_HeHII:		/* HeH+ */
    case IO_HD:		/* HD */
    case IO_DI:		/* D */
    case IO_DII:		/* D plus  */
      break;
#endif

    case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; n++)
	P[offset + n].Hsml = *fp++;
      break;

    case IO_DELAYTIME:
#ifdef WINDS
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].DelayTime = *fp++;
#endif
      break;

    case IO_AGE:		/* Age of stars */
#ifdef STELLARAGE
      for(n = 0; n < pc; n++)
	{
	  if(P[offset + n].Type == 4)
	    MPP(offset + n).StellarAge = *fp++;
#if defined(BLACK_HOLES)
	  if(P[offset + n].Type == 5)
	    BPP(offset + n).StellarAge = *fp++;
#endif
	}
#endif
      break;

    case IO_Z:			/* Gas and star metallicity */
#ifdef METALS
      for(n = 0; n < pc; n++)
	P[offset + n].Metallicity = *fp++;
#endif
      break;

    case IO_EGYPROM:		/* SN Energy Reservoir */
      break;

    case IO_EGYCOLD:		/* Cold  SN Energy Reservoir */
      break;

    case IO_VRMS:		/* Turbulence on kernel scale */
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Vrms = *fp++;
#endif
      break;
    case IO_VBULK:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  SphP[offset_gas + n].Vbulk[k] = *fp++;
#endif
      break;
    case IO_VTAN:
#ifdef JD_DECOMPOSE_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Vtan = *fp++;
#endif
      break;
    case IO_VRAD:
#ifdef JD_DECOMPOSE_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Vrad = *fp++;
#endif
      break;
    case IO_VDIV:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].DivVel = *fp++;
#endif
      break;
    case IO_VROT:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].r.CurlVel = *fp++;
#endif
      break;
    case IO_TRUENGB:
#ifdef JD_VTURB
      for(n = 0; n < pc; n++)
	P[offset + n].TrueNGB = *fp++;
#endif
      break;
    case IO_DPP:
#ifdef JD_DPP
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Dpp = *fp++;
#endif
      break;

    case IO_BFLD:		/* Magnetic field */
#ifdef MAGNETIC
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  SphP[offset_gas + n].BPred[k] = *fp++;
#ifdef TRACEDIVB
      SphP[offset_gas + n].divB = 0;
#endif
#ifdef MAGNETICZERO
      for(k = 0; k < 3; k++)
	SphP[offset_gas + n].BPred[k] = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
      SphP[offset_gas + n].Phi = 0;
      SphP[offset_gas + n].PhiPred = 0;
#endif
#endif
      break;

    case IO_LMBCR_pNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
      for(n = 0; n < pc; n++)
	for(k = 0; k < LMB_CR_PROTONS; k++)
	  {
	    SphP[offset_gas + n].CRpNorm[k] = *fp++;
	    SphP[offset_gas + n].CRpNorm[k] = pow(10.0, SphP[offset_gas + n].CRpNorm[k]);
	  }
#endif
      break;

    case IO_LMBCR_eNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
      for(n = 0; n < pc; n++)
	for(k = 0; k < LMB_CR_ELECTRONS; k++)
	  {
	    SphP[offset_gas + n].CReNorm[k] = *fp++;
	    SphP[offset_gas + n].CReNorm[k] = pow(10.0, SphP[offset_gas + n].CReNorm[k]);
	  }
#endif
      break;

    case IO_LMBCR_pSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
      for(n = 0; n < pc; n++)
	for(k = 0; k < LMB_CR_PROTONS; k++)
	  SphP[offset_gas + n].CRpSlope[k] = *fp++;
#endif
      break;

    case IO_LMBCR_eSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
      for(n = 0; n < pc; n++)
	for(k = 0; k < LMB_CR_ELECTRONS; k++)
	  SphP[offset_gas + n].CReSlope[k] = *fp++;
#endif
      break;

    case IO_LMBCR_pCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].CRpCut = *fp++;
#endif
      break;

    case IO_LMBCR_eCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].CReCut = *fp++;
#endif
      break;

    case IO_LMBCR_pPRESSURE:
#if defined LMB_SPECTRAL_CRs_PRESSURE_FROM_ICS
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].CRpPressure = *fp++;
#endif
      break;

    case IO_LMBCR_ePRESSURE:
#if defined LMB_SPECTRALCRs_PRESSURE_FROM_ICS
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].CRePressure = *fp++;
#endif
      break;

    case IO_RHO_OLD:
#if defined READ_RHO_OLD
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].DensityOld = *fp++;
#endif
      break;

    case IO_MACH:
#if defined READ_MACH
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Shock_Mach = *fp++;
#endif
      break;

    case IO_BHMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	BPP(offset + n).BH_Mass = *fp++;
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	BPP(offset + n).BH_Mdot = *fp++;
#endif
      break;

    case IO_BHPROGS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	BPP(offset + n).BH_CountProgs = *fp++;
#endif
      break;

    case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
      for(n = 0; n < pc; n++)
	BPP(offset + n).BH_Mass_radio = *fp++;
#endif
      break;

    case IO_EOSXNUC:
      break;
    case IO_DUSTL:

#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
      if(type == 0)
	for(n = 0; n < pc; n++)
	  for(k=0; k<LT_NMetP; k++)
	    SphP[offset_gas+n].DustL[k] = (float) *fp++;
#endif
      break;
    case IO_DUSTS:
#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
      if(type == 0)
	for(n = 0; n < pc; n++)
 	  for(k=0; k<LT_NMetP; k++)
	    SphP[offset_gas+n].DustS[k] = (float) *fp++;
#endif
      break;

    case IO_Zs:
#ifdef LT_STELLAREVOLUTION
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    for(k = 0; k < LT_NMetP; k++)
	      MPP(offset + n).Metals[k] = (float) *fp++;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  for(k = 0; k < LT_NMetP; k++)
	    SphP[offset_gas + n].Metals[k] = (float) *fp++;
#endif
      break;

    case IO_ZAGE:
#ifdef LT_ZAGE
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    MPP(offset + n).ZAge = *fp++;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    /* note this is not the weight that was used when the snapshot has been written */
	    SphP[offset_gas + n].ZAgeW = get_metalmass(SphP[offset_gas + n].Metals);
#ifndef LT_LOGZAGE
	    SphP[offset_gas + n].ZAge = *fp++ * SphP[offset_gas + n].ZAgeW;
#else
	    if(SphP[offset_gas + n].ZAgeW > 0)
	      SphP[offset_gas + n].ZAge = log10(*fp++ * SphP[offset_gas + n].ZAgeW);
	    else
	      SphP[offset_gas + n].ZAge = 0;
#endif
	  }
#endif
      break;

    case IO_ZAGE_LLV:
#ifdef LT_ZAGE_LLV
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    MPP(offset + n).ZAge_llv = *fp++;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    /* note this is not the weight that was used when the snapshot has been written */
	    SphP[offset_gas + n].ZAgeW_llv = SphP[offset_gas + n].Metals[Iron];
#ifndef LT_LOGZAGE
	    SphP[offset_gas + n].ZAge_llv = *fp++ * SphP[offset_gas + n].ZAgeW_llv;
#else
	    if(SphP[offset_gas + n].ZAgeW_llv > 0)
	      SphP[offset_gas + n].ZAge_llv = log10(*fp++ * SphP[offset_gas + n].ZAgeW_llv);
	    else
	      SphP[offset_gas + n].ZAge_llv = 0;
#endif
	  }
#endif
      break;

    case IO_iMass:
#ifdef LT_STELLAREVOLUTION
      if(type == 4)
	for(n = 0; n < pc; n++)
	  MPP(offset + n).iMass = *fp++;
#endif
      break;

    case IO_CONTRIB:
#if defined(LT_STELLAREVOLUTION) && defined(LT_TRACK_CONTRIBUTES)
      if(type == 4)
	{
	  for(n = 0; n < pc; n++)
	    MPP(offset + n).contrib = *contribp++;
	}
      else if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset_gas + n].contrib = *contribp++;
#endif
      break;

    case IO_RADGAMMA:
      break;

    case IO_DMHSML:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].DM_Hsml = *fp_single++;
#endif
      break;

    case IO_DMDENSITY:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].u.DM_Density = *fp_single++;
#endif
      break;

    case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].v.DM_VelDisp = *fp_single++;
#endif
      break;

    case IO_EULERA:
#ifdef READ_EULER
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].EulerA = *fp++;
#endif
      break;

    case IO_EULERB:
#ifdef READ_EULER
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].EulerB = *fp++;
#endif
      break;

    case IO_ALFA2_DYN:
      break;

    case IO_ETA2_DYN:
      break;


    case IO_FB_M_H:		/* particle hot mass */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset_gas + n].M_h = *fp_single++;
#endif
      break;

    case IO_FB_M_C:		/* particle cold mass */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset_gas + n].M_c = *fp_single++;
#endif
      break;

    case IO_FB_M_MO:		/* particle molecular mass */
#if defined (GM_MUPPI) && defined(MV_OUTPUT_MMOL)
      if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset + n].Fcoll = *fp_single++;
#endif
      break;

    case IO_FB_E_H:		/* particle hot energy */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    SphP[offset_gas + n].E_h = *fp_single++;
	    SphP[offset_gas + n].E_h *= 1.e50 / All.UnitEnergy_in_cgs;
	  }
#endif
      break;


    case IO_FB_M_SF:		/* particle SF mass */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset_gas + n].M_sf = *fp_single++;
#endif
      break;

    case IO_FB_MF:		/* number of time steps in multiphase */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset_gas + n].MultiPhase = *ip_single++;
#endif
      break;

    case IO_FB_NMF:		/* number of multiphase stages */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset_gas + n].NMF = *ip_single++;
#endif
      break;

    case IO_FB_EOUT:		/* Energy output from MP particles */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    SphP[offset_gas + n].E_out = *fp_single++;
	    SphP[offset_gas + n].E_out *= 1.e50 / All.UnitEnergy_in_cgs;
	  }
#endif
      break;
    case IO_FB_EREC:		/* Energy received from  MP particles */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    SphP[offset_gas + n].E_rec = *fp_single++;
	    SphP[offset_gas + n].E_rec *= 1.e50 / All.UnitEnergy_in_cgs;
	  }
#endif
      break;

    case IO_FB_CLOCK:		/* time to still pass in MP */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  SphP[offset_gas + n].clock = *fp_single++;
#endif
      break;

    case IO_FB_E_TOT_0:	/* multi phase hot energy time i-1 */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    SphP[offset_gas + n].Egy_tot_0 = *fp_single++;
	    SphP[offset_gas + n].Egy_tot_0 *= 1.e50 / All.UnitEnergy_in_cgs;
	  }
#endif
    case IO_FB_TDYN:		/* cold phase dynamical time */
#ifdef GM_MUPPI
      if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    SphP[offset_gas + n].tdyn = *fp_single++;
	  }
      break;

    case IO_FB_TSTARTMP:
      if(type == 0)
	for(n = 0; n < pc; n++)
	  {
	    SphP[offset_gas + n].t_startMP = *fp_single++;
	  }
      break;
#endif


    case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].XColdCloud = *fp++;
#endif
      break;

    case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Temperature = *fp++;
#endif
      break;

    case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
      for(n = 0; n < pc; n++)
	SphP[offset_gas + n].Temperature = *fp++;
#endif
      break;

      /* the other input fields (if present) are not needed to define the
         initial conditions of the code */

    case IO_SFR:
    case IO_ZSMOOTH:
    case IO_allZSMOOTH:
    case IO_POT:
    case IO_ACCEL:
    case IO_DTENTR:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_DBDT:
    case IO_DIVB:
    case IO_ABVC:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_BSMTH:
    case IO_DENN:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_AMDC:
    case IO_PHI:
    case IO_XPHI:
    case IO_GRADPHI:
    case IO_TIDALTENSORPS:
    case IO_ROTB:
    case IO_SROTB:
    case IO_FLOW_DETERMINANT:
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_ANNIHILATION_RADIATION:
    case IO_EOSTEMP:
    case IO_PRESSURE:
    case IO_PRESHOCK_CSND:
    case IO_EDDINGTON_TENSOR:
    case IO_LAST_CAUSTIC:
    case IO_HSMS:
    case IO_ACRS:
    case IO_ACRB:
    case IO_PSUM:
    case IO_SIDMNUMNGB:
    case IO_NUMTOTALSCATTER:
    case IO_SIDMHSML:
    case IO_SIDMDENSITY:
    case IO_SIDMVELDISP:
    case IO_AGS_SOFT:
    case IO_AGS_DENS:
    case IO_AGS_ZETA:
    case IO_AGS_OMEGA:
    case IO_AGS_CORR:
    case IO_AGS_NGBS:
    case IO_MG_PHI:
    case IO_MG_GRAD_PHI:
    case IO_MG_ACCEL:
      break;

    case IO_LASTENTRY:
      endrun(220);
      break;
    }
}



/*! This function reads a snapshot file and distributes the data it contains
 *  to tasks 'readTask' to 'lastTask'.
 */
#ifdef GADGET3_IO_LIB
void read_file(char *fname, int readTask, int lastTask, int num_files)
#else
void read_file(char *fname, int readTask, int lastTask)
#endif
{
  size_t blockmaxlen;
  int i, n_in_file, n_for_this_task, ntask, pc, offset = 0, offset_gas = 0, offset_stars = 0, offset_bhs =
    0, task;
  int blksize1, blksize2;
  MPI_Status status;
  FILE *fd = 0;
  int nall, nread;
  int type, bnr;
  char label[4], buf[500];
  int bytes_per_blockelement, npart, nextblock, typelist[6];
  enum iofields blocknr;
  size_t bytes;

#ifdef HAVE_HDF5
  int rank, pcsum;
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_dataspace_in_file;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory, hdf5_dataset;
  hsize_t dims[2], count[2], start[2];
#endif

#ifdef GM_MUPPI
  FILE *FdMuppi = 0x0, *fd_muppi_back;
  char fmuppiname[500];
#endif


#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	{
	  if(!(fd = fopen(fname, "rb")))
	    {
	      printf("can't open file `%s' for reading initial conditions.\n", fname);
	      endrun(123);
	    }

#ifdef GM_MUPPI
	  sprintf(fmuppiname, "%s.muppi", fname);
	  if((FdMuppi = fopen(fmuppiname, "rb")) != 0x0)
	    {
	      fd_muppi_back = fd;
	    }
	  else
	    printf("can't open auxiliary file %s\n", fmuppiname);
	  fflush(stdout);
#endif

	  if(All.ICFormat == 2)
	    {
	      SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_file = blksize1;
#endif
	      my_fread(&label, sizeof(char), 4, fd);
	      my_fread(&nextblock, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_Nbyte((char *) &nextblock, 1, 4);
#endif
	      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3],
		     nextblock);
	      SKIP2;
	    }

	  SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  if(All.ICFormat == 1)
	    {
	      if(blksize1 != 256)
		swap_file = 1;
	    }
#endif
	  my_fread(&header, sizeof(header), 1, fd);
	  SKIP2;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_Nbyte((char *) &blksize1, 1, 4);
	  swap_Nbyte((char *) &blksize2, 1, 4);
#endif

	  if(blksize1 != 256 || blksize2 != 256)
	    {
	      printf("incorrect header format\n");
	      fflush(stdout);
	      endrun(890);
	      /* Probable error is wrong size of fill[] in header file. Needs to be 256 bytes in total. */
	    }
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_header();
#endif
#ifdef PATCH_IO
	  printf("Task %d: WARNING, PATCHING HEADER FOR MagneticumDM simulation !!!\n", ThisTask);
	  header.mass[2] = header.mass[1];
	  header.mass[1] = header.mass[0];
	  header.mass[0] = 0;
	  header.npart[2] = header.npart[1];
	  header.npart[1] = header.npart[0];
	  header.npart[0] = 0;
	  header.npartTotal[2] = header.npartTotal[1];
	  header.npartTotal[1] = header.npartTotal[0];
	  header.npartTotal[0] = 0;
	  header.npartTotalHighWord[2] = header.npartTotalHighWord[1];
	  header.npartTotalHighWord[1] = header.npartTotalHighWord[0];
	  header.npartTotalHighWord[0] = 0;
#endif
	}


#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
	{
	  read_header_attributes_in_hdf5(fname);

	  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
#ifdef OLD_HDF5
		  hdf5_grp[type] = H5Gopen(hdf5_file, buf);
#else
		  hdf5_grp[type] = H5Gopen(hdf5_file, buf, H5P_DEFAULT);
#endif
		}
	    }
	}
#endif

      for(task = readTask + 1; task <= lastTask; task++)
	{
	  MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MYMPI_COMM_WORLD);
#ifdef AUTO_SWAP_ENDIAN_READIC
	  MPI_Ssend(&swap_file, 1, MPI_INT, task, TAG_SWAP, MYMPI_COMM_WORLD);
#endif
#ifdef GM_MUPPI
	  MPI_Ssend(&FdMuppi, sizeof(FILE *), MPI_BYTE, task, TAG_MPPI, MYMPI_COMM_WORLD);
#endif

	}

    }
  else
    {
      MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MYMPI_COMM_WORLD, &status);
#ifdef AUTO_SWAP_ENDIAN_READIC
      MPI_Recv(&swap_file, 1, MPI_INT, readTask, TAG_SWAP, MYMPI_COMM_WORLD, &status);
#endif
#ifdef GM_MUPPI
      MPI_Recv(&FdMuppi, sizeof(FILE *), MPI_BYTE, readTask, TAG_MPPI, MYMPI_COMM_WORLD, &status);
#endif

    }

#ifdef INPUT_IN_DOUBLEPRECISION
  if(header.flag_doubleprecision == 0)
    {
      if(ThisTask == 0)
	printf
	  ("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
      endrun(11);
    }
#else
  if(header.flag_doubleprecision)
    {
      if(ThisTask == 0)
	printf
	  ("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
      endrun(10);
    }
#endif

  if(All.TotNumPart == 0)
    {
#ifdef GADGET3_IO_LIB		
      if(header.num_files <= 1 || num_files <= 1)
#else
      if(header.num_files <= 1)
#endif	  
	for(i = 0; i < 6; i++)
	  {
	    header.npartTotal[i] = header.npart[i];
#ifdef SFR
	    header.npartTotalHighWord[i] = 0;
#endif
	  }

      All.TotN_gas = header.npartTotal[0] + (((long long) header.npartTotalHighWord[0]) << 32);
#ifdef LT_STELLAREVOLUTION
      All.TotN_stars = header.npartTotal[4];
#endif

      for(i = 0, All.TotNumPart = 0; i < 6; i++)
	{
	  All.TotNumPart += header.npartTotal[i];
	  All.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
	}


#ifdef GENERATE_GAS_IN_ICS
      if(RestartFlag == 0)
	{
	  if(ThisTask == 0)
	    {
	      printf("  Generating space for GAS particles ...\n");
	      fflush(stdout);
	    }
	  All.TotN_gas += header.npartTotal[1];
	  All.TotN_gas += (((long long) header.npartTotalHighWord[1]) << 32);
	  All.TotNumPart += header.npartTotal[1];
	  All.TotNumPart += (((long long) header.npartTotalHighWord[1]) << 32);
	}
#endif


      for(i = 0; i < 6; i++)
	All.MassTable[i] = header.mass[i];

      All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));	/* sets the maximum number of particles that may */
      All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));	/* sets the maximum number of particles that may
										   reside on a processor */
      if(ThisTask == 0)
	{
	  printf("Task %d/%d/%d: MaxPart %d, MaxSphP %d (%lld,%lld)\n", ThisTask, readTask, lastTask,
		 All.MaxPart, All.MaxPartSph, All.TotNumPart, All.TotN_gas);
	  fflush(stdout);
	}

#ifdef INHOMOG_GASDISTR_HINT
      All.MaxPartSph = All.MaxPart;
#endif

#ifdef LT_STELLAREVOLUTION
      if(All.TotN_stars == 0)
	All.MaxPartMet = All.PartAllocFactor * (All.TotN_gas / NTask) * All.SFfactor * All.Generations;
      else
	All.MaxPartMet =
	  All.PartAllocFactor * (All.TotN_stars / NTask +
				 (All.TotN_gas / NTask) * All.SFfactor * All.Generations);
      if(ThisTask == 0)
	{
	  printf("       MaxStarPart %d (%lld,%lld,%g,%d)\n", All.MaxPartMet, All.TotN_gas, All.TotN_stars,
		 All.SFfactor, All.Generations);
	  fflush(stdout);
	}
#endif

#if defined(BLACK_HOLES)
      if(All.TotBHs == 0)
	All.MaxPartBH = All.PartAllocFactor * (All.TotN_gas / NTask) * All.BHfactor;
      else
	All.MaxPartBH = All.PartAllocFactor * (All.TotBHs / NTask + (All.TotN_gas / NTask) * All.BHfactor);
      if(ThisTask == 0)
	{
	  printf("       MaxBlackHolePart %d (%lld,%d,%g)\n", All.MaxPartBH, All.TotN_gas, All.TotBHs,
		 All.BHfactor);
	  fflush(stdout);
	}
#endif

      allocate_memory();

      if(!(CommBuffer = mymalloc("CommBuffer", bytes = All.BufferSize * 1024 * 1024)))
	{
	  printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}

      if(RestartFlag >= 2)
	{
	  All.Time = All.TimeBegin = header.time;
	  set_cosmo_factors_for_current_time();
	}

#ifdef LT_STELLAREVOLUTION

#ifndef GADGET3_IO_LIB
      time_age = get_age(All.Time);
#endif

#endif
    }

  if(ThisTask == readTask)
    {
      for(i = 0, n_in_file = 0; i < 6; i++)
	n_in_file += header.npart[i];

      if(ThisTask == 0)
	printf("Total Number of Particles: %llu\n",
	       header.npartTotal[0] + header.npartTotal[1] + header.npartTotal[2] +
	       header.npartTotal[3] + header.npartTotal[4] + header.npartTotal[5] +
	       (((long long) header.npartTotalHighWord[0]) << 32) +
	       (((long long) header.npartTotalHighWord[1]) << 32) +
	       (((long long) header.npartTotalHighWord[2]) << 32) +
	       (((long long) header.npartTotalHighWord[3]) << 32) +
	       (((long long) header.npartTotalHighWord[4]) << 32) +
	       (((long long) header.npartTotalHighWord[5]) << 32));
      printf("\nreading file `%s' on task=%d (contains %d particles.)\n"
	     "distributing this file to tasks %d-%d\n"
	     "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 1 (halo):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 2 (disk):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 3 (bulge): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 4 (stars): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 5 (bndry): %8d  (tot=%6d%09d) masstab=%g\n\n", fname, ThisTask, n_in_file, readTask,
	     lastTask, header.npart[0], (int) (header.npartTotal[0] / 1000000000),
	     (int) (header.npartTotal[0] % 1000000000), All.MassTable[0], header.npart[1],
	     (int) (header.npartTotal[1] / 1000000000), (int) (header.npartTotal[1] % 1000000000),
	     All.MassTable[1], header.npart[2], (int) (header.npartTotal[2] / 1000000000),
	     (int) (header.npartTotal[2] % 1000000000), All.MassTable[2], header.npart[3],
	     (int) (header.npartTotal[3] / 1000000000), (int) (header.npartTotal[3] % 1000000000),
	     All.MassTable[3], header.npart[4], (int) (header.npartTotal[4] / 1000000000),
	     (int) (header.npartTotal[4] % 1000000000), All.MassTable[4], header.npart[5],
	     (int) (header.npartTotal[5] / 1000000000), (int) (header.npartTotal[5] % 1000000000),
	     All.MassTable[5]);
      fflush(stdout);
    }


  ntask = lastTask - readTask + 1;


  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  for(type = 0, nall = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;


      if(type == 0)
	{
	  if(N_gas + n_for_this_task > All.MaxPartSph)
	    {
	      printf("Not enough space on task=%d for SPH particles (space for %d, need at least %d)\n",
		     ThisTask, All.MaxPartSph, N_gas + n_for_this_task);
	      fflush(stdout);
	      endrun(172);
	    }
	}

      nall += n_for_this_task;
    }

  if(NumPart + nall > All.MaxPart)
    {
      printf("Not enough space on task=%d (space for %d, need at least %d)\n",
	     ThisTask, All.MaxPart, NumPart + nall);
      fflush(stdout);
      endrun(173);
    }

  // This invalidates the pointers where particles are linked to their extra physics data field
  // Therefore this breaks down when reading multiple files to one MPI rank.
  // This is now done by bookkeeping offsets for all components and sorting gas particles back
  // to the front at the end of teh reading
  //  memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(struct particle_data));


  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

#ifdef GADGET_MAX_HSML
	  if (blocknr == IO_HSML)
		  break;
#endif

      if(blocknr == IO_LASTENTRY)
          break;

#ifndef GADGET_READ_ID
	  if (blocknr == IO_ID) /* skip for LONGIDS depends*/
		  continue;
#endif

#ifdef NO_UTHERM_IN_IC_FILE
      if(RestartFlag == 0 && blocknr == IO_U)
	continue;
#endif

#ifdef GM_STARDENSITY		/* no need to read these in input */
      if(blocknr == IO_RHOSTAR || blocknr == IO_HSMLSTAR || blocknr == IO_NEIGHSTAR)
	continue;
#endif


      if(RestartFlag == 5 && blocknr > IO_MASS)	/* if we only do power spectra, we don't need to read other blocks beyond the mass */
	continue;

      if(blockpresent(blocknr))
	{

#ifndef CHEMISTRY
	  if(RestartFlag == 0 && blocknr > IO_U && blocknr != IO_BFLD
#ifdef READ_HSML
	     && blocknr != IO_HSML
#endif
#ifdef LMB_SPECTRAL_CRs_FROM_ICS
	     && blocknr != IO_LMBCR_pNORM && blocknr != IO_LMBCR_eNORM
	     && blocknr != IO_LMBCR_pSLOPE && blocknr != IO_LMBCR_eSLOPE
	     && blocknr != IO_LMBCR_pCUT && blocknr != IO_LMBCR_eCUT
#endif // LMB_SPECTRAL_CRs_FROM_ICS
#ifdef READ_RHO_OLD
	     && blocknr != IO_RHO_OLD
#endif
#ifdef READ_MACH
	     && blocknr != IO_MACH
#endif // READ_MACH
#ifdef READ_EULER
	     && blocknr != IO_EULERB && blocknr != IO_EULERA
#endif
#ifdef LT_Zs_IN_IC
#ifdef LT_STELLAREVOLUTION
	     && blocknr != IO_Zs
#endif
#endif
	    )
#else
	  if(RestartFlag == 0 && blocknr > IO_HM)
#endif

#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
	    if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_SHEET_ORIENTATION))
	      if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_INIT_DENSITY))
		if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_CAUSTIC_COUNTER))
#endif
		  continue;	/* ignore all other blocks in initial conditions */


#ifdef SUBFIND_RESHUFFLE_AND_POTENTIAL
	  if(blocknr == IO_POT)
	    continue;
#endif


#ifdef BINISET
	  if(RestartFlag == 0 && blocknr == IO_BFLD)
	    continue;
#endif

#ifdef OUTPUTACCELERATION
	  if(RestartFlag == 2 && blocknr == IO_ACCEL)
	    continue;
#endif

	  if(RestartFlag == 2 && blocknr == IO_DELAYTIME)
	    continue;

#ifdef OUTPUTCOOLRATE
	  if(RestartFlag == 2 && blocknr == IO_COOLRATE)
	    continue;
#endif

#ifdef SUBFIND
	  if(RestartFlag == 2 && blocknr == IO_HSMS)
	    continue;
#endif

#ifdef LT_STELLAREVOLUTION
	  if(blocknr == IO_Zs)
	    Zs_present = 1;
#ifdef LT_ZAGE
	  if(blocknr == IO_ZAGE)
	    ZAge_present = 1;
#endif
#ifdef LT_ZAGE_LLV
	  if(blocknr == IO_ZAGE_LLV)
	    ZAge_llv_present = 1;
#endif
#ifdef LT_TRACK_CONTRIBUTES
	  if(blocknr == IO_CONTRIB && FdTrck == 0x0)
	    {
	      if(ThisTask == 0)
		printf("auxiliary file for block %d not present\n", blocknr);
	      continue;
	    }
#endif
#endif

#ifdef GM_MUPPI
	  if((blocknr == IO_FB_EOUT || blocknr == IO_FB_EREC || blocknr == IO_FB_CLOCK ||
	      blocknr == IO_FB_E_TOT_0 || blocknr == IO_FB_TDYN ||
	      blocknr == IO_FB_TCOOL || blocknr == IO_FB_TSTARTMP || blocknr == IO_FB_EKINREC)
	     && FdMuppi == 0x0)
	    {
	      if(ThisTask == 0)
		printf("auxiliary file for block %d not present\n", blocknr);
	      fflush(stdout);
	      continue;
	    }
#endif


#ifdef ADAPTGRAVSOFT
	  if(blocknr == IO_AGS_SOFT)
	    continue;
	  if(blocknr == IO_AGS_DENS)
	    continue;
	  if(blocknr == IO_AGS_ZETA)
	    continue;
	  if(blocknr == IO_AGS_OMEGA)
	    continue;
	  if(blocknr == IO_AGS_NGBS)
	    continue;
	  if(blocknr == IO_AGS_CORR)
	    continue;
#endif

#if defined (UM_CHEMISTRY) && defined (UM_CHEMISTRY_INISET)
	  if(RestartFlag == 0 && blocknr == IO_NE)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_NH)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_HII)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_HM)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_HeI)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_HeII)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_HeIII)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_H2I)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_H2II)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_HD)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_DI)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_DII)
	    continue;
	  if(RestartFlag == 0 && blocknr == IO_HeHII)
	    continue;
#endif

	  if(blocknr == IO_HSMS)
	    continue;

	  if(ThisTask == readTask && (ThisTask == 0 || VERBOSE_LEVEL > 2))
	    {
	      get_dataset_name(blocknr, buf);
	      printf("Task %d: reading block %d (%s) \n", ThisTask, bnr, buf);
	      fflush(stdout);
	    }

	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 1);

	  blockmaxlen = (size_t) ((All.BufferSize * 1024 * 1024) / bytes_per_blockelement);

	  npart = get_particles_in_block(blocknr, &typelist[0]);

#ifdef GADGET3_IO_LIB
	  if(npart > 0)
	    {
            // 6 particle types, max. 1000 block types
            all_particles_blocks[6 * blocknr + 0] += typelist[0];
            all_particles_blocks[6 * blocknr + 1] += typelist[1];
            all_particles_blocks[6 * blocknr + 2] += typelist[2];
            all_particles_blocks[6 * blocknr + 3] += typelist[3];
            all_particles_blocks[6 * blocknr + 4] += typelist[4];
            all_particles_blocks[6 * blocknr + 5] += typelist[5];
		}
#endif

	  if(npart > 0)
	    {
	      if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP)
		if(ThisTask == readTask)
		  {
#if defined(LT_TRACK_CONTRIBUTES)
		    if(blocknr == IO_CONTRIB)
		      fd = FdTrck;
#endif
#ifdef GM_MUPPI
		    if(blocknr == IO_FB_EOUT || blocknr == IO_FB_EREC || blocknr == IO_FB_CLOCK ||
		       blocknr == IO_FB_E_TOT_0 || blocknr == IO_FB_TDYN ||
		       blocknr == IO_FB_TCOOL || blocknr == IO_FB_TSTARTMP || blocknr == IO_FB_EKINREC)
		      fd = FdMuppi;
#endif


		    if(All.ICFormat == 2)
		      {
			get_Tab_IO_Label(blocknr, label);
			find_block(label, fd);
		      }

		    if(All.ICFormat == 1 || All.ICFormat == 2)
		      SKIP;
		  }

	      // Here we re-set the different offsets to the previously read number of particles for each block
	      for(type = 0, offset = NumPart, offset_gas = N_gas, offset_stars = N_stars, offset_bhs =
		  N_BHs, nread = 0; type < 6; type++)
		{
		  n_in_file = header.npart[type];
#ifdef HAVE_HDF5
		  pcsum = 0;
#endif
		  if(typelist[type] == 0)
		    {
		      n_for_this_task = n_in_file / ntask;
		      if((ThisTask - readTask) < (n_in_file % ntask))
			n_for_this_task++;

		      offset += n_for_this_task;
		    }
		  else
		    {
		      for(task = readTask; task <= lastTask; task++)
			{
			  n_for_this_task = n_in_file / ntask;
			  if((task - readTask) < (n_in_file % ntask))
			    n_for_this_task++;

			  if(task == ThisTask)
			    if(NumPart + n_for_this_task > All.MaxPart)
			      {
				printf("too many particles. %d %d %d\n", NumPart, n_for_this_task,
				       All.MaxPart);
				endrun(1313);
			      }

			  if(ThisTask == readTask && blocknr == IO_POS
			     && (ThisTask == 0 || VERBOSE_LEVEL > 2))
			    {
			      printf("Task %d:    type = %d, offset = %d, n_for_this_task = %d\n",
				     ThisTask, type, offset, n_for_this_task);
			      if(type == 0)
				printf("Task %d:    offset_gas = %d, N_gas = %d\n",
				       ThisTask, offset_gas, N_gas);
			      if(type == 1)
				printf("Task %d:    offset = %d, NumPart = %d\n", ThisTask, offset, NumPart);
#ifdef LT_STELLAREVOLUTION
			      if(type == 4)
				printf("Task %d:    offset_stars = %d, N_stars = %d\n",
				       ThisTask, offset_stars, N_stars);
#endif
#ifdef BLACK_HOLES
			      if(type == 5)
				printf("Task %d:    offset_bhs = %d, N_BHs = %d\n",
				       ThisTask, offset_bhs, N_BHs);
#endif
			    }
			  // Note that we shop here the reading due to teh buffer size, so offsets have to be increased
			  // after each segment which was read
			  do
			    {
			      pc = n_for_this_task;

			      if(pc > (int) blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == readTask)
				{
				  if(All.ICFormat == 1 || All.ICFormat == 2)
				    {
				      if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY
					 && blocknr != IO_DMVELDISP)
					{
					  my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
					  nread += pc;
					}
				      else
					{
					  nread += pc;
					}
				    }

#ifdef HAVE_HDF5
				  if(All.ICFormat == 3 && pc > 0)
				    {
				      get_dataset_name(blocknr, buf);
#ifdef OLD_HDF5
				      hdf5_dataset = H5Dopen(hdf5_grp[type], buf);
#else
				      hdf5_dataset = H5Dopen(hdf5_grp[type], buf, H5P_DEFAULT);
#endif
				      dims[0] = header.npart[type];
				      dims[1] = get_values_per_blockelement(blocknr);
				      if(dims[1] == 1)
					rank = 1;
				      else
					rank = 2;

				      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

				      dims[0] = pc;
				      hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

				      start[0] = pcsum;
				      start[1] = 0;

				      count[0] = pc;
				      count[1] = get_values_per_blockelement(blocknr);
				      pcsum += pc;

				      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							  start, NULL, count, NULL);

				      switch (get_datatype_in_block(blocknr))
					{
					case 0:
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
					  break;
					case 1:
#ifdef INPUT_IN_DOUBLEPRECISION
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
					  break;
					case 2:
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
					  break;
					}

				      H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory,
					      hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

				      H5Tclose(hdf5_datatype);
				      H5Sclose(hdf5_dataspace_in_memory);
				      H5Sclose(hdf5_dataspace_in_file);
				      H5Dclose(hdf5_dataset);
				    }
#endif
				}

			      if(ThisTask == readTask && task != readTask && pc > 0)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					  TAG_PDATA, MYMPI_COMM_WORLD);

			      if(ThisTask != readTask && task == ThisTask && pc > 0)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask,
					 TAG_PDATA, MYMPI_COMM_WORLD, &status);

			      if(ThisTask == task)
				{
				  empty_read_buffer(blocknr, offset, pc, type, offset_gas, offset_stars,
						    offset_bhs);
				  // Now updating the offsets for the different particle types after (partially) reading a block
				  offset += pc;
				  if(type == 0)
				    offset_gas += pc;
				  if(type == 4)
				    offset_stars += pc;
				  if(type == 5)
				    offset_bhs += pc;
				}

			      n_for_this_task -= pc;
			    }
			  while(n_for_this_task > 0);
			}
		    }
		}

	      if(ThisTask == readTask)
		{
		  if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP)
		    if(All.ICFormat == 1 || All.ICFormat == 2)
		      {
			SKIP2;

#ifdef GM_MUPPI
			if(blocknr == IO_FB_EOUT || blocknr == IO_FB_EREC || blocknr == IO_FB_CLOCK ||
			   blocknr == IO_FB_E_TOT_0 || blocknr == IO_FB_TDYN ||
			   blocknr == IO_FB_TCOOL || blocknr == IO_FB_TSTARTMP || blocknr == IO_FB_EKINREC)
			  {
			    fd = fd_muppi_back;
			    if(FdMuppi != 0x0)
			      fclose(FdMuppi);
			  }
#endif


#ifdef AUTO_SWAP_ENDIAN_READIC
			swap_Nbyte((char *) &blksize1, 1, 4);
			swap_Nbyte((char *) &blksize2, 1, 4);
#endif
			if(blksize1 != blksize2)
			  {
			    printf("incorrect block-sizes detected!\n");
			    printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, bnr,
				   blksize1, blksize2);
			    if(blocknr == IO_ID)
			      {
				printf
				  ("Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
			      }
			    fflush(stdout);
			    endrun(1889);
			  }
		      }
		}
	    }
	}
    }

  // Now updating the current number of read particles of teh different species
  for(type = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;

      NumPart += n_for_this_task;

      if(type == 0)
	N_gas += n_for_this_task;
#ifdef LT_STELLAREVOLUTION
      if(type == 4)
	N_stars += n_for_this_task;
#endif
#if defined(BLACK_HOLES)
      if(type == 5)
	N_BHs += n_for_this_task;
#endif
    }

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	fclose(fd);
#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
	{
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Fclose(hdf5_file);
	}
#endif
    }
}



/*! This function determines on how many files a given snapshot is distributed.
 */
int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if(All.ICFormat == 3)
    {
      sprintf(buf, "%s.%d.hdf5", fname, 0);
      sprintf(buf1, "%s.hdf5", fname);
    }

  if(ThisTask == 0)
    printf("Trying to read file %s in format %d\n", buf, All.ICFormat);

#ifndef  HAVE_HDF5
  if(All.ICFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }
#endif

  header.num_files = 0;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "rb")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.ICFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = 1;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
#ifdef PATCH_IO
	      printf("Task %d: WARNING, PATCHING HEADER FOR MagneticumDM simulation !!!\n", ThisTask);
	      header.mass[2] = header.mass[1];
	      header.mass[1] = header.mass[0];
	      header.mass[0] = 0;
	      header.npart[2] = header.npart[1];
	      header.npart[1] = header.npart[0];
	      header.npart[0] = 0;
	      header.npartTotal[2] = header.npartTotal[1];
	      header.npartTotal[1] = header.npartTotal[0];
	      header.npartTotal[0] = 0;
	      header.npartTotalHighWord[2] = header.npartTotalHighWord[1];
	      header.npartTotalHighWord[1] = header.npartTotalHighWord[0];
	      header.npartTotalHighWord[0] = 0;
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(All.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf);
#endif
	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, 1, MPI_INT, 0, MYMPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf1, "rb")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.ICFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = 1;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
#ifdef PATCH_IO
	      printf("Task %d: WARNING, PATCHING HEADER FOR MagneticumDM simulation !!!\n", ThisTask);
	      header.mass[2] = header.mass[1];
	      header.mass[1] = header.mass[0];
	      header.mass[0] = 0;
	      header.npart[2] = header.npart[1];
	      header.npart[1] = header.npart[0];
	      header.npart[0] = 0;
	      header.npartTotal[2] = header.npartTotal[1];
	      header.npartTotal[1] = header.npartTotal[0];
	      header.npartTotal[0] = 0;
	      header.npartTotalHighWord[2] = header.npartTotalHighWord[1];
	      header.npartTotalHighWord[1] = header.npartTotalHighWord[0];
	      header.npartTotalHighWord[0] = 0;
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(All.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf1);
#endif

	  header.num_files = 1;
	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, 1, MPI_INT, 0, MYMPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      printf("\nCan't find initial conditions file.");
      printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
      fflush(stdout);
    }

  endrun(0);
  return 0;
}


/*! This function assigns a certain number of files to processors, such that
 *  each processor is exactly assigned to one file, and the number of cpus per
 *  file is as homogenous as possible. The number of files may at most be
 *  equal to the number of processors.
 */

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last)
{
  int ntask, filesleft, filesright, tasksleft;

  if(nfiles > 1)
    {
      ntask = lasttask - firsttask + 1;

      filesleft = (int) ((((double) (ntask / 2)) / ntask) * nfiles);
      if(filesleft <= 0)
	filesleft = 1;
      if(filesleft >= nfiles)
	filesleft = nfiles - 1;

      filesright = nfiles - filesleft;

      tasksleft = ntask / 2;
      distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
      distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr, master,
		      last);
    }
  else
    {
      if(ThisTask >= firsttask && ThisTask <= lasttask)
	{
	  *filenr = firstfile;
	  *master = firsttask;
	  *last = lasttask;
	}
    }
}


#ifdef HAVE_HDF5
void read_header_attributes_in_hdf5(char *fname)
{
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
  int i;

  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#ifdef OLD_HDF5
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
#else
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header", H5P_DEFAULT);
#endif
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);

#ifdef LT_STELLAREVOLUTION
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
#endif

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_IC_Info");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
  H5Aclose(hdf5_attribute);


  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);
}
#endif





#ifdef AUTO_SWAP_ENDIAN_READIC
/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data, int n, int m)
{
  int i, j;
  char old_data[16];

  if(swap_file != 8)
    {
      for(j = 0; j < n; j++)
	{
	  memcpy(&old_data[0], &data[j * m], m);
	  for(i = 0; i < m; i++)
	    {
	      data[j * m + i] = old_data[m - i - 1];
	    }
	}
    }
}

/*------------------------------------------------------------------*/
/*----------- procedure to swap header if needed -------------------*/
/*------------------------------------------------------------------*/

void swap_header()
{
  swap_Nbyte((char *) &header.npart, 6, 4);
  swap_Nbyte((char *) &header.mass, 6, 8);
  swap_Nbyte((char *) &header.time, 1, 8);
  swap_Nbyte((char *) &header.redshift, 1, 8);
  swap_Nbyte((char *) &header.flag_sfr, 1, 4);
  swap_Nbyte((char *) &header.flag_feedback, 1, 4);
  swap_Nbyte((char *) &header.npartTotal, 6, 4);
  swap_Nbyte((char *) &header.flag_cooling, 1, 4);
  swap_Nbyte((char *) &header.num_files, 1, 4);
  swap_Nbyte((char *) &header.BoxSize, 1, 8);
  swap_Nbyte((char *) &header.Omega0, 1, 8);
  swap_Nbyte((char *) &header.OmegaLambda, 1, 8);
  swap_Nbyte((char *) &header.HubbleParam, 1, 8);
  swap_Nbyte((char *) &header.flag_stellarage, 1, 4);
  swap_Nbyte((char *) &header.flag_metals, 1, 4);
  swap_Nbyte((char *) &header.npartTotalHighWord, 6, 4);
  swap_Nbyte((char *) &header.flag_entropy_instead_u, 1, 4);
  swap_Nbyte((char *) &header.flag_doubleprecision, 1, 4);
  swap_Nbyte((char *) &header.flag_ic_info, 1, 4);
  swap_Nbyte((char *) &header.lpt_scalingfactor, 1, 4);
}

#endif

/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
void find_block(char *label, FILE *fd)
{
  unsigned int blocksize = 0, blksize;
  char blocklabel[5] = { "    " };

#define FBSKIP  {my_fread(&blksize,sizeof(int),1,fd);}

  rewind(fd);

#if defined(LT_TRACK_CONTRIBUTES)
  int blcksize;
  if(strcmp(label, "TRCK") == 0)
    {
	  my_fread(&blcksize, sizeof(int), 1, fd);
      blcksize += sizeof(int);
#ifdef _WIN32
	  _fseeki64(fd, (__int64)blcksize, SEEK_CUR);
#else
      fseek(fd, blcksize, SEEK_CUR);
#endif
    }
#endif

  while(!feof(fd) && blocksize == 0)
    {
      FBSKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
      swap_file = blksize;
      swap_Nbyte((char *) &blksize, 1, 4);
#endif
      if(blksize != 8)
	{
	  printf("Incorrect Format (blksize=%u)!\n", blksize);
	  exit(1891);
	}
      else
	{
	  my_fread(blocklabel, 4 * sizeof(char), 1, fd);
	  my_fread(&blocksize, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_Nbyte((char *) &blocksize, 1, 4);
#endif
	  /*
	     printf("Searching <%c%c%c%c>, found Block <%s> with %d bytes\n",
	     label[0],label[1],label[2],label[3],blocklabel,blocksize);
	   */
	  FBSKIP;
	  if(strncmp(label, blocklabel, 4) != 0)
	    {
#ifdef _WIN32
		  _fseeki64(fd, (__int64)blocksize, 1);
#else
	      fseek(fd, blocksize, 1);
#endif
	      blocksize = 0;
	    }
	}
    }
  if(feof(fd))
    {
      printf("Block '%c%c%c%c' not found !\n", label[0], label[1], label[2], label[3]);
      fflush(stdout);
      endrun(1890);
    }
}
