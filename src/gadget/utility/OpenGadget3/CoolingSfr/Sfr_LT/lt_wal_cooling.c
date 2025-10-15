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
#include <dirent.h>
#include <math.h>
#include <ctype.h>
#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"

#ifdef LT_METAL_COOLING_WAL
#define H5Dopen_vers 2

#define H5Acreate_vers 2

#define H5Dcreate_vers 2
#include <hdf5.h>

#include "lt_wal_cooling.h"	/* this includes the variables
				   needed to use Wiersma et al
				   cooling */
#include "lt_wal_cooling_nh.h"	/* this includes variables
				   needed to calculate H ionization
				   using a pure H+He ioniz. equ.
				   with the UV background */
#include "lt_error_codes.h"


/* :: -------------------------------------------- ::  */
/*    SEGMENT .CODE                                    */



/* :: ---------------------------- PUBLIC ROUTINES ::  */
/*                                 MAIN ROUTINES       */
/*                                                     */
/* +  WalCoolInitialize(void)                          */
/*     initializes tables and memory                   */
/*                                                     */
/* +  WalCool_tables_load()                            */
/*     checks whether or not a new table(s) must be    */
/*     loaded and loads them.                          */
/*                                                     */
/* +  *set_metallicities(int, double *, double)        */
/*     set up an array with the metallicities of those */
/*     elements that are in the cooling matrix         */
/*                                                     */
/*                                                     */

/*                                 USEFUL ROUTINES     */
/*                                                     */
/* +  set_cooltable_index(int)                         */
/*     set by hands the index of the current low-z     */
/*     cooling table                                   */
/*                                                     */
/* +  get_max_cool_redshift()                          */
/*     returns the highest cooling redshift            */
/*                                                     */
/* +  GetCoolingTime(..)                               */
/*     returnd the cooling time                        */
/*                                                     */
/* +  DoCooling(..)                                    */
/*     performs the cooling                            */
/*                                                     */
/* +  convert_u_to_temp(double, double, double*)       */
/*     returnd the temperature                         */
/*                                                     */
/* 
*/

int DBG = 0;
/* static */

void read_cooling_tables_dummy()
{
  int i;

  ZBins = 1;
  TBins = 1;
  ZMin = ZMax = -4.0;
//  TMin = 4.0;
//  TMax = 8.0;
  cl_set_Tmin(4.0);
  cl_set_Tmax(8.0);

  if(ThisTask == 0)
    printf("Setting dummy cooling tables to zero for WAL cooling ...\n");

  CoolTvalue = (double *) mymalloc("CoolTvalue", (TBins + ZBins + TBins * ZBins) * sizeof(double));
  memset(CoolTvalue, 0, (TBins + ZBins + TBins * ZBins) * sizeof(double));

  CoolingTables = (double **) mymalloc("CoolingTables", ZBins * sizeof(double *));
  memset(CoolingTables, 0, ZBins * sizeof(double *));

  CoolZvalue = CoolTvalue + TBins;
  CoolingTables[0] = CoolZvalue + ZBins;
  for(i = 1; i < ZBins; i++)
    CoolingTables[i] = CoolingTables[i - 1] + TBins;

}

void WalCool_Initialize()
/* initialize the cooling: loads the common tables and allocate memory to store 2 cooling tables */
{

  if(ThisTask == 0)
    {
      WalCool_get_redshift_table();
      printf("\n[Wiersma et al. Cooling] %d tables found, spanning redshift from %g downto %g %g\n\n",
	     WalCool_CoolTables_num, WalCool_CoolTables_redshifts[WalCool_CoolTables_num - 1].redshift,
	     WalCool_CoolTables_redshifts[0].redshift, max_cool_redshift);
    }

  MPI_Bcast(&WalCool_CoolTables_num, sizeof(int), MPI_BYTE, 0, MYMPI_COMM_WORLD);
  MPI_Bcast(&max_cool_redshift, sizeof(double), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  if(ThisTask != 0)
    WalCool_CoolTables_redshifts =
      (struct my_direntry *) mymalloc("WalCool_CoolTables_redshifts",
				      sizeof(struct my_direntry) * WalCool_CoolTables_num);

  MPI_Bcast(WalCool_CoolTables_redshifts, sizeof(struct my_direntry) * WalCool_CoolTables_num, MPI_BYTE, 0,
	    MYMPI_COMM_WORLD);

  /* get general informations from headers */
  WalCool_Initialize_get_header();

  /* allocate space for tables */
  WalCoolMfreeSize = WalCool_n_Hef * WalCool_n_Rho * WalCool_n_T;
  WalCoolSize = WalCool_n_El * WalCool_n_Rho * WalCool_n_T;

  WALCOOLTABLES = (float *) mymalloc("WALCOOLTABLES", sizeof(float) * 2 * (WalCoolMfreeSize + WalCoolSize));

  /* set pointers to H+He cooling */
  WalMfreeCoolTables = WALCOOLTABLES;

  /* set pointer to metal cooling tables */
  WalCoolTables = WALCOOLTABLES + 2 * WalCoolMfreeSize;


  WalCool_redshift_index = -1;
  HighZTable = 0;
  LowZTable = 1;

  return;
}


void WalCool_tables_load(double redshift)
{
  int change = 0;

  if(redshift > max_cool_redshift)
    {
      if(collisional_table_loaded)
	return;

      WalCool_redshift_index = -1;
      WalCool_get_collis_table();
      collisional_table_loaded = 1;
      if(ThisTask == 0)
	printf("[WalCool] getting collisional table\n");
      fflush(stdout);

    }
  else
    {
      if(WalCool_redshift_index == -1)
	{
	  WalCool_redshift_index = WalCool_CoolTables_num - 1;
	  change = 2;
	}
      else if(WalCool_redshift_index == 0)
	return;

      for(;
	  (redshift < WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift)
	  && (WalCool_redshift_index >= 0); WalCool_redshift_index--)
	change++;

      if(change == 1)
	{

	  /*
	     HighZTable ^= 1;
	     LowZTable ^= 1;
	   */

	  if(ThisTask == 0)
	    printf("[WalCool] sup = %d, inf = %d, getting %d in inf, z=%f,having %d in sup, z=%f\n",
		   HighZTable, LowZTable, WalCool_redshift_index,
		   WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift, WalCool_redshift_index + 1,
		   WalCool_CoolTables_redshifts[WalCool_redshift_index + 1].redshift);

	  WalCool_get_table(WalCool_redshift_index, LowZTable);
	  WalCool_get_table(WalCool_redshift_index + 1, HighZTable);

	  MPI_Bcast(&WalMfreeCoolTables[MFREE_COOLRATE_z_IDX(LowZTable)], WalCoolMfreeSize * sizeof(float),
		    MPI_BYTE, 0, MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCoolTables[COOLRATE_z_IDX(LowZTable)], WalCoolSize * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_enHS[ENHS_z_IDX(LowZTable)], WalCool_arraysS_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_enH[ENH_z_IDX(LowZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_mu[MU_z_IDX(LowZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_UtoT[UtoT_z_IDX(LowZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);

	  MPI_Bcast(&WalMfreeCoolTables[MFREE_COOLRATE_z_IDX(HighZTable)], WalCoolMfreeSize * sizeof(float),
		    MPI_BYTE, 0, MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCoolTables[COOLRATE_z_IDX(HighZTable)], WalCoolSize * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_enHS[ENHS_z_IDX(HighZTable)], WalCool_arraysS_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_enH[ENH_z_IDX(HighZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_mu[MU_z_IDX(HighZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);
	  MPI_Bcast(&WalCool_UtoT[UtoT_z_IDX(HighZTable)], WalCool_arrays_size * sizeof(float), MPI_BYTE, 0,
		    MYMPI_COMM_WORLD);

	  return;
	}
      else if(change > 1)
	{
	  if(ThisTask == 0)
	    printf("[WalCool] sup = %d, inf = %d, getting %d(z=%f) in inf, %d(z=%f) in sup\n", HighZTable,
		   LowZTable, WalCool_redshift_index,
		   WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift, WalCool_redshift_index + 1,
		   WalCool_CoolTables_redshifts[WalCool_redshift_index + 1].redshift);

	  WalCool_get_table(WalCool_redshift_index, LowZTable);
	  WalCool_get_table(WalCool_redshift_index + 1, HighZTable);
	}
    }

  MPI_Bcast(WALCOOLTABLES, (2 * (WalCoolMfreeSize + WalCoolSize)) * sizeof(float), MPI_BYTE, 0,
	    MYMPI_COMM_WORLD);
  MPI_Bcast(WalCool_enHS, (3 * 2 * WalCool_arrays_size + 2 * WalCool_arraysS_size) * sizeof(float), MPI_BYTE,
	    0, MYMPI_COMM_WORLD);

  return;
}


// legacy version
int get_cool_redshift(double Redshift, double *DZ)
{
  /* find index for redshift */
  if(WalCool_redshift_index >= 0)
    *DZ =
      (Redshift -
       WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift) /
      (WalCool_CoolTables_redshifts[WalCool_redshift_index + 1].redshift -
       WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift);
  else
    *DZ = 0;
  return 0;
}


// Updated version. Returns value directly for readability
double lt_cl_GetCoolRedshift(double Redshift)
{
  double DZ;			// TODO what is DZ, anyway?
  if(WalCool_redshift_index >= 0)
    {
      int z0 = WalCool_CoolTables_redshifts[WalCool_redshift_index].redshift;
      int z1 = WalCool_CoolTables_redshifts[WalCool_redshift_index + 1].redshift;

      DZ = (Redshift - z0) / (z1 - z0);
    }
  else
    {
      DZ = 0;
    }

  return DZ;
}

int get_cool_n_el()
{
  return WalCool_n_Ab;
}

int Is_a_Coolant(int i)
{
  return Specie_is_coolant[i];
}


void set_cooltable_index(int idx)
{
  WalCool_redshift_index = idx;
  return;
}


double get_max_cool_redshift()
{
  return max_cool_redshift;
}


double *set_metallicities(int i, double *Metallicities, double dens_to_phys_fact)
/*
 * dens_to_phys_fact must contain all factors needed to convert densities
 * from code units to physical code units (also cosmological factors must be
 * included)
 */
{
  double nonHydmass;
  MyDouble mass;
  int pos, j, k;

  dens_to_phys_fact *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  mass = P[i].Mass;
  if(P[i].Type != 0)
    {
      printf("[%d] : error : request to calculate metallicities for particle"
	     " %d that is not a gas particle (type %d)\n", ThisTask, i, P[i].Type);
      endrun(LT_ERR_WALCOOL_Z_FOR_NOTGAS);
    }

  for(pos = 0, j = 0; j < WalCool_n_Ab; j++, pos++)
    if(SpeciesIdx[j] >= 0)
      {
	if(j == myHyd)
	  {
	    for(nonHydmass = 0, k = 0; k < LT_NMetP; k++)
#ifdef GL_CR_DUST
	      nonHydmass += SphP[i].Metals[k] + SphP[i].DustL[k] + SphP[i].DustS[k];
#else
	      nonHydmass += SphP[i].Metals[k];
#endif
	    Metallicities[SpeciesPos[j]] =
	      (double) (mass - nonHydmass) / mass * SphP[i].Density * dens_to_phys_fact;
	  }
	else
#if defined(LT_METAL_COOLING_on_SMOOTH_Z)
	  Metallicities[SpeciesPos[j]] = (double) SphP[i].Zsmooth[SpeciesIdx[j]];
#else
	  Metallicities[SpeciesPos[j]] = (double) SphP[i].Metals[SpeciesIdx[j]] / mass;
#endif

	if(j == myHel && UseHeNumberRatio)
	  Metallicities[SpeciesPos[j]] *= SphP[i].Density * dens_to_phys_fact;
      }
    else
      printf("%d element %d has no idx\n", ThisTask, j);

  Metallicities[HydPos] /= PROTONMASS;
  if(UseHeNumberRatio)
    Metallicities[HelPos] /= (4 * PROTONMASS * Metallicities[HydPos]);	/* it would be better to use number ratio, however for    
									 * some unknown reason, the hdf5 lib does not read in the 
									 * Helium_number_ratio_bins data set                      
									 */
  return Metallicities;
}

#ifdef SUBFIND

double *set_metallicities_subfind(int i, double *Metallicities, double dens_to_phys_fact)
/*
 * dens_to_phys_fact must contain all factors needed to convert densities
 * from code units to physical code units (also cosmological factors must be
 * included)
 */
{
  double nonHydmass;
  MyDouble mass;
  int pos, j, k;

  dens_to_phys_fact *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  mass = P[i].Mass;
  if(P[i].Type != 0)
    {
      printf("[%d] : error : request to calculate metallicities for particle"
	     " %d that is not a gas particle\n", ThisTask, i);
      endrun(LT_ERR_WALCOOL_Z_FOR_NOTGAS);
    }

  for(pos = 0, j = 0; j < WalCool_n_Ab; j++, pos++)
    if(SpeciesIdx[j] >= 0)
      {
	if(j == myHyd)
	  {
	    for(nonHydmass = 0, k = 0; k < LT_NMetP; k++)
#ifdef GM_EG_GL_DUST
	      nonHydmass +=
		SphP[P[i].origindex].Metals[k] + SphP[P[i].origindex].DustL[k] +
		SphP[P[i].origindex].DustS[k];
#else
	      nonHydmass += SphP[P[i].origindex].Metals[k];
#endif
	    Metallicities[SpeciesPos[j]] =
	      (double) (mass - nonHydmass) / mass * SphP[P[i].origindex].Density * dens_to_phys_fact;
	  }
	else
#if defined(LT_METAL_COOLING_on_SMOOTH_Z)
	  Metallicities[SpeciesPos[j]] = (double) SphP[P[i].origindex].Zsmooth[SpeciesIdx[j]];
#else
	  Metallicities[SpeciesPos[j]] = (double) SphP[P[i].origindex].Metals[SpeciesIdx[j]] / mass;
#endif

	if(j == myHel && UseHeNumberRatio)
	  Metallicities[SpeciesPos[j]] *= SphP[P[i].origindex].Density * dens_to_phys_fact;
      }
    else
      printf("%d element %d has no idx\n", ThisTask, j);

  Metallicities[HydPos] /= PROTONMASS;
  if(UseHeNumberRatio)
    Metallicities[HelPos] /= (4 * PROTONMASS * Metallicities[HydPos]);	/* it would be better to use number ratio, however for    
									 * some unknown reason, the hdf5 lib does not read in the 
									 * Helium_number_ratio_bins data set                      
									 */
  return Metallicities;
}

/* return the neutral hydrogen mass of particle i in internal units */
double get_HImass_subfind(int i, double Redshift, double DZ)
{
  double u, ne, mH, nh0 = 0, nHeII = 0;
  double fac_entr_to_u, a3inv;
  double rho, temp;
  double Metallicities[LT_NMet];

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    {
      a3inv = 1;
    }
  mH = P[i].Mass * HYDROGEN_MASSFRAC;
  //mH *= All.UnitMass_in_g / All.HubbleParam;

  rho = SphP[P[i].origindex].Density * a3inv;

  fac_entr_to_u = pow(rho, get_gamma_minus1(i)) / get_gamma_minus1(i);

  u = DMAX(All.MinEgySpec, SphP[P[i].origindex].Entropy * fac_entr_to_u);

  ne = SphP[P[i].origindex].elec;

  set_metallicities_subfind(i, &Metallicities[0], a3inv);

  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  temp = convert_u_to_tempMET(u, Redshift, DZ, &Metallicities[0]);

  nh0 = AbundanceRatios(temp, rho, &ne, &nh0, &nHeII);

  return nh0 * mH;

}

#endif

int set_indexes(double U, double DZ, double *Metallicities, double *dX, int *IDX)
/*
 * quantities must be given in physical units (Metallicities, for instance,
 * as returned by set_metallicities() )
 *
 * NOTE : the index U_IDX here is set referring to the INTERNAL ENERGY, not to
 *        the TEMPERATURE (you must do it after calling convert_u_to_temp, like
 *        in CoolingRateFromU, and in any case before calling CoolingRate
 */
{
  /* find index for redshift */
  dz = DZ;
  IDX[z_IDX] = WalCool_redshift_index;

  /* find index for hydrogen number density */
  IDX[H_IDX] = find_index(LOG, WalCool_Rho, Rmin, Rrange, WalCool_n_Rho, Metallicities[HydPos], &dX[H_IDX]);

  /* find index for the ration between he number density and H number density */
  IDX[He_IDX] =
    find_index(LIN, WalCool_Hef, WalCool_Hef[0], WalCool_Hef[WalCool_n_Hef - 1] - WalCool_Hef[0],
	       WalCool_n_Hef, Metallicities[HelPos], &dX[He_IDX]);

  /* find index for internal energy */
  IDX[U_IDX] = find_index(LOG, WalCool_U, Umin, Urange, WalCool_n_T, U, &dX[U_IDX]);

  return 0;
}

int set_indexes_T(double temp, double DZ, double *Metallicities, double *dX, int *IDX)
/*
 * quantities must be given in physical units (Metallicities, for instance,
 * as returned by set_metallicities() )
 *
 * NOTE : the index U_IDX here is set referring to the INTERNAL ENERGY, not to
 *        the TEMPERATURE (you must do it after calling convert_u_to_temp, like
 *        in CoolingRateFromU, and in any case before calling CoolingRate
 */
{
  /* find index for redshift */
  dz = DZ;
  IDX[z_IDX] = WalCool_redshift_index;

  /* find index for hydrogen number density */
  IDX[H_IDX] = find_index(LOG, WalCool_Rho, Rmin, Rrange, WalCool_n_Rho, Metallicities[HydPos], &dX[H_IDX]);

  /* find index for the ration between he number density and H number density */
  IDX[He_IDX] =
    find_index(LIN, WalCool_Hef, WalCool_Hef[0], WalCool_Hef[WalCool_n_Hef - 1] - WalCool_Hef[0],
	       WalCool_n_Hef, Metallicities[HelPos], &dX[He_IDX]);

  /* find index for Temperature */
  IDX[U_IDX] = find_index(LOG, WalCool_T, Tmin_wal, Trange, WalCool_n_T, (float) temp, &dX[U_IDX]);

  return 0;
}


#ifdef GL_DUST_COOLING
double GetCoolingTimeMET(double U, double Rho, double Redshift, double DZ, double DL, double DS,
			 double *Metallicities, double *temp)
/* quantities must be in physical code units
 */
{
  double Lambda, CTime, ratefact;
  /* convert in physical cgs units */
  Rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  U *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  ratefact = Metallicities[HydPos] * Metallicities[HydPos] / Rho;
  if((Lambda = CoolingRateFromU(U, Redshift, DZ, DL, DS, Metallicities, temp)) > 0)
    CTime = 0;
  else
    {
      /* note that if Lambda=0, CTime is formally Inf */
      CTime = U / (-ratefact * Lambda);
      CTime *= All.HubbleParam / All.UnitTime_in_s;
    }

  return CTime;
}
#else
double GetCoolingTimeMET(double U, double Rho, double Redshift, double DZ, double *Metallicities,
			 double *temp)
/* quantities must be in physical code units
 */
{
  double Lambda, CTime, ratefact;
  /* convert in physical cgs units */
  Rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  U *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  ratefact = Metallicities[HydPos] * Metallicities[HydPos] / Rho;

  if((Lambda = CoolingRateFromU(U, Redshift, DZ, Metallicities, temp)) > 0)
    CTime = 0;
  else
    {
      /* note that if Lambda=0, CTime is formally Inf */
      CTime = U / (-ratefact * Lambda);
      CTime *= All.HubbleParam / All.UnitTime_in_s;
    }

  return CTime;
}
#endif // GL_DUST_COOLING


double convert_u_to_tempMET(double U, double Redshift, double DZ, double *Metallicities)
/*
 * quantities must be given in physical units (Metallicities, for instance,
 * as returned by set_metallicities)
 *
 * Metallicities contains densitiy of metals as well as of H and He
 */
{
  double T, dX[4];
  int IDX[4];

  set_indexes(U, DZ, Metallicities, &dX[0], &IDX[0]);

  if(WalCool_redshift_index >= 0)
    {
/*       for(i = 0; i < (1 << 4); i++) */
/*         { */
/*           for(i = 0; i < 4; i++) */
/*             S[i] = j & (1 << i); */

/*           V[i] = WalCool_UtoT[UtoT_IDX(WalCool_redshift_index + S[0], */
/*                                        He_idx + S[1], */
/*                                        T_idx + S[2], */
/*                                        R_idx + S[3])]; */
/*         } */

/*       T = get_Ndim_interp(4, V, dC) */


      T = (1 - dz) * (1 - dHe) * (1 - dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU, inH)] +
	(dz) * (1 - dHe) * (1 - dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU, inH)] +
	(1 - dz) * (dHe) * (1 - dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU, inH)] +
	(dz) * (dHe) * (1 - dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU, inH)] +
	(1 - dz) * (1 - dHe) * (dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU + 1, inH)] +
	(dz) * (1 - dHe) * (dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU + 1, inH)] +
	(1 - dz) * (dHe) * (dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU + 1, inH)] +
	(dz) * (dHe) * (dU) * (1 - dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU + 1, inH)] +
	(1 - dz) * (1 - dHe) * (1 - dU) * (dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU, inH + 1)] +
	(dz) * (1 - dHe) * (1 - dU) * (dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU, inH + 1)] +
	(1 - dz) * (dHe) * (1 - dU) * (dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU, inH + 1)] +
	(dz) * (dHe) * (1 - dU) * (dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU, inH + 1)] +
	(1 - dz) * (1 - dHe) * (dU) * (dnH) * WalCool_UtoT[UtoT_IDX(0, iHe, iU + 1, inH + 1)] +
	(dz) * (1 - dHe) * (dU) * (dnH) * WalCool_UtoT[UtoT_IDX(1, iHe, iU + 1, inH + 1)] +
	(1 - dz) * (dHe) * (dU) * (dnH) * WalCool_UtoT[UtoT_IDX(0, iHe + 1, iU + 1, inH + 1)] +
	(dz) * (dHe) * (dU) * (dnH) * WalCool_UtoT[UtoT_IDX(1, iHe + 1, iU + 1, inH + 1)];

    }
  else
    {
      //iz  = 0;
      //inH = 0;

      T = (1 - dHe) * (1 - dU) * WalCool_UtoT[UtoT_collis_IDX(iHe, iU)] +
	dHe * (1 - dU) * WalCool_UtoT[UtoT_collis_IDX(iHe + 1, iU)] +
	(1 - dHe) * dU * WalCool_UtoT[UtoT_collis_IDX(iHe, iU + 1)] +
	dHe * dU * WalCool_UtoT[UtoT_collis_IDX(iHe + 1, iU + 1)];

    }

  return T;
}

#ifdef GL_DUST_COOLING
double CoolingRateFromU(double U, double Redshift, double DZ, double DL, double DS, double *Metallicities,
			double *temp)
{
  *temp = convert_u_to_tempMET(U, Redshift, DZ, Metallicities);

  return CoolingRateMET(*temp, Redshift, DZ, DL, DS, Metallicities);
#else
double CoolingRateFromU(double U, double Redshift, double DZ, double *Metallicities, double *temp)
{
  *temp = convert_u_to_tempMET(U, Redshift, DZ, Metallicities);

  return CoolingRateMET(*temp, Redshift, DZ, Metallicities);
#endif // GL_DUST_COOLING


}



double InnerInterpolation(float *zhigh, float *zlow, double *dX, int *IDX)
{
  double Lambda;
  double LambdaH, LambdaL;

  LambdaH = (1 - dU) * (1 - dnH) * zhigh[COOLRATE_IIDX(iU, inH)] +
    (dU) * (1 - dnH) * zhigh[COOLRATE_IIDX(iU + 1, inH)] +
    (1 - dU) * (dnH) * zhigh[COOLRATE_IIDX(iU, inH + 1)] +
    (dU) * (dnH) * zhigh[COOLRATE_IIDX(iU + 1, inH + 1)];

  LambdaL = (1 - dU) * (1 - dnH) * zlow[COOLRATE_IIDX(iU, inH)] +
    (dU) * (1 - dnH) * zlow[COOLRATE_IIDX(iU + 1, inH)] +
    (1 - dU) * (dnH) * zlow[COOLRATE_IIDX(iU, inH + 1)] + (dU) * (dnH) * zlow[COOLRATE_IIDX(iU + 1, inH + 1)];

  Lambda = (1 - dz) * LambdaH + dz * LambdaL;

  return Lambda;
}

#ifdef GL_DUST_COOLING
double cooling_rate_dust(double, double, double);
double CoolingRateMET(double Temp, double Redshift, double DZ, double DL, double DS, double *Metallicities)
#else
double CoolingRateMET(double Temp, double Redshift, double DZ, double *Metallicities)
#endif				// GL_DUST_COOLING
/*
  NOTE : when calling this routine directly, take care that the index U_IDX
         has been set referring to temperature and not to internal energy.
         set_indexes do the job for internal energy, then you should call
         find_index using the temperature value and specifying the temperature
         table, like in CoolingRateFromU.

         if you do not use this routine directly, but only through DoCooling and
         GetCoolingTime, do not mind about this note.
 */
{
  double Lambda1, Lambda2, Lambda, LambdaCmpt;
  double ne_He1, ne_He2, ne, neS, fact_Z, fact_ne;
  double redshift_fact;
  int el;

  double dX[4];
  int IDX[4];

  set_indexes_T(Temp, DZ, Metallicities, &dX[0], &IDX[0]);

  /* H + He cooling first */

  if(Temp < Tmin_wal)
    {
      return 0;

      /*
         if(All.ComovingIntegrationOn)
         {
         ne = 1.0; // everything is neutral
         redshift_fact = (1 + Redshift) * (1 + Redshift);
         redshift_fact *= redshift_fact;

         LambdaCmpt = 5.65e-36 * ne * (Temp - TCMB(Redshift)) * redshift_fact / Metallicities[HydPos];  // / Metallicities[HydPos];
         return -LambdaCmpt;
         }
         else
         return 0;
       */
    }


  LambdaCmpt = 0;
  if(iz >= 0)
    /* there is a photo-ionizing background */
    {
      Lambda1 = InnerInterpolation(&WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(HighZTable, iHe)],
				   &WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(LowZTable, iHe)], &dX[0],
				   &IDX[0]);
      Lambda2 =
	InnerInterpolation(&WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(HighZTable, iHe + 1)],
			   &WalMfreeCoolTables[MFREE_COOLRATE_Hef_IDX(LowZTable, iHe + 1)], &dX[0], &IDX[0]);

      ne_He1 = InnerInterpolation(&WalCool_enH[ENH_IDX(HighZTable, iHe, 0, 0)],
				  &WalCool_enH[ENH_IDX(LowZTable, iHe, 0, 0)], &dX[0], &IDX[0]);
      ne_He2 = InnerInterpolation(&WalCool_enH[ENH_IDX(HighZTable, iHe + 1, 0, 0)],
				  &WalCool_enH[ENH_IDX(LowZTable, iHe + 1, 0, 0)], &dX[0], &IDX[0]);

      ne = dHe * ne_He1 + (1 - dHe) * ne_He2;

      neS =
	InnerInterpolation(&WalCool_enHS[ENHS_z_IDX(HighZTable)], &WalCool_enHS[ENHS_z_IDX(LowZTable)],
			   &dX[0], &IDX[0]);
    }
  else
    /* no UV background */
    {
      // interpolates ne in T separately for each He
      ne_He1 =
	dU * WalCool_enH[ENH_collis_IDX(iHe, iU)] + (1 - dU) * WalCool_enH[ENH_collis_IDX(iHe, iU + 1)];
      ne_He2 =
	dU * WalCool_enH[ENH_collis_IDX(iHe + 1, iU)] + (1 -
							 dU) * WalCool_enH[ENH_collis_IDX(iHe + 1, iU + 1)];
      // interpolates ne between the two He
      ne = dHe * ne_He1 + (1 - dHe) * ne_He2;

      neS = dU * WalCool_enHS[ENHS_collis_IDX(iU)] + (1 - dU) * WalCool_enHS[ENHS_collis_IDX(iU + 1)];

      // adding inverse Comp. Scatt. off the CMB
      if(All.ComovingIntegrationOn)
	{
	  redshift_fact = (1 + Redshift) * (1 + Redshift);
	  redshift_fact *= redshift_fact;

	  LambdaCmpt = 5.65e-36 * ne * (Temp - TCMB(Redshift)) * redshift_fact / Metallicities[HydPos];	// / Metallicities[HydPos];
	}

      // interpolates lambda separately for each He
      Lambda1 = dU * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe, iU)] +
	(1 - dU) * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe, iU + 1)];
      Lambda2 = dU * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe + 1, iU)] +
	(1 - dU) * WalMfreeCoolTables[MFREE_COOLRATE_collis_IDX(iHe + 1, iU + 1)];

    }

  // interpolates lambda between the two He and adds Compton
  Lambda = LambdaCmpt + dHe * Lambda1 + (1 - dHe) * Lambda2;

  fact_ne = ne / neS;

  /* add all the elements */

  for(el = 0; el < WalCool_n_El; el++)
    {
      fact_Z = fact_ne * Metallicities[el] / WalCool_Ab[el + 2];

      if(iz >= 0)
	Lambda += fact_Z * InnerInterpolation(&WalCoolTables[COOLRATE_el_IDX(HighZTable, el)],
					      &WalCoolTables[COOLRATE_el_IDX(LowZTable, el)], &dX[0],
					      &IDX[0]);
      else
	{
	  Lambda += fact_Z * (dU * WalCoolTables[COOLRATE_collis_IDX(iU, el)] +
			      (1 - dU) * WalCoolTables[COOLRATE_collis_IDX(iU + 1, el)]);
	}
    }
#ifdef GL_DUST_COOLING

/*  double before;
  before = Lambda;
  Lambda +=  cooling_rate_dust(Temp, DL, 0.5); //
  if (Temp > 1e7)
    printf("%g %g %g %g %g %g  \n", Temp, DL, DS, before, Lambda, Lambda-before); fflush(stdout);
  #else
*/
  Lambda += cooling_rate_dust(Temp, DL, 0.05);	//
  Lambda += cooling_rate_dust(Temp, DS, 0.005);	//

#endif // GL_DUST_COOLING




  if(DBG)
    DBG = 0;

  return -Lambda;
}

#ifdef LB_PRESSURE_DEPENDENT_ACCRETION
double GetUFromLambda(double Lambda_in, double Rho, double *Metallicities, double Redshift, double DZ,
		      double dt, double *temp, double U_min, double U_max)
/*
 * note: all quantitites must be given in physical code units !
 * note: assumes that the Metallicities array has been given by set_metallicities()
 */
{
  double ratefact, Lambda, dLambda, mu;
  int iter = 0;
  double u, u_bottom, u_top, du;

  mu = (1 + 4 * yhelium) / (1 + yhelium + 1);

  u_bottom = U_min;
  u_top = U_max;

  do				// iterate to convergence
    {
      u = 0.5 * (u_bottom + u_top);

      Lambda = CoolingRateFromU(u, Redshift, DZ, Metallicities, temp);

      dLambda = Lambda - Lambda_in;

      if(Lambda_in > Lambda)
	u_top = u;
      else
	u_bottom = u;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("dLambda = %g | Lambda_in= %g\n", dLambda, Lambda_in);
    }
  while(fabs(dLambda / Lambda_in) > 1.0e-6 && iter < MAXITER);

  return u;

}
#endif

#ifdef GL_DUST_COOLING
double DoCoolingMET(double U_in, double Rho, double *Metallicities, double DL, double DS, double Redshift,
		    double DZ, double dt, double *temp)
#else // GL_DUST_COOLING
double DoCoolingMET(double U_in, double Rho, double *Metallicities, double Redshift, double DZ, double dt,
		    double *temp)
#endif				// GL_DUST_COOLING
/*
 * note: all quantitites must be given in physical code units !
 * note: assumes that the Metallicities array has been given by set_metallicities()
 */
{

  double u, u_bottom, u_top, du;
  double ratefact, Lambda;
  int iter, flag = 0;
  double fac_entr_to_u, rho_in, u_out, t_out;
#ifdef LT_WAL_DEBUG
  int i, do_heat = 0;
  double DEBUG_ARRAY_h[201][7];
  double DEBUG_ARRAY_c[201][7];
#endif

  rho_in = Rho;

  /* convert in physical cgs units */
  Rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  U_in *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

  /* in Metallicities H abundance must be in number */
  /* abundance and He's one in number ratio over H  */
  u = U_in;
  u_bottom = u_top = u;

  ratefact = Metallicities[HydPos] * Metallicities[HydPos] / Rho;

#ifdef GL_DUST_COOLING
  Lambda = CoolingRateFromU(u, Redshift, DZ, DL, DS, Metallicities, temp);
#else // GL_DUST_COOLING
  Lambda = CoolingRateFromU(u, Redshift, DZ, Metallicities, temp);
#endif // GL_DUST_COOLING


  /* cooling/heating are too weak to calculate */
  du = ratefact * Lambda * dt;
  if(fabs(du / u) < 1.0e-6)
    {
      u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;	/* to internal units */
      return u;
    }
  if(U_in + du < 0)
    {
      flag = 1;
      fac_entr_to_u = pow(rho_in, GAMMA_MINUS1) / GAMMA_MINUS1;
      u_out = All.MinEgySpec / fac_entr_to_u;
      t_out =
	convert_u_to_tempMET(u * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, Redshift, DZ,
			     Metallicities);
    }

  /* bracketing  */
  iter = 0;
  //  if(u - U_in - ratefact * Lambda * dt < 0)
  if(Lambda > 0)		// heating 
    {
      u_top *= sqrt(1.15);
      u_bottom /= sqrt(1.15);
#ifdef GL_DUST_COOLING
      while(u_top - U_in -
	    ratefact * CoolingRateFromU(u_top, Redshift, DZ, DL, DS, Metallicities, temp) * dt < 0)
#else // GL_DUST_COOLING
      while(u_top - U_in - ratefact * CoolingRateFromU(u_top, Redshift, DZ, Metallicities, temp) * dt < 0)
#endif // GL_DUST_COOLING

	{
#ifdef LT_WAL_DEBUG
	  DEBUG_ARRAY_h[iter][0] = u_bottom;
	  DEBUG_ARRAY_h[iter][1] = u_top;
	  DEBUG_ARRAY_h[iter][2] = CoolingRateFromU(u_top, Redshift, DZ, Metallicities, temp);
	  DEBUG_ARRAY_h[iter][3] = *temp;
	  DEBUG_ARRAY_h[iter][4] = (double) inH;
	  DEBUG_ARRAY_h[iter][5] = (double) iHe;
	  DEBUG_ARRAY_h[iter][6] = (double) iU;
	  do_heat++;
#endif
	  iter++;
	  u_top *= 1.15;
	  u_bottom *= 1.15;

	  if(iter > 180)
#ifdef GL_DUST_COOLING
	    printf("\t\t[WalCool][Hb][%d][%d]  %g %g %g %g %g %g %g\n", ThisTask, iter, u_bottom,
		   u_top, U_in, Rho, dt, CoolingRateFromU(u_top, Redshift, DZ, DL, DS, Metallicities, temp),
		   ratefact * CoolingRateFromU(u_top, Redshift, DZ, DL, DS, Metallicities, temp) * dt);
	  fflush(stdout);
#else
	    printf("\t\t[WalCool][Hb][%d][%d]  %g %g %g %g %g %g %g\n", ThisTask, iter, u_bottom,
		   u_top, U_in, Rho, dt, CoolingRateFromU(u_top, Redshift, DZ, Metallicities, temp),
		   ratefact * CoolingRateFromU(u_top, Redshift, DZ, Metallicities, temp) * dt);
	  fflush(stdout);
#endif // GL_DUST_COOLING

	  if(iter > 200)
	    {
#ifdef LT_WAL_DEBUG
	      printf("\n\n EPIC FAIL in bracketing (heating). here come the details:\n"
		     "U in was : %g\n"
		     "Rho is   : %g\n"
		     "Z[HydPos] is : %g\n"
		     "ratefact is : %g\n"
		     "lambda is : %g\n"
		     "dt is : %g\n\n", U_in, Rho, Metallicities[HydPos], ratefact, Lambda, dt);

	      for(i = 0; i <= 200; i++)
		printf("\t[%d] %g %g %g %g %g %g %g\n", i,
		       DEBUG_ARRAY_h[i][0],
		       DEBUG_ARRAY_h[i][1],
		       DEBUG_ARRAY_h[i][2],
		       DEBUG_ARRAY_h[i][3], DEBUG_ARRAY_h[i][4], DEBUG_ARRAY_h[i][5], DEBUG_ARRAY_h[i][6]);
#endif
	      endrun(110011001);
	    }
	}
    }

  iter = 0;
  //  if(u - U_in - ratefact * Lambda * dt > 0)
  if(Lambda < 0)		// cooling
    {
      u_top *= sqrt(1.15);
      u_bottom /= sqrt(1.15);
#ifdef GL_DUST_COOLING
      du = ratefact * CoolingRateFromU(u_bottom, Redshift, DZ, DL, DS, Metallicities, temp) * dt;
#else // GL_DUST_COOLING
      du = ratefact * CoolingRateFromU(u_bottom, Redshift, DZ, Metallicities, temp) * dt;
#endif // GL_DUST_COOLING
      while(u_bottom - U_in - du > 0)
	{
#ifdef LT_WAL_DEBUG
	  DEBUG_ARRAY_c[iter][0] = u_bottom;
	  DEBUG_ARRAY_c[iter][1] = u_top;
#ifdef GL_DUST_COOLING
	  DEBUG_ARRAY_c[iter][2] = CoolingRateFromU(u_bottom, Redshift, DZ, DL, DS, Metallicities, temp);
#else // GL_DUST_COOLING
	  DEBUG_ARRAY_c[iter][2] = CoolingRateFromU(u_bottom, Redshift, DZ, Metallicities, temp);
#endif // GL_DUST_COOLING
	  DEBUG_ARRAY_c[iter][3] = *temp;
	  DEBUG_ARRAY_c[iter][4] = (double) inH;
	  DEBUG_ARRAY_c[iter][5] = (double) iHe;
	  DEBUG_ARRAY_c[iter][6] = (double) iU;
#endif
	  iter++;
	  if(U_in + du < 0 && flag == 0)
	    {
	      flag = iter;
	      fac_entr_to_u = pow(rho_in, GAMMA_MINUS1) / GAMMA_MINUS1;
	      u_out = All.MinEgySpec / fac_entr_to_u;
	      t_out =
		convert_u_to_tempMET(u * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, Redshift, DZ,
				     Metallicities);
	    }
	  u_top /= 1.15;
	  u_bottom /= 1.15;
	  if((iter > 180) && flag == 0)
	    {


#ifdef GL_DUST_COOLING
	      printf("\t\t[WalCool][Cb][%d][%d] %g %g %g %g %g %g %g\n",
		     ThisTask, iter,
		     u_bottom, u_top, U_in, Rho, dt,
		     CoolingRateFromU(u_bottom, Redshift, DZ, DL, DS, Metallicities, temp),
		     ratefact * CoolingRateFromU(u_bottom, Redshift, DZ, DL, DS, Metallicities, temp) * dt);
#else // GL_DUST_COOLING
	      printf("\t\t[WalCool][Cb][%d][%d] %g %g %g %g %g %g %g\n",
		     ThisTask, iter,
		     u_bottom, u_top, U_in, Rho, dt,
		     CoolingRateFromU(u_bottom, Redshift, DZ, Metallicities, temp),
		     ratefact * CoolingRateFromU(u_bottom, Redshift, DZ, Metallicities, temp) * dt);
#endif // GL_DUST_COOLING

	      fflush(stdout);
	      DBG = 1;
	    }
	  if(iter > 200)
	    {
	      if(flag != 0)
		{
		  u = u_out;
		  *temp = t_out;
		  printf("WARNING: Task %d paricle on too long timestep for cooling (iter=%d)!!!\n", ThisTask,
			 flag);
		  return u;
		}
#ifdef LT_WAL_DEBUG
	      printf("\n\n EPIC FAIL in bracketing (cooling). here come the details:\n"
		     "U in was : %g\n"
		     "Rho is   : %g\n"
		     "Z[HydPos] is : %g\n"
		     "ratefact is : %g\n"
		     "lambda is : %g\n"
		     "dt is : %g\n\n", U_in, Rho, Metallicities[HydPos], ratefact, Lambda, dt);
	      if(do_heat > 0)
		{
		  printf("Did heating first ...\n");
		  for(i = 0; i < do_heat; i++)
		    printf("\t[%d] %g %g %g %g %g %g %g\n", i,
			   DEBUG_ARRAY_h[i][0],
			   DEBUG_ARRAY_h[i][1],
			   DEBUG_ARRAY_h[i][2],
			   DEBUG_ARRAY_h[i][3],
			   DEBUG_ARRAY_h[i][4], DEBUG_ARRAY_h[i][5], DEBUG_ARRAY_h[i][6]);

		}
	      for(i = 0; i <= 200; i++)
		printf("\t[%d] %g %g %g %g %g %g %g\n", i,
		       DEBUG_ARRAY_c[i][0],
		       DEBUG_ARRAY_c[i][1],
		       DEBUG_ARRAY_c[i][2],
		       DEBUG_ARRAY_c[i][3], DEBUG_ARRAY_c[i][4], DEBUG_ARRAY_c[i][5], DEBUG_ARRAY_c[i][6]);
#endif
	      endrun(110011002);
	    }
#ifdef GL_DUST_COOLING
	  du = ratefact * CoolingRateFromU(u_bottom, Redshift, DZ, DL, DS, Metallicities, temp) * dt;
#else // GL_DUST_COOLING
	  du = ratefact * CoolingRateFromU(u_bottom, Redshift, DZ, Metallicities, temp) * dt;
#endif // GL_DUST_COOLING
	}
    }


  iter = 0;

  do				// iterate to convergence
    {
      u = 0.5 * (u_bottom + u_top);
#ifdef GL_DUST_COOLING
      Lambda = CoolingRateFromU(u, Redshift, DZ, DL, DS, Metallicities, temp);
#else // GL_DUST_COOLING
      Lambda = CoolingRateFromU(u, Redshift, DZ, Metallicities, temp);
#endif // GL_DUST_COOLING

      if(u - U_in - ratefact * Lambda * dt > 0)
	{
	  u_top = u;
	}
      else
	{
	  u_bottom = u;
	}

      du = u_top - u_bottom;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("u= %g\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    printf("failed to converge in DoCoolingMET()\n");

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;	/* to internal units */

  return u;

}









/* :: -------------------------- INTERNAL ROUTINES ::  */


#ifdef GL_DUST_COOLING
double cooling_rate_dust(double temp, double delta, double a_grain)
/*
returns cooling rate as a function of gas temperature, electron number density [cgs], contributed by grains
giving a dust to gas delta and radius a_grain in micron
*/
{
  double x, lambda, h_tilde, a3;


  x = 2.71e8 * pow(a_grain, 0.66666) / temp;
  a3 = a_grain * a_grain * a_grain;

  if(x >= 4.5)
    {
      h_tilde = 5.38e-18 * pow(temp, 1.5) * a_grain * a_grain;
    }
  else if(x >= 1.5)
    {
      h_tilde = 3.37e-13 * pow(a_grain, 2.41) * pow(temp, 0.88);
    }
  else
    {
      h_tilde = 6.48e-6 * a3;
    }

  lambda = 3.3 * 7.8531e-14 * h_tilde * delta / a3;
  // the factor 3.3 originally wrong in Vogelsberger+ is included here
  // Quoting from my comment to them
  // 3.3 is "n_e*n/n_H^2 (where n is the density including all ions and
  // electrons). For a fully ionized mix of H and He with X=0.75 and
  // Y=0.25, this amounts to 30/9 (or about 3.33)"
  // 7.8531e-14 (mean molecular mass) x (proton mass) / (grain mass)
  // for a 1 mircon grain, assuming a material density of 3 gr/cm^3 and for mu=0.59
  //cgs: IDL> 0.59*1.6726219d-24 /(4.*!pi/3*3*1e-4^3)
  //si:  IDL> 0.59*1.6726219d-27 /(4.*!pi/3*3e3*1e-6^3)

  return lambda;
}
#endif // GL_DUST_COOLING

int WalCool_get_redshift_table(void)
/* this routine gets from the tables' directory what are the redshifts at
   which the tables have been defined */
{

  int n;
  DIR *WalCool_CoolTables_dir;
  struct dirent *dentry;


  /* as first, task 0 gets the content of the directory */
  /* we DO NOT use scandir because it uses the system's */
  /* malloc routine and we want to stick onto mymalloc  */
  /* (isnt'it ?)                                        */

  if((WalCool_CoolTables_dir = opendir(All.WalCool_CoolTables_path)) == NULL)
    {
      printf("The directory %s that should contain the cooling tables from Wiersma et al. does not exist\n",
	     All.WalCool_CoolTables_path);
      exit(9191000);
    }

  n = 0;			/* this first cycle is needed to determine how many   */
  /* tables are in the directory                        */
  while((dentry = readdir(WalCool_CoolTables_dir)))
    if(is_it_a_tablefile(dentry->d_name, -1))
      n++;

  WalCool_CoolTables_num = n;

  closedir(WalCool_CoolTables_dir);
  /* re-open the directory */
  WalCool_CoolTables_dir = opendir(All.WalCool_CoolTables_path);

  WalCool_CoolTables_redshifts = (struct my_direntry *) mymalloc("WalCool_CoolTables_redshifts", sizeof(struct my_direntry) * n);	/* allocate memory to store the entries  */

  n = 0;
  while((dentry = readdir(WalCool_CoolTables_dir)))
    if(is_it_a_tablefile(dentry->d_name, n))
      n++;

  closedir(WalCool_CoolTables_dir);

  qsort(WalCool_CoolTables_redshifts, WalCool_CoolTables_num, sizeof(struct my_direntry),
	compare_dir_redshifts);

  max_cool_redshift = WalCool_CoolTables_redshifts[WalCool_CoolTables_num - 1].redshift;
  return WalCool_CoolTables_num;
}



int is_it_a_tablefile(char *name, int n)
/* a file is a table file whether its name is made by the following pattern */
/*      z_[:digit:].[:digit:].hdf5                                          */
/* where [:digit:] means un unknown number od digits                        */
{
  int pos, dotpos, ndots;

  pos = strlen(name);

  if(pos < 10)			/* check for the minimum expected lenght */
    return 0;

  if(strcmp(&name[pos - 5], ".hdf5") != 0x0)	/* check whether the suffix is ".hdf5" */
    return 0;

  if(strncmp(name, "z_", 2) != 0x0)	/* check whether the prefix is "z_" */
    return 0;

  dotpos = pos - 5;		/* points to the '.' */

  for(ndots = 0, pos = 2; pos < dotpos; pos++)
    {
      if(!isdigit(name[pos]))
	{
	  if(name[pos] == '.')
	    {
	      if(ndots == 1)
		return 0;
	      else
		ndots = 1;
	    }
	  else
	    return 0;
	}
    }

  if(n >= 0)
    {
      strncpy(WalCool_CoolTables_redshifts[n].redshift_str, &name[2], dotpos - 2);
      WalCool_CoolTables_redshifts[n].redshift_str[dotpos - 2] = '\0';
      WalCool_CoolTables_redshifts[n].redshift = atof(WalCool_CoolTables_redshifts[n].redshift_str);
    }

  return 1;
}


void WalCool_get_table(int redshift_index, int tab_index)
/* load the cooling table corresponding to a given entry in the redshift table;
   the table is loaded in the tab_index position */
{
  char filename[300], dname[300];
  hid_t fid, dset;
  herr_t hd_report;

  int element;

  if(ThisTask == 0)
    {
      sprintf(filename, "%s/z_%s.hdf5", All.WalCool_CoolTables_path,
	      WalCool_CoolTables_redshifts[redshift_index].redshift_str);
      fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* get cooling from metal-free: H + He */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Metal_free/Net_Cooling");
#else
      dset = H5Dopen(fid, "/Metal_free/Net_Cooling", H5P_DEFAULT);
#endif
      hd_report =
	H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		&WalMfreeCoolTables[MFREE_COOLRATE_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);



      /* get cooling for metals */
      for(element = 0; element < WalCool_n_El; element++)
	{
	  sprintf(dname, "/%s/Net_Cooling", WalCool_El_names[element]);
#ifdef OLD_HDF5
	  dset = H5Dopen(fid, dname);
#else
	  dset = H5Dopen(fid, dname, H5P_DEFAULT);
#endif
	  hd_report =
	    H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &WalCoolTables[COOLRATE_el_IDX(tab_index, element)]);
	  hd_report = H5Dclose(dset);
	}

      /* get ne / nH table */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h");
#else
      dset = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report =
	H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_enH[ENH_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);


      /* get ne / nH solar table */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Solar/Electron_density_over_n_h");
#else
      dset = H5Dopen(fid, "/Solar/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report =
	H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_enHS[ENHS_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);

      /* get mu table */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Metal_free/Mean_particle_mass");
#else
      dset = H5Dopen(fid, "/Metal_free/Mean_particle_mass", H5P_DEFAULT);
#endif
      hd_report =
	H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_mu[MU_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);


      /* get U to T table */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Metal_free/Temperature/Temperature");
#else
      dset = H5Dopen(fid, "/Metal_free/Temperature/Temperature", H5P_DEFAULT);
#endif
      hd_report =
	H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_UtoT[UtoT_z_IDX(tab_index)]);
      hd_report = H5Dclose(dset);

      hd_report = H5Fclose(fid);
    }

  MPI_Barrier(MYMPI_COMM_WORLD);

  return;
}


void WalCool_get_collis_table()
/* loads the purely collisional cooling table, to be used at high redshift */
{
  char filename[300], dname[300];
  hid_t fid, dset;
  herr_t hd_report;

  float *rtable;
  int element, Tidx;

  if(ThisTask == 0)
    {
      sprintf(filename, "%s/z_collis.hdf5", All.WalCool_CoolTables_path);
      fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* - elecron density over n_h  NO BCKGRND */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h");
#else
      dset = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_enH);
      hd_report = H5Dclose(dset);


      /* - solar elecron density over n_h  NO BCKGRND */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Solar/Electron_density_over_n_h");
#else
      dset = H5Dopen(fid, "/Solar/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_enHS);
      hd_report = H5Dclose(dset);


      /* - mean particle mass  NO BCKGRND */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Metal_free/Mean_particle_mass");
#else
      dset = H5Dopen(fid, "/Metal_free/Mean_particle_mass", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_mu);
      hd_report = H5Dclose(dset);


      /* - temperature conversion  NO BCKGRND */
#ifdef OLD_HDF5
      dset = H5Dopen(fid, "/Metal_free/Temperature/Temperature");
#else
      dset = H5Dopen(fid, "/Metal_free/Temperature/Temperature", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_UtoT);
      hd_report = H5Dclose(dset);


      /* get cooling from metal-free: H + He */
      sprintf(dname, "/Metal_free/Net_Cooling");
#ifdef OLD_HDF5
      dset = H5Dopen(fid, dname);
#else
      dset = H5Dopen(fid, dname, H5P_DEFAULT);
#endif
      hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalMfreeCoolTables);
      hd_report = H5Dclose(dset);


      rtable = (float *) mymalloc("temporary_table", sizeof(float) * WalCool_n_T);
      /* get cooling from metals */
      for(element = 0; element < WalCool_n_El; element++)
	{
	  sprintf(dname, "/%s/Net_Cooling", WalCool_El_names[element]);
#ifdef OLD_HDF5
	  dset = H5Dopen(fid, dname);
#else
	  dset = H5Dopen(fid, dname, H5P_DEFAULT);
#endif
	  hd_report = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rtable);
	  hd_report = H5Dclose(dset);

	  for(Tidx = 0; Tidx < WalCool_n_T; Tidx++)
	    WalCoolTables[COOLRATE_collis_IDX(Tidx, element)] = rtable[Tidx];
	}

      myfree(rtable);
    }

  MPI_Barrier(MYMPI_COMM_WORLD);

  return;

}

void WalCool_Initialize_get_header(void)
/* loads everything is needed to use the tables: Rho, t, He tabs etc. */
{
#define NAME_LEN 20
  char filename[300];
  hid_t fid, dataset;
  herr_t hd_report;

  int dimensions[6], k = 0, i, pos, size;

  if(ThisTask == 0)
    {
      sprintf(filename, "%s/z_%s.hdf5", All.WalCool_CoolTables_path,
	      WalCool_CoolTables_redshifts[WalCool_CoolTables_num - 1].redshift_str);
      fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

      /* get the arrays dimensions */

#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Number_of_density_bins");
#else
      dataset = H5Dopen(fid, "/Header/Number_of_density_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_Rho);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_Rho;

#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Number_of_temperature_bins");
#else
      dataset = H5Dopen(fid, "/Header/Number_of_temperature_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_T);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_T;

      /* NOTE :: /Header/Number_of_helium_fractions SEEMS to BE ZERO!
         We Use an alternative way, however.

         #ifdef OLD_HDF5
         dataset = H5Dopen(fid, "/Header/Number_of_helium_fractions");
         #else
         dataset = H5Dopen(fid, "/Header/Number_of_helium_fractions", H5P_DEFAULT);
         #endif
         hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_Hef);
         hd_report = H5Dclose(dataset);
         dimensions[k++] = WalCool_n_Hef;
       */

      if(UseHeNumberRatio)
#ifdef OLD_HDF5
	dataset = H5Dopen(fid, "/Metal_free/Helium_number_ratio_bins");
#else
	dataset = H5Dopen(fid, "/Metal_free/Helium_number_ratio_bins", H5P_DEFAULT);
#endif
      else
#ifdef OLD_HDF5
	dataset = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins");
#else
	dataset = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins", H5P_DEFAULT);
#endif
      WalCool_n_Hef = H5Dget_storage_size(dataset) / sizeof(float);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_Hef;

#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Number_of_species");
#else
      dataset = H5Dopen(fid, "/Header/Number_of_species", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_El);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_El;

#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Abundances/Number_of_abundances");
#else
      dataset = H5Dopen(fid, "/Header/Abundances/Number_of_abundances", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &WalCool_n_Ab);
      hd_report = H5Dclose(dataset);
      dimensions[k++] = WalCool_n_Ab;
    }

  MPI_Bcast(dimensions, 6, MPI_INT, 0, MYMPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      k = 0;
      WalCool_n_Rho = dimensions[k++];
      WalCool_n_T = dimensions[k++];
      WalCool_n_Hef = dimensions[k++];
      WalCool_n_El = dimensions[k++];
      WalCool_n_Ab = dimensions[k++];
    }

  enH_mu_Size = WalCool_n_Rho * WalCool_n_T * WalCool_n_Hef;
  UtoT_Size = WalCool_n_Rho * WalCool_n_T * WalCool_n_Hef;

  /* allocate memory for array */

  WalCool_arrays_size = WalCool_n_Hef * WalCool_n_T * WalCool_n_Rho;
  WalCool_arraysS_size = WalCool_n_T * WalCool_n_Rho;

  size = (WalCool_n_Rho +	/* density bins */
	  2 * WalCool_n_T +	/* T and U bins */
	  WalCool_n_Hef +	/* Helium fractions */
	  WalCool_n_Ab +	/* solar abundances */
	  3 * 2 * WalCool_arrays_size +	/* mean molecular mass, electrons over n_H and temp. conversion @ 2 redshifts */
	  2 * WalCool_arraysS_size);	/* ne over nH for solar composition */


  WalCool_indexes_and_arrays = (float *) mymalloc("WalCool_indexe_and_arrays", sizeof(float) * size);

  WalCool_Rho = WalCool_indexes_and_arrays;
  WalCool_T = WalCool_Rho + WalCool_n_Rho;
  WalCool_U = WalCool_T + WalCool_n_T;
  WalCool_Hef = WalCool_U + WalCool_n_T;
  WalCool_Ab = WalCool_Hef + WalCool_n_Hef;

  WalCool_enHS = WalCool_Ab + WalCool_n_Ab;
  WalCool_enH = WalCool_enHS + 2 * WalCool_arraysS_size;
  WalCool_mu = WalCool_enH + 2 * WalCool_arrays_size;
  WalCool_UtoT = WalCool_mu + 2 * WalCool_arrays_size;

  /* get the arrays */
  if(ThisTask == 0)
    {
      /* - solar abundances by mass */
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Header/Abundances/Solar_mass_fractions");
#else
      dataset = H5Dopen(fid, "/Header/Abundances/Solar_mass_fractions", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_Ab);
      hd_report = H5Dclose(dataset);


      /* - elecron density over n_h */
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h");
#else
      dataset = H5Dopen(fid, "/Metal_free/Electron_density_over_n_h", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_enH);
      hd_report = H5Dclose(dataset);


      /* - helium fraction bins */
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins");
#else
      dataset = H5Dopen(fid, "/Metal_free/Helium_mass_fraction_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_Hef);
      hd_report = H5Dclose(dataset);


      /* - density bins */
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Metal_free/Hydrogen_density_bins");
#else
      dataset = H5Dopen(fid, "/Metal_free/Hydrogen_density_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_Rho);
      hd_report = H5Dclose(dataset);


      /* - mean particle mass */
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Metal_free/Mean_particle_mass");
#else
      dataset = H5Dopen(fid, "/Metal_free/Mean_particle_mass", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_mu);
      hd_report = H5Dclose(dataset);


      /* - energy density bins */
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Metal_free/Temperature/Energy_density_bins");
#else
      dataset = H5Dopen(fid, "/Metal_free/Temperature/Energy_density_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_U);
      hd_report = H5Dclose(dataset);


      /* - temperature bins */
#ifdef OLD_HDF5
      dataset = H5Dopen(fid, "/Metal_free/Temperature_bins");
#else
      dataset = H5Dopen(fid, "/Metal_free/Temperature_bins", H5P_DEFAULT);
#endif
      hd_report = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, WalCool_T);
      hd_report = H5Dclose(dataset);

      hd_report = H5Fclose(fid);
    }

  MPI_Bcast(&WalCool_indexes_and_arrays[0], size * sizeof(float), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  Trange = log10(WalCool_T[WalCool_n_T - 1] / WalCool_T[0]);
  Tmin_wal = log10(WalCool_T[0]);
  Urange = log10(WalCool_U[WalCool_n_T - 1] / WalCool_U[0]);
  Umin = log10(WalCool_U[0]);
  Rrange = log10(WalCool_Rho[WalCool_n_Rho - 1] / WalCool_Rho[0]);
  Rmin = log10(WalCool_Rho[0]);
  Herange = WalCool_Hef[WalCool_n_Hef - 1] - WalCool_Hef[0];
  Hemin = WalCool_Hef[0];

  /* allocate memory for names */

  WalCool_names_tab =
    (char **) mymalloc("WalCool_names_tab", (WalCool_n_El + 2 * WalCool_n_Ab) * sizeof(char *));
  WalCool_names =
    (char *) mymalloc("WalCool_names", (WalCool_n_El + 2 * WalCool_n_Ab) * NAME_LEN * sizeof(char));

  WalCool_El_names = WalCool_names_tab;
  for(pos = i = 0; i < WalCool_n_El; i++, pos += NAME_LEN)
    WalCool_El_names[i] = &WalCool_names[pos];

  k = pos;
  WalCool_Ab_names = WalCool_El_names + WalCool_n_El;
  for(i = 0; i < WalCool_n_Ab; i++, pos += NAME_LEN)
    WalCool_Ab_names[i] = &WalCool_names[pos];

  k = pos;
  WalCool_Ab_symbols = WalCool_Ab_names + WalCool_n_Ab;
  for(i = 0; i < WalCool_n_Ab; i++, pos += NAME_LEN)
    WalCool_Ab_symbols[i] = &WalCool_names[pos];

  /* get the names */
  if(ThisTask == 0)
    {
      char line[200];
      FILE *myfile;

      sprintf(filename, "%s/Species_names.txt", All.WalCool_CoolTables_path);
      if((myfile = fopen(filename, "r")) == NULL)
	{
	  printf("Can't open file %s\n", filename);
	  endrun(225533);
	}
      for(i = 0; i < WalCool_n_El; i++)
	{
	  fgets(WalCool_El_names[i], NAME_LEN, myfile);
	  if(WalCool_El_names[i][strlen(WalCool_El_names[i]) - 1] == '\n')
	    WalCool_El_names[i][strlen(WalCool_El_names[i]) - 1] = '\0';
	}

      fclose(myfile);
      sprintf(filename, "%s/Abundances_names.txt", All.WalCool_CoolTables_path);
      if((myfile = fopen(filename, "r")) == NULL)
	{
	  printf("Can't open file %s\n", filename);
	  endrun(225544);
	}
      for(i = 0; i < WalCool_n_Ab; i++)
	{
	  fgets(line, 200, myfile);
	  if(sscanf(line, "%s %s\n", WalCool_Ab_names[i], WalCool_Ab_symbols[i]) < 2)
	    {
	      printf("I've got some problem in reading %s/Abundances_names.txt @ line %d : %s\n",
		     All.WalCool_CoolTables_path, i, line);
	      endrun(LT_ERR_WALCOOL_IMPOSSIBLE_TO_READ_ABUNDANCE);
	    }
	  if(WalCool_Ab_names[i][strlen(WalCool_Ab_names[i]) - 1] == '\n')
	    WalCool_Ab_names[i][strlen(WalCool_Ab_names[i]) - 1] = '\0';
	  if(WalCool_Ab_symbols[i][strlen(WalCool_Ab_symbols[i]) - 1] == '\n')
	    WalCool_Ab_names[i][strlen(WalCool_Ab_symbols[i]) - 1] = '\0';
	}
      fclose(myfile);
    }

  MPI_Bcast(WalCool_names, (WalCool_n_El + 2 * WalCool_n_Ab) * NAME_LEN, MPI_BYTE, 0, MYMPI_COMM_WORLD);

  MPI_Barrier(MYMPI_COMM_WORLD);

#ifdef LT_STELLAREVOLUTION

  SpeciesIdx = (int *) mymalloc("SpeciesIdx", LT_NMet * sizeof(int));	/* defines which position a coolant has in the normal chemical array */
  Specie_is_coolant = (int *) mymalloc("Specie_is_coolant", LT_NMet * sizeof(int));	/* defined whether a specie is a coolant */
  SpeciesPos = (int *) mymalloc("SpeciesPos", WalCool_n_Ab * sizeof(int));	/* defines the ordinal position in the tables for a coolant */


  for(i = 0; i < LT_NMet; i++)
    SpeciesIdx[i] = Specie_is_coolant[i] = -1;

  for(pos = 0, i = 0; i < WalCool_n_Ab; i++)
    {
      if(strcmp(WalCool_Ab_names[i], "Hydrogen") == 0)
	{
	  myHyd = i;
	  HydPos = WalCool_n_El;
	  SpeciesPos[i] = HydPos;
	  SpeciesIdx[i] = Hyd;
	}
      else if(strcmp(WalCool_Ab_names[i], "Helium") == 0)
	{
	  myHel = i;
	  HelPos = WalCool_n_El + 1;
	  SpeciesPos[i] = HelPos;
	  SpeciesIdx[i] = Hel;
	  Specie_is_coolant[Hel] = SpeciesPos[i];
	}
      else
	{
	  for(k = 0; k < LT_NMet; k++)
	    if(strcmp(MetNames[k], WalCool_Ab_symbols[i]) == 0)
	      break;
	  if(k == LT_NMet)
	    {
	      if(ThisTask == 0)
		printf("error initializing cooling: %s (%s) is absent in chemical evolution\n",
		       WalCool_Ab_names[i], WalCool_Ab_symbols[i]);
	      fflush(stdout);
	      endrun(LT_ERR_WALCOOL_EL_MISSED);
	    }
	  SpeciesIdx[i] = k;
	  SpeciesPos[i] = pos++;
	  Specie_is_coolant[k] = SpeciesPos[i];
	}

    }

  MPI_Barrier(MYMPI_COMM_WORLD);
#endif

  return;
}


int compare_dir_redshifts(const void *A, const void *B)
/* used to sort the tables' redshifts */
{
  if((*(struct my_direntry *) A).redshift > (*(struct my_direntry *) B).redshift)
    return 1;
  if((*(struct my_direntry *) A).redshift < (*(struct my_direntry *) B).redshift)
    return -1;
  return 0;
}


int compare_redshifts(const void *A, const void *B)
/* used to sort the tables' redshifts */
{
  if((*(float *) A) > (*(struct my_direntry *) B).redshift)
    return 1;
  if((*(float *) A) < (*(struct my_direntry *) B).redshift)
    return -1;
  return 0;
}


int INLINE_FUNC find_index(int mode, float *base, double min, double range, int N, float value, double *d)
/* mode = LOG means logarithmic spacing
   mode = LIN means linear spacing
   where LOG and LIN are defined in the .data segment
*/
{
  int idx;


  if(value > base[N - 1])
    {
      *d = 0;
      return N - 2;
    }
  else if(value <= base[0])
    {
      *d = 1;
      return 0;
    }
  /* guess */

  if(mode == LOG)
    idx = (int) (log10(value) - min) / range * N;
  else
    idx = (int) (value - min) / range * N;

  if(idx < 0)
    idx = 0;
  if(idx > N - 1)
    idx = N - 1;
  /* check the guess */
  for(; value > base[idx] && idx < N - 1; idx++);

  for(; value < base[idx] && idx > 0; idx--);

  if(mode == LOG)		/* note: no pathological behaviour should be possible here */
    *d = log10(value / base[idx]) / log10(base[idx + 1] / base[idx]);	/* Hence, we avoid checks in order to speed-up             */
  else
    *d = (value - base[idx]) / (base[idx + 1] - base[idx]);

  return idx;
}


// ----------------------------------------------------------------------------------------------------------------------------

/*
  The routines here below are taken from the original cooling.c on the aim of
  providing the calculation for ionization state of H in presence of a UV
  background.
  this way, the HI field in the snapshot could be meaningful.

 */

// ----------------------------------------------------------------------------------------------------------------------------


double AbundanceRatios(double temp, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer)
{
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units
									   caller should care to convert
									   to physical from comoving */

  find_abundances_and_rates(log10(temp), rho, ne_guess);
  if(nH0_pointer != 0x0)
    *nH0_pointer = Cool.nH0;
  if(nHeII_pointer != 0x0)
    *nHeII_pointer = Cool.nHepp;

  return Cool.nH0;
}


/* this function computes the actual abundance ratios 
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess)
/*
  note: arguments are expected to be in physical units; use AbundanceRatios instead.
 */
{
  double neold, nenew;
  int j, niter;
  double Tlow, flow, fhi, t;

  double logT_input, rho_input, ne_input;

  logT_input = logT;
  rho_input = rho;
  ne_input = *ne_guess;

  if(logT <= Cool.Tmin)		/* everything neutral */
    {
      cl_set_nH0(1);
//      nH0 = 1.0;
      cl_set_nHe0(yhelium);
//      nHe0 = yhelium;
      cl_set_nHp(0);
//      nHp = 0;
      cl_set_nHep(0);
//      nHep = 0;
      cl_set_nHepp(0);
//      nHepp = 0;
      cl_set_ne(0);
//      ne = 0;
      *ne_guess = 0;
      return;
    }

  if(logT >= Cool.Tmax)		/* everything is ionized */
    {
      cl_set_nH0(0);
//      nH0 = 0;
      cl_set_nHe0(0);
//      nHe0 = 0;
      cl_set_nHp(1.0);
//      nHp = 1.0;
      cl_set_nHep(0);
//      nHep = 0;
      cl_set_nHepp(yhelium);
//      nHepp = yhelium;
      cl_set_ne(Cool.nHp + 2.0 * Cool.nHepp);
//      ne = nHp + 2.0 * nHepp;
      *ne_guess = Cool.ne;	/* note: in units of the hydrogen number density */
      return;
    }

  t = (logT - Cool.Tmin) / Cool.deltaT;
  j = (int) t;
  Tlow = Cool.Tmin + Cool.deltaT * j;
  //  Thi = Tlow + deltaT;
  fhi = t - j;
  flow = 1 - fhi;

  if(*ne_guess == 0)
    *ne_guess = 1.0;

  cl_set_nHcgs(XH * rho / PROTONMASS);
//  nHcgs = XH * rho / PROTONMASS;      /* hydrogen number dens in cgs units */

  cl_set_ne(*ne_guess);
//  ne = *ne_guess;
  neold = Cool.ne;
  niter = 0;
  cl_set_necgs(Cool.ne * Cool.nHcgs);

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do
    {
      niter++;

      cl_set_aHp(flow * Cool.AlphaHp[j] + fhi * Cool.AlphaHp[j + 1]);
      cl_set_aHep(flow * Cool.AlphaHep[j] + fhi * Cool.AlphaHep[j + 1]);
      cl_set_aHepp(flow * Cool.AlphaHepp[j] + fhi * Cool.AlphaHepp[j + 1]);
      cl_set_ad(flow * Cool.Alphad[j] + fhi * Cool.Alphad[j + 1]);
      cl_set_geH0(flow * Cool.GammaeH0[j] + fhi * Cool.GammaeH0[j + 1]);
      cl_set_geHe0(flow * Cool.GammaeHe0[j] + fhi * Cool.GammaeHe0[j + 1]);
      cl_set_geHep(flow * Cool.GammaeHep[j] + fhi * Cool.GammaeHep[j + 1]);
//      aHp = flow * Cool.AlphaHp[j] + fhi * Cool.AlphaHp[j + 1];
//      aHep = flow * Cool.AlphaHep[j] + fhi * Cool.AlphaHep[j + 1];
//      aHepp = flow * Cool.AlphaHepp[j] + fhi * Cool.AlphaHepp[j + 1];
//      ad = flow * Cool.Alphad[j] + fhi * Cool.Alphad[j + 1];
//      geH0 = flow * Cool.GammaeH0[j] + fhi * Cool.GammaeH0[j + 1];
//      geHe0 = flow * Cool.GammaeHe0[j] + fhi * Cool.GammaeHe0[j + 1];
//      geHep = flow * Cool.GammaeHep[j] + fhi * Cool.GammaeHep[j + 1];

      if(Cool.necgs <= 1.e-25 || Cool.J_UV == 0)
	{
	  cl_set_gJH0ne(0);
	  cl_set_gJHe0ne(0);
	  cl_set_gJHepne(0);
//        gJH0ne = gJHe0ne = gJHepne = 0;
	}
      else
	{
	  cl_set_gJH0ne(Cool.gJH0 / Cool.necgs);
	  cl_set_gJHe0ne(Cool.gJHe0 / Cool.necgs);
	  cl_set_gJHepne(Cool.gJHep / Cool.necgs);
//        gJH0ne = Cool.gJH0 / Cool.necgs;
//        gJHe0ne = Cool.gJHe0 / Cool.necgs;
//        gJHepne = Cool.gJHep / Cool.necgs;
	}

      cl_set_nH0(Cool.aHp / (Cool.aHp + Cool.geH0 + Cool.gJH0ne));
//      nH0 = Cool.aHp / (Cool.aHp + Cool.geH0 + Cool.gJH0ne);  /* eqn (33) */
      cl_set_nHp(1.0 - Cool.nH0);
//      nHp = 1.0 - Cool.nH0;           /* eqn (34) */

      if((Cool.gJHe0ne + Cool.geHe0) <= SMALLNUM)	/* no ionization at all */
	{
	  cl_set_nHep(0.0);
//        nHep = 0.0;
	  cl_set_nHepp(0.0);
//        nHepp = 0.0;
	  cl_set_nHe0(yhelium);
//        nHe0 = yhelium;
	}
      else
	{
	  cl_set_nHep(yhelium /
		      (1.0 + (Cool.aHep + Cool.ad) / (Cool.geHe0 + Cool.gJHe0ne) +
		       (Cool.geHep + Cool.gJHepne) / Cool.aHepp));
//        nHep = yhelium / (1.0 + (Cool.aHep + Cool.ad) / (Cool.geHe0 + Cool.gJHe0ne) + (Cool.geHep + Cool.gJHepne) / Cool.aHepp);      /* eqn (35) */
	  cl_set_nHe0(Cool.nHep * (Cool.aHep + Cool.ad) / (Cool.geHe0 + Cool.gJHe0ne));
//        nHe0 = nHep * (Cool.aHep + Cool.ad) / (Cool.geHe0 + Cool.gJHe0ne);    /* eqn (36) */
	  cl_set_nHepp(Cool.nHep * (Cool.geHep + Cool.gJHepne) / Cool.aHepp);
//        nHepp = Cool.nHep * (Cool.geHep + Cool.gJHepne) / Cool.aHepp; /* eqn (37) */
	}

      neold = Cool.ne;

      cl_set_ne(Cool.nHp + Cool.nHep + 2 * Cool.nHepp);
//      ne = nHp + nHep + 2 * nHepp;    /* eqn (38) */
      cl_set_necgs(Cool.ne * Cool.nHcgs);

      if(Cool.J_UV == 0)
	break;

      nenew = 0.5 * (Cool.ne + neold);
      cl_set_ne(nenew);
//      ne = nenew;
      cl_set_necgs(Cool.ne * Cool.nHcgs);

      if(fabs(Cool.ne - neold) < 1.0e-4)
	break;

      if(niter > (MAXITER - 10))
	printf("ne= %g  niter=%d init=(%g<%g<%g?,rho=%g,ne=%g)\n", Cool.ne, niter, Cool.Tmin, logT_input,
	       Cool.Tmax, rho_input, ne_input);
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("no convergence reached in find_abundances_and_rates()\n");
      printf("logT_input= %g  rho_input= %g  ne_input= %g\n", logT_input, rho_input, ne_input);
      //printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
      //     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      endrun(13);
    }

  cl_set_bH0(flow * Cool.BetaH0[j] + fhi * Cool.BetaH0[j + 1]);
  cl_set_bHep(flow * Cool.BetaHep[j] + fhi * Cool.BetaHep[j + 1]);
  cl_set_bff(flow * Cool.Betaff[j] + fhi * Cool.Betaff[j + 1]);
//  bH0 = flow * Cool.BetaH0[j] + fhi * Cool.BetaH0[j + 1];
//  bHep = flow * Cool.BetaHep[j] + fhi * Cool.BetaHep[j + 1];
//  bff = flow * Cool.Betaff[j] + fhi * Cool.Betaff[j + 1];

  *ne_guess = Cool.ne;
}


/* double binterp(float x1, float x2, float y1, float y2, float vx1y1, float vx2y1, float vx1y2, float vx2y2, float x, float y) */
/* /\* calculate the bilinear interpolation on a square like the one here below *\/ */
/* /\* */
/*         ^ */
/*         |   (Vx1y2)          (Vx2y2) */
/*     y2  |     +---------------+ */
/*         |     |               | */
/*         |     |               | */
/*         |     |               | */
/*         |     |        (x,y)  | */
/*         |     |..........+    | */
/*         |     |          :    | */
/*         |     |          :    | */
/*     y1  |     +---------------+ */
/*         |   (Vx1y1)          (Vx2y1) */
/*         | */
/*         +---------------------------------->     */
/*              x1               x2 */

/*  *\/ */
/* { */
/*   double value; */
/*   double nx, ny; */

/*   if(x1 == x2) */
/*     nx = 0; */
/*   else */
/*     nx = (x - x1) / (x2 - x1); */

/*   if(y1 == y2) */
/*     ny = 0; */
/*   else */
/*     ny = (y - y1) / (y2 - y1); */

/*   value = (vx1y1 * (1 - nx)*(1 - ny) + */
/*            vx2y1 * nx * (1 - ny) + */
/*            vx1y2 * ny * (1 - nx) + */
/*            vx2y2 * nx * ny); */

/*   return value; */
/* } */


/* int INLINE_FUNC find_1D_index_in_ND_matrix(int N, int *L, int *gridpoints) */
/* { */
/*   int index, i, partial; */

/*   /\* find the index in the uni-dimensional array that represents the */
/*      N-dimensional matrix */

/*      L is an array that stores the length of each dimension of the */
/*      matrix; */
/*      gridpoint is an array that specifies, in N-dimensional coordinates. */
/*      the grid point whose 1-dimensional index we want to find out. */
/*    *\/ */
/*   for(index = 0, i = N-1; i >= 0; i--) */
/*     { */
/*       for(partial = gridpoint[i], k = i-1; k >= 0; k--) */
/*         partial *= L[k]; */
/*       index += partial; */
/*     } */

/*   return index; */
/* } */



/* float get_Ndim_interp(int N, double *P, double *Dx) */
/* /\* */
/*   This routine gives the linear interpolation in a N-dimensional */
/*   function defined on a grid. */

/*   --==[ INPUT VALUES ]==-- */

/*   N  is the number of dimensions */
/*   P  is a 2^N long array that stores the values of the function */
/*      at the grid points that encompasse the interpolation point */
/*   Dx is an N-dimensional array that stores the normalized */
/*      displacement of the interpolation point with respect to the */
/*      lower grid point: */
/*         Dx_i =  (x_i - x0_i) / (x1_i - x0_i) */

/*         in two dimensions, for instance: */

/*         ^ */
/*         |                             */
/*    x1_1 |     +---------------+ */
/*         |     |               | */
/*         |     |               | */
/*         |     |               | */
/*         |     |      (x_0,x_1)| */
/*         |     |..........+    | */
/*         |     |          :    | */
/*         |     |          :    | */
/*    x0_1 |     +---------------+ */
/*         |                            */
/*         | */
/*         +---------------------------------->     */
/*              x0_0             x1_0 */


/*   --==[ HOW TO FILL THE ARRAYS ]==-- */

/*   As first you must define what is the order of the variables x,i */
/*   (that is irrelevant: interpolation is commutative ini this respect). */

/*   Then you fill the Dx array with that order. */

/*   The P array is filled in the same "walking order" of the N-dimensional */
/*   cube - which has 2^N verteces - that encloses the interpolation point. */

/*  *\/ */
/* { */
/*   double INTERP = 0, partial; */
/*   int    i, j; */

/*   for(j = 0; j < (1<<N); j++) */
/*     { */
/*       for(i = 0; i < N; i++) */
/*         S[i] = j & (1 << i); */

/*       for(partial = 1, i = 0; (i < N) && (partial> 0); i++) */
/*         { */
/*           if(S[i] == 0) */
/*             partial *= (1.0 - Dx[i]); */
/*           else */
/*             partial *= Dx[i]; */

/*           partial *= P[i]; */
/*         } */

/*       INTERP += partial; */
/*     } */

/*   return (float)INTERP; */


/*   /\* */
/*     The above routine is the algorithmical translation of the following */
/*     formula: */
/*                                                _                       _ */
/*                                  ----         /     _____               \                       */
/*                                  \            |      | |      [m,j]     |                       */
/*      interpolated value =        /            |      | |     U      [X] | F(m1 + s1,..,mn+sn)   */
/*                                  ----         \_ j in {1,..n}  s_j     _/                       */
/*                             s1,..,sn in 0,1} */

/*     where: */

/*     n are the dimensions */
/*     m1,..,mn are the grid points */
/*     x1,..,xn are the coordinates of the interpolation point */
/*     s1,..,sn are either 0 or 1 and means the leap from a grid point to */
/*              the subsequent one along the same coordinate */

/*      [m,j] */
/*     U     [X] = (1 - s_j) + (x_j - m_j)(2s_j - 1) */
/*      s_j */


/*    *\/ */
/* } */

#ifdef GL_DUST_COOLING
double get_cooling_time_fromT(double temp, double u,
			      double rho, double *ne_guess, double Redshift, double DZ, double DL, double DS,
			      double *Metallicities)
#else
double get_cooling_time_fromT(double temp, double u,
			      double rho, double *ne_guess, double Redshift, double DZ, double *Metallicities)
#endif				// GL_DUST_COOLING
{
  double u_old;
  double ratefact, XH;
  double LambdaNet, coolingtime;
  double DoCool_rho_input, DoCool_ne_guess_input, nHcgs;
  FILE *f;


  DoCool_rho_input = rho;
  DoCool_ne_guess_input = *ne_guess;

  XH = HYDROGEN_MASSFRAC;

  nHcgs = XH * rho / PROTONMASS;

  ratefact = nHcgs * nHcgs / rho;

  u_old = u;
#ifdef GL_DUST_COOLING
  LambdaNet = CoolingRateMET(temp, Redshift, DZ, DL, DS, Metallicities);
#else
  LambdaNet = CoolingRateMET(temp, Redshift, DZ, Metallicities);
#endif // GL_DUST_COOLING
  if(LambdaNet >= 0)
    return 0;

  coolingtime = u_old / (-ratefact * LambdaNet);

  return coolingtime;
}

#endif
