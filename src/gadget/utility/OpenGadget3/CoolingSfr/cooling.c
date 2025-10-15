/*
 * @file
 * This file is part of the developer version of GADGET3 and contains
 * the license conditions for its usage.
 *
 * @author GADGET-development team, led by Volker Springel and Klaus Dolag.
 *
 * @section LICENSE
 * Copyright (c) 2016, Volker Springel, Klaus Dolag, and all contributing
 * authors (see change logs). All rights reserved.
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
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"
//#include "cooling.h"

/*! \file cooling.c
 *  \brief implements gas cooling routines
 *
 *  Relevant papers:
 *  Katz 95, https://arxiv.org/abs/astro-ph/9509107, Section 3
 *
 *  As part of the effort to properly document the cooling code please add relevant papers if you find them.
 */




#ifdef COOLING


// TODO does this construction serve a purpose?
#define __STR(M) #M
#define __XSTR(N) __STR(N)
#pragma message "GADGET_COOLING set to  " __XSTR(GADGET_COOLING)
#undef __STR
#undef __XSTR

#ifdef LT_METAL_COOLING_WAL
extern double DZ, Redshift;
#ifdef GL_DUST_COOLING
double DL,DS=0; // here defined only to allow compilation of cl_CoolingOnly
#endif
#endif

double XH = HYDROGEN_MASSFRAC;	/* hydrogen abundance by mass */
double yhelium;

struct global_cooling_data Cool;

/* table input (from file TREECOOL) for ionizing parameters */

float inlogz[TABLESIZE];
float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
int nheattab;			/* length of table */

/**********************************************************************************************************************/
/*                                                      InitCool                                                      */
/**********************************************************************************************************************/
// Everything in this block is strictly executed by calling InitCool()

// called once in begrun.c and lt_sn_init.c


void InitCool(void)
{
  if(ThisTask == 0)
    {
      printf("Initializing cooling ...\n");
    }

  InitCoolMemory();
  cl_set_Tmin(4.0);
  cl_set_Tmax(8.0);
  MakeCoolingTable();

  if(ThisTask == 0)
    printf("   reading TREECOOL tables ...\n");
  ReadIonizeParams("TREECOOL");

  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  IonizeParams();

#if GADGET_COOLING == COOLING_SUTHERLAND
  if(ThisTask == 0)
    printf("   reading S&D tables ...\n");
  read_cooling_tables("coolingtable");

#elif GADGET_COOLING == COOLING_WAL
  read_cooling_tables_dummy();
  if(ThisTask == 0)
    printf("   reading WAL tables ...\n");
  WalCool_Initialize();

#elif GADGET_COOLING == COOLING_LUEDERS
  if(ThisTask == 0)
    printf("   reading delauni sampled CLOUDI tables ...\n");
#endif
}



static void InitCoolMemory(void)
{
  Cool.BetaH0 = (double *) mymalloc("Cool.BetaH0", (NCOOLTAB + 1) * sizeof(double));
  Cool.BetaHep = (double *) mymalloc("Cool.BetaHep", (NCOOLTAB + 1) * sizeof(double));
  Cool.AlphaHp = (double *) mymalloc("Cool.AlphaHp", (NCOOLTAB + 1) * sizeof(double));
  Cool.AlphaHep = (double *) mymalloc("Cool.AlphaHep", (NCOOLTAB + 1) * sizeof(double));
  Cool.Alphad = (double *) mymalloc("Cool.Alphad", (NCOOLTAB + 1) * sizeof(double));
  Cool.AlphaHepp = (double *) mymalloc("Cool.AlphaHepp", (NCOOLTAB + 1) * sizeof(double));
  Cool.GammaeH0 = (double *) mymalloc("Cool.GammaeH0", (NCOOLTAB + 1) * sizeof(double));
  Cool.GammaeHe0 = (double *) mymalloc("Cool.GammaeHe0", (NCOOLTAB + 1) * sizeof(double));
  Cool.GammaeHep = (double *) mymalloc("Cool.GammaeHep", (NCOOLTAB + 1) * sizeof(double));
  Cool.Betaff = (double *) mymalloc("Cool.Betaff", (NCOOLTAB + 1) * sizeof(double));
}


/*!
 * Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105,
 * 19 Hydrogen, Helium III recombination rates and collisional ionization
 * cross-sections are updated
 */
 // TODO FUnction has side effects, but only within cooling scope (except maybe XH, yhelium, not sure)
 // TODO Function overrides HYDROGEN_ONLY (XH is not HYDROGEN_MASSFRAC)
 // No dependencies
static void MakeCoolingTable(void)
{
#ifdef NEW_RATES
  double dE, P, A, X, K, U, T_eV;
  double b0, b1, b2, b3, b4, b5, c0, c1, c2, c3, c4, c5, y;	/* used in Scholz-Walter fit */
  double E1s_2, Gamma1s_2s, Gamma1s_2p;
#endif

  XH = 0.76;
  yhelium = (1 - XH) / (4 * XH);

  if(All.MinGasTemp > 0.0)
    {
      Cool.Tmin = log10(0.1 * All.MinGasTemp);
    }
  else
    {
      Cool.Tmin = 1.0;
    }

  Cool.deltaT = (Cool.Tmax - Cool.Tmin) / NCOOLTAB;

  for(int i = 0; i <= NCOOLTAB; i++)
    {
      Cool.BetaH0[i] = Cool.BetaHep[i] = Cool.Betaff[i] = Cool.AlphaHp[i] = Cool.AlphaHep[i] =
	Cool.AlphaHepp[i] = Cool.Alphad[i] = Cool.GammaeH0[i] = Cool.GammaeHe0[i] = Cool.GammaeHep[i] = 0;

      double T = pow(10.0, Cool.Tmin + Cool.deltaT * i);
      double Tfact = 1.0 / (1 + sqrt(T / 1.0e5));

      // TODO figure out who wrote this and slap them
      if(118348 / T < 70)
	{
	  Cool.BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;
	}

#ifdef NEW_RATES
      // TODO turn into struct. magic numbers
      /* Scholtz-Walters 91 fit */
      if(T >= 2.0e3 && T < 1e8)
	{

	  if(T >= 2.0e3 && T < 6.0e4)
	    {
	      b0 = -3.299613e1;
	      b1 = 1.858848e1;
	      b2 = -6.052265;
	      b3 = 8.603783e-1;
	      b4 = -5.717760e-2;
	      b5 = 1.451330e-3;

	      c0 = -1.630155e2;
	      c1 = 8.795711e1;
	      c2 = -2.057117e1;
	      c3 = 2.359573;
	      c4 = -1.339059e-1;
	      c5 = 3.021507e-3;
	    }
	  else
	    {
	      if(T >= 6.0e4 && T < 6.0e6)
		{
		  b0 = 2.869759e2;
		  b1 = -1.077956e2;
		  b2 = 1.524107e1;
		  b3 = -1.080538;
		  b4 = 3.836975e-2;
		  b5 = -5.467273e-4;

		  c0 = 5.279996e2;
		  c1 = -1.939399e2;
		  c2 = 2.718982e1;
		  c3 = -1.883399;
		  c4 = 6.462462e-2;
		  c5 = -8.811076e-4;
		}
	      else
		{
		  b0 = -2.7604708e3;
		  b1 = 7.9339351e2;
		  b2 = -9.1198462e1;
		  b3 = 5.1993362;
		  b4 = -1.4685343e-1;
		  b5 = 1.6404093e-3;

		  c0 = -2.8133632e3;
		  c1 = 8.1509685e2;
		  c2 = -9.4418414e1;
		  c3 = 5.4280565;
		  c4 = -1.5467120e-1;
		  c5 = 1.7439112e-3;
		}

	      y = log(T);
	      E1s_2 = 10.2;	/* eV */

	      Gamma1s_2s = exp(b0 + b1 * y + b2 * y * y + b3 * y * y * y +
			       b4 * y * y * y * y + b5 * y * y * y * y * y);
	      Gamma1s_2p = exp(c0 + c1 * y + c2 * y * y + c3 * y * y * y +
			       c4 * y * y * y * y + c5 * y * y * y * y * y);

	      T_eV = T / eV_to_K;
	      Cool.BetaH0[i] = E1s_2 * eV_to_erg * (Gamma1s_2s + Gamma1s_2p) * exp(-E1s_2 / T_eV);
	    }
	}
#endif

      // TODO magic numbers...
      if(473638 / T < 70)
	{
	  Cool.BetaHep[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;
	}

      Cool.Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));

#ifdef NEW_RATES
      Cool.AlphaHp[i] = 6.28e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);
      Cool.AlphaHepp[i] = 3.36e-10 * pow(T / 1000, -0.2) / (1. + pow(T / 4.0e6, 0.7)) / sqrt(T);
#else
      Cool.AlphaHp[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);	/* old Cen92 fit */
      Cool.AlphaHepp[i] = 4. * Cool.AlphaHp[i];	/* old Cen92 fit */
#endif

      Cool.AlphaHep[i] = 1.5e-10 * pow(T, -0.6353);

      if(470000 / T < 70)
	{
	  Cool.Alphad[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));
	}

#ifdef NEW_RATES
      T_eV = T / eV_to_K;

      /* Voronov 97 fit */
      /* hydrogen */
      dE = 13.6;
      P = 0.0;
      A = 0.291e-7;
      X = 0.232;
      K = 0.39;

      U = dE / T_eV;
      Cool.GammaeH0[i] = A * (1.0 + P * sqrt(U)) * pow(U, K) * exp(-U) / (X + U);

      /* Helium */
      dE = 24.6;
      P = 0.0;
      A = 0.175e-7;
      X = 0.18;
      K = 0.35;

      U = dE / T_eV;
      Cool.GammaeHe0[i] = A * (1.0 + P * sqrt(U)) * pow(U, K) * exp(-U) / (X + U);

      /* Hellium II */
      dE = 54.4;
      P = 1.0;
      A = 0.205e-8;
      X = 0.265;
      K = 0.25;

      U = dE / T_eV;
      Cool.GammaeHep[i] = A * (1.0 + P * sqrt(U)) * pow(U, K) * exp(-U) / (X + U);

#else // NEW_RATES
      if(157809.1 / T < 70)
	{
	  Cool.GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;
	}

      if(285335.4 / T < 70)
	{
	  Cool.GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;
	}

      if(631515.0 / T < 70)
	{
	  Cool.GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
	}
#endif
    }				// end for
}


/*!
 * table input (from file TREECOOL) for ionizing parameters
 */
static void ReadIonizeParams(char *fname)
{
  FILE *fdcool;

  if(!(fdcool = fopen(fname, "r")))
    {
      printf(" Cannot read ionization table in file `%s'\n", fname);
      endrun(456);
    }

  for(int i = 0; i < TABLESIZE; i++)
    {
      gH0[i] = 0;
    }

  for(int i = 0; i < TABLESIZE; i++)
    {
      if(fscanf
	 (fdcool, "%g %g %g %g %g %g %g", &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i],
	  &eHep[i]) == EOF)
	{
	  break;
	}
    }

  fclose(fdcool);

  /*  nheattab is the number of entries in the table */
  nheattab = 0;
  for(int i = 0; i < TABLESIZE; i++)
    {
      if(gH0[i] != 0.0)
	{
	  nheattab++;
	}
      else
	{
	  break;
	}
    }

  if(ThisTask == 0)
    {
      printf("\n\nread ionization table with %d entries in file `%s'.\n\n", nheattab, fname);
    }
}


/**********************************************************************************************************************/
/*                                                  End of InitCool                                                   */
/**********************************************************************************************************************/





/**********************************************************************************************************************/
/*                                              Setter/Getter Functions                                               */
/**********************************************************************************************************************/

// As part of the refactoring of the cooling code all cooling code globals are consolidated into struct Cool
// Modifying Cool from outside this file is STRONGLY discouraged
// For backwards compatibility these setter and getter functions are provided

// Currently only used once, in lt_wal_cooling.c/read_cooling_tables_dummy()
void cl_set_Tmin(double value)
{
  Cool.Tmin = value;
}


// Currently only used once, in lt_wal_cooling.c/read_cooling_tables_dummy()
void cl_set_Tmax(double value)
{
  Cool.Tmax = value;
}


// These are used by lt_wal_cooling.c, mostly
void cl_set_ne(double value)
{
  Cool.ne = value;
}

void cl_set_necgs(double value)
{
  Cool.necgs = value;
}

void cl_set_nHcgs(double value)
{
  Cool.nHcgs = value;
}


// The follwoing functions are called exactly once, in lt_wal_cooling.c
// the variable Cool.bH0 is also only used once, btw! TODO can probably be made non global somehow
void cl_set_bH0(double value)
{
  Cool.bH0 = value;
}

void cl_set_bHep(double value)
{
  Cool.bHep = value;
}

void cl_set_bff(double value)
{
  Cool.bff = value;
}


void cl_set_aHp(double value)
{
  Cool.aHp = value;
}


void cl_set_aHep(double value)
{
  Cool.aHep = value;
}

void cl_set_aHepp(double value)
{
  Cool.aHepp = value;
}

void cl_set_ad(double value)
{
  Cool.ad = value;
}

void cl_set_geH0(double value)
{
  Cool.geH0 = value;
}

void cl_set_geHe0(double value)
{
  Cool.geHe0 = value;
}

void cl_set_geHep(double value)
{
  Cool.geHep = value;
}

void cl_set_gJH0ne(double value)
{
  Cool.gJH0ne = value;
}

void cl_set_gJHe0ne(double value)
{
  Cool.gJHe0ne = value;
}

void cl_set_gJHepne(double value)
{
  Cool.gJHepne = value;
}


void cl_set_nH0(double value)
{
  Cool.nH0 = value;
}

void cl_set_nHp(double value)
{
  Cool.nHp = value;
}

void cl_set_nHep(double value)
{
  Cool.nHep = value;
}

void cl_set_nHe0(double value)
{
  Cool.nHe0 = value;
}

void cl_set_nHepp(double value)
{
  Cool.nHepp = value;
}


/**********************************************************************************************************************/
/*                                           End of Setter/Getter Functions                                           */
/**********************************************************************************************************************/


/**********************************************************************************************************************/
/*                                                  Static Functions                                                  */
/**********************************************************************************************************************/

// These functions are only used in cooling.c

/*!
 * Calculates (heating rate-cooling rate)/n_h^2 in cgs units
 */
//#ifdef LT_METAL_COOLING
//double CoolingRateZ(double logT, double rho, double *nelec, double Z)
//#else
//double CoolingRateSTD(double logT, double rho, double *nelec)
//#endif
//{
static double CoolingRate(CoolingRateArgs *args)
{
//   unpack arguments
//   TODO using directly would be faster
  double logT = args->logT;
  double rho = args->rho;
  double *nelec = args->nelec;

#ifdef LT_METAL_COOLING
  double Z = args->Z;
#endif


  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double redshift;
  double T;

#ifdef LT_METAL_COOLING
  double LambdaMet;
#endif
#if defined(UM_METAL_COOLING) && defined(UM_MET_IN_LT_COOLING)
  double LambdaMetalLines;
#endif

  if(logT <= Cool.Tmin)
    {
      logT = Cool.Tmin + 0.5 * Cool.deltaT;	/* floor at Tmin */
    }

  Cool.nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

  if(logT < Cool.Tmax)
    {
      find_abundances_and_rates(logT, rho, nelec);

      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
#ifdef LT_METAL_COOLING
      /* get cooling rate from sutherland and dopita tables */
      /* keeping abundances and rates calculations */

#if defined(UM_METAL_COOLING) && defined(UM_MET_IN_LT_COOLING)
      if((Z >= ZMin) && (logT >= log10(UM_LT_INF_LIMIT)))
#else
      if((Z > ZMin) && (logT >= Cool.Tmin))
#endif // defined...
	{
	  if((LambdaMet =
	      (GetMetalLambda(logT, Z) - GetMetalLambda(logT, ZMin)) * Cool.ne * (Cool.nHp + Cool.nHep +
										  Cool.nHepp)) < 0)
	    LambdaMet = 0;
	}
      else
	{
	  LambdaMet = 0;
	}

#endif // LT_METAL_COOLING

      T = pow(10.0, logT);

#if defined(UM_CHEMISTRY) && defined(UM_METAL_COOLING) && defined(UM_MET_IN_LT_COOLING)
#ifdef UM_e_MET_IMPACTS
      LambdaMetalLines = um_metal_line_cooling(T, rho, LT_NMet, ne * Cool.nHcgs) / (Cool.nHcgs * Cool.nHcgs);
#else
      LambdaMetalLines = um_metal_line_cooling(T, rho, LT_NMet) / (Cool.nHcgs * Cool.nHcgs);
#endif
#endif

      LambdaExcH0 = Cool.bH0 * Cool.ne * Cool.nH0;
      LambdaExcHep = Cool.bHep * Cool.ne * Cool.nHep;
      LambdaExc = LambdaExcH0 + LambdaExcHep;	/* excitation */

      LambdaIonH0 = 2.18e-11 * Cool.geH0 * Cool.ne * Cool.nH0;
      LambdaIonHe0 = 3.94e-11 * Cool.geHe0 * Cool.ne * Cool.nHe0;
      LambdaIonHep = 8.72e-11 * Cool.geHep * Cool.ne * Cool.nHep;
      LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* ionization */

      LambdaRecHp = 1.036e-16 * T * Cool.ne * (Cool.aHp * Cool.nHp);
      LambdaRecHep = 1.036e-16 * T * Cool.ne * (Cool.aHep * Cool.nHep);
      LambdaRecHepp = 1.036e-16 * T * Cool.ne * (Cool.aHepp * Cool.nHepp);
      LambdaRecHepd = 6.526e-11 * Cool.ad * Cool.ne * Cool.nHep;
      LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;

      LambdaFF = Cool.bff * (Cool.nHp + Cool.nHep + 4 * Cool.nHepp) * Cool.ne;

      Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

#if defined(UM_CHEMISTRY) && defined(UM_METAL_COOLING) && defined(UM_MET_IN_LT_COOLING)
      Lambda += LambdaMetalLines;
#endif

      if(All.ComovingIntegrationOn)
	{
	  redshift = 1 / All.Time - 1;
	  LambdaCmptn = 5.65e-36 * Cool.ne * (T - 2.73 * (1. + redshift)) *
	    pow(1. + redshift, 4.) / Cool.nHcgs;

	  Lambda += LambdaCmptn;
	}
      else
	{
	  LambdaCmptn = 0;
	}
      Heat = 0;
      if(Cool.J_UV != 0)
	Heat += (Cool.nH0 * Cool.epsH0 + Cool.nHe0 * Cool.epsHe0 + Cool.nHep * Cool.epsHep) / Cool.nHcgs;

#ifdef LT_METAL_COOLING
      Lambda += LambdaMet;
#endif
    }
  else				/* here we're outside of tabulated rates, T>Cool.Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are
         present. Assumes no heating. */

      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep =
	LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

      /* very hot: H and He both fully ionized */
      Cool.nHp = 1.0;
      Cool.nHep = 0;
      Cool.nHepp = yhelium;
      Cool.ne = Cool.nHp + 2.0 * Cool.nHepp;
      *nelec = Cool.ne;		/* note: in units of the hydrogen number density */

      T = pow(10.0, logT);
      LambdaFF = 1.42e-27 * sqrt(T) *
	(1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (Cool.nHp + 4 * Cool.nHepp) * Cool.ne;

      if(All.ComovingIntegrationOn)
	{
	  redshift = 1 / All.Time - 1;
	  /* add inverse Compton cooling off the microwave background */
	  LambdaCmptn = 5.65e-36 * Cool.ne * (T - 2.73 * (1. + redshift)) *
	    pow(1. + redshift, 4.) / Cool.nHcgs;
	}
      else
	LambdaCmptn = 0;

      Lambda = LambdaFF + LambdaCmptn;
    }

  return (Heat - Lambda);
}


// Called only once in cl_GetCoolingTimeAll
// function is fine, TODO check dependencies, write documentation
static double cl_GetCoolingRateFromUAll(coolingdata_in in, coolingdata_out *out)
{
  double lambda = 0;

#if GADGET_COOLING == COOLING_WAL
  out->Temperature = convert_u_to_tempMET(in.u_old, in.rho_cgs, in.DZ, &in.Metals[0]);
#ifdef GL_DUST_COOLING
  lambda = CoolingRateMET(out->Temperature, in.Redshift, in.DZ, in.DL, in.DS, &in.Metals[0]);
#else
  lambda = CoolingRateMET(out->Temperature, in.Redshift, in.DZ, &in.Metals[0]);
#endif
  return lambda;
#else
  out->Temperature = convert_u_to_tempSTD(in.u_old, in.rho_cgs, &out->ne_guess);
#endif

  CoolingRateArgs args;
  args.logT = log10(out->Temperature);
  args.rho = in.rho_cgs;
  args.nelec = &out->ne_guess;

#if GADGET_COOLING == COOLING_KATZ
  //  lambda = CoolingRateSTD(log10(out->Temperature), in.rho_cgs, &out->ne_guess);
  lambda = CoolingRate(&args);
#elif GADGET_COOLING == COOLING_SUTHERLAND
  //  lambda = CoolingRateZ(log10(out->Temperature), in.rho_cgs, &out->ne_guess, in.Z);
  args.Z = in.Z;
  lambda = CoolingRate(&args);
#elif GADGET_COOLING == COOLING_LUEDERS

#endif
  return lambda;
}


/**********************************************************************************************************************/
/*                                              End of Static Functions                                               */
/**********************************************************************************************************************/


/**********************************************************************************************************************/
/*                                                  Public Functions                                                  */
/**********************************************************************************************************************/

// Goal: These functions should serve as the public interface of the cooling code
// There is one function for each purpose that works differently under different circumstances

double cl_GetCooling()
{
  return 5;
}


double cl_GetCoolingRate()
{
  return 5;
}


double cl_get_CoolingTime()
{
  return 5;
}


/**********************************************************************************************************************/
/*                                              End of Public Functions                                               */
/**********************************************************************************************************************/



// called only once, in io.c
// No side effects, TODO check dependencies, what is DZ?
void particle2in_cooling(coolingdata_in *in, int i)
{
  double a3inv;
  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    {
      a3inv = 1.;
    }

  in->rho_cgs = SphP[i].Density * a3inv * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs;
  in->u_old = SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1) *
    All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

#if defined(METALS)
  in->Z = P[i].Metallicity;
#elif defined(LT_STELLAREVOLUTION)
  in->Z = get_metallicity_solarunits(get_metallicity(i, Iron));
#endif

#if GADGET_COOLING == COOLING_WAL
  double Redshift;
  if(All.ComovingIntegrationOn)
    {
      Redshift = 1.0 / All.Time - 1;	/* note that ComovingIntegration is alway on with wal cooling */
    }
  else
    {
      Redshift = 0.0;
    }

//  get_cool_redshift(Redshift, &in->DZ);
  in->DZ = lt_cl_GetCoolRedshift(Redshift);

  for(int k = 0; k < LT_NMet; k++)
    {
      in->Metals[k] = SphP[i].Metals[k];
    }
#endif
}


// Called only once, in io.c
// function is fine, TODO remove this comment, write documentation
double cl_GetCoolingTimeAll(coolingdata_in in, coolingdata_out *out)
{
  out->LambdaNet = cl_GetCoolingRateFromUAll(in, out);

  // Do we actually have heating due to UV background?
  if(out->LambdaNet >= 0)
    {
      return 0;
    }

  double nHcgs = XH * in.rho_cgs / PROTONMASS;	/* hydrogen number dens in cgs units */
  double ratefact = nHcgs * nHcgs / in.rho_cgs;
  double coolingtime = in.u_old / (-ratefact * out->LambdaNet);
  coolingtime *= All.HubbleParam / All.UnitTime_in_s;

  return coolingtime;
}


// Called in lots of places, mostly cooling.c
// function is fine, TODO check dependencies, write documentation
double DoCooling(coolingdata_in in, coolingdata_out *out, double dt)
{
  double unew = in.u_old;
#if GADGET_COOLING == COOLING_KATZ
  unew = DoCoolingSTD(in.u_old, in.rho_cgs, dt, &out->ne_guess);
#elif GADGET_COOLING == COOLING_SUTHERLAND
  unew = DoCoolingZ(in.u_old, in.rho_cgs, dt, &out->ne_guess, in.Z, &out->Temperature);
#elif GADGET_COOLING == COOLING_WAL
#ifdef GL_DUST_COOLING
  unew = DoCoolingMET(in.u_old, in.rho_cgs, &in.Metals[0], in.DL, in.DS, in.Redshift, in.DZ, dt, &out->Temperature);
#else
  unew = DoCoolingMET(in.u_old, in.rho_cgs, &in.Metals[0], in.Redshift, in.DZ, dt, &out->Temperature);
#endif
#elif GADGET_COOLING == COOLING_LUEDERS

#endif
  return unew;
}


// TODO Has side effects to chemistry noneq code, but no dependencies. Mostly fine.
void IonizeParams(void)
{
  if(!All.ComovingIntegrationOn)
    {
      SetZeroIonization();
      return;
    }

  double redshift = 1 / All.Time - 1;
  double logz = log10(redshift + 1.0);

  int ilow = 0;
  for(int i = 0; i < nheattab; i++)
    {
      if(inlogz[i] < logz)
	ilow = i;
      else
	break;
    }

  double dzlow = logz - inlogz[ilow];
  double dzhi = inlogz[ilow + 1] - logz;

  if(ThisTask == 0)
    {
      printf(" UM check: at redshift %g/%g: %d/%g/%g\n", redshift,
	     pow(10., inlogz[nheattab - 1]), nheattab, gH0[ilow], gH0[ilow + 1]);
    }

  if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
      SetZeroIonization();
#ifdef UM_CHEMISTRY
      gJH0uvb = gJHe0uvb = gJHepuvb = 0;
      epsH0uvb = epsHe0uvb = epsHepuvb = 0;
#endif
      return;
    }
  else
    {
      Cool.J_UV = 1.e-21;	/* irrelevant as long as it's not 0 */
    }

#ifdef UM_CHEMISTRY_NO_UV_BACKGROUND_IN_COSMO_EVOLUTION
  Cool.J_UV = 0.;		/* neglect UV background */
#endif

  Cool.gJH0 = JAMPL * pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ilow + 1])) / (dzlow + dzhi));
  Cool.gJHe0 = JAMPL * pow(10., (dzhi * log10(gHe[ilow]) + dzlow * log10(gHe[ilow + 1])) / (dzlow + dzhi));
  Cool.gJHep = JAMPL * pow(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ilow + 1])) / (dzlow + dzhi));
  Cool.epsH0 = JAMPL * pow(10., (dzhi * log10(eH0[ilow]) + dzlow * log10(eH0[ilow + 1])) / (dzlow + dzhi));
  Cool.epsHe0 = JAMPL * pow(10., (dzhi * log10(eHe[ilow]) + dzlow * log10(eHe[ilow + 1])) / (dzlow + dzhi));
  Cool.epsHep = JAMPL * pow(10., (dzhi * log10(eHep[ilow]) + dzlow * log10(eHep[ilow + 1])) / (dzlow + dzhi));

#ifdef UM_CHEMISTRY
  gJH0uvb = Cool.gJH0;
  gJHe0uvb = Cool.gJHe0;
  gJHepuvb = Cool.gJHep;
  epsH0uvb = Cool.epsH0;
  epsHe0uvb = Cool.epsHe0;
  epsHepuvb = Cool.epsHep;

  if(ThisTask == 0)
    {
      printf(" **** interpolated gJH0uvb = %g, gJHe0uvb = %g, gJHepuvb = %g \n", gJH0uvb, gJHe0uvb, gJHepuvb);
      printf(" **** interpolated epsH0uvb = %g, epsHe0uvb = %g, epsHepuvb = %g\n\n", epsH0uvb, epsHe0uvb,
	     epsHepuvb);
    }
#endif

  return;
}


void SetZeroIonization(void)
{
  Cool.gJHe0 = Cool.gJHep = Cool.gJH0 = 0;
  Cool.epsHe0 = Cool.epsHep = Cool.epsH0 = 0;
  Cool.J_UV = 0;
}


// called only from run.c
/*!
 * normal cooling routine when star formation is disabled
 */
 // TODO Needs more work.
void cl_CoolingOnly(void)
{
  double ascale = 1;
#ifdef LT_METAL_COOLING_WAL
  if(All.ComovingIntegrationOn)
    {
      Redshift = 1.0 / ascale - 1;
    }
  else
    {
      Redshift = 0.0;
    }

//  get_cool_redshift(Redshift, &DZ);
  DZ = lt_cl_GetCoolRedshift(Redshift);

  if(LT_NMet < get_cool_n_el())
    {
      endrun(1009);
    }
  WalCool_tables_load(Redshift);
#endif

  double a3inv, hubble_a = 0, time_hubble_a;
  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
    }
  else
    {
      a3inv = time_hubble_a = hubble_a = ascale = 1;
    }

  for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef BLACK_HOLES
      if(P[i].Type == 0 && P[i].Mass > 0)
#else
      if(P[i].Type == 0)
#endif
	{
	  double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */

	  double dtime = dt;
	  if(All.ComovingIntegrationOn)
	    {
	      dtime = All.Time * dt / time_hubble_a;
	    }

	  // electron abundance (gives ionization state and mean molecular weight)
	  double ne = SphP[i].elec;

	  double fac_entr_to_u = pow(SphP[i].Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
	  double uold = DMAX(All.MinEgySpec, SphP[i].Entropy * fac_entr_to_u);
	  double Zcool;

#ifdef LT_METAL_COOLING_WAL
	  double Metallicities[LT_NMet];

	  set_metallicities(i, &Metallicities[0], a3inv);	/* note: sets the physical density too! */
	  Zcool = get_metallicity_solarunits(get_metallicity(i, Iron));
#elif defined(LT_METAL_COOLING)
	  Zcool = get_metallicity_solarunits(get_metallicity(i, Iron));
#endif

	  double unew;

#ifdef UM_CHEMISTRY

	  double ti_step = P[i].TimeBin ? (1 << P[i].TimeBin) : 0;
	  double tstart = P[i].Ti_begstep + ti_step / 2;	/* midpoint of old step */
	  double tend = P[i].Ti_begstep + ti_step + ti_step / 2;	/* midpoint of new step */
	  double a_start, a_end;
	  if(All.ComovingIntegrationOn)
	    {
	      a_start = All.TimeBegin * exp(tstart * All.Timebase_interval);
	      a_end = All.TimeBegin * exp(tend * All.Timebase_interval);
	    }
	  else
	    {
	      a_start = tstart * All.Timebase_interval;
	      a_end = tend * All.Timebase_interval;
	    }

	  /* flag=1: cooling */
#ifdef LT_METAL_COOLING_WAL
	  unew =
	    Um_DoCooling(uold, SphP[i].Density * a3inv, dtime, &ne, &Metallicities[0], Redshift, DZ, i, 1);
#elif defined(LT_METAL_COOLING)
	  unew = Um_DoCooling(uold, SphP[i].Density * a3inv, dtime, &ne, Zcool, i, flag);
#else
	  unew = Um_DoCooling(uold, SphP[i].Density * a3inv, dtime, &ne, i, 1);
#endif // LT_METAL_COOLING_WAL


#else // UM_CHEMISTRY


	  double Temperature;
#ifndef RT_COOLING_PHOTOHEATING
#ifdef LT_METAL_COOLING_WAL
#ifdef GL_DUST_COOLING
	  unew =
	    DoCoolingMET(uold, SphP[i].Density * a3inv, &Metallicities[0], DL, DS, Redshift, DZ, dtime, &Temperature);
#else
	  unew =
	    DoCoolingMET(uold, SphP[i].Density * a3inv, &Metallicities[0], Redshift, DZ, dtime, &Temperature);
#endif
#elif defined(LT_METAL_COOLING)
	  unew = DoCoolingZ(uold, SphP[i].Density * a3inv, dtime, &ne, Zcool, &Temperature);
#else
	  unew = DoCoolingSTD(uold, SphP[i].Density * a3inv, dtime, &ne);
#endif // LT_METAL_COOLING_WAL
#else // RT_COOLING_PHOTOHEATING
	  unew = uold + dt * fac_entr_to_u * (rt_DoHeating(i, dt) + rt_DoCooling(i, dt));
#endif // RT_COOLING_PHOTOHEATING

#endif // UM_CHEMISTRY

#if GADGET_HYDRO == HYDRO_MFM
	  SphP[i].InternalEnergy = unew;
#endif
	  SphP[i].Entropy = unew / fac_entr_to_u;
	  SphP[i].elec = ne;

#if GADGET_HYDRO == HYDRO_MFM
	  SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
	  mfm_synchronize_entropy(i);
#endif
	  SphP[i].EntropyPred = SphP[i].Entropy;
	  SphP[i].Pressure = get_pressure(i);
	}
    }
}


#ifndef LT_METAL_COOLING_WAL

#define eV_to_K 11606.0
#define eV_to_erg 1.60184e-12


#ifdef LT_METAL_COOLING
double GetMetalLambda(double, double);
#endif

/*!
 * returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
 // TODO relatively little side effects? Needs work. Unify function header.
#ifdef LT_METAL_COOLING
double DoCoolingZ(double u_old, double rho, double dt, double *ne_guess, double Z, double *temp)
#else
double DoCoolingSTD(double u_old, double rho, double dt, double *ne_guess)
#endif				// LT_METAL_COOLING
{
  double du;
  int iter = 0;

  double DoCool_u_old_input = u_old;
  double DoCool_rho_input = rho;
  double DoCool_dt_input = dt;
  double DoCool_ne_guess_input = *ne_guess;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

  Cool.nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  double ratefact = Cool.nHcgs * Cool.nHcgs / rho;

  double u = u_old;
  double u_lower = u;
  double u_upper = u;

#ifdef LT_METAL_COOLING
  double LambdaNet = CoolingRateFromU(u, rho, ne_guess, Z, temp);
#else
  double LambdaNet = CoolingRateFromU(u, rho, ne_guess);
#endif

  /* bracketing */

  // heating
  if(u - u_old - ratefact * LambdaNet * dt < 0)
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
#ifdef LT_METAL_COOLING
      while(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess, Z, temp) * dt < 0)
#else
      while(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess) * dt < 0)
#endif
	{
	  u_upper *= 1.1;
	  u_lower *= 1.1;
	}
    }

  if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
#ifdef LT_METAL_COOLING
      while(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess, Z, temp) * dt > 0)
#else
      while(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess) * dt > 0)
#endif
	{
	  u_upper /= 1.1;
	  u_lower /= 1.1;
	}
    }

  do
    {
      u = 0.5 * (u_lower + u_upper);

#ifdef LT_METAL_COOLING
      LambdaNet = CoolingRateFromU(u, rho, ne_guess, Z, temp);
#else
      LambdaNet = CoolingRateFromU(u, rho, ne_guess);
#endif

      if(u - u_old - ratefact * LambdaNet * dt > 0)
	{
	  u_upper = u;
	}
      else
	{
	  u_lower = u;
	}

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("u= %g\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in DoCooling()\n");
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= "
	     "%g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      endrun(10);
    }

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;	/* to internal units */

  return u;
}

/*!
 * returns cooling time.
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
 // TODO fix header, figure out dependencies
#ifdef LT_METAL_COOLING
double GetCoolingTimeZ(double u_old, double rho, double *ne_guess, double Z, double *temp)
#else
double GetCoolingTimeSTD(double u_old, double rho, double *ne_guess)
#endif
{
  double coolingtime;

  double DoCool_u_old_input = u_old;
  double DoCool_rho_input = rho;
  double DoCool_ne_guess_input = *ne_guess;

  // convert to physical cgs units
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  Cool.nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  double ratefact = Cool.nHcgs * Cool.nHcgs / rho;

  double u = u_old;

#ifdef LT_METAL_COOLING
  double LambdaNet = CoolingRateFromU(u, rho, ne_guess, Z, temp);
#else
  double LambdaNet = CoolingRateFromU(u, rho, ne_guess);
#endif

  /* bracketing */

  if(LambdaNet >= 0)		/* ups, we have actually heating due to UV background */
    {
      return 0;
    }

  coolingtime = u_old / (-ratefact * LambdaNet);
  coolingtime *= All.HubbleParam / All.UnitTime_in_s;

  return coolingtime;
}


// called in too many places
/*!
 * this function determines the electron fraction, and hence the mean
 * molecular weight. With it arrives at a self-consistent temperature.
 * Element abundances and the rates for the emission are also computed
 */
double convert_u_to_tempSTD(double u, double rho, double *ne_guess)
{
  // TODO i think these intermediate variables can be safely removed
  double u_input = u;
  double rho_input = rho;
  double ne_input = *ne_guess;

  double mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
  double temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  double max = 0;
  double temp_old;
  int iter = 0;
  do
    {
      double ne_old = *ne_guess;

      find_abundances_and_rates(log10(temp), rho, ne_guess);
      temp_old = temp;

      mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);

      double temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      max =
	DMAX(max,
	     temp_new / (1 + yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

      temp = temp_old + (temp_new - temp_old) / (1 + max);
      iter++;

      if(iter > (MAXITER - 10))
#if defined(LT_SMOOTH_Z) || defined(LT_METAL_COOLING)
	if(!ignore_failure_in_convert_u)
#endif
	  printf("-> temp= %g ne=%g\n", temp, *ne_guess);
    }
  while(fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);

  if(iter >= MAXITER)
    {
#if defined(LT_SMOOTH_Z) || defined(LT_METAL_COOLING)
      if(!ignore_failure_in_convert_u)
#endif
	{
	  printf("failed to converge in convert_u_to_tempSTD()\n");
	  printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input, ne_input);
//      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= "
//             "%g\nDoCool_ne_guess_input= %g\n",
//             DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input,
//             DoCool_ne_guess_input);

	  endrun(12);
	}
    }

  return temp;
}


/*!
 * this function computes the actual abundance ratios
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess)
{
  double logT_input = logT;
  double rho_input = rho;
  double ne_input = *ne_guess;

  if(logT <= Cool.Tmin)		/* everything neutral */
    {
      Cool.nH0 = 1.0;
      Cool.nHe0 = yhelium;
      Cool.nHp = 0;
      Cool.nHep = 0;
      Cool.nHepp = 0;
      Cool.ne = 0;
      *ne_guess = 0;
      return;
    }

  if(logT >= Cool.Tmax)		/* everything is ionized */
    {
      Cool.nH0 = 0;
      Cool.nHe0 = 0;
      Cool.nHp = 1.0;
      Cool.nHep = 0;
      Cool.nHepp = yhelium;
      Cool.ne = Cool.nHp + 2.0 * Cool.nHepp;
      *ne_guess = Cool.ne;	/* note: in units of the hydrogen number density */
      return;
    }

  double t = (logT - Cool.Tmin) / Cool.deltaT;
  int j = (int) t;
  double fhi = t - j;
  double flow = 1 - fhi;

  if(*ne_guess == 0)
    {
      *ne_guess = 1.0;
    }

  Cool.nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

  Cool.ne = *ne_guess;
  double neold = Cool.ne;
  Cool.necgs = Cool.ne * Cool.nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  int niter = 0;
  do
    {
      niter++;

      Cool.aHp = flow * Cool.AlphaHp[j] + fhi * Cool.AlphaHp[j + 1];
      Cool.aHep = flow * Cool.AlphaHep[j] + fhi * Cool.AlphaHep[j + 1];
      Cool.aHepp = flow * Cool.AlphaHepp[j] + fhi * Cool.AlphaHepp[j + 1];
      Cool.ad = flow * Cool.Alphad[j] + fhi * Cool.Alphad[j + 1];
      Cool.geH0 = flow * Cool.GammaeH0[j] + fhi * Cool.GammaeH0[j + 1];
      Cool.geHe0 = flow * Cool.GammaeHe0[j] + fhi * Cool.GammaeHe0[j + 1];
      Cool.geHep = flow * Cool.GammaeHep[j] + fhi * Cool.GammaeHep[j + 1];

      if(Cool.necgs <= 1.e-25 || Cool.J_UV == 0)
	{
	  Cool.gJH0ne = Cool.gJHe0ne = Cool.gJHepne = 0;
	}
      else
	{
	  Cool.gJH0ne = Cool.gJH0 / Cool.necgs;
	  Cool.gJHe0ne = Cool.gJHe0 / Cool.necgs;
	  Cool.gJHepne = Cool.gJHep / Cool.necgs;
	}

      // Katz95 eqs. 33 and 34
      Cool.nH0 = Cool.aHp / (Cool.aHp + Cool.geH0 + Cool.gJH0ne);
      Cool.nHp = 1.0 - Cool.nH0;	/* eqn (34) */

      if((Cool.gJHe0ne + Cool.geHe0) <= SMALLNUM)	/* no ionization at all */
	{
	  Cool.nHep = 0.0;
	  Cool.nHepp = 0.0;
	  Cool.nHe0 = yhelium;
	}
      else
	{
	  // Katz95 eqs. 35-37
	  Cool.nHep = yhelium / (1.0 + (Cool.aHep + Cool.ad) / (Cool.geHe0 + Cool.gJHe0ne) +
				 (Cool.geHep + Cool.gJHepne) / Cool.aHepp);
	  Cool.nHe0 = Cool.nHep * (Cool.aHep + Cool.ad) / (Cool.geHe0 + Cool.gJHe0ne);
	  Cool.nHepp = Cool.nHep * (Cool.geHep + Cool.gJHepne) / Cool.aHepp;
	}

      neold = Cool.ne;
      // Katz95 eq. 38
      Cool.ne = Cool.nHp + Cool.nHep + 2 * Cool.nHepp;

      Cool.necgs = Cool.ne * Cool.nHcgs;

      if(Cool.J_UV == 0)
	break;

      Cool.ne = 0.5 * (Cool.ne + neold);
      Cool.necgs = Cool.ne * Cool.nHcgs;

      if(fabs(Cool.ne - neold) < 1.0e-4)
	break;

      if(niter > (MAXITER - 10))
	printf("ne= %g  niter=%d\n", Cool.ne, niter);
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("no convergence reached in find_abundances_and_rates()\n");
      printf("logT_input= %g  rho_input= %g  ne_input= %g\n", logT_input, rho_input, ne_input);
//    printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= "
//           "%g\nDoCool_ne_guess_input= %g\n",
//           DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input,
//           DoCool_ne_guess_input);
      endrun(13);
    }

  Cool.bH0 = flow * Cool.BetaH0[j] + fhi * Cool.BetaH0[j + 1];
  Cool.bHep = flow * Cool.BetaHep[j] + fhi * Cool.BetaHep[j + 1];
  Cool.bff = flow * Cool.Betaff[j] + fhi * Cool.Betaff[j + 1];

  *ne_guess = Cool.ne;
}


/*!
 * this function first computes the self-consistent temperature
 * and abundance ratios, and then it calculates
 * (heating rate-cooling rate)/n_h^2 in cgs units
 * TODO is that description even accurate?
 */
 // TODO Unify function headers
#ifdef LT_METAL_COOLING
double CoolingRateFromU(double u, double rho, double *ne_guess, double Z, double *temp)
#else
double CoolingRateFromU(double u, double rho, double *ne_guess)
#endif
{
#ifndef LT_METAL_COOLING
  double temp = convert_u_to_tempSTD(u, rho, ne_guess);
#endif
  CoolingRateArgs args;

#ifndef LT_METAL_COOLING
  args.logT = log10(temp);
#else
  args.logT = log10(*temp);
#endif
  args.rho = rho;
  args.nelec = ne_guess;

#ifdef LT_METAL_COOLING
  args.Z = Z;
#endif

  return CoolingRate(&args);

//#ifdef LT_METAL_COOLING
//  return CoolingRateZ(log10(temp), rho, ne_guess, Z);
//#else
//  return CoolingRateSTD(log10(temp), rho, ne_guess);
//#endif
}


/*!
 * this function computes the self-consistent temperature
 * and abundance ratios
 */
// This function is fine, mostly just calls taht other function - TODO remove comment
double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer)
{
  double temp;

//  DoCool_u_old_input = u;
//  DoCool_rho_input = rho;
//  DoCool_ne_guess_input = *ne_guess;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  temp = convert_u_to_tempSTD(u, rho, ne_guess);

  *nH0_pointer = Cool.nH0;
  *nHeII_pointer = Cool.nHep;
  *ne_guess = Cool.ne;

  return temp;
}




#ifdef LT_METAL_COOLING
/*  M E T A L   C O O L I N G  */

double GetMetalLambda(double logT, double Zsol)
{
  double t, u, log_temp = logT;
  double x_1y, x_1y_1, xy_1;
  int Zi, Ti;

  /* note: Zsol is supposed to be the Iron abundance in solar units
     (anders and grevesse 1989: 2.62e-3)
   */

  if(getindex(&CoolZvalue[0], 0, ZBins - 1, &Zsol, &Zi) == -2)
    return 0;
  if(Zi == ZBins - 1)
    t = 0;
  else
    t = (Zsol - CoolZvalue[Zi]) / (CoolZvalue[Zi + 1] - CoolZvalue[Zi]);

  if(getindex(&CoolTvalue[0], 0, TBins - 1, &log_temp, &Ti) == -2)
    return 0;
  if(Ti == TBins - 1)
    u = 0;
  else
    u = (log_temp - CoolTvalue[Ti]) / (CoolTvalue[Ti + 1] - CoolTvalue[Ti]);

  if(t * u > 0)
    x_1y_1 = CoolingTables[Zi + 1][Ti + 1];
  else
    x_1y_1 = 0;

  if(t == 0)
    x_1y = 0;
  else
    x_1y = CoolingTables[Zi + 1][Ti];

  if(u == 0)
    xy_1 = 0;
  else
    xy_1 = CoolingTables[Zi][Ti + 1];

  return pow(10,
	     (1 - t) * (1 - u) * CoolingTables[Zi][Ti] + t * (1 - u) * x_1y + t * u * x_1y_1 + (1 -
												t) * u *
	     xy_1);
}
#endif

#endif /* closes LT_METAL_COOLING_WAL */

#endif // COOLING
