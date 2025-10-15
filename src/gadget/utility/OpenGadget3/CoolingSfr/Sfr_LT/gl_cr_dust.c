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


// GLG 240117 "La programmazione non e' una corsa ma un tiro al bersaglio" (parafrasi da Susanna Tamaro)
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"
#ifdef GM_MUPPI
/*#include "../../CoolingSfr/Muppi/muppi.h"*/
#include "../../CoolingSfr/Muppi/muppi_proto.h"
/* #include "../CoolingSfr/Muppi/muppi_kicks.h" GM: to be implemented */
#endif


#if defined(GL_GIVENACC) && !defined(GL_GIVENGRAIN)	// GL150218
#error nonsense to use GL_GIVENACC without GL_GIVENGRAIN
#endif

#if defined(GL_CR_ALMOSTGIVENACC) && !defined(GL_GIVENACC)	// GL150218
#error nonsense to use GL_CR_ALMOSTGIVENACC without GL_GIVENACC
#endif


int Xind_C = 1;			// indices in dcondeff*
// GL080218 introduced definitions below one upon a time only for GIVENACC form sometime in 2020 allways needed
/* XXXXXXXXXXX NOTE NOTE NOTE !!!!!!!!!!!!!!!!
 in lt_sn.h defined same quantities with same name but without initial X.
 they are used for grain production in lt_sn.c
 if change here any of these quantities, DO THE SAME THERE AND VICEVERSA!!!!
*/
int Xind_O = 3;			// indices o oxigen in dcondeff*
int Xind_Mg = 6, Xind_Fe = 9;	// indices in dcondeff required only for ALMOSTGIVENACC*
int Xind_Si = 8;
// here "olivine" version Mg, Si, Fe
#define n_sil 3			// n_sil is the number of elements other than O entering in assumed silicate grain
int Xind_sil[n_sil] = { 6, 8, 9 };	// indeces in dcondeff* of elements other than O entering in assumed silicate grain
double Xmu_sil[n_sil] = { 24., 28., 56. };	// atomic weights
double Aacc_C = 2502., Aacc_sil_O = 1565., Aacc_sil[n_sil] = { 733., 911., 2523. };	// coefficients for accretion timescale of gas oto grains

// Aacc computed woth the following piece of IDL code. Uses Eq. 20 of Hirashita & Kuo 2011. with nx=n*(M_x/M)*<mu>/mu_x
/*pro acc_fac, metal
  if (metal eq 'C') then begin
      fx=1.0   & mu=12. & s_sma=2.26
  endif
  if (metal eq 'O') then begin
      fx=0.371   & mu=16.    & s_sma=3.3
  endif
  if (metal eq 'Mg') then begin
      fx=0.141   & mu=24.3  & s_sma=3.3 ; Mg
  endif
  if (metal eq 'Si') then begin
      fx=0.163   & mu=28.11 & s_sma=3.3 ; Si
  endif
  if (metal eq 'Fe') then begin
      fx=0.32    & mu=55.9 & s_sma=3.3 ; Fe
  endif
  a=5e-3 * 1e-4
  HYDROGEN_MASSFRAC = 0.7
  mu_c = (4./ (3. * HYDROGEN_MASSFRAC + 1.) )
  S_big=0.3
  n=1e3
  T=50
  k=1.38D-16
  mh=1.66d-24

  mx=mu*mh
  coef=a*fx*s_sma*mu/(sqrt(mx)*n*mu_c*S_big) * sqrt(2*!Pi/k/T)

  coef_yr=coef/(365.25*24*3600)

  print, metal, coef_yr, coef_yr/3
end
*/
int Xn_ato_sil[n_sil] = { 1, 1, 1 }, Xn_ato_O = 4.;	// number of atoms in assumed silicate c. e.g. 1 1 1 4 for O in standard olivine

  // NOTE that if GL_GIVENGRAIN is defined, only the C elements of the arrays dcondeffX are used.
double Xmu_O = 16.0, Xmu_C = 12.0;	// mus are atomic weights


#define mu_h  ((double)4./ (5. * HYDROGEN_MASSFRAC + 3.) )	/* molecular weight of the hot phase */
#define mu_c  ((double)4./ (3. * HYDROGEN_MASSFRAC + 1.) )	/* molecular weight of the cold phase */
#define SOLAR_MASS_IN_CODE_UNITS  (SOLAR_MASS / All.UnitMass_in_g)

#if !(defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST))
#error Dust evolution GL_DUST requires LT_STELLAREVOLUTION!
#endif


#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)

#define YEAR_IN_CODE_UNITS        (SEC_PER_YEAR / All.UnitTime_in_s)

/* START OF DUST EVOLUTION FUNCTION */
void evolve_dust(void)
{
  int i, j;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, time_hubble_a, rhomean, tstart, tend;
  double taumax = 1e11;		// GLG 270219 to improve speed, I will set to 0 timescales above a maximum, and computations are then avoided
  double facsn;			// GLG 270219
  double Ztot, Zgas, Mgc, D, Dl, Ds, taush, tauco, tau, tauac, tauspu_s, tauspu_l, eta, Nsne, msw;
  double XCold, M_h, M_c, Vol_tot, T_h, E_h, n_g, xn, n_hot, dum;
  double Mdl[LT_NMetP], Mds[LT_NMetP], dtDl[LT_NMetP], dtDs[LT_NMetP];

  double T_d = 2.e6, w_d = 2.5;	//constants entering into the formula for sputtering

  double fdense = 0.5;		// fiducial for fraction of gas particles dense enough for accretion and coagulation.
  // this constant value never used in practise.

// GLG & CR 160217, revised 210820, previously base on wrong interpretation by Rosa. PD.
// following quantities are used to compute an fdense depending on particle density (NOT with Muppi, only effective model)
// based on Wada & Norman 07 results
  double sig_pdf = 2.5;		// sigma of assumed PDF for unresolved gas densities
  double n_th = 1e3;		// assumed threshold number density for accretion and coagulation (unresolved! used in PDF)
  double n_0 = 2.12, arg_erf;	// normalization of PDF and auxiliary argument to be put in erf

  double Metals_accretion;
  double mwk, T_wk;

  tstart = second();

  if(All.ComovingIntegrationOn)
    {				/* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
      rhomean = 3 * hubble_a * hubble_a / 8 / M_PI / All.G;	/* mean density at time a in unit's code */
    }
  else
    {
      a3inv = ascale = time_hubble_a = 1;
      rhomean = 1;
    }

  /* this it used to get the hot phase numerical density in the effective model */
  T_h = All.TempSupernova / All.FactorEVP;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      /* SphP[i].d.Density -> density
         SphP[i].Temperature -> temperature
         P[i].Mass -> mass
         SphP[i].Metals[1...15] -> mass of metals (code units)
         SphP[i].numSnIa, SphP[i].numSnII -> number of supernovae
         SphP[i].DustL[1...15], SphP[i].DustS[1...15] -> large and small dust grains
       */

      if(P[i].Type == 0)
	{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */
	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

	  XCold = 0;
	  Ztot = 0.;
	  Zgas = 0.;
	  D = 0.;
	  Dl = 0.;
	  Ds = 0.;
	  tauac = 0.;
	  taush = 0.;
	  tauco = 0.;
	  tauspu_l = 0.;
	  tauspu_s = 0.;
	  Metals_accretion = 0.;
	  Mgc = P[i].Mass;

#ifdef GM_MUPPI			// GL 190219 we are running with MUPPI
	  /* Converting relevant quantities into h true */
// XXXXXX  GLG & CR 240620 ci siamo accorti che in compute_ism_properties le cose vanno passate con h calcolato e nel calcolo del Vol_tot va sottratto massa andata in stelle.
// l'errore sembra ereditato da caso modello effettivo, 040920, inoltre fino a 190920 a3inv missing. Non ho parole

	  Vol_tot = Mgc / SphP[i].Density / a3inv / All.HubbleParam / All.HubbleParam / All.HubbleParam;	//GL 0423
	  M_c = SphP[i].M_c / All.HubbleParam;	// tutte ste cose all'uscita di Muppi_main sono rimesse in unità con h, quindi qua le ripuliamo di nuovo
	  M_h = SphP[i].M_h / All.HubbleParam;
	  E_h = SphP[i].E_h / All.HubbleParam;
	  T_wk = SphP[i].Temperature;	// GL 0423, funalmente ora in OG MUppi si usa questo, non un suo campo -Temp
	  if(Mgc > 0)
	    {
	      XCold = M_c / Mgc * All.HubbleParam;
	    }
	  else
	    {
	      XCold = 0;
	    }
#else
	  /* these are used to get the hot phase numerical density in the effective model */
	  XCold = SphP[i].XColdCloud;
	  Vol_tot = Mgc / SphP[i].Density / a3inv / All.HubbleParam / All.HubbleParam / All.HubbleParam;	//GLG before 190920 a3inv missing
	  M_c = XCold * Mgc / All.HubbleParam;
	  M_h = (1. - XCold) * Mgc / All.HubbleParam;
	  T_wk = SphP[i].Temperature;	// GLG 190118 but was inside sputtering computation before 250219
#endif
	  n_g = SphP[i].Density * a3inv * All.UnitDensity_in_cgs / mu_h / PROTONMASS * All.HubbleParam * All.HubbleParam;	// GLG before 270620 h was missing!!!! and before 190920 also a3inv


#ifdef GL_GIVENACC		// GL 080218
	  // find the element that constraints accretion to preserve given grain atom proportions
	  double wk, nummolwk;
	  int ind_acc, ind_acc_sil;
	  nummolwk = SphP[i].Metals[Xind_O] / pow(Xmu_O, 1.5) / Xn_ato_O;
	  ind_acc = Xind_O;
	  ind_acc_sil = -1;
	  for(j = 0; j < n_sil; j++)
	    {
#ifdef GL_CR_ALMOSTGIVENACC	// GL 170719
	      if(Xind_sil[j] == Xind_Mg || Xind_sil[j] == Xind_Fe)
		{
		  wk = (SphP[i].Metals[Xind_Mg] / Xmu_sil[Xind_Mg] * A_acc_sil[Xind_Mg];
			+SphP[i].Metals[Xind_Fe] / Xmu_sil[Xind_Fe] * A_acc_sil[Xind_Fe]) / (Xn_ato_sil[Xind_Mg] + Xn_ato_sil[Xind_Fe]);	// XXX To be tested...
		  // not sure the factor A_acc_sil[] used like this works perfectly in this case ALMOSTGIVENACC
		  // still probably the best we can do.
		}
	      else
		{
		  wk = SphP[i].Metals[Xind_sil[j]] / Xmu_sil[j] / Xn_ato_sil[j] * A_acc_sil[j];;
		}
#else
	      wk = SphP[i].Metals[Xind_sil[j]] / pow(Xmu_sil[j], 1.5) / Xn_ato_sil[j];
#endif /* GL_CR_ALMOSTGIVENACC */
	      if(wk < nummolwk)
		{
		  nummolwk = wk;
		  ind_acc = Xind_sil[j];
		  ind_acc_sil = j;
		}
	    }
#endif /* GL_GIVENACC */


	  for(j = 1; j < LT_NMetP; j++)
	    {
	      dtDl[j] = 0.;
	      dtDs[j] = 0.;
	      Ztot += SphP[i].Metals[j] + SphP[i].DustL[j] + SphP[i].DustS[j];
	      Zgas += SphP[i].Metals[j];

	      Mdl[j] = SphP[i].DustL[j];
	      Mds[j] = SphP[i].DustS[j];
	      if(P[i].Mass != 0)
		{
		  D += (Mdl[j] + Mds[j]);
		  Dl += Mdl[j];
		  Ds += Mds[j];
		}		// CLOSES if (P[i].Mass != 0)
	    }			// CLOSES for(j=1;j<LT_NMet; j++)


	  if(P[i].Mass != 0)
	    {
	      Ztot /= Mgc;	//this is a METALLICITY not the total metal mass
	      Zgas /= Mgc;	//this is a METALLICITY not the total metal mass
	      Dl /= Mgc;
	      Ds /= Mgc;
	      D /= Mgc;
	    }

	  if(D > 1. || Dl > 1. || Ds > 1. || Ztot > 1. || Zgas > 1)
	    {
	      printf(" DUST ERROR. Dl=%e, Ds=%e, D=%e Ztot=%g Zgas=%g Mg=%g Task %d part %d\n", Dl, Ds, D,
		     Ztot, Zgas, Mgc, ThisTask, i);
	      fflush(stdout);
	    }

	  //========================= Calculation of timescales

	  // GLG 140217 moved outside following if (Dl...)
	  if(XCold > 0)		//particle is multi-phase
	    {

#ifdef GM_MUPPI			// GL 190219 we are running with MUPPI
	      compute_ism_properties(M_h, E_h, M_c, Vol_tot, &T_h, &dum, &dum, &dum, &n_hot, &xn);
	      //double Fcoll = molecular_fraction(n_hot*T_h);
	      fdense = SphP[i].Fcoll;	/* BR or Kr or H2_from_Dust+He *///GL 0423
//                fdense = SphP[i].M_H2/All.HubbleParam/HYDROGEN_MASSFRAC/M_c; /* BR or Kr or H2_from_Dust+He*/

#else
	      compute_ism_properties_effective_model(M_h, M_c, Vol_tot, T_h, &dum, &dum, &dum, &n_hot, &xn);
// GLG 200820 proper interpretation of Wada and Norman 07
	      arg_erf = (log(n_th / n_0) / sig_pdf - sig_pdf) / 1.4142;
	      fdense = 0.5 * (1.0 - erf(arg_erf));	// GLG & CR 160217
	      arg_erf = (log(xn / n_0) / sig_pdf - sig_pdf) / 1.4142;
	      fdense /= 0.5 * (1.0 - erf(arg_erf));

// old version, based on Rosa's misenterpretation of Wada and Norman 07 results!
//                n_0 = xn/exp( 0.5* sig_pdf * sig_pdf );
//                arg_erf = ( log(n_th/n_0) / sig_pdf - sig_pdf )  / 1.4142;
//                fdense = 0.5 * (1.0 -erf(arg_erf)); // GLG & CR 160217
#endif




//printf(" %e %e %e %e %e %e %e %e %e GLG DEBUG M_c, M_h, T_wk, XCold, n_hot, xn, fdense, Fcoll, rat\n",M_c, M_h, T_wk, XCold, n_hot, xn, fdense, Fcoll, fdense/Fcoll);
//fflush(stdout);
	      T_wk = T_h;
	    }
	  else			//numerical density of the whole particle
	    {
	      xn = n_g;
	    }

	  if(Dl > 0 && xn < 1e3)
	    {
	      /* GLG July 2017 eliminated here additional condition  && xn <1. Replaced with a smooth transition
	         of the normalization (equation B4 from Aoyama+16) from v=10 km/s for xn<1 to v=0.1 km/s for xn/1000. This gives
	         the factor pow(xn,0.6666) at xn>1 just below. While this chage seems to be irrelevant in LR Dianoga runs, it is
	         fundamental with Muppi galaxy runs to get enough small grains */
	      taush = 5.41e7 * 0.01 / Dl / xn;	//Aoyama+2016 (23) //GL 240117
	      if(xn > 1)
		{
		  taush = taush * pow(xn, 0.6666);
		}

	      //xn: n_gas if not multiphase, n_h if multiphase
	      if(taush < taumax)	// GLG 280219 convert to code unit or nullify if too long
		{
		  taush = taush * YEAR_IN_CODE_UNITS;
		}
	      else
		{
		  taush = 0;
		}
	    }

	  if(Ds > 0 && XCold > 0. && fdense > 1e-8)	//only for MP particles
	    {			/* for 'dense' clouds.
				   vco is the typical velocity dispersion of small grains in [km/s]. */
// GLG 100217 now XCold included whenever fdense is used
	      double vco = 0.2;
	      tauco = (0.1 / vco) * 2.71e5 * (0.01 / Ds) / fdense / XCold;	//
	      if(tauco < taumax)	// GLG 280219 convert to code unit or nullify if too long
		{
		  tauco = tauco * YEAR_IN_CODE_UNITS;
		}
	      else
		{
		  tauco = 0;
		}
	    }



// GLG & CR from 240620 we apply sputtering also to multiphase particles, of course only to hot phase
//              {
	  // GLG 090217 inspired by a combination of Hirashita 2017 (eq 10) and Tsai & Mathews 1995 eq 14 and 15
	  // GLG 180118 adopting 0.05 micron as typical size of large grains, the mean for a power law -3.5
	  // distribution between 0.03 and 0.25 micron
	  // GLG 200417 both factors corrected. Before (run16) the former too small by 2, the latter too large by 100
// GLG 180118 corrected factor using a more proper mean molecular weight mu and )
	  if(T_wk > 3e7)
	    T_wk = 3e7;		// GLG190118
	  double n_forspu;
	  if(XCold <= 0.)
	    {
	      n_forspu = xn;	// no multiphase
	    }
	  else
	    {
	      n_forspu = n_hot;	// multiphase
	    }
	  tauspu_l = 5 * 1.65e7 / 3. * (0.05 / 0.1) * 0.01 / n_forspu * (pow(T_d / T_wk, w_d) + 1);	// factor 3 is because Tsai & Mathews give timescale for radius not mass
	  tauspu_s = 0.1 * tauspu_l;	// GLG 090217 adopting for small grains a size 10 times smaller than for large ones
	  if(tauspu_l < taumax)	// GLG 280219 convert to code unit or nullify if too long
	    {
	      tauspu_l = tauspu_l * YEAR_IN_CODE_UNITS;
	    }
	  else
	    {
	      tauspu_l = 0;
	    }

	  if(tauspu_s < taumax)	// GLG 280219 convert to code unit or nullify if too long
	    {
	      tauspu_s = tauspu_s * YEAR_IN_CODE_UNITS;
	    }
	  else
	    {
	      tauspu_s = 0;
	    }
//              }

	  //printf("      DDDDDDDDD Task %d part %d dtime %e taush %e taucc %e\n",
	  //ThisTask, i, dtime, taush, tauco); fflush(stdout);


	  //SN shocks
//SphP[i].numSnIa =0 ;
//#ifdef GM_MUPPI //for this hot fix we adopt coeff of kroupa93 with MUPPI and of Chabrier_pow with effective model
//SphP[i].numSnII = 0.00524 * SphP[i].Sfr * dtime/YEAR_IN_CODE_UNITS; // GLG 141020 (updated) hot fix aprox coefff valido per IMF papero MUPPI+dust
//#else
//SphP[i].numSnII = 0.01195 * SphP[i].Sfr * dtime/YEAR_IN_CODE_UNITS; // GLG 141020 hot fix aprox coefff valido per IMF papero MUPPI+dust
//#endif
	  Nsne = SphP[i].numSnIa + SphP[i].numSnII;
	  msw = 1360 * SOLAR_MASS_IN_CODE_UNITS / 4.;	// mass swept by a single blast in code units Aoymama+2016 (15)
	  // factor 4 because we take v_s=200 km/s
	  eta = 0.1 * msw / (Mgc / All.HubbleParam);	// we expect msw always << than mgas
	  facsn = (1 - pow(1 - eta, Nsne));

//if (Nsne !=0 ) printf(" %e %e %e %e %e %e %e GLG DEBUG dt SphP[i].numSnIa, SphP[i].numSnII, XCold, Mgc, eta, facsn \n",dtime/YEAR_IN_CODE_UNITS, SphP[i].numSnIa, SphP[i].numSnII, XCold, Mgc, eta, facsn);
//if (Nsne !=0 )
//if (SphP[i].numSnIa!=0 )
//printf(" %e %e %e %e %e %e %e %e %e %e %e GLG DEBUG ascale, X, Y, Z, dt, SnIa, SnII, XCold, Mgc, SFR, D\n",ascale, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], dtime/YEAR_IN_CODE_UNITS, SphP[i].numSnIa, SphP[i].numSnII, XCold, Mgc, SphP[i].Sfr, D);
//fflush(stdout);
//

//GLG 270219 now compute variations due to the various channels.
// Computation are avoided when timescales arrive here still 0, or set again to 0 because too long

	  for(j = 1; j < LT_NMetP; j++)
	    {

	      if(taush > 0)
		{
		  mwk = dtime * Mdl[j] / taush;	//GLG 230117 se non mi muovo famo mezzanotte e mi tocca cambiare data:-)
		  if(mwk > SphP[i].DustL[j])
		    mwk = SphP[i].DustL[j];	// GLG 230117 not enough large grains available to transfer to small
		  dtDs[j] += mwk;
		  dtDl[j] -= mwk;
		}

	      if(tauco > 0)
		{
		  mwk = dtime * Mds[j] / tauco;	//GLG 230117 se non mi muovo famo mezzanotte e mi tocca cambiare data:-)
		  if(mwk > SphP[i].DustS[j])
		    mwk = SphP[i].DustS[j];	// GLG 230117 not enough small grains available to transfer large
		  dtDl[j] += mwk;
		  dtDs[j] -= mwk;
		}


	      if(XCold > 0. && SphP[i].Metals[j] > 0. && fdense > 1e-8)	//GLG 230117 GLG & CR 270620 corrected misinterpretation of Aoyama formula
		// now rest directly on Uses Eq. 20 of Hirashita & Kuo 2011
		{
// GLG 100217 now XCold included whenever fdense is used
		  tauac = 1.0 * Mgc / fdense / XCold;	// this factors alsways present XXXX 3.33 fake GLG
#ifdef GL_GIVENACC		// GL 080218
		  if(j != Xind_C)
		    {
		      if(ind_acc_sil > -1)
			{
			  tauac = Aacc_sil[ind_acc_sil] * tauac / (SphP[i].Metals[ind_acc]);
			}
		      else
			{
			  tauac = Aacc_sil_O * tauac / (SphP[i].Metals[ind_acc]);
			}
		    }
		  else
		    {
		      tauac = Aacc_C * tauac / (SphP[i].Metals[j]);
		    }
#else
		  tauac = tauac / (SphP[i].Metals[j]);
		  if(j != Xind_C)
		    {
		      tauac = Aacc_sil[Xind_Si] * tauac;
		    }
		  else
		    {
		      tauac = Aacc_C * tauac;
		    }

#endif
		  if(tauac < taumax)	// GLG 280219 convert to code unit or nullify if too long
		    {
		      tauac = tauac * YEAR_IN_CODE_UNITS;
		    }
		  else
		    {
		      tauac = 0;
		    }
		}
	      else
		{
		  tauac = 0.;
		}

	      if(tauac > 0)
		{
		  Metals_accretion = dtime * Mds[j] / tauac;
		  if(Metals_accretion > XCold * fdense * SphP[i].Metals[j])	// GLG modified 24.01.22 includimg fdense & Xcold (thanks MP)
		    //GLG230115 check to be done here, to avoid creating small grains with metals that do not exist
		    {
		      printf
			(" DUST INCONSISTENCY Task %d Part %d wants to accrete a mass %e of metal %d, but only has %e\n",
			 ThisTask, i, Metals_accretion, j, SphP[i].Metals[j]);
		      fflush(stdout);
		      Metals_accretion = XCold * fdense * SphP[i].Metals[j];
		    }
		  /* decreasing metals amount for metals accretion on small dust grain */
		  SphP[i].Metals[j] -= Metals_accretion;
		  dtDs[j] += Metals_accretion;
		}




	      if(SphP[i].DustL[j] + dtDl[j] < 0.)
		{
		  printf
		    (" DUST INCONSISTENCY Task %d Part %d wants to decrease large dust %d to less than zero.\n",
		     ThisTask, i, j);
		  printf(" Present=%e decrement=%e dt=%e n=%e Temp=%e XCold=%e\n", SphP[i].DustL[j], dtDl[j],
			 dtime / YEAR_IN_CODE_UNITS, xn, T_wk, XCold);
		  fflush(stdout);
		  SphP[i].DustL[j] = 0.0;
		}
	      else
		SphP[i].DustL[j] += dtDl[j];	// *** CHECK INCREASING TIMESTEP ***

	      if(SphP[i].DustS[j] + dtDs[j] < 0.)
		{
		  printf
		    (" DUST INCONSISTENCY Task %d Part %d wants to decrease small dust %d to less than zero.\n",
		     ThisTask, i, j);
		  printf(" Present=%e decrement=%e dt=%e n=%e Temp=%e XCold=%e\n", SphP[i].DustS[j], dtDs[j],
			 dtime / YEAR_IN_CODE_UNITS, xn, T_wk, XCold);
		  fflush(stdout);
		  SphP[i].DustS[j] = 0.0;
		}
	      else
		SphP[i].DustS[j] += dtDs[j];	// *** float or double ??? ***



	      // CRF 171221. We make SN destruction and sputtering after the rest of the processes. 
	      // Since these two processes retun metals to the gas phase we need to be sure that the
	      // amount of metals returned to the gas phase is exactly the same than the destroyed dust.
	      // We also add here the metals restitution by SN (was not present)
	      mwk = facsn * SphP[i].DustL[j];
	      if(mwk > SphP[i].DustL[j])
		mwk = SphP[i].DustL[j];
	      SphP[i].DustL[j] -= mwk;
	      SphP[i].Metals[j] += mwk;

	      mwk = facsn * SphP[i].DustS[j];
	      if(mwk > SphP[i].DustS[j])
		mwk = SphP[i].DustS[j];
	      SphP[i].DustS[j] -= mwk;
	      SphP[i].Metals[j] += mwk;


	      double Xhot;
	      Xhot = 1. - XCold;
	      if(tauspu_l > 0)	// GLG 070217
		{
		  mwk = dtime * SphP[i].DustL[j] * Xhot / tauspu_l;
		  if(mwk > Xhot * SphP[i].DustL[j])
		    mwk = Xhot * SphP[i].DustL[j];	//  enough large grains available
		  SphP[i].Metals[j] += mwk;
		  SphP[i].DustL[j] -= mwk;
		}


	      if(tauspu_s > 0)	// GLG 070217
		{
		  mwk = dtime * SphP[i].DustS[j] * Xhot / tauspu_s;
		  if(mwk > Xhot * SphP[i].DustS[j])
		    mwk = Xhot * SphP[i].DustS[j];	//  enough large grains available
		  SphP[i].Metals[j] += mwk;
		  SphP[i].DustS[j] -= mwk;
		}



	      if(SphP[i].Metals[j] < 0.0)
		{
		  /* note, this should not happen any more */
		  printf
		    (" EG-WARNING, hot fix on metals! Task %d Time %f Metal %d Value %g change %g dtime %g  Mds[j] %g \n",
		     ThisTask, All.Time, j, SphP[i].Metals[j], Metals_accretion, dtime, Mds[j]);
		  fflush(stdout);
		  printf("             Mdl %g Mgc %g Ztot %g  D %g Dl %g Ds %g\n", Mdl[j], Mgc, Ztot, D, Dl,
			 Ds);
		  fflush(stdout);

		  SphP[i].Metals[j] = 0.0;
		}


	      if(D > 1. || Dl > 1. || Ds > 1. || Ztot > 1. || Zgas > 1)
		{
		  printf
		    (" DUST ERROR. Dl=%e, Ds=%e, dtDl=%e dtDs=%g Met=%g Ztot=%g Zgas=%g Mg=%g Xhot=%g tauS=%g tauL=%g Task %d part %d j %d\n",
		     SphP[i].DustL[j], SphP[i].DustS[j], dtDl[j], dtDs[j], SphP[i].Metals[j], Ztot, Zgas, Mgc,
		     Xhot, tauspu_s, tauspu_l, ThisTask, i, j);
		  fflush(stdout);
		}

	    }			// CLOSES for(j=1;j<LT_NMetP; j++)

/*
//if (taush/YEAR_IN_CODE_UNITS > 0 && taush/YEAR_IN_CODE_UNITS <1e9)
//{
printf(" n=%e Dl=%g Ds=%g Temp=%e XCold=%e\n", xn, Dl, Ds, T_wk, XCold); fflush(stdout);
printf("            dtime %g taush %g tauco %g tauac %g  tausp_s %g tausp_l %g\n", dtime/YEAR_IN_CODE_UNITS, taush/YEAR_IN_CODE_UNITS,
       tauco/YEAR_IN_CODE_UNITS, tauac/YEAR_IN_CODE_UNITS, tauspu_s/YEAR_IN_CODE_UNITS, tauspu_l/YEAR_IN_CODE_UNITS); fflush(stdout);
//}
*/
	  /* zeroing the sn count */
	  SphP[i].numSnIa = 0.;
	  SphP[i].numSnII = 0.;
	}			// CLOSES if(P[i].Type==0)
    }				// CLOSES for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])

  tend = second();
  All.CPU_Sum[CPU_COOLINGSFR] = timediff(tstart, tend);
  /* XXXX Maybe we want to add a CPU_DUST */

}				// CLOSES evolve_dust()


void compute_ism_properties_effective_model(double M_h, double M_c, double Vol_tot,
					    double T_h, double *fill_h, double *Rho_h, double *Rho_c,
					    double *n_h, double *n_c)
{
  double Vol_h, Vol_c, Frac_h, Frac_c;

  Frac_h = M_h / (M_c + M_h);
  Frac_c = 1.0 - Frac_h;
  *fill_h = 1.0 / (1.0 + Frac_c / Frac_h * mu_h / mu_c * All.TempClouds / T_h);
  Vol_h = Vol_tot * *fill_h;
  Vol_c = Vol_tot * (1.0 - *fill_h);
  *Rho_h = M_h / Vol_h;
  *n_h = *Rho_h * All.UnitDensity_in_cgs / mu_h / PROTONMASS;
  *Rho_c = M_c / Vol_c;
  *n_c = *Rho_c * All.UnitDensity_in_cgs / mu_c / PROTONMASS;

  return;
}


#endif // CLOSES defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
