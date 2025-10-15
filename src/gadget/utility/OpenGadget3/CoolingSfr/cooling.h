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
#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

// Number of steps in COOLing TABle
#define NCOOLTAB  2000

#define SMALLNUM 1.0e-60
#define COOLLIM  0.1
#define HEATLIM	 20.0

#define JAMPL	1.0		/* amplitude factor relative to input table */
#define TABLESIZE 200		/* Max # of lines in TREECOOL */

extern double XH;		/* hydrogen abundance by mass */
extern double yhelium;

//extern double Tmin;   /* in log10 */
//extern double Tmax;
//extern double deltaT;

//extern double *BetaH0, *BetaHep, *Betaff;
//extern double *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp;
//extern double *GammaeH0, *GammaeHe0, *GammaeHep;

//extern double J_UV, gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;

//extern double ne, necgs, nHcgs;
//extern double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
//extern double gJH0ne, gJHe0ne, gJHepne;
//extern double nH0, nHp, nHep, nHe0, nHepp;

extern float inlogz[TABLESIZE];
extern float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
extern float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
extern int nheattab;		/* length of table */

void cl_set_Tmin(double value);
void cl_set_Tmax(double value);
void cl_set_deltaT(double value);
void cl_set_ne(double value);
void cl_set_necgs(double value);
void cl_set_nHcgs(double value);
void cl_set_bH0(double value);
void cl_set_bHep(double value);
void cl_set_bff(double value);
void cl_set_aHp(double value);
void cl_set_aHep(double value);
void cl_set_aHepp(double value);
void cl_set_ad(double value);
void cl_set_geH0(double value);
void cl_set_geHe0(double value);
void cl_set_geHep(double value);
void cl_set_gJH0ne(double value);
void cl_set_gJHe0ne(double value);
void cl_set_gJHepne(double value);
void cl_set_nH0(double value);
void cl_set_nHp(double value);
void cl_set_nHep(double value);
void cl_set_nHe0(double value);
void cl_set_nHepp(double value);


// This struct makes cooling data globally available. Do not write to it outside of cooling.c
// note that "global" in this case means within each process - not all processes!
extern struct global_cooling_data
{
  // in log10
  double Tmin;
  double Tmax;
  double deltaT;

  // Recombination rates
  double *AlphaHp;
  double *AlphaHep;
  double *Alphad;		// likely dielectronic recombination of He
  double *AlphaHepp;

  double *BetaH0;
  double *BetaHep;
  double *Betaff;

  // Collisional ionization rates
  double *GammaeH0;
  double *GammaeHe0;
  double *GammaeHep;

#ifdef _WIN32
  double J_UV;
  double gJH0;
  double gJHep;
  double gJHe0;
  double epsH0;
  double epsHep;
  double epsHe0;
#else
  double J_UV = 0;
  double gJH0 = 0;
  double gJHep = 0;
  double gJHe0 = 0;
  double epsH0 = 0;
  double epsHep = 0;
  double epsHe0 = 0;
#endif

  double ne;
  double necgs;
  double nHcgs;

  double bH0;
  double bHep;
  double bff;
  double aHp;
  double aHep;
  double aHepp;
  double ad;
  double geH0;
  double geHe0;
  double geHep;

  double gJH0ne;
  double gJHe0ne;
  double gJHepne;

  double nH0;
  double nHp;
  double nHep;
  double nHe0;
  double nHepp;
} Cool;


typedef struct
{
  MyFloat u_old;
  MyFloat rho_cgs;
#if defined(METALS) || defined(LT_STELLAREVOLUTION)
  MyFloat Z;
#endif
#ifdef LT_STELLAREVOLUTION
  double Metals[LT_NMetP];
  double DZ, Redshift;
#ifdef GL_DUST_COOLING
  double DL, DS;
#endif
#endif
} coolingdata_in;

typedef struct
{
  double ne_guess;
  double LambdaNet;
  double Temperature;
} coolingdata_out;


typedef struct
{
  double logT;
  double rho;
  double *nelec;
#ifdef LT_METAL_COOLING
  double Z;
#endif
} CoolingRateArgs;		// TODO better name. Possible to merge with something else?



extern void particle2in_cooling(coolingdata_in * in, int i);
extern double cl_GetCoolingTimeAll(coolingdata_in in, coolingdata_out * out);
//extern double cl_GetCoolingRateFromUAll(coolingdata_in in, coolingdata_out *out);
extern double DoCoolingAll(coolingdata_in in, coolingdata_out * out);


double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer);
void find_abundances_and_rates(double logT, double rho, double *ne_guess);
void InitCool(void);
static void InitCoolMemory(void);
void IonizeParams(void);
void IonizeParamsFunction(void);
double INLINE_FUNC LogTemp(double u, double ne);
static void MakeCoolingTable(void);
static void ReadIonizeParams(char *fname);
void SetZeroIonization(void);
void TestCool(void);

#if !defined(LT_METAL_COOLING) && !defined(LT_METAL_COOLING_WAL)

double convert_u_to_tempSTD(double u, double rho, double *ne_guess);
//double CoolingRateSTD(double logT, double rho, double *nelec);
//double CoolingRate(CoolingRateArgs * args);
double CoolingRateFromU(double u, double rho, double *ne_guess);
double DoCoolingSTD(double u_old, double rho, double dt, double *ne_guess);
double GetCoolingTimeSTD(double u_old, double rho, double *ne_guess);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess);

#else

#if defined(LT_METAL_COOLING)
double convert_u_to_tempSTD(double u, double rho, double *ne_guess);
double GetCoolingTimeZ(double u_old, double rho, double *ne_guess, double Z, double *temp);
double CoolingRateFromU(double u, double rho, double *ne_guess, double Z, double *temp);
double DoCoolingZ(double u_old, double rho, double dt, double *ne_guess, double Z, double *temp);
//double CoolingRateZ(double logT, double rho, double *nelec, double Z);
//double CoolingRate(CoolingRateArgs * args);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess,
			    double Z, double *temp);
#endif

#if defined(LT_METAL_COOLING_WAL)
#ifdef GL_DUST_COOLING
double CoolingRateMET(double Temp, double Redshift, double DZ, double DL, double DS, double *Metallicities);
#else
double CoolingRateMET(double Temp, double Redshift, double DZ, double *Metallicities);
#endif // GL_DUST_COOLING
#ifdef GL_DUST_COOLING
double CoolingRateFromU(double U, double Redshift, double DZ, double DL, double DS, double *Metallicities, 
double *temp);
#else
double CoolingRateFromU(double U, double Redshift, double DZ, double *Metallicities, double *temp);
#endif // GL_DUST_COOLING
#ifdef GL_DUST_COOLING
double DoCoolingMET(double U_in, double Rho, double *Metallicities, double DL, double DS, double Redshift, 
        double DZ, double dt, double *temp);
#else
double DoCoolingMET(double U_in, double Rho, double *Metallicities, double Redshift, double DZ, double dt,
		    double *temp);
#endif // GL_DUST_COOLING
double GetUFromLambda(double Lambda_in, double Rho, double *Metallicities, double Redshift, double DZ,
		      double dt, double *temp, double u_min, double u_max);
double convert_u_to_tempMET(double U, double Redshift, double DZ, double *Metallicities);
#ifdef GL_DUST_COOLING
double GetCoolingTimeMET(double U, double Rho, double Redshift, double DZ, double DL, double DS, 
        double *Metallicities, double *temp);
#else
double GetCoolingTimeMET(double U, double Rho, double Redshift, double DZ, double *Metallicities,
			 double *temp);
#endif // GL_DUST_COOLING
int get_cool_redshift(double, double *);
double lt_cl_GetCoolRedshift(double);
double get_max_cool_redshift();
int get_cool_n_el();
int Is_a_Coolant(int);
void WalCool_set_PID(MyIDType);
void WalCool_tables_load(double);
void WalCool_get_collis_table();
void set_cooltable_index(int);
double *set_metallicities(int, double *, double);
#if defined(SUBFIND)
double *set_metallicities_subfind(int, double *, double);
double get_HImass_subfind(int i, double, double);
#endif
void read_cooling_tables_dummy();
void WalCool_Initialize();
#endif

#endif
