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

* Changes by NVIDIA Corp. (see change logs) are licensed as follows:
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*/


/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#ifdef _WIN32
#   include <windows.h>
#   define __restrict__
#   define __thread
#   define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef GADGET3_IO_LIB
#include <mpi.h>
#endif
#include <stdio.h>
#include <math.h>

#ifndef GADGET3_IO_LIB
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../gadgetconfig.h"
#include "../CodeBase/switches.h"
#include "../CodeBase/precision.h"

#ifdef FPE_TRAPPING
#define _GNU_SOURCE
#include <fenv.h>
#endif


#ifdef MPISENDRECV_CHECKSUM
#define MPI_Sendrecv MPI_Check_Sendrecv
#endif

#ifdef MPISENDRECV_SIZELIMIT
#define MPI_Sendrecv MPI_Sizelimited_Sendrecv
#endif

#include "../System/tags.h"
#if defined(CHEMISTRY) || defined (UM_CHEMISTRY)
#include "../Chemistry/chemistry.h"
#endif

#ifdef GADGET3_IO_LIB
#   include <assert.h>
#else
#   include "../System/assert.h"
#endif

#if GADGET_HYDRO == HYDRO_MFM
#include "../Mfm/Eos.hpp"
#include "../Mfm/FluidVector.hpp"
#include "../Mfm/FluxSolver.hpp"
#include "../Mfm/Matrix.hpp"
#include "../Mfm/Vector.hpp"
using namespace OpenGadget3;
extern Eos < NUMDIMS > *eos;
extern FluxSolver < NUMDIMS > *fluxSolver;
//extern SlopeLimiter<NUMDIMS> *slopeLimiter;
#endif


//see https://stackoverflow.com/questions/4744437/testing-if-float-value-is-nan
#ifndef ISNAN
#define ISNAN(x) (x != x)
#endif

#ifdef BLACK_HOLES
#define IS_SPH_ALIVE(Pi) ( ((Pi).Type == 0) &&  ((Pi).Mass > 0) )
#define IS_BLACKHOLE(Pi) ( (Pi).Type == 5 )
#else
#define IS_SPH_ALIVE(Pi) ((Pi).Type==0)
#define IS_BLACKHOLE(Pi) 0
#endif

#if defined(LT_STELLAREVOLUTION) || defined(SNIA_HEATING)
#if defined(LT_STELLAREVOLUTION)
#define IS_STAR_ACTIVE(Pp, n) ( ((Pp)[(n)].Type == 4) && (is_chemically_active((n))) )
#else
#define IS_STAR_ACTIVE(Pp, n) ((Pp)[(n)].Type == 4)
#endif
#else // neither LT_STELLAREVOLUTION nor SNIA_HEATING
#define IS_STAR_ACTIVE(Pp, n) 0
#endif

#ifndef VERBOSE_LEVEL
#define VERBOSE_LEVEL 0
#endif
#ifndef VERBOSE_DEFAULT_CONDITION
#define VERBOSE_DEFAULT_CONDITION ThisTask==0
#endif
#define VERBOSE_OUTPUT(f_, ...) do{  printf((f_), ##__VA_ARGS__); printf("\n"); fflush(stdout); }while(0)
#define VERBOSE_BELOW_LEVEL(level, f_, ...)  if(level<=VERBOSE_LEVEL && VERBOSE_DEFAULT_CONDITION)  VERBOSE_OUTPUT(f_,  ##__VA_ARGS__)
#define VERBOSE(level, f_, ...) VERBOSE_BELOW_LEVEL(level, f_,  ##__VA_ARGS__)
#define VERBOSE_IF(condition,  f_, ...) if(condition) VERBOSE_OUTPUT(f_, ##__VA_ARGS__)
#define VERBOSE_CONDITION(condition, level, f_, ...) if(condition) VERBOSE(level, f_, ##__VA_ARGS__)
#define VERBOSE_BELOW_LEVEL_ALL(level, f_, ...)  if(level<=VERBOSE_LEVEL)  VERBOSE_OUTPUT(f_,  ##__VA_ARGS__)
#define VERBOSE_ALL(level, f_, ...) VERBOSE_BELOW_LEVEL_ALL(level, f_,  ##__VA_ARGS__)
#define VERBOSE_CONDITION_ALL(condition, level, f_, ...) if(condition) VERBOSE_ALL(level, f_, ##__VA_ARGS__)
#define PANIC_W_OUTPUT_CONDITION(output_condition, f_, ...) EndRunString(output_condition, __FUNCTION__, __FILE__,__LINE__, (f_), ##__VA_ARGS__)
#define PANIC(f_, ...) PANIC_W_OUTPUT_CONDITION(ThisTask==0, f_,  ##__VA_ARGS__)
#define PANIC_IF(condition, f_, ...) if(condition) PANIC_W_OUTPUT_CONDITION(1, f_, ##__VA_ARGS__)
#define PANIC_ANY(f_, ...) PANIC_W_OUTPUT_CONDITION(1, f_, ##__VA_ARGS__)
// see https://gcc.gnu.org/onlinedocs/cpp/Swallowing-the-Semicolon.html
#ifdef _OPENMP
#define THREAD_ID_OR_ZERO() omp_get_thread_num()
#else
#define THREAD_ID_OR_ZERO() 0
#endif



#define OMP_SORT_THRESH 1000

#ifdef MYSORT
#ifdef _OPENMP
#define MYSORT_DATAINDEX omp_mysort
#else
#define MYSORT_DATAINDEX mysort_dataindex
#endif
#else // MYSORT
#ifdef _OPENMP
#define MYSORT_DATAINDEX omp_qsort
#else
#define MYSORT_DATAINDEX qsort
#endif
#endif

#ifdef _OPENMP
#define STD_SORT omp_qsort
#else
#define STD_SORT qsort
#endif

// compiler specific data alignment hints
// XLC compiler
#if defined(__xlC__)
#define ALIGN(n) __attribute__((__aligned__(n)))
// GNU compiler
#elif defined(__GNUC__)
#define ALIGN(n) __attribute__((__aligned__(n)))
// Intel Compiler
#elif defined(__INTEL_COMPILER)
// GNU Intel Compiler
#define ALIGN(n) __declspec(align(n))
// Unknown Compiler
#else
#define ALIGN(n)
#endif



#ifdef PERIODIC
#ifdef POWER6
#define NEAREST_X(x) ((x)  - boxSize_X * __frin ( (x) * inverse_boxSize_X))
#define NEAREST_Y(y) ((y)  - boxSize_Y * __frin ( (y) * inverse_boxSize_Y))
#define NEAREST_Z(z) ((z)  - boxSize_Z * __frin ( (z) * inverse_boxSize_Z))
#define NEAREST(x) ((x)  - boxSize * __frin ( (x) * inverse_boxSize))
#else
#define NEAREST_X(x) (((x)>boxHalf_X)?((x)-boxSize_X):(((x)<-boxHalf_X)?((x)+boxSize_X):(x)))
#define NEAREST_Y(y) (((y)>boxHalf_Y)?((y)-boxSize_Y):(((y)<-boxHalf_Y)?((y)+boxSize_Y):(y)))
#define NEAREST_Z(z) (((z)>boxHalf_Z)?((z)-boxSize_Z):(((z)<-boxHalf_Z)?((z)+boxSize_Z):(z)))
#define NEAREST(x) (((x)>boxHalf)?((x)-boxSize):(((x)<-boxHalf)?((x)+boxSize):(x)))
#define __fsel(crit,age,alt) (((crit) >= 0.0) ? (age) : (alt))
#endif
#else
#define NEAREST_X(x) (x)
#define NEAREST_Y(x) (x)
#define NEAREST_Z(x) (x)
#define NEAREST(x) (x)
#endif

#define ASSIGN_ADD(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
#define ASSIGN_MIN(x,y,mode) (mode == 0 ? (x=y) : (x=std::min(x,y)))
#define ASSIGN_MAX(x,y,mode) (mode == 0 ? (x=y) : (x=std::max(x,y)))

#define  GADGETVERSION   "0.1"	/*!< code version string */

#ifndef  GENERATIONS
#define  GENERATIONS     2	/*!< Number of star particles that may be created per gas particle */
#endif

#ifdef   ENLARGE_DYNAMIC_RANGE_IN_TIME
typedef long long integertime;
#define  TIMEBINS        60
#define  TIMEBASE        (((long long)1)<<TIMEBINS)	/*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
							 *   where TIMESPAN needs to be a power of 2.
							 */
#else

#ifdef GADGET3_IO_LIB
#define integertime int
#else
typedef int integertime;
#endif

#define  TIMEBINS        29
#define  TIMEBASE        (1<<TIMEBINS)	/*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
					 *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
					 *   to 2^29
					 */
#endif


#ifndef  MULTIPLEDOMAINS
#define  MULTIPLEDOMAINS     1
#endif

#ifdef ONEDIM
#define DIMS 1
#else
#ifdef TWODIMS			/* will only be compiled in 2D case */
#define DIMS 2
#else
#define DIMS 3
#endif
#endif

#ifndef  TOPNODEFACTOR
#define  TOPNODEFACTOR       2.5
#endif

#ifndef  GRAVCOSTLEVELS
#define  GRAVCOSTLEVELS      6
#endif

#define  NUMBER_OF_MEASUREMENTS_TO_RECORD  6	/* this is the number of past executions of a timebin that the reported average CPU-times average over */


#if !defined(fSIDM) && !defined(rSIDM)
#define NODELISTLENGTH       8
#else
#define  NODELISTLENGTH      NODELISTLENGTH_default	//8
#endif


typedef unsigned long long peanokey;


#define  BITS_PER_DIMENSION 21	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))

#define  BITS_PER_DIMENSION_SAVE_KEYS 10
#define  PEANOCELLS_SAVE_KEYS (((peanokey)1)<<(3*BITS_PER_DIMENSION_SAVE_KEYS))

#if defined(fSIDM) || defined(rSIDM)
#include "../mSIDM/msidm.hpp"
#endif

#ifdef LT_STELLAREVOLUTION
#include "../CoolingSfr/Sfr_LT/lt.h"
#endif

#ifdef LT_ADD_GAL_TO_SUB
extern float *tempiAS, *CB07, *Filters_Effective_Wavelenght;
#ifdef OBSERVER_FRAME
extern float *CB07obs;
#ifdef INTERP_OBSERVER_FRAME
extern float CB07obs_redshifts[INTERP_OBSERVER_FRAME];
extern int flag_readCB07obs_redshifts;
#endif
#endif
#endif

#ifdef GM_MUPPI
#define ConeAngle_FB ((MyAtLeastDouble)1.047195)	/* NB: 60 degrees of cone semi-aperture -> 2.094395 rad */
					  /* NB: 45 degrees of cone semi-aperture -> 1.570796 rad */
					  /* NB: 30 degrees of cone semi-aperture -> 1.047195 rad */
					  /* NB: 15 degrees of cone semi-aperture -> 0.523360 rad */

#endif


#define  check_particles()          check_particles_info( __FUNCTION__, __FILE__, __LINE__)
#define  check_node()               check_node_info( __FUNCTION__, __FILE__, __LINE__)

#ifndef DISABLE_MEMORY_MANAGER
#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)
#else

#if (DISABLE_MEMORY_MANAGER==2)
#define mymalloc(x,y) calloc(y,1)
#define mymalloc_movable(x,y,z) calloc(z,1)
#else /* DISABLE_MEMORY_MANAGER==1 */
#define mymalloc(x,y) malloc(y)
#define mymalloc_movable(x,y,z) malloc(z)
#endif /* DISABLE_MEMORY_MANAGER==1 */

#define  myrealloc(x, y)           realloc(x, y)
#define  myrealloc_movable(x, y)   realloc(x, y)

#define  myfree(x)                 free(x)
#define  myfree_movable(x)         free(x)

#define  report_memory_usage(x, y) printf("Memory manager disabled.\n")
#endif

#ifdef LT_STELLAREVOLUTION
#define endrun(x) EndRun(x, __FUNCTION__, __FILE__, __LINE__)
#endif

#ifndef  GAMMA
#define  GAMMA         (5.0/3.0)	/*!< adiabatic index of simulated gas */
#endif

#define  GAMMA_INV     (1./GAMMA)
#define  GAMMA_MINUS1  (GAMMA-1)
#define  GAMMA_MINUS1_INV  (1./(GAMMA-1))

#ifndef  HYDROGEN_ONLY
#define  HYDROGEN_MASSFRAC 0.76	/*!< mass fraction of hydrogen, relevant only for radiative cooling */
#else
#define  HYDROGEN_MASSFRAC 1	/*!< mass fraction of hydrogen, relevant only for radiative cooling */
#endif

#define  METAL_YIELD       0.02	/*!< effective metal yield for star formation */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  RNDTABLE 8192

/* ... often used physical constants (cgs units) */

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.38066e-16
#define  GAS_CONST   8.31425e7
#define  C_LIGHT     2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#define  LYMAN_ALPHA      1215.6e-8	/* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8	/* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615
#define  ELECTRONVOLT_IN_ERGS      1.60217733e-12


#ifdef MAGNETIC
#ifdef SFR
#define POW_CC 1./3.
#endif

#ifdef MU0_UNITY
#define MU0 1.0
#define MU0_1 1.0
#else
#define MU0  12.566370614	/* 4pi */
#define MU0_1 0.079577472	/* 1/4pi */
#endif
#endif

#ifdef NAVIERSTOKES
#define  LOG_LAMBDA      37.8	/* logarithmic Coulomb factor */
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY) || defined(INCLUDE_RADIATION) || defined(KSPACE_NEUTRINOS_2)
#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */
#endif

#if defined(KSPACE_NEUTRINOS_2) || defined(INCLUDE_RADIATION)
/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 *The correction is (3.046/3)^(1/4), for N_eff = 3.046 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00381)	/* Neutrino + antineutrino background temperature in Kelvin */
#endif

/*With slightly relativistic massive neutrinos, for consistency we need to include radiation.
 * A note on normalisation (as of 08/02/2012):
 * CAMB appears to set Omega_Lambda + Omega_Matter+Omega_K = 1,
 * calculating Omega_K in the code and specifying Omega_Lambda and Omega_Matter in the paramfile.
 * This means that Omega_tot = 1+ Omega_r + Omega_g, effectively
 * making h0 (very) slightly larger than specified.
 */
#ifdef INCLUDE_RADIATION
/*Stefan-Boltzmann constant in cgs units*/
#define STEFAN_BOLTZMANN 5.670373e-5
/* Omega_g = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H^2)*/
#define OMEGAG (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*C*C*C*HUBBLE*HUBBLE)*pow(T_CMB0,4)/All.HubbleParam/All.HubbleParam)
#if (defined KSPACE_NEUTRINOS_2) || (defined NEUTRINOS)
    /*Neutrinos are included elsewhere */
#define OMEGAR OMEGAG
#else
    /*Neutrinos are included in the radiation */
    /*For massless neutrinos, rho_nu/rho_g = 7/8 (T_nu/T_cmb)^4 *N_eff, but we absorbed N_eff into T_nu above */
#define OMEGANU (OMEGAG*7/8.*pow(TNU/T_CMB0,4)*3)
    /*With massless neutrinos only, add the neutrinos to the radiation */
#define OMEGAR (OMEGAG+OMEGANU)
#endif
#else
	/*Default is no radiation */
#define OMEGAR 0.
#endif

/* For convenience define OMEGAK*/
#define OMEGAK (1-All.Omega0 - All.OmegaLambda)

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

/*Determines the maximum size of arrays related to the number of CR populations */
#ifndef NUMCRPOP		/*!< Number of CR populations pressent in parameter file */
#define NUMCRPOP 1
#endif



#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif


/* some flags for the field "flag_ic_info" in the file header */
#define FLAG_ZELDOVICH_ICS     1
#define FLAG_SECOND_ORDER_ICS  2
#define FLAG_EVOLVED_ZELDOVICH 3
#define FLAG_EVOLVED_2LPT      4
#define FLAG_NORMALICS_2LPT    5


#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25
#endif
#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5
#endif

#define COND_TIMESTEP_PARAMETER 0.25
#define VISC_TIMESTEP_PARAMETER 0.25

#define MAXLEN_OUTPUTLIST 14000	/*!< maxmimum number of entries in output list */

enum
{
  /*!< length of the lookup table used to hold the drift and kick factors */
  DRIFT_TABLE_LENGTH = 1000
};

#define MAXITER 150

#ifndef LINKLENGTH
#define LINKLENGTH 0.2
#endif

#ifndef FOF_GROUP_MIN_LEN
#define FOF_GROUP_MIN_LEN 32
#endif

#define MINRESTFAC 0.05

#ifdef SUBFIND_DENSITY_AND_POTENTIAL	/*!< activate needed options */
#define ONLY_PRODUCE_HSML_FILES
#define COMPUTE_POTENTIAL_ENERGY
#define SUBFIND_RESHUFFLE_AND_POTENTIAL
#define SUBFIND_RESHUFFLE_CATALOGUE
#endif
#ifdef ADAPTGRAVSOFT
#ifndef AGS_MAXSOFT
/*!< AGS_MAXSOFT gives the maximum allowed gravitational softening when using the TreePM method.
 *  The quantity is given in units of the scale used for the force split (ASMTH) */
#define AGS_MAXSOFT 0.5
#endif
#ifndef AGS_MINSOFT
/*!< AGS_MINSOFT gives the minimum allowed gravitational softening when using the TreePM method.
 *  The quantity is given in units of the scale used for the force split (ASMTH) */
#define AGS_MINSOFT 0.001
#endif
#endif

#ifndef GDE_TYPES
#define GDE_TYPES 2
#endif


struct unbind_data
{
  int index;
};


#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif

#define FLT(x) (x)

#ifdef GDE_BIGFLOAT
#include "bigfloat.h"
typedef BigFloat MyBigFloat;
extern int xErrNo;
#else
typedef MyDouble MyBigFloat;
#endif

#ifdef GDE_BIGFLOAT
#define GDE_ABS(x) ((x).abs())
#define GDE_SQRT(x) ((x).sqrt())
#define GDE_LOG(x) ((x).log())
#else
#define GDE_ABS(x) (fabs(x))
#define GDE_SQRT(x) (sqrt(x))
#define GDE_LOG(x) (log(x))
#endif



enum cpufields
{

  CPU_ALL,
  CPU_TREEWALK1,
  CPU_TREEWALK2,
  CPU_TREEWAIT1,
  CPU_TREEWAIT2,
  CPU_TREESEND,
  CPU_TREERECV,
  CPU_TREEMISC,
  CPU_TREEBUILD,
  CPU_TREEUPDATE,
  CPU_TREEHMAXUPDATE,
  CPU_DOMAIN,
  CPU_DENSCOMPUTE,
  CPU_DENSWAIT,
  CPU_DENSCOMM,
  CPU_DENSMISC,
  CPU_HYDCOMPUTE,
  CPU_HYDWAIT,
  CPU_HYDCOMM,
  CPU_HYDMISC,
  CPU_DRIFT,
  CPU_TIMELINE,
  CPU_POTENTIAL,
  CPU_MESH,
  CPU_PEANO,
  CPU_COOLINGSFR,
  CPU_SNAPSHOT,
  CPU_FOF,
  CPU_BLACKHOLES,
  CPU_MISC,
#ifdef GM_STARDENSITY
  CPU_STARDENSMISC,
  CPU_STARCOMPUTE,
  CPU_STARWAIT,
  CPU_STARCOMM,
#endif
  CPU_SMTHCOMPUTE,
  CPU_SMTHWAIT,
  CPU_SMTHCOMM,
  CPU_SMTHMISC,
  CPU_HOTNGBS,
  CPU_WEIGHTS_HOT,
  CPU_ENRICH_HOT,
  CPU_WEIGHTS_COLD,
  CPU_ENRICH_COLD,
  CPU_CSMISC,
  CPU_HYDNETWORK,
  CPU_AGSDENSCOMPUTE,
  CPU_AGSDENSWAIT,
  CPU_AGSDENSCOMM,
  CPU_AGSDENSMISC,
  CPU_AGSTREEHMAXUPD,
  CPU_MG_CIC,
  CPU_MG_FIELDSOLVE,
  CPU_MG_EFF_MASS,
  CPU_CONDUCTION,
#if GADGET_HYDRO == HYDRO_MFM
  CPU_MFMGRADCOMPUTE,
  CPU_MFMGRADWAIT,
  CPU_MFMGRADCOMM,
  CPU_MFMGRADMISC,
  CPU_MFMLIMCOMPUTE,
  CPU_MFMLIMWAIT,
  CPU_MFMLIMCOMM,
  CPU_MFMLIMMISC,
  CPU_MFMFLUXCOMPUTE,
  CPU_MFMFLUXWAIT,
  CPU_MFMFLUXCOMM,
  CPU_MFMFLUXMISC,
#endif

  CPU_MUPPI,
  CPU_MUPPINET,

#if defined(fSIDM) || defined(rSIDM)
  CPU_MSIDMCOMPUTE,
  CPU_MSIDMPREPP,
  CPU_MSIDMCOMM,
  CPU_MSIDMWAIT,
  CPU_MSIDMNETWORK,
  CPU_MSIDMMISC,
#endif

// this gives the number of timers that have to be added to get the total time.
// must be just after last of them.
  CPU_PARTS4SUM,
  CPU_PARTS4SUMM1 = CPU_PARTS4SUM - 1,	// ensure counting continues correctly
// elements below this line are excluded from the sum

#ifdef SUBTIMERS
  CPU_SUBFIND_TREEBUILD_SPECIES,
  CPU_SUBFIND_TREEBUILD,
  CPU_SUBFIND_SMOOTHINGLENGTH,
  CPU_SUBFIND_DENSITY,
  CPU_SUBFIND_DMDENSITY,
  CPU_SUBFIND_SAVE_DENSITY,
  CPU_SUBFIND_EXCHANGE,
  CPU_SUBFIND_COLLHALOS,
  CPU_SUBFIND_SORTLOCAL,
  CPU_SUBFIND_LOCALGROUPS,
  CPU_SUBFIND_UNSORTLOCAL,
  CPU_SUBFIND_EXCHANGERETURN,
  CPU_SUBFIND_DOMAIN,
  CPU_SUBFIND_MASSES,
  CPU_SUBFIND_CONT,
  CPU_SUBFIND_SAVEFINAL,
  CPU_SUBFIND_TOT,
  CPU_FOF_TOT,
#endif

#if defined(fSIDM) || defined(rSIDM)	/* these should be excluded in CPU_PARTS4SUM */
  CPU_FSIDMSCATTER,
  CPU_FSIDMSCATTER_DF,
  CPU_FSIDMSCATTER_RC,
  CPU_RSIDMSCATTER,
  CPU_MSIDMPARALLEL,
  CPU_MSIDMLOCAL,
#endif

  CPU_PARTS			/* this gives the number of parts above (must be last) */
};

#define CPU_STRING_LEN 120

#if !defined(QUINTIC_KERNEL)
#if !defined(TWODIMS) && !defined(ONEDIM)
//#define  NUMDIMS 3            /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786	/*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#ifndef  ONEDIM
//#define  NUMDIMS 2            /*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI	/*!< Coefficient for kernel normalization. */
#else
//#define  NUMDIMS 1             /*!< For 1D-normalized kernel */
#define  KERNEL_COEFF_1  (4.0/3)
#define  KERNEL_COEFF_2  (8.0)
#define  KERNEL_COEFF_3  (24.0)
#define  KERNEL_COEFF_4  (16.0)
#define  KERNEL_COEFF_5  (8.0/3)
#define  KERNEL_COEFF_6  (-8.0)
#define  NORM_COEFF      2.0
#endif
#endif
#else /* here comes the QUINTIC kernel */
#if !defined(TWODIMS) && !defined(ONEDIM)
//#define  NUMDIMS 3
#define  NORM_COEFF      4.188790204786
#else
#ifndef  ONEDIM
//#define  NUMDIMS 2
#define  NORM_COEFF      M_PI
#else
//#define  NUMDIMS 1
#define  NORM_COEFF      2.0
#endif
#endif
#endif /* end of !QUINTIC */

#define SPP SphP

#ifdef PERIODIC
extern MyDouble boxSize, boxHalf, inverse_boxSize;


#ifdef LONG_X
extern MyDouble boxSize_X, boxHalf_X, inverse_boxSize_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#define inverse_boxSize_X inverse_boxSize
#endif
#ifdef LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y, inverse_boxSize_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#define inverse_boxSize_Y inverse_boxSize
#endif
#ifdef LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z, inverse_boxSize_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#define inverse_boxSize_Z inverse_boxSize
#endif
#endif

#ifdef PERIODIC
#ifndef POWER6
#define NGB_PERIODIC_LONG(x,box,hbox) ((fabs(x)>hbox)?(box-fabs(x)):fabs(x))
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp)
#else
#ifdef DOUBLEPRECISION
#define NGB_PERIODIC_LONG(x,box,hbox) __fsel(hbox-fabs(x),fabs(x),box-fabs(x))
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),__fsel(boxHalf_X-xtmp,xtmp,boxSize_X-xtmp))
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),__fsel(boxHalf_Y-xtmp,xtmp,boxSize_Y-xtmp))
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),__fsel(boxHalf_Z-xtmp,xtmp,boxSize_Z-xtmp))
#else
#define NGB_PERIODIC_LONG(x,box,hbox) __fsel(hbox-fabsf(x),fabsf(x),box-fabsf(x))
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabsf(x),__fsels(boxHalf_X-xtmp,xtmp,boxSize_X-xtmp))
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabsf(x),__fsels(boxHalf_Y-xtmp,xtmp,boxSize_Y-xtmp))
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabsf(x),__fsels(boxHalf_Z-xtmp,xtmp,boxSize_Z-xtmp))
#endif
#endif
#else
#define NGB_PERIODIC_LONG(x,box,hbox) fabs(x)
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#endif

#define FACT1 0.366025403785	/* FACT1 = 0.5 * (sqrt(3)-1) */
#define FACT2 0.86602540	/* FACT2 = 0.5 * sqrt(3) */



/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

#ifdef UM_SEN_H2_SOLOMON_PROCESS
#define JLWn 23
extern float JLW[JLWn];
extern float z_JLW[JLWn];
#endif

#ifdef _OPENMP
extern int NActivePart;
extern int *ActiveParticleList;
#if defined(LT_STELLAREVOLUTION)
extern int NActivePartMerk;
#endif
#endif

extern int FirstActiveParticle;
extern int *NextActiveParticle;
extern unsigned char *ProcessedFlag;

extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountSph[TIMEBINS];
extern int TimeBinActive[TIMEBINS];

extern int FirstInTimeBin[TIMEBINS];
extern int LastInTimeBin[TIMEBINS];
extern int *NextInTimeBin;
extern int *PrevInTimeBin;

#ifdef SFR
extern double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
extern double TimeBin_BH_mass[TIMEBINS];
extern double TimeBin_BH_dynamicalmass[TIMEBINS];
extern double TimeBin_BH_Mdot[TIMEBINS];
extern double TimeBin_BH_Medd[TIMEBINS];
#endif

#ifdef KD_MONITOR_PERFORMANCE
extern int count_performance;
#endif

extern __thread int ThisThread;	/* the Number of the local thread */
extern int NThread;		/* number of threads */
extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern int PTask;		/*!< note: NTask = 2^PTask */

#if !defined(GADGET3_IO_LIB) || defined(GADGET3_IO_LIB_MPI)
extern MPI_Comm MYMPI_COMM_WORLD;
#endif

#ifdef INVARIANCETEST
extern int World_ThisTask;
extern int World_NTask;
extern int Color;
extern MPI_Comm MPI_CommLocal;
#ifndef DO_NOT_REDEFINE_MPI_COMM_WORLD
#undef  MPI_COMM_WORLD
#define MPI_COMM_WORLD MPI_CommLocal
#endif
#endif


extern double CPUThisRun;	/*!< Sums CPU time of current process */

#ifdef RECOMPOSE_DOMAIN_FACTOR
extern int LastDomainUpdate;
#endif

extern int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;

extern int NumSphUpdate;	/*!< number of active SPH particles on local processor in current timestep  */

extern int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
extern int RestartSnapNum;

extern double PostProcessDt;	/*!< For the post-processing FKP the timestep is passed as a command line argument */

extern int SelRnd;

extern int TakeLevel;

extern int *Exportflag;		/*!< Buffer used for flagging whether a particle needs to be exported to another process */
extern int *Exportnodecount;
extern int *Exportindex;

extern int *Send_offset, *Send_count, *Recv_count, *Recv_offset;

#ifdef IMPORT_ALLtoALLv_FROM_G4
extern MPI_Win win;
#endif

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;

extern double CPU_Step[CPU_PARTS];
extern char CPU_Symbol[CPU_PARTS];
extern char CPU_SymbolImbalance[CPU_PARTS];
extern char CPU_String[CPU_STRING_LEN + 1];

extern double WallclockTime;	/*!< This holds the last wallclock time measurement for timings measurements */

extern int Flag_FullStep;	/*!< Flag used to signal that the current step involves all particles */

extern size_t HighMark_run, HighMark_domain, HighMark_gravtree, HighMark_pmperiodic,
  HighMark_pmnonperiodic, HighMark_sphdensity, HighMark_sphhydro, HighMark_addSPH;


#ifdef GM_MUPPI
extern size_t HighMark_mppin;
#endif


#ifdef ADAPTGRAVSOFT
extern size_t HighMark_agsdensity;
#endif

#if defined(fSIDM) || defined(rSIDM)
extern size_t HighMark_msidm;
#endif

extern int TreeReconstructFlag;
extern int GlobFlag;

extern char DumpFlag;


extern int NumPart;		/*!< number of particles on the LOCAL processor */

extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
extern int N_stars;		/*!< number of star particles in the LOCAL processor */
extern int N_BHs;

extern long long Ntype[6];	/*!< total number of particles of each type */
extern int NtypeLocal[6];	/*!< local number of particles of each type */

#ifndef GADGET3_IO_LIB
extern gsl_rng *random_generator;	/*!< the random number generator used */
#endif


#ifdef SFR
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
#endif

#define d_fftw_real fftw_real

extern double TimeOfLastTreeConstruction;	/*!< holds what it says */

extern int *Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */
extern int *RecyclingBuffer;
extern int HasFisrtParticleBeenProcessed;

extern double *R2ngblist;

extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int *DomainStartList, *DomainEndList;

extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern int *DomainList, DomainNumChanged;

#ifndef GADGET3_IO_LIB
extern peanokey *__restrict__ Key, *__restrict__ KeySorted;
#endif


extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  MyFloat GravCost;
  int Daughter;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;

extern int NTopnodes, NTopleaves;

extern double RndTable[RNDTABLE];


#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
#endif


#ifdef WRITE_KEY_FILES
extern peanokey *KeyIndex;
extern int *NPartPerKey, *PartKeyOffset;
extern int NKeys[6];
extern long long NKeysTotal[6];
#endif


/* variables for input/output , usually only used on process 0 */


extern char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,		/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU,			/*!< file handle for cpu.txt log-file. */
 *FdTimebin, *FdParamChangeLog;

#ifdef SCFPOTENTIAL
extern FILE *FdSCF;
#endif

#ifdef SFR
extern FILE *FdSfr;		/*!< file handle for sfr.txt log-file. */
#endif

#ifdef GM_MUPPI
extern FILE *FdExit;
#ifdef GM_COUNT_PARTICLES_IN_CONE
extern FILE *FdCone;
#endif
#ifdef MV_GM_STELLAR_KIN_FB2_OUTPUT
extern FILE *FdWind;
#endif
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
extern FILE *FdMolFrac;
#endif
#ifdef MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z_OUTPUT
extern FILE *FdHypernovae;
#endif
#ifdef MV_GM_AGNMUPPI
extern FILE *FdAGNmuppi;
#endif
#ifdef MV_GM_AGNMUPPI_OUTPUT
extern FILE *FdEgyAGN;
extern FILE *FdEgyAGNtot;
#ifdef MV_AGNMUPPI_COOLING_OFF
extern FILE *FdCoolOffAGNmuppi;
#endif
#endif
#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
extern FILE *FdCovFacAGN;
#endif
#endif

#ifdef DISTORTIONTENSORPS
#ifdef PMGRID
extern FILE *FdTidaltensor;	/*!< file handle for tidaltensor.txt log-file. */
#endif
#endif

#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;	/*!< file handle for blackholes.txt log-file. */
extern FILE *FdBlackHolesDetails;
#endif

#ifdef LMB_SPECTRAL_CRs
extern FILE *FdCRs;		/*!< file handle for CRs.txt log-file. */
extern FILE *FdCRsDetailsAdiabatic;
extern FILE *FdCRsDetailsRadiative;
extern FILE *FdCRsDetailsDpp;
extern FILE *FdCRsDetailsShockInjection;
extern FILE *FdCRsDetailsSfrInjection;
#endif



#ifdef FORCETEST
extern FILE *FdForceTest;	/*!< file handle for forcetest.txt log-file. */
#endif

#ifdef DARKENERGY
extern FILE *FdDE;		/*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
extern FILE *FdXXL;		/*!< file handle for xxl.txt log-file. */

#ifdef MAGNETIC
extern double MeanB;

#ifdef TRACEDIVB
extern double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
extern double MeanAlpha;
#endif
#endif

/*! table for the cosmological drift factors */
extern double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
extern double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];

#ifdef MAGNETIC
/*! table for the cosmological kick factor for induction equation */
extern double MagKickTable[DRIFT_TABLE_LENGTH];
#endif

#if defined(fSIDM) || defined(rSIDM)
/*! table for the cosmological kick factor for mSIDM */
extern double MSIDMKickTable[DRIFT_TABLE_LENGTH];
#endif


extern void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */



extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

#ifdef LT_STELLAREVOLUTION
  long long TotN_stars;		/*!<  total star particle number (global value) */
#endif

#ifdef NEUTRINOS
  long long TotNumNeutrinos;
  int tree_nu_on;
  double Time_tree_on_nu;
#endif

#ifdef BLACK_HOLES
  int TotBHs;
  int MaxPartBH;
#endif

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				   processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				   processor. */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int DoDynamicUpdate;

  int NumFilesPerSnapshot;	/*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/*!< maximum number of files that may be written simultaneously when
					   writing/reading restart-files, or when writing snapshot files */

  double BufferSize;		/*!< size of communication buffer in MB */
  int BunchSize;		/*!< number of particles fitting into the buffer in the parallel tree algorithm  */


  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  int DesNumNgb;		/*!< Desired number of SPH neighbours */
#ifdef SUBFIND
  int DesLinkNgb;
  double ErrTolThetaSubfind;
#endif

  double MaxNumNgbDeviation;	/*!< Maximum allowed deviation neighbour number */
#ifdef START_WITH_EXTRA_NGBDEV
  double MaxNumNgbDeviationStart;	/*!< Maximum allowed deviation neighbour number to start with */
#endif

  double ArtBulkViscConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;		/*!< the minimum allowed temperature expressed as energy per unit mass */
#ifdef ARTIFICIAL_CONDUCTIVITY
  double ArtCondConstant;
#ifdef TIME_DEP_ART_COND
  double ArtCondMin;
#endif
#endif

#ifdef KSPACE_NEUTRINOS
  double OmegaNu;
  int KspaceNeutrinoSeed;
  int Nsample;
  int SphereMode;
  char KspaceDirWithTransferfunctions[500];
  char KspaceBaseNameTransferfunctions[500];
  double PrimordialIndex;
  double Sigma8;
  double InputSpectrum_UnitLength_in_cm;
#endif

#if defined  KSPACE_NEUTRINOS_2
  double OmegaNu;
  char KspaceTransferFunction[500];
  double TimeTransfer;
  double OmegaBaryonCAMB;
  double InputSpectrum_UnitLength_in_cm;
  /*Allow for three massive neutrino species:
   * Could be made configurable at some point
   * Neutrino masses are in eV*/
#define NUSPECIES 3
  double MNu[NUSPECIES];
#endif


  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;	/*!< count particle updates since last domain decomposition */

  /* some variable for dynamic work-load adjustment based on CPU measurements */

  double cf_atime, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a;	/* various cosmological factors that are only a function of the current scale factor, and in Newtonian runs are set to 1 */

  /* system of units  */

  double UnitTime_in_s,		/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,	/*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,		/*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/*!< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,	/*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,		/*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,	/*!< factor to convert internal time to megayears/h */
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION
    UnitSurfDens_in_Msunpc2,	/*!< factor to convert internal surface density to Msun/pc^2 units */
#endif
    GravityConstantInternal,	/*!< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;				/*!< Gravity-constant in internal units */
  double UnitDensity_in_Gev_per_cm3;	/*!< factor to convert internal density unit to GeV/c^2 / cm^3 */
  /* Cosmology */

  double Hubble;		/*!< Hubble-constant in internal units */
  double Omega0,		/*!< matter density in units of the critical density (at z=0) */
    OmegaLambda,		/*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,		/*!< baryon density in units of the critical density (at z=0) */
    HubbleParam;		/*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */

  double BoxSize;		/*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;	/*!< flags that periodic boundaries are enabled */
  int ResubmitOn;		/*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;	/*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */
  int TypeOfTimestepCriterion;	/*!< gives type of timestep criterion (only 0 supported right now - unlike
				   gadget-1.1) */
  int OutputListOn;		/*!< flags that output times are listed in a specified file */
  int CoolingOn;		/*!< flags that cooling is enabled */
  int StarformationOn;		/*!< flags that star formation is enabled */

  int HighestActiveTimeBin;
  int HighestOccupiedTimeBin;

  /* parameters determining output frequency */

  int SnapshotFileCount;	/*!< number of snapshot that is written next */
  double TimeBetSnapshot,	/*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,	/*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,	/*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,	/*!< cpu-time when last restart-file was written */
    TimeBetStatistics,		/*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;		/*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;		/*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,			/*!< current time of the simulation */
    TimeBegin,			/*!< time of initial conditions of the simulation */
    TimeStep,			/*!< difference between current times of previous and current timestep */
    TimeMax;			/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval;	/*!< factor to convert from floating point time interval to integer timeline */
  integertime Ti_Current;	/*!< current time on integer timeline */
  integertime Previous_Ti_Current;
  integertime Ti_nextoutput;	/*!< next output time on integer timeline */
  integertime Ti_lastoutput;

#if defined(OUTPUT_LONGRANGE_POTENTIAL)
  int dumped_after_snap;	/* for whether the potential has been dumped after the snapshot! */
  integertime Ti_nextpotdump;	/*!< next potential dump time on integer timeline */
  int PotentialFileCount;	/*!< number of potential that is written next */
#endif

#ifdef PMGRID
  integertime PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[3][3], UpperCorner[3][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[3];
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
  double Epsilon;
#endif

  integertime Ti_nextlineofsight;

  int CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];

  int LevelToTimeBin[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_PARTS];	/*!< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,	/*!< minimum allowed timestep. Normally, the simulation terminates if the
				   timestep determined by the timestep criteria falls below this limit. */
    MaxSizeTimestep;		/*!< maximum allowed timestep */

  double MaxRMSDisplacementFac;	/*!< this determines a global timestep criterion for cosmological simulations
				   in comoving coordinates.  To this end, the code computes the rms velocity
				   of all particles, and limits the timestep such that the rms displacement
				   is a fraction of the mean particle separation (determined from the
				   particle mass and the cosmological parameters). This parameter specifies
				   this fraction. */

  int MaxMemSize;

  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency;	/*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional,	/*!< minimum allowed SPH smoothing length in units of SPH gravitational
				   softening length */
    MinGasHsml;			/*!< minimum allowed SPH smoothing length */

#ifdef MAXHSML
  double MaxHsml;
#endif

  double SofteningGas,		/*!< for type 0 */
    SofteningHalo,		/*!< for type 1 */
    SofteningDisk,		/*!< for type 2 */
    SofteningBulge,		/*!< for type 3 */
    SofteningStars,		/*!< for type 4 */
    SofteningBndry;		/*!< for type 5 */

  double SofteningGasMaxPhys,	/*!< for type 0 */
    SofteningHaloMaxPhys,	/*!< for type 1 */
    SofteningDiskMaxPhys,	/*!< for type 2 */
    SofteningBulgeMaxPhys,	/*!< for type 3 */
    SofteningStarsMaxPhys,	/*!< for type 4 */
    SofteningBndryMaxPhys;	/*!< for type 5 */

  double SofteningTable[6];	/*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];	/*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];


  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100],
    InfoFile[100], TimingsFile[100], TimebinFile[100], RestartFile[100], ResubmitCommand[100],
    OutputListFilename[100];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength;		/*!< number of times stored in table of desired output times */



#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
  double ReferenceGasMass;
#endif


#ifdef ADAPTIVE_FORCE_ACC
  double ErrTolForceAccParam;
#endif


#ifdef DISTORTIONTENSORPS
  /* present day velocity dispersion of DM particle in cm/s (e.g. Neutralino = 0.03 cm/s) */
  double DM_velocity_dispersion;
  double TidalCorrection;
#ifdef GDE_LEAN
  double GDEInitStreamDensity;
#endif
#endif

#ifdef SFR			/* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
#ifdef LT_STELLAREVOLUTION
  double OrigGasMass;
#endif
  double EgySpecSN;
  double FactorSN;
  double EgySpecCold;
  double FactorEVP;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double MaxSfrTimescale;
  double FactorForSofterEQS;
#ifdef WINDS
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelMaxTimeFactor;	/* maximum free travel time in units of the Hubble time at the current simulation redshift */
  double WindFreeTravelDensFac;
#ifdef VARIABLE_WINDS
  double VariableWindVelFactor;	/* wind velocity in units of the halo escape velcoity */
  double VariableWindSpecMomentum;	/* momentum available for wind per unit mass of stars formed, in internal velocity units */
  double HaloConcentrationNorm;	/* concentration c0 of a halo of unit mass */
  double HaloConcentrationSlope;	/* slope n of mass concentration relation, namely c = c0 * M_200,crit^n */
#endif
#endif
#endif

#if defined(SFR) || defined(BLACK_HOLES)
  unsigned int bits;
#endif

#ifdef DARKENERGY
  double DarkEnergyParam;	/*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */
#ifdef TIMEDEPGRAV
  double Gini;
#endif
#endif
#endif

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
#endif

#if defined(SNIA_HEATING)
  double SnIaHeatingRate;
#endif

#ifdef TIME_DEP_ART_VISC
  double ViscSource0;		/*!< Given sourceterm in viscosity evolution */
  double DecayLength;		/*!< Number of h for the viscosity decay */
  double ViscSource;		/*!< Reduced sourceterm in viscosity evolution */
  double DecayTime;		/*!< Calculated decaytimescale */
  double AlphaMin;		/*!< Minimum of allowed viscosity parameter */
#endif

#ifdef CONDUCTION
  double ConductionCoeff;	/*!< Thermal Conductivity */
#ifdef CONDUCTION_SATURATION
  double ElectronFreePathFactor;	/*!< Factor to get electron mean free path */
#endif

  integertime Conduction_Ti_begstep, Conduction_Ti_endstep;
  double MaxSizeConductionStep;

#ifdef CONDUCTION_INCLUDEMAGNETIC_PERPENDICULAR
  double ConductionPerpendicularConstants;
#endif
#endif				// CONDUCTION

#ifdef LMB_SPECTRAL_CRs_DIFFUSION
  integertime CRDiffusion_Ti_begstep, CRDiffusion_Ti_endstep;
  double MaxSizeCRDiffusionStep;

  double CR_Kappa_10k;
  double CR_AlphaKappa;
#ifdef LMB_SPECTRAL_CRs_DIFFUSION_INCLUDEMAGNETIC_PERPENDICULAR
  double CRDiffusionPerpendicularConstants;
#endif
#endif

#ifdef VSMOOTH
  double VSmoothScale;
#endif

#ifdef MAGNETIC
#ifdef BINISET
  double BiniX, BiniY, BiniZ;	/*!< Initial values for B */
#endif

#if defined(BSMOOTH) || defined(BSMOOTH_TIME)
  int BSmoothInt;
  double BSmoothFrac;
  int MainTimestepCounts;
#ifdef SETMAINTIMESTEPCOUNT
  int MainTimestepCountIni;
#endif
#endif				/* BSMOOTH || BSMOOTH_TIME */

#if defined(MAGNETIC_DISSIPATION)
#define MAGNETIC_SIGNALVEL	/*!< Make explicit the dependence with MAGNETIC_SIGNALVEL as described in Dolag2009 */
  double ArtMagDispConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial magnetic disipation */
#ifdef TIME_DEP_MAGN_DISP
  double ArtMagDispMin;
  double ArtMagDispSource;
  double ArtMagDispTime;
#endif
#endif

#ifdef MAGNETIC_SN_SEEDING
  double SnSeedRadius;
  double SnSeedBubble;
  double SnSeedField;
  double SnSeedSoftening;
#endif

#ifdef DIVBCLEANING_DEDNER
  double DivBcleanParabolicSigma;
  double DivBcleanHyperbolicSigma;
  double DivBcleanQ;
#endif


#ifdef MAGNETIC_DIFFUSION
  double MagneticEta;
#endif
#ifdef MAGNETIC_DIFFUSION_LIMIT
  double MagneticDiffSpeed;
#endif
#endif				/* MAGNETIC */

  int LevelOfStrickness;	// Defines how strict the code is

  int BlackHolesOn;		// Switches black hole module on, see Springel&DiMatheo 2005
#ifdef BLACK_HOLES
  int BlackHoleThermalFeedbackOn;	// Switches on Thermal Feedback
  int BlackHoleKineticFeedbackOn;	// Switches on Kinetic Feedback
  double massDMpart;		// DM particle mass to use, unless a fixed mass is given in the header
  double TimeNextOnTheFlyFoF;	//
  double TimeBetOnTheFlyFoF;	// Frequency of Black Hole seeding
  double MinFoFMassForNewSeed;	// Minmum FoF mass for seeding Black Holes (internal units)
  double SeedBlackHoleMass;	// Mass of seeded Black Hole
  double BlackHoleSeedStarMassFraction;	// If > 0, seeding is on stellar body with reduced mass
  double BlackHoleSeedDMFraction;	// Star to DM fraction to exceed for seeding BH
  double BlackHoleSeedGasFraction;	// Gas fraction needede for seeding BH
  double BHfactor;		// Maximum number of black holes per gas particle to form between restarts
  double BlackHoleNgbFactor;	// Boost of Black Hole neighbour number
  double BlackHoleMaxAccretionRadius;	// Maximum Black Hole accreation radius (internal units)
  int BlackHoleSwallowGasOn;	// Removal of gas particles through black hole accretion (0=Off, 1=All, 2=Cold gas only)
  int BlackHoleIgnoreMomentum;	// Switches On/Off the transfer of momentum while swallowing
  int BlackHoleAccretionSlicesOn;	// Quantization for swallowing gas particles with GENERATIONS
  int BlackHoleDetails;		// Switches On/Off writing of black hole details by each MPI rank
  int BlackHoleRepositioningOn;	// Switches On/Off black hole repositioning, see Springel&DiMatheo 2005

  int BlackHoleHotColdAccretionOn;	// Switches On/Off hot/cold accreation, see Steinborn+ 2015
  double BlackHoleAccretionFactor;	// Used as alpha Factor in Bondi for total (or hot gas) accreation
  double BlackHoleColdAccretionFactor;	// Used for alpha for the cold gas accreation, values see Gaspari+ 2013
  int BlackHoleVelocityInBondiOn;	// Use the relative velocity in Bondi Formulae
  int BlackHoleVariableAccretionFactorOn;	// Switches on the scaled accreation factor, see Booth&Schaye 2013
  double BlackHoleVariableAccretionSlope;	// Slope for variable accreation factor

  int BlackHoleVariableEfficiencyOn;	// Switches On/Off mass dependent  accreation, see Steinborn+ 2015
  double BlackHoleRadiativeEfficiency;	// Radiative efficiency for standard model
  double BlackHoleRadiativeEfficiencySlope;	// Slope for mass dependent efficiency (see Steinborn+ 2015)
  double BlackHoleRadiativeEfficiencyNorm;	// Norm for mass dependent efficiency (see Steinborn+ 2015)
  double BlackHoleRadiativeEfficiencySlopeMdot;	// Slope for mdot-mbh dependent efficiency
  double BlackHoleRadiativeEfficiencySlopeMbh;	// Slope for mdot-mbh dependent efficiency
  double BlackHolePhysMdotUpperLimit;	// Maximum limit for the accretion rate in Msun/yr
  double BlackHoleRadiativeEfficiencyMax;	// Slope for mass dependent efficiency (see Steinborn+ 2015)
  double BlackHoleRadiativeEfficiencyMin;	// Norm for mass dependent efficiency (see Steinborn+ 2015)
  double BlackHoleEddingtonFactor;	// Maximum allowd accreation in terms of Eddington accreation

  int BlackHoleOutflowModelOn;	// Switches On/Off radiation outflow model, see Steinborn+ 2015
  double BlackHoleFeedbackFactor;	// Fixed feedback factor for standard model
  double BlackHoleRadioModeFeedbackBoost;	// Boosts for FeedbackFactor in radi mode, see Fabjan+ 2010
  double BlackHoleRadioModeTreshold;	// Treshhold for switching to radio mode
  int BlackHoleLimitFeedbackOn;	// Limits feedback (0=Off, 1=to hot particles only, 2=maximum pressure given by halo, see Vogelsberger+ 2014, 3=1+2)
  double BlackHoleLimitMaxTemp;	// Limits feedback to heat individual particle to maximum temperature
  double BlackHoleColdGasTemperatureTresh;	// Temperature Limit to seperate cold gas from hot particles (e.g. 5e4)

  int BlackHoleFrictionForceOn;	// Switches On/Off friction force for black holes, see Hirschmann+ 2014
  int BlackHoleFrictionForceDynaicsOn;	// Switches On/Off friction force using subfind details for black holes
  double BlackHoleFrictionForceMaxCorrect;	// Sets a maximum to the friction correction

  double BlackHoleMergeDistFrac;	// Fraction of relative velocity compared to sound velocity for merging
  double BlackHoleMergeCsndFrac;	// Distance compared to softening length for merging, see Hirschmann+ 2014
  double BlackHoleMergeBindingFrac;	// Fraction of binding energy for merging, see Hirschmann+ 2014
#endif



#ifdef REIONIZATION
  int not_yet_reionized;	/*!< flag that makes sure that there is only one reionization */
#endif


#ifdef NAVIERSTOKES
  double NavierStokes_ShearViscosity;
  double FractionSpitzerViscosity;
  double ShearViscosityTemperature;
#endif

#ifdef NAVIERSTOKES_BULK
  double NavierStokes_BulkViscosity;
#endif

#ifdef NAVIERSTOKES_VISCOSITY_SATURATION
  double IonMeanFreePath;
#endif


#ifdef RELAXOBJECT
  double RelaxBaseFac;
  double RelaxFac;
#endif

  int SpectralCRsOn;

#ifdef LMB_SPECTRAL_CRs

  int CR_ProtonsOn;		/*!< evolve CR protons */
  int CR_ElectronsOn;		/*!< evolve CR electrons */
  int CR_Details;		/*!< output CR properties by every MPI rank */
  int CR_TrackParticle;		/*!< if != 0 tracks particle of that ID */
  int CR_AdiabaticHighUpdateRate;	/*!< Flag to update spectrum on every step */
  int CR_SpektrumHighUpdateRate;	/*!< Flag to update spectrum on every step */

#ifdef JD_SHOCK
  int CR_DSA_Model;		/*!< Flag for diffuse shock acceleration model */
#ifdef LMB_CRs_CALCULATE_P_INJ
  double CR_Chi_pinj;		/*!< Multiplication factor for Q * p_thermal */
#endif
#endif

#ifdef LMB_CR_PROTONS
  double CRp_pmin, CRp_pmax;	/*!< min and max momenta of cosmic ray protons */
  double CRp_bound[LMB_CR_PROTONS + 1];	/*!< boundaries of proton momentum bins */
#endif				// LMB_CR_PROTONS
#ifdef LMB_CR_ELECTRONS
  double CRe_pmin, CRe_pmax;	/*!< min and max momenta of cosmic ray electrons */
  double CRe_bound[LMB_CR_ELECTRONS + 1];	/*!< boundaries of electron momentum bins */
#endif				// LMB_CR_ELECTRONS

  double CR_SlopeBig;		/*!< maximum slope of momentum bin */
  double CR_Dpp;		/*!< constant momentum diffusion coefficient */
  double CR_EpsilonB;		/*!< ratio of the general magn. field along the shock normal to the amplitude
				   of the downstream MHD wave turbulence (Kang & Ryu 2011) */
  double CR_DSlope;		/*!< Safety net to avoid energy from diverging */
  double CRe_Ratio;		/*!< Ratio of injected e^- per inj. p */
  double CR_T0Inj;		/*!< Starting time of CR injection */
  double CR_bound[LMB_SPECTRAL_CRs + 1];	/*!< boundaries of cosmic rays momentum bins */
  double CR_Gamma;		/*!< adiabatic index of CRs  */

#ifdef LMB_SPECTRAL_CRs_SEED
#ifdef LMB_CR_PROTONS
  double CRp_InitSlope;		/*!< init slopes of spectra of protons */
#endif				// LMB_CR_PROTONS
#ifdef LMB_CR_ELECTRONS
  double CRe_InitSlope;		/*!< init slopes of spectra of electrons */
#endif				// LMB_CR_ELECTRONS
#endif				// LMB_SPECTRAL_CRs_SEED

#ifdef LMB_SPECTRAL_CRs_ARTIFICIAL_CONDUCTIVITY
  double CR_ArtCondConstant;
#endif
#endif

#ifdef LT_STELLAREVOLUTION
  int MaxPartMet;		/*!< This gives the maxmimum number of STAR particles that can be stored on one
				   processor. */

  double Time_Age;		/*!< current cosmic time in Gyrs */

  int TestSuite, StarBits;
  double SFfactor;		/*!< the expected maximum factor of proportionality between TotN_gas and TotN_star */
  int Generations;

  int NeighInfNum;		/*!< minimum number of neighbours for metal and egy spreading */
  int DesNumNgbSN,		/*!< desired number of neighbours for sn spreading */
    SpreadNumNgbDev,		/*!< maximum deviation in number of neighbours for sn spreading */
    LeftNumNgbSN, RightNumNgbSN;	/*!< range of neighbours for sn spreading */
  double SpreadNeighCoeff;	/*!< store DesNumNgbSN / DesNumNgb */

  double MinChemSpreadL;
  double MinChemTimeStep;
  double Enrich_SFGas_Th;

  double LocalSpreadFactor;

#ifdef GADGET3_IO_LIB
  long double MetalsCheckSum[LT_NMetP];
#else
  long double MetalsCheckSum[LT_NMetP] = { 0 };
#endif  
  char SFfilename[300], IMFfilename[300], SnIaDataFile[300], SnIIDataFile[300], AGBDataFile[300];
  int Ia_Nset_ofYields, II_Nset_ofYields, AGB_Nset_ofYields;

  double Mup,			/*!< SnII threshold (e.g. 8 Msun) */
    MBm,			/*!< min bin. system mass */
    MBM,			/*!< max bin. system mass */
    BinFrac,			/*!< bin. system fraction */
    MBms;			/*!< min star mass in SnIa binary systems */
  double SnIaRemn;		/*!< SnIa Remn (1.4Msun ?) */
  /* lifetimes */
  double inf_lifetime, mean_lifetime, sup_lifetime;
  /* energy provided by sn explosions */
  double SnIaEgy, SnIIEgy;
  /* define IRA range */
  double metIRA_ThMass, egyIRA_ThMass;
  double SnII_Step_Prec, LLv_Step_Prec;
  double referenceZ_toset_SF_DensTh;
  int SFTh_Zdep, referenceZbin_SFTh;

#ifdef LT_POPIII
  double PopIII_Zlimit;
  int PopIII_IMF_idx;
#endif

#ifdef LT_HOT_EJECTA
  double EgySpecEjecta;
#endif

#ifdef LT_STARBURSTS
  int StarBurstCondition;
  double SB_Density_Thresh;
  double SB_DEntropy_Thresh;
#endif

#ifdef LT_STOP_COOL_BELOW_Z
  double Below_this_redshift_stop_cooling;
#endif

#ifdef LT_SMOOTH_Z
#if defined(LT_SMOOTH_SIZE) && !defined(LT_SMOOTH_NGB)
  double SmoothRegionSize, SmoothRegionSizeMax;
#endif
#if defined(LT_SMOOTH_NGB) && !defined(LT_SMOOTH_SIZE)
  int DesNumNgbSmooth;
#endif
#endif

  double MaxChemSpreadL;

#endif

#ifdef LT_METAL_COOLING_WAL
  char WalCool_CoolTables_path[200];
#endif

#ifdef LT_ADD_GAL_TO_SUB
  char BC_SED_Path[200];
#endif

#ifdef GENERATE_GAS_IN_ICS
#ifdef GENERATE_GAS_TG
  int GenGasRefFac;
#endif
#endif

#if defined (UM_CHEMISTRY) && defined (UM_CHEMISTRY_INISET)
  /* used if read initial composition from the parameter file */
  double Startelec;
  double StartHI, StartHII, StartHM;
  double StartHeI, StartHeII, StartHeIII;
  double StartH2I, StartH2II;
  double StartHD, StartDI, StartDII;
  double StartHeHII;
#endif

#ifdef GM_MUPPI			/*GM: NOTE: currently unused */
  int MuppiOn;
  int MuppiDebugOn;
  int CountParticlesInConeOn;
  int StellarKinFB2On;
  int StellarKinFB2OutputOn;

#ifdef LT_STELLAREVOLUTION
  double beta_sf[10];
  double f_rest[10];
  double E_SN_51[10];
#endif
  MyAtLeastDouble FracEgyIn;
  MyAtLeastDouble FracEgyOut;
  MyAtLeastDouble FracEgyKin;
  MyAtLeastDouble StarCM[3];
  MyAtLeastDouble Mstar;	/* used in run for print out star pos */
  int Nstar;
#ifdef GM_MUPPI_DEBUG
  MyAtLeastDouble Egy_out, Egy_in_thermal, Egy_in_kinetic;
  MyAtLeastDouble Egy_outCum, Egy_in_thermalCum, Egy_in_kineticCum;
#endif
#if defined (MV_GM_AGNMUPPI) && !defined(MV_GM_COVER_FACT_MCLS)
  double CouplingBHEnHot;
  double CouplingBHEnCold;
#endif
#ifdef MV_GM_AGNMUPPI_OUTPUT
  double TotE_AGNM_tot, TotE_AGNM_cool, TotE_AGNM_sfr_c, TotE_AGNM_sfr_h;
#endif
#endif


#ifdef MOL_CLOUDS
  unsigned int MOL_CLOUDS_NumMCs;
  unsigned int MOL_CLOUDS_NumStars;
  double MOL_CLOUDS_MassMCs;
  double MOL_CLOUDS_MassStars;

  MyFloat MOL_CLOUDS_m1_tilde, MOL_CLOUDS_m2_tilde;
#endif

#ifdef ADAPTGRAVSOFT
  int AGS_DesNumNgb;
  double AGS_MaxNumNgbDeviation;
#ifdef PMGRID
  double AGS_MaxSoft[2], AGS_MinSoft[2];
#endif
#endif

#ifdef fSIDM
  double FSIDM_CrossSectionPerMass_in_cgs;	/*!< Read from parameter file */
#ifdef mSIDM_TIMESTEP
  double FSIDM_tfact;
#endif
#endif
#ifdef rSIDM
  double RSIDM_CrossSectionPerMass_in_cgs;	/*!< Read from parameter file */
#ifdef mSIDM_TIMESTEP
  double RSIDM_tfact;
#endif
#ifdef rSIDM_DISSIPATION
  double RSIDM_f_diss;
#endif
#endif
#ifdef mSIDM_TIMESTEP
#if !defined(mSIDM_VDEP0)
  double MSIDM_OmegaInit;
#endif
#endif
#if defined(fSIDM) || defined(rSIDM)
  unsigned int MSIDM_SEED;	/*!< For random angles in scattering */
  unsigned int MSIDM_MaxIntOrder;
  unsigned int MSIDM_IntOrder;
#ifdef mSIDM_VDEP0
  double MSIDM_vdep_w;
  double MSIDM_vdep_alpha;
  double MSIDM_vdep_beta;
#endif
#ifndef ADAPTGRAVSOFT
  int DesNumNgb_pt1;
  double MaxNumNgbDeviation_pt1;
#endif
#endif
#if GADGET_HYDRO == HYDRO_MFM
  double EkinSwitchFraction;
  double EpotSwitchFraction;
#endif
#ifdef SPHERICAL_BOUNDARY
  struct
  {
    double radius;
    int type;
  } spherical_boundary;
#endif

}

All;


#if defined(fSIDM) || defined(rSIDM)
extern class LambdaNumerical calc_lambda;
extern class GaussLegendre gauss_legendre;
extern class MSIDM_Scattering msidm_scatter;
#endif


/*! This structure contains data which is the SAME for all tasks as the All data structure
 * but its values are computed on the fly thus it can change in shape betwen restarts
 */

#ifndef GADGET3_IO_LIB

extern struct derived_global_data_all_processes
{
#ifdef DOMAIN_DECOMPOSITION_ON_TB
  int DoDomainDecompositionTBFull, DoDomainDecompositionTBDiffuse, DoDomainDecompositionTBFix;
#endif
} AllDerived;


#endif


/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern ALIGN(32)
     struct particle_data
     {
       MyLongDouble Pos[3];	/*!< particle position at its current time */
       short int Type;		/*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
       short int TimeBin;
       MyIDType ID;

       integertime Ti_begstep;	/*!< marks start of current timestep of particle on integer timeline */
       integertime Ti_current;	/*!< current time of the particle */

       MyFloat Mass;		/*!< particle mass */

       MyFloat Vel[3];		/*!< particle velocity at its current time */
       MyFloat dp[3];

       MyFloat GravAccel[3];	/*!< particle acceleration due to gravity */
#ifdef PMGRID
       MyFloat GravPM[3];	/*!< particle acceleration due to long-range PM gravity force */
#endif

#ifdef FORCETEST
       MyFloat GravAccelDirect[3];	/*!< particle acceleration calculated by direct summation */
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
       MyFloat Potential;	/*!< gravitational potential */
#ifdef KD_USE_IPOT_FOR_BH_MERGER
       MyFloat IPotential;	/* Integrated potential for BH mergers  */
#endif
#endif

#ifdef DISTORTIONTENSORPS
       MyBigFloat distortion_tensorps[6][6];	/*!< phase space distortion tensor */
       MyBigFloat last_determinant;	/*!< last real space distortion tensor determinant */
       MyBigFloat stream_density;	/*!< physical stream density that is going to be integrated */
       double tidal_tensorps[3][3];	/*!< tidal tensor (=second derivatives of grav. potential) */
       float caustic_counter;	/*!< caustic counter */
#ifndef GDE_LEAN
       MyBigFloat annihilation;	/*!< integrated annihilation rate */
       MyBigFloat analytic_annihilation;	/*!< analytically integrated annihilation rate */
       MyBigFloat rho_normed_cutoff_current;	/*!< current and last normed_cutoff density in rho_max/rho_init * sqrt(sigma) */
       MyBigFloat rho_normed_cutoff_last;
       MyBigFloat s_1_current, s_2_current, s_3_current;	/*! < current and last stretching factor */
       MyBigFloat s_1_last, s_2_last, s_3_last;
       MyBigFloat second_deriv_current;	/*! < current and last second derivative */
       MyBigFloat second_deriv_last;
       double V_matrix[3][3];	/*!< initial orientation of CDM sheet the particle is embedded in */
       float init_density;	/*!< initial stream density */
       float analytic_caustics;	/*!< number of caustics that were integrated analytically */
       float a0;
#endif


#ifdef OUTPUT_LAST_CAUSTIC
       MyFloat lc_Time;		/*!< time of caustic passage */
       MyFloat lc_Pos[3];	/*!< position of caustic */
       MyFloat lc_Vel[3];	/*!< particle velocity when passing through caustic */
       MyFloat lc_rho_normed_cutoff;	/*!< normed_cutoff density at caustic */
       MyFloat lc_Dir_x[3];	/*!< principal axis frame of smear out */
       MyFloat lc_Dir_y[3];
       MyFloat lc_Dir_z[3];
       MyFloat lc_smear_x;	/*!< smear out length */
       MyFloat lc_smear_y;
       MyFloat lc_smear_z;
#endif
#ifdef PMGRID
       double tidal_tensorpsPM[3][3];	/*!< for TreePM simulations, long range tidal field */
#endif
#endif
#ifdef AR_FIX_ASSIGN
//particles with only ngb from other CPUs will 
//assign wrong values on non summed values in out2particle density (e.g. .BH_TimeBinGasNeighbor)
       int assigned;
#endif


       MyFloat OldAcc;		/*!< magnitude of old gravitational force. Used in relative opening
				   criterion */
#if defined(EVALPOTENTIAL) && defined(PMGRID)
       MyFloat PM_Potential;
#endif

#if (defined(fSIDM) || defined(rSIDM)) && defined(mSIDM_TIMESTEP)
#if !defined(mSIDM_VDEP0)
#ifdef fSIDM
       MyFloat omega_frequent;
#endif
#ifdef rSIDM
       MyFloat omega_rare;
#endif
#endif
#endif

#ifdef STELLARAGE
#if !defined(LT_STELLAREVOLUTION)
       MyFloat StellarAge;	/*!< formation time of star particle */
#endif
#endif
#ifdef METALS
       MyFloat Metallicity;	/*!< metallicity of gas or star particle */
#endif				/* closes METALS */

       MyFloat Hsml;
       MyFloat NumNgb;
#if  defined(SNIA_HEATING)
       MyFloat DensAroundStar;
#endif

#if defined(DISKPOT) && !defined(SUBFIND) && !defined(ORDER_SNAPSHOTS_BY_ID)
       int GrNr;
#endif



#if defined(SUBFIND)
       int GrNr;
       int SubNr;
       int DM_NumNgb;
#if defined(KD_MAIN_HALO)
       int SubLevel;
#endif
       unsigned short targettask, origintask2;
       int origintask, submark, origindex;
       MyFloat DM_Hsml;
       union
       {
	 MyFloat DM_Density;
	 MyFloat DM_Potential;
       } u;
       union
       {
	 MyFloat DM_VelDisp;
	 MyFloat DM_BindingEnergy;
       } v;
#ifdef DENSITY_SPLIT_BY_TYPE
       union
       {
	 MyFloat int_energy;
	 MyFloat density_sum;
       } w;
#endif

#ifdef SAVE_HSML_IN_IC_ORDER
       MyIDType ID_ic_order;
#endif
#ifdef SUBFIND_ALTERNATIVE_COLLECTIVE
       peanokey Key;
#endif
#endif

#if defined(ORDER_SNAPSHOTS_BY_ID) && !defined(SUBFIND)
       int GrNr;
       int SubNr;
#endif



       float GravCost[GRAVCOSTLEVELS];	/*!< weight factor used for balancing the work-load */

int dt_step;

       union
       {
	 unsigned int BHID;
	 unsigned int MetID;
       } pt;

#ifdef SCF_HYBRID
       MyDouble GravAccelSum[3];
       MyFloat MassBackup;
#endif

#ifdef MOL_CLOUDS
       MyFloat MOL_CLOUDS_TimeBorn;
       MyFloat MOL_CLOUDS_LifeTime;

       unsigned int MOL_CLOUDS_index;
#endif

       int TrueNGB;		/*!< Number of neighbours inside hsml */

#ifdef ADAPTGRAVSOFT
       MyFloat AGS_Density;	/* !< (mass/number) density of particle */
#ifdef AGS_OUTPUTNGBS
       int AGS_REALNumNgb;
#endif
       MyFloat AGS_NumNgb;
       MyFloat AGS_Hsml;
       MyFloat AGS_zeta;	/*!< factor in the correction term */
       MyFloat AGS_omega;	/*!< factor in the correction term */
#ifdef AGS_OUTPUTCORR
       MyFloat AGS_corr;	/*!< factor in the correction term */
#endif
#endif

#if defined(FOF_EXTENDED_PROPERTIES) &&  !defined(SUBFIND)
       int GrNr;
       unsigned short targettask, origintask;
#endif

#if defined(GM_STARDENSITY) || defined(GM_USE_STARDENS_AND_POTMIN)
       MyFloat StarHsml;

       MyFloat StarNumNgb;
       MyFloat StarDensity;	/* SPH-smoothed stellar density of star particles */

#ifdef EVALPOTENTIAL
       MyFloat SmoothPot;	/* SPH-smoothed grav potential of star  particles 
				   WARNING, smoothed on  neighbouring star particles ONLY */
#endif

#endif


     }
 *__restrict__ P,		/*!< holds particle data on local processor */
*__restrict__ DomainPartBuf;	/*!< buffer for particle data used in domain decomposition */

#ifndef GDE_LEAN
#define GDE_TIMEBEGIN(i) (P[i].a0)
#define GDE_VMATRIX(i, a, b) (P[i].V_matrix[a][b])
#define GDE_INITDENSITY(i) (P[i].init_density)
#else
#define GDE_TIMEBEGIN(i) (All.TimeBegin)
#define GDE_VMATRIX(i, a, b) (0.0)
#define GDE_INITDENSITY(i) (All.GDEInitStreamDensity)
#endif

#ifdef LT_STELLAREVOLUTION
#define MPP(i) MetP[P[(i)].pt.MetID]
#else
#define MPP(i) P[(i)]
#endif

#if defined(BLACK_HOLES)

#define BPP(i) BHP[P[(i)].pt.BHID]

     extern struct bh_particle_data
     {
       unsigned int PID;

       MyFloat StellarAge;	/*!< formation time of star particle */

       MyIDType SwallowID;

       MyFloat SwallowPot;

       int BH_CountProgs;

       MyFloat BH_Mass;
       MyFloat BH_Mdot;
       int BH_TimeBinGasNeighbor;
       MyFloat BH_Density;
       MyFloat BH_Entropy;
       MyFloat BH_SurroundingGasVel[3];

       MyFloat BH_ColdDensity;
       MyFloat BH_ColdEntropy;
       MyFloat BH_SurroundingColdGasVel[3];

       MyFloat BH_HotDensity;
       MyFloat BH_HotEntropy;
       MyFloat BH_SurroundingHotGasVel[3];

#ifdef LB_PRESSURE_DEPENDENT_ACCRETION
       MyFloat BH_SurroundingMetallicity;
       MyFloat BH_SurroundingGasPressure;
       MyFloat BH_SurroundingGasMass;
       MyFloat BH_SurroundingTemperature;
       MyFloat BH_SfrThreshold;
#endif

       MyFloat BH_accreted_Mass;
       MyFloat BH_accreted_BHMass;

       MyFloat BH_accreted_momentum[3];

       MyFloat BH_SurroundingVel[3];
       MyFloat BH_SurroundingDensity;
       MyFloat BH_sigma;
       MyFloat BH_bmax;

       MyFloat BH_TotalFeedbackEfficiency;

       MyFloat mean_hsml;
       MyFloat mean_rho;

       MyLongDouble BH_MinPotPos[3];
       MyFloat BH_MinPot;

#if defined(GM_MUPPI) && defined(GM_MV_ANGMOM)	/* GM note: this will be used in BH+MUPPI */
       MyFloat Vphi[3], VphiMod;
#endif
#ifdef AD_DYNFRIC
       MyFloat BHDF_DynFric[3];
       int     BHDF_ngb;
#endif

#if defined(GM_REPOSITION_ON_STARDENSITY_MAX) || defined(GM_USE_STARDENS_AND_POTMIN)
       MyLongDouble BH_StarDensPos[3];
       MyFloat BH_StarDens;
#ifdef GM_USE_VELS_ON_REPOSITIONING
       MyLongDouble BH_StarDensVel[3];
#endif
#endif

     }
 *__restrict__ BHP, *__restrict__ DomainBHBuf;

#endif // BLACK_HOLES
					       /* [----------- start LT block ------------- ] */
#ifdef LT_STELLAREVOLUTION	/* [LT] define the structure which hosts the stellar data] */

     extern struct met_particle_data
     {
       float LastChemTime;	/*!< the last and next next time of evolution for Ia and II */
       float iMass;		/*!< initial mass of this SSP                               */
       float Metals[LT_NMetP];	/*!< the metal array (H is not stored here)                 */
       double weight;		/*!< used in spreading                                      */
       unsigned int PID;

#ifdef STELLARAGE
       MyFloat StellarAge;	/*!< formation time of star particle */
#endif

#ifdef LT_TRACK_CONTRIBUTES
       Contrib contrib;
#endif
#ifdef LT_ZAGE
       float ZAge;
#endif
#ifdef LT_ZAGE_LLV
       float ZAge_llv;
#endif
       int ChemTimeBin;
#ifdef LT_STARS_GUESSHSML
       MyFloat mean_hsml;
       MyFloat mean_rho;
#endif
#ifdef LT_SMOOTH_Z
       float Zsmooth;
#endif

     }
 *__restrict__ MetP,		/*!< holds metal particle data on local processor */
*__restrict__ DomainMetBuf;	/*!< buffer for metal data in domain decomposition */
/* [----------- end LT block ------------- ]*/
#endif


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
     extern struct sph_particle_data
     {
       MyLongDouble Entropy;	/*!< entropy (actually entropic function) of particle */
       MyLongDouble EntropyPred;	/*!< predicted value of the entropy at the current time */
       MyFloat VelPred[3];	/*!< predicted SPH particle velocity at the current time */
       MyFloat MaxSignalVel;	/*!< maximum signal velocity */

       MyLongDouble Density;
       MyLongDouble DtEntropy;
       MyLongDouble DensityNorm;

       MyLongDouble HydroAccel[3];
       MyFloat DivVel;		/*!< local velocity divergence */

#if GADGET_HYDRO == HYDRO_SPH
       MyLongDouble DhsmlDensityFactor;
       MyFloat Pressure;	/*!< current pressure */
#elif GADGET_HYDRO == HYDRO_PESPH
       MyFloat NumDens;		/* SPH number density */
       MyFloat Pressure;	/* Thermal pressure (summation variable in Pressure-Entropy SPH) */
       MyFloat DhsmlNumDensFactor;	/*!< correction factor needed in Pressure-Entropy formulation of SPH */
       MyFloat DhsmlPressureFactor;	/*!< correction factor needed in Pressure-Entropy formulation of SPH */
#elif GADGET_HYDRO == HYDRO_MFM
       MyFloat NumDens;		/*!< SPH number density */
       MyFloat DhsmlNumDensFactor;	/*!< correction factor needed in Pressure-Entropy formulation of SPH */
       MyFloat Pressure;	/*!< Thermal pressure (TODO : Only added as stop-gap to prevent compiltion errors) */
       MyFloat DistNgbSqdMax;	/*!< Maximum distance squared to all neighbours (Used in Scalar slope limiter) */
       MyFloat AlphaSlope[5];	/*!< Slope limiter weight factor */
       QFluidVector < NUMDIMS > dQ;	/*!< Partial sums of conserved fluid vector */
       QFluidVector < NUMDIMS > dQdt;	/*!< Rate of change of conserved quantities */
       WFluidVector < NUMDIMS > Wprim;	/*!< Primitive fluid vector */
       WFluidVector < NUMDIMS > Wmax;	/*!< Max. values of W for all neighbours (used in slope limiter) */
       WFluidVector < NUMDIMS > Wmin;	/*!< Min. values of W for all neighbours (used in slope limiter) */
       SquareMatrix < NUMDIMS > Bgrad;	/*!< Matrix required for computing accurate gradients */
       Vector < NUMDIMS > Vel_smooth;	/*!< averaged velocity over neighbours */
       Vector < NUMDIMS > gradW[NUMDIMS + 2];	/*!< Gradients of primitive quantities */
       MyFloat ConditionNumber;	/*!< Condition number of gradient calculation */
       MyFloat MaxVelSquareNgb;	/*!< Maximum velocity squared (w.r.t. neighbours), used for energy/entropy switch */
       MyFloat InternalEnergy;	/*!< internal energy of particle, used for evolution instead of entropy */
       MyFloat InternalEnergyPred;	/*!< prediced internal energy at the current time */
       MyFloat DtInternalEnergy;	/*!< rate fo change of the internal energy */
#endif

       union
       {
	 MyFloat CurlVel;	/*!< local velocity curl */
	 MyFloat Rot[3];	/*!< local velocity curl */
       } r;

#ifdef NAVIERSTOKES
       union
       {
	 MyFloat dvel[3][3];
	 struct
	 {
	   MyFloat DivVel;
	   MyFloat CurlVel;
	   MyFloat StressDiag[3];
	   MyFloat StressOffDiag[3];
#ifdef NAVIERSTOKES_BULK
	   MyFloat StressBulk;
#endif
	 } s;
       } u;
#endif


       MyIDType SwallowID;

       MyFloat Injected_BH_Energy;

#ifdef COOLING
       MyFloat elec;		/*!< electron fraction, expressed as local electron number
				   density normalized to the hydrogen number density. Gives
				   indirectly ionization state and mean molecular weight. */
#endif

#ifdef SFR
       MyFloat Sfr;
#endif
#ifdef WINDS
       MyFloat DelayTime;	/*!< remaining maximum decoupling time of wind particle */
#endif

#ifdef JD_VTURB
       MyFloat Vrms;		/*!< RMS velocity inside kernel around Vbulk */
       MyFloat Vbulk[3];	/*!< Mean velocity inside kernel */
#ifdef JD_DECOMPOSE_VTURB
       MyFloat Vtan;		/*!< RMS of tangential component of velocity inside kernel around Vbulk */
       MyFloat Vrad;		/*!< RMS of radial component of velocity inside kernel around Vbulk */
#endif
#if (defined(JD_DPP) || defined(LMB_DPP))
       MyFloat Dpp;		/*!< Reacceleration Coefficient as (Cassano+ '04) */
#endif
#endif

#ifdef JD_SHOCK
       MyFloat Shock_N1[3];	/*!< Parallel to shock > */
       MyFloat Shock_Up_V1[3];	/*!< Upstream parallel velocity > */
       MyFloat Shock_Down_V1[3];	/*!< Downstream parallel velocity > */
       MyFloat Shock_Up_Signal;	/*!< Upstream signal velocity > */
       MyFloat Shock_Mach;	/*!< Machnumber > */
       MyFloat Shock_Speed;	/*!< Shock speed > */
       MyFloat Shock_Compress;	/*!< Compression ratio > */
       MyFloat Shock_Up_Weight;	/*!< Upstream weight > */
       MyFloat Shock_Down_Weight;	/*!< Downstream weight > */
       MyFloat Shock_Up_Rho;	/*!< Upstream density > */
       MyFloat Shock_Down_Rho;	/*!< Downstream density > */
       MyFloat Shock_Up_Pressure;	/*!< Upstream pressure > */
       MyFloat Shock_Down_Pressure;	/*!< Downstream pressure > */
#ifdef JD_SHOCK_VELDIV
       MyFloat Shock_N2[3];	/*!< Perpendicular to shock > */
       MyFloat Shock_N3[3];	/*!< Perpendicular to shock > */
       MyFloat Shock_Up_V2[3];	/*!< Upstream perpendicular velocity > */
       MyFloat Shock_Down_V2[3];	/*!< Downstream perpendicular velocity > */
       MyFloat Shock_Up_V3[3];	/*!< Upstream perpendicular velocity > */
       MyFloat Shock_Down_V3[3];	/*!< Downstream perpendicular velocity > */
       MyFloat Shock_Up_Weight23[2];	/*!< Upstream perpendicular weights > */
       MyFloat Shock_Down_Weight23[2];	/*!< Downstream perpendicular weights > */
#endif				// JD_SHOCK_VELDIV
#ifdef LMB_CRs_CALCULATE_P_INJ
       MyFloat Shock_Down_Temperature;
#endif				// LMB_CRs_CALCULATE_P_INJ
#ifdef JD_SHOCK_MAGNETIC
       MyFloat Shock_Mach_Alfven;	/*!< Alfven Machnumber > */
       MyFloat Shock_Up_Alfven;	/*!< Upstream Alfven velocity > */
       MyFloat Shock_Down_Alfven;	/*!< Downstream Alfven velocity > */
#endif				// JD_SHOCK_MAGNETIC
#ifdef LMB_SHOCK_OBLIQUITY
       MyFloat Shock_Obliquity;
#endif
#endif				// JD_SHOCK

#ifdef TIME_DEP_ART_COND
       MyFloat Calpha, GradA[3];

#ifndef NOGRAVITY
       MyFloat Climit, Cgrav[3];
#endif
#endif

#ifdef AB_SHEAR
       MyFloat Chi[3][3];
       MyFloat Xix[3], Xiy[3], Xiz[3];
#endif

#ifdef MAGNETIC
       MyFloat B[3];

#ifdef WINDS
       MyFloat WindB[3];
#endif

       MyFloat BPred[3];

#ifdef MAGNETIC_SN_SEEDING
       MyFloat MagSeed[3];
#endif

#ifdef DIVBFORCE3
       MyFloat magacc[3];
       MyFloat magcorr[3];
#endif
       MyFloat DtB[3];
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP) || defined(DIVBCLEANING_DEDNER)
       MyFloat divB;
       MyFloat reldivB;
#endif
#if defined(BSMOOTH) || defined(BSMOOTH_TIME)
       MyFloat BSmooth[3];
#endif
#ifdef TIME_DEP_MAGN_DISP
       MyFloat Balpha, DtBalpha;
#ifdef AB_ART_DISP
       MyFloat ArtMagDispMatrix[3][3];
#endif
#endif
#ifdef DIVBCLEANING_DEDNER
       MyFloat Phi, PhiPred, DtPhi;
       MyFloat GradPhi[3];
#ifdef SMOOTH_PHI
       MyFloat SmoothPhi;
#endif
#endif
#if defined(SMOOTH_DIVB) || defined(DIVBCLEANING_DEDNER)
       MyFloat SmoothDivB;
#endif

#if defined(ROT_IN_MAG_DIS) || defined(OUTPUT_ROTB)
       MyFloat RotB[3];
#ifdef SMOOTH_ROTB
       MyFloat SmoothedRotB[3];
#endif
#endif

#endif				/* MAGNETIC */

#ifdef VSMOOTH
       MyFloat VSmooth[3];
       int SmoothNgb;
#endif


#ifdef TIME_DEP_ART_VISC
       MyFloat alpha, Dtalpha;
#ifdef AB_ART_VISC
       MyFloat visc_R, visc_OldDvel;
#endif
#endif
#ifdef NAVIERSTOKES_TIMESTEP
       MyFloat ViscEntropyChange;
#endif
#ifdef CONDUCTION_SATURATION
       MyFloat GradEntr[3];
#endif
#if defined(CONDUCTION_INCLUDEMAGNETIC) || defined(AA_OUTPUT_GRADS)
       MyFloat GradEnergy[3];
#endif
#ifdef AA_OUTPUT_GRADS
       MyFloat GradPressure[3];
#endif





#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
       MyFloat HI;
       MyFloat HII;

       MyFloat HeI;
       MyFloat HeII;
       MyFloat HeIII;

       MyFloat H2I;
       MyFloat H2II;

       MyFloat HM;

       MyFloat Gamma;
       MyFloat t_elec, t_cool;

#ifdef UM_CHEMISTRY
       MyFloat Um_MeanMolecularWeight;
       MyFloat HD;
       MyFloat DI;
       MyFloat DII;
       MyFloat HeHII;
#endif

#endif

       //if we restart from a run w.o. soa, we still declare it in the all vars.
#if !defined(AR_SOA_LITE) || defined(AR_SOA_LITE_FROM_RESTART_AOS)
       short int wakeup;	/*!< flag to wake up particle */
#endif


#ifdef LMB_SPECTRAL_CRs
       // Protons
       MyFloat CRpPressure;	/*!< pressure of CR p -> always present! */
#ifdef LMB_CR_PROTONS
       MyAtLeastDouble CRpNorm[LMB_CR_PROTONS];	/*!< normalization of CR protons spectrum */
       MyFloat CRpSlope[LMB_CR_PROTONS];	/*!< slope of CR protons spectrum */
       MyFloat CRpCut;		/*!< cutoff of CR protons spectrum  */
       MyFloat CRpN[LMB_CR_PROTONS];	/*!< number of CR p */
       MyFloat CRpE[LMB_CR_PROTONS];	/*!< energy of CR p */
#endif				// LMB_CR_PROTONS

       // Electrons
       MyFloat CRePressure;	/*!< pressure of CR e -> always present! */
#ifdef LMB_CR_ELECTRONS
       MyAtLeastDouble CReNorm[LMB_CR_ELECTRONS];	/*!< normalization of CR electrons spectrum */
       MyFloat CReSlope[LMB_CR_ELECTRONS];	/*!< slope of CR electrons spectrum */
       MyFloat CReCut;		/*!< cutoff of CR electrons spectrum  */
       MyFloat CReN[LMB_CR_ELECTRONS];	/*!< number of CR e */
       MyFloat CReE[LMB_CR_ELECTRONS];	/*!< energy of CR e */
#ifdef LMB_CR_OUTPUT_SYNCHROTRON
       MyAtLeastDouble SynchEmissivity[8];	/*!< Synchrotron emissivity in erg / cm^3 / s / Hz 
                                                for 60MHz, 150MHz, 235MHz, 325MHz, 610MHz, 944MHz, 1.4GHz and 1.656GHz 
                                            */
#endif				// LMB_CR_OUTPUT_SYNCHROTRON
#endif				// LMB_CR_ELECTRONS

#ifdef LMB_SPECTRAL_CRs_ARTIFICIAL_CONDUCTIVITY
#ifdef LMB_CR_PROTONS
       MyFloat DtCRpE[LMB_CR_PROTONS];	/*!< time derivative of CR p energy */
       MyFloat DtCRpN[LMB_CR_PROTONS];	/*!< time derivative of CR p number */
#endif				// LMB_CR_PROTONS
#ifdef LMB_CR_ELECTRONS
       MyFloat DtCReE[LMB_CR_ELECTRONS];	/*!< time derivative of CR e energy */
       MyFloat DtCReN[LMB_CR_ELECTRONS];	/*!< time derivative of CR e number */
#endif				// LMB_CR_ELECTRONS
#endif				// LMB_SPECTRAL_CRs_ARTIFICIAL_CONDUCTIVITY

#ifdef LMB_SPECTRAL_CRs_DIFFUSION_INCLUDEMAGNETIC
       MyFloat GradCRpPressure[3];	// CR proton pressure gradient 
       MyFloat GradCRePressure[3];	// CR electron pressure gradient 
#endif

       MyFloat DensityOld;

#ifdef LMB_SPECTRAL_CRs_SN_SEEDING
       MyFloat CRsSNe;		/*! Energy from SN to be injected into CRs */
#endif				// LMB_SPECTRAL_CRs_SN_SEEDING

#endif

#if (defined(SFR) && defined(MAGNETIC)) || defined(LT_STELLAREVOLUTION)
       float XColdCloud;
#endif
#if defined(LT_STELLAREVOLUTION) || defined(GM_MUPPI)
       float Temperature;
#endif

#ifdef LT_STELLAREVOLUTION
       MyFloat MassRes;
       MyFloat EgyRes;		/*!< external (Sn) energy resorvoir */
       float Metals[LT_NMetP];	/*!< H is not stored here */

#ifdef GL_CR_DUST
       float DustL[LT_NMetP];  /*!< the metal in large dust grains array (H is not stored here) */
       float DustS[LT_NMetP];  /*!< the metal in small dust grains array (H is not stored here) */
       float numSnIa, numSnII;  /* number of supernovae from the last dust evolution calculation (spreaded!) */
#endif //GL_CR_DUST

#ifdef GM_MUPPI
       MyAtLeastDouble EgyStep;
#endif
#ifndef LT_LOCAL_IRA
       double mstar;
#endif
#ifdef LT_TRACK_CONTRIBUTES
       Contrib contrib;
#endif
#ifdef LT_ZAGE
       MyFloat ZAge, ZAgeW;
#endif
#ifdef LT_ZAGE_LLV
       MyFloat ZAge_llv, ZAgeW_llv;
#endif
#ifdef LT_TRACK_WINDS
       MyFloat AvgHsml;
#endif

#endif

#if defined(LT_SMOOTH_Z)
#if defined(LT_SMOOTH_ALLMETALS)
       MyFloat Zsmooth[LT_NMetP];
#else
       MyFloat Zsmooth;
       MyFloat Zsmooth_a;
       MyFloat Zsmooth_b;
#if defined(LT_SMOOTH_SIZE)
       MyFloat SmoothDens;
       MyFloat SmoothDens_b;
       int SmoothNgb;
#endif
#endif				/* closes LT_SMOOTH_ALLMETALS  */
#endif				/* closes LT_SMOOTH_Z */

#ifdef LT_SMOOTH_XCLD		/* smooth the cloud fraction */
       float XCLDsmooth;
#endif

#ifdef PHIDOT
       MyFloat PhiDot;
#endif

#ifdef GM_MUPPI			/* GM: this variables will be detached from SphP */
       MyAtLeastDouble M_sf;
       MyAtLeastDouble M_h;
       MyAtLeastDouble M_c;
       MyAtLeastDouble Fcoll;
       MyAtLeastDouble E_h;
       MyAtLeastDouble Egy_tot_0;
       MyAtLeastDouble tdyn;
       MyAtLeastDouble clock;
       int MultiPhase;
       int out_part;
       int NMF;
       MyAtLeastDouble Eout_norm;
       MyAtLeastDouble Ekin_rec[3];
       MyAtLeastDouble Ekin_norm;
#ifdef MV_GM_STELLAR_KIN_FB2
       MyAtLeastDouble Ekin_norm_fountain;
       int wait;
#endif
#ifdef GM_COUNT_PARTICLES_IN_CONE
       MyAtLeastDouble Ekin_total;
       int nnorm;
#endif
       MyIDType Nnp;
       MyAtLeastDouble NnpD;
#ifdef WINDS
       MyAtLeastDouble DelayTime0;
#endif
#ifdef LT_STELLAREVOLUTION
       double Mgas0;
#endif
       MyAtLeastDouble GradDens[3];
       MyAtLeastDouble E_out;
       MyAtLeastDouble E_rec;
       MyAtLeastDouble E_kin, xkin, ykin, zkin;
       MyAtLeastDouble t_startMP;
       int NoArtCond;
       MyAtLeastDouble AveS;
       MyAtLeastDouble AveT;
#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
       double CouplingBHEnHot;
       double CouplingBHEnCold;
#endif
#endif

#ifdef EB_SFR_MAGNETIC
       MyFloat DensityThreshold;
#endif

     }
 *__restrict__ SphP,		/*!< holds SPH particle data on local processor */
*__restrict__ DomainSphBuf;	/*!< buffer for SPH particle data in domain decomposition */


#if defined(AR_PARALLEL_LOCAL_REFINE) || defined(AR_NODE_RECURSIVE_PARALLEL)
     extern int NumUsedThreads;	//when parallelise recursive region we check/set the number of free threads with atomic-capture
     extern int NumberOfAllocatedNodeLocks;
#endif


#ifdef AR_SOA_LITE
     extern struct sph_soa
     {
       short int *wakeup;	/*!< flag to wake up particle */
     } SphPSoa;
//#define AR_ALLOC(name, typ, size) {int bytes; name = (typ*)mymalloc(#name, bytes = size*sizeof(typ)); PANIC_IF(name==NULL, "Failed to allocate `" #name  "`");bytes;}
#define AR_ALLOC(name, typ, size, counter) {int bytes; name = (typ*)mymalloc(#name, bytes = size*sizeof(typ)); counter+=bytes; bzero((void*)name, size*sizeof(typ));}
#define SPHP(i, prop) (SphPSoa.prop[i])
#else
#define SPHP(i, prop) (SphP[i].prop)
#endif

     extern peanokey *DomainKeyBuf;

/* global state of system
*/
     extern struct state_of_system
     {
       double Mass, EnergyKin, EnergyPot, EnergyInt,
#if defined(MAGNETIC_STATISTICS)
	 EnergyMag,
#endif
       EnergyTot,
	 Momentum[4],
	 AngMomentum[4],
	 CenterOfMass[4],
	 MassComp[6],
	 EnergyKinComp[6],
	 EnergyPotComp[6],
	 EnergyIntComp[6],
	 EnergyTotComp[6], MomentumComp[6][4], AngMomentumComp[6][4], CenterOfMassComp[6][4];
#if defined(MAGNETIC_STATISTICS) && defined(TRACEDIVB)
       double DivBerr;
#endif
     }
SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

     extern struct data_index
     {
       int Task;
       int Index;
       int IndexGet;
     }
 *DataIndexTable;		/*!< the particles to be exported are grouped
				   by task-number. This table allows the
				   results to be disentangled again and to be
				   assigned to the correct particle */

     extern struct data_nodelist
     {
       int NodeList[NODELISTLENGTH];
     }
 *DataNodeList, *DataNodeListIn,	/* ancestral modules (as FoF) will keep their DataNodeList in their exchange data structres */
*DataNodeListGet;		/* new modules (where nodelist is separated), for readability purposes will use DataNodeListIn and DataNodeListGet */

#if defined(fSIDM) || defined(rSIDM)
     extern struct data_nodelist_msidm
     {
       int NodeList[NODELISTLENGTH_mSIDM];
     }
 *DataNodeList_mSIDM, *DataNodeListIn_mSIDM,	/* ancestral modules (as FoF) will keep their DataNodeList in their exchange data structres */
*DataNodeListGet_mSIDM;		/* new modules (where nodelist is separated), for readability purposes will use DataNodeListIn and DataNodeListGet */
#endif

     extern struct gravdata_in
     {
       MyLongDouble Pos[3];
       int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
       MyFloat Soft;
#endif
       MyFloat OldAcc;
       int NodeList[NODELISTLENGTH];
#ifdef ADAPTGRAVSOFT
       MyFloat AGS_zeta, AGS_omega, AGS_Hsml, Mass;
#endif
     }
 *GravDataIn,			/*!< holds particle data to be exported to other processors */
*GravDataGet;			/*!< holds particle data imported from other processors */


     extern struct gravdata_out
     {
       MyLongDouble Acc[3];
#ifdef EVALPOTENTIAL
       MyLongDouble Potential;
#endif
#ifdef DISTORTIONTENSORPS
       MyLongDouble tidal_tensorps[3][3];
#endif
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
       MyFloat AGS_corr;
#endif
     }
 *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
*GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


     extern struct potdata_out
     {
       MyLongDouble Potential;
     }
 *PotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
*PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

     extern struct phidotdata_out
     {
       MyDouble PhiDot;
     }
 *PhiDotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
*PhiDotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


     extern struct info_block
     {
       char label[4];
       char type[8];
       int ndim;
       int is_present[6];
     }
 *InfoBlock;


/*! Header for the standard file format.
 */
     extern struct io_header
     {
       int npart[6];		/*!< number of particles of each type in this file */
       double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
       double time;		/*!< time of snapshot file */
       double redshift;		/*!< redshift of snapshot file */
       int flag_sfr;		/*!< flags whether the simulation was including star formation */
       int flag_feedback;	/*!< flags whether feedback was included (obsolete) */
       unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
					   different from npart if one is dealing with a multi-file snapshot. */
       int flag_cooling;	/*!< flags whether cooling was included  */
       int num_files;		/*!< number of files in multi-file snapshot */
       double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
       double Omega0;		/*!< matter density in units of critical density */
       double OmegaLambda;	/*!< cosmological constant parameter */
       double HubbleParam;	/*!< Hubble parameter in units of 100 km/sec/Mpc */
       int flag_stellarage;	/*!< flags whether the file contains formation times of star particles */
       int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
				   particles */
       unsigned int npartTotalHighWord[6];	/*!< High word of the total number of particles of each type */
       int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
       int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */

       int flag_ic_info;	/*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
				   or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
				   For snapshots files, the value informs whether the simulation was evolved from
				   Zeldoch or 2lpt ICs. Encoding is as follows:
				   FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
				   FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
				   FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
				   FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
				   FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
				   All other values, including 0 are interpreted as "don't know" for backwards compatability.
				 */
       float lpt_scalingfactor;	/*!< scaling factor for 2lpt initial conditions */

       char fill[18];		/*!< fills to 256 Bytes */

       char names[15][2];
     }
header;				/*!< holds header for snapshot files */


     enum iofields
     { IO_POS,
       IO_VEL,
       IO_ID,
       IO_MASS,
       IO_SECONDORDERMASS,
       IO_U,
       IO_RHO,
       IO_NE,
       IO_NH,
       IO_HSML,
       IO_SFR,
       IO_AGE,
       IO_HSMS,
       IO_ACRS,
       IO_Z,
       IO_BHMASS,
       IO_BHMDOT,
       IO_BHPROGS,
       IO_BHMBUB,
       IO_BHMINI,
       IO_BHMRAD,
       IO_ACRB,
       IO_POT,
       IO_IPOT,
       IO_PHIDOT,
       IO_ACCEL,
       IO_CR_C0,
       IO_CR_Q0,
       IO_CR_P0,
       IO_CR_E0,
       IO_CR_n0,
       IO_CR_ThermalizationTime,
       IO_CR_DissipationTime,
       IO_HII,
       IO_HeI,
       IO_HeII,
       IO_HeIII,
       IO_H2I,
       IO_H2II,
       IO_HM,
       IO_HD,
       IO_DI,
       IO_DII,
       IO_HeHII,
       IO_DTENTR,
       IO_STRESSDIAG,
       IO_STRESSOFFDIAG,
       IO_STRESSBULK,
       IO_SHEARCOEFF,
       IO_TSTP,
       IO_BFLD,
       IO_BSMTH,
       IO_DBDT,
       IO_DIVB,
       IO_RELDIVB,
       IO_ABVC,
       IO_ACVC,
       IO_AMDC,
       IO_VTURB,
       IO_LTURB,
       IO_ALFA2_DYN,
       IO_ETA2_DYN,
       IO_PHI,
       IO_XPHI,
       IO_GRADPHI,
       IO_ROTB,
       IO_SROTB,
       IO_COOLRATE,
       IO_CONDRATE,
       IO_DENN,
       IO_EGYPROM,
       IO_EGYCOLD,
       IO_MAXSIGNALVEL,
       IO_MACH,
       IO_SHSPEED,
       IO_SHCOMPRESS,
       IO_SHNORMAL,
       IO_SHRHOUP,
       IO_SHPRESUP,
       IO_SHPRESDOWN,
       IO_SHVELUP,
       IO_SHVELDOWN,
       IO_MACH_ALFVEN,
       IO_SHALFVENUP,
       IO_SHALFVENDOWN,
       IO_SHOBLIQUITY,
       IO_DTENERGY,
       IO_PRESHOCK_CSND,
       IO_PRESHOCK_DENSITY,
       IO_PRESHOCK_ENERGY,
       IO_PRESHOCK_XCR,
       IO_DENSITY_JUMP,
       IO_ENERGY_JUMP,
       IO_CRINJECT,
       IO_TIDALTENSORPS,
       IO_DISTORTIONTENSORPS,
       IO_EULERA,
       IO_EULERB,
       IO_FLOW_DETERMINANT,
       IO_PHASE_SPACE_DETERMINANT,
       IO_ANNIHILATION_RADIATION,
       IO_STREAM_DENSITY,
       IO_EOSTEMP,
       IO_EOSXNUC,
       IO_PRESSURE,
       IO_RADGAMMA,
       IO_RAD_ACCEL,
       IO_EDDINGTON_TENSOR,
       IO_LAST_CAUSTIC,
       IO_SHEET_ORIENTATION,
       IO_INIT_DENSITY,
       IO_CAUSTIC_COUNTER,
       IO_DMHSML,		/* for 'SUBFIND_RESHUFFLE_CATALOGUE' option */
       IO_DMDENSITY,
       IO_DMVELDISP,
       IO_VRMS,
       IO_VBULK,
       IO_VRAD,
       IO_VTAN,
       IO_TRUENGB,
       IO_WNGB,
       IO_VDIV,
       IO_VROT,
       IO_VORT,
       IO_DPP,
       IO_LMBCR_pNORM,
       IO_LMBCR_eNORM,
       IO_LMBCR_pSLOPE,
       IO_LMBCR_eSLOPE,
       IO_LMBCR_pCUT,
       IO_LMBCR_eCUT,
       IO_LMBCR_ePRESSURE,
       IO_LMBCR_pPRESSURE,
       IO_LMBCR_eSYNCHROTRON,
       IO_RHO_OLD,
       IO_FB_M_H,
       IO_FB_M_C,
       IO_FB_M_MO,
       IO_FB_E_H,
       IO_FB_M_SF,
       IO_FB_MF,
       IO_FB_NMF,
       IO_FB_EOUT,
       IO_FB_EREC,
       IO_FB_CLOCK,
       IO_FB_E_TOT_0,
       IO_FB_TDYN,
       IO_FB_GRADRHO,
       IO_FB_TCOOL,
       IO_FB_EKINREC,
       IO_FB_TSTARTMP,

       IO_iMass,
       IO_Zs,
       IO_ContribII,
       IO_ContribIa,
       IO_ContribAGB,
       IO_ZAGE,
       IO_ZAGE_LLV,
       IO_CLDX,
       IO_HTEMP,
       IO_TEMP,
       IO_CONTRIB,
       IO_ZSMOOTH,
       IO_allZSMOOTH,
       IO_CHEM,
       IO_DELAYTIME,

       IO_ABUNDANCE,
       IO_ABND_GRAD,
       IO_ABND_SPH,
       IO_DIFFUSING_CB,
       IO_DIFFUSING_CA,
       IO_DIFFUSING_CD,
       IO_DIFF_TSTEP,

       IO_PSUM,
       IO_SIDMNUMNGB,
       IO_NUMTOTALSCATTER,
       IO_SIDMHSML,
       IO_SIDMDENSITY,
       IO_SIDMVELDISP,

       IO_AGS_SOFT,
       IO_AGS_DENS,
       IO_AGS_ZETA,
       IO_AGS_OMEGA,
       IO_AGS_CORR,
       IO_AGS_NGBS,

       IO_MG_PHI,
       IO_MG_GRAD_PHI,
       IO_MG_ACCEL,

       IO_GRADENERGY,
       IO_GRADPRESSURE,

       IO_QP,
       IO_DQP,

       IO_MFM_GRADIENTS_VX,
       IO_MFM_GRADIENTS_VY,
       IO_MFM_GRADIENTS_VZ,
       IO_MFM_GRADIENTS_RHO,
       IO_MFM_GRADIENTS_U,

       /* GM if GM_STARDENSITY  is on, records stellar density, hsml and grav potential if EVALPOTENTIAL is also on */
       IO_RHOSTAR,
       IO_HSMLSTAR,		/* used if GM_STARDENSITY_FIXAPERTURE is NOT defined */
       IO_NEIGHSTAR,		/* used if GM_STARDENSITY_FIXAPERTURE IS defined */

       IO_MAINHALO,		/* Rank in the last subfind output */

       IO_DUSTL, /* used for GL_CR_DUST */
       IO_DUSTS,

       IO_LASTENTRY		/* This should be kept - it signals the end of the list */
     };

     enum siofields
     { SIO_GLEN,
       SIO_GOFF,
       SIO_MTOT,
       SIO_GPOS,
       SIO_MMEA,
       SIO_RMEA,
       SIO_MCRI,
       SIO_RCRI,
       SIO_MTOP,
       SIO_RTOP,
       SIO_M500,
       SIO_R500,
       SIO_M5CC,
       SIO_R5CC,
       SIO_M25K,
       SIO_R25K,
       SIO_DMEA,
       SIO_DCRI,
       SIO_DTOP,
       SIO_D500,
       SIO_D5CC,
       SIO_D25K,
       SIO_MGAS,
       SIO_MSTR,
       SIO_TGAS,
       SIO_LGAS,
       SIO_YGAS,
       SIO_NCON,
       SIO_MCON,
       SIO_BGPOS,
       SIO_BGMTOP,
       SIO_BGRTOP,
       SIO_NSUB,
       SIO_FSUB,
       SIO_SLEN,
       SIO_SOFF,
       SIO_PFOF,
       SIO_MSUB,
       SIO_SPOS,
       SIO_SVEL,
       SIO_SCM,
       SIO_SPIN,
       SIO_DSUB,
       SIO_VMAX,
       SIO_RVMAX,
       SIO_RHMS,
       SIO_SHMR,
       SIO_MBID,
       SIO_GRNR,
       SIO_SMST,
       SIO_SLUM,
       SIO_SLATT,
       SIO_SLOBS,
       SIO_DUST,
       SIO_SAGE,
       SIO_SZ,
       SIO_SSFR,
       SIO_SSYNCHROTRON,
       SIO_MHI,
       SIO_SSANG,
       SIO_SBVAL,
       SIO_SZGAS,
       SIO_SZCGAS,
       SIO_PPOS,
       SIO_PVEL,
       SIO_PTYP,
       SIO_PMAS,
       SIO_PAGE,
       SIO_PID,

       SIO_LASTENTRY
     };

/*
 * Variables for Tree
 * ------------------
 */

     extern int Nexport, Nimport;
     extern int BufferFullFlag;
     extern int NextParticle;
     extern int NextGroup;
     extern int NextJ;
     extern int TimerFlag;

// note, the ALIGN(32) directive will effectively pad the structure size
// to a multiple of 32 bytes
     extern ALIGN(32)
     struct NODE
     {
       MyFloat center[3];	/*!< geometrical center of node */
       MyFloat len;		/*!< sidelength of treenode */
#ifdef AR_XMAS_TREE
       unsigned char allocated;
#endif

       union
       {
	 int suns[8];		/*!< temporary pointers to daughter nodes */
	 struct
	 {
	   MyFloat s[3];	/*!< center of mass of node */
	   MyFloat mass;	/*!< mass of node */
	   int sibling;		/*!< this gives the next node in the walk in case the current node can be used */
	   int nextnode;	/*!< this gives the next node in case the current node needs to be opened */
	   int father;		/*!< this gives the parent node of each node (or -1 if we have the root node) */
	   unsigned int bitflags;	/*!< flags certain node properties */
	 }
	 d;
       }
       u;

       integertime Ti_current;
       double GravCost;


#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTGRAVSOFT) || defined(fSIDM) || defined(rSIDM)
       MyFloat maxsoft;		/*!< hold the maximum gravitational softening of particle in the
				   node if either the ADAPTIVE_GRAVSOFT_FORGAS or the ADAPTGRAVSOFT option is selected */
#endif

#if defined(ADAPTGRAVSOFT) || defined(fSIDM) || defined(rSIDM)
       int Flag;
#endif
     }
 *Nodes_base,			/*!< points to the actual memory allocted for the nodes */
*Nodes;				/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
				   gives the first allocated node */
#ifdef AR_TREEBUILD_PARALLEL
     extern omp_lock_t *NodesLockList;	//locks for the tree build
#endif


#ifdef ACC_FORCETREE
/*! auxialiary variable used to set-up non-recursive walk */
//extern int last;
#define NTAB 1000
/*! variables for short-range lookup table */
     extern float shortrange_table[NTAB], shortrange_table_potential[NTAB];
#ifdef DISTORTIONTENSORPS
     extern float shortrange_table_tidal[NTAB];
#endif
     extern double fac_intp;
#define EN 64
     extern MyFloat potcorr[EN + 1][EN + 1][EN + 1];
     extern MyFloat fcorrx[EN + 1][EN + 1][EN + 1];
     extern MyFloat fcorry[EN + 1][EN + 1][EN + 1];
     extern MyFloat fcorrz[EN + 1][EN + 1][EN + 1];

#endif

     extern struct extNODE
     {
       MyLongDouble dp[3];
#ifdef GRAVITY_CENTROID
       int suns[8];
#endif
       MyFloat vs[3];
       MyFloat vmax;
       MyFloat hmax;		/*!< maximum SPH smoothing length in node. Only used for gas particles */
       integertime Ti_lastkicked;
       int Flag;

     }
 *Extnodes, *Extnodes_base;


     extern int MaxNodes;	/*!< maximum allowed number of internal nodes */
     extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


     extern int *Nextnode;	/*!< gives next node in tree walk  (nodes array) */
#ifdef GROUP_LEAVES
     extern int *NextGroupedNode;	/*!< gives next node in tree walk  (nodes array) */
#endif

     extern int *Father;	/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
     extern double Rs, R200;
     extern double Dc;
     extern double RhoCrit, V200;
     extern double fac;
#endif

#if defined  (UM_METAL_COOLING)
/* --==[ link with LT_ stuffs]==-- */
     extern float *um_ZsPoint, um_FillEl_mu, um_mass;
#endif


#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
/* ----- chemistry part ------- */

 // UM: H number fraction is 0.93, mass fraction is 0.76  
#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
     extern double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
     extern double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T],
       k10a[N_T], k11a[N_T];
     extern double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T],
       k20a[N_T], k21a[N_T];
     extern double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
     extern double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
     extern double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
       sigma30[N_nu], sigma31[N_nu];
#endif
#endif


#ifdef UM_CHEMISTRY
     extern double kHeHII1a[N_T], kHeHII2a[N_T], kHeHII3a[N_T];
     extern double kHD1a[N_T], kHD2a[N_T], kHD3a[N_T], kHD4a[N_T], kHD5a[N_T], kHD6a[N_T], kHD7a[N_T];
#endif

#ifdef UM_CHEMISTRY
     extern double gJH0uvb;
     extern double gJHe0uvb;
     extern double gJHepuvb;
     extern double epsH0uvb;
     extern double epsHe0uvb;
     extern double epsHepuvb;

     extern double HIshielding;
#endif

#if defined (UM_CHEMISTRY) && defined (CHEMISTRY)
#error you cannot define both UM_CHEMISTRY and CHEMISTRY !
#endif


#ifdef UM_CHEMISTRY
#define T_SUP_INTERPOL_LIMIT        1.e4
#endif


#ifdef SCFPOTENTIAL
     extern long scf_seed;
     extern MyDouble *Anltilde, *coeflm, *twoalpha, *c1, *c2, *c3;
     extern MyDouble *cosmphi, *sinmphi;
     extern MyDouble *ultrasp, *ultraspt, *ultrasp1;
     extern MyDouble *dblfact, *plm, *dplm;
     extern MyDouble *sinsum, *cossum;
     extern MyDouble *sinsum_all, *cossum_all;
#ifdef SCF_SCALEFAC
     extern float *scalefac_nlm;
#endif
#endif


     extern int maxThreads;

#ifdef KD_ALTERNATIVE_GROUP_SORT
     extern int MaxNgroups;
#endif


     extern __thread int mythreadid;
     extern __thread integertime last_time0_drift, last_time1_drift,
       last_time0_hydrokick, last_time1_hydrokick, last_time0_gravkick, last_time1_gravkick;
     extern __thread double last_value_drift, last_value_hydrokick, last_value_gravkick;
#ifdef MAGNETIC
     extern __thread integertime last_time0_magkick, last_time1_magkick;
     extern __thread double last_value_magkick;
#endif
#if defined(fSIDM) || defined(rSIDM)
     extern __thread integertime last_time0_mSIDMkick, last_time1_mSIDMkick;
     extern __thread double last_value_mSIDMkick;
#endif

//#ifdef _OPENMP
//#pragma omp threadprivate(last_time0_drift, last_time1_drift, last_value_drift, last_time0_hydrokick, last_time1_hydrokick, last_value_hydrokick, last_time0_gravkick, last_time1_gravkick, last_value_gravkick)
//#ifdef MAGNETIC
//#pragma omp threadprivate(last_time0_magkick, last_time1_magkick, last_value_magkick)
//#endif
//#endif

#ifdef GADGET3_IO_LIB
     extern int all_particles_blocks[6 * 1000]; // 6 particle types, max. 1000 block types
#endif

#endif /* ALLVARS_H  - please do not put anything below this line */

#ifdef PMGRID
#ifdef KD_FFTW_MEMORY_BALANCE
     extern int dofftw;
#endif
#endif


#ifdef AA_RELAX_ATMOSPHERE
     extern int relaxCounter;
#endif


#ifdef USE_NVTX
#include "nvToolsExt.h"

     extern const uint32_t colors[];	//= { 0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff, 0x0000ffff, 0x00ff0000, 0x00ffffff };
     extern const int num_colors;	// = sizeof(colors)/sizeof(uint32_t);
//we generate a color according to the first letter of the name
/* nvtxEventAttributes_t eventAttrib = {0}; \#define PUSH_T(name) do { \ */
// we force the user to add a semicolon after!

#define PUSH_T(name) do {     int color_id = (name[2]*4+name[1]*2+name[0])%num_colors;     nvtxEventAttributes_t eventAttrib = {0};     eventAttrib.version = NVTX_VERSION;     eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;     eventAttrib.colorType = NVTX_COLOR_ARGB;     eventAttrib.color = colors[color_id];     eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;     eventAttrib.message.ascii = name;     nvtxRangePushEx(&eventAttrib); } while(0)
#define POP_T() do{nvtxRangePop();}while(0)
#else
#define PUSH_T(name) do{}while(0)
#define POP_T() do{}while(0)
#endif

#ifdef ACC_GRAVITY
#define PRAGMA_GRAVITY_DATA_REGION() _Pragma("acc data				      \
  copy(NodeGravCosts[0:NodesSize], PGravCosts[0:NumPart], LitePOut[0:NActivePart]) \
         if(acc_all.go_gpu)")

#else
#define PRAGMA_GRAVITY_DATA_REGION()
#endif

#if defined(ACC) && !defined(ACC_DECLARED)
//since allvars.h can be included in many places multiple times, it happens that we 
//redeclare openacc variables, with this ifdef we check that we run #pragma acc declare only once
#define ACC_DECLARED
     extern int AccErrorCode;
     extern int NodesSize;
#pragma acc declare create  (All,  TakeLevel,  MaxNodes,  Numnodestree, NActivePart, AccErrorCode, NumPart, ActiveParticleList, NodesSize)

#if defined( LT_METAL_COOLING)|| defined(LT_METAL_COOLING_WAL)
     extern int WalCool_n[6];
#define WalCool_n_Rho WalCool_n[0]	/*!< number of density bins      */
#define WalCool_n_Ab  WalCool_n[4]	/*!< number of solar abundances  */
     extern int *SpeciesIdx;
     extern int *SpeciesPos;
     extern int myHyd, myHel, HydPos, HelPos;
     extern int UseHeNumberRatio;
#pragma acc declare create(SFs, CoolZvalue, ZBins)
#pragma  acc declare create (WalCool_n, UseHeNumberRatio, myHyd,myHel,HydPos,HelPos, SpeciesIdx,SpeciesPos)
#endif

#ifdef PERIODIC
#pragma acc declare create (boxHalf, boxSize,fac_intp, fcorrx, fcorry, fcorrz, potcorr)
#ifdef LONG_X
#pragma acc declare create  ( boxHalf_X, boxSize_X)
#endif
#ifdef LONG_Y
#pragma acc declare create  ( boxHalf_Y, boxSize_Y)
#endif
#ifdef LONG_Z
#pragma acc declare create  ( boxHalf_Z, boxSize_Z)
#endif
#endif
#pragma acc declare create(Nodes_base)
#pragma acc declare create(Nextnode)
#pragma acc declare create (Extnodes_base)
#pragma acc declare create (shortrange_table, shortrange_table_potential)
#ifdef DISTORTIONTENSORPS
#pragma acc declare create(shortrange_table_tidal)
#endif



#endif
