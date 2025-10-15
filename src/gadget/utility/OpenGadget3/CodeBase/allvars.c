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

#ifdef GADGET3_IO_LIB
#include <mpi.h>
#endif

#include "allvars.h"



#ifdef PERIODIC
MyDouble boxSize, boxHalf, inverse_boxSize;

#ifdef LONG_X
MyDouble boxSize_X, boxHalf_X, inverse_boxSize_X;
#else
#endif
#ifdef LONG_Y
MyDouble boxSize_Y, boxHalf_Y, inverse_boxSize_Y;
#else
#endif
#ifdef LONG_Z
MyDouble boxSize_Z, boxHalf_Z, inverse_boxSize_Z;
#else
#endif
#endif

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
MPI_Status mpistat;
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

#ifdef KD_MONITOR_PERFORMANCE
int count_performance = 0;
#endif

__thread int mythreadid;
__thread integertime last_time0_drift, last_time1_drift,
  last_time0_hydrokick, last_time1_hydrokick, last_time0_gravkick, last_time1_gravkick;
__thread double last_value_drift, last_value_hydrokick, last_value_gravkick;
#ifdef MAGNETIC
__thread integertime last_time0_magkick, last_time1_magkick;
__thread double last_value_magkick;
#endif

#if defined(fSIDM) || defined(rSIDM)
__thread integertime last_time0_mSIDMkick, last_time1_mSIDMkick;
__thread double last_value_mSIDMkick;
#endif

#ifdef AXION_DM
__thread integertime last_time0_axkick, last_time1_axkick;
__thread double last_value_axkick;
#endif

__thread int ThisThread;	/* the Number of the local thread */
int maxThreads = 1;		/* number of threads */
int ThisTask;			/*!< the number of the local processor  */
int NTask;			/*!< number of processors */
int PTask;			/*!< note: NTask = 2^PTask */

#if !defined(GADGET3_IO_LIB) || defined(GADGET3_IO_LIB_MPI)
MPI_Comm MYMPI_COMM_WORLD;
#endif

double CPUThisRun;		/*!< Sums CPU time of current process */

#ifdef RECOMPOSE_DOMAIN_FACTOR
int LastDomainUpdate = 0;
#endif

int NumForceUpdate;		/*!< number of active particles on local processor in current timestep  */
long long GlobNumForceUpdate;
int NumSphUpdate;		/*!< number of active SPH particles on local processor in current timestep  */

int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
int RestartSnapNum;

double PostProcessDt;

int SelRnd;

int *Exportflag;		/*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount;

#ifdef IMPORT_ALLtoALLv_FROM_G4
MPI_Win win;
#endif

int TakeLevel;

#ifdef _OPENMP
int NActivePart;
int *ActiveParticleList;
#if defined(LT_STELLAREVOLUTION)
int NActivePartMerk;
#endif
#endif

int FirstActiveParticle;
int *NextActiveParticle;
unsigned char *ProcessedFlag;

int TimeBinCount[TIMEBINS];
int TimeBinCountSph[TIMEBINS];
int TimeBinActive[TIMEBINS];

int FirstInTimeBin[TIMEBINS];
int LastInTimeBin[TIMEBINS];
int *NextInTimeBin;
int *PrevInTimeBin;

#ifdef AXION_DM
int TimeBinCountAx[TIMEBINS];
#endif

size_t HighMark_run, HighMark_domain, HighMark_gravtree,
  HighMark_pmperiodic, HighMark_pmnonperiodic, HighMark_sphdensity, HighMark_sphhydro, HighMark_addSPH;

#ifdef GM_MUPPI
size_t HighMark_mppin;
#endif


#ifdef ADAPTGRAVSOFT
size_t HighMark_agsdensity;
#endif

#if defined(fSIDM) || defined(rSIDM)
size_t HighMark_msidm;
#endif

#ifdef SFR
double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
double TimeBin_BH_mass[TIMEBINS];
double TimeBin_BH_dynamicalmass[TIMEBINS];
double TimeBin_BH_Mdot[TIMEBINS];
double TimeBin_BH_Medd[TIMEBINS];
#endif



char DumpFlag = 1;

size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;

double CPU_Step[CPU_PARTS];
char CPU_Symbol[CPU_PARTS] =
  { '-', '*', '=', ';', '<', '[', '^', ':', '.', '~', '|', '+', '"', '/', '`', ',', '>', '@', '#', '&', '$',
  ']', '(', '?', ')', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '\\', '%', '{', '}'
};

char CPU_SymbolImbalance[CPU_PARTS] =
  { 'a', 't', 'u', 'v', 'b', 'w', 'd', 'r', 'h', 'm', 'n', 'l', 'o', 'p', 's', 'f', 'i', 'g', 'c', 'e', 'x',
  'y', 'z', 'A', 'I', 'W', 'T', 'V', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'
};

char CPU_String[CPU_STRING_LEN + 1];

double WallclockTime;		/*!< This holds the last wallclock time measurement for timings measurements */

int Flag_FullStep;		/*!< Flag used to signal that the current step involves all particles */


int TreeReconstructFlag;
int GlobFlag;

int NumPart;			/*!< number of particles on the LOCAL processor */

int N_gas;			/*!< number of gas particles on the LOCAL processor  */
int N_stars;			/*!< number of star particles in the LOCAL processor */
int N_BHs;

long long Ntype[6];		/*!< total number of particles of each type */
int NtypeLocal[6];		/*!< local number of particles of each type */

#ifndef GADGET3_IO_LIB
gsl_rng *random_generator;	/*!< the random number generator used */
#endif


#ifdef SFR
int Stars_converted;		/*!< current number of star particles in gas particle block */
#endif

#ifdef OUTPUT_LIGHTCONES
int lc_step_data_present;
int lc_step_data_len;
double lc_offset, lc_addistmin;

struct output_data *step_data;
#endif

double TimeOfLastTreeConstruction;	/*!< holds what it says */

int *Ngblist;			/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */

int *RecyclingBuffer;
int HasFisrtParticleBeenProcessed;

double *R2ngblist;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int *DomainStartList, *DomainEndList;



double *DomainWork;
int *DomainCount;
int *DomainCountSph;
int *DomainTask;
int *DomainNodeIndex;
int *DomainList, DomainNumChanged;

peanokey *__restrict__ Key, *__restrict__ KeySorted;

struct topnode_data *TopNodes;


int NTopnodes, NTopleaves;


double RndTable[RNDTABLE];

#ifdef SUBFIND
int GrNr;
int NumPartGroup;
#endif


#ifdef WRITE_KEY_FILES
peanokey *KeyIndex;
int *NPartPerKey, *PartKeyOffset;
int NKeys[6];
long long NKeysTotal[6];
#endif

#ifdef LT_ADD_GAL_TO_SUB
float *tempiAS, *CB07, *Filters_Effective_Wavelenght;
#ifdef OBSERVER_FRAME
float *CB07obs;
#ifdef INTERP_OBSERVER_FRAME
float CB07obs_redshifts[INTERP_OBSERVER_FRAME];
int flag_readCB07obs_redshifts;
#endif
#endif
#endif


/* variables for input/output , usually only used on process 0 */


char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,			/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU,			/*!< file handle for cpu.txt log-file. */
 *FdTimebin, *FdParamChangeLog;

#ifdef SCFPOTENTIAL
FILE *FdSCF;
#endif

#ifdef SFR
FILE *FdSfr;			/*!< file handle for sfr.txt log-file. */

#ifdef GM_MUPPI

FILE *FdExit;
#ifdef GM_COUNT_PARTICLES_IN_CONE
FILE *FdCone;
#endif
#ifdef MV_GM_STELLAR_KIN_FB2_OUTPUT
FILE *FdWind;
#endif
#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION_OUTPUT
FILE *FdMolFrac;
#endif
#ifdef MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z_OUTPUT
FILE *FdHypernovae;
#endif
#ifdef MV_GM_AGNMUPPI
FILE *FdAGNmuppi;
#endif
#ifdef MV_GM_AGNMUPPI_OUTPUT
FILE *FdEgyAGN;
FILE *FdEgyAGNtot;
#endif
#if defined(MV_GM_AGNMUPPI) && defined(MV_GM_COVER_FACT_MCLS)
FILE *FdCovFacAGN;
#endif
#ifdef MV_AGNMUPPI_COOLING_OFF
FILE *FdCoolOffAGNmuppi;
#endif

#endif // closes GM_MUPPI

#endif //closes SFR


#ifdef DISTORTIONTENSORPS
#ifdef PMGRID
FILE *FdTidaltensor;		/*!< file handle for Tidaltensor.txt log-file. */
#endif
#endif

#ifdef BLACK_HOLES
FILE *FdBlackHoles;		/*!< file handle for blackholes.txt log-file. */
FILE *FdBlackHolesDetails;
#endif

#ifdef LMB_SPECTRAL_CRs
FILE *FdCRs;			/*!< file handle for CRs.txt log-file. */
FILE *FdCRsDetailsAdiabatic;
FILE *FdCRsDetailsRadiative;
FILE *FdCRsDetailsDpp;
FILE *FdCRsDetailsShockInjection;
FILE *FdCRsDetailsSfrInjection;
#endif

#ifdef FORCETEST
FILE *FdForceTest;		/*!< file handle for forcetest.txt log-file. */
#endif

#ifdef DARKENERGY
FILE *FdDE;			/*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
FILE *FdXXL;			/*!< file handle for xxl.txt log-file. */

#ifdef MAGNETIC
double MeanB;

#ifdef TRACEDIVB
double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
double MeanAlpha;
#endif
#endif

/*! table for the cosmological drift factors */
double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
double HydroKickTable[DRIFT_TABLE_LENGTH];

#ifdef AXION_DM
/*! table for the cosmological kick factor for axion hydrodynmical forces */
double AxHydroKickTable[DRIFT_TABLE_LENGTH];
#endif

#ifdef MAGNETIC
/*! table for the cosmological kick factor for induction equation */
double MagKickTable[DRIFT_TABLE_LENGTH];
#endif

#if defined(fSIDM) || defined(rSIDM)
/*! table for the cosmological kick factor for mSIDM */
double MSIDMKickTable[DRIFT_TABLE_LENGTH];
#endif


void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;


/*! This structure contains data which is the SAME for all tasks as the All data structure
 * but its values are computed on the fly thus it can change in shape betwen restarts
 */

#ifndef GADGET3_IO_LIB
struct derived_global_data_all_processes AllDerived;
#endif

/*! Structure for global variables needed for frequent self-interacting dark matter */
#if defined(fSIDM) || defined(rSIDM)
LambdaNumerical calc_lambda;
GaussLegendre gauss_legendre;
MSIDM_Scattering msidm_scatter;
#endif


/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *__restrict__ P,	/*!< holds particle data on local processor */
 *__restrict__ DomainPartBuf;	/*!< buffer for particle data used in domain decomposition */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *__restrict__ SphP,	/*!< holds SPH particle data on local processor */
 *__restrict__ DomainSphBuf;	/*!< buffer for SPH particle data in domain decomposition */

#if defined(AR_PARALLEL_LOCAL_REFINE) || defined(AR_NODE_RECURSIVE_PARALLEL)
int NumUsedThreads;		//when parallelise recursive region we check/set the number of free threads with atomic-capture
int NumberOfAllocatedNodeLocks;
#endif

#ifdef AR_SOA_LITE
sph_soa SphPSoa;
#endif


#ifdef AXION_DM
struct ax_particle_data *AxP,	/*!< holds AX particle data on local processor */
 *__restrict__ DomainAxBuf;	/*!< buffer for AX particle data in domain decomposition */
#endif

#ifdef LT_STELLAREVOLUTION
struct met_particle_data *__restrict__ MetP, *__restrict__ DomainMetBuf;
#endif

#if defined(BLACK_HOLES)
struct bh_particle_data *__restrict__ BHP;
#endif

peanokey *DomainKeyBuf;

/* global state of system
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

struct data_index *DataIndexTable;	/*!< the particles to be exported are grouped
					   by task-number. This table allows the
					   results to be disentangled again and to be
					   assigned to the correct particle */

struct data_nodelist *DataNodeList, *DataNodeListIn, *DataNodeListGet;

#if defined(fSIDM) || defined(rSIDM)
struct data_nodelist_msidm *DataNodeList_mSIDM, *DataNodeListIn_mSIDM, *DataNodeListGet_mSIDM;
#endif

struct gravdata_in *GravDataIn,	/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


struct gravdata_out *GravDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct potdata_out *PotDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

struct phidotdata_out *PhiDotDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PhiDotDataOut;		/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct info_block *InfoBlock;

/*! Header for the standard file format.
 */
struct io_header header;	/*!< holds header for snapshot files */





/*
 * Variables for Tree
 * ------------------
 */

int Nexport, Nimport;
int BufferFullFlag;
int NextParticle;
int NextGroup;
int NextJ;
int TimerFlag;

struct NODE *Nodes_base,	/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
				   gives the first allocated node */

#ifdef AR_TREEBUILD_PARALLEL
omp_lock_t *NodesLockList;	//locks for the tree build
#endif

struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/*!< maximum allowed number of internal nodes */
int Numnodestree;		/*!< number of (internal) nodes in each tree */


int *Nextnode;			/*!< gives next node in tree walk  (nodes array) */
#ifdef GROUP_LEAVES
int *NextGroupedNode;		/*!< gives next node in tree walk  (nodes array) */
#endif

int *Father;			/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif

#if defined (UM_CHEMISTRY) || defined (UM_METAL_COOLING)
/* --==[ link with LT_ stuffs]==-- */
float *um_ZsPoint, um_FillEl_mu, um_mass;

/* char *PT_Symbols; */
/* double *PT_Masses; */
/* int NPT; */

/* double **II_AvgFillNDens, **IIShLv_AvgFillNDens, **Ia_AvgFillNDens, **AGB_AvgFillNDens; */
#endif

#if defined (CHEMISTRY) || defined (UM_CHEMISTRY)
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T],
  k10a[N_T], k11a[N_T];
double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T],
  k20a[N_T], k21a[N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
  sigma30[N_nu], sigma31[N_nu];
#endif
#endif

#ifdef UM_CHEMISTRY
double gJH0uvb;
double gJHe0uvb;
double gJHepuvb;
double epsH0uvb;
double epsHe0uvb;
double epsHepuvb;

double HIshielding;
#endif

#if defined (UM_CHEMISTRY)
double kHeHII1a[N_T], kHeHII2a[N_T], kHeHII3a[N_T];
double kHD1a[N_T], kHD2a[N_T], kHD3a[N_T], kHD4a[N_T], kHD5a[N_T], kHD6a[N_T], kHD7a[N_T];
#endif

#ifdef LT_STELLAREVOLUTION
#include "../CoolingSfr/Sfr_LT/lt.c"
#endif



#ifdef SCFPOTENTIAL
long scf_seed;
MyDouble *Anltilde, *coeflm, *twoalpha, *c1, *c2, *c3;
MyDouble *cosmphi, *sinmphi;
MyDouble *ultrasp, *ultraspt, *ultrasp1;
MyDouble *dblfact, *plm, *dplm;
MyDouble *sinsum, *cossum;
MyDouble *sinsum_all, *cossum_all;
#ifdef SCF_SCALEFAC
float *scalefac_nlm;
#endif
#endif

#ifdef KD_ALTERNATIVE_GROUP_SORT
int MaxNgroups;
#endif

#ifdef PMGRID
#ifdef KD_FFTW_MEMORY_BALANCE
int dofftw = 0;
#endif
#endif

#ifdef AA_RELAX_ATMOSPHERE
int relaxCounter = 0;
#endif

#ifdef USE_NVTX
#include "nvToolsExt.h"

const uint32_t colors[] = { 0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff,
  0x0000ffff, 0x00ff0000,
  0x00800000, 0x00000080, 0x00008000,
  0x00808080, 0x00808000, 0x00008080, 0x00800080
};


const int num_colors = sizeof(colors) / sizeof(uint32_t);
#endif



#if GADGET_HYDRO == HYDRO_MFM
Eos < NUMDIMS > *eos;
FluxSolver < NUMDIMS > *fluxSolver;
//SlopeLimiter<NUMDIMS> *slopeLimiter;
#endif

#ifdef GADGET3_IO_LIB
int all_particles_blocks[6 * 1000] = {0}; // 6 particle types, max. 1000 block types
#endif