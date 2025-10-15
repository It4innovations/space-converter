/*
 * Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
 *
 * This program is free software : you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "gadget_extract_iolib.h"

#ifdef WITH_NANOVDB

#ifdef WITH_TBB
#   define NANOVDB_USE_TBB
#   define NANOVDB_USE_INTRINSICS
#endif

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/GridBuilder.h>
#	include <nanovdb/util/CreateNanoGrid.h>
#	include <nanovdb/util/IO.h>
#else
#	include <nanovdb/tools/GridBuilder.h>
#	include <nanovdb/tools/CreateNanoGrid.h>
#	include <nanovdb/io/IO.h>
#endif

#endif

#include <float.h>
#include <iostream>

#ifdef WITH_OPENMP
#	include <omp.h>
#endif

#include "gadget3_io_lib.h"
#include "convert_common.h"

#include <mpi.h>

extern "C"
{

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"

}

#include <string.h>

#define RETURN_NORM_EMPTY return 0; //return std::numeric_limits<float>::quiet_NaN();//0;
#define RETURN_NORM_VALUE(v) return (float)(v);
#define RETURN_NORM_VECTOR3(v) return (float)common::calculate_dmagnitude3(v[0], v[1], v[2]);
#define RETURN_NORM_DVECTORN(v,n) return (float)common::calculate_dmagnituden(v,n);
#define RETURN_NORM_FVECTORN(v,n) return (float)common::calculate_fmagnituden(v,n);

#define RETURN_COMP_EMPTY return 0;
#define RETURN_COMP_VALUE(v) return 1;
#define RETURN_COMP_VECTOR3(v) return 3;
#define RETURN_COMP_DVECTORN(v,n) return n;
#define RETURN_COMP_FVECTORN(v,n) return n;

#define RETURN_ORIG_EMPTY RETURN_COMP_EMPTY
#define RETURN_ORIG_VALUE(v) out_value[0] = (float)v; RETURN_COMP_VALUE(v)
#define RETURN_ORIG_VECTOR3(v) out_value[0] = (float)v[0]; out_value[1] = (float)v[1]; out_value[2] = (float)v[2]; RETURN_COMP_VECTOR3(v)
#define RETURN_ORIG_DVECTORN(v,n) for(int iv=0;iv<n;iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_DVECTORN(v,n)
#define RETURN_ORIG_FVECTORN(v,n) for(int iv=0;iv<n;iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_FVECTORN(v,n)

//using namespace nanovdb;

namespace gadget {
	namespace io {

		size_t get_local_num_particles()
		{
			return NumPart;
		}

		size_t get_global_num_particles()
		{
			return All.TotNumPart;
		}

		const char* get_cpufields_name(int field)
		{
			switch (field)
			{
			case CPU_ALL: return "CPU_ALL";
			case CPU_TREEWALK1: return "CPU_TREEWALK1";
			case CPU_TREEWALK2: return "CPU_TREEWALK2";
			case CPU_TREEWAIT1: return "CPU_TREEWAIT1";
			case CPU_TREEWAIT2: return "CPU_TREEWAIT2";
			case CPU_TREESEND: return "CPU_TREESEND";
			case CPU_TREERECV: return "CPU_TREERECV";
			case CPU_TREEMISC: return "CPU_TREEMISC";
			case CPU_TREEBUILD: return "CPU_TREEBUILD";
			case CPU_TREEUPDATE: return "CPU_TREEUPDATE";
			case CPU_TREEHMAXUPDATE: return "CPU_TREEHMAXUPDATE";
			case CPU_DOMAIN: return "CPU_DOMAIN";
			case CPU_DENSCOMPUTE: return "CPU_DENSCOMPUTE";
			case CPU_DENSWAIT: return "CPU_DENSWAIT";
			case CPU_DENSCOMM: return "CPU_DENSCOMM";
			case CPU_DENSMISC: return "CPU_DENSMISC";
			case CPU_HYDCOMPUTE: return "CPU_HYDCOMPUTE";
			case CPU_HYDWAIT: return "CPU_HYDWAIT";
			case CPU_HYDCOMM: return "CPU_HYDCOMM";
			case CPU_HYDMISC: return "CPU_HYDMISC";
			case CPU_DRIFT: return "CPU_DRIFT";
			case CPU_TIMELINE: return "CPU_TIMELINE";
			case CPU_POTENTIAL: return "CPU_POTENTIAL";
			case CPU_MESH: return "CPU_MESH";
			case CPU_PEANO: return "CPU_PEANO";
			case CPU_COOLINGSFR: return "CPU_COOLINGSFR";
			case CPU_SNAPSHOT: return "CPU_SNAPSHOT";
			case CPU_FOF: return "CPU_FOF";
			case CPU_BLACKHOLES: return "CPU_BLACKHOLES";
			case CPU_MISC: return "CPU_MISC";
#ifdef GM_STARDENSITY		
			case CPU_STARDENSMISC: return "CPU_STARDENSMISC";
			case CPU_STARCOMPUTE: return "CPU_STARCOMPUTE";
			case CPU_STARWAIT: return "CPU_STARWAIT";
			case CPU_STARCOMM: return "CPU_STARCOMM";
#endif		
			case CPU_SMTHCOMPUTE: return "CPU_SMTHCOMPUTE";
			case CPU_SMTHWAIT: return "CPU_SMTHWAIT";
			case CPU_SMTHCOMM: return "CPU_SMTHCOMM";
			case CPU_SMTHMISC: return "CPU_SMTHMISC";
			case CPU_HOTNGBS: return "CPU_HOTNGBS";
			case CPU_WEIGHTS_HOT: return "CPU_WEIGHTS_HOT";
			case CPU_ENRICH_HOT: return "CPU_ENRICH_HOT";
			case CPU_WEIGHTS_COLD: return "CPU_WEIGHTS_COLD";
			case CPU_ENRICH_COLD: return "CPU_ENRICH_COLD";
			case CPU_CSMISC: return "CPU_CSMISC";
			case CPU_HYDNETWORK: return "CPU_HYDNETWORK";
			case CPU_AGSDENSCOMPUTE: return "CPU_AGSDENSCOMPUTE";
			case CPU_AGSDENSWAIT: return "CPU_AGSDENSWAIT";
			case CPU_AGSDENSCOMM: return "CPU_AGSDENSCOMM";
			case CPU_AGSDENSMISC: return "CPU_AGSDENSMISC";
			case CPU_AGSTREEHMAXUPD: return "CPU_AGSTREEHMAXUPD";
			case CPU_MG_CIC: return "CPU_MG_CIC";
			case CPU_MG_FIELDSOLVE: return "CPU_MG_FIELDSOLVE";
			case CPU_MG_EFF_MASS: return "CPU_MG_EFF_MASS";
			case CPU_CONDUCTION: return "CPU_CONDUCTION";
#if GADGET_HYDRO == HYDRO_MFM		
			case CPU_MFMGRADCOMPUTE: return "CPU_MFMGRADCOMPUTE";
			case CPU_MFMGRADWAIT: return "CPU_MFMGRADWAIT";
			case CPU_MFMGRADCOMM: return "CPU_MFMGRADCOMM";
			case CPU_MFMGRADMISC: return "CPU_MFMGRADMISC";
			case CPU_MFMLIMCOMPUTE: return "CPU_MFMLIMCOMPUTE";
			case CPU_MFMLIMWAIT: return "CPU_MFMLIMWAIT";
			case CPU_MFMLIMCOMM: return "CPU_MFMLIMCOMM";
			case CPU_MFMLIMMISC: return "CPU_MFMLIMMISC";
			case CPU_MFMFLUXCOMPUTE: return "CPU_MFMFLUXCOMPUTE";
			case CPU_MFMFLUXWAIT: return "CPU_MFMFLUXWAIT";
			case CPU_MFMFLUXCOMM: return "CPU_MFMFLUXCOMM";
			case CPU_MFMFLUXMISC: return "CPU_MFMFLUXMISC";
#endif		
			case CPU_MUPPI: return "CPU_MUPPI";
			case CPU_MUPPINET: return "CPU_MUPPINET";
#if defined(fSIDM) || defined(rSIDM)		
			case CPU_MSIDMCOMPUTE: return "CPU_MSIDMCOMPUTE";
			case CPU_MSIDMPREPP: return "CPU_MSIDMPREPP";
			case CPU_MSIDMCOMM: return "CPU_MSIDMCOMM";
			case CPU_MSIDMWAIT: return "CPU_MSIDMWAIT";
			case CPU_MSIDMNETWORK: return "CPU_MSIDMNETWORK";
			case CPU_MSIDMMISC: return "CPU_MSIDMMISC";
			case CPU_PARTS4SUM: return "CPU_PARTS4SUM";
			case CPU_PARTS4SUMM1: return "CPU_PARTS4SUMM1";
			case CPU_SUBFIND_TREEBUILD_SPECIES: return "CPU_SUBFIND_TREEBUILD_SPECIES";
			case CPU_SUBFIND_TREEBUILD: return "CPU_SUBFIND_TREEBUILD";
			case CPU_SUBFIND_SMOOTHINGLENGTH: return "CPU_SUBFIND_SMOOTHINGLENGTH";
			case CPU_SUBFIND_DENSITY: return "CPU_SUBFIND_DENSITY";
			case CPU_SUBFIND_DMDENSITY: return "CPU_SUBFIND_DMDENSITY";
			case CPU_SUBFIND_SAVE_DENSITY: return "CPU_SUBFIND_SAVE_DENSITY";
			case CPU_SUBFIND_EXCHANGE: return "CPU_SUBFIND_EXCHANGE";
			case CPU_SUBFIND_COLLHALOS: return "CPU_SUBFIND_COLLHALOS";
			case CPU_SUBFIND_SORTLOCAL: return "CPU_SUBFIND_SORTLOCAL";
			case CPU_SUBFIND_LOCALGROUPS: return "CPU_SUBFIND_LOCALGROUPS";
			case CPU_SUBFIND_UNSORTLOCAL: return "CPU_SUBFIND_UNSORTLOCAL";
			case CPU_SUBFIND_EXCHANGERETURN: return "CPU_SUBFIND_EXCHANGERETURN";
			case CPU_SUBFIND_DOMAIN: return "CPU_SUBFIND_DOMAIN";
			case CPU_SUBFIND_MASSES: return "CPU_SUBFIND_MASSES";
			case CPU_SUBFIND_CONT: return "CPU_SUBFIND_CONT";
			case CPU_SUBFIND_SAVEFINAL: return "CPU_SUBFIND_SAVEFINAL";
			case CPU_SUBFIND_TOT: return "CPU_SUBFIND_TOT";
			case CPU_FOF_TOT: return "CPU_FOF_TOT";
#endif
#if defined(fSIDM) || defined(rSIDM)
			case CPU_FSIDMSCATTER: return "CPU_FSIDMSCATTER";
			case CPU_FSIDMSCATTER_DF: return "CPU_FSIDMSCATTER_DF";
			case CPU_FSIDMSCATTER_RC: return "CPU_FSIDMSCATTER_RC";
			case CPU_RSIDMSCATTER: return "CPU_RSIDMSCATTER";
			case CPU_MSIDMPARALLEL: return "CPU_MSIDMPARALLEL";
			case CPU_MSIDMLOCAL: return "CPU_MSIDMLOCAL";
#endif		
			}

			return "UNKNOWN";
		}

		void print_CPU_steps()
		{
#if 0	
			for (int s = 0; s < CPU_PARTS; s++) {
				printf("%d: %s: %f\n", s, get_cpufields_name(s), CPU_Step[s]);
			}
#endif
			// only CPU_SNAPSHOT
			printf("%d: %s: %f\n", CPU_SNAPSHOT, get_cpufields_name(CPU_SNAPSHOT), CPU_Step[CPU_SNAPSHOT]);
		}

		int get_particle_type(uint64_t id) {
			return P[id].Type;
		}

		void get_particle_position(uint64_t id, double* pos)
		{
			pos[0] = P[id].Pos[0];
			pos[1] = P[id].Pos[1];
			pos[2] = P[id].Pos[2];
		}

		float get_particle_norm_value(int blocknr, uint64_t id)
		{
			//if (blocknr == IO_LASTENTRY) //counter
			//    return 1.0f;

			int type = P[id].Type;
			//enum iofields blocknr

			switch (blocknr)
			{
			case IO_POS:		/* positions */
				RETURN_NORM_VECTOR3(P[id].Pos);

			case IO_VEL:		/* velocities */
				RETURN_NORM_VECTOR3(P[id].Vel);

			case IO_ID:		/* particle ID */
				RETURN_NORM_VALUE(P[id].ID);

			case IO_MASS:		/* particle mass */
				//RETURN_NORM_VALUE(((double)P[id].Mass * (double)1e10)); // fix value 1^10
				RETURN_NORM_VALUE(P[id].Mass);

			case IO_SHEET_ORIENTATION:	/* initial particle sheet orientation */
				RETURN_NORM_EMPTY;

			case IO_INIT_DENSITY:	/* initial stream density */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
				RETURN_NORM_VALUE(GDE_INITDENSITY(id));
#endif

			case IO_CAUSTIC_COUNTER:	/* initial caustic counter */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
				RETURN_NORM_VALUE(P[id].caustic_counter);
#endif
				RETURN_NORM_EMPTY;

			case IO_SECONDORDERMASS:
				//RETURN_NORM_VALUE(P[id].Mass);
				RETURN_NORM_VALUE(P[id].OldAcc);

			case IO_U:			/* temperature */
#ifdef JD_RELAX_CLUSTERS
				RETURN_NORM_VALUE(SphP[id].U);
#else
				RETURN_NORM_VALUE(SphP[id].Entropy);
#endif

			case IO_RHO:		/* density */
				RETURN_NORM_VALUE(SphP[id].Density);

			case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(CHEMISTRY) || defined(UM_CHEMISTRY)
				RETURN_NORM_VALUE(SphP[id].elec);
#endif
				RETURN_NORM_EMPTY;

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)	//now that's weird. But I am afraid of removing this #if, I am not sure where it ends. AR.
			case IO_NH:		/* neutral hydrogen abundance */
				RETURN_NORM_VALUE(SphP[id].HI);

			case IO_HII:		/* ionized hydrogen abundance */
				RETURN_NORM_VALUE(SphP[id].HII);

			case IO_HeI:		/* neutral Helium */
				RETURN_NORM_VALUE(SphP[id].HeI);

			case IO_HeII:		/* ionized Heluum */
				RETURN_NORM_VALUE(SphP[id].HeII);

			case IO_HeIII:		/* double ionised Helium */
				RETURN_NORM_VALUE(SphP[id].HeIII);
#endif

			case IO_H2I:		/* H2 molecule */
				RETURN_NORM_VALUE(SphP[id].H2I);

			case IO_H2II:		/* ionised H2 molecule */
				RETURN_NORM_VALUE(SphP[id].H2II);

			case IO_HM:		/* H minus */
				RETURN_NORM_VALUE(SphP[id].HM);

			case IO_HeHII:		/* HeH+ */
#if defined (UM_CHEMISTRY)
				RETURN_NORM_VALUE(SphP[id].HeHII);
#endif
				RETURN_NORM_EMPTY;

			case IO_HD:		/* HD */
#if defined (UM_CHEMISTRY)
				RETURN_NORM_VALUE(SphP[id].HD);
#endif
				RETURN_NORM_EMPTY;

			case IO_DI:		/* D */
#if defined (UM_CHEMISTRY)
				RETURN_NORM_VALUE(SphP[id].DI);
#endif
				RETURN_NORM_EMPTY;

			case IO_DII:		/* D plus */
#if defined (UM_CHEMISTRY)
				RETURN_NORM_VALUE(SphP[id].DII);
#endif
				RETURN_NORM_EMPTY;

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
				RETURN_NORM_EMPTY;
#endif

			case IO_HSML:		/* SPH smoothing length */
				RETURN_NORM_VALUE(P[id].Hsml);

			case IO_DELAYTIME:
#ifdef WINDS
				RETURN_NORM_VALUE(SphP[id].DelayTime);
#endif
				RETURN_NORM_EMPTY;

			case IO_AGE:		/* Age of stars */
#ifdef STELLARAGE
				RETURN_NORM_VALUE(MPP(id).StellarAge);
#if defined(BLACK_HOLES)
				//RETURN_NORM_VALUE(BPP(id).StellarAge); //TODO
#endif

#endif
				RETURN_NORM_EMPTY;

			case IO_Z:			/* Gas and star metallicity */
#ifdef METALS
				RETURN_NORM_VALUE(P[id].Metallicity);
#endif
				RETURN_NORM_EMPTY;

			case IO_EGYPROM:		/* SN Energy Reservoir */
				RETURN_NORM_EMPTY;

			case IO_EGYCOLD:		/* Cold  SN Energy Reservoir */
				RETURN_NORM_EMPTY;

			case IO_VRMS:		/* Turbulence on kernel scale */
#ifdef JD_VTURB
				RETURN_NORM_VALUE(SphP[id].Vrms);
#endif
				RETURN_NORM_EMPTY;
			case IO_VBULK:
#ifdef JD_VTURB
				RETURN_NORM_VECTOR3(SphP[id].Vbulk);
#endif
				RETURN_NORM_EMPTY;
			case IO_VTAN:
#ifdef JD_DECOMPOSE_VTURB
				RETURN_NORM_VALUE(SphP[id].Vtan);
#endif
				RETURN_NORM_EMPTY;
			case IO_VRAD:
#ifdef JD_DECOMPOSE_VTURB
				RETURN_NORM_VALUE(SphP[id].Vrad);
#endif
				RETURN_NORM_EMPTY;
			case IO_VDIV:
#ifdef JD_VTURB
				RETURN_NORM_VALUE(SphP[id].DivVel);
#endif
				RETURN_NORM_EMPTY;
			case IO_VROT:
#ifdef JD_VTURB
				RETURN_NORM_VALUE(SphP[id].r.CurlVel);
#endif
				RETURN_NORM_EMPTY;
			case IO_TRUENGB:
#ifdef JD_VTURB
				RETURN_NORM_VALUE(P[id].TrueNGB);
#endif
				RETURN_NORM_EMPTY;
			case IO_DPP:
#ifdef JD_DPP
				RETURN_NORM_VALUE(SphP[id].Dpp);
#endif
				RETURN_NORM_EMPTY;

			case IO_BFLD:		/* Magnetic field */
#ifdef MAGNETIC
				RETURN_NORM_VECTOR3(SphP[id].BPred);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_pNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_NORM_DVECTORN(SphP[id].CRpNorm, LMB_CR_PROTONS);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_eNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_NORM_DVECTORN(SphP[id].CReNorm, LMB_CR_PROTONS);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_pSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_NORM_DVECTORN(SphP[id].CRpSlope, LMB_CR_PROTONS);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_eSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_NORM_DVECTORN(SphP[id].CReSlope, LMB_CR_PROTONS);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_pCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_NORM_VALUE(SphP[id].CRpCut);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_eCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_NORM_VALUE(SphP[id].CReCut);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_pPRESSURE:
#if defined LMB_SPECTRAL_CRs_PRESSURE_FROM_ICS
				RETURN_NORM_VALUE(SphP[id].CRpPressure);
#endif
				RETURN_NORM_EMPTY;

			case IO_LMBCR_ePRESSURE:
#if defined LMB_SPECTRALCRs_PRESSURE_FROM_ICS
				RETURN_NORM_VALUE(SphP[id].CRePressure);
#endif
				RETURN_NORM_EMPTY;

			case IO_RHO_OLD:
#if defined READ_RHO_OLD
				RETURN_NORM_VALUE(SphP[id].DensityOld);
#endif
				RETURN_NORM_EMPTY;

			case IO_MACH:
#if defined READ_MACH
				RETURN_NORM_VALUE(SphP[id].Shock_Mach);
#endif
				RETURN_NORM_EMPTY;

			case IO_BHMASS:
#ifdef BLACK_HOLES
				RETURN_NORM_VALUE(BPP(id).BH_Mass);
#endif
				RETURN_NORM_EMPTY;

			case IO_BHMDOT:
#ifdef BLACK_HOLES
				RETURN_NORM_VALUE(BPP(id).BH_Mdot);
#endif
				RETURN_NORM_EMPTY;

			case IO_BHPROGS:
#ifdef BLACK_HOLES
				RETURN_NORM_VALUE(BPP(id).BH_CountProgs);
#endif
				RETURN_NORM_EMPTY;

			case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
				RETURN_NORM_VALUE(BPP(id).BH_Mass_radio);
#endif
				RETURN_NORM_EMPTY;

			case IO_EOSXNUC:
				RETURN_NORM_EMPTY;

			case IO_DUSTL:
#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
				RETURN_NORM_FVECTORN(SphP[id].DustL, LT_NMetP);
#endif
				RETURN_NORM_EMPTY;

			case IO_DUSTS:
#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
				if (type == 0)
					RETURN_NORM_FVECTORN(SphP[id].DustS, LT_NMetP);
#endif
				RETURN_NORM_EMPTY;

			case IO_Zs:
#ifdef LT_STELLAREVOLUTION
				if (type == 4) {
					RETURN_NORM_FVECTORN(MPP(id).Metals, LT_NMetP);
				}
				else if (type == 0) {
					RETURN_NORM_FVECTORN(SphP[id].Metals, LT_NMetP);
				}
#endif
				RETURN_NORM_EMPTY;

			case IO_ZAGE:
#ifdef LT_ZAGE
				if (type == 4) {
					RETURN_NORM_VALUE(MPP(id).ZAge);
				}
				else if (type == 0) {
					RETURN_NORM_VALUE(SphP[id].ZAge);
				}
#endif
				RETURN_NORM_EMPTY;

			case IO_ZAGE_LLV:
#ifdef LT_ZAGE_LLV
				if (type == 4) {
					RETURN_NORM_VALUE(MPP(id).ZAge_llv);
				}
				else if (type == 0) {
					RETURN_NORM_VALUE(SphP[id].ZAge_llv);
				}
#endif
				RETURN_NORM_EMPTY;

			case IO_iMass:
#ifdef LT_STELLAREVOLUTION
				if (type == 4)
					RETURN_NORM_VALUE(MPP(id).iMass);
#endif
				RETURN_NORM_EMPTY;

			case IO_CONTRIB:
#if defined(LT_STELLAREVOLUTION) && defined(LT_TRACK_CONTRIBUTES)
				if (type == 4) {
					RETURN_NORM_VALUE(MPP(id).contrib);
				}
				else if (type == 0) {
					RETURN_NORM_VALUE(SphP[id].contrib);
				}
#endif
				RETURN_NORM_EMPTY;

			case IO_RADGAMMA:
				RETURN_NORM_EMPTY;

			case IO_DMHSML:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_NORM_VALUE(P[id].DM_Hsml);
#endif
				RETURN_NORM_EMPTY;

			case IO_DMDENSITY:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_NORM_VALUE(P[id].u.DM_Density);
#endif
				RETURN_NORM_EMPTY;

			case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_NORM_VALUE(P[id].v.DM_VelDisp);
#endif
				RETURN_NORM_EMPTY;

			case IO_EULERA:
#ifdef READ_EULER
				RETURN_NORM_VALUE(SphP[id].EulerA);
#endif
				RETURN_NORM_EMPTY;

			case IO_EULERB:
#ifdef READ_EULER
				RETURN_NORM_VALUE(SphP[id].EulerB);
#endif
				RETURN_NORM_EMPTY;

			case IO_ALFA2_DYN:
				RETURN_NORM_EMPTY;

			case IO_ETA2_DYN:
				RETURN_NORM_EMPTY;


			case IO_FB_M_H:		/* particle hot mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].M_h);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_M_C:		/* particle cold mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].M_c);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_M_MO:		/* particle molecular mass */
#if defined (GM_MUPPI) && defined(MV_OUTPUT_MMOL)
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].Fcoll);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_E_H:		/* particle hot energy */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].E_h);
#endif
				RETURN_NORM_EMPTY;


			case IO_FB_M_SF:		/* particle SF mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].M_sf);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_MF:		/* number of time steps in multiphase */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].MultiPhase);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_NMF:		/* number of multiphase stages */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].NMF);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_EOUT:		/* Energy output from MP particles */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].E_out);
#endif
				RETURN_NORM_EMPTY;
			case IO_FB_EREC:		/* Energy received from  MP particles */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].E_rec);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_CLOCK:		/* time to still pass in MP */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].clock);
#endif
				RETURN_NORM_EMPTY;

			case IO_FB_E_TOT_0:	/* multi phase hot energy time i-1 */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].Egy_tot_0);
#endif
			case IO_FB_TDYN:		/* cold phase dynamical time */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].tdyn);

				RETURN_NORM_EMPTY;
#endif

			case IO_FB_TSTARTMP:
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_NORM_VALUE(SphP[id].t_startMP);
#endif
				RETURN_NORM_EMPTY;

			case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
				RETURN_NORM_VALUE(SphP[id].XColdCloud);
#endif
				RETURN_NORM_EMPTY;

			case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
				RETURN_NORM_VALUE(SphP[id].Temperature);
#endif
				RETURN_NORM_EMPTY;

			case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
				RETURN_NORM_VALUE(SphP[id].Temperature);
#endif
				RETURN_NORM_EMPTY;

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
				RETURN_NORM_EMPTY;

			case IO_LASTENTRY:
				//endrun(220);
				RETURN_NORM_EMPTY;
			}

			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value)
		{
			//if (blocknr == IO_LASTENTRY) //counter
			//    return 1.0f;

			int type = P[id].Type;
			//enum iofields blocknr

			switch (blocknr)
			{
			case IO_POS:		/* positions */
				RETURN_ORIG_VECTOR3(P[id].Pos);

			case IO_VEL:		/* velocities */
				RETURN_ORIG_VECTOR3(P[id].Vel);

			case IO_ID:		/* particle ID */
				RETURN_ORIG_VALUE(P[id].ID);

			case IO_MASS:		/* particle mass */
				//RETURN_ORIG_VALUE(((double)P[id].Mass * (double)1e10)); // fix value 1^10
				RETURN_ORIG_VALUE(P[id].Mass);

			case IO_SHEET_ORIENTATION:	/* initial particle sheet orientation */
				RETURN_ORIG_EMPTY;

			case IO_INIT_DENSITY:	/* initial stream density */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
				RETURN_ORIG_VALUE(GDE_INITDENSITY(id));
#endif

			case IO_CAUSTIC_COUNTER:	/* initial caustic counter */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
				RETURN_ORIG_VALUE(P[id].caustic_counter);
#endif
				RETURN_ORIG_EMPTY;

			case IO_SECONDORDERMASS:
				//RETURN_ORIG_VALUE(P[id].Mass);
				RETURN_ORIG_VALUE(P[id].OldAcc);

			case IO_U:			/* temperature */
#ifdef JD_RELAX_CLUSTERS
				RETURN_ORIG_VALUE(SphP[id].U);
#else
				RETURN_ORIG_VALUE(SphP[id].Entropy);
#endif

			case IO_RHO:		/* density */
				RETURN_ORIG_VALUE(SphP[id].Density);

			case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(CHEMISTRY) || defined(UM_CHEMISTRY)
				RETURN_ORIG_VALUE(SphP[id].elec);
#endif
				RETURN_ORIG_EMPTY;

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)	//now that's weird. But I am afraid of removing this #if, I am not sure where it ends. AR.
			case IO_NH:		/* neutral hydrogen abundance */
				RETURN_ORIG_VALUE(SphP[id].HI);

			case IO_HII:		/* ionized hydrogen abundance */
				RETURN_ORIG_VALUE(SphP[id].HII);

			case IO_HeI:		/* neutral Helium */
				RETURN_ORIG_VALUE(SphP[id].HeI);

			case IO_HeII:		/* ionized Heluum */
				RETURN_ORIG_VALUE(SphP[id].HeII);

			case IO_HeIII:		/* double ionised Helium */
				RETURN_ORIG_VALUE(SphP[id].HeIII);
#endif

			case IO_H2I:		/* H2 molecule */
				RETURN_ORIG_VALUE(SphP[id].H2I);

			case IO_H2II:		/* ionised H2 molecule */
				RETURN_ORIG_VALUE(SphP[id].H2II);

			case IO_HM:		/* H minus */
				RETURN_ORIG_VALUE(SphP[id].HM);

			case IO_HeHII:		/* HeH+ */
#if defined (UM_CHEMISTRY)
				RETURN_ORIG_VALUE(SphP[id].HeHII);
#endif
				RETURN_ORIG_EMPTY;

			case IO_HD:		/* HD */
#if defined (UM_CHEMISTRY)
				RETURN_ORIG_VALUE(SphP[id].HD);
#endif
				RETURN_ORIG_EMPTY;

			case IO_DI:		/* D */
#if defined (UM_CHEMISTRY)
				RETURN_ORIG_VALUE(SphP[id].DI);
#endif
				RETURN_ORIG_EMPTY;

			case IO_DII:		/* D plus */
#if defined (UM_CHEMISTRY)
				RETURN_ORIG_VALUE(SphP[id].DII);
#endif
				RETURN_ORIG_EMPTY;

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
				RETURN_ORIG_EMPTY;
#endif

			case IO_HSML:		/* SPH smoothing length */
				RETURN_ORIG_VALUE(P[id].Hsml);

			case IO_DELAYTIME:
#ifdef WINDS
				RETURN_ORIG_VALUE(SphP[id].DelayTime);
#endif
				RETURN_ORIG_EMPTY;

			case IO_AGE:		/* Age of stars */
#ifdef STELLARAGE
				RETURN_ORIG_VALUE(MPP(id).StellarAge);
#if defined(BLACK_HOLES)
				//RETURN_ORIG_VALUE(BPP(id).StellarAge); //TODO
#endif

#endif
				RETURN_ORIG_EMPTY;

			case IO_Z:			/* Gas and star metallicity */
#ifdef METALS
				RETURN_ORIG_VALUE(P[id].Metallicity);
#endif
				RETURN_ORIG_EMPTY;

			case IO_EGYPROM:		/* SN Energy Reservoir */
				RETURN_ORIG_EMPTY;

			case IO_EGYCOLD:		/* Cold  SN Energy Reservoir */
				RETURN_ORIG_EMPTY;

			case IO_VRMS:		/* Turbulence on kernel scale */
#ifdef JD_VTURB
				RETURN_ORIG_VALUE(SphP[id].Vrms);
#endif
				RETURN_ORIG_EMPTY;
			case IO_VBULK:
#ifdef JD_VTURB
				RETURN_ORIG_VECTOR3(SphP[id].Vbulk);
#endif
				RETURN_ORIG_EMPTY;
			case IO_VTAN:
#ifdef JD_DECOMPOSE_VTURB
				RETURN_ORIG_VALUE(SphP[id].Vtan);
#endif
				RETURN_ORIG_EMPTY;
			case IO_VRAD:
#ifdef JD_DECOMPOSE_VTURB
				RETURN_ORIG_VALUE(SphP[id].Vrad);
#endif
				RETURN_ORIG_EMPTY;
			case IO_VDIV:
#ifdef JD_VTURB
				RETURN_ORIG_VALUE(SphP[id].DivVel);
#endif
				RETURN_ORIG_EMPTY;
			case IO_VROT:
#ifdef JD_VTURB
				RETURN_ORIG_VALUE(SphP[id].r.CurlVel);
#endif
				RETURN_ORIG_EMPTY;
			case IO_TRUENGB:
#ifdef JD_VTURB
				RETURN_ORIG_VALUE(P[id].TrueNGB);
#endif
				RETURN_ORIG_EMPTY;
			case IO_DPP:
#ifdef JD_DPP
				RETURN_ORIG_VALUE(SphP[id].Dpp);
#endif
				RETURN_ORIG_EMPTY;

			case IO_BFLD:		/* Magnetic field */
#ifdef MAGNETIC
				RETURN_ORIG_VECTOR3(SphP[id].BPred);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_pNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_ORIG_DVECTORN(SphP[id].CRpNorm, LMB_CR_PROTONS);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_eNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_ORIG_DVECTORN(SphP[id].CReNorm, LMB_CR_PROTONS);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_pSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_ORIG_DVECTORN(SphP[id].CRpSlope, LMB_CR_PROTONS);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_eSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_ORIG_DVECTORN(SphP[id].CReSlope, LMB_CR_PROTONS);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_pCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_ORIG_VALUE(SphP[id].CRpCut);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_eCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_ORIG_VALUE(SphP[id].CReCut);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_pPRESSURE:
#if defined LMB_SPECTRAL_CRs_PRESSURE_FROM_ICS
				RETURN_ORIG_VALUE(SphP[id].CRpPressure);
#endif
				RETURN_ORIG_EMPTY;

			case IO_LMBCR_ePRESSURE:
#if defined LMB_SPECTRALCRs_PRESSURE_FROM_ICS
				RETURN_ORIG_VALUE(SphP[id].CRePressure);
#endif
				RETURN_ORIG_EMPTY;

			case IO_RHO_OLD:
#if defined READ_RHO_OLD
				RETURN_ORIG_VALUE(SphP[id].DensityOld);
#endif
				RETURN_ORIG_EMPTY;

			case IO_MACH:
#if defined READ_MACH
				RETURN_ORIG_VALUE(SphP[id].Shock_Mach);
#endif
				RETURN_ORIG_EMPTY;

			case IO_BHMASS:
#ifdef BLACK_HOLES
				RETURN_ORIG_VALUE(BPP(id).BH_Mass);
#endif
				RETURN_ORIG_EMPTY;

			case IO_BHMDOT:
#ifdef BLACK_HOLES
				RETURN_ORIG_VALUE(BPP(id).BH_Mdot);
#endif
				RETURN_ORIG_EMPTY;

			case IO_BHPROGS:
#ifdef BLACK_HOLES
				RETURN_ORIG_VALUE(BPP(id).BH_CountProgs);
#endif
				RETURN_ORIG_EMPTY;

			case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
				RETURN_ORIG_VALUE(BPP(id).BH_Mass_radio);
#endif
				RETURN_ORIG_EMPTY;

			case IO_EOSXNUC:
				RETURN_ORIG_EMPTY;

			case IO_DUSTL:
#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
				RETURN_ORIG_FVECTORN(SphP[id].DustL, LT_NMetP);
#endif
				RETURN_ORIG_EMPTY;

			case IO_DUSTS:
#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
				if (type == 0)
					RETURN_ORIG_FVECTORN(SphP[id].DustS, LT_NMetP);
#endif
				RETURN_ORIG_EMPTY;

			case IO_Zs:
#ifdef LT_STELLAREVOLUTION
				if (type == 4) {
					RETURN_ORIG_FVECTORN(MPP(id).Metals, LT_NMetP);
				}
				else if (type == 0) {
					RETURN_ORIG_FVECTORN(SphP[id].Metals, LT_NMetP);
				}
#endif
				RETURN_ORIG_EMPTY;

			case IO_ZAGE:
#ifdef LT_ZAGE
				if (type == 4) {
					RETURN_ORIG_VALUE(MPP(id).ZAge);
				}
				else if (type == 0) {
					RETURN_ORIG_VALUE(SphP[id].ZAge);
				}
#endif
				RETURN_ORIG_EMPTY;

			case IO_ZAGE_LLV:
#ifdef LT_ZAGE_LLV
				if (type == 4) {
					RETURN_ORIG_VALUE(MPP(id).ZAge_llv);
				}
				else if (type == 0) {
					RETURN_ORIG_VALUE(SphP[id].ZAge_llv);
				}
#endif
				RETURN_ORIG_EMPTY;

			case IO_iMass:
#ifdef LT_STELLAREVOLUTION
				if (type == 4)
					RETURN_ORIG_VALUE(MPP(id).iMass);
#endif
				RETURN_ORIG_EMPTY;

			case IO_CONTRIB:
#if defined(LT_STELLAREVOLUTION) && defined(LT_TRACK_CONTRIBUTES)
				if (type == 4) {
					RETURN_ORIG_VALUE(MPP(id).contrib);
				}
				else if (type == 0) {
					RETURN_ORIG_VALUE(SphP[id].contrib);
				}
#endif
				RETURN_ORIG_EMPTY;

			case IO_RADGAMMA:
				RETURN_ORIG_EMPTY;

			case IO_DMHSML:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_ORIG_VALUE(P[id].DM_Hsml);
#endif
				RETURN_ORIG_EMPTY;

			case IO_DMDENSITY:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_ORIG_VALUE(P[id].u.DM_Density);
#endif
				RETURN_ORIG_EMPTY;

			case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_ORIG_VALUE(P[id].v.DM_VelDisp);
#endif
				RETURN_ORIG_EMPTY;

			case IO_EULERA:
#ifdef READ_EULER
				RETURN_ORIG_VALUE(SphP[id].EulerA);
#endif
				RETURN_ORIG_EMPTY;

			case IO_EULERB:
#ifdef READ_EULER
				RETURN_ORIG_VALUE(SphP[id].EulerB);
#endif
				RETURN_ORIG_EMPTY;

			case IO_ALFA2_DYN:
				RETURN_ORIG_EMPTY;

			case IO_ETA2_DYN:
				RETURN_ORIG_EMPTY;


			case IO_FB_M_H:		/* particle hot mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].M_h);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_M_C:		/* particle cold mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].M_c);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_M_MO:		/* particle molecular mass */
#if defined (GM_MUPPI) && defined(MV_OUTPUT_MMOL)
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].Fcoll);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_E_H:		/* particle hot energy */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].E_h);
#endif
				RETURN_ORIG_EMPTY;


			case IO_FB_M_SF:		/* particle SF mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].M_sf);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_MF:		/* number of time steps in multiphase */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].MultiPhase);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_NMF:		/* number of multiphase stages */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].NMF);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_EOUT:		/* Energy output from MP particles */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].E_out);
#endif
				RETURN_ORIG_EMPTY;
			case IO_FB_EREC:		/* Energy received from  MP particles */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].E_rec);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_CLOCK:		/* time to still pass in MP */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].clock);
#endif
				RETURN_ORIG_EMPTY;

			case IO_FB_E_TOT_0:	/* multi phase hot energy time i-1 */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].Egy_tot_0);
#endif
			case IO_FB_TDYN:		/* cold phase dynamical time */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].tdyn);

				RETURN_ORIG_EMPTY;
#endif

			case IO_FB_TSTARTMP:
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_ORIG_VALUE(SphP[id].t_startMP);
#endif
				RETURN_ORIG_EMPTY;

			case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
				RETURN_ORIG_VALUE(SphP[id].XColdCloud);
#endif
				RETURN_ORIG_EMPTY;

			case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
				RETURN_ORIG_VALUE(SphP[id].Temperature);
#endif
				RETURN_ORIG_EMPTY;

			case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
				RETURN_ORIG_VALUE(SphP[id].Temperature);
#endif
				RETURN_ORIG_EMPTY;

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
				RETURN_ORIG_EMPTY;

			case IO_LASTENTRY:
				//endrun(220);
				RETURN_ORIG_EMPTY;
			}

			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id)
		{
			//if (blocknr == IO_LASTENTRY) //counter
			//    return 1.0f;

			int type = P[id].Type;
			//enum iofields blocknr

			switch (blocknr)
			{
			case IO_POS:		/* positions */
				RETURN_COMP_VECTOR3(P[id].Pos);

			case IO_VEL:		/* velocities */
				RETURN_COMP_VECTOR3(P[id].Vel);

			case IO_ID:		/* particle ID */
				RETURN_COMP_VALUE(P[id].ID);

			case IO_MASS:		/* particle mass */
				//RETURN_COMP_VALUE(((double)P[id].Mass * (double)1e10)); // fix value 1^10
				RETURN_COMP_VALUE(P[id].Mass);

			case IO_SHEET_ORIENTATION:	/* initial particle sheet orientation */
				RETURN_COMP_EMPTY;

			case IO_INIT_DENSITY:	/* initial stream density */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
				RETURN_COMP_VALUE(GDE_INITDENSITY(id));
#endif

			case IO_CAUSTIC_COUNTER:	/* initial caustic counter */
#if defined(DISTORTIONTENSORPS) && defined(GDE_READIC)
				RETURN_COMP_VALUE(P[id].caustic_counter);
#endif
				RETURN_COMP_EMPTY;

			case IO_SECONDORDERMASS:
				//RETURN_COMP_VALUE(P[id].Mass);
				RETURN_COMP_VALUE(P[id].OldAcc);

			case IO_U:			/* temperature */
#ifdef JD_RELAX_CLUSTERS
				RETURN_COMP_VALUE(SphP[id].U);
#else
				RETURN_COMP_VALUE(SphP[id].Entropy);
#endif

			case IO_RHO:		/* density */
				RETURN_COMP_VALUE(SphP[id].Density);

			case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(CHEMISTRY) || defined(UM_CHEMISTRY)
				RETURN_COMP_VALUE(SphP[id].elec);
#endif
				RETURN_COMP_EMPTY;

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)	//now that's weird. But I am afraid of removing this #if, I am not sure where it ends. AR.
			case IO_NH:		/* neutral hydrogen abundance */
				RETURN_COMP_VALUE(SphP[id].HI);

			case IO_HII:		/* ionized hydrogen abundance */
				RETURN_COMP_VALUE(SphP[id].HII);

			case IO_HeI:		/* neutral Helium */
				RETURN_COMP_VALUE(SphP[id].HeI);

			case IO_HeII:		/* ionized Heluum */
				RETURN_COMP_VALUE(SphP[id].HeII);

			case IO_HeIII:		/* double ionised Helium */
				RETURN_COMP_VALUE(SphP[id].HeIII);
#endif

			case IO_H2I:		/* H2 molecule */
				RETURN_COMP_VALUE(SphP[id].H2I);

			case IO_H2II:		/* ionised H2 molecule */
				RETURN_COMP_VALUE(SphP[id].H2II);

			case IO_HM:		/* H minus */
				RETURN_COMP_VALUE(SphP[id].HM);

			case IO_HeHII:		/* HeH+ */
#if defined (UM_CHEMISTRY)
				RETURN_COMP_VALUE(SphP[id].HeHII);
#endif
				RETURN_COMP_EMPTY;

			case IO_HD:		/* HD */
#if defined (UM_CHEMISTRY)
				RETURN_COMP_VALUE(SphP[id].HD);
#endif
				RETURN_COMP_EMPTY;

			case IO_DI:		/* D */
#if defined (UM_CHEMISTRY)
				RETURN_COMP_VALUE(SphP[id].DI);
#endif
				RETURN_COMP_EMPTY;

			case IO_DII:		/* D plus */
#if defined (UM_CHEMISTRY)
				RETURN_COMP_VALUE(SphP[id].DII);
#endif
				RETURN_COMP_EMPTY;

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
				RETURN_COMP_EMPTY;
#endif

			case IO_HSML:		/* SPH smoothing length */
				RETURN_COMP_VALUE(P[id].Hsml);

			case IO_DELAYTIME:
#ifdef WINDS
				RETURN_COMP_VALUE(SphP[id].DelayTime);
#endif
				RETURN_COMP_EMPTY;

			case IO_AGE:		/* Age of stars */
#ifdef STELLARAGE
				RETURN_COMP_VALUE(MPP(id).StellarAge);
#if defined(BLACK_HOLES)
				//RETURN_COMP_VALUE(BPP(id).StellarAge); //TODO
#endif

#endif
				RETURN_COMP_EMPTY;

			case IO_Z:			/* Gas and star metallicity */
#ifdef METALS
				RETURN_COMP_VALUE(P[id].Metallicity);
#endif
				RETURN_COMP_EMPTY;

			case IO_EGYPROM:		/* SN Energy Reservoir */
				RETURN_COMP_EMPTY;

			case IO_EGYCOLD:		/* Cold  SN Energy Reservoir */
				RETURN_COMP_EMPTY;

			case IO_VRMS:		/* Turbulence on kernel scale */
#ifdef JD_VTURB
				RETURN_COMP_VALUE(SphP[id].Vrms);
#endif
				RETURN_COMP_EMPTY;
			case IO_VBULK:
#ifdef JD_VTURB
				RETURN_COMP_VECTOR3(SphP[id].Vbulk);
#endif
				RETURN_COMP_EMPTY;
			case IO_VTAN:
#ifdef JD_DECOMPOSE_VTURB
				RETURN_COMP_VALUE(SphP[id].Vtan);
#endif
				RETURN_COMP_EMPTY;
			case IO_VRAD:
#ifdef JD_DECOMPOSE_VTURB
				RETURN_COMP_VALUE(SphP[id].Vrad);
#endif
				RETURN_COMP_EMPTY;
			case IO_VDIV:
#ifdef JD_VTURB
				RETURN_COMP_VALUE(SphP[id].DivVel);
#endif
				RETURN_COMP_EMPTY;
			case IO_VROT:
#ifdef JD_VTURB
				RETURN_COMP_VALUE(SphP[id].r.CurlVel);
#endif
				RETURN_COMP_EMPTY;
			case IO_TRUENGB:
#ifdef JD_VTURB
				RETURN_COMP_VALUE(P[id].TrueNGB);
#endif
				RETURN_COMP_EMPTY;
			case IO_DPP:
#ifdef JD_DPP
				RETURN_COMP_VALUE(SphP[id].Dpp);
#endif
				RETURN_COMP_EMPTY;

			case IO_BFLD:		/* Magnetic field */
#ifdef MAGNETIC
				RETURN_COMP_VECTOR3(SphP[id].BPred);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_pNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_COMP_DVECTORN(SphP[id].CRpNorm, LMB_CR_PROTONS);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_eNORM:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_COMP_DVECTORN(SphP[id].CReNorm, LMB_CR_PROTONS);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_pSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_COMP_DVECTORN(SphP[id].CRpSlope, LMB_CR_PROTONS);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_eSLOPE:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_COMP_DVECTORN(SphP[id].CReSlope, LMB_CR_PROTONS);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_pCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_PROTONS))
				RETURN_COMP_VALUE(SphP[id].CRpCut);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_eCUT:
#if (defined(LMB_SPECTRAL_CRs_FROM_ICS) && defined(LMB_CR_ELECTRONS))
				RETURN_COMP_VALUE(SphP[id].CReCut);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_pPRESSURE:
#if defined LMB_SPECTRAL_CRs_PRESSURE_FROM_ICS
				RETURN_COMP_VALUE(SphP[id].CRpPressure);
#endif
				RETURN_COMP_EMPTY;

			case IO_LMBCR_ePRESSURE:
#if defined LMB_SPECTRALCRs_PRESSURE_FROM_ICS
				RETURN_COMP_VALUE(SphP[id].CRePressure);
#endif
				RETURN_COMP_EMPTY;

			case IO_RHO_OLD:
#if defined READ_RHO_OLD
				RETURN_COMP_VALUE(SphP[id].DensityOld);
#endif
				RETURN_COMP_EMPTY;

			case IO_MACH:
#if defined READ_MACH
				RETURN_COMP_VALUE(SphP[id].Shock_Mach);
#endif
				RETURN_COMP_EMPTY;

			case IO_BHMASS:
#ifdef BLACK_HOLES
				RETURN_COMP_VALUE(BPP(id).BH_Mass);
#endif
				RETURN_COMP_EMPTY;

			case IO_BHMDOT:
#ifdef BLACK_HOLES
				RETURN_COMP_VALUE(BPP(id).BH_Mdot);
#endif
				RETURN_COMP_EMPTY;

			case IO_BHPROGS:
#ifdef BLACK_HOLES
				RETURN_COMP_VALUE(BPP(id).BH_CountProgs);
#endif
				RETURN_COMP_EMPTY;

			case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
				RETURN_COMP_VALUE(BPP(id).BH_Mass_radio);
#endif
				RETURN_COMP_EMPTY;

			case IO_EOSXNUC:
				RETURN_COMP_EMPTY;

			case IO_DUSTL:
#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
				RETURN_COMP_FVECTORN(SphP[id].DustL, LT_NMetP);
#endif
				RETURN_COMP_EMPTY;

			case IO_DUSTS:
#if defined(LT_STELLAREVOLUTION) && defined(GL_CR_DUST)
				if (type == 0)
					RETURN_COMP_FVECTORN(SphP[id].DustS, LT_NMetP);
#endif
				RETURN_COMP_EMPTY;

			case IO_Zs:
#ifdef LT_STELLAREVOLUTION
				if (type == 4) {
					RETURN_COMP_FVECTORN(MPP(id).Metals, LT_NMetP);
				}
				else if (type == 0) {
					RETURN_COMP_FVECTORN(SphP[id].Metals, LT_NMetP);
				}
#endif
				RETURN_COMP_EMPTY;

			case IO_ZAGE:
#ifdef LT_ZAGE
				if (type == 4) {
					RETURN_COMP_VALUE(MPP(id).ZAge);
				}
				else if (type == 0) {
					RETURN_COMP_VALUE(SphP[id].ZAge);
				}
#endif
				RETURN_COMP_EMPTY;

			case IO_ZAGE_LLV:
#ifdef LT_ZAGE_LLV
				if (type == 4) {
					RETURN_COMP_VALUE(MPP(id).ZAge_llv);
				}
				else if (type == 0) {
					RETURN_COMP_VALUE(SphP[id].ZAge_llv);
				}
#endif
				RETURN_COMP_EMPTY;

			case IO_iMass:
#ifdef LT_STELLAREVOLUTION
				if (type == 4)
					RETURN_COMP_VALUE(MPP(id).iMass);
#endif
				RETURN_COMP_EMPTY;

			case IO_CONTRIB:
#if defined(LT_STELLAREVOLUTION) && defined(LT_TRACK_CONTRIBUTES)
				if (type == 4) {
					RETURN_COMP_VALUE(MPP(id).contrib);
				}
				else if (type == 0) {
					RETURN_COMP_VALUE(SphP[id].contrib);
				}
#endif
				RETURN_COMP_EMPTY;

			case IO_RADGAMMA:
				RETURN_COMP_EMPTY;

			case IO_DMHSML:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_COMP_VALUE(P[id].DM_Hsml);
#endif
				RETURN_COMP_EMPTY;

			case IO_DMDENSITY:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_COMP_VALUE(P[id].u.DM_Density);
#endif
				RETURN_COMP_EMPTY;

			case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
				RETURN_COMP_VALUE(P[id].v.DM_VelDisp);
#endif
				RETURN_COMP_EMPTY;

			case IO_EULERA:
#ifdef READ_EULER
				RETURN_COMP_VALUE(SphP[id].EulerA);
#endif
				RETURN_COMP_EMPTY;

			case IO_EULERB:
#ifdef READ_EULER
				RETURN_COMP_VALUE(SphP[id].EulerB);
#endif
				RETURN_COMP_EMPTY;

			case IO_ALFA2_DYN:
				RETURN_COMP_EMPTY;

			case IO_ETA2_DYN:
				RETURN_COMP_EMPTY;


			case IO_FB_M_H:		/* particle hot mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].M_h);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_M_C:		/* particle cold mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].M_c);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_M_MO:		/* particle molecular mass */
#if defined (GM_MUPPI) && defined(MV_OUTPUT_MMOL)
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].Fcoll);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_E_H:		/* particle hot energy */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].E_h);
#endif
				RETURN_COMP_EMPTY;


			case IO_FB_M_SF:		/* particle SF mass */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].M_sf);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_MF:		/* number of time steps in multiphase */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].MultiPhase);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_NMF:		/* number of multiphase stages */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].NMF);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_EOUT:		/* Energy output from MP particles */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].E_out);
#endif
				RETURN_COMP_EMPTY;
			case IO_FB_EREC:		/* Energy received from  MP particles */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].E_rec);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_CLOCK:		/* time to still pass in MP */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].clock);
#endif
				RETURN_COMP_EMPTY;

			case IO_FB_E_TOT_0:	/* multi phase hot energy time i-1 */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].Egy_tot_0);
#endif
			case IO_FB_TDYN:		/* cold phase dynamical time */
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].tdyn);

				RETURN_COMP_EMPTY;
#endif

			case IO_FB_TSTARTMP:
#ifdef GM_MUPPI
				if (type == 0)
					RETURN_COMP_VALUE(SphP[id].t_startMP);
#endif
				RETURN_COMP_EMPTY;

			case IO_CLDX:
#ifdef LT_STELLAREVOLUTION
				RETURN_COMP_VALUE(SphP[id].XColdCloud);
#endif
				RETURN_COMP_EMPTY;

			case IO_HTEMP:
#ifdef LT_STELLAREVOLUTION
				RETURN_COMP_VALUE(SphP[id].Temperature);
#endif
				RETURN_COMP_EMPTY;

			case IO_TEMP:
#ifdef LT_METAL_COOLING_WAL
				RETURN_COMP_VALUE(SphP[id].Temperature);
#endif
				RETURN_COMP_EMPTY;

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
				RETURN_COMP_EMPTY;

			case IO_LASTENTRY:
				//endrun(220);
				RETURN_COMP_EMPTY;
			}

			RETURN_COMP_EMPTY;
		}

		//double get_particle_radius(uint64_t id) {
		//	double hsml = get_particle_hsml(id);
		//	double radius = 2.0 * hsml;

		//	if (std::isnan(radius))
		//		radius = 0;

		//	return radius;
		//}

		double get_particle_hsml(uint64_t id) {
			int type = P[id].Type;
			if (type != 0) {
				return 0;
			}

			double volume = 0;

			double hsml = get_particle_norm_value(IO_HSML, id);

#if 0
			if (false/*std::isnan(hsml)*/)
			{
				double mass = get_particle_mass(id);// *All.UnitMass_in_g; //SOLAR_MASS;// / (double)1e10;
				double density = get_particle_rho(id);// *All.UnitDensity_in_cgs;

				volume = (mass / density);
			}
			else {
				double hsml_in_cm = hsml;// *All.UnitLength_in_cm;
				volume = All.HubbleParam * hsml_in_cm * hsml_in_cm * hsml_in_cm;
			}

			double radius = cbrt(volume * 3.0 * M_PI / 4.0);// / All.UnitLength_in_cm;
			hsml = radius / 2.0;
#endif

#if 0
			if (true/*std::isnan(hsml)*/) {
				double mass = get_particle_mass(id);
				double density = get_particle_rho(id);
				double cn = 0.1f; // All.HubbleParam; //4.1f; //TODO

				hsml = cn * cbrt(mass / density);
			}
#endif

			return hsml;
		}

		double get_particle_mass(uint64_t id) {
			int type = P[id].Type;
			if (type != 0) {
				return 0;
			}			

			return get_particle_norm_value(IO_MASS, id);
		}

		int get_particle_rho_blocknr() {
			return IO_RHO;
		}

		double get_particle_rho(uint64_t id) {
			int type = P[id].Type;
			if (type != 0) {
				return 0;
			}			

			return get_particle_norm_value(IO_RHO, id);
		}

		std::string get_particle_unit(int blocknr)
		{
			switch (blocknr)
			{
			case IO_POS:		/* positions */
				return "Mpc/h";

			case IO_VEL:		/* velocities */
				return "km/s";

			case IO_MASS:		/* particle mass */
				return "1e10 Msun/h";
			}

			return "";
		}

		void gadget_init_lib(int world_rank, int world_size)
		{
			init_lib(world_rank, world_size);
		}
		void gadget_finish_lib()
		{
			finish_lib();
		}
		void gadget_set_parameter(int ICFormat, int SnapFormat, int NumFilesWrittenInParallel, int MaxMemSize, double BufferSize, double PartAllocFactor, int TotBHs)
		{
			set_parameter(ICFormat, SnapFormat, NumFilesWrittenInParallel, MaxMemSize, BufferSize, PartAllocFactor, TotBHs);
		}
		void gadget_read_parameter_file(char* fname, char* tag[], void** addr, int* id, int nt)
		{
			read_parameter_file(fname, tag, addr, id, nt);
		}
		void gadget_read_ic(char* fname)
		{
			read_ic(fname);
		}
		void gadget_mymalloc_init()
		{
			mymalloc_init();
		}

		void gadget_set_units()
		{
			set_units();
		}

		std::string gadget_get_type_name(int type)
		{
			switch (type) {
			case 0:
				return "Gas";
			case 1:
				return "DM Halos";
			case 2:
				return "Disk";
			case 3:
				return "Bulge";
			case 4:
				return "Stars";
#ifdef BLACK_HOLES
			case 5:
				return "Black Holes";
#else
			case 5:
				return "Bndry";
#endif
			}

			return "Unknown";
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks)
		{
			types_and_blocks.resize(sizeof(all_particles_blocks) / sizeof(int));
			memcpy(types_and_blocks.data(), all_particles_blocks, sizeof(all_particles_blocks));
		}

		std::string gadget_get_dataset_name(int blocknr)
		{
			char buf[500];
			get_dataset_name((iofields)blocknr, buf);

			return std::string(buf);
		}

		void print_types_and_blocks_local()
		{
			printf("\n");

			//char buf[500];
			for (int type = 0; type < 6; type++) {
				printf("Type: %s (%d)\n", gadget_get_type_name(type).c_str(), type);
				for (int blocknr = 0; blocknr < IO_LASTENTRY; blocknr++) {

					if (all_particles_blocks[6 * blocknr + type] > 0) {
						std::string buf = gadget_get_dataset_name((iofields)blocknr);
						printf("\t%s (%d)\n", buf.c_str(), blocknr);
					}
				}
			}
		}

		void print_types_and_blocks(std::vector<int>& types_and_blocks)
		{
			printf("\nAll snapshots contain:\n");

			//char buf[500];
			for (int type = 0; type < 6; type++) {
				printf("Type: %s (%d)\n", gadget_get_type_name(type).c_str(), type);
				for (int blocknr = 0; blocknr < IO_LASTENTRY; blocknr++) {

					if (types_and_blocks[6 * blocknr + type] > 0) {
						std::string buf = gadget_get_dataset_name((iofields)blocknr);
						printf("\t%s (%d)\n", buf.c_str(), blocknr);
					}
				}
			}
		}

		int get_num_types() {
			return 6;
		}
		int get_num_blocks() {
			return IO_LASTENTRY;
		}
	}//io
}//gadget