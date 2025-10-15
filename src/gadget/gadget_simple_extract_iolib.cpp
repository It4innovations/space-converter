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

#include "gadget_simple_extract_iolib.h"

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

 //#include "gadget_simple.h"
#include "convert_common.h"

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <string>
 //#include <format>

#ifdef _WIN32
#define FSEEKO(s,off,orig) _fseeki64(s,(__int64)(off),orig)
#define FTELLO(s) _ftelli64(s)
#else
#define FSEEKO fseeko
#define FTELLO ftello
#endif

#ifdef WITH_OPENMP
# include <omp.h>
#endif

//#include "gadget_simple.h"
#include "convert_vdb.h"

#include <mpi.h>

extern "C"
{

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"

}

//#define VALUE_CONVERT_DEFAULT	0
//#define VALUE_CONVERT_LOG		1

//#include <vector>
//#include <string.h>

#define BNAMESIZE 4                    // the size of a block name in char

#define NTYPES 6                       // how many types of particles there are
#define NALL   NTYPES                  //

#define ALL_PARTICLES  0               // used in get_particles_num
#define FILE_PARTICLES 1

#define GAS   0
#define DM1   1
#define HALO  1 
#define DM2   2 
#define DISK  2
#define DM3   3
#define BULGE 3
#define STARS 4
#define BH    5
#define BNDRY 5


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

#ifdef GADGET_MAX_HSML
#	define GADGET_MAX_BLOCKS IO_SFR //(IO_HSML+1)
#else
#	define GADGET_MAX_BLOCKS IO_LASTENTRY
#endif

namespace gadget_simple {


	typedef float float_t;                 // the floating type for reading the data

	// shortcuts for some integer types
	typedef unsigned int uint;
	typedef unsigned long long int ull;
	typedef ull nparts_t[NTYPES + 1];

	typedef char blockname_t[BNAMESIZE + 1];
	typedef blockname_t blockname_list[];

	//typedef struct {
	//	int     type;
	//	int     nQ;
	//	std::vector<std::string> Q;
	//} componentlist_t;

	//struct block_names {
	//	enum names {
	//		HEAD,
	//		POS,
	//		VEL,
	//		ID,
	//		MASS,
	//		U,
	//		TEMP,
	//		RHO,
	//		NE,
	//		NH,
	//		HSML,
	//		SFR,
	//		AGE,
	//		Z,
	//		Zs,
	//		iM,
	//		ZAGE,
	//		ZALV,
	//		CLDX,
	//		TSTP,
	//		POT,
	//		ACCE,
	//		ENDT,
	//		IDU,
	//		HOTT,
	//		MHOT,
	//		MCLD,
	//		EHOT,
	//		MSF,
	//		MFST,
	//		NMF,
	//		EOUT,
	//		EREC,
	//		EOLD,
	//		TDYN,
	//		SFRo,
	//		CLCK,
	//		Egy0,
	//		GRAD,
	//		BHMA,
	//		BHMD,
	//		BHPC,
	//		ACRB
	//	};
	//};

	//extern char** envp;

	//typedef struct { int tag1; char name[BNAMESIZE]; int len; int tag2; } head_t;
	//typedef struct {
	//	int    npart[NTYPES];			// number of particles of each type in this file
	//	double mass[NTYPES];		        // the mass table
	//	double time;                          // the expansion factor of the universe
	//	double redshift;                      // the redshift of the universe
	//	int    flag_sfr;                      //     ( signals whether the star formation was on )
	//	int    flag_feedback;                 //     ( signals whether the feedback was on )
	//	uint   npartTotal[NTYPES];	        // the total number of particles of each type in this snapshot.
	//	// this differs from npart if this is a  multi-file snapshot
	//	int    flag_cooling;                  //     ( signals whether the cooling was on )
	//	int    num_files;                     // how many files for this snapshot
	//	double BoxSize;                       // the size of the computational box
	//	double Omega0;                        //     ( the content of matter in the universe )
	//	double OmegaLambda;                   //     ( the content of dark energy in the universe )
	//	double HubbleParam;                   //     ( the hubble parameter )
	//	int    flag_stellarage;               //     ( whether the stars age is saved )
	//	int    flag_metals;                   //     ( signals whether the stellar evolution is on )
	//	uint   npartTotalHighWord[NTYPES];    // the highest 4 bytes of the total number of particles
	//	int    flag_entropy_instead_u;        //     ( signals that the entropy is integrated instead of the energy )
	//	int    flag_doubleprecision;          // signals that data are saved in double precision
	//	int    flag_ic_info;
	//	float  lpt_scalingfactor;

	//	char   fill[18];
	//	char   names[15][2];
	//} io_header;


	typedef struct
	{
		blockname_t name;
		int  active[NTYPES];
		size_t  fact;
		size_t  datasize;
	} blockdata_t;

	int check_number_of_files(std::string& basename);              // check the consistency of the files number
	int  seek_block(FILE*, const blockname_t);                         // seek a specific block in a specific file
	int  get_particles_num(io_header* snaphead, nparts_t* nparts_file, nparts_t* nparts_all);              // get how many particles are in a snapshot
	//int  get_known_block_list(const char*, blockdata_t**);
	int  get_knownblock_idx(blockdata_t*, int, const blockname_t);
	int position_and_skip_particles(blockdata_t*, int, const blockname_t, int, nparts_t*, FILE*);

	//int extract_data(std::string& basename);
	void gadget_simple_read_snaphead(FILE* fd, io_header* snaphead);

	class GadgetData {

	public:
		GadgetData() :data_name(""), num_elements(0), size_of_elem(0), data_size(0) {
			//data.resize(0);
		}

		float get_fvalue(size_t i) {
			if (num_elements != 1) {
				if (num_elements != 0)
				printf("Error: get_fvalue called with num_elements != 1 (num_elements: %lld)\n", num_elements);
				return 0;
			}
			if (size_of_elem == 8) {
				return ((double*)data.data())[i];
			}
			else if (size_of_elem == 4) {
				return ((float*)data.data())[i];
			}
			return 0;
		}

		int get_ivalue(size_t i) {
			if (num_elements != 1) {
				if (num_elements != 0)
				printf("Error: get_ivalue called with num_elements != 1 (num_elements: %lld)\n", num_elements);
				return 0;
			}
			if (size_of_elem == 8) {
				return ((long long int*)data.data())[i];
			}
			else if (size_of_elem == 4) {
				return ((int*)data.data())[i];
			}
			return 0;
		}

		double get_dvalue(size_t i) {
			if (num_elements != 1) {
				if (num_elements != 0)
				printf("Error: get_dvalue called with num_elements != 1 (num_elements: %lld)\n", num_elements);
				return 0;
			}
			if (size_of_elem == 8) {
				return ((double*)data.data())[i];
			}
			else if (size_of_elem == 4) {
				return ((float*)data.data())[i];
			}
			return 0;
		}

		void get_fvalue3(size_t i, float v[3]) {
			if (num_elements != 3) {
				if (i == 0)
					printf("Error: get_fvalue3 called with num_elements != 3 (%s: %d)\n", data_name.c_str(), num_elements);

				v[0] = v[1] = v[2] = 0;
				return;
			}
			if (size_of_elem == 8) {
				v[0] = ((double*)data.data())[i * num_elements + 0];
				v[1] = ((double*)data.data())[i * num_elements + 1];
				v[2] = ((double*)data.data())[i * num_elements + 2];
			}
			else if (size_of_elem == 4) {
				v[0] = ((float*)data.data())[i * num_elements + 0];
				v[1] = ((float*)data.data())[i * num_elements + 1];
				v[2] = ((float*)data.data())[i * num_elements + 2];
			}
			else {
				v[0] = v[1] = v[2] = 0; // Invalid size_of_elem	
			}
		}

		void get_dvalue3(size_t i, double v[3]) {
			if (num_elements != 3) {
				if (i == 0)
					printf("Error: get_dvalue3 called with num_elements != 3 (%s: %d)\n", data_name.c_str(), num_elements);

				v[0] = v[1] = v[2] = 0;
				return;
			}
			if (size_of_elem == 8) {
				v[0] = ((double*)data.data())[i * num_elements + 0];
				v[1] = ((double*)data.data())[i * num_elements + 1];
				v[2] = ((double*)data.data())[i * num_elements + 2];
			}
			else if (size_of_elem == 4) {
				v[0] = ((float*)data.data())[i * num_elements + 0];
				v[1] = ((float*)data.data())[i * num_elements + 1];
				v[2] = ((float*)data.data())[i * num_elements + 2];
			}
			else {
				v[0] = v[1] = v[2] = 0; // Invalid size_of_elem	
			}
		}

		std::string get_name() {
			return data_name;
		}

		size_t get_count() {
			if (num_elements == 0 || size_of_elem == 0) {
				//printf("Error: get_count called with num_elements == 0 or size_of_elem == 0\n");
				return 0;
			}
			return (data_size) / (num_elements * size_of_elem);
		}

		size_t get_nelem() {
			return num_elements;
		}

		// Function to merge data from another GadgetData object
		void merge_data(const GadgetData& other, size_t count) {
			if (data_size > 0) {
				if (this->data_name != other.data_name) {
					printf("Cannot merge: data names do not match\n");
					return;
				}

				if (this->num_elements != other.num_elements) {
					printf("Cannot merge: Number of elements per entry does not match\n");
					return;
				}

				if (this->size_of_elem != other.size_of_elem) {
					printf("Cannot merge: Type of elements per entry does not match\n");
					return;
				}
			}
			else {
				this->data_name = other.data_name;
				this->num_elements = other.num_elements;
				this->size_of_elem = other.size_of_elem;
				//this->count = other.count;
			}

			size_t alloc_size = count * other.num_elements * other.size_of_elem;

			if (data_size + alloc_size > data.size()) {
				printf("Cannot merge: data_size + alloc_size > data.size()\n");
				return;
			}

			memcpy(&data[data_size], other.data.data(), alloc_size);
			data_size += alloc_size;
		}

		void reserve_memory(size_t size) {
			data.resize(size);
			memset(data.data(), 0, size);
			data_size = 0;
		}

	public:
		//int ptype;
		std::string data_name;
		size_t num_elements;
		//size_t count;
		int size_of_elem;
		std::vector<uint8_t> data;
		size_t data_size;
		//size_t extra_data_space;
	};

	GadgetData gadget_datas[NTYPES][GADGET_MAX_BLOCKS];
	nparts_t npart_all;
	nparts_t npart_local_rank;

	int nfiles = 0;
	std::vector<blockdata_t> knownblocks;

	//std::vector<componentlist_t> components_list;
	//io_header  snapheader;
	double redshift = 0.0;		
	double hubble_param = 1.0;

	double steps_time[2];
	//
	//float g_bbox_pos[3] = { 0,0,0 };
	//float g_bbox_size[3] = { 0,0,0 };
	//int g_vdb_size[3] = { 0,0,0 };
	//float g_vdb_transform = 0;
	//std::string g_vdb_out_dir;
	//float g_min_value = 0;
	//float g_max_value = 0;
	//int g_value_convert = 0;
	//extract_data_t g_extract_data;

	//char** envp;

	/* Size of each input chunk to be
	   read and allocate for. */
#ifndef  READALL_CHUNK
#define  READALL_CHUNK  262144
#endif

#define  READALL_OK          0  /* Success */
#define  READALL_INVALID    -1  /* Invalid parameters */
#define  READALL_ERROR      -2  /* Stream error */
#define  READALL_TOOMUCH    -3  /* Too much input */
#define  READALL_NOMEM      -4  /* Out of memory */

	int readall(FILE* in, char** dataptr, size_t* sizeptr)
	{
		char* data = NULL, * temp;
		size_t size = 0;
		size_t used = 0;
		size_t n;

		/* None of the parameters can be NULL. */
		if (in == NULL || dataptr == NULL || sizeptr == NULL)
			return READALL_INVALID;

		/* A read error already occurred? */
		if (ferror(in))
			return READALL_ERROR;

		while (1) {

			if (used + READALL_CHUNK + 1 > size) {
				size = used + READALL_CHUNK + 1;

				/* Overflow check. Some ANSI C compilers
				   may optimize this away, though. */
				if (size <= used) {
					free(data);
					return READALL_TOOMUCH;
				}

				temp = (char*)realloc(data, size);
				if (temp == NULL) {
					free(data);
					return READALL_NOMEM;
				}
				data = temp;
			}

			n = fread(data + used, 1, READALL_CHUNK, in);
			if (n == 0)
				break;

			used += n;
		}

		if (ferror(in)) {
			free(data);
			return READALL_ERROR;
		}

		temp = (char*)realloc(data, used + 1);
		if (temp == NULL) {
			free(data);
			return READALL_NOMEM;
		}
		data = temp;
		data[used] = '\0';

		*dataptr = data;
		*sizeptr = used;

		return READALL_OK;
	}

	void fill_known_blocks() {
		knownblocks.clear();
		for (int bl = 0; bl < GADGET_MAX_BLOCKS; bl++) {
			blockdata_t block;
			char label[4];
			get_Tab_IO_Label((iofields)bl, label);
			block.name[0] = label[0];
			block.name[1] = label[1];
			block.name[2] = label[2];
			block.name[3] = label[3];
			block.name[4] = '\0';

			block.active[0] = 1;
			block.active[1] = 1;
			block.active[2] = 1;
			block.active[3] = 1;
			block.active[4] = 1;
			block.active[5] = 1;
			//block.active[6] = 1;

			block.fact = get_values_per_blockelement((iofields)bl);
			if (block.fact == 0) {
				//printf("WARNING: block %s has no values per block element, skipping it\n", block.name);
				continue;
			}
			block.datasize = (block.fact == 0) ? 0 : get_bytes_per_blockelement((iofields)bl, 1) / block.fact;

			knownblocks.push_back(block);
		}
	}

	int read_part_file_per_rank(std::string& basename, int world_rank, int world_size)
	{
		// More ranks than files: split each file among  world_size / nfiles ranks
		int ranks_per_file = world_size / nfiles;
		int file_rank_id = world_rank / ranks_per_file;
		int local_rank = world_rank % ranks_per_file;

		//size_t tot_count = 0;
		nparts_t orig_count;// gio_data.readNumElems(gio_data_rank_id);
		nparts_t offset;

		// find npart_local_rank
		int f = file_rank_id;
		//nparts_t npart_local;                     // will contains the nr. of particles for this file
		//std::vector<char>     name_vec(strlen(basename) + 10);
		std::string name;     // the name of this file

		// open the fth sub-file
		//
		if (nfiles > 1)
			//name = std::format("{}.{}", basename, f);
			name = std::string(basename) + "." + std::to_string(f);
		else
			//name = std::format("{}", basename);
			name = std::string(basename);

		FILE* file = fopen(name.c_str(), "rb");
		if (file == NULL)
		{
			// that should not happen, but deal with it
			//
			printf("I cannot open the file %s\n", name.c_str());
			exit(-1);
		}

		// get how many particles are in this file
		//
		io_header snaphead;
		gadget_simple_read_snaphead(file, &snaphead);
		get_particles_num(&snaphead, &orig_count, &npart_all);

		npart_local_rank[NALL] = 0;
		for (int i = 0; i < NTYPES; i++) {
			std::size_t chunk_size = orig_count[i] / ranks_per_file;
			offset[i] = local_rank * chunk_size;
			if (local_rank == ranks_per_file - 1) {
				// Last rank gets the rest of the particles
				chunk_size = orig_count[i] - offset[i];
			}

			npart_local_rank[NALL] += (npart_local_rank[i] = chunk_size);
		}

		printf("rank %d: processing file %s (%d/%d), particles: %lld\n", world_rank, name.c_str(), f + 1, nfiles, npart_local_rank[NALL]);
		printf("rank %d: particles per type: ", world_rank);
		for (int i = 0; i < NTYPES; i++) {
			printf("%lld ", npart_local_rank[i]);
		}
		printf("\n");

		// ---------------------------------------------------------
		//
		// FIRST LOOP
		// loop over the sub-files
		//
		// ---------------------------------------------------------   
		//{

		//nparts_t npart_local;                     // will contains the nr. of particles for this file
		//std::vector<char>     name_vec(strlen(basename) + 10);
		//std::string name;     // the name of this file


		// open the fth sub-file
		//
		//if (nfiles > 1)
		//	//name = std::format("{}.{}", basename, f);
		//	name = std::string(basename) + "." + std::to_string(f);
		//else
		//	//name = std::format("{}", basename);
		//	name = std::string(basename);

		file = fopen(name.c_str(), "rb");
		if (file == NULL)
		{
			// that should not happen, but deal with it
			//
			printf("I cannot open the file %s\n", name.c_str());
			exit(-1);
		}

		// get how many particles are in this file
		//
		//get_particles_num(file, &npart_local, &npart_all);
		//printf("processing file %d over %d:\n", f, nfiles);


		// -----------------------------------------------------------------
		//
		// SECOND LOOP
		// loop over the "components"
		//
		// ...........................      
		// A "component" is a list of physical quantities that you want
		// to map onto the grid.
		// Each component refers to particles of one specific type
		// ( there are 6 possible types: 0,..,5 that are gas, dark-matter
		// (3 types of dm), stars and black holes.
		//------------------------------------------------------------------

		for (int type = 0; type < NTYPES; type++) //ptypes
		{
			if (npart_local_rank[type] == 0) {
				continue;
			}

			//pos_t* pos;
			int      block_idx;
			//int      type;

			// get the type to which this component refers to
			//
			//type = c;// g_extract_data.components_list[c].type;

			// -----------------------------------------------------------
			//
			// THIRD LOOP
			// loop over the list of quantities for this component
			//
			// -----------------------------------------------------------			

			for (int q = 0; q < knownblocks.size(); q++)
			{
				//float_t* data = NULL;
				float_t  mass = 0;

				char label[4];
				get_Tab_IO_Label((iofields)q, label);

				//printf("\t\tprocessing field %s\n", label);

				if (strncmp(label, "MASS", 4) == 0)
				{
					// special treatment for the mass table
					//
					mass = snaphead.mass[type];

				}

				if (mass == 0)
				{
					// set the file position at the right place
					//
					block_idx = position_and_skip_particles(knownblocks.data(), knownblocks.size(), label, type, &orig_count, file);

					if (block_idx < 0) { //not found
						continue;
					}

					GadgetData gd;
					gd.data_name = std::string(label, 4);
					gd.num_elements = knownblocks[block_idx].fact;
					gd.size_of_elem = knownblocks[block_idx].datasize;
					gd.data.resize(knownblocks[block_idx].fact * knownblocks[block_idx].datasize * npart_local_rank[type]);

					//allocate the memory needed for the positions
					//data = (float_t*)malloc(npart_local[type] * knownblocks[block_idx].fact * knownblocks[block_idx].datasize);
					//memset(data, -1, npart_local[type] * knownblocks[block_idx].fact * knownblocks[block_idx].datasize);
					//load the positions for this type

					FSEEKO(file, knownblocks[block_idx].fact * knownblocks[block_idx].datasize * offset[type], SEEK_CUR);
					fread((void*)gd.data.data(), knownblocks[block_idx].fact * knownblocks[block_idx].datasize, npart_local_rank[type], file);

					//printf("data type: %d - %d\n", knownblocks[block_idx].fact, knownblocks[block_idx].datasize);

					if (gadget_datas[type][q].data.size() == 0) {
						gadget_datas[type][q].reserve_memory(knownblocks[block_idx].fact * knownblocks[block_idx].datasize * npart_local_rank[type]);
					}

					gadget_datas[type][q].merge_data(gd, npart_local_rank[type]);
				}
				else {
					// allocate the memory needed for the positions
					//data = (float_t*)malloc(npart_local[type] * sizeof(float_t));

					GadgetData gd;
					gd.data_name = std::string(label, 4);
					gd.num_elements = 1;
					gd.size_of_elem = sizeof(float_t);
					gd.data.resize(sizeof(float_t) * npart_local_rank[type]);

					float_t* data = (float_t*)gd.data.data();

#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
					for (int i = 0; i < npart_local_rank[type]; i++) {
						data[i] = mass;
					}

					if (gadget_datas[type][q].data.size() == 0) {
						gadget_datas[type][q].reserve_memory(knownblocks[block_idx].fact * knownblocks[block_idx].datasize * npart_local_rank[type]);
					}

					gadget_datas[type][q].merge_data(gd, npart_local_rank[type]);
				}

			}  // closes the loop over the list of quantities for the component c

			//printf("\n");
			//free(pos);

		}  // closes the loop over the array of components

		//printf("\n");
		fclose(file);
		//}  // closes the loop over the sub-files    

		return 0;
	}

	int read_multiple_files_per_rank(std::string& basename, int world_rank, int world_size)
	{
		int ranks_per_process = nfiles / world_size;
		int start_ranks = world_rank * ranks_per_process;
		int num_ranks = (world_rank == world_size - 1) ? nfiles - start_ranks : ranks_per_process;

		// find npart_local_rank
		npart_local_rank[NALL] = 0;
		for (int f = start_ranks; f < start_ranks + num_ranks; f++)
		{
			nparts_t npart_local;                     // will contains the nr. of particles for this file
			//std::vector<char>     name_vec(strlen(basename) + 10);
			std::string name;     // the name of this file

			// open the fth sub-file
			//
			if (nfiles > 1)
				//name = std::format("{}.{}", basename, f);
				name = std::string(basename) + "." + std::to_string(f);
			else
				//name = std::format("{}", basename);
				name = std::string(basename);

			FILE* file = fopen(name.c_str(), "rb");
			if (file == NULL)
			{
				// that should not happen, but deal with it
				//
				printf("I cannot open the file %s\n", name.c_str());
				exit(-1);
			}

			// get how many particles are in this file
			//
			//get_particles_num(file, &npart_local, &npart_all);
			io_header snaphead;
			gadget_simple_read_snaphead(file, &snaphead);
			get_particles_num(&snaphead, &npart_local, &npart_all);

			for (int i = 0; i < NTYPES; i++) {
				npart_local_rank[i] += npart_local[i];
				npart_local_rank[NALL] += npart_local[i];
			}

			fclose(file);
		}

		// ---------------------------------------------------------
		//
		// FIRST LOOP
		// loop over the sub-files
		//
		// ---------------------------------------------------------   
		for (int f = start_ranks; f < start_ranks + num_ranks; f++)
		{

			nparts_t npart_local;                     // will contains the nr. of particles for this file
			//std::vector<char>     name_vec(strlen(basename) + 10);
			std::string name;     // the name of this file


			// open the fth sub-file
			//
			if (nfiles > 1)
				//name = std::format("{}.{}", basename, f);
				name = std::string(basename) + "." + std::to_string(f);
			else
				//name = std::format("{}", basename);
				name = std::string(basename);

			FILE* file = fopen(name.c_str(), "rb");
			if (file == NULL)
			{
				// that should not happen, but deal with it
				//
				printf("I cannot open the file %s\n", name.c_str());
				exit(-1);
			}

			// get how many particles are in this file
			//
			//get_particles_num(file, &npart_local, &npart_all);
			io_header snaphead;
			gadget_simple_read_snaphead(file, &snaphead);
			get_particles_num(&snaphead, &npart_local, &npart_all);


			printf("rank %d: processing file %s (%d/%d), particles: %lld\n", world_rank, name.c_str(), f + 1, nfiles, npart_local[NALL]);
			printf("rank %d: particles per type: ", world_rank);
			for (int i = 0; i < NTYPES; i++) {
				printf("%lld ", npart_local[i]);
			}
			printf("\n");


			// -----------------------------------------------------------------
			//
			// SECOND LOOP
			// loop over the "components"
			//
			// ...........................      
			// A "component" is a list of physical quantities that you want
			// to map onto the grid.
			// Each component refers to particles of one specific type
			// ( there are 6 possible types: 0,..,5 that are gas, dark-matter
			// (3 types of dm), stars and black holes.
			//------------------------------------------------------------------

			for (int type = 0; type < NTYPES; type++) //ptypes
			{
				if (npart_local[type] == 0) {
					continue;
				}

				//pos_t* pos;
				int      block_idx;
				//int      type;

				// get the type to which this component refers to
				//
				//type = c;// g_extract_data.components_list[c].type;

				// -----------------------------------------------------------
				//
				// THIRD LOOP
				// loop over the list of quantities for this component
				//
				// -----------------------------------------------------------			

				for (int q = 0; q < knownblocks.size(); q++)
				{
					//float_t* data = NULL;
					float_t  mass = 0;

					char label[4];
					get_Tab_IO_Label((iofields)q, label);

					//printf("\t\tprocessing field %s\n", label);

					if (strncmp(label, "MASS", 4) == 0)
					{
						// special treatment for the mass table
						//
						mass = snaphead.mass[type];

					}

					if (mass == 0)
					{
						// set the file position at the right place
						//
						block_idx = position_and_skip_particles(knownblocks.data(), knownblocks.size(), label, type, &npart_local, file);

						if (block_idx < 0) { //not found
							continue;
						}

						GadgetData gd;
						gd.data_name = std::string(label, 4);
						gd.num_elements = knownblocks[block_idx].fact;
						gd.size_of_elem = knownblocks[block_idx].datasize;
						gd.data.resize(knownblocks[block_idx].fact * knownblocks[block_idx].datasize * npart_local[type]);

						//allocate the memory needed for the positions
						//data = (float_t*)malloc(npart_local[type] * knownblocks[block_idx].fact * knownblocks[block_idx].datasize);
						//memset(data, -1, npart_local[type] * knownblocks[block_idx].fact * knownblocks[block_idx].datasize);
						//load the positions for this type

						//format 2 - skip 4 bytes
						int blksize = 0;
						fread(&blksize, sizeof(int), 1, file);

						fread((void*)gd.data.data(), knownblocks[block_idx].fact * knownblocks[block_idx].datasize, npart_local[type], file);

						//printf("data type: %d - %d\n", knownblocks[block_idx].fact, knownblocks[block_idx].datasize);

						if (gadget_datas[type][q].data.size() == 0) {
							gadget_datas[type][q].reserve_memory(knownblocks[block_idx].fact * knownblocks[block_idx].datasize * npart_local_rank[type]);
						}

						gadget_datas[type][q].merge_data(gd, npart_local[type]);
					}
					else {
						// allocate the memory needed for the positions
						//data = (float_t*)malloc(npart_local[type] * sizeof(float_t));

						GadgetData gd;
						gd.data_name = std::string(label, 4);
						gd.num_elements = 1;
						gd.size_of_elem = sizeof(float_t);
						gd.data.resize(sizeof(float_t) * npart_local[type]);

						float_t* data = (float_t*)gd.data.data();

#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
						for (int i = 0; i < npart_local[type]; i++) {
							data[i] = mass;
						}

						if (gadget_datas[type][q].data.size() == 0) {
							gadget_datas[type][q].reserve_memory(knownblocks[block_idx].fact * knownblocks[block_idx].datasize * npart_local_rank[type]);
						}

						gadget_datas[type][q].merge_data(gd, npart_local[type]);
					}

				}  // closes the loop over the list of quantities for the component c

				//printf("\n");
				//free(pos);

			}  // closes the loop over the array of components

			//printf("\n");
			fclose(file);
		}  // closes the loop over the sub-files    

		return 0;
	}

	int gadget_simple_swap_file = 8;
	void gadget_simple_swap_Nbyte(char* data, int n, int m)
	{
		int i, j;
		char old_data[16];

		if (gadget_simple_swap_file != 8)
		{
			for (j = 0; j < n; j++)
			{
				memcpy(&old_data[0], &data[j * m], m);
				for (i = 0; i < m; i++)
				{
					data[j * m + i] = old_data[m - i - 1];
				}
			}
		}
	}

	void gadget_simple_swap_snaphead(io_header* snaphead)
	{
		gadget_simple_swap_Nbyte((char*)&snaphead->npart, 6, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->mass, 6, 8);
		gadget_simple_swap_Nbyte((char*)&snaphead->time, 1, 8);
		gadget_simple_swap_Nbyte((char*)&snaphead->redshift, 1, 8);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_sfr, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_feedback, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->npartTotal, 6, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_cooling, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->num_files, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->BoxSize, 1, 8);
		gadget_simple_swap_Nbyte((char*)&snaphead->Omega0, 1, 8);
		gadget_simple_swap_Nbyte((char*)&snaphead->OmegaLambda, 1, 8);
		gadget_simple_swap_Nbyte((char*)&snaphead->HubbleParam, 1, 8);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_stellarage, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_metals, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->npartTotalHighWord, 6, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_entropy_instead_u, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_doubleprecision, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->flag_ic_info, 1, 4);
		gadget_simple_swap_Nbyte((char*)&snaphead->lpt_scalingfactor, 1, 4);
	}

	void gadget_simple_read_snaphead(FILE* fd, io_header* snaphead)
	{		
		//char buf[200], buf1[200];
		int dummy;

		//if (All.ICFormat == 1 || All.ICFormat == 2)
		//{
		//	if (All.ICFormat == 2)
		//	{
		fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		gadget_simple_swap_file = dummy;
#endif
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		//}

		fread(&dummy, sizeof(dummy), 1, fd);
		//#ifdef AUTO_SWAP_ENDIAN_READIC
		//					if (All.ICFormat == 1)
		//					{
		//						if (dummy == 256)
		//							gadget_simple_swap_file = 8;
		//						else
		//							gadget_simple_swap_file = 1;
		//					}
		//#endif
		fread(snaphead, sizeof(io_header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		gadget_simple_swap_snaphead(snaphead);
#endif
		//#ifdef PATCH_IO
		//					printf("Task %d: WARNING, PATCHING HEADER FOR MagneticumDM simulation !!!\n", ThisTask);
		//					snaphead->mass[2] = snaphead->mass[1];
		//					snaphead->mass[1] = snaphead->mass[0];
		//					snaphead->mass[0] = 0;
		//					snaphead->npart[2] = snaphead->npart[1];
		//					snaphead->npart[1] = snaphead->npart[0];
		//					snaphead->npart[0] = 0;
		//					snaphead->npartTotal[2] = snaphead->npartTotal[1];
		//					snaphead->npartTotal[1] = snaphead->npartTotal[0];
		//					snaphead->npartTotal[0] = 0;
		//					snaphead->npartTotalHighWord[2] = snaphead->npartTotalHighWord[1];
		//					snaphead->npartTotalHighWord[1] = snaphead->npartTotalHighWord[0];
		//					snaphead->npartTotalHighWord[0] = 0;
		//#endif
		fread(&dummy, sizeof(dummy), 1, fd);
		//}
		//fclose(fd);

		//#ifdef HAVE_HDF5
		//				if (All.ICFormat == 3)
		//					read_header_attributes_in_hdf5(buf);
		//#endif
	}

	int check_number_of_files(std::string& basename)
	{
#if 0
		int   nfiles = 1;
		FILE* filein = NULL;
		std::string fname;// = (char*)malloc(strlen(wdir) + strlen(basename) + 10);
		io_header* snaphead = (snaphead == NULL ? (io_header*)calloc(1, sizeof(io_header)) : snaphead);

		/*
		 * first, determine whether it seems that we are in a multi-file case;
		 * we check the existence of a files named as basename
		 */
		 //fname = std::format("{}{}", wdir, basename);
		fname = basename;
		if ((filein = fopen(fname.c_str(), "rb")) == NULL)
		{
			// a file named as basename does NOT exist;
			// that meaans that we are in a multi-file case
			// and hence we check the existance of
			// $basename.0 file
			//
			//fname = std::format("{}{}.0", wdir, basename);
			fname = basename + std::string(".0");
			if ((filein = fopen(fname.c_str(), "rb")) == NULL)
			{
				// ooops: $basename.0 does not exist either;
				// let's signal an error and exit
				//
				printf("[io] unable to find the file under %s/ "
					"both as %s and as %s.0\n",
					basename.c_str(), basename.c_str(), basename.c_str());
				//free(fname);
				if (snaphead == NULL) free(snaphead);
				return 0;
			}

			// get the snapshot snaphead
			//
			seek_block(filein, "HEAD");
			fread(snaphead, sizeof(io_header), 1, filein);
			gadget_simple_swap_snaphead(snaphead);

			// now let's check how many files $basename.x we find
			//
			do
			{
				// close the file open before
				fclose(filein);
				// check for and additional file
				//fname = std::format("{}{}.{}", wdir, basename, nfiles);
				fname = basename + std::string(".") + std::to_string(nfiles);
				nfiles += ((filein = fopen(fname.c_str(), "rb")) != NULL);
			} while (filein != NULL);

		}
		else
		{
			//
			// the file $basename exists; hence this snapshot is
			// made of just a single file
			//

			seek_block(filein, "HEAD");
			fread(snaphead, sizeof(io_header), 1, filein);
			gadget_simple_swap_snaphead(snaphead);

			fclose(filein);
		}

		if (snaphead->num_files != nfiles)
		{
			printf("There is a mismatch between the nr. of files found (%d) "
				"and the number of files specified in the snaphead (%d)\n",
				nfiles, snaphead->num_files);
			//nfiles = -1;
		}

		//free(fname);
		if (snaphead == NULL) free(snaphead);

		return nfiles;
#endif

		FILE* fd;
		char buf[200], buf1[200];
		int dummy;
		//int ICFormat = 2;

		io_header snaphead;

		sprintf(buf, "%s.%d", basename.c_str(), 0);
		sprintf(buf1, "%s", basename.c_str());

		//if (All.ICFormat == 3)
		//{
		//	sprintf(buf, "%s.%d.hdf5", basename.c_str(), 0);
		//	sprintf(buf1, "%s.hdf5", basename.c_str());
		//}

		//if (world_rank == 0)
		//printf("rank %d: Trying to read file %s in format %d\n", world_rank, buf, 2);

		//#ifndef  HAVE_HDF5
		//		if (All.ICFormat == 3)
		//		{
		//			if (ThisTask == 0)
		//				printf("Code wasn't compiled with HDF5 support enabled!\n");
		//			exit(0);
		//		}
		//#endif

		snaphead.num_files = 0;

		//if (ThisTask == 0)
		//{
		if ((fd = fopen(buf, "rb")))
		{
//			//if (All.ICFormat == 1 || All.ICFormat == 2)
//			//{
//			//	if (All.ICFormat == 2)
//			//	{
//			fread(&dummy, sizeof(dummy), 1, fd);
//#ifdef AUTO_SWAP_ENDIAN_READIC
//			gadget_simple_swap_file = dummy;
//#endif
//			fread(&dummy, sizeof(dummy), 1, fd);
//			fread(&dummy, sizeof(dummy), 1, fd);
//			fread(&dummy, sizeof(dummy), 1, fd);
//			//}
//
//			fread(&dummy, sizeof(dummy), 1, fd);
//			//#ifdef AUTO_SWAP_ENDIAN_READIC
//			//					if (All.ICFormat == 1)
//			//					{
//			//						if (dummy == 256)
//			//							gadget_simple_swap_file = 8;
//			//						else
//			//							gadget_simple_swap_file = 1;
//			//					}
//			//#endif
//			fread(&snaphead, sizeof(io_header), 1, fd);
//#ifdef AUTO_SWAP_ENDIAN_READIC
//			swap_header();
//#endif
//			//#ifdef PATCH_IO
//			//					printf("Task %d: WARNING, PATCHING HEADER FOR MagneticumDM simulation !!!\n", ThisTask);
//			//					snaphead->mass[2] = snaphead->mass[1];
//			//					snaphead->mass[1] = snaphead->mass[0];
//			//					snaphead->mass[0] = 0;
//			//					snaphead->npart[2] = snaphead->npart[1];
//			//					snaphead->npart[1] = snaphead->npart[0];
//			//					snaphead->npart[0] = 0;
//			//					snaphead->npartTotal[2] = snaphead->npartTotal[1];
//			//					snaphead->npartTotal[1] = snaphead->npartTotal[0];
//			//					snaphead->npartTotal[0] = 0;
//			//					snaphead->npartTotalHighWord[2] = snaphead->npartTotalHighWord[1];
//			//					snaphead->npartTotalHighWord[1] = snaphead->npartTotalHighWord[0];
//			//					snaphead->npartTotalHighWord[0] = 0;
//			//#endif
//			fread(&dummy, sizeof(dummy), 1, fd);
//			//}
			gadget_simple_read_snaphead(fd, &snaphead);
			printf("snaphead.redshift: %f, snaphead.HubbleParam: %f\n", snaphead.redshift, snaphead.HubbleParam);

			fclose(fd);
//			//#ifdef HAVE_HDF5
//			//				if (All.ICFormat == 3)
//			//					read_header_attributes_in_hdf5(buf);
//			//#endif
			
		}
		//}

//#ifdef AUTO_SWAP_ENDIAN_READIC
//		MPI_Bcast(&gadget_simple_swap_file, 1, MPI_INT, 0, MYMPI_COMM_WORLD);
//#endif
//		MPI_Bcast(&snaphead, sizeof(snaphead), MPI_BYTE, 0, MYMPI_COMM_WORLD);

		if (snaphead.num_files > 0)
			return snaphead.num_files;

		//if (ThisTask == 0)
		//{
		if ((fd = fopen(buf1, "rb")))
		{
//			//if (All.ICFormat == 1 || All.ICFormat == 2)
//			//{
//			//	if (All.ICFormat == 2)
//			//	{
//			fread(&dummy, sizeof(dummy), 1, fd);
//#ifdef AUTO_SWAP_ENDIAN_READIC
//			gadget_simple_swap_file = dummy;
//#endif
//			fread(&dummy, sizeof(dummy), 1, fd);
//			fread(&dummy, sizeof(dummy), 1, fd);
//			fread(&dummy, sizeof(dummy), 1, fd);
//			//}
//
//			fread(&dummy, sizeof(dummy), 1, fd);
//			//#ifdef AUTO_SWAP_ENDIAN_READIC
//			//					if (All.ICFormat == 1)
//			//					{
//			//						if (dummy == 256)
//			//							gadget_simple_swap_file = 8;
//			//						else
//			//							gadget_simple_swap_file = 1;
//			//					}
//			//#endif
//			fread(&snaphead, sizeof(io_header), 1, fd);
//#ifdef AUTO_SWAP_ENDIAN_READIC
//			swap_header();
//#endif
//			//#ifdef PATCH_IO
//			//					printf("Task %d: WARNING, PATCHING HEADER FOR MagneticumDM simulation !!!\n", ThisTask);
//			//					snaphead->mass[2] = snaphead->mass[1];
//			//					snaphead->mass[1] = snaphead->mass[0];
//			//					snaphead->mass[0] = 0;
//			//					snaphead->npart[2] = snaphead->npart[1];
//			//					snaphead->npart[1] = snaphead->npart[0];
//			//					snaphead->npart[0] = 0;
//			//					snaphead->npartTotal[2] = snaphead->npartTotal[1];
//			//					snaphead->npartTotal[1] = snaphead->npartTotal[0];
//			//					snaphead->npartTotal[0] = 0;
//			//					snaphead->npartTotalHighWord[2] = snaphead->npartTotalHighWord[1];
//			//					snaphead->npartTotalHighWord[1] = snaphead->npartTotalHighWord[0];
//			//					snaphead->npartTotalHighWord[0] = 0;
//			//#endif
//			fread(&dummy, sizeof(dummy), 1, fd);
//			//}
			gadget_simple_read_snaphead(fd, &snaphead);
			printf("snaphead.redshift: %f, snaphead.HubbleParam: %f\n", snaphead.redshift, snaphead.HubbleParam);

			fclose(fd);

			//#ifdef HAVE_HDF5
			//				if (All.ICFormat == 3)
			//					read_header_attributes_in_hdf5(buf1);
			//#endif

			snaphead.num_files = 1;
		}
		//}

		//#ifdef AUTO_SWAP_ENDIAN_READIC
		//		MPI_Bcast(&gadget_simple_swap_file, 1, MPI_INT, 0, MYMPI_COMM_WORLD);
		//#endif
		//		MPI_Bcast(&snaphead, sizeof(snaphead), MPI_BYTE, 0, MYMPI_COMM_WORLD);

		if (snaphead.num_files > 0)
			return snaphead.num_files;

		//if (ThisTask == 0)
		//{
		//	printf("\nCan't find initial conditions file.");
		//	printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
		//	fflush(stdout);
		//}

		//exit(0);
		return 0;
	}

	int seek_block(FILE* file, const blockname_t name)
		/*
		 * This function get to the block specified by name
		 * in a format-2 file pointed by *file.
		 *
		 * RETURN VALUE:
		 * 0 if the seeking is not successful, otherwise
		 * the 4-bytes long tag at the begin of the data block
		 * (i.e. the block's length in bytes).
		 * The file position is set to the begin of the data
		 * block, right after the length tag.
		 *
		 */
	{
#if 0
#define FBSKIP  {ret = fread(&blksize,sizeof(int),1,file);}

		if (file == NULL)
			return -1;

		//head_t head;
		size_t ret = 0;
		unsigned int blocksize = 0, blksize;
		char blocklabel[5] = { "    " };

		rewind(file);

		//ret = fread(&head, sizeof(head_t), 1, file);
		{
			FBSKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
			gadget_simple_swap_file = blksize;
			gadget_simple_swap_Nbyte((char*)&blksize, 1, 4);
#endif
			if (blksize != 8)
			{
				printf("Incorrect Format (blksize=%u)!\n", blksize);
				exit(1891);
			}
			ret = fread(blocklabel, 4 * sizeof(char), 1, file);
			ret = fread(&blocksize, sizeof(int), 1, file);
#ifdef AUTO_SWAP_ENDIAN_READIC
			gadget_simple_swap_Nbyte((char*)&blocksize, 1, 4);
#endif
			FBSKIP;
		}

		while ((!feof(file)) && (strncmp(blocklabel, name, 4) != 0) && (ret == 1))
		{
			FSEEKO(file, (off_t)blocksize, SEEK_CUR);
			if (!feof(file)) {
				//ret = fread(&head, sizeof(head_t), 1, file);
				FBSKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
				gadget_simple_swap_file = blksize;
				gadget_simple_swap_Nbyte((char*)&blksize, 1, 4);
#endif
				//if (blksize != 8)
				//{
				//	printf("Incorrect Format (blksize=%u)!\n", blksize);
				//	exit(1891);
				//}
				ret = fread(blocklabel, 4 * sizeof(char), 1, file);
				ret = fread(&blocksize, sizeof(int), 1, file);
#ifdef AUTO_SWAP_ENDIAN_READIC
				gadget_simple_swap_Nbyte((char*)&blocksize, 1, 4);
#endif
				FBSKIP;
			}
		}

		//printf("Block pos: %lld, %lld\n", ftell(file), _ftelli64(file));

		if (!feof(file) && (ret == 1))
		{
			int   tag1, tag2;
			ret = fread(&tag1, sizeof(int), 1, file);
			//off_t pos = FTELLO(file);
			//ret = FSEEKO(file, tag1, SEEK_CUR);
			//ret = fread(&tag2, sizeof(int), 1, file);
			//ret = FSEEKO(file, pos, SEEK_SET);
			//if ((tag1 != tag2)) {
			//	printf("[SEEKBLOCK] error in tags of block %s: %d vs %d\n",
			//		name, tag1, tag2);
			//	return 0;
			//}

			//if (tag1 == 0)
			//	printf("[SEEKBLOCK] error 0-valued tags for block %s\n", name);

			////printf("Block pos: %lld, %lld\n", ftell(file), _ftelli64(file));

			//return tag1;
//			FBSKIP;
//#ifdef AUTO_SWAP_ENDIAN_READIC
//			//gadget_simple_swap_file = blksize;
//			gadget_simple_swap_Nbyte((char*)&blksize, 1, 4);
//#endif
			return gadget_simple_swap_file;
		}
		else
			return -1;
#else

#define FBSKIP  {ret = fread(&blksize,sizeof(int),1,file);}

		unsigned int blocksize = 0, blksize;
		char blocklabel[5] = { "    " };
		size_t ret = 1;

		rewind(file);

		while (!feof(file) && blocksize == 0 && ret == 1)
		{
			FBSKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
			gadget_simple_swap_file = blksize;
			gadget_simple_swap_Nbyte((char*)&blksize, 1, 4);
#endif
			if (blksize != 8)
			{
				printf("Incorrect Format (blksize=%u)!\n", blksize);
				exit(1891);
			}
			else
			{
				ret = fread(blocklabel, 4 * sizeof(char), 1, file);
				ret = fread(&blocksize, sizeof(int), 1, file);
#ifdef AUTO_SWAP_ENDIAN_READIC
				gadget_simple_swap_Nbyte((char*)&blocksize, 1, 4);
#endif
				/*
					printf("Searching <%c%c%c%c>, found Block <%s> with %d bytes\n",
					label[0],label[1],label[2],label[3],blocklabel,blocksize);
				*/
				FBSKIP;
				if (strncmp(name, blocklabel, 4) != 0)
				{
#ifdef _WIN32
					_fseeki64(file, (__int64)blocksize, 1);
#else
					fseek(file, blocksize, 1);
#endif
					blocksize = 0;
				}
			}
		}
		if (feof(file))
		{
			printf("Block '%c%c%c%c' not found !\n", name[0], name[1], name[2], name[3]);
			fflush(stdout);
			exit(1890);
		}

		return (blocksize == 0) ? -1 : 1;
#endif

	}


	int get_particles_num(io_header *snaphead, nparts_t* nparts_file, nparts_t* nparts_all)
	{

		//int ret = seek_block(file, "HEAD");

		//if (ret == 0)
		//{
		//	printf("unable to find the HEAD block in snapshot files\n");
		//	return -1;
		//}

		//io_header snaphead;

		//ret = fread(&snaphead, sizeof(io_header), 1, file);
		//gadget_simple_swap_snaphead(&snaphead);

		if (nparts_file != NULL)
		{
			(*nparts_file)[NALL] = 0;
			for (int i = 0; i < NTYPES; i++)
				(*nparts_file)[NALL] += ((*nparts_file)[i] = snaphead->npart[i]);
		}

		if (nparts_all != NULL)
		{
			(*nparts_all)[NALL] = 0;
			for (int i = 0; i < NTYPES; i++)
				(*nparts_all)[NALL] += ((*nparts_all)[i] =
					(((long long)snaphead->npartTotalHighWord[i]) << 32) + snaphead->npartTotal[i]);
		}

		return snaphead->num_files;
	}


	int position_and_skip_particles(blockdata_t* knownblocks, int n_knownblocks, const blockname_t name, int current_type, nparts_t* npart, FILE* file)
		//
		// this routine sets the file position so that to skip the particles of type < current_type,
		// for the block with name <name>
		// it returns the block index in the array of blocks metadata
		//
	{
		// retrieve the index of "POS " block in the array of blocks metadata
		//
		int idx = get_knownblock_idx(knownblocks, n_knownblocks, name);
		if (idx < 0)
		{
			// some error occurred (probably the wrong name has been passed)
			//
			return idx;
		}

		// positioning at the begin of the block
		int ret = seek_block(file, name);
		if (ret < 0)
		{
			// some error occurred (either a block with that name has not been found
			// or the block is corrupted )
			//
			return ret;
		}

		// skip the positions for types < current type;
		// to_skip contains the cumulative nr. of particles
		// having types < current type
		uint to_skip = 0;
		for (int t = 0; t < current_type; t++)
			to_skip += (*npart)[t];

		// actual positioning in the file
		//
		int64_t size = to_skip *         // how many particles are to be skipped     
			knownblocks[idx].fact *       // how many data per particle
			knownblocks[idx].datasize;    // the single data size

		// repositioning
		//
		FSEEKO(file, size, SEEK_CUR);

		// return the block index
		//
		return idx;
	}

	int get_knownblock_idx(blockdata_t* knownblocks, int n_knownblocks, const blockname_t name)
	{
		// identify the b-th component among the known blocks
		//
		int idx = 0;
		while ((idx < n_knownblocks) &&
			(strncmp(name, knownblocks[idx].name, strlen(knownblocks[idx].name)) != 0))
			idx++;
		if (idx == n_knownblocks)
			idx = -1;

		return idx;
	}

	namespace io {

		size_t get_local_num_particles()
		{
			return npart_local_rank[NALL];
		}

		size_t get_global_num_particles()
		{
			return npart_all[NALL];
		}

		void print_CPU_steps()
		{
			printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		int get_particle_type(uint64_t id) {
			size_t temp_count = 0;
			for (int pt = 0; pt < NTYPES; pt++) {
				temp_count += npart_local_rank[pt];
				if (id < temp_count) {
					return pt;
				}
			}

			return 0;
		}

		uint64_t get_particle_type_offset(int pt) {
			uint64_t offset = 0;
			for (int o = 0; o < pt; o++) {
				offset += npart_local_rank[o];
			}

			return offset;
		}

		void get_particle_position(uint64_t id, double* pos)
		{
			int pt = get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);
			gadget_datas[pt][IO_POS].get_dvalue3(id - offset, pos);
		}

		float get_particle_norm_value(int blocknr, uint64_t id)
		{
			int pt = get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			//int
			if (blocknr == IO_ID ||
				blocknr == IO_TRUENGB ||
				blocknr == IO_MAINHALO ||
				blocknr == IO_BHPROGS)
			{
				RETURN_NORM_VALUE(gadget_datas[pt][blocknr].get_ivalue(id - offset));
			}


			if (gadget_datas[pt][blocknr].get_nelem() == 1) {
				RETURN_NORM_VALUE(gadget_datas[pt][blocknr].get_fvalue(id - offset));
			}

			if (gadget_datas[pt][blocknr].get_nelem() == 3) {
				float v[3];
				gadget_datas[pt][blocknr].get_fvalue3(id - offset, v);
				RETURN_NORM_VECTOR3(v);
			}

			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value)
		{
			int pt = get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			//int
			if (blocknr == IO_ID ||
				blocknr == IO_TRUENGB ||
				blocknr == IO_MAINHALO ||
				blocknr == IO_BHPROGS)
			{
				RETURN_ORIG_VALUE(gadget_datas[pt][blocknr].get_ivalue(id - offset));
			}

			if (gadget_datas[pt][blocknr].get_nelem() == 1) {
				RETURN_ORIG_VALUE(gadget_datas[pt][blocknr].get_fvalue(id - offset));
			}

			if (gadget_datas[pt][blocknr].get_nelem() == 3) {
				float v[3];
				gadget_datas[pt][blocknr].get_fvalue3(id - offset, v);
				RETURN_ORIG_VECTOR3(v);
			}

			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id)
		{
			int pt = get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			return gadget_datas[pt][blocknr].get_nelem();
		}

		//double get_particle_radius(uint64_t id) {
		//	double hsml = get_particle_hsml(id);
		//	double radius = 2.0 * hsml;

		//	if (std::isnan(radius))
		//		radius = 0;

		//	return radius;
		//}

		double get_particle_hsml(uint64_t id) {
			int pt = get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			return gadget_datas[pt][IO_HSML].get_dvalue(id - offset);
		}

		double get_particle_mass(uint64_t id) {
			int pt = get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			return gadget_datas[pt][IO_MASS].get_dvalue(id - offset);
		}

		int get_particle_rho_blocknr() {
			return IO_RHO;
		}

		double get_particle_rho(uint64_t id) {
			int pt = get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			return gadget_datas[pt][IO_RHO].get_dvalue(id - offset);
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

		void gadget_init_lib(std::string& fname, int world_rank, int world_size)
		{
			steps_time[0] = omp_get_wtime();

			fill_known_blocks();

			nfiles = check_number_of_files(fname);

			// working dir is the current one
			if (nfiles <= 0)
			{
				if (nfiles == 0)
					printf("I'm sorry to be unable to find any snapshot "
						"neither as ./%1$s nor as ./%1$s.x", fname);

				//free(basename);
				exit(-1);
			}

			if (world_rank == 0) {
				printf("Total number of files: %d\n", nfiles);
			}

			memset(npart_local_rank, 0, sizeof(nparts_t)); // initialize the local particles count to 0

			if (nfiles < world_size) {
				read_part_file_per_rank(fname, world_rank, world_size);
			}
			else
			{
				read_multiple_files_per_rank(fname, world_rank, world_size);
			}

			steps_time[1] = omp_get_wtime();
		}

		void gadget_finish_lib()
		{
		}
		//void gadget_set_parameter(int ICFormat, int SnapFormat, int NumFilesWrittenInParallel, int MaxMemSize, double BufferSize, double PartAllocFactor, int TotBHs)
		//{
		//}
		//void gadget_read_parameter_file(char* fname, char* tag[], void** addr, int* id, int nt)
		//{
		//}
		//void gadget_read_ic(std::string& fname)
		//{
		//	gadget_simple::extract(fname);
		//}
		//void gadget_mymalloc_init()
		//{
		//}

		//void gadget_set_units()
		//{
		//}

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
			types_and_blocks.resize(NTYPES * GADGET_MAX_BLOCKS, 0);
			for (int type = 0; type < NTYPES; type++) {
				for (int blocknr = 0; blocknr < GADGET_MAX_BLOCKS; blocknr++) {
					if (gadget_datas[type][blocknr].get_count() > 0) {
						types_and_blocks[NTYPES * blocknr + type] = 1;
					}
				}
			}
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
			for (int type = 0; type < NTYPES; type++) {
				printf("Type: %s (%d)\n", gadget_get_type_name(type).c_str(), type);
				for (int blocknr = 0; blocknr < GADGET_MAX_BLOCKS; blocknr++) {

					if (gadget_datas[type][blocknr].get_count() > 0) {
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
			for (int type = 0; type < NTYPES; type++) {
				printf("Type: %s (%d)\n", gadget_get_type_name(type).c_str(), type);
				for (int blocknr = 0; blocknr < GADGET_MAX_BLOCKS; blocknr++) {

					if (types_and_blocks[NTYPES * blocknr + type] > 0) {
						std::string buf = gadget_get_dataset_name((iofields)blocknr);
						printf("\t%s (%d)\n", buf.c_str(), blocknr);
					}
				}
			}
		}

		int get_num_types() {
			return NTYPES;
		}
		int get_num_blocks() {
			return GADGET_MAX_BLOCKS;
		}

		double get_redshift()
		{
			return redshift;
		}

		double get_hubble_param()
		{
			return hubble_param;
		}
	}//io
}//gadget