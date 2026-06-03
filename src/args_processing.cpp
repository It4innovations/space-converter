/*
 * Copyright(C) 2023-2026 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
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

#include "args_processing.h"
#include <iostream>

namespace space_converter {

	/**
	 * @brief Display usage information and command-line options for the space_converter tool.
	 * 
	 * This function prints all available command-line arguments, their required parameters,
	 * and format-specific options. It exits the program after displaying the information.
	 */
	void usage()
	{
		std::cout << "./space_converter --data-type [GADGET, GADGET_SIMPLE, CHANGA_TIPSY, CHANGA_NCHILADA, HACC_GENERICIO, HACC_BIN] <options> <args>" << std::endl;

		// === General Options ===
		std::cout << "options:" << std::endl;
		std::cout << "\t--grid-dim X" << std::endl;
		std::cout << "\t--output-path X" << std::endl;
		std::cout << "\t--server X" << std::endl;
		std::cout << "\t--port X" << std::endl;
		std::cout << "\t--info" << std::endl;
		std::cout << "\t--nanovdb" << std::endl;
		//std::cout << "\t--dense" << std::endl;
		std::cout << "\t--dense-file X" << std::endl;
		std::cout << "\t--anim START END STEP" << std::endl;
		std::cout << "\t--anim-merge" << std::endl;		
		std::cout << "\t--raw-particles X" << std::endl;
		// std::cout << "\t--rawpart2vdb" << std::endl;
		std::cout << "\t--export-data TYPE DATASET" << std::endl;
		std::cout << "\t--dense-type X" << std::endl;
		std::cout << "\t--bbox-sphere x y z r" << std::endl;
		std::cout << "\t--simple-density" << std::endl;
		std::cout << "\t--offset-position X Y Z" << std::endl;

		// === Neighbor Search Options ===
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)		
		std::cout << "\t--calc-radius-neigh N              : Calculate radius for N nearest neighbors" << std::endl;
		std::cout << "\t--calc-radius-neigh-rho-kernel X   : Specify density kernel type for neighbor search" << std::endl;
#endif
		std::cout << "\t--calc-radius-neigh-file X         : Load neighbor radii from file" << std::endl;

#ifdef WITH_CUDAKDTREE
		std::cout << "\t--cudakdtree                       : Use CUDA-accelerated KDTree for neighbor search" << std::endl;
		std::cout << "\t--cudakdtree-cpu                   : Use CPU-based CUDA KDTree implementation" << std::endl;
#endif

#ifdef WITH_NANOFLANN
		std::cout << "\t--nanoflann                        : Use nanoflann library for neighbor search" << std::endl;
#endif
	
		// === Advanced Options ===
		std::cout << "\t--sort-by-radius                   : Sort particles by radius" << std::endl;
		std::cout << "\t--sort-by-non-overlap              : Sort particles spatially to avoid atomic operations" << std::endl;		std::cout << "\t--bbox x1 y1 z1 x2 y2 z2           : Define axis-aligned bounding box" << std::endl;
		std::cout << "\t--save-mpi-rank                    : Save MPI rank information                    : Save MPI rank information" << std::endl;

#ifdef WITH_GPU_CUDA
		std::cout << "\t--gpu-cuda                         : Enable GPU acceleration for CUDA-based computations" << std::endl;
#endif		

		// === Format-Specific Arguments ===
		std::cout << "\nGADGET args:" << std::endl;
		std::cout << "\t--param-file X              : GADGET parameter file" << std::endl;
		std::cout << "\t--gadget-file X             : GADGET snapshot file path" << std::endl;
		std::cout << "\t--max-mem-size X            : Maximum memory size (MB)" << std::endl;
		std::cout << "\t--buffer-size X             : Buffer size for particle loading" << std::endl;
		std::cout << "\t--part-alloc-factor X       : Particle allocation factor" << std::endl;
		std::cout << "\t--bh-count X                : Number of black holes in simulation" << std::endl;

		std::cout << "\nGADGET_SIMPLE args:" << std::endl;
		std::cout << "\t--gadget-file X             : GADGET snapshot file path (simplified format)" << std::endl;		
		
		std::cout << "\nCHANGA_TIPSY args:" << std::endl;
		std::cout << "\t--tipsy-file X              : TIPSY format file path" << std::endl;

		std::cout << "\nCHANGA_NCHILADA args:" << std::endl;
		std::cout << "\t--nc-dir X                  : NChilada data directory" << std::endl;

		std::cout << "\nHACC_GENERICIO args:" << std::endl;
		std::cout << "\t--genericio-file X          : GenericIO format file path" << std::endl;
		std::cout << "\t--pos-names x y z           : Position field names in GenericIO file" << std::endl;
		std::cout << "\t--vel-names vx vy vz        : Velocity field names in GenericIO file" << std::endl;

		std::cout << "\nHACC_BIN args:" << std::endl;
		std::cout << "\t--haccbin-file X            : HACC binary format file path" << std::endl;

		exit(0);
	}

	/**
	 * @brief Parse command-line arguments and populate configuration structures.
	 * 
	 * This function processes all command-line arguments provided to the application
	 * and fills the FromCL and SpaceData structures with the parsed values. It handles
	 * various data formats, visualization options, and computation parameters.
	 * 
	 * @param from_cl Reference to FromCL structure for command-line configuration
	 * @param space_data Reference to SpaceData structure for spatial/simulation data
	 * @param argc Number of command-line arguments
	 * @param argv Array of command-line argument strings
	 */
	void parse_args(FromCL& from_cl, common::SpaceData &space_data, int argc, char** argv)
	{
		// Require at least 3 arguments (program name + --data-type + format)
		if (argc < 3) {
			usage();
		}

		// Parse each command-line argument
		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			
			// === Core Configuration ===
			if (arg == "--data-type") {
				// Specify input data format (GADGET, HACC, etc.)
				from_cl.data_type = argv[++i];
			}
			else if (arg == "--grid-dim") {
				// Set the grid resolution/dimension for voxelization
				space_data.bbox_dim = std::stoi(argv[++i]);
			}
			else if (arg == "-o" || arg == "--output-path") {
				// Output directory for generated files
				from_cl.output_path = argv[++i];
			}
			else if (arg == "--server") {
				// Server address for remote communication
				from_cl.server = argv[++i];
			}
			else if (arg == "--port") {
				// Server port for remote communication
				from_cl.port = std::stoi(argv[++i]);
			}
			else if (arg == "--info") {
				// Display dataset information only (no conversion)
				from_cl.info = true;
				from_cl.remote = false;
			}
			// === Output Format Options ===
#ifdef WITH_OPENVDB
			else if (arg == "--nanovdb") {
				// Use NanoVDB format (GPU-friendly VDB variant)
				from_cl.use_nanovdb = true;
			}
#endif
			else if (arg == "--dense-file") {
				// Export dense matrix representation to file
				//from_cl.use_dense2file = true;
				space_data.extracted_dense_type = (common::SpaceData::ExtractedDenseType)std::stoi(argv[++i]);
			}
			
			// === Animation Processing ===
			else if (arg == "--anim") {
				// Process animation sequence (start and end frame numbers)
				if(space_data.anim_type == common::SpaceData::AnimType::eNone)
					space_data.anim_type = common::SpaceData::AnimType::eAllPath;

				space_data.anim_start = std::stoi(argv[++i]);
				space_data.anim_end = std::stoi(argv[++i]);
				space_data.anim_step = std::stoi(argv[++i]);
			}
			else if (arg == "--anim-merge") {
				// Merge all animation frames into a single output
				space_data.anim_type = common::SpaceData::AnimType::eAllMerge;
			}
			
			// === Particle Export Options ===
			else if (arg == "--raw-particles") {
				// Export raw particle data without VDB conversion
				space_data.extracted_type = common::SpaceData::ExtractedType::eParticle; // eRawParticles
				space_data.extracted_particle_type = (common::SpaceData::ExtractedParticleType)std::stoi(argv[++i]);
			}
			else if (arg == "--save-mpi-rank") {
				from_cl.use_save_mpirank = true;
			}			
			// else if (arg == "--rawpart2vdb") {
			// 	from_cl.use_rawpart2vdb = true;
			// }			
			else if (arg == "--export-data") {
				space_data.particle_type = std::stoi(argv[++i]);
				space_data.block_name_id = std::stoi(argv[++i]);
				//from_cl.export_dense_type = std::stoi(argv[++i]);				
				from_cl.remote = false;
			}
			// === Density Computation Options ===
			else if (arg == "--dense-type") {
				// Set density computation type (kernel function)
				space_data.dense_type = (common::SpaceData::DenseType)std::stoi(argv[++i]);
				if (space_data.dense_type == common::SpaceData::DenseType::eNone)
					space_data.extracted_type = common::SpaceData::ExtractedType::eSparse;
				else
					space_data.extracted_type = common::SpaceData::ExtractedType::eDense;
			}
			else if (arg == "--dense-norm") {
				// Set density normalization method
				space_data.dense_norm = (common::SpaceData::DenseNorm)std::stoi(argv[++i]);
			}
			else if (arg == "--simple-density") {
				// Use simplified density calculation method
				space_data.use_simple_density = true;
				if (space_data.dense_type == common::SpaceData::DenseType::eNone)
					space_data.dense_type = common::SpaceData::DenseType::eCubic;

				space_data.extracted_type = common::SpaceData::ExtractedType::eDense;
			}
			
			// === Spatial Filtering Options ===
			else if (arg == "--bbox-sphere") {
				// Define spherical bounding region (center xyz + radius)
				space_data.use_bbox_sphere = true;
				space_data.bbox_sphere_pos[0] = std::stof(argv[++i]);
				space_data.bbox_sphere_pos[1] = std::stof(argv[++i]);
				space_data.bbox_sphere_pos[2] = std::stof(argv[++i]);
				space_data.bbox_sphere_r = std::stof(argv[++i]);
			}
			else if (arg == "--bbox") {
				// Define axis-aligned bounding box (min xyz, max xyz)
				space_data.bbox_min[0] = std::stof(argv[++i]);
				space_data.bbox_min[1] = std::stof(argv[++i]);
				space_data.bbox_min[2] = std::stof(argv[++i]);

				space_data.bbox_max[0] = std::stof(argv[++i]);
				space_data.bbox_max[1] = std::stof(argv[++i]);
				space_data.bbox_max[2] = std::stof(argv[++i]);
			}
			else if (arg == "--offset-position") {
				// Apply position offset to all particles (translation)
				space_data.offset_position[0] = std::stof(argv[++i]);
				space_data.offset_position[1] = std::stof(argv[++i]);
				space_data.offset_position[2] = std::stof(argv[++i]);
			}
			// === Neighbor Search Configuration ===
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
			else if (arg == "--calc-radius-neigh") {
				// Calculate smoothing radius based on N nearest neighbors
				space_data.calc_radius_neigh = std::stoi(argv[++i]);
			}
			else if (arg == "--calc-radius-neigh-rho-kernel") {
				// Specify kernel for density calculation with neighbor search
				space_data.calc_radius_neigh_rho_kernel = (common::SpaceData::DenseType)std::stoi(argv[++i]);
			}			
#endif			
			else if (arg == "--calc-radius-neigh-file") {
				// Load pre-computed neighbor radii from file
				space_data.calc_radius_neigh_file = argv[++i];
			}
			
			// === Neighbor Search Implementation ===
#ifdef WITH_CUDAKDTREE			
			else if (arg == "--cudakdtree") {
				// Use CUDA-accelerated KDTree for neighbor search
				from_cl.use_cudakdtree = true;
			}
			else if (arg == "--cudakdtree-cpu") {
				// Use CPU implementation of CUDA KDTree
				from_cl.use_cudakdtree_cpu = true;
			}
#endif
#ifdef WITH_NANOFLANN
			else if (arg == "--nanoflann") {
				// Use nanoflann library for neighbor search
				from_cl.use_nanoflann = true;
			}
#endif		
#ifdef WITH_GPU_CUDA
			else if (arg == "--gpu-cuda") {
				// Enable GPU acceleration for CUDA-based computations
				from_cl.use_gpu_cuda = true;
			}
#endif
			else if (arg == "--sort-by-radius") {
				from_cl.use_sort_by_radius = true;
			}
			else if (arg == "--sort-by-non-overlap") {
				from_cl.use_sort_by_non_overlap = true;
				}
			else if (arg == "-h" || arg == "--help") {
				// Display usage information
				usage();
			}
		}
	}

} // namespace space_converter

/*
 * ============================================================================
 * EXAMPLE USAGE COMMANDS
 * ============================================================================
 * 
 * Below are example command-line invocations demonstrating various use cases
 * and workflow patterns. These examples are organized by data format and
 * feature type for easy reference.
 * ============================================================================
 */

// === CHANGA Format Examples ===

// Basic TIPSY file processing with data export
//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\galaxy_merger.dat --export-data 0 1 0 --output-path f:\temp\

// NCHILADA format processing
//--data-type CHANGA_NCHILADA --grid-dim 1000 --output-path e:\temp\changa\tmp\ --nc-dir e:\temp\changa\galaxy_merger.dat.data

// TIPSY with info display only
//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --export-data 0 1 0 --output-path f:\temp\ --info

// TIPSY with remote server communication
//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --output-path f:\temp\ --port 5000

// TIPSY with multi-resolution processing
//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --output-path f:\temp\ --port 5000 --multires

// Accretion disk simulation
//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\accretiondisklowresstd --output-path f:\temp\ --port 5000

//--data-type CHANGA_TIPSY --tipsy-file e:\temp\21\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.00995 --output-path e:\temp\21\ --grid-dim 100 --export-data 0 1 --bbox 452.542480 506.425049 478.175903 510.356018 564.238525 535.989380 --gpu-cuda --dense-type 6 # --cudakdtree --calc-radius-neigh 10 

// === GADGET Format Examples ===

// Basic GADGET processing (small dataset)
//--data-type GADGET --gadget-file f:\temp\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\

// GADGET with remote server (larger dataset)
//--data-type GADGET --gadget-file f:\temp\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --port 5000

// GADGET with data export
//--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --export-data 0 1 0

// GADGET with black hole tracking
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110

// GADGET with info display
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --info

// GADGET galaxy simulation
//--data-type GADGET --gadget-file e:\temp\galaxy\snap_SQ_mitCGM_mbar5e5 --max-mem-size 10000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 100 --output-path f:\temp\ --export-data 0 1 0 --bh-count 4 --info

// GADGET with NanoVDB output
//--data-type GADGET --gadget-file F:\work\blender\SPHtoGrid.jl\snap_sedov --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --port 5005 --bh-count 1 --nanovdb

// GADGET with multi-resolution
//--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --port 5005 --bh-count 1 --multires


// === GADGET Animation Examples ===

// Black hole simulation animation (frames 0-1)
//--data-type GADGET --gadget-file e:\temp\black_hole\13_1e4_B1e9_z19\snapdir_{}\snap_{} --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 500 --output-path e:\temp\ --port 5000 --anim 0 1 --bh-count 1

// Black hole simulation animation (frames 1000-1001)
//--data-type GADGET --gadget-file e:\temp\black_hole\13_1e4_B1e9_z19\snapdir_{}\snap_{} --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 500 --output-path e:\temp\ --port 5000 --anim 1000 1001 --bh-count 1

// Animation with merged frames and raw particles
//--data-type GADGET --gadget-file e:\temp\black_hole\13_1e4_B1e9_z19\snapdir_{}\snap_{} --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 500 --output-path e:\temp\ --port 5000 --anim 1000 1001 --bh-count 1 --anim-merge --raw-particles


// === GADGET with Neighbor Search ===

// Neighbor radius calculation with KDTree
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh 2 --bbox 386.881104 590.247009 613.632996 393.591431 596.957336 620.343323 --export-data 0 1 1

// Using pre-computed neighbor radii from file
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh-file e:\temp\gadget_rad10.bin

// Advanced neighbor search with bounding box
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh 10 --bbox 386.881104 590.247009 613.632996 393.591431 596.957336 620.343323 --export-data 0 1 1

// Multi-resolution with neighbor calculation
//--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --export-data 0 1 1 --bh-count 1 --calc-radius-neigh 10 --multires

// Large-scale neighbor search
//--data-type GADGET --gadget-file e:/temp/galaxy/snap_SQ_mitCGM_mbar5e4 --max-mem-size 100000 --buffer-size 500.0 --part-alloc-factor 1.2 --grid-dim 1000 --output-path f:\temp\ --port 5000 --bh-count 4 --calc-radius-neigh 50


// === GADGET_SIMPLE Format Examples ===

// Basic GADGET_SIMPLE processing
//--data-type GADGET_SIMPLE --gadget-file e:\temp\alice\data\snapdir_1000\snap_1000 --output-path e:\temp\ --port 5000

// GADGET_SIMPLE with neighbor search
//--data-type GADGET_SIMPLE --gadget-file e:\temp\alice\data\snapdir_1000\snap_1000 --output-path e:\temp\ --port 5000 --calc-radius-neigh 32

// GADGET_SIMPLE animation sequence
//--data-type GADGET_SIMPLE --gadget-file e:\temp\alice\data\H1e5_DF_p03\snapdir_{}\snap_{} --output-path e:\temp\ --port 5000 --anim 0 3

// GADGET_SIMPLE with CUDA KDTree
//--data-type GADGET_SIMPLE --gadget-file e:\temp\gadget\very_small_example\snap_081 --output-path e:\temp\ --port 5000 --calc-radius-neigh 32 --cudakdtree

//--data-type GADGET_SIMPLE --gadget-file e:\temp\gadget\very_small_example\snap_081 --output-path e:\temp\20_new\ --export-data 0 1 --dense-type 6 --gpu-cuda  --sort-by-radius --nanovdb

// === HACC_GenericIO Format Examples ===

// GenericIO with info display
//--data-type HACC_GENERICIO --genericio-file e:\data\jar091\down\m000.mpicosmo.499 --grid-dim 1000 --output-path f:\temp\ --info

// GenericIO with remote server
//--data-type HACC_GENERICIO --genericio-file e:\temp\hacc\farpoint\m000p-499.select.mpicosmo --grid-dim 1000 --output-path f:\temp\ --port 5000


// === HACCBIN Format Examples ===

// HACC binary format with info
//--data-type HACC_BIN --haccbin-file e:\temp\hacc\Full.cosmo.0 --grid-dim 1000 --output-path f:\temp\ --info