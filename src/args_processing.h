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

#pragma once

#include <string>
#include "data_common.h"

/**
 * @namespace space_converter
 * @brief Provides functionality for converting simulation data formats and
 *        managing command-line argument parsing for configuration.
 * 
 * This namespace contains tools for parsing command-line arguments and
 * converting various scientific simulation data formats (GADGET, TIPSY, HDF5, etc.)
 * into volumetric representations suitable for visualization.
 */
namespace space_converter {

	/**
	 * @struct FromCL
	 * @brief Stores configuration options parsed from command-line arguments.
	 * 
	 * This structure holds all the parsed command-line options that control
	 * the behavior of the space converter tool, including I/O settings,
	 * processing options, and feature flags.
	 */
	struct FromCL {
		// === Core Configuration ===
		
		/// Type of input data format (e.g., "GADGET", "TIPSY")
		std::string data_type = "GADGET";
		
		/// Output directory path for generated files
		std::string output_path;
		
		/// Server address for remote operations (default: "localhost")
		std::string server = "localhost";
		
		/// Port number for server communication (default: 7000)
		int port = 7000;
		
		/// Flag to display dataset information only without processing
		bool info = false;

		/// Flag to enable remote processing mode
		bool remote = true;

		// === Output Format Options ===
		
#ifdef WITH_OPENVDB
		/// Flag to use NanoVDB format (GPU-friendly VDB variant)
		bool use_nanovdb = false;
#else
		/// Flag to use NanoVDB format (default when OpenVDB not available)
		bool use_nanovdb = true;
#endif

		/// Flag to save MPI rank information for parallel processing
		bool use_save_mpirank = false;
		
		/// Flag to convert raw particle data to VDB file format
		// bool use_rawpart2vdb = false;

		/// Flag to export dense matrix representation directly to file
		// bool use_dense2file = false;

		// === Neighbor Search Options ===
		
#ifdef WITH_CUDAKDTREE
		/// Flag to use CUDA-accelerated KDTree for neighbor search
		bool use_cudakdtree = false;
		
		/// Flag to use CPU-based implementation of CUDA KDTree
		bool use_cudakdtree_cpu = false;
#endif

#ifdef WITH_NANOFLANN
		/// Flag to use nanoflann library for neighbor search
		bool use_nanoflann = false;
#endif

#ifdef WITH_GPU_CUDA
		/// Flag to enable GPU acceleration for CUDA-based computations
		bool use_gpu_cuda = false;
#endif

		bool use_sort_by_radius = false;
		bool use_sort_by_non_overlap = false;

		// === MPI Configuration ===
		
		/// Rank of the current process in MPI (default: 0 for single process)
		int world_rank = 0;
		
		/// Total number of processes in MPI communicator (default: 1 for single process)
		int world_size = 1;

		/**
		 * @brief Print all current configuration values.
		 * 
		 * This function outputs all member variables of the FromCL structure
		 * to standard output, useful for debugging and verification of parsed
		 * command-line arguments.
		 */
		void print_info() const;
	};

	/**
	 * @brief Parse command-line arguments and populate configuration structures.
	 * 
	 * This function processes command-line arguments passed to the application and
	 * fills both the FromCL structure (for general configuration) and SpaceData
	 * structure (for spatial/simulation-specific parameters). It validates arguments
	 * and calls usage() if insufficient arguments are provided.
	 * 
	 * @param from_cl Reference to FromCL structure to populate with parsed configuration
	 * @param space_data Reference to SpaceData structure for spatial/simulation parameters
	 * @param argc Number of command-line arguments (including program name)
	 * @param argv Array of command-line argument strings
	 * 
	 * @note The function expects at least 3 arguments: program name, --data-type, and format.
	 *       If fewer arguments are provided, it displays usage information and exits.
	 * 
	 * @see FromCL
	 * @see common::SpaceData
	 * @see usage()
	 */
	void parse_args(space_converter::FromCL& from_cl, common::SpaceData& space_data, int argc, char** argv);

} // namespace space_converter
