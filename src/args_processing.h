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

#pragma once

#include <string>
#include "data_common.h"

 // Namespace space_converter provides functionality for converting simulation data
 // and managing command-line argument parsing for configuration.
namespace space_converter {

	// Struct FromCL stores configuration options parsed from the command line.
	struct FromCL {
		std::string data_type = "GADGET";   // Type of data (e.g., "GADGET", "TIPSY").
		//int grid_dim = 100;               // Dimension of the computational grid.
		std::string output_path;           // Path to save the output data.
		std::string server = "localhost";  // Server address for remote operations.
		int port = 7000;                   // Port for server communication.
		bool info = false;                 // Flag to display information.

		bool remote = true;                // Flag for remote processing.

		//int export_type = -1;              // Export type (-1 indicates not set).
		//int export_dataset = -1;           // Export dataset (-1 indicates not set).
		//int export_extracted_type = 0;     // Export dataset (0 indicates Sparse).
		//int export_dense_type = 0;         // Export dense_type (0 indicates not dense)
		//int export_dense_norm = 0;         // Export dense_norm (0 indicates no normalization)

#ifdef WITH_OPENVDB
		bool use_nanovdb = false;          // Flag to indicate use of NanoVDB.
#else
		bool use_nanovdb = true;           // Flag to indicate use of NanoVDB.
#endif
		//bool use_raw_particles = false;    // Flag to indicate use of particles (no VDB).
		bool use_save_mpirank = false;
		bool use_rawpart2vdb = false;      // Flag to export RAW particles to vdb file.

		//bool use_dense = false;            // Flag to indicate use of Dense format.
		bool use_dense2file = false;       // Flag to export RAW dense matrix to file.
		//bool use_anim = false;             // Flag to read timesteps.
		//int anim_type = 0;                 // Flag to read timesteps and merge to one.
		//int anim_start = -1;               // Start frame (-1 indicates not set).
		//int anim_end = -1;                 // End frame (-1 indicates not set).

#ifdef WITH_CUDAKDTREE
		bool use_cudakdtree = false;
		bool use_cudakdtree_cpu = false;
#endif

#ifdef WITH_NANOFLANN
		bool use_nanoflann = false;
#endif

		//#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
		//        int calc_radius_neigh = -1;
		//#endif
		//        std::string calc_radius_neigh_file = "";
		//
		//        bool use_bbox_sphere = false;
		//        float bbox_sphere_pos[3] = {0.0f,0.0f,0.0f};
		//        float bbox_sphere_r = 0.0f;
		//        bool use_simple_density = false;
		//        float offset_position[3] = {0.0f,0.0f,0.0f};

		bool use_multires = false;
		//bool use_bbox = false;
		//float bbox_pos[6] = { 0.0f,0.0f,0.0f, 0.0f,0.0f,0.0f };

		int world_rank = 0;                // Rank of the current process in MPI.
		int world_size = 1;                // Total number of processes in MPI.
	};

	// Function parse_args parses command-line arguments and populates the FromCL structure.
	// @param from_cl: Reference to the FromCL struct to populate with parsed values.
	// @param argc: Number of command-line arguments.
	// @param argv: Array of command-line argument strings.
	void parse_args(space_converter::FromCL& from_cl, common::SpaceData& space_data, int argc, char** argv);

} // namespace space_converter
