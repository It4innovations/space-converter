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

#include "args_processing.h"
#include "convert_vdb.h"
#include "data_common.h"

/**
 * @namespace space_converter
 * @brief Provides the core data processing pipeline for converting simulation data to VDB format.
 * 
 * This namespace contains functions for the complete conversion workflow:
 * - Initialization and cleanup of format-specific converters
 * - Bounding box calculation and spatial partitioning
 * - Grid creation and particle-to-volume conversion
 * - MPI reduction and aggregation across distributed processes
 * - Final grid optimization and output (VDB, NanoVDB, or raw formats)
 */
namespace space_converter {

	// ========================================================================
	// Converter Lifecycle Management
	// ========================================================================

	/**
	 * @brief Initialize the VDB converter for a specific simulation format.
	 * 
	 * Creates a format-specific converter (GADGET, TIPSY, etc.) based on
	 * command-line arguments. Initializes I/O libraries, allocates particle buffers,
	 * and optionally computes particle radii using neighbor search algorithms.
	 * 
	 * @param argc Number of command-line arguments
	 * @param argv Array of command-line argument strings
	 * @param from_cl Reference to FromCL struct with parsed configuration
	 * @param space_data Reference to SpaceData for storing dataset metadata
	 * @return Pointer to initialized ConvertVDBBase-derived object
	 * 
	 * @note Caller is responsible for calling deinit_converter() to free resources.
	 */
	common::vdb::ConvertVDBBase* init_converter(int argc, char** argv, space_converter::FromCL& from_cl, common::SpaceData& space_data);

	/**
	 * @brief Clean up and destroy the VDB converter.
	 * 
	 * Finalizes I/O library resources and deallocates the converter object.
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object to destroy
	 */
	void deinit_converter(common::vdb::ConvertVDBBase* convert_vdb_base);

	// ========================================================================
	// Dataset Information and Analysis
	// ========================================================================

	/**
	 * @brief Query and print dataset information across all MPI ranks.
	 * 
	 * Gathers particle type counts and data block information from all processes
	 * using MPI_Allreduce. Prints summary to stdout on rank 0 if info mode enabled.
	 * 
	 * @param convert_vdb_base Pointer to initialized ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param types_and_blocks_global Output vector containing global particle type and block counts
	 */
	void print_info(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, std::vector<int>& types_and_blocks_global);

	/**
	 * @brief Calculate spatial bounding box for particle data.
	 * 
	 * Computes min/max extents across all particles using MPI reduction. Optionally
	 * makes the bounding box cubic (symmetric in X/Y/Z) for uniform voxel spacing.
	 * Result stored in space_data.bbox_min_orig and bbox_max_orig.
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param space_data Reference to SpaceData to store bounding box results
	 * @param particle_type Particle type to analyze (-1 for all types)
	 */
	void find_bbox(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, int particle_type = -1);

	// ========================================================================
	// Grid Creation and Conversion
	// ========================================================================

	/**
	 * @brief Create an empty VDB grid with specified resolution and transform.
	 * 
	 * Allocates grid structure based on extraction type (sparse VDB or dense volume).
	 * Sets up coordinate transformation from world space to index space.
	 * 
	 * @param grid_main Reference to VDBParticles object to initialize
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param space_data Reference to SpaceData with grid dimensions and transform
	 */
	void create_grid(common::vdb::ConvertVDBBase* convert_vdb_base, common::vdb::VDBParticles& grid_main, space_converter::FromCL& from_cl, common::SpaceData& space_data);

	/**
	 * @brief Convert particle data to volumetric grid using density estimation.
	 * 
	 * Reads particles from I/O library and splatts them into the VDB grid using
	 * various density kernels (simple, Gaussian, etc.). Handles spatial filtering,
	 * radius-based smoothing, and optional neighbor search for adaptive radii.
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param space_data Reference to SpaceData with conversion parameters
	 * @param grid_main Reference to VDBParticles object to populate
	 */
	void convert_to_grid(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main);

	// ========================================================================
	// Min/Max Value Analysis
	// ========================================================================

	/**
	 * @brief Find minimum and maximum density values across all MPI ranks.
	 * 
	 * Uses MPI_Allreduce to compute global min/max from local values computed
	 * during grid conversion. Also aggregates total particle count.
	 * 
	 * @param from_cl Reference to FromCL struct with MPI configuration
	 * @param space_data Reference to SpaceData to store min/max results
	 */
	void find_minmax_value(space_converter::FromCL& from_cl, common::SpaceData& space_data);

	/**
	 * @brief Find min/max values after MPI reduction for animation sequences.
	 * 
	 * For animation modes, computes min/max across all frames to ensure consistent
	 * value ranges throughout the sequence.
	 * 
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param space_data Reference to SpaceData to store reduced min/max
	 */
	void find_minmax_reduced_value(space_converter::FromCL& from_cl, common::SpaceData& space_data);

	// ========================================================================
	// MPI Reduction and Grid Aggregation
	// ========================================================================

	/**
	 * @brief Combine grids from all MPI processes using reduction operations.
	 * 
	 * For dense grids: Uses mpi_reduce() to sum values element-wise.
	 * For sparse grids: Merges active voxel topology from all ranks.
	 * Supports animation modes with per-frame or merged output.
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with MPI configuration
	 * @param space_data Reference to SpaceData with metadata
	 * @param grid_main Reference to local VDBParticles from this rank
	 * @param grid_main_sum Output VDBParticles containing reduced result (valid on rank 0)
	 */
	void reduction(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, common::vdb::VDBParticles& grid_main_sum);

	// ========================================================================
	// Grid Finalization and Optimization
	// ========================================================================

	/**
	 * @brief Apply final transformations and optimizations to the grid.
	 * 
	 * Operations may include:
	 * - Converting NanoVDB to OpenVDB format	 
	 * - Grid pruning and compression
	 * - Normalization and value remapping
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param space_data Reference to SpaceData with metadata
	 * @param grid_main_sum Reference to reduced VDBParticles from reduction phase
	 * @param grid_main_final Output VDBParticles containing finalized grid
	 */
	void finalize_grid(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main_sum, common::vdb::VDBParticles& grid_main_final);

	// ========================================================================
	// Output and Serialization
	// ========================================================================

	/**
	 * @brief Save VDB grid to file in various formats.
	 * 
	 * Supports OpenVDB (.vdb), NanoVDB, and raw particle formats. Handles both
	 * single-file and per-frame output for animation sequences.
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with output configuration
	 * @param space_data Reference to SpaceData with file path and metadata
	 * @param grid_main_final Reference to finalized VDBParticles to save
	 * @param particle_type Type of VDB data (FloatGrid, NanoVDB, RawParticles)
	 * @param only_rank0 If true, only rank 0 writes output (default: true)
	 */
	void save_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main_final, common::vdb::VDBParticleType particle_type, bool only_rank0 = true);

	/**
	 * @brief Save dense volume as raw binary file.
	 * 
	 * Exports grid as uncompressed 3D array for tools that cannot read VDB format.
	 * Includes companion metadata file with dimensions and value ranges.
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param space_data Reference to SpaceData with output path
	 * @param grid_main Reference to VDBParticles containing dense data
	 * @param only_rank0 If true, only rank 0 writes output (default: true)
	 */
	void save_raw_volume(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, bool only_rank0 = true);

	/**
	 * @brief Save raw particle positions to OpenVDB PointDataGrid.
	 * 
	 * Converts raw particles to VDB's native point representation for efficient
	 * storage and rendering of particle systems.
	 * 
	 * @param convert_vdb_base Pointer to ConvertVDBBase object
	 * @param from_cl Reference to FromCL struct with configuration
	 * @param space_data Reference to SpaceData with metadata
	 * @param grid_main Reference to VDBParticles containing particle data
	 */
	void save_raw_particles_to_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main);

	// ========================================================================
	// Testing and Debugging
	// ========================================================================

	/**
	 * @brief Test harness for converter functionality (debug builds only).
	 * 
	 * Runs validation tests on the converter implementation. Not used in production.
	 * 
	 * @param argc Number of command-line arguments
	 * @param argv Array of command-line argument strings
	 * @param from_cl Reference to FromCL struct with configuration
	 */
	void test_converter(int argc, char** argv, space_converter::FromCL& from_cl);

} // namespace space_converter
