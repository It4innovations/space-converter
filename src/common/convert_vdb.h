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

#include <map>
#include <string>
#include <vector>

#include "convert_common.h"
#include "data_common.h"
#include "data_cache.h"

#include "dense_common.h"
#include "sparse_common.h"
#include "raw_common.h"

#ifdef WITH_NANOVDB

#ifdef WITH_TBB
#define NANOVDB_USE_TBB
#define NANOVDB_USE_INTRINSICS
#endif

#include <nanovdb/tools/GridBuilder.h>
#include <nanovdb/tools/CreateNanoGrid.h>
#include <nanovdb/io/IO.h>

#endif

#ifdef WITH_OPENVDB
#	include <openvdb/openvdb.h>
#	include <openvdb/points/PointDataGrid.h>
#endif

namespace space_converter {
	namespace common {
		namespace vdb {

			/// Normalize accumulated dense-grid values in place (voxel-volume or
			/// per-voxel weight normalization; shared by the dense -> NanoVDB /
			/// OpenVDB / CUB converters). Implemented in convert_vdb.cpp.
			void normalize_dense_values(VoxelDenseManager* dense_manager, double transform_scale, common::SpaceData::DenseNorm dense_norm);

			/**
			 * @brief Base class for converting particle data to VDB grids
			 *
			 * Provides functionality to convert various particle formats into VDB grid
			 * representations (dense, sparse, or raw). Handles spatial indexing, filtering,
			 * and normalization. Derived classes implement particle I/O library specifics.
			 */
			class ConvertVDBBase {
			public:
				cache::CacheManager cache_manager; ///< Cache manager for storing intermediate data

			public:
				virtual ~ConvertVDBBase() = default;

				/**
				 * @brief Convert particle data from I/O library format to VDB grid
				 *
				 * Main conversion function that reads particles and rasterizes them into the
				 * specified grid type (dense / sparse / raw, decided by grid.type). Handles
				 * filtering and spatial indexing; dispatches to the CPU or GPU kernels.
				 *
				 * @param particle_type Type of particles to convert
				 * @param particle_radius_multiplier Particle radius multiplier (0 = off)
				 * @param bbox_min/bbox_max Requested zoom box in object space [0, object_size]
				 * @param bbox_dim Grid dimension (resolution)
				 * @param bbox_min_orig/bbox_size_orig Cube-symmetrized dataset bbox (world space)
				 * @param grid Output VDB grid container
				 */
				void convert_iolib_to_grid(
					int particle_type,
					float particle_radius_multiplier,
					float* bbox_min,
					float* bbox_max,
					int bbox_dim,
					int* bbox_min_orig,
					double bbox_size_orig,
					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					int block_name_id,
					float object_size,
					float& min_value,
					float& max_value,
					size_t& particles_count,
					VDBParticles& grid,
					double& transform_scale,
					float filter_min,
					float filter_max,
					common::SpaceData::AnimType anim_type,
					int frame_req,
					int frame,
					float* bbox_sphere_pos,
					float bbox_sphere_r,
					bool use_simple_density,
					bool use_norm_value,
					float* offset_position
				);

				/**
				 * @brief Merge one VDB grid into another
				 *
				 * Combines grid data from grid_recv into grid_dst, handling different
				 * grid types (dense, sparse, serialized).
				 */
				 /**
				  * @brief Merge one VDB grid into another
				  * @param grid_dst Destination grid to merge into
				  * @param grid_recv Source grid to merge from
				  *
				  * Combines grid data from grid_recv into grid_dst, handling different
				  * grid types (dense, sparse, serialized).
				  */
				void merge_grid(
					VDBParticles& grid_dst,
					VDBParticles& grid_recv
				);


#ifdef WITH_NANOVDB

				/**
				 * @brief Convert dense grid to NanoVDB format
				 */
				std::shared_ptr<nanovdb::tools::build::FloatGrid> dense_to_nanovdb(VoxelDenseManager* dense_manager, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);
				std::shared_ptr<nanovdb::tools::build::FloatGrid> sparse_to_nanovdb(VoxelSparseManager* voxel_manager);

#endif


#ifdef WITH_OPENVDB
				/**
				 * @brief Convert dense grid to OpenVDB format
				 */
				 /**
				  * @brief Convert dense grid to OpenVDB format
				  * @param particles Dense particle grid
				  * @param transform_scale Grid transformation scale
				  * @param dense_type Type of density calculation
				  * @param dense_norm Normalization type
				  * @return OpenVDB grid pointer
				  */
				openvdb::FloatGrid::Ptr dense_to_openvdb(VoxelDenseManager* dense_manager, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);
				openvdb::FloatGrid::Ptr sparse_to_openvdb(VoxelSparseManager* voxel_manager);

				/**
				 * @brief Serialize OpenVDB grid to binary buffer
				 * @param grid OpenVDB grid to serialize
				 * @param file_content Output buffer for serialized data
				 */
				void openvdb_to_vector(openvdb::FloatGrid::Ptr grid, std::vector<uint8_t>& file_content);

				/**
				 * @brief Serialize two OpenVDB grids to binary buffer
				 */
				void openvdb_to_vector2(openvdb::FloatGrid::Ptr grid1, openvdb::FloatGrid::Ptr grid2, std::vector<uint8_t>& file_content);

				/**
				 * @brief Deserialize OpenVDB grid from binary buffer
				 */
				openvdb::FloatGrid::Ptr vector_to_openvdb(std::vector<uint8_t>& file_content);
#endif

				/**
				 * @brief Read particle radius data from file
				 * @param calc_radius_neigh_file Path to radius data file
				 */
				void read_radius_from_file(std::string& calc_radius_neigh_file);

				/**
				 * @brief Write particle radius data to file
				 * @param calc_radius_neigh_file Path to output file
				 */
				void write_radius_from_file(std::string& calc_radius_neigh_file);

				/**
				 * @brief Find particle positions
				 */
				void find_particle_positions();

				/**
				 * @brief Find particle values for a given particle type
				 * @param ptype Particle type
				 * @param block_name_id Data block ID
				 * @return True if values were found, false otherwise
				 */
				bool find_particle_values(int ptype, int block_name_id);

#ifdef WITH_CUDAKDTREE
				/**
				 * @brief Build per-type float4 KD-trees and store in cache_manager.
				 *
				 * Called once during init_converter() after find_particle_positions().
				 * Each tree node stores the offset-adjusted world position (x,y,z) and the
				 * original particle index bit-cast into w.  The tree is later used by
				 * convert_to_dense_grid_loop_over_voxels_cpu/gpu instead of building it
				 * on every conversion call.
				 *
				 * @param offset_position  Per-axis position offset (from --offset-position).
				 */
				void build_voxel_kdtree(float* offset_position);

				/**
				 * @brief Calculate particle radii using CUDA KD-Tree neighbor search
				 * @param calc_radius_neigh Number of neighbors to consider
				 * @param calc_radius_neigh_file Output file for radius data
				 * @param use_cycling Enable cycling boundary conditions
				 * @param use_cudakdtree_cpu Use CPU fallback instead of GPU
				 * @param maxRadius Maximum search radius
				 * @param rho_kernel Density kernel type for SPH calculations
				 */
				void calculate_radius_by_cudakdtree(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, bool use_cudakdtree_cpu, float maxRadius, common::SpaceData::DenseType& rho_kernel);
#endif

#ifdef WITH_NANOFLANN
				/**
				 * @brief Calculate particle radii using Nanoflann KD-Tree neighbor search
				 * @param calc_radius_neigh Number of neighbors to consider
				 * @param calc_radius_neigh_file Output file for radius data
				 * @param use_cycling Enable cycling boundary conditions
				 * @param rho_kernel Density kernel type for SPH calculations
				 */
				void calculate_radius_by_nanoflann(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, common::SpaceData::DenseType& rho_kernel);
#endif

				/**
				 * @brief Find bounding box for particles of given type
				 * @param particle_type Type of particles to consider
				 * @param bbox_min Output minimum coordinates
				 * @param bbox_max Output maximum coordinates
				 * @param offset_position Optional position offset
				 * @param use_gpu Flag to use GPU (true) or CPU (false) implementation
				 */
				virtual void iolib_find_bbox(
					int particle_type,
					float* bbox_min,
					float* bbox_max,
					float* offset_position
				);

				/**
				 * @brief Find minimum and maximum values for data block
				 * @param particle_type Type of particles
				 * @param block_nr Data block number
				 * @param v_min Output minimum value
				 * @param v_max Output maximum value
				 */
				 //virtual void iolib_find_minmax(
				 //	int particle_type,
				 //	int block_nr,
				 //	float& v_min,
				 //	float& v_max
				 //);

				 /**
				  * @brief Get density value for particle
				  * @param id Particle ID
				  * @return Density value
				  */
				virtual double get_particle_rho(uint64_t id);

				/**
				 * @brief Get available particle types and data blocks
				 * @param types_and_blocks Output vector of type/block pairs
				 */
				virtual void get_types_and_blocks(std::vector<int>& types_and_blocks);

				/**
				 * @brief Get normalized value for particle in data block
				 * @param blocknr Data block number
				 * @param id Particle ID
				 * @return Normalized value
				 */
				virtual float get_particle_norm_value(int blocknr, uint64_t id);

				/**
				 * @brief Get value for particle in data block
				 * @param blocknr Data block number
				 * @param id Particle ID
				 * @param out_value Output buffer for value(s)
				 * @return Number of components
				 */
				virtual int get_particle_value(int blocknr, uint64_t id, float* out_value);

				/**
				 * @brief Get number of components for particle value
				 * @param blocknr Data block number
				 * @param id Particle ID
				 * @return Number of components
				 */
				virtual int get_particle_value_comp(int blocknr, uint64_t id);

				/**
				 * @brief Format filename with number substitution
				 * @param pattern Filename pattern with format specifiers
				 * @param number Number to substitute into pattern
				 * @param zero_pad Number of digits for zero-padding
				 * @return Formatted filename
				 */
				std::string format_filename(const std::string& pattern, int number, int zero_pad = 0);

			public:
				// ============================================================
				// Pure Virtual Methods - Must be implemented by derived classes
				// ============================================================

				/** @brief Print CPU processing steps for debugging */
				virtual void print_CPU_steps() = 0;

				/** @brief Get normalized particle value (internal implementation) */
				virtual float get_particle_norm_value_internal(int blocknr, uint64_t id) = 0;

				/** @brief Get particle value (internal implementation) */
				virtual int get_particle_value_internal(int blocknr, uint64_t id, float* out_value) = 0;

				/** @brief Get particle value component count (internal implementation) */
				virtual int get_particle_value_comp_internal(int blocknr, uint64_t id) = 0;

				/** @brief Get particle type ID */
				virtual int get_particle_type(uint64_t id) = 0;

				/** @brief Get particle position coordinates */
				virtual void get_particle_position(uint64_t id, double* pos) const = 0;

				/** @brief Get number of particles on this MPI rank */
				virtual size_t get_local_num_particles() const = 0;

				/** @brief Get total number of particles across all ranks */
				virtual size_t get_global_num_particles() const = 0;

				/** @brief Get particle smoothing length (SPH) */
				virtual double get_particle_hsml(uint64_t id) = 0;

				/** @brief Get particle mass */
				virtual double get_particle_mass(uint64_t id) = 0;

				/** @brief Get particle density (internal implementation) */
				virtual double get_particle_rho_internal(uint64_t id) = 0;

				/** @brief Get data block number containing density values */
				virtual int get_particle_rho_blocknr() = 0;

				/** @brief Initialize I/O library with MPI parameters */
				virtual void init_lib(int argc, char** argv, int world_rank, int world_size) = 0;

				/** @brief Finalize and cleanup I/O library */
				virtual void finish_lib() = 0;

				/** @brief Get particle types and data blocks (internal implementation) */
				virtual void get_types_and_blocks_internal(std::vector<int>& types_and_blocks) = 0;

				/** @brief Print local particle types and blocks */
				virtual void print_types_and_blocks_local() = 0;

				/** @brief Print particle types and blocks across all ranks */
				virtual void print_types_and_blocks(std::vector<int>& types_and_blocks) = 0;

				/** @brief Get human-readable name for particle type */
				virtual std::string get_type_name(int type) = 0;

				/** @brief Get human-readable name for data block */
				virtual std::string get_dataset_name(int blocknr) = 0;

				/** @brief Get formatted string of all particle data type names */
				virtual std::string get_particle_data_type_names(std::vector<int>& types_and_blocks) = 0;

				/** @brief Get number of particle types */
				virtual int get_num_types() = 0;

				/** @brief Get number of data blocks */
				virtual int get_num_blocks() = 0;
			};
		}// vdb
	} //common
} // space_converter