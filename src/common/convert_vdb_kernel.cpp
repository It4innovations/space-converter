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

#include "convert_vdb_kernel.h"
#include "dense_common.h"
#include "raw_common.h"
#include "sparse_common.h"
#include "convert_vdb.h"
#include <float.h>
#include <unordered_map>
#include <stdexcept>

#ifdef WITH_OPENMP
# include <omp.h>
#endif

#include "data_common.h"

namespace common {
	namespace vdb {
		namespace kernel {

			// Using declarations for cleaner code
			using common::vdb::sparse::packCoord3;

			/**
			 * @brief CPU implementation of bounding box calculation using OpenMP
			 *
			 * Parallel reduction to compute min/max coordinates across all particles.
			 */
			void find_bbox_cpu(
				const float* pos_particles,
				size_t num_particles,
				const float* offset_position,
				float* bbox_min,
				float* bbox_max
			) {
				if (num_particles == 0) {
					bbox_min[0] = bbox_min[1] = bbox_min[2] = 0.0f;
					bbox_max[0] = bbox_max[1] = bbox_max[2] = 0.0f;
					return;
				}

				float min_x = bbox_min[0], min_y = bbox_min[1], min_z = bbox_min[2];
				float max_x = bbox_max[0], max_y = bbox_max[1], max_z = bbox_max[2];

				// Parallel reduction to find min/max coordinates
#ifdef WITH_OPENMP
#pragma omp parallel for reduction(min : min_x, min_y, min_z) reduction(max : max_x, max_y, max_z)
#endif
				for (size_t i = 0; i < num_particles; ++i) {
					// Load position (x, y, z are interleaved: [x0, y0, z0, x1, y1, z1, ...])
					float px = pos_particles[i * 3 + 0] - offset_position[0];
					float py = pos_particles[i * 3 + 1] - offset_position[1];
					float pz = pos_particles[i * 3 + 2] - offset_position[2];

					if (px < min_x) min_x = px;
					if (py < min_y) min_y = py;
					if (pz < min_z) min_z = pz;

					if (px > max_x) max_x = px;
					if (py > max_y) max_y = py;
					if (pz > max_z) max_z = pz;
				}

				bbox_min[0] = min_x;
				bbox_min[1] = min_y;
				bbox_min[2] = min_z;

				bbox_max[0] = max_x;
				bbox_max[1] = max_y;
				bbox_max[2] = max_z;
			}


			/**
			 * @brief Process particles and convert to sparse grid format (NanoVDB/OpenVDB)
			 *
			 * Particles are rasterized into voxels using thread-local hash maps to avoid race conditions.
			 * Thread-local maps are later merged into a global voxel manager.
			 *
			 * @param pos_particles Particle positions (interleaved xyz)
			 * @param radius_particles Per-particle radii (can be nullptr)
			 * @param num_particles Total number of particles to process
			 * @param particle_radius_multiplier Multiplier for particle radius
			 * @param bbox_min_orig Bounding box minimum corner
			 * @param bbox_size_orig Bounding box size
			 * @param bbox_dim Grid resolution
			 * @param offset_position Position offset to apply
			 * @param filter_min Minimum value filter threshold
			 * @param filter_max Maximum value filter threshold
			 * @param bbox_sphere_pos Spherical filter center
			 * @param bbox_sphere_r Spherical filter radius
			 * @param anim_type Animation extraction type
			 * @param frame_req Requested frame number
			 * @param frame Current frame number
			 * @param use_simple_density Use uniform density (1.0) for all particles
			 * @param voxel_manager Output voxel manager for storing sparse data
			 * @param min_value [out] Minimum particle value encountered
			 * @param max_value [out] Maximum particle value encountered
			 * @param particles_count [out] Number of particles processed
			 */
			void convert_to_sparse_grid_cpu(
				const float* pos_particles,
				const size_t* particle_ids,
				const float* radius_particles,
				const float* value_particles,
				size_t num_particles,
				float particle_radius_multiplier,
				int* bbox_min_orig,
				double bbox_size_orig,
				int bbox_dim,
				float* offset_position,
				float filter_min,
				float filter_max,
				float* bbox_sphere_pos,
				float bbox_sphere_r,
				common::SpaceData::AnimType anim_type,
				int frame_req,
				int frame,
				bool use_simple_density,
				double particle_radius_const,

				double bbox_x_min_norm,
				double bbox_x_max_norm,
				double bbox_y_min_norm,
				double bbox_y_max_norm,
				double bbox_z_min_norm,
				double bbox_z_max_norm,
				double scale_space_diagonal,

				VoxelSparseManager* voxel_manager,
				float& min_value,
				float& max_value,
				size_t& particles_count
			) {
				common::vdb::sparse::VoxelOpenMPManager* omp_manager = dynamic_cast<common::vdb::sparse::VoxelOpenMPManager*>(voxel_manager);
				if (omp_manager) {
					// Initialize output tracking variables
					float min = FLT_MAX;
					float max = -FLT_MAX;
					size_t particles_count_temp = 0;

					// Create thread-local hash maps to avoid race conditions during parallel processing
#ifdef WITH_OPENMP
					const int thread_count = omp_get_max_threads();
#else
					const int thread_count = 1;
#endif
					std::vector<std::unordered_map<uint64_t, float>> sparse_thread_maps(thread_count);

#ifdef WITH_OPENMP
#pragma omp parallel 
					{
						// Get thread-local hash map for accumulating voxel values
						int tid = omp_get_thread_num();
						auto& sparse_local_map = sparse_thread_maps[tid];
						sparse_local_map.reserve((size_t)num_particles / (size_t)thread_count + 64);

#pragma omp for reduction(min : min) reduction(max : max) reduction(+ : particles_count_temp) schedule(dynamic, 64) //TODO
#else
						{
							// For non-OpenMP builds, use single thread map
							auto& sparse_local_map = sparse_thread_maps[0];
							sparse_local_map.reserve(num_particles);
#endif
							for (size_t ic = 0; ic < num_particles; ic++) {
								// Cached particle index
								size_t cached_idx = particle_ids[ic];

								// Particle processing outputs
								double Pos[3];
								double px_norm, py_norm, pz_norm;
								int px, py, pz;
								float v_orig = 0.0f;

								// Check if particle passes filters and compute voxel coordinates
								bool should_process = process_particle(
									cached_idx,
									pos_particles,
									radius_particles,
									value_particles,
									offset_position,
									bbox_min_orig,
									bbox_size_orig,
									bbox_dim,
									scale_space_diagonal,
									particle_radius_multiplier,
									bbox_x_min_norm, bbox_x_max_norm,
									bbox_y_min_norm, bbox_y_max_norm,
									bbox_z_min_norm, bbox_z_max_norm,
									filter_min, filter_max,
									bbox_sphere_pos, bbox_sphere_r,
									anim_type, frame_req, frame,
									use_simple_density,
									particle_radius_const,

									false, // is_dense_grid = false for sparse
									// Outputs
									Pos, px_norm, py_norm, pz_norm,
									px, py, pz, v_orig
								);

								if (!should_process) {
									continue;
								}

								// Count processed particles
								particles_count_temp++;

								// Track value range for normalization
								if (v_orig < min) min = v_orig;
								if (v_orig > max) max = v_orig;

								// Pack voxel coordinates into 64-bit key and accumulate value
								const uint64_t key = packCoord3(px, py, pz);
								sparse_local_map[key] += v_orig;
							}

					}

					// Store output statistics
					particles_count = particles_count_temp;
					min_value = min;
					max_value = max;

					// Merge thread-local maps into global voxel manager
					// First, estimate total unique voxels across all threads
					size_t total_local_unique = 0;
					for (const auto& m : sparse_thread_maps) {
						total_local_unique += m.size();
					}

					// Merge all thread-local maps (accumulating values for duplicate voxels)
					std::unordered_map<uint64_t, float> merged;
					merged.reserve(total_local_unique + total_local_unique / 4 + 1);

					for (const auto& m : sparse_thread_maps) {
						for (const auto& kv : m) {
							merged[kv.first] += kv.second;
						}
					}

					// Transfer merged voxels to voxel manager
					for (const auto& kv : merged) {
						voxel_manager->insertOrUpdatePackedSequential(kv.first, kv.second);
					}
				} else {
					// OpenVDB, NanoVDB
					// Initialize output tracking variables
					float min = FLT_MAX;
					float max = -FLT_MAX;
					size_t particles_count_temp = 0;

					for (size_t ic = 0; ic < num_particles; ic++) {
						// Cached particle index
						size_t cached_idx = particle_ids[ic];

						// Particle processing outputs
						double Pos[3];
						double px_norm, py_norm, pz_norm;
						int px, py, pz;
						float v_orig = 0.0f;

						// Check if particle passes filters and compute voxel coordinates
						bool should_process = process_particle(
							cached_idx,
							pos_particles,
							radius_particles,
							value_particles,
							offset_position,
							bbox_min_orig,
							bbox_size_orig,
							bbox_dim,
							scale_space_diagonal,
							particle_radius_multiplier,
							bbox_x_min_norm, bbox_x_max_norm,
							bbox_y_min_norm, bbox_y_max_norm,
							bbox_z_min_norm, bbox_z_max_norm,
							filter_min, filter_max,
							bbox_sphere_pos, bbox_sphere_r,
							anim_type, frame_req, frame,
							use_simple_density,
							particle_radius_const,

							false, // is_dense_grid = false for sparse
							// Outputs
							Pos, px_norm, py_norm, pz_norm,
							px, py, pz, v_orig
						);

						if (!should_process) {
							continue;
						}

						// Count processed particles
						particles_count_temp++;

						// Track value range for normalization
						if (v_orig < min) min = v_orig;
						if (v_orig > max) max = v_orig;

						common::vdb::sparse::Voxel voxel(px, py, pz, v_orig);
						voxel_manager->insertOrUpdate(&voxel, 1);
					}
				}
			}

		/**
		 * @brief Process particles and convert to dense grid format
		 *
		 * Particles are splatted into a dense 3D grid using SPH-like kernel smoothing.
		 * Each particle affects multiple neighboring voxels based on its radius.
		 *
		 * @param pos_particles Particle positions (interleaved xyz)
		 * @param radius_particles Per-particle radii (can be nullptr)
		 * @param num_particles Total number of particles to process
		 * @param particle_radius_multiplier Multiplier for particle radius
		 * @param bbox_min_orig Bounding box minimum corner
		 * @param bbox_size_orig Bounding box size
		 * @param bbox_dim Grid resolution
		 * @param scale_space_diagonal Space diagonal scaling factor
		 * @param dense_type SPH kernel type to use for splatting
		 * @param dense_norm Normalization method for density
		 * @param block_name_id Data block identifier
		 * @param offset_position Position offset to apply
		 * @param filter_min Minimum value filter threshold
		 * @param filter_max Maximum value filter threshold
		 * @param bbox_sphere_pos Spherical filter center
		 * @param bbox_sphere_r Spherical filter radius
		 * @param anim_type Animation extraction type
		 * @param frame_req Requested frame number
		 * @param frame Current frame number
		 * @param use_simple_density Use uniform density (1.0) for all particles
		 * @param grid Dense grid data structure to fill
		 * @param min_value [out] Minimum particle value encountered
		 * @param max_value [out] Maximum particle value encountered
		 * @param particles_count [out] Number of particles processed
		 */
		void convert_to_dense_grid_cpu(
			const float* pos_particles,
			const size_t* particle_ids,
			const float* radius_particles,
			const float* value_particles,
			size_t num_particles,
			float particle_radius_multiplier,
			int* bbox_min_orig,
			double bbox_size_orig,
			int bbox_dim,

			common::SpaceData::DenseType dense_type,
			common::SpaceData::DenseNorm dense_norm,
			int block_name_id,
			float* offset_position,
			float filter_min,
			float filter_max,
			float* bbox_sphere_pos,
			float bbox_sphere_r,
			common::SpaceData::AnimType anim_type,
			int frame_req,
			int frame,
			bool use_simple_density,
			double particle_radius_const,

			double bbox_x_min_norm,
			double bbox_x_max_norm,
			double bbox_y_min_norm,
			double bbox_y_max_norm,
			double bbox_z_min_norm,
			double bbox_z_max_norm,
			double scale_space_diagonal,

			VoxelDenseManager* grid,
			float& min_value,
			float& max_value,
			size_t& particles_count
		) {
			// Initialize output tracking variables
			float min = FLT_MAX;
			float max = -FLT_MAX;
			size_t particles_count_temp = 0;

#ifdef WITH_OPENMP
#pragma omp parallel for reduction(min : min) reduction(max : max) reduction(+ : particles_count_temp) schedule(dynamic, 64) //TODO
#endif
			for (size_t ic = 0; ic < num_particles; ic++) {
				// Cached particle index
				size_t cached_idx = particle_ids[ic];

				// Particle processing outputs
				double Pos[3];
				double px_norm, py_norm, pz_norm;
				int px, py, pz;
				float v_orig = 0.0f;

				// Check if particle passes filters and compute voxel coordinates
				bool should_process = process_particle(
					cached_idx,
					pos_particles,
					radius_particles,
					value_particles,
					offset_position,
					bbox_min_orig,
					bbox_size_orig,
					bbox_dim,
					scale_space_diagonal,
					particle_radius_multiplier,
					bbox_x_min_norm, bbox_x_max_norm,
					bbox_y_min_norm, bbox_y_max_norm,
					bbox_z_min_norm, bbox_z_max_norm,
					filter_min, filter_max,
					bbox_sphere_pos, bbox_sphere_r,
					anim_type, frame_req, frame,
					use_simple_density,
					particle_radius_const, 

					true, // is_dense_grid = true for dense
					// Outputs
					Pos, px_norm, py_norm, pz_norm,
					px, py, pz, v_orig
				);

				if (!should_process) {
					continue;
				}

				// Count processed particles
				particles_count_temp++;

				// Track value range for normalization
				if (v_orig < min) min = v_orig;
				if (v_orig > max) max = v_orig;

				// Splat particle into dense grid using SPH kernel function
				// This will affect multiple neighboring voxels based on particle radius
				fill_voxels(
					grid->data_density.data(),
					grid->data_temp.size() > 0 ? grid->data_temp.data() : nullptr,
					grid->offset,
					grid->dims,
					cached_idx,
					v_orig,
					bbox_dim,
					bbox_min_orig,
					bbox_size_orig,
					scale_space_diagonal,
					dense_type,
					dense_norm,
					particle_radius_multiplier,
					radius_particles,
					particle_radius_const,
					block_name_id,
					Pos
				);
			}

			// Store output statistics
			particles_count = particles_count_temp;
			min_value = min;
			max_value = max;
		}

		/**
		 * @brief Process particles and store as raw particle data (no rasterization)
		 *
		 * Particles are stored directly with their positions, values, radii, and frame numbers.
		 * No grid conversion is performed - this preserves exact particle data.
		 *
		 * @param pos_particles Particle positions (interleaved xyz)
		 * @param radius_particles Per-particle radii (can be nullptr)
		 * @param num_particles Total number of particles to process
		 * @param particle_radius_multiplier Multiplier for particle radius
		 * @param bbox_min_orig Bounding box minimum corner
		 * @param bbox_size_orig Bounding box size
		 * @param bbox_dim Grid resolution
		 * @param object_size Scaling factor for output positions
		 * @param block_name_id Data block identifier for value extraction
		 * @param offset_position Position offset to apply
		 * @param filter_min Minimum value filter threshold
		 * @param filter_max Maximum value filter threshold
		 * @param bbox_sphere_pos Spherical filter center
		 * @param bbox_sphere_r Spherical filter radius
		 * @param anim_type Animation extraction type
		 * @param frame_req Requested frame number
		 * @param frame Current frame number
		 * @param use_simple_density Use uniform density (1.0) for all particles
		 * @param grid Raw particle data structure to fill
		 * @param min_value [out] Minimum particle value encountered
		 * @param max_value [out] Maximum particle value encountered
		 * @param particles_count [out] Number of particles processed
		 */
		void convert_to_raw_particles_cpu(
			const float* pos_particles,
			const size_t* particle_ids,
			const float* radius_particles,
			const float* value_particles,
			size_t num_particles,
			float particle_radius_multiplier,
			int* bbox_min_orig,
			double bbox_size_orig,
			int bbox_dim,
			float object_size,
			int block_name_id,
			float* offset_position,
			float filter_min,
			float filter_max,
			float* bbox_sphere_pos,
			float bbox_sphere_r,
			common::SpaceData::AnimType anim_type,
			int frame_req,
			int frame,
			bool use_simple_density,
			double particle_radius_const,

			double bbox_x_min_norm,
			double bbox_x_max_norm,
			double bbox_y_min_norm,
			double bbox_y_max_norm,
			double bbox_z_min_norm,
			double bbox_z_max_norm,
			double scale_space_diagonal,

			RawParticles& raw_particles,
			float& min_value,
			float& max_value,
			size_t& particles_count
		) {
			// Initialize output tracking variables
			float min = FLT_MAX;
			float max = -FLT_MAX;
			size_t particles_count_temp = 0;

			// Temporary arrays to store particle data (will be appended to output at end)
			RawParticles::ParticleData raw_positions;
			raw_positions.name = "position";
			raw_positions.num_comp = 3; // xyz

			RawParticles::ParticleData raw_values;
			raw_values.name = "value";

			RawParticles::ParticleData raw_radius;
			raw_radius.name = "radius";
			raw_radius.num_comp = 1;

			RawParticles::ParticleData raw_frame;
			raw_frame.name = "frame";
			raw_frame.num_comp = 1;

			// Note: Raw particle storage cannot be easily parallelized due to push_back operations
			// Consider using thread-local vectors and merging if performance is critical
			for (size_t ic = 0; ic < num_particles; ic++) {
				// Cached particle index
				size_t cached_idx = particle_ids[ic];

				// Particle processing outputs
				double Pos[3];
				double px_norm, py_norm, pz_norm;
				int px, py, pz;
				float v_orig = 0.0f;

				// Check if particle passes filters
				bool should_process = process_particle(
					cached_idx,
					pos_particles,
					radius_particles,
					value_particles,
					offset_position,
					bbox_min_orig,
					bbox_size_orig,
					bbox_dim,
					scale_space_diagonal,
					particle_radius_multiplier,
					bbox_x_min_norm, bbox_x_max_norm,
					bbox_y_min_norm, bbox_y_max_norm,
					bbox_z_min_norm, bbox_z_max_norm,
					filter_min, filter_max,
					bbox_sphere_pos, bbox_sphere_r,
					anim_type, frame_req, frame,
					use_simple_density,
					particle_radius_const,

					false, // is_dense_grid = false for raw
					// Outputs
					Pos, px_norm, py_norm, pz_norm,
					px, py, pz, v_orig
				);

				if (!should_process) {
					continue;
				}

				// Count processed particles
				particles_count_temp++;

				// Track value range for normalization
				if (v_orig < min) min = v_orig;
				if (v_orig > max) max = v_orig;

				// Store particle position in world space coordinates
				raw_positions.values.push_back((double)px_norm * (double)object_size);
				raw_positions.values.push_back((double)py_norm * (double)object_size);
				raw_positions.values.push_back((double)pz_norm * (double)object_size);

				// Get particle attribute values (may be multi-component, e.g., velocity vector)
				// TODO: Implement get_particle_value_comp and get_particle_value
				// int n_comp = get_particle_value_comp(block_name_id, cached_idx);
				// std::vector<float> values(n_comp);
				// get_particle_value(block_name_id, cached_idx, values.data());
				// raw_values.num_comp = n_comp;
				// raw_values.values.insert(raw_values.values.end(), values.begin(), values.end());

				// Placeholder: store simple value
				raw_values.num_comp = 1;
				raw_values.values.push_back(v_orig);

				// Compute and store particle radius
				//double pr = get_particle_radius(
				//	cached_idx,
				//	bbox_dim,
				//	bbox_min_orig,
				//	bbox_size_orig,
				//	scale_space_diagonal,
				//	particle_radius_multiplier,
				//	radius_particles,
				//	particle_radius_const
				//);
				double pr = get_particle_radius(
					cached_idx,
					1,
					bbox_min_orig,
					bbox_size_orig,
					1.0/object_size,
					particle_radius_multiplier,
					radius_particles,
					particle_radius_const
				);
				raw_radius.values.push_back(pr);

				// Store animation frame number for temporal data
				raw_frame.values.push_back(frame);
			}

			// Store output statistics
			particles_count = particles_count_temp;
			min_value = min;
			max_value = max;

			// Append particle data arrays to output structure
			raw_particles.data.push_back(raw_positions);
			raw_particles.data.push_back(raw_values);
			raw_particles.data.push_back(raw_radius);
			raw_particles.data.push_back(raw_frame);
		}

		/**
		 * @brief Main dispatcher function to convert particle data to grid format (CPU implementation)
		 *
		 * Routes to appropriate grid conversion function based on grid type:
		 * - Sparse grids (NanoVDB/OpenVDB): Hash-based voxel accumulation
		 * - Dense grids: SPH kernel splatting into 3D array
		 * - Raw particles: Direct storage without rasterization
		 *
		 * This function serves as a unified interface for all grid conversion types.
		 */
		void convert_iolib_to_grid_cpu(
			const float* pos_particles,
			const size_t* particle_ids,
			const float* radius_particles,
			const float* value_particles,
			size_t num_particles,
			float particle_radius_multiplier,
			std::string grid_name,
			float grid_transform,
			float* bbox_min,
			float* bbox_max,
			int bbox_dim,
			int* bbox_min_orig,
			double bbox_size_orig,
			common::SpaceData::ExtractedType extracted_type,
			common::SpaceData::ExtractedParticleType extracted_particle_type,
			common::SpaceData::DenseType dense_type,
			common::SpaceData::DenseNorm dense_norm,
			int block_name_id,
			float object_size,
			float& min_value,
			float& max_value,
			float min_value_global,
			float max_value_global,
			size_t& particles_count,
			VDBParticles& grid,
			double& transform_scale,
			float filter_min,
			float filter_max,
			float min_rho,
			float max_rho,
			common::SpaceData::AnimType anim_type,
			int frame_req,
			int frame,
			float* bbox_sphere_pos,
			float bbox_sphere_r,
			bool use_simple_density,
			double particle_radius_const,

			double bbox_x_min_norm,
			double bbox_x_max_norm,
			double bbox_y_min_norm,
			double bbox_y_max_norm,
			double bbox_z_min_norm,
			double bbox_z_max_norm,
			double scale_space_diagonal,

			float* offset_position,

			bool use_dense_loop_over_voxels,
			int calc_radius_neigh,
			float4* kdtree_nodes,
			int     kdtree_N
		) {
			// Route to appropriate conversion function based on grid type
			if (grid.type == VDBParticleType::eNanoVDB || grid.type == VDBParticleType::eOpenVDB) {
				// Sparse grid conversion using hash-based voxel accumulation
				convert_to_sparse_grid_cpu(
					pos_particles,
					particle_ids,
					radius_particles,
					value_particles,
					num_particles,
					particle_radius_multiplier,
					bbox_min_orig,
					bbox_size_orig,
					bbox_dim,
					offset_position,
					filter_min,
					filter_max,
					bbox_sphere_pos,
					bbox_sphere_r,
					anim_type,
					frame_req,
					frame,
					use_simple_density,
					particle_radius_const,

					bbox_x_min_norm,
					bbox_x_max_norm,
					bbox_y_min_norm,
					bbox_y_max_norm,
					bbox_z_min_norm,
					bbox_z_max_norm,
					scale_space_diagonal,

					grid.sparse_grid.get(),
					min_value,
					max_value,
					particles_count
				);
			}
			else if (grid.type == VDBParticleType::eDense) {
#ifdef WITH_CUDAKDTREE
				if (use_dense_loop_over_voxels && calc_radius_neigh > 0
					    && kdtree_nodes != nullptr && kdtree_N > 0) {
					convert_to_dense_grid_loop_over_voxels_cpu(
						kdtree_nodes,
						kdtree_N,
						radius_particles,
						value_particles,
						particle_radius_multiplier,
						bbox_min_orig,
						bbox_size_orig,
						bbox_dim,
						dense_type,
						dense_norm,
						block_name_id,
						offset_position,
						use_simple_density,
						particle_radius_const,
						scale_space_diagonal,
						calc_radius_neigh,
						grid.dense_grid.get(),
						min_value,
						max_value,
						particles_count
					);
				}
				else
#endif
				// Dense grid conversion using SPH kernel splatting
				convert_to_dense_grid_cpu(
					pos_particles,
					particle_ids,
					radius_particles,
					value_particles,
					num_particles,
					particle_radius_multiplier,
					bbox_min_orig,
					bbox_size_orig,
					bbox_dim,

					dense_type,
					dense_norm,
					block_name_id,
					offset_position,
					filter_min,
					filter_max,
					bbox_sphere_pos,
					bbox_sphere_r,
					anim_type,
					frame_req,
					frame,
					use_simple_density,
					particle_radius_const,

					bbox_x_min_norm,
					bbox_x_max_norm,
					bbox_y_min_norm,
					bbox_y_max_norm,
					bbox_z_min_norm,
					bbox_z_max_norm,
					scale_space_diagonal,

					grid.dense_grid.get(),
					min_value,
					max_value,
					particles_count
				);
			}
			else if (grid.type == VDBParticleType::eRawParticles) {
				// Raw particle storage without rasterization
				convert_to_raw_particles_cpu(
					pos_particles,
					particle_ids,
					radius_particles,
					value_particles,
					num_particles,
					particle_radius_multiplier,
					bbox_min_orig,
					bbox_size_orig,
					bbox_dim,
					object_size,
					block_name_id,
					offset_position,
					filter_min,
					filter_max,
					bbox_sphere_pos,
					bbox_sphere_r,
					anim_type,
					frame_req,
					frame,
					use_simple_density,
					particle_radius_const,

					bbox_x_min_norm,
					bbox_x_max_norm,
					bbox_y_min_norm,
					bbox_y_max_norm,
					bbox_z_min_norm,
					bbox_z_max_norm,
					scale_space_diagonal,

					grid.raw_particles,
					min_value,
					max_value,
					particles_count
				);
			}
		}

		} // namespace kernel
	} // namespace vdb
} // namespace common