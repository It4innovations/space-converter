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

#include <float.h>
#include <vector>
#include <cmath>
#include "data_common.h"
#include "../utility/dense_utility.h"

namespace common {
	namespace vdb {

		class VoxelSparseManager;
		class VDBParticles;
		class VoxelDenseManager;

		namespace kernel {

#ifdef __CUDACC__
			__host__ __device__
#endif            
				inline double get_particle_radius(
					uint64_t pid,
					int bbox_dim,
					int* bbox_min_orig,
					double bbox_size_orig,
					double scale_space_diagonal,
					float particle_fix_size,
					//int particle_type,
					const float* radius_particles,
					double radius_particle_const
				) {
				if (radius_particle_const > 0.0) {
					return radius_particle_const;
				}

				// Conversion factor from world space to voxel space
				double norm_fac = (double)bbox_dim / scale_space_diagonal;

				// Get SPH properties for this particle
				//double hsml = get_particle_hsml(pid);  // Smoothing length
				//double mass = get_particle_mass(pid);
				//double rho = get_particle_rho(pid);    // Density

				double radius = 0.0;

				// Use precomputed radius from k-NN search if available (overrides hsml)
				if (radius_particles != nullptr) {
					radius = radius_particles[pid];
				}

				if (particle_fix_size != 0.0f) {
					radius *= particle_fix_size;
				}

				// Convert radius to voxel units
				double len_to_pix = norm_fac / bbox_size_orig;
				double radiusxyz_max = radius * len_to_pix;

				return radiusxyz_max;
			}

#ifdef __CUDACC__
			__host__ __device__
#endif            
				inline void fill_voxels(
					//common::vdb::VoxelDenseManager& grid,
					float* data_density,
					float* data_temp,
					size_t* offset,
					size_t* dims,

					size_t pid,
					float value,
					int bbox_dim,
					int* bbox_min_orig,
					double bbox_size_orig,
					double scale_space_diagonal,
					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					float particle_fix_size,
					const float* radius_particles,
					int block_name_id,
					double* pos
				) {

				double norm_fac = (double)bbox_dim / scale_space_diagonal;
				double len2pix = norm_fac / bbox_size_orig;

				if (block_name_id == 0) { //only pos
					value = 1.0f;
				}

				double hsml = get_particle_radius(
					pid,
					bbox_dim,
					bbox_min_orig,
					bbox_size_orig,
					scale_space_diagonal,
					particle_fix_size,
					radius_particles,
					1.0
				);

				int iradiusx = static_cast<int>(hsml);
				int iradiusy = static_cast<int>(hsml);
				int iradiusz = static_cast<int>(hsml);

				///////////////////////////////////////////////////////////////////////////////
				double dpx = ((double)pos[0] - (double)bbox_min_orig[0]) * len2pix;
				double dpy = ((double)pos[1] - (double)bbox_min_orig[1]) * len2pix;
				double dpz = ((double)pos[2] - (double)bbox_min_orig[2]) * len2pix;

				int px = static_cast<int>(dpx);
				int py = static_cast<int>(dpy);
				int pz = static_cast<int>(dpz);

				///////////////////////////////////////////////////////////////////////////////

				for (int sx = px - iradiusx; sx <= px + iradiusx; sx++) {
					for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
						for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

							int osx = sx - offset[0];
							int osy = sy - offset[1];
							int osz = sz - offset[2];

							if (osx < 0 || osy < 0 || osz < 0)
								continue;

							if (osx > dims[0] - 1 || osy > dims[1] - 1 || osz > dims[2] - 1)
								continue;

							double density = 1.0;

							double dx = dpx - sx;
							double dy = dpy - sy;
							double dz = dpz - sz;
#ifdef __CUDA_ARCH__
							double distance_norm = sqrt(dx * dx + dy * dy + dz * dz);
#else                        
							double distance_norm = std::sqrt(dx * dx + dy * dy + dz * dz);
#endif

							double W = 0.0;
							double h = hsml;

							if (iradiusx == 0 && iradiusy == 0 && iradiusz == 0) {
								W = 1.0; // full value
							}
							else {
								double q = distance_norm / h;

								W = utility::dense::sph_kernel::W(dense_type, q, 1.0 / h);
							}

							density = W;

							//final density
							float d = density * value;

							double norm = 0.0;
							if (dense_norm == common::SpaceData::DenseNorm::eNone) {
								norm = 0.0;
							}
							else if (dense_norm == common::SpaceData::DenseNorm::eCount) {
								norm = 1.0; // count
							}
							else if (dense_norm == common::SpaceData::DenseNorm::eSPHInterpolation) {
								norm = density;
							}

							//final norm
							float n = norm;

							//x + y * dims[0] + z * dims[0] * dims[1];
							size_t gindex = osx + osy * dims[0] + osz * dims[0] * dims[1];
#ifdef __CUDA_ARCH__
							atomicAdd(data_density + gindex, d);
#else
#pragma omp atomic
							data_density[gindex] += d;
#endif


							if (data_temp != nullptr) {
#ifdef __CUDA_ARCH__
								atomicAdd(data_temp + gindex, n);
#else
#pragma omp atomic
								data_temp[gindex] += n;
#endif

							}
						}
					}
				}
			}

			/**
			 * @brief Process a single particle for grid conversion
			 *
			 * This function can be called from both CPU and CUDA kernels.
			 * Returns true if particle should be processed, false if it should be skipped.
			 */
#ifdef __CUDACC__
			__host__ __device__
#endif
				inline bool process_particle(
					size_t cached_idx,
					const float* pos_particles,
					const float* radius_particles,
					const float* value_particles,
					const float* offset_position,
					int* bbox_min_orig,
					double bbox_size_orig,
					int bbox_dim,
					double scale_space_diagonal,
					float particle_fix_size,
					double bbox_x_min_norm,
					double bbox_x_max_norm,
					double bbox_y_min_norm,
					double bbox_y_max_norm,
					double bbox_z_min_norm,
					double bbox_z_max_norm,
					float filter_min,
					float filter_max,
					const float* bbox_sphere_pos,
					float bbox_sphere_r,
					common::SpaceData::AnimType anim_type,
					int frame_req,
					int frame,
					bool use_simple_density,
					bool is_dense_grid,
					// Output parameters
					double* out_Pos,
					double& out_px_norm,
					double& out_py_norm,
					double& out_pz_norm,
					int& out_px,
					int& out_py,
					int& out_pz,
					float& out_v_orig
				) {
				// Filter by animation frame if extracting specific frame
				if (anim_type == common::SpaceData::AnimType::eFrameExtract) {
					if (frame_req != frame)
						return false;
				}

				// Get position from cached data
				out_Pos[0] = pos_particles[cached_idx * 3 + 0];
				out_Pos[1] = pos_particles[cached_idx * 3 + 1];
				out_Pos[2] = pos_particles[cached_idx * 3 + 2];

				// Apply position offset if specified
				if (offset_position[0] != 0.0f || offset_position[1] != 0.0f || offset_position[2] != 0.0f) {
					out_Pos[0] -= offset_position[0];
					out_Pos[1] -= offset_position[1];
					out_Pos[2] -= offset_position[2];
				}

				// Normalize position to [0,1] range based on bounding box
				out_px_norm = ((double)out_Pos[0] - (double)bbox_min_orig[0]) / bbox_size_orig;
				out_py_norm = ((double)out_Pos[1] - (double)bbox_min_orig[1]) / bbox_size_orig;
				out_pz_norm = ((double)out_Pos[2] - (double)bbox_min_orig[2]) / bbox_size_orig;

				// Check if particle is within bounding box
				// Dense grids need extended bounds (particle radius) for proper splatting
				if (is_dense_grid) {
					// Get particle radius for boundary check with tolerance
					double radiusxyz_max = get_particle_radius(
						cached_idx,
						bbox_dim,
						bbox_min_orig,
						bbox_size_orig,
						scale_space_diagonal,
						particle_fix_size,
						radius_particles,
						1.0
					);

					if (out_px_norm + radiusxyz_max < bbox_x_min_norm || out_px_norm - radiusxyz_max > bbox_x_max_norm)
						return false;

					if (out_py_norm + radiusxyz_max < bbox_y_min_norm || out_py_norm - radiusxyz_max > bbox_y_max_norm)
						return false;

					if (out_pz_norm + radiusxyz_max < bbox_z_min_norm || out_pz_norm - radiusxyz_max > bbox_z_max_norm)
						return false;
				}
				else {
					// Simple point-in-box test for sparse grids
					if (out_px_norm < bbox_x_min_norm || out_px_norm > bbox_x_max_norm)
						return false;

					if (out_py_norm < bbox_y_min_norm || out_py_norm > bbox_y_max_norm)
						return false;

					if (out_pz_norm < bbox_z_min_norm || out_pz_norm > bbox_z_max_norm)
						return false;
				}

				// Convert normalized position to voxel coordinates
				out_px = static_cast<int>((double)out_px_norm * (double)bbox_dim / scale_space_diagonal);
				out_py = static_cast<int>((double)out_py_norm * (double)bbox_dim / scale_space_diagonal);
				out_pz = static_cast<int>((double)out_pz_norm * (double)bbox_dim / scale_space_diagonal);

				// Get particle value (density, temperature, etc.)
				// Note: This will need to be provided by the caller since it requires data access
				out_v_orig = value_particles[cached_idx];

				if (use_simple_density) {
					// Simple density mode: treat all particles with uniform weight
					out_v_orig = 1.0f;
					return true;  // Skip remaining checks if using simple density
				}

				// Apply value range filter (caller must set out_v_orig first)
				if (out_v_orig < filter_min || out_v_orig > filter_max)
					return false;

				// Apply spherical bounding filter if specified (radius > 0)
				if (bbox_sphere_r > 0.0f)
				{
					// Compute squared distance from particle to sphere center
					// Using squared distance avoids expensive sqrt operation
					double distSquared =
						(out_Pos[0] - bbox_sphere_pos[0]) * (out_Pos[0] - bbox_sphere_pos[0]) +
						(out_Pos[1] - bbox_sphere_pos[1]) * (out_Pos[1] - bbox_sphere_pos[1]) +
						(out_Pos[2] - bbox_sphere_pos[2]) * (out_Pos[2] - bbox_sphere_pos[2]);

					// Skip particles outside the spherical region
					if (distSquared > (bbox_sphere_r * bbox_sphere_r)) {
						return false;
					}
				}

				return true;
			};

			//////////////////

#ifdef WITH_GPU_CUDA

			/**
			 * @brief GPU implementation of bounding box calculation
			 *
			 * Launches CUDA kernel to compute min/max coordinates across all particles.
			 * Implemented in convert_vdb_kernel.cu
			 */
			void find_bbox_gpu(
				const float* d_pos_particles,
				size_t num_particles,
				const float* offset_position,
				float* bbox_min,
				float* bbox_max
			);

			///**
			// * @brief GPU implementation of sparse grid conversion
			// *
			// * Converts particles to sparse grid on GPU using CUDA.
			// * Implemented in convert_vdb_kernel.cu
			// */
			//void convert_to_sparse_grid_gpu(
			//	const float* pos_particles,
			//	const float* radius_particles,
			//	size_t num_particles,
			//	float particle_fix_size,
			//	int* bbox_min_orig,
			//	double bbox_size_orig,
			//	int bbox_dim,
			//	float* offset_position,
			//	float filter_min,
			//	float filter_max,
			//	float* bbox_sphere_pos,
			//	float bbox_sphere_r,
			//	common::SpaceData::AnimType anim_type,
			//	int frame_req,
			//	int frame,
			//	bool use_simple_density,
			//	VoxelSparseManager* voxel_manager,
			//	float& min_value,
			//	float& max_value,
			//	size_t& particles_count
			//);

			///**
			// * @brief GPU implementation of dense grid conversion
			// *
			// * Converts particles to dense grid on GPU using CUDA with SPH splatting.
			// * Implemented in convert_vdb_kernel.cu
			// */
			//void convert_to_dense_grid_gpu(
			//	const float* pos_particles,
			//	const float* radius_particles,
			//	size_t num_particles,
			//	float particle_fix_size,
			//	int* bbox_min_orig,
			//	double bbox_size_orig,
			//	int bbox_dim,
			//	double scale_space_diagonal,
			//	common::SpaceData::DenseType dense_type,
			//	common::SpaceData::DenseNorm dense_norm,
			//	int block_name_id,
			//	float* offset_position,
			//	float filter_min,
			//	float filter_max,
			//	float* bbox_sphere_pos,
			//	float bbox_sphere_r,
			//	common::SpaceData::AnimType anim_type,
			//	int frame_req,
			//	int frame,
			//	bool use_simple_density,
			//	VoxelDenseManager& grid,
			//	float& min_value,
			//	float& max_value,
			//	size_t& particles_count
			//);

			/**
			 * @brief Main GPU dispatcher for converting I/O library data to grids
			 *
			 * Routes to appropriate GPU conversion function based on grid type.
			 * Implemented in convert_vdb_kernel.cu
			 */
			void convert_iolib_to_grid_gpu(
				const float* pos_particles,
				const size_t* particle_ids,
				const float* radius_particles,
				const float* value_particles,
				size_t num_particles,
				float particle_fix_size,
				std::string grid_name,
				float grid_transform,
				float* bbox_min,
				float* bbox_max,
				int bbox_dim,
				int* bbox_min_orig,
				double bbox_size_orig,
				common::SpaceData::ExtractedType extracted_type,
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

				double bbox_x_min_norm,
				double bbox_x_max_norm,
				double bbox_y_min_norm,
				double bbox_y_max_norm,
				double bbox_z_min_norm,
				double bbox_z_max_norm,
				double scale_space_diagonal,

				float* offset_position
			);


#endif // WITH_GPU_CUDA

			/**
			 * @brief CPU implementation of bounding box calculation using OpenMP
			 *
			 * Parallel reduction to compute min/max coordinates across all particles.
			 * Implemented in convert_vdb_kernel.cpp
			 */
			void find_bbox_cpu(
				const float* pos_particles,
				size_t num_particles,
				const float* offset_position,
				float* bbox_min,
				float* bbox_max
			);

			void convert_iolib_to_grid_cpu(
				const float* pos_particles,
				const size_t* particle_ids,
				const float* radius_particles,
				const float* value_particles,
				size_t num_particles,
				float particle_fix_size,
				std::string grid_name,
				float grid_transform,
				float* bbox_min,
				float* bbox_max,
				int bbox_dim,
				int* bbox_min_orig,
				double bbox_size_orig,
				common::SpaceData::ExtractedType extracted_type,
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

				double bbox_x_min_norm,
				double bbox_x_max_norm,
				double bbox_y_min_norm,
				double bbox_y_max_norm,
				double bbox_z_min_norm,
				double bbox_z_max_norm,
				double scale_space_diagonal,

				float* offset_position
			);



		} // namespace kernel
	} // namespace vdb
} // namespace common