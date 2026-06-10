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
#include "../utility/gpu_utility.h"
#include "../utility/gpu_logging.h"
#include "sparse_common.h"
#include "dense_common.h"
#include "data_common.h"
#include <float.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <limits>
#include <algorithm>

#ifdef WITH_CUDAKDTREE
#include "cukd/traverse-stack-free.h"
#include "cukd/builder.h"
#include "cukd/knn.h"
#include "cukd/cukd-math.h"
using cukd::divRoundUp;
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace space_converter {
	namespace common {
		namespace vdb {
			namespace kernel {

				// Using declarations for cleaner code
				using common::vdb::sparse::packCoord3;

				/**
				 * @brief CUDA kernel to find bounding box on GPU
				 *
				 * Each thread processes particles and atomically updates min/max values.
				 */
				__global__ void find_bbox_kernel_cuda(
					const float* pos_particles,
					size_t num_particles,
					const float* offset_position,
					float* bbox_min,
					float* bbox_max
				) {
					// Thread index
					size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

					// Thread-local min/max
					float local_min_x = FLT_MAX, local_min_y = FLT_MAX, local_min_z = FLT_MAX;
					float local_max_x = -FLT_MAX, local_max_y = -FLT_MAX, local_max_z = -FLT_MAX;

					// Process particles with stride
					for (size_t i = idx; i < num_particles; i += blockDim.x * gridDim.x) {
						// Load position (x, y, z are interleaved: [x0, y0, z0, x1, y1, z1, ...])
						float px = pos_particles[i * 3 + 0] - offset_position[0];
						float py = pos_particles[i * 3 + 1] - offset_position[1];
						float pz = pos_particles[i * 3 + 2] - offset_position[2];

						// Update local min/max
						if (px < local_min_x) local_min_x = px;
						if (py < local_min_y) local_min_y = py;
						if (pz < local_min_z) local_min_z = pz;

						if (px > local_max_x) local_max_x = px;
						if (py > local_max_y) local_max_y = py;
						if (pz > local_max_z) local_max_z = pz;
					}

					// Atomic min/max updates using compare-and-swap pattern for floats
					// This ensures correct floating-point atomic operations
					if (local_min_x < FLT_MAX) {
						int* addr_as_int = (int*)&bbox_min[0];
						int old = *addr_as_int, assumed;
						do {
							assumed = old;
							old = atomicCAS(addr_as_int, assumed,
								__float_as_int(fminf(__int_as_float(assumed), local_min_x)));
						} while (assumed != old);
					}

					if (local_min_y < FLT_MAX) {
						int* addr_as_int = (int*)&bbox_min[1];
						int old = *addr_as_int, assumed;
						do {
							assumed = old;
							old = atomicCAS(addr_as_int, assumed,
								__float_as_int(fminf(__int_as_float(assumed), local_min_y)));
						} while (assumed != old);
					}

					if (local_min_z < FLT_MAX) {
						int* addr_as_int = (int*)&bbox_min[2];
						int old = *addr_as_int, assumed;
						do {
							assumed = old;
							old = atomicCAS(addr_as_int, assumed,
								__float_as_int(fminf(__int_as_float(assumed), local_min_z)));
						} while (assumed != old);
					}

					if (local_max_x > -FLT_MAX) {
						int* addr_as_int = (int*)&bbox_max[0];
						int old = *addr_as_int, assumed;
						do {
							assumed = old;
							old = atomicCAS(addr_as_int, assumed,
								__float_as_int(fmaxf(__int_as_float(assumed), local_max_x)));
						} while (assumed != old);
					}

					if (local_max_y > -FLT_MAX) {
						int* addr_as_int = (int*)&bbox_max[1];
						int old = *addr_as_int, assumed;
						do {
							assumed = old;
							old = atomicCAS(addr_as_int, assumed,
								__float_as_int(fmaxf(__int_as_float(assumed), local_max_y)));
						} while (assumed != old);
					}

					if (local_max_z > -FLT_MAX) {
						int* addr_as_int = (int*)&bbox_max[2];
						int old = *addr_as_int, assumed;
						do {
							assumed = old;
							old = atomicCAS(addr_as_int, assumed,
								__float_as_int(fmaxf(__int_as_float(assumed), local_max_z)));
						} while (assumed != old);
					}
				}

				/**
				 * @brief GPU implementation of bounding box calculation
				 *
				 * Launches CUDA kernel to compute min/max coordinates across all particles.
				 */
				void find_bbox_gpu(
					const float* d_pos_particles,
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

					// Allocate device memory for offset and results
					float* d_offset, * d_bbox_min, * d_bbox_max;
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_offset, 3 * sizeof(float)));
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_bbox_min, 3 * sizeof(float)));
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_bbox_max, 3 * sizeof(float)));

					// Initialize GPU results with extreme values
					float init_min[3] = { bbox_min[0], bbox_min[1], bbox_min[2] };
					float init_max[3] = { bbox_max[0], bbox_max[1], bbox_max[2] };
					CUDA_CHECK_ERROR(cudaMemcpy(d_offset, offset_position, 3 * sizeof(float), cudaMemcpyHostToDevice));
					CUDA_CHECK_ERROR(cudaMemcpy(d_bbox_min, init_min, 3 * sizeof(float), cudaMemcpyHostToDevice));
					CUDA_CHECK_ERROR(cudaMemcpy(d_bbox_max, init_max, 3 * sizeof(float), cudaMemcpyHostToDevice));

					// Launch kernel
					int blockSize = 256;
					int numBlocks = (num_particles + blockSize - 1) / blockSize;
					numBlocks = min(numBlocks, 1024); // Limit grid size

					GPU_KERNEL_TIME_START(find_bbox_kernel_cuda);
					find_bbox_kernel_cuda << <numBlocks, blockSize >> > (
						d_pos_particles,
						num_particles,
						d_offset,
						d_bbox_min,
						d_bbox_max
						);
					GPU_KERNEL_TIME_END(find_bbox_kernel_cuda);
					//CUDA_CHECK_LAST_ERROR();
					CUDA_SYNC_CHECK();

					// Copy results back to host
					CUDA_CHECK_ERROR(cudaMemcpy(bbox_min, d_bbox_min, 3 * sizeof(float), cudaMemcpyDeviceToHost));
					CUDA_CHECK_ERROR(cudaMemcpy(bbox_max, d_bbox_max, 3 * sizeof(float), cudaMemcpyDeviceToHost));

					// Free device memory
					CUDA_CHECK_ERROR(cudaFree(d_offset));
					CUDA_CHECK_ERROR(cudaFree(d_bbox_min));
					CUDA_CHECK_ERROR(cudaFree(d_bbox_max));

					// Wait for GPU to finish
					//CUDA_SYNC_CHECK();
				}

				/**
				 * @brief CUDA kernel for sparse grid conversion
				 *
				 * Each thread processes particles and atomically accumulates values into voxels.
				 * Uses a hash map approach with packed voxel coordinates.
				 */
				__global__ void convert_to_sparse_grid_kernel_cuda(
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

					uint64_t* voxel_keys,      // Output: packed voxel coordinates
					float* voxel_values,       // Output: voxel values
					uint64_t* particle_count   // Output: number of processed particles
				) {
					size_t idx_c = blockIdx.x * blockDim.x + threadIdx.x;

					if (idx_c >= num_particles) return;

					// Cached particle ID for indexing
					size_t idx = particle_ids[idx_c];

					// Particle processing outputs
					double Pos[3];
					double px_norm, py_norm, pz_norm;
					int px, py, pz;
					float v_orig = 0.0f;

					//    float filter_max,
					//    const float* bbox_sphere_pos,
					//    float bbox_sphere_r,
					//    common::SpaceData::AnimType anim_type,
					//    int frame_req,
					//    int frame,
					//    bool use_simple_density,
					//    double particle_radius_const,

					//    bool is_dense_grid,
					//    // Output parameters
					//    double* out_Pos,
					//    double& out_px_norm,
					//    double& out_py_norm,
					//    double& out_pz_norm,
					//    int& out_px,
					//    int& out_py,
					//    int& out_pz,
					//    float& out_v_orig

					// Check if particle passes filters and compute voxel coordinates
					bool should_process = process_particle(
						idx, //size_t cached_idx,
						pos_particles, //    const float* pos_particles,
						radius_particles,//    const float* radius_particles,
						value_particles,//    const float* value_particles,
						offset_position,//    const float* offset_position,
						bbox_min_orig,//    int* bbox_min_orig,
						bbox_size_orig,//    double bbox_size_orig,
						bbox_dim,//    int bbox_dim,
						scale_space_diagonal,//    double scale_space_diagonal,
						particle_radius_multiplier,//    float particle_radius_multiplier,
						bbox_x_min_norm,//    double bbox_x_min_norm,
						bbox_x_max_norm,//    double bbox_x_max_norm,
						bbox_y_min_norm,//    double bbox_y_min_norm,
						bbox_y_max_norm,//    double bbox_y_max_norm,
						bbox_z_min_norm, //    double bbox_z_min_norm,
						bbox_z_max_norm,//    double bbox_z_max_norm,
						filter_min,//    float filter_min,
						filter_max,
						bbox_sphere_pos,
						bbox_sphere_r,
						anim_type,
						frame_req,
						frame,
						use_simple_density,
						particle_radius_const,
						false, // is_dense_grid = false for sparse
						// Outputs
						Pos, px_norm, py_norm, pz_norm,
						px, py, pz, v_orig
					);

					// Store results
					if (should_process) {
						voxel_keys[idx] = packCoord3(px, py, pz);
						voxel_values[idx] = v_orig;
						atomicAdd((unsigned long long*)particle_count, 1ULL);
					}
				}

				/**
				 * @brief CUDA kernel for dense grid conversion with SPH splatting
				 *
				 * Each thread processes particles and splats them into the dense grid using
				 * SPH kernel functions. Uses atomic operations to accumulate density values.
				 */
				__global__ void convert_to_dense_grid_kernel_cuda(
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

					float* grid_data_density,
					float* grid_data_temp,
					size_t* grid_offset,
					size_t* grid_dims,

					uint64_t* particle_count  // Output: number of processed particles
				) {
					size_t idx_c = blockIdx.x * blockDim.x + threadIdx.x;

					if (idx_c >= num_particles) return;

					// Cached particle ID for indexing
					size_t idx = particle_ids[idx_c];

					// Particle processing outputs
					double Pos[3];
					double px_norm, py_norm, pz_norm;
					int px, py, pz;
					float v_orig = 0.0f;

					// Check if particle passes filters and compute voxel coordinates
					bool should_process = process_particle(
						idx,
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

					if (should_process) {
						atomicAdd((unsigned long long*)particle_count, 1ULL);
						// Splat particle into dense grid using SPH kernel function
						fill_voxels(
							grid_data_density,
							grid_data_temp,
							grid_offset,
							grid_dims,
							idx,
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
				}

				/**
				 * @brief GPU implementation of sparse grid conversion
				 *
				 * Converts particles to sparse grid on GPU using CUDA.
				 */
				void convert_to_sparse_grid_gpu(
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
					if (num_particles == 0) {
						min_value = 0.0f;
						max_value = 0.0f;
						particles_count = 0;
						return;
					}

					// Attempt to dynamic_cast to VoxelGPUManagerSortReduce
					common::vdb::sparse::VoxelGPUManagerSortReduce* voxel_manager_gpu = dynamic_cast<common::vdb::sparse::VoxelGPUManagerSortReduce*>(voxel_manager);
					if (!voxel_manager_gpu) {
						printf("[convert_to_sparse_grid_gpu] Error: Incompatible manager type\n");
						return;
					}

					//// Allocate device memory
					//uint64_t* d_voxel_keys;
					//float* d_voxel_values;
					//size_t* d_particle_valid;
					//float* d_particle_values;
					//
					//CUDA_CHECK_ERROR(CUDA_MALLOC(&d_voxel_keys, num_particles * sizeof(uint64_t)));
					//CUDA_CHECK_ERROR(CUDA_MALLOC(&d_voxel_values, num_particles * sizeof(float)));
					//CUDA_CHECK_ERROR(CUDA_MALLOC(&d_particle_valid, num_particles * sizeof(size_t)));
					//CUDA_CHECK_ERROR(CUDA_MALLOC(&d_particle_values, num_particles * sizeof(float)));

					// Reset particle counter before kernel launch
					CUDA_CHECK_ERROR(cudaMemset(voxel_manager_gpu->d_particle_count, 0, sizeof(uint64_t)));

					// Upload per-call host parameters into persistent device buffers
					CUDA_CHECK_ERROR(cudaMemcpy(voxel_manager_gpu->d_bbox_min_orig, bbox_min_orig, 3 * sizeof(int), cudaMemcpyHostToDevice));
					CUDA_CHECK_ERROR(cudaMemcpy(voxel_manager_gpu->d_offset_position, offset_position, 3 * sizeof(float), cudaMemcpyHostToDevice));
					CUDA_CHECK_ERROR(cudaMemcpy(voxel_manager_gpu->d_bbox_sphere_pos, bbox_sphere_pos, 3 * sizeof(float), cudaMemcpyHostToDevice));

					// Launch kernel
					int blockSize = 256;
					int numBlocks = (num_particles + blockSize - 1) / blockSize;

					GPU_KERNEL_TIME_START(convert_to_sparse_grid_kernel_cuda);
					convert_to_sparse_grid_kernel_cuda << <numBlocks, blockSize >> > (
						pos_particles,
						particle_ids,
						radius_particles,
						value_particles,
						num_particles,
						particle_radius_multiplier,
						voxel_manager_gpu->d_bbox_min_orig,
						bbox_size_orig,
						bbox_dim,
						voxel_manager_gpu->d_offset_position,
						filter_min,
						filter_max,
						voxel_manager_gpu->d_bbox_sphere_pos,
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

						voxel_manager_gpu->d_keys,
						voxel_manager_gpu->d_vals,
						voxel_manager_gpu->d_particle_count
						//d_particle_valid,
						//d_particle_values
						);
					GPU_KERNEL_TIME_END(convert_to_sparse_grid_kernel_cuda);
					//CUDA_CHECK_LAST_ERROR();
					CUDA_SYNC_CHECK();

					// Read back the number of particles that passed the should_process filter
					uint64_t raw_count = 0;
					CUDA_CHECK_ERROR(cudaMemcpy(&raw_count, voxel_manager_gpu->d_particle_count, sizeof(uint64_t), cudaMemcpyDeviceToHost));
					particles_count = static_cast<size_t>(raw_count);
					voxel_manager_gpu->update(raw_count);

					// Get min/max using CUB (operates on d_vals_out which contains reduced values)
					voxel_manager_gpu->find_min_max(min_value, max_value);
				}

				/**
				 * @brief GPU implementation of dense grid conversion
				 *
				 * Converts particles to dense grid on GPU using CUDA with SPH splatting.
				 */
				void convert_to_dense_grid_gpu(
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
					if (num_particles == 0) {
						min_value = 0.0f;
						max_value = 0.0f;
						particles_count = 0;
						return;
					}

					// Attempt to dynamic_cast to VoxelGPUManagerSortReduce
					common::vdb::dense::VoxelGPUDenseManager* voxel_manager_gpu = dynamic_cast<common::vdb::dense::VoxelGPUDenseManager*>(grid);
					if (!voxel_manager_gpu) {
						printf("[convert_to_dense_grid_gpu] Error: Incompatible manager type\n");
						return;
					}

					//                // Allocate device memory for grid and particle data
					//                float* d_grid_data_density;
					//                float* d_grid_data_temp;
					//                size_t* d_grid_offset;
					//                size_t* d_grid_dims;
					//                size_t* d_particle_valid;
					//                float* d_particle_values;
					//                
					//                size_t grid_size = grid.size();
					//                
					//                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_grid_data_density, grid_size * sizeof(float)));
					//                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_grid_data_temp, grid_size * sizeof(float)));
					//                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_grid_offset, 3 * sizeof(size_t)));
					//                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_grid_dims, 3 * sizeof(size_t)));
					//                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_particle_valid, num_particles * sizeof(size_t)));
					//                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_particle_values, num_particles * sizeof(float)));
					//                
					//                // Initialize grid on device
					//                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_data_density, grid.data_density.data(), grid_size * sizeof(float), cudaMemcpyHostToDevice));
					//#ifndef WITH_NO_DATA_TEMP
					//                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_data_temp, grid.data_temp.data(), grid_size * sizeof(float), cudaMemcpyHostToDevice));
					//#else
					//                d_grid_data_temp = nullptr;
					//#endif
					//                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_offset, grid.offset, 3 * sizeof(size_t), cudaMemcpyHostToDevice));
					//                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_dims, grid.dims, 3 * sizeof(size_t), cudaMemcpyHostToDevice));

									// Copy host pointers to device memory
					int* d_bbox_min_orig;
					float* d_offset_position;
					float* d_bbox_sphere_pos;
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_bbox_min_orig, 3 * sizeof(int)));
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_offset_position, 3 * sizeof(float)));
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_bbox_sphere_pos, 3 * sizeof(float)));
					CUDA_CHECK_ERROR(cudaMemcpy(d_bbox_min_orig, bbox_min_orig, 3 * sizeof(int), cudaMemcpyHostToDevice));
					CUDA_CHECK_ERROR(cudaMemcpy(d_offset_position, offset_position, 3 * sizeof(float), cudaMemcpyHostToDevice));
					CUDA_CHECK_ERROR(cudaMemcpy(d_bbox_sphere_pos, bbox_sphere_pos, 3 * sizeof(float), cudaMemcpyHostToDevice));

					//sync d_offset
					CUDA_CHECK_ERROR(cudaMemcpy(voxel_manager_gpu->d_offset, voxel_manager_gpu->offset, 3 * sizeof(size_t), cudaMemcpyHostToDevice));

					// Launch kernel
					int blockSize = 256;
					int numBlocks = (num_particles + blockSize - 1) / blockSize;

					GPU_KERNEL_TIME_START(convert_to_dense_grid_kernel_cuda);
					convert_to_dense_grid_kernel_cuda << <numBlocks, blockSize >> > (
						pos_particles, //    const float* pos_particles,
						particle_ids, //const size_t* particle_ids
						radius_particles, //    const float* radius_particles,
						value_particles, //    const float* value_particles,
						num_particles, //    size_t num_particles,
						particle_radius_multiplier, //    float particle_radius_multiplier,
						d_bbox_min_orig, //    int* bbox_min_orig,
						bbox_size_orig,//    double bbox_size_orig,
						bbox_dim, //    int bbox_dim,

						dense_type, //    common::SpaceData::DenseType dense_type,
						dense_norm, //    common::SpaceData::DenseNorm dense_norm,
						block_name_id, //    int block_name_id,
						d_offset_position, //    float* offset_position,
						filter_min, //    float filter_min,
						filter_max,//    float filter_max,
						d_bbox_sphere_pos, //    float* bbox_sphere_pos,
						bbox_sphere_r, //float bbox_sphere_r,
						anim_type, //    common::SpaceData::AnimType anim_type,
						frame_req, //    int frame_req,
						frame, //    int frame,
						use_simple_density, //    bool use_simple_density,
						particle_radius_const, //    double particle_radius_const,

						bbox_x_min_norm, //    double bbox_x_min_norm,
						bbox_x_max_norm,  //    double bbox_x_max_norm,
						bbox_y_min_norm,  //    double bbox_y_min_norm,
						bbox_y_max_norm, //    double bbox_y_max_norm,
						bbox_z_min_norm,  //    double bbox_z_min_norm,
						bbox_z_max_norm, //    double bbox_z_max_norm,
						scale_space_diagonal, //    double scale_space_diagonal,

						voxel_manager_gpu->d_data_density, //    float* grid_data_density,
						voxel_manager_gpu->d_data_temp, //    float* grid_data_temp,
						voxel_manager_gpu->d_offset, //    size_t* grid_offset,
						voxel_manager_gpu->d_dims, //size_t* grid_dims,
						voxel_manager_gpu->d_particle_count //    uint64_t* particle_count
						);
					GPU_KERNEL_TIME_END(convert_to_dense_grid_kernel_cuda);
					//CUDA_CHECK_LAST_ERROR();
					CUDA_SYNC_CHECK();

					// Free temporary device memory
					CUDA_CHECK_ERROR(cudaFree(d_bbox_min_orig));
					CUDA_CHECK_ERROR(cudaFree(d_offset_position));
					CUDA_CHECK_ERROR(cudaFree(d_bbox_sphere_pos));

					//                // Copy results back to host
					//                CUDA_CHECK_ERROR(cudaMemcpy(grid.data_density.data(), d_grid_data_density, grid_size * sizeof(float), cudaMemcpyDeviceToHost));
					//#ifndef WITH_NO_DATA_TEMP
					//                CUDA_CHECK_ERROR(cudaMemcpy(grid.data_temp.data(), d_grid_data_temp, grid_size * sizeof(float), cudaMemcpyDeviceToHost));
					//#endif
					//                
					//                std::vector<size_t> h_particle_valid(num_particles);
					//                std::vector<float> h_particle_values(num_particles);
					//                CUDA_CHECK_ERROR(cudaMemcpy(h_particle_valid.data(), d_particle_valid, num_particles * sizeof(size_t), cudaMemcpyDeviceToHost));
					//                CUDA_CHECK_ERROR(cudaMemcpy(h_particle_values.data(), d_particle_values, num_particles * sizeof(float), cudaMemcpyDeviceToHost));
					//                
					//                // Free device memory
					//                CUDA_CHECK_ERROR(cudaFree(d_grid_data_density));
					//                if (d_grid_data_temp) CUDA_CHECK_ERROR(cudaFree(d_grid_data_temp));
					//                CUDA_CHECK_ERROR(cudaFree(d_grid_offset));
					//                CUDA_CHECK_ERROR(cudaFree(d_grid_dims));
					//                CUDA_CHECK_ERROR(cudaFree(d_particle_valid));
					//                CUDA_CHECK_ERROR(cudaFree(d_particle_values));
					//                
					//                // Compute statistics on CPU
					//                float min = FLT_MAX;
					//                float max = -FLT_MAX;
					//                size_t count = 0;
					//                
					//                for (size_t i = 0; i < num_particles; ++i) {
					//                    if (h_particle_valid[i]) {
					//                        count++;
					//                        float val = h_particle_values[i];
					//                        if (val < min) min = val;
					//                        if (val > max) max = val;
					//                    }
					//                }

									// Read back the number of particles that passed the should_process filter
					uint64_t raw_count = 0;
					CUDA_CHECK_ERROR(cudaMemcpy(&raw_count, voxel_manager_gpu->d_particle_count, sizeof(uint64_t), cudaMemcpyDeviceToHost));
					particles_count = static_cast<size_t>(raw_count);
					//min_value = (count > 0) ? min : 0.0f;
					//max_value = (count > 0) ? max : 0.0f;
					voxel_manager_gpu->find_min_max(min_value, max_value);
				}

				/**
				 * @brief Main GPU dispatcher for converting I/O library data to grids
				 *
				 * Routes to appropriate GPU conversion function based on grid type.
				 */
				void convert_iolib_to_grid_gpu(
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
					void* d_kdtree_nodes,
					int     kdtree_N
				) {
					// Route to appropriate conversion function based on grid type
					if (grid.type == VDBParticleType::eNanoVDB || grid.type == VDBParticleType::eOpenVDB || grid.type == VDBParticleType::eCUB) {
						// Sparse grid conversion using GPU hash-based voxel accumulation                    
						convert_to_sparse_grid_gpu(
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
							&& d_kdtree_nodes != nullptr && kdtree_N > 0) {
							convert_to_dense_grid_loop_over_voxels_gpu(
								d_kdtree_nodes,
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
							// Dense grid conversion using GPU SPH kernel splatting
							convert_to_dense_grid_gpu(
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
						// Raw particle storage - not yet implemented for GPU
						printf("Error: Raw particle storage not yet implemented for GPU\n");
					}
				}

				// ─────────────────────────────────────────────────────────────────────────────
				// Voxel-centric dense grid conversion – GPU kernels and launcher.
				// CPU implementation lives in cudakdtree_tool.cu (compiled by WITH_CUDAKDTREE).
				// Both require WITH_CUDAKDTREE; the GPU launcher additionally requires WITH_GPU_CUDA.
				// ─────────────────────────────────────────────────────────────────────────────
#ifdef WITH_CUDAKDTREE

		// cukd traits for float4 nodes: (x,y,z)=voxel-unit position, w=encoded particle idx.
				struct float4_particle_data_traits_gpu {
					using point_t = float3;
					using data_t = float4;
					enum { has_explicit_dim = false };

					static inline __both__ point_t get_point(const float4& n)
					{
						return make_float3(n.x, n.y, n.z);
					}

					static inline __both__ float get_coord(const float4& n, int d)
					{
						return d == 0 ? n.x : (d == 1 ? n.y : n.z);
					}

					static inline __both__ int  get_dim(const float4&) { return -1; }
					static inline __both__ void set_dim(float4&, int) {}
				};

				// ── Helper kernel: fill voxel query positions for a batch ─────────────────
				// Produces offset-adjusted world-space positions matching the tree coordinate system.
				__global__ void fill_voxel_queries_kernel(
					float3* d_queries,
					size_t  voxel_start,
					size_t  batch_size,
					size_t  dim0,
					size_t  dim1,
					size_t  off0,
					size_t  off1,
					size_t  off2,
					// Voxel index → offset-adjusted world position:  world = bbox_min + wx/len2pix
					float   bbox_min_x,
					float   bbox_min_y,
					float   bbox_min_z,
					double  inv_len2pix   // = bbox_size_orig * scale_space_diagonal / bbox_dim
				) {
					size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
					if (tid >= batch_size) return;
					size_t v = voxel_start + tid;
					size_t wx = (v % dim0) + off0;
					size_t wy = ((v / dim0) % dim1) + off1;
					size_t wz = (v / (dim0 * dim1)) + off2;
					float3 q;
					q.x = (float)((double)bbox_min_x + (double)wx * inv_len2pix);
					q.y = (float)((double)bbox_min_y + (double)wy * inv_len2pix);
					q.z = (float)((double)bbox_min_z + (double)wz * inv_len2pix);
					d_queries[tid] = q;
				}

				// ── KNN query kernel for float4 tree ──────────────────────────────────────
				__global__ void runQuery_float4_kernel(
					float4* tree,
					int       tree_N,
					uint64_t* candidateLists,
					int       k,
					float     maxRadius,
					float3* queries,
					int       numQueries
				) {
					size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
					if (tid >= (size_t)numQueries) return;
					float3 qp = queries[tid];
					cukd::FlexHeapCandidateList cl(candidateLists + (size_t)k * tid, k, maxRadius);
					cukd::stackFree::knn<cukd::FlexHeapCandidateList,
						float4,
						float4_particle_data_traits_gpu>(cl, qp, tree, tree_N);
				}

				// ── Accumulation kernel: SPH contributions from candidates ────────────────
				__global__ void accumulate_voxel_sph_kernel(
					float4* d_tree,
					uint64_t* d_cand,
					int       k,
					size_t    voxel_start,
					size_t    batch_size,
					size_t    dim0,
					size_t    dim1,
					size_t    dim2,
					const float* radius_particles,
					const float* value_particles,
					float     particle_radius_multiplier,
					double    particle_radius_const,
					int* d_bbox_min_orig,
					double    bbox_size_orig,
					int       bbox_dim,
					double    scale_space_diagonal,
					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					int       block_name_id,
					bool      use_simple_density,
					float* grid_data_density,
					float* grid_data_temp,
					uint64_t* d_particle_count
				) {
					size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
					if (tid >= batch_size) return;

					size_t v = voxel_start + tid;
					size_t lx = v % dim0;
					size_t ly = (v / dim0) % dim1;
					size_t lz = v / (dim0 * dim1);
					if (lz >= dim2) return;

					uint64_t* my_cand = d_cand + tid * (size_t)k;

					double density_sum = 0.0;
					double norm_sum = 0.0;
					bool   contributed = false;

					for (int i = 0; i < k; ++i) {
						// Decode packed entry: upper 32 bits = dist2 (float bits), lower 32 = tree index
						int   tree_j = (int)(uint32_t(my_cand[i]));
						float dist2 = cukd::uint_as_float((uint32_t)(my_cand[i] >> 32));
						if (tree_j < 0) continue;
						if (isinf(dist2)) continue;

						// Recover original particle index from w component (bit-cast)
						uint32_t orig_idx32;
						memcpy(&orig_idx32, &d_tree[tree_j].w, sizeof(uint32_t));
						size_t orig_idx = (size_t)orig_idx32;

						float val = use_simple_density ? 1.0f : value_particles[orig_idx];
						if (block_name_id == 0) val = 1.0f;

						// h in world units; W normalised in voxel units for consistency with
						// the particle-centric path (q is dimensionless, same either way).
						double h_world = 0.0;
						if (particle_radius_const != 0.0) {
							h_world = particle_radius_const;
						}
						else if (particle_radius_multiplier != 0.0f) {
							h_world = (double)radius_particles[orig_idx] * (double)particle_radius_multiplier;
						}
						else {
							h_world = (double)radius_particles[orig_idx];
						}

						if (h_world <= 0.0) continue;

						double norm_fac = (double)bbox_dim / scale_space_diagonal;
						double len2pix = norm_fac / bbox_size_orig;
						double h_voxel = h_world * len2pix;

						double dist_world = sqrt((double)dist2);
						double q_ratio = dist_world / h_world;  // = dist_voxel / h_voxel
						if (q_ratio > 1.0 + 1e-6) continue;

						double W = utility::dense::sph_kernel::W(dense_type, q_ratio, 1.0 / h_voxel);

						density_sum += W * (double)val;
						if (dense_norm == common::SpaceData::DenseNorm::eCount)
							norm_sum += 1.0;
						else if (dense_norm == common::SpaceData::DenseNorm::eSPHInterpolation)
							norm_sum += W;
						contributed = true;
					}

					if (density_sum != 0.0 || norm_sum != 0.0) {
						size_t gidx = lx + ly * dim0 + lz * dim0 * dim1;
						atomicAdd(&grid_data_density[gidx], (float)density_sum);
						if (grid_data_temp && norm_sum != 0.0)
							atomicAdd(&grid_data_temp[gidx], (float)norm_sum);
						if (contributed)
							atomicAdd((unsigned long long*)d_particle_count, 1ULL);
					}
				}

#ifdef WITH_GPU_CUDA
				/**
				 * @brief GPU voxel-centric dense conversion using a pre-built device KD-tree.
				 *
				 * The tree was built once by build_voxel_kdtree() in init_converter() and is
				 * stored in cache_manager.d_voxel_kdtree_per_ptype.  No rebuild happens here.
				 */
				void convert_to_dense_grid_loop_over_voxels_gpu(
					float4* d_tree,             // pre-built device float4 KD-tree
					int           tree_N,
					const float* d_radius_particles,
					const float* d_value_particles,
					float  particle_radius_multiplier,
					int* bbox_min_orig,
					double bbox_size_orig,
					int    bbox_dim,

					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					int    block_name_id,
					float* /*offset_position*/,       // not needed (tree already in offset-adj. space)
					bool   use_simple_density,
					double particle_radius_const,
					double scale_space_diagonal,

					int calc_radius_neigh,

					VoxelDenseManager* grid,
					float& min_value,
					float& max_value,
					size_t& particles_count
				) {
					if (tree_N <= 0 || d_tree == nullptr || !grid) {
						min_value = max_value = 0.0f;
						particles_count = 0;
						return;
					}

					common::vdb::dense::VoxelGPUDenseManager* voxel_gpu =
						dynamic_cast<common::vdb::dense::VoxelGPUDenseManager*>(grid);
					if (!voxel_gpu) {
						printf("[convert_to_dense_grid_loop_over_voxels_gpu] Error: Incompatible manager\n");
						return;
					}

					const size_t dim0 = grid->dims[0];
					const size_t dim1 = grid->dims[1];
					const size_t dim2 = grid->dims[2];
					const size_t off0 = grid->offset[0];
					const size_t off1 = grid->offset[1];
					const size_t off2 = grid->offset[2];

					// ── Upload bbox_min_orig for the accumulation kernel ──────────────────
					int* d_bbox_min_orig_dev = nullptr;
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_bbox_min_orig_dev, 3 * sizeof(int)));
					CUDA_CHECK_ERROR(cudaMemcpy(d_bbox_min_orig_dev, bbox_min_orig, 3 * sizeof(int), cudaMemcpyHostToDevice));

					// Sync dense manager metadata
					CUDA_CHECK_ERROR(cudaMemcpy(voxel_gpu->d_offset, voxel_gpu->offset, 3 * sizeof(size_t), cudaMemcpyHostToDevice));
					CUDA_CHECK_ERROR(cudaMemcpy(voxel_gpu->d_dims, voxel_gpu->dims, 3 * sizeof(size_t), cudaMemcpyHostToDevice));

					// ── Process voxels in batches ─────────────────────────────────────────
					const size_t total_voxels = dim0 * dim1 * dim2;
					const size_t batch_size = std::min((size_t)65536, total_voxels);
					const int    k = calc_radius_neigh;

					float3* d_queries = nullptr;
					uint64_t* d_cand = nullptr;
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_queries, batch_size * sizeof(float3)));
					CUDA_CHECK_ERROR(CUDA_MALLOC(&d_cand, batch_size * (size_t)k * sizeof(uint64_t)));

					CUDA_CHECK_ERROR(cudaMemset(voxel_gpu->d_particle_count, 0, sizeof(uint64_t)));

					for (size_t vox_start = 0; vox_start < total_voxels; vox_start += batch_size) {
						size_t this_batch = std::min(batch_size, total_voxels - vox_start);

						double inv_l2p = bbox_size_orig * scale_space_diagonal / (double)bbox_dim;
						fill_voxel_queries_kernel << <divRoundUp(this_batch, 256ULL), 256 >> > (
							d_queries, vox_start, this_batch, dim0, dim1, off0, off1, off2,
							(float)bbox_min_orig[0], (float)bbox_min_orig[1], (float)bbox_min_orig[2],
							inv_l2p);
						CUDA_SYNC_CHECK();

						runQuery_float4_kernel << <divRoundUp(this_batch, 256ULL), 256 >> > (
							d_tree, tree_N,
							d_cand, k, std::numeric_limits<float>::infinity(),
							d_queries, (int)this_batch);
						CUDA_SYNC_CHECK();

						accumulate_voxel_sph_kernel << <divRoundUp(this_batch, 256ULL), 256 >> > (
							d_tree, d_cand, k,
							vox_start, this_batch,
							dim0, dim1, dim2,
							d_radius_particles, d_value_particles,
							particle_radius_multiplier, particle_radius_const,
							d_bbox_min_orig_dev, bbox_size_orig, bbox_dim, scale_space_diagonal,
							dense_type, dense_norm, block_name_id, use_simple_density,
							voxel_gpu->d_data_density, voxel_gpu->d_data_temp,
							voxel_gpu->d_particle_count);
						CUDA_SYNC_CHECK();
					}

					// ── Gather statistics ─────────────────────────────────────────────────
					uint64_t raw_count = 0;
					CUDA_CHECK_ERROR(cudaMemcpy(&raw_count, voxel_gpu->d_particle_count, sizeof(uint64_t), cudaMemcpyDeviceToHost));
					particles_count = (size_t)raw_count;
					voxel_gpu->find_min_max(min_value, max_value);

					CUDA_CHECK_ERROR(cudaFree(d_bbox_min_orig_dev));
					CUDA_CHECK_ERROR(cudaFree(d_queries));
					CUDA_CHECK_ERROR(cudaFree(d_cand));
				}
#endif // WITH_GPU_CUDA

#endif // WITH_CUDAKDTREE

			} // namespace kernel
		} // namespace vdb
	} // namespace common
} // namespace vdbconverter