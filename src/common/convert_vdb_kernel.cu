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
#include "sparse_common.h"
#include "dense_common.h"
#include "data_common.h"
#include <float.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <unordered_map>
#include <vector>
#include <stdexcept>

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
                float *d_offset, *d_bbox_min, *d_bbox_max;
                CUDA_CHECK_ERROR(cudaMalloc(&d_offset, 3 * sizeof(float)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_bbox_min, 3 * sizeof(float)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_bbox_max, 3 * sizeof(float)));
                
                // Initialize GPU results with extreme values
                float init_min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
                float init_max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
                CUDA_CHECK_ERROR(cudaMemcpy(d_offset, offset_position, 3 * sizeof(float), cudaMemcpyHostToDevice));
                CUDA_CHECK_ERROR(cudaMemcpy(d_bbox_min, init_min, 3 * sizeof(float), cudaMemcpyHostToDevice));
                CUDA_CHECK_ERROR(cudaMemcpy(d_bbox_max, init_max, 3 * sizeof(float), cudaMemcpyHostToDevice));
                
                // Launch kernel
                int blockSize = 256;
                int numBlocks = (num_particles + blockSize - 1) / blockSize;
                numBlocks = min(numBlocks, 1024); // Limit grid size
                
                find_bbox_kernel_cuda<<<numBlocks, blockSize>>>(
                    d_pos_particles,
                    num_particles,
                    d_offset,
                    d_bbox_min,
                    d_bbox_max
                );
                CUDA_CHECK_LAST_ERROR();
                
                // Copy results back to host
                CUDA_CHECK_ERROR(cudaMemcpy(bbox_min, d_bbox_min, 3 * sizeof(float), cudaMemcpyDeviceToHost));
                CUDA_CHECK_ERROR(cudaMemcpy(bbox_max, d_bbox_max, 3 * sizeof(float), cudaMemcpyDeviceToHost));
                
                // Free device memory
                CUDA_CHECK_ERROR(cudaFree(d_offset));
                CUDA_CHECK_ERROR(cudaFree(d_bbox_min));
                CUDA_CHECK_ERROR(cudaFree(d_bbox_max));
                
                // Wait for GPU to finish
                CUDA_SYNC_CHECK();
            }

            /**
             * @brief CUDA kernel for sparse grid conversion
             * 
             * Each thread processes particles and atomically accumulates values into voxels.
             * Uses a hash map approach with packed voxel coordinates.
             */
            __global__ void convert_to_sparse_grid_kernel_cuda(
                const float* pos_particles,
                const float* radius_particles,
                const float* value_particles,
                size_t num_particles,
                float particle_fix_size,
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

                double bbox_x_min_norm,
                double bbox_x_max_norm,
                double bbox_y_min_norm,
                double bbox_y_max_norm,
                double bbox_z_min_norm,
                double bbox_z_max_norm,
                double scale_space_diagonal,

                uint64_t* voxel_keys,      // Output: packed voxel coordinates
                float* voxel_values,       // Output: voxel values
                size_t* particle_valid,    // Output: per-particle validity flag
                float* particle_values     // Output: per-particle value for reduction
            ) {
                size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
                
                if (idx >= num_particles) return;
                                
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
                    particle_fix_size,
                    bbox_x_min_norm, bbox_x_max_norm,
                    bbox_y_min_norm, bbox_y_max_norm,
                    bbox_z_min_norm, bbox_z_max_norm,
                    filter_min, filter_max,
                    bbox_sphere_pos, bbox_sphere_r,
                    anim_type, frame_req, frame,
                    use_simple_density,
                    false, // is_dense_grid = false for sparse
                    // Outputs
                    Pos, px_norm, py_norm, pz_norm,
                    px, py, pz, v_orig
                );
                
                // Store results
                if (should_process) {
                    particle_valid[idx] = 1;
                    voxel_keys[idx] = packCoord3(px, py, pz);
                    voxel_values[idx] = v_orig;
                    particle_values[idx] = v_orig;
                } else {
                    particle_valid[idx] = 0;
                    voxel_values[idx] = 0.0f;
                    particle_values[idx] = 0.0f;
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
                const float* radius_particles,
				const float* value_particles,
                size_t num_particles,
                float particle_fix_size,
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
                size_t* particle_valid,
                float* particle_values
            ) {
                size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
                
                if (idx >= num_particles) return;               
                
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
                    particle_fix_size,
                    bbox_x_min_norm, bbox_x_max_norm,
                    bbox_y_min_norm, bbox_y_max_norm,
                    bbox_z_min_norm, bbox_z_max_norm,
                    filter_min, filter_max,
                    bbox_sphere_pos, bbox_sphere_r,
                    anim_type, frame_req, frame,
                    use_simple_density,
                    true, // is_dense_grid = true for dense
                    // Outputs
                    Pos, px_norm, py_norm, pz_norm,
                    px, py, pz, v_orig
                );
                
                if (!should_process) {
                    particle_valid[idx] = 0;
                    particle_values[idx] = 0.0f;
                    return;
                }
                
                particle_valid[idx] = 1;
                particle_values[idx] = v_orig;
                
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
                    particle_fix_size,
                    radius_particles,
                    block_name_id,
                    Pos
                );
            }
            
            /**
             * @brief GPU implementation of sparse grid conversion
             * 
             * Converts particles to sparse grid on GPU using CUDA.
             */
            void convert_to_sparse_grid_gpu(
                const float* pos_particles,
                const float* radius_particles,
				const float* value_particles,
                size_t num_particles,
                float particle_fix_size,
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

                double bbox_x_min_norm,
                double bbox_x_max_norm,
                double bbox_y_min_norm,
                double bbox_y_max_norm,
                double bbox_z_min_norm,
                double bbox_z_max_norm,
                double scale_space_diagonal,

                VoxelManager* voxel_manager,
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
                
                // Allocate device memory
                uint64_t* d_voxel_keys;
                float* d_voxel_values;
                size_t* d_particle_valid;
                float* d_particle_values;
                
                CUDA_CHECK_ERROR(cudaMalloc(&d_voxel_keys, num_particles * sizeof(uint64_t)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_voxel_values, num_particles * sizeof(float)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_particle_valid, num_particles * sizeof(size_t)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_particle_values, num_particles * sizeof(float)));
                
                // Launch kernel
                int blockSize = 256;
                int numBlocks = (num_particles + blockSize - 1) / blockSize;
                
                convert_to_sparse_grid_kernel_cuda<<<numBlocks, blockSize>>>(
                    pos_particles,
                    radius_particles,
                    value_particles,
                    num_particles,
                    particle_fix_size,
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

                    bbox_x_min_norm,
                    bbox_x_max_norm,
                    bbox_y_min_norm,
                    bbox_y_max_norm,
                    bbox_z_min_norm,
                    bbox_z_max_norm,
                    scale_space_diagonal,

                    d_voxel_keys,
                    d_voxel_values,
                    d_particle_valid,
                    d_particle_values
                );
                CUDA_CHECK_LAST_ERROR();
                
                // Copy results back to host
                std::vector<uint64_t> h_voxel_keys(num_particles);
                std::vector<float> h_voxel_values(num_particles);
                std::vector<size_t> h_particle_valid(num_particles);
                std::vector<float> h_particle_values(num_particles);
                
                CUDA_CHECK_ERROR(cudaMemcpy(h_voxel_keys.data(), d_voxel_keys, num_particles * sizeof(uint64_t), cudaMemcpyDeviceToHost));
                CUDA_CHECK_ERROR(cudaMemcpy(h_voxel_values.data(), d_voxel_values, num_particles * sizeof(float), cudaMemcpyDeviceToHost));
                CUDA_CHECK_ERROR(cudaMemcpy(h_particle_valid.data(), d_particle_valid, num_particles * sizeof(size_t), cudaMemcpyDeviceToHost));
                CUDA_CHECK_ERROR(cudaMemcpy(h_particle_values.data(), d_particle_values, num_particles * sizeof(float), cudaMemcpyDeviceToHost));
                
                // Free device memory
                CUDA_CHECK_ERROR(cudaFree(d_voxel_keys));
                CUDA_CHECK_ERROR(cudaFree(d_voxel_values));
                CUDA_CHECK_ERROR(cudaFree(d_particle_valid));
                CUDA_CHECK_ERROR(cudaFree(d_particle_values));
                
                // Process results on CPU
                float min = FLT_MAX;
                float max = -FLT_MAX;
                size_t count = 0;
                
                // Accumulate voxels into hash map (could be optimized with GPU hash map)
                std::unordered_map<uint64_t, float> voxel_map;
                for (size_t i = 0; i < num_particles; ++i) {
                    if (h_particle_valid[i]) {
                        count++;
                        float val = h_particle_values[i];
                        if (val < min) min = val;
                        if (val > max) max = val;
                        
                        voxel_map[h_voxel_keys[i]] += h_voxel_values[i];
                    }
                }
                
                particles_count = count;
                min_value = (count > 0) ? min : 0.0f;
                max_value = (count > 0) ? max : 0.0f;
                
                // TODO: Transfer voxel_map to voxel_manager
                // This requires implementing VoxelManager GPU support
            }
            
            /**
             * @brief GPU implementation of dense grid conversion
             * 
             * Converts particles to dense grid on GPU using CUDA with SPH splatting.
             */
            void convert_to_dense_grid_gpu(
                const float* pos_particles,
                const float* radius_particles,
				const float* value_particles,
                size_t num_particles,
                float particle_fix_size,
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

                double bbox_x_min_norm,
                double bbox_x_max_norm,
                double bbox_y_min_norm,
                double bbox_y_max_norm,
                double bbox_z_min_norm,
                double bbox_z_max_norm,
                double scale_space_diagonal,

                DenseParticles& grid,
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
                
                // Allocate device memory for grid and particle data
                float* d_grid_data_density;
                float* d_grid_data_temp;
                size_t* d_grid_offset;
                size_t* d_grid_dims;
                size_t* d_particle_valid;
                float* d_particle_values;
                
                size_t grid_size = grid.size();
                
                CUDA_CHECK_ERROR(cudaMalloc(&d_grid_data_density, grid_size * sizeof(float)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_grid_data_temp, grid_size * sizeof(float)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_grid_offset, 3 * sizeof(size_t)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_grid_dims, 3 * sizeof(size_t)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_particle_valid, num_particles * sizeof(size_t)));
                CUDA_CHECK_ERROR(cudaMalloc(&d_particle_values, num_particles * sizeof(float)));
                
                // Initialize grid on device
                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_data_density, grid.data_density.data(), grid_size * sizeof(float), cudaMemcpyHostToDevice));
#ifndef WITH_NO_DATA_TEMP
                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_data_temp, grid.data_temp.data(), grid_size * sizeof(float), cudaMemcpyHostToDevice));
#else
                d_grid_data_temp = nullptr;
#endif
                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_offset, grid.offset, 3 * sizeof(size_t), cudaMemcpyHostToDevice));
                CUDA_CHECK_ERROR(cudaMemcpy(d_grid_dims, grid.dims, 3 * sizeof(size_t), cudaMemcpyHostToDevice));
                
                // Launch kernel
                int blockSize = 256;
                int numBlocks = (num_particles + blockSize - 1) / blockSize;
                
                convert_to_dense_grid_kernel_cuda<<<numBlocks, blockSize>>>(
                    pos_particles,
                    radius_particles,
					value_particles,
                    num_particles,
                    particle_fix_size,
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

                    bbox_x_min_norm,
                    bbox_x_max_norm,
                    bbox_y_min_norm,
                    bbox_y_max_norm,
                    bbox_z_min_norm,
                    bbox_z_max_norm,
                    scale_space_diagonal,

                    d_grid_data_density,
                    d_grid_data_temp,
                    d_grid_offset,
                    d_grid_dims,
                    d_particle_valid,
                    d_particle_values
                );
                CUDA_CHECK_LAST_ERROR();
                
                // Copy results back to host
                CUDA_CHECK_ERROR(cudaMemcpy(grid.data_density.data(), d_grid_data_density, grid_size * sizeof(float), cudaMemcpyDeviceToHost));
#ifndef WITH_NO_DATA_TEMP
                CUDA_CHECK_ERROR(cudaMemcpy(grid.data_temp.data(), d_grid_data_temp, grid_size * sizeof(float), cudaMemcpyDeviceToHost));
#endif
                
                std::vector<size_t> h_particle_valid(num_particles);
                std::vector<float> h_particle_values(num_particles);
                CUDA_CHECK_ERROR(cudaMemcpy(h_particle_valid.data(), d_particle_valid, num_particles * sizeof(size_t), cudaMemcpyDeviceToHost));
                CUDA_CHECK_ERROR(cudaMemcpy(h_particle_values.data(), d_particle_values, num_particles * sizeof(float), cudaMemcpyDeviceToHost));
                
                // Free device memory
                CUDA_CHECK_ERROR(cudaFree(d_grid_data_density));
                if (d_grid_data_temp) CUDA_CHECK_ERROR(cudaFree(d_grid_data_temp));
                CUDA_CHECK_ERROR(cudaFree(d_grid_offset));
                CUDA_CHECK_ERROR(cudaFree(d_grid_dims));
                CUDA_CHECK_ERROR(cudaFree(d_particle_valid));
                CUDA_CHECK_ERROR(cudaFree(d_particle_values));
                
                // Compute statistics on CPU
                float min = FLT_MAX;
                float max = -FLT_MAX;
                size_t count = 0;
                
                for (size_t i = 0; i < num_particles; ++i) {
                    if (h_particle_valid[i]) {
                        count++;
                        float val = h_particle_values[i];
                        if (val < min) min = val;
                        if (val > max) max = val;
                    }
                }
                
                particles_count = count;
                min_value = (count > 0) ? min : 0.0f;
                max_value = (count > 0) ? max : 0.0f;
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
            ) {
                // Route to appropriate conversion function based on grid type
                if (grid.type == VDBParticleType::eNanoVDB || grid.type == VDBParticleType::eOpenVDB) {
                    // Sparse grid conversion using GPU hash-based voxel accumulation                    
                    convert_to_sparse_grid_gpu(
                        pos_particles,
                        radius_particles,
						value_particles,
                        num_particles,
                        particle_fix_size,
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
                        bbox_x_min_norm,
                        bbox_x_max_norm,
                        bbox_y_min_norm,
                        bbox_y_max_norm,
                        bbox_z_min_norm,
                        bbox_z_max_norm,
                        scale_space_diagonal,
                        grid.sparse_particles.get(),
                        min_value,
                        max_value,
                        particles_count
                    );
                }
                else if (grid.type == VDBParticleType::eDense) {
                    // Dense grid conversion using GPU SPH kernel splatting                    
                    convert_to_dense_grid_gpu(
                        pos_particles,
                        radius_particles,
						value_particles,
                        num_particles,
                        particle_fix_size,
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

                        bbox_x_min_norm,
                        bbox_x_max_norm,
                        bbox_y_min_norm,
                        bbox_y_max_norm,
                        bbox_z_min_norm,
                        bbox_z_max_norm,
                        scale_space_diagonal,

                        grid.dense_grid,
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
            
        } // namespace kernel
    } // namespace vdb
} // namespace common