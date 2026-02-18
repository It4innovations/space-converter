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
#include <float.h>
#include <cuda_runtime.h>

namespace common {
    namespace vdb {
        namespace kernel {
            
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
            
        } // namespace kernel
    } // namespace vdb
} // namespace common