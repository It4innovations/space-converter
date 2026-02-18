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

 #include "data_cache.h"
 #include "../utility/gpu_utility.h"
 #include <cuda_runtime.h>
 #include <thrust/device_vector.h>
 #include <thrust/sort.h>
 #include <thrust/copy.h>
 #include <thrust/execution_policy.h>
 #include <iostream>

namespace common {
    namespace cache {     
        
        void CacheManager::copy_particles_to_gpu() {
            // Free existing GPU memory if any
            free_gpu_memory();

            // Allocate and copy pos_particles_per_ptype
            d_pos_particles_per_ptype.resize(pos_particles_per_ptype.size());
            for (size_t i = 0; i < pos_particles_per_ptype.size(); ++i) {
                size_t size_bytes = pos_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(cudaMalloc(&d_pos_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_pos_particles_per_ptype[i], pos_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_pos_particles_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy particles_id_ordered_per_ptype
            d_particles_id_ordered_per_ptype.resize(particles_id_ordered_per_ptype.size());
            for (size_t i = 0; i < particles_id_ordered_per_ptype.size(); ++i) {
                size_t size_bytes = particles_id_ordered_per_ptype[i].size() * sizeof(size_t);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(cudaMalloc(&d_particles_id_ordered_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_particles_id_ordered_per_ptype[i], particles_id_ordered_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_particles_id_ordered_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy radius_particles_per_ptype
            d_radius_particles_per_ptype.resize(radius_particles_per_ptype.size());
            for (size_t i = 0; i < radius_particles_per_ptype.size(); ++i) {
                size_t size_bytes = radius_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(cudaMalloc(&d_radius_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_radius_particles_per_ptype[i], radius_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_radius_particles_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy particles_ptype_offset
            d_particles_ptype_offset = nullptr;
            if (!particles_ptype_offset.empty()) {
                size_t size_bytes = particles_ptype_offset.size() * sizeof(size_t);
                CUDA_CHECK_ERROR(cudaMalloc(&d_particles_ptype_offset, size_bytes));
                CUDA_CHECK_ERROR(cudaMemcpy(d_particles_ptype_offset, particles_ptype_offset.data(), 
                         size_bytes, cudaMemcpyHostToDevice));
            }

            // Allocate and copy rho_particles_per_ptype
            d_rho_particles_per_ptype.resize(rho_particles_per_ptype.size());
            for (size_t i = 0; i < rho_particles_per_ptype.size(); ++i) {
                size_t size_bytes = rho_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(cudaMalloc(&d_rho_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_rho_particles_per_ptype[i], rho_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_rho_particles_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy mass_particles_per_ptype
            d_mass_particles_per_ptype.resize(mass_particles_per_ptype.size());
            for (size_t i = 0; i < mass_particles_per_ptype.size(); ++i) {
                size_t size_bytes = mass_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(cudaMalloc(&d_mass_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_mass_particles_per_ptype[i], mass_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_mass_particles_per_ptype[i] = nullptr;
                }
            }

            // Check for CUDA errors
            CUDA_CHECK_LAST_ERROR();
        }

        void CacheManager::free_gpu_memory() {
            // Free pos_particles_per_ptype
            for (auto ptr : d_pos_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_pos_particles_per_ptype.clear();

            // Free particles_id_ordered_per_ptype
            for (auto ptr : d_particles_id_ordered_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_particles_id_ordered_per_ptype.clear();

            // Free radius_particles_per_ptype
            for (auto ptr : d_radius_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_radius_particles_per_ptype.clear();

            // Free particles_ptype_offset
            if (d_particles_ptype_offset != nullptr) {
                CUDA_CHECK_ERROR(cudaFree(d_particles_ptype_offset));
                d_particles_ptype_offset = nullptr;
            }

            // Free rho_particles_per_ptype
            for (auto ptr : d_rho_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_rho_particles_per_ptype.clear();

            // Free mass_particles_per_ptype
            for (auto ptr : d_mass_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_mass_particles_per_ptype.clear();

            // Check for CUDA errors
            CUDA_CHECK_LAST_ERROR();
        }

        void CacheManager::sort_particles_by_radius_gpu() {
            // Sort particle IDs by radius on GPU for each particle type using Thrust
            for (size_t ptype = 0; ptype < radius_particles_per_ptype.size(); ++ptype) {
                if (ptype >= particles_id_ordered_per_ptype.size()) {
                    continue;
                }

                auto& radii = radius_particles_per_ptype[ptype];
                auto& ids = particles_id_ordered_per_ptype[ptype];

                if (radii.size() != ids.size() || radii.empty()) {
                    continue;
                }

                size_t n = radii.size();

                // Create device vectors from host data
                thrust::device_vector<float> d_radii(radii.begin(), radii.end());
                thrust::device_vector<size_t> d_ids(ids.begin(), ids.end());

                // Sort by key: radii are keys, ids are values
                // This sorts both arrays so that radii is in ascending order
                // and ids are reordered accordingly
                thrust::sort_by_key(d_radii.begin(), d_radii.end(), d_ids.begin());

                // Copy sorted data back to host
                thrust::copy(d_radii.begin(), d_radii.end(), radii.begin());
                thrust::copy(d_ids.begin(), d_ids.end(), ids.begin());
            }

            // Check for CUDA errors
            CUDA_CHECK_LAST_ERROR();
        }

        void CacheManager::sort_particles_by_radius_gpu_inplace() {
            // Sort particle IDs by radius directly on GPU using device pointers
            // No CPU<->GPU transfers, operates on d_radius_particles_per_ptype and d_particles_id_ordered_per_ptype
            
            for (size_t ptype = 0; ptype < d_radius_particles_per_ptype.size(); ++ptype) {
                if (ptype >= d_particles_id_ordered_per_ptype.size()) {
                    continue;
                }

                float* d_radii_ptr = d_radius_particles_per_ptype[ptype];
                size_t* d_ids_ptr = d_particles_id_ordered_per_ptype[ptype];

                if (d_radii_ptr == nullptr || d_ids_ptr == nullptr) {
                    continue;
                }

                // Get the size from the CPU vectors
                if (ptype >= radius_particles_per_ptype.size() || 
                    ptype >= particles_id_ordered_per_ptype.size()) {
                    continue;
                }

                size_t n = radius_particles_per_ptype[ptype].size();
                if (n != particles_id_ordered_per_ptype[ptype].size() || n == 0) {
                    continue;
                }

                // Wrap raw device pointers with thrust::device_ptr
                thrust::device_ptr<float> d_radii_thrust(d_radii_ptr);
                thrust::device_ptr<size_t> d_ids_thrust(d_ids_ptr);

                // Sort by key: radii are keys, ids are values
                // This sorts both arrays in-place on the GPU
                thrust::sort_by_key(d_radii_thrust, d_radii_thrust + n, d_ids_thrust);
            }

            // Check for CUDA errors
            CUDA_CHECK_LAST_ERROR();
        }
    }
}