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

namespace common {
    namespace vdb {
        namespace kernel {

#ifdef WITH_GPU_CUDA
            
            /**
             * @brief GPU implementation of bounding box calculation
             * 
             * Launches CUDA kernel to compute min/max coordinates across all particles.
             * Implemented in convert_vdb.cu
             */
            void find_bbox_gpu(
                const float* d_pos_particles,
                size_t num_particles,
                const float* offset_position,
                float* bbox_min,
                float* bbox_max
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
            
        } // namespace kernel
    } // namespace vdb
} // namespace common