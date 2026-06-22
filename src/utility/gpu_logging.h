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

#ifdef WITH_GPU_CUDA
#include <cuda_runtime.h>
#endif

#include <stdio.h>

//#define SPACE_CONVERTER_GPU_LOGGING_ENABLED

namespace space_converter {
    namespace gpu_logging {
        // Global rank for GPU logging (can be set from MPI)
        extern int g_gpu_log_rank;

        /**
         * @brief Set the MPI rank for GPU logging output
         *
         * @param rank The MPI rank of the current process
         */
        inline void set_gpu_log_rank(int rank) {
            g_gpu_log_rank = rank;
        }

        /**
         * @brief Get the current GPU logging rank
         *
         * @return Current MPI rank
         */
        inline int get_gpu_log_rank() {
            return g_gpu_log_rank;
        }
	} // namespace gpu_logging
} // namespace space_converter

#ifdef SPACE_CONVERTER_GPU_LOGGING_ENABLED
/**
 * @brief Macro for timing CUDA kernel launches
 * 
 * Usage:
 *   GPU_KERNEL_TIME_START("kernel_name");
 *   my_kernel<<<grid, block>>>(...);
 *   GPU_KERNEL_TIME_END("kernel_name");
 * 
 * Prints: "[Rank X] GPU Kernel 'kernel_name': Y.YYY ms"
 */
#define GPU_KERNEL_TIME_START(kernel_name) \
    cudaEvent_t gpu_start_##kernel_name, gpu_stop_##kernel_name; \
    cudaEventCreate(&gpu_start_##kernel_name); \
    cudaEventCreate(&gpu_stop_##kernel_name); \
    cudaEventRecord(gpu_start_##kernel_name, 0);

#define GPU_KERNEL_TIME_END(kernel_name) \
    cudaEventRecord(gpu_stop_##kernel_name, 0); \
    cudaEventSynchronize(gpu_stop_##kernel_name); \
    float gpu_elapsed_##kernel_name = 0.0f; \
    cudaEventElapsedTime(&gpu_elapsed_##kernel_name, gpu_start_##kernel_name, gpu_stop_##kernel_name); \
    printf("[rank #%d]: GPU Kernel '%s': %.3f ms\n", gpu_logging::g_gpu_log_rank, #kernel_name, gpu_elapsed_##kernel_name); \
    cudaEventDestroy(gpu_start_##kernel_name); \
    cudaEventDestroy(gpu_stop_##kernel_name);

/**
 * @brief Macro for timing CUB operations (like DeviceRadixSort, DeviceReduce)
 * 
 * Usage:
 *   GPU_CUB_TIME_START("sort_operation");
 *   cub::DeviceRadixSort::SortPairs(...);
 *   GPU_CUB_TIME_END("sort_operation");
 * 
 * Prints: "[Rank X] GPU CUB 'sort_operation': Y.YYY ms"
 */
#define GPU_CUB_TIME_START(operation_name) \
    cudaEvent_t gpu_cub_start_##operation_name, gpu_cub_stop_##operation_name; \
    cudaEventCreate(&gpu_cub_start_##operation_name); \
    cudaEventCreate(&gpu_cub_stop_##operation_name); \
    cudaEventRecord(gpu_cub_start_##operation_name, 0);

#define GPU_CUB_TIME_END(operation_name) \
    cudaEventRecord(gpu_cub_stop_##operation_name, 0); \
    cudaEventSynchronize(gpu_cub_stop_##operation_name); \
    float gpu_cub_elapsed_##operation_name = 0.0f; \
    cudaEventElapsedTime(&gpu_cub_elapsed_##operation_name, gpu_cub_start_##operation_name, gpu_cub_stop_##operation_name); \
    printf("[rank #%d]: GPU CUB '%s': %.3f ms\n", gpu_logging::g_gpu_log_rank, #operation_name, gpu_cub_elapsed_##operation_name); \
    cudaEventDestroy(gpu_cub_start_##operation_name); \
    cudaEventDestroy(gpu_cub_stop_##operation_name);

#else
#define GPU_KERNEL_TIME_START(kernel_name)
#define GPU_KERNEL_TIME_END(kernel_name)
#define GPU_CUB_TIME_START(operation_name)
#define GPU_CUB_TIME_END(operation_name)
#endif
