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
 * Source: SPHtoGrid: https://github.com/LudwigBoess/SPHtoGrid.jl, SPHKernels: https://github.com/LudwigBoess/SPHKernels.jl
 */

#pragma once

#include <iostream>

/**
 * @brief Check CUDA error and print detailed message if error occurred
 * 
 * This macro wraps CUDA API calls and checks their return value.
 * If an error occurs, it prints the error message with file and line information.
 * 
 * Usage: CUDA_CHECK_ERROR(cudaMalloc(&d_ptr, size));
 */
#define CUDA_CHECK_ERROR(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error in " << __FILE__ << " at line " << __LINE__ << ": " \
                      << cudaGetErrorString(err) << " (error code: " << err << ")" << std::endl; \
        } \
    } while (0)

/**
 * @brief Check for CUDA errors after kernel launch or async operations
 * 
 * This macro checks cudaGetLastError() to catch errors from kernel launches
 * or other asynchronous CUDA operations.
 * 
 * Usage: CUDA_CHECK_LAST_ERROR();
 */
#define CUDA_CHECK_LAST_ERROR() \
    do { \
        cudaError_t err = cudaGetLastError(); \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error (last error) in " << __FILE__ << " at line " << __LINE__ << ": " \
                      << cudaGetErrorString(err) << " (error code: " << err << ")" << std::endl; \
        } \
    } while (0)

/**
 * @brief Synchronize device and check for errors
 * 
 * This macro calls cudaDeviceSynchronize() and checks for any errors.
 * Useful for ensuring all previous CUDA operations completed successfully.
 * 
 * Usage: CUDA_SYNC_CHECK();
 */
#define CUDA_SYNC_CHECK() \
    do { \
        cudaError_t err = cudaDeviceSynchronize(); \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Sync Error in " << __FILE__ << " at line " << __LINE__ << ": " \
                      << cudaGetErrorString(err) << " (error code: " << err << ")" << std::endl; \
        } \
    } while (0)

/**
 * @brief Debug macro to copy GPU array to host and print first 10 values
 * @details Useful for logging/debugging GPU computations
 * @param gpu_ptr Pointer to GPU memory
 * @param size Total size of the array
 * @param label Descriptive label for the log output
 */
#ifdef WITH_GPU_CUDA
#define DEBUG_PRINT_GPU_ARRAY(gpu_ptr, size, label, type) \
	{ \
		size_t debug_count = std::min((size_t)10, (size_t)size); \
		std::vector<type> host_debug_data(size); \
		CUDA_CHECK_ERROR(cudaMemcpy(host_debug_data.data(), gpu_ptr, sizeof(type) * size, cudaMemcpyDeviceToHost)); \
		printf("DEBUG [%s], size: %zu, values: ", label, size); \
        size_t c = 0; \
		for (size_t i = 0; i < size; i++) { \
            if (host_debug_data[i] == 0.0f) continue; \
			printf("%zu: %.6f ", i, (float)host_debug_data[i]); \
            if (++c >= debug_count) break; \
		} \
		printf("\n"); \
	}
#else
#define DEBUG_PRINT_GPU_ARRAY(gpu_ptr, size, label)
#endif    

#ifdef WITH_CUDA_MALLOCMANAGED
#   define CUDA_MALLOC cudaMallocManaged
#else
#   define CUDA_MALLOC cudaMalloc
#endif 
