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
 * GPU Device Compatibility Layer for CUDA and HIP (AMD), modeled after
 * https://github.com/It4innovations/braas-hpc-gpujpeg/blob/braas-hpc-dev/src/gpujpeg_device_compat.h
 *
 * Include this header instead of <cuda_runtime.h> / <cub/cub.cuh>.
 *
 * Backend selection:
 *  - WITH_GPU_HIP  (CMake) -> HIP/ROCm backend (AMD GPUs, e.g. gfx90a on LUMI-G).
 *  - WITH_GPU_CUDA (CMake) -> CUDA backend.
 *  - neither               -> CPU-only build, this header is a no-op.
 *
 * Note: a HIP build also defines the legacy WITH_GPU_CUDA guard macro (see
 * src/CMakeLists.txt) so that the existing "#ifdef WITH_GPU_CUDA" GPU code
 * paths stay enabled; the CUDA runtime API used by them (cudaMalloc,
 * cudaMemcpy, cudaEvent_*, ...) is remapped to HIP below.
 *
 * CUB equivalent on AMD: hipCUB (backed by rocPRIM) provides the same API
 * (hipcub::DeviceRadixSort::SortPairs, hipcub::DeviceReduce::{Min,Max,Reduce,
 * ReduceByKey}, ...); "namespace cub = hipcub" makes existing cub:: calls
 * compile unchanged. Thrust code needs no change: rocThrust keeps the
 * <thrust/...> headers and the thrust:: namespace.
 */

#ifndef SPACE_CONVERTER_GPU_DEVICE_COMPAT_H
#define SPACE_CONVERTER_GPU_DEVICE_COMPAT_H

#include <stddef.h>  // for size_t

// ==============================================================================
// Backend selection (WITH_GPU_HIP must win: HIP builds also define WITH_GPU_CUDA)
// ==============================================================================
#if defined(WITH_GPU_HIP) || defined(__HIPCC__) || defined(__HIP_PLATFORM_AMD__)
#    define SPACE_GPU_USE_HIP 1
#elif defined(WITH_GPU_CUDA) || defined(__CUDACC__)
#    define SPACE_GPU_USE_CUDA 1
#endif

// ==============================================================================
// CUDA Backend
// ==============================================================================
#if defined(SPACE_GPU_USE_CUDA)

#include <cuda_runtime.h>

#ifdef __CUDACC__
#include <device_launch_parameters.h>
#include <cub/cub.cuh>
#endif

// Warp size
#define GPU_WARP_SIZE                       32

// Unified gpu* aliases (for new code; existing code may keep the cuda* names)
#define gpuError_t                          cudaError_t
#define gpuSuccess                          cudaSuccess
#define gpuGetErrorString                   cudaGetErrorString
#define gpuGetLastError                     cudaGetLastError
#define gpuDeviceSynchronize                cudaDeviceSynchronize

#define gpuMalloc                           cudaMalloc
#define gpuMallocManaged                    cudaMallocManaged
#define gpuFree                             cudaFree
#define gpuMemcpy                           cudaMemcpy
#define gpuMemcpyAsync                      cudaMemcpyAsync
#define gpuMemset                           cudaMemset
#define gpuMemsetAsync                      cudaMemsetAsync

#define gpuMemcpyKind                       cudaMemcpyKind
#define gpuMemcpyHostToDevice               cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost               cudaMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice             cudaMemcpyDeviceToDevice
#define gpuMemcpyDefault                    cudaMemcpyDefault

#define gpuStream_t                         cudaStream_t
#define gpuStreamSynchronize                cudaStreamSynchronize

#define gpuEvent_t                          cudaEvent_t
#define gpuEventCreate                      cudaEventCreate
#define gpuEventDestroy                     cudaEventDestroy
#define gpuEventRecord                      cudaEventRecord
#define gpuEventSynchronize                 cudaEventSynchronize
#define gpuEventElapsedTime                 cudaEventElapsedTime

#define gpuSetDevice                        cudaSetDevice
#define gpuGetDevice                        cudaGetDevice
#define gpuGetDeviceCount                   cudaGetDeviceCount

// ==============================================================================
// HIP Backend (AMD)
// ==============================================================================
#elif defined(SPACE_GPU_USE_HIP)

// Define HIP platform for AMD GPUs
#ifndef __HIP_PLATFORM_AMD__
#define __HIP_PLATFORM_AMD__
#endif

#ifdef __HIPCC__
// HIP device-capable compiler (hipcc / amdclang++): full runtime + kernel language
#include <hip/hip_runtime.h>
// hipCUB: CUB-compatible primitives for AMD, backed by rocPRIM
#include <hipcub/hipcub.hpp>
namespace cub = hipcub;
#else
// Plain host compiler (g++ / CC): host-side HIP runtime API only
#include <hip/hip_runtime_api.h>
#endif

// Warp (wavefront) size: 64 on AMD CDNA (gfx90a); use the warpSize builtin in
// device code when the exact value matters.
#define GPU_WARP_SIZE                       64

#ifdef __HIPCC__
// CUDA warp-sync collectives -> HIP collectives (mask argument is dropped;
// the XOR-butterfly reductions used in this project are wave64-safe)
#define __shfl_xor_sync(mask, ...)          __shfl_xor(__VA_ARGS__)
#define __shfl_sync(mask, ...)              __shfl(__VA_ARGS__)
#define __shfl_down_sync(mask, ...)         __shfl_down(__VA_ARGS__)
#define __shfl_up_sync(mask, ...)           __shfl_up(__VA_ARGS__)
#endif

// ------------------------------------------------------------------------------
// Legacy CUDA runtime API -> HIP mapping, so the existing .cu sources compile
// unchanged (hipify-style macro table, restricted to the APIs used here).
// ------------------------------------------------------------------------------
#define cudaError_t                         hipError_t
#define cudaSuccess                         hipSuccess
#define cudaGetErrorString                  hipGetErrorString
#define cudaGetLastError                    hipGetLastError
#define cudaDeviceSynchronize               hipDeviceSynchronize

#define cudaMalloc                          hipMalloc
#define cudaMallocManaged                   hipMallocManaged
#define cudaFree                            hipFree
#define cudaMemcpy                          hipMemcpy
#define cudaMemcpyAsync                     hipMemcpyAsync
#define cudaMemset                          hipMemset
#define cudaMemsetAsync                     hipMemsetAsync

#define cudaMemcpyKind                      hipMemcpyKind
#define cudaMemcpyHostToDevice              hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost              hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice            hipMemcpyDeviceToDevice
#define cudaMemcpyDefault                   hipMemcpyDefault

#define cudaStream_t                        hipStream_t
#define cudaStreamSynchronize               hipStreamSynchronize

#define cudaEvent_t                         hipEvent_t
#define cudaEventCreate                     hipEventCreate
#define cudaEventDestroy                    hipEventDestroy
#define cudaEventRecord                     hipEventRecord
#define cudaEventSynchronize                hipEventSynchronize
#define cudaEventElapsedTime                hipEventElapsedTime

#define cudaSetDevice                       hipSetDevice
#define cudaGetDevice                       hipGetDevice
#define cudaGetDeviceCount                  hipGetDeviceCount
#define cudaGetDeviceProperties             hipGetDeviceProperties
#define cudaDeviceProp                      hipDeviceProp_t

// Unified gpu* aliases (for new code)
#define gpuError_t                          hipError_t
#define gpuSuccess                          hipSuccess
#define gpuGetErrorString                   hipGetErrorString
#define gpuGetLastError                     hipGetLastError
#define gpuDeviceSynchronize                hipDeviceSynchronize

#define gpuMalloc                           hipMalloc
#define gpuMallocManaged                    hipMallocManaged
#define gpuFree                             hipFree
#define gpuMemcpy                           hipMemcpy
#define gpuMemcpyAsync                      hipMemcpyAsync
#define gpuMemset                           hipMemset
#define gpuMemsetAsync                      hipMemsetAsync

#define gpuMemcpyKind                       hipMemcpyKind
#define gpuMemcpyHostToDevice               hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost               hipMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice             hipMemcpyDeviceToDevice
#define gpuMemcpyDefault                    hipMemcpyDefault

#define gpuStream_t                         hipStream_t
#define gpuStreamSynchronize                hipStreamSynchronize

#define gpuEvent_t                          hipEvent_t
#define gpuEventCreate                      hipEventCreate
#define gpuEventDestroy                     hipEventDestroy
#define gpuEventRecord                      hipEventRecord
#define gpuEventSynchronize                 hipEventSynchronize
#define gpuEventElapsedTime                 hipEventElapsedTime

#define gpuSetDevice                        hipSetDevice
#define gpuGetDevice                        hipGetDevice
#define gpuGetDeviceCount                   hipGetDeviceCount

#endif // SPACE_GPU_USE_HIP / SPACE_GPU_USE_CUDA

#endif // SPACE_CONVERTER_GPU_DEVICE_COMPAT_H
