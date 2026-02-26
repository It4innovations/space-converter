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

#include "sparse_common.h"
#include <cuda_runtime.h>
#include <cub/cub.cuh>

#include "data_common.h"
#include "../utility/gpu_utility.h"

namespace common {
	namespace vdb {
		namespace sparse {

			// ---------------------------------------------
			// Kernel: Voxels -> (keys, values)
			// ---------------------------------------------
			__global__ void voxelsToKeyValue(const Voxel* __restrict__ in,
				int n,
				uint64_t* __restrict__ keys,
				float* __restrict__ vals)
			{
				int tid = blockIdx.x * blockDim.x + threadIdx.x;
				if (tid >= n) return;

				Voxel v = in[tid];
				keys[tid] = packCoord3(v.i, v.j, v.k);
				vals[tid] = v.value;
			}

			// ---------------------------------------------
			// Kernel: Pixels -> (keys, values)
			// ---------------------------------------------			
			__global__ void pixelsToKeyValue(const Voxel* __restrict__ in,
				int n,
				uint64_t* __restrict__ keys,
				float* __restrict__ vals)
			{
				int tid = blockIdx.x * blockDim.x + threadIdx.x;
				if (tid >= n) return;

				Voxel v = in[tid];
				keys[tid] = packCoord3(v.i, v.j, v.k);
				vals[tid] = v.value;
			}

			// ---------------------------------------------
			// Kernel: (keys, values) -> Voxels (unique output)
			// ---------------------------------------------
			__global__ void keyValueToVoxels(const uint64_t* __restrict__ keys,
				const float* __restrict__ vals,
				int n,
				Voxel* __restrict__ out)
			{
				int tid = blockIdx.x * blockDim.x + threadIdx.x;
				if (tid >= n) return;

				int x, y, z;
				unpackCoord3(keys[tid], x, y, z);
				out[tid] = Voxel(x, y, z, vals[tid]);
			}

			// CustomSum functor
			struct CustomSum
			{
				template <typename T>
				__device__ __forceinline__
				T operator()(const T &a, const T &b) const {
					return a + b;
				}
			};

			// ---------------------------------------------
			// Fast manager: per-batch accumulate via sort+reduce
			// ---------------------------------------------
			VoxelGPUManagerSortReduce::VoxelGPUManagerSortReduce(): 
				m_max(0), m_last_count(0),
				d_keys(nullptr), d_vals(nullptr),
				d_keys_alt(nullptr), d_vals_alt(nullptr),
				d_keys_out(nullptr), d_vals_out(nullptr),
				d_num_out(nullptr),
				d_sort_temp(nullptr), m_sort_temp_bytes(0),
				d_reduce_temp(nullptr), m_reduce_temp_bytes(0),
				d_minmax_temp(nullptr), m_minmax_temp_bytes(0),
				d_min_out(nullptr), d_max_out(nullptr),
				d_particle_count(nullptr),
				d_bbox_min_orig(nullptr), d_offset_position(nullptr), d_bbox_sphere_pos(nullptr),
				common::vdb::VoxelSparseManager()
			{

			}

			void VoxelGPUManagerSortReduce::init(unsigned int expected_voxels)
			{
				    m_max = expected_voxels;
					//cudaMalloc(&d_inVoxels, m_max * sizeof(Voxel));

					CUDA_CHECK_ERROR(cudaMalloc(&d_keys, m_max * sizeof(uint64_t)));
					CUDA_CHECK_ERROR(cudaMalloc(&d_vals, m_max * sizeof(float)));

					CUDA_CHECK_ERROR(cudaMalloc(&d_keys_alt, m_max * sizeof(uint64_t)));
					CUDA_CHECK_ERROR(cudaMalloc(&d_vals_alt, m_max * sizeof(float)));

					// Worst-case: all unique => output size == input size
					CUDA_CHECK_ERROR(cudaMalloc(&d_keys_out, m_max * sizeof(uint64_t)));
					CUDA_CHECK_ERROR(cudaMalloc(&d_vals_out, m_max * sizeof(float)));
					CUDA_CHECK_ERROR(cudaMalloc(&d_num_out, sizeof(int)));

					// Precompute temp storage sizes for CUB primitives (max)
					size_t sort_bytes = 0;
					cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes,
						d_keys, d_keys_alt,
						d_vals, d_vals_alt,
						(int)m_max);
					m_sort_temp_bytes = sort_bytes;
					CUDA_CHECK_ERROR(cudaMalloc(&d_sort_temp, m_sort_temp_bytes));

					size_t reduce_bytes = 0;
					CustomSum op_sum;
					cub::DeviceReduce::ReduceByKey(nullptr, reduce_bytes,
						d_keys_alt, d_keys_out,
						d_vals_alt, d_vals_out,
						d_num_out,
						op_sum,
						(int)m_max);
				m_reduce_temp_bytes = reduce_bytes;
				CUDA_CHECK_ERROR(cudaMalloc(&d_reduce_temp, m_reduce_temp_bytes));

				// Allocate min/max output buffers
				CUDA_CHECK_ERROR(cudaMalloc(&d_min_out, sizeof(float)));
				CUDA_CHECK_ERROR(cudaMalloc(&d_max_out, sizeof(float)));

				// Allocate and zero the processed-particle counter
				CUDA_CHECK_ERROR(cudaMalloc(&d_particle_count, sizeof(uint64_t)));
				CUDA_CHECK_ERROR(cudaMemset(d_particle_count, 0, sizeof(uint64_t)));

				// Persistent device buffers for per-call host parameters
				CUDA_CHECK_ERROR(cudaMalloc(&d_bbox_min_orig,   3 * sizeof(int)));
				CUDA_CHECK_ERROR(cudaMalloc(&d_offset_position, 3 * sizeof(float)));
				CUDA_CHECK_ERROR(cudaMalloc(&d_bbox_sphere_pos, 3 * sizeof(float)));

				// Precompute temp storage for min/max reductions
				size_t min_bytes = 0;
				cub::DeviceReduce::Min(nullptr, min_bytes, d_vals_out, d_min_out, (int)m_max);
				
				size_t max_bytes = 0;
				cub::DeviceReduce::Max(nullptr, max_bytes, d_vals_out, d_max_out, (int)m_max);
				
				// Use the larger of the two
				m_minmax_temp_bytes = (min_bytes > max_bytes) ? min_bytes : max_bytes;
				CUDA_CHECK_ERROR(cudaMalloc(&d_minmax_temp, m_minmax_temp_bytes));
			}

			VoxelGPUManagerSortReduce::~VoxelGPUManagerSortReduce()
			{
					// cudaFree(d_inVoxels);

					if (d_keys) {
						CUDA_CHECK_ERROR(cudaFree(d_keys));
					}

					if (d_vals) {
						CUDA_CHECK_ERROR(cudaFree(d_vals));
					}					

					if (d_keys_alt) {
						CUDA_CHECK_ERROR(cudaFree(d_keys_alt));
					}
					if (d_vals_alt) {
						CUDA_CHECK_ERROR(cudaFree(d_vals_alt));
					}

					if (d_keys_out) {
						CUDA_CHECK_ERROR(cudaFree(d_keys_out));
					}
					if (d_vals_out) {
						CUDA_CHECK_ERROR(cudaFree(d_vals_out));
					}

					if (d_num_out) {
						CUDA_CHECK_ERROR(cudaFree(d_num_out));
					}

					if (d_sort_temp) {
						CUDA_CHECK_ERROR(cudaFree(d_sort_temp));
					}
					if (d_reduce_temp) {
						CUDA_CHECK_ERROR(cudaFree(d_reduce_temp));
					}
					if (d_minmax_temp) {
						CUDA_CHECK_ERROR(cudaFree(d_minmax_temp));
					}
					if (d_min_out) {
						CUDA_CHECK_ERROR(cudaFree(d_min_out));
					}
					if (d_max_out) {
						CUDA_CHECK_ERROR(cudaFree(d_max_out));
					}
					if (d_particle_count) {
						CUDA_CHECK_ERROR(cudaFree(d_particle_count));
					}
					if (d_bbox_min_orig) {
						CUDA_CHECK_ERROR(cudaFree(d_bbox_min_orig));
					}
					if (d_offset_position) {
						CUDA_CHECK_ERROR(cudaFree(d_offset_position));
					}
					if (d_bbox_sphere_pos) {
						CUDA_CHECK_ERROR(cudaFree(d_bbox_sphere_pos));
					}
			}

			// Accumulate a batch of voxels: output becomes unique (i,j,k) with summed values
			// Returns number of unique voxels in the batch.
			int VoxelGPUManagerSortReduce::insertOrUpdate(const Voxel* h_voxels, int num_voxels)
			{
					if (num_voxels <= 0) {
						m_last_count = 0;
						return 0;
					}
					if ((size_t)num_voxels > m_max) {
						printf("[VoxelGPUManagerSortReduce] ERROR: num_voxels=%d exceeds max=%zu\n",
							num_voxels, m_max);
						return 0;
					}

					// H2D
					Voxel* d_inVoxels = nullptr;
					CUDA_CHECK_ERROR(cudaMalloc(&d_inVoxels, m_max * sizeof(Voxel)));
					CUDA_CHECK_ERROR(cudaMemcpy(d_inVoxels, h_voxels, num_voxels * sizeof(Voxel), cudaMemcpyHostToDevice));

					// map -> key/value
					{
						int block = 256;
						int grid = (num_voxels + block - 1) / block;
						voxelsToKeyValue << <grid, block >> > (d_inVoxels, num_voxels, d_keys, d_vals);
					}

					CUDA_CHECK_ERROR(cudaFree(d_inVoxels));

					// sort by key (pairs)
					cub::DeviceRadixSort::SortPairs(d_sort_temp, m_sort_temp_bytes,
					d_keys, d_keys_alt,
						d_vals, d_vals_alt,
						num_voxels);

					// reduce-by-key (sum values for identical keys)
					CustomSum op_sum;
					cub::DeviceReduce::ReduceByKey(d_reduce_temp, m_reduce_temp_bytes,
						d_keys_alt, d_keys_out,
						d_vals_alt, d_vals_out,
						d_num_out,
						op_sum,
						num_voxels);

					int h_num_out = 0;
					CUDA_CHECK_ERROR(cudaMemcpy(&h_num_out, d_num_out, sizeof(int), cudaMemcpyDeviceToHost));

				m_last_count = h_num_out;
				return h_num_out;
			}

			int VoxelGPUManagerSortReduce::update(size_t count)
			{
				// sort by key (pairs)
				cub::DeviceRadixSort::SortPairs(d_sort_temp, m_sort_temp_bytes,
					d_keys, d_keys_alt,
					d_vals, d_vals_alt,
					count);

				// reduce-by-key (sum values for identical keys)
				CustomSum op_sum;
				cub::DeviceReduce::ReduceByKey(d_reduce_temp, m_reduce_temp_bytes,
					d_keys_alt, d_keys_out,
					d_vals_alt, d_vals_out,
					d_num_out,
					op_sum,
					count);

				int h_num_out = 0;
				CUDA_CHECK_ERROR(cudaMemcpy(&h_num_out, d_num_out, sizeof(int), cudaMemcpyDeviceToHost));

				m_last_count = h_num_out;
				return h_num_out;
			}

			// Extract the last accumulated unique voxels back to host
			int VoxelGPUManagerSortReduce::extractAll(Voxel** h_output_voxels)
			{
					const int n = m_last_count;
					if (n <= 0) {
						*h_output_voxels = nullptr;
						return 0;
					}

					Voxel* d_out_voxels = nullptr;
					CUDA_CHECK_ERROR(cudaMalloc(&d_out_voxels, n * sizeof(Voxel)));

					{
						int block = 256;
						int grid = (n + block - 1) / block;
						keyValueToVoxels << <grid, block >> > (d_keys_out, d_vals_out, n, d_out_voxels);
					}

					*h_output_voxels = new Voxel[n];
					CUDA_CHECK_ERROR(cudaMemcpy(*h_output_voxels, d_out_voxels, n * sizeof(Voxel), cudaMemcpyDeviceToHost));

					CUDA_CHECK_ERROR(cudaDeviceSynchronize());

					CUDA_CHECK_ERROR(cudaFree(d_out_voxels));

				return n;
			}

			// ================================================================================
			// CPU-side methods
			// ================================================================================

			// CPU-side serialization: copy from device to host, then serialize
			void VoxelGPUManagerSortReduce::serializeCPU(uint8_t* bin_data) {
				// int    m_last_count = 0;
				// uint64_t* d_keys_out = nullptr;
				// float* d_vals_out = nullptr;

				// Allocate host memory for keys and values
				uint64_t* h_keys = new uint64_t[m_last_count];
				float* h_vals = new float[m_last_count];

				// Copy directly from device
				CUDA_CHECK_ERROR(cudaMemcpy(h_keys, d_keys_out, m_last_count * sizeof(uint64_t), cudaMemcpyDeviceToHost));
				CUDA_CHECK_ERROR(cudaMemcpy(h_vals, d_vals_out, m_last_count * sizeof(float), cudaMemcpyDeviceToHost));

				uint8_t* ptr = bin_data;				
				memcpy(ptr, &m_last_count, sizeof(int));
				ptr += sizeof(int);
				memcpy(ptr, h_keys, m_last_count * sizeof(uint64_t));
				ptr += m_last_count * sizeof(uint64_t);
				memcpy(ptr, h_vals, m_last_count * sizeof(float));

				delete[] h_keys;
				delete[] h_vals;
			}

			// CPU-side deserialization: deserialize then copy to device
			void VoxelGPUManagerSortReduce::deserializeCPU(uint8_t* bin_data) {
				if (bin_data == nullptr) {
					return; // Invalid data
				}

				// int    m_last_count = 0;
				// uint64_t* d_keys_out = nullptr;
				// float* d_vals_out = nullptr;
				uint8_t* ptr = bin_data;				
				memcpy(&m_last_count, ptr, sizeof(int));
				m_max = m_last_count;
				ptr += sizeof(int);

				CUDA_CHECK_ERROR(cudaMalloc(&d_keys_out, m_last_count * sizeof(uint64_t)));
				CUDA_CHECK_ERROR(cudaMalloc(&d_vals_out, m_last_count * sizeof(float)));

				CUDA_CHECK_ERROR(cudaMemcpy(d_keys_out, ptr, m_last_count * sizeof(uint64_t), cudaMemcpyHostToDevice));
				ptr += m_last_count * sizeof(uint64_t);
				CUDA_CHECK_ERROR(cudaMemcpy(d_vals_out, ptr, m_last_count * sizeof(float), cudaMemcpyHostToDevice));
			}

			// CPU-side merge: combine key-value pairs from other manager
			void VoxelGPUManagerSortReduce::mergeCPU(common::vdb::VoxelSparseManager* _other) {
				// Attempt to dynamic_cast to VoxelGPUManagerSortReduce
				VoxelGPUManagerSortReduce* other = dynamic_cast<VoxelGPUManagerSortReduce*>(_other);
				if (!other) {
					printf("[VoxelGPUManagerSortReduce::mergeCPU] Error: Incompatible manager type\n");
					return;
				}

				if (other->m_last_count <= 0) {
					return; // Nothing to merge
				}

				const int n1 = m_last_count;
				const int n2 = other->m_last_count;
				const int n_total = n1 + n2;

				if ((size_t)n_total > m_max) {
					printf("[VoxelGPUManagerSortReduce::mergeCPU] ERROR: merged count=%d exceeds max=%zu\n",
						n_total, m_max);
					return;
				}

				// Concatenate key-value pairs
				CUDA_CHECK_ERROR(cudaMemcpy(d_keys, d_keys_out, n1 * sizeof(uint64_t), cudaMemcpyDeviceToDevice));
				CUDA_CHECK_ERROR(cudaMemcpy(d_vals, d_vals_out, n1 * sizeof(float), cudaMemcpyDeviceToDevice));
				CUDA_CHECK_ERROR(cudaMemcpy(d_keys + n1, other->d_keys_out, n2 * sizeof(uint64_t), cudaMemcpyDeviceToDevice));
				CUDA_CHECK_ERROR(cudaMemcpy(d_vals + n1, other->d_vals_out, n2 * sizeof(float), cudaMemcpyDeviceToDevice));

				// Sort
				cub::DeviceRadixSort::SortPairs(d_sort_temp, m_sort_temp_bytes,
					d_keys, d_keys_alt,
					d_vals, d_vals_alt,
					n_total);

				// Reduce-by-key
				CustomSum op_sum;
				cub::DeviceReduce::ReduceByKey(d_reduce_temp, m_reduce_temp_bytes,
					d_keys_alt, d_keys_out,
					d_vals_alt, d_vals_out,
					d_num_out,
					op_sum,
					n_total);

				int h_num_out = 0;
				CUDA_CHECK_ERROR(cudaMemcpy(&h_num_out, d_num_out, sizeof(int), cudaMemcpyDeviceToHost));
				m_last_count = h_num_out;
			}

			// Get min/max values from reduced voxel data using CUB
			void VoxelGPUManagerSortReduce::find_min_max(float& min_value, float& max_value)
			{
				if (m_last_count <= 0) {
					min_value = 0.0f;
					max_value = 0.0f;
					return;
				}

				// Use CUB's DeviceReduce::Min and Max on the reduced values
				cub::DeviceReduce::Min(d_minmax_temp, m_minmax_temp_bytes, 
					d_vals_out, d_min_out, m_last_count);
				
				// Reuse the same temp buffer for max (CUB allows this)
				cub::DeviceReduce::Max(d_minmax_temp, m_minmax_temp_bytes, 
					d_vals_out, d_max_out, m_last_count);

				// Copy results to host
				CUDA_CHECK_ERROR(cudaMemcpy(&min_value, d_min_out, sizeof(float), cudaMemcpyDeviceToHost));
				CUDA_CHECK_ERROR(cudaMemcpy(&max_value, d_max_out, sizeof(float), cudaMemcpyDeviceToHost));
			}

			void VoxelGPUManagerSortReduce::get_keys_values_from_device(uint64_t* h_keys, float* h_vals) {
				if (h_keys && d_keys_out) {
					CUDA_CHECK_ERROR(cudaMemcpy(h_keys, d_keys_out, m_last_count * sizeof(uint64_t), cudaMemcpyDeviceToHost));
				}
				if (h_vals && d_vals_out) {
					CUDA_CHECK_ERROR(cudaMemcpy(h_vals, d_vals_out, m_last_count * sizeof(float), cudaMemcpyDeviceToHost));
				}
			}

			// ================================================================================
			// GPU-side methods (no host memory transfers except for control)
			// ================================================================================

			// GPU-side serialization: allocate device memory and write serialized data
			void VoxelGPUManagerSortReduce::serialize(uint8_t* d_data) {

				uint8_t* ptr = d_data;				
				CUDA_CHECK_ERROR(cudaMemcpy(ptr, &m_last_count, sizeof(int), cudaMemcpyHostToDevice));
				ptr += sizeof(int);
				CUDA_CHECK_ERROR(cudaMemcpy(ptr, d_keys_out, m_last_count * sizeof(uint64_t), cudaMemcpyDeviceToDevice));
				ptr += m_last_count * sizeof(uint64_t);
				CUDA_CHECK_ERROR(cudaMemcpy(ptr, d_vals_out, m_last_count * sizeof(float), cudaMemcpyDeviceToDevice));

			}

			// GPU-side deserialization: read from device memory
			void VoxelGPUManagerSortReduce::deserialize(uint8_t* d_data) {

				if (d_data == nullptr) {
					return; // Invalid data
				}

				// int    m_last_count = 0;
				// uint64_t* d_keys_out = nullptr;
				// float* d_vals_out = nullptr;
				uint8_t* ptr = d_data;				
				CUDA_CHECK_ERROR(cudaMemcpy(&m_last_count, ptr, sizeof(int), cudaMemcpyDeviceToHost));
				m_max = m_last_count;
				ptr += sizeof(int);

				CUDA_CHECK_ERROR(cudaMalloc(&d_keys_out, m_last_count * sizeof(uint64_t)));
				CUDA_CHECK_ERROR(cudaMalloc(&d_vals_out, m_last_count * sizeof(float)));

				CUDA_CHECK_ERROR(cudaMemcpy(d_keys_out, ptr, m_last_count * sizeof(uint64_t), cudaMemcpyDeviceToDevice));
				ptr += m_last_count * sizeof(uint64_t);
				CUDA_CHECK_ERROR(cudaMemcpy(d_vals_out, ptr, m_last_count * sizeof(float), cudaMemcpyDeviceToDevice));
				
			}

			// GPU-side merge: combine with another manager without using host memory
			void VoxelGPUManagerSortReduce::merge(common::vdb::VoxelSparseManager* _other) {
				VoxelGPUManagerSortReduce* other = dynamic_cast<VoxelGPUManagerSortReduce*>(_other);
				if (!other) {
					printf("[VoxelGPUManagerSortReduce::merge] Error: Incompatible manager type\n");
					return;
				}

				if (other->m_last_count <= 0) {
					return; // Nothing to merge
				}

				const int n1 = m_last_count;
				const int n2 = other->m_last_count;
				const int n_total = n1 + n2;

				if ((size_t)n_total > m_max) {
					printf("[VoxelGPUManagerSortReduce::merge] ERROR: merged count=%d exceeds max=%zu\n",
						n_total, m_max);
					return;
				}

				// Concatenate key-value pairs from both managers
				// Copy this manager's data to beginning
				CUDA_CHECK_ERROR(cudaMemcpy(d_keys, d_keys_out, n1 * sizeof(uint64_t), cudaMemcpyDeviceToDevice));
				CUDA_CHECK_ERROR(cudaMemcpy(d_vals, d_vals_out, n1 * sizeof(float), cudaMemcpyDeviceToDevice));

				// Copy other manager's data after this manager's data
				CUDA_CHECK_ERROR(cudaMemcpy(d_keys + n1, other->d_keys_out, n2 * sizeof(uint64_t), cudaMemcpyDeviceToDevice));
				CUDA_CHECK_ERROR(cudaMemcpy(d_vals + n1, other->d_vals_out, n2 * sizeof(float), cudaMemcpyDeviceToDevice));

				// Sort combined data
				cub::DeviceRadixSort::SortPairs(d_sort_temp, m_sort_temp_bytes,
					d_keys, d_keys_alt,
					d_vals, d_vals_alt,
					n_total);

				// Reduce-by-key
				CustomSum op_sum;
				cub::DeviceReduce::ReduceByKey(d_reduce_temp, m_reduce_temp_bytes,
					d_keys_alt, d_keys_out,
					d_vals_alt, d_vals_out,
					d_num_out,
					op_sum,
					n_total);

				int h_num_out = 0;
				CUDA_CHECK_ERROR(cudaMemcpy(&h_num_out, d_num_out, sizeof(int), cudaMemcpyDeviceToHost));
				m_last_count = h_num_out;
			}

		}
	}
}
