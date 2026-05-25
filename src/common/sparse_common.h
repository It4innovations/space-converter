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

 //#include "convert_vdb.h"

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <memory>
#include "data_common.h"

// CUDA includes
#ifdef WITH_GPU_CUDA
#include <cuda_runtime.h>
#endif

#define COORD_BITS 21
#define COORD_BIAS (1 << (COORD_BITS - 1)) // 2^20
#define COORD_MASK ((1ull << COORD_BITS) - 1ull)

// Constants for hash table management
#define EMPTY_KEY -999999
#define HASH_TABLE_LOAD_FACTOR 0.5f

namespace common {
	namespace vdb {
		namespace sparse {
			/**
			 * @brief Single voxel entry in the sparse grid
			 */
			struct Voxel {
				int i;         ///< X coordinate (voxel index)
				int j;         ///< Y coordinate (voxel index)
				int k;         ///< Z coordinate (voxel index)
				float value;   ///< Voxel value (density, temperature, etc.)

#ifdef __CUDACC__
				__host__ __device__
#endif
					Voxel() : i(0), j(0), k(0), value(0.0f) {}

#ifdef __CUDACC__
				__host__ __device__
#endif					
					Voxel(int _i, int _j, int _k, float _value)
					: i(_i), j(_j), k(_k), value(_value) {
				}
			};

#ifdef __CUDACC__
			__host__ __device__
#endif			
				inline uint64_t packCoord3(int x, int y, int z)
			{
				// Note: for production, you may want to clamp or assert in debug.
				uint64_t ux = (uint64_t)(x + COORD_BIAS) & COORD_MASK;
				uint64_t uy = (uint64_t)(y + COORD_BIAS) & COORD_MASK;
				uint64_t uz = (uint64_t)(z + COORD_BIAS) & COORD_MASK;

				return (ux) | (uy << COORD_BITS) | (uz << (2 * COORD_BITS));
			}

#ifdef __CUDACC__			
			__host__ __device__
#endif			
				inline void unpackCoord3(uint64_t key, int& x, int& y, int& z)
			{
				uint64_t ux = (key)&COORD_MASK;
				uint64_t uy = (key >> COORD_BITS) & COORD_MASK;
				uint64_t uz = (key >> (2 * COORD_BITS)) & COORD_MASK;

				x = (int)ux - COORD_BIAS;
				y = (int)uy - COORD_BIAS;
				z = (int)uz - COORD_BIAS;
			}

			// Hash map entry for voxel storage
			struct VoxelHashEntry {
				int i, j, k;
				float value;
				int occupied;  // 0 = empty, -1 = being written, 1 = fully written
			};

			// ---------------------------------------------
			// OpenMP-based CPU voxel manager
			// ---------------------------------------------
			class VoxelOpenMPManager : public common::vdb::VoxelSparseManager {
			public:
				VoxelHashEntry* hash_table;
				unsigned int table_size;
				int insert_count;

				// Hash function for 3D coordinates (CPU version)
				inline unsigned int hash3D_cpu(int i, int j, int k) const;

			public:
				VoxelOpenMPManager();
				~VoxelOpenMPManager();

			public:
				void init(unsigned int expected_voxels) override;

				// Helper for sequential insertion
				void insertOrUpdatePackedSequential(uint64_t key, float value) override;

				// Serialization: write current voxel data to binary buffer
				void serialize(uint8_t* bin_data) override;

				// Deserialization: read voxel data from binary buffer
				void deserialize(uint8_t* bin_data) override;

				// Merge: combine voxels from another manager (accumulate values)
				void merge(common::vdb::VoxelSparseManager* other) override;

				size_t mem_size() const override {
					return sizeof(int) + sizeof(unsigned int) + sizeof(VoxelHashEntry) * table_size;
				}

				void merge(uint8_t* bin_data) override {
					// Deserialize the incoming data into a temporary manager and then merge
					VoxelOpenMPManager temp_manager;
					temp_manager.deserialize(bin_data);
					this->merge(&temp_manager);
				}

			public:

				// Insert or update voxels using OpenMP parallelization with thread-local pre-aggregation
				void insertOrUpdate(Voxel* h_voxels, int num_voxels);

				// Extract all voxels from hash table using OpenMP
				int extractAll(Voxel** h_output_voxels);

				// Clear hash table
				void clear();
			};

#ifdef WITH_GPU_CUDA
			// ---------------------------------------------
			// GPU-based voxel manager using sort+reduce
			// ---------------------------------------------
			class VoxelGPUManagerSortReduce : public common::vdb::VoxelSparseManager {
			public:
				VoxelGPUManagerSortReduce();
				~VoxelGPUManagerSortReduce();

			public:
				void init(unsigned int expected_voxels) override;

				// Helper for sequential insertion
				//void insertOrUpdatePackedSequential(uint64_t key, float value) override;

				// Serialization: write current voxel data to binary buffer
				void serialize(uint8_t* bin_data) override;

				// Deserialization: read voxel data from binary buffer
				void deserialize(uint8_t* bin_data) override;

				// Merge: combine voxels from another manager (accumulate values)
				void merge(common::vdb::VoxelSparseManager* other) override;

				size_t mem_size() const override {
					// TODO: this is a rough estimate. For accurate memory usage, we would need to track the actual number of voxels after reduction.
					//return sizeof(size_t) + sizeof(size_t) + sizeof(uint64_t) * m_max + sizeof(float) * m_max;

					// int    m_last_count = 0;
					// uint64_t* d_keys_out = nullptr;
					// float* d_vals_out = nullptr;
					return sizeof(int) + sizeof(uint64_t) * m_last_count + sizeof(float) * m_last_count;
				}

				void merge(uint8_t* bin_data) override {
					// Deserialize the incoming data into a temporary manager and then merge
					// Note: bin_data is assumed to be host memory (from std::vector in MPI operations)
					VoxelGPUManagerSortReduce temp_manager;
					temp_manager.deserializeCPU(bin_data);
					this->mergeCPU(&temp_manager);
				}


			public:
				// Accumulate a batch of voxels: output becomes unique (i,j,k) with summed values
				// Returns number of unique voxels in the batch.
				int insertOrUpdate(const Voxel* h_voxels, int num_voxels);

				// Extract the last accumulated unique voxels back to host
				int extractAll(Voxel** h_output_voxels);

				// CPU-side serialization: write current voxel data to binary buffer
				void serializeCPU(uint8_t* bin_data);

				// CPU-side deserialization: read voxel data from binary buffer
				void deserializeCPU(uint8_t* bin_data);

				// CPU-side merge: combine voxels from another manager (uses host memory)
				void mergeCPU(common::vdb::VoxelSparseManager* other);

				void mergeCPU(uint8_t* bin_data) {
					// Deserialize the incoming data into a temporary manager and then merge
					VoxelGPUManagerSortReduce temp_manager;
					temp_manager.deserializeCPU(bin_data);
					this->mergeCPU(&temp_manager);
				}

				// GPU-side serialization: write to device memory (device-to-device)
				void serializeGPU(uint8_t* d_data);

				// GPU-side deserialization: read from device memory (device-to-device)
				void deserializeGPU(uint8_t* d_data);

				void get_keys_values_from_device(uint64_t* h_keys, float* h_vals);

				// Get min/max values from reduced voxel data using CUB
				void find_min_max(float& min_value, float& max_value);

				int update(size_t count);

			public:
				size_t m_max = 0;
				int    m_last_count = 0;

				uint64_t* d_keys_out = nullptr;
				float* d_vals_out = nullptr;

				//Voxel* d_inVoxels = nullptr;

				uint64_t* d_keys = nullptr;
				float* d_vals = nullptr;

				uint64_t* d_keys_alt = nullptr;
				float* d_vals_alt = nullptr;

				int* d_num_out = nullptr;

				void* d_sort_temp = nullptr;
				size_t    m_sort_temp_bytes = 0;

				void* d_reduce_temp = nullptr;
				size_t    m_reduce_temp_bytes = 0;

				// Min/max reduction support
				void* d_minmax_temp = nullptr;
				size_t m_minmax_temp_bytes = 0;
				float* d_min_out = nullptr;
				float* d_max_out = nullptr;

				// Atomic counter: number of particles that passed should_process
				uint64_t* d_particle_count = nullptr;

				// Persistent device copies of per-call host parameters
				int* d_bbox_min_orig = nullptr;
				float* d_offset_position = nullptr;
				float* d_bbox_sphere_pos = nullptr;
			};
#endif // WITH_GPU_CUDA

		} //sparse
	}// vdb
} //common