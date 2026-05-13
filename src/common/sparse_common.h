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

// NanoVDB includes
#include <nanovdb/NanoVDB.h>
#include <nanovdb/GridHandle.h>
#include <nanovdb/tools/GridBuilder.h>
#include <nanovdb/tools/CreateNanoGrid.h>

// OpenVDB includes
#ifdef WITH_OPENVDB
#	include <openvdb/openvdb.h>
#	include <openvdb/Grid.h>
#	include <openvdb/tree/Tree.h>
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
					: i(_i), j(_j), k(_k), value(_value) {}
			};

			// /**
			//  * @brief Sparse grid representation storing only non-empty voxels
			//  * 
			//  * Stores sparse voxel data as (i, j, k, value) tuples, replacing
			//  * the NanoVDB and OpenVDB grid representations with a simpler format.
			//  * Used for VDBParticleType::eSparse representation.
			//  */
			// struct SparseParticles
			// {
			// 	std::vector<Voxel> voxels;  ///< Collection of non-empty voxels

			// 	/**
			// 	 * @brief Add a voxel to the sparse grid
			// 	 * @param i X coordinate
			// 	 * @param j Y coordinate
			// 	 * @param k Z coordinate
			// 	 * @param value Voxel value
			// 	 */
			// 	void add_voxel(int i, int j, int k, float value) {
			// 		voxels.emplace_back(i, j, k, value);
			// 	}

			// 	/**
			// 	 * @brief Clear all voxel data
			// 	 */
			// 	void clear() {
			// 		voxels.clear();
			// 	}

			// 	/**
			// 	 * @brief Get the number of active voxels
			// 	 * @return Number of voxels in the sparse grid
			// 	 */
			// 	size_t size() const {
			// 		return voxels.size();
			// 	}

			// 	/**
			// 	 * @brief Reserve space for voxels to avoid reallocation
			// 	 * @param count Expected number of voxels
			// 	 */
			// 	void reserve(size_t count) {
			// 		voxels.reserve(count);
			// 	}

			// 	/**
			// 	 * @brief Serialize sparse grid data to a binary file
			// 	 * @param filename Path to output file
			// 	 */
			// 	void serialize(const std::string& filename) const {
			// 		std::ofstream out(filename, std::ios::binary);
			// 		if (!out) {
			// 			throw std::runtime_error("Failed to open file for writing");
			// 		}

			// 		// Write the number of voxels
			// 		size_t voxel_count = voxels.size();
			// 		out.write(reinterpret_cast<const char*>(&voxel_count), sizeof(voxel_count));

			// 		// Write all voxel data
			// 		out.write(reinterpret_cast<const char*>(voxels.data()), voxel_count * sizeof(Voxel));

			// 		out.close();
			// 	}

			// 	/**
			// 	 * @brief Deserialize sparse grid data from a binary file
			// 	 * @param filename Path to input file
			// 	 */
			// 	void deserialize(const std::string& filename) {
			// 		std::ifstream in(filename, std::ios::binary);
			// 		if (!in) {
			// 			throw std::runtime_error("Failed to open file for reading");
			// 		}

			// 		// Read the number of voxels
			// 		size_t voxel_count;
			// 		in.read(reinterpret_cast<char*>(&voxel_count), sizeof(voxel_count));

			// 		// Read all voxel data
			// 		voxels.resize(voxel_count);
			// 		in.read(reinterpret_cast<char*>(voxels.data()), voxel_count * sizeof(Voxel));

			// 		in.close();
			// 	}

			// 	/**
			// 	 * @brief Serialize sparse grid data to a binary buffer (for MPI/network transfer)
			// 	 * @param bin_data Output buffer to store serialized data
			// 	 */
			// 	void serialize(std::vector<uint8_t>& bin_data) const {
			// 		size_t voxel_count = voxels.size();
			// 		size_t total_size = sizeof(voxel_count) + voxel_count * sizeof(Voxel);
					
			// 		bin_data.resize(total_size);
			// 		uint8_t* ptr = bin_data.data();

			// 		// Write voxel count
			// 		memcpy(ptr, &voxel_count, sizeof(voxel_count));
			// 		ptr += sizeof(voxel_count);

			// 		// Write voxel data
			// 		memcpy(ptr, voxels.data(), voxel_count * sizeof(Voxel));
			// 	}

			// 	/**
			// 	 * @brief Deserialize sparse grid data from a binary buffer
			// 	 * @param bin_data Input buffer containing serialized data
			// 	 */
			// 	void deserialize(const std::vector<uint8_t>& bin_data) {
			// 		const uint8_t* ptr = bin_data.data();

			// 		// Read voxel count
			// 		size_t voxel_count;
			// 		memcpy(&voxel_count, ptr, sizeof(voxel_count));
			// 		ptr += sizeof(voxel_count);

			// 		// Read voxel data
			// 		voxels.resize(voxel_count);
			// 		memcpy(voxels.data(), ptr, voxel_count * sizeof(Voxel));
			// 	}

			// 	/**
			// 	 * @brief Merge another SparseParticles object into this one
			// 	 * @param other Source sparse grid to merge
			// 	 * 
			// 	 * Appends all voxels from the other sparse grid into this one.
			// 	 */
			// 	void merge(const SparseParticles& other) {
			// 		voxels.insert(voxels.end(), other.voxels.begin(), other.voxels.end());
			// 	}
			// };

			// ---------------------------------------------
			// Key packing (lossless within chosen range)
			// 21 bits per axis => range [-2^20, 2^20-1]
			// i,j,k must be within [-1,048,576, 1,048,575]
			// ---------------------------------------------
			//static constexpr int   COORD_BITS = 21;
			//static constexpr int   COORD_BIAS = 1 << (COORD_BITS - 1); // 2^20
			//static constexpr uint64_t COORD_MASK = (1ull << COORD_BITS) - 1ull;

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
			inline void unpackCoord3(uint64_t key, int &x, int &y, int &z)
			{
				uint64_t ux = (key) & COORD_MASK;
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
			class VoxelOpenMPManager: public common::vdb::VoxelSparseManager {
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
				void serialize(uint8_t *bin_data) override;
				
				// Deserialization: read voxel data from binary buffer
				void deserialize(uint8_t *bin_data) override;
				
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

			// ---------------------------------------------
			// CPU-based voxel manager using OpenMP for parallelization
			// ---------------------------------------------
			class VoxelCPUManager : public common::vdb::VoxelSparseManager {
			private:
				VoxelHashEntry* hash_table;
				unsigned int table_size;
				int insert_count;

				// Hash function for 3D coordinates (CPU version)
				inline unsigned int hash3D_cpu(int i, int j, int k) const {
					unsigned int hash = 73856093u * i ^ 19349663u * j ^ 83492791u * k;
					return hash % table_size;
				}

				inline void insertOrUpdatePackedSequential(uint64_t key, float value) {
					int x, y, z;
					unpackCoord3(key, x, y, z);
					unsigned int slot = hash3D_cpu(x, y, z);

					for (unsigned int probe = 0; probe < table_size; probe++) {
						if (hash_table[slot].occupied == 0) {
							hash_table[slot].i = x;
							hash_table[slot].j = y;
							hash_table[slot].k = z;
							hash_table[slot].value = value;
							hash_table[slot].occupied = 1;
							insert_count++;
							return;
						}

						if (hash_table[slot].i == x &&
							hash_table[slot].j == y &&
							hash_table[slot].k == z) {
							hash_table[slot].value += value;
							return;
						}

						slot++;
						if (slot == table_size) slot = 0;
					}
				}

			public:
				VoxelCPUManager() : table_size(0), insert_count(0), hash_table(nullptr), common::vdb::VoxelSparseManager() {}

				VoxelCPUManager(unsigned int expected_voxels) : common::vdb::VoxelSparseManager() {
					init(expected_voxels);
				}

				~VoxelCPUManager() {
					if (hash_table) delete[] hash_table;
				}

			public:
				void init(unsigned int expected_voxels) override {
					// Size hash table with load factor consideration
					table_size = (unsigned int)(expected_voxels / HASH_TABLE_LOAD_FACTOR);
					insert_count = 0;

					// Allocate hash table on CPU
					hash_table = new VoxelHashEntry[table_size];

					// Initialize hash table
#pragma omp parallel for
					for (unsigned int i = 0; i < table_size; i++) {
						hash_table[i].i = EMPTY_KEY;
						hash_table[i].j = EMPTY_KEY;
						hash_table[i].k = EMPTY_KEY;
						hash_table[i].value = 0.0f;
						hash_table[i].occupied = 0;
					}
				}

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
					VoxelCPUManager temp_manager;
					temp_manager.deserialize(bin_data);
					this->merge(&temp_manager);
				}

				// Insert or update voxels using OpenMP parallelization with thread-local pre-aggregation
				void insertOrUpdate(Voxel* h_voxels, int num_voxels);

				// Extract all voxels from hash table using OpenMP
				int extractAll(Voxel** h_output_voxels);
			};


#ifdef WITH_GPU_CUDA
			// ---------------------------------------------
			// GPU-based voxel manager with hash table
			// ---------------------------------------------
			class VoxelGPUManager : public common::vdb::VoxelSparseManager {
			private:
				VoxelHashEntry* d_hash_table;
				unsigned int table_size;
				int* d_insert_count;
				
				// Device buffers for sorting/aggregation
				Voxel* d_voxels;
				size_t d_voxels_capacity;
				uint64_t* d_keys;
				float* d_vals;
				uint64_t* d_keys_out;
				float* d_vals_out;
				size_t d_work_capacity;
				
				// Pinned host memory for faster transfers
				Voxel* h_pinned_voxels;
				size_t h_pinned_capacity;
				
				// CUB temporary storage
				void* d_temp_storage;
				size_t temp_storage_bytes;

				cudaEvent_t insert_start_event;
				cudaEvent_t insert_stop_event;

				void ensureVoxelBuffer(size_t required_voxels);
				void ensureWorkBuffers(size_t required_size);
				void ensurePinnedBuffer(size_t required_size);
				
			public:
				VoxelGPUManager() : table_size(0), d_hash_table(nullptr), d_insert_count(nullptr),
					d_voxels(nullptr), d_voxels_capacity(0),
					d_keys(nullptr), d_vals(nullptr), d_keys_out(nullptr), d_vals_out(nullptr), d_work_capacity(0),
					h_pinned_voxels(nullptr), h_pinned_capacity(0),
					d_temp_storage(nullptr), temp_storage_bytes(0),
					common::vdb::VoxelSparseManager() {}
				
				VoxelGPUManager(unsigned int expected_voxels);
				~VoxelGPUManager();

			public:
				void init(unsigned int expected_voxels) override;
				void insertOrUpdatePackedSequential(uint64_t key, float value) override;
				void serialize(uint8_t* bin_data) override;
				void deserialize(uint8_t* bin_data) override;
				void merge(common::vdb::VoxelSparseManager* other) override;
				
				void merge(uint8_t* bin_data) override {
					VoxelGPUManager temp_manager;
					temp_manager.deserialize(bin_data);
					this->merge(&temp_manager);
				}
				
				size_t mem_size() const override {
					return sizeof(unsigned int) + sizeof(int) + 
						   table_size * sizeof(VoxelHashEntry) +
						   d_voxels_capacity * sizeof(Voxel) +
						   d_work_capacity * (2 * sizeof(uint64_t) + 2 * sizeof(float)) +
						   h_pinned_capacity * sizeof(Voxel) +
						   temp_storage_bytes;
				}
				
			public:
				// Optimized insert with pre-aggregation using sorting + reduce-by-key
				void insertOrUpdate(Voxel* h_voxels, int num_voxels);
				
				// Extract all voxels from hash table (optimized with pinned memory)
				int extractAll(Voxel** h_output_voxels);
				
				// Clear hash table
				void clear();
			};

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
				void serialize(uint8_t *bin_data) override;
				
				// Deserialization: read voxel data from binary buffer
				void deserialize(uint8_t *bin_data) override;
				
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
					VoxelGPUManagerSortReduce temp_manager;
					temp_manager.deserialize(bin_data);
					this->merge(&temp_manager);
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
				int*   d_bbox_min_orig   = nullptr;
				float* d_offset_position = nullptr;
				float* d_bbox_sphere_pos = nullptr;
			};
#endif // WITH_GPU_CUDA

			// ---------------------------------------------
			// NanoVDB-based voxel manager (CPU construction)
			// ---------------------------------------------
			class VoxelNanoVDBManager : public common::vdb::VoxelSparseManager {
			private:
				nanovdb::GridHandle<> gridHandle;
#if OPENVDB_VERSION == 11
				std::shared_ptr<nanovdb::build::FloatGrid> persistentGrid;
#else
				std::shared_ptr<nanovdb::tools::build::FloatGrid> persistentGrid;
#endif
				
			public:
				VoxelNanoVDBManager() : common::vdb::VoxelSparseManager() {
					// Initialize persistent destination grid
#if OPENVDB_VERSION == 11
					persistentGrid = std::make_shared<nanovdb::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#else
					persistentGrid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#endif
				}
				
				~VoxelNanoVDBManager() {}

			public:
				void init(unsigned int expected_voxels) override {
					// NanoVDB doesn't need pre-allocation
				}
				
				void insertOrUpdatePackedSequential(uint64_t key, float value) override {
					// Not implemented for NanoVDB (uses batch insertion)
				}
				
				void serialize(uint8_t* bin_data) override;
				void deserialize(uint8_t* bin_data) override;
				void merge(common::vdb::VoxelSparseManager* other) override;
				
				void merge(uint8_t* bin_data) override {
					VoxelNanoVDBManager temp_manager;
					temp_manager.deserialize(bin_data);
					this->merge(&temp_manager);
				}
				
				size_t mem_size() const override {
					return gridHandle.size();
				}
				
			public:
				// Insert or update voxels (like VoxelGPUManager::insertOrUpdate)
				void buildFromVoxels(Voxel* h_voxels, int num_voxels);
				
				// Extract voxels back from NanoVDB grid
				int extractAll(Voxel** h_output_voxels);
				
				// Clear persistent grid
				void clear() {
#if OPENVDB_VERSION == 11
					persistentGrid = std::make_shared<nanovdb::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#else
					persistentGrid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#endif
					gridHandle = nanovdb::GridHandle<>();
				}
				
				// Alias for consistency with VoxelGPUManager
				void insertOrUpdate(Voxel* h_voxels, int num_voxels) {
					buildFromVoxels(h_voxels, num_voxels);
				}
				
				// Get grid handle for GPU operations
				const nanovdb::GridHandle<>& getGridHandle() const {
					return gridHandle;
				}
			};

#ifdef WITH_OPENVDB
			// ---------------------------------------------
			// OpenVDB-based voxel manager (CPU construction)
			// ---------------------------------------------
			class VoxelOpenVDBManager : public common::vdb::VoxelSparseManager {
			private:
				openvdb::FloatGrid::Ptr persistentGrid;
				
			public:
				VoxelOpenVDBManager() : common::vdb::VoxelSparseManager() {
					// Initialize persistent destination grid
					persistentGrid = openvdb::FloatGrid::create(0.0f);
					persistentGrid->setName("density");
					persistentGrid->setGridClass(openvdb::GRID_FOG_VOLUME);
				}
				
				~VoxelOpenVDBManager() {}

			public:
				void init(unsigned int expected_voxels) override {
					// OpenVDB doesn't need pre-allocation
				}
				
				void insertOrUpdatePackedSequential(uint64_t key, float value) override {
					// Not implemented for OpenVDB (uses batch insertion)
				}
				
				void serialize(uint8_t* bin_data) override;
				void deserialize(uint8_t* bin_data) override;
				void merge(common::vdb::VoxelSparseManager* other) override;
				
				void merge(uint8_t* bin_data) override {
					VoxelOpenVDBManager temp_manager;
					temp_manager.deserialize(bin_data);
					this->merge(&temp_manager);
				}
				
				size_t mem_size() const override {
					return persistentGrid->memUsage();
				}
				
			public:
				// Insert or update voxels (like VoxelNanoVDBManager::buildFromVoxels)
				void buildFromVoxels(Voxel* h_voxels, int num_voxels);
				
				// Extract voxels back from OpenVDB grid
				int extractAll(Voxel** h_output_voxels);
				
				// Clear persistent grid
				void clear() {
					persistentGrid = openvdb::FloatGrid::create(0.0f);
					persistentGrid->setName("density");
					persistentGrid->setGridClass(openvdb::GRID_FOG_VOLUME);
				}
				
				// Alias for consistency with other managers
				void insertOrUpdate(Voxel* h_voxels, int num_voxels) {
					buildFromVoxels(h_voxels, num_voxels);
				}
				
				// Get grid pointer for external operations
				openvdb::FloatGrid::Ptr getGrid() const {
					return persistentGrid;
				}
			};
#endif // WITH_OPENVDB
		} //sparse
	}// vdb
} //common