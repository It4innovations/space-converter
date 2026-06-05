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
#include <cstdint>
#include <unordered_map>
#include "convert_vdb.h"

#include <omp.h>

#ifdef WITH_OPENVDB
#	include <openvdb/openvdb.h>
#	include <openvdb/tools/Composite.h>
#endif

 // Constants
#define EMPTY_KEY -999999
#define HASH_TABLE_LOAD_FACTOR 0.5f

namespace common {
	namespace vdb {
		namespace sparse {

			// Hash function for 3D coordinates (CPU version)
			inline unsigned int VoxelOpenMPManager::hash3D_cpu(int i, int j, int k) const {
				unsigned int hash = 73856093u * i ^ 19349663u * j ^ 83492791u * k;
				return hash % table_size;
			}

			VoxelOpenMPManager::VoxelOpenMPManager() : table_size(0), insert_count(0), hash_table(nullptr), common::vdb::VoxelSparseManager() {
			}

			VoxelOpenMPManager::~VoxelOpenMPManager() {
				if (hash_table != nullptr) {
					delete[] hash_table;
				}
			}

			void VoxelOpenMPManager::init(unsigned int expected_voxels) {
				// Size hash table with load factor consideration
				table_size = (unsigned int)(expected_voxels / HASH_TABLE_LOAD_FACTOR);
				insert_count = 0;

				// Allocate hash table on CPU
				if (hash_table != nullptr) {
					delete[] hash_table;
				}
				hash_table = new VoxelHashEntry[table_size];

				// Initialize hash table using OpenMP
#pragma omp parallel for
				for (unsigned int i = 0; i < table_size; i++) {
					hash_table[i].i = EMPTY_KEY;
					hash_table[i].j = EMPTY_KEY;
					hash_table[i].k = EMPTY_KEY;
					hash_table[i].value = 0.0f;
					hash_table[i].occupied = 0;
				}
			}

			// Helper for sequential insertion
			void VoxelOpenMPManager::insertOrUpdatePackedSequential(uint64_t key, float value) {
				int x, y, z;
				unpackCoord3(key, x, y, z);
				unsigned int slot = hash3D_cpu(x, y, z);

				for (unsigned int probe = 0; probe < table_size; probe++) {
					unsigned int idx = (slot + probe) % table_size;
					if (hash_table[idx].occupied == 0) {
						hash_table[idx].i = x;
						hash_table[idx].j = y;
						hash_table[idx].k = z;
						hash_table[idx].value = value;
						hash_table[idx].occupied = 1;
						insert_count++;
						return;
					}

					if (hash_table[idx].i == x &&
						hash_table[idx].j == y &&
						hash_table[idx].k == z) {
						hash_table[idx].value += value;
						return;
					}

				}
			}

			// Insert or update voxels using OpenMP parallelization with thread-local pre-aggregation
			void VoxelOpenMPManager::insertOrUpdate(void* _voxels, size_t num_voxels) {
				if (num_voxels <= 0) return;

				Voxel* h_voxels = (Voxel*)_voxels;

				auto start_total = omp_get_wtime();

				// OPTIMIZATION: Use thread-local hash maps to pre-aggregate voxels
				// This drastically reduces hash table contention
				auto start_preagg = omp_get_wtime();
				const int thread_count = omp_get_max_threads();
				std::vector<std::unordered_map<uint64_t, float>> thread_maps(thread_count);

#pragma omp parallel
				{
					int tid = omp_get_thread_num();
					auto& local_map = thread_maps[tid];
					local_map.reserve((size_t)num_voxels / (size_t)thread_count + 64);

#pragma omp for schedule(static)
					for (size_t idx = 0; idx < num_voxels; idx++) {
						const Voxel& v = h_voxels[idx];
						const uint64_t key = packCoord3(v.i, v.j, v.k);
						local_map[key] += v.value;
					}
				}
				auto end_preagg = omp_get_wtime();

				// Merge thread-local maps into single map
				auto start_merge = omp_get_wtime();
				size_t total_local_unique = 0;
				for (const auto& m : thread_maps) {
					total_local_unique += m.size();
				}
				std::unordered_map<uint64_t, float> merged;
				merged.reserve(total_local_unique + total_local_unique / 4 + 1);

				for (const auto& m : thread_maps) {
					for (const auto& kv : m) {
						merged[kv.first] += kv.second;
					}
				}
				auto end_merge = omp_get_wtime();

				// Sequential insertion into hash table (no contention)
				auto start_hash_insert = omp_get_wtime();
				for (const auto& kv : merged) {
					insertOrUpdatePackedSequential(kv.first, kv.second);
				}
				auto end_hash_insert = omp_get_wtime();
				auto end_total = omp_get_wtime();
				auto preagg_duration = (end_preagg - start_preagg) * 1000.0;
				auto merge_duration = (end_merge - start_merge) * 1000.0;
				auto hash_insert_duration = (end_hash_insert - start_hash_insert) * 1000.0;
				auto total_duration = (end_total - start_total) * 1000.0;

				printf("[insertOrUpdate] Pre-agg: %.3f ms, Merge: %.3f ms, Hash insert: %.3f ms, Total: %.3f ms\n",
					preagg_duration, merge_duration,
					hash_insert_duration, total_duration);
			}

			// Extract all voxels from hash table using OpenMP
//			int VoxelOpenMPManager::extractAll(Voxel** h_output_voxels) {
//				auto start_total = omp_get_wtime();
//
//				// Thread-local storage for voxel collection
//				std::vector<std::vector<Voxel>> thread_voxels(omp_get_max_threads());
//
//				// Parallel extraction - each thread collects its voxels
//#pragma omp parallel
//				{
//					int thread_id = omp_get_thread_num();
//
//#pragma omp for schedule(static)
//					for (unsigned int idx = 0; idx < table_size; idx++) {
//						if (hash_table[idx].occupied == 1) {
//							thread_voxels[thread_id].push_back(Voxel(
//								hash_table[idx].i,
//								hash_table[idx].j,
//								hash_table[idx].k,
//								hash_table[idx].value
//							));
//						}
//					}
//				}
//
//				// Combine thread-local results
//				int total_count = 0;
//				for (const auto& vec : thread_voxels) {
//					total_count += vec.size();
//				}
//
//				*h_output_voxels = new Voxel[total_count];
//				int offset = 0;
//				for (const auto& vec : thread_voxels) {
//					std::copy(vec.begin(), vec.end(), *h_output_voxels + offset);
//					offset += vec.size();
//				}
//
//				auto end_total = omp_get_wtime();
//				auto total_duration = (end_total - start_total) * 1000.0;
//				printf("[extractAll] Total: %.3f ms\n", total_duration);
//
//				return total_count;
//			}

			// Clear hash table
			void VoxelOpenMPManager::clear() {
#pragma omp parallel for
				for (unsigned int i = 0; i < table_size; i++) {
					hash_table[i].i = EMPTY_KEY;
					hash_table[i].j = EMPTY_KEY;
					hash_table[i].k = EMPTY_KEY;
					hash_table[i].value = 0.0f;
					hash_table[i].occupied = 0;
				}
				insert_count = 0;
			}

			// Serialization: write current voxel data to binary buffer
			void VoxelOpenMPManager::serialize(uint8_t* bin_data) {
				// int insert_count;
				// unsigned int table_size;
				// VoxelHashEntry* hash_table;

				//size_t expected_size = sizeof(int) + sizeof(unsigned int) + sizeof(VoxelHashEntry) * table_size;
				uint8_t* ptr = bin_data;

				memcpy(ptr, &insert_count, sizeof(int));
				ptr += sizeof(int);
				memcpy(ptr, &table_size, sizeof(unsigned int));
				ptr += sizeof(unsigned int);
				memcpy(ptr, hash_table, sizeof(VoxelHashEntry) * table_size);
			}

			// Deserialization: read voxel data from binary buffer
			void VoxelOpenMPManager::deserialize(uint8_t* bin_data) {
				// int insert_count;
				// unsigned int table_size;
				// VoxelHashEntry* hash_table;

				// size_t expected_size = sizeof(int) + sizeof(unsigned int) + sizeof(VoxelHashEntry) * table_size;
				// if (bin_data.size() != expected_size) {
				// 	printf("[VoxelOpenMPManager::deserialize] Error: Expected size %zu, got %zu\n", expected_size, bin_data.size());
				// 	return; // Invalid data
				// }

				const uint8_t* ptr = bin_data;
				memcpy(&insert_count, ptr, sizeof(int));
				ptr += sizeof(int);
				memcpy(&table_size, ptr, sizeof(unsigned int));
				ptr += sizeof(unsigned int);

				if (hash_table) {
					delete[] hash_table;
				}
				hash_table = new VoxelHashEntry[table_size];
				memcpy(hash_table, ptr, sizeof(VoxelHashEntry) * table_size);
			}

			// Merge: combine voxels from another manager (accumulate values)
			void VoxelOpenMPManager::merge(common::vdb::VoxelSparseManager* _other) {

				// Attempt to dynamic_cast to VoxelOpenMPManager
				VoxelOpenMPManager* other = dynamic_cast<VoxelOpenMPManager*>(_other);
				if (!other) {
					printf("[VoxelOpenMPManager::merge] Error: Incompatible manager type\n");
					return;
				}

				// Parallel insertion directly from other's hash table using thread-local maps
				const int thread_count = omp_get_max_threads();
				std::vector<std::unordered_map<uint64_t, float>> thread_maps(thread_count);

#pragma omp parallel
				{
					int tid = omp_get_thread_num();
					auto& local_map = thread_maps[tid];

#pragma omp for schedule(static)
					for (unsigned int i = 0; i < other->table_size; i++) {
						if (other->hash_table[i].occupied == 1) {
							uint64_t key = packCoord3(other->hash_table[i].i, other->hash_table[i].j, other->hash_table[i].k);
							local_map[key] += other->hash_table[i].value;
						}
					}
				}

				// Merge thread-local maps
				size_t total_local_unique = 0;
				for (const auto& m : thread_maps) {
					total_local_unique += m.size();
				}
				std::unordered_map<uint64_t, float> merged;
				merged.reserve(total_local_unique + total_local_unique / 4 + 1);

				for (const auto& m : thread_maps) {
					for (const auto& kv : m) {
						merged[kv.first] += kv.second;
					}
				}

				// Insert into this manager's hash table
				for (const auto& kv : merged) {
					insertOrUpdatePackedSequential(kv.first, kv.second);
				}
			}

#ifdef WITH_NANOVDB
			// =====================================================================
			// VoxelNanoVDBManager Implementation
			// =====================================================================

			// Helper to get the BuildType from a Grid
			template<typename GridT>
			struct GridBuildType;

			template<typename BuildT>
			struct GridBuildType<nanovdb::tools::build::Grid<BuildT>> {
				using type = BuildT;
			};

			template<typename GridT>
			using GridBuildType_t = typename GridBuildType<GridT>::type;

			// Helper to convert float to the appropriate value type
			template<typename ValueT>
			inline ValueT floatToValue(float f) {
				return ValueT(f);
			}

			template<>
			inline float floatToValue<float>(float f) {
				return f;
			}

			template<>
			inline nanovdb::math::Vec3<float> floatToValue<nanovdb::math::Vec3<float>>(float f) {
				return nanovdb::math::Vec3<float>(f, f, f);
			}

			template<typename T>
			VoxelNanoVDBManager<T>::VoxelNanoVDBManager() : insert_count(0), common::vdb::VoxelSparseManager() {
				using BuildType = GridBuildType_t<T>;
				nano_grid = std::make_shared<T>(floatToValue<BuildType>(0.0f), "density", nanovdb::GridClass::FogVolume);
			}

			template<typename T>
			VoxelNanoVDBManager<T>::~VoxelNanoVDBManager() {
				// NanoVDB grid is managed by shared_ptr, no manual cleanup needed
			}

			template<typename T>
			void VoxelNanoVDBManager<T>::init(unsigned int expected_voxels) {
				// Clear and re-initialize the NanoVDB grid
				using BuildType = GridBuildType_t<T>;
				nano_grid = std::make_shared<T>(floatToValue<BuildType>(0.0f), "density", nanovdb::GridClass::FogVolume);
				insert_count = 0;

				// Set transform if available
				if (transform_scale > 0.0) {
					nano_grid->setTransform(transform_scale, nanovdb::Vec3d(0, 0, 0));
				}
			}

			template<typename T>
			void VoxelNanoVDBManager<T>::insertOrUpdatePackedSequential(uint64_t key, float value) {
				using BuildType = GridBuildType_t<T>;
				int x, y, z;
				unpackCoord3(key, x, y, z);

				auto accessor = nano_grid->getAccessor();
				nanovdb::Coord coord(x, y, z);

				// NanoVDB GridBuilder supports accumulation: getValue + setValue
				auto current_value = accessor.getValue(coord);
				accessor.setValue(coord, current_value + floatToValue<BuildType>(value));
				insert_count++;
			}

			template<typename T>
			void VoxelNanoVDBManager<T>::insertOrUpdate(void* _voxels, size_t num_voxels) {
				if (num_voxels <= 0) return;

				using BuildType = GridBuildType_t<T>;
				Voxel* h_voxels = (Voxel*)_voxels;

				auto accessor = nano_grid->getAccessor();
				for (size_t idx = 0; idx < num_voxels; idx++) {
					const Voxel& v = h_voxels[idx];
					nanovdb::Coord coord(v.i, v.j, v.k);

					// NanoVDB GridBuilder supports accumulation: getValue + setValue
					auto current_value = accessor.getValue(coord);
					accessor.setValue(coord, current_value + floatToValue<BuildType>(v.value));
					insert_count++;
				}
			}

			//int VoxelNanoVDBManager::extractAll(Voxel** h_output_voxels) {
			//	auto start_total = omp_get_wtime();

			//	// Collect all active voxels from NanoVDB grid
			//	std::vector<Voxel> voxels;

			//	// NanoVDB accessor to iterate through active voxels
			//	auto accessor = nano_grid->getAccessor();

			//	// Note: GridBuilder doesn't provide a convenient iterator over active values
			//	// We need to iterate through the bounding box
			//	// This is less efficient than hash table extraction but necessary with NanoVDB

			//	// Get bounding box
			//	auto bbox = nano_grid->getBBox();
			//	if (bbox.empty()) {
			//		*h_output_voxels = nullptr;
			//		auto end_total = omp_get_wtime();
			//		printf("[VoxelNanoVDBManager::extractAll] Total: %.3f ms (empty grid)\n", (end_total - start_total) * 1000.0);
			//		return 0;
			//	}

			//	// Iterate through bounding box and collect non-zero values
			//	for (auto iter = bbox.begin(); iter; ++iter) {
			//		nanovdb::Coord coord = *iter;
			//		float value = accessor.getValue(coord);
			//		if (value != 0.0f) {
			//			voxels.push_back(Voxel(coord[0], coord[1], coord[2], value));
			//		}
			//	}

			//	// Allocate output array
			//	int total_count = static_cast<int>(voxels.size());
			//	*h_output_voxels = new Voxel[total_count];
			//	std::copy(voxels.begin(), voxels.end(), *h_output_voxels);

			//	auto end_total = omp_get_wtime();
			//	auto total_duration = (end_total - start_total) * 1000.0;
			//	printf("[VoxelNanoVDBManager::extractAll] Total: %.3f ms\n", total_duration);

			//	return total_count;
			//}

			template<typename T>
			void VoxelNanoVDBManager<T>::clear() {
				// Create a new empty grid
				using BuildType = GridBuildType_t<T>;
				nano_grid = std::make_shared<T>(floatToValue<BuildType>(0.0f), "density", nanovdb::GridClass::FogVolume);
				if (transform_scale > 0.0) {
					nano_grid->setTransform(transform_scale, nanovdb::Vec3d(0, 0, 0));
				}
				insert_count = 0;
			}

			template<typename T>
			void VoxelNanoVDBManager<T>::serialize(uint8_t* bin_data) {
				nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle = nanovdb::tools::createNanoGrid(*nano_grid);
				size_t buffer_size = grid_handle.bufferSize();
				memcpy(bin_data, grid_handle.data(), buffer_size);
			}

			template<typename T>
			void VoxelNanoVDBManager<T>::deserialize(uint8_t* bin_data) {
				// Deserialize from a finalized NanoGrid into our builder grid
				// First, clear the current grid
				using BuildType = GridBuildType_t<T>;
				nano_grid = std::make_shared<T>(floatToValue<BuildType>(0.0f), "density", nanovdb::GridClass::FogVolume);
				if (transform_scale > 0.0) {
					nano_grid->setTransform(transform_scale, nanovdb::Vec3d(0, 0, 0));
				}
				insert_count = 0;

				// Now populate it from the serialized data
				merge(bin_data);
			}

			template<typename T>
			void VoxelNanoVDBManager<T>::merge(uint8_t* bin_data) {
				using BuildType = GridBuildType_t<T>;
				auto acc_dst = nano_grid->getAccessor();
				auto* grid_src_float = (nanovdb::NanoGrid<float>*)bin_data;

				// loop over child nodes of the root node
				for (auto it2 = grid_src_float->tree().root().cbeginChild(); it2; ++it2) {
					// loop over child nodes of the upper internal node
					for (auto it1 = it2->cbeginChild(); it1; ++it1) {
						// loop over child nodes of the lower internal node
						for (auto it0 = it1->cbeginChild(); it0; ++it0) {
							// loop over values
							for (auto it = it0->cbeginValueOn(); it; ++it) {
								float f = *it;
								BuildType v = floatToValue<BuildType>(f);

								nanovdb::Coord xyz = it.getCoord();

								if (acc_dst.isValueOn(xyz)) {
									v = v + acc_dst.getValue(xyz); //ADD
								}
								acc_dst.setValue(xyz, v);
							}
						}
					}
				}

			}

			template<typename T>
			void VoxelNanoVDBManager<T>::merge(common::vdb::VoxelSparseManager* _other) {
				// Attempt to dynamic_cast to VoxelNanoVDBManager
				VoxelNanoVDBManager* other = dynamic_cast<VoxelNanoVDBManager*>(_other);
				if (!other) {
					printf("[VoxelNanoVDBManager::merge] Error: Incompatible manager type\n");
					return;
				}

				if (!other->nano_grid) {
					return;
				}

				//// Get accessor for destination grid
				//auto acc_dst = nano_grid->getAccessor();

				//// Iterate directly over active voxels in source grid using tree structure
				//// This avoids the overhead of extractAll() and intermediate allocations
				//auto bbox = other->nano_grid->getBBox();
				//if (bbox.empty()) {
				//	return;
				//}

				//auto acc_src = other->nano_grid->getAccessor();
				//int merged_count = 0;

				//// Iterate through bounding box and merge non-zero values
				//for (auto iter = bbox.begin(); iter; ++iter) {
				//	nanovdb::Coord xyz = *iter;
				//	float v = acc_src.getValue(xyz);

				//	if (v != 0.0f) {
				//		// Accumulate: add source value to destination value
				//		float current_value = acc_dst.getValue(xyz);
				//		acc_dst.setValue(xyz, current_value + v);
				//		merged_count++;
				//	}
				//}

				//insert_count += merged_count;

				auto acc_dst = nano_grid->getAccessor();
				nanovdb::GridHandle<nanovdb::HostBuffer> other_grid_handle = nanovdb::tools::createNanoGrid(*other->nano_grid);
				auto* grid_src_float = (nanovdb::NanoGrid<float>*)other_grid_handle.data();
				merge((uint8_t*)grid_src_float);
			}

			template<typename T>
			size_t VoxelNanoVDBManager<T>::mem_size() const {
				// TODO
				nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle = nanovdb::tools::createNanoGrid(*nano_grid);
				return grid_handle.bufferSize();
			}

			// Explicit template instantiations
			template class VoxelNanoVDBManager<nanovdb::tools::build::FloatGrid>;
			template class VoxelNanoVDBManager<nanovdb::tools::build::Vec3fGrid>;
#endif // WITH_NANOVDB

#ifdef WITH_OPENVDB
			// =====================================================================
			// VoxelOpenVDBManager Implementation
			// =====================================================================

			// Helper to get the ValueType from an OpenVDB Grid
			template<typename GridT>
			struct GridValueType;

			template<typename ValueT>
			struct GridValueType<openvdb::Grid<openvdb::tree::Tree4<ValueT, 5, 4, 3>::Type>> {
				using type = ValueT;
			};

			template<typename GridT>
			using GridValueType_t = typename GridValueType<GridT>::type;

			// Helper to convert float to the appropriate value type for OpenVDB
			template<typename ValueT>
			inline ValueT floatToValueVDB(float f) {
				return ValueT(f);
			}

			template<>
			inline float floatToValueVDB<float>(float f) {
				return f;
			}

			template<>
			inline openvdb::Vec3f floatToValueVDB<openvdb::Vec3f>(float f) {
				return openvdb::Vec3f(f, f, f);
			}

			// Helper to get default value for type
			template<typename ValueT>
			inline ValueT getDefaultValue() {
				return ValueT(0.0f);
			}

			template<>
			inline float getDefaultValue<float>() {
				return 0.0f;
			}

			template<>
			inline openvdb::Vec3f getDefaultValue<openvdb::Vec3f>() {
				return openvdb::Vec3f(0.0f, 0.0f, 0.0f);
			}

			template<typename T>
			VoxelOpenVDBManager<T>::VoxelOpenVDBManager() : insert_count(0), common::vdb::VoxelSparseManager() {
				// Initialize OpenVDB library (thread-safe, only initializes once)
				openvdb::initialize();

				// Create a new OpenVDB grid with the appropriate type
				using ValueType = typename T::ValueType;
				vdb_grid = std::make_shared<T>(getDefaultValue<ValueType>());
				vdb_grid->setGridClass(openvdb::GRID_FOG_VOLUME);
				vdb_grid->setName("density");
			}

			template<typename T>
			VoxelOpenVDBManager<T>::~VoxelOpenVDBManager() {
				// Grid is managed by shared_ptr, no manual cleanup needed
			}

			template<typename T>
			void VoxelOpenVDBManager<T>::init(unsigned int expected_voxels) {
				// Clear and re-initialize the OpenVDB grid
				using ValueType = typename T::ValueType;
				vdb_grid = std::make_shared<T>(getDefaultValue<ValueType>());
				vdb_grid->setGridClass(openvdb::GRID_FOG_VOLUME);
				vdb_grid->setName("density");
				insert_count = 0;

				// Set transform if available
				if (transform_scale > 0.0) {
					openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
					vdb_grid->setTransform(transform);
				}
			}

			template<typename T>
			void VoxelOpenVDBManager<T>::insertOrUpdatePackedSequential(uint64_t key, float value) {
				using ValueType = typename T::ValueType;
				int x, y, z;
				unpackCoord3(key, x, y, z);

				auto accessor = vdb_grid->getAccessor();
				openvdb::Coord coord(x, y, z);

				// OpenVDB supports accumulation: getValue + setValue
				auto current_value = accessor.getValue(coord);
				accessor.setValue(coord, current_value + floatToValueVDB<ValueType>(value));
				insert_count++;
			}

			template<typename T>
			void VoxelOpenVDBManager<T>::insertOrUpdate(void* _voxels, size_t num_voxels) {
				if (num_voxels <= 0) return;

				using ValueType = typename T::ValueType;
				Voxel* h_voxels = (Voxel*)_voxels;

				auto accessor = vdb_grid->getAccessor();
				for (size_t idx = 0; idx < num_voxels; idx++) {
					const Voxel& v = h_voxels[idx];
					openvdb::Coord coord(v.i, v.j, v.k);

					// OpenVDB supports accumulation: getValue + setValue
					auto current_value = accessor.getValue(coord);
					accessor.setValue(coord, current_value + floatToValueVDB<ValueType>(v.value));
					insert_count++;
				}
			}

			//int VoxelOpenVDBManager::extractAll(Voxel** h_output_voxels) {
			//	auto start_total = omp_get_wtime();

			//	// Collect all active voxels from OpenVDB grid
			//	std::vector<Voxel> voxels;

			//	// Use OpenVDB's value iterator to efficiently iterate over active voxels
			//	for (auto iter = vdb_grid->cbeginValueOn(); iter; ++iter) {
			//		openvdb::Coord coord = iter.getCoord();
			//		float value = *iter;
			//		if (value != 0.0f) {
			//			voxels.push_back(Voxel(coord[0], coord[1], coord[2], value));
			//		}
			//	}

			//	// Allocate output array
			//	int total_count = static_cast<int>(voxels.size());
			//	if (total_count > 0) {
			//		*h_output_voxels = new Voxel[total_count];
			//		std::copy(voxels.begin(), voxels.end(), *h_output_voxels);
			//	}
			//	else {
			//		*h_output_voxels = nullptr;
			//	}

			//	auto end_total = omp_get_wtime();
			//	auto total_duration = (end_total - start_total) * 1000.0;
			//	printf("[VoxelOpenVDBManager::extractAll] Total: %.3f ms\n", total_duration);

			//	return total_count;
			//}

			template<typename T>
			void VoxelOpenVDBManager<T>::clear() {
				// Create a new empty grid
				using ValueType = typename T::ValueType;
				vdb_grid = std::make_shared<T>(getDefaultValue<ValueType>());
				vdb_grid->setGridClass(openvdb::GRID_FOG_VOLUME);
				vdb_grid->setName("density");
				if (transform_scale > 0.0) {
					openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
					vdb_grid->setTransform(transform);
				}
				insert_count = 0;
			}

			template<typename T>
			void VoxelOpenVDBManager<T>::serialize(uint8_t* bin_data) {
				// Serialize OpenVDB grid to buffer
				// Create a grid container and write to stream
				openvdb::GridPtrVec grids;
				grids.push_back(vdb_grid);

				std::ostringstream ostr(std::ios_base::binary);
				openvdb::io::Stream(ostr).write(grids);
				std::string str = ostr.str();

				// Write size and data
				size_t size = str.size();
				memcpy(bin_data, &size, sizeof(size_t));
				memcpy(bin_data + sizeof(size_t), str.data(), size);
			}

			template<typename T>
			void VoxelOpenVDBManager<T>::deserialize(uint8_t* bin_data) {
				// Deserialize from a finalized OpenVDB grid
				// First, clear the current grid
				using ValueType = typename T::ValueType;
				vdb_grid = std::make_shared<T>(getDefaultValue<ValueType>());
				vdb_grid->setGridClass(openvdb::GRID_FOG_VOLUME);
				vdb_grid->setName("density");
				if (transform_scale > 0.0) {
					openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
					vdb_grid->setTransform(transform);
				}
				insert_count = 0;

				// Now populate it from the serialized data
				merge(bin_data);
			}

			template<typename T>
			void VoxelOpenVDBManager<T>::merge(uint8_t* bin_data) {
				// Read serialized grid and merge it
				size_t size;
				memcpy(&size, bin_data, sizeof(size_t));
				
				std::string str(reinterpret_cast<const char*>(bin_data + sizeof(size_t)), size);
				std::istringstream istr(str, std::ios_base::binary);
				
				openvdb::io::Stream stream(istr);
				openvdb::GridPtrVecPtr grids_src = stream.getGrids();
				
				if (grids_src && !grids_src->empty()) {
					typename T::Ptr grid_src = openvdb::gridPtrCast<T>((*grids_src)[0]);
					if (grid_src) {
						// Use OpenVDB's optimized composite operation for efficient parallel merging
						openvdb::tools::compSum(*vdb_grid, *grid_src);
					}
				}
			}

			template<typename T>
			void VoxelOpenVDBManager<T>::merge(common::vdb::VoxelSparseManager* _other) {
				// Attempt to dynamic_cast to VoxelOpenVDBManager
				VoxelOpenVDBManager* other = dynamic_cast<VoxelOpenVDBManager*>(_other);
				if (!other) {
					printf("[VoxelOpenVDBManager::merge] Error: Incompatible manager type\n");
					return;
				}

				if (!other->vdb_grid) {
					return;
				}

				// Use OpenVDB's optimized composite operation for efficient parallel merging
				// This is significantly faster than manual iteration
				openvdb::tools::compSum(*vdb_grid, *other->vdb_grid);

				insert_count += other->insert_count;
			}

			template<typename T>
			size_t VoxelOpenVDBManager<T>::mem_size() const {
				// Estimate memory size of the grid
				return vdb_grid->memUsage();
			}

			// Explicit template instantiations
			template class VoxelOpenVDBManager<openvdb::FloatGrid>;
			template class VoxelOpenVDBManager<openvdb::Vec3fGrid>;
#endif // WITH_OPENVDB
		}
	}
}