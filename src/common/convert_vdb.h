/*
 * Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
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

#include <map>
#include <string>
#include <vector>

#include "convert_common.h"
#include "data_common.h"

#ifdef WITH_NANOVDB

#ifdef WITH_TBB
#define NANOVDB_USE_TBB
#define NANOVDB_USE_INTRINSICS
#endif

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/GridBuilder.h>
#	include <nanovdb/util/CreateNanoGrid.h>
#	include <nanovdb/util/IO.h>
#else
#	include <nanovdb/tools/GridBuilder.h>
#	include <nanovdb/tools/CreateNanoGrid.h>
#	include <nanovdb/io/IO.h>
#endif

#endif

#ifdef WITH_OPENVDB
#	include <openvdb/openvdb.h>
#	include <openvdb/points/PointDataGrid.h>
#endif

namespace common {
	namespace vdb {
		struct RawParticles {
			struct ParticleData {
				std::string name;
				int num_comp = 0;
				std::vector<float> values;

				//// Assignment operator for ParticleData
				//ParticleData& operator=(const ParticleData& other) {
				//	if (this != &other) {
				//		name = other.name;
				//		num_comp = other.num_comp;
				//		values = other.values;
				//	}
				//	return *this;
				//}

				//// Move assignment operator for ParticleData
				//ParticleData& operator=(ParticleData&& other) noexcept {
				//	if (this != &other) {
				//		name = std::move(other.name);
				//		num_comp = other.num_comp;
				//		values = std::move(other.values);
				//	}
				//	return *this;
				//}
			};

			std::vector<ParticleData> data;

			//// Assignment operator for RawParticles
			//RawParticles& operator=(const RawParticles& other) {
			//	if (this != &other) {
			//		data = other.data; // std::vector has its own assignment operator
			//	}
			//	return *this;
			//}

			//// Move assignment operator for RawParticles
			//RawParticles& operator=(RawParticles&& other) noexcept {
			//	if (this != &other) {
			//		data = std::move(other.data);
			//	}
			//	return *this;
			//}

			void serialize(const std::string& filename) const {
				std::ofstream out(filename, std::ios::binary);
				if (!out) {
					throw std::runtime_error("Failed to open file for writing");
				}

				// Write the size of the data vector
				size_t dataSize = data.size();
				out.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));

				// Write each ParticleData
				for (const auto& particle : data) {
					// Write the size of the name string
					size_t nameSize = particle.name.size();
					out.write(reinterpret_cast<const char*>(&nameSize), sizeof(nameSize));

					// Write the characters of the name string
					out.write(particle.name.data(), nameSize);

					// Write num_comp
					out.write(reinterpret_cast<const char*>(&particle.num_comp), sizeof(particle.num_comp));

					// Write the size of the values vector
					size_t valuesSize = particle.values.size();
					out.write(reinterpret_cast<const char*>(&valuesSize), sizeof(valuesSize));

					// Write the values
					out.write(reinterpret_cast<const char*>(particle.values.data()), valuesSize * sizeof(float));
				}

				out.close();
			}

			void deserialize(const std::string& filename) {
				std::ifstream in(filename, std::ios::binary);
				if (!in) {
					throw std::runtime_error("Failed to open file for reading");
				}

				// Read the size of the data vector
				size_t dataSize;
				in.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));

				data.resize(dataSize);

				// Read each ParticleData
				for (auto& particle : data) {
					// Read the size of the name string
					size_t nameSize;
					in.read(reinterpret_cast<char*>(&nameSize), sizeof(nameSize));

					// Read the characters of the name string
					particle.name.resize(nameSize);
					in.read(&particle.name[0], nameSize);

					// Read num_comp
					in.read(reinterpret_cast<char*>(&particle.num_comp), sizeof(particle.num_comp));

					// Read the size of the values vector
					size_t valuesSize;
					in.read(reinterpret_cast<char*>(&valuesSize), sizeof(valuesSize));

					// Read the values
					particle.values.resize(valuesSize);
					in.read(reinterpret_cast<char*>(particle.values.data()), valuesSize * sizeof(float));
				}

				in.close();
			}

			void serialize(std::vector<uint8_t>& bin_data) const {
				std::ostringstream oss(std::ios::binary);

				// Write the size of the data vector
				size_t dataSize = data.size();
				oss.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));

				// Write each ParticleData
				for (const auto& particle : data) {
					// Write the size of the name string
					size_t nameSize = particle.name.size();
					oss.write(reinterpret_cast<const char*>(&nameSize), sizeof(nameSize));

					// Write the characters of the name string
					oss.write(particle.name.data(), nameSize);

					// Write num_comp
					oss.write(reinterpret_cast<const char*>(&particle.num_comp), sizeof(particle.num_comp));

					// Write the size of the values vector
					size_t valuesSize = particle.values.size();
					oss.write(reinterpret_cast<const char*>(&valuesSize), sizeof(valuesSize));

					// Write the values
					oss.write(reinterpret_cast<const char*>(particle.values.data()), valuesSize * sizeof(float));
				}

				const std::string& str = oss.str();
				bin_data.assign(str.begin(), str.end());
			}

			void deserialize(const std::vector<uint8_t>& bin_data) {
				std::istringstream iss(std::string(bin_data.begin(), bin_data.end()), std::ios::binary);

				// Read the size of the data vector
				size_t dataSize;
				iss.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));

				data.resize(dataSize);

				// Read each ParticleData
				for (auto& particle : data) {
					// Read the size of the name string
					size_t nameSize;
					iss.read(reinterpret_cast<char*>(&nameSize), sizeof(nameSize));

					// Read the characters of the name string
					particle.name.resize(nameSize);
					iss.read(&particle.name[0], nameSize);

					// Read num_comp
					iss.read(reinterpret_cast<char*>(&particle.num_comp), sizeof(particle.num_comp));

					// Read the size of the values vector
					size_t valuesSize;
					iss.read(reinterpret_cast<char*>(&valuesSize), sizeof(valuesSize));

					// Read the values
					particle.values.resize(valuesSize);
					iss.read(reinterpret_cast<char*>(particle.values.data()), valuesSize * sizeof(float));
				}
			}

			void merge(const RawParticles& other) {
				for (const auto& otherParticle : other.data) {
					auto it = std::find_if(data.begin(), data.end(), [&](const ParticleData& particle) {
						return particle.name == otherParticle.name;
						});

					if (it != data.end()) {
						// If a particle with the same name exists, append values
						it->values.insert(it->values.end(), otherParticle.values.begin(), otherParticle.values.end());
					}
					else {
						data.push_back(otherParticle);
					}
				}
			}
		};

		struct DenseParticles
		{
			std::vector<float> data_density;
#ifndef WITH_NO_DATA_TEMP
			std::vector<float> data_temp;
#endif			
			size_t dims[3] = { 0,0,0 };
			//int type = 0;
			size_t offset[3] = { 0,0,0 };

			void clear() {
				data_density.clear();
#ifndef WITH_NO_DATA_TEMP				
				data_temp.clear();
#endif				

				memset(dims, 0, 3 * sizeof(size_t));
				memset(offset, 0, 3 * sizeof(size_t));
			}

			void create(size_t x, size_t y, size_t z) {
				dims[0] = x;
				dims[1] = y;
				dims[2] = z;

				data_density.resize(size());
				memset(data_density.data(), 0, memsize());
#ifndef WITH_NO_DATA_TEMP
				data_temp.resize(size());
				memset(data_temp.data(), 0, memsize());
#endif				
			}

			size_t x() {
				return dims[0];
			}

			size_t y() {
				return dims[1];
			}

			size_t z() {
				return dims[2];
			}

			size_t size() {
				return dims[0] * dims[1] * dims[2];
			}

			size_t memsize() {
				return dims[0] * dims[1] * dims[2] * sizeof(float);
			}

			size_t get_index(size_t x, size_t y, size_t z) {
				return x + y * dims[0] + z * dims[0] * dims[1];
			}
		};

		class VDBParticles
		{
		public:
			DenseParticles dense_grid;
#if OPENVDB_VERSION == 11
			std::shared_ptr<nanovdb::build::FloatGrid> nano_grid;
#else
			std::shared_ptr<nanovdb::tools::build::FloatGrid> nano_grid;
#endif

#ifdef WITH_OPENVDB
			openvdb::FloatGrid::Ptr vdb_grid;
#endif
			std::vector<uint8_t> vector_grid;

			RawParticles raw_particles;

			enum VDBParticleType
			{
				eDense,
				eVector,
				eNanoVDB,
				eOpenVDB,
				eRawParticles
			};

			VDBParticleType type;
		};

		class ConvertVDBBase {
		public:

			void convert_iolib_to_grid(
				int particle_type,
				float particle_fix_size,
				std::string grid_name,
				float grid_transform,
				float* bbox_min,
				float* bbox_max,
				int bbox_dim,
				int* bbox_min_orig,
				double bbox_size_orig,
				common::SpaceData::ExtractedType extracted_type,
				common::SpaceData::DenseType dense_type,
				common::SpaceData::DenseNorm dense_norm,
				int block_name_id,
				float object_size,
				float& min_value,
				float& max_value,
				float min_value_global,
				float max_value_global,
				size_t& particles_count,
				VDBParticles& grid,
				double& transform_scale,
				float filter_min,
				float filter_max,
				float min_rho,
				float max_rho,
				common::SpaceData::AnimType anim_type,
				int frame_req,
				int frame,
				float *bbox_sphere_pos,
				float bbox_sphere_r,
				bool use_simple_density,
				float *offset_position
			);

			void merge_grid(
				VDBParticles& grid_dst,
				VDBParticles& grid_recv	
			);		
#if 0
			void fill_voxels_v1(common::vdb::DenseParticles& grid,
				size_t pid, int px, int py, int pz, float v,
				int bbox_dim, int* bbox_min_orig, int* bbox_max_orig, double scale_space_diagonal);

			void fill_voxels_v2(common::vdb::DenseParticles& grid,
				size_t pid, int px, int py, int pz, float value,
				int bbox_dim, int* bbox_min_orig, int* bbox_max_orig, double scale_space_diagonal, float particle_fix_size);

			void fill_voxels_v3(common::vdb::DenseParticles& grid,
				size_t pid, float value,
				int bbox_dim, int* bbox_min_orig, int* bbox_max_orig,
				double scale_space_diagonal,
				int dense_type, float particle_fix_size);

			void fill_voxels_v4(
				common::vdb::DenseParticles& grid,
				size_t pid, 
				float value,
				int bbox_dim, 
				int* bbox_min_orig, 
				int* bbox_max_orig,
				double scale_space_diagonal,
				common::SpaceData::DenseType dense_type, 
				float particle_fix_size,
				int particle_type,
				int block_name_id,
				double *pos		
				);
#endif

			void fill_voxels_v5(
				common::vdb::DenseParticles& grid,
				size_t pid, 
				float value,
				int bbox_dim, 
				int* bbox_min_orig, 
				double bbox_size_orig,
				double scale_space_diagonal,
				common::SpaceData::DenseType dense_type, 
				common::SpaceData::DenseNorm dense_norm,
				float particle_fix_size,
				int particle_type,
				int block_name_id,
				double *pos		
				);				

#ifdef WITH_NANOVDB
			

#if OPENVDB_VERSION == 11
			std::shared_ptr<nanovdb::build::FloatGrid> dense_to_nanovdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);
#else
			std::shared_ptr<nanovdb::tools::build::FloatGrid> dense_to_nanovdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);
#endif

#endif


#ifdef WITH_OPENVDB
			openvdb::FloatGrid::Ptr dense_to_openvdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);
			//openvdb::FloatGrid::Ptr dense_to_openvdbE(DenseParticles& particles, double transform_scale);
			//openvdb::FloatGrid::Ptr dense_to_openvdbI(DenseParticles& particles, double transform_scale);

			void openvdb_to_vector(openvdb::FloatGrid::Ptr grid, std::vector<uint8_t>& file_content);
			void openvdb_to_vector2(openvdb::FloatGrid::Ptr grid1, openvdb::FloatGrid::Ptr grid2, std::vector<uint8_t>& file_content);
			openvdb::FloatGrid::Ptr vector_to_openvdb(std::vector<uint8_t>& file_content);
#endif

			std::vector< std::vector<float> > radius_particles_per_ptype;
			std::vector<size_t> particles_ptype_offset;
			std::vector< std::vector<float> > rho_particles_per_ptype;

			double redshift = 0.0;
			double hubble_param = 1.0;

#ifdef WITH_EMBREE
			void* rtc_device;
			void* rtc_scene;
			void create_embree_scene(int particle_type);
#endif

			void read_radius_from_file(std::string &calc_radius_neigh_file);
			void write_radius_from_file(std::string& calc_radius_neigh_file);

#ifdef WITH_CUDAKDTREE
			void calculate_radius_by_cudakdtree(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, bool use_cudakdtree_cpu, float maxRadius, common::SpaceData::DenseType& rho_kernel);
#endif

#ifdef WITH_NANOFLANN
			void calculate_radius_by_nanoflann(int calc_radius_neigh, std::string &calc_radius_neigh_file, bool use_cycling, common::SpaceData::DenseType& rho_kernel);
#endif

			virtual double get_particle_radius(
				uint64_t pid,
				int bbox_dim,
				int* bbox_min_orig,
				double bbox_size_orig,
				double scale_space_diagonal,
				float particle_fix_size,
				int particle_type
			);

			virtual void iolib_find_bbox(
				int particle_type,
				float* bbox_min,
				float* bbox_max,
				float* offset_position
			);

			virtual void iolib_find_minmax(
				int particle_type,
				int block_nr,
				float& v_min,
				float& v_max
			);

			virtual double get_particle_rho(uint64_t id);

			virtual void get_types_and_blocks(std::vector<int>& types_and_blocks);
			virtual float get_particle_norm_value(int blocknr, uint64_t id);
			virtual int get_particle_value(int blocknr, uint64_t id, float* out_value);
			virtual int get_particle_value_comp(int blocknr, uint64_t id);

		public:
			// Abstract methods
			virtual void print_CPU_steps() = 0;
			virtual float get_particle_norm_value_internal(int blocknr, uint64_t id) = 0;
			virtual int get_particle_value_internal(int blocknr, uint64_t id, float* out_value) = 0;
			virtual int get_particle_value_comp_internal(int blocknr, uint64_t id) = 0;
			virtual int get_particle_type(uint64_t id) = 0;
			virtual void get_particle_position(uint64_t id, double* pos) const = 0;
			virtual size_t get_local_num_particles() const = 0;
			virtual size_t get_global_num_particles() const = 0;

			virtual double get_particle_hsml(uint64_t id) = 0;
			virtual double get_particle_mass(uint64_t id) = 0;
			virtual double get_particle_rho_internal(uint64_t id) = 0;
			virtual int get_particle_rho_blocknr() = 0;

			virtual void init_lib(int argc, char** argv, int world_rank, int world_size) = 0;
			virtual void finish_lib() = 0;

			virtual void get_types_and_blocks_internal(std::vector<int>& types_and_blocks) = 0;
			virtual void print_types_and_blocks_local() = 0;
			virtual void print_types_and_blocks(std::vector<int>& types_and_blocks) = 0;

			virtual std::string get_type_name(int type) = 0;
			virtual std::string get_dataset_name(int blocknr) = 0;

			virtual std::string get_particle_data_type_names(std::vector<int>& types_and_blocks) = 0;

			virtual int get_num_types() = 0;
			virtual int get_num_blocks() = 0;
		};
	}// vdb
} //common