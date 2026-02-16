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
		/**
		 * @brief Structure to hold raw particle data for serialization and transmission
		 * 
		 * RawParticles stores particle data in a compact format suitable for serialization,
		 * network transmission, and file I/O operations.
		 */
		struct RawParticles {
			/**
			 * @brief Single particle data attribute (position, velocity, etc.)
			 */
			struct ParticleData {
				std::string name;              ///< Name of the data attribute (e.g., "position", "velocity")
				int num_comp = 0;              ///< Number of components per particle (3 for position, 1 for scalar)
				std::vector<float> values;     ///< Flat array of values (size = num_particles * num_comp)
			};

			std::vector<ParticleData> data;  ///< Collection of all particle data attributes

			/**
			 * @brief Serialize particle data to a binary file
			 * @param filename Path to output file
			 */
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

			/**
			 * @brief Deserialize particle data from a binary file
			 * @param filename Path to input file
			 */
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

			/**
			 * @brief Serialize particle data to a binary buffer (for MPI/network transfer)
			 * @param bin_data Output buffer to store serialized data
			 */
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
			/**
			 * @brief Deserialize particle data from a binary buffer
			 * @param bin_data Input buffer containing serialized data
			 */
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

			/**
			 * @brief Merge another RawParticles object into this one
			 * @param other Source particle data to merge
			 * 
			 * Appends values for matching particle data names, or adds new data attributes.
			 */
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

		/**
		 * @brief Dense regular grid representation of particle data
		 * 
		 * Stores particle data rasterized onto a uniform 3D grid.
		 * Used for volume rendering and analysis.
		 */
		struct DenseParticles
		{
			std::vector<float> data_density;   ///< Primary density/value data
#ifndef WITH_NO_DATA_TEMP
			std::vector<float> data_temp;      ///< Temporary accumulation buffer for normalization
#endif
			size_t dims[3] = { 0,0,0 };        ///< Grid dimensions [x, y, z]
			size_t offset[3] = { 0,0,0 };      ///< Grid offset in global coordinate space

			/**
			 * @brief Clear all grid data and reset dimensions
			 */
			void clear() {
				data_density.clear();
#ifndef WITH_NO_DATA_TEMP
				data_temp.clear();
#endif
				memset(dims, 0, 3 * sizeof(size_t));
				memset(offset, 0, 3 * sizeof(size_t));
			}

			/**
			 * @brief Create and initialize grid with specified dimensions
			 * @param x Width of the grid
			 * @param y Height of the grid
			 * @param z Depth of the grid
			 */
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

			/** @brief Get grid width */
			size_t x() {
				return dims[0];
			}

			/** @brief Get grid height */
			size_t y() {
				return dims[1];
			}

			/** @brief Get grid depth */
			size_t z() {
				return dims[2];
			}

			/** @brief Get total number of voxels */
			size_t size() {
				return dims[0] * dims[1] * dims[2];
			}

			/** @brief Get total memory size in bytes */
			size_t memsize() {
				return dims[0] * dims[1] * dims[2] * sizeof(float);
			}

			/**
			 * @brief Convert 3D coordinates to linear index
			 * @param x X coordinate
			 * @param y Y coordinate
			 * @param z Z coordinate
			 * @return Linear array index
			 */
			size_t get_index(size_t x, size_t y, size_t z) {
				return x + y * dims[0] + z * dims[0] * dims[1];
			}
		};

		/**
		 * @brief Container for different VDB grid representations
		 * 
		 * Can hold particle data in various formats: dense grid, sparse NanoVDB,
		 * OpenVDB, serialized binary, or raw particle data.
		 */
		class VDBParticles
		{
		public:
			DenseParticles dense_grid;  ///< Dense regular grid representation
#if OPENVDB_VERSION == 11
			std::shared_ptr<nanovdb::build::FloatGrid> nano_grid;  ///< NanoVDB sparse grid (v11)
#else
			std::shared_ptr<nanovdb::tools::build::FloatGrid> nano_grid;  ///< NanoVDB sparse grid
#endif

#ifdef WITH_OPENVDB
			openvdb::FloatGrid::Ptr vdb_grid;  ///< OpenVDB sparse grid
#endif
			std::vector<uint8_t> vector_grid;  ///< Serialized binary grid data (for MPI transfer)

			RawParticles raw_particles;  ///< Raw particle point cloud data

			/**
			 * @brief Type of VDB particle representation
			 */
			enum VDBParticleType
			{
				eDense,        ///< Dense regular grid
				eVector,       ///< Serialized binary format
				eNanoVDB,      ///< NanoVDB sparse grid
				eOpenVDB,      ///< OpenVDB sparse grid
				eRawParticles  ///< Raw point cloud
			};

			VDBParticleType type;  ///< Current representation type
		};

		/**
		 * @brief Base class for converting particle data to VDB grids
		 * 
		 * Provides functionality to convert various particle formats into VDB grid
		 * representations (dense, sparse, or raw). Handles spatial indexing, filtering,
		 * and normalization. Derived classes implement particle I/O library specifics.
		 */
		class ConvertVDBBase {
		public:

			/**
			 * @brief Convert particle data from I/O library format to VDB grid
			 * 
			 * Main conversion function that reads particles and rasterizes them into the
			 * specified grid type. Handles filtering, normalization, and spatial indexing.
			 * 
			 * @param particle_type Type of particles to convert
			 * @param particle_fix_size Fixed particle size multiplier
			 * @param grid_name Name for the output grid
			 * @param grid_transform Grid transformation scale
			 * @param bbox_min Bounding box minimum coordinates
			 * @param bbox_max Bounding box maximum coordinates
			 * @param bbox_dim Grid dimension (resolution)
			 * @param dense_type Type of density calculation
			 * @param grid Output VDB grid container
			 */
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

			/**
			 * @brief Merge one VDB grid into another
			 * 
			 * Combines grid data from grid_recv into grid_dst, handling different
			 * grid types (dense, sparse, serialized).
			 */
			/**
			 * @brief Merge one VDB grid into another
			 * @param grid_dst Destination grid to merge into
			 * @param grid_recv Source grid to merge from
			 * 
			 * Combines grid data from grid_recv into grid_dst, handling different
			 * grid types (dense, sparse, serialized).
			 */
			void merge_grid(
				VDBParticles& grid_dst,
				VDBParticles& grid_recv
			);

			/**
			 * @brief Fill voxels with particle contribution using advanced splatting
			 * @param grid Dense grid to fill
			 * @param pid Particle ID
			 * @param value Particle value to splat
			 * @param bbox_dim Bounding box dimension
			 * @param bbox_min_orig Bounding box minimum coordinates
			 * @param bbox_size_orig Original bounding box size
			 * @param scale_space_diagonal Diagonal scaling factor
			 * @param dense_type Type of density calculation
			 * @param dense_norm Normalization type
			 * @param particle_fix_size Fixed particle size multiplier
			 * @param particle_type Type of particle
			 * @param block_name_id Data block identifier
			 * @param pos Particle position
			 */
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
			/**
			 * @brief Convert dense grid to NanoVDB format (OpenVDB v11)
			 */
			std::shared_ptr<nanovdb::build::FloatGrid> dense_to_nanovdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);
#else
			/**
			 * @brief Convert dense grid to NanoVDB format
			 */
			std::shared_ptr<nanovdb::tools::build::FloatGrid> dense_to_nanovdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);
#endif

#endif


#ifdef WITH_OPENVDB
			/**
			 * @brief Convert dense grid to OpenVDB format
			 */
			/**
			 * @brief Convert dense grid to OpenVDB format
			 * @param particles Dense particle grid
			 * @param transform_scale Grid transformation scale
			 * @param dense_type Type of density calculation
			 * @param dense_norm Normalization type
			 * @return OpenVDB grid pointer
			 */
			openvdb::FloatGrid::Ptr dense_to_openvdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm);

			/**
			 * @brief Serialize OpenVDB grid to binary buffer
			 * @param grid OpenVDB grid to serialize
			 * @param file_content Output buffer for serialized data
			 */
			void openvdb_to_vector(openvdb::FloatGrid::Ptr grid, std::vector<uint8_t>& file_content);
			
			/**
			 * @brief Serialize two OpenVDB grids to binary buffer
			 */
			void openvdb_to_vector2(openvdb::FloatGrid::Ptr grid1, openvdb::FloatGrid::Ptr grid2, std::vector<uint8_t>& file_content);
			
			/**
			 * @brief Deserialize OpenVDB grid from binary buffer
			 */
			openvdb::FloatGrid::Ptr vector_to_openvdb(std::vector<uint8_t>& file_content);
#endif

			// Particle radius data organized by particle type
			std::vector< std::vector<float> > radius_particles_per_ptype;  ///< Particle radii per type
			std::vector<size_t> particles_ptype_offset;                    ///< Particle count offsets per type
			std::vector< std::vector<float> > rho_particles_per_ptype;     ///< Density values per type

			// Cosmological simulation parameters
			double redshift = 0.0;              ///< Cosmological redshift value
			double hubble_param = 1.0;          ///< Hubble parameter (h)
			double radius_particle_const = 0.0; ///< Constant particle radius value

#ifdef WITH_EMBREE
			void* rtc_device;  ///< Embree ray tracing device
			void* rtc_scene;   ///< Embree ray tracing scene
			
			/**
			 * @brief Create Embree ray tracing scene for particle intersection
			 * @param particle_type Type of particles to include in scene
			 */
			void create_embree_scene(int particle_type);
#endif

			/**
			 * @brief Read particle radius data from file
			 * @param calc_radius_neigh_file Path to radius data file
			 */
			void read_radius_from_file(std::string &calc_radius_neigh_file);
			
			/**
			 * @brief Write particle radius data to file
			 * @param calc_radius_neigh_file Path to output file
			 */
			void write_radius_from_file(std::string& calc_radius_neigh_file);

#ifdef WITH_CUDAKDTREE
			/**
			 * @brief Calculate particle radii using CUDA KD-Tree neighbor search
			 * @param calc_radius_neigh Number of neighbors to consider
			 * @param calc_radius_neigh_file Output file for radius data
			 * @param use_cycling Enable cycling boundary conditions
			 * @param use_cudakdtree_cpu Use CPU fallback instead of GPU
			 * @param maxRadius Maximum search radius
			 * @param rho_kernel Density kernel type for SPH calculations
			 */
			void calculate_radius_by_cudakdtree(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, bool use_cudakdtree_cpu, float maxRadius, common::SpaceData::DenseType& rho_kernel);
#endif

#ifdef WITH_NANOFLANN
			/**
			 * @brief Calculate particle radii using Nanoflann KD-Tree neighbor search
			 * @param calc_radius_neigh Number of neighbors to consider
			 * @param calc_radius_neigh_file Output file for radius data
			 * @param use_cycling Enable cycling boundary conditions
			 * @param rho_kernel Density kernel type for SPH calculations
			 */
			void calculate_radius_by_nanoflann(int calc_radius_neigh, std::string &calc_radius_neigh_file, bool use_cycling, common::SpaceData::DenseType& rho_kernel);
#endif

			/**
			 * @brief Get effective radius for a particle
			 * @param pid Particle ID
			 * @param bbox_dim Bounding box dimension
			 * @param bbox_min_orig Original bounding box minimum
			 * @param bbox_size_orig Original bounding box size
			 * @param scale_space_diagonal Diagonal scaling factor
			 * @param particle_fix_size Fixed size multiplier
			 * @param particle_type Type of particle
			 * @return Particle radius in grid coordinates
			 */
			virtual double get_particle_radius(
				uint64_t pid,
				int bbox_dim,
				int* bbox_min_orig,
				double bbox_size_orig,
				double scale_space_diagonal,
				float particle_fix_size,
				int particle_type
			);

			/**
			 * @brief Find bounding box for particles of given type
			 * @param particle_type Type of particles to consider
			 * @param bbox_min Output minimum coordinates
			 * @param bbox_max Output maximum coordinates
			 * @param offset_position Optional position offset
			 */
			virtual void iolib_find_bbox(
				int particle_type,
				float* bbox_min,
				float* bbox_max,
				float* offset_position
			);

			/**
			 * @brief Find minimum and maximum values for data block
			 * @param particle_type Type of particles
			 * @param block_nr Data block number
			 * @param v_min Output minimum value
			 * @param v_max Output maximum value
			 */
			virtual void iolib_find_minmax(
				int particle_type,
				int block_nr,
				float& v_min,
				float& v_max
			);

			/**
			 * @brief Get density value for particle
			 * @param id Particle ID
			 * @return Density value
			 */
			virtual double get_particle_rho(uint64_t id);

			/**
			 * @brief Get available particle types and data blocks
			 * @param types_and_blocks Output vector of type/block pairs
			 */
			virtual void get_types_and_blocks(std::vector<int>& types_and_blocks);
			
			/**
			 * @brief Get normalized value for particle in data block
			 * @param blocknr Data block number
			 * @param id Particle ID
			 * @return Normalized value
			 */
			virtual float get_particle_norm_value(int blocknr, uint64_t id);
			
			/**
			 * @brief Get value for particle in data block
			 * @param blocknr Data block number
			 * @param id Particle ID
			 * @param out_value Output buffer for value(s)
			 * @return Number of components
			 */
			virtual int get_particle_value(int blocknr, uint64_t id, float* out_value);
			
			/**
			 * @brief Get number of components for particle value
			 * @param blocknr Data block number
			 * @param id Particle ID
			 * @return Number of components
			 */
			virtual int get_particle_value_comp(int blocknr, uint64_t id);

			/**
			 * @brief Format filename with number substitution
			 * @param pattern Filename pattern with format specifiers
			 * @param number Number to substitute into pattern
			 * @return Formatted filename
			 */
			std::string format_filename(const std::string& pattern, int number);

		public:
			// ============================================================
			// Pure Virtual Methods - Must be implemented by derived classes
			// ============================================================
			
			/** @brief Print CPU processing steps for debugging */
			virtual void print_CPU_steps() = 0;
			
			/** @brief Get normalized particle value (internal implementation) */
			virtual float get_particle_norm_value_internal(int blocknr, uint64_t id) = 0;
			
			/** @brief Get particle value (internal implementation) */
			virtual int get_particle_value_internal(int blocknr, uint64_t id, float* out_value) = 0;
			
			/** @brief Get particle value component count (internal implementation) */
			virtual int get_particle_value_comp_internal(int blocknr, uint64_t id) = 0;
			
			/** @brief Get particle type ID */
			virtual int get_particle_type(uint64_t id) = 0;
			
			/** @brief Get particle position coordinates */
			virtual void get_particle_position(uint64_t id, double* pos) const = 0;
			
			/** @brief Get number of particles on this MPI rank */
			virtual size_t get_local_num_particles() const = 0;
			
			/** @brief Get total number of particles across all ranks */
			virtual size_t get_global_num_particles() const = 0;

			/** @brief Get particle smoothing length (SPH) */
			virtual double get_particle_hsml(uint64_t id) = 0;
			
			/** @brief Get particle mass */
			virtual double get_particle_mass(uint64_t id) = 0;
			
			/** @brief Get particle density (internal implementation) */
			virtual double get_particle_rho_internal(uint64_t id) = 0;
			
			/** @brief Get data block number containing density values */
			virtual int get_particle_rho_blocknr() = 0;

			/** @brief Initialize I/O library with MPI parameters */
			virtual void init_lib(int argc, char** argv, int world_rank, int world_size) = 0;
			
			/** @brief Finalize and cleanup I/O library */
			virtual void finish_lib() = 0;

			/** @brief Get particle types and data blocks (internal implementation) */
			virtual void get_types_and_blocks_internal(std::vector<int>& types_and_blocks) = 0;
			
			/** @brief Print local particle types and blocks */
			virtual void print_types_and_blocks_local() = 0;
			
			/** @brief Print particle types and blocks across all ranks */
			virtual void print_types_and_blocks(std::vector<int>& types_and_blocks) = 0;

			/** @brief Get human-readable name for particle type */
			virtual std::string get_type_name(int type) = 0;
			
			/** @brief Get human-readable name for data block */
			virtual std::string get_dataset_name(int blocknr) = 0;

			/** @brief Get formatted string of all particle data type names */
			virtual std::string get_particle_data_type_names(std::vector<int>& types_and_blocks) = 0;

			/** @brief Get number of particle types */
			virtual int get_num_types() = 0;
			
			/** @brief Get number of data blocks */
			virtual int get_num_blocks() = 0;
		};
	}// vdb
} //common