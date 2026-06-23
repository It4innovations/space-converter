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

#include "ipic3d_hdf5_extract_iolib.h"

#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include <hdf5.h>

#include "convert_common.h"

#define RETURN_NORM_EMPTY return 0;
#define RETURN_NORM_VALUE(v) return (float)(v);
#define RETURN_NORM_VECTOR3(v) return (float)space_converter::common::calculate_dmagnitude3(v[0], v[1], v[2]);

#define RETURN_COMP_EMPTY return 0;
#define RETURN_COMP_VALUE(v) return 1;
#define RETURN_COMP_VECTOR3(v) return 3;

#define RETURN_ORIG_EMPTY RETURN_COMP_EMPTY
#define RETURN_ORIG_VALUE(v) out_value[0] = (float)v; RETURN_COMP_VALUE(v)
#define RETURN_ORIG_VECTOR3(v) out_value[0] = (float)v[0]; out_value[1] = (float)v[1]; out_value[2] = (float)v[2]; RETURN_COMP_VECTOR3(v)

namespace ipic3d {
    namespace io {

        // Data structure to hold particle data
        struct ParticleData {
            std::vector<double> x, y, z;    // Position
            std::vector<double> u, v, w;    // Velocity
            std::vector<double> q;          // Charge
            std::vector<uint64_t> id;       // Particle ID
            size_t count;                   // Number of particles
            int species_type;               // Species type (0-3)
        };

        ParticleData particle_data;
        int64_t g_total_particles = 0;
        double steps_time[2];
        int current_species_id = 0;
        int current_cycle_id = 0;

        void print_CPU_steps() {
            printf("iPIC3D HDF5 reading time: %f seconds\n", steps_time[1] - steps_time[0]);
        }

        float get_particle_norm_value(int blocknr, uint64_t id) {
            if (id >= particle_data.count) RETURN_NORM_EMPTY;

            switch (blocknr) {
            case IPIC3DBlockType::Pos: {
                double pos[3] = { particle_data.x[id], particle_data.y[id], particle_data.z[id] };
                RETURN_NORM_VECTOR3(pos);
            }
            case IPIC3DBlockType::Vel: {
                double vel[3] = { particle_data.u[id], particle_data.v[id], particle_data.w[id] };
                RETURN_NORM_VECTOR3(vel);
            }
            case IPIC3DBlockType::Charge:
                RETURN_NORM_VALUE(particle_data.q[id]);
            case IPIC3DBlockType::ID:
                RETURN_NORM_VALUE(particle_data.id[id]);
            default:
                RETURN_NORM_EMPTY;
            }
        }

        int get_particle_value(int blocknr, uint64_t id, float* out_value) {
            if (id >= particle_data.count) RETURN_ORIG_EMPTY;

            switch (blocknr) {
            case IPIC3DBlockType::Pos: {
                double pos[3] = { particle_data.x[id], particle_data.y[id], particle_data.z[id] };
                RETURN_ORIG_VECTOR3(pos);
            }
            case IPIC3DBlockType::Vel: {
                double vel[3] = { particle_data.u[id], particle_data.v[id], particle_data.w[id] };
                RETURN_ORIG_VECTOR3(vel);
            }
            case IPIC3DBlockType::Charge:
                RETURN_ORIG_VALUE(particle_data.q[id]);
            case IPIC3DBlockType::ID:
                RETURN_ORIG_VALUE(particle_data.id[id]);
            default:
                RETURN_ORIG_EMPTY;
            }
        }

        int get_particle_value_comp(int blocknr, uint64_t id) {
            if (id >= particle_data.count) RETURN_COMP_EMPTY;

            switch (blocknr) {
            case IPIC3DBlockType::Pos:
            case IPIC3DBlockType::Vel:
                RETURN_COMP_VECTOR3(nullptr);
            case IPIC3DBlockType::Charge:
            case IPIC3DBlockType::ID:
                RETURN_COMP_VALUE(nullptr);
            default:
                RETURN_COMP_EMPTY;
            }
        }

        int get_particle_type(uint64_t id) {
            return particle_data.species_type;
        }

        void get_particle_position(uint64_t id, double* pos) {
            if (id < particle_data.count) {
                pos[0] = particle_data.x[id];
                pos[1] = particle_data.y[id];
                pos[2] = particle_data.z[id];
            }
        }

        size_t get_local_num_particles() {
            return particle_data.count;
        }

        size_t get_global_num_particles() {
            return g_total_particles;
        }

        double get_particle_hsml(uint64_t id) {
            // iPIC3D doesn't have smoothing length, return a default value
            return 1.0;
        }

        double get_particle_mass(uint64_t id) {
            // Use charge as a proxy for mass in iPIC3D
            if (id < particle_data.count) {
                return particle_data.q[id];
            }
            return 1.0;
        }

        double get_particle_rho(uint64_t id) {
            // Density not directly available, return default
            return 1.0;
        }

        int get_particle_rho_blocknr() {
            return 0; // No density block
        }

        std::string get_dataset_name(int blocknr) {
            switch (blocknr) {
            case IPIC3DBlockType::Pos:
                return "Position";
            case IPIC3DBlockType::Vel:
                return "Velocity";
            case IPIC3DBlockType::Charge:
                return "Charge";
            case IPIC3DBlockType::ID:
                return "ID";
            default:
                return "Unknown";
            }
        }

        void get_types_and_blocks(std::vector<int>& types_and_blocks) {
            types_and_blocks.resize(IPIC3DParticleType::PTMax * IPIC3DBlockType::BTMax, 0);
            
            // Mark available blocks for the current species
            types_and_blocks[IPIC3DBlockType::BTMax * IPIC3DBlockType::Pos + particle_data.species_type] = 1;
            types_and_blocks[IPIC3DBlockType::BTMax * IPIC3DBlockType::Vel + particle_data.species_type] = 1;
            types_and_blocks[IPIC3DBlockType::BTMax * IPIC3DBlockType::Charge + particle_data.species_type] = 1;
            types_and_blocks[IPIC3DBlockType::BTMax * IPIC3DBlockType::ID + particle_data.species_type] = 1;
        }

        void print_types_and_blocks_local() {
            std::vector<int> types_and_blocks;
            get_types_and_blocks(types_and_blocks);

            printf("\n");
            for (int type = 0; type < IPIC3DParticleType::PTMax; type++) {
                printf("Type: Species_%d (%d)\n", type, type);
                for (int blocknr = 0; blocknr < IPIC3DBlockType::BTMax; blocknr++) {
                    if (types_and_blocks[IPIC3DBlockType::BTMax * blocknr + type] > 0) {
                        printf("\t%s (%d)\n", get_dataset_name(blocknr).c_str(), blocknr);
                    }
                }
            }
        }

        void print_types_and_blocks(std::vector<int>& types_and_blocks) {
            printf("\nAll snapshots contain:\n");
            for (int type = 0; type < IPIC3DParticleType::PTMax; type++) {
                printf("Type: Species_%d (%d)\n", type, type);
                for (int blocknr = 0; blocknr < IPIC3DBlockType::BTMax; blocknr++) {
                    if (types_and_blocks[IPIC3DBlockType::BTMax * blocknr + type] > 0) {
                        printf("\t%s (%d)\n", get_dataset_name(blocknr).c_str(), blocknr);
                    }
                }
            }
        }

        // Replace "{}" placeholders in pattern with the given file number.
        std::string format_filename(const std::string& pattern, int number) {
            std::string result = pattern;
            const std::string placeholder = "{}";
            size_t pos = result.find(placeholder);
            while (pos != std::string::npos) {
                std::string number_str = std::to_string(number);
                result.replace(pos, placeholder.length(), number_str);
                pos = result.find(placeholder, pos + number_str.length());
            }
            return result;
        }

        // Split `total` items into `parts` contiguous chunks, returning the start
        // index and size of chunk `idx` (remainder distributed to the first chunks).
        void partition_range(size_t total, size_t parts, size_t idx, size_t& start, size_t& count) {
            size_t base = total / parts;
            size_t rem = total % parts;
            count = base + (idx < rem ? 1 : 0);
            start = idx * base + std::min(idx, rem);
        }

        // Append a [start, start + count) slice of src's fields onto dst's fields.
        void append_particle_range(const ParticleData& src, ParticleData& dst, size_t start, size_t count) {
            dst.x.insert(dst.x.end(), src.x.begin() + start, src.x.begin() + start + count);
            dst.y.insert(dst.y.end(), src.y.begin() + start, src.y.begin() + start + count);
            dst.z.insert(dst.z.end(), src.z.begin() + start, src.z.begin() + start + count);
            dst.u.insert(dst.u.end(), src.u.begin() + start, src.u.begin() + start + count);
            dst.v.insert(dst.v.end(), src.v.begin() + start, src.v.begin() + start + count);
            dst.w.insert(dst.w.end(), src.w.begin() + start, src.w.begin() + start + count);
            dst.q.insert(dst.q.end(), src.q.begin() + start, src.q.begin() + start + count);
            dst.id.insert(dst.id.end(), src.id.begin() + start, src.id.begin() + start + count);
        }

        // Helper function to read a dataset from HDF5
        template<typename T>
        void read_hdf5_dataset(hid_t file_id, const std::string& dataset_path, std::vector<T>& data) {
            hid_t dataset_id = H5Dopen(file_id, dataset_path.c_str(), H5P_DEFAULT);
            if (dataset_id < 0) {
                throw std::runtime_error("Failed to open dataset: " + dataset_path);
            }

            hid_t dataspace_id = H5Dget_space(dataset_id);
            hssize_t num_elements = H5Sget_simple_extent_npoints(dataspace_id);
            
            data.resize(num_elements);
            
            hid_t memtype;
            if (std::is_same<T, double>::value) {
                memtype = H5T_NATIVE_DOUBLE;
            } else if (std::is_same<T, float>::value) {
                memtype = H5T_NATIVE_FLOAT;
            } else if (std::is_same<T, uint64_t>::value) {
                memtype = H5T_NATIVE_UINT64;
            } else if (std::is_same<T, int64_t>::value) {
                memtype = H5T_NATIVE_INT64;
            } else {
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                throw std::runtime_error("Unsupported data type for HDF5 reading");
            }
            
            herr_t status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
            
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            
            if (status < 0) {
                throw std::runtime_error("Failed to read dataset: " + dataset_path);
            }
        }

        // Open one HDF5 file (file_index substituted into the "{}" placeholder of
        // hdf5_file_pattern) and read all particle datasets for species/cycle into out.
        void read_particle_file(const std::string& hdf5_file_pattern, int file_index,
                                 const std::string& species_path, const std::string& cycle_suffix,
                                 int species_id, ParticleData& out) {
            std::string file_path = format_filename(hdf5_file_pattern, file_index);
            hid_t file_id = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (file_id < 0) {
                throw std::runtime_error("Failed to open HDF5 file: " + file_path);
            }

            try {
                out.species_type = species_id;
                read_hdf5_dataset(file_id, species_path + "/x" + cycle_suffix, out.x);
                read_hdf5_dataset(file_id, species_path + "/y" + cycle_suffix, out.y);
                read_hdf5_dataset(file_id, species_path + "/z" + cycle_suffix, out.z);
                read_hdf5_dataset(file_id, species_path + "/u" + cycle_suffix, out.u);
                read_hdf5_dataset(file_id, species_path + "/v" + cycle_suffix, out.v);
                read_hdf5_dataset(file_id, species_path + "/w" + cycle_suffix, out.w);
                read_hdf5_dataset(file_id, species_path + "/q" + cycle_suffix, out.q);
                read_hdf5_dataset(file_id, species_path + "/ID" + cycle_suffix, out.id);
                out.count = out.x.size();
            } catch (const std::exception& e) {
                H5Fclose(file_id);
                throw;
            }

            H5Fclose(file_id);
        }

        void init_lib(std::string hdf5_file, int world_rank, int world_size,
                      int species_id, int cycle_id, int num_files) {
#ifdef WITH_OPENMP
            steps_time[0] = omp_get_wtime();
#endif
            current_species_id = species_id;
            current_cycle_id = cycle_id;
            particle_data.species_type = species_id;

            std::string species_path = "/particles/species_" + std::to_string(species_id);
            std::string cycle_suffix = "/cycle_" + std::to_string(cycle_id);

            if (world_rank == 0) {
                printf("Reading iPIC3D HDF5 file(s): %s\n", hdf5_file.c_str());
                printf("Species: %d, Cycle: %d, num_files: %d, world_size: %d\n", species_id, cycle_id, num_files, world_size);
            }

            ParticleData result;
            result.species_type = species_id;
            result.count = 0;

            if (num_files >= world_size) {
                // Each rank owns one or more whole files exclusively; concatenate them.
                size_t file_start, file_count;
                partition_range((size_t)num_files, (size_t)world_size, (size_t)world_rank, file_start, file_count);

                for (size_t i = 0; i < file_count; i++) {
                    ParticleData file_data;
                    read_particle_file(hdf5_file, (int)(file_start + i), species_path, cycle_suffix, species_id, file_data);
                    append_particle_range(file_data, result, 0, file_data.count);
                }
            } else {
                // Fewer files than ranks: a group of ranks shares each file and
                // splits its particles among the group.
                size_t rank_start = 0, rank_count = 0;
                int file_index = -1;
                size_t local_index = 0, group_size = 1;

                for (int f = 0; f < num_files; f++) {
                    partition_range((size_t)world_size, (size_t)num_files, (size_t)f, rank_start, rank_count);
                    if ((size_t)world_rank >= rank_start && (size_t)world_rank < rank_start + rank_count) {
                        file_index = f;
                        local_index = (size_t)world_rank - rank_start;
                        group_size = rank_count;
                        break;
                    }
                }

                ParticleData file_data;
                read_particle_file(hdf5_file, file_index, species_path, cycle_suffix, species_id, file_data);

                size_t local_start, local_count;
                partition_range(file_data.count, group_size, local_index, local_start, local_count);
                append_particle_range(file_data, result, local_start, local_count);
            }

            result.count = result.x.size();
            particle_data = result;
            g_total_particles = particle_data.count;

            printf("Rank %d has %lld particles\n", world_rank, (long long)particle_data.count);

#ifdef WITH_OPENMP
            steps_time[1] = omp_get_wtime();
#endif
        }

        void finish_lib() {
            particle_data.x.clear();
            particle_data.y.clear();
            particle_data.z.clear();
            particle_data.u.clear();
            particle_data.v.clear();
            particle_data.w.clear();
            particle_data.q.clear();
            particle_data.id.clear();
            particle_data.count = 0;
        }

    } // namespace io
} // namespace ipic3d
