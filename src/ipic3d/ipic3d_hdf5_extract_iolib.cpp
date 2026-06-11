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

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include <omp.h>
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

        void init_lib(std::string hdf5_file, int world_rank, int world_size, 
                      int species_id, int cycle_id) {
            steps_time[0] = omp_get_wtime();

            current_species_id = species_id;
            current_cycle_id = cycle_id;
            particle_data.species_type = species_id;

            if (world_rank == 0) {
                printf("Reading iPIC3D HDF5 file: %s\n", hdf5_file.c_str());
                printf("Species: %d, Cycle: %d\n", species_id, cycle_id);
            }

            // Open HDF5 file
            hid_t file_id = H5Fopen(hdf5_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (file_id < 0) {
                throw std::runtime_error("Failed to open HDF5 file: " + hdf5_file);
            }

            try {
                // Construct dataset paths
                std::string species_path = "/particles/species_" + std::to_string(species_id);
                std::string cycle_suffix = "/cycle_" + std::to_string(cycle_id);

                // Read particle data
                read_hdf5_dataset(file_id, species_path + "/x" + cycle_suffix, particle_data.x);
                read_hdf5_dataset(file_id, species_path + "/y" + cycle_suffix, particle_data.y);
                read_hdf5_dataset(file_id, species_path + "/z" + cycle_suffix, particle_data.z);
                read_hdf5_dataset(file_id, species_path + "/u" + cycle_suffix, particle_data.u);
                read_hdf5_dataset(file_id, species_path + "/v" + cycle_suffix, particle_data.v);
                read_hdf5_dataset(file_id, species_path + "/w" + cycle_suffix, particle_data.w);
                read_hdf5_dataset(file_id, species_path + "/q" + cycle_suffix, particle_data.q);
                read_hdf5_dataset(file_id, species_path + "/ID" + cycle_suffix, particle_data.id);

                particle_data.count = particle_data.x.size();
                
                // For simplicity, assume each rank gets equal particles
                // In a real implementation, you'd want proper domain decomposition
                size_t particles_per_rank = particle_data.count / world_size;
                size_t start_idx = world_rank * particles_per_rank;
                size_t end_idx = (world_rank == world_size - 1) ? particle_data.count : (world_rank + 1) * particles_per_rank;
                size_t local_count = end_idx - start_idx;

                // Extract local portion
                if (world_size > 1) {
                    ParticleData local_data;
                    local_data.species_type = species_id;
                    local_data.count = local_count;
                    local_data.x.assign(particle_data.x.begin() + start_idx, particle_data.x.begin() + end_idx);
                    local_data.y.assign(particle_data.y.begin() + start_idx, particle_data.y.begin() + end_idx);
                    local_data.z.assign(particle_data.z.begin() + start_idx, particle_data.z.begin() + end_idx);
                    local_data.u.assign(particle_data.u.begin() + start_idx, particle_data.u.begin() + end_idx);
                    local_data.v.assign(particle_data.v.begin() + start_idx, particle_data.v.begin() + end_idx);
                    local_data.w.assign(particle_data.w.begin() + start_idx, particle_data.w.begin() + end_idx);
                    local_data.q.assign(particle_data.q.begin() + start_idx, particle_data.q.begin() + end_idx);
                    local_data.id.assign(particle_data.id.begin() + start_idx, particle_data.id.begin() + end_idx);
                    
                    g_total_particles = particle_data.count;
                    particle_data = local_data;
                } else {
                    g_total_particles = particle_data.count;
                }

                if (world_rank == 0) {
                    printf("Total particles read: %lld\n", (long long)g_total_particles);
                }
                printf("Rank %d has %lld particles\n", world_rank, (long long)particle_data.count);

            } catch (const std::exception& e) {
                H5Fclose(file_id);
                throw;
            }

            H5Fclose(file_id);
            steps_time[1] = omp_get_wtime();
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
