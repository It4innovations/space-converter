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

#include <cstdint>
#include <string>
#include <vector>

 // Namespace ipic3d provides functionalities for handling iPIC3D HDF5 data.
namespace ipic3d {

    // Enum defining particle types in the iPIC3D format.
    enum IPIC3DParticleType {
        Species_0 = 0,
        Species_1 = 1,
        Species_2 = 2,
        Species_3 = 3,

        PTMax = 4           // Maximum value for ParticleType (used for validation or iteration)
    };

    // Enum defining block types, representing various attributes of iPIC3D particles.
    enum IPIC3DBlockType {
        Pos = 0,        // Position (x, y, z)
        Vel = 1,        // Velocity (u, v, w)
        Charge = 2,     // Charge (q)
        ID = 3,         // Particle ID

        BTMax           // Maximum value for BlockType (used for validation or iteration)
    };

    // Namespace io contains I/O operations and utilities for iPIC3D HDF5 data processing.
    namespace io {

        // Print the steps executed on the CPU during the iPIC3D process.
        void print_CPU_steps();

        // Retrieve a float value associated with a particle in a specific block.
        // @param blocknr: The block number.
        // @param id: The unique identifier of the particle.
        // @return A float value specific to the particle and block.
        float get_particle_norm_value(int blocknr, uint64_t id);

        int get_particle_value(int blocknr, uint64_t id, float* out_value);
        int get_particle_value_comp(int blocknr, uint64_t id);

        // Get the type of a particle based on its unique identifier.
        // @param id: The unique identifier of the particle.
        // @return An integer representing the particle type.
        int get_particle_type(uint64_t id);

        // Retrieve the position of a particle in 3D space.
        // @param id: The unique identifier of the particle.
        // @param pos: A pointer to an array to store the particle's position (x, y, z).
        void get_particle_position(uint64_t id, double* pos);

        // Get the number of particles on the local process.
        // @return The local particle count.
        size_t get_local_num_particles();

        // Get the total number of particles across all processes.
        // @return The global particle count.
        size_t get_global_num_particles();

        // Initialize the iPIC3D library.
        void init_lib(std::string hdf5_file, int world_rank, int world_size,
            int species_id, int cycle_id, int num_files);

        // Finalize and clean up the iPIC3D library.
        void finish_lib();

        // Retrieve the particle types and blocks for processing.
        // @param types_and_blocks: A vector to store the particle types and block information.
        void get_types_and_blocks(std::vector<int>& types_and_blocks);

        // Print particle types and blocks information for the local process.
        void print_types_and_blocks_local();

        // Print the particle types and blocks information globally.
        // @param types_and_blocks: A vector containing types and block information.
        void print_types_and_blocks(std::vector<int>& types_and_blocks);

        // Get the name of a dataset block.
        // @param blocknr: The block number.
        // @return The name of the dataset as a string.
        std::string get_dataset_name(int blocknr);

        // Retrieve the smoothing length (HSML) of a particle.
        // @param id: The unique identifier of the particle.
        // @return The smoothing length of the particle (default value for iPIC3D).
        double get_particle_hsml(uint64_t id);

        // Retrieve the mass of a particle.
        // @param id: The unique identifier of the particle.
        // @return The mass of the particle (derived from charge for iPIC3D).
        double get_particle_mass(uint64_t id);

        // Retrieve the density (rho) of a particle.
        // @param id: The unique identifier of the particle.
        // @return The density of the particle.
        double get_particle_rho(uint64_t id);

        // Get the block number for density data.
        // @return The block number for density.
        int get_particle_rho_blocknr();

    } // namespace io
} // namespace ipic3d
