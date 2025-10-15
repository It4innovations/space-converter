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

#include <cstdint>
#include <string>
#include <vector>

 // Namespace changa contains functionality for handling Tipsy format particle data in Changa simulations.
namespace changa {
    // Namespace tipsy defines particle and block types, along with I/O utilities for working with Tipsy data.
    namespace tipsy {

        // Enum defining particle types in the Tipsy format.
        enum ParticleType {
            Gas = 0,        // Gas particles
            Dark,           // Dark matter particles
            Stars,          // Star particles

            PTMax           // Maximum value for ParticleType (used for validation or iteration)
        };

        // Enum defining block types, representing various attributes of Tipsy particles.
        enum BlockType {
            Pos = 0,        // Position
            Mass,           // Mass
            Vel,            // Velocity
            Soft,           // Softening length
            Phi,            // Gravitational potential
            HSmooth,        // Smoothing length
            Rho,            // Density
            Temp,           // Temperature
            Metals,         // Metallicity
            TForm,          // Formation time

            BTMax           // Maximum value for BlockType (used for validation or iteration)
        };

        // Namespace io provides utility functions for I/O operations on Tipsy particle data.
        namespace io {

            // Print CPU-specific processing steps for Tipsy operations.
            void print_CPU_steps();

            // Retrieve a float value from a specific block and particle.
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

            // Retrieve the radius of a particle.
            // @param id: The unique identifier of the particle.
            // @return The radius of the particle.
            //double get_particle_radius(uint64_t id);

            // Retrieve the smoothing length (HSML) of a particle.
            // @param id: The unique identifier of the particle.
            // @return The smoothing length of the particle.
            double get_particle_hsml(uint64_t id);

            // Retrieve the mass of a particle.
            // @param id: The unique identifier of the particle.
            // @return The mass of the particle.
            double get_particle_mass(uint64_t id);

            // Retrieve the density (rho) of a particle.
            // @param id: The unique identifier of the particle.
            // @return The density of the particle.
            double get_particle_rho(uint64_t id);
            int get_particle_rho_blocknr();

            // Retrieve the unit information of a specific block.
            // @param blocknr: The block number.
            // @return A string representing the unit of the block.
            std::string get_particle_unit(int blocknr);

            // Initialize the Tipsy library for reading particle data.
            // @param basedir: The base directory containing the Tipsy files.
            // @param world_rank: The rank of the current process in the communicator.
            // @param world_size: The total number of processes in the communicator.
            void init_lib(std::string basedir, int world_rank, int world_size);

            // Finalize and clean up the Tipsy library.
            void finish_lib();

            // Retrieve the types and blocks associated with the particle data.
            // @param types_and_blocks: A vector to store the retrieved types and blocks.
            void get_types_and_blocks(std::vector<int>& types_and_blocks);

            // Retrieve the names of blocks.
            // @param blocknr: The block number.
            // @return The name of block.
            std::string get_dataset_name(int blocknr);

        } // namespace io

    } // namespace tipsy

} // namespace changa
