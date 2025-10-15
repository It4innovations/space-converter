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

#include "convert_vdb.h"

 // Namespace csv contains functionality for handling VDB data conversion specific to the CSV format (mainly for testing).
namespace csv {

    // ConvertVDBCSV is a derived class from common::vdb::ConvertVDBBase.
    // This class implements methods for converting particle data from CSV/Excel format into VDB format.
    class ConvertVDBCSV : public common::vdb::ConvertVDBBase {
    public:
        // Print the steps executed on the CPU during the conversion process.
        void print_CPU_steps() override;

        // Retrieve a float value associated with a particle in a specific block.
        // @param blocknr: The block number.
        // @param id: The unique identifier of the particle.
        // @return A float value specific to the particle.
        float get_particle_norm_value_internal(int blocknr, uint64_t id) override;

        int get_particle_value_internal(int blocknr, uint64_t id, float* out_value) override;
        int get_particle_value_comp_internal(int blocknr, uint64_t id) override;

        // Get the type of a particle based on its unique identifier.
        // @param id: The unique identifier of the particle.
        // @return An integer representing the particle type.
        int get_particle_type(uint64_t id) override;

        // Retrieve the position of a particle in 3D space.
        // @param id: The unique identifier of the particle.
        // @param pos: A pointer to an array to store the particle's position (x, y, z).
        void get_particle_position(uint64_t id, double* pos) const override;

        // Get the number of particles on the local process.
        // @return The local particle count.
        size_t get_local_num_particles() const override;

        // Get the total number of particles across all processes.
        // @return The global particle count.
        size_t get_global_num_particles() const override;

        // Retrieve the radius of a particle.
        // @param id: The unique identifier of the particle.
        // @return The radius of the particle.
        //double get_particle_radius(uint64_t id) override;

        // Retrieve the smoothing length (HSML) of a particle.
        // @param id: The unique identifier of the particle.
        // @return The smoothing length of the particle.
        double get_particle_hsml(uint64_t id) override;

        // Retrieve the mass of a particle.
        // @param id: The unique identifier of the particle.
        // @return The mass of the particle.
        double get_particle_mass(uint64_t id) override;

        // Retrieve the density (rho) of a particle.
        // @param id: The unique identifier of the particle.
        // @return The density of the particle.
        double get_particle_rho_internal(uint64_t id) override;
        int get_particle_rho_blocknr() override;

        // Initialize the library with command-line arguments and parallel processing configuration.
        // @param argc: The number of command-line arguments.
        // @param argv: The array of command-line arguments.
        // @param world_rank: The rank of the process in the communicator.
        // @param world_size: The total number of processes in the communicator.
        void init_lib(int argc, char** argv, int world_rank, int world_size) override;

        // Finalize and clean up the library resources.
        void finish_lib() override;

        // Retrieve the particle types and blocks for processing.
        // @param types_and_blocks: A vector to store the particle types and block information.
        void get_types_and_blocks_internal(std::vector<int>& types_and_blocks) override;

        // Print particle types and blocks information for the local process.
        void print_types_and_blocks_local() override;

        // Print the particle types and blocks information globally.
        // @param types_and_blocks: A vector containing types and block information.
        void print_types_and_blocks(std::vector<int>& types_and_blocks) override;

        // Get the name of a particle type.
        // @param type: The integer representation of the type.
        // @return The name of the type as a string.
        std::string get_type_name(int type) override;

        // Get the name of a dataset block.
        // @param blocknr: The block number.
        // @return The name of the dataset as a string.
        std::string get_dataset_name(int blocknr) override;

        // Retrieve a string description of the particle data types based on types and blocks.
        // @param types_and_blocks: A vector containing the particle types and block information.
        // @return A string representing the data types of the particles.
        std::string get_particle_data_type_names(std::vector<int>& types_and_blocks) override;

        int get_num_types() override;
        int get_num_blocks() override;

    }; // class ConvertVDBCSV

} // namespace csv
