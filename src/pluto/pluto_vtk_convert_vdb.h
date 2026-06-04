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

#include "convert_vdb.h"

 // Namespace plutovtk contains functionality for handling VDB data conversion specific to PLUTO VTK format.
namespace plutovtk {

    // ConvertVDBPlutoVTK is a derived class from common::vdb::ConvertVDBBase.
    // This class implements methods for converting VTK structured grid data from PLUTO simulations into VDB format.
    class ConvertVDBPlutoVTK : public common::vdb::ConvertVDBBase {
    public:
        // Print the steps executed on the CPU during the conversion process.
        void print_CPU_steps() override;

        // Retrieve a float value associated with a voxel in a specific block.
        // @param blocknr: The block number.
        // @param id: The unique identifier of the voxel.
        // @return A float value specific to the voxel.
        float get_particle_norm_value_internal(int blocknr, uint64_t id) override;

        int get_particle_value_internal(int blocknr, uint64_t id, float* out_value) override;
        int get_particle_value_comp_internal(int blocknr, uint64_t id) override;

        // Get the type of a voxel based on its unique identifier.
        // @param id: The unique identifier of the voxel.
        // @return An integer representing the voxel type.
        int get_particle_type(uint64_t id) override;

        // Retrieve the position of a voxel in 3D space.
        // @param id: The unique identifier of the voxel.
        // @param pos: A pointer to an array to store the voxel's position (x, y, z).
        void get_particle_position(uint64_t id, double* pos) const override;

        // Get the number of voxels on the local process.
        // @return The local voxel count.
        size_t get_local_num_particles() const override;

        // Get the total number of voxels across all processes.
        // @return The global voxel count.
        size_t get_global_num_particles() const override;

        // Retrieve the smoothing length (HSML) of a voxel (cell size).
        // @param id: The unique identifier of the voxel.
        // @return The smoothing length of the voxel.
        double get_particle_hsml(uint64_t id) override;

        // Retrieve the mass of a voxel.
        // @param id: The unique identifier of the voxel.
        // @return The mass of the voxel.
        double get_particle_mass(uint64_t id) override;

        // Retrieve the density (rho) of a voxel.
        // @param id: The unique identifier of the voxel.
        // @return The density of the voxel.
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

        // Retrieve the voxel types and blocks for processing.
        // @param types_and_blocks: A vector to store the voxel types and block information.
        void get_types_and_blocks_internal(std::vector<int>& types_and_blocks) override;

        // Print voxel types and blocks information for the local process.
        void print_types_and_blocks_local() override;

        // Print the voxel types and blocks information globally.
        // @param types_and_blocks: A vector containing types and block information.
        void print_types_and_blocks(std::vector<int>& types_and_blocks) override;

        // Get the name of a voxel type.
        // @param type: The integer representation of the type.
        // @return The name of the type as a string.
        std::string get_type_name(int type) override;

        // Get the name of a dataset block.
        // @param blocknr: The block number.
        // @return The name of the dataset as a string.
        std::string get_dataset_name(int blocknr) override;

        // Retrieve a string description of the voxel data types based on types and blocks.
        // @param types_and_blocks: A vector containing the voxel types and block information.
        // @return A string representing the data types of the voxels.
        std::string get_particle_data_type_names(std::vector<int>& types_and_blocks) override;

        int get_num_types() override;
        int get_num_blocks() override;

    }; // class ConvertVDBPlutoVTK

} // namespace plutovtk
