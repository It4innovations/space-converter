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

 // Namespace plutovtk provides functionalities for handling VTK data from PLUTO simulations.
namespace plutovtk {

    // Enum defining particle types (for uniformity with other formats, treat grid cells as particles).
    enum PlutoVTKParticleType {
        PLUTO = 0,

        PTMax           // Maximum value for ParticleType (used for validation or iteration)
    };

    // Enum defining block types, representing various attributes of the structured grid cells.
    enum PlutoVTKBlockType {
        Rho = 0,        // Density scalar field

        BTMax           // Maximum value for BlockType (used for validation or iteration)
    };

    // Namespace io contains I/O operations and utilities for VTK data processing.
    namespace io {

        // Print the steps executed on the CPU during the VTK process.
        void print_CPU_steps();

        // Retrieve a float value associated with a voxel in a specific block.
        // @param blocknr: The block number.
        // @param id: The unique identifier of the voxel (cell).
        // @return A float value specific to the voxel and block.
        float get_particle_norm_value(int blocknr, uint64_t id);

        int get_particle_value(int blocknr, uint64_t id, float* out_value);
        int get_particle_value_comp(int blocknr, uint64_t id);

        // Get the type of a voxel based on its unique identifier.
        // @param id: The unique identifier of the voxel.
        // @return An integer representing the voxel type.
        int get_particle_type(uint64_t id);

        // Retrieve the position of a voxel in 3D space.
        // @param id: The unique identifier of the voxel (cell).
        // @param pos: A pointer to an array to store the voxel's position (x, y, z).
        void get_particle_position(uint64_t id, double* pos);

        // Get the number of voxels on the local process.
        // @return The local voxel count.
        size_t get_local_num_particles();

        // Get the total number of voxels across all processes.
        // @return The global voxel count.
        size_t get_global_num_particles();

        // Initialize the VTK library.
        void init_lib(std::string vtk_file, int world_rank, int world_size, 
            std::vector<std::string> &scalar_names_vec);

        // Finalize and clean up the VTK library.
        void finish_lib();

        // Retrieve the voxel types and blocks for processing.
        // @param types_and_blocks: A vector to store the voxel types and block information.
        void get_types_and_blocks(std::vector<int>& types_and_blocks);

        // Print voxel types and blocks information for the local process.
        void print_types_and_blocks_local();

        // Print the voxel types and blocks information globally.
        // @param types_and_blocks: A vector containing types and block information.
        void print_types_and_blocks(std::vector<int>& types_and_blocks);

        // Get the name of a dataset block.
        // @param blocknr: The block number.
        // @return The name of the dataset as a string.
        std::string get_dataset_name(int blocknr);

        // Retrieve the smoothing length (HSML) of a voxel (cell size).
        // @param id: The unique identifier of the voxel.
        // @return The smoothing length of the voxel.
        double get_particle_hsml(uint64_t id);

        // Retrieve the mass of a voxel (using density and cell volume).
        // @param id: The unique identifier of the voxel.
        // @return The mass of the voxel.
        double get_particle_mass(uint64_t id);

        // Retrieve the density (rho) of a voxel.
        // @param id: The unique identifier of the voxel.
        // @return The density of the voxel.
        double get_particle_rho(uint64_t id);
        int get_particle_rho_blocknr();

        // Retrieve the unit information of a specific block.
        // @param blocknr: The block number.
        // @return A string representing the unit of the block.
        std::string get_particle_unit(int blocknr);

    } // namespace io

} // namespace plutovtk
