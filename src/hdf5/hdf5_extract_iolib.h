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

#include <string>
#include <vector>

 // Namespace hdf5 provides functionalities for handling HDF5 data.
namespace hdf5 {

    // Enum defining particle types in the Tipsy format.
    enum HDF5ParticleType {
        Gas = 0,        // Gas particles
        Dark,           // Dark matter particles
        Stars,          // Star particles

        PTMax           // Maximum value for ParticleType (used for validation or iteration)
    };

    // Enum defining block types, representing various attributes of Tipsy particles.
    enum HDF5BlockType {
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

    // Namespace io contains I/O operations and utilities for HDF5 data processing.
    namespace io {

        // Print the steps executed on the CPU during the HDF5 process.
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

        // Initialize the HDF5 library.
        void init_lib(std::string basefile, int world_rank, int world_size);

        // Finalize and clean up the HDF5 library.
        void finish_lib();

        // Set parameters for the HDF5.
        // @param ICFormat: Initial condition file format.
        // @param SnapFormat: Snapshot file format.
        // @param NumFilesWrittenInParallel: Number of files written in parallel.
        // @param MaxMemSize: Maximum memory size allowed.
        // @param BufferSize: Size of the buffer.
        // @param PartAllocFactor: Particle allocation factor.
        // @param TotBHs: Count of BHs.
        //void hdf5_set_parameter(int ICFormat, int SnapFormat, int NumFilesWrittenInParallel, int MaxMemSize, double BufferSize, double PartAllocFactor, int TotBHs);

        // Read a parameter file for HDF5.
        // @param fname: Name of the parameter file.
        // @param tag: Array of parameter tags.
        // @param addr: Array of addresses for the parameters.
        // @param id: Array of parameter IDs.
        // @param nt: Number of tags.
        //void hdf5_read_parameter_file(char* fname, char* tag[], void** addr, int* id, int nt);

        // Read initial condition data from a file.
        // @param fname: Name of the file containing initial conditions.
        //void hdf5_read_ic(std::string fname);

        // Initialize memory allocation for HDF5.
        //void hdf5_mymalloc_init();

        // Set units for HDF5 data.
        //void hdf5_set_units();

        // Retrieve the particle types and blocks for processing.
        // @param types_and_blocks: A vector to store the particle types and block information.
        void get_types_and_blocks(std::vector<int>& types_and_blocks);

        // Print particle types and blocks information for the local process.
        void print_types_and_blocks_local();

        // Print the particle types and blocks information globally.
        // @param types_and_blocks: A vector containing types and block information.
        void print_types_and_blocks(std::vector<int>& types_and_blocks);

        // Get the name of a particle type.
        // @param type: The integer representation of the type.
        // @return The name of the type as a string.
        //std::string get_type_name(int type);

        // Get the name of a dataset block.
        // @param blocknr: The block number.
        // @return The name of the dataset as a string.
        std::string get_dataset_name(int blocknr);

        // Retrieve the radius of a particle.
        // @param id: The unique identifier of the particle.
        // @return The radius of the particle.
        double get_particle_radius(uint64_t id);

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

        // Retrieve the unit information of a specific block.
        // @param blocknr: The block number.
        // @return A string representing the unit of the block.
        std::string get_particle_unit(int blocknr);

    } // namespace io

} // namespace hdf5
