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

 // Namespace ramses provides functionalities for handling RAMSES AMR outputs
 // (output_NNNNN directories: info, amr, hydro and part files).
namespace ramses {

    // Particle types. Type 0 are the leaf cells of the AMR grid (treated as
    // particles, like the grid readers of the other formats); the others are
    // the RAMSES particle families.
    enum RamsesParticleType {
        Gas = 0,        // AMR leaf cells
        DM = 1,         // family 1
        Star = 2,       // family 2
        Cloud = 3,      // family 3 (sink cloud particles)
        Debris = 4,     // family 4
        Other = 5,      // family 5, undefined (127) and gas tracers (0)
        Sink = 6,       // sink particles (sink_NNNNN.csv, read by rank 0)

        PTMax           // Maximum value for ParticleType (used for validation or iteration)
    };

    // Fixed blocks; the variables of hydro_file_descriptor.txt (gas) and the
    // remaining fields of part_file_descriptor.txt (particles) follow from BTMax on.
    enum RamsesBlockType {
        Pos = 0,        // positions only (weight 1)
        Mass = 1,       // gas: density x cell volume, particles: mass
        Rho = 2,        // gas: density, particles: none
        Vel = 3,        // velocity vector (gas: velocity_x/y/z, particles: velocity_x/y/z)
        Level = 4,      // AMR level of the cell / of the grid the particle is attached to

        BTMax           // Maximum value for BlockType (used for validation or iteration)
    };

    // Namespace io contains I/O operations and utilities for RAMSES data processing.
    namespace io {

        // Print the steps executed on the CPU during the reading process.
        void print_CPU_steps();

        // Scalar ("norm") value of a particle in a block (vectors: magnitude).
        float get_particle_norm_value(int blocknr, uint64_t id);

        // Original components of a particle in a block; returns their count.
        int get_particle_value(int blocknr, uint64_t id, float* out_value);

        // Component count of a particle in a block (0 = not available).
        int get_particle_value_comp(int blocknr, uint64_t id);

        // Type of a particle (RamsesParticleType).
        int get_particle_type(uint64_t id);

        // Position of a particle / cell centre in code length units ([0, boxlen]).
        void get_particle_position(uint64_t id, double* pos);

        // Number of particles + leaf cells read by this rank.
        size_t get_local_num_particles();

        // Number of particles + leaf cells over all ranks.
        size_t get_global_num_particles();

        // Read-time selection: keep only the cells / particles whose field
        // `name` (a hydro variable or a particle field, e.g. "density" or
        // "birth_time") lies in [min, max]. Types without that field are
        // not affected.
        struct FieldFilter {
            std::string name;
            double min;
            double max;
        };

        // Initialize: read the output directory. Each rank reads a contiguous
        // range of the ncpu RAMSES files (amr/hydro/part_NNNNN.outCCCCC).
        // @param output_dir: path of the output_NNNNN directory
        // @param read_gas / read_particles: skip the grid or the particles entirely
        // @param level_max: cells above this level are dropped and their
        //        parent cells kept as leaves (0 = all levels)
        void init_lib(const std::string& output_dir, int world_rank, int world_size,
            bool read_gas, bool read_particles, int level_max,
            const std::vector<FieldFilter>& filters);

        // Finalize and release the data.
        void finish_lib();

        // Available (type, block) pairs: types_and_blocks[PTMax * block + type] > 0.
        void get_types_and_blocks(std::vector<int>& types_and_blocks);

        // Print types and blocks.
        void print_types_and_blocks_local();
        void print_types_and_blocks(std::vector<int>& types_and_blocks);

        // Name of a block.
        std::string get_dataset_name(int blocknr);

        // Smoothing length: gas cell size; particles: size of the cells of
        // the AMR level they are attached to (levelp); sinks: size of the
        // cells of their level (the "level" column of the sink file).
        double get_particle_hsml(uint64_t id);

        // Mass: gas density x cell volume; particles: mass (code units).
        double get_particle_mass(uint64_t id);

        // Density: gas density (code units); particles: 0.
        double get_particle_rho(uint64_t id);
        int get_particle_rho_blocknr();

        // Physical units of the run (cgs), from the info file.
        double get_unit_l();
        double get_unit_d();
        double get_unit_t();
        double get_time();

    } // namespace io

} // namespace ramses
