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

 // Namespace fil provides functionalities for handling the 3D output of FIL
 // (Frankfurt/IllinoisGRMHD) and other Einstein Toolkit codes: Carpet HDF5
 // files (CarpetIOHDF5 out3D / out_vars, one file per variable or group, or
 // per process, *.file_N.h5) holding the box-in-box mesh-refinement hierarchy.
namespace fil {

	// Particle types. The grid points of the refinement hierarchy that are not
	// covered by a finer level are the only type (treated as particles, like the
	// grid readers of the other formats).
	enum FilParticleType {
		Cell = 0,       // grid points (finest level wins)

		PTMax           // Maximum value for ParticleType (used for validation or iteration)
	};

	// Fixed blocks; the grid functions of the files follow from BTMax on (one
	// scalar block per variable, then one vector block per vector group: names
	// ending in [0] [1] [2], e.g. HydroBase::vel, or <base>x <base>y <base>z,
	// e.g. IllinoisGRMHD Bx By Bz).
	enum FilBlockType {
		Pos = 0,        // positions only (weight 1)
		Mass = 1,       // Rho x cell volume (coordinate volume delta^3)
		Rho = 2,        // the density variable (rho, rho_b, else the first variable)
		Level = 3,      // refinement level of the point

		BTMax           // Maximum value for BlockType (used for validation or iteration)
	};

	// Namespace io contains I/O operations and utilities for Carpet HDF5 data processing.
	namespace io {

		struct Options {
			std::vector<std::string> files;     // Carpet HDF5 files (--fil-file, repeatable)
			std::vector<std::string> dirs;      // directories: every *.h5 in them (--fil-dir, repeatable)
			std::vector<std::string> vars;      // variables to load (empty = all)
			long long iteration = -1;           // -1 = the latest iteration present in the files
			int level_min = 0;                  // refinement levels read: level_min .. level_max
			int level_max = -1;                 // -1 = all levels
			bool keep_ghosts = false;           // keep the inter-process ghost zones
		};

		// Print the steps executed on the CPU during the reading process.
		void print_CPU_steps();

		// Scalar ("norm") value of a particle in a block (vectors: magnitude).
		float get_particle_norm_value(int blocknr, uint64_t id);

		// Original components of a particle in a block; returns their count.
		int get_particle_value(int blocknr, uint64_t id, float* out_value);

		// Component count of a particle in a block (0 = not available).
		int get_particle_value_comp(int blocknr, uint64_t id);

		// Type of a particle (FilParticleType).
		int get_particle_type(uint64_t id);

		// Grid point position, Cartesian, in code length units (M for GR runs).
		void get_particle_position(uint64_t id, double* pos);

		// Number of grid points read by this rank.
		size_t get_local_num_particles();

		// Number of grid points over all ranks.
		size_t get_global_num_particles();

		// Initialize: index the files, pick the iteration, and read it. The
		// patches (refinement level, component) are split into contiguous
		// ranges of about equal point counts over the ranks; every rank reads
		// its own patches only (all ranks read the patch geometry).
		void init_lib(const Options& options, int world_rank, int world_size);

		// Finalize and release the data.
		void finish_lib();

		// Available (type, block) pairs: types_and_blocks[PTMax * block + type] > 0.
		void get_types_and_blocks(std::vector<int>& types_and_blocks);

		// Print types and blocks.
		void print_types_and_blocks_local();
		void print_types_and_blocks(std::vector<int>& types_and_blocks);

		// Name of a block.
		std::string get_dataset_name(int blocknr);

		// Smoothing length: the grid spacing of the point's refinement level.
		double get_particle_hsml(uint64_t id);

		// Mass: Rho x delta^3 (code units).
		double get_particle_mass(uint64_t id);

		// Density: Rho (code units).
		double get_particle_rho(uint64_t id);
		int get_particle_rho_blocknr();

		// Simulation time and iteration of the data read.
		double get_time();
		long long get_iteration();

	} // namespace io

} // namespace fil
