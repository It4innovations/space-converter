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

 // Namespace fil_grace provides functionalities for handling the 3D volume output of
 // GRACE (General Relativistic Astrophysics Code for Exascale, the Kokkos / p4est
 // successor of FIL; https://github.com/GRACE-astro/grace): one HDF5 file per
 // output (volume_out_NNNNNN.h5) holding the cells of all blocks (p4est leaves)
 // of the block-structured AMR mesh.
namespace fil_grace {

	// Particle types. The cells of the AMR blocks are the only type (treated as
	// particles, like the grid readers of the other formats).
	enum FilGraceParticleType {
		Cell = 0,       // cells of the leaf blocks

		PTMax           // Maximum value for ParticleType (used for validation or iteration)
	};

	// Fixed blocks; the datasets of the file follow from BTMax on (one block per
	// scalar dataset, then one block per vector dataset).
	enum FilGraceBlockType {
		Pos = 0,        // positions only (weight 1)
		Mass = 1,       // Rho x cell volume (coordinate volume dx^3)
		Rho = 2,        // the density variable (rho, dens, else the first scalar)
		Level = 3,      // refinement level of the cell's block

		BTMax           // Maximum value for BlockType (used for validation or iteration)
	};

	// Namespace io contains I/O operations and utilities for GRACE HDF5 data processing.
	namespace io {

		struct Options {
			std::string file;                   // volume output file (--fil-grace-file)
			std::vector<std::string> vars;      // datasets to load (empty = all)
			int level_min = 0;                  // refinement levels read: level_min .. level_max
			int level_max = -1;                 // -1 = all levels
			bool has_region = false;            // read only the blocks that intersect region
			double region[6] = { 0, 0, 0, 0, 0, 0 };   // xmin ymin zmin xmax ymax zmax (code units)
			std::string mirror;                 // axes to mirror about the plane coordinate = 0 ("z", "xy", ...)
			int block_size = 0;                 // cells per block edge (0 = derived from the file)
		};

		// Print the steps executed on the CPU during the reading process.
		void print_CPU_steps();

		// Scalar ("norm") value of a particle in a block (vectors: magnitude).
		float get_particle_norm_value(int blocknr, uint64_t id);

		// Original components of a particle in a block; returns their count.
		int get_particle_value(int blocknr, uint64_t id, float* out_value);

		// Component count of a particle in a block (0 = not available).
		int get_particle_value_comp(int blocknr, uint64_t id);

		// Type of a particle (FilGraceParticleType).
		int get_particle_type(uint64_t id);

		// Cell centre, Cartesian, in code length units (M for GR runs).
		void get_particle_position(uint64_t id, double* pos);

		// Number of cells read by this rank.
		size_t get_local_num_particles();

		// Number of cells over all ranks.
		size_t get_global_num_particles();

		// Initialize: read the block geometry and the selected datasets. The
		// blocks are split into contiguous ranges of equal block counts over the
		// ranks; every rank reads its own blocks only.
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

		// Smoothing length: the cell size of the cell's block.
		double get_particle_hsml(uint64_t id);

		// Mass: Rho x dx^3 (code units).
		double get_particle_mass(uint64_t id);

		// Density: Rho (code units).
		double get_particle_rho(uint64_t id);
		int get_particle_rho_blocknr();

		// Simulation time and iteration of the data read.
		double get_time();
		long long get_iteration();

	} // namespace io

} // namespace fil_grace
