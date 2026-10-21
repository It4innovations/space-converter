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

 // Namespace bhac provides functionalities for handling BHAC (MPI-AMRVAC family)
 // native snapshots (dataNNNN.dat: block-AMR forest + per-block cell data).
namespace bhac {

	// Particle types. The leaf cells of the block-AMR grid are the only type
	// (treated as particles, like the grid readers of the other formats).
	enum BhacParticleType {
		Cell = 0,       // leaf cells

		PTMax           // Maximum value for ParticleType (used for validation or iteration)
	};

	// Fixed blocks; the stored variables of the snapshot (wnames of the .par
	// file, e.g. d s1 s2 s3 tau b1 b2 b3 Ds dtr1 lfac xi) follow from BTMax on.
	enum BhacBlockType {
		Pos = 0,        // positions only (weight 1)
		Mass = 1,       // Rho x cell volume (flat-space coordinate volume)
		Rho = 2,        // rest-mass density d / lfac if both are stored, else the first variable
		BVec = 3,       // magnetic field b1 b2 b3 (contravariant, code coordinate basis) as a
		                // Cartesian vector through the flat-space Jacobian of the coordinates
		Level = 4,      // AMR level of the cell

		BTMax           // Maximum value for BlockType (used for validation or iteration)
	};

	// Namespace io contains I/O operations and utilities for BHAC data processing.
	namespace io {

		// Mapping of the code coordinates (x1, x2, x3) to Cartesian positions.
		enum class Coord {
			Cart,       // x, y, z
			Sph,        // r, theta, phi (Kerr-Schild / Boyer-Lindquist / flat spherical)
			MKS,        // modified Kerr-Schild: r = R0 + exp(x1), theta = x2 + h/2 sin(2 x2), phi = x3
		};

		struct Options {
			std::string dat_file;       // dataNNNN.dat
			std::string par_file;       // the run's .par file (base grid, domain, typeaxial, wnames)
			Coord coord = Coord::MKS;   // used when typeaxial = spherical; slab is always Cart
			bool coord_set = false;
			double mks_h = 0.0;         // MKS theta squeeze (coordpar(h_) of the run)
			double mks_r0 = 0.0;        // MKS radial offset (coordpar(R0_) of the run)
			double rmin = 0.0;          // keep cells with rmin <= r <= rmax (r = |x| for slab)
			double rmax = 0.0;          // 0 = no limit
		};

		// Print the steps executed on the CPU during the reading process.
		void print_CPU_steps();

		// Scalar ("norm") value of a particle in a block (vectors: magnitude).
		float get_particle_norm_value(int blocknr, uint64_t id);

		// Original components of a particle in a block; returns their count.
		int get_particle_value(int blocknr, uint64_t id, float* out_value);

		// Component count of a particle in a block (0 = not available).
		int get_particle_value_comp(int blocknr, uint64_t id);

		// Type of a particle (BhacParticleType).
		int get_particle_type(uint64_t id);

		// Cell centre, Cartesian, in code length units (M for GRMHD runs).
		void get_particle_position(uint64_t id, double* pos);

		// Number of leaf cells read by this rank.
		size_t get_local_num_particles();

		// Number of leaf cells over all ranks.
		size_t get_global_num_particles();

		// Initialize: read the snapshot. The leaf blocks (in file order) are
		// split into contiguous ranges over the ranks; every rank reads its own
		// blocks only.
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

		// Smoothing length: the largest Cartesian extent of the cell
		// (max(dr, r dtheta, r sin(theta) dphi) for spherical grids).
		double get_particle_hsml(uint64_t id);

		// Mass: Rho x cell volume (code units).
		double get_particle_mass(uint64_t id);

		// Density: Rho (code units).
		double get_particle_rho(uint64_t id);
		int get_particle_rho_blocknr();

		// Simulation time and iteration of the snapshot.
		double get_time();
		int get_iteration();

	} // namespace io

} // namespace bhac
