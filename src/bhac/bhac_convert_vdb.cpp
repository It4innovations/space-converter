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

#include "bhac_convert_vdb.h"
#include "bhac_extract_iolib.h"

#include <iostream>
#include <stdexcept>

namespace bhac {
	void ConvertVDBBhac::print_CPU_steps() {
		bhac::io::print_CPU_steps();
	}

	float ConvertVDBBhac::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return bhac::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBBhac::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return bhac::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBBhac::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return bhac::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBBhac::get_particle_type(uint64_t id) {
		return bhac::io::get_particle_type(id);
	}
	void ConvertVDBBhac::get_particle_position(uint64_t id, double* pos) const {
		bhac::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBBhac::get_local_num_particles() const {
		return bhac::io::get_local_num_particles();
	}
	size_t ConvertVDBBhac::get_global_num_particles() const {
		return bhac::io::get_global_num_particles();
	}

	double ConvertVDBBhac::get_particle_hsml(uint64_t id) {
		return bhac::io::get_particle_hsml(id);
	}

	double ConvertVDBBhac::get_particle_mass(uint64_t id) {
		return bhac::io::get_particle_mass(id);
	}

	double ConvertVDBBhac::get_particle_rho_internal(uint64_t id) {
		return bhac::io::get_particle_rho(id);
	}

	int ConvertVDBBhac::get_particle_rho_blocknr() {
		return bhac::io::get_particle_rho_blocknr();
	}

	void ConvertVDBBhac::init_lib(int argc, char** argv, int world_rank, int world_size) {
		bhac::io::Options options;

		bool use_anim = false;
		int anim_start = -1;
		int anim_end = -1;
		int anim_step = -1;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--bhac-file") {
				options.dat_file = argv[++i];
			}
			else if (arg == "--bhac-par") {
				options.par_file = argv[++i];
			}
			else if (arg == "--bhac-coord") {
				const std::string c = argv[++i];
				if (c == "mks") options.coord = bhac::io::Coord::MKS;
				else if (c == "sph" || c == "ks") options.coord = bhac::io::Coord::Sph;
				else if (c == "cart") options.coord = bhac::io::Coord::Cart;
				else throw std::runtime_error("--bhac-coord: expected mks, sph (ks) or cart, got " + c);
				options.coord_set = true;
			}
			else if (arg == "--bhac-mks") {
				options.mks_h = std::stod(argv[++i]);
				options.mks_r0 = std::stod(argv[++i]);
			}
			else if (arg == "--bhac-rrange") {
				options.rmin = std::stod(argv[++i]);
				options.rmax = std::stod(argv[++i]);
			}
			else if (arg == "--anim") {
				use_anim = true;
				anim_start = std::stoi(argv[++i]);
				anim_end = std::stoi(argv[++i]);
				anim_step = std::stoi(argv[++i]);
			}
		}

		// Anim: one snapshot per rank (the pattern holds the snapshot number)
		if (use_anim) {
			// Clamp to anim_end so trailing ranks do not address snapshots past the animation range
			int anim_frame = anim_start + anim_step * world_rank;
			if (anim_frame > anim_end)
				anim_frame = anim_end;
			options.dat_file = format_filename(options.dat_file, anim_frame);
			std::cout << "Reading BHAC snapshot: " << options.dat_file << std::endl;
			world_rank = 0;
			world_size = 1;
		}

		bhac::io::init_lib(options, world_rank, world_size);

		print_CPU_steps();
	}

	void ConvertVDBBhac::finish_lib()
	{
		bhac::io::finish_lib();
	}

	void ConvertVDBBhac::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		bhac::io::get_types_and_blocks(types_and_blocks);
	}

	void ConvertVDBBhac::print_types_and_blocks_local() {
		bhac::io::print_types_and_blocks_local();
	}

	void ConvertVDBBhac::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");
		bhac::io::print_types_and_blocks(types_and_blocks);
	}

	std::string ConvertVDBBhac::get_type_name(int type) {
		switch ((bhac::BhacParticleType)type) {
		case bhac::BhacParticleType::Cell: return "Cell";
		default: break;
		}
		return "Unknown";
	}

	std::string ConvertVDBBhac::get_dataset_name(int blocknr) {
		return bhac::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBBhac::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < bhac::BhacParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < (int)types_and_blocks.size() / bhac::BhacParticleType::PTMax; bnr++) {
				if (types_and_blocks[bhac::BhacParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBBhac::get_num_types() {
		return bhac::BhacParticleType::PTMax;
	}

	int ConvertVDBBhac::get_num_blocks() {
		std::vector<int> types_and_blocks;
		bhac::io::get_types_and_blocks(types_and_blocks);
		return (int)types_and_blocks.size() / bhac::BhacParticleType::PTMax;
	}
}
