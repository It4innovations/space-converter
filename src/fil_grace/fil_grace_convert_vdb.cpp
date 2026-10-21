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

#include "fil_grace_convert_vdb.h"
#include "fil_grace_extract_iolib.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace fil_grace {
	void ConvertVDBFilGrace::print_CPU_steps() {
		fil_grace::io::print_CPU_steps();
	}

	float ConvertVDBFilGrace::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return fil_grace::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBFilGrace::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return fil_grace::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBFilGrace::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return fil_grace::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBFilGrace::get_particle_type(uint64_t id) {
		return fil_grace::io::get_particle_type(id);
	}
	void ConvertVDBFilGrace::get_particle_position(uint64_t id, double* pos) const {
		fil_grace::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBFilGrace::get_local_num_particles() const {
		return fil_grace::io::get_local_num_particles();
	}
	size_t ConvertVDBFilGrace::get_global_num_particles() const {
		return fil_grace::io::get_global_num_particles();
	}

	double ConvertVDBFilGrace::get_particle_hsml(uint64_t id) {
		return fil_grace::io::get_particle_hsml(id);
	}

	double ConvertVDBFilGrace::get_particle_mass(uint64_t id) {
		return fil_grace::io::get_particle_mass(id);
	}

	double ConvertVDBFilGrace::get_particle_rho_internal(uint64_t id) {
		return fil_grace::io::get_particle_rho(id);
	}

	int ConvertVDBFilGrace::get_particle_rho_blocknr() {
		return fil_grace::io::get_particle_rho_blocknr();
	}

	void ConvertVDBFilGrace::init_lib(int argc, char** argv, int world_rank, int world_size) {
		fil_grace::io::Options options;

		bool use_anim = false;
		int anim_start = -1;
		int anim_end = -1;
		int anim_step = -1;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--fil-grace-file") {
				options.file = argv[++i];
			}
			else if (arg == "--fil-grace-vars") {
				// comma- or space-separated list in one argument
				std::string list = argv[++i];
				for (char& c : list)
					if (c == ',') c = ' ';
				std::stringstream ss(list);
				std::string v;
				while (ss >> v)
					options.vars.push_back(v);
			}
			else if (arg == "--fil-grace-levels") {
				options.level_min = std::stoi(argv[++i]);
				options.level_max = std::stoi(argv[++i]);
			}
			else if (arg == "--fil-grace-region") {
				options.has_region = true;
				for (int a = 0; a < 6; a++)
					options.region[a] = std::stod(argv[++i]);
			}
			else if (arg == "--fil-grace-mirror") {
				options.mirror = argv[++i];
			}
			else if (arg == "--fil-grace-block-size") {
				options.block_size = std::stoi(argv[++i]);
			}
			else if (arg == "--anim") {
				use_anim = true;
				anim_start = std::stoi(argv[++i]);
				anim_end = std::stoi(argv[++i]);
				anim_step = std::stoi(argv[++i]);
			}
		}

		// Anim: one volume output per rank (the frame number is the GRACE iteration,
		// zero-padded to 6 digits as in volume_out_NNNNNN.h5; it replaces "{}" in --fil-grace-file)
		if (use_anim) {
			// Clamp to anim_end so trailing ranks do not address iterations past the animation range
			int anim_frame = anim_start + anim_step * world_rank;
			if (anim_frame > anim_end)
				anim_frame = anim_end;
			options.file = format_filename(options.file, anim_frame, 6);
			std::cout << "Reading FIL_GRACE volume output: " << options.file << std::endl;
			world_rank = 0;
			world_size = 1;
		}

		fil_grace::io::init_lib(options, world_rank, world_size);

		print_CPU_steps();
	}

	void ConvertVDBFilGrace::finish_lib()
	{
		fil_grace::io::finish_lib();
	}

	void ConvertVDBFilGrace::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		fil_grace::io::get_types_and_blocks(types_and_blocks);
	}

	void ConvertVDBFilGrace::print_types_and_blocks_local() {
		fil_grace::io::print_types_and_blocks_local();
	}

	void ConvertVDBFilGrace::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");
		fil_grace::io::print_types_and_blocks(types_and_blocks);
	}

	std::string ConvertVDBFilGrace::get_type_name(int type) {
		switch ((fil_grace::FilGraceParticleType)type) {
		case fil_grace::FilGraceParticleType::Cell: return "Cell";
		default: break;
		}
		return "Unknown";
	}

	std::string ConvertVDBFilGrace::get_dataset_name(int blocknr) {
		return fil_grace::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBFilGrace::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < fil_grace::FilGraceParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < (int)types_and_blocks.size() / fil_grace::FilGraceParticleType::PTMax; bnr++) {
				if (types_and_blocks[fil_grace::FilGraceParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBFilGrace::get_num_types() {
		return fil_grace::FilGraceParticleType::PTMax;
	}

	int ConvertVDBFilGrace::get_num_blocks() {
		std::vector<int> types_and_blocks;
		fil_grace::io::get_types_and_blocks(types_and_blocks);
		return (int)types_and_blocks.size() / fil_grace::FilGraceParticleType::PTMax;
	}
}
