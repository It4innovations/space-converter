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

#include "haccbin_convert_vdb.h"
#include "haccbin_extract_iolib.h"

namespace haccbin {
	std::string format_filename(const std::string& pattern, int number) {
		// Create the formatted number as a string
		std::ostringstream formattedNumber;
		// if (number < 1000) {
		// 	formattedNumber << std::setw(3) << std::setfill('0') << number;
		// }
		// else {
		formattedNumber << number;
		//}

		std::string result = pattern;
		std::string placeholder = "{}";
		size_t pos = result.find(placeholder);

		while (pos != std::string::npos) {
			// Replace the first occurrence of the placeholder with the formatted number
			result.replace(pos, placeholder.length(), formattedNumber.str());

			// Find the next occurrence of the placeholder
			pos = result.find(placeholder, pos + formattedNumber.str().length());
		}

		return result;
	}

	void ConvertVDBHACCBin::print_CPU_steps() {
		haccbin::io::print_CPU_steps();
	}

	float ConvertVDBHACCBin::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return haccbin::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBHACCBin::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return haccbin::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBHACCBin::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return haccbin::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBHACCBin::get_particle_type(uint64_t id) {
		return haccbin::io::get_particle_type(id);
	}
	void ConvertVDBHACCBin::get_particle_position(uint64_t id, double* pos) const {
		haccbin::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBHACCBin::get_local_num_particles() const {
		return haccbin::io::get_local_num_particles();
	}
	size_t ConvertVDBHACCBin::get_global_num_particles() const {
		return haccbin::io::get_global_num_particles();
	}

	//double ConvertVDBHACCBin::get_particle_radius(uint64_t id) {
	//	return haccbin::io::get_particle_radius(id);
	//}

	double ConvertVDBHACCBin::get_particle_hsml(uint64_t id) {
		return haccbin::io::get_particle_hsml(id);
	}

	double ConvertVDBHACCBin::get_particle_mass(uint64_t id) {
		return haccbin::io::get_particle_mass(id);
	}

	double ConvertVDBHACCBin::get_particle_rho_internal(uint64_t id) {
		return haccbin::io::get_particle_rho(id);
	}

	int ConvertVDBHACCBin::get_particle_rho_blocknr() {
		return haccbin::io::get_particle_rho_blocknr();
	}

	void ConvertVDBHACCBin::init_lib(int argc, char** argv, int world_rank, int world_size) {	
		std::string haccbin_file;
		bool use_anim = false;
		int anim_start = -1;
		int anim_end = -1;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--haccbin-file") {
				haccbin_file = argv[++i];
			}
			else if (arg == "--anim") {
				use_anim = true;
				anim_start = std::stoi(argv[++i]);
				anim_end = std::stoi(argv[++i]);
			}
		}

		// Anim
		if (use_anim) {
			haccbin_file = format_filename(haccbin_file, anim_start + world_rank);
			std::cout << "Reading timestep file: " << haccbin_file << std::endl;
		}

		if (use_anim) {
			world_rank = 0;
			world_size = 1;
		}

		haccbin::io::init_lib(haccbin_file, world_rank, world_size);

		print_CPU_steps();

		//std::string param_file;
		//std::string snap_file;
		//int MaxMemSize = 1000;
		//double BufferSize = 100;
		//double PartAllocFactor = 1.2;
		//bool use_anim = false;
		//int anim_start = -1;
		//int anim_end = -1;
		//int TotBHs = 0;

		//for (int i = 1; i < argc; i++) {
		//	const std::string arg = argv[i];
		//	if (arg == "--param-file") {
		//		param_file = argv[++i]; //std::stof(av[++i]);
		//	}
		//	else if (arg == "--gadget-file") {
		//		snap_file = argv[++i]; //std::stoi(argv[++i]);
		//	}
		//	else if (arg == "--max-mem-size") {
		//		MaxMemSize = std::stoi(argv[++i]);
		//	}
		//	else if (arg == "--buffer-size") {
		//		BufferSize = std::stof(argv[++i]);
		//	}
		//	else if (arg == "--part-alloc-factor") {
		//		PartAllocFactor = std::stof(argv[++i]);
		//	}
		//	else if (arg == "--anim") {
		//		use_anim = true;
		//		anim_start = std::stoi(argv[++i]);
		//		anim_end = std::stoi(argv[++i]);
		//	}
		//	else if (arg == "--bh-count") {
		//		TotBHs = std::stoi(argv[++i]);
		//	}
		//}

		//// Anim
		//if (use_anim) {
		//	snap_file = format_filename(snap_file, anim_start + world_rank);
		//	std::cout << "Reading step file: " << snap_file << std::endl;
		//}

		//if (use_anim) {
		//	world_rank = 0;
		//	world_size = 1;
		//}

		//gadget::io::gadget_init_lib(world_rank, world_size);

		////read_parameter_file(param_file, NULL, NULL, NULL, 0);
		////void set_parameter(int ICFormat, int SnapFormat, int NumFilesWrittenInParallel, int MaxMemSize, double BufferSize, double PartAllocFactor)
		//gadget::io::gadget_set_parameter(2, 2, world_size, MaxMemSize, BufferSize, PartAllocFactor, TotBHs);
		//gadget::io::gadget_mymalloc_init();
		//gadget::io::gadget_set_units();

		//// Allocate memory for the char array
		//std::vector<char> csnap_file(snap_file.length() + 1);

		//// Copy the contents of the std::string to the char array
		//std::strcpy(csnap_file.data(), snap_file.c_str());

		//gadget::io::gadget_read_ic(csnap_file.data());

		//print_CPU_steps();
	}

	void ConvertVDBHACCBin::finish_lib()
	{
		haccbin::io::finish_lib();
	}

	void ConvertVDBHACCBin::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		haccbin::io::get_types_and_blocks(types_and_blocks);
	}
	void ConvertVDBHACCBin::print_types_and_blocks_local() {
		std::vector<int> types_and_blocks;
		haccbin::io::get_types_and_blocks(types_and_blocks);

		printf("\n");

		//char buf[500];
		for (int type = 0; type < haccbin::HACCBinParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / haccbin::HACCBinParticleType::PTMax; blocknr++) {

				if (types_and_blocks[haccbin::HACCBinParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	void ConvertVDBHACCBin::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");

		//char buf[500];
		for (int type = 0; type < haccbin::HACCBinParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / haccbin::HACCBinParticleType::PTMax; blocknr++) {

				if (types_and_blocks[haccbin::HACCBinParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	std::string ConvertVDBHACCBin::get_type_name(int type) {		

		haccbin::HACCBinParticleType pt = (haccbin::HACCBinParticleType)type;
		switch (pt) {
		case haccbin::HACCBinParticleType::DarkMatter:
			return "DarkMatter";
		case haccbin::HACCBinParticleType::Baryon:
			return "Baryon";
		case haccbin::HACCBinParticleType::BaryonStart:
			return "BaryonStar";
		case haccbin::HACCBinParticleType::BaryonWind:
			return "BaryonWind";
		case haccbin::HACCBinParticleType::BaryonGas:
			return "BaryonGas";
		case haccbin::HACCBinParticleType::DarkMatterAGN:
			return "DarkMatterAGN";
		}
		return "Unknown";
	}
	std::string ConvertVDBHACCBin::get_dataset_name(int blocknr) {
		return haccbin::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBHACCBin::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < haccbin::HACCBinParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < types_and_blocks.size() / haccbin::HACCBinParticleType::PTMax; bnr++) {
				if (types_and_blocks[haccbin::HACCBinParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}
	int ConvertVDBHACCBin::get_num_types() {
		return haccbin::HACCBinParticleType::PTMax;
	}
	int ConvertVDBHACCBin::get_num_blocks() {
		std::vector<int> types_and_blocks;
		haccbin::io::get_types_and_blocks(types_and_blocks);
		return types_and_blocks.size() / haccbin::HACCBinParticleType::PTMax;
	}
}