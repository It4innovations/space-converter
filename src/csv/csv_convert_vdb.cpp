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

#include "csv_convert_vdb.h"
#include "csv_extract_iolib.h"

namespace csv {
	void ConvertVDBCSV::print_CPU_steps() {
		csv::io::print_CPU_steps();
	}

	float ConvertVDBCSV::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return csv::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBCSV::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return csv::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBCSV::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return csv::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBCSV::get_particle_type(uint64_t id) {
		return csv::io::get_particle_type(id);
	}
	void ConvertVDBCSV::get_particle_position(uint64_t id, double* pos) const {
		csv::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBCSV::get_local_num_particles() const {
		return csv::io::get_local_num_particles();
	}
	size_t ConvertVDBCSV::get_global_num_particles() const {
		return csv::io::get_global_num_particles();
	}

	//double ConvertVDBCSV::get_particle_radius(uint64_t id) {
	//	return csv::io::get_particle_radius(id);
	//}

	double ConvertVDBCSV::get_particle_hsml(uint64_t id) {
		return csv::io::get_particle_hsml(id);
	}

	double ConvertVDBCSV::get_particle_mass(uint64_t id) {
		return csv::io::get_particle_mass(id);
	}

	double ConvertVDBCSV::get_particle_rho_internal(uint64_t id) {
		return csv::io::get_particle_rho(id);
	}

	int ConvertVDBCSV::get_particle_rho_blocknr() {
		return csv::io::get_particle_rho_blocknr();
	}

	void ConvertVDBCSV::init_lib(int argc, char** argv, int world_rank, int world_size) {	
		std::string csv_file;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--csv-file") {
				csv_file = argv[++i]; //std::stoi(argv[++i]);
			}
		}

		csv::io::init_lib(csv_file, world_rank, world_size);

		print_CPU_steps();
	}

	void ConvertVDBCSV::finish_lib()
	{
		csv::io::finish_lib();
	}

	void ConvertVDBCSV::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		csv::io::get_types_and_blocks(types_and_blocks);
	}
	void ConvertVDBCSV::print_types_and_blocks_local() {
		std::vector<int> types_and_blocks;
		csv::io::get_types_and_blocks(types_and_blocks);

		printf("\n");

		//char buf[500];
		for (int type = 0; type < csv::CSVParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / csv::CSVParticleType::PTMax; blocknr++) {

				if (types_and_blocks[csv::CSVParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	void ConvertVDBCSV::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");

		//char buf[500];
		for (int type = 0; type < csv::CSVParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / csv::CSVParticleType::PTMax; blocknr++) {

				if (types_and_blocks[csv::CSVParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	std::string ConvertVDBCSV::get_type_name(int type) {		

		csv::CSVParticleType pt = (csv::CSVParticleType)type;
		switch (pt) {
		case csv::CSVParticleType::Gas:
			return "Gas";
		case csv::CSVParticleType::Dark:
			return "Dark";
		case csv::CSVParticleType::Stars:
			return "Stars";
		}
		return "Unknown";
	}
	std::string ConvertVDBCSV::get_dataset_name(int blocknr) {
		return csv::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBCSV::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < csv::CSVParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < types_and_blocks.size() / csv::CSVParticleType::PTMax; bnr++) {
				if (types_and_blocks[csv::CSVParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBCSV::get_num_types() {
		return csv::CSVParticleType::PTMax;
	}
	int ConvertVDBCSV::get_num_blocks() {
		std::vector<int> types_and_blocks;
		csv::io::get_types_and_blocks(types_and_blocks);
		return types_and_blocks.size() / csv::CSVParticleType::PTMax;
	}
}