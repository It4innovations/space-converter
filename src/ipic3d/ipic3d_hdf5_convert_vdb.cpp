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

#include "ipic3d_hdf5_convert_vdb.h"
#include "ipic3d_hdf5_extract_iolib.h"

namespace ipic3d {
    void ConvertVDBIPIC3DHDF5::print_CPU_steps() {
        ipic3d::io::print_CPU_steps();
    }

    float ConvertVDBIPIC3DHDF5::get_particle_norm_value_internal(int blocknr, uint64_t id) {
        return ipic3d::io::get_particle_norm_value(blocknr, id);
    }

    int ConvertVDBIPIC3DHDF5::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
        return ipic3d::io::get_particle_value(blocknr, id, value);
    }

    int ConvertVDBIPIC3DHDF5::get_particle_value_comp_internal(int blocknr, uint64_t id) {
        return ipic3d::io::get_particle_value_comp(blocknr, id);
    }

    int ConvertVDBIPIC3DHDF5::get_particle_type(uint64_t id) {
        return ipic3d::io::get_particle_type(id);
    }

    void ConvertVDBIPIC3DHDF5::get_particle_position(uint64_t id, double* pos) const {
        ipic3d::io::get_particle_position(id, pos);
    }

    size_t ConvertVDBIPIC3DHDF5::get_local_num_particles() const {
        return ipic3d::io::get_local_num_particles();
    }

    size_t ConvertVDBIPIC3DHDF5::get_global_num_particles() const {
        return ipic3d::io::get_global_num_particles();
    }

    double ConvertVDBIPIC3DHDF5::get_particle_hsml(uint64_t id) {
        return ipic3d::io::get_particle_hsml(id);
    }

    double ConvertVDBIPIC3DHDF5::get_particle_mass(uint64_t id) {
        return ipic3d::io::get_particle_mass(id);
    }

    double ConvertVDBIPIC3DHDF5::get_particle_rho_internal(uint64_t id) {
        return ipic3d::io::get_particle_rho(id);
    }

    int ConvertVDBIPIC3DHDF5::get_particle_rho_blocknr() {
        return ipic3d::io::get_particle_rho_blocknr();
    }

    void ConvertVDBIPIC3DHDF5::init_lib(int argc, char** argv, int world_rank, int world_size) {
        std::string hdf5_file;
        int species_id = 0;
        int cycle_id = 0;
        int num_files = 1;

        bool use_anim = false;
        int anim_start = -1;
        int anim_end = -1;
        int anim_step = -1;

        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--hdf5-file") {
                hdf5_file = argv[++i];
            }
            else if (arg == "--species-id") {
                species_id = std::stoi(argv[++i]);
            }
            else if (arg == "--cycle-id") {
                cycle_id = std::stoi(argv[++i]);
            }
            else if (arg == "--anim") {
                use_anim = true;
                anim_start = std::stoi(argv[++i]);
                anim_end = std::stoi(argv[++i]);
                anim_step = std::stoi(argv[++i]);
            }
            else if (arg == "--num-files") {
                num_files = std::stoi(argv[++i]);
            }
        }

        // TODO
        // Anim: Use world_rank to select which cycle to load
        // if (use_anim) {
        //     cycle_id = anim_start + anim_step * world_rank;
        //     std::cout << "Reading timestep cycle: " << cycle_id << std::endl;
        //     world_rank = 0;
        //     world_size = 1;
        // }

        ipic3d::io::init_lib(hdf5_file, world_rank, world_size, species_id, cycle_id, num_files);

        print_CPU_steps();
    }

    void ConvertVDBIPIC3DHDF5::finish_lib() {
        ipic3d::io::finish_lib();
    }

    void ConvertVDBIPIC3DHDF5::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
        ipic3d::io::get_types_and_blocks(types_and_blocks);
    }

    void ConvertVDBIPIC3DHDF5::print_types_and_blocks_local() {
        std::vector<int> types_and_blocks;
        ipic3d::io::get_types_and_blocks(types_and_blocks);

        printf("\n");

        for (int type = 0; type < ipic3d::IPIC3DParticleType::PTMax; type++) {
            printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
            for (int blocknr = 0; blocknr < types_and_blocks.size() / ipic3d::IPIC3DParticleType::PTMax; blocknr++) {
                if (types_and_blocks[ipic3d::IPIC3DParticleType::PTMax * blocknr + type] > 0) {
                    std::string buf = get_dataset_name(blocknr);
                    printf("\t%s (%d)\n", buf.c_str(), blocknr);
                }
            }
        }
    }

    void ConvertVDBIPIC3DHDF5::print_types_and_blocks(std::vector<int>& types_and_blocks) {
        printf("\nAll snapshots contain:\n");

        for (int type = 0; type < ipic3d::IPIC3DParticleType::PTMax; type++) {
            printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
            for (int blocknr = 0; blocknr < types_and_blocks.size() / ipic3d::IPIC3DParticleType::PTMax; blocknr++) {
                if (types_and_blocks[ipic3d::IPIC3DParticleType::PTMax * blocknr + type] > 0) {
                    std::string buf = get_dataset_name(blocknr);
                    printf("\t%s (%d)\n", buf.c_str(), blocknr);
                }
            }
        }
    }

    std::string ConvertVDBIPIC3DHDF5::get_type_name(int type) {
        ipic3d::IPIC3DParticleType pt = (ipic3d::IPIC3DParticleType)type;
        switch (pt) {
        case ipic3d::IPIC3DParticleType::Species_0:
            return "Species_0";
        case ipic3d::IPIC3DParticleType::Species_1:
            return "Species_1";
        case ipic3d::IPIC3DParticleType::Species_2:
            return "Species_2";
        case ipic3d::IPIC3DParticleType::Species_3:
            return "Species_3";
        }
        return "Unknown";
    }

    std::string ConvertVDBIPIC3DHDF5::get_dataset_name(int blocknr) {
        return ipic3d::io::get_dataset_name(blocknr);
    }

    std::string ConvertVDBIPIC3DHDF5::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
        std::string particle_data_types = "";

        for (int t = 0; t < ipic3d::IPIC3DParticleType::PTMax; t++) {
            for (int bnr = 0; bnr < types_and_blocks.size() / ipic3d::IPIC3DParticleType::PTMax; bnr++) {
                if (types_and_blocks[ipic3d::IPIC3DParticleType::PTMax * bnr + t] == 0)
                    continue;

                //std::string type_name = get_type_name(t);
                //std::string block_name = get_dataset_name(bnr);

                //if (particle_data_types.length() > 0)
                //    particle_data_types += ",";

                //particle_data_types += type_name + "_" + block_name;
                particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
            }
        }

        return particle_data_types;
    }

    int ConvertVDBIPIC3DHDF5::get_num_types() {
        return ipic3d::IPIC3DParticleType::PTMax;
    }

    int ConvertVDBIPIC3DHDF5::get_num_blocks() {
        std::vector<int> types_and_blocks;
        ipic3d::io::get_types_and_blocks(types_and_blocks);
        return types_and_blocks.size() / ipic3d::IPIC3DParticleType::PTMax;
    }

} // namespace ipic3d
