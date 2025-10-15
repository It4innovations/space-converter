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

#include "changa_tipsy_extract_iolib.h"

#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

 //#include "tree_xdr.h"
#include "TipsyFile.h"

#include <omp.h>

#include "convert_common.h"

#define RETURN_NORM_EMPTY return 0; //return std::numeric_limits<float>::quiet_NaN();//0;
#define RETURN_NORM_VALUE(v) return (float)(v);
#define RETURN_NORM_VECTOR3(v) return (float)common::calculate_dmagnitude3(v[0], v[1], v[2]);
#define RETURN_NORM_DVECTORN(v,n) return (float)common::calculate_dmagnituden(v,n);
#define RETURN_NORM_FVECTORN(v,n) return (float)common::calculate_fmagnituden(v,n);

#define RETURN_COMP_EMPTY return 0;
#define RETURN_COMP_VALUE(v) return 1;
#define RETURN_COMP_VECTOR3(v) return 3;
#define RETURN_COMP_DVECTORN(v,n) return n;
#define RETURN_COMP_FVECTORN(v,n) return n;

#define RETURN_ORIG_EMPTY RETURN_COMP_EMPTY
#define RETURN_ORIG_VALUE(v) out_value[0] = (float)v; RETURN_COMP_VALUE(v)
#define RETURN_ORIG_VECTOR3(v) out_value[0] = (float)v[0]; out_value[1] = (float)v[1]; out_value[2] = (float)v[2]; RETURN_COMP_VECTOR3(v)
#define RETURN_ORIG_DVECTORN(v,n) for(int iv=0;iv<n;iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_DVECTORN(v,n)
#define RETURN_ORIG_FVECTORN(v,n) for(int iv=0;iv<n;iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_FVECTORN(v,n)

#ifdef _WIN32
#	define STDMAX max
#	define STDMIN min
#else
#	define STDMAX std::max
#	define STDMIN std::min
#endif

namespace changa {
	namespace tipsy {
		namespace io {

			//static double fh_time; // gross, but quick way to get time
			double steps_time[2];

			//Tipsy::TipsyFile* tf_data = nullptr;
			Tipsy::PartialTipsyFile* tf_data = nullptr;
			size_t g_total_particles = 0;
			int g_start_particles = 0;
			int g_num_particles = 0;

			int g_smoothlength_blocknr = -1;

			void init_lib(std::string basefile, int world_rank, int world_size) {
				steps_time[0] = omp_get_wtime();

				//tf_data = new Tipsy::TipsyFile(basefile);
				Tipsy::TipsyReader treader(basefile);
				Tipsy::header th = treader.getHeader();

				g_total_particles = th.nbodies;//STDMAX(th.nsph, STDMAX(th.ndark, th.nstar));

				// Each process calculates its range of particles to read
				int particles_per_process = g_total_particles / world_size;
				g_start_particles = world_rank * particles_per_process;
				g_num_particles = (world_rank == world_size - 1) ? g_total_particles - g_start_particles : particles_per_process;

				tf_data = new Tipsy::PartialTipsyFile(basefile, g_start_particles, g_start_particles + g_num_particles);

				//check smoothlength
				for (int i = 0; i < tf_data->aux_values.size(); i++) {
					if (tf_data->aux_values[i].name == "smoothlength") {
						g_smoothlength_blocknr = i + changa::tipsy::BlockType::BTMax;
						break;
					}
				}

				steps_time[1] = omp_get_wtime();
			}

			void finish_lib() {
				if (tf_data != nullptr)
					delete tf_data;

				tf_data = nullptr;
			}

			uint64_t get_particle_type_offset(ParticleType pt) {
				switch (pt) {
				case ParticleType::Gas:
					return 0;
				case ParticleType::Dark:
					return tf_data->gas.size();
				case ParticleType::Stars:
					return tf_data->gas.size() + tf_data->darks.size();
				}

				return 0;
			}

			void print_CPU_steps() {
				printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
			}

			float get_particle_norm_value(int blocknr, uint64_t id) {
				ParticleType pt = (ParticleType)get_particle_type(id);
				uint64_t offset = get_particle_type_offset(pt);

				BlockType bt = (BlockType)blocknr;

				switch (pt) {
				case ParticleType::Gas: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_NORM_VECTOR3(tf_data->gas[id - offset].pos);
					case BlockType::Mass:
						RETURN_NORM_VALUE(tf_data->gas[id - offset].mass);
					case BlockType::Vel:
						RETURN_NORM_VECTOR3(tf_data->gas[id - offset].vel);
					case BlockType::Phi:
						RETURN_NORM_VALUE(tf_data->gas[id - offset].phi);
					case BlockType::HSmooth:
						RETURN_NORM_VALUE(tf_data->gas[id - offset].hsmooth);
					case BlockType::Rho:
						RETURN_NORM_VALUE(tf_data->gas[id - offset].rho);
					case BlockType::Temp:
						RETURN_NORM_VALUE(tf_data->gas[id - offset].temp);
					case BlockType::Metals:
						RETURN_NORM_VALUE(tf_data->gas[id - offset].metals);
					}
					break;
				}
				case ParticleType::Dark: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_NORM_VECTOR3(tf_data->darks[id - offset].pos);
					case BlockType::Mass:
						RETURN_NORM_VALUE(tf_data->darks[id - offset].mass);
					case BlockType::Vel:
						RETURN_NORM_VECTOR3(tf_data->darks[id - offset].vel);
					case BlockType::Soft:
						RETURN_NORM_VALUE(tf_data->darks[id - offset].eps);
					case BlockType::Phi:
						RETURN_NORM_VALUE(tf_data->darks[id - offset].phi);
					}
					break;
				}
				case ParticleType::Stars: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_NORM_VECTOR3(tf_data->stars[id - offset].pos);
					case BlockType::Mass:
						RETURN_NORM_VALUE(tf_data->stars[id - offset].mass);
					case BlockType::Vel:
						RETURN_NORM_VECTOR3(tf_data->stars[id - offset].vel);
					case BlockType::Soft:
						RETURN_NORM_VALUE(tf_data->stars[id - offset].eps);
					case BlockType::Phi:
						RETURN_NORM_VALUE(tf_data->stars[id - offset].phi);
					case BlockType::Metals:
						RETURN_NORM_VALUE(tf_data->stars[id - offset].metals);
					case BlockType::TForm:
						RETURN_NORM_VALUE(tf_data->stars[id - offset].tform);
					}
					break;
				}
				}

				if (blocknr >= changa::tipsy::BlockType::BTMax && tf_data->aux_values.size() > blocknr - changa::tipsy::BlockType::BTMax) {

					Tipsy::PartialTipsyFile::AuxValue& av = tf_data->aux_values[blocknr - changa::tipsy::BlockType::BTMax];

					if (pt == ParticleType::Gas) {
						if (av.num_components == 1) {
							RETURN_NORM_VALUE(av.gas[id - offset]);
						}

						if (av.num_components == 3) {
							float v[3];
							v[0] = av.gas[(id - offset) * 3 + 0];
							v[1] = av.gas[(id - offset) * 3 + 1];
							v[2] = av.gas[(id - offset) * 3 + 2];
							RETURN_NORM_VECTOR3(v);
						}
					}

					if (pt == ParticleType::Dark) {
						if (av.num_components == 1) {
							RETURN_NORM_VALUE(av.darks[id - offset]);
						}

						if (av.num_components == 3) {
							float v[3];
							v[0] = av.darks[(id - offset) * 3 + 0];
							v[1] = av.darks[(id - offset) * 3 + 1];
							v[2] = av.darks[(id - offset) * 3 + 2];
							RETURN_NORM_VECTOR3(v);
						}
					}

					if (pt == ParticleType::Stars) {
						if (av.num_components == 1) {
							RETURN_NORM_VALUE(av.stars[id - offset]);
						}

						if (av.num_components == 3) {
							float v[3];
							v[0] = av.stars[(id - offset) * 3 + 0];
							v[1] = av.stars[(id - offset) * 3 + 1];
							v[2] = av.stars[(id - offset) * 3 + 2];
							RETURN_NORM_VECTOR3(v);
						}
					}
				}

				RETURN_NORM_EMPTY;
			}

			int get_particle_value(int blocknr, uint64_t id, float* out_value) {
				ParticleType pt = (ParticleType)get_particle_type(id);
				uint64_t offset = get_particle_type_offset(pt);

				BlockType bt = (BlockType)blocknr;

				switch (pt) {
				case ParticleType::Gas: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_ORIG_VECTOR3(tf_data->gas[id - offset].pos);
					case BlockType::Mass:
						RETURN_ORIG_VALUE(tf_data->gas[id - offset].mass);
					case BlockType::Vel:
						RETURN_ORIG_VECTOR3(tf_data->gas[id - offset].vel);
					case BlockType::Phi:
						RETURN_ORIG_VALUE(tf_data->gas[id - offset].phi);
					case BlockType::HSmooth:
						RETURN_ORIG_VALUE(tf_data->gas[id - offset].hsmooth);
					case BlockType::Rho:
						RETURN_ORIG_VALUE(tf_data->gas[id - offset].rho);
					case BlockType::Temp:
						RETURN_ORIG_VALUE(tf_data->gas[id - offset].temp);
					case BlockType::Metals:
						RETURN_ORIG_VALUE(tf_data->gas[id - offset].metals);
					}
					break;
				}
				case ParticleType::Dark: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_ORIG_VECTOR3(tf_data->darks[id - offset].pos);
					case BlockType::Mass:
						RETURN_ORIG_VALUE(tf_data->darks[id - offset].mass);
					case BlockType::Vel:
						RETURN_ORIG_VECTOR3(tf_data->darks[id - offset].vel);
					case BlockType::Soft:
						RETURN_ORIG_VALUE(tf_data->darks[id - offset].eps);
					case BlockType::Phi:
						RETURN_ORIG_VALUE(tf_data->darks[id - offset].phi);
					}
					break;
				}
				case ParticleType::Stars: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_ORIG_VECTOR3(tf_data->stars[id - offset].pos);
					case BlockType::Mass:
						RETURN_ORIG_VALUE(tf_data->stars[id - offset].mass);
					case BlockType::Vel:
						RETURN_ORIG_VECTOR3(tf_data->stars[id - offset].vel);
					case BlockType::Soft:
						RETURN_ORIG_VALUE(tf_data->stars[id - offset].eps);
					case BlockType::Phi:
						RETURN_ORIG_VALUE(tf_data->stars[id - offset].phi);
					case BlockType::Metals:
						RETURN_ORIG_VALUE(tf_data->stars[id - offset].metals);
					case BlockType::TForm:
						RETURN_ORIG_VALUE(tf_data->stars[id - offset].tform);
					}
					break;
				}
				}

				if (blocknr >= changa::tipsy::BlockType::BTMax && tf_data->aux_values.size() > blocknr - changa::tipsy::BlockType::BTMax) {

					Tipsy::PartialTipsyFile::AuxValue& av = tf_data->aux_values[blocknr - changa::tipsy::BlockType::BTMax];

					if (pt == ParticleType::Gas) {
						if (av.num_components == 1) {
							RETURN_ORIG_VALUE(av.gas[id - offset]);
						}

						if (av.num_components == 3) {
							float v[3];
							v[0] = av.gas[(id - offset) * 3 + 0];
							v[1] = av.gas[(id - offset) * 3 + 1];
							v[2] = av.gas[(id - offset) * 3 + 2];
							RETURN_ORIG_VECTOR3(v);
						}
					}

					if (pt == ParticleType::Dark) {
						if (av.num_components == 1) {
							RETURN_ORIG_VALUE(av.darks[id - offset]);
						}

						if (av.num_components == 3) {
							float v[3];
							v[0] = av.darks[(id - offset) * 3 + 0];
							v[1] = av.darks[(id - offset) * 3 + 1];
							v[2] = av.darks[(id - offset) * 3 + 2];
							RETURN_ORIG_VECTOR3(v);
						}
					}

					if (pt == ParticleType::Stars) {
						if (av.num_components == 1) {
							RETURN_ORIG_VALUE(av.stars[id - offset]);
						}

						if (av.num_components == 3) {
							float v[3];
							v[0] = av.stars[(id - offset) * 3 + 0];
							v[1] = av.stars[(id - offset) * 3 + 1];
							v[2] = av.stars[(id - offset) * 3 + 2];
							RETURN_ORIG_VECTOR3(v);
						}
					}
				}

				RETURN_ORIG_EMPTY;
			}

			int get_particle_value_comp(int blocknr, uint64_t id) {
				ParticleType pt = (ParticleType)get_particle_type(id);
				uint64_t offset = get_particle_type_offset(pt);

				BlockType bt = (BlockType)blocknr;

				switch (pt) {
				case ParticleType::Gas: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_COMP_VECTOR3(tf_data->gas[id - offset].pos);
					case BlockType::Mass:
						RETURN_COMP_VALUE(tf_data->gas[id - offset].mass);
					case BlockType::Vel:
						RETURN_COMP_VECTOR3(tf_data->gas[id - offset].vel);
					case BlockType::Phi:
						RETURN_COMP_VALUE(tf_data->gas[id - offset].phi);
					case BlockType::HSmooth:
						RETURN_COMP_VALUE(tf_data->gas[id - offset].hsmooth);
					case BlockType::Rho:
						RETURN_COMP_VALUE(tf_data->gas[id - offset].rho);
					case BlockType::Temp:
						RETURN_COMP_VALUE(tf_data->gas[id - offset].temp);
					case BlockType::Metals:
						RETURN_COMP_VALUE(tf_data->gas[id - offset].metals);
					}
					break;
				}
				case ParticleType::Dark: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_COMP_VECTOR3(tf_data->darks[id - offset].pos);
					case BlockType::Mass:
						RETURN_COMP_VALUE(tf_data->darks[id - offset].mass);
					case BlockType::Vel:
						RETURN_COMP_VECTOR3(tf_data->darks[id - offset].vel);
					case BlockType::Soft:
						RETURN_COMP_VALUE(tf_data->darks[id - offset].eps);
					case BlockType::Phi:
						RETURN_COMP_VALUE(tf_data->darks[id - offset].phi);
					}
					break;
				}
				case ParticleType::Stars: {

					switch (bt) {
					case BlockType::Pos:
						RETURN_COMP_VECTOR3(tf_data->stars[id - offset].pos);
					case BlockType::Mass:
						RETURN_COMP_VALUE(tf_data->stars[id - offset].mass);
					case BlockType::Vel:
						RETURN_COMP_VECTOR3(tf_data->stars[id - offset].vel);
					case BlockType::Soft:
						RETURN_COMP_VALUE(tf_data->stars[id - offset].eps);
					case BlockType::Phi:
						RETURN_COMP_VALUE(tf_data->stars[id - offset].phi);
					case BlockType::Metals:
						RETURN_COMP_VALUE(tf_data->stars[id - offset].metals);
					case BlockType::TForm:
						RETURN_COMP_VALUE(tf_data->stars[id - offset].tform);
					}
					break;
				}
				}

				if (blocknr >= changa::tipsy::BlockType::BTMax && tf_data->aux_values.size() > blocknr - changa::tipsy::BlockType::BTMax) {

					Tipsy::PartialTipsyFile::AuxValue& av = tf_data->aux_values[blocknr - changa::tipsy::BlockType::BTMax];

					if (pt == ParticleType::Gas) {
						if (av.num_components == 1) {
							RETURN_COMP_VALUE(av.gas[id - offset]);
						}

						if (av.num_components == 3) {
							RETURN_COMP_VECTOR3(v);
						}
					}

					if (pt == ParticleType::Dark) {
						if (av.num_components == 1) {
							RETURN_COMP_VALUE(av.darks[id - offset]);
						}

						if (av.num_components == 3) {
							RETURN_COMP_VECTOR3(v);
						}
					}

					if (pt == ParticleType::Stars) {
						if (av.num_components == 1) {
							RETURN_COMP_VALUE(av.stars[id - offset]);
						}

						if (av.num_components == 3) {
							RETURN_COMP_VECTOR3(v);
						}
					}
				}

				RETURN_COMP_EMPTY;
			}

			void get_types_and_blocks(std::vector<int>& types_and_blocks) {
				// 3 * 10
				types_and_blocks.resize(ParticleType::PTMax * (BlockType::BTMax + tf_data->aux_values.size()));
				memset(types_and_blocks.data(), 0, types_and_blocks.size() * sizeof(int));
				//memcpy(types_and_blocks.data(), all_particles_blocks, sizeof(all_particles_blocks));

				int* tab = types_and_blocks.data();

				for (int ipt = 0; ipt < ParticleType::PTMax; ipt++) {
					// tipsy default
					for (int ibt = 0; ibt < BlockType::BTMax; ibt++) {

						ParticleType pt = (ParticleType)ipt;
						BlockType bt = (BlockType)ibt;

						switch (pt) {
						case ParticleType::Gas: {

							switch (bt) {
							case BlockType::Pos:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Mass:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Vel:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Phi:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::HSmooth:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Rho:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Temp:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Metals:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							}
							break;
						}
						case ParticleType::Dark: {

							switch (bt) {
							case BlockType::Pos:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Mass:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Vel:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Soft:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Phi:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							}
							break;
						}
						case ParticleType::Stars: {

							switch (bt) {
							case BlockType::Pos:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Mass:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Vel:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Soft:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Phi:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::Metals:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							case BlockType::TForm:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							}
							break;
						}
						}
					}

					// tipsy aux
					for (int ibt = BlockType::BTMax; ibt < BlockType::BTMax + tf_data->aux_values.size(); ibt++) {
						ParticleType pt = (ParticleType)ipt;

						switch (pt) {
						case ParticleType::Gas: {
							tab[ibt * ParticleType::PTMax + ipt]++;
							break;
						}
						case ParticleType::Dark: {
							tab[ibt * ParticleType::PTMax + ipt]++;
							break;
						}
						case ParticleType::Stars: {
							tab[ibt * ParticleType::PTMax + ipt]++;
							break;
						}
						}
					}
				}
			}

			std::string get_dataset_name(int blocknr) {
				changa::tipsy::BlockType bt = (changa::tipsy::BlockType)blocknr;
				switch (bt) {
				case changa::tipsy::BlockType::Pos:
					return "Pos";
				case changa::tipsy::BlockType::Mass:
					return "Mass";
				case changa::tipsy::BlockType::Vel:
					return "Vel";
				case changa::tipsy::BlockType::Soft:
					return "Soft";
				case changa::tipsy::BlockType::Phi:
					return "Phi";
				case changa::tipsy::BlockType::HSmooth:
					return "HSmooth";
				case changa::tipsy::BlockType::Rho:
					return "Rho";
				case changa::tipsy::BlockType::Temp:
					return "Temp";
				case changa::tipsy::BlockType::Metals:
					return "Metals";
				case changa::tipsy::BlockType::TForm:
					return "TForm";
				}

				if (blocknr >= changa::tipsy::BlockType::BTMax && tf_data->aux_values.size() > blocknr - changa::tipsy::BlockType::BTMax) {
					return tf_data->aux_values[blocknr - changa::tipsy::BlockType::BTMax].name;
				}

				return "Unknown";
			}

			int get_particle_type(uint64_t id) {
				if (id < tf_data->gas.size()) {
					return ParticleType::Gas;
				}
				else if (id < tf_data->gas.size() + tf_data->darks.size()) {
					return ParticleType::Dark;
				}
				else {
					return ParticleType::Stars;
				}
			}

			void get_particle_position(uint64_t id, double* pos) {
				ParticleType pt = (ParticleType)get_particle_type(id);
				uint64_t offset = get_particle_type_offset(pt);

				switch (pt) {
				case ParticleType::Gas: {
					pos[0] = tf_data->gas[id - offset].pos[0];
					pos[1] = tf_data->gas[id - offset].pos[1];
					pos[2] = tf_data->gas[id - offset].pos[2];
					break;
				}
				case ParticleType::Dark: {
					pos[0] = tf_data->darks[id - offset].pos[0];
					pos[1] = tf_data->darks[id - offset].pos[1];
					pos[2] = tf_data->darks[id - offset].pos[2];
					break;
				}
				case ParticleType::Stars: {
					pos[0] = tf_data->stars[id - offset].pos[0];
					pos[1] = tf_data->stars[id - offset].pos[1];
					pos[2] = tf_data->stars[id - offset].pos[2];
					break;
				}
				}
			}

			size_t get_local_num_particles() {
				return tf_data->h.nbodies; //tf_data->gas.size() + tf_data->darks.size() + tf_data->stars.size();
			}
			size_t get_global_num_particles() {
				return tf_data->fullHeader.nbodies; //tf_data->gas.size() + tf_data->darks.size() + tf_data->stars.size();
			}

			//double get_particle_radius(uint64_t id) {
			//	return get_particle_hsml(id) * 2.0;
			//}

			double get_particle_hsml(uint64_t id) {
				if (g_smoothlength_blocknr != -1)
					return get_particle_norm_value(g_smoothlength_blocknr, id);

				return 0.0;
			}

			double get_particle_mass(uint64_t id) {
				return get_particle_norm_value(BlockType::Mass, id);
			}

			int get_particle_rho_blocknr() {
				return BlockType::Rho;
			}

			double get_particle_rho(uint64_t id) {
				return get_particle_norm_value(BlockType::Rho, id);
			}

			std::string get_particle_unit(int blocknr) {
				return "";
			}
		}
	}
}