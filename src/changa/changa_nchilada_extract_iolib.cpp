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

#include "changa_nchilada_extract_iolib.h"

#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "tree_xdr.h"
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

namespace changa {
	namespace nchilada {
		namespace io {

			static double fh_time; // gross, but quick way to get time
			double steps_time[2];

			//Tipsy64::TipsyFile64* tf_data = nullptr;
			struct NChiladaBlockData {
				void* data = nullptr;
				FieldHeader fh = {};

				uint64_t start_particle_rank = 0;
				uint64_t num_particles_per_rank = 0;
			};

			struct NChiladaData {
				//std::map<ParticleType, std::map<BlockType, NChiladaBlockData>> data_map;
				std::vector<std::vector<NChiladaBlockData>> data_map;

				//uint64_t nSph = 0, nDark = 0, nStar = 0;
				NChiladaData() {
					data_map.resize(ParticleType::PTMax);

					for (int ipt = 0; ipt < ParticleType::PTMax; ipt++) {
						data_map[ipt].resize(BlockType::BTMax);
					}
				}
				~NChiladaData() {
					for (auto& particlePair : data_map) {
						for (auto& blockPair : particlePair) {
							deleteField(blockPair.fh, blockPair.data);
							blockPair.data = nullptr;
						}
						particlePair.clear();
					}
					data_map.clear();
				}

				uint64_t get_size_global(ParticleType pt, BlockType bt) const {
					if(data_map.size() <= pt || data_map[pt].size() <= bt)
						return 0;

					return data_map[pt][bt].fh.numParticles;
				}
				uint64_t get_size(ParticleType pt, BlockType bt) const {
					if(data_map.size() <= pt || data_map[pt].size() <= bt)
						return 0;

					return data_map[pt][bt].num_particles_per_rank;
				}
				bool get_value(ParticleType pt, BlockType bt, uint64_t id, float& value) const {
					if(data_map.size() <= pt || data_map[pt].size() <= bt)
						return false;

					if (data_map[pt][bt].data != nullptr) {
						switch (data_map[pt][bt].fh.code) {
						case TypeHandling::float32:
							value = static_cast<float*>(data_map[pt][bt].data)[data_map[pt][bt].fh.dimensions * id];
							break;
						case TypeHandling::float64:
							value = static_cast<double*>(data_map[pt][bt].data)[data_map[pt][bt].fh.dimensions * id];
							break;
						default:
							throw XDRException("I don't recognize the type of this field!");
						}

						return true;

					}
					return false;
				}
				bool get_vector3(ParticleType pt, BlockType bt, uint64_t id, float* value) const {
					if(data_map.size() <= pt || data_map[pt].size() <= bt)
						return false;
											
					if (data_map[pt][bt].data != nullptr) {

						if (data_map[pt][bt].fh.dimensions != 3)
							throw XDRException("Wrong dimensions in get_vector3");

						switch (data_map[pt][bt].fh.code) {
						case TypeHandling::float32:
							for (unsigned int j = 0; j < data_map[pt][bt].fh.dimensions; ++j)
								value[j] = static_cast<float*>(data_map[pt][bt].data)[data_map[pt][bt].fh.dimensions * id + j];

							break;
						case TypeHandling::float64:
							for (unsigned int j = 0; j < data_map[pt][bt].fh.dimensions; ++j)
								value[j] = static_cast<double*>(data_map[pt][bt].data)[data_map[pt][bt].fh.dimensions * id + j];
							break;
						default:
							throw XDRException("I don't recognize the type of this field!");
						}

						return true;

					}
					return false;
				}

			};
			NChiladaData* nch_data = nullptr;

			void readFieldData(std::string filename, NChiladaBlockData& bdata, int world_rank, int world_size)
			{
				FILE* infile = fopen(filename.c_str(), "rb");
				if (!infile) {
					std::string smess("Couldn't open field file: ");
					smess += filename;
					throw XDRException(smess);
				}

				XDR xdrs;
				xdrstdio_create(&xdrs, infile, XDR_DECODE);

				if (!xdr_template(&xdrs, &bdata.fh)) {
					throw XDRException("Couldn't read header from file!");
				}
				if (bdata.fh.magic != FieldHeader::MagicNumber) {
					throw XDRException("This file does not appear to be a field file (magic number doesn't match).");
				}
				//if (fh.dimensions != dim) {
				//	throw XDRException("Wrong dimension of positions.");
				//}

				//u_int64_t numParticles = 0; u_int64_t startParticle = 0; numParticles = fh.numParticles;

				// Calculate the range of particles each rank will process
				uint64_t particles_per_rank = bdata.fh.numParticles / world_size;
				uint64_t start_index = world_rank * particles_per_rank;
				uint64_t end_index = (world_rank + 1) * particles_per_rank;

				// Handle any remaining particles
				if (world_rank == world_size - 1) {
					end_index = bdata.fh.numParticles;
				}

				bdata.start_particle_rank = start_index;
				bdata.num_particles_per_rank = end_index - start_index;

				bdata.data = readField(bdata.fh, &xdrs, bdata.num_particles_per_rank, bdata.start_particle_rank);

				if (bdata.data == 0) {
					throw XDRException("Had problems reading in the field");
				}
				xdr_destroy(&xdrs);
				fclose(infile);
				//return data;
			}

			void getData(NChiladaBlockData& bdata, std::string basedir, std::string typedir, std::string type, int world_rank, int world_size)
			{
				std::string filename = basedir + "/" + typedir + "/" + type;

#ifdef _WIN32
				if (_access(filename.c_str(), 0) == 0)
#else
				if (access(filename.c_str(), F_OK) == 0)
#endif
					readFieldData(filename.c_str(), bdata, world_rank, world_size);
				//deleteField(fh, data);
			}


			void init_lib(std::string basedir, int world_rank, int world_size) {
				steps_time[0] = omp_get_wtime();

				if (nch_data != nullptr)
					delete nch_data;

				//tf_data = new Tipsy64::TipsyFile64("tst", nSph, nDark, nStar);
				nch_data = new NChiladaData();

				//tf_data->h.time = fh_time;
				//Pos
				getData(nch_data->data_map[ParticleType::Dark][BlockType::Pos], basedir, "/dark", "/pos", world_rank, world_size);
				//getMass(bdata_dark, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Dark][BlockType::Mass], basedir, "/dark", "/mass", world_rank, world_size);
				//getVel(bdata_dark, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Dark][BlockType::Vel], basedir, "/dark", "/vel", world_rank, world_size);
				//getSoft(bdata_dark, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Dark][BlockType::Soft], basedir, "/dark", "/soft", world_rank, world_size);
				//getPhi(bdata_dark, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Dark][BlockType::Phi], basedir, "/dark", "/phi", world_rank, world_size);


				//if (nSph > 0) {
				//strncpy(filename, basedir.c_str(), FILELEN, world_rank, world_size);
				//strcat(filename, "/gas", world_rank, world_size);
				//getPos(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::Pos], basedir, "/dark", "/pos", world_rank, world_size);
				//getMass(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::Mass], basedir, "/dark", "/mass", world_rank, world_size);
				//getVel(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::Vel], basedir, "/dark", "/vel", world_rank, world_size);
				//getPhi(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::Phi], basedir, "/dark", "/phi", world_rank, world_size);
				//getHSmooth(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::HSmooth], basedir, "/dark", "/soft", world_rank, world_size);
				//getRho(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::Rho], basedir, "/dark", "/GasDensity", world_rank, world_size);
				//getTemp(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::Temp], basedir, "/dark", "/temperature", world_rank, world_size);
				//getMetalsOx(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::MetalsOx], basedir, "/dark", "/OxMassFrac", world_rank, world_size);
				//getMetalsFe(tf_data->gas, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Gas][BlockType::MetalsFe], basedir, "/dark", "/FeMassFrac", world_rank, world_size);
				//}

				//if (nStar > 0) {
				//	strncpy(filename, basedir.c_str(), FILELEN, world_rank, world_size);
				//	strcat(filename, "/star", world_rank, world_size);
				//getPos(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::Pos], basedir, "/star", "/pos", world_rank, world_size);
				//	getMass(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::Mass], basedir, "/star", "/mass", world_rank, world_size);
				//	getVel(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::Vel], basedir, "/star", "/vel", world_rank, world_size);
				//	getSoft(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::Soft], basedir, "/star", "/soft", world_rank, world_size);
				//	getPhi(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::Phi], basedir, "/star", "/phi", world_rank, world_size);
				//	getMetalsOx(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::MetalsOx], basedir, "/star", "/OxMassFrac", world_rank, world_size);
				//	getMetalsFe(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::MetalsFe], basedir, "/star", "/FeMassFrac", world_rank, world_size);
				//getTForm(tf_data->stars, filename, world_rank, world_size);
				getData(nch_data->data_map[ParticleType::Stars][BlockType::TForm], basedir, "/star", "/timeform", world_rank, world_size);

				//}

				steps_time[1] = omp_get_wtime();
			}

			void finish_lib() {
				if (nch_data != nullptr)
					delete nch_data;

				nch_data = nullptr;
			}

			uint64_t get_particle_type_offset(ParticleType pt) {
				switch (pt) {
				case ParticleType::Gas:
					return 0;
				case ParticleType::Dark:
					return nch_data->get_size(ParticleType::Gas, BlockType::Pos);
				case ParticleType::Stars:
					return nch_data->get_size(ParticleType::Gas, BlockType::Pos) + nch_data->get_size(ParticleType::Dark, BlockType::Pos);
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
				float value;
				float vector3[3];

				switch (pt) {
				case ParticleType::Gas: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						nch_data->get_vector3(pt, bt, id - offset, vector3);
						RETURN_NORM_VECTOR3(vector3);
					case BlockType::Mass:
					case BlockType::Phi:
					case BlockType::HSmooth:
					case BlockType::Rho:
					case BlockType::Temp:
					case BlockType::MetalsOx:
					case BlockType::MetalsFe:
						nch_data->get_value(pt, bt, id - offset, value);
						RETURN_NORM_VALUE(value);

					}
					break;
				}
				case ParticleType::Dark: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						nch_data->get_vector3(pt, bt, id - offset, vector3);
						RETURN_NORM_VECTOR3(vector3);
					case BlockType::Soft:
					case BlockType::Phi:
						nch_data->get_value(pt, bt, id - offset, value);
						RETURN_NORM_VALUE(value);
					}
					break;
				}
				case ParticleType::Stars: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						nch_data->get_vector3(pt, bt, id - offset, vector3);
						RETURN_NORM_VECTOR3(vector3);
					case BlockType::Mass:
					case BlockType::Soft:
					case BlockType::Phi:
					case BlockType::MetalsOx:
					case BlockType::MetalsFe:
					case BlockType::TForm:
						nch_data->get_value(pt, bt, id - offset, value);
						RETURN_NORM_VALUE(value);
					}
					break;
				}
				}

				RETURN_NORM_EMPTY;
			}

			int get_particle_value(int blocknr, uint64_t id, float* out_value) {
				ParticleType pt = (ParticleType)get_particle_type(id);
				uint64_t offset = get_particle_type_offset(pt);

				BlockType bt = (BlockType)blocknr;
				float value;
				float vector3[3];

				switch (pt) {
				case ParticleType::Gas: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						nch_data->get_vector3(pt, bt, id - offset, vector3);
						RETURN_ORIG_VECTOR3(vector3);
					case BlockType::Mass:
					case BlockType::Phi:
					case BlockType::HSmooth:
					case BlockType::Rho:
					case BlockType::Temp:
					case BlockType::MetalsOx:
					case BlockType::MetalsFe:
						nch_data->get_value(pt, bt, id - offset, value);
						RETURN_ORIG_VALUE(value);

					}
					break;
				}
				case ParticleType::Dark: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						nch_data->get_vector3(pt, bt, id - offset, vector3);
						RETURN_ORIG_VECTOR3(vector3);
					case BlockType::Soft:
					case BlockType::Phi:
						nch_data->get_value(pt, bt, id - offset, value);
						RETURN_ORIG_VALUE(value);
					}
					break;
				}
				case ParticleType::Stars: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						nch_data->get_vector3(pt, bt, id - offset, vector3);
						RETURN_ORIG_VECTOR3(vector3);
					case BlockType::Mass:
					case BlockType::Soft:
					case BlockType::Phi:
					case BlockType::MetalsOx:
					case BlockType::MetalsFe:
					case BlockType::TForm:
						nch_data->get_value(pt, bt, id - offset, value);
						RETURN_ORIG_VALUE(value);
					}
					break;
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
					case BlockType::Vel:
						RETURN_COMP_VECTOR3(vector3);
					case BlockType::Mass:
					case BlockType::Phi:
					case BlockType::HSmooth:
					case BlockType::Rho:
					case BlockType::Temp:
					case BlockType::MetalsOx:
					case BlockType::MetalsFe:
						RETURN_COMP_VALUE(value);

					}
					break;
				}
				case ParticleType::Dark: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						RETURN_COMP_VECTOR3(vector3);
					case BlockType::Soft:
					case BlockType::Phi:
						RETURN_COMP_VALUE(value);
					}
					break;
				}
				case ParticleType::Stars: {

					switch (bt) {
					case BlockType::Pos:
					case BlockType::Vel:
						RETURN_COMP_VECTOR3(vector3);
					case BlockType::Mass:
					case BlockType::Soft:
					case BlockType::Phi:
					case BlockType::MetalsOx:
					case BlockType::MetalsFe:
					case BlockType::TForm:
						RETURN_COMP_VALUE(value);
					}
					break;
				}
				}

				RETURN_COMP_EMPTY;
			}

			void get_types_and_blocks(std::vector<int>& types_and_blocks) {
				// 3 * 10
				types_and_blocks.resize(ParticleType::PTMax * BlockType::BTMax);
				memset(types_and_blocks.data(), 0, types_and_blocks.size() * sizeof(int));
				//memcpy(types_and_blocks.data(), all_particles_blocks, sizeof(all_particles_blocks));

				int* tab = types_and_blocks.data();

				for (int ipt = 0; ipt < ParticleType::PTMax; ipt++) {
					for (int ibt = 0; ibt < BlockType::BTMax; ibt++) {

						ParticleType pt = (ParticleType)ipt;
						BlockType bt = (BlockType)ibt;

						switch (pt) {
						case ParticleType::Gas: {

							switch (bt) {
							case BlockType::Pos:
							case BlockType::Mass:
							case BlockType::Vel:
							case BlockType::Phi:
							case BlockType::HSmooth:
							case BlockType::Rho:
							case BlockType::Temp:
							case BlockType::MetalsOx:
							case BlockType::MetalsFe:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							}

							break;
						}
						case ParticleType::Dark: {

							switch (bt) {
							case BlockType::Pos:
							case BlockType::Mass:
							case BlockType::Vel:
							case BlockType::Soft:
							case BlockType::Phi:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							}

							break;
						}
						case ParticleType::Stars: {

							switch (bt) {
							case BlockType::Pos:
							case BlockType::Mass:
							case BlockType::Vel:
							case BlockType::Soft:
							case BlockType::Phi:
							case BlockType::MetalsOx:
							case BlockType::MetalsFe:
							case BlockType::TForm:
								tab[ibt * ParticleType::PTMax + ipt]++;
								break;
							}

							break;
						}
						}
					}
				}
			}

			int get_particle_type(uint64_t id) {
				if (id < nch_data->get_size(ParticleType::Gas, BlockType::Pos)) {
					return ParticleType::Gas;
				}
				else if (id < nch_data->get_size(ParticleType::Gas, BlockType::Pos) + nch_data->get_size(ParticleType::Dark, BlockType::Pos)) {
					return ParticleType::Dark;
				}
				else {
					return ParticleType::Stars;
				}
			}

			void get_particle_position(uint64_t id, double* pos) {
				ParticleType pt = (ParticleType)get_particle_type(id);
				uint64_t offset = get_particle_type_offset(pt);

				float vector3[3];
				nch_data->get_vector3(pt, BlockType::Pos, id - offset, vector3);

				pos[0] = vector3[0];
				pos[1] = vector3[1];
				pos[2] = vector3[2];
			}

			size_t get_local_num_particles() {
				return nch_data->get_size(ParticleType::Gas, BlockType::Pos) + nch_data->get_size(ParticleType::Dark, BlockType::Pos) + nch_data->get_size(ParticleType::Stars, BlockType::Pos);
			}
			size_t get_global_num_particles() {
				return nch_data->get_size_global(ParticleType::Gas, BlockType::Pos) + nch_data->get_size_global(ParticleType::Dark, BlockType::Pos) + nch_data->get_size_global(ParticleType::Stars, BlockType::Pos);
			}

			//double get_particle_radius(uint64_t id) {
			//	return 0;
			//}

			double get_particle_hsml(uint64_t id) {
				return 0;
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