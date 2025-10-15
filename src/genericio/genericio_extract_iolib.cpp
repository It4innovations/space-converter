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

#include "genericio_extract_iolib.h"

#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>

#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include <omp.h>

#include "convert_common.h"
#include "GenericIO.h"

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

namespace genericio {
	namespace io {

		//class GIODataBase {
		//public:
		//	virtual ~GIODataBase() {}
		//	virtual float get_fvalue(size_t i) = 0;
		//	virtual std::string get_name() = 0;
		//	virtual size_t get_count() = 0;
		//	virtual size_t get_nelem() = 0;
		//};

		//template <class T>
		class GIOData /* : public GIODataBase*/ {
		public:
			enum GIODataType {
				Float,
				Double,
				UnsignedChar,
				SignedChar,
				Int16,
				UInt16,
				Int32,
				UInt32,
				Int64,
				UInt64
			};

		public:
			GIOData(): num_elements(0), type_element(GIODataType::Float), data_name(""), data_size(0), extra_data_space(0) {}
			GIOData(gio::GenericIO& G, size_t c, size_t nelem, GIODataType dt, const std::string& n)
				: num_elements(nelem), type_element(dt), data_name(n), extra_data_space(0)
			{
				extra_data_space = G.requestedExtraSpace();
				size_t alloc_size = c * nelem * get_size_of_elem(dt) + extra_data_space;
				//printf("alloc_size: %lld\n", alloc_size);
				data.resize(alloc_size);
				data_size = c * nelem * get_size_of_elem(dt);

				switch (type_element) {
				case GIODataType::Float:
					G.addScalarizedVariable(n, (float*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::Double:
					G.addScalarizedVariable(n, (double*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::UnsignedChar:
					G.addScalarizedVariable(n, (unsigned char*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::SignedChar:
					G.addScalarizedVariable(n, (signed char*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::Int16:
					G.addScalarizedVariable(n, (int16_t*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::UInt16:
					G.addScalarizedVariable(n, (uint16_t*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::Int32:
					G.addScalarizedVariable(n, (int32_t*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::UInt32:
					G.addScalarizedVariable(n, (uint32_t*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::Int64:
					G.addScalarizedVariable(n, (int64_t*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				case GIODataType::UInt64:
					G.addScalarizedVariable(n, (uint64_t*)data.data(), nelem, gio::GenericIO::VarHasExtraSpace);
					break;
				}				
			}

			// ~GIOData() {
			// 	// Clear the data vector to release memory
			// 	data.clear();
			
			// 	// Reset fundamental data members
			// 	data_name.clear();
			// 	num_elements = 0;
			// 	type_element = GIODataType::Float; // Reset to a default type
			
			// 	// You can add debug logs if necessary
			// 	// std::cout << "GIOData object destroyed." << std::endl;
			// }			

			float get_fvalue(size_t i) {			
				switch (type_element) {
				case GIODataType::Float:
					return ((float*)data.data())[i];
				case GIODataType::Double:
					return ((double*)data.data())[i];
				case GIODataType::UnsignedChar:
					return ((unsigned char*)data.data())[i];
				case GIODataType::SignedChar:
					return ((signed char*)data.data())[i];
				case GIODataType::Int16:
					return ((int16_t*)data.data())[i];
				case GIODataType::UInt16:
					return ((uint16_t*)data.data())[i];
				case GIODataType::Int32:
					return ((int32_t*)data.data())[i];
				case GIODataType::UInt32:
					return ((uint32_t*)data.data())[i];
				case GIODataType::Int64:
					return ((int64_t*)data.data())[i];
				case GIODataType::UInt64:
					return ((uint64_t*)data.data())[i];
				default:
					return 0; // Return 0 for unknown types
				}
			}

			std::string get_name() {
				return data_name;
			}

			size_t get_count() {
				return (data_size)  / (num_elements * get_size_of_elem(type_element));
			}

			size_t get_nelem() {
				return num_elements;
			}

			static size_t get_size_of_elem(GIODataType dt) {
				switch (dt) {
				case GIODataType::Float:
					return sizeof(float);
				case GIODataType::Double:
					return sizeof(double);
				case GIODataType::UnsignedChar:
					return sizeof(unsigned char);
				case GIODataType::SignedChar:
					return sizeof(signed char);
				case GIODataType::Int16:
					return sizeof(int16_t);
				case GIODataType::UInt16:
					return sizeof(uint16_t);
				case GIODataType::Int32:
					return sizeof(int32_t);
				case GIODataType::UInt32:
					return sizeof(uint32_t);
				case GIODataType::Int64:
					return sizeof(int64_t);
				case GIODataType::UInt64:
					return sizeof(uint64_t);
				default:
					printf("Cannot get_size_of_elem: Return 0 for unknown types\n");
					return 0; // Return 0 for unknown types
				}
			}

			// Function to merge data from another GIOData object
			void merge_data(const GIOData& other, size_t count) {
				if (data_size > 0) {
					if (this->data_name != other.data_name) {
						printf("Cannot merge: data names do not match\n");
						return;
					}

					if (this->num_elements != other.num_elements) {
						printf("Cannot merge: Number of elements per entry does not match\n");
						return;
					}

					if (this->type_element != other.type_element) {
						printf("Cannot merge: Type of elements per entry does not match\n");
						return;
					}
				}
				else {
					this->data_name = other.data_name;
					this->num_elements = other.num_elements;
					this->type_element = other.type_element;
					//this->count = other.count;
				}

				size_t alloc_size = count * other.num_elements * get_size_of_elem(other.type_element);

				if (data_size + alloc_size > data.size()) {
					printf("Cannot merge: data_size + alloc_size > data.size()\n");
					return;
				}

				//if (offset == 0) {
					// Append the data from `other` to this object
					//data.insert(data.end(), other.data.begin(), other.data.begin() + alloc_size);
				memcpy(&data[data_size], other.data.data(), alloc_size);
//				}
//				else {
//#pragma omp parallel for
//					for (size_t i = 0; i < count; i++) {
//						switch (type_element) {
//						case GIODataType::Float:
//							((float*)data.data())[i] = ((float*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::Double:
//							((double*)data.data())[i] = ((double*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::UnsignedChar:
//							((unsigned char*)data.data())[i] = ((unsigned char*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::SignedChar:
//							((signed char*)data.data())[i] = ((signed char*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::Int16:
//							((int16_t*)data.data())[i] = ((int16_t*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::UInt16:
//							((uint16_t*)data.data())[i] = ((uint16_t*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::Int32:
//							((int32_t*)data.data())[i] = ((int32_t*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::UInt32:
//							((uint32_t*)data.data())[i] = ((uint32_t*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::Int64:
//							((int64_t*)data.data())[i] = ((int64_t*)other.data.data())[i] + offset;
//							break;
//						case GIODataType::UInt64:
//							((uint64_t*)data.data())[i] = ((uint64_t*)other.data.data())[i] + offset;
//							break;
//						}
//					}
//				}
				data_size += alloc_size;
			}

			void reserve_memory(size_t size) {
				data.resize(size);
				data_size = 0;
			}

		public:
			std::string data_name;
			size_t num_elements;
			//size_t count;
			GIODataType type_element;
			std::vector<uint8_t> data;
			size_t data_size;
			size_t extra_data_space;
		};

		//template <typename T>
		GIOData* addGIOData(gio::GenericIO::VariableInfo& V, gio::GenericIO& GIO, size_t count) {

			GIOData::GIODataType DT;

			if (V.IsFloat && V.IsSigned && V.ElementSize == sizeof(float)) {
				DT = GIOData::GIODataType::Float;
			}
			else if (V.IsFloat && V.IsSigned && V.ElementSize == sizeof(double)) {
				DT = GIOData::GIODataType::Double;
			}
			else if (!V.IsFloat && V.IsSigned && V.ElementSize == sizeof(signed char)) {
				DT = GIOData::GIODataType::SignedChar;
			}
			else if (!V.IsFloat && !V.IsSigned && V.ElementSize == sizeof(unsigned char)) {
				DT = GIOData::GIODataType::UnsignedChar;
			}
			else if (!V.IsFloat && V.IsSigned && V.ElementSize == sizeof(int16_t)) {
				DT = GIOData::GIODataType::Int16;
			}
			else if (!V.IsFloat && !V.IsSigned && V.ElementSize == sizeof(uint16_t)) {
				DT = GIOData::GIODataType::UInt16;
			}
			else if (!V.IsFloat && V.IsSigned && V.ElementSize == sizeof(int32_t)) {
				DT = GIOData::GIODataType::Int32;
			}
			else if (!V.IsFloat && !V.IsSigned && V.ElementSize == sizeof(uint32_t)) {
				DT = GIOData::GIODataType::UInt32;
			}
			else if (!V.IsFloat && V.IsSigned && V.ElementSize == sizeof(int64_t)) {
				DT = GIOData::GIODataType::Int64;
			}
			else if (!V.IsFloat && !V.IsSigned && V.ElementSize == sizeof(uint64_t)) {
				DT = GIOData::GIODataType::UInt64;
			}else{
				printf("Unrecognized type for %s\n", V.Name.c_str());
				return nullptr;
			}

			return new GIOData(GIO, count, V.Size / V.ElementSize, DT, V.Name);
		}

		std::vector<GIOData> gio_datas;
		//std::map<std::string, size_t> gio_datas_map;
		//std::vector<GIOData*> gio_pos_idx(3, nullptr);
		int64_t gio_pos_idx[3] = { -1,-1,-1 };
		//std::vector<GIOData*> gio_vel_idx(3, nullptr);
		int64_t gio_vel_idx[3] = { -1,-1,-1 };

		int64_t gio_mass_idx = -1;
		int64_t gio_rho_idx = -1;
		int64_t gio_hsml_idx = -1;

		//static double fh_time; // gross, but quick way to get time
		double steps_time[2];

		//Tipsy::TipsyFile* gio_data = nullptr;
		//gio::GenericIO* gio_data = nullptr;
		int64_t g_total_particles = 0;
		//size_t g_start_particles = 0;
		//size_t g_num_particles = 0;

		//int g_smoothlength_blocknr = -1;

		std::vector<std::string> split_string(const std::string& input, char delimiter) {
			std::vector<std::string> tokens;
			std::stringstream ss(input);
			std::string token;

			while (std::getline(ss, token, delimiter)) {
				tokens.push_back(token);
			}

			return tokens;
		}

		void read_multiple_files_per_rank(
			std::string basefile, 
			int world_rank, 
			int world_size,
			std::vector<std::string>& pos_names_vec,
			std::vector<std::string>& vel_names_vec,
			std::string& vel_name_mass,
			std::string& vel_name_rho,
			std::string& vel_name_hsml,
			gio::GenericIO& gio_data,
			int NR
		) {
			//int NR = gio_data.readNRanks();

			std::vector<gio::GenericIO::VariableInfo> variable_info;
			gio_data.getVariableInfo(variable_info);

			g_total_particles = gio_data.readTotalNumElems();

			//gio_datas.resize(variable_info.size());

			int ranks_per_process = NR / world_size;
			int start_ranks = world_rank * ranks_per_process;
			int num_ranks = (world_rank == world_size - 1) ? NR - start_ranks : ranks_per_process;

			////size_t max_nelem = gio_data.readNumElems(-1);
			size_t max_count = 0;
			size_t tot_count = 0;
			for (int i = start_ranks; i < start_ranks + num_ranks; ++i) {
				size_t count = gio_data.readNumElems(i);

				if (max_count < count)
					max_count = count;

				tot_count += count;
			}

			std::vector<GIOData*> gio_datas_temp;
			//std::vector<GIOData*> gio_pos_data_temp(3, nullptr);
			//std::vector<GIOData*> gio_vel_data_temp(3, nullptr);

			for (size_t v = 0; v < variable_info.size(); ++v) {
				GIOData* P = addGIOData(variable_info[v], gio_data, max_count);
				if (!P) continue;

				if (P->get_name() == pos_names_vec[0]) {
					gio_pos_idx[0] = v;
				}
				else if (P->get_name() == pos_names_vec[1]) {
					gio_pos_idx[1] = v;
				}
				else if (P->get_name() == pos_names_vec[2]) {
					gio_pos_idx[2] = v;
				}
				else if (P->get_name() == vel_names_vec[0]) {
					gio_vel_idx[0] = v;
				}
				else if (P->get_name() == vel_names_vec[1]) {
					gio_vel_idx[1] = v;
				}
				else if (P->get_name() == vel_names_vec[2]) {
					gio_vel_idx[2] = v;
				}
				else if (P->get_name() == vel_name_mass) {
					gio_mass_idx = v;
				}
				else if (P->get_name() == vel_name_rho) {
					gio_rho_idx = v;
				}
				else if (P->get_name() == vel_name_hsml) {
					gio_hsml_idx = v;
				}

				gio_datas_temp.push_back(P);
			}

			//int dims[3];
			//gio_data.readDims(dims);


			//gio_data.readPhysOrigin(PhysOrigin);
			//gio_data.readPhysScale(PhysScale);
			//gio_data.readCoords();

			gio_datas.resize(gio_datas_temp.size());

			for (size_t v = 0; v < gio_datas_temp.size(); ++v) {
				gio_datas[v].reserve_memory(
					/*gio_datas_temp[v]->data.size() * num_ranks*/
					gio_datas_temp[v]->num_elements * gio_datas_temp[v]->get_size_of_elem(gio_datas_temp[v]->type_element) * tot_count);
			}

			for (int i = start_ranks; i < start_ranks + num_ranks; ++i) {
				size_t count = gio_data.readNumElems(i);
				//printf("rank %d: reading count: %lld\n", world_rank, count);

				gio_data.readData(i, false);

				for (size_t v = 0; v < gio_datas_temp.size(); ++v) {
					gio_datas[v].merge_data(*gio_datas_temp[v], count/*, offset*/);
				}
			}

			printf("rank %d: has total count: %lld\n", world_rank, gio_datas[0].get_count());

			for (size_t v = 0; v < gio_datas_temp.size(); ++v) {
				delete gio_datas_temp[v];
			}

			//gio_data.readData(-1, true);
			//for (size_t k = 0; k < gio_datas.size(); ++k) {
			//	//gio_datas[k]->print(std::cout, j);
			//}

			//steps_time[1] = omp_get_wtime();
		}

		void read_part_file_per_rank(
			std::string basefile, 
			int world_rank, 
			int world_size,
			std::vector<std::string> &pos_names_vec,
			std::vector<std::string> &vel_names_vec,
			std::string &vel_name_mass,
			std::string &vel_name_rho,
			std::string &vel_name_hsml,
			gio::GenericIO &gio_data,
			int NR
		) {	
			//world_rank = 1; //TODO
			//world_size = 2;

			//int NR = gio_data.readNRanks();

			std::vector<gio::GenericIO::VariableInfo> variable_info;
			gio_data.getVariableInfo(variable_info);

			g_total_particles = gio_data.readTotalNumElems();

			//gio_datas.resize(variable_info.size());

			//int ranks_per_process = NR / world_size;
			//int start_ranks = world_rank * ranks_per_process;
			//int num_ranks = (world_rank == world_size - 1) ? NR - start_ranks : ranks_per_process;

			// More ranks than files: split each file among  world_size / NR ranks
			int ranks_per_file = world_size / NR;
			int gio_data_rank_id = world_rank / ranks_per_file;
			int local_rank = world_rank % ranks_per_file;

			////size_t max_nelem = gio_data.readNumElems(-1);
			size_t max_count = 0;
			size_t tot_count = 0;
			//for (int i = start_ranks; i < start_ranks + num_ranks; ++i) 
			//{
			size_t orig_count = gio_data.readNumElems(gio_data_rank_id);

			std::size_t chunk_size = orig_count / ranks_per_file;
			std::size_t offset = local_rank * chunk_size;
			if (local_rank == ranks_per_file - 1) {
				// Last rank gets the rest of the particles
				chunk_size = orig_count - offset;
			}

			max_count = chunk_size;
			tot_count = chunk_size;
			//}

			std::vector<GIOData*> gio_datas_temp;

			for (size_t v = 0; v < variable_info.size(); ++v) {
				GIOData* P = addGIOData(variable_info[v], gio_data, max_count);
				if (!P) continue;

				if (P->get_name() == pos_names_vec[0]) {
					gio_pos_idx[0] = v;
				}
				else if (P->get_name() == pos_names_vec[1]) {
					gio_pos_idx[1] = v;
				}
				else if (P->get_name() == pos_names_vec[2]) {
					gio_pos_idx[2] = v;
				}
				else if (P->get_name() == vel_names_vec[0]) {
					gio_vel_idx[0] = v;
				}
				else if (P->get_name() == vel_names_vec[1]) {
					gio_vel_idx[1] = v;
				}
				else if (P->get_name() == vel_names_vec[2]) {
					gio_vel_idx[2] = v;
				}
				else if (P->get_name() == vel_name_mass) {
					gio_mass_idx = v;
				}
				else if (P->get_name() == vel_name_rho) {
					gio_rho_idx = v;
				}
				else if (P->get_name() == vel_name_hsml) {
					gio_hsml_idx = v;
				}				

				gio_datas_temp.push_back(P);
			}

			//int dims[3];
			//gio_data.readDims(dims);			
			//gio_data.readPhysOrigin(PhysOrigin);
			//gio_data.readPhysScale(PhysScale);
			//gio_data.readCoords();

			gio_datas.resize(gio_datas_temp.size());

			for (size_t v = 0; v < gio_datas_temp.size(); ++v) {
				gio_datas[v].reserve_memory(
					/*gio_datas_temp[v]->data.size() * num_ranks*/ 
					gio_datas_temp[v]->num_elements * gio_datas_temp[v]->get_size_of_elem(gio_datas_temp[v]->type_element) * tot_count);
			}

			//for (int i = start_ranks; i < start_ranks + num_ranks; ++i) 
			{
				//size_t count = gio_data.readNumElems(i);
				//printf("rank %d: reading count: %lld\n", world_rank, count);

				//gio_data.readData(i, false);
				gio_data.readData(gio_data_rank_id, offset, tot_count, false);

//				int NErrs[3] = { 0, 0, 0 };
//				uint64_t TotalReadSize = 0;
//				gio_data.readData(gio_data_rank_id, offset, 0, TotalReadSize, NErrs);
//
//				int AllNErrs[3];
//#ifndef GENERICIO_NO_MPI
//				MPI_Allreduce(NErrs, AllNErrs, 3, MPI_INT, MPI_SUM, Comm);
//#else
//				AllNErrs[0] = NErrs[0]; AllNErrs[1] = NErrs[1]; AllNErrs[2] = NErrs[2];
//#endif
//
//				if (AllNErrs[0] > 0 || AllNErrs[1] > 0 || AllNErrs[2] > 0) {
//					std::stringstream ss;
//					ss << "Experienced " << AllNErrs[0] << " I/O error(s), " <<
//						AllNErrs[1] << " CRC error(s) and " << AllNErrs[2] <<
//						" decompression CRC error(s) reading: " << "OpenFileName";
//					throw std::runtime_error(ss.str());
//				}

				for (size_t v = 0; v < gio_datas_temp.size(); ++v) {
					gio_datas[v].merge_data(*gio_datas_temp[v], tot_count/*, offset*/);
				}				
			}

			printf("rank %d: has total count: %lld\n", world_rank, gio_datas[0].get_count());

			for (size_t v = 0; v < gio_datas_temp.size(); ++v) {
				delete gio_datas_temp[v];
			}			
		}

		void init_lib(
			std::string basefile, 
			int world_rank, 
			int world_size,
			std::vector<std::string>& pos_names_vec,
			std::vector<std::string>& vel_names_vec,
			std::string& vel_name_mass,
			std::string& vel_name_rho,
			std::string& vel_name_hsml
		) {
			steps_time[0] = omp_get_wtime();

			if (pos_names_vec.size() != 3) {
				pos_names_vec.resize(3);
				pos_names_vec[0] = "x";
				pos_names_vec[1] = "y";
				pos_names_vec[2] = "z";
			}

			if (vel_names_vec.size() != 3) {
				//vel_names_vec = split_string("vx,vy,vz", ',');
				vel_names_vec.resize(3);
				vel_names_vec[0] = "vx";
				vel_names_vec[1] = "vy";
				vel_names_vec[2] = "vz";
			}

#ifdef GENERICIO_NO_MPI
			unsigned int method = gio::GenericIO::FileIOPOSIX;  //gio::GenericIO::FileIOPOSIX; //gio::GenericIO::FileIOMPI;
			gio::GenericIO gio_data(basefile, method); // = new gio::GenericIO
#else
			unsigned int method = gio::GenericIO::FileIOMPI;  //gio::GenericIO::FileIOPOSIX; //gio::GenericIO::FileIOMPI;
			gio::GenericIO gio_data(MPI_COMM_SELF, basefile, method); // = new gio::GenericIO
#endif			

			gio_data.openAndReadHeader(gio::GenericIO::MismatchAllowed/*gio::GenericIO::MismatchRedistribute*/, -1, true);

			int NR = gio_data.readNRanks();

#if 0 //DEBUG
			NR = 1;
			world_rank = 0;
			world_size = 2;
#endif

			if (world_rank == 0) {
				printf("Total number of ranks: %d\n", NR);
			}

			if (NR < world_size) {
				read_part_file_per_rank(basefile, world_rank, world_size, pos_names_vec, vel_names_vec, vel_name_mass, vel_name_rho, vel_name_hsml, gio_data, NR);
			}
			else 
			{
				read_multiple_files_per_rank(basefile, world_rank, world_size, pos_names_vec, vel_names_vec, vel_name_mass, vel_name_rho, vel_name_hsml, gio_data, NR);
			}

			steps_time[1] = omp_get_wtime();
		}


		void finish_lib() {
			// if (gio_data != nullptr)
			// 	delete gio_data;

			// gio_data = nullptr;
		}

		uint64_t get_particle_type_offset(GenericIOParticleType pt) {
			switch (pt) {
			case GenericIOParticleType::HACC:
				return 0;
			}

			return 0;
		}

		void print_CPU_steps() {
			printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		float get_particle_norm_value(int blocknr, uint64_t id) {
			GenericIOParticleType pt = (GenericIOParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			GenericIOBlockType bt = (GenericIOBlockType)blocknr;

			switch (pt) {
			case GenericIOParticleType::HACC: {

				switch (bt) {
				case GenericIOBlockType::Pos:
					float pos[3];
					pos[0] = gio_datas[gio_pos_idx[0]].get_fvalue(id);
					pos[1] = gio_datas[gio_pos_idx[1]].get_fvalue(id);
					pos[2] = gio_datas[gio_pos_idx[2]].get_fvalue(id);

					RETURN_NORM_VECTOR3(pos);
				case GenericIOBlockType::Vel:
					float vel[3];
					vel[0] = gio_datas[gio_vel_idx[0]].get_fvalue(id);
					vel[1] = gio_datas[gio_vel_idx[1]].get_fvalue(id);
					vel[2] = gio_datas[gio_vel_idx[2]].get_fvalue(id);

					RETURN_NORM_VECTOR3(vel);
				}
				break;
			}
			//case GenericIOParticleType::Dark: {

			//	switch (bt) {
			//	case GenericIOBlockType::Pos:
			//		RETURN_NORM_VECTOR3(gio_data->darks[id - offset].pos);
			//	case GenericIOBlockType::Mass:
			//		RETURN_NORM_VALUE(gio_data->darks[id - offset].mass);
			//	case GenericIOBlockType::Vel:
			//		RETURN_NORM_VECTOR3(gio_data->darks[id - offset].vel);
			//	case GenericIOBlockType::Soft:
			//		RETURN_NORM_VALUE(gio_data->darks[id - offset].eps);
			//	case GenericIOBlockType::Phi:
			//		RETURN_NORM_VALUE(gio_data->darks[id - offset].phi);
			//	}
			//	break;
			//}
			//case GenericIOParticleType::Stars: {

			//	switch (bt) {
			//	case GenericIOBlockType::Pos:
			//		RETURN_NORM_VECTOR3(gio_data->stars[id - offset].pos);
			//	case GenericIOBlockType::Mass:
			//		RETURN_NORM_VALUE(gio_data->stars[id - offset].mass);
			//	case GenericIOBlockType::Vel:
			//		RETURN_NORM_VECTOR3(gio_data->stars[id - offset].vel);
			//	case GenericIOBlockType::Soft:
			//		RETURN_NORM_VALUE(gio_data->stars[id - offset].eps);
			//	case GenericIOBlockType::Phi:
			//		RETURN_NORM_VALUE(gio_data->stars[id - offset].phi);
			//	case GenericIOBlockType::Metals:
			//		RETURN_NORM_VALUE(gio_data->stars[id - offset].metals);
			//	case GenericIOBlockType::TForm:
			//		RETURN_NORM_VALUE(gio_data->stars[id - offset].tform);
			//	}
			//	break;
			//}
			}

			if (blocknr >= GenericIOBlockType::BTMax && gio_datas.size() > blocknr - GenericIOBlockType::BTMax) {

				GIOData* av = &gio_datas[blocknr - GenericIOBlockType::BTMax];

				if (pt == GenericIOParticleType::HACC) {
					if (av->get_nelem() == 1) {
						RETURN_NORM_VALUE(av->get_fvalue(id - offset));
					}

					if (av->get_nelem() == 3) {
						float v[3];
						v[0] = av->get_fvalue((id - offset) * 3 + 0);
						v[1] = av->get_fvalue((id - offset) * 3 + 1);
						v[2] = av->get_fvalue((id - offset) * 3 + 2);
						RETURN_NORM_VECTOR3(v);
					}
				}
			}

			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value) {
			GenericIOParticleType pt = (GenericIOParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			GenericIOBlockType bt = (GenericIOBlockType)blocknr;

			switch (pt) {
			case GenericIOParticleType::HACC: {

				switch (bt) {
				case GenericIOBlockType::Pos:
					float pos[3];
					pos[0] = gio_datas[gio_pos_idx[0]].get_fvalue(id);
					pos[1] = gio_datas[gio_pos_idx[1]].get_fvalue(id);
					pos[2] = gio_datas[gio_pos_idx[2]].get_fvalue(id);

					RETURN_ORIG_VECTOR3(pos);
				case GenericIOBlockType::Vel:
					float vel[3];
					vel[0] = gio_datas[gio_vel_idx[0]].get_fvalue(id);
					vel[1] = gio_datas[gio_vel_idx[1]].get_fvalue(id);
					vel[2] = gio_datas[gio_vel_idx[2]].get_fvalue(id);

					RETURN_ORIG_VECTOR3(vel);

					//	case GenericIOBlockType::Mass:
					//		RETURN_ORIG_VALUE(gio_data->gas[id - offset].mass);
					//	case GenericIOBlockType::Vel:
					//		RETURN_ORIG_VECTOR3(gio_data->gas[id - offset].vel);
					//	case GenericIOBlockType::Phi:
					//		RETURN_ORIG_VALUE(gio_data->gas[id - offset].phi);
					//	case GenericIOBlockType::HSmooth:
					//		RETURN_ORIG_VALUE(gio_data->gas[id - offset].hsmooth);
					//	case GenericIOBlockType::Rho:
					//		RETURN_ORIG_VALUE(gio_data->gas[id - offset].rho);
					//	case GenericIOBlockType::Temp:
					//		RETURN_ORIG_VALUE(gio_data->gas[id - offset].temp);
					//	case GenericIOBlockType::Metals:
					//		RETURN_ORIG_VALUE(gio_data->gas[id - offset].metals);
				}
				break;
			}
			//case GenericIOParticleType::Dark: {

			//	switch (bt) {
			//	case GenericIOBlockType::Pos:
			//		RETURN_ORIG_VECTOR3(gio_data->darks[id - offset].pos);
			//	case GenericIOBlockType::Mass:
			//		RETURN_ORIG_VALUE(gio_data->darks[id - offset].mass);
			//	case GenericIOBlockType::Vel:
			//		RETURN_ORIG_VECTOR3(gio_data->darks[id - offset].vel);
			//	case GenericIOBlockType::Soft:
			//		RETURN_ORIG_VALUE(gio_data->darks[id - offset].eps);
			//	case GenericIOBlockType::Phi:
			//		RETURN_ORIG_VALUE(gio_data->darks[id - offset].phi);
			//	}
			//	break;
			//}
			//case GenericIOParticleType::Stars: {

			//	switch (bt) {
			//	case GenericIOBlockType::Pos:
			//		RETURN_ORIG_VECTOR3(gio_data->stars[id - offset].pos);
			//	case GenericIOBlockType::Mass:
			//		RETURN_ORIG_VALUE(gio_data->stars[id - offset].mass);
			//	case GenericIOBlockType::Vel:
			//		RETURN_ORIG_VECTOR3(gio_data->stars[id - offset].vel);
			//	case GenericIOBlockType::Soft:
			//		RETURN_ORIG_VALUE(gio_data->stars[id - offset].eps);
			//	case GenericIOBlockType::Phi:
			//		RETURN_ORIG_VALUE(gio_data->stars[id - offset].phi);
			//	case GenericIOBlockType::Metals:
			//		RETURN_ORIG_VALUE(gio_data->stars[id - offset].metals);
			//	case GenericIOBlockType::TForm:
			//		RETURN_ORIG_VALUE(gio_data->stars[id - offset].tform);
			//	}
			//	break;
			//}
			}

			if (blocknr >= GenericIOBlockType::BTMax && gio_datas.size() > blocknr - GenericIOBlockType::BTMax) {

				GIOData* av = &gio_datas[blocknr - GenericIOBlockType::BTMax];

				if (pt == GenericIOParticleType::HACC) {
					if (av->get_nelem() == 1) {
						RETURN_ORIG_VALUE(av->get_fvalue(id - offset));
					}

					if (av->get_nelem() == 3) {
						float v[3];
						v[0] = av->get_fvalue((id - offset) * 3 + 0);
						v[1] = av->get_fvalue((id - offset) * 3 + 1);
						v[2] = av->get_fvalue((id - offset) * 3 + 2);
						RETURN_ORIG_VECTOR3(v);
					}
				}
			}


			//if (blocknr >= GenericIOBlockType::BTMax && get_gio_values() > blocknr - GenericIOBlockType::BTMax) {

			//	GenericIOAuxValue& av = gio_data->aux_values[blocknr - GenericIOBlockType::BTMax];

			//	if (pt == GenericIOParticleType::HACC) {
			//		if (av.num_components == 1) {
			//			RETURN_ORIG_VALUE(av.gas[id - offset]);
			//		}

			//		if (av.num_components == 3) {
			//			float v[3];
			//			v[0] = av.gas[(id - offset) * 3 + 0];
			//			v[1] = av.gas[(id - offset) * 3 + 1];
			//			v[2] = av.gas[(id - offset) * 3 + 2];
			//			RETURN_ORIG_VECTOR3(v);
			//		}
			//	}

			//	if (pt == GenericIOParticleType::Dark) {
			//		if (av.num_components == 1) {
			//			RETURN_ORIG_VALUE(av.darks[id - offset]);
			//		}

			//		if (av.num_components == 3) {
			//			float v[3];
			//			v[0] = av.darks[(id - offset) * 3 + 0];
			//			v[1] = av.darks[(id - offset) * 3 + 1];
			//			v[2] = av.darks[(id - offset) * 3 + 2];
			//			RETURN_ORIG_VECTOR3(v);
			//		}
			//	}

			//	if (pt == GenericIOParticleType::Stars) {
			//		if (av.num_components == 1) {
			//			RETURN_ORIG_VALUE(av.stars[id - offset]);
			//		}

			//		if (av.num_components == 3) {
			//			float v[3];
			//			v[0] = av.stars[(id - offset) * 3 + 0];
			//			v[1] = av.stars[(id - offset) * 3 + 1];
			//			v[2] = av.stars[(id - offset) * 3 + 2];
			//			RETURN_ORIG_VECTOR3(v);
			//		}
			//	}
			//}

			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id) {
			GenericIOParticleType pt = (GenericIOParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			GenericIOBlockType bt = (GenericIOBlockType)blocknr;

			switch (pt) {
			case GenericIOParticleType::HACC: {

				switch (bt) {
				case GenericIOBlockType::Pos:
					float pos[3];
					RETURN_COMP_VECTOR3(pos);

				case GenericIOBlockType::Vel:
					float vel[3];
					RETURN_COMP_VECTOR3(vel);
					//	case GenericIOBlockType::Mass:
					//		RETURN_COMP_VALUE(gio_data->gas[id - offset].mass);
					//	case GenericIOBlockType::Vel:
					//		RETURN_COMP_VECTOR3(gio_data->gas[id - offset].vel);
					//	case GenericIOBlockType::Phi:
					//		RETURN_COMP_VALUE(gio_data->gas[id - offset].phi);
					//	case GenericIOBlockType::HSmooth:
					//		RETURN_COMP_VALUE(gio_data->gas[id - offset].hsmooth);
					//	case GenericIOBlockType::Rho:
					//		RETURN_COMP_VALUE(gio_data->gas[id - offset].rho);
					//	case GenericIOBlockType::Temp:
					//		RETURN_COMP_VALUE(gio_data->gas[id - offset].temp);
					//	case GenericIOBlockType::Metals:
					//		RETURN_COMP_VALUE(gio_data->gas[id - offset].metals);
				}
				break;
			}
			//case GenericIOParticleType::Dark: {

			//	switch (bt) {
			//	case GenericIOBlockType::Pos:
			//		RETURN_COMP_VECTOR3(gio_data->darks[id - offset].pos);
			//	case GenericIOBlockType::Mass:
			//		RETURN_COMP_VALUE(gio_data->darks[id - offset].mass);
			//	case GenericIOBlockType::Vel:
			//		RETURN_COMP_VECTOR3(gio_data->darks[id - offset].vel);
			//	case GenericIOBlockType::Soft:
			//		RETURN_COMP_VALUE(gio_data->darks[id - offset].eps);
			//	case GenericIOBlockType::Phi:
			//		RETURN_COMP_VALUE(gio_data->darks[id - offset].phi);
			//	}
			//	break;
			//}
			//case GenericIOParticleType::Stars: {

			//	switch (bt) {
			//	case GenericIOBlockType::Pos:
			//		RETURN_COMP_VECTOR3(gio_data->stars[id - offset].pos);
			//	case GenericIOBlockType::Mass:
			//		RETURN_COMP_VALUE(gio_data->stars[id - offset].mass);
			//	case GenericIOBlockType::Vel:
			//		RETURN_COMP_VECTOR3(gio_data->stars[id - offset].vel);
			//	case GenericIOBlockType::Soft:
			//		RETURN_COMP_VALUE(gio_data->stars[id - offset].eps);
			//	case GenericIOBlockType::Phi:
			//		RETURN_COMP_VALUE(gio_data->stars[id - offset].phi);
			//	case GenericIOBlockType::Metals:
			//		RETURN_COMP_VALUE(gio_data->stars[id - offset].metals);
			//	case GenericIOBlockType::TForm:
			//		RETURN_COMP_VALUE(gio_data->stars[id - offset].tform);
			//	}
			//	break;
			//}
			}

			if (blocknr >= GenericIOBlockType::BTMax && gio_datas.size() > blocknr - GenericIOBlockType::BTMax) {

				GIOData* av = &gio_datas[blocknr - GenericIOBlockType::BTMax];

				if (pt == GenericIOParticleType::HACC) {
					return av->get_nelem();
				}
			}

			//if (blocknr >= GenericIOBlockType::BTMax && get_gio_values() > blocknr - GenericIOBlockType::BTMax) {

			//	GenericIOAuxValue& av = gio_data->aux_values[blocknr - GenericIOBlockType::BTMax];

			//	if (pt == GenericIOParticleType::HACC) {
			//		if (av.num_components == 1) {
			//			RETURN_COMP_VALUE(av.gas[id - offset]);
			//		}

			//		if (av.num_components == 3) {
			//			RETURN_COMP_VECTOR3(v);
			//		}
			//	}

			//	if (pt == GenericIOParticleType::Dark) {
			//		if (av.num_components == 1) {
			//			RETURN_COMP_VALUE(av.darks[id - offset]);
			//		}

			//		if (av.num_components == 3) {
			//			RETURN_COMP_VECTOR3(v);
			//		}
			//	}

			//	if (pt == GenericIOParticleType::Stars) {
			//		if (av.num_components == 1) {
			//			RETURN_COMP_VALUE(av.stars[id - offset]);
			//		}

			//		if (av.num_components == 3) {
			//			RETURN_COMP_VECTOR3(v);
			//		}
			//	}
			//}

			RETURN_COMP_EMPTY;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			// 3 * 10
			types_and_blocks.resize(GenericIOParticleType::PTMax * (GenericIOBlockType::BTMax + gio_datas.size()/* + get_gio_values()*/));
			memset(types_and_blocks.data(), 0, types_and_blocks.size() * sizeof(int));
			//memcpy(types_and_blocks.data(), all_particles_blocks, sizeof(all_particles_blocks));

			int* tab = types_and_blocks.data();

			for (int ipt = 0; ipt < GenericIOParticleType::PTMax; ipt++) {
				// tipsy default
				for (int ibt = 0; ibt < GenericIOBlockType::BTMax; ibt++) {

					GenericIOParticleType pt = (GenericIOParticleType)ipt;
					GenericIOBlockType bt = (GenericIOBlockType)ibt;

					switch (pt) {
					case GenericIOParticleType::HACC: {

						switch (bt) {
						case GenericIOBlockType::Pos:
							if (gio_pos_idx[0] != -1 && gio_pos_idx[1] != -1 && gio_pos_idx[2] != -1)
								tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							break;
						case GenericIOBlockType::Vel:
							if (gio_vel_idx[0] != -1 && gio_vel_idx[1] != -1 && gio_vel_idx[2] != -1)
								tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							break;
							//case GenericIOBlockType::Mass:
							//	tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							//	break;
							//case GenericIOBlockType::Vel:
							//	tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							//	break;
							//case GenericIOBlockType::Phi:
							//	tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							//	break;
							//case GenericIOBlockType::HSmooth:
							//	tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							//	break;
							//case GenericIOBlockType::Rho:
							//	tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							//	break;
							//case GenericIOBlockType::Temp:
							//	tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							//	break;
							//case GenericIOBlockType::Metals:
							//	tab[ibt * GenericIOParticleType::PTMax + ipt]++;
							//	break;
						}
						break;
					}
					//case GenericIOParticleType::Dark: {

					//	switch (bt) {
					//	case GenericIOBlockType::Pos:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Mass:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Vel:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Soft:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Phi:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	}
					//	break;
					//}
					//case GenericIOParticleType::Stars: {

					//	switch (bt) {
					//	case GenericIOBlockType::Pos:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Mass:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Vel:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Soft:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Phi:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::Metals:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	case GenericIOBlockType::TForm:
					//		tab[ibt * GenericIOParticleType::PTMax + ipt]++;
					//		break;
					//	}
					//	break;
					//}
					}
				}

				//GIODataBase* av = gio_datas[blocknr - GenericIOBlockType::BTMax];

				// tipsy aux
				for (int ibt = GenericIOBlockType::BTMax; ibt < GenericIOBlockType::BTMax + gio_datas.size() /*get_gio_values()*/; ibt++) {
					GenericIOParticleType pt = (GenericIOParticleType)ipt;

					switch (pt) {
					case GenericIOParticleType::HACC: {
						tab[ibt * GenericIOParticleType::PTMax + ipt]++;
						break;
					}
					}
				}
			}
		}

		std::string get_dataset_name(int blocknr) {
			GenericIOBlockType bt = (GenericIOBlockType)blocknr;
			switch (bt) {
			case GenericIOBlockType::Pos:
				return "Pos[x,y,z]";
				//case GenericIOBlockType::Mass:
				//	return "Mass";
			case GenericIOBlockType::Vel:
				return "Vel[vx,vy,vz]";
				//case GenericIOBlockType::Soft:
				//	return "Soft";
				//case GenericIOBlockType::Phi:
				//	return "Phi";
				//case GenericIOBlockType::HSmooth:
				//	return "HSmooth";
				//case GenericIOBlockType::Rho:
				//	return "Rho";
				//case GenericIOBlockType::Temp:
				//	return "Temp";
				//case GenericIOBlockType::Metals:
				//	return "Metals";
				//case GenericIOBlockType::TForm:
				//	return "TForm";
			}

			if (blocknr >= GenericIOBlockType::BTMax && gio_datas.size() > blocknr - GenericIOBlockType::BTMax) {
				return gio_datas[blocknr - GenericIOBlockType::BTMax].get_name();
			}

			return "Unknown";
		}

		int get_particle_type(uint64_t id) {
			//if (id < gio_data->gas.size()) {
			//	return GenericIOParticleType::HACC;
			//}
			//else if (id < gio_data->gas.size() + gio_data->darks.size()) {
			//	return GenericIOParticleType::Dark;
			//}
			//else 
			//{
			return GenericIOParticleType::HACC;
			//}
		}

		void get_particle_position(uint64_t id, double* pos) {
			GenericIOParticleType pt = (GenericIOParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			switch (pt) {
			case GenericIOParticleType::HACC: {
				pos[0] = gio_datas[gio_pos_idx[0]].get_fvalue(id);
				pos[1] = gio_datas[gio_pos_idx[1]].get_fvalue(id);
				pos[2] = gio_datas[gio_pos_idx[2]].get_fvalue(id);
				break;
			}
			//case GenericIOParticleType::Dark: {
			//	pos[0] = gio_data->darks[id - offset].pos[0];
			//	pos[1] = gio_data->darks[id - offset].pos[1];
			//	pos[2] = gio_data->darks[id - offset].pos[2];
			//	break;
			//}
			//case GenericIOParticleType::Stars: {
			//	pos[0] = gio_data->stars[id - offset].pos[0];
			//	pos[1] = gio_data->stars[id - offset].pos[1];
			//	pos[2] = gio_data->stars[id - offset].pos[2];
			//	break;
			//}
			}
		}

		size_t get_local_num_particles() {
			if (gio_pos_idx[0] == -1)
				return 0;

			return gio_datas[gio_pos_idx[0]].get_count();// gio_data->h.nbodies; //gio_data->gas.size() + gio_data->darks.size() + gio_data->stars.size();
		}
		size_t get_global_num_particles() {
			//return gio_data->readTotalNumElems();  //gio_datas[0]->get_count();//gio_data->nbodies; //gio_data->gas.size() + gio_data->darks.size() + gio_data->stars.size();
			return g_total_particles;
		}

		//double get_particle_radius(uint64_t id) {
		//	return get_particle_hsml(id) * 2.0;
		//}

		double get_particle_hsml(uint64_t id) {
			// if (g_smoothlength_blocknr != -1)
			// 	return get_particle_norm_value(g_smoothlength_blocknr, id);

			if (gio_hsml_idx != -1)
				return gio_datas[gio_hsml_idx].get_fvalue(id);

			return 0.0;
		}

		double get_particle_mass(uint64_t id) {
			if (gio_mass_idx != -1)
				return gio_datas[gio_mass_idx].get_fvalue(id);

			return 0;// get_particle_norm_value(GenericIOBlockType::Mass, id);
		}

		int get_particle_rho_blocknr() {
			return GenericIOBlockType::BTMax + gio_rho_idx;
		}

		double get_particle_rho(uint64_t id) {
			if (gio_rho_idx != -1)
				return gio_datas[gio_rho_idx].get_fvalue(id);

			return 0;//get_particle_norm_value(GenericIOBlockType::Rho, id);
		}

		std::string get_particle_unit(int blocknr) {
			return "";
		}
	}
}