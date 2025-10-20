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

#include "data_processing.h"
#include "data_communication.h"

#include "gadget/gadget_convert_vdb.h"
#include "gadget/gadget_simple_convert_vdb.h"
#include "changa/changa_nchilada_convert_vdb.h"
#include "changa/changa_tipsy_convert_vdb.h"
#include "csv/csv_convert_vdb.h"
#include "common/convert_vdb.h"
#include "haccbin/haccbin_convert_vdb.h"

#ifdef WITH_GENERICIO
#	include "genericio/genericio_convert_vdb.h"
#endif


#ifdef WITH_HDF5
#	include "hdf5/hdf5_convert_vdb.h"
#endif

#ifdef WITH_OPENMP
#	include <omp.h>
#endif

#ifdef WITH_OPENVDB

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/NanoToOpenVDB.h>
#else
#	include <nanovdb/tools/NanoToOpenVDB.h>
#endif

#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointDataGrid.h>
#include <openvdb/tools/Count.h>

#include <fstream>
#include <filesystem>

#endif

#ifdef _WIN32
#	undef max              // disable the macro
#	undef min              // disable the macro
#	include <algorithm>    // for std::max
#endif


namespace space_converter {
#if 0
	class PrinterBase {
	public:
		virtual ~PrinterBase() {}
		virtual void print(std::ostream& os, size_t i) = 0;
	};

	template <class T>
	class Printer : public PrinterBase {
	public:
		Printer(gio::GenericIO& G, size_t MNE, size_t NE, const std::string& N)
			: NumElements(NE), Data(MNE* NE + G.requestedExtraSpace() / sizeof(T)) {
			G.addScalarizedVariable(N, Data, NE, gio::GenericIO::VarHasExtraSpace);
		}

		virtual void print(std::ostream& os, size_t i) {
			for (size_t j = 0; j < NumElements; ++j) {
				//os << scientific << setprecision(std::numeric_limits<T>::digits10) << Data[i * NumElements + j];
				os << Data[i * NumElements + j];

				if (j != NumElements - 1)
					os << "\t";
			}
		}

	protected:
		size_t NumElements;
		std::vector<T> Data;
	};

	template <typename T>
	PrinterBase* addPrinter(gio::GenericIO::VariableInfo& V,
		gio::GenericIO& GIO, size_t MNE) {
		if (sizeof(T) != V.ElementSize)
			return 0;

		if (V.IsFloat == std::numeric_limits<T>::is_integer)
			return 0;
		if (V.IsSigned != std::numeric_limits<T>::is_signed)
			return 0;

		return new Printer<T>(GIO, MNE, V.Size / V.ElementSize, V.Name);
	}
#endif

	void test_converter(int argc, char** argv, space_converter::FromCL& from_cl)
	{
#if 0
		//convert raw to vdb
		common::vdb::VDBParticles grid_main;
		common::SpaceData space_data_temp(0, 0, 512);

		space_data_temp.dense_type = common::SpaceData::DenseType::eType1;
		space_converter::create_grid(grid_main, from_cl, space_data_temp);

		FILE* file = fopen("e:\\temp\\cycles_anari\\magnetic_reconnection_512x512x512_float32.raw", "rb");
		size_t size = grid_main.dense_grid.memsize();
		fread((char*)grid_main.dense_grid.data_density.data(), size, 1, file);
		fclose(file);

		changa::ConvertVDBChangaTipsy* convert_vdb_base_temp = new changa::ConvertVDBChangaTipsy();
		common::vdb::VDBParticles grid_main_final;

		space_data_temp.transform_scale = 1.0;

		finalize_grid(convert_vdb_base_temp, from_cl, space_data_temp, grid_main, grid_main_final);
		save_vdb(convert_vdb_base_temp, from_cl, space_data_temp, grid_main_final);
#endif

#if 0
		std::string inFileName = "e:\\temp\\wdas\\wdas_cloud.vdb";
		//vec3i dims(0);
		std::string outFileName = "e:\\temp\\wdas\\wdas_cloud.raw";
		std::string outputFormat = "sparse";
		std::string gridName = "density";

		// Create a file object and open the VDB file
		openvdb::io::File file(inFileName);
		file.open();

		openvdb::GridBase::Ptr baseGrid;
		for (openvdb::io::File::NameIterator nameIter = file.beginName();
			nameIter != file.endName(); ++nameIter)
		{
			std::cout << "grid " << nameIter.gridName() << std::endl;
		}

		// Retrieve grid from the file
		//openvdb::GridPtrVec grids = file.getGrids();
		baseGrid = file.readGrid(gridName);

		// Close the file
		file.close();

		openvdb::FloatGrid::Ptr floatGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

		// Compute the active voxel bounding box
		openvdb::CoordBBox bbox = floatGrid->evalActiveVoxelBoundingBox();
		//if (floatGrid->evalActiveVoxelBoundingBox(bbox)) {
			// Get the dimensions of the bounding box
		openvdb::Coord dim = bbox.dim();
		std::cout << "dims: " << dim.x() << "," << dim.y() << "," << dim.z() << std::endl;

		// Create a 3D vector to store the grid values
		//std::vector<std::vector<std::vector<float>>> denseMatrix(dim.x(), std::vector<std::vector<float>>(dim.y(), std::vector<float>(dim.z(), 0.0f)));
		std::vector<float> denseMatrix((size_t)dim.x() * (size_t)dim.y() * (size_t)dim.z(), 0.0f);
		std::vector<int> denseMatrixCount((size_t)dim.x() * (size_t)dim.y() * (size_t)dim.z(), 0);

		float min = FLT_MAX;
		float max = -FLT_MAX;

		// Iterate over the active voxels and copy values to the dense matrix
		for (openvdb::FloatGrid::ValueOnCIter iter = floatGrid->cbeginValueOn(); iter; ++iter) {
			float v = *iter;
			if (v < min)
				min = v;

			if (v > max)
				max = v;
		}

		std::cout << "min: " << min << ", max:" << max << std::endl;

		// Iterate over the active voxels and copy values to the dense matrix
		for (openvdb::FloatGrid::ValueOnCIter iter = floatGrid->cbeginValueOn(); iter; ++iter) {
			openvdb::Coord xyz = iter.getCoord() - bbox.min();
			float v = *iter;// (((float)*iter - min) / (max - min));
			denseMatrix[xyz.x() + xyz.y() * (size_t)dim.x() + xyz.z() * (size_t)dim.x() * (size_t)dim.y()] = v;
			denseMatrixCount[xyz.x() + xyz.y() * (size_t)dim.x() + xyz.z() * (size_t)dim.x() * (size_t)dim.y()] += 1;
		}

		//for (size_t s = 0; s < denseMatrixCount.size(); s++) {
		//	if (denseMatrixCount[s] > 1) {
		//		printf("%lld: %d\n", denseMatrixCount[s]);
		//	}
		//}

		std::ofstream out(outFileName.c_str(), std::ios::binary);
		if (!out) {
			throw std::runtime_error("Failed to open file for writing");
		}
		out.write((const char*)denseMatrix.data(), denseMatrix.size() * sizeof(denseMatrix[0]));
#endif

#if 0
		bool ShowMap = false;
		bool NoData = false;
		bool PrintRankInfo = true;
		std::string FileName = "e:\\data\\jar091\\down\\m000.mpicosmo.499";

		int Rank, NRanks;
		MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
		MPI_Comm_size(MPI_COMM_WORLD, &NRanks);

		if (Rank == 0) {
			unsigned Method = gio::GenericIO::FileIOPOSIX;
			const char* EnvStr = getenv("GENERICIO_USE_MPIIO");
			if (EnvStr && std::string(EnvStr) == "1")
				Method = gio::GenericIO::FileIOMPI;

			gio::GenericIO GIO(MPI_COMM_SELF, FileName, Method);

			GIO.openAndReadHeader(gio::GenericIO::MismatchAllowed, -1, !ShowMap);

			int NR = GIO.readNRanks();

			std::vector<gio::GenericIO::VariableInfo> VI;
			GIO.getVariableInfo(VI);

			size_t MaxNElem = 0;
			for (int i = 0; i < NR; ++i) {
				size_t NElem = GIO.readNumElems(i);
				MaxNElem = std::max(MaxNElem, NElem);
			}

			std::vector<PrinterBase*> Printers;
			for (size_t i = 0; i < VI.size(); ++i) {
				PrinterBase* P = 0;

#define ADD_PRINTER(T) \
      if (!P) P = addPrinter<T>(VI[i], GIO, MaxNElem)
				ADD_PRINTER(float);
				ADD_PRINTER(double);
				ADD_PRINTER(unsigned char);
				ADD_PRINTER(signed char);
				ADD_PRINTER(int16_t);
				ADD_PRINTER(uint16_t);
				ADD_PRINTER(int32_t);
				ADD_PRINTER(uint32_t);
				ADD_PRINTER(int64_t);
				ADD_PRINTER(uint64_t);
#undef ADD_PRINTER 

				if (!P) throw std::runtime_error("Don't know how to print variable: " + VI[i].Name);
				Printers.push_back(P);
			}

			int Dims[3];
			GIO.readDims(Dims);

			std::cout << "# " << FileName << ": " << NR << " rank(s): " <<
				Dims[0] << "x" << Dims[1] << "x" << Dims[2];

			uint64_t TotalNumElems = GIO.readTotalNumElems();
			if (TotalNumElems != (uint64_t)-1)
				std::cout << ": " << GIO.readTotalNumElems() << " row(s)";
			std::cout << std::endl;

			double PhysOrigin[3], PhysScale[3];
			GIO.readPhysOrigin(PhysOrigin);
			GIO.readPhysScale(PhysScale);
			if (PhysScale[0] != 0.0 || PhysScale[1] != 0.0 || PhysScale[2] != 0.0) {
				std::cout << "# physical coordinates: (" << PhysOrigin[0] << "," <<
					PhysOrigin[1] << "," << PhysOrigin[2] << ") -> (" <<
					PhysScale[0] << "," << PhysScale[1] << "," <<
					PhysScale[2] << ")" << std::endl;


				std::vector<gio::GenericIO::VariableInfo> VIX, VIY, VIZ;
				for (size_t i = 0; i < VI.size(); ++i) {
					if (VI[i].IsPhysCoordX) VIX.push_back(VI[i]);
					if (VI[i].IsPhysCoordY) VIY.push_back(VI[i]);
					if (VI[i].IsPhysCoordZ) VIZ.push_back(VI[i]);
				}

				if (!VIX.empty()) {
					std::cout << "# x variables: ";
					for (size_t i = 0; i < VIX.size(); ++i) {
						if (i != 0) std::cout << ", ";
						std::cout << VIX[i].Name;
						if (VIX[i].MaybePhysGhost)
							std::cout << " [maybe ghost]";
					}
					std::cout << std::endl;
				}
				if (!VIY.empty()) {
					std::cout << "# y variables: ";
					for (size_t i = 0; i < VIY.size(); ++i) {
						if (i != 0) std::cout << ", ";
						std::cout << VIY[i].Name;
						if (VIY[i].MaybePhysGhost)
							std::cout << " [maybe ghost]";
					}
					std::cout << std::endl;
				}
				if (!VIZ.empty()) {
					std::cout << "# z variables: ";
					for (size_t i = 0; i < VIZ.size(); ++i) {
						if (i != 0) std::cout << ", ";
						std::cout << VIZ[i].Name;
						if (VIZ[i].MaybePhysGhost)
							std::cout << " [maybe ghost]";
					}
					std::cout << std::endl;
				}
			}

			std::cout << "# ";
			for (size_t i = 0; i < VI.size(); ++i) {
				if (VI[i].Size == VI[i].ElementSize) {
					std::cout << VI[i].Name;
				}
				else {
					size_t NumElements = VI[i].Size / VI[i].ElementSize;
					for (size_t j = 0; j < NumElements; ++j) {
						std::cout << VI[i].Name;
						if (j == 0) {
							std::cout << ".x";
						}
						else if (j == 1) {
							std::cout << ".y";
						}
						else if (j == 2) {
							std::cout << ".z";
						}
						else if (j == 3) {
							std::cout << ".w";
						}
						else {
							std::cout << ".w" << (j - 3);
						}

						if (j != NumElements - 1)
							std::cout << "\t";
					}
				}

				if (i != VI.size() - 1)
					std::cout << "\t";
			}
			std::cout << std::endl;

			std::cout << "# ";
			for (size_t i = 0; i < VI.size(); ++i) {
				if (VI[i].Size == VI[i].ElementSize) {
					std::cout << "(" << (VI[i].IsFloat ? "f" :
						(VI[i].IsSigned ? "s" : "u")) << 8 * VI[i].Size << ")";
				}
				else {
					size_t NumElements = VI[i].Size / VI[i].ElementSize;
					for (size_t j = 0; j < NumElements; ++j) {
						std::cout << "(" << (VI[i].IsFloat ? "f" :
							(VI[i].IsSigned ? "s" : "u")) <<
							8 * VI[i].ElementSize << ")";

						if (j != NumElements - 1)
							std::cout << "\t";
					}
				}

				if (i != VI.size() - 1)
					std::cout << "\t";
			}
			std::cout << std::endl;

			for (int i = 0; i < NR; ++i) {
				size_t NElem = GIO.readNumElems(i);
				int Coords[3];
				GIO.readCoords(Coords, i);
				if (PrintRankInfo)
					std::cout << "# rank " << GIO.readGlobalRankNumber(i) << ": " <<
					Coords[0] << "," << Coords[1] << "," <<
					Coords[2] << ": " << NElem << " row(s)" << std::endl;
				if (NoData)
					continue;

				GIO.readData(i, false);
				for (size_t j = 0; j < NElem; ++j) {
					for (size_t k = 0; k < Printers.size(); ++k) {
						Printers[k]->print(std::cout, j);
						if (k != Printers.size() - 1)
							std::cout << "\t";
					}
					std::cout << std::endl;
				}
			}

			for (size_t i = 0; i < Printers.size(); ++i) {
				delete Printers[i];
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();

#endif
		exit(0);
	}

	common::vdb::ConvertVDBBase* init_converter(int argc, char** argv, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		common::vdb::ConvertVDBBase* convert_vdb_base = nullptr;

		if (from_cl.world_rank == 0)
			printf("Data: %s\n", from_cl.data_type.c_str());

		if (from_cl.data_type == "GADGET") {
			convert_vdb_base = new gadget::ConvertVDBGadget();
		}
		else if (from_cl.data_type == "GADGET_SIMPLE") {
			convert_vdb_base = new gadget_simple::ConvertVDBGadgetSimple();
		}
		else if (from_cl.data_type == "CHANGA_TIPSY") {
			convert_vdb_base = new changa::ConvertVDBChangaTipsy();
		}
		else if (from_cl.data_type == "CHANGA_NCHILADA") {
			convert_vdb_base = new changa::ConvertVDBChangaNchilada();
		}
		else if (from_cl.data_type == "CSV") {
			convert_vdb_base = new csv::ConvertVDBCSV();
		}
#ifdef WITH_GENERICIO
		else if (from_cl.data_type == "GENERICIO") {
			convert_vdb_base = new genericio::ConvertVDBGenericIO();
		}
#endif
#ifdef WITH_HDF5
		else if (from_cl.data_type == "HDF5") {
			convert_vdb_base = new hdf5::ConvertVDBHDF5();
		}
#endif
		else if (from_cl.data_type == "HACCBIN") {
			convert_vdb_base = new haccbin::ConvertVDBHACCBin();
		}
		else {
			throw std::runtime_error("Unknown data type [GADGET, CHANGA_TIPSY, CHANGA_NCHILADA, CSV, GENERICIO, HDF5, HACCBIN]): " + from_cl.data_type);
		}

//#if 1 //test

		const char* converter_split_init_count = getenv("CONVERTER_SPLIT_INIT_COUNT");
		if (converter_split_init_count) {
			printf("WARNING: This is a test version of the converter, it will not work with real data!\n");

			// Ensure only one process enters the function at a time
			// for (int i = 0; i < from_cl.world_size; ++i) {			
			// 	if (from_cl.world_rank == i) {
			// 		printf("rank: %d: init_lib\n", from_cl.world_rank);
			// 		convert_vdb_base->init_lib(argc, argv, from_cl.world_rank, from_cl.world_size);
			// 	}
			// 	// Synchronize all processes before moving to the next
			// 	MPI_Barrier(MPI_COMM_WORLD);
			// }

			MPI_Barrier(MPI_COMM_WORLD);
			double t_init = omp_get_wtime();

			int group_size = atoi(converter_split_init_count);
			int num_groups = (from_cl.world_size + group_size - 1) / group_size;

			for (int g = 0; g < num_groups; ++g) {
				int start = g * group_size;
				int end = std::min(start + group_size, from_cl.world_size);

				if (from_cl.world_rank >= start && from_cl.world_rank < end) {
					//printf("rank: %d: init_lib\n", from_cl.world_rank);
					convert_vdb_base->init_lib(argc, argv, from_cl.world_rank, from_cl.world_size);
				}

				// Synchronize all processes before next group proceeds
				MPI_Barrier(MPI_COMM_WORLD);
			}

			MPI_Barrier(MPI_COMM_WORLD);
			double t_init_end = omp_get_wtime();
			if (from_cl.world_rank == 0) {
				printf("rank: %d: CONVERTER_SPLIT_INIT time: %f\n", from_cl.world_rank, t_init_end - t_init);
			}
					
			//printf("rank: %d: init_lib done\n", from_cl.world_rank);
		}
	else {
		convert_vdb_base->init_lib(argc, argv, from_cl.world_rank, from_cl.world_size);
	}

#ifdef WITH_EMBREE
		double t_embree = omp_get_wtime();
		convert_vdb_base->create_embree_scene(/*particle_type*/ 0); //TODO
		printf("rank: %d: create_embree_scene: %f\n", world_rank, omp_get_wtime() - t_embree);
#endif

// #if 0
// 		convert_vdb_base->calculate_radius_by_cudakdtree(from_cl.calc_radius_neigh, from_cl.calc_radius_neigh_file);
// 		std::vector<float> radius_particles_cuda(convert_vdb_base->radius_particles.size());
// 		memcpy(radius_particles_cuda.data(), convert_vdb_base->radius_particles.data(), radius_particles_cuda.size() * sizeof(float));

// 		from_cl.calc_radius_neigh_file = "e:\\temp\\gadget_rad10.bin";
// 		from_cl.calc_radius_neigh = -1;
// 		convert_vdb_base->calculate_radius_by_nanoflann(from_cl.calc_radius_neigh, from_cl.calc_radius_neigh_file);
// 		std::vector<float> radius_particles_nano(convert_vdb_base->radius_particles.size());
// 		memcpy(radius_particles_nano.data(), convert_vdb_base->radius_particles.data(), radius_particles_nano.size() * sizeof(float));

// 		// Calculate MSE between nanoflann and cudakdtree radius calculations
// 		double mse = 0.0;
// 		if (radius_particles_nano.size() == radius_particles_cuda.size() && radius_particles_nano.size() > 0) {
// 			double sum_squared_diff = 0.0;

// #pragma omp parallel for reduction(+:sum_squared_diff)
// 			for (size_t i = 0; i < radius_particles_nano.size(); i++) {
// 				double diff = static_cast<double>(radius_particles_nano[i]) - static_cast<double>(radius_particles_cuda[i]);
// 				sum_squared_diff += diff * diff;
// 			}
// 			mse = sum_squared_diff / radius_particles_nano.size();
			
// 			if (from_cl.world_rank == 0) {
// 				printf("MSE between nanoflann and cudakdtree radius calculations: %e\n", mse);
// 			}
// 		} else if (from_cl.world_rank == 0) {
// 			printf("Cannot calculate MSE: arrays have different sizes or are empty\n");
// 		}

// 		fflush(0);
// 		MPI_Finalize();
// 		exit(0);

// #endif

		// Calculate radius for particles
		bool use_cycling = space_data.anim_type == common::SpaceData::AnimType::eNone;// || from_cl.use_anim_merge;		

		std::string calc_radius_neigh_file = space_data.calc_radius_neigh_file;
		// if (calc_radius_neigh_file.length() > 0) {
		// 	calc_radius_neigh_file = calc_radius_neigh_file + "." + std::to_string(from_cl.world_rank);
		// }

#ifdef WITH_CUDAKDTREE
		if ((from_cl.use_cudakdtree || from_cl.use_cudakdtree_cpu) && (space_data.calc_radius_neigh > 0 || calc_radius_neigh_file.length() > 0)) {
			float max_radius = std::numeric_limits<float>::infinity();
			const char* converter_knn_max_radius = getenv("CONVERTER_KNN_MAX_RADIUS");
			if (converter_knn_max_radius) {
				max_radius = static_cast<float>(atof(converter_knn_max_radius));
				printf("Using max_radius: %f\n", max_radius);
			}
			common::SpaceData::DenseType rho_kernel = space_data.calc_radius_neigh_rho_kernel;
			convert_vdb_base->calculate_radius_by_cudakdtree(space_data.calc_radius_neigh, calc_radius_neigh_file, use_cycling, from_cl.use_cudakdtree_cpu, max_radius, rho_kernel);
		}
		else
#endif

#ifdef WITH_NANOFLANN		
		if (from_cl.use_nanoflann && (space_data.calc_radius_neigh > 0 || calc_radius_neigh_file.length() > 0)) {
			common::SpaceData::DenseType rho_kernel = space_data.calc_radius_neigh_rho_kernel;
			convert_vdb_base->calculate_radius_by_nanoflann(space_data.calc_radius_neigh, calc_radius_neigh_file, use_cycling, rho_kernel);
		}		
		else
#endif
		if (calc_radius_neigh_file.length() > 0) {
			//calc_radius_neigh_file = calc_radius_neigh_file + "." + std::to_string(from_cl.world_rank + from_cl.anim_start);
			convert_vdb_base->read_radius_from_file(calc_radius_neigh_file);
		}

		return convert_vdb_base;
	}

	void deinit_converter(common::vdb::ConvertVDBBase* convert_vdb_base)
	{
		convert_vdb_base->finish_lib();
		delete convert_vdb_base;
	}

#ifdef WITH_OPENVDB
	void eval_min_max(openvdb::FloatGrid::Ptr grid, float &min_val, float &max_val)
	{
		openvdb::math::MinMax<float> extrema = openvdb::tools::minMax(grid->tree());
		min_val = extrema.min();
		max_val = extrema.max();
	}
#endif

	void print_info(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, std::vector<int>& types_and_blocks_global)
	{
		std::vector<int> types_and_blocks_local;
		convert_vdb_base->get_types_and_blocks(types_and_blocks_local);

		types_and_blocks_global.resize(types_and_blocks_local.size());

		// Use MPI_Allreduce to find the sum across all arrays for each element
		MPI_Allreduce(types_and_blocks_local.data(), types_and_blocks_global.data(), types_and_blocks_local.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		if (from_cl.world_rank == 0) {
			convert_vdb_base->print_types_and_blocks(types_and_blocks_global);
		}

		if (from_cl.info) {
			convert_vdb_base->finish_lib();

			MPI_Finalize();
			exit(0);
		}

	}

	void find_bbox(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, int particle_type)
	{
		double t_message_type2 = omp_get_wtime();
		double t_find_bbox = omp_get_wtime();

		float bbox_min_orig_local[3], bbox_max_orig_local[3], bbox_min_orig[3], bbox_max_orig[3];

		convert_vdb_base->iolib_find_bbox(
			particle_type,
			bbox_min_orig_local,
			bbox_max_orig_local,
			space_data.offset_position
		);

		if (from_cl.world_rank == 0)
			printf("rank: %d: find bbox local: %f\n", from_cl.world_rank, omp_get_wtime() - t_message_type2);

		t_find_bbox = omp_get_wtime();

		// Use MPI_Allreduce to find the minimum across all arrays for each element
		MPI_Allreduce(bbox_min_orig_local, bbox_min_orig, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(bbox_max_orig_local, bbox_max_orig, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

		space_data.bbox_min_orig[0] = static_cast<int>(bbox_min_orig[0] - 1.0f);
		space_data.bbox_min_orig[1] = static_cast<int>(bbox_min_orig[1] - 1.0f);
		space_data.bbox_min_orig[2] = static_cast<int>(bbox_min_orig[2] - 1.0f);

		space_data.bbox_max_orig[0] = static_cast<int>(bbox_max_orig[0] + 1.0f);
		space_data.bbox_max_orig[1] = static_cast<int>(bbox_max_orig[1] + 1.0f);
		space_data.bbox_max_orig[2] = static_cast<int>(bbox_max_orig[2] + 1.0f);

		space_data.bbox_size_orig = std::max((double)space_data.bbox_max_orig[0] - (double)space_data.bbox_min_orig[0], (double)space_data.bbox_max_orig[1] - (double)space_data.bbox_min_orig[1]);
		space_data.bbox_size_orig = std::max(space_data.bbox_size_orig, (double)space_data.bbox_max_orig[2] - (double)space_data.bbox_min_orig[2]);

#if 1 // symmetric by x,y,z
		space_data.bbox_min_orig[0] = (space_data.bbox_min_orig[0] + space_data.bbox_max_orig[0]) / 2.0 - space_data.bbox_size_orig / 2.0;
		space_data.bbox_min_orig[1] = (space_data.bbox_min_orig[1] + space_data.bbox_max_orig[1]) / 2.0 - space_data.bbox_size_orig / 2.0;
		space_data.bbox_min_orig[2] = (space_data.bbox_min_orig[2] + space_data.bbox_max_orig[2]) / 2.0 - space_data.bbox_size_orig / 2.0;

		space_data.bbox_max_orig[0] = (space_data.bbox_min_orig[0] + space_data.bbox_max_orig[0]) / 2.0 + space_data.bbox_size_orig / 2.0;
		space_data.bbox_max_orig[1] = (space_data.bbox_min_orig[1] + space_data.bbox_max_orig[1]) / 2.0 + space_data.bbox_size_orig / 2.0;
		space_data.bbox_max_orig[2] = (space_data.bbox_min_orig[2] + space_data.bbox_max_orig[2]) / 2.0 + space_data.bbox_size_orig / 2.0;
#endif

		if (from_cl.world_rank == 0)
			printf("rank: %d: find bbox mpi: %f, box_size: %f\n", from_cl.world_rank, omp_get_wtime() - t_message_type2, space_data.bbox_size_orig);
	}

	void create_grid(common::vdb::VDBParticles& grid_main, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		//Dense
		if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {
			grid_main.type = common::vdb::VDBParticles::VDBParticleType::eDense;
			grid_main.dense_grid.create(space_data.bbox_dim, space_data.bbox_dim, space_data.bbox_dim);
		}
		//Raw
		else if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
			grid_main.type = common::vdb::VDBParticles::VDBParticleType::eRawParticles;
		}
		else if (from_cl.use_nanovdb) {
			grid_main.type = common::vdb::VDBParticles::VDBParticleType::eNanoVDB;
#if OPENVDB_VERSION == 11
			grid_main.nano_grid = std::make_shared<nanovdb::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#else
			grid_main.nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#endif
		}
		else {
#ifdef WITH_OPENVDB
			grid_main.type = common::vdb::VDBParticles::VDBParticleType::eOpenVDB;
			grid_main.vdb_grid = openvdb::FloatGrid::create(0.0f);
			grid_main.vdb_grid->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);
			grid_main.vdb_grid->setName("density");
#endif
		}
	}

	void convert_to_grid(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main)
	{
		double t_convert = omp_get_wtime();

		convert_vdb_base->convert_iolib_to_grid(
			space_data.particle_type,
			space_data.particle_fix_size,
			"density",
			//grid_dims,
			space_data.grid_transform,
			space_data.bbox_min,
			space_data.bbox_max,
			space_data.bbox_dim,
			space_data.bbox_min_orig,
			space_data.bbox_size_orig,
			space_data.extracted_type,
			space_data.dense_type,
			space_data.dense_norm,
			space_data.block_name_id,
			space_data.object_size,
			space_data.min_value_local,
			space_data.max_value_local,
			space_data.min_value,
			space_data.max_value,
			space_data.particles_count_local,
			grid_main,
			space_data.transform_scale,
			space_data.filter_min,
			space_data.filter_max,
			space_data.min_rho,
			space_data.max_rho,
			space_data.anim_type,
			space_data.frame,
			space_data.anim_start + from_cl.world_rank,
			space_data.bbox_sphere_pos,
			space_data.bbox_sphere_r,
			space_data.use_simple_density,
			space_data.offset_position
		);

		if (from_cl.world_rank == 0) {
			printf("rank: %d: grid convert: %f\n", from_cl.world_rank, omp_get_wtime() - t_convert);

			float bbox_x_mid_orig = ((((space_data.bbox_min[0] / (float)space_data.object_size) + 
				(space_data.bbox_max[0] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[0] - (double)space_data.bbox_min_orig[0])) + (double)space_data.bbox_min_orig[0];

			float bbox_y_mid_orig = ((((space_data.bbox_min[1] / (float)space_data.object_size) +
				(space_data.bbox_max[1] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[1] - (double)space_data.bbox_min_orig[1])) + (double)space_data.bbox_min_orig[1];

			float bbox_z_mid_orig = ((((space_data.bbox_min[2] / (float)space_data.object_size) +
				(space_data.bbox_max[2] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[2] - (double)space_data.bbox_min_orig[2])) + (double)space_data.bbox_min_orig[2];

			float bbox_x_r_orig = ((((-space_data.bbox_min[0] / (float)space_data.object_size) +
				(space_data.bbox_max[0] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[0] - (double)space_data.bbox_min_orig[0]));// +(double)space_data.bbox_min_orig[0];

			float bbox_y_r_orig = ((((-space_data.bbox_min[1] / (float)space_data.object_size) +
				(space_data.bbox_max[1] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[1] - (double)space_data.bbox_min_orig[1]));// + (double)space_data.bbox_min_orig[1];

			float bbox_z_r_orig = ((((-space_data.bbox_min[2] / (float)space_data.object_size) +
				(space_data.bbox_max[2] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[2] - (double)space_data.bbox_min_orig[2]));// + (double)space_data.bbox_min_orig[2];

			printf("rank: %d: bbox-sphere coord: %f, %f, %f, bbox-sphere radius: %f\n", 
				from_cl.world_rank, bbox_x_mid_orig, bbox_y_mid_orig, bbox_z_mid_orig,
				std::max(bbox_x_r_orig, std::max(bbox_y_r_orig, bbox_z_r_orig)));

			printf("rank: %d: bbox coord: %f, %f, %f, %f, %f, %f\n",
				from_cl.world_rank, space_data.bbox_min[0], space_data.bbox_min[1], space_data.bbox_min[2],
				space_data.bbox_max[0], space_data.bbox_max[1], space_data.bbox_max[2]);
		}
	}

	void find_minmax_value(space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		double t_find_minmax = omp_get_wtime();

		// Use MPI_Allreduce to find the minimum across all arrays for each element
		MPI_Allreduce(&space_data.min_value_local, &space_data.min_value, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&space_data.max_value_local, &space_data.max_value, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

		MPI_Allreduce(&space_data.particles_count_local, &space_data.particles_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		if (from_cl.world_rank == 0)
			printf("rank: %d: find minmax mpi: %f, particles_count: %lld\n", from_cl.world_rank, omp_get_wtime() - t_find_minmax, space_data.particles_count);
	}

	void find_minmax_reduced_value(space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
			double t_find_minmax = omp_get_wtime();

			float min_value_reduced = space_data.min_value_reduced;
			float max_value_reduced = space_data.max_value_reduced;

			// Use MPI_Allreduce to find the minimum across all arrays for each element
			MPI_Allreduce(&min_value_reduced, &space_data.min_value_reduced, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
			MPI_Allreduce(&max_value_reduced, &space_data.max_value_reduced, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

			if (from_cl.world_rank == 0)
				printf("rank: %d: find reduced minmax mpi: %f\n", from_cl.world_rank, omp_get_wtime() - t_find_minmax);
		}
	}

	void reduction(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, common::vdb::VDBParticles& grid_main_sum)
	{
		double t_grid = omp_get_wtime();

#if 0 //TODO
		if (space_data.anim_type == common::SpaceData::AnimType::eFrameCache) {
			int frame_req = space_data.frame;
			int frame = space_data.anim_start + from_cl.world_rank;
			int rank_req = frame_req - space_data.anim_start;

			if (rank_req < 0 || rank_req >= from_cl.world_size) {
				return; //out of range
			}

			//Dense
			if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {
				grid_main_sum.type = common::vdb::VDBParticles::VDBParticleType::eDense;

				if (from_cl.world_rank == 0) {
					grid_main_sum.dense_grid.create(space_data.bbox_dim, space_data.bbox_dim, space_data.bbox_dim);
					memcpy(grid_main_sum.dense_grid.offset, grid_main.dense_grid.offset, sizeof(grid_main.dense_grid.offset));
					//grid_main_sum.dense_grid.type = grid_main.dense_grid.type;

					if (frame_req == frame) {
						memcpy(grid_main_sum.dense_grid.data_density.data(), grid_main.dense_grid.data_density.data(), grid_main.dense_grid.memsize());
#ifndef WITH_NO_DATA_TEMP						
						memcpy(grid_main_sum.dense_grid.data_temp.data(), grid_main.dense_grid.data_temp.data(), grid_main.dense_grid.memsize());
#endif						
					}
					else {
						mpi_recv(grid_main_sum.dense_grid.data_density.data(), grid_main_sum.dense_grid.data_density.size(), MPI_BYTE, from_cl.world_rank, 0);
#ifndef WITH_NO_DATA_TEMP						
						mpi_recv(grid_main_sum.dense_grid.data_temp.data(), grid_main_sum.dense_grid.data_temp.size(), MPI_BYTE, from_cl.world_rank, 0);
#endif						
					}
				}
				else if (frame_req == frame) {
					mpi_recv(grid_main.dense_grid.data_density.data(), grid_main.dense_grid.data_density.size(), MPI_BYTE, from_cl.world_rank, 0);
#ifndef WITH_NO_DATA_TEMP					
					mpi_recv(grid_main.dense_grid.data_temp.data(), grid_main.dense_grid.data_temp.size(), MPI_BYTE, from_cl.world_rank, 0);
#endif					
				}
			}
			else {
				grid_main_sum.type = common::vdb::VDBParticles::VDBParticleType::eVector;

				if (from_cl.world_rank == 0) {
					size_t ns = 0;
					mpi_recv(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank, 0);
					grid_main_sum.vector_grid.resize(ns);

					mpi_recv(grid_main_sum.vector_grid.data(), ns, MPI_BYTE, from_cl.world_rank, 0);
				}
				else if (frame_req == frame) {
					if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
						std::vector<uint8_t> grid_handle_main;
						grid_main_sum.raw_particles.serialize(grid_handle_main);
						size_t ns = grid_handle_main.size();
						mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank, 0);
						mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank, 0);
					}
					else if (from_cl.use_nanovdb) {
#if OPENVDB_VERSION == 11
						nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_main = nanovdb::createNanoGrid(*grid_main_sum.nano_grid);
#else
						nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_main = nanovdb::tools::createNanoGrid(*grid_main_sum.nano_grid);
#endif
						size_t ns = grid_handle_main.size();
						mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank, 0);
						mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank, 0);
					}
					else {
						std::vector<uint8_t> grid_handle_main;
#ifdef WITH_OPENVDB
						convert_vdb_base->openvdb_to_vector(grid_main_sum.vdb_grid, grid_handle_main);
#endif
						size_t ns = grid_handle_main.size();
						mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank, 0);
						mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank, 0);
					}
				}
			}
		}
		else
#endif
		// Dense		
		if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {
			//common::vdb::DenseParticles grid_main_sum;
			grid_main_sum.type = common::vdb::VDBParticles::VDBParticleType::eDense;

			if (from_cl.world_rank == 0 || space_data.anim_type != common::SpaceData::AnimType::eNone) {
				grid_main_sum.dense_grid.create(space_data.bbox_dim, space_data.bbox_dim, space_data.bbox_dim);
				memcpy(grid_main_sum.dense_grid.offset, grid_main.dense_grid.offset, sizeof(grid_main.dense_grid.offset));
				//grid_main_sum.dense_grid.type = grid_main.dense_grid.type;
			}

			if (space_data.anim_type != common::SpaceData::AnimType::eNone) {
				// Copy
				memcpy(grid_main_sum.dense_grid.data_density.data(), grid_main.dense_grid.data_density.data(), grid_main.dense_grid.memsize());
#ifndef WITH_NO_DATA_TEMP				
				memcpy(grid_main_sum.dense_grid.data_temp.data(), grid_main.dense_grid.data_temp.data(), grid_main.dense_grid.memsize());
#endif				
			}
			else {
				// Reduce all the local sums to a global sum in the root process (rank 0)			
				mpi_reduce(grid_main.dense_grid.data_density.data(), grid_main_sum.dense_grid.data_density.data(), grid_main.dense_grid.size());
#ifndef WITH_NO_DATA_TEMP				
				mpi_reduce(grid_main.dense_grid.data_temp.data(), grid_main_sum.dense_grid.data_temp.data(), grid_main.dense_grid.size());
#endif				
			}

			grid_main.dense_grid.clear();
		}
		else {
			grid_main_sum.type = grid_main.type;

			if (grid_main_sum.type == common::vdb::VDBParticles::VDBParticleType::eOpenVDB) {
#ifdef WITH_OPENVDB
				grid_main_sum.vdb_grid = grid_main.vdb_grid;
#endif
			}
			else if (grid_main_sum.type == common::vdb::VDBParticles::VDBParticleType::eNanoVDB)
				grid_main_sum.nano_grid = grid_main.nano_grid;
			else if (grid_main_sum.type == common::vdb::VDBParticles::VDBParticleType::eRawParticles)
				grid_main_sum.raw_particles = grid_main.raw_particles;

			if (space_data.anim_type == common::SpaceData::AnimType::eNone || space_data.anim_type == common::SpaceData::AnimType::eAllMerge || space_data.anim_type == common::SpaceData::AnimType::eFrameExtract) {
				// Logarithmic reduction steps
				for (int step = 1; step < from_cl.world_size; step *= 2) {
					if (from_cl.world_rank % (2 * step) == 0) {
						if (from_cl.world_rank + step < from_cl.world_size) {

							size_t ns = 0;
							mpi_recv(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank + step, 0);
							//std::vector<uint8_t> nanogrid_recv(ns);
							common::vdb::VDBParticles nanogrid_recv;
							nanogrid_recv.type = common::vdb::VDBParticles::VDBParticleType::eVector;
							nanogrid_recv.vector_grid.resize(ns);

							mpi_recv(nanogrid_recv.vector_grid.data(), ns, MPI_BYTE, from_cl.world_rank + step, 0);
							convert_vdb_base->merge_grid(grid_main_sum, nanogrid_recv);
						}
					}
					else if (from_cl.world_rank % (2 * step) == step) {
						if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
							std::vector<uint8_t> grid_handle_main;
							grid_main_sum.raw_particles.serialize(grid_handle_main);
							size_t ns = grid_handle_main.size();
							mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank - step, 0);
							mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank - step, 0);
						}
						else if (from_cl.use_nanovdb) {
#if OPENVDB_VERSION == 11
							nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_main = nanovdb::createNanoGrid(*grid_main_sum.nano_grid);
#else
							nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_main = nanovdb::tools::createNanoGrid(*grid_main_sum.nano_grid);
#endif
							size_t ns = grid_handle_main.size();
							mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank - step, 0);
							mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank - step, 0);
						}
						else {
							std::vector<uint8_t> grid_handle_main;
#ifdef WITH_OPENVDB
							convert_vdb_base->openvdb_to_vector(grid_main_sum.vdb_grid, grid_handle_main);
#endif
							size_t ns = grid_handle_main.size();
							mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank - step, 0);
							mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank - step, 0);
						}

						break; // This rank is done participating
					}
				}
			}
		}

		if (from_cl.world_rank == 0) {
			printf("rank: %d: merged time: %f\n", from_cl.world_rank, omp_get_wtime() - t_grid);
		}
	}
#ifdef WITH_MULTIRES
	class Openvdb_Multi_Res_Grids{    
	public:
		std::vector<openvdb::FloatGrid::Ptr> grids_array;
		std::vector<openvdb::FloatGrid::Accessor> grids_accessors;
		std::vector<float> data_range;
		std::vector<std::vector<float>> data_bbox;
		std::vector<openvdb::math::Transform::Ptr> point_to_inx_transforms;

#if 0
		// error variables
        double sum_diff_abs;
        double sum_diff_squared;
		double sum_diff_squared_norm;
		float global_max_value;
		int num_voxels_for_metrics; // non zero values
#endif

		Openvdb_Multi_Res_Grids(std::vector<std::vector<float>> data_bb, int max_division){
#if 0
			sum_diff_abs = 0.0;
			sum_diff_squared = 0.0;
			sum_diff_squared_norm = 0.0;
			num_voxels_for_metrics = 0;
#endif

			for (size_t i = 0; i < data_bb[0].size(); ++i) {
				data_range.push_back(data_bb[1][i] - data_bb[0][i]);
			}
			data_bbox = data_bb;

			// create empty grid for each resolution
			int num_vox = 1;
			for (int i = 0; i < (max_division + 1); i++){
				append_grid(num_vox, i);
				num_vox *= 2;
			}
				
		}

		void append_grid(int grid_num_voxels, int grid_inx) {
			openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);
			grids_array.push_back(grid);
			
			double x = data_range[0] / grid_num_voxels;
			double y = data_range[1] / grid_num_voxels;
			double z = data_range[2] / grid_num_voxels;

			// for houdini - move by half a voxel size
			// openvdb::math::Mat4d transform_mat(
			// 	x, 0.0, 0.0, 0.0,
			// 	0.0, y, 0.0, 0.0,
			// 	0.0, 0.0, z, 0.0,
			// 	data_bbox[0][0] + 0.5 * x, (data_bbox[0][1] + 0.5 * y), (data_bbox[0][2] + 0.5 * z), 1.0
			// );

			// for blender - do not move by half a voxel size
			openvdb::math::Mat4d transform_mat(
				x, 0.0, 0.0, 0.0,
				0.0, y, 0.0, 0.0,
				0.0, 0.0, z, 0.0,
				double(data_bbox[0][0]), double(data_bbox[0][1]), double(data_bbox[0][2]), 1.0
			);

			openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_mat);
			grids_array.back()->setTransform(transform);
			grids_array.back()->setName("density" + std::to_string(grid_inx));
			grids_array.back()->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);

			// create accessor
			grids_accessors.push_back(grids_array.back()->getAccessor());

			// create transform mat for conversion of point to index
			openvdb::math::Mat4d point_to_inx_transform_mat(
				x, 0.0, 0.0, 0.0,
				0.0, y, 0.0, 0.0,
				0.0, 0.0, z, 0.0,
				double(data_bbox[0][0]), double(data_bbox[0][1]), double(data_bbox[0][2]), 1.0
			);
			openvdb::math::Transform::Ptr transform_point_to_inx = openvdb::math::Transform::createLinearTransform(point_to_inx_transform_mat);
			point_to_inx_transforms.push_back(transform_point_to_inx);
		}

		void set_value_in_grid(int grid_inx, std::vector<float> position_world, float value){
			openvdb::Vec3d position_index = (point_to_inx_transforms[grid_inx])->worldToIndex(openvdb::Vec3d(position_world[0], position_world[1], position_world[2])); // get index position
			grids_accessors[grid_inx].setValue(openvdb::Coord(round(position_index[0]), round(position_index[1]), round(position_index[2])), value); // set value
		}

		void save_vdb_grids(std::string output_dir_path, int division_start, openvdb::math::Mat4d &grid_transform){
			// print num voxels/max/min in each level
			size_t activeVoxelCountSum = 0;
			for(int lvl = 0; lvl < grids_array.size(); lvl++){
				float min_val, max_val;
				//grids_array[lvl]->evalMinMax(min_val, max_val);
				eval_min_max(grids_array[lvl], min_val, max_val);
				size_t activeVoxelCount = grids_array[lvl]->activeVoxelCount();
				printf("lvl: %d\t voxels:%lld\t values: %f - %f \n", lvl, activeVoxelCount, min_val, max_val);

				activeVoxelCountSum += activeVoxelCount;

				//openvdb::math::Mat4d mat4d = grids_array[lvl]->transform().baseMap()->getAffineMap()->getMat4();
				//openvdb::math::Mat4d matResult = mat4d * grid_transform;

				// Create a new transform from the resulting matrix
				//openvdb::math::Transform::Ptr newTransform = openvdb::math::Transform::createLinearTransform(matResult);
				//grids_array[lvl]->setTransform(newTransform);
			}

			std::cout << "activeVoxelCountSum: " << activeVoxelCountSum << std::endl;
			
			// save each level of resolution (starting from coarsest resolution given by division_start) into one vdb
			openvdb::GridPtrVec grids;
			std::string file_path = output_dir_path + "/grids_" + std::to_string(omp_get_wtime()) + ".vdb";
			printf("Save multi res grids: %s\n", file_path.c_str());
			openvdb::io::File file(file_path);
			for (int i = division_start; i < grids_array.size(); i++){
				size_t activeVoxelCount = grids_array[i]->activeVoxelCount();
				if (activeVoxelCount > 0)
					grids.push_back(grids_array[i]);
			}
			file.write(grids);
    		file.close();	
		}

	};

	void prepare_input_multi_res(openvdb::FloatGrid::Ptr moto_grid, 
		common::vdb::VDBParticles& grid_main_sum, 
		common::SpaceData::DenseType dense_type, 
		std::vector<std::vector<float>>& points,
		std::vector<std::vector<float>>& extended_points, 
		std::vector<std::vector<float>>& data_bbox, 
		std::vector<float>& global_min_precision, 
		float value_precision, 
		int dim_size) {

		// --- create array of 4D points (coordinates in 3D space and value) + extend points to have at most max difference for neighbouring values
#if 0
		if (dense_type != common::SpaceData::DenseType::eNone) { // dense grid
			std::vector<float> input_data_values = grid_main_sum.dense_grid.data_density;
			int max_inx = dim_size - 1;

			std::cout << "Data dim: (" << grid_main_sum.dense_grid.x() << ", " << grid_main_sum.dense_grid.y() << ", " << grid_main_sum.dense_grid.z() << ")\n";

			data_bbox.push_back(std::vector<float>({ 0.0, 0.0, 0.0 }));
			data_bbox.push_back(std::vector<float>({ float(max_inx), float(max_inx), float(max_inx) }));
			global_min_precision = { 1.0 * 2.0, 1.0 * 2.0, 1.0 * 2.0, value_precision };

			float max_value_diff = value_precision / 2.0;

			for (int x = 0; x < dim_size; x++) {
				for (int y = 0; y < dim_size; y++) {
					for (int z = 0; z < dim_size; z++) {
						float point_value = input_data_values[grid_main_sum.dense_grid.get_index(x, y, z)];

						// add 4D point to points
						std::vector<float> point = { float(x), float(y), float(z), point_value };
						points.push_back(point);
						//extended_points.push_back(point);

						// get neighbours
						std::vector<std::vector<float>> neighbours;
						if (z < max_inx) { // +z, up
							float neighbour_voxel_value = input_data_values[grid_main_sum.dense_grid.get_index(x, y, z + 1)];
							std::vector<float> neighbour = { float(x), float(y), float(z + 1), neighbour_voxel_value };
							neighbours.push_back(neighbour);
						}
						if (y < max_inx) { // +y, forward
							float neighbour_voxel_value = input_data_values[grid_main_sum.dense_grid.get_index(x, y + 1, z)];
							std::vector<float> neighbour = { float(x), float(y + 1), float(z), neighbour_voxel_value };
							neighbours.push_back(neighbour);
						}
						if (x < max_inx) { // +x, right
							float neighbour_voxel_value = input_data_values[grid_main_sum.dense_grid.get_index(x + 1, y, z)];
							std::vector<float> neighbour = { float(x + 1), float(y), float(z), neighbour_voxel_value };
							neighbours.push_back(neighbour);
						}

						// for each neighbour: if the difference between the point and the neighbour is bigger than the max allowed difference, create points between the point and the neighbour and add them to extended points
						// TODO: create extended points between the point and neighbour also in other dimensions, not only 4th dimension ("move" them also in other dimensions)
						for (auto neighbour : neighbours) {
							float diff_value = neighbour[3] - point[3];
							float diff_value_abs = fabs(diff_value);

							if (diff_value_abs > max_value_diff) {
								int num_added_points = int(diff_value_abs / max_value_diff);
								if (diff_value > 0) {
									for (int added_inx = 1; added_inx < num_added_points + 1; added_inx++) {
										extended_points.push_back({ point[0], point[1], point[2], point[3] + added_inx * max_value_diff });
									}
								}
								else {
									for (int added_inx = 1; added_inx < num_added_points + 1; added_inx++) {
										extended_points.push_back({ point[0], point[1], point[2], point[3] - added_inx * max_value_diff });
									}
								}
							}
						}
					}
				}
			}

		}
		else 
#endif
		{ // vdg grid
			openvdb::FloatGrid::Ptr input_vdb_grid = grid_main_sum.vdb_grid;

			// ---motorbike
			//openvdb::FloatGrid::Ptr input_vdb_grid = moto_grid;

			// openvdb::CoordBBox inx_bbox = input_vdb_grid->evalActiveVoxelBoundingBox();		
			// openvdb::Coord inx_bbox_min = inx_bbox.min();
			// openvdb::Coord inx_bbox_max = inx_bbox.max();

			openvdb::Coord inx_bbox_min = openvdb::Coord(0, 0, 0);
			openvdb::Coord inx_bbox_max = openvdb::Coord(inx_bbox_min[0] + dim_size - 1, inx_bbox_min[1] + dim_size - 1, inx_bbox_min[2] + dim_size - 1);

			std::cout << "Min Index: (" << inx_bbox_min.x() << ", " << inx_bbox_min.y() << ", " << inx_bbox_min.z() << ")\n";
			std::cout << "Max Index: (" << inx_bbox_max.x() << ", " << inx_bbox_max.y() << ", " << inx_bbox_max.z() << ")\n";

			openvdb::FloatGrid::Accessor vdb_grid_accessor = input_vdb_grid->getAccessor();

			openvdb::Vec3d bbox_world_start = input_vdb_grid->transform().indexToWorld(openvdb::Vec3d(inx_bbox_min.x(), inx_bbox_min.y(), inx_bbox_min.z()));
			openvdb::Vec3d bbox_world_end = input_vdb_grid->transform().indexToWorld(openvdb::Vec3d(inx_bbox_max.x(), inx_bbox_max.y(), inx_bbox_max.z()));
			data_bbox.push_back(std::vector<float>({ float(bbox_world_start.x()), float(bbox_world_start.y()), float(bbox_world_start.z()) }));
			data_bbox.push_back(std::vector<float>({ float(bbox_world_end.x()), float(bbox_world_end.y()), float(bbox_world_end.z()) }));

			openvdb::Vec3d bbox_world_start_neighbour = input_vdb_grid->transform().indexToWorld(openvdb::Vec3d(inx_bbox_min.x() + 1, inx_bbox_min.y() + 1, inx_bbox_min.z() + 1));
			global_min_precision = { float((bbox_world_start_neighbour[0] - bbox_world_start[0]) * 2.0), float((bbox_world_start_neighbour[1] - bbox_world_start[1]) * 2.0), float((bbox_world_start_neighbour[2] - bbox_world_start[2]) * 2.0), value_precision };

			float max_value_diff = value_precision / 2.0;

			for (int x = inx_bbox_min.x(); x <= inx_bbox_max.x(); x++) {
				for (int y = inx_bbox_min.y(); y <= inx_bbox_max.y(); y++) {
					for (int z = inx_bbox_min.z(); z <= inx_bbox_max.z(); z++) {
						openvdb::Coord voxel_inx(x, y, z);
						float voxel_value = vdb_grid_accessor.getValue(voxel_inx);
						openvdb::Vec3d voxel_world_center = input_vdb_grid->transform().indexToWorld(openvdb::Vec3d(voxel_inx.x(), voxel_inx.y(), voxel_inx.z()));

						//printf("inx: <%3d %3d %3d>\tworld center: <%.3f %.3f %.3f> == %.3f\n", voxel_inx.x(), voxel_inx.y(), voxel_inx.z(), voxel_world_center.x(), voxel_world_center.y(), voxel_world_center.z(), voxel_value);

						// add 4D point to points
						std::vector<float> point = { float(voxel_world_center.x()), float(voxel_world_center.y()), float(voxel_world_center.z()), voxel_value };
						//printf("point <%.3f, %.3f, %.3f, %.3f>\n", point[0], point[1], point[2], point[3]);
						points.push_back(point);
						//extended_points.push_back(point);

						// get neighbours
						std::vector<std::vector<float>> neighbours;
						if (z < inx_bbox_max.z()) { // +z, up
							openvdb::Coord neighbour_voxel_inx(x, y, z + 1);
							float neighbour_voxel_value = vdb_grid_accessor.getValue(neighbour_voxel_inx);
							openvdb::Vec3d neighbour_voxel_world_center = input_vdb_grid->transform().indexToWorld(openvdb::Vec3d(neighbour_voxel_inx.x(), neighbour_voxel_inx.y(), neighbour_voxel_inx.z()));
							std::vector<float> neighbour = { float(neighbour_voxel_world_center.x()), float(neighbour_voxel_world_center.y()), float(neighbour_voxel_world_center.z()), neighbour_voxel_value };
							//printf("neighbour <%.3f, %.3f, %.3f, %.3f>\n", neighbour[0], neighbour[1], neighbour[2], neighbour[3]);
							neighbours.push_back(neighbour);
						}
						if (y < inx_bbox_max.y()) { // +y, forward
							openvdb::Coord neighbour_voxel_inx(x, y + 1, z);
							float neighbour_voxel_value = vdb_grid_accessor.getValue(neighbour_voxel_inx);
							openvdb::Vec3d neighbour_voxel_world_center = input_vdb_grid->transform().indexToWorld(openvdb::Vec3d(neighbour_voxel_inx.x(), neighbour_voxel_inx.y(), neighbour_voxel_inx.z()));
							std::vector<float> neighbour = { float(neighbour_voxel_world_center.x()), float(neighbour_voxel_world_center.y()), float(neighbour_voxel_world_center.z()), neighbour_voxel_value };
							//printf("neighbour <%.3f, %.3f, %.3f, %.3f>\n", neighbour[0], neighbour[1], neighbour[2], neighbour[3]);
							neighbours.push_back(neighbour);
						}
						if (x < inx_bbox_max.x()) { // +x, right
							openvdb::Coord neighbour_voxel_inx(x + 1, y, z);
							float neighbour_voxel_value = vdb_grid_accessor.getValue(neighbour_voxel_inx);
							openvdb::Vec3d neighbour_voxel_world_center = input_vdb_grid->transform().indexToWorld(openvdb::Vec3d(neighbour_voxel_inx.x(), neighbour_voxel_inx.y(), neighbour_voxel_inx.z()));
							std::vector<float> neighbour = { float(neighbour_voxel_world_center.x()), float(neighbour_voxel_world_center.y()), float(neighbour_voxel_world_center.z()), neighbour_voxel_value };
							//printf("neighbour <%.3f, %.3f, %.3f, %.3f>\n", neighbour[0], neighbour[1], neighbour[2], neighbour[3]);
							neighbours.push_back(neighbour);
						}

						// for each neighbour: if the difference between the point and the neighbour is bigger than the max allowed difference, create points between the point and the neighbour and add them to extended points
						// TODO: create extended points between the point and neighbour also in other dimensions, not only 4th dimension ("move" them also in other dimensions)
						for (auto neighbour : neighbours) {
							float diff_value = neighbour[3] - point[3];
							float diff_value_abs = fabs(diff_value);
							//printf("diff: <%.3f>\n", diff_value);
							//printf("diff_abs: <%.3f>\n", diff_value_abs);

							if (diff_value_abs > max_value_diff) {
								//printf("diff: <%.3f>\n", diff_value);
								int num_added_points = int(diff_value_abs / max_value_diff);
								if (diff_value > 0) {
									for (int added_inx = 1; added_inx < num_added_points + 1; added_inx++) {
										extended_points.push_back({ point[0], point[1], point[2], point[3] + added_inx * max_value_diff });
									}
								}
								else {
									for (int added_inx = 1; added_inx < num_added_points + 1; added_inx++) {
										extended_points.push_back({ point[0], point[1], point[2], point[3] - added_inx * max_value_diff });
									}
								}
							}
						}
					}
				}
			}
		}

		printf("extended points size: %lld\t points size: %lld\n", extended_points.size(), points.size());

	}

	void save_voxel(Openvdb_Multi_Res_Grids& openvdb_grids, std::vector<float> voxel_world_position, std::vector<std::vector<float>> points, int division_inx){
		float average_value = 0.0;
		if(points.size() != 0){
			float sum_value = 0.0;
			for(int i = 0; i < points.size(); i++){
				sum_value += points[i][3];
			}
			average_value = sum_value / points.size();
		}

		if(average_value > FDATA_EPSILON){
#if 0
			// difference computation
			for(auto point: points){
				float diff = point[3] - average_value;
				float diff_norm = diff/openvdb_grids.global_max_value; // normalize
				float abs_diff = abs(diff);
				float squared_diff = diff*diff;
				float squared_diff_norm = diff_norm*diff_norm;
				openvdb_grids.sum_diff_abs += abs_diff;
				openvdb_grids.sum_diff_squared += squared_diff;
				openvdb_grids.sum_diff_squared_norm += squared_diff_norm;
				openvdb_grids.num_voxels_for_metrics++;
			}
#endif	
			openvdb_grids.set_value_in_grid(division_inx, voxel_world_position, average_value);
		}

		// else {
		// 	printf("Voxel value is 0.0.");
		// }	

	}

	// number of boxes needed to cover the data (data min in each dim is 0.0)
	int box_count(std::vector<std::vector<float>> points, std::vector<float> box_size, int size_inx){
		// set max inx of the box in each dim
		int max_box_inx = std::pow(2, size_inx) - 1;
		// create set of boxes - for each point compute inx of the box and add it to the set
		std::set<std::array<int, 4>> boxes;
		for(auto point: points){
			std::array<int, 4> point_box;
			for(int dim_i = 0; dim_i < point.size(); dim_i++){
				int box_inx = int((point[dim_i]) / box_size[dim_i]); // get inx of the box for the point in given dimension
				point_box[dim_i] = (box_inx <= max_box_inx) ?  box_inx: max_box_inx; // add box inx or maxium box inx for given dimension
			}

			boxes.insert(point_box);
		}
		
		return boxes.size(); // return number of boxes
	}

	// number of boxes needed to cover the data (data min in each dim is 0.0)  const std::vector<int>&
	int box_count_2D(const std::vector<std::vector<float>>& points, const std::vector<std::vector<float>>& extended_points, std::vector<double> box_size, std::vector<float> voxel_size, int size_inx, double time_now, std::string folder_name, std::string copy_folder){
		// set max inx of the box in each dim
		int max_box_inx = std::pow(2, size_inx) - 1;
		// create set of boxes - for each point compute inx of the box and add it to the set
		std::set<std::array<int, 2>> boxes;
		for(int i = 0; i < points.size(); i++){
			std::array<int, 2> point_box;
			int box_inx_x = int(points[i][0] / box_size[0]); // get inx of the box for the point in given dimension
			point_box[0] = (box_inx_x <= max_box_inx) ?  box_inx_x: max_box_inx; // add box inx or maxium box inx for given dimension

			int box_inx_y = int(float(points[i][3]) / box_size[1]); // get inx of the box for the point in given dimension
			point_box[1] = (box_inx_y <= max_box_inx) ?  box_inx_y: max_box_inx; // add box inx or maxium box inx for given dimension

			boxes.insert(point_box);
		}

		for(int i = 0; i < extended_points.size(); i++){
			std::array<int, 2> point_box;
			int box_inx_x = int(extended_points[i][0] / box_size[0]); // get inx of the box for the point in given dimension
			point_box[0] = (box_inx_x <= max_box_inx) ?  box_inx_x: max_box_inx; // add box inx or maxium box inx for given dimension

			int box_inx_y = int(float(extended_points[i][3]) / box_size[1]); // get inx of the box for the point in given dimension
			point_box[1] = (box_inx_y <= max_box_inx) ?  box_inx_y: max_box_inx; // add box inx or maxium box inx for given dimension

			boxes.insert(point_box);
		}
		
		// // boxes output file
		// std::string file_name = "boxes_" + std::to_string(size_inx) + ".txt";
		// std::ofstream rectangle_file(folder_name + "/" + file_name);
		// for (auto box: boxes) {
		// 	float start_x = box[0] * box_size[0];
		// 	float start_y = box[1] * box_size[1];

		// 	rectangle_file << start_x << " " << start_y << " " << box_size[0] << " " << box_size[1] << "\n";
		// }
		// rectangle_file.close();
		// // copy to folder "last"
		// std::filesystem::copy_file(folder_name + "/" + file_name, copy_folder + "/" + file_name, std::filesystem::copy_options::overwrite_existing);


		return boxes.size(); // return number of boxes
	}

	double compute_slope(const std::vector<double>& x, const std::vector<double>& y) {
		int n = x.size();
		double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

		for (int i = 0; i < n; i++) {
			sumX += x[i];
			sumY += y[i];
			sumXY += x[i] * y[i];
			sumXX += x[i] * x[i];
		}

		// Compute slope (m)
		double slope = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
		return slope;
	}


	double fract_dim(std::vector<std::vector<float>> points, std::vector<std::vector<float>> extended_points, std::vector<float> global_min_precision, int points_orig_size, std::string out_full_filepath, std::vector<float> voxel_size, std::vector<int> num_voxels, int division_current){
		// --- move points to have min values [0,0,0,0] + convert to 2D
		int points_dim = points[0].size();

		// find min and max value in each dimension
		std::vector<float> min_values({points[0][0], points[0][1], points[0][2], points[0][3]});
		std::vector<float> max_values({points[0][0], points[0][1], points[0][2], points[0][3]});
		for (auto point: points) {
			for (int dim = 0; dim < points_dim; dim++) {
				min_values[dim] = std::min(min_values[dim], point[dim]);
				max_values[dim] = std::max(max_values[dim], point[dim]);
			}
		}
		// move points and extended points based on found min value + save in x/[0] value the position in 1D
		for(int i = 0; i < points.size(); i++){
			points[i] = {points[i][0] - min_values[0], points[i][1] - min_values[1], points[i][2] - min_values[2], points[i][3] - min_values[3]};
			points[i][0] = (points[i][2]/voxel_size[2]) * num_voxels[0] * num_voxels[1] + (points[i][1]/voxel_size[1]) * num_voxels[0]  + (points[i][0]/voxel_size[0]);
			//points[i][0] = i; //test
		}
		for(int i = 0; i < extended_points.size(); i++){
			extended_points[i] = {extended_points[i][0] - min_values[0], extended_points[i][1] - min_values[1], extended_points[i][2] - min_values[2], extended_points[i][3] - min_values[3]};
			extended_points[i][0] = (extended_points[i][2]/voxel_size[2]) * num_voxels[0] * num_voxels[1] + (extended_points[i][1]/voxel_size[1]) * num_voxels[0]  + (extended_points[i][0]/voxel_size[0]);
			//extended_points[i][0] = i; //test
		}

		// --- box count output files
		double time_now = omp_get_wtime();
		std::string file_name = "points_" + std::to_string(time_now) + ".txt";
		std::string folder_name = out_full_filepath + "/out_multi_res/docs/box_count/" + std::to_string(division_current) + "/";
		std::string copy_folder = out_full_filepath + "/out_multi_res/docs/box_count/last";

		// std::filesystem::create_directory(folder_name);
		// // create file with points (the ones used for box count computation)
		// std::ofstream points_file(folder_name + "/" + file_name);
		// for (int i = 0; i < points.size(); i++) {
		// 	points_file << points[i][0] << " " << points[i][3] << "\n";
		// }
		// points_file.close();
		// // copy to folder "last"
		
		// std::filesystem::create_directory(copy_folder);
		// std::filesystem::copy_file(folder_name + "/" + file_name, copy_folder + "/" + file_name, std::filesystem::copy_options::overwrite_existing);

		// // --- check if data range (max) is bigger than global min precision
		// for(int dim = 0; dim < points_dim; dim++){
		// 	// move max values based on found min value
		// 	max_values[dim] = max_values[dim] - min_values[dim];
		// 	// if max value is smaller then global min precision -> data are smaller than desired data resolution -> f dim is 1/2/3 (line/surface)
		// 	if(max_values[dim] < global_min_precision[dim]){
		// 		return 3.0;
		// 	}
		// }

		max_values[3] = max_values[3] - min_values[3];
		float size_1d = float((num_voxels[0] * num_voxels[1] * num_voxels[2]) - 1);
		if(max_values[3] < global_min_precision[3] || size_1d < 2.0){
			return 1.0;
		}

		// // --- compute sizes - take max sizes and divide them by 2 to gain smaller sizes (stop when box size is smaller then precision)
		// std::vector<std::vector<double>> box_sizes;
		// box_sizes.push_back(max_values);
		// while(true){
		// 	std::vector<double> new_size;
		// 	for(int dim = 0; dim < points_dim; dim++){
		// 		new_size.push_back(box_sizes.back()[dim]/2.0);
		// 	}
		// 	if (new_size[3] < global_min_precision[3]/2.0 || new_size[2] < global_min_precision[2]/2.0 || new_size[1] < global_min_precision[1]/2.0 || new_size[0] < global_min_precision[0]/2.0){
		// 		if(box_sizes.size() >= 2){
		// 			break;
		// 		} else{
		// 			// TODO odstranit, pokud je otestovano, ze se to nedeje
		// 			printf("Err: this should not happen");
		// 			return -1.0;
		// 		}
		// 	} else{
		// 		box_sizes.push_back(new_size);
		// 	}
		// }

		// --- compute sizes - take max sizes and divide them by 2 to gain smaller sizes (stop when box size is smaller then precision)
		std::vector<std::vector<double>> box_sizes;
		box_sizes.push_back({size_1d, max_values[3]});
		//box_sizes.push_back({double(points.size() - 1), max_values[3]}); //test
		while(true){
			std::vector<double> new_size;
			for(int dim = 0; dim < box_sizes[0].size(); dim++){
				new_size.push_back(box_sizes.back()[dim]/2.0);
			}
			if (new_size[1] < global_min_precision[3]/2.0 || new_size[0] < 1.0){
				if(box_sizes.size() >= 2){
					break;
				} else{
					// TODO odstranit, pokud je otestovano, ze se to nedeje
					printf("Err: this should not happen.\n");
					return -1.0;
				}
			} else{
				box_sizes.push_back(new_size);
			}
		}

		// --- box count
		std::vector<int> counts;
		for(int size_i = 0; size_i < box_sizes.size(); size_i++){ // TODO skip size 0, which covers all data (box count is 1)
			//counts.push_back(box_count(points, box_sizes[size_i], size_i));
			counts.push_back(box_count_2D(points, extended_points, box_sizes[size_i], voxel_size, size_i, time_now, folder_name, copy_folder));
		}

		// --- set box sizes for fract dim computation
		std::vector<double> sizes = {1.0};
		for(int i = 0; i < (box_sizes.size() - 1); i++){
			sizes.push_back(sizes.back() / 2.0);
		}
			

		// --- Linear fit to determine the fractal dimension
		std::vector<double> x;
		std::vector<double> log_counts;
		for(int i = 0; i < sizes.size(); i++){
			x.push_back(std::log(1/sizes[i]));
			log_counts.push_back(std::log(counts[i]));
		}

		double fractal_dimension = compute_slope(x, log_counts); // TODO - nahradit funkci z knihovny

		return fractal_dimension;
	}
		
	void divide_grid_multi_res(Openvdb_Multi_Res_Grids& openvdb_grids, std::vector<std::vector<float>> grid_bbox, std::vector<std::vector<float>> points, std::vector<std::vector<float>> extended_points, int division_current, int division_max, int division_start, float f_dim_threshold, std::vector<float> global_min_precision, std::string out_full_filepath, std::vector<float> voxel_size, std::vector<int> num_points_3D){

		// stop if max division reached
		if (division_current == division_max){
			save_voxel(openvdb_grids, grid_bbox[0], points, division_current);
			return;
		}
			
		// stop if points is empty
		if(points.size() == 0){
			printf("ERR: points is empty\n");
			save_voxel(openvdb_grids, grid_bbox[0], points, division_current);
			return;
		}

		// stop if points size is one
		if(points.size() == 1){
			//printf("div %d: save one point\n", division_current);
			save_voxel(openvdb_grids, grid_bbox[0], points, division_current);
			return;
		}
			
		// compute fract dim for div start and higher
		double f_dim;
		if (division_current >= division_start){
			// double x = ((grid_bbox[1][0] - grid_bbox[0][0])/voxel_size[0]) + 1;
			// double y = ((grid_bbox[1][1] - grid_bbox[0][1])/voxel_size[1]) + 1;
			// double z = ((grid_bbox[1][2] - grid_bbox[0][2])/voxel_size[2]) + 1;
			// std::vector<int> num_voxels({int(x), int(y), int(z)});
			if(num_points_3D[0] * num_points_3D[1] * num_points_3D[2] != points.size()){
				printf("Err: wrong number of points 3D for fd computation.");
				return;
			}
			f_dim = fract_dim(points, extended_points, global_min_precision, points.size(), out_full_filepath, voxel_size, num_points_3D, division_current);
			//printf("div: %d\tfd: %f\n", division_current, f_dim);
			if(division_current == division_start){
				printf("%f\n", f_dim);
				//printf("div: %d\tfd: %f\n", division_current, f_dim);
			}
		} else {
			f_dim = 3.9;
		}
			

		// divide if f_dim is above threshold
		if (f_dim >= f_dim_threshold){
		//if (true){
			division_current += 1;

			// divide points
			std::vector<float> grid_mid;
			std::vector<float> new_grid_size;
			for(int dim_i = 0; dim_i < grid_bbox[0].size(); dim_i++){
				grid_mid.push_back((grid_bbox[0][dim_i] + grid_bbox[1][dim_i]) / 2.0);
				new_grid_size.push_back(grid_mid[dim_i] - grid_bbox[0][dim_i]);
			}

			std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> new_grids_points(2, std::vector<std::vector<std::vector<std::vector<float>>>>(2, std::vector<std::vector<std::vector<float>>>(2)));  // 2x2x2 array of arrays
			std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> new_grids_points_extended(2, std::vector<std::vector<std::vector<std::vector<float>>>>(2, std::vector<std::vector<std::vector<float>>>(2)));  // 2x2x2 array of arrays
			std::vector<std::vector<std::vector<std::vector<std::set<float>>>>> new_grids_num_points_3D_set(2, std::vector<std::vector<std::vector<std::set<float>>>>(2, std::vector<std::vector<std::set<float>>>(2, std::vector<std::set<float>>(3))));  // 2x2x2 array of sets
			std::vector<int> new_grid_inx({0,0,0});

			for(auto point: points){
				for(int dim_i = 0; dim_i < point.size() - 1; dim_i++){
					if(point[dim_i] < grid_mid[dim_i]){
						new_grid_inx[dim_i] = 0;
					}	
					else{
						new_grid_inx[dim_i] = 1;
					}	
				}
				new_grids_points[new_grid_inx[0]][new_grid_inx[1]][new_grid_inx[2]].push_back(point);

				// fill num voxels in each dimension in the new regions
				new_grids_num_points_3D_set[new_grid_inx[0]][new_grid_inx[1]][new_grid_inx[2]][0].insert(point[0]);
				new_grids_num_points_3D_set[new_grid_inx[0]][new_grid_inx[1]][new_grid_inx[2]][1].insert(point[1]);
				new_grids_num_points_3D_set[new_grid_inx[0]][new_grid_inx[1]][new_grid_inx[2]][2].insert(point[2]);					

			}

			for(auto point: extended_points){
				for(int dim_i = 0; dim_i < point.size() - 1; dim_i++){
					if(point[dim_i] < grid_mid[dim_i]){
						new_grid_inx[dim_i] = 0;
					}	
					else{
						new_grid_inx[dim_i] = 1;
					}	
				}
				new_grids_points_extended[new_grid_inx[0]][new_grid_inx[1]][new_grid_inx[2]].push_back(point);
			}

			// clear whole points and extended points
			std::vector<std::vector<float>>().swap(points);
			std::vector<std::vector<float>>().swap(extended_points);

			for(int i = 0; i < new_grids_points.size(); i++){
				for(int j = 0; j < new_grids_points[i].size(); j++){
					for(int k = 0; k < new_grids_points[i][j].size(); k++){

						std::vector<std::vector<float>> new_grid_bb({{grid_bbox[0][0] + i * new_grid_size[0], grid_bbox[0][1] + j * new_grid_size[1], grid_bbox[0][2] + k * new_grid_size[2]}, 
						{grid_mid[0] + i * new_grid_size[0], grid_mid[1] + j * new_grid_size[1], grid_mid[2] + k * new_grid_size[2]}});
						
						std::vector<int> new_grids_num_points_3D({int(new_grids_num_points_3D_set[i][j][k][0].size()), int(new_grids_num_points_3D_set[i][j][k][1].size()), int(new_grids_num_points_3D_set[i][j][k][2].size())});

						divide_grid_multi_res(openvdb_grids, new_grid_bb, new_grids_points[i][j][k], new_grids_points_extended[i][j][k], division_current, division_max, division_start, f_dim_threshold, global_min_precision, out_full_filepath, voxel_size, new_grids_num_points_3D);
					}
				}	
			}
				
		} else{
			save_voxel(openvdb_grids, grid_bbox[0], points, division_current);
			return;
		}
	}

	std::vector<double> uniform_filter(const std::vector<double>& image, int win_size, int dim_size) {
		//int rows = image.size(), cols = image[0].size();
		std::vector<double> result(image.size(), 0.0);
		int half_win = win_size / 2;

#ifdef WITH_OPENMP
# pragma omp parallel for
#endif    
		for (int i = half_win; i < dim_size - half_win; ++i) {
			for (int j = half_win; j < dim_size - half_win; ++j) {
				for (int k = half_win; k < dim_size - half_win; ++k) {
					double sum = 0.0;
					for (int m = -half_win; m <= half_win; ++m) {
						for (int n = -half_win; n <= half_win; ++n) {
							for (int u = -half_win; u <= half_win; ++u) {
								size_t im_id = size_t(i + m) + size_t(j + n) * dim_size + size_t(k + u) * dim_size * dim_size;
								sum += image[im_id];
							}
						}
					}

					size_t res_id = size_t(i) + size_t(j) * dim_size + size_t(k) * dim_size * dim_size;
					result[res_id] = sum / (win_size * win_size);
				}
			}
		}
		return result;
	}

	// Helper function to set up cubic interpolation weights
	inline void set_cubic_spline_weights(float w[4], float t) {
		float t2 = t * t;
		float t3 = t2 * t;
		
		w[0] = -1.0f/6.0f * t3 + 0.5f * t2 - 0.5f * t + 1.0f/6.0f;
		w[1] = 0.5f * t3 - t2 + 2.0f/3.0f;
		w[2] = -0.5f * t3 + 0.5f * t2 + 0.5f * t + 1.0f/6.0f;
		w[3] = 1.0f/6.0f * t3;
	}		

	double get_mres_value3(openvdb::FloatGrid::Ptr grid_reference, std::vector<openvdb::FloatGrid::Ptr>& multi_grids, openvdb::Coord coord) {
		double acc_value_sum = 0.0;
	
		// Convert reference coordinate to world space
		openvdb::Vec3d coord_w = grid_reference->transform().indexToWorld(coord);
	
		// For each grid in the multi-resolution hierarchy
		for (size_t i = 0; i < multi_grids.size(); ++i) {
			openvdb::FloatGrid::Ptr grid2 = multi_grids[i];
			
			// Convert world space to target grid index space
			openvdb::Vec3d coord_grid2 = grid2->transform().worldToIndex(coord_w);
			
			// Extract the integer and fractional parts
			int ix = int(std::floor(coord_grid2.x()));
			int iy = int(std::floor(coord_grid2.y()));
			int iz = int(std::floor(coord_grid2.z()));
			
			float tx = coord_grid2.x() - ix;
			float ty = coord_grid2.y() - iy;
			float tz = coord_grid2.z() - iz;
			
			// Get grid dimensions
			openvdb::Coord dim = grid2->evalActiveVoxelDim();
			int width = dim.x();
			int height = dim.y();
			int depth = dim.z();
			
			// Check if any cubic samples would be outside the grid
			if (ix < -1 || ix > width || iy < -1 || iy > height || iz < -1 || iz > depth) {
				// Point is outside or too close to the edge for cubic interpolation
				continue;
			}
			
			// Compute cubic interpolation sample positions
			int pix = ix - 1;
			int nix = ix + 1;
			int nnix = ix + 2;
			
			int piy = iy - 1;
			int niy = iy + 1;
			int nniy = iy + 2;
			
			int piz = iz - 1;
			int niz = iz + 1;
			int nniz = iz + 2;
			
			// If all cubic samples are inside the grid, perform full cubic interpolation
			if (pix >= 0 && nnix < width && piy >= 0 && nniy < height && piz >= 0 && nniz < depth) {
				// Create accessor for efficient grid access
				openvdb::FloatGrid::Accessor acc = grid2->getAccessor();
				
				// Setup cubic spline weights
				float u[4], v[4], w[4];
				set_cubic_spline_weights(u, tx);
				set_cubic_spline_weights(v, ty);
				set_cubic_spline_weights(w, tz);
				
				// Perform tricubic interpolation
				float value = 0.0f;
				
				// Triple nested loop for 4x4x4 neighborhood
				for (int k = 0; k < 4; k++) {
					int z_pos = (k == 0) ? piz : (k == 1) ? iz : (k == 2) ? niz : nniz;
					
					for (int j = 0; j < 4; j++) {
						int y_pos = (j == 0) ? piy : (j == 1) ? iy : (j == 2) ? niy : nniy;
						float v_weight = v[j];
						
						for (int i = 0; i < 4; i++) {
							int x_pos = (i == 0) ? pix : (i == 1) ? ix : (i == 2) ? nix : nnix;
							
							float sample = acc.getValue(openvdb::Coord(x_pos, y_pos, z_pos));
							value += u[i] * v_weight * w[k] * sample;
						}
					}
				}
				
				acc_value_sum += value;
			}
			else {
				// Edge case: Use OpenVDB's built-in interpolation
				openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*grid2);
				acc_value_sum += sampler.wsSample(coord_w);
			}
		}
	
		return acc_value_sum;
	}

	double get_mres_value(openvdb::FloatGrid::Ptr grid_reference, std::vector<openvdb::FloatGrid::Ptr>& multi_grids, openvdb::Coord coord) {
		double acc_value_sum = 0.0;

		//printf("%d: (", num_voxels);
		for (size_t i = 0; i < multi_grids.size(); ++i) {
			openvdb::FloatGrid::Ptr grid2 = multi_grids[i];
			openvdb::FloatGrid::Accessor acc2 = grid2->getAccessor();

			openvdb::Vec3d coord_w = grid_reference->transform().indexToWorld(coord);
			openvdb::Vec3d coord_grid2 = grid2->transform().worldToIndex(coord_w);
			openvdb::Coord coord_grid2_int = openvdb::Coord(int(coord_grid2.x()), int(coord_grid2.y()), int(coord_grid2.z()));
			acc_value_sum += acc2.getValue(coord_grid2_int);
		}

		return acc_value_sum;
	}

	double trilinear_interpolation(
		float V000, float V100, float V010, float V110,
		float V001, float V101, float V011, float V111,
		float x, float y, float z
	) {
		// Interpolace pod�l osy x
		double V00 = (1.0 - x) * V000 + x * V100;
		double V10 = (1.0 - x) * V010 + x * V110;
		double V01 = (1.0 - x) * V001 + x * V101;
		double V11 = (1.0 - x) * V011 + x * V111;

		// Interpolace pod�l osy y
		double V0 = (1.0 - y) * V00 + y * V10;
		double V1 = (1.0 - y) * V01 + y * V11;

		// Interpolace pod�l osy z
		double V = (1.0 - z) * V0 + z * V1;

		return V;
	}

	double get_mres_value4(openvdb::FloatGrid::Ptr grid_reference, std::vector<openvdb::FloatGrid::Ptr>& multi_grids, openvdb::Coord coord) {
		double acc_value_sum = 0.0;
		int count = 0;

		for (int dx = -1; dx <= 1; ++dx) {
			for (int dy = -1; dy <= 1; ++dy) {
				for (int dz = -1; dz <= 1; ++dz) {
					// TODO
					//if (dx == 1 && dy == 1 && dz == 1) continue;
					//acc_value_sum += get_mres_value2(grid_reference, multi_grids, openvdb::Coord(coord.x() + dx, coord.y() + dy, coord.z() + dz));
					++count;
				}
			}
		}

		return acc_value_sum / (double)count;
	}

	double compute_ssim(openvdb::FloatGrid::Ptr grid_reference, std::vector<openvdb::FloatGrid::Ptr>& multi_grids, int dim_size, float min_value, float max_value, float min_value_mres, float max_value_mres,
		double data_range = 1.0, int win_size = 7, double K1 = 0.01, double K2 = 0.03) {
		size_t N = dim_size * size_t(dim_size) * dim_size;
		std::vector<double> im1(N);
		std::vector<double> im2(N);

		//float minValue = FLT_MAX;
		//float maxValue = -FLT_MAX;
		//float minVDB = FLT_MAX;
		//float maxVDB = -FLT_MAX;

		openvdb::FloatGrid::Accessor reference_grid_accessor = grid_reference->getAccessor();

		// compute min/max's
//#ifdef WITH_OPENMP
//# pragma omp parallel for reduction(min: minValue, minVDB) reduction(max: maxValue, maxVDB)
//#endif  
//		for (int z = 0; z < dim_size; ++z) {
//			for (int y = 0; y < dim_size; ++y) {
//				for (int x = 0; x < dim_size; ++x) {
//					openvdb::Coord coord(x, y, z);
//					float value0 = reference_grid_accessor.getValue(coord);
//
//					minValue = fminf(minValue, value0);
//					maxValue = fmaxf(maxValue, value0);
//
//					float value1 = get_mres_value(grid_reference, multi_grids, coord);
//					minVDB = fminf(minVDB, value1);
//					maxVDB = fmaxf(maxVDB, value1);
//				}
//			}
//		}

#ifdef WITH_OPENMP
# pragma omp parallel for
#endif  
		for (int k = 0; k < dim_size; ++k) {
			for (int j = 0; j < dim_size; ++j) {
				for (int i = 0; i < dim_size; ++i) {
					// float value0 = getValue(input,i,j,k);
					// float value1 = getValue(comp,i,j,k);

					openvdb::Coord coord(i, j, k);
					float value0 = reference_grid_accessor.getValue(coord);
					value0 = (value0 - min_value) / (max_value - min_value);

					float value1 = get_mres_value(grid_reference, multi_grids, coord);
					value1 = (value1 - min_value_mres) / (max_value_mres - min_value_mres);

					size_t id = size_t(i) + size_t(j) * dim_size + size_t(k) * dim_size * dim_size;
					im1[id] = value0;
					im2[id] = value1;
				}
			}
		}

		//int rows = im1.size(), cols = im1[0].size();
		double C1 = (K1 * data_range) * (K1 * data_range);
		double C2 = (K2 * data_range) * (K2 * data_range);

		auto mu1 = uniform_filter(im1, win_size, dim_size);
		auto mu2 = uniform_filter(im2, win_size, dim_size);

		std::vector<double> mu1_sq(im1.size(), 0.0), mu2_sq(im2.size(), 0.0), mu1_mu2(im1.size(), 0.0);

#ifdef WITH_OPENMP
# pragma omp parallel for
#endif    
		for (int i = 0; i < dim_size; ++i) {
			for (int j = 0; j < dim_size; ++j) {
				for (int k = 0; k < dim_size; ++k) {
					size_t id = size_t(i) + size_t(j) * dim_size + size_t(k) * dim_size * dim_size;
					mu1_sq[id] = mu1[id] * mu1[id];
					mu2_sq[id] = mu2[id] * mu2[id];
					mu1_mu2[id] = mu1[id] * mu2[id];
				}
			}
		}

		auto sigma1_sq = uniform_filter(im1, win_size, dim_size);
		auto sigma2_sq = uniform_filter(im2, win_size, dim_size);
		auto sigma12 = uniform_filter(im1, win_size, dim_size);

#ifdef WITH_OPENMP
# pragma omp parallel for
#endif    
		for (int i = 0; i < dim_size; ++i) {
			for (int j = 0; j < dim_size; ++j) {
				for (int k = 0; k < dim_size; ++k) {
					size_t id = size_t(i) + size_t(j) * dim_size + size_t(k) * dim_size * dim_size;
					sigma1_sq[id] -= mu1_sq[id];
					sigma2_sq[id] -= mu2_sq[id];
					sigma12[id] -= mu1_mu2[id];
				}
			}
		}

		double ssim_sum = 0.0;
		size_t count = 0;

#ifdef WITH_OPENMP
# pragma omp parallel for reduction(+: ssim_sum, count)
#endif    
		for (int i = win_size / 2; i < dim_size - win_size / 2; ++i) {
			for (int j = win_size / 2; j < dim_size - win_size / 2; ++j) {
				for (int k = win_size / 2; k < dim_size - win_size / 2; ++k) {
					size_t id = size_t(i) + size_t(j) * dim_size + size_t(k) * dim_size * dim_size;

					double numerator = (2 * mu1_mu2[id] + C1) * (2 * sigma12[id] + C2);
					double denominator = (mu1_sq[id] + mu2_sq[id] + C1) * (sigma1_sq[id] + sigma2_sq[id] + C2);
					ssim_sum += numerator / denominator;
					count++;
				}
			}
		}

		double ssim = ssim_sum / count;
		if (ssim > 1.0) // TODO
			ssim = 2.0 - ssim;

		return ssim;
	}

	struct Stats
	{
		float minValue{ FLT_MAX }, maxValue{ -FLT_MAX };
		float minVDB{ FLT_MAX }, maxVDB{ -FLT_MAX };

		double mse{ 0.0 };
		double snr{ 0.0 };
		double psnr{ 0.0 };
	};

	Stats compute_stats(openvdb::FloatGrid::Ptr grid_reference, std::vector<openvdb::FloatGrid::Ptr>& multi_grids, int dim_size)
	{
		Stats res;

		float minValue = FLT_MAX;
		float maxValue = -FLT_MAX;
		float minVDB = FLT_MAX;
		float maxVDB = -FLT_MAX;

		openvdb::FloatGrid::Accessor reference_grid_accessor = grid_reference->getAccessor();

		// compute min/max's
#ifdef WITH_OPENMP
# pragma omp parallel for reduction(min: minValue, minVDB) reduction(max: maxValue, maxVDB)
#endif  
		for (int z = 0; z < dim_size; ++z) {
			for (int y = 0; y < dim_size; ++y) {
				for (int x = 0; x < dim_size; ++x) {
					openvdb::Coord coord(x, y, z);
					float value0 = reference_grid_accessor.getValue(coord);

					minValue = fminf(minValue, value0);
					maxValue = fmaxf(maxValue, value0);

					float value1 = get_mres_value(grid_reference, multi_grids, coord);
					minVDB = fminf(minVDB, value1);
					maxVDB = fmaxf(maxVDB, value1);
				}
			}
		}

		res.minValue = minValue;
		res.maxValue = maxValue;
		res.minVDB = minVDB;
		res.maxVDB = maxVDB;

		double sumSquared{ 0.0 };
		double sumSquaredErr{ 0.0 };


#ifdef WITH_OPENMP
# pragma omp parallel for reduction(+: sumSquared, sumSquaredErr)
#endif  
		for (int z = 0; z < dim_size; ++z) {
			for (int y = 0; y < dim_size; ++y) {
				for (int x = 0; x < dim_size; ++x) {
					openvdb::Coord coord(x, y, z);
					float value0 = reference_grid_accessor.getValue(coord);
					value0 = (value0 - minValue) / (maxValue - minValue);

					float value1 = get_mres_value(grid_reference, multi_grids, coord);
					value1 = (value1 - minVDB) / (maxVDB - minVDB);

					double sqr = double(value0) * double(value0);
					double diff = double(value0) - double(value1);
					sumSquared += sqr;
					sumSquaredErr += diff * diff;
				}
			}
		}

		size_t N = dim_size * size_t(dim_size) * dim_size;
		
		res.mse = sumSquaredErr / N;
		double signalMean = sumSquared / N;
		double noiseMean = res.mse;
		if (noiseMean == 0.0) {
			res.snr = INFINITY;
			res.psnr = INFINITY;
		}
		else {
			res.snr = 20 * log10(sqrt(signalMean) / sqrt(noiseMean));
			res.psnr = 10 * log10(1.0 / noiseMean);
		}		

		return res;
	}


	void multires_calc_statistisc(openvdb::FloatGrid::Ptr grid_reference, std::vector<openvdb::FloatGrid::Ptr> &multi_grids, int dim_size, std::string output_dir_path) {

		//openvdb::math::MinMax<float> min_max_ref = openvdb::tools::minMax(grid_reference->tree());
		//float min_value = min_max_ref.min(); // TODO pouzit i min value na norm
		//float max_value = min_max_ref.max();		

		//float min_value_mres = FLT_MAX; // TODO pouzit i min value na norm
		//float max_value_mres = -FLT_MAX;

		//for (size_t i = 0; i < multi_grids.size(); ++i) {
		//	openvdb::math::MinMax<float> min_max = openvdb::tools::minMax(multi_grids[i]->tree());
		//	if (min_value_mres != 0.0f)
		//		min_value_mres = std::min(min_value_mres, min_max.min());

		//	max_value_mres = std::max(max_value_mres, min_max.max());
		//}

		float min_value = FLT_MAX;
		float max_value = -FLT_MAX;
		float min_value_mres = FLT_MAX;
		float max_value_mres = -FLT_MAX;

		openvdb::FloatGrid::Accessor reference_grid_accessor = grid_reference->getAccessor();

		// compute min/max's
#ifdef WITH_OPENMP
# pragma omp parallel for reduction(min: min_value, min_value_mres) reduction(max: max_value, max_value_mres)
#endif  
		for (int z = 0; z < dim_size; ++z) {
			for (int y = 0; y < dim_size; ++y) {
				for (int x = 0; x < dim_size; ++x) {
					openvdb::Coord coord(x, y, z);
					float value0 = reference_grid_accessor.getValue(coord);

					min_value = fminf(min_value, value0);
					max_value = fmaxf(max_value, value0);

					float value1 = get_mres_value(grid_reference, multi_grids, coord);
					min_value_mres = fminf(min_value_mres, value1);
					max_value_mres = fmaxf(max_value_mres, value1);
				}
			}
		}

		printf("Min/Max value: %f/%f, mres: %f/%f\n", min_value, max_value, min_value_mres, max_value_mres);

		if (fabs(min_value - min_value_mres) > FDATA_EPSILON || fabs(max_value - max_value_mres) > FDATA_EPSILON) {
			printf("WARNING: Min/Max values are not equal, using only %f/%f.\n", min_value, max_value);
			min_value_mres = min_value;
			max_value_mres = max_value;			
		}

		double sum_squared = 0.0;
		double sum_squared_err = 0.0;
//		double sum_squared_err_zerro = 0.0;
		double num_voxels = 0.0;
		//int num_diff1 = 0;
		//int num_diff2 = 0;
		size_t num_zero_match1 = 0;
		size_t num_zero_match2 = 0;

		//openvdb::FloatGrid::Accessor reference_grid_accessor = grid_reference->getAccessor();
		openvdb::Coord inx_bbox_min = openvdb::Coord(0, 0, 0);
		openvdb::Coord inx_bbox_max = openvdb::Coord(inx_bbox_min[0] + dim_size - 1, inx_bbox_min[1] + dim_size - 1, inx_bbox_min[2] + dim_size - 1);

		openvdb::FloatGrid::Ptr vdb_grid = openvdb::FloatGrid::create(0.0f);
		vdb_grid->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);
		vdb_grid->setName("density");
		vdb_grid->setTransform(grid_reference->transformPtr());
		auto acc_vdb_grid = vdb_grid->getAccessor();

		//for (auto iter = grid_reference->cbeginValueOn(); iter; ++iter) { 
		for (int x = inx_bbox_min.x(); x <= inx_bbox_max.x(); x++) {
			for (int y = inx_bbox_min.y(); y <= inx_bbox_max.y(); y++) {
				for (int z = inx_bbox_min.z(); z <= inx_bbox_max.z(); z++) {				

					//openvdb::Coord coord = iter.getCoord();
					openvdb::Coord coord(x, y, z);

					//float value_reference = *iter;                
					double value_reference = reference_grid_accessor.getValue(coord);
					double value_reference_norm = (value_reference - min_value) / (max_value - min_value); // map to 0-1

					//openvdb::Vec3d world_pos = grid_reference->transform().indexToWorld(openvdb::Vec3d(coord.x(), coord.y(), coord.z()));

					//bool non_zero_found = false;
					//int num_none_zero = 0;

					double acc_value_sum = get_mres_value(grid_reference, multi_grids, coord);
					double acc_value_sum_norm = (acc_value_sum - min_value_mres) / (max_value_mres - min_value_mres);

					//printf("%d: (", num_voxels);
					//for (size_t i = 0; i < multi_grids.size(); ++i) {
					//	openvdb::FloatGrid::Ptr grid2 = multi_grids[i];
					//	openvdb::FloatGrid::Accessor acc2 = grid2->getAccessor();

					//	// sampler
					//	// openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::PointSampler> sampler(grid2->getConstAccessor(), grid2->transform());
					//	// float sampledValue = sampler.wsSample(world_pos);
					//	// sampledValue = sampledValue/max_value; // map to 0-1

					//	openvdb::Vec3d coord_w = grid_reference->transform().indexToWorld(coord);
					//	openvdb::Vec3d coord_grid2 = grid2->transform().worldToIndex(coord_w);
					//	openvdb::Coord coord_grid2_int = openvdb::Coord(int(coord_grid2.x()), int(coord_grid2.y()), int(coord_grid2.z()));
					//	acc_value_sum += acc2.getValue(coord_grid2_int);

					//}
					acc_vdb_grid.setValue(coord, acc_value_sum);

					double diff = value_reference_norm - acc_value_sum_norm;
					sum_squared_err += diff * diff;
					sum_squared += value_reference_norm * value_reference_norm;

					num_voxels = num_voxels + 1.0;

					if (value_reference_norm < FDATA_EPSILON)
						num_zero_match1++;
					if (acc_value_sum_norm < FDATA_EPSILON)
						num_zero_match2++;
				}
			}
		}
		printf("num_zero_match1: %lld\n", num_zero_match1);
		printf("num_zero_match2: %lld\n", num_zero_match2);

		std::string file_path = output_dir_path + "/grids_" + std::to_string(omp_get_wtime()) + "_merge.vdb";
		printf("Save multi res grids (merged): %s\n", file_path.c_str());
		openvdb::io::File file(file_path);

		file.write({ vdb_grid });
		file.close();

		//double MSE = sum_squared_err / (num_diff1 + num_diff2);
		double MSE = sum_squared_err / num_voxels;
		//double MSE_zero = sum_squared_err / (num_diff1 + num_diff2 + num_zero_match);
		double PSNR = 10.0 * std::log10(1.0 / MSE); // max value is 1.0 if MSE is normalized
		//double PSNR = 10.0 * std::log10(1.0 * max_value * max_value / MSE); // max value is 1.0 if MSE is normalized
		//double PSNR = 10.0 * std::log10((max_value * max_value)/ MSE); //

		double SNR = 20.0 * std::log10(sqrt(sum_squared / num_voxels) / sqrt(MSE));

		double SSIM = compute_ssim(grid_reference, multi_grids, dim_size, min_value, max_value, min_value_mres, max_value_mres);

		printf("PSNR: %e\n", PSNR);
		printf("SNR: %e\n", SNR);
		printf("SSIM: %e\n", SSIM);
		//printf("PSNR_zero: %f\n", PSNR_zero);
		printf("MSE: %e\n", MSE);
		//printf("MSE_voxels: %f\n", MSE_voxels);
		//printf("MSE_zero: %f\n", MSE_zero);
		printf("sum_squared_err: %f\n", sum_squared_err);
		printf("num_voxels: %lld\n", (size_t)num_voxels);
		//printf("num_diff1: %d\n", num_diff1);
		//printf("num_diff2: %d\n", num_diff2);
		//printf("diff: %d\n", num_diff1 + num_diff2);
		//printf("num_zero_match: %d\n", num_zero_match);
		//printf("diff + zero_match %d\n", num_diff1 + num_diff2 + num_zero_match);

		//Stats s = compute_stats(grid_reference, multi_grids, dim_size);
		//std::cout << "min/max (in) ....: [" << s.minValue << ',' << s.maxValue << "]\n";
		//std::cout << "min/max (out) ...: [" << s.minVDB << ',' << s.maxVDB << "]\n";
		//std::cout << "MSE .............: " << s.mse << '\n';
		//std::cout << "SNR .............: " << s.snr << '\n';
		//std::cout << "PSNR ............: " << s.psnr << '\n';
	}

	void create_multi_resolution_grid(common::vdb::ConvertVDBBase* convert_vdb_base,
		openvdb::FloatGrid::Ptr moto_grid, common::vdb::VDBParticles& grid_main_sum, 
		common::SpaceData& space_data,
		common::SpaceData::DenseType dense_type,
		std::string out_full_filepath){
		int dim_size = space_data.bbox_dim;

		// --- create array of points from vdb grid
		printf("Prepare input data\n");
		std::vector<std::vector<float>> points;
		std::vector<std::vector<float>> extended_points;
		std::vector<std::vector<float>> data_bbox;
		std::vector<float> global_min_precision; // max difference for neighbouring values (twice the resolution - because at least two sizes are needed to compute fractal dimension)

		float prec_factor = 0.001f;
		const char* env_prec_factor = std::getenv("PREC_FACTOR");
		if (env_prec_factor != nullptr) {
			prec_factor = std::stof(env_prec_factor);
			std::cout << "PREC_FACTOR: " << prec_factor << std::endl;
		}
		
		float value_precision = space_data.max_value_reduced * prec_factor; // maximal precision in which we are interested (eg 0.2)
		std::cout << "Value precision: " << value_precision << std::endl;

		float f_dim_threshold = 1.2; // 1.20
		const char* env_f_dim_threshold = std::getenv("F_DIM_THRESHOLD");
		if (env_f_dim_threshold != nullptr) {
			f_dim_threshold = std::stof(env_f_dim_threshold);
			std::cout << "F_DIM_THRESHOLD: " << f_dim_threshold << std::endl;
		}		

		prepare_input_multi_res(moto_grid, grid_main_sum, dense_type, points, extended_points, data_bbox, global_min_precision, value_precision, dim_size);
		std::vector<float> voxel_size({(data_bbox[1][0] - data_bbox[0][0])/(dim_size - 1), (data_bbox[1][1] - data_bbox[0][1])/(dim_size - 1), (data_bbox[1][2] - data_bbox[0][2])/(dim_size - 1)});

		// --- save input data
		// std::ofstream dataFile(out_full_filepath + "/out_multi_res/docs/input_values.txt");
		// std::ofstream dataFile_extended(out_full_filepath + "/out_multi_res/docs/input_values_extended.txt");
		// for (int i = 0; i < points.size(); i++) {
		// 	float x = (points[i][2]/voxel_size[2]) * dim_size * dim_size + (points[i][1]/voxel_size[1]) * dim_size  + (points[i][0]/voxel_size[0]); // TODO pokud nezacinaji body v nule, potreba prepsat
		// 	//float x = points[i][0] = i; //test
		// 	float y = points[i][3]; 
		// 	dataFile << x << " " << y << "\n";
		// 	dataFile_extended << x << " " << y << "\n";
		// }
		// dataFile.close();

		// for (int i = 0; i < extended_points.size(); i++) {
		// 	float x = (extended_points[i][2]/voxel_size[2]) * dim_size * dim_size + (extended_points[i][1]/voxel_size[1]) * dim_size  + (extended_points[i][0]/voxel_size[0]); // TODO pokud nezacinaji body v nule, potreba prepsat
		// 	//float x = points[i][0] = i; //test
		// 	float y = extended_points[i][3];
		// 	dataFile_extended << x << " " << y << "\n";
		// }
		// dataFile_extended.close();

		// --- move values above zero (openvdb takes zero as background)
		// find min value and max value (for normalization)
#if 0
		float min_value = points[0][3];
		float max_value = points[0][3];
		for (auto point: points){
			min_value = std::min(min_value, point[3]);
			max_value = std::max(max_value, point[3]);
		}
		std::cout << "Min value: " << min_value << std::endl;
		std::cout << "Max value: " << max_value << std::endl;
		
		// move values if min value is < 0
		if(min_value < 0){
			float min_target_value = 1.0;
			for (int i = 0; i < points.size(); i++){
				points[i][3] = min_target_value + points[i][3] - min_value;
			}
			for (int i = 0; i < extended_points.size(); i++){
				extended_points[i][3] = min_target_value + extended_points[i][3] - min_value;
			}
		}
#endif		
		// --- create multi resolution grid from input points
		printf("Divide input grid\n");
		int division_current = 0;
		int division_max = 10;
		int division_start = 1;		
		printf("f_dim_threshold: %f\n", f_dim_threshold);

		Openvdb_Multi_Res_Grids openvdb_grids(data_bbox, division_max);
#if 0
		openvdb_grids.global_max_value = max_value;
		printf("Global max: %f\n", openvdb_grids.global_max_value);
#endif
		std::vector<int> num_points_3D({dim_size, dim_size, dim_size});

		//std::string out_folder = // TODO nastavit vystupni slozku podle casu a pak tuhle slozku poslat dal

		divide_grid_multi_res(openvdb_grids, data_bbox, points, extended_points, division_current, division_max, division_start, f_dim_threshold, global_min_precision, out_full_filepath, voxel_size, num_points_3D);

		openvdb::FloatGrid::Ptr grid_reference = grid_main_sum.vdb_grid;
		
		// --- save grids to openvdb files
		std::string output_dir = out_full_filepath + "/out_multi_res/last";
		printf("Save multi res grids: %s\n", output_dir.c_str());
		openvdb_grids.save_vdb_grids(output_dir, division_start, grid_reference->transform().baseMap()->getAffineMap()->getMat4());

#if 0
		// --- print error values
		float norm_value = points.size() * (max_value - min_value);
		float NRMS = (std::sqrt(openvdb_grids.sum_diff_squared))/norm_value;
		float NAVG = openvdb_grids.sum_diff_abs/norm_value;
		float MSE = openvdb_grids.sum_diff_squared_norm/openvdb_grids.num_voxels_for_metrics;
		printf("MSE: %f\n", MSE);
		printf("sum_diff_squared_norm: %f\n", openvdb_grids.sum_diff_squared_norm);
		printf("voxels for metrics: %d\n", openvdb_grids.num_voxels_for_metrics);
		// printf("RMS: %f\n", NRMS);
		// printf("AVG: %f\n", NAVG);
#endif
		
#if 0
		if (dense_type != common::SpaceData::DenseType::eNone) { // dense grid
			//std::vector<float> input_data_values = grid_main_sum.dense_grid.data_density;
			grid_reference = convert_vdb_base->dense_to_openvdb(grid_main_sum.dense_grid, space_data.transform_scale);
		}
		else {
			grid_reference = grid_main_sum.vdb_grid;
		}
#endif

		multires_calc_statistisc(grid_reference, openvdb_grids.grids_array, dim_size, output_dir);

	}
#endif
	void finalize_grid(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main_sum, common::vdb::VDBParticles& grid_main_final)
	{
		if (from_cl.world_rank == 0 || space_data.anim_type != common::SpaceData::AnimType::eNone) {
			double t = omp_get_wtime();

			//Dense
			if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {

				grid_main_final.type = common::vdb::VDBParticles::VDBParticleType::eVector;

				if (from_cl.use_nanovdb) {
					auto nanovdb_handleI = convert_vdb_base->dense_to_nanovdb(grid_main_sum.dense_grid, space_data.transform_scale, space_data.dense_type, space_data.dense_norm);
					//nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*nanovdb_handleI);

#if OPENVDB_VERSION == 11
					nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::createNanoGrid(*nanovdb_handleI);
#else
					nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*nanovdb_handleI);
#endif					

					grid_main_final.vector_grid.resize(grid_handle_final.size());
					memcpy(grid_main_final.vector_grid.data(), grid_handle_final.data(), grid_handle_final.size());

					nanovdb::NanoGrid<float>* nanogrid = (nanovdb::NanoGrid<float>*) grid_handle_final.data();

#if OPENVDB_VERSION == 11					
					space_data.min_value_reduced = nanogrid->tree().root().data()->getMin();
					space_data.max_value_reduced = nanogrid->tree().root().data()->getMax();
#else
					nanogrid->tree().extrema(space_data.min_value_reduced, space_data.max_value_reduced);
#endif
				}
				else {

#ifdef WITH_OPENVDB
					//auto openvdb_handleE = convert_vdb_base->dense_to_openvdbE(grid_main_sum.dense_grid, space_data.transform_scale);
					auto openvdb_handleI = convert_vdb_base->dense_to_openvdb(grid_main_sum.dense_grid, space_data.transform_scale, space_data.dense_type, space_data.dense_norm);
					
					//convert_vdb_base->openvdb_to_vector2(openvdb_handleE, openvdb_handleI, grid_main_final.vector_grid);
					convert_vdb_base->openvdb_to_vector(openvdb_handleI, grid_main_final.vector_grid);
					std::cout << "Active voxels in input grid: " << openvdb_handleI->activeVoxelCount() << std::endl;

					//printf("rank: %d: minI: %e, maxI: %e\n", world_rank, min_value_reduced, max_value_reduced);

					//openvdb_handleI->evalMinMax(space_data.min_value_reduced, space_data.max_value_reduced);
					eval_min_max(openvdb_handleI, space_data.min_value_reduced, space_data.max_value_reduced);

#ifdef WITH_MULTIRES
					if (from_cl.use_multires) {
						// ----- MULTI-RES
						printf("START MULTI-RES\n");
						std::string out_full_filepath = from_cl.output_path;
// #if 0
// 						//////////////////////////////////////////////////////
// 						// Input and output file names
// 						const std::string inputFile = "f:\\temp\\multires_final3\\aneurism_256x256x256_float32.raw";

// 						// Open input file
// 						std::ifstream input(inputFile, std::ios::binary);
// 						if (!input) {
// 							throw std::runtime_error("Failed to open input file.");
// 						}

// 						constexpr size_t INPUT_DIM_X = 256;
// 						constexpr size_t INPUT_DIM_Y = 256;
// 						constexpr size_t INPUT_DIM_Z = 256;

// 						// Process data one value at a time
// 						float inputValue;
// 						float outputValue;

// 						common::vdb::DenseParticles dense_grid;
// 						dense_grid.create(INPUT_DIM_X, INPUT_DIM_Y, INPUT_DIM_Z);

// 						input.read(reinterpret_cast<char*>(dense_grid.data_density.data()), dense_grid.memsize());

// 						//for (size_t z = 0; z < INPUT_DIM_Z; ++z) {
// 						//	for (size_t y = 0; y < INPUT_DIM_Y; ++y) {
// 						//		for (size_t x = 0; x < INPUT_DIM_X; ++x) {
// 						//			// Calculate the corresponding input index
// 						//			size_t inputIndex = dense_grid.get_index(x, y, z); //x + y * (size_t)INPUT_DIM_X + z * (size_t)INPUT_DIM_X * INPUT_DIM_Y;

// 						//			// Seek to the appropriate position in the input file
// 						//			input.seekg(inputIndex * sizeof(float), std::ios::beg);
// 						//			if (!input.read(reinterpret_cast<char*>(&inputValue), sizeof(float))) {
// 						//				throw std::runtime_error("Failed to read value from input file.");
// 						//			}

// 						//			// Convert the value and write to the output file
// 						//			//outputValue = static_cast<float>(inputValue);
// 						//			//outputValue = inputValue;// == 0.f ? 0.f : log10(inputValue);
// 						//			dense_grid.data_density[inputIndex] = inputValue;
// 						//		}
// 						//	}
// 						//}

// 						// Close files
// 						input.close();

// 						std::cout << "Reading complete! Output written to " << inputFile << std::endl;

// 						openvdb_handleI = convert_vdb_base->dense_to_openvdb(dense_grid, 1.0f, common::SpaceData::DenseType::eNone);
// 						convert_vdb_base->openvdb_to_vector(openvdb_handleI, grid_main_final.vector_grid);
// 						eval_min_max(openvdb_handleI, space_data.min_value_reduced, space_data.max_value_reduced);
// 						//////////////////////////////////////////////////////
// #endif

						//create_multi_resolution_grid(convert_vdb_base, nullptr, grid_main_sum, space_data, out_full_filepath);

						common::vdb::VDBParticles grid_main_sum2;
						grid_main_sum2.type = common::vdb::VDBParticles::VDBParticleType::eOpenVDB;
						grid_main_sum2.vdb_grid = openvdb_handleI;

						common::SpaceData::DenseType dense_type = common::SpaceData::DenseType::eNone;

						create_multi_resolution_grid(convert_vdb_base, nullptr, grid_main_sum2, space_data, dense_type, out_full_filepath);

						printf("END MULTI-RES\n");
						// ----- END MULTI-RES
					}
#endif

#endif
				}

				printf("rank: %d: minI: %e, maxI: %e, reduced: minI: %e, maxI: %e\n", from_cl.world_rank, space_data.min_value, space_data.max_value, space_data.min_value_reduced, space_data.max_value_reduced);
				printf("rank: %d: final grid time: %f\n", from_cl.world_rank, omp_get_wtime() - t);

				save_raw_volume(convert_vdb_base, from_cl, space_data, grid_main_sum);

			}
			else {
				grid_main_final.type = common::vdb::VDBParticles::VDBParticleType::eVector;

				if (grid_main_sum.type == common::vdb::VDBParticles::VDBParticleType::eVector) {
					grid_main_final.vector_grid = std::move(grid_main_sum.vector_grid);

					if (from_cl.use_nanovdb) {
						nanovdb::NanoGrid<float>* nanogrid = (nanovdb::NanoGrid<float>*) grid_main_final.vector_grid.data();

#if OPENVDB_VERSION == 11					
						space_data.min_value_reduced = nanogrid->tree().root().data()->getMin();
						space_data.max_value_reduced = nanogrid->tree().root().data()->getMax();
#else
						nanogrid->tree().extrema(space_data.min_value_reduced, space_data.max_value_reduced);
#endif
					}
				}
				else if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
					grid_main_sum.raw_particles.serialize(grid_main_final.vector_grid);

					space_data.min_value_reduced = space_data.min_value;
					space_data.max_value_reduced = space_data.max_value;

					save_raw_particles_to_vdb(convert_vdb_base, from_cl, space_data, grid_main_sum);

					printf("rank: %d: Particles count: %lld\n", from_cl.world_rank, (size_t)(grid_main_sum.raw_particles.data[0].values.size() / 3));
				}
				else if (from_cl.use_nanovdb) {
#if OPENVDB_VERSION == 11
					nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::createNanoGrid(*grid_main_sum.nano_grid);
#else
					nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*grid_main_sum.nano_grid);
#endif

#if 0				//TODO
					openvdb::GridBase::Ptr openvdb_base_handle = nanovdb::nanoToOpenVDB(grid_handle_final);
					openvdb::FloatGrid::Ptr openvdb_handle = openvdb::gridPtrCast<openvdb::FloatGrid>(openvdb_base_handle);

					convert_vdb_base->openvdb_to_vector(openvdb_handle, grid_main_final.vector_grid);
					openvdb_handle->evalMinMax(space_data.min_value_reduced, space_data.max_value_reduced);
					eval_min_max(openvdb_handle, space_data.min_value_reduced, space_data.max_value_reduced);
#else
					grid_main_final.vector_grid.resize(grid_handle_final.size());
					memcpy(grid_main_final.vector_grid.data(), grid_handle_final.data(), grid_handle_final.size());

					nanovdb::NanoGrid<float>* nanogrid = (nanovdb::NanoGrid<float>*) grid_handle_final.data();

#if OPENVDB_VERSION == 11					
					space_data.min_value_reduced = nanogrid->tree().root().data()->getMin();
					space_data.max_value_reduced = nanogrid->tree().root().data()->getMax();
#else
					nanogrid->tree().extrema(space_data.min_value_reduced, space_data.max_value_reduced);
#endif

#endif
				}
				else {

#ifdef WITH_OPENVDB
					convert_vdb_base->openvdb_to_vector(grid_main_sum.vdb_grid, grid_main_final.vector_grid);
					//grid_main_sum.vdb_grid->evalMinMax(space_data.min_value_reduced, space_data.max_value_reduced);
					eval_min_max(grid_main_sum.vdb_grid, space_data.min_value_reduced, space_data.max_value_reduced);
					std::cout << "Active voxels in input grid: " << grid_main_sum.vdb_grid->activeVoxelCount() << std::endl;

#ifdef WITH_MULTIRES					
					if (from_cl.use_multires) {
						// ----- MULTI-RES
						printf("START MULTI-RES\n");
						std::string out_full_filepath = from_cl.output_path;
						float min_val, max_val;
						//grid_main_sum.vdb_grid->evalMinMax(min_val, max_val);
						eval_min_max(grid_main_sum.vdb_grid, min_val, max_val);
						printf("Input value range: %f, %f \n", min_val, max_val);
						create_multi_resolution_grid(convert_vdb_base, nullptr, grid_main_sum, space_data, space_data.dense_type, out_full_filepath);

						// // --test plane
						// openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);
						// openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
						// //int z = 50;
						// for (int x = 0; x < 100; ++x) {
						// 	for (int y = 0; y < 100; ++y) {
						// 		for (int z = 0; z < 100; ++z) {
						// 			accessor.setValue(openvdb::Coord(x, y, z), float(x+y+z));
						// 		}
						// 	}
						// }
						// create_multi_resolution_grid(grid, out_full_filepath);

						// // --test motorbike
						// openvdb::io::File file("/home/hra/IT4I/data/navy_multi_grid/input_data/navy_input_motorbike_bigger_clip_200.vdb");
						// file.open();
						// openvdb::GridBase::Ptr baseGrid;
						// baseGrid = file.readGrid("p");
						// file.close();
						// openvdb::FloatGrid::Ptr originalGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

						//  // Create a new grid to store the subregion [50,150]
						// openvdb::FloatGrid::Ptr subGrid = openvdb::FloatGrid::create(0.0f);

						// openvdb::math::Mat4d matrix = originalGrid->transform().baseMap()->getAffineMap()->getMat4();
						// openvdb::math::Transform::Ptr newTransform = openvdb::math::Transform::createLinearTransform(matrix);
						// subGrid->setTransform(newTransform);

						// openvdb::FloatGrid::Accessor subAccessor = subGrid->getAccessor();

						// // Copy only voxels within the range [50, 150]
						// int inx_start = 50;
						// int inx_end = 150;
						// for (int x = inx_start; x < inx_end; ++x) {
						// 	for (int y = inx_start; y < inx_end; ++y) {
						// 		for (int z = inx_start; z < inx_end; ++z) {
						// 			openvdb::Coord coord(x, y, z);

						// 			subAccessor.setValue(coord, originalGrid->tree().getValue(coord));

						// 		}
						// 	}
						// }

						// // Print number of active voxels in the subgrid
						// std::cout << "Active voxels in subgrid: " << subGrid->activeVoxelCount() << std::endl;

						// subGrid->evalMinMax(min_val, max_val);
						// eval_min_max(subGrid.vdb_grid, min_val, max_val);
						// printf("Input value range: %f, %f \n", min_val, max_val);

						// create_multi_resolution_grid(subGrid, nullptr, space_data.dense_type, space_data.bbox_dim, out_full_filepath);
						// // end test motorbike

						printf("END MULTI-RES\n");
						// ----- END MULTI-RES
					}
#endif

#endif

// #if 0 //print values
// 					// Iterate over the active voxels and copy values to the dense matrix
// 					for (openvdb::FloatGrid::ValueOnCIter iter = grid_main_sum.vdb_grid->cbeginValueOn(); iter; ++iter) {
// 						openvdb::Coord xyz = iter.getCoord();
// 						float v = *iter;

// 						printf("VDB Values: %d, %d, %d : %f\n", xyz.x(), xyz.y(), xyz.z(), v);
// 					}
// #endif
				}

				printf("rank: %d: minI: %e, maxI: %e, reduced: minI: %e, maxI: %e\n", from_cl.world_rank, space_data.min_value, space_data.max_value, space_data.min_value_reduced, space_data.max_value_reduced);
				printf("rank: %d: grid_handle_merged_size: %lld\n", from_cl.world_rank, grid_main_final.vector_grid.size());
				printf("rank: %d: final grid time: %f\n", from_cl.world_rank, omp_get_wtime() - t);
			}

			//printf("rank: %d: normalize_values: %f, voxels_count: %lld, voxels_count_zero: %lld\n", from_cl.world_rank, omp_get_wtime() - t, voxels_count, voxels_count_zero);

			if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {
				grid_main_sum.dense_grid.clear();
			}
		}
	}


	void save_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main_final, common::vdb::VDBParticles::VDBParticleType particle_type, bool only_rank0)
	{
		if (!only_rank0 || from_cl.world_rank == 0 || (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) &&
			(space_data.anim_type == common::SpaceData::AnimType::eFrameExtract && space_data.frame == space_data.anim_start + from_cl.world_rank || space_data.anim_type == common::SpaceData::AnimType::eAllPath)) {
			std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);

			if (!only_rank0 || space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
				char temp[1024];
				sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
				full_filepath = full_filepath + "_" + std::string(temp);
				//space_data.anim_task_counter++;
			}

			if (particle_type == common::vdb::VDBParticles::VDBParticleType::eRawParticles) {
				full_filepath = full_filepath + std::string(".part");
			}
			else if (particle_type == common::vdb::VDBParticles::VDBParticleType::eVector) {
				if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
					full_filepath = full_filepath + std::string(".bin");
				}
				else if (from_cl.use_nanovdb) {
					full_filepath = full_filepath + std::string(".nvdb");
				}
				else {
					full_filepath = full_filepath + std::string(".vdb");
				}				
			}
			else if (particle_type == common::vdb::VDBParticles::VDBParticleType::eNanoVDB) {
				full_filepath = full_filepath + std::string(".nvdb");
			}
			else {
				full_filepath = full_filepath + std::string(".vdb");
			}
			space_data.full_filepath = full_filepath;

			//	nanovdb::io::writeGrid(space_data.full_filepath, grid_handle_final);

			if (particle_type == common::vdb::VDBParticles::VDBParticleType::eOpenVDB) {
#ifdef WITH_OPENVDB
				openvdb::io::File(full_filepath).write({ grid_main_final.vdb_grid });
#endif

				printf("finished: %s\n", full_filepath.c_str());
			}
			else if (particle_type == common::vdb::VDBParticles::VDBParticleType::eVector) {

				// Open a file in binary mode
				std::ofstream output_file(full_filepath, std::ios::binary);
				if (!output_file) {
					printf("Unable to open file for writing: %s\n", full_filepath.c_str());
					return;
				}

				// Write the content of the vector to the file
				output_file.write((char*)grid_main_final.vector_grid.data(), grid_main_final.vector_grid.size());

				// Close the file
				output_file.close();

				printf("finished: %s\n", full_filepath.c_str());

			}
			else if (particle_type == common::vdb::VDBParticles::VDBParticleType::eNanoVDB) {
				//nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::createNanoGrid(*grid_main_final.nano_grid);

#if OPENVDB_VERSION == 11
				nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::createNanoGrid(*grid_main_final.nano_grid);
#else
				nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*grid_main_final.nano_grid);
#endif				

				nanovdb::io::writeGrid(space_data.full_filepath, grid_handle_final);
				printf("finished: %s\n", full_filepath.c_str());
			}
			else if (particle_type == common::vdb::VDBParticles::VDBParticleType::eRawParticles) {
				std::vector<uint8_t> grid_handle_main;
				grid_main_final.raw_particles.serialize(grid_handle_main);
				
				// Open a file in binary mode
				std::ofstream output_file(full_filepath, std::ios::binary);
				if (!output_file) {
					printf("Unable to open file for writing: %s\n", full_filepath.c_str());
					return;
				}

				// Write the content of the vector to the file
				output_file.write((char*)grid_handle_main.data(), grid_handle_main.size());

				// Close the file
				output_file.close();

				printf("finished: %s\n", full_filepath.c_str());
			}
			else if (particle_type == common::vdb::VDBParticles::VDBParticleType::eDense) {
				save_raw_volume(convert_vdb_base, from_cl, space_data, grid_main_final, only_rank0);
			}
			else {
				printf("Unknown Type for Saving\n");
			}
		}
	}

	void save_raw_volume(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, bool only_rank0)
	{
		if (from_cl.use_dense2file && grid_main.dense_grid.data_density.size() > 0) {
			if (!only_rank0 || from_cl.world_rank == 0) {
				std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);
				if (!only_rank0 || space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
					char temp[1024];
					sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
					full_filepath = full_filepath + "_" + std::string(temp);
					//space_data.anim_task_counter++;
				}

				full_filepath = full_filepath
					+ std::string("_") + std::to_string(grid_main.dense_grid.x())
					+ std::string("_") + std::to_string(grid_main.dense_grid.y())
					+ std::string("_") + std::to_string(grid_main.dense_grid.z())
					+ std::string("_float.raw");

				space_data.full_filepath = full_filepath;

				std::ofstream out(full_filepath.c_str(), std::ios::binary);
				out.write((const char*)grid_main.dense_grid.data_density.data(), grid_main.dense_grid.data_density.size() * sizeof(grid_main.dense_grid.data_density[0]));

				printf("finished: %s\n", full_filepath.c_str());
			}
		}
	}

	void save_raw_particles_to_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main)
	{
#ifdef WITH_OPENVDB
		if (from_cl.use_rawpart2vdb && grid_main.raw_particles.data.size() > 0) {
			if (from_cl.world_rank == 0) {
				std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);

				if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
					char temp[1024];
					sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
					full_filepath = full_filepath + "_" + std::string(temp);
					//space_data.anim_task_counter++;
				}

				full_filepath = full_filepath + std::string(".part.vdb");
				space_data.full_filepath = full_filepath;

				// positions
				std::vector<openvdb::Vec3f> points;
				for (size_t type = 0; type < grid_main.raw_particles.data.size(); type++) {
					const std::string attrib_name = grid_main.raw_particles.data[type].name;
					if (attrib_name == "position") {
						const size_t num_points = grid_main.raw_particles.data[type].values.size() / 3;
						points.resize(num_points);
#pragma omp parallel for
						for (size_t i = 0; i < num_points; i++) {
							points[i] = openvdb::Vec3f(grid_main.raw_particles.data[type].values[i * 3 + 0], grid_main.raw_particles.data[type].values[i * 3 + 1], grid_main.raw_particles.data[type].values[i * 3 + 2]);
						}
						break;
					}
				}

				// Compute voxel size and create the transform				
				openvdb::points::PointAttributeVector positions_wrapper(points);
				float vSize = openvdb::points::computeVoxelSize(positions_wrapper, 8);
				openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(vSize);

				// Create a PointIndexGrid. This can be done automatically on creation of the grid, however as this index grid is
				// required for the position and radius attributes, we create one we can use for both attribute creations.
				openvdb::tools::PointIndexGrid::Ptr pt_index = openvdb::tools::createPointIndexGrid<openvdb::tools::PointIndexGrid>(positions_wrapper, *transform);

				// Create a PointDataGrid containing these four points and using the point index grid. This requires the positions
				// wrapper.
				openvdb::points::PointDataGrid::Ptr pt_grid = openvdb::points::createPointDataGrid<openvdb::points::NullCodec, openvdb::points::PointDataGrid>(
					*pt_index, positions_wrapper, *transform);

				pt_grid->setName("points");

				// Attribute definition
				const openvdb::NamePair vdb_var = openvdb::points::TypedAttributeArray<float>::attributeType();
				const openvdb::NamePair vdb_var3 = openvdb::points::TypedAttributeArray<openvdb::Vec3f>::attributeType();

				for (size_t type = 0; type < grid_main.raw_particles.data.size(); type++) {
					const std::string attrib_name = grid_main.raw_particles.data[type].name;
					if (attrib_name == "position") {
						continue;
					}

					if (grid_main.raw_particles.data[type].num_comp == 1) {
						// Add the density attribute to the grid
						openvdb::points::appendAttribute(pt_grid->tree(), attrib_name, vdb_var);
						openvdb::points::populateAttribute(pt_grid->tree(), pt_index->tree(), attrib_name, openvdb::points::PointAttributeVector(grid_main.raw_particles.data[type].values));
					}
					else if (grid_main.raw_particles.data[type].num_comp == 3) {
						const size_t num_points = grid_main.raw_particles.data[type].values.size() / 3;
						std::vector<openvdb::Vec3f> values_vector(num_points);

#pragma omp parallel for
						for (size_t i = 0; i < num_points; i++) {
							values_vector[i] = openvdb::Vec3f(grid_main.raw_particles.data[type].values[i * 3 + 0], grid_main.raw_particles.data[type].values[i * 3 + 1], grid_main.raw_particles.data[type].values[i * 3 + 2]);
						}

						// Add the density attribute to the grid
						openvdb::points::appendAttribute(pt_grid->tree(), attrib_name, vdb_var3);
						openvdb::points::populateAttribute(pt_grid->tree(), pt_index->tree(), attrib_name, openvdb::points::PointAttributeVector(values_vector));
					}
				}

				openvdb::io::File file(full_filepath);
				file.setCompression(openvdb::io::COMPRESS_BLOSC | openvdb::io::COMPRESS_ACTIVE_MASK);
				file.write({ pt_grid });
				file.close();

				printf("finished: %s\n", full_filepath.c_str());
			}
		}
#endif
	}
}