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

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>

//#ifdef WITH_NANOVDB

#ifdef WITH_TBB
#define NANOVDB_USE_TBB
#define NANOVDB_USE_INTRINSICS
#endif

#ifdef WITH_OPENVDB
#  include <openvdb/tools/Dense.h>
#endif
#ifdef WITH_NANOVDB
#  define NANOVDB_USE_OPENVDB
#  include <nanovdb/NanoVDB.h>
#  if NANOVDB_MAJOR_VERSION_NUMBER > 32 || \
      (NANOVDB_MAJOR_VERSION_NUMBER == 32 && NANOVDB_MINOR_VERSION_NUMBER >= 7)
#    include <nanovdb/tools/CreateNanoGrid.h>
#  else
#    include <nanovdb/util/OpenToNanoVDB.h>
#  endif

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/GridBuilder.h>
#	include <nanovdb/util/IO.h>
#else
#	include <nanovdb/tools/GridBuilder.h>
#	include <nanovdb/io/IO.h>
#endif

#endif

//#if OPENVDB_VERSION == 11
//#	include <nanovdb/util/GridBuilder.h>
//#	include <nanovdb/util/IO.h>
//#else
//#	include <nanovdb/tools/GridBuilder.h>
//#	include <nanovdb/io/IO.h>
//#endif
//
//#include <nanovdb/NanoVDB.h>
//#endif
//
//#if defined(WITH_OPENVDB) && defined(WITH_NANOVDB)
//
//#if OPENVDB_VERSION == 11
//#	include <nanovdb/util/OpenToNanoVDB.h>
//#else
//#	include <nanovdb/tools/OpenToNanoVDB.h>
//#endif
//
//#	include <openvdb/openvdb.h>
//#	include <openvdb/io/Stream.h>
//#endif

int main(int argc, char** argv) {
	
	if (argc != 3 && argc != 9 && argc != 6) {
		printf("usage: vdb2nano file.vdb -f 0.0 1.0 -g grid-name -o out.nvdb\n");
		printf("usage: vdb2nano file.vdb -g grid-name -o out.nvdb\n");
		printf("usage: vdb2nano file.vdb --print\n");		
		exit(0);
	}

	std::string out_file_nvdb = "";
	std::string grid_name = "";
	std::string vdb_file = "";
	float filter_min = 0.0f;
	float filter_max = 1.0f;
	bool use_filter = false;
	bool only_print = false;

	for (int i = 1; i < argc; i++) {
		const std::string arg = argv[i];
		if (arg == "-o") {
			out_file_nvdb = argv[++i];
		}
		else if (arg == "-g") {
			grid_name = argv[++i];
		}
		else if (arg == "-f") {
			use_filter = true;
			filter_min = atof(argv[++i]);
			filter_max = atof(argv[++i]);
		}
		else if (arg == "--print") {
			only_print = true;
		}
		else {
			vdb_file = arg;
		}
	}

#if defined(WITH_OPENVDB) && defined(WITH_NANOVDB)

	if (!only_print) {
		if (out_file_nvdb.empty()) {
			std::cerr << "Error: Output filename not specified. Use -o option.\n";
			return 1;
		}

		if (grid_name.empty()) {
			std::cerr << "Error: Grid name not specified. Use -g option.\n";
			return 1;
		}
	}

	if (vdb_file.empty()) {
		std::cerr << "Error: No input file specified.\n";
		return 1;
	}

	openvdb::initialize();  // Initialize OpenVDB if available.

	// Process each input file
	try {
		// Open the VDB file
		openvdb::io::File input_file(vdb_file);
		input_file.open();

		// Get all grids from the file
		openvdb::GridPtrVecPtr file_grids = input_file.getGrids();
		std::cout << "Processing file: " << vdb_file << " (" << file_grids->size() << " grids)\n";

		if (only_print) {
			std::cout << "Grids in file: " << vdb_file << "\n";
			for (auto& base_grid : *file_grids) {

				openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);

				if (grid) {
					std::cout << " - " << grid->getName() << "\n";
					openvdb::math::MinMax<float> extrema = openvdb::tools::minMax(grid->tree());
					std::cout << "   - min: " << extrema.min() << ", max: " << extrema.max() << "\n";
				}
				else {
					std::cout << " - " << grid->getName() << "\n";
				}
			}
		}
		else {
			for (auto& base_grid : *file_grids) {
				if (base_grid->getName() == grid_name) {

					if (use_filter) {
						std::cout << "Filtering values between " << filter_min << " and " << filter_max << std::endl;

#if OPENVDB_VERSION == 11
						nanovdb::build::FloatGrid nanogrid(0.0f, grid_name, nanovdb::GridClass::FogVolume);
#else
						nanovdb::tools::build::FloatGrid nanogrid(0.0f, grid_name, nanovdb::GridClass::FogVolume);
#endif
						auto acc_dst = nanogrid.getAccessor();

						openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);

						size_t orig_count = 0;
						size_t filter_count = 0;
						// Filtering
						for (auto iter = grid->cbeginValueOn(); iter; ++iter) {
							float value = *iter;
							if (value >= filter_min && value <= filter_max) {
								openvdb::Coord coord = iter.getCoord();
								acc_dst.setValue(nanovdb::Coord(coord.x(), coord.y(), coord.z()), value);
								filter_count++;
							}
							orig_count++;
						}

						std::cout << "After filtering: " << filter_count << "/" << orig_count << " active voxels (saved: " << int(100.0f - (float)filter_count * 100.0f / (float)orig_count) << "% from " << grid->activeVoxelCount() << ")" << std::endl;

#if OPENVDB_VERSION == 11
						nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle = nanovdb::createNanoGrid(nanogrid);
#else
						nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle = nanovdb::tools::createNanoGrid(nanogrid);
#endif

						nanovdb::io::writeGrid(out_file_nvdb, grid_handle);
						printf("finished: %s\n", out_file_nvdb.c_str());

						break;
					}
					else {

						//openvdb::FloatGrid::Ptr fgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
						//openvdb::CoordBBox bbox;
						//fgrid->tree().evalActiveVoxelBoundingBox(bbox);

						//std::vector<float> data_density(bbox.dim().x()* bbox.dim().y()* bbox.dim().z());

						//openvdb::tools::Dense<float, openvdb::tools::LayoutXYZ> denseIn(bbox, data_density.data());
						//openvdb::tools::copyToDense(fgrid->tree(), denseIn, false);

						//openvdb::FloatGrid::Ptr grid_new = openvdb::FloatGrid::create(0.0f);
						//grid_new->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);
						//grid_new->setName(fgrid->getName());
						//grid_new->setTransform(fgrid->transformPtr());

						//openvdb::tools::Dense<const float, openvdb::tools::LayoutXYZ> denseOut(bbox, data_density.data());
						//openvdb::tools::copyFromDense(denseOut, grid_new->tree(), 0.0f);

#if OPENVDB_VERSION == 11									

#  if NANOVDB_MAJOR_VERSION_NUMBER > 32 || \
      (NANOVDB_MAJOR_VERSION_NUMBER == 32 && NANOVDB_MINOR_VERSION_NUMBER >= 7)
						const openvdb::FloatGrid grid(*openvdb::gridConstPtrCast<openvdb::FloatGrid>(base_grid));
						nanovdb::GridHandle<nanovdb::HostBuffer> nanogrid = nanovdb::createNanoGrid<openvdb::FloatGrid, float>(grid);
#  else
						openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
						nanovdb::GridHandle<nanovdb::HostBuffer> nanogrid = nanovdb::openToNanoVDB(grid);
#  endif

#else
						const openvdb::FloatGrid grid(*openvdb::gridConstPtrCast<openvdb::FloatGrid>(base_grid));
						nanovdb::GridHandle<nanovdb::HostBuffer> nanogrid = nanovdb::tools::createNanoGrid<openvdb::FloatGrid, float>(grid);
#endif				

						nanovdb::io::writeGrid(out_file_nvdb, nanogrid);
						printf("finished: %s\n", out_file_nvdb.c_str());

						break;
					}
				}
			}
		}

		input_file.close();		
	}
	catch (const std::exception& e) {
		std::cerr << "Error processing file " << vdb_file << ": " << e.what() << std::endl;
	}
#endif

	return 0;
}