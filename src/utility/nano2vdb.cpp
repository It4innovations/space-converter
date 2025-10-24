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

#ifdef WITH_NANOVDB

#ifdef WITH_TBB
#define NANOVDB_USE_TBB
#define NANOVDB_USE_INTRINSICS
#endif

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/GridBuilder.h>
#	include <nanovdb/util/IO.h>
#else
#	include <nanovdb/tools/GridBuilder.h>
#	include <nanovdb/io/IO.h>
#endif

#include <nanovdb/NanoVDB.h>
#endif

#if defined(WITH_OPENVDB) && defined(WITH_NANOVDB)

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/NanoToOpenVDB.h>
#else
#	include <nanovdb/tools/NanoToOpenVDB.h>
#endif

#	include <openvdb/openvdb.h>
#	include <openvdb/io/Stream.h>
#endif

int main(int argc, char** argv) {
	
	if (argc < 3) {
		printf("usage: nano2vdb file.nvdb file.vdb\n");
		exit(0);
	}

	std::string file_nvdb = argv[1];
	std::string file_vdb = argv[2];

#if defined(WITH_OPENVDB) && defined(WITH_NANOVDB)
	openvdb::initialize();  // Initialize OpenVDB if available.
	auto handle = nanovdb::io::readGrid<nanovdb::HostBuffer>(file_nvdb);

	auto openvdb_handle = nanovdb::tools::nanoToOpenVDB(handle);

	openvdb::io::File(file_vdb).write({ openvdb_handle });
	printf("finished: %s\n", file_vdb.c_str());
#endif

	return 0;
}