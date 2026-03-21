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

#include <cstdint>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

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

#	include <nanovdb/NanoVDB.h>
#endif

namespace {

enum class CombineMode {
	Mesh,
	Zip,
};

struct Args {
	std::string input_nvdb;
	std::string output_bin;
	CombineMode mode = CombineMode::Mesh;
	uint32_t grid_index = 0;
};

void print_usage() {
	printf("usage: nano2float3 file.nvdb file.bin [--mode mesh|zip] [--grid-index N]\n");
	printf("  mesh: dense export over indexBBox (includes inactive voxels/background values)\n");
	printf("  zip : export only active voxels\n");
}

bool parse_u32(const std::string& s, uint32_t& out) {
	char* end_ptr = nullptr;
	const unsigned long value = std::strtoul(s.c_str(), &end_ptr, 10);
	if (end_ptr == nullptr || *end_ptr != '\0') {
		return false;
	}
	if (value > static_cast<unsigned long>(UINT32_MAX)) {
		return false;
	}
	out = static_cast<uint32_t>(value);
	return true;
}

bool parse_args(int argc, char** argv, Args& out) {
	if (argc < 3) {
		print_usage();
		return false;
	}

	out.input_nvdb = argv[1];
	out.output_bin = argv[2];

	for (int i = 3; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "--mode") {
			if (i + 1 >= argc) {
				fprintf(stderr, "Error: --mode expects mesh or zip\n");
				return false;
			}
			const std::string mode = argv[++i];
			if (mode == "mesh") {
				out.mode = CombineMode::Mesh;
			} else if (mode == "zip") {
				out.mode = CombineMode::Zip;
			} else {
				fprintf(stderr, "Error: unknown mode '%s'\n", mode.c_str());
				return false;
			}
		} else if (arg == "--grid-index") {
			if (i + 1 >= argc) {
				fprintf(stderr, "Error: --grid-index expects an unsigned integer\n");
				return false;
			}
			if (!parse_u32(argv[++i], out.grid_index)) {
				fprintf(stderr, "Error: invalid --grid-index value\n");
				return false;
			}
		} else if (arg == "-h" || arg == "--help") {
			print_usage();
			exit(0);
		} else {
			fprintf(stderr, "Error: unknown argument '%s'\n", arg.c_str());
			print_usage();
			return false;
		}
	}

	return true;
}

uint32_t bswap32(uint32_t x) {
	return ((x & 0x000000FFu) << 24u) |
			   ((x & 0x0000FF00u) << 8u) |
			   ((x & 0x00FF0000u) >> 8u) |
			   ((x & 0xFF000000u) >> 24u);
}

bool host_is_little_endian() {
	const uint32_t value = 1u;
	return *reinterpret_cast<const uint8_t*>(&value) == 1u;
}

void write_float3_le(std::ofstream& out, float x, float y, float z) {
	float packed[3] = { x, y, z };
	if (!host_is_little_endian()) {
		for (int i = 0; i < 3; ++i) {
			uint32_t raw = 0u;
			std::memcpy(&raw, &packed[i], sizeof(uint32_t));
			raw = bswap32(raw);
			std::memcpy(&packed[i], &raw, sizeof(uint32_t));
		}
	}
	out.write(reinterpret_cast<const char*>(packed), sizeof(packed));
}

} // namespace

int main(int argc, char** argv) {
#ifdef WITH_NANOVDB
	Args args;
	if (!parse_args(argc, argv, args)) {
		return 1;
	}

	try {
		auto handle = nanovdb::io::readGrid<nanovdb::HostBuffer>(args.input_nvdb);

		if (handle.gridCount() == 0) {
			fprintf(stderr, "Error: no grids found in '%s'\n", args.input_nvdb.c_str());
			return 1;
		}
		if (args.grid_index >= handle.gridCount()) {
			fprintf(stderr, "Error: --grid-index %u out of range [0, %u)\n", args.grid_index, handle.gridCount());
			return 1;
		}

		const auto* meta = handle.gridMetaData(args.grid_index);
		if (!meta) {
			fprintf(stderr, "Error: failed to access grid metadata\n");
			return 1;
		}

		const auto* grid = handle.grid<nanovdb::Vec3f>(args.grid_index);
		if (!grid) {
			fprintf(stderr,
				"Error: selected grid #%u ('%s') is not a Vec3f grid\n",
				args.grid_index,
				meta->shortGridName());
			return 1;
		}

		if (meta->isEmpty()) {
			fprintf(stderr, "Error: selected grid is empty\n");
			return 1;
		}

		std::ofstream out(args.output_bin, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!out.is_open()) {
			fprintf(stderr, "Error: cannot open output file '%s'\n", args.output_bin.c_str());
			return 1;
		}

		const auto& bbox = meta->indexBBox();
		const nanovdb::Coord min = bbox.min();
		const nanovdb::Coord max = bbox.max();

		const auto& tree = grid->tree();
		uint64_t entries = 0u;

		for (int z = min[2]; z <= max[2]; ++z) {
			for (int y = min[1]; y <= max[1]; ++y) {
				for (int x = min[0]; x <= max[0]; ++x) {
					const nanovdb::Coord ijk(x, y, z);
					if (args.mode == CombineMode::Zip && !tree.isActive(ijk)) {
						continue;
					}
					const auto v = tree.getValue(ijk);
					write_float3_le(out, v[0], v[1], v[2]);
					++entries;
				}
			}
		}

		if (!out.good()) {
			fprintf(stderr, "Error: failed while writing '%s'\n", args.output_bin.c_str());
			return 1;
		}

		const uint64_t byte_size = entries * 3u * sizeof(float);
		printf("Grid: #%u (%s)\n", args.grid_index, meta->shortGridName());
		printf("Mode: %s\n", args.mode == CombineMode::Mesh ? "mesh" : "zip");
		printf("Index bbox: min=(%d,%d,%d) max=(%d,%d,%d)\n", min[0], min[1], min[2], max[0], max[1], max[2]);
		printf("Wrote %llu float3 entries to: %s\n", static_cast<unsigned long long>(entries), args.output_bin.c_str());
		printf("Byte size: %llu\n", static_cast<unsigned long long>(byte_size));
		return 0;
	} catch (const std::exception& exc) {
		fprintf(stderr, "Error: %s\n", exc.what());
		return 1;
	}
#else
	(void)argc;
	(void)argv;
	fprintf(stderr, "Error: nano2float3 was built without WITH_NANOVDB\n");
	return 1;
#endif
}