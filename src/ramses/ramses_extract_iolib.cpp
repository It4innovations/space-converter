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

#include "ramses_extract_iolib.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "convert_common.h"
#include "reader_return_macros.h"

// Reader data contract (see docs/SpaceConverter_Code_Analysis_2026-08.md §7):
//   positions             -> code length units, [0, boxlen] (cell centres for the gas)
//   get_particle_mass(id) -> gas: density x cell volume; particles: mass (code units)
//   get_particle_rho(id)  -> gas: density (code units); particles: 0
//   get_particle_hsml(id) -> gas: cell size; particles: size of the cells of the AMR
//                            level the particle is attached to (levelp), i.e. an
//                            adaptive smoothing length that follows the refinement
//   particle id space     -> per rank, [gas leaf cells][DM][Star][Cloud][Debris][Other],
//                            contiguous per type in type order (the base class relies on it)
//   MPI split             -> the ncpu RAMSES files are split into contiguous ranges
//                            over the ranks; every rank reads its own files only
//
// File formats (RAMSES >= 2017, "# version: 1" field descriptors): Fortran
// unformatted sequential records with 4-byte length markers, native endianness.

namespace ramses {
	namespace io {

		// ---------------------------------------------------------------------
		// Fortran unformatted sequential file
		// ---------------------------------------------------------------------
		class FortranFile {
		public:
			// The whole file is read at once: skipping the many small records of
			// the other CPUs' grids with seekg() would refill the stream buffer
			// from the (Lustre) file system on every record.
			explicit FortranFile(const std::string& path) : path_(path) {
				std::ifstream f(path, std::ios::binary | std::ios::ate);
				if (!f) {
					throw std::runtime_error("RAMSES: cannot open " + path);
				}
				data_.resize((size_t)f.tellg());
				f.seekg(0);
				if (!data_.empty() && !f.read(data_.data(), (std::streamsize)data_.size())) {
					throw std::runtime_error("RAMSES: cannot read " + path);
				}
			}

			// Read one record into buf.
			void read_record(std::vector<char>& buf) {
				size_t n = record_begin();
				buf.assign(data_.begin() + pos_, data_.begin() + pos_ + n);
				record_end(n);
			}

			void skip_record() {
				record_end(record_begin());
			}

			void skip_records(int n) {
				for (int i = 0; i < n; i++)
					skip_record();
			}

			template <typename T>
			std::vector<T> read_array() {
				size_t n = record_begin();
				std::vector<T> out(n / sizeof(T));
				if (!out.empty())
					std::memcpy(out.data(), &data_[pos_], out.size() * sizeof(T));
				record_end(n);
				return out;
			}

			template <typename T>
			T read_scalar() {
				std::vector<T> v = read_array<T>();
				if (v.empty()) {
					throw std::runtime_error("RAMSES: empty record in " + path_);
				}
				return v[0];
			}

			std::string read_string() {
				std::vector<char> buf;
				read_record(buf);
				std::string s(buf.begin(), buf.end());
				s.erase(s.find_last_not_of(" \0", std::string::npos, 2) + 1);
				return s;
			}

		private:
			int32_t read_marker() {
				if (pos_ + 4 > data_.size()) {
					throw std::runtime_error("RAMSES: unexpected end of " + path_);
				}
				int32_t n;
				std::memcpy(&n, &data_[pos_], 4);
				pos_ += 4;
				if (n < 0) {
					// gfortran splits records > 2 GB into sub-records (negative
					// markers); per-CPU RAMSES records never get that large
					throw std::runtime_error("RAMSES: sub-records (records > 2 GB) are not supported: " + path_);
				}
				return n;
			}

			// Opening marker; returns the payload size, pos_ at the payload
			size_t record_begin() {
				size_t n = (size_t)read_marker();
				if (pos_ + n + 4 > data_.size()) {
					throw std::runtime_error("RAMSES: corrupt record in " + path_);
				}
				return n;
			}

			// Skip the payload and check the closing marker
			void record_end(size_t n) {
				pos_ += n;
				if ((size_t)read_marker() != n) {
					throw std::runtime_error("RAMSES: corrupt record in " + path_);
				}
			}

			std::string path_;
			std::vector<char> data_;
			size_t pos_ = 0;
		};

		// One field of a *_file_descriptor.txt: name and Fortran kind character
		struct FieldDesc {
			std::string name;
			char kind = 'd';
		};

		// A block beyond BTMax: one hydro variable and/or one particle field
		// of the same name (e.g. "metallicity" exists in both)
		struct ExtraBlock {
			std::string name;
			int gas_var = -1;   // index into gas_vars, -1 = not a gas block
			int part_var = -1;  // index into part_extra, -1 = not a particle block
		};

		// ---------------------------------------------------------------------
		// Reader state
		// ---------------------------------------------------------------------
		std::string output_dir_path;
		int iout = 0;
		int ncpu = 0;
		int ndim = 3;
		int nlevelmax = 0;
		int levelmin = 0;
		double boxlen = 1.0;
		double time_code = 0.0;
		double unit_l = 1.0, unit_d = 1.0, unit_t = 1.0;

		std::vector<FieldDesc> hydro_desc;
		std::vector<FieldDesc> part_desc;
		std::vector<ExtraBlock> extra_blocks;
		int hydro_idx_density = -1;
		int hydro_idx_vel[3] = { -1, -1, -1 };

		// Gas leaf cells
		size_t ngas = 0;
		std::vector<double> gas_pos;                 // 3 per cell
		std::vector<float> gas_dx;                   // cell size (code length)
		std::vector<uint8_t> gas_level;
		std::vector<std::vector<float>> gas_vars;    // [hydro var][cell]

		// Particles (ordered by type after init)
		size_t npart_local = 0;
		std::vector<double> part_pos;                // 3 per particle
		std::vector<float> part_vel;                 // 3 per particle
		std::vector<double> part_mass;
		std::vector<uint8_t> part_level;
		std::vector<std::vector<float>> part_extra;  // [part_extra_desc][particle]
		std::vector<std::string> part_extra_names;

		// id ranges of the types: [type_offset[t], type_offset[t + 1])
		uint64_t type_offset[PTMax + 1] = { 0 };
		size_t global_num = 0;

		// Read-time filters resolved to field indices: (index, min, max)
		struct ResolvedFilter {
			int index;
			double min;
			double max;
		};
		std::vector<ResolvedFilter> gas_filters;   // index into gas_vars
		std::vector<ResolvedFilter> part_filters;  // index into the extra particle fields, -1 = mass

		double steps_time[2];

		// ---------------------------------------------------------------------
		// Helpers
		// ---------------------------------------------------------------------
		static std::string cpu_file(const std::string& kind, int icpu) {
			char name[64];
			snprintf(name, sizeof(name), "/%s_%05d.out%05d", kind.c_str(), iout, icpu);
			return output_dir_path + name;
		}

		static std::vector<FieldDesc> read_descriptor(const std::string& path) {
			std::vector<FieldDesc> desc;
			std::ifstream f(path);
			if (!f)
				return desc;
			std::string line;
			while (std::getline(f, line)) {
				if (line.empty() || line[0] == '#')
					continue;
				// "ivar, name, kind"
				std::stringstream ss(line);
				std::string ivar, name, kind;
				std::getline(ss, ivar, ',');
				std::getline(ss, name, ',');
				std::getline(ss, kind, ',');
				auto trim = [](std::string s) {
					s.erase(0, s.find_first_not_of(" \t"));
					s.erase(s.find_last_not_of(" \t\r\n") + 1);
					return s;
				};
				FieldDesc d;
				d.name = trim(name);
				kind = trim(kind);
				d.kind = kind.empty() ? 'd' : kind[0];
				if (!d.name.empty())
					desc.push_back(d);
			}
			return desc;
		}

		static void read_info(const std::string& path) {
			std::ifstream f(path);
			if (!f) {
				throw std::runtime_error("RAMSES: cannot open " + path);
			}
			std::string line;
			while (std::getline(f, line)) {
				size_t eq = line.find('=');
				if (eq == std::string::npos)
					continue;
				std::string key = line.substr(0, eq);
				key.erase(key.find_last_not_of(" \t") + 1);
				std::string val = line.substr(eq + 1);
				if (key == "ncpu") ncpu = std::stoi(val);
				else if (key == "ndim") ndim = std::stoi(val);
				else if (key == "levelmin") levelmin = std::stoi(val);
				else if (key == "levelmax") nlevelmax = std::stoi(val);
				else if (key == "boxlen") boxlen = std::stod(val);
				else if (key == "time") time_code = std::stod(val);
				else if (key == "unit_l") unit_l = std::stod(val);
				else if (key == "unit_d") unit_d = std::stod(val);
				else if (key == "unit_t") unit_t = std::stod(val);
			}
		}

		// Record payload -> float values of one field kind
		static void convert_field(const std::vector<char>& buf, char kind, std::vector<float>& out) {
			size_t n = 0;
			switch (kind) {
			case 'd': n = buf.size() / 8; out.resize(n); for (size_t i = 0; i < n; i++) { double v; std::memcpy(&v, &buf[i * 8], 8); out[i] = (float)v; } break;
			case 'f': n = buf.size() / 4; out.resize(n); std::memcpy(out.data(), buf.data(), n * 4); break;
			case 'q': n = buf.size() / 8; out.resize(n); for (size_t i = 0; i < n; i++) { int64_t v; std::memcpy(&v, &buf[i * 8], 8); out[i] = (float)v; } break;
			case 'i': n = buf.size() / 4; out.resize(n); for (size_t i = 0; i < n; i++) { int32_t v; std::memcpy(&v, &buf[i * 4], 4); out[i] = (float)v; } break;
			case 'h': n = buf.size() / 2; out.resize(n); for (size_t i = 0; i < n; i++) { int16_t v; std::memcpy(&v, &buf[i * 2], 2); out[i] = (float)v; } break;
			case 'b': n = buf.size(); out.resize(n); for (size_t i = 0; i < n; i++) out[i] = (float)(int8_t)buf[i]; break;
			default: n = buf.size() / 4; out.resize(n); for (size_t i = 0; i < n; i++) { int32_t v; std::memcpy(&v, &buf[i * 4], 4); out[i] = (float)(v != 0); } break;
			}
		}

		static void partition_range(size_t total, size_t parts, size_t idx, size_t& start, size_t& count) {
			size_t base = total / parts;
			size_t rem = total % parts;
			count = base + (idx < rem ? 1 : 0);
			start = idx * base + std::min(idx, rem);
		}

		// ---------------------------------------------------------------------
		// AMR + hydro of one CPU file: append the leaf cells
		// ---------------------------------------------------------------------
		static void read_cpu_gas(int icpu, int level_max) {
			FortranFile amr(cpu_file("amr", icpu));
			FortranFile hyd(cpu_file("hydro", icpu));

			// --- AMR header
			amr.skip_record();                                   // ncpu
			amr.skip_record();                                   // ndim
			std::vector<int32_t> nxyz = amr.read_array<int32_t>();
			int nlev = amr.read_scalar<int32_t>();
			amr.skip_record();                                   // ngridmax
			int nboundary = amr.read_scalar<int32_t>();
			amr.skip_record();                                   // ngrid_current
			amr.skip_record();                                   // boxlen
			amr.skip_records(11);                                // time variables .. mass_sph
			amr.skip_records(2);                                 // headl, taill
			std::vector<int32_t> numbl = amr.read_array<int32_t>();
			amr.skip_record();                                   // numbtot
			std::vector<int32_t> numbb;
			if (nboundary > 0) {
				amr.skip_records(2);                             // headb, tailb
				numbb = amr.read_array<int32_t>();
			}
			amr.skip_record();                                   // headf, tailf, ...
			std::string ordering = amr.read_string();
			amr.skip_records(ordering.rfind("bisection", 0) == 0 ? 5 : 1);
			amr.skip_records(3);                                 // coarse son, flag1, cpu_map

			// --- hydro header
			hyd.skip_record();                                   // ncpu
			int nvar = hyd.read_scalar<int32_t>();
			hyd.skip_records(4);                                 // ndim, nlevelmax, nboundary, gamma
			if (nvar != (int)hydro_desc.size()) {
				throw std::runtime_error("RAMSES: hydro file has " + std::to_string(nvar) +
					" variables, hydro_file_descriptor.txt " + std::to_string(hydro_desc.size()));
			}

			const int twotondim = 1 << ndim;
			const int nskip_grid = 3 + ndim + 1 + 2 * ndim + 3 * twotondim;
			// Grid centres xg are in coarse-cell units, with the (unit) domain
			// shifted by the boundary coarse cells: nx = 3 for non-periodic
			// boundaries, domain [1, 2]; nx = 1 when periodic, domain [0, 1]
			// (same convention as utils/f90/amr2cube.f90)
			double xbound[3];
			for (int d = 0; d < 3; d++)
				xbound[d] = (double)((nxyz.size() > (size_t)d ? nxyz[d] : 1) / 2);

			std::vector<char> rec;
			std::vector<std::vector<double>> xg(ndim);
			std::vector<std::vector<int32_t>> son(twotondim);
			std::vector<std::vector<double>> vals(twotondim * nvar);

			for (int ilevel = 1; ilevel <= nlev; ilevel++) {
				const double dx = std::ldexp(1.0, -ilevel);
				for (int ibound = 1; ibound <= ncpu + nboundary; ibound++) {
					int ncache = (ibound <= ncpu)
						? numbl[(size_t)(ibound - 1) + (size_t)ncpu * (ilevel - 1)]
						: numbb[(size_t)(ibound - ncpu - 1) + (size_t)nboundary * (ilevel - 1)];

					hyd.skip_record();                           // ilevel
					int ncache_h = hyd.read_scalar<int32_t>();
					if (ncache_h != ncache) {
						throw std::runtime_error("RAMSES: amr/hydro grid counts differ in cpu " + std::to_string(icpu));
					}
					if (ncache <= 0)
						continue;

					const bool own = (ibound == icpu) && (level_max <= 0 || ilevel <= level_max);
					if (!own) {
						amr.skip_records(nskip_grid);
						hyd.skip_records(twotondim * nvar);
						continue;
					}

					amr.skip_records(3);                         // ind_grid, next, prev
					for (int d = 0; d < ndim; d++)
						xg[d] = amr.read_array<double>();
					amr.skip_records(1 + 2 * ndim);              // father, nbor
					for (int c = 0; c < twotondim; c++)
						son[c] = amr.read_array<int32_t>();
					amr.skip_records(2 * twotondim);             // cpu_map, flag1

					for (int c = 0; c < twotondim; c++)
						for (int v = 0; v < nvar; v++)
							vals[(size_t)c * nvar + v] = hyd.read_array<double>();

					const bool at_cut = (level_max > 0 && ilevel == level_max);
					for (int c = 0; c < twotondim; c++) {
						for (int g = 0; g < ncache; g++) {
							if (son[c][g] != 0 && !at_cut)
								continue;                        // refined further: not a leaf
							bool keep = true;
							for (const ResolvedFilter& f : gas_filters) {
								double v = vals[(size_t)c * nvar + f.index][g];
								keep = keep && v >= f.min && v <= f.max;
							}
							if (!keep)
								continue;
							for (int d = 0; d < 3; d++) {
								double off = (d < ndim) ? ((double)((c >> d) & 1) - 0.5) * dx : 0.0;
								double x = (d < ndim) ? xg[d][g] + off - xbound[d] : 0.5;
								gas_pos.push_back(x * boxlen);
							}
							gas_dx.push_back((float)(dx * boxlen));
							gas_level.push_back((uint8_t)ilevel);
							for (int v = 0; v < nvar; v++)
								gas_vars[v].push_back((float)vals[(size_t)c * nvar + v][g]);
						}
					}
				}
			}
		}

		// ---------------------------------------------------------------------
		// Particles of one CPU file, appended to per-type temporary buffers
		// ---------------------------------------------------------------------
		struct PartBuffers {
			std::vector<double> pos;
			std::vector<float> vel;
			std::vector<double> mass;
			std::vector<uint8_t> level;
			std::vector<std::vector<float>> extra;
		};

		static int family_to_type(int family) {
			switch (family) {
			case 1: return DM;
			case 2: return Star;
			case 3: return Cloud;
			case 4: return Debris;
			default: return (family <= 0) ? -1 : Other;  // tracers (<= 0) are skipped
			}
		}

		static void read_cpu_particles(int icpu, std::vector<PartBuffers>& buffers) {
			FortranFile f(cpu_file("part", icpu));
			f.skip_record();                                     // ncpu
			f.skip_record();                                     // ndim
			int npart = f.read_scalar<int32_t>();
			f.skip_records(5);                                   // localseed, nstar_tot, mstar_tot, mstar_lost, nsink
			if (npart <= 0)
				return;

			std::vector<char> rec;
			std::vector<std::vector<double>> pos(3), vel(3);
			std::vector<double> mass;
			std::vector<float> level, family;
			std::vector<std::vector<float>> extra(part_extra_names.size());

			size_t iextra = 0;
			for (const FieldDesc& d : part_desc) {
				const std::string& n = d.name;
				if (n.rfind("position_", 0) == 0 || n.rfind("velocity_", 0) == 0 || n == "mass") {
					f.read_record(rec);
					std::vector<double> v(rec.size() / 8);
					if (d.kind == 'd') {
						std::memcpy(v.data(), rec.data(), v.size() * 8);
					}
					else {
						std::vector<float> tmp;
						convert_field(rec, d.kind, tmp);
						v.assign(tmp.begin(), tmp.end());
					}
					if (n == "mass") mass = std::move(v);
					else {
						int axis = n.back() - 'x';
						if (axis >= 0 && axis < 3)
							(n[0] == 'p' ? pos : vel)[axis] = std::move(v);
					}
				}
				else if (n == "levelp") {
					f.read_record(rec);
					convert_field(rec, d.kind, level);
				}
				else {
					f.read_record(rec);
					if (n == "family")
						convert_field(rec, d.kind, family);
					convert_field(rec, d.kind, extra[iextra]);
					iextra++;
				}
			}

			for (int i = 0; i < npart; i++) {
				int t = family.empty() ? (int)DM : family_to_type((int)family[i]);
				if (t < 0)
					continue;
				bool keep = true;
				for (const ResolvedFilter& f : part_filters) {
					double v;
					if (f.index < 0)
						v = mass.empty() ? 0.0 : mass[i];
					else
						v = extra[f.index].empty() ? 0.0 : (double)extra[f.index][i];
					keep = keep && v >= f.min && v <= f.max;
				}
				if (!keep)
					continue;
				PartBuffers& b = buffers[t];
				for (int a = 0; a < 3; a++) {
					b.pos.push_back(pos[a].empty() ? 0.0 : pos[a][i]);
					b.vel.push_back(vel[a].empty() ? 0.0f : (float)vel[a][i]);
				}
				b.mass.push_back(mass.empty() ? 0.0 : mass[i]);
				b.level.push_back(level.empty() ? (uint8_t)levelmin : (uint8_t)level[i]);
				b.extra.resize(extra.size());
				for (size_t e = 0; e < extra.size(); e++)
					b.extra[e].push_back(extra[e].empty() ? 0.0f : extra[e][i]);
			}
		}

		// ---------------------------------------------------------------------
		// Public API
		// ---------------------------------------------------------------------
		void init_lib(const std::string& output_dir, int world_rank, int world_size,
			bool read_gas, bool read_particles, int level_max,
			const std::vector<FieldFilter>& filters) {
#ifdef WITH_OPENMP
			steps_time[0] = omp_get_wtime();
#endif
			finish_lib();

			output_dir_path = output_dir;
			while (output_dir_path.size() > 1 && output_dir_path.back() == '/')
				output_dir_path.pop_back();

			std::smatch m;
			std::string base = output_dir_path.substr(output_dir_path.find_last_of('/') + 1);
			if (!std::regex_search(base, m, std::regex("output_([0-9]+)"))) {
				throw std::runtime_error("RAMSES: expected an output_NNNNN directory, got " + output_dir);
			}
			iout = std::stoi(m[1].str());

			char info_name[64];
			snprintf(info_name, sizeof(info_name), "/info_%05d.txt", iout);
			read_info(output_dir_path + info_name);

			hydro_desc = read_descriptor(output_dir_path + "/hydro_file_descriptor.txt");
			part_desc = read_descriptor(output_dir_path + "/part_file_descriptor.txt");
			read_gas = read_gas && !hydro_desc.empty();
			read_particles = read_particles && !part_desc.empty();

			// Hydro variables with a fixed block of their own
			for (size_t v = 0; v < hydro_desc.size(); v++) {
				const std::string& n = hydro_desc[v].name;
				if (n == "density") hydro_idx_density = (int)v;
				if (n == "velocity_x") hydro_idx_vel[0] = (int)v;
				if (n == "velocity_y") hydro_idx_vel[1] = (int)v;
				if (n == "velocity_z") hydro_idx_vel[2] = (int)v;
			}
			if (read_gas && hydro_idx_density < 0) {
				throw std::runtime_error("RAMSES: no 'density' in hydro_file_descriptor.txt (conservative output is not supported)");
			}

			// Extra blocks: every hydro variable, then the particle fields not
			// covered by the fixed blocks, merged by name
			for (size_t v = 0; v < hydro_desc.size(); v++) {
				ExtraBlock b;
				b.name = hydro_desc[v].name;
				b.gas_var = (int)v;
				extra_blocks.push_back(b);
			}
			for (const FieldDesc& d : part_desc) {
				const std::string& n = d.name;
				if (n.rfind("position_", 0) == 0 || n.rfind("velocity_", 0) == 0 || n == "mass" || n == "levelp")
					continue;
				int pidx = (int)part_extra_names.size();
				part_extra_names.push_back(n);
				auto it = std::find_if(extra_blocks.begin(), extra_blocks.end(),
					[&](const ExtraBlock& b) { return b.name == n; });
				if (it != extra_blocks.end()) {
					it->part_var = pidx;
				}
				else {
					ExtraBlock b;
					b.name = n;
					b.part_var = pidx;
					extra_blocks.push_back(b);
				}
			}

			// Filters -> field indices
			for (const FieldFilter& f : filters) {
				bool used = false;
				for (size_t v = 0; v < hydro_desc.size(); v++) {
					if (hydro_desc[v].name == f.name) {
						gas_filters.push_back({ (int)v, f.min, f.max });
						used = true;
					}
				}
				if (f.name == "mass") {
					part_filters.push_back({ -1, f.min, f.max });
					used = true;
				}
				for (size_t e = 0; e < part_extra_names.size(); e++) {
					if (part_extra_names[e] == f.name) {
						part_filters.push_back({ (int)e, f.min, f.max });
						used = true;
					}
				}
				if (!used && world_rank == 0) {
					printf("Warning: --ramses-filter field '%s' is neither a hydro variable nor a particle field\n", f.name.c_str());
				}
			}

			// This rank's CPU files
			size_t first = 0, count = 0;
			partition_range((size_t)ncpu, (size_t)world_size, (size_t)world_rank, first, count);

			gas_vars.assign(hydro_desc.size(), std::vector<float>());
			std::vector<PartBuffers> buffers(PTMax);
			for (size_t k = 0; k < count; k++) {
				int icpu = (int)(first + k) + 1;
				if (read_gas)
					read_cpu_gas(icpu, level_max);
				if (read_particles)
					read_cpu_particles(icpu, buffers);
			}
			ngas = gas_dx.size();

			// Particles in type order
			type_offset[0] = 0;
			type_offset[1] = ngas;
			part_extra.assign(part_extra_names.size(), std::vector<float>());
			for (int t = 1; t < PTMax; t++) {
				PartBuffers& b = buffers[t];
				part_pos.insert(part_pos.end(), b.pos.begin(), b.pos.end());
				part_vel.insert(part_vel.end(), b.vel.begin(), b.vel.end());
				part_mass.insert(part_mass.end(), b.mass.begin(), b.mass.end());
				part_level.insert(part_level.end(), b.level.begin(), b.level.end());
				for (size_t e = 0; e < part_extra.size() && e < b.extra.size(); e++)
					part_extra[e].insert(part_extra[e].end(), b.extra[e].begin(), b.extra[e].end());
				type_offset[t + 1] = type_offset[t] + b.mass.size();
				b = PartBuffers();  // release the buffer of this type
			}
			npart_local = part_mass.size();

			size_t local = ngas + npart_local;
			global_num = local;
#ifdef WITH_MPI
			unsigned long long l = local, g = 0;
			MPI_Allreduce(&l, &g, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
			global_num = g;
#endif

			if (world_rank == 0) {
				printf("RAMSES output %05d: ncpu %d, levels %d..%d, boxlen %g, time %g (unit_t %g s)\n",
					iout, ncpu, levelmin, nlevelmax, boxlen, time_code, unit_t);
				printf("RAMSES hydro variables:");
				for (const FieldDesc& d : hydro_desc) printf(" %s", d.name.c_str());
				printf("\nRAMSES particle fields:");
				for (const FieldDesc& d : part_desc) printf(" %s", d.name.c_str());
				printf("\n");
			}
			printf("Rank %d: RAMSES files %zu..%zu, leaf cells %zu, particles %zu (DM %llu, stars %llu)\n",
				world_rank, first + 1, first + count, ngas, npart_local,
				(unsigned long long)(type_offset[DM + 1] - type_offset[DM]),
				(unsigned long long)(type_offset[Star + 1] - type_offset[Star]));

#ifdef WITH_OPENMP
			steps_time[1] = omp_get_wtime();
#endif
		}

		void finish_lib() {
			hydro_desc.clear();
			part_desc.clear();
			extra_blocks.clear();
			hydro_idx_density = -1;
			hydro_idx_vel[0] = hydro_idx_vel[1] = hydro_idx_vel[2] = -1;
			ngas = 0;
			std::vector<double>().swap(gas_pos);
			std::vector<float>().swap(gas_dx);
			std::vector<uint8_t>().swap(gas_level);
			gas_vars.clear();
			npart_local = 0;
			std::vector<double>().swap(part_pos);
			std::vector<float>().swap(part_vel);
			std::vector<double>().swap(part_mass);
			std::vector<uint8_t>().swap(part_level);
			part_extra.clear();
			part_extra_names.clear();
			gas_filters.clear();
			part_filters.clear();
			for (int t = 0; t <= PTMax; t++)
				type_offset[t] = 0;
			global_num = 0;
		}

		void print_CPU_steps() {
			//printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		int get_particle_type(uint64_t id) {
			for (int t = 0; t < PTMax; t++) {
				if (id < type_offset[t + 1])
					return t;
			}
			return Other;
		}

		static inline bool is_gas(uint64_t id) {
			return id < ngas;
		}

		static inline size_t part_index(uint64_t id) {
			return (size_t)(id - ngas);
		}

		// Components of a block for one particle; fills v (up to 3), returns the count
		static int fetch(int blocknr, uint64_t id, float* v) {
			const bool gas = is_gas(id);
			switch (blocknr) {
			case Pos:
				v[0] = 1.0f;
				return 1;
			case Mass:
				v[0] = (float)get_particle_mass(id);
				return 1;
			case Rho:
				if (!gas) return 0;
				v[0] = gas_vars[hydro_idx_density][id];
				return 1;
			case Vel:
				if (gas) {
					for (int a = 0; a < 3; a++)
						v[a] = hydro_idx_vel[a] >= 0 ? gas_vars[hydro_idx_vel[a]][id] : 0.0f;
				}
				else {
					for (int a = 0; a < 3; a++)
						v[a] = part_vel[part_index(id) * 3 + a];
				}
				return 3;
			case Level:
				v[0] = gas ? (float)gas_level[id] : (float)part_level[part_index(id)];
				return 1;
			default:
				break;
			}
			int e = blocknr - BTMax;
			if (e < 0 || e >= (int)extra_blocks.size())
				return 0;
			const ExtraBlock& b = extra_blocks[e];
			if (gas) {
				if (b.gas_var < 0) return 0;
				v[0] = gas_vars[b.gas_var][id];
				return 1;
			}
			if (b.part_var < 0 || part_extra[b.part_var].empty()) return 0;
			v[0] = part_extra[b.part_var][part_index(id)];
			return 1;
		}

		float get_particle_norm_value(int blocknr, uint64_t id) {
			float v[3];
			int n = fetch(blocknr, id, v);
			if (n == 1) {
				RETURN_NORM_VALUE(v[0]);
			}
			if (n == 3) {
				RETURN_NORM_VECTOR3(v);
			}
			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value) {
			float v[3];
			int n = fetch(blocknr, id, v);
			if (n == 1) {
				RETURN_ORIG_VALUE(v[0]);
			}
			if (n == 3) {
				RETURN_ORIG_VECTOR3(v);
			}
			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id) {
			float v[3];
			return fetch(blocknr, id, v);
		}

		void get_particle_position(uint64_t id, double* pos) {
			const double* p = is_gas(id) ? &gas_pos[id * 3] : &part_pos[part_index(id) * 3];
			pos[0] = p[0];
			pos[1] = p[1];
			pos[2] = p[2];
		}

		size_t get_local_num_particles() {
			return ngas + npart_local;
		}

		size_t get_global_num_particles() {
			return global_num;
		}

		double get_particle_hsml(uint64_t id) {
			if (is_gas(id))
				return gas_dx[id];
			return boxlen * std::ldexp(1.0, -(int)part_level[part_index(id)]);
		}

		double get_particle_mass(uint64_t id) {
			if (is_gas(id)) {
				double dx = gas_dx[id];
				return (double)gas_vars[hydro_idx_density][id] * dx * dx * dx;
			}
			return part_mass[part_index(id)];
		}

		double get_particle_rho(uint64_t id) {
			if (is_gas(id))
				return gas_vars[hydro_idx_density][id];
			return 0.0;
		}

		int get_particle_rho_blocknr() {
			return Rho;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			const int nblocks = BTMax + (int)extra_blocks.size();
			types_and_blocks.assign((size_t)PTMax * nblocks, 0);

			for (int t = 0; t < PTMax; t++) {
				if (type_offset[t + 1] == type_offset[t])
					continue;
				const bool gas = (t == Gas);
				types_and_blocks[PTMax * Pos + t] = 1;
				types_and_blocks[PTMax * Mass + t] = 1;
				types_and_blocks[PTMax * Vel + t] = 1;
				types_and_blocks[PTMax * Level + t] = 1;
				if (gas)
					types_and_blocks[PTMax * Rho + t] = 1;
				for (size_t e = 0; e < extra_blocks.size(); e++) {
					bool ok = gas ? extra_blocks[e].gas_var >= 0 : extra_blocks[e].part_var >= 0;
					if (ok)
						types_and_blocks[PTMax * (BTMax + (int)e) + t] = 1;
				}
			}
		}

		static const char* type_name(int t) {
			static const char* names[PTMax] = { "Gas", "DM", "Star", "Cloud", "Debris", "Other" };
			return (t >= 0 && t < PTMax) ? names[t] : "Unknown";
		}

		void print_types_and_blocks(std::vector<int>& types_and_blocks) {
			int nblocks = (int)types_and_blocks.size() / PTMax;
			for (int t = 0; t < PTMax; t++) {
				bool any = false;
				for (int b = 0; b < nblocks; b++)
					any = any || types_and_blocks[PTMax * b + t] > 0;
				if (!any)
					continue;
				printf("Type: %s (%d)\n", type_name(t), t);
				for (int b = 0; b < nblocks; b++) {
					if (types_and_blocks[PTMax * b + t] > 0)
						printf("\t%s (%d)\n", get_dataset_name(b).c_str(), b);
				}
			}
		}

		void print_types_and_blocks_local() {
			std::vector<int> types_and_blocks;
			get_types_and_blocks(types_and_blocks);
			printf("\n");
			print_types_and_blocks(types_and_blocks);
		}

		std::string get_dataset_name(int blocknr) {
			switch (blocknr) {
			case Pos: return "Pos";
			case Mass: return "Mass";
			case Rho: return "Rho";
			case Vel: return "Vel";
			case Level: return "Level";
			default: break;
			}
			int e = blocknr - BTMax;
			if (e >= 0 && e < (int)extra_blocks.size())
				return extra_blocks[e].name;
			return "unknown";
		}

		double get_unit_l() { return unit_l; }
		double get_unit_d() { return unit_d; }
		double get_unit_t() { return unit_t; }
		double get_time() { return time_code; }

	} // namespace io
} // namespace ramses
