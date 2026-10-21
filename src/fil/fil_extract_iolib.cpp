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

#include "fil_extract_iolib.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>

#include <hdf5.h>
#include <unistd.h>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "convert_common.h"
#include "reader_return_macros.h"

// Reader data contract (see docs/SpaceConverter_Code_Analysis_2026-08.md §7):
//   positions             -> Cartesian grid point coordinates in code length units
//                            (M for GR runs), x = origin + i delta
//   get_particle_mass(id) -> Rho x delta^3 (coordinate volume of the point's cell)
//   get_particle_rho(id)  -> the density variable: rho (HydroBase), else rho_b
//                            (IllinoisGRMHD), else the first variable read
//   get_particle_hsml(id) -> the grid spacing of the point's refinement level, an
//                            adaptive smoothing length that follows the refinement
//   blocks                -> the grid functions as they are in the files (one scalar
//                            block per variable, plus the magnitude / components of
//                            vector groups); derived physics (temperature, b^2,
//                            magnetisation, ...) belongs in the consumer (shader)
//   particle id space     -> per rank, [grid points], one type
//   refinement            -> finest level wins: a point is dropped when its cell
//                            (the delta-cube around it) lies inside a patch of a
//                            finer selected level; points at the edges of finer
//                            patches are kept on both levels (overlap, no gaps)
//   ghost zones           -> the cctk_nghostzones points on the faces a patch
//                            shares with another process (cctk_bbox = 0) are
//                            dropped (they duplicate the neighbour's interior);
//                            outer and refinement boundaries are kept
//   MPI split             -> the patches (level, component) are split into
//                            contiguous ranges of about equal point counts over
//                            the ranks; every rank reads its own patches only
//
// File format (CarpetIOHDF5, Carpet/CarpetIOHDF5/src/Output.cc): one dataset per
// grid function, iteration, time level, refinement level and component, named
//   "<THORN>::<var> it=<iteration> tl=<timelevel>[ m=<map>] rl=<level>[ c=<component>]"
// (vector elements "<THORN>::<var>[i] ..."), 3D arrays [nz][ny][nx] (x fastest,
// float or double) with the attributes origin, delta (double[3]), iorigin
// (int[3]), level, time, cctk_nghostzones (int[3]) and cctk_bbox (int[6],
// lower/upper face per direction: 1 = outer or refinement boundary, 0 = shared
// with another process). One file per variable, per group (one_file_per_group)
// or per process (<name>.file_<N>.h5); 1D/2D output (<name>.x.h5, <name>.xy.h5,
// ...) uses the same dataset names with lower-rank arrays and is skipped.
// The dataset names of every file are cached in "<file>.names" (see list_names).

namespace fil {
	namespace io {

		// ---------------------------------------------------------------------
		// Reader state
		// ---------------------------------------------------------------------
		Options opt;

		struct DsetRef {
			int file = -1;
			std::string name;
		};

		// One component of one refinement level (the geometry of the reference variable)
		struct Patch {
			int level = 0;
			int comp = 0;
			double origin[3] = { 0.0, 0.0, 0.0 };
			double delta[3] = { 1.0, 1.0, 1.0 };
			int n[3] = { 1, 1, 1 };            // x, y, z point counts of the dataset
			int lo[3] = { 0, 0, 0 };           // kept points: lo <= i < hi
			int hi[3] = { 1, 1, 1 };
			double time = 0.0;
			std::vector<DsetRef> dsets;         // per selected variable (file -1: missing)
			size_t weight = 0;                  // kept point count (before the refinement mask)
		};

		std::vector<std::string> file_names;
		std::vector<std::string> var_full;      // selected variables, "THORN::name" as in the files
		std::vector<std::string> var_names;     // block names (short where unique)

		struct Vector {
			std::string name;
			int var[3];
		};
		std::vector<Vector> vectors;

		int idx_rho = -1;
		long long iteration = -1;
		double time_code = 0.0;
		int level_lo = 0, level_hi = 0;

		// Grid points of this rank
		size_t npoints = 0;
		std::vector<double> point_pos;                // 3 per point
		std::vector<float> point_hsml;
		std::vector<uint8_t> point_level;
		std::vector<std::vector<float>> point_vars;   // [variable][point]

		size_t global_num = 0;

		double steps_time[2];

		// ---------------------------------------------------------------------
		// Helpers
		// ---------------------------------------------------------------------
		static std::string lower(std::string s) {
			std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
			return s;
		}

		// "THORN::name" -> "name"
		static std::string short_name(const std::string& full) {
			size_t p = full.rfind("::");
			return p == std::string::npos ? full : full.substr(p + 2);
		}

		static bool ends_with(const std::string& s, const std::string& suffix) {
			return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
		}

		// Carpet 1D/2D output, checkpoints and other non-3D files of a directory
		static bool is_3d_output_name(const std::string& path) {
			const std::string f = std::filesystem::path(path).filename().string();
			if (!ends_with(f, ".h5"))
				return false;
			if (f.rfind("checkpoint", 0) == 0)
				return false;
			static const std::regex lowdim("\\.(x|y|z|d|xy|xz|yz)(\\.file_[0-9]+)?\\.h5$");
			return !std::regex_search(f, lowdim);
		}

		// One entry of the dataset index
		struct DsetName {
			std::string var;
			long long it;
			int tl;
			int map;
			int rl;
			int comp;
		};

		static bool parse_dset_name(const std::string& s, DsetName& d) {
			static const std::regex re("^(\\S+) it=(-?[0-9]+) tl=([0-9]+)(?: ml=[0-9]+)?(?: m=([0-9]+))?(?: rl=([0-9]+))?(?: c=([0-9]+))?$");
			std::smatch m;
			if (!std::regex_match(s, m, re))
				return false;
			d.var = m[1].str();
			d.it = std::stoll(m[2].str());
			d.tl = std::stoi(m[3].str());
			d.map = m[4].matched ? std::stoi(m[4].str()) : 0;
			d.rl = m[5].matched ? std::stoi(m[5].str()) : 0;
			d.comp = m[6].matched ? std::stoi(m[6].str()) : 0;
			return true;
		}

		static herr_t collect_link(hid_t, const char* name, const H5L_info_t*, void* op_data) {
			static_cast<std::vector<std::string>*>(op_data)->push_back(name);
			return 0;
		}

		// Dataset names of a file. Listing a large Carpet file (tens of thousands of
		// datasets) walks its whole group B-tree with small random reads, which takes
		// minutes on Lustre; the names are therefore cached in "<file>.names" next to
		// the file (used while it is not older than the file, written when possible).
		static void list_names(hid_t fid, const std::string& path, std::vector<std::string>& names) {
			namespace fs = std::filesystem;
			const std::string cache = path + ".names";
			std::error_code ec_file, ec_cache;
			const auto t_file = fs::last_write_time(path, ec_file);
			const auto t_cache = fs::last_write_time(cache, ec_cache);
			if (!ec_file && !ec_cache && t_cache >= t_file) {
				std::ifstream in(cache);
				std::string line;
				while (std::getline(in, line)) {
					if (!line.empty())
						names.push_back(line);
				}
				if (!names.empty())
					return;
			}
			hsize_t idx = 0;
			H5Literate(fid, H5_INDEX_NAME, H5_ITER_NATIVE, &idx, collect_link, &names);

			// write atomically (other ranks / processes may do the same)
			const std::string tmp = cache + ".tmp" + std::to_string((long long)getpid());
			{
				std::ofstream out(tmp);
				if (!out)
					return;
				for (const std::string& n : names)
					out << n << '\n';
				if (!out)
					return;
			}
			std::error_code ec;
			fs::rename(tmp, cache, ec);
			if (ec)
				fs::remove(tmp, ec);
		}

		// Attribute with exactly n elements, converted to memtype
		static bool read_attr(hid_t obj, const char* name, hid_t memtype, void* buf, hssize_t n) {
			if (H5Aexists(obj, name) <= 0)
				return false;
			hid_t a = H5Aopen(obj, name, H5P_DEFAULT);
			if (a < 0)
				return false;
			hid_t sp = H5Aget_space(a);
			const bool ok = H5Sget_simple_extent_npoints(sp) == n && H5Aread(a, memtype, buf) >= 0;
			H5Sclose(sp);
			H5Aclose(a);
			return ok;
		}

		// Components per point along x (the fastest index) of a patch read into a buffer
		static void read_patch_var(hid_t file, const Patch& p, const DsetRef& ref, std::vector<float>& buf) {
			hid_t ds = H5Dopen2(file, ref.name.c_str(), H5P_DEFAULT);
			if (ds < 0) {
				throw std::runtime_error("FIL: cannot open dataset '" + ref.name + "' in " + file_names[ref.file]);
			}
			hid_t sp = H5Dget_space(ds);
			hsize_t dims[3] = { 0, 0, 0 };
			const int rank = H5Sget_simple_extent_ndims(sp);
			if (rank == 3)
				H5Sget_simple_extent_dims(sp, dims, nullptr);
			H5Sclose(sp);
			if (rank != 3 || (int)dims[2] != p.n[0] || (int)dims[1] != p.n[1] || (int)dims[0] != p.n[2]) {
				H5Dclose(ds);
				throw std::runtime_error("FIL: dataset '" + ref.name + "' does not match the patch shape of level " +
					std::to_string(p.level) + " component " + std::to_string(p.comp));
			}
			buf.resize((size_t)p.n[0] * p.n[1] * p.n[2]);
			const herr_t st = H5Dread(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
			H5Dclose(ds);
			if (st < 0) {
				throw std::runtime_error("FIL: cannot read dataset '" + ref.name + "' of " + file_names[ref.file]);
			}
		}

		// ---------------------------------------------------------------------
		// Public API
		// ---------------------------------------------------------------------
		void init_lib(const Options& options, int world_rank, int world_size) {
#ifdef WITH_OPENMP
			steps_time[0] = omp_get_wtime();
#endif
			finish_lib();
			opt = options;

			// Silence the HDF5 error stack; failures are reported by the reader
			H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

			// Files: --fil-file as given, --fil-dir all 3D output files (sorted)
			file_names = opt.files;
			for (const std::string& dir : opt.dirs) {
				std::vector<std::string> found;
				std::error_code ec;
				for (const auto& e : std::filesystem::directory_iterator(dir, ec)) {
					if (e.is_regular_file() && is_3d_output_name(e.path().string()))
						found.push_back(e.path().string());
				}
				if (ec) {
					throw std::runtime_error("FIL: cannot list the directory " + dir + ": " + ec.message());
				}
				std::sort(found.begin(), found.end());
				file_names.insert(file_names.end(), found.begin(), found.end());
			}
			if (file_names.empty()) {
				throw std::runtime_error("FIL: no input files (use --fil-file FILE or --fil-dir DIR)");
			}

			// Index: every tl = 0 dataset of every file
			std::vector<hid_t> fids(file_names.size(), -1);
			struct Entry {
				DsetName d;
				int file;
				std::string name;
			};
			std::vector<Entry> entries;
			for (size_t f = 0; f < file_names.size(); f++) {
				fids[f] = H5Fopen(file_names[f].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
				if (fids[f] < 0) {
					throw std::runtime_error("FIL: cannot open the HDF5 file " + file_names[f]);
				}
				std::vector<std::string> names;
				list_names(fids[f], file_names[f], names);
				for (const std::string& n : names) {
					Entry e;
					if (!parse_dset_name(n, e.d) || e.d.tl != 0 || e.d.map != 0)
						continue;
					e.file = (int)f;
					e.name = n;
					entries.push_back(std::move(e));
				}
			}
			if (entries.empty()) {
				throw std::runtime_error("FIL: no Carpet datasets (\"<THORN>::<var> it=... tl=0 rl=... c=...\") in the input files");
			}

			// Variables present, and the ones selected
			std::set<std::string> all_vars;
			for (const Entry& e : entries)
				all_vars.insert(e.d.var);
			if (opt.vars.empty()) {
				var_full.assign(all_vars.begin(), all_vars.end());
			}
			else {
				for (const std::string& want : opt.vars) {
					const std::string w = lower(want);
					std::vector<std::string> hit;
					for (const std::string& v : all_vars) {
						// "vel" also selects the vector elements vel[0], vel[1], vel[2]
						const std::string lv = lower(v), ls = lower(short_name(v));
						if (lv == w || ls == w || ls.rfind(w + "[", 0) == 0 || lv.rfind(w + "[", 0) == 0)
							hit.push_back(v);
					}
					if (hit.empty()) {
						std::string list;
						for (const std::string& v : all_vars) list += " " + v;
						throw std::runtime_error("FIL: variable '" + want + "' is not in the files; present:" + list);
					}
					for (const std::string& v : hit) {
						if (std::find(var_full.begin(), var_full.end(), v) == var_full.end())
							var_full.push_back(v);
					}
				}
			}
			std::map<std::string, int> var_index;
			for (size_t v = 0; v < var_full.size(); v++)
				var_index[var_full[v]] = (int)v;

			// Iteration: requested, else the latest one that holds every selected variable
			std::map<long long, std::set<int>> vars_at_it;
			for (const Entry& e : entries) {
				auto vi = var_index.find(e.d.var);
				if (vi != var_index.end())
					vars_at_it[e.d.it].insert(vi->second);
			}
			if (opt.iteration >= 0) {
				if (!vars_at_it.count(opt.iteration)) {
					std::string near;
					auto it = vars_at_it.lower_bound(opt.iteration);
					if (it != vars_at_it.end()) near += " " + std::to_string(it->first);
					if (it != vars_at_it.begin()) near = " " + std::to_string(std::prev(it)->first) + near;
					throw std::runtime_error("FIL: iteration " + std::to_string(opt.iteration) +
						" is not in the files (nearest:" + near + ")");
				}
				iteration = opt.iteration;
			}
			else {
				iteration = vars_at_it.rbegin()->first;
				for (auto it = vars_at_it.rbegin(); it != vars_at_it.rend(); ++it) {
					if (it->second.size() == var_full.size()) {
						iteration = it->first;
						break;
					}
				}
			}

			// Datasets of the iteration: (level, component) -> per variable; the
			// geometry of a patch comes from the first 3D dataset found for it
			std::map<std::pair<int, int>, Patch> patch_map;
			std::set<std::string> duplicates;
			for (const Entry& e : entries) {
				if (e.d.it != iteration)
					continue;
				auto vi = var_index.find(e.d.var);
				if (vi == var_index.end())
					continue;
				Patch& p = patch_map[{ e.d.rl, e.d.comp }];
				if (p.dsets.empty()) {
					p.level = e.d.rl;
					p.comp = e.d.comp;
					p.dsets.assign(var_full.size(), DsetRef());
					p.n[0] = -1;
				}
				if (p.dsets[vi->second].file >= 0) {
					duplicates.insert(e.d.var);
					continue;
				}
				// 3D datasets only (1D/2D output shares the dataset names)
				hid_t ds = H5Dopen2(fids[e.file], e.name.c_str(), H5P_DEFAULT);
				if (ds < 0)
					continue;
				hid_t sp = H5Dget_space(ds);
				const int rank = H5Sget_simple_extent_ndims(sp);
				hsize_t dims[3] = { 0, 0, 0 };
				if (rank == 3)
					H5Sget_simple_extent_dims(sp, dims, nullptr);
				H5Sclose(sp);
				if (rank != 3) {
					H5Dclose(ds);
					continue;
				}
				if (p.n[0] < 0) {
					p.n[0] = (int)dims[2];
					p.n[1] = (int)dims[1];
					p.n[2] = (int)dims[0];
					if (!read_attr(ds, "origin", H5T_NATIVE_DOUBLE, p.origin, 3) ||
						!read_attr(ds, "delta", H5T_NATIVE_DOUBLE, p.delta, 3)) {
						H5Dclose(ds);
						throw std::runtime_error("FIL: dataset '" + e.name + "' lacks the origin / delta attributes");
					}
					read_attr(ds, "time", H5T_NATIVE_DOUBLE, &p.time, 1);
					int ngh[3] = { 0, 0, 0 };
					int bbox[6] = { 1, 1, 1, 1, 1, 1 };
					read_attr(ds, "cctk_nghostzones", H5T_NATIVE_INT, ngh, 3);
					read_attr(ds, "cctk_bbox", H5T_NATIVE_INT, bbox, 6);
					for (int d = 0; d < 3; d++) {
						p.lo[d] = 0;
						p.hi[d] = p.n[d];
						if (!opt.keep_ghosts) {
							if (bbox[2 * d] == 0) p.lo[d] = std::min(ngh[d], p.n[d]);
							if (bbox[2 * d + 1] == 0) p.hi[d] = std::max(p.n[d] - ngh[d], p.lo[d]);
						}
					}
				}
				H5Dclose(ds);
				p.dsets[vi->second] = { e.file, e.name };
			}
			std::vector<Patch> patches;
			int level_top = -1;
			for (auto& kv : patch_map) {
				if (kv.second.n[0] < 0)
					continue;
				level_top = std::max(level_top, kv.second.level);
				patches.push_back(std::move(kv.second));
			}
			if (patches.empty()) {
				throw std::runtime_error("FIL: no 3D datasets of the selected variables at iteration " + std::to_string(iteration));
			}
			level_lo = std::max(opt.level_min, 0);
			level_hi = (opt.level_max >= 0) ? std::min(opt.level_max, level_top) : level_top;
			patches.erase(std::remove_if(patches.begin(), patches.end(),
				[](const Patch& p) { return p.level < level_lo || p.level > level_hi; }), patches.end());
			if (patches.empty()) {
				throw std::runtime_error("FIL: no patches on the levels " + std::to_string(level_lo) + ".." + std::to_string(level_hi));
			}
			std::sort(patches.begin(), patches.end(), [](const Patch& a, const Patch& b) {
				return a.level != b.level ? a.level < b.level : a.comp < b.comp;
				});
			time_code = patches.back().time;

			// Block names: the short variable name where it is unique
			std::map<std::string, int> short_count;
			for (const std::string& v : var_full)
				short_count[short_name(v)]++;
			var_names.clear();
			for (const std::string& v : var_full)
				var_names.push_back(short_count[short_name(v)] > 1 ? v : short_name(v));

			// Density variable
			idx_rho = -1;
			for (const char* cand : { "rho", "rho_b" }) {
				for (size_t v = 0; v < var_full.size() && idx_rho < 0; v++) {
					if (lower(short_name(var_full[v])) == cand)
						idx_rho = (int)v;
				}
			}
			if (idx_rho < 0)
				idx_rho = 0;

			// Vector groups: <base>[0..2] or <base>x/y/z (base not ending in x, y, z,
			// so that the metric components gxx, gxy, ... do not form vectors)
			vectors.clear();
			{
				std::map<std::string, std::array<int, 3>> groups;
				static const std::regex re_idx("^(.+)\\[([0-2])\\]$");
				static const std::regex re_xyz("^(.*[^xyzXYZ:])([xyz])$");
				for (size_t v = 0; v < var_full.size(); v++) {
					std::smatch m;
					std::string base;
					int c = -1;
					if (std::regex_match(var_full[v], m, re_idx)) {
						base = m[1].str();
						c = m[2].str()[0] - '0';
					}
					else if (std::regex_match(var_full[v], m, re_xyz)) {
						base = m[1].str();
						c = m[2].str()[0] - 'x';
					}
					if (c < 0)
						continue;
					auto g = groups.find(base);
					if (g == groups.end())
						g = groups.emplace(base, std::array<int, 3>{ -1, -1, -1 }).first;
					g->second[c] = (int)v;
				}
				std::map<std::string, int> vec_count;
				for (const auto& g : groups)
					vec_count[short_name(g.first)]++;
				for (const auto& g : groups) {
					if (g.second[0] < 0 || g.second[1] < 0 || g.second[2] < 0)
						continue;
					Vector vec;
					vec.name = vec_count[short_name(g.first)] > 1 ? g.first : short_name(g.first);
					for (int c = 0; c < 3; c++)
						vec.var[c] = g.second[c];
					vectors.push_back(vec);
				}
			}

			// Refinement mask needs, for every patch, the finer patches that overlap it
			struct Box {
				int level;
				double lo[3], hi[3];    // the region of the kept points, extended by delta / 2
			};
			std::vector<Box> boxes(patches.size());
			for (size_t i = 0; i < patches.size(); i++) {
				const Patch& p = patches[i];
				boxes[i].level = p.level;
				for (int d = 0; d < 3; d++) {
					boxes[i].lo[d] = p.origin[d] + (p.lo[d] - 0.5) * p.delta[d];
					boxes[i].hi[d] = p.origin[d] + (p.hi[d] - 0.5) * p.delta[d];
				}
				size_t w = 1;
				for (int d = 0; d < 3; d++)
					w *= (size_t)std::max(p.hi[d] - p.lo[d], 0);
				patches[i].weight = w;
			}

			// This rank's patches: contiguous ranges of about equal point counts
			size_t total = 0;
			for (const Patch& p : patches)
				total += p.weight;
			size_t first = patches.size(), last = patches.size();
			{
				size_t acc = 0;
				for (size_t i = 0; i < patches.size(); i++) {
					// the rank that owns the patch's first point
					const size_t owner = total > 0 ? std::min((size_t)world_size - 1, (acc * (size_t)world_size) / total) : 0;
					if (owner == (size_t)world_rank) {
						if (first == patches.size()) first = i;
						last = i + 1;
					}
					acc += patches[i].weight;
				}
			}

			// Read
			point_vars.assign(var_full.size(), std::vector<float>());
			std::vector<std::vector<float>> buf(var_full.size());
			std::vector<uint8_t> covered;
			size_t points_level[64] = { 0 };
			for (size_t i = first; i < last; i++) {
				const Patch& p = patches[i];
				const size_t nxy = (size_t)p.n[0] * p.n[1];

				// cells covered by a finer patch
				std::vector<size_t> finer;
				for (size_t j = 0; j < patches.size(); j++) {
					if (boxes[j].level <= p.level)
						continue;
					bool overlap = true;
					for (int d = 0; d < 3; d++) {
						const double plo = p.origin[d] + (p.lo[d] - 0.5) * p.delta[d];
						const double phi = p.origin[d] + (p.hi[d] - 0.5) * p.delta[d];
						overlap = overlap && boxes[j].lo[d] < phi && boxes[j].hi[d] > plo;
					}
					if (overlap)
						finer.push_back(j);
				}
				covered.assign(nxy * p.n[2], 0);
				if (!finer.empty()) {
					const double tol = 1e-6 * p.delta[0];
#pragma omp parallel for schedule(static)
					for (int k = p.lo[2]; k < p.hi[2]; k++) {
						double x[3];
						x[2] = p.origin[2] + k * p.delta[2];
						for (int jj = p.lo[1]; jj < p.hi[1]; jj++) {
							x[1] = p.origin[1] + jj * p.delta[1];
							for (int ii = p.lo[0]; ii < p.hi[0]; ii++) {
								x[0] = p.origin[0] + ii * p.delta[0];
								for (size_t f : finer) {
									bool in = true;
									for (int d = 0; d < 3 && in; d++) {
										const double h = 0.5 * p.delta[d];
										in = x[d] - h >= boxes[f].lo[d] - tol && x[d] + h <= boxes[f].hi[d] + tol;
									}
									if (in) {
										covered[(size_t)ii + (size_t)p.n[0] * jj + nxy * k] = 1;
										break;
									}
								}
							}
						}
					}
				}

				for (size_t v = 0; v < var_full.size(); v++) {
					if (p.dsets[v].file < 0) {
						if (world_rank == 0 || world_size == 1) {
							printf("Warning: FIL variable %s has no data on level %d component %d at iteration %lld (zero)\n",
								var_full[v].c_str(), p.level, p.comp, iteration);
						}
						buf[v].assign(nxy * p.n[2], 0.0f);
						continue;
					}
					read_patch_var(fids[p.dsets[v].file], p, p.dsets[v], buf[v]);
				}

				const float hsml = (float)std::max(p.delta[0], std::max(p.delta[1], p.delta[2]));
				for (int k = p.lo[2]; k < p.hi[2]; k++) {
					for (int jj = p.lo[1]; jj < p.hi[1]; jj++) {
						for (int ii = p.lo[0]; ii < p.hi[0]; ii++) {
							const size_t c = (size_t)ii + (size_t)p.n[0] * jj + nxy * k;
							if (covered[c])
								continue;
							point_pos.push_back(p.origin[0] + ii * p.delta[0]);
							point_pos.push_back(p.origin[1] + jj * p.delta[1]);
							point_pos.push_back(p.origin[2] + k * p.delta[2]);
							point_hsml.push_back(hsml);
							point_level.push_back((uint8_t)p.level);
							for (size_t v = 0; v < var_full.size(); v++)
								point_vars[v].push_back(buf[v][c]);
							points_level[std::min(p.level, 63)]++;
						}
					}
				}
			}
			npoints = point_hsml.size();

			for (hid_t f : fids) {
				if (f >= 0) H5Fclose(f);
			}

			global_num = npoints;
#ifdef WITH_MPI
			unsigned long long l = npoints, g = 0;
			MPI_Allreduce(&l, &g, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
			global_num = g;
#endif

			if (world_rank == 0) {
				printf("FIL (Carpet HDF5): %zu file(s), iteration %lld, t %g, levels %d..%d of 0..%d, %zu patches\n",
					file_names.size(), iteration, time_code, level_lo, level_hi, level_top, patches.size());
				int lev = -1;
				for (const Patch& p : patches) {
					if (p.level == lev)
						continue;
					lev = p.level;
					int ncomp = 0;
					double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
					for (size_t j = 0; j < patches.size(); j++) {
						if (patches[j].level != lev)
							continue;
						ncomp++;
						for (int d = 0; d < 3; d++) {
							lo[d] = std::min(lo[d], patches[j].origin[d] + patches[j].lo[d] * patches[j].delta[d]);
							hi[d] = std::max(hi[d], patches[j].origin[d] + (patches[j].hi[d] - 1) * patches[j].delta[d]);
						}
					}
					printf("FIL level %d: delta %g, %d component(s), extent [%g %g %g] .. [%g %g %g]\n",
						lev, p.delta[0], ncomp, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
				}
				printf("FIL variables:");
				for (size_t v = 0; v < var_full.size(); v++) printf(" %s", var_full[v].c_str());
				printf("\nFIL density variable: %s\n", var_full[idx_rho].c_str());
				if (!vectors.empty()) {
					printf("FIL vectors:");
					for (const Vector& vec : vectors) printf(" %s", vec.name.c_str());
					printf("\n");
				}
				for (const std::string& d : duplicates)
					printf("Warning: FIL variable %s has several datasets for the same patch (restarts?); the first file wins\n", d.c_str());
			}
			printf("Rank %d: FIL patches %zu..%zu of %zu, grid points %zu (per level:", world_rank,
				first, last, patches.size(), npoints);
			for (int l = level_lo; l <= level_hi && l < 64; l++) printf(" %zu", points_level[l]);
			printf(")\n");

#ifdef WITH_OPENMP
			steps_time[1] = omp_get_wtime();
#endif
		}

		void finish_lib() {
			file_names.clear();
			var_full.clear();
			var_names.clear();
			vectors.clear();
			idx_rho = -1;
			npoints = 0;
			std::vector<double>().swap(point_pos);
			std::vector<float>().swap(point_hsml);
			std::vector<uint8_t>().swap(point_level);
			point_vars.clear();
			global_num = 0;
		}

		void print_CPU_steps() {
			//printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		int get_particle_type(uint64_t id) {
			return Cell;
		}

		// Components of a block for one point; fills v (up to 3), returns the count
		static int fetch(int blocknr, uint64_t id, float* v) {
			switch (blocknr) {
			case Pos:
				v[0] = 1.0f;
				return 1;
			case Mass:
				v[0] = (float)get_particle_mass(id);
				return 1;
			case Rho:
				if (idx_rho < 0) return 0;
				v[0] = point_vars[idx_rho][id];
				return 1;
			case Level:
				v[0] = (float)point_level[id];
				return 1;
			default:
				break;
			}
			int e = blocknr - BTMax;
			if (e >= 0 && e < (int)point_vars.size()) {
				v[0] = point_vars[e][id];
				return 1;
			}
			e -= (int)point_vars.size();
			if (e >= 0 && e < (int)vectors.size()) {
				for (int a = 0; a < 3; a++)
					v[a] = point_vars[vectors[e].var[a]][id];
				return 3;
			}
			return 0;
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
			pos[0] = point_pos[id * 3 + 0];
			pos[1] = point_pos[id * 3 + 1];
			pos[2] = point_pos[id * 3 + 2];
		}

		size_t get_local_num_particles() {
			return npoints;
		}

		size_t get_global_num_particles() {
			return global_num;
		}

		double get_particle_hsml(uint64_t id) {
			return point_hsml[id];
		}

		double get_particle_mass(uint64_t id) {
			const double h = point_hsml[id];
			return idx_rho < 0 ? 0.0 : (double)point_vars[idx_rho][id] * h * h * h;
		}

		double get_particle_rho(uint64_t id) {
			return idx_rho < 0 ? 0.0 : point_vars[idx_rho][id];
		}

		int get_particle_rho_blocknr() {
			return Rho;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			const int nblocks = BTMax + (int)point_vars.size() + (int)vectors.size();
			types_and_blocks.assign((size_t)PTMax * nblocks, 0);
			if (global_num == 0)
				return;
			for (int b = 0; b < nblocks; b++)
				types_and_blocks[PTMax * b + Cell] = 1;
		}

		void print_types_and_blocks(std::vector<int>& types_and_blocks) {
			int nblocks = (int)types_and_blocks.size() / PTMax;
			for (int t = 0; t < PTMax; t++) {
				bool any = false;
				for (int b = 0; b < nblocks; b++)
					any = any || types_and_blocks[PTMax * b + t] > 0;
				if (!any)
					continue;
				printf("Type: Cell (%d)\n", t);
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
			case Level: return "Level";
			default: break;
			}
			int e = blocknr - BTMax;
			if (e >= 0 && e < (int)var_names.size())
				return var_names[e];
			e -= (int)var_names.size();
			if (e >= 0 && e < (int)vectors.size())
				return vectors[e].name;
			return "unknown";
		}

		double get_time() { return time_code; }
		long long get_iteration() { return iteration; }

	} // namespace io
} // namespace fil
