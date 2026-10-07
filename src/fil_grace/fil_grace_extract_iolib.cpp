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

#include "fil_grace_extract_iolib.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include <hdf5.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "convert_common.h"
#include "reader_return_macros.h"

// Reader data contract (see docs/SpaceConverter_Code_Analysis_2026-08.md §7):
//   positions             -> Cartesian cell centres in code length units (M for GR
//                            runs), x = block corner + (i + 1/2) dx
//   get_particle_mass(id) -> Rho x dx dy dz (coordinate volume of the cell)
//   get_particle_rho(id)  -> the density variable: rho, else dens, else the first
//                            scalar dataset read
//   get_particle_hsml(id) -> the cell size of the cell's block, an adaptive
//                            smoothing length that follows the refinement
//   blocks                -> the datasets as they are in the file (one block per
//                            scalar dataset, then one per vector dataset: magnitude /
//                            components); derived physics (temperature, b^2,
//                            magnetisation, ...) belongs in the consumer (shader)
//   particle id space     -> per rank, [cells] x mirror copies, one type
//   refinement            -> the blocks are the leaves of the p4est forest: they do
//                            not overlap, every cell is kept (--fil-grace-levels drops
//                            whole blocks and leaves their volume empty)
//   ghost zones           -> not in the files
//   mirror                -> --fil-grace-mirror AXES adds the mirror images about the
//                            planes coordinate = 0 (runs with reflection symmetry
//                            store one half only); the mirrored vector component
//                            changes sign (polar vectors: velocity; for the axial
//                            magnetic field only the magnitude is right)
//   MPI split             -> the blocks are split into contiguous ranges of equal
//                            block counts over the ranks; every rank reads its own
//                            blocks only
//
// File format (GRACE src/IO/hdf5_output.cpp, volume_out_<iteration, 6 digits>.h5):
// attributes Time (double) and Iteration (uint) of the file; with n^3 cells per
// block (GRACE blocks are cubic) and nq blocks
//   /Points   double [nq (n+1)^3][3]  block vertices, x fastest; the first and the
//                                     last vertex of a block are its corners
//   /Cells    uint64 [nq n^3][8]      hexahedron connectivity (not read; 8 vertices
//                                     = 3D Cartesian mesh, the only one supported)
//   /<var>    double [nq n^3]         scalar per cell, index i + n (j + n (k + n q))
//   /<var>    double [nq n^3][3]      vector per cell
//   /Level /Rank /Tree_ID /Quad_ID    per cell, written with output_extra_quantities

namespace fil_grace {
	namespace io {

		// ---------------------------------------------------------------------
		// Reader state
		// ---------------------------------------------------------------------
		Options opt;

		// One selected block (p4est leaf) of this rank
		struct Block {
			double lo[3] = { 0.0, 0.0, 0.0 };
			double dx[3] = { 1.0, 1.0, 1.0 };
			int level = 0;
		};

		std::vector<Block> blocks;
		size_t nb = 0;                                  // cells per block edge
		size_t nb3 = 0;                                 // cells per block

		std::vector<std::string> scalar_names;
		std::vector<std::vector<float>> scalar_data;    // [variable][cell]
		std::vector<std::string> vector_names;
		std::vector<std::vector<float>> vector_data;    // [variable][3 x cell]

		int idx_rho = -1;
		long long iteration = -1;
		double time_code = 0.0;

		size_t ncells = 0;                              // cells read (one copy)
		int mirror_axis[3] = { 0, 0, 0 };
		int n_mirror = 0;                               // mirrored axes
		size_t n_copies = 1;                            // 2^n_mirror

		size_t global_num = 0;

		// ---------------------------------------------------------------------
		// Helpers
		// ---------------------------------------------------------------------
		static std::string lower(std::string s) {
			std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
			return s;
		}

		static herr_t collect_link(hid_t, const char* name, const H5L_info_t*, void* op_data) {
			static_cast<std::vector<std::string>*>(op_data)->push_back(name);
			return 0;
		}

		// Dimensions of a dataset (empty when it does not exist)
		static std::vector<hsize_t> dset_dims(hid_t fid, const std::string& name) {
			std::vector<hsize_t> dims;
			if (H5Lexists(fid, name.c_str(), H5P_DEFAULT) <= 0)
				return dims;
			hid_t d = H5Dopen2(fid, name.c_str(), H5P_DEFAULT);
			if (d < 0)
				return dims;
			hid_t sp = H5Dget_space(d);
			const int nd = H5Sget_simple_extent_ndims(sp);
			if (nd > 0) {
				dims.resize(nd);
				H5Sget_simple_extent_dims(sp, dims.data(), nullptr);
			}
			H5Sclose(sp);
			H5Dclose(d);
			return dims;
		}

		static bool read_attr(hid_t obj, const char* name, hid_t memtype, void* buf) {
			if (H5Aexists(obj, name) <= 0)
				return false;
			hid_t a = H5Aopen(obj, name, H5P_DEFAULT);
			if (a < 0)
				return false;
			const bool ok = H5Aread(a, memtype, buf) >= 0;
			H5Aclose(a);
			return ok;
		}

		// Contiguous range of selected blocks: file blocks q .. q + n - 1 -> local blocks dst ..
		struct Run {
			hsize_t q;
			hsize_t n;
			size_t dst;
		};

		// Read the rows of the runs (ncomp values per cell) of a dataset as float
		static void read_runs(hid_t fid, const std::string& name, const std::vector<Run>& runs, int ncomp, std::vector<float>& out) {
			out.assign(ncells * ncomp, 0.0f);
			if (runs.empty())
				return;
			hid_t d = H5Dopen2(fid, name.c_str(), H5P_DEFAULT);
			if (d < 0)
				throw std::runtime_error("FIL_GRACE: cannot open the dataset " + name);
			hid_t fsp = H5Dget_space(d);
			for (const Run& r : runs) {
				hsize_t start[2] = { r.q * nb3, 0 };
				hsize_t count[2] = { r.n * nb3, (hsize_t)ncomp };
				H5Sselect_hyperslab(fsp, H5S_SELECT_SET, start, nullptr, count, nullptr);
				hid_t msp = H5Screate_simple(ncomp > 1 ? 2 : 1, count, nullptr);
				const herr_t err = H5Dread(d, H5T_NATIVE_FLOAT, msp, fsp, H5P_DEFAULT, out.data() + r.dst * nb3 * ncomp);
				H5Sclose(msp);
				if (err < 0) {
					H5Sclose(fsp);
					H5Dclose(d);
					throw std::runtime_error("FIL_GRACE: cannot read the dataset " + name);
				}
			}
			H5Sclose(fsp);
			H5Dclose(d);
		}

		// ---------------------------------------------------------------------
		// Init
		// ---------------------------------------------------------------------
		void init_lib(const Options& options, int world_rank, int world_size) {
			finish_lib();
			opt = options;

			if (opt.file.empty())
				throw std::runtime_error("FIL_GRACE: --fil-grace-file is required");

			// Mirror axes
			n_mirror = 0;
			for (char c : lower(opt.mirror)) {
				if (c < 'x' || c > 'z')
					throw std::runtime_error("FIL_GRACE: --fil-grace-mirror takes the axes x, y, z (e.g. z or xy)");
				const int a = c - 'x';
				if (std::find(mirror_axis, mirror_axis + n_mirror, a) == mirror_axis + n_mirror)
					mirror_axis[n_mirror++] = a;
			}
			n_copies = (size_t)1 << n_mirror;

			H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
			hid_t fid = H5Fopen(opt.file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (fid < 0)
				throw std::runtime_error("FIL_GRACE: cannot open " + opt.file);

			unsigned int it = 0;
			if (read_attr(fid, "Iteration", H5T_NATIVE_UINT, &it))
				iteration = it;
			read_attr(fid, "Time", H5T_NATIVE_DOUBLE, &time_code);

			// Mesh: nq blocks of nb^3 cells
			const std::vector<hsize_t> pdims = dset_dims(fid, "Points");
			const std::vector<hsize_t> cdims = dset_dims(fid, "Cells");
			if (pdims.size() != 2 || cdims.size() != 2 || pdims[1] != 3) {
				H5Fclose(fid);
				throw std::runtime_error("FIL_GRACE: " + opt.file + " has no /Points [N][3] and /Cells datasets (not a GRACE volume output)");
			}
			if (cdims[1] != 8) {
				H5Fclose(fid);
				throw std::runtime_error("FIL_GRACE: only 3D Cartesian meshes are supported (/Cells has " +
					std::to_string(cdims[1]) + " vertices per cell, expected 8)");
			}
			nb = 0;
			if (opt.block_size > 0) {
				nb = (size_t)opt.block_size;
			}
			else {
				// cells / points = nb^3 / (nb + 1)^3
				for (size_t n = 1; n <= 4096 && nb == 0; n++) {
					const double c = (double)n * n * n, p = (double)(n + 1) * (n + 1) * (n + 1);
					if (std::fabs((double)cdims[0] * p - (double)pdims[0] * c) < 0.5 * c)
						nb = n;
				}
			}
			nb3 = nb * nb * nb;
			const size_t npq = (nb + 1) * (nb + 1) * (nb + 1);
			if (nb == 0 || cdims[0] % nb3 != 0 || pdims[0] != (cdims[0] / nb3) * npq) {
				H5Fclose(fid);
				throw std::runtime_error("FIL_GRACE: cannot derive the block size from /Points and /Cells (blocks not cubic?); use --fil-grace-block-size");
			}
			const hsize_t ncells_glob = cdims[0];
			const hsize_t nq_glob = ncells_glob / nb3;

			// This rank's blocks
			const hsize_t q0 = nq_glob * (hsize_t)world_rank / (hsize_t)world_size;
			const hsize_t q1 = nq_glob * (hsize_t)(world_rank + 1) / (hsize_t)world_size;
			const hsize_t nq = q1 - q0;

			// Block corners: the first and the last vertex of every block
			std::vector<double> corners((size_t)nq * 6);
			std::vector<unsigned int> levels((size_t)nq, 0);
			bool have_level = false;
			if (nq > 0) {
				hid_t d = H5Dopen2(fid, "Points", H5P_DEFAULT);
				hid_t fsp = H5Dget_space(d);
				hsize_t start[2] = { q0 * npq, 0 };
				hsize_t stride[2] = { npq, 1 };
				hsize_t count[2] = { nq, 1 };
				hsize_t block[2] = { 1, 3 };
				H5Sselect_hyperslab(fsp, H5S_SELECT_SET, start, stride, count, block);
				start[0] = q0 * npq + npq - 1;
				H5Sselect_hyperslab(fsp, H5S_SELECT_OR, start, stride, count, block);
				hsize_t mdims[2] = { 2 * nq, 3 };
				hid_t msp = H5Screate_simple(2, mdims, nullptr);
				const herr_t err = H5Dread(d, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, corners.data());
				H5Sclose(msp);
				H5Sclose(fsp);
				H5Dclose(d);
				if (err < 0) {
					H5Fclose(fid);
					throw std::runtime_error("FIL_GRACE: cannot read /Points of " + opt.file);
				}

				// Refinement level of the blocks (the first cell of each)
				const std::vector<hsize_t> ldims = dset_dims(fid, "Level");
				if (ldims.size() == 1 && ldims[0] == ncells_glob) {
					hid_t dl = H5Dopen2(fid, "Level", H5P_DEFAULT);
					hid_t lsp = H5Dget_space(dl);
					hsize_t lstart[1] = { q0 * nb3 };
					hsize_t lstride[1] = { nb3 };
					hsize_t lcount[1] = { nq };
					H5Sselect_hyperslab(lsp, H5S_SELECT_SET, lstart, lstride, lcount, nullptr);
					hid_t lmsp = H5Screate_simple(1, lcount, nullptr);
					have_level = H5Dread(dl, H5T_NATIVE_UINT, lmsp, lsp, H5P_DEFAULT, levels.data()) >= 0;
					H5Sclose(lmsp);
					H5Sclose(lsp);
					H5Dclose(dl);
				}
			}

			// Without /Level: levels relative to the coarsest block of the file
			std::vector<double> dx0((size_t)nq);
			double dx_max = 0.0;
			for (hsize_t i = 0; i < nq; i++) {
				dx0[i] = (corners[6 * i + 3] - corners[6 * i]) / (double)nb;
				dx_max = std::max(dx_max, dx0[i]);
			}
			int any_level = have_level ? 1 : 0;
#ifdef WITH_MPI
			{
				double g = 0.0;
				MPI_Allreduce(&dx_max, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				dx_max = g;
				int ga = 0, la = (nq == 0 || have_level) ? 1 : 0;
				MPI_Allreduce(&la, &ga, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
				any_level = ga;
			}
#endif

			// Selection: level range and region (a block is kept when it or one of its
			// mirror images intersects the region)
			std::vector<Run> runs;
			blocks.clear();
			for (hsize_t i = 0; i < nq; i++) {
				Block b;
				bool ok = true;
				for (int a = 0; a < 3; a++) {
					const double lo = corners[6 * i + a], hi = corners[6 * i + 3 + a];
					b.lo[a] = lo;
					b.dx[a] = (hi - lo) / (double)nb;
					if (!(b.dx[a] > 0.0)) {
						H5Fclose(fid);
						throw std::runtime_error("FIL_GRACE: block " + std::to_string(q0 + i) + " is not an axis-aligned box (non-Cartesian coordinates?)");
					}
					if (opt.has_region) {
						const double rlo = opt.region[a], rhi = opt.region[3 + a];
						bool in = lo < rhi && hi > rlo;
						if (!in && std::find(mirror_axis, mirror_axis + n_mirror, a) != mirror_axis + n_mirror)
							in = -hi < rhi && -lo > rlo;
						ok = ok && in;
					}
				}
				b.level = any_level ? (int)levels[i] : (int)std::lround(std::log2(dx_max / dx0[i]));
				if (b.level < opt.level_min || (opt.level_max >= 0 && b.level > opt.level_max))
					ok = false;
				if (!ok)
					continue;
				if (!runs.empty() && runs.back().q + runs.back().n == q0 + i)
					runs.back().n++;
				else
					runs.push_back({ q0 + i, 1, blocks.size() });
				blocks.push_back(b);
			}
			ncells = blocks.size() * nb3;

			// Datasets: the requested ones in the order given, else all (sorted by name)
			std::vector<std::string> names;
			{
				hsize_t idx = 0;
				H5Literate(fid, H5_INDEX_NAME, H5_ITER_INC, &idx, collect_link, &names);
			}
			static const char* mesh_sets[] = { "points", "cells" };
			static const char* extra_sets[] = { "level", "rank", "tree_id", "quad_id" };
			std::vector<std::string> selected;
			if (opt.vars.empty()) {
				for (const std::string& n : names) {
					const std::string l = lower(n);
					bool skip = false;
					for (const char* m : mesh_sets) skip = skip || l == m;
					for (const char* m : extra_sets) skip = skip || l == m;
					if (!skip)
						selected.push_back(n);
				}
			}
			else {
				for (const std::string& v : opt.vars) {
					std::string found;
					for (const std::string& n : names) {
						if (lower(n) == lower(v))
							found = n;
					}
					if (found.empty()) {
						std::string avail;
						for (const std::string& n : names) avail += " " + n;
						H5Fclose(fid);
						throw std::runtime_error("FIL_GRACE: no dataset '" + v + "' in " + opt.file + " (available:" + avail + ")");
					}
					selected.push_back(found);
				}
			}

			scalar_names.clear();
			vector_names.clear();
			for (const std::string& n : selected) {
				const std::vector<hsize_t> dims = dset_dims(fid, n);
				if (dims.size() == 1 && dims[0] == ncells_glob)
					scalar_names.push_back(n);
				else if (dims.size() == 2 && dims[0] == ncells_glob && dims[1] == 3)
					vector_names.push_back(n);
				else if (world_rank == 0)
					printf("Warning: FIL_GRACE dataset %s is not a per-cell scalar or vector (skipped)\n", n.c_str());
			}
			if (scalar_names.empty() && vector_names.empty()) {
				H5Fclose(fid);
				throw std::runtime_error("FIL_GRACE: no per-cell datasets selected in " + opt.file);
			}

			scalar_data.assign(scalar_names.size(), std::vector<float>());
			vector_data.assign(vector_names.size(), std::vector<float>());
			try {
				for (size_t v = 0; v < scalar_names.size(); v++)
					read_runs(fid, scalar_names[v], runs, 1, scalar_data[v]);
				for (size_t v = 0; v < vector_names.size(); v++)
					read_runs(fid, vector_names[v], runs, 3, vector_data[v]);
			}
			catch (...) {
				H5Fclose(fid);
				throw;
			}
			H5Fclose(fid);

			// Density variable
			idx_rho = -1;
			for (const char* cand : { "rho", "dens" }) {
				for (size_t v = 0; v < scalar_names.size() && idx_rho < 0; v++) {
					if (lower(scalar_names[v]) == cand)
						idx_rho = (int)v;
				}
			}
			if (idx_rho < 0 && !scalar_names.empty())
				idx_rho = 0;

			// Summary
			unsigned long long lcount[64] = { 0 }, gcount[64] = { 0 };
			double ldx[64] = { 0.0 }, gdx[64] = { 0.0 };
			for (const Block& b : blocks) {
				const int l = std::min(std::max(b.level, 0), 63);
				lcount[l]++;
				ldx[l] = std::max(ldx[l], b.dx[0]);
			}
			global_num = ncells * n_copies;
#ifdef WITH_MPI
			{
				unsigned long long l = global_num, g = 0;
				MPI_Allreduce(&l, &g, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
				global_num = g;
				MPI_Allreduce(lcount, gcount, 64, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(ldx, gdx, 64, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			}
#else
			std::memcpy(gcount, lcount, sizeof(gcount));
			std::memcpy(gdx, ldx, sizeof(gdx));
#endif

			if (world_rank == 0) {
				printf("FIL_GRACE: %s, iteration %lld, t %g, %llu blocks of %zu^3 cells%s\n", opt.file.c_str(), iteration, time_code,
					(unsigned long long)nq_glob, nb, any_level ? "" : " (no /Level dataset: levels relative to the coarsest block)");
				for (int l = 0; l < 64; l++) {
					if (gcount[l] > 0)
						printf("FIL_GRACE level %d: dx %g, %llu block(s) selected\n", l, gdx[l], gcount[l]);
				}
				printf("FIL_GRACE scalars:");
				for (const std::string& n : scalar_names) printf(" %s", n.c_str());
				if (idx_rho >= 0)
					printf("\nFIL_GRACE density variable: %s", scalar_names[idx_rho].c_str());
				printf("\nFIL_GRACE vectors:");
				for (const std::string& n : vector_names) printf(" %s", n.c_str());
				printf("\n");
				if (n_mirror > 0)
					printf("FIL_GRACE mirror: %zu copies (axes %s)\n", n_copies, opt.mirror.c_str());
			}
			printf("Rank %d: FIL_GRACE blocks %llu..%llu, selected %zu, cells %zu\n", world_rank,
				(unsigned long long)q0, (unsigned long long)q1, blocks.size(), ncells * n_copies);
		}

		void finish_lib() {
			std::vector<Block>().swap(blocks);
			scalar_names.clear();
			vector_names.clear();
			scalar_data.clear();
			vector_data.clear();
			idx_rho = -1;
			ncells = 0;
			n_mirror = 0;
			n_copies = 1;
			global_num = 0;
		}

		void print_CPU_steps() {
		}

		int get_particle_type(uint64_t id) {
			return Cell;
		}

		// Cell of a particle id (ids past the cells are the mirror copies) and its copy number
		static inline size_t cell_of(uint64_t id, unsigned int& copy) {
			if (n_copies == 1) {
				copy = 0;
				return (size_t)id;
			}
			copy = (unsigned int)(id / ncells);
			return (size_t)(id % ncells);
		}

		// Components of a block for one cell; fills v (up to 3), returns the count
		static int fetch(int blocknr, uint64_t id, float* v) {
			unsigned int copy;
			const size_t c = cell_of(id, copy);
			switch (blocknr) {
			case Pos:
				v[0] = 1.0f;
				return 1;
			case Mass:
				v[0] = (float)get_particle_mass(id);
				return 1;
			case Rho:
				if (idx_rho < 0) return 0;
				v[0] = scalar_data[idx_rho][c];
				return 1;
			case Level:
				v[0] = (float)blocks[c / nb3].level;
				return 1;
			default:
				break;
			}
			int e = blocknr - BTMax;
			if (e >= 0 && e < (int)scalar_data.size()) {
				v[0] = scalar_data[e][c];
				return 1;
			}
			e -= (int)scalar_data.size();
			if (e >= 0 && e < (int)vector_data.size()) {
				for (int a = 0; a < 3; a++)
					v[a] = vector_data[e][3 * c + a];
				for (int m = 0; m < n_mirror; m++) {
					if (copy & (1u << m))
						v[mirror_axis[m]] = -v[mirror_axis[m]];
				}
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
			unsigned int copy;
			const size_t c = cell_of(id, copy);
			const Block& b = blocks[c / nb3];
			const size_t r = c % nb3;
			const size_t ijk[3] = { r % nb, (r / nb) % nb, r / (nb * nb) };
			for (int a = 0; a < 3; a++)
				pos[a] = b.lo[a] + ((double)ijk[a] + 0.5) * b.dx[a];
			for (int m = 0; m < n_mirror; m++) {
				if (copy & (1u << m))
					pos[mirror_axis[m]] = -pos[mirror_axis[m]];
			}
		}

		size_t get_local_num_particles() {
			return ncells * n_copies;
		}

		size_t get_global_num_particles() {
			return global_num;
		}

		double get_particle_hsml(uint64_t id) {
			unsigned int copy;
			const Block& b = blocks[cell_of(id, copy) / nb3];
			return std::max(b.dx[0], std::max(b.dx[1], b.dx[2]));
		}

		double get_particle_mass(uint64_t id) {
			if (idx_rho < 0)
				return 0.0;
			unsigned int copy;
			const size_t c = cell_of(id, copy);
			const Block& b = blocks[c / nb3];
			return (double)scalar_data[idx_rho][c] * b.dx[0] * b.dx[1] * b.dx[2];
		}

		double get_particle_rho(uint64_t id) {
			if (idx_rho < 0)
				return 0.0;
			unsigned int copy;
			return scalar_data[idx_rho][cell_of(id, copy)];
		}

		int get_particle_rho_blocknr() {
			return Rho;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			const int nblocks = BTMax + (int)scalar_names.size() + (int)vector_names.size();
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
			if (e >= 0 && e < (int)scalar_names.size())
				return scalar_names[e];
			e -= (int)scalar_names.size();
			if (e >= 0 && e < (int)vector_names.size())
				return vector_names[e];
			return "unknown";
		}

		double get_time() { return time_code; }
		long long get_iteration() { return iteration; }

	} // namespace io
} // namespace fil_grace
