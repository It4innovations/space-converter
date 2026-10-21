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

#include "bhac_extract_iolib.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
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
//   positions             -> Cartesian cell centres in code length units (M for GRMHD);
//                            spherical grids: x = r sin(theta) cos(phi), y = r sin(theta)
//                            sin(phi), z = r cos(theta) with (r, theta, phi) from the code
//                            coordinates (--bhac-coord, flat-space embedding)
//   get_particle_mass(id) -> Rho x cell volume (flat-space coordinate volume)
//   get_particle_rho(id)  -> Rho: d / lfac when both are stored (the rest-mass density
//                            of the conserved D = lfac rho), else the first variable
//   get_particle_hsml(id) -> largest Cartesian extent of the cell (max(dr, r dtheta,
//                            r sin(theta) dphi)), an adaptive smoothing length that
//                            follows the refinement and the stretching of the grid
//   blocks                -> the stored variables as they are in the file (conserved
//                            variables + auxiliaries, names from wnames of the .par);
//                            derived physics (pressure, b^2, sigma, ...) belongs in the
//                            consumer (shader), not here
//   particle id space     -> per rank, [leaf cells], one type
//   MPI split             -> the leaf blocks (file order) are split into contiguous
//                            ranges over the ranks; every rank reads its own blocks only
//
// File format (BHAC amrio.t write_snapshot*, native endianness, no record markers):
//   nleafs x { w(nx1,nx2,nx3,1:nw) [+ ws(0:nx1,0:nx2,0:nx3,1:nws) if staggered] }
//       blocks in Morton order = the order of the leaves in the forest below,
//       each variable a Fortran-ordered (x1 fastest) array of doubles without ghost cells
//   forest: one 4-byte logical per tree node, depth first; roots in ig1-fastest order
//       over ng1 x ng2 x ng3 level-1 blocks, children in ic1-fastest order
//   nx(1:ndim) int32, eqpar(1:neqpar) double,
//   nleafs, levmax, ndim, ndir, nw, nws, neqpar, it  int32, t double
// The base grid and the domain are not in the file: they come from the .par
// file (nxlone^D, xprobmin^D, xprobmax^D, typeaxial; spherical: x2, x3 in units of 2 pi).

namespace bhac {
	namespace io {

		// ---------------------------------------------------------------------
		// Reader state
		// ---------------------------------------------------------------------
		Options opt;

		// .par settings
		int nxlone[3] = { 1, 1, 1 };
		double xprobmin[3] = { 0.0, 0.0, 0.0 };
		double xprobmax[3] = { 1.0, 1.0, 1.0 };
		std::string typeaxial = "slab";
		std::vector<std::string> var_names;

		// .dat header
		int ndim = 3, ndir = 3, nw = 0, nws = 0, neqpar = 0, nleafs = 0, levmax = 1, iteration = 0;
		int nx[3] = { 1, 1, 1 };
		double time_code = 0.0;
		std::vector<double> eqpar;

		// variables with a fixed block of their own
		int idx_d = -1, idx_lfac = -1;
		int idx_b[3] = { -1, -1, -1 };

		// Leaf cells of this rank
		size_t ncells = 0;
		std::vector<double> cell_pos;                // 3 per cell
		std::vector<float> cell_hsml;
		std::vector<float> cell_vol;
		std::vector<float> cell_rho;
		std::vector<float> cell_b;                   // 3 per cell (empty: no b1 b2 b3)
		std::vector<uint8_t> cell_level;
		std::vector<std::vector<float>> cell_vars;   // [variable][cell]

		size_t global_num = 0;

		double steps_time[2];

		// ---------------------------------------------------------------------
		// Helpers
		// ---------------------------------------------------------------------
		static void partition_range(size_t total, size_t parts, size_t idx, size_t& start, size_t& count) {
			size_t base = total / parts;
			size_t rem = total % parts;
			count = base + (idx < rem ? 1 : 0);
			start = idx * base + std::min(idx, rem);
		}

		static std::string lower(std::string s) {
			std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
			return s;
		}

		// Fortran real literal ("0.5d0", "1.0D-3") -> double
		static double fortran_double(std::string s) {
			for (char& c : s)
				if (c == 'd' || c == 'D') c = 'e';
			return std::stod(s);
		}

		// The key = value pairs of a Fortran namelist file (keys lower case,
		// quotes and '!' comments removed). Values spread over several lines
		// (e.g. typeB) are not needed and only their first line is kept.
		static std::map<std::string, std::string> read_namelist(const std::string& path) {
			std::ifstream f(path);
			if (!f) {
				throw std::runtime_error("BHAC: cannot open the .par file " + path);
			}
			std::map<std::string, std::string> kv;
			const std::regex pair_re("([A-Za-z_][A-Za-z0-9_]*(?:\\([0-9, ]+\\))?)\\s*=\\s*('[^']*'|\"[^\"]*\"|[^,\\s]+)");
			std::string line;
			while (std::getline(f, line)) {
				// strip a '!' comment that is not inside quotes
				bool quoted = false;
				char q = 0;
				for (size_t i = 0; i < line.size(); i++) {
					char c = line[i];
					if (quoted) {
						if (c == q) quoted = false;
					}
					else if (c == '\'' || c == '"') {
						quoted = true;
						q = c;
					}
					else if (c == '!') {
						line.erase(i);
						break;
					}
				}
				for (std::sregex_iterator it(line.begin(), line.end(), pair_re), end; it != end; ++it) {
					std::string key = lower((*it)[1].str());
					key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
					std::string val = (*it)[2].str();
					if (val.size() >= 2 && (val[0] == '\'' || val[0] == '"'))
						val = val.substr(1, val.size() - 2);
					kv[key] = val;
				}
			}
			return kv;
		}

		static void read_par(const std::string& path) {
			std::map<std::string, std::string> kv = read_namelist(path);
			for (int d = 0; d < ndim; d++) {
				const std::string n = std::to_string(d + 1);
				if (!kv.count("nxlone" + n) || !kv.count("xprobmin" + n) || !kv.count("xprobmax" + n)) {
					throw std::runtime_error("BHAC: " + path + " lacks nxlone" + n + ", xprobmin" + n + " or xprobmax" + n);
				}
				nxlone[d] = std::stoi(kv["nxlone" + n]);
				xprobmin[d] = fortran_double(kv["xprobmin" + n]);
				xprobmax[d] = fortran_double(kv["xprobmax" + n]);
			}
			typeaxial = kv.count("typeaxial") ? lower(kv["typeaxial"]) : std::string("slab");
			if (typeaxial == "spherical") {
				// amrio.t readparameters: xprob^LIM^DE = xprob^LIM^DE * two * dpi
				for (int d = 1; d < ndim; d++) {
					xprobmin[d] *= 2.0 * M_PI;
					xprobmax[d] *= 2.0 * M_PI;
				}
			}
			else if (typeaxial != "slab") {
				throw std::runtime_error("BHAC: typeaxial '" + typeaxial + "' is not supported (slab, spherical)");
			}

			var_names.clear();
			if (kv.count("wnames")) {
				std::stringstream ss(kv["wnames"]);
				std::string n;
				while (ss >> n)
					var_names.push_back(n);
			}
		}

		// Tail of the .dat file (see the format note above)
		static void read_header(std::ifstream& f, size_t file_size, size_t& tail_size) {
			const size_t tail_fixed = 8 * sizeof(int32_t) + sizeof(double);
			if (file_size < tail_fixed) {
				throw std::runtime_error("BHAC: " + opt.dat_file + " is too short for a snapshot");
			}
			int32_t v[8];
			f.seekg((std::streamoff)(file_size - tail_fixed));
			f.read(reinterpret_cast<char*>(v), sizeof(v));
			f.read(reinterpret_cast<char*>(&time_code), sizeof(double));
			nleafs = v[0];
			levmax = v[1];
			ndim = v[2];
			ndir = v[3];
			nw = v[4];
			nws = v[5];
			neqpar = v[6];
			iteration = v[7];
			if (!f || ndim < 1 || ndim > 3 || ndir < ndim || ndir > 3 || nw < 1 || nws < 0 ||
				neqpar < 0 || nleafs < 1 || levmax < 1 || levmax > 30) {
				throw std::runtime_error("BHAC: " + opt.dat_file + " has no valid snapshot tail (ndim " + std::to_string(ndim) +
					", nw " + std::to_string(nw) + ", nleafs " + std::to_string(nleafs) + "); only native-endian .dat files are supported");
			}

			tail_size = tail_fixed + (size_t)ndim * sizeof(int32_t) + (size_t)neqpar * sizeof(double);
			if (file_size < tail_size) {
				throw std::runtime_error("BHAC: " + opt.dat_file + " is too short for its header");
			}
			f.seekg((std::streamoff)(file_size - tail_size));
			int32_t n[3] = { 1, 1, 1 };
			f.read(reinterpret_cast<char*>(n), (std::streamsize)(ndim * sizeof(int32_t)));
			eqpar.assign((size_t)neqpar, 0.0);
			if (neqpar > 0)
				f.read(reinterpret_cast<char*>(eqpar.data()), (std::streamsize)(neqpar * sizeof(double)));
			if (!f) {
				throw std::runtime_error("BHAC: cannot read the header of " + opt.dat_file);
			}
			for (int d = 0; d < 3; d++)
				nx[d] = (d < ndim) ? n[d] : 1;
		}

		// One leaf block of the forest
		struct Leaf {
			int level;
			int ig[3];
		};

		// Depth-first walk of the forest (forest.t read_node): leaves in file order
		static void walk_node(const std::vector<int32_t>& nodes, size_t& pos, int level, const int ig[3], std::vector<Leaf>& leaves) {
			if (pos >= nodes.size()) {
				throw std::runtime_error("BHAC: the forest of " + opt.dat_file + " ends early");
			}
			const bool leaf = nodes[pos++] != 0;
			if (leaf) {
				leaves.push_back({ level, { ig[0], ig[1], ig[2] } });
				return;
			}
			const int n3 = ndim > 2 ? 2 : 1, n2 = ndim > 1 ? 2 : 1;
			for (int ic3 = 1; ic3 <= n3; ic3++) {
				for (int ic2 = 1; ic2 <= n2; ic2++) {
					for (int ic1 = 1; ic1 <= 2; ic1++) {
						const int ic[3] = { ic1, ic2, ic3 };
						int child[3];
						for (int d = 0; d < 3; d++)
							child[d] = (d < ndim) ? 2 * (ig[d] - 1) + ic[d] : 1;
						walk_node(nodes, pos, level + 1, child, leaves);
					}
				}
			}
		}

		// Cell centre (code coordinates) -> Cartesian position, cell extents
		// and the flat-space Jacobian columns used for vectors
		struct CellGeom {
			double pos[3];
			double ext[3];      // Cartesian extent along the three code directions
			double r;           // distance from the origin
			double e[3][3];     // e[k] = Cartesian unit vector of code direction k
			double scale[3];    // length per code coordinate unit along direction k
		};

		static void cell_geometry(const double x[3], const double dx[3], CellGeom& g) {
			if (typeaxial == "slab" || opt.coord == Coord::Cart) {
				for (int k = 0; k < 3; k++) {
					g.pos[k] = x[k];
					g.ext[k] = (k < ndim) ? dx[k] : 0.0;
					g.scale[k] = 1.0;
					for (int a = 0; a < 3; a++)
						g.e[k][a] = (a == k) ? 1.0 : 0.0;
				}
				g.r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
				return;
			}

			// spherical: (x1, x2, x3) -> (r, theta, phi); 1D: equatorial, 2D: phi = 0
			double r, drdx1, th, dthdx2;
			const double x2 = (ndim > 1) ? x[1] : 0.5 * M_PI;
			const double ph = (ndim > 2) ? x[2] : 0.0;
			if (opt.coord == Coord::MKS) {
				drdx1 = std::exp(x[0]);
				r = opt.mks_r0 + drdx1;
				th = x2 + 0.5 * opt.mks_h * std::sin(2.0 * x2);
				dthdx2 = 1.0 + opt.mks_h * std::cos(2.0 * x2);
			}
			else {
				r = x[0];
				drdx1 = 1.0;
				th = x2;
				dthdx2 = 1.0;
			}
			const double st = std::sin(th), ct = std::cos(th), sp = std::sin(ph), cp = std::cos(ph);
			g.r = r;
			g.pos[0] = r * st * cp;
			g.pos[1] = r * st * sp;
			g.pos[2] = r * ct;
			g.scale[0] = drdx1;
			g.scale[1] = r * dthdx2;
			g.scale[2] = r * st;
			const double e[3][3] = {
				{ st * cp, st * sp, ct },       // r hat
				{ ct * cp, ct * sp, -st },      // theta hat
				{ -sp, cp, 0.0 },               // phi hat
			};
			std::memcpy(g.e, e, sizeof(e));
			g.ext[0] = g.scale[0] * dx[0];
			g.ext[1] = (ndim > 1) ? g.scale[1] * dx[1] : 0.0;
			// 2D (axisymmetric) cells stand for the whole ring
			g.ext[2] = (ndim > 2) ? g.scale[2] * dx[2] : 2.0 * M_PI * r * st;
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

			std::ifstream f(opt.dat_file, std::ios::binary);
			if (!f) {
				throw std::runtime_error("BHAC: cannot open " + opt.dat_file);
			}
			f.seekg(0, std::ios::end);
			const size_t file_size = (size_t)f.tellg();

			size_t tail_size = 0;
			read_header(f, file_size, tail_size);
			if (opt.par_file.empty()) {
				throw std::runtime_error("BHAC: --bhac-par is required (the base grid and the domain are not stored in the .dat file)");
			}
			read_par(opt.par_file);
			if (!opt.coord_set)
				opt.coord = (typeaxial == "spherical") ? Coord::MKS : Coord::Cart;

			// Variable names: wnames of the .par, else w1..wN
			if ((int)var_names.size() != nw) {
				if (world_rank == 0 && !var_names.empty()) {
					printf("Warning: BHAC wnames has %zu names, the snapshot %d variables; using w1..w%d\n",
						var_names.size(), nw, nw);
				}
				var_names.clear();
				for (int v = 0; v < nw; v++)
					var_names.push_back("w" + std::to_string(v + 1));
			}
			for (int v = 0; v < nw; v++) {
				const std::string n = lower(var_names[v]);
				if (n == "d") idx_d = v;
				if (n == "lfac") idx_lfac = v;
				if (n == "b1") idx_b[0] = v;
				if (n == "b2") idx_b[1] = v;
				if (n == "b3") idx_b[2] = v;
			}
			const bool have_b = idx_b[0] >= 0 && idx_b[1] >= 0 && idx_b[2] >= 0;

			// Block layout
			size_t ncell_block = 1, nstg_block = 1;
			for (int d = 0; d < ndim; d++) {
				ncell_block *= (size_t)nx[d];
				nstg_block *= (size_t)nx[d] + 1;
			}
			const size_t block_bytes = ncell_block * (size_t)nw * sizeof(double) +
				(nws > 0 ? nstg_block * (size_t)nws * sizeof(double) : 0);
			const size_t forest_start = block_bytes * (size_t)nleafs;
			if (forest_start + (size_t)nleafs * sizeof(int32_t) + tail_size > file_size) {
				throw std::runtime_error("BHAC: " + opt.dat_file + " is shorter than its " + std::to_string(nleafs) + " blocks");
			}
			const size_t nnodes = (file_size - tail_size - forest_start) / sizeof(int32_t);

			// Base grid in blocks
			int ng[3] = { 1, 1, 1 };
			double dxlone[3] = { 1.0, 1.0, 1.0 };
			for (int d = 0; d < ndim; d++) {
				if (nxlone[d] % nx[d] != 0) {
					throw std::runtime_error("BHAC: nxlone" + std::to_string(d + 1) + " = " + std::to_string(nxlone[d]) +
						" is not a multiple of the block size " + std::to_string(nx[d]));
				}
				ng[d] = nxlone[d] / nx[d];
				dxlone[d] = (xprobmax[d] - xprobmin[d]) / (double)nxlone[d];
			}

			// Forest -> leaves in file order
			std::vector<int32_t> nodes(nnodes);
			f.seekg((std::streamoff)forest_start);
			f.read(reinterpret_cast<char*>(nodes.data()), (std::streamsize)(nnodes * sizeof(int32_t)));
			if (!f) {
				throw std::runtime_error("BHAC: cannot read the forest of " + opt.dat_file);
			}
			std::vector<Leaf> leaves;
			leaves.reserve((size_t)nleafs);
			size_t pos = 0;
			for (int ig3 = 1; ig3 <= ng[2]; ig3++) {
				for (int ig2 = 1; ig2 <= ng[1]; ig2++) {
					for (int ig1 = 1; ig1 <= ng[0]; ig1++) {
						const int ig[3] = { ig1, ig2, ig3 };
						walk_node(nodes, pos, 1, ig, leaves);
					}
				}
			}
			if (pos != nnodes || (int)leaves.size() != nleafs) {
				throw std::runtime_error("BHAC: the forest of " + opt.dat_file + " does not match the base grid of the .par file (" +
					std::to_string(leaves.size()) + " leaves of " + std::to_string(nleafs) + ", " + std::to_string(pos) +
					" of " + std::to_string(nnodes) + " nodes read)");
			}

			// Blocks that can hold cells within --bhac-rrange (r is monotonic in x1
			// on spherical grids; slab blocks are all kept), so that the ranges of
			// the ranks below are all non-empty whenever there are enough blocks
			std::vector<size_t> selected;
			selected.reserve(leaves.size());
			for (size_t b = 0; b < leaves.size(); b++) {
				if (typeaxial == "spherical" && (opt.rmin > 0.0 || opt.rmax > 0.0)) {
					const double dx1 = std::ldexp(dxlone[0], 1 - leaves[b].level);
					double x1[2] = { xprobmin[0] + (double)(leaves[b].ig[0] - 1) * nx[0] * dx1, 0.0 };
					x1[1] = x1[0] + nx[0] * dx1;
					double r[2];
					for (int s = 0; s < 2; s++)
						r[s] = (opt.coord == Coord::MKS) ? opt.mks_r0 + std::exp(x1[s]) : x1[s];
					if (r[1] < opt.rmin || (opt.rmax > 0.0 && r[0] > opt.rmax))
						continue;
				}
				selected.push_back(b);
			}

			// This rank's blocks
			size_t first = 0, count = 0;
			partition_range(selected.size(), (size_t)world_size, (size_t)world_rank, first, count);

			cell_vars.assign((size_t)nw, std::vector<float>());
			std::vector<double> w(ncell_block * (size_t)nw);
			for (size_t k = first; k < first + count; k++) {
				const size_t b = selected[k];
				const Leaf& leaf = leaves[b];
				f.seekg((std::streamoff)(block_bytes * b));
				f.read(reinterpret_cast<char*>(w.data()), (std::streamsize)(w.size() * sizeof(double)));
				if (!f) {
					throw std::runtime_error("BHAC: cannot read block " + std::to_string(b) + " of " + opt.dat_file);
				}

				double dx[3] = { 0.0, 0.0, 0.0 }, xmin[3] = { 0.0, 0.0, 0.0 };
				for (int d = 0; d < ndim; d++) {
					dx[d] = std::ldexp(dxlone[d], 1 - leaf.level);
					xmin[d] = xprobmin[d] + (double)(leaf.ig[d] - 1) * nx[d] * dx[d];
				}

				for (int k = 0; k < nx[2]; k++) {
					for (int j = 0; j < nx[1]; j++) {
						for (int i = 0; i < nx[0]; i++) {
							const size_t c = (size_t)i + (size_t)nx[0] * ((size_t)j + (size_t)nx[1] * (size_t)k);
							const int ijk[3] = { i, j, k };
							double x[3] = { 0.0, 0.0, 0.0 };
							for (int d = 0; d < ndim; d++)
								x[d] = xmin[d] + ((double)ijk[d] + 0.5) * dx[d];

							CellGeom g;
							cell_geometry(x, dx, g);
							if (g.r < opt.rmin || (opt.rmax > 0.0 && g.r > opt.rmax))
								continue;

							for (int a = 0; a < 3; a++)
								cell_pos.push_back(g.pos[a]);
							cell_hsml.push_back((float)std::max(g.ext[0], std::max(g.ext[1], g.ext[2])));
							// missing slab dimensions count as unit length
							double vol = 1.0;
							for (int a = 0; a < 3; a++)
								vol *= g.ext[a] > 0.0 ? g.ext[a] : 1.0;
							cell_vol.push_back((float)vol);
							cell_level.push_back((uint8_t)leaf.level);
							for (int v = 0; v < nw; v++)
								cell_vars[v].push_back((float)w[(size_t)v * ncell_block + c]);

							double rho = w[c];
							if (idx_d >= 0 && idx_lfac >= 0) {
								const double lfac = w[(size_t)idx_lfac * ncell_block + c];
								rho = w[(size_t)idx_d * ncell_block + c] / (lfac > 0.0 ? lfac : 1.0);
							}
							cell_rho.push_back((float)rho);

							if (have_b) {
								double bc[3] = { 0.0, 0.0, 0.0 };
								for (int kdir = 0; kdir < 3; kdir++) {
									const double bk = w[(size_t)idx_b[kdir] * ncell_block + c] * g.scale[kdir];
									for (int a = 0; a < 3; a++)
										bc[a] += bk * g.e[kdir][a];
								}
								for (int a = 0; a < 3; a++)
									cell_b.push_back((float)bc[a]);
							}
						}
					}
				}
			}
			ncells = cell_hsml.size();

			global_num = ncells;
#ifdef WITH_MPI
			unsigned long long l = ncells, g = 0;
			MPI_Allreduce(&l, &g, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
			global_num = g;
#endif

			if (world_rank == 0) {
				const char* coord_name = (opt.coord == Coord::MKS) ? "mks" : (opt.coord == Coord::Sph) ? "sph" : "cart";
				printf("BHAC snapshot %s: it %d, t %g, ndim %d, ndir %d, nw %d, nws %d, block %dx%dx%d, "
					"base grid %dx%dx%d blocks, levels 1..%d, %d leaf blocks\n",
					opt.dat_file.c_str(), iteration, time_code, ndim, ndir, nw, nws, nx[0], nx[1], nx[2],
					ng[0], ng[1], ng[2], levmax, nleafs);
				printf("BHAC domain (%s, coord %s", typeaxial.c_str(), coord_name);
				if (opt.coord == Coord::MKS && typeaxial == "spherical")
					printf(", h %g, R0 %g", opt.mks_h, opt.mks_r0);
				printf("):");
				for (int d = 0; d < ndim; d++)
					printf(" x%d [%g, %g]", d + 1, xprobmin[d], xprobmax[d]);
				printf("\nBHAC variables:");
				for (const std::string& n : var_names) printf(" %s", n.c_str());
				printf("\nBHAC eqpar:");
				for (double e : eqpar) printf(" %g", e);
				printf("\n");
			}
			printf("Rank %d: BHAC blocks %zu of %zu selected (%zu..%zu), leaf cells %zu\n", world_rank, count, selected.size(),
				count ? selected[first] : 0, count ? selected[first + count - 1] : 0, ncells);

#ifdef WITH_OPENMP
			steps_time[1] = omp_get_wtime();
#endif
		}

		void finish_lib() {
			var_names.clear();
			eqpar.clear();
			idx_d = idx_lfac = -1;
			idx_b[0] = idx_b[1] = idx_b[2] = -1;
			ncells = 0;
			std::vector<double>().swap(cell_pos);
			std::vector<float>().swap(cell_hsml);
			std::vector<float>().swap(cell_vol);
			std::vector<float>().swap(cell_rho);
			std::vector<float>().swap(cell_b);
			std::vector<uint8_t>().swap(cell_level);
			cell_vars.clear();
			global_num = 0;
		}

		void print_CPU_steps() {
			//printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		int get_particle_type(uint64_t id) {
			return Cell;
		}

		// Components of a block for one cell; fills v (up to 3), returns the count
		static int fetch(int blocknr, uint64_t id, float* v) {
			switch (blocknr) {
			case Pos:
				v[0] = 1.0f;
				return 1;
			case Mass:
				v[0] = (float)get_particle_mass(id);
				return 1;
			case Rho:
				v[0] = cell_rho[id];
				return 1;
			case BVec:
				if (cell_b.empty()) return 0;
				for (int a = 0; a < 3; a++)
					v[a] = cell_b[id * 3 + a];
				return 3;
			case Level:
				v[0] = (float)cell_level[id];
				return 1;
			default:
				break;
			}
			int e = blocknr - BTMax;
			if (e < 0 || e >= (int)cell_vars.size())
				return 0;
			v[0] = cell_vars[e][id];
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
			pos[0] = cell_pos[id * 3 + 0];
			pos[1] = cell_pos[id * 3 + 1];
			pos[2] = cell_pos[id * 3 + 2];
		}

		size_t get_local_num_particles() {
			return ncells;
		}

		size_t get_global_num_particles() {
			return global_num;
		}

		double get_particle_hsml(uint64_t id) {
			return cell_hsml[id];
		}

		double get_particle_mass(uint64_t id) {
			return (double)cell_rho[id] * (double)cell_vol[id];
		}

		double get_particle_rho(uint64_t id) {
			return cell_rho[id];
		}

		int get_particle_rho_blocknr() {
			return Rho;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			const int nblocks = BTMax + (int)cell_vars.size();
			types_and_blocks.assign((size_t)PTMax * nblocks, 0);
			if (global_num == 0)
				return;
			types_and_blocks[PTMax * Pos + Cell] = 1;
			types_and_blocks[PTMax * Mass + Cell] = 1;
			types_and_blocks[PTMax * Rho + Cell] = 1;
			types_and_blocks[PTMax * Level + Cell] = 1;
			if (idx_b[0] >= 0 && idx_b[1] >= 0 && idx_b[2] >= 0)
				types_and_blocks[PTMax * BVec + Cell] = 1;
			for (size_t e = 0; e < cell_vars.size(); e++)
				types_and_blocks[PTMax * (BTMax + (int)e) + Cell] = 1;
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
			case BVec: return "B";
			case Level: return "Level";
			default: break;
			}
			int e = blocknr - BTMax;
			if (e >= 0 && e < (int)var_names.size())
				return var_names[e];
			return "unknown";
		}

		double get_time() { return time_code; }
		int get_iteration() { return iteration; }

	} // namespace io
} // namespace bhac
