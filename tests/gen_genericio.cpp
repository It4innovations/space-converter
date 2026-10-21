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

/*
 * Synthetic HACC_GENERICIO dataset writer for the smoke tests.
 *
 * Uses the bundled GenericIO library (MPI build — the writer path is only
 * compiled with MPI enabled) to produce a small deterministic dataset the
 * HACC_GENERICIO reader can consume:
 *
 *   gen_genericio <output-file>
 *
 * 900 particles on a jittered lattice inside a NON-CUBIC box (64 x 128 x 256),
 * variables: x, y, z, vx, vy, vz, mass, rho, hh (float32) and id (int64).
 * Matching converter flags:
 *   --pos-names x y z --vel-names vx vy vz
 *   --mass-name mass --rho-name rho --hsml-name hh
 *
 * Deterministic: a tiny LCG instead of rand() so re-runs are byte-identical.
 * Run with a single MPI rank (one GenericIO rank in the file).
 */

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include <mpi.h>

#include "GenericIO.h"

namespace {

	const size_t kNumParticles = 900;
	const double kBox[3] = { 64.0, 128.0, 256.0 };  // Mpc/h, non-cubic on purpose

	// Deterministic little LCG (do not use rand(): implementation-defined)
	struct Lcg {
		uint64_t state = 20260831u;
		double uniform(double lo, double hi) {
			state = state * 6364136223846793005ULL + 1442695040888963407ULL;
			double u = (double)(state >> 11) / (double)(1ULL << 53);
			return lo + u * (hi - lo);
		}
	};

	template <typename T>
	std::vector<T> with_extra_space(size_t n) {
		return std::vector<T>(n + gio::GenericIO::requestedExtraSpace() / sizeof(T));
	}

} // namespace

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

	int world_size = 1;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	if (world_size != 1) {
		fprintf(stderr, "gen_genericio: run with exactly 1 MPI rank\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	std::string out_path = (argc > 1) ? argv[1] : "genericio_test";

	Lcg rng;
	const size_t n = kNumParticles;
	int side = 1;
	while ((size_t)side * side * side < n)
		side++;

	std::vector<float> x = with_extra_space<float>(n);
	std::vector<float> y = with_extra_space<float>(n);
	std::vector<float> z = with_extra_space<float>(n);
	std::vector<float> vx = with_extra_space<float>(n);
	std::vector<float> vy = with_extra_space<float>(n);
	std::vector<float> vz = with_extra_space<float>(n);
	std::vector<float> mass = with_extra_space<float>(n);
	std::vector<float> rho = with_extra_space<float>(n);
	std::vector<float> hh = with_extra_space<float>(n);
	std::vector<int64_t> id = with_extra_space<int64_t>(n);

	for (size_t i = 0; i < n; i++) {
		int ix = (int)(i % side);
		int iy = (int)((i / side) % side);
		int iz = (int)(i / ((size_t)side * side));
		x[i] = (float)((ix + 0.5 + rng.uniform(-0.3, 0.3)) / side * kBox[0]);
		y[i] = (float)((iy + 0.5 + rng.uniform(-0.3, 0.3)) / side * kBox[1]);
		z[i] = (float)((iz + 0.5 + rng.uniform(-0.3, 0.3)) / side * kBox[2]);
		vx[i] = (float)rng.uniform(-100.0, 100.0);
		vy[i] = (float)rng.uniform(-100.0, 100.0);
		vz[i] = (float)rng.uniform(-100.0, 100.0);
		mass[i] = 1.0f;
		rho[i] = 0.5f + (float)(i % 100) * 0.01f;
		hh[i] = 2.0f;
		id[i] = (int64_t)(i + 1);
	}

	{
		gio::GenericIO gio_out(MPI_COMM_WORLD, out_path);
		gio_out.setNumElems(n);
		for (int d = 0; d < 3; d++) {
			gio_out.setPhysOrigin(0.0, d);
			gio_out.setPhysScale(kBox[d], d);
		}

		const unsigned extra = gio::GenericIO::VarHasExtraSpace;
		gio_out.addVariable("x", x, extra | gio::GenericIO::VarIsPhysCoordX);
		gio_out.addVariable("y", y, extra | gio::GenericIO::VarIsPhysCoordY);
		gio_out.addVariable("z", z, extra | gio::GenericIO::VarIsPhysCoordZ);
		gio_out.addVariable("vx", vx, extra);
		gio_out.addVariable("vy", vy, extra);
		gio_out.addVariable("vz", vz, extra);
		gio_out.addVariable("mass", mass, extra);
		gio_out.addVariable("rho", rho, extra);
		gio_out.addVariable("hh", hh, extra);
		gio_out.addVariable("id", id, extra);

		gio_out.write();
	}

	printf("wrote %s: %zu particles, box (%g, %g, %g)\n",
		out_path.c_str(), n, kBox[0], kBox[1], kBox[2]);

	MPI_Finalize();
	return 0;
}
