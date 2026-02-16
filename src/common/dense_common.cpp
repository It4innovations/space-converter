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

#include "dense_common.h"

#ifdef WITH_OPENMP
# include <omp.h>
#endif

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "utility/dense_utility.h"

namespace common {
	namespace vdb {
		void ConvertVDBBase::fill_voxels_v5(
			common::vdb::DenseParticles& grid,
			size_t pid,
			float value,
			int bbox_dim,
			int* bbox_min_orig,
			double bbox_size_orig,
			double scale_space_diagonal,
			common::SpaceData::DenseType dense_type,
			common::SpaceData::DenseNorm dense_norm,
			float particle_fix_size,
			int particle_type,
			int block_name_id,
			double* pos
		) {

			double norm_fac = (double)bbox_dim / scale_space_diagonal;
			double len2pix = norm_fac / bbox_size_orig;

			if (block_name_id == 0) { //only pos
				value = 1.0f;
			}

			double hsml = get_particle_radius(
				pid,
				bbox_dim,
				bbox_min_orig,
				bbox_size_orig,
				scale_space_diagonal,
				particle_fix_size,
				particle_type
			);

			int iradiusx = static_cast<int>(hsml);
			int iradiusy = static_cast<int>(hsml);
			int iradiusz = static_cast<int>(hsml);

			///////////////////////////////////////////////////////////////////////////////
			double dpx = ((double)pos[0] - (double)bbox_min_orig[0]) * len2pix;
			double dpy = ((double)pos[1] - (double)bbox_min_orig[1]) * len2pix;
			double dpz = ((double)pos[2] - (double)bbox_min_orig[2]) * len2pix;

			int px = static_cast<int>(dpx);
			int py = static_cast<int>(dpy);
			int pz = static_cast<int>(dpz);

			///////////////////////////////////////////////////////////////////////////////

			for (int sx = px - iradiusx; sx <= px + iradiusx; sx++) {
				for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
					for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

						int osx = sx - grid.offset[0];
						int osy = sy - grid.offset[1];
						int osz = sz - grid.offset[2];

						if (osx < 0 || osy < 0 || osz < 0)
							continue;

						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
							continue;

						double density = 1.0;

						double dx = dpx - sx;
						double dy = dpy - sy;
						double dz = dpz - sz;
						double distance_norm = std::sqrt(dx * dx + dy * dy + dz * dz);

						double W = 0.0;
						double h = hsml;

						if (iradiusx == 0 && iradiusy == 0 && iradiusz == 0) {
							W = 1.0; // full value
						}
						else {
							double q = distance_norm / h;

							W = utility::dense::sph_kernel::W(dense_type, q, 1.0 / h);
						}

						density = W;

						//final density
						float d = density * value;

						double norm = 0.0;
						if (dense_norm == common::SpaceData::DenseNorm::eNone) {
							norm = 0.0;
						}
						else if (dense_norm == common::SpaceData::DenseNorm::eCount) {
							norm = 1.0; // count
						}
						else if (dense_norm == common::SpaceData::DenseNorm::eSPHInterpolation) {
							norm = density;
						}

						//final norm
						float n = norm;

						size_t gindex = grid.get_index(osx, osy, osz);

#pragma omp atomic
						grid.data_density[gindex] += d;

#ifndef WITH_NO_DATA_TEMP
#pragma omp atomic
						grid.data_temp[gindex] += n;
#endif
					}
				}
			}
		}

	}// dense
} //common