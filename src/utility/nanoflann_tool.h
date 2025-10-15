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

#pragma once

#include <cmath>
#include <vector>
#include "data_common.h"

namespace utility {
	namespace nanoflann_tool {

        void run_knn(float* points, size_t N, int k, std::vector<float>& radius_particles, std::vector<float>& rho_particles, std::vector<float>& mass_particles, bool use_cycling, common::SpaceData::DenseType& rho_kernel);

	}// nanoflann
} //utility