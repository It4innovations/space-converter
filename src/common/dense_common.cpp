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
#include <cstring>

namespace space_converter {
	namespace common {
		namespace vdb {
			namespace dense {

				// ─────────────────────────────────────────────────────────────────────
				// VoxelCPUDenseManager  (CPU-only concrete implementation)
				// ─────────────────────────────────────────────────────────────────────

				void VoxelCPUDenseManager::clear() {
					data_density.clear();
					data_temp.clear();
					memset(dims, 0, 3 * sizeof(size_t));
					memset(offset, 0, 3 * sizeof(size_t));
				}

				void VoxelCPUDenseManager::create(size_t x, size_t y, size_t z, bool allocate_data_temp) {
					dims[0] = x;  dims[1] = y;  dims[2] = z;

					data_density.resize(size());
					memset(data_density.data(), 0, memsize());
					if (allocate_data_temp) {
						data_temp.resize(size());
						memset(data_temp.data(), 0, memsize());
					}
				}

			} // dense

		} // vdb
	} // common
} // namespace space_converter
