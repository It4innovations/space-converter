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

#pragma once

//#include "convert_vdb.h"

#include <vector>

namespace common {
	namespace vdb {
		// /**
		//  * Structure for storing particle simulation parameters
		//  * Used in density projection and rendering algorithms
		//  */
		// struct particle_sim_v2
		// {
		// 	float e;          // Attribute value (e.g., temperature, energy)
		// 	float x, y, z;    // Particle position
		// 	float r;          // Particle radius/smoothing length
		// 	float I;          // Particle intensity/mass
		// };

		// double custom_asinh(double val);
		// //float custom_xexp(float v);

		// /**
		//  * Template class for fast exponential function evaluation using lookup tables
		//  * Provides significant speedup for repeated exponential calculations
		//  * @tparam T Floating point type (float or double)
		//  */
        // template<typename T> 
        // class exptable
        // {
        // private:
        //     T expfac, taylorlimit;
        //     std::vector<T> tab1, tab2;  // Two-level lookup table for exponential values
        //     enum {
        //         nbits = 10,                    // Number of bits for table indexing
        //         dim1 = 1 << nbits,             // First dimension size (1024)
        //         mask1 = dim1 - 1,              // Mask for first level (1023)
        //         dim2 = (1 << nbits) << nbits,  // Second dimension size
        //         mask2 = dim2 - 1,              // Mask for second level
        //         mask3 = ~mask2                 // Mask for out-of-range check
        //     };

        // public:
        //     /**
        //      * Constructor: Build exponential lookup tables
        //      * @param maxexp Maximum exponent value to support
        //      */
        //     exptable(T maxexp)
        //         : expfac(dim2 / maxexp)
        //     {
        //         tab1.resize(dim1);                    
        //         tab2.resize(dim1);

        //         using namespace std;
        //         // Pre-compute exponential values for two-level lookup
        //         for (int m = 0; m < dim1; ++m)
        //         {
        //             tab1[m] = exp(m * dim1 / expfac);  // First level: larger steps
        //             tab2[m] = exp(m / expfac);          // Second level: finer steps
        //         }
        //         // Compute threshold for Taylor series approximation
        //         taylorlimit = sqrt(T(2) * abs(maxexp) / dim2);
        //     }

        //     /**
        //      * Get Taylor series threshold for small arguments
        //      * For |arg| < taylorLimit(), use Taylor series instead of table lookup
        //      * @return Threshold value for Taylor series approximation
        //      */
        //     T taylorLimit() const { return taylorlimit; }

        //     /**
        //      * Compute exp(arg) using lookup table
        //      * @param arg Exponent value
        //      * @return Approximate value of exp(arg)
        //      */
        //     T operator() (T arg) const
        //     {
        //         int iarg = int(arg * expfac);
        //         if (iarg & mask3)  // Out of range check
        //             return (iarg < 0) ? T(1) : T(0);
        //         // Two-level lookup: multiply coarse and fine values
        //         return tab1[iarg >> nbits] * tab2[iarg & mask1];
        //     }

        //     /**
        //      * Compute exp(arg) - 1 with better precision for small arguments
        //      * For small |arg|, uses Taylor series: exp(x) - 1 ≈ x
        //      * @param arg Exponent value
        //      * @return Approximate value of exp(arg) - 1
        //      */
        //     T expm1(T arg) const
        //     {
        //         if (std::abs(arg) < taylorlimit) 
        //             return arg;  // Taylor series: exp(x) - 1 ≈ x for small x

        //         return operator()(arg) - T(1);
        //     }
        // };

		//class VoxelDenseManager : public VoxelManager {

			/**
			 * @brief Dense regular grid representation of particle data
			 * 
			 * Stores particle data rasterized onto a uniform 3D grid.
			 * Used for volume rendering and analysis.
			 */
			struct DenseParticles
			{
				std::vector<float> data_density;   ///< Primary density/value data
	#ifndef WITH_NO_DATA_TEMP
				std::vector<float> data_temp;      ///< Temporary accumulation buffer for normalization
	#endif
				size_t dims[3] = { 0,0,0 };        ///< Grid dimensions [x, y, z]
				size_t offset[3] = { 0,0,0 };      ///< Grid offset in global coordinate space

				/**
				 * @brief Clear all grid data and reset dimensions
				 */
				void clear() {
					data_density.clear();
	#ifndef WITH_NO_DATA_TEMP
					data_temp.clear();
	#endif
					memset(dims, 0, 3 * sizeof(size_t));
					memset(offset, 0, 3 * sizeof(size_t));
				}

				/**
				 * @brief Create and initialize grid with specified dimensions
				 * @param x Width of the grid
				 * @param y Height of the grid
				 * @param z Depth of the grid
				 */
				void create(size_t x, size_t y, size_t z) {
					dims[0] = x;
					dims[1] = y;
					dims[2] = z;

					data_density.resize(size());
					memset(data_density.data(), 0, memsize());
	#ifndef WITH_NO_DATA_TEMP
					data_temp.resize(size());
					memset(data_temp.data(), 0, memsize());
	#endif
				}

				/** @brief Get grid width */
				size_t x() {
					return dims[0];
				}

				/** @brief Get grid height */
				size_t y() {
					return dims[1];
				}

				/** @brief Get grid depth */
				size_t z() {
					return dims[2];
				}

				/** @brief Get total number of voxels */
				size_t size() {
					return dims[0] * dims[1] * dims[2];
				}

				/** @brief Get total memory size in bytes */
				size_t memsize() {
					return dims[0] * dims[1] * dims[2] * sizeof(float);
				}

				/**
				 * @brief Convert 3D coordinates to linear index
				 * @param x X coordinate
				 * @param y Y coordinate
				 * @param z Z coordinate
				 * @return Linear array index
				 */
				size_t get_index(size_t x, size_t y, size_t z) {
					return x + y * dims[0] + z * dims[0] * dims[1];
				}
			};
		//};     
	}// vdb
} //common