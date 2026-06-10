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

#include "convert_common.h"

#include <iostream>
#include <string>
#include <vector>
#include <float.h>
#include <cmath>

const std::string WHITESPACE = " \n\r\t\f\v";

namespace space_converter {
	namespace common {

		/**
		 * Calculate the magnitude (length) of a 3D vector
		 * @param x X component of the vector
		 * @param y Y component of the vector
		 * @param z Z component of the vector
		 * @return The Euclidean magnitude of the vector
		 */
		double calculate_dmagnitude3(double x, double y, double z) {
			return (float)sqrt(x * x + y * y + z * z);
		}

		/**
		 * Calculate the magnitude (length) of an N-dimensional vector (double precision)
		 * @param v Pointer to array of double values
		 * @param n Number of dimensions
		 * @return The Euclidean magnitude of the vector
		 */
		double calculate_dmagnituden(double* v, int n) {
			double res = 0;

			for (int i = 0; i < n; i++) {
				res += v[i] * v[i];
			}

			return (float)sqrt(res);
		}

		/**
		 * Calculate the magnitude (length) of an N-dimensional vector (single precision)
		 * @param v Pointer to array of float values
		 * @param n Number of dimensions
		 * @return The Euclidean magnitude of the vector
		 */
		double calculate_fmagnituden(float* v, int n) {
			double res = 0;

			for (int i = 0; i < n; i++) {
				res += v[i] * v[i];
			}

			return (float)sqrt(res);
		}

		/**
		 * Remove leading whitespace from a string
		 * @param s Input string to trim
		 * @return String with leading whitespace removed
		 */
		std::string ltrim(const std::string& s)
		{
			size_t start = s.find_first_not_of(WHITESPACE);
			return (start == std::string::npos) ? "" : s.substr(start);
		}

		/**
		 * Remove trailing whitespace from a string
		 * @param s Input string to trim
		 * @return String with trailing whitespace removed
		 */
		std::string rtrim(const std::string& s)
		{
			size_t end = s.find_last_not_of(WHITESPACE);
			return (end == std::string::npos) ? "" : s.substr(0, end + 1);
		}

		/**
		 * Remove leading and trailing whitespace from a string
		 * @param s Input string to trim
		 * @return String with both leading and trailing whitespace removed
		 */
		std::string trim(const std::string& s) {
			return rtrim(ltrim(s));
		}

		/**
		 * Gaussian kernel function for smoothing and interpolation
		 * @param x Distance from the kernel center
		 * @param h Bandwidth parameter controlling the spread of the Gaussian function
		 * @return The weighted kernel value at distance x
		 */
		double gaussian_kernel(double x, double h) {
			return std::exp(-(x * x) / (2 * h * h)) / std::sqrt(2 * 3.14159265358979323846 * h * h);
		}

		/**
		 * Generate normalized Gaussian-distributed density values
		 * The densities are normalized so their sum equals 1.0
		 * @param densities Vector to fill with normalized Gaussian values (size determines number of samples)
		 * @param h Bandwidth parameter controlling the spread of the Gaussian function
		 */
		void generate_normalized_gaussian(std::vector<double>& densities, double h) {
			// Generate N Gaussian-distributed values and compute their densities
			double total_sum = 0;
			for (int i = 0; i < densities.size(); i++) {
				densities[i] = gaussian_kernel((float)i, h);
				total_sum += densities[i];
			}

			// Normalize the densities so they sum to 1.0
			for (int i = 0; i < densities.size(); i++) {
				densities[i] /= total_sum;
			}
		}
	} // namespace common
} // namespace space_converter
