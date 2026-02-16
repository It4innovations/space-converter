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

#include <map>
#include <string>
#include <vector>

#ifdef _WIN32
#	define STDMAX max
#	define STDMIN min
#else
#	define STDMAX std::max
#	define STDMIN std::min
#endif

namespace common {

	// String manipulation utilities
	std::string ltrim(const std::string& s);        // Remove leading whitespace
	std::string rtrim(const std::string& s);        // Remove trailing whitespace
	std::string trim(const std::string& s);         // Remove leading and trailing whitespace

	// Gaussian kernel functions for smoothing and density estimation
	double gaussian_kernel(double x, double h = 1.0);                            // Compute Gaussian kernel value
	void generate_normalized_gaussian(std::vector<double>& densities, double h = 1.0);  // Generate normalized Gaussian distribution

	// Vector magnitude calculation utilities
	double calculate_dmagnitude3(double x, double y, double z);  // Calculate 3D vector magnitude
	double calculate_dmagnituden(double* v, int n);              // Calculate N-D vector magnitude (double)
	double calculate_fmagnituden(float* v, int n);               // Calculate N-D vector magnitude (float)
} //common