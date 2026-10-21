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

namespace space_converter {

    /**
     * @brief Initialize logging system based on command-line configuration.
     * 
     * This function sets up the logging system, including log levels, output
     * destinations (console, file), and formatting. It can be extended to support
     * more advanced features like log rotation or remote logging.
     * 
     * @param world_rank Rank of the MPI process (used to differentiate logs in parallel runs)
     */
    void LOG_Init(int world_rank);

    /**
     * @brief Start a named performance measurement section.
     * 
     * This function marks the beginning of a code section for performance measurement.
     * It can be used in conjunction with LOG_MeasureStop() to measure execution time
     * and other metrics for specific code blocks.
     * 
     * @param name Name of the measurement section (e.g., "Data Loading", "Grid Conversion")
     */
    void LOG_MeasureStart(const char* name);

    /**
     * @brief Stop a named performance measurement section and log results.
     *
     * This function marks the end of a code section for performance measurement.
     * It calculates the elapsed time and logs it along with the section name.
     *
     * @param name Name of the measurement section (should match the name used in LOG_MeasureStart)
     */
    void LOG_MeasureStop(const char* name);

    /**
     * @brief Whether verbose diagnostic output is enabled.
     *
     * Controlled by the SPACE_CONVERTER_VERBOSE environment variable (any
     * non-empty value other than "0" enables it). Per-chunk/per-merge timing
     * printfs are guarded with this so large MPI runs stay quiet by default.
     * The value is read once and cached.
     */
    bool LOG_Verbose();

} // namespace space_converter