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

/*
 * Shared RETURN_* helper macros for the per-solver readers.
 *
 * Every *_extract_iolib.cpp implements three accessor families:
 *   get_particle_norm_value(blocknr, id)  -> one float per particle ("norm" value:
 *                                            scalars pass through, vectors collapse
 *                                            to their Euclidean magnitude)
 *   get_particle_value(blocknr, id, out)  -> original components, returns count
 *   get_particle_value_comp(blocknr, id)  -> component count only
 *
 * These macros used to be copy-pasted (and slowly diverging) in each reader;
 * they are defined once here.  Note: "norm" does NOT mean normalized to [0,1] —
 * it means scalar-reduced (magnitude for vectors).
 */

#include "convert_common.h"

// -- get_particle_norm_value helpers ----------------------------------------
#define RETURN_NORM_EMPTY return 0;
#define RETURN_NORM_VALUE(v) return (float)(v);
#define RETURN_NORM_VECTOR3(v) return (float)space_converter::common::calculate_dmagnitude3(v[0], v[1], v[2]);
#define RETURN_NORM_DVECTORN(v, n) return (float)space_converter::common::calculate_dmagnituden(v, n);
#define RETURN_NORM_FVECTORN(v, n) return (float)space_converter::common::calculate_fmagnituden(v, n);

// -- get_particle_value_comp helpers ----------------------------------------
#define RETURN_COMP_EMPTY return 0;
#define RETURN_COMP_VALUE(v) return 1;
#define RETURN_COMP_VECTOR3(v) return 3;
#define RETURN_COMP_DVECTORN(v, n) return n;
#define RETURN_COMP_FVECTORN(v, n) return n;

// -- get_particle_value helpers (fill out_value[], return component count) --
#define RETURN_ORIG_EMPTY RETURN_COMP_EMPTY
#define RETURN_ORIG_VALUE(v) out_value[0] = (float)v; RETURN_COMP_VALUE(v)
#define RETURN_ORIG_VECTOR3(v) out_value[0] = (float)v[0]; out_value[1] = (float)v[1]; out_value[2] = (float)v[2]; RETURN_COMP_VECTOR3(v)
#define RETURN_ORIG_DVECTORN(v, n) for (int iv = 0; iv < n; iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_DVECTORN(v, n)
#define RETURN_ORIG_FVECTORN(v, n) for (int iv = 0; iv < n; iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_FVECTORN(v, n)
