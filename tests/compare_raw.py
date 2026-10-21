#!/usr/bin/env python3

#####################################################################################################################
# Copyright(C) 2023-2026 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
#
# This program is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#####################################################################################################################
"""Compare two raw float32 volume dumps with a tolerance.

Bit-identity cannot be expected across MPI rank counts or OpenMP schedules:
float accumulation order differs (atomics, partial-sum reduction), so voxel
values legitimately differ in the last ulps. This checks a relative/absolute
tolerance instead.

Usage: compare_raw.py A.raw B.raw [rel_tol] [abs_tol] [max_bad_fraction]
max_bad_fraction (default 0) allows a small share of out-of-tolerance voxels —
useful for k-NN comparisons where equal-distance ties legitimately pick
different neighbor sets per backend or rank split.
Exit 0 when equal within tolerance, 1 otherwise.
"""

import struct
import sys


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    rel = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-4
    abs_tol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6
    max_bad_fraction = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0

    with open(a_path, "rb") as f:
        a = f.read()
    with open(b_path, "rb") as f:
        b = f.read()
    if len(a) != len(b):
        print(f"size mismatch: {len(a)} vs {len(b)} bytes")
        return 1

    n = len(a) // 4
    va = struct.unpack(f"<{n}f", a)
    vb = struct.unpack(f"<{n}f", b)

    worst = 0.0
    worst_i = -1
    bad = 0
    suma = 0.0
    sumb = 0.0
    for i in range(n):
        x, y = va[i], vb[i]
        suma += x
        sumb += y
        d = abs(x - y)
        if d > abs_tol + rel * max(abs(x), abs(y)):
            bad += 1
            if d > worst:
                worst, worst_i = d, i
    print(f"voxels: {n}, out-of-tolerance: {bad} ({100.0 * bad / n:.4f}%), "
          f"sum(A)={suma:.6g}, sum(B)={sumb:.6g}")
    if bad:
        print(f"worst diff {worst:g} at voxel {worst_i}: {va[worst_i]:g} vs {vb[worst_i]:g}")
    if bad > max_bad_fraction * n:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
