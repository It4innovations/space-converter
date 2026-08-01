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
"""Cut a small, still-valid NChilada snapshot out of a large one.

Full ChaNGa snapshots do not fit in memory on a handful of nodes (the reader keeps
every block resident), which makes them impractical for testing. This copies the
first N particles of every field file, keeping the on-disk layout intact:

    FieldHeader (28 bytes, XDR big-endian) | min | max | numParticles * dimensions values

min/max are recomputed over the particles that are kept, so the subset stays
self-consistent. Fields whose stored min equals max hold no data array at all
(the reader expands the single value) and are copied through unchanged.

Usage: nchilada_subset.py <src-dir> <dst-dir> [num-particles]
"""

import os
import shutil
import struct
import sys

HEADER_FMT = ">i d Q I i"          # magic, time, numParticles, dimensions, code
HEADER_SIZE = struct.calcsize(HEADER_FMT)
MAGIC = 1062053

# TypeHandling::DataTypeCode -> (struct char, size)
CODES = {
    1: ("b", 1), 2: ("B", 1), 3: ("h", 2), 4: ("H", 2), 5: ("i", 4),
    6: ("I", 4), 7: ("q", 8), 8: ("Q", 8), 9: ("f", 4), 10: ("d", 8),
}


def subset_field(src, dst, want):
    with open(src, "rb") as f:
        raw = f.read(HEADER_SIZE)
        if len(raw) < HEADER_SIZE:
            return None
        magic, time, nparts, dims, code = struct.unpack(HEADER_FMT, raw)
        if magic != MAGIC:
            return None
        if code not in CODES:
            return None

        ch, esize = CODES[code]
        keep = min(want, nparts)
        rec = ch * dims                       # one particle = dims values

        lo = f.read(esize * dims)
        hi = f.read(esize * dims)
        if len(lo) < esize * dims or len(hi) < esize * dims:
            return None

        # min == max marks a constant field stored without a data array
        if lo == hi:
            with open(dst, "wb") as g:
                g.write(struct.pack(HEADER_FMT, magic, time, keep, dims, code))
                g.write(lo)
                g.write(hi)
            return keep

        data = f.read(esize * dims * keep)
        if len(data) < esize * dims * keep:
            return None

        vals = struct.unpack(">" + rec * keep, data)
        per_dim_min = [min(vals[d::dims]) for d in range(dims)]
        per_dim_max = [max(vals[d::dims]) for d in range(dims)]

        with open(dst, "wb") as g:
            g.write(struct.pack(HEADER_FMT, magic, time, keep, dims, code))
            g.write(struct.pack(">" + rec, *per_dim_min))
            g.write(struct.pack(">" + rec, *per_dim_max))
            g.write(data)
        return keep


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1

    src_root, dst_root = sys.argv[1], sys.argv[2]
    want = int(sys.argv[3]) if len(sys.argv) > 3 else 100000

    os.makedirs(dst_root, exist_ok=True)
    for name in os.listdir(src_root):
        src = os.path.join(src_root, name)
        dst = os.path.join(dst_root, name)
        if os.path.isfile(src):
            shutil.copy2(src, dst)          # description.xml and friends
            continue
        os.makedirs(dst, exist_ok=True)
        for field in sorted(os.listdir(src)):
            sf, df = os.path.join(src, field), os.path.join(dst, field)
            if not os.path.isfile(sf):
                continue
            kept = subset_field(sf, df, want)
            print(f"  {name}/{field}: {'skipped (not a field file)' if kept is None else str(kept) + ' particles'}")
            if kept is None and os.path.exists(df):
                os.remove(df)
    return 0


if __name__ == "__main__":
    sys.exit(main())
