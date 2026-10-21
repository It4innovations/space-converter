#!/usr/bin/env python3
"""Generate a small synthetic HACC_BIN file.

Record layout (src/hacc/haccbin_extract_iolib.cpp, #pragma pack(1), little-endian,
58 bytes/particle):
  float x, vx, y, vy, z, vz, mass, uu, hh, mu, rho, phi   (12 * 4 B)
  int64 id                                                (8 B)
  uint16 mask                                             (2 B)

Deterministic: 800 particles on a jittered lattice in a non-cubic box.
"""

import struct
import sys
import random

N = 800
BOX = (64.0, 128.0, 256.0)   # Mpc/h, non-cubic on purpose
SEED = 20260830


def main(out_path):
    rng = random.Random(SEED)
    side = max(1, int(round(N ** (1.0 / 3.0))))
    rec = struct.Struct("<12fqH")
    data = bytearray()
    for i in range(N):
        ix, iy, iz = i % side, (i // side) % side, i // (side * side)
        x = (ix + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[0]
        y = (iy + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[1]
        z = (iz + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[2]
        vx, vy, vz = rng.uniform(-100, 100), rng.uniform(-100, 100), rng.uniform(-100, 100)
        mass = 1.0
        uu = 100.0 + i * 0.1
        hh = 2.0                       # smoothing length (Mpc/h)
        mu = 0.6
        rho = 0.5 + (i % 100) * 0.01
        phi = -1.0
        data += rec.pack(x, vx, y, vy, z, vz, mass, uu, hh, mu, rho, phi, i + 1, 0)

    with open(out_path, "wb") as f:
        f.write(data)
    print(f"wrote {out_path}: {N} particles, {rec.size} B/record, box {BOX}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "haccbin_test")
