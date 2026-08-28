#!/usr/bin/env python3
"""Generate a small synthetic CHANGA_TIPSY snapshot (standard/XDR, big-endian).

Layout (matches src/changa/utility/structures/TipsyReader.cpp):
  header (32 B): double time, uint nbodies, int ndim=3, uint nsph, uint ndark,
                 uint nstar, 4-byte pad — all big-endian (XDR)
  gas particles  (nsph):  mass, pos[3], vel[3], rho, temp, hsmooth, metals, phi
                          -> 12 big-endian float32 (48 B)
  dark particles (ndark): mass, pos[3], vel[3], eps, phi -> 9 floats (36 B)
  star particles (nstar): mass, pos[3], vel[3], metals, tform, eps, phi
                          -> 11 floats (44 B)

Deterministic; 700 gas + 350 dark + 0 stars on jittered lattices in a
non-cubic box, unequal masses per family (gas 1.0, dark 2.5).
"""

import struct
import sys
import random

NSPH = 700
NDARK = 350
NSTAR = 0
BOX = (10.0, 20.0, 40.0)   # non-cubic on purpose
SEED = 20260829


def jittered_positions(n, rng):
    side = max(1, int(round(n ** (1.0 / 3.0))))
    pts = []
    i = 0
    while len(pts) < n:
        ix, iy, iz = i % side, (i // side) % side, i // (side * side)
        # Tipsy positions are conventionally centered around 0
        x = ((ix + 0.5 + rng.uniform(-0.3, 0.3)) / side - 0.5) * BOX[0]
        y = ((iy + 0.5 + rng.uniform(-0.3, 0.3)) / side - 0.5) * BOX[1]
        z = ((iz + 0.5 + rng.uniform(-0.3, 0.3)) / side - 0.5) * BOX[2]
        pts.append((x, y, z))
        i += 1
    return pts


def main(out_path):
    rng = random.Random(SEED)
    nbodies = NSPH + NDARK + NSTAR

    data = struct.pack(">d5ii", 0.25, nbodies, 3, NSPH, NDARK, NSTAR, 0)
    # ">d5ii": time, nbodies, ndim, nsph, ndark, nstar, pad -> 8 + 6*4 = 32 B

    for p in jittered_positions(NSPH, rng):
        vel = (rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1))
        rho = 0.5 + rng.uniform(0.0, 1.0)
        temp = 1e4 + rng.uniform(0.0, 1e3)
        hsmooth = 0.4
        metals = 0.02
        phi = -1.0
        data += struct.pack(">12f", 1.0, *p, *vel, rho, temp, hsmooth, metals, phi)

    for p in jittered_positions(NDARK, rng):
        vel = (rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1))
        eps = 0.6
        phi = -1.0
        data += struct.pack(">9f", 2.5, *p, *vel, eps, phi)

    with open(out_path, "wb") as f:
        f.write(data)
    print(f"wrote {out_path}: {NSPH} gas + {NDARK} dark + {NSTAR} star, box {BOX}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tipsy_test")
