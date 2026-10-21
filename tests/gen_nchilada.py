#!/usr/bin/env python3
"""Generate a small synthetic CHANGA_NCHILADA dataset for smoke tests.

Layout written (see src/changa/changa_nchilada_extract_iolib.cpp init_lib()):
  <dir>/gas/{pos,mass,vel,phi,smoothlength,GasDensity,temperature,OxMassFrac,FeMassFrac}
  <dir>/dark/{pos,mass,vel,soft,phi}
  <dir>/star/{pos,mass,vel,soft,phi,OxMassFrac,FeMassFrac,timeform}   (0 particles)

Each file is an NChilada "field" file (tree_xdr.h, XDR = big-endian):
  FieldHeader (28 bytes): int magic (1062053), double time, uint64 numParticles,
                          uint32 dimensions (1 or 3), int32 code (float32 = 9)
  then minValue, maxValue (dimensions values each), then numParticles*dimensions
  values. If min == max the reader fills the array from min without reading the
  data section, so constant fields still round-trip correctly.

Contents:
  - gas:  600 particles on a jittered lattice in a NON-CUBIC box, mass 1.0,
          varying GasDensity/temperature/phi/metals, smoothlength 2.0
  - dark: 300 particles, mass 2.5 (unequal to gas on purpose), soft 1.5
  - star: every requested file present but empty (avoids "block not found"
          warnings and exercises the 0-particle path)
  All per-family files share the same numParticles (init_lib drops blocks
  whose count differs from the family's pos block).

Deterministic (fixed seed): re-running produces byte-identical files.
"""

import os
import struct
import sys
import random

NGAS = 600
NDARK = 300
BOX = (100.0, 200.0, 400.0)   # non-cubic on purpose (exercises bbox symmetrization)
SEED = 20260828
TIME = 0.5

FIELD_MAGIC = 1062053
CODE_FLOAT32 = 9              # TypeHandling::DataTypeCode::float32


def field_file(values, dimensions):
    """Serialize one XDR field file: header, min, max, data (all big-endian float32).

    values: flat list of floats, length = numParticles * dimensions.
    """
    assert len(values) % dimensions == 0
    n = len(values) // dimensions
    if n > 0:
        vmin = [min(values[d::dimensions]) for d in range(dimensions)]
        vmax = [max(values[d::dimensions]) for d in range(dimensions)]
    else:
        vmin = [0.0] * dimensions   # empty file still carries min/max
        vmax = [0.0] * dimensions
    out = struct.pack(">idQIi", FIELD_MAGIC, TIME, n, dimensions, CODE_FLOAT32)
    out += struct.pack(">%df" % dimensions, *vmin)
    out += struct.pack(">%df" % dimensions, *vmax)
    out += struct.pack(">%df" % len(values), *values)
    return out


def write_field(basedir, family, name, values, dimensions):
    fam_dir = os.path.join(basedir, family)
    os.makedirs(fam_dir, exist_ok=True)
    with open(os.path.join(fam_dir, name), "wb") as f:
        f.write(field_file(values, dimensions))


def jittered_positions(n, rng):
    """n points on a jittered lattice inside BOX (flat [x,y,z,x,y,z,...] list)."""
    side = max(1, int(round(n ** (1.0 / 3.0))))
    flat = []
    i = 0
    while len(flat) < 3 * n:
        ix, iy, iz = i % side, (i // side) % side, i // (side * side)
        flat.append((ix + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[0])
        flat.append((iy + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[1])
        flat.append((iz + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[2])
        i += 1
    return flat


def velocities(n, rng):
    return [rng.uniform(-1.0, 1.0) for _ in range(3 * n)]


def main(out_dir):
    rng = random.Random(SEED)

    # --- gas family ---------------------------------------------------------
    write_field(out_dir, "gas", "pos", jittered_positions(NGAS, rng), 3)
    write_field(out_dir, "gas", "vel", velocities(NGAS, rng), 3)
    # Constant mass 1.0 (min == max shortcut in the reader)
    write_field(out_dir, "gas", "mass", [1.0] * NGAS, 1)
    write_field(out_dir, "gas", "phi", [-1.0 - i * 0.001 for i in range(NGAS)], 1)
    write_field(out_dir, "gas", "smoothlength", [2.0] * NGAS, 1)
    # Same easily assertable ranges as gen_gadget_simple.py
    write_field(out_dir, "gas", "GasDensity", [0.5 + (i % 100) * 0.01 for i in range(NGAS)], 1)
    write_field(out_dir, "gas", "temperature", [10.0 + i * 0.01 for i in range(NGAS)], 1)
    write_field(out_dir, "gas", "OxMassFrac", [0.01 + (i % 10) * 0.001 for i in range(NGAS)], 1)
    write_field(out_dir, "gas", "FeMassFrac", [0.002 + (i % 10) * 0.0005 for i in range(NGAS)], 1)

    # --- dark family --------------------------------------------------------
    write_field(out_dir, "dark", "pos", jittered_positions(NDARK, rng), 3)
    write_field(out_dir, "dark", "vel", velocities(NDARK, rng), 3)
    # Unequal masses per family: gas 1.0, dark 2.5 (mass-weighted tests rely on this)
    write_field(out_dir, "dark", "mass", [2.5] * NDARK, 1)
    write_field(out_dir, "dark", "soft", [1.5] * NDARK, 1)
    write_field(out_dir, "dark", "phi", [-2.0 - i * 0.001 for i in range(NDARK)], 1)

    # --- star family: requested by init_lib but empty in this dataset -------
    for name, dims in (("pos", 3), ("vel", 3), ("mass", 1), ("soft", 1),
                       ("phi", 1), ("OxMassFrac", 1), ("FeMassFrac", 1),
                       ("timeform", 1)):
        write_field(out_dir, "star", name, [], dims)

    print(f"wrote {out_dir}: {NGAS} gas + {NDARK} dark + 0 star particles, box {BOX}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "nchilada_test")
