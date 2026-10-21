#!/usr/bin/env python3
"""Generate a small synthetic GADGET_SIMPLE (format-2) snapshot for smoke tests.

Layout written (little-endian, format 2):
  For every block: [int 8]["XXXX"][int payload+8][int 8] [int payload][payload][int payload]
  HEAD payload is the 256-byte io_header (see OpenGadget3 allvars.h).

Contents:
  - 2 particle types: type 0 ("gas", N0 particles) and type 1 ("halo", N1),
    on jittered lattices inside a deliberately NON-CUBIC box.
  - POS/VEL (float3, all types), ID (uint32, all types),
    MASS (float, all types; header mass table zeroed -> per-particle masses,
    UNEQUAL between the two types), U/RHO/HSML (float, gas only).

Deterministic (fixed seed): re-running produces byte-identical files.
"""

import struct
import sys
import random

N0 = 1000   # gas
N1 = 500    # halo / dark matter
BOX = (100.0, 200.0, 400.0)   # non-cubic on purpose (exercises bbox symmetrization)
SEED = 20260828


def block(name, payload):
    assert len(name) == 4
    head = struct.pack("<i4si i", 8, name.encode(), len(payload) + 8, 8)
    rec = struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))
    return head + rec


def make_header():
    npart = [N0, N1, 0, 0, 0, 0]
    mass = [0.0] * 6          # zero -> masses come from the MASS block
    header = b""
    header += struct.pack("<6i", *npart)
    header += struct.pack("<6d", *mass)
    header += struct.pack("<d", 0.5)          # time
    header += struct.pack("<d", 1.0)          # redshift
    header += struct.pack("<i", 0)            # flag_sfr
    header += struct.pack("<i", 0)            # flag_feedback
    header += struct.pack("<6I", *npart)      # npartTotal
    header += struct.pack("<i", 0)            # flag_cooling
    header += struct.pack("<i", 1)            # num_files
    header += struct.pack("<d", max(BOX))     # BoxSize
    header += struct.pack("<d", 0.3)          # Omega0
    header += struct.pack("<d", 0.7)          # OmegaLambda
    header += struct.pack("<d", 0.7)          # HubbleParam
    header += struct.pack("<i", 0)            # flag_stellarage
    header += struct.pack("<i", 0)            # flag_metals
    header += struct.pack("<6I", 0, 0, 0, 0, 0, 0)  # npartTotalHighWord
    header += struct.pack("<i", 0)            # flag_entropy_instead_u
    header += struct.pack("<i", 0)            # flag_doubleprecision
    header += struct.pack("<i", 0)            # flag_ic_info
    header += struct.pack("<f", 0.0)          # lpt_scalingfactor
    header += b"\0" * 18                      # fill
    header += b"\0" * 30                      # names[15][2]
    assert len(header) == 256, len(header)
    return header


def jittered_positions(n, rng):
    """n points on a jittered lattice inside BOX."""
    import math
    side = max(1, int(round(n ** (1.0 / 3.0))))
    pts = []
    i = 0
    while len(pts) < n:
        ix, iy, iz = i % side, (i // side) % side, i // (side * side)
        x = (ix + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[0]
        y = (iy + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[1]
        z = (iz + 0.5 + rng.uniform(-0.3, 0.3)) / side * BOX[2]
        pts.append((x, y, z))
        i += 1
    return pts


def main(out_path):
    rng = random.Random(SEED)
    pos = jittered_positions(N0, rng) + jittered_positions(N1, rng)
    ntot = N0 + N1

    pos_payload = b"".join(struct.pack("<3f", *p) for p in pos)
    vel_payload = b"".join(struct.pack("<3f", rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(ntot))
    id_payload = b"".join(struct.pack("<I", i + 1) for i in range(ntot))
    # Unequal masses per type: gas 1.0, halo 2.5 (mass-weighted tests rely on this)
    mass_payload = b"".join(struct.pack("<f", 1.0) for _ in range(N0)) + \
                   b"".join(struct.pack("<f", 2.5) for _ in range(N1))
    # Gas-only blocks with known, easily assertable ranges
    u_payload = b"".join(struct.pack("<f", 10.0 + i * 0.01) for i in range(N0))
    rho_payload = b"".join(struct.pack("<f", 0.5 + (i % 100) * 0.01) for i in range(N0))
    hsml_payload = b"".join(struct.pack("<f", 2.0) for _ in range(N0))

    data = b""
    data += block("HEAD", make_header())
    data += block("POS ", pos_payload)
    data += block("VEL ", vel_payload)
    data += block("ID  ", id_payload)
    data += block("MASS", mass_payload)
    data += block("U   ", u_payload)
    data += block("RHO ", rho_payload)
    data += block("HSML", hsml_payload)

    with open(out_path, "wb") as f:
        f.write(data)
    print(f"wrote {out_path}: {ntot} particles ({N0} gas + {N1} halo), box {BOX}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "snap_test")
