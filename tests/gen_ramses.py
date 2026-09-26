#!/usr/bin/env python3
"""Write a tiny synthetic RAMSES output directory (output_00001) for the smoke tests.

Two CPU files, periodic box (nx = 1, nboundary = 0), boxlen 10, 3 variables
(density, velocity_x..z, pressure):
  cpu 1: the level-1 grid (8 cells, cell 0 refined -> 7 leaves)
  cpu 2: the level-2 grid inside cell 0 (8 leaves)
  -> 15 leaf cells.
Particles: cpu 1 30 DM + 20 stars (birth_time -9.5 .. 9.5, 10 of them < 0),
cpu 2 10 DM -> 40 DM, 20 stars.

usage: gen_ramses.py OUT_PARENT_DIR   (creates OUT_PARENT_DIR/output_00001)
"""
import os
import struct
import sys

NCPU = 2
NLEV = 2
BOXLEN = 10.0
HYDRO = ["density", "velocity_x", "velocity_y", "velocity_z", "pressure"]
PART = [("position_x", "d"), ("position_y", "d"), ("position_z", "d"),
        ("velocity_x", "d"), ("velocity_y", "d"), ("velocity_z", "d"),
        ("mass", "d"), ("identity", "i"), ("levelp", "i"), ("family", "b"),
        ("tag", "b"), ("birth_time", "d"), ("metallicity", "d")]


def rec(f, fmt, values):
    data = struct.pack("<" + fmt, *values)
    f.write(struct.pack("<i", len(data)) + data + struct.pack("<i", len(data)))


def rec_str(f, s, n=128):
    data = s.ljust(n).encode()
    f.write(struct.pack("<i", len(data)) + data + struct.pack("<i", len(data)))


# grids[level] = (owner cpu, centre, son flags of the 8 cells)
GRIDS = {1: (1, (0.5, 0.5, 0.5), [1, 0, 0, 0, 0, 0, 0, 0]),
         2: (2, (0.25, 0.25, 0.25), [0] * 8)}


def numbl(icpu):
    """numbl(icpu, ilevel) of the file of CPU icpu, Fortran order (icpu fastest).
    Each file holds its own grids only (no ghost grids of the other CPU)."""
    v = []
    for lev in range(1, NLEV + 1):
        for cpu in range(1, NCPU + 1):
            v.append(1 if GRIDS[lev][0] == cpu == icpu else 0)
    return v


def write_amr(path, icpu):
    with open(path, "wb") as f:
        rec(f, "i", [NCPU]); rec(f, "i", [3]); rec(f, "3i", [1, 1, 1])
        rec(f, "i", [NLEV]); rec(f, "i", [100]); rec(f, "i", [0]); rec(f, "i", [2])
        rec(f, "d", [BOXLEN])
        rec(f, "3i", [1, 1, 1]); rec(f, "d", [1.0]); rec(f, "d", [1.0]); rec(f, "d", [0.5])
        rec(f, "2d", [0.1, 0.1]); rec(f, "2d", [0.1, 0.1]); rec(f, "2i", [10, 10])
        rec(f, "3d", [0, 0, 0]); rec(f, "7d", [0] * 7); rec(f, "5d", [1, 0, 1, 0, 0]); rec(f, "d", [0])
        nb = numbl(icpu)
        rec(f, f"{len(nb)}i", nb); rec(f, f"{len(nb)}i", nb); rec(f, f"{len(nb)}i", nb)
        rec(f, "20i", [0] * 20)
        rec(f, "5i", [0] * 5)
        rec_str(f, "hilbert")
        rec(f, "3d", [0.0, 1.0, 2.0])
        rec(f, "i", [1]); rec(f, "i", [0]); rec(f, "i", [1])
        for lev in range(1, NLEV + 1):
            owner, centre, sons = GRIDS[lev]
            for ib in range(1, NCPU + 1):
                if ib != owner or ib != icpu:
                    continue
                rec(f, "i", [lev]); rec(f, "i", [0]); rec(f, "i", [0])      # ind_grid, next, prev
                for d in range(3):
                    rec(f, "d", [centre[d]])
                rec(f, "i", [0])                                            # father
                for _ in range(6):
                    rec(f, "i", [0])                                        # nbor
                for c in range(8):
                    rec(f, "i", [sons[c]])                                  # son
                for _ in range(8):
                    rec(f, "i", [owner])                                    # cpu_map
                for _ in range(8):
                    rec(f, "i", [0])                                        # flag1


def write_hydro(path, icpu):
    with open(path, "wb") as f:
        rec(f, "i", [NCPU]); rec(f, "i", [len(HYDRO)]); rec(f, "i", [3])
        rec(f, "i", [NLEV]); rec(f, "i", [0]); rec(f, "d", [5.0 / 3.0])
        for lev in range(1, NLEV + 1):
            owner = GRIDS[lev][0]
            for ib in range(1, NCPU + 1):
                ncache = 1 if ib == owner == icpu else 0
                rec(f, "i", [lev]); rec(f, "i", [ncache])
                if ncache:
                    for c in range(8):
                        for v in range(len(HYDRO)):
                            rec(f, "d", [lev + 0.1 * c + 0.01 * v])


def particles(icpu):
    p = []
    if icpu == 1:
        for i in range(30):
            p.append(dict(x=(0.3 * i) % BOXLEN, fam=1, bt=0.0))
        for i in range(20):
            p.append(dict(x=(0.45 * i + 0.1) % BOXLEN, fam=2, bt=i - 9.5))
    else:
        for i in range(10):
            p.append(dict(x=(0.9 * i + 0.2) % BOXLEN, fam=1, bt=0.0))
    return p


def write_part(path, icpu):
    ps = particles(icpu)
    n = len(ps)
    with open(path, "wb") as f:
        rec(f, "i", [NCPU]); rec(f, "i", [3]); rec(f, "i", [n]); rec(f, "4i", [1, 2, 3, 4])
        rec(f, "i", [20]); rec(f, "d", [0.0]); rec(f, "d", [0.0]); rec(f, "i", [0])
        cols = {
            "position_x": [q["x"] for q in ps], "position_y": [(q["x"] * 1.7) % BOXLEN for q in ps],
            "position_z": [(q["x"] * 2.3) % BOXLEN for q in ps],
            "velocity_x": [1.0] * n, "velocity_y": [0.0] * n, "velocity_z": [0.0] * n,
            "mass": [0.5] * n, "identity": list(range(n)), "levelp": [2] * n,
            "family": [q["fam"] for q in ps], "tag": [0] * n,
            "birth_time": [q["bt"] for q in ps], "metallicity": [0.02] * n,
        }
        for name, kind in PART:
            rec(f, f"{n}{kind}", cols[name])


def main():
    out = os.path.join(sys.argv[1], "output_00001")
    os.makedirs(out, exist_ok=True)
    for icpu in range(1, NCPU + 1):
        write_amr(f"{out}/amr_00001.out{icpu:05d}", icpu)
        write_hydro(f"{out}/hydro_00001.out{icpu:05d}", icpu)
        write_part(f"{out}/part_00001.out{icpu:05d}", icpu)
    with open(f"{out}/info_00001.txt", "w") as f:
        f.write(f"ncpu        ={NCPU:11d}\nndim        ={3:11d}\nlevelmin    ={1:11d}\n"
                f"levelmax    ={NLEV:11d}\nngridmax    ={100:11d}\nnstep_coarse={10:11d}\n\n"
                f"boxlen      ={BOXLEN:23.15E}\ntime        ={0.5:23.15E}\naexp        ={1.0:23.15E}\n"
                f"unit_l      ={3.0857e21:23.15E}\nunit_d      ={6.77e-23:23.15E}\nunit_t      ={4.7e14:23.15E}\n")
    with open(f"{out}/hydro_file_descriptor.txt", "w") as f:
        f.write("# version:  1\n# ivar, variable_name, variable_type\n")
        for i, name in enumerate(HYDRO):
            f.write(f"{i + 1:2d}, {name}, d\n")
    with open(f"{out}/part_file_descriptor.txt", "w") as f:
        f.write("# version:  1\n# ivar, variable_name, variable_type\n")
        for i, (name, kind) in enumerate(PART):
            f.write(f"{i + 1:2d}, {name}, {kind}\n")
    print(out)


if __name__ == "__main__":
    main()
