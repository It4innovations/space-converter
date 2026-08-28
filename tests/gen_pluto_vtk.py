#!/usr/bin/env python3
"""Generate a small synthetic PLUTO_VTK file (legacy VTK, RECTILINEAR_GRID).

Matches what src/pluto/pluto_vtk_extract_iolib.cpp reads: a legacy .vtk
rectilinear grid with CELL_DATA scalar arrays (PLUTO's standard VTK output).
ASCII data mode keeps the generator trivial (vtkRectilinearGridReader handles
both ASCII and binary).

Grid: 9 x 13 x 17 points -> 8 x 12 x 16 = 1536 cells, deliberately non-cubic
and with NON-UNIFORM spacing on the x axis (exercises the per-cell spacing
fix). Two cell arrays: rho (first array = density) and prs.
"""

import sys

NX, NY, NZ = 9, 13, 17          # point dimensions
LEN = (8.0, 24.0, 64.0)         # domain extents, non-cubic


def coords(n, length, stretch):
    # Mildly stretched spacing: x_i = length * (i/n)^stretch
    return [length * (i / (n - 1)) ** stretch for i in range(n)]


def main(out_path):
    xs = coords(NX, LEN[0], 1.3)   # non-uniform on purpose
    ys = coords(NY, LEN[1], 1.0)
    zs = coords(NZ, LEN[2], 1.0)
    ncells = (NX - 1) * (NY - 1) * (NZ - 1)

    lines = []
    lines.append("# vtk DataFile Version 2.0")
    lines.append("space-converter synthetic PLUTO grid")
    lines.append("ASCII")
    lines.append("DATASET RECTILINEAR_GRID")
    lines.append(f"DIMENSIONS {NX} {NY} {NZ}")
    for name, vals in (("X_COORDINATES", xs), ("Y_COORDINATES", ys), ("Z_COORDINATES", zs)):
        lines.append(f"{name} {len(vals)} float")
        lines.append(" ".join(f"{v:.6g}" for v in vals))
    lines.append(f"CELL_DATA {ncells}")

    def cell_values(fn):
        vals = []
        for k in range(NZ - 1):
            for j in range(NY - 1):
                for i in range(NX - 1):
                    vals.append(fn(i, j, k))
        return vals

    # rho: a smooth, strictly positive deterministic field
    rho = cell_values(lambda i, j, k: 1.0 + 0.1 * i + 0.01 * j + 0.001 * k)
    lines.append("SCALARS rho float")
    lines.append("LOOKUP_TABLE default")
    lines.append(" ".join(f"{v:.6g}" for v in rho))

    prs = cell_values(lambda i, j, k: 0.5 + 0.05 * k)
    lines.append("SCALARS prs float")
    lines.append("LOOKUP_TABLE default")
    lines.append(" ".join(f"{v:.6g}" for v in prs))

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out_path}: {ncells} cells ({NX}x{NY}x{NZ} points), extents {LEN}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "pluto_test.vtk")
