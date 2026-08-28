# space-converter smoke tests

Small end-to-end tests against a deterministic synthetic dataset. They map to
the test plan in `docs/SpaceConverter_Code_Analysis_2026-08.md` §11.4 and act as
regression tests for the bug fixes described there (§2, §7).

## What is here

| File | Purpose |
|---|---|
| `gen_gadget_simple.py` | Format-2 GADGET_SIMPLE snapshot: 1000 gas + 500 halo particles (unequal masses) on jittered lattices in a **non-cubic** box; POS/VEL/ID/MASS for all types, U/RHO/HSML gas-only. Deterministic. |
| `gen_tipsy.py` | CHANGA_TIPSY snapshot (standard/XDR big-endian): 700 gas + 350 dark, non-cubic box, unequal masses. |
| `gen_haccbin.py` | HACC_BIN file: 800 particles, 58-byte packed records. |
| `gen_pluto_vtk.py` | PLUTO legacy VTK rectilinear grid: 8x12x16 cells, non-uniform x spacing, `rho` + `prs` CELL_DATA arrays. |
| `gen_nchilada.py` | CHANGA_NCHILADA directory (XDR field files): 600 gas + 300 dark (+ empty star family), non-cubic box. |
| `gen_ipic3d.c` | IPIC3D_HDF5 restart + settings files (C, compile with `h5cc` from the HDF5 module): 2 species (400 e⁻ + 200 ions, opposite charges), 8³ grid with ghost layer, non-cubic box. |
| `gen_genericio.cpp` | HACC_GENERICIO writer (C++, built by CMake as `gen_genericio` in WITH_HACC builds — the self-describing format is written via the bundled GenericIO library): 900 particles, x/y/z/vx/vy/vz/mass/rho/hh/id. |
| `protocol_client.py` | Minimal TCP client speaking the bspace wire protocol (info request + one sparse extraction). |
| `compare_raw.py` | Tolerant comparison of raw float32 volume dumps. |
| `run_smoke_tests.sh` | The test driver (T1–T10). Exit code = number of failures. |

## Tests

- **T1** – argument validation: unknown flag, truncated `--bbox`, `--grid-dim 0` must be rejected.
- **T2** – single-rank sparse extraction of particle type 0: particle count and a `.vdb` output.
- **T3** – same for particle **type 1** (regression for the global-vs-compact particle-id bug — garbage/OOB before the fix).
- **T4** – dense grid (WendlandC6, `--radius-const`, RAW dump) must be equal **within tolerance** (`compare_raw.py`) for 1..N MPI ranks — bit-identity is not expected because float accumulation order differs with rank count and OpenMP scheduling (regression for the MPI reduction/decomposition class of bugs; also caught the `LOG_MeasureStart/Stop` barrier deadlock).
- **T5** – k-NN radius sanity at 1 and 2 ranks (regression for the destroyed-candidate-list cycling bug: radius = 0 / rho = NaN). Skipped when the build has no cudakdtree/nanoflann.
- **T6** – `--skip-cache-manager` streaming run must produce a grid equal within tolerance to the cached run (type 1 — also exercises the id fix).
- **T7** – (not scripted) CPU vs GPU comparison; requires a GPU build.
- **T8** – TCP protocol round-trip with `protocol_client.py` (guards the wire format shared with the addon).
- **T9** – the non-cubic dataset's bbox must be symmetrized to the longest axis (regression for the mixed min/max symmetrization bug).
- **T10** – per-reader extraction against the synthetic datasets above (TIPSY gas+dark, HACC_BIN, PLUTO_VTK, NCHILADA gas+dark, IPIC3D both species — note iPIC3D counts include the reader's synthetic grid points: 400+1000 and 200+1000 — HACC_GENERICIO via the `gen_genericio` tool, and GADGET/CodeBase, which reads the same format-2 snapshot as GADGET_SIMPLE with no parameter file needed; its dark-mass case additionally asserts min=max=2.5, a direct regression check for the "non-gas mass is 0" fix. Readers missing from the build are skipped with a hint. **Every reader is covered.**

To test ALL readers, build the "full" CPU variant first (adds HACC, PLUTO/VTK
via the ParaView module, and nanoflann so T5 activates):

```bash
srun --jobid=<JOBID> --overlap -N1 -n1 -c64 \
    bash /mnt/proj1/open-36-34/milanjaros/projects/blender/scripts/build_spaceconverter_barng_cpu_full.sh
SC_BIN=/mnt/proj1/open-36-34/milanjaros/projects/blender/install/space_converter_barng_cpu_full/bin/space_converter \
    srun --jobid=<JOBID> --overlap -N1 -n1 -c16 bash tests/run_smoke_tests.sh
```

For the GPU-side test plan (to be executed on a GPU cluster), see
`docs/GPU_Test_Plan.md`.

## Running (cluster, never on a login node)

```bash
# inside an existing Slurm job:
srun --jobid=<JOBID> --overlap -N1 -n1 -c16 bash tests/run_smoke_tests.sh
```

Configuration via environment: `SC_BIN` (binary path), `SC_OUT` (scratch dir,
default `<repo>/temp/tests`), `SC_NP` (max ranks for T4, default 3). MPI cases
use `mpirun` when present in the environment, otherwise `srun --overlap` inside
the surrounding job.

Suggested order when tests fail after a change: T1/T2 first (they gate
everything), then T3/T6 (id conventions), then T4/T5 (MPI behavior).
