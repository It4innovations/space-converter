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
| `gen_ramses.py` | RAMSES output directory (Fortran unformatted amr/hydro/part files + info and field descriptors): 2 CPU files, 15 AMR leaf cells over two levels, 40 DM + 20 stars (10 born before t = 0), 5 sinks in `sink_00001.csv` (3 born at t >= 0.015). |
| `gen_nchilada.py` | CHANGA_NCHILADA directory (XDR field files): 600 gas + 300 dark (+ empty star family), non-cubic box. |
| `gen_ipic3d.c` | IPIC3D_HDF5 restart + settings files (C, compile with `h5cc` from the HDF5 module): 2 species (400 e⁻ + 200 ions, opposite charges), 8³ grid with ghost layer, non-cubic box. |
| `gen_fil.c` | FIL / Einstein Toolkit Carpet HDF5 3D output (C, compile with `h5cc`): two refinement levels (level 0 split into two components with 2 inter-process ghost zones, level 1 a 9³ refined box), iterations 0 and 4, `HYDROBASE::rho` in per-process files (`*.file_N.h5`), `HYDROBASE::vel[0..2]` in a per-group file, `ADMBASE::alp`, plus a 2D slice file (`*.xy.h5`) with the same dataset names that must be ignored. Finest level wins: 702 + 729 = 1431 points. |
| `gen_fil_grace.c` | FIL_GRACE / GRACE HDF5 volume output (C, compile with `h5cc`): 11 blocks of 4³ cells on two levels, z >= 0 only (reflection symmetry), `rho`, `alp` and the vector `Bvec` (\|Bvec\| = 0.5); `volume_out_000004.h5` with `/Level`, `volume_out_000008.h5` without it (levels derived from the cell sizes). |
| `gen_genericio.cpp` | HACC_GENERICIO writer (C++, built by CMake as `gen_genericio` in WITH_HACC builds — the self-describing format is written via the bundled GenericIO library): 900 particles, x/y/z/vx/vy/vz/mass/rho/hh/id. |
| `protocol_client.py` | Minimal TCP client speaking the bspace wire protocol (info request + one sparse extraction). |
| `compare_raw.py` | Tolerant comparison of raw float32 volume dumps. |
| `run_smoke_tests.sh` | The test driver (T1–T11). Exit code = number of failures. |

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
- **T11** – `--dense-norm 2` (sum(qW)/sum(W)) on 1 and 2 ranks must finish and agree within tolerance (regression for the weight-grid MPI reduction that only rank 0 entered → every multi-rank run with `--dense-norm 1/2` hung; each rank runs under `timeout 300`).
- **T10** – per-reader extraction against the synthetic datasets above (TIPSY gas+dark, HACC_BIN, PLUTO_VTK, RAMSES (leaf cells, DM, `--ramses-filter`, sinks), NCHILADA gas+dark, IPIC3D both species — note iPIC3D counts include the reader's synthetic grid points without the tiles' ghost layer: 400+512 and 200+512 — FIL (Carpet HDF5: finest-level-wins point count 1431, `--fil-levels 0 0` → 729, `--fil-ghosts` → 1728, explicit `--fil-file` list including a 2D file, the vector group `vel` with magnitude 0.5, the default = latest iteration and `--fil-iteration 0`, and a 2-rank run with the same global count), FIL_GRACE (704 cells, iteration/time from the file attributes, `--fil-grace-mirror z` → 1408, `--fil-grace-levels 0 0` → 192, `--fil-grace-region` → 64, \|Bvec\| = 0.5, levels from cell sizes when `/Level` is missing, and a 2-rank run), HACC_GENERICIO via the `gen_genericio` tool, and GADGET/CodeBase, which reads the same format-2 snapshot as GADGET_SIMPLE with no parameter file needed; its dark-mass case additionally asserts min=max=2.5, a direct regression check for the "non-gas mass is 0" fix. Readers missing from the build are skipped with a hint. Every reader except BHAC has a synthetic case (BHAC was validated against real snapshots only).

To test all readers, build a "full" CPU variant first: enable `WITH_HACC`,
`WITH_IPIC3D`, `WITH_PLUTO` (turns on `WITH_VTK`), `WITH_RAMSES`, `WITH_FIL`,
`WITH_FIL_GRACE` (the HDF5 readers turn on `WITH_HDF5`) and `WITH_NANOFLANN` (so
T5 activates). The `build_spaceconverter_*.sh` scripts named below are the
developers' cluster scripts kept outside this repository (the per-cluster
scripts in `scripts/space_converter/` are the public starting point); paths are
examples. Barbora, `build_spaceconverter_barng_cpu_full.sh` (all 27 tests at
that time pass, none skipped, 2026-09-25):

```bash
srun --jobid=<JOBID> --overlap -N1 -n1 -c64 \
    bash <workspace>/scripts/build_spaceconverter_barng_cpu_full.sh
SC_BIN=<workspace>/install/space_converter_barng_cpu_full/bin/space_converter \
    srun --jobid=<JOBID> --overlap -N1 -n1 -c16 bash tests/run_smoke_tests.sh
```

The GPU paths (CUDA and HIP) are tested by `tests/run_gpu_tests.sh` (G1–G12, every
GPU result compared against the same binary run on the CPU), on a GPU node.

## Running (cluster, never on a login node)

```bash
# inside an existing Slurm job:
srun --jobid=<JOBID> --overlap -N1 -n1 -c16 bash tests/run_smoke_tests.sh
```

Configuration via environment: `SC_BIN` (binary path), `SC_OUT` (scratch dir,
default `<repo>/temp/tests`), `SC_NP` (max ranks for T4, default 3), `SC_NODE`
(pin the MPI ranks to one node of the job, e.g. when the other nodes are busy),
`SC_LUMI_PARTITION` (LUMI only: `G` for the HIP build (default), `C` for the
LUMI-C CPU builds; `cray-hdf5` is loaded for the `h5cc`-built generators).

Barbora CPU nodes with the FIL and FIL_GRACE readers (build:
`build_spaceconverter_bar_cpu_fil_grace.sh`; 44 pass, 0 fail, 2026-10-06):

```bash
ssh <node of the job> "cd <repo>; export SLURM_JOB_ID=<JOBID>; \
    SC_BIN=<workspace>/install/space_converter_bar_cpu_fil_grace/bin/space_converter \
    SC_NODE=<node> bash tests/run_smoke_tests.sh"
```

LUMI-C with the FIL reader (build: `build_spaceconverter_lumic_fil.sh`;
32 pass, 0 fail, 2026-09-29):

```bash
SC_BIN=<workspace>/install/space_converter_lumic_fil/bin/space_converter \
SC_LUMI_PARTITION=C SC_NODE=<node> \
    srun --jobid=<JOBID> --overlap -w <node> -N1 -n1 -c16 bash tests/run_smoke_tests.sh
```

MPI cases use `mpirun` when present in the environment, otherwise `srun --overlap` inside
the surrounding job.

Suggested order when tests fail after a change: T1/T2 first (they gate
everything), then T3/T6 (id conventions), then T4/T5 (MPI behavior).
