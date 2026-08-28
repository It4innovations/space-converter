# SpaceConverter — Full Code Analysis & Bug Hunt (August 2026)

Scope: complete review of the project's own sources (`src/`, `blender_addons/bspace/`,
`scripts/`), with special focus on the **voxelization path** — how particle datasets from
the different solvers (GADGET, GADGET_SIMPLE, CHANGA TIPSY/NChilada, HACC bin/GenericIO,
iPIC3D HDF5, PLUTO VTK) are turned into dense/sparse VDB grids, and whether the result is
consistent between the command line and the bspace Blender addon, between CPU and GPU, and
between 1 and N MPI ranks. Third-party submodules were only read where the project calls
into them (cudaKDTree `FlexHeapCandidateList`, OpenVDB/NanoVDB APIs).

Severity legend: **HIGH** = produces wrong output or crashes in realistic configurations,
**MED** = wrong/inconsistent output in specific configurations, **LOW** = quality/usability.

---

## 1. Pipeline overview (as implemented)

```
main.cpp
 ├─ parse_args()                      args_processing.cpp
 ├─ init_converter()                  data_processing.cpp
 │   ├─ <solver>::init_lib()          reads snapshot (per-rank subset)
 │   ├─ find_particle_positions()     caches pos/mass/rho/hsml per particle type
 │   ├─ calculate_radius_by_{cudakdtree,nanoflann}()   optional k-NN radius+rho
 │   └─ sorts (radius / Morton), optional voxel KD-tree
 ├─ find_bbox()                       global bbox via MPI_Allreduce, cube-symmetrized
 └─ per request (TCP from bspace, or once in --export mode):
     ├─ recv_requested_data()         bbox_min/max, grid_dim, dense type/norm, filters…
     ├─ create_grid()                 dense dim³ buffer | sparse manager (OpenVDB/NanoVDB/GPU)
     ├─ convert_to_grid()             convert_iolib_to_grid → CPU/GPU kernels
     │   ├─ sparse: 1 particle → 1 voxel (no kernel splat), value accumulated
     │   └─ dense:  SPH splat via fill_voxels() (Cubic…WendlandC8), or voxel-centric k-NN loop
     ├─ reduction()                   MPI_Reduce (dense) / log-tree serialize+merge (sparse)
     ├─ finalize_grid()               normalization (Count/SPH/VoxelVolume), dense→VDB
     └─ send_vdb()/save_vdb()
```

Key coordinate conventions:

- `bbox_min_orig`/`bbox_size_orig` (ints/double): cube-symmetrized world-space bbox of the
  whole dataset (or of one type for the BBOX query).
- `bbox_min/bbox_max` (floats in `[0, object_size]`): the zoom sub-box requested by the
  client; `object_size` defaults to 1000.
- `scale_space_diagonal = |bbox_max-bbox_min| / (object_size·√3)`: zoom factor (1 = full box).
- voxel coordinate = `(pos − bbox_min_orig) / bbox_size_orig · bbox_dim / scale_space_diagonal`;
  the dense grid is a `bbox_dim³` window at `offset[] = bbox_min_norm · bbox_dim / scale`.
- `transform_scale = scale · object_size / bbox_dim` = voxel size in "object" units.

---

## 2. HIGH-severity findings (core voxelization)

### 2.1 Per-type particle caches are indexed with *global* reader IDs → any particle type > 0 is broken
- Where:
  - `src/common/convert_vdb_kernel.h:352-355, 412` (`process_particle`:
    `pos_particles[cached_idx*3]`, `value_particles[cached_idx]`),
    `convert_vdb_kernel.h:101-103` (`radius_particles[pid]`)
  - IDs come from `particles_id_ordered_per_ptype`, filled with the **global** reader index
    `i` in `find_particle_positions` (`src/common/convert_vdb.cpp:255`).
- The cached arrays (`pos_particles_per_ptype[ptype]`, `values_particles`,
  `radius_particles_per_ptype[ptype]`) are **compact per-type arrays** (size = particles of
  that type), but the kernels index them with global IDs. All readers use a global ID space
  grouped by type (see e.g. `gadget_simple_extract_iolib.cpp:1681-1707`:
  `get_particle_type(id)` + `get_particle_type_offset`). Hence:
  - ptype 0: global id == within-type id → works (this is why the bug goes unnoticed).
  - ptype ≥ 1: `cached_idx ∈ [offset, offset+N_t)` indexes arrays of size `N_t` →
    **out-of-bounds reads**; positions/values/radii are garbage or crash.
- The rest of the code *does* rebase correctly (`get_particle_rho`:
  `convert_vdb.cpp:1626` uses `id − particles_ptype_offset[type]`), which confirms the
  intended convention and that the kernels miss the rebase.
- Same pattern on GPU (`convert_vdb_kernel.cu:258, 381`: `idx = particle_ids[idx_c]`).
- Note: the `--skip-cache-manager` streaming path builds identity ids per chunk
  (`convert_vdb.cpp:813-814`) and is **correct** — cached vs. streaming runs of a
  multi-type dataset therefore disagree.
- Fix: store within-type indices in `particles_id_ordered_per_ptype` (0..N_t−1); the sorts
  in `data_cache.cpp` then also become self-consistent (see 3.6), and reader-facing calls
  must add `particles_ptype_offset[ptype]`.

### 2.2 `find_particle_values()` fetches values of the wrong particles for ptype > 0
- Where: `src/common/convert_vdb.cpp:316-318`:
  ```cpp
  for (size_t idx = 0; idx < num_particles; idx++) {
      float value = get_particle_norm_value(block_name_id, idx);   // idx = within-type
  ```
  `get_particle_norm_value(blocknr, id)` expects a **global** id
  (`convert_vdb.cpp:1659-1665` → `get_particle_rho(id)` → `get_particle_type(id)`).
- For ptype ≥ 1 the values cached are those of the *first N_t particles in the file*
  (i.e. particles of type 0), not of the selected type. Combined with 2.1 the extracted
  block values for any non-first particle type are doubly wrong.
- Fix: `get_particle_norm_value(block_name_id, particle_ids[idx])` (with global ids), or
  keep `idx` but rebase inside once 2.1 is fixed to within-type convention.

### 2.3 k-NN radius/density with MPI cycling is destroyed every round (cudaKDTree, CPU & GPU)
- Where: `src/utility/cudakdtree_tool.cu:338` (GPU: `cudaMemset(d_cand, 0, …)` inside the
  round loop) and `:635` (CPU: `memset(cand.data(), 0, …)`), combined with the continuation
  constructor `FlexHeapCandidateList(list, k, round == 0 ? maxRadius : -1.f)`
  (`cudakdtree_tool.cu:83-85, 451-453`).
- cukd semantics (`submodules/cudaKDTree/cukd/knn.h:361-372, 391-397`): with
  `cutOffRadius < 0` the list is *used as-is*; `push()` inserts only if
  `entry < entry[0]`. A zeroed list has `entry[0] = 0` (dist2 = 0), so **no candidate can
  ever be inserted in rounds ≥ 1**, and `returnValue()` = 0.
- Consequence: with `--calc-radius-neigh` and `use_cycling` (the default whenever
  `anim_type == eNone`) on **more than 1 MPI rank**, the last round overwrites the results
  with `radius = 0` and `rho = NaN` (h_inv = ∞ → `q = NaN` → `W = NaN`) for *every*
  particle. Radius 0 silently degrades dense SPH splatting to single-voxel deposits.
  Single-rank runs (only round 0) are unaffected — another CLI-vs-cluster inconsistency.
- Additional design flaw: with `CUDAKDTREE_NUM_SPLITS > 1` the candidate buffer only holds
  one batch, so cross-round continuation is impossible even without the memset; the
  extract step also runs (and overwrites outputs) every round instead of once at the end.
- The nanoflann implementation (`nanoflann_tool.cpp:196-360`) does cycling correctly
  (accumulates all rounds' candidates, sorts at the end) — but see 3.7/3.8.
- Fix: allocate the candidate list for *all* queries, initialize it once before the round
  loop, never clear between rounds, extract once after the last round (i.e. restore the
  `#if 0 // nosplit` structure, or per-batch persistent slices `cand[start*k .. )`).

### 2.4 `particle_radius_multiplier` is received over TCP but never MPI-broadcast
- Where: `src/data_communication.cpp:458` (recv) vs. broadcast list `:472-486` — the
  multiplier is missing from the `MPI_Bcast` block; every other received field is there.
- On a multi-rank server driven from bspace, "Particle Size Multiplier" applies only to
  rank 0's particle subset; all other ranks keep `0.0f` (= multiplier disabled, since
  `fill_voxels` treats 0 as "don't scale"). Grids differ per rank and the merged output is
  a mixture. The user's production script `run_spaceconverter_barng_cpu.sh` (n=6144) is
  exactly this configuration.
- Fix: add `MPI_Bcast(&space_data.particle_radius_multiplier, sizeof(float), MPI_BYTE, 0, …)`.

### 2.5 Cube-symmetrization of the global bbox mixes new min with old max
- Where: `src/data_processing.cpp:498-504`:
  ```cpp
  space_data.bbox_min_orig[0] = (min+max)/2 - S/2;          // overwrites min
  space_data.bbox_max_orig[0] = (bbox_min_orig[0]+max)/2 + S/2;  // uses NEW min + OLD max
  ```
  With extent `e` on an axis and `S` = max extent, the resulting extent is `(e+3S)/4`
  instead of `S` — the "cube" is not a cube for any non-cubic dataset.
- `bbox_min_orig`/`bbox_size_orig` (used by the conversion kernels) stay self-consistent,
  but `bbox_max_orig` is wrong. It is consumed by `send_bbox`
  (`data_communication.cpp:400-414`), which normalizes **per axis** by
  `dims_orig[i] = bbox_max_orig[i] − bbox_min_orig[i]`. The per-type sub-bbox shown in the
  addon (and any sub-box the user then selects from it) is therefore scaled by the wrong
  factor per axis → zoom regions land displaced relative to what the conversion kernels
  cut out. Only exact cubes are safe.
- Also: the assignments truncate `double → int` toward zero instead of `floor` (negative
  coordinates round the wrong way; masked by the ±1 padding at `:485-491`).
- Fix: compute the center once from the *original* min/max, then
  `min = c − S/2; max = c + S/2` per axis (with `floor`/`ceil`).

---

## 3. MED-severity findings

### 3.1 CPU and GPU report different `min/max` semantics
- CPU sparse/dense kernels track min/max of the **per-particle values** that passed the
  filters (`convert_vdb_kernel.cpp:228-229, 452-453`); the GPU paths compute min/max of
  the **accumulated voxel values** after reduce (`convert_vdb_kernel.cu:565, 758, 1213`
  `find_min_max` over `d_vals` / `d_data_density` — the dense variant includes all the
  zero voxels, so `min` is almost always 0).
- These values are sent to bspace (`send_vdb`, `data_communication.cpp:535-538`) and used
  for shader/value mapping → the same dataset produces differently normalized volumes with
  and without `--gpu`.

### 3.2 `--simple-density` silently disables the spherical filter (and value filters)
- `process_particle` (`convert_vdb_kernel.h:414-417`): `use_simple_density` returns `true`
  *before* the `filter_min/max` and `bbox_sphere` checks. Skipping the value filter is
  arguably intended (value ≡ 1), but skipping `--bbox-sphere` is a spatial-filter loss the
  user cannot see.

### 3.3 Dense extraction with `dense_type = eNone` writes (almost) nothing
- `fill_voxels` computes `W_norm(dense_type=eNone) = 0` (`dense_utility.h:593-620` has no
  eNone branch) so every splat with radius ≥ 1 voxel contributes 0; only sub-voxel-radius
  particles hit the `iradius == 0` fast path. The CLI allows this state only transiently
  (`--dense-type 0` forces sparse mode, `args_processing.cpp:224-227`), but the TCP path
  does not validate it — the addon *can* send `extracted_type=eDense, dense_type=0` and
  get an empty grid with no diagnostic.

### 3.4 SPH density estimate uses the query particle's mass for all neighbors
- `cudakdtree_tool.cu:117, 486` and `nanoflann_tool.cpp:349`:
  `rho[i] = Σ_j m_i · W(r_ij/h_i)` — SPHtoGrid's estimator is `Σ_j m_j · W(...)`.
  Correct only for equal-mass simulations; systematically wrong densities for variable
  particle masses (Gadget star formation, zoom-in simulations, multi-species runs).
  (Neighbor masses are admittedly not available for remote-rank trees — that limitation
  should at least be documented.)

### 3.5 nanoflann rho accumulates on top of stale reader densities
- `nanoflann_tool.cpp:326-349`: `rho_particles.resize(N)` then `rho_particles[q] += …`.
  The vector passed in is `cache_manager.rho_particles_per_ptype[ptype]`, already filled
  with the reader's own `get_particle_rho` values by `find_particle_positions`
  (`convert_vdb.cpp:249-250, 270`). `resize` to the same size keeps them → the k-NN rho is
  *added to* the solver's rho instead of replacing it. The cudakdtree path memsets its
  device buffer and is unaffected. `--nanoflann` vs `--cudakdtree-cpu` therefore give
  different densities on datasets that provide a rho block.

### 3.6 The two particle sorts assume mutually incompatible id conventions
- `sort_particles_by_radius_cpu` (`data_cache.cpp:98-113`) sorts positions-in-the-array by
  `radii[i]` (within-type indexing) — consistent *before* any other sort.
- `sort_particles_by_nonoverlap_cpu` (`data_cache.cpp:177-193`) reads
  `positions[ids[i]*3]` — i.e. treats the stored ids as *within-type* indices, but they
  are global reader ids (see 2.1): out-of-bounds/garbage Morton keys for ptype ≥ 1, and
  after a radius sort even for ptype 0 the radius/`ids` pairing used by a subsequent
  Morton sort is stale. Fixing 2.1 (within-type ids everywhere) resolves this too.

### 3.7 NanoVDB Float3 grids: MPI merge reinterprets Vec3f grids as float grids
- `sparse_common.cpp:522-546` (`VoxelNanoVDBManager<T>::merge(bin_data)`): the received
  buffer is cast to `nanovdb::NanoGrid<float>*` regardless of `T`; for the
  `Vec3fGrid` instantiation (used when `--no-norm-value`) the traversal reads garbage and
  `floatToValue<Vec3f>(f)` splats scalars into all three components. The OpenVDB manager
  does this correctly with a typed `gridPtrCast<T>` + `compSum` (`sparse_common.cpp:895`).
- Wider design note: the sparse insert interface is scalar-only
  (`insertOrUpdatePackedSequential(key, float)`), and vector blocks are reduced to their
  magnitude by the readers' `RETURN_NORM_VECTOR3` macro anyway — so the whole Float3
  ("no-norm-value") pipeline currently stores `(|v|,|v|,|v|)`, never a real vector.

### 3.8 Radius `∞` / rho `∞` sentinels are propagated as real values
- When fewer than k neighbors are found the code stores `inf` in radius and rho
  (`nanoflann_tool.cpp:356-358`, `cudakdtree_tool.cu:133-135`). `fill_voxels` then
  computes `iradius = (int)inf` (UB; in practice INT_MIN/huge) and loops over the whole
  grid row range, or `get_particle_radius` returns `inf · len2pix`. No downstream
  filtering of non-finite radii exists.

### 3.9 Dense-grid offset cast: negative zoom box wraps to huge `size_t`
- `convert_vdb.cpp:685-687`:
  `offset[i] = (size_t)(bbox_min_norm · bbox_dim / scale)` — if the client sends a bbox
  with a negative component (the addon's BBOX empty can extend below 0), the negative
  double converts to a huge unsigned value; the dense window then never intersects any
  particle (empty output) or indexes wildly. Sparse mode has the complementary problem:
  voxel coords are computed with `static_cast<int>` truncation (`convert_vdb_kernel.h:406-408`),
  so the two voxels adjacent to a zero-plane collapse into voxel 0 (`ifloor` exists in the
  same header but is not used here).

### 3.10 `packCoord3` silently wraps beyond ±2²⁰ voxels
- `sparse_common.h:86-96`: 21-bit biased packing with `& COORD_MASK`, no clamp. Voxel
  indices grow as `bbox_dim / scale_space_diagonal`; a deep zoom (e.g. `--grid-dim 1000`
  and a sub-box smaller than ~0.1 % of the domain) exceeds 2²⁰ and aliases to unrelated
  coordinates — corrupted key space with no warning.

### 3.11 Dense radius-tolerance test mixes voxel units with normalized units
- `process_particle` (`convert_vdb_kernel.h:373-391`): `radiusxyz_max` is in **voxel**
  units (from `get_particle_radius`) but is added to `out_px_norm` (normalized [0,1]) and
  compared against `bbox_*_norm`. The tolerance is `bbox_dim/scale` times too large →
  effectively no rejection for dense grids. Correctness is saved by clipping inside
  `fill_voxels`, so the impact is wasted work per particle, not wrong output.

### 3.12 GPU vs CPU: `kernel_norm` uses `powf` on doubles in device code
- `dense_utility.h:94-98` etc.: device branch computes `powf(h_inv, (float)dim)` in single
  precision while the host uses `std::pow` in double — one more small CPU/GPU numeric
  divergence inside every splat (h_inv³ for small h in voxel units can be large).

---

## 4. LOW-severity / robustness notes (core)

- `args_processing.cpp` has no unknown-argument error and no bounds checks on `argv[++i]`
  — typos are silently ignored (this converts several script bugs below into silent
  behavior changes); a flag at the end of the line dereferences `argv[argc]`.
- `usage()` omits `--radius-const`, `--no-norm-value`, `--cub`; the example commands at
  the bottom of the file use retired arities (`--anim S E`, 3-arg `--export-data`).
- `--dense-type` / `--raw-particles` / `--simple-density` all overwrite
  `extracted_type` — order of flags on the command line changes the result.
- `find_bbox` truncation `static_cast<int>(min − 1.0f)` rounds toward zero (masked by the
  1-unit padding).
- `save_vdb` condition `data_processing.cpp:1178-1179` relies on `&&`/`||` precedence
  without parentheses; currently equivalent to the (probable) intent but fragile.
- `VoxelOpenMPManager::insertOrUpdatePackedSequential` drops the voxel silently when the
  table is full (`sparse_common.cpp:86-105`); serialize ships the whole (mostly empty)
  table over MPI.
- `find_min_max` GPU dense scans all `dim³` voxels including implicit zeros — reported
  `min_value_reduced` is ~always 0 regardless of data.
- `mpi_reduce`d `size_t` particle counts use `MPI_LONG_LONG` (`data_processing.cpp:718`) —
  fine on LP64, wrong on Windows/LLP64 only in spirit (both 64-bit).
- OpenMP `#pragma omp parallel for` over `size_t` loop indices appears in several places;
  fine for GCC/Clang OpenMP ≥ 3.0, breaks MSVC builds.
- `save_raw_particles_to_vtp`: `points->SetPoint` inside `#pragma omp parallel for` is
  safe only because `SetNumberOfPoints` preallocated — but `vtkCellArray::InsertNextCell`
  right after is serial; OK. The parallel `SetPoint` on vtkPoints is not guaranteed
  thread-safe by VTK, though.
- `--server` is parsed but dead — the converter always listens
  (`data_communication.cpp:263` + `renderengine_tcp.h` default `is_server = true`).

---

## 5. CLI vs. bspace addon consistency

The addon does **not** build a command line; it speaks the TCP protocol
(`recv_requested_data`). Wire order and sizes match on both sides (verified field by
field). Differences that make "the same dataset" voxelize differently:

| Parameter | CLI (batch) | addon (TCP) | Consequence |
|---|---|---|---|
| `particle_radius_multiplier` | **cannot be set at all** | property, default 0.0 | addon results with multiplier ≠ 0 are unreproducible from the CLI; multi-rank servers apply it on rank 0 only (bug 2.4) |
| `grid_dim` | default 100, scripts use 100–4000 | default 100 | different resolutions per launch path |
| port | default 7000 | pref default 5000 | server unreachable with defaults on both sides |
| `filter_min/max` | ±FLT_MAX | ±1e39 → stored as **±inf** (float32 overflow in Blender property) | passes-all today, but diverges from CLI defaults and from any `isfinite` handling |
| `dense_type` | 0–6 (0 forces sparse) | enum "1"–"6" only | addon cannot express eNone (good); CLI order-dependence (see §4) |
| `--bbox-sphere`, `--offset-position`, `--no-norm-value`, `--radius-const`, `--simple-density` | CLI only | not exposed | addon cannot reproduce CLI runs using them |
| `object_size` | fixed 1000 default | `bbox_size` property | must match or coordinates shift |
| CUB output (`--cub`) | supported | **missing** from `file_type_items` (`bspace_panel.py:81-87`) | `IndexError` when a `--cub` server responds |

Other addon findings (from the addon/scripts review):

- `bspace_pref.py:29`: port `max=65565` (typo for 65535).
- `bspace_panel.py:1363, 1754`: index guard `> len(list)` instead of `>=` → `IndexError`
  at `index == len`.
- `recvall` returns short buffers on disconnect → downstream `struct.error` instead of a
  clear "not connected" (`bspace_panel.py:2353-2372`).
- NanoVDB reception with empty `nvdb_converter_path` silently imports a `.vdb` that was
  never created (`bspace_panel.py:1901-1945`).
- Volume sequences configured only for `anim_type == '1'`, but the server sends
  `frames = world_size` for FrameExtract (3) too.
- Temp path built by string concatenation from an (empty by default) pref → writes to
  filesystem root; no `os.path.join`.
- One shared `BSPACE()` socket as a class attribute; never closed on unregister; draw
  handler leaks; `register_export` checkbox state not restored on file load.
- `grid_dim`/`bbox_size` properties lack `min=` clamps; a 0/negative grid dim is sent to
  the server unvalidated (the server does not validate either).

---

## 6. Cluster / launch script findings

- **`--gpu-cuda` no longer exists** — `run_spaceconverter_dev_gpu.sh:65` and
  `run_spaceconverter_p12_gpu.sh:62` silently run the CPU pipeline (parser only knows
  `--gpu`; unknown args are ignored, see §4). HIGH for anyone comparing GPU numbers.
- **`--data-type GENERICIO`** in `aurora/run_aurora.sh:46` and `polaris/run_polaris.sh:48`
  → immediate `runtime_error` (valid name is `HACC_GENERICIO`).
- `barbora/run_bar.sh:53-58`: relative `--gadget-file` (the `data=` variable is never
  used) and stale 3-argument `--export-data 1 1 ${out_d}` — the output path is silently
  swallowed.
- `karolina/run_kar.sh:168`: `--grid-dim 1000 … --grid-dim 100` — last one wins (100) on a
  4096-rank run.
- Old arities (`--anim START END` without STEP, bare `--raw-particles`) in many commented
  example lines — re-enabling them throws in `std::stoi` or reads `argv[argc]`.
- `polaris/run_polaris.sh:48` calls `gpu_bind_kar2.sh` which does not exist in the tree;
  Polaris has 4 GPUs, kar-style binding assumes 8.
- `gpu_bind_bar.sh:33`: `CUDA_VISIBLE_DEVICES=$LOCAL_RANK` without `% 4` — breaks with
  more than 4 tasks/node on Barbora.
- `leonardo/run_leo.sh:39`: uses an unset `${cosmo25_nc_path}`, no output path/grid dim,
  default port 7000 (which no addon default matches).
- `karolina/gpu.sh` and `gpu_bind_kar.sh` end with unquoted `$@` and no `exec`; only
  `lumi/gpu.sh` does `exec "$@"` correctly.
- The user's `run_spaceconverter_barng_cpu.sh` (GADGET_SIMPLE, n=6144, port 5000) parses
  cleanly; it is however the configuration most affected by bugs 2.3 (if
  `--calc-radius-neigh` is ever added) and 2.4 (radius multiplier on 6144 ranks).

---

## 7. Solver readers (per-format data flow and issues)

### 7.1 What each reader hands to the voxelizer

| Solver | Position | "mass" | "rho" | hsml (radius) | Endianness | MPI split |
|---|---|---|---|---|---|---|
| ChaNGa Tipsy | `gas/darks/stars[].pos` | tipsy `mass` (all 3 types) | gas `rho` | gas `hsmooth`, DM/stars `eps`; aux `smoothlength` wins | XDR in TipsyReader; **aux swap keyed on host, not file** | contiguous slice per rank, remainder → last rank |
| ChaNGa NChilada | per-type XDR field files | `mass` field — **Dark mass unreadable (7.2-H3)** | gas `GasDensity` | gas `smoothlength`/`soft`, DM/stars `soft` | XDR | each field file sliced independently |
| GADGET (CodeBase) | `P[id].Pos` via `read_ic` | gas only, **0 for other types** | `SphP` gas only | `P[id].Hsml` gas only, else 0 | handled by gadget lib | gadget domain decomposition |
| GADGET_SIMPLE | POS block per type | MASS block or header table | RHO block | HSML block | header swapped, **payloads never swapped** | whole files or per-file chunks |
| HACC bin | packed 58-B records | `mass` field | `rho` field | `hh` | assumed LE | contiguous range, remainder → last rank |
| HACC GenericIO | user-named vars | `--mass-name` var, else 0 | `--rho-name` var, else 0 | `--hsml-name` var, else 0 | GenericIO lib | GIO ranks over MPI ranks |
| iPIC3D HDF5 | particles + synthetic grid nodes | **signed charge `q`** (grid: 1.0) | grid `rho` moment (particles: 1.0) | constant 1.0 | HDF5 | proper `partition_range` |
| PLUTO VTK | rectilinear node coords as "particles" | `rho · Δx³` | first cell array (assumed rho) | `Δx` (x-spacing only) | VTK lib | **computed but never applied (7.2-H1)** |

### 7.2 HIGH-severity reader bugs

- **H1 — PLUTO: MPI decomposition never applied.**
  `pluto_vtk_extract_iolib.cpp:124-126, 294-296`: `start_cell = world_rank * cells_per_rank`
  is computed and discarded; `get_particle_position(id)` decodes `id` as a *global* cell
  index. With >1 rank every rank voxelizes the same first chunk (counted multiple times in
  the reduction) and the middle of the grid is only covered by the last rank.
- **H2 — PLUTO: point-dims vs cell-dims confusion → sheared volume.**
  `pluto_vtk_extract_iolib.cpp:109-121, 268-281`: the `-1` on `cell_dims` was commented
  out. Cell arrays have `(nx−1)(ny−1)(nz−1)` tuples with x-stride `nx−1`, but positions
  are decoded with stride `nx`; the id→value vs id→position mappings diverge by one cell
  per row, progressively shearing the whole volume; ids past the array end silently read
  value 0 / rho 1.0.
- **H3 — NChilada: dark-matter mass always 0.**
  `changa_nchilada_extract_iolib.cpp:357-368` (and 423-434, 485-494): `BlockType::Mass`
  missing from the `Dark` switch in all three accessors although the file is loaded and
  the block advertised. Mass-weighted voxelization of DM contributes nothing.
- **H4 — GADGET_SIMPLE (`nfiles < world_size` path): missing 4-byte Fortran record tag.**
  `gadget_simple_extract_iolib.cpp:715-716` vs the correct `:934-938`: the shared-file
  path never consumes the block's leading length tag, so every value is shifted by one
  32-bit word (x/y/z rotate through wrong slots). The multi-file path is correct — the
  same snapshot read with a different rank/file ratio gives different data.
- **H5 — GADGET_SIMPLE: MASS block offset ignores the header mass table.**
  `gadget_simple_extract_iolib.cpp:1631-1639`: skips `npart` of *all* lower types although
  the MASS block only contains types with `header.mass[type] == 0` (the common case: DM in
  the table, gas per-particle) → masses read from the wrong place.
- **H6 — GADGET_SIMPLE: gas-only blocks (RHO/HSML/U…) read for every type.**
  `fill_known_blocks` marks all blocks active for all 6 types
  (`gadget_simple_extract_iolib.cpp:525-529`); reading them for types 1–5 runs past the
  block into tags/next blocks, and the garbage is then *advertised as available data*
  (`:1913`, `get_count() > 0`).

### 7.3 MEDIUM-severity reader bugs

- **M1 — HACC GenericIO: rho block id underflow.** `hacc_genericio_extract_iolib.cpp:1301-1303`:
  without `--rho-name`, `get_particle_rho_blocknr()` returns `BTMax + (−1) = 1` = the
  **Vel** block; the base class then hijacks velocity requests and returns 0 for all.
- **M2 — HACC GenericIO: `gio_*_idx` indexes `variable_info` but `gio_datas` is compacted**
  (`:422-455`): one unrecognized variable type shifts every later index (latent).
- **M3 — iPIC3D: types/blocks table written block-major (`BTMax*s + block`) but read
  type-major (`PTMax*blocknr + type`)** — `ipic3d_hdf5_extract_iolib.cpp:378-392` vs
  `ipic3d_hdf5_convert_vdb.cpp:111,127,157`. The (type, block) availability shown in
  bspace is scrambled.
- **M4 — iPIC3D: mass = signed charge** (`:316-323`) — electron species get negative
  weights (voxel subtraction, /≈0 normalization). Also `get_global_num_particles()`
  returns the local count without reduction (`:957-958`).
- **M5 — iPIC3D: real particles and synthetic grid nodes share one type** — value-0
  phantoms dilute count-normalized output; `/fields/E*,B*` duplicated into every species.
- **M6 — GADGET_SIMPLE & HACC GenericIO: `world_size` not a multiple of `nfiles` →
  trailing ranks open nonexistent files** (`gadget_simple…:546-548`, `hacc_genericio…:527-529`);
  remainder never distributed (iPIC3D's `partition_range` is the correct pattern).
- **M7 — GADGET_SIMPLE: no endianness conversion of payloads** — header detection exists,
  `gadget_simple_swap_Nbyte` is never applied to the data (`:716, 938`); a big-endian
  snapshot reads as garbage after a valid header.
- **M8 — Tipsy aux files: byte swap keyed on host endianness** (`TipsyFile.cpp:715, 789-805`)
  — assumes aux is always big-endian; a native little-endian tipsy + aux dataset gets its
  `smoothlength` byte-swapped into garbage radii. Double-precision aux silently passes the
  component check as "2 components" and reads as zeros.
- *(M9 — the bbox symmetrization bug — is the same finding as §2.5, confirmed independently.)*
- **M10 — GADGET (CodeBase): mass/rho/hsml zero for all non-gas types**
  (`gadget_extract_iolib.cpp:1957-2016`) — inconsistent with tipsy/nchilada/gadget-simple:
  the same physics ingested via a different reader voxelizes differently. `SphP[id]`
  dereferenced with no type guard for RHO/U requests on non-gas ids (OOB).
- **M11 — PLUTO: block-id space off by one** — advertised blocks `0..n−1`, accessors map
  `blocknr − BTMax` with `BTMax = 1`; every non-rho scalar shifted by one, last scalar
  unreachable (`pluto_vtk_extract_iolib.cpp:330-338` vs `:190-192, 362-364`).

### 7.4 LOW-severity reader notes

- HACC bin: `int` rank arithmetic overflows at ≥2³¹ particles (`haccbin…:328-330`);
  same pattern in tipsy (masked by 32-bit header counts).
- GADGET_SIMPLE: stale/possibly negative `block_idx` in the mass-table branch
  (`:746, 968`, `knownblocks[-1]` UB); file handle opened twice, closed once
  (`:569, 627`); every `fread` return value ignored (truncated files → silent garbage).
- PLUTO: `cell_size` from `x_coords[1]−x_coords[0]` only — wrong hsml and `rho·Δx³` mass
  on stretched/anisotropic grids; `--scalar-names` parsed but never used.
- iPIC3D: grid points placed at nodes though the comment says cell centers (half-cell
  offset); hsml hard-coded 1.0 regardless of `Dx`.
- NChilada: each field file sliced by its own `numParticles` — count mismatch between
  block files silently misaligns them.
- All wrappers: anim mode uses `anim_start + anim_step·world_rank` and ignores
  `anim_end` — ranks beyond the range read timesteps past the end.

### 7.5 Cross-solver semantics ("what is density?")

The voxelizer treats `mass`, `rho`, `hsml` and the block value as uniform quantities, but
their meaning differs per reader (see 7.1): masses are code units (gadget's `1e10 M⊙/h`
conversion is deliberately commented out, `gadget_extract_iolib.cpp:248`) vs physical
(`rho·Δx³` in PLUTO) vs charge (iPIC3D) vs zero (HACC GenericIO without `--mass-name`,
gadget non-gas). `rho` is an SPH density, a grid moment, "whatever array 0 is", or 0.
`hsml` is a smoothing length, a softening, a constant 1.0, or an x-spacing. Consequences:

- `DenseNorm::eSPHInterpolation` and any rho-dependent path behave differently or
  degenerately per solver;
- particles with hsml = 0 silently collapse to single-voxel deposits (gadget non-gas,
  GenericIO without `--hsml-name`) — visually a different "grain" than SPH-smoothed gas;
- `eVoxelVolume` produces a "physical" density only when the deposited value is actually
  a mass in consistent units — true for some readers only.

A per-reader unit/semantic contract (documented and asserted at init) is the systemic fix.

---

## 8. Software-engineering review (beyond the numerics)

An assessment of the codebase as software, independent of the domain bugs above. Several
of the HIGH bugs in §2/§7 are direct consequences of these structural issues.

### 8.1 Architecture & API design

- **The 40-parameter function.** `convert_iolib_to_grid` and both kernels take ~40 positional
  parameters (`convert_vdb.h`, `convert_vdb_kernel.h:529-656`), most copied 4× (base → CPU
  dispatcher → CPU kernels; base → GPU dispatcher → GPU kernels). Adding one parameter means
  editing 8+ signatures; forgetting one call site compiles fine because nearly everything is
  `float`/`int`/`double` (this is exactly how bug 2.4 — the un-broadcast multiplier — and the
  unused `min_rho/max_rho/grid_transform` parameters happened). A single `ConversionRequest`
  struct passed by const reference (and `MPI_Bcast`-able as one block) removes the class of bug.
- **No single source of truth for the id convention.** "Particle id" means *global reader id*
  in the readers and `get_particle_rho`, but *within-type index* in `find_particle_values`
  and implicitly in the kernels (bugs 2.1/2.2, 3.6). Nothing in the type system distinguishes
  them — both are `size_t`/`uint64_t`. Strong typedefs (`GlobalId`, `TypeLocalId`) or accessor
  methods on the cache would have made the bug uncompilable.
- **Duplicated protocol definitions.** The TCP wire format exists as hand-ordered
  `send_data_data`/`recv_data_data` call sequences in C++ *and* as a mirrored sequence of
  `struct.pack/unpack` calls in `bspace_panel.py`, with no shared schema, no version field,
  no magic/handshake, and native endianness. Every new field must be added in 4 places
  (send, recv, Bcast, Python) in exactly the same position; bug 2.4 is one missed place.
  The same applies to the enum duplication (`file_type_items` in Python vs
  `FileTypeIdentifier` in C++ — the missing CUB entry is the Python copy lagging).
- **Dispatch by `dynamic_cast` chains.** `finalize_grid`, `save_vdb`, `sparse_to_*`,
  `merge_grid` each probe the concrete manager type with 2–4 `dynamic_cast`s and duplicate
  the body per type. The `VoxelSparseManager` interface exists but is bypassed; making
  `serialize/merge/extrema/save` proper virtuals would collapse ~400 lines and remove the
  type-punning NanoVDB Float3 merge bug (3.7), which is precisely a case where one branch
  of a copy-pasted chain was edited and another was not.
- **CPU/GPU kernel duplication.** `process_particle`/`fill_voxels` are shared (good), but the
  surrounding orchestration (batching, min/max, counters) is written twice with different
  semantics (bug 3.1). The nosplit/batched `#if 0` forks in `cudakdtree_tool.cu` left the
  live path broken (bug 2.3) — dead alternatives kept in-tree with `#if 0` instead of git
  history are actively dangerous here.

### 8.2 Error handling & robustness

- **Errors are `printf` + `return`.** Readers and managers report failure to stdout and
  return normally (`read_radius_from_file`, gadget-simple `fread`s, VoxelOpenMP table-full
  drop). The pipeline then continues with partial/zero data — the worst failure mode for a
  scientific tool, because the output *looks* valid. There is no error propagation path at
  all between reader → converter → client.
- **`exit()` deep inside library code** (readers call `exit(-1)` on open failure; `print_info`
  calls `exit(0)`) — kills all MPI ranks without `MPI_Abort`, risking hangs on partial exits.
- **No input validation layer.** Neither `parse_args` (unknown flags ignored, `argv[++i]`
  unchecked) nor `recv_requested_data` (grid dim ≤ 0, negative bbox, enum ranges) validates
  anything; the addon does not clamp either. Every malformed input becomes a silent numeric
  bug instead of a diagnostic (§4, 3.3, 3.9).
- **Sentinel values instead of status.** `radius = ∞`, `rho = ∞/NaN`, `EMPTY_KEY = -999999`,
  `calc_radius_neigh = -1` — all flow into arithmetic unfiltered (3.8).

### 8.3 Resource & concurrency hygiene

- Raw `new[]`/`delete[]` and `cudaMalloc/cudaFree` pairs throughout (`sparse_common.cpp`,
  `convert_vdb_kernel.cu`, `cudakdtree_tool.cu`) with early returns between them — leaks on
  every error path; no RAII wrappers despite C++17.
- `#pragma omp critical` around a per-neighbor `emplace_back` (`nanoflann_tool.cpp:308-313`)
  serializes the hottest loop; thread-local buffers are used correctly elsewhere
  (`convert_to_sparse_grid_cpu`), so the idiom is known but inconsistently applied.
- OpenMP loops over `size_t` indices and unparenthesized `&&/||` conditions (§4) —
  MSVC-incompatible and precedence-fragile respectively.
- The remote server loop holds all caches forever; `clear()`ing dense grids but sharing
  `sparse_grid` pointers between `grid_main`/`grid_main_sum`/`grid_main_final` makes
  ownership untraceable (works only because everything is `shared_ptr`).

### 8.4 Testing, build, observability

- **There are no tests.** Not one unit test in the project's own code (`test_converter` is
  an empty TODO). Every bug in §2 would be caught by: a 2-type synthetic dataset round-trip
  test (2.1/2.2), a 2-rank cycling test against the 1-rank result (2.3, 2.4), a non-cubic
  bbox assertion (2.5), and a reader golden-file test per format (§7). These are cheap to
  write: the kernels are pure functions of arrays.
- **Commented-out code as version control.** Hundreds of lines of `#if 0`/`//`-blocks,
  including entire alternate implementations, retired protocols, and stale example command
  lines that no longer parse (§4). They actively mislead (the correct PLUTO cell-center
  code exists *only* as a comment next to the wrong live code, 7.2-H2).
- **Debug printf as logging.** `LOG_Measure*` exists but most diagnostics are bare `printf`
  with `rank #%d` prefixes; no verbosity levels, no way to silence per-chunk spam
  (`[insertOrUpdate] Pre-agg: …` on every merge) on 6144 ranks.
- Magic numbers inline: `COORD_BITS 21`, `HASH_TABLE_LOAD_FACTOR 0.5`, chunk `1<<20`,
  `MAX_MPI_BYTES`, per-file duplicated `RETURN_NORM_*` macros (8 copies, already divergent
  in their NaN handling comments).
- The addon: one 2700-line `bspace_panel.py` mixing UI, protocol, shader setup and file
  I/O; module-level socket singleton; no reconnect/timeout strategy beyond one retry.

### 8.5 What would give the most engineering leverage

1. A `ConversionRequest`/`GridSpec` struct + one `MPI_Bcast(&req, sizeof req)` (kills the
   §2.4 class), and a versioned, length-prefixed TCP message with a shared field list.
2. Strong id types or a `CacheView` accessor that hides the global↔per-type rebase
   (kills the §2.1 class).
3. `Result`/exception-based error propagation from readers; replace `exit()` with
   `MPI_Abort`; check every `fread`.
4. A minimal test harness: synthetic 2-type particle set + golden VDB hashes, run at
   1 and 2 ranks, CPU and GPU — this single fixture would have caught 2.1–2.5, 3.1, 3.5.
5. Delete `#if 0` archaeology (git remembers), promote the sparse-manager interface to
   real virtuals, unify the duplicated GPU/CPU orchestration.

---

## 9. Consistency matrix — "same dataset, different result" causes

| Axis | Cause | Ref |
|---|---|---|
| 1 vs N MPI ranks | k-NN cycling destroys candidates (radius=0, rho=NaN) | 2.3 |
| 1 vs N MPI ranks | radius multiplier applied on rank 0 only | 2.4 |
| CPU vs GPU | min/max semantics (particle values vs voxel values) | 3.1 |
| CPU vs GPU | `powf` vs `pow` in kernel norm | 3.12 |
| cached vs `--skip-cache-manager` | global-vs-compact id mismatch (ptype > 0) | 2.1 |
| `--nanoflann` vs `--cudakdtree-cpu` | rho accumulates on stale reader rho | 3.5 |
| CLI vs addon | radius multiplier not settable from CLI; defaults tables §5 | §5 |
| ptype 0 vs ptype ≥ 1 | id-space bug + wrong value block | 2.1/2.2 |
| cubic vs non-cubic dataset | bbox symmetrization + per-axis send_bbox scaling | 2.5 |
| zoom depth | packCoord3 21-bit wrap; negative-offset size_t wrap | 3.9/3.10 |
| flag order on CLI | `--dense-type`/`--raw-particles`/`--simple-density` overwrite each other | §4 |
| 1 vs N MPI ranks (PLUTO) | decomposition computed but not applied — chunks duplicated/missing | 7.2-H1 |
| rank/file ratio (GADGET_SIMPLE) | shared-file path misaligned by the missing record tag | 7.2-H4 |
| reader choice (GADGET vs GADGET_SIMPLE vs Tipsy) | non-gas mass/hsml = 0 vs real values | 7.3-M10 |
| DM vs gas (NChilada) | dark mass always 0 | 7.2-H3 |
| species sign (iPIC3D) | mass = signed charge; electrons subtract | 7.3-M4 |
| file endianness (GADGET_SIMPLE, Tipsy aux) | payloads not/incorrectly swapped | 7.3-M7/M8 |

## 10. Priority fix list

Correctness (voxelization output changes):

1. **2.1 + 2.2 + 3.6** — switch `particles_id_ordered_per_ptype` to within-type indices,
   rebase reader calls with `particles_ptype_offset`; fixes every ptype ≥ 1 extraction and
   both sorts. (One coordinated change; touch CPU + GPU kernels + `find_particle_values`.)
2. **2.3** — restore persistent candidate lists across cycling rounds in
   `cudakdtree_tool.cu` (GPU and CPU paths); extract once after the last round.
3. **2.4** — add the missing `MPI_Bcast` of `particle_radius_multiplier`.
4. **7.2-H1/H2** — PLUTO: apply `start_cell`, restore `cell_dims = dims − 1` and
   cell-center positions.
5. **7.2-H4/H5/H6 + 7.3-M7** — GADGET_SIMPLE block reading (record tag, mass-table offset,
   per-type block presence, payload endianness).
6. **7.2-H3** — NChilada: add `BlockType::Mass` to the Dark switches.
7. **2.5** — correct the bbox symmetrization (compute center before overwriting) and use
   `floor`; make `send_bbox` use `bbox_size_orig` for all axes.
8. **3.1** — unify min/max semantics (recommend: voxel-value extrema everywhere, they are
   what the shader mapping actually needs).
9. **3.5** — zero the rho buffer before accumulation in nanoflann; **7.3-M1** GenericIO
   rho-block underflow; **7.3-M3** iPIC3D table convention; **7.3-M4** `fabs(q)`.

Tooling / integration:

10. Scripts: replace `--gpu-cuda` → `--gpu`, `GENERICIO` → `HACC_GENERICIO`; add an
    unknown-argument error in `parse_args` so this class of bug becomes loud.
11. Addon: add CUB file type, fix port max, `>=` guards, port-default alignment.
12. Guard rails: reject non-finite radii (3.8), clamp/validate client bbox and grid dim,
    clamp `packCoord3` inputs with a warning (3.10); distribute the file remainder in
    GADGET_SIMPLE/GenericIO rank mapping (7.3-M6).

Engineering (prevents recurrence — see §8.5):

13. `ConversionRequest` struct + single-`Bcast` + versioned TCP schema.
14. Strong id types / cache accessor hiding the global↔per-type rebase.
15. Error propagation instead of `printf`+continue; `MPI_Abort` instead of `exit`.
16. Minimal test fixture: synthetic 2-type dataset, golden hashes, 1 vs 2 ranks,
    CPU vs GPU — would have caught 2.1–2.5, 3.1, 3.5 and most of §7.
17. Remove `#if 0` archaeology; make sparse managers real virtual interfaces.

---

## 11. Cleanup & hardening plan (actionable checklists)

Concrete follow-up work items complementing §8/§10 — suitable for splitting into small,
reviewable commits. Each item is independent unless noted.

### 11.1 Code cleanup checklist

- [ ] **Delete all `#if 0` / commented-out implementation blocks** (git history keeps them):
      `cudakdtree_tool.cu` (nosplit forks), `sparse_common.cpp` (extractAll variants,
      temp-file serialize), `convert_vdb.cpp` (commented thread-collection blocks in both
      `calculate_radius_by_*`, old OpenVDB grid paths), `data_processing.cpp`
      (nano_grid/vdb_grid remnants, commented timing code), `convert_vdb_kernel.h`
      (commented GPU prototypes). Exception: first *revive* the correct PLUTO cell-center
      code (7.2-H2) before deleting its commented version.
- [ ] **Remove dead parameters** threaded through the 40-arg signatures: `min_rho`,
      `max_rho`, `grid_transform`, `grid_name` (unused in CPU kernels), `bbox_min/bbox_max`
      past the norm computation; remove unused `dims_orig_max` (`data_communication.cpp:405`).
- [ ] **Remove or implement dead flags**: `--server` (parsed, never used — converter always
      listens); `--dense-file` naming (it takes an enum int, not a file).
- [ ] **Deduplicate macros**: one shared header for `RETURN_NORM_*` / `RETURN_ORIG_*`
      (8 divergent copies across readers), `EMPTY_KEY`/`HASH_TABLE_LOAD_FACTOR` (defined in
      both `sparse_common.h` and `.cpp`), `FBSKIP`-style record-tag helpers in the gadget
      readers.
- [ ] **Delete stale example command lines** at the bottom of `args_processing.cpp`
      (~120 lines, several with retired argument arities) — move 2–3 *working* examples
      into `README.md`/`usage()` instead.
- [ ] **Drop per-call debug spam** on hot paths (`[insertOrUpdate] Pre-agg…` per merge,
      `Starting round…` per rank) behind a verbosity flag; route the rest through
      `utility/logging.h` instead of bare `printf`.
- [ ] **Naming/style pass**: `sended:` → `sent:`, "does exist" → "does not exist"
      (`convert_vdb.cpp:87`), mixed tab/space indentation in `sparse_common.cpp` and
      `nanoflann_tool.cpp`; add `.clang-format` and format `src/` (one-time commit).
- [ ] Parenthesize the compound conditions flagged in §4 (`save_vdb`,
      `save_raw_volume`, `save_vti_volume`) even where currently equivalent.

### 11.2 Comment & in-code documentation plan

- [ ] **Document the particle-id convention at its source**: a block comment on
      `CacheManager` (`data_cache.h`) defining *global reader id* vs *within-type index*,
      which arrays use which, and that `particles_ptype_offset` is the bridge — plus
      one-line comments at every rebase site. (Do this together with the fix for 2.1.)
- [ ] **Write the §1 coordinate math into the code**: the derivation of
      `scale_space_diagonal`, `transform_scale`, `len2pix`, `offset[]` currently lives in
      nobody's head — put it as a header comment in `convert_vdb_kernel.h` above
      `get_particle_radius`/`fill_voxels`, with units (world / normalized / voxel) for
      each variable.
- [ ] **Protocol table**: one comment block (mirrored verbatim in `data_communication.cpp`
      and `bspace_panel.py`) listing every field of every message in order, with C type,
      Python struct code and size — until a real schema (§8.5-1) replaces it.
- [ ] **Per-reader unit contract** (§7.5): a standard header comment in each
      `*_extract_iolib.cpp` stating what `get_particle_mass/rho/hsml` return, in which
      units, for which types, and what is 0/unsupported — so the voxelizer's assumptions
      are checkable per solver.
- [ ] **Fix comments that are wrong**: iPIC3D "cell center" comment vs node positions
      (§7.4); `find_bbox` "Expand by 1 unit" (it also truncates); doc-comments on
      `usage()` examples with retired arities; `VoxelDenseManager` header notes referring
      to removed members.
- [ ] Doxygen hygiene: many `@param` lists in `data_processing.cpp` /
      `convert_vdb_kernel.cpp` are copy-pasted and list parameters that the function does
      not have (or miss new ones) — regenerate or trim to short prose.

### 11.3 Argument parsing & help overhaul

- [ ] **Fail loudly on unknown arguments**: final `else { fprintf(stderr, "Unknown
      argument: %s\n", arg); usage(); }` in `parse_args` — this single change turns the
      silent script bugs of §6 (`--gpu-cuda`, stale `--export-data` arity) into visible
      errors. Readers that re-parse `argc/argv` for their own flags must be given a shared
      whitelist (or parse first and strip), otherwise solver-specific flags would now trip
      the error — simplest: collect all known flags in one table used by both.
- [ ] **Bounds-check every `argv[++i]`** (helper `next_arg(i, argc, argv, "--flag")` that
      errors instead of dereferencing `argv[argc]`), and wrap `std::stoi/stof` to report
      *which* flag had a bad value.
- [ ] **Validate values**: `--grid-dim > 0` (and warn above e.g. 4096); enum ranges for
      `--dense-type [0..6]`, `--dense-norm [0..3]`, `--raw-particles [0..3]`,
      `--dense-file [0..2]`; `--port [1..65535]`; `--anim start<=end, step>0`;
      `--bbox min<max` per axis. Apply the same checks to values arriving via TCP in
      `recv_requested_data` (the addon must not be able to request `grid_dim = 0`).
- [ ] **Remove order-dependence**: `--dense-type`, `--simple-density`, `--raw-particles`
      each overwrite `extracted_type`; parse flags first, derive `extracted_type` once at
      the end with a defined precedence, and error on contradictions
      (e.g. `--raw-particles` + `--dense-type 6`).
- [ ] **Complete `usage()`**: add the missing `--radius-const`, `--no-norm-value`,
      `--cub`, `--port` default, `--export` arity note; state defaults for every option
      (grid-dim 100, port 7000, object_size 1000); one working example per data type
      (replacing the dead ones, see 11.1); mention which options are TCP-only
      (`particle_radius_multiplier`, filters from the addon) so batch users know what they
      cannot reproduce (§5).
- [ ] Align defaults with the addon while touching this: either change CLI default port
      to 5000 or the addon pref to 7000 (§5), and document it in `usage()`.

### 11.4 Smoke-test plan (minimal but catching the real bugs)

Infrastructure: a `tests/` directory with CTest; a tiny generator
`tests/gen_synthetic.py` writing (a) a GADGET_SIMPLE snapshot and (b) a Tipsy file with
**two particle types** (e.g. 1000 gas + 500 dark), known masses (unequal!), positions on a
jittered lattice inside a deliberately **non-cubic** box, and a known rho/hsml block.
Assertions compare against golden values with tolerance, or between two runs of the
converter itself (self-consistency needs no golden data at all).

- [ ] **T1 — arg parsing**: unknown flag exits non-zero; truncated `--bbox 1 2 3` exits
      non-zero; `--grid-dim 0` exits non-zero. (Guards 11.3; trivial shell tests.)
- [ ] **T2 — single-rank golden**: `--export 0 <rho> --grid-dim 32` sparse and dense
      (WendlandC6, each `--dense-norm`); assert particle count, min/max, and a hash (or
      sum + a few probed voxels) of the output grid. Catches regressions in the whole
      §2/§3 area once fixed.
- [ ] **T3 — ptype > 0**: same as T2 with `--export 1 <mass>`; assert the *dark* count and
      that the deposited total equals the dark mass sum. Fails today (2.1/2.2); becomes
      the regression test for the id-convention fix.
- [ ] **T4 — rank invariance**: run T2 with `mpirun -n 1` and `-n 2` (and `-n 3` to hit
      remainder paths, 7.3-M6); the merged grids must be bit-identical (dense reduce is
      order-independent for these counts; sparse compare after sort). Catches 2.3, 2.4,
      7.2-H1, H4 today.
- [ ] **T5 — k-NN sanity**: regular lattice with spacing `d`, `--calc-radius-neigh 6
      --cudakdtree-cpu`: assert radius ≈ d for interior particles and rho > 0 finite;
      repeat under `mpirun -n 2` (fails today, 2.3); repeat with `--nanoflann` and assert
      both backends agree (fails today, 3.5).
- [ ] **T6 — streaming equivalence**: T2 with and without `--skip-cache-manager` must
      match exactly (the code claims it; the test enforces it — and fails for ptype>0
      today).
- [ ] **T7 — CPU vs GPU** (only where a GPU is present; tag as optional in CTest): T2
      with `--gpu`, tolerance-compare grid values, and assert the *reported* min/max match
      CPU semantics after the 3.1 unification.
- [ ] **T8 — protocol round-trip**: a ~100-line python client (reusing the addon's pack
      order, or the addon module imported headless) that connects to a `--port` server,
      requests info + bbox + one extraction, and checks it receives a parseable VDB and
      sane min/max. Guards the wire format (§5) without Blender.
- [ ] **T9 — non-cubic bbox**: dataset with extents 1:2:4; assert the symmetrized bbox is
      a cube containing all particles and that a sub-box extraction of the dense grid
      equals the corresponding window of the full-box grid (fails today, 2.5).
- [ ] **T10 — reader golden files**: commit tiny (<100 kB) fixture snapshots per format
      (both endiannesses for GADGET_SIMPLE/Tipsy-aux once M7/M8 are fixed; a 2-type
      NChilada set for H3; a small PLUTO VTK with CELL_DATA for H1/H2) and assert
      positions/masses/values of a few known particles via `--info` + a debug dump flag.

Suggested order: T1/T2 first (they gate everything), then T4+T5 (the MPI bugs), then the
rest as the corresponding fixes land. All tests run in seconds; wire them into CI
(`ctest` after build) so the §6 class of silent drift cannot recur.
