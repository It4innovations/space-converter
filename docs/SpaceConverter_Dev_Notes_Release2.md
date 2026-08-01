# SpaceConverter Development Branch Release2 Notes

## Comparison Basis

| Main reference (Release1) | `main` at `b19c7a7` |
| Development reference (Release2) | `origin/dev` at `1081ec2` |

## Release2 Overview

The development branch is a major architectural and performance update. Particle conversion has been reorganized around reusable dense and sparse voxel managers, a shared particle cache, and CPU/CUDA conversion kernels. The release adds an end-to-end GPU path, a CUB-backed sparse representation, new iPIC3D HDF5 and PLUTO VTK readers, expanded particle and dense export formats, and a substantially updated BSpace add-on.

The build is also reorganized into an installable CMake package. HACC-related readers now live under a unified module, generic CSV and generic HDF5 readers have been removed, and the former multiresolution and conditional temporary-buffer options are no longer part of the active build configuration.

## Highlights

- Added CUDA execution for dense and sparse particle-to-voxel conversion.
- Added optional CUDA-aware MPI for direct GPU-buffer reduction and transfer.
- Added optional CUDA managed-memory allocation.
- Added CPU, OpenMP, NanoVDB, OpenVDB, and CUB voxel-manager implementations.
- Added a shared particle cache with CPU/GPU storage and sorting operations.
- Added particle-radius and non-overlapping spatial sorting modes.
- Added a voxel-centric dense conversion mode backed by KD-tree queries.
- Added iPIC3D HDF5 input with particle species, electromagnetic fields, charge density, and pressure moments.
- Added PLUTO legacy VTK rectilinear-grid input.
- Added VTK PolyData (`.vtp`) particle output and VTK ImageData (`.vti`) dense output.
- Added native Blender point-cloud import alongside Geometry Nodes particle visualization.
- Updated BSpace to version 1.1.0 and Blender 5.0.
- Added an exported CMake target set and package configuration.
- Updated OpenVDB and added the BRAAS HPC render-engine submodule.

## 1. Conversion Architecture

### 1.1 Unified namespace

Core code and format adapters now use the top-level `space_converter` namespace. Common components are nested under `space_converter::common`, with VDB and cache functionality below it. This removes previously mixed global, `common`, and format-specific ownership.

Downstream C++ code using internal headers must update qualified names, for example:

```cpp
space_converter::common::SpaceData
space_converter::common::vdb::ConvertVDBBase
space_converter::common::vdb::VDBParticles
```

### 1.2 Dense and sparse manager abstractions

The monolithic grid container has been replaced by manager interfaces:

- `VoxelDenseManager` owns regular-grid dimensions, offsets, primary values, and the optional normalization buffer.
- `VoxelSparseManager` defines sparse initialization, insertion, transformation, serialization, deserialization, merging, memory reporting, and cleanup.
- `VDBParticles` now owns managers through `std::shared_ptr` rather than embedding every possible representation directly.

The representation enum is now:

```text
eDense
eVector
eNanoVDB
eOpenVDB
eCUB
eRawParticles
```

This design allows the processing pipeline to choose CPU or GPU storage without duplicating high-level conversion and reduction logic.

### 1.3 New source modules

The following core modules were added:

- `src/common/data_cache.{h,cpp,cu}`
- `src/common/convert_vdb_kernel.{h,cpp,cu}`
- `src/common/sparse_common.{h,cpp,cu}`
- `src/common/dense_common.cu`
- `src/common/raw_common.h`

Responsibilities are now separated as follows:

| Module | Responsibility |
|---|---|
| `data_cache` | Particle positions, values, radii, densities, masses, ordering, and GPU mirrors |
| `convert_vdb_kernel` | Shared CPU and CUDA particle-to-voxel conversion kernels |
| `sparse_common` | Sparse storage back ends and serialized sparse formats |
| `dense_common` | CPU dense manager and common dense behavior |
| `dense_common.cu` | GPU dense allocation, copies, and CUDA lifecycle |
| `raw_common` | Raw particle field containers and binary serialization |

### 1.4 Removed TCP implementation duplication

`src/common/connection_tcp.{h,cpp}` was removed. TCP transport is supplied through the newly added BRAAS HPC render-engine submodule, consolidating network functionality with the external rendering infrastructure.

## 2. Particle Cache and Preprocessing

### 2.1 `CacheManager`

A new `CacheManager` centralizes reusable particle data:

- positions grouped by particle type;
- original and reordered particle identifiers;
- radii grouped by type;
- per-type particle offsets;
- density and mass arrays;
- selected field values;
- CPU and GPU working buffers;
- CUDA KD-tree data where enabled;
- execution flags such as GPU use and sorting mode.

The manager provides explicit initialization and cleanup and replaces several independent arrays previously owned by `ConvertVDBBase`.

### 2.2 Position and value caching

Particle positions and selected field values are collected once and reused by:

- bbox calculation;
- radius calculation;
- sorting;
- dense conversion;
- sparse conversion;
- CPU/GPU transfer;
- raw-particle export.

Field lookup was optimized to reduce repeated virtual reader calls during voxelization.

### 2.3 Radius sorting

The new `--sort-by-radius` option sorts particle identifiers by effective particle radius. CPU and GPU implementations preserve the mapping back to original particle data.

`--skip-cache-manager` turns the particle cache off on the CPU path. Instead of materialising
positions, radii, values and identifiers for every particle (about 32 bytes each), the
conversion pulls particles from the reader in chunks of 1 M and splats each chunk with the
same kernel, so results are unchanged. It costs a second pass over the snapshot and is what
makes snapshots larger than memory tractable — `cosmo25.2304g.nc` (547 GB, 12.2 G particles
per type) needs 4 LUMI-G nodes with the cache but fits on 2 without it. The option is ignored,
with a note on rank 0, when something needs random access to all particles at once: `--gpu`,
`--sort-by-radius`, `--sort-by-non-overlap`, `--calc-radius-neigh`, or raw particle export.

This ordering is intended to improve memory locality and conversion behavior when support radii differ significantly.

### 2.4 Non-overlap sorting

The new `--sort-by-non-overlap` option spatially reorganizes particles into groups whose voxel contributions do not overlap. This permits conversion passes with reduced or eliminated atomic-update contention.

Both CPU and CUDA code paths are available. The algorithm operates on cached positions and radii and produces an ordered particle-ID sequence consumed by conversion kernels.

## 3. GPU and CUDA Support

### 3.1 Build options

The release adds:

| CMake option | Default | Purpose |
|---|---:|---|
| `WITH_GPU_CUDA` | OFF | Compile CUDA dense/sparse conversion kernels |
| `WITH_CUDA_AWARE_MPI` | OFF | Use CUDA-aware MPI for direct GPU-buffer communication |
| `WITH_CUDA_MALLOCMANAGED` | OFF | Use CUDA unified/managed memory |
| `WITH_GPU_HIP` | OFF | Compile the same kernels for AMD GPUs via HIP/ROCm |
| `WITH_HIP_AWARE_MPI` | OFF | Use GPU-aware MPI for direct device-buffer communication |
| `WITH_HIP_MALLOCMANAGED` | OFF | Use HIP unified/managed memory |

`WITH_CUDAKDTREE` remains available for neighbor search. It is separate from `WITH_GPU_CUDA`: KD-tree acceleration can be built without enabling the complete GPU voxelization pipeline.

`WITH_GPU_HIP` and `WITH_GPU_CUDA` are mutually exclusive. Enabling the HIP backend also
defines `WITH_GPU_CUDA` internally, because the GPU code paths are guarded by that macro;
`src/utility/gpu_device_compat.h` remaps the CUDA runtime API onto HIP and aliases
`cub` to `hipcub`, so the `.cu` sources compile unchanged for both vendors. The AMD build
compiles them through thin `src/common/hip/*.cpp` wrappers that only `#include` the
matching `.cu` file, so no source is duplicated.

### 3.2 CUDA language level and compiler options

- CUDA targets use CUDA C++17.
- Extended lambdas are enabled for NanoVDB CUDA functionality.
- Windows builds enable the conforming MSVC preprocessor through NVCC.
- Selected NVCC diagnostics are suppressed for host/device call and deep kernel stack warnings.
- CUDA debug builds enable host and device debug information.
- CUDA toolkit include paths are propagated to C++ compilation where shared headers expose CUDA types.

### 3.3 CUDA initialization order

When GPU conversion is compiled, the process initializes a CUDA context before `MPI_Init`. This avoids UCX/MPI implementations probing CUDA without an established context.

The initialization:

1. queries device count;
2. selects device zero as the initial device;
3. forces context creation;
4. reports a warning when initialization fails.

Rank-aware device selection and platform launch scripts remain responsible for final device binding.

### 3.4 GPU dense manager

`VoxelGPUDenseManager` adds device-side density and normalization buffers. It supports:

- allocation and zero initialization;
- optional managed memory;
- host-to-device and device-to-host synchronization;
- direct use by CUDA conversion kernels;
- cleanup through the manager lifecycle;
- MPI reduction either through CUDA-aware MPI or host staging.

### 3.5 GPU sparse manager

`VoxelGPUManagerSortReduce` implements sparse insertion as a GPU sort-and-reduce pipeline:

1. conversion kernels emit packed voxel keys and values;
2. CUB radix sort orders key/value pairs;
3. CUB reduction combines duplicate voxel keys;
4. the compact result is retained on device or serialized;
5. final data can be converted to OpenVDB or NanoVDB.

Capacity, temporary sort storage, reduced counts, and serialization alignment are managed explicitly.

### 3.6 CUB output representation

A new `eCUB` internal grid type and `FTI_CUB = 5` transport identifier were added.

The CUB representation is designed for GPU-friendly sparse data. Its serialized form contains metadata followed by aligned packed keys and float values. The implementation includes fixes for:

- serialized padding and 64-bit key alignment;
- transform preservation;
- min/max evaluation;
- sparse merges;
- CUB temporary-storage use;
- Linux and Windows compilation.

The `--cub` CLI switch selects this representation.

### 3.7 CUDA-aware MPI

When both `WITH_GPU_CUDA` and `WITH_CUDA_AWARE_MPI` are enabled, dense and supported sparse reductions can operate directly on device buffers. Otherwise, the same manager stages data through host memory.

This option requires an MPI stack built for CUDA/UCX or equivalent GPU-direct transport. It is disabled by default because unsupported MPI installations can fail at runtime.

### 3.8 GPU diagnostics

New files `gpu_logging.{h,cpp}` and `gpu_utility.h` provide:

- rank-aware GPU logging;
- optional CUDA synchronization around measured sections;
- CUDA error checking;
- debug printing of device arrays;
- managed-memory allocation macros;
- consistent GPU resource diagnostics.

Verbose GPU logging remains compile-time controlled.

## 4. Conversion Algorithms

### 4.1 Shared CPU/CUDA kernels

Sparse and dense conversion logic was moved into matching `.cpp` and `.cu` implementations with a common declaration layer. Both paths receive the same:

- cached position and value arrays;
- bbox and transform metadata;
- particle radii and masses;
- normalization mode;
- selected SPH kernel;
- filter and selection settings.

This reduces behavioral drift between CPU and CUDA implementations.

### 4.2 Sparse conversion managers

The release provides:

- an OpenMP sparse manager;
- a templated NanoVDB manager for scalar and `Vec3f` grids;
- a templated OpenVDB manager for scalar and vector grids;
- a GPU CUB sort/reduce manager.

Each manager implements its own insertion, merge, serialization, and transform behavior behind `VoxelSparseManager`.

### 4.3 Dense conversion modes

The default dense conversion continues to loop over particles and splat each particle into its support region.

A new optional mode, `--dense-loop-over-voxels`, reverses the traversal:

1. build a KD-tree over particle positions;
2. iterate output voxels;
3. query nearby particles for each voxel;
4. evaluate SPH contributions;
5. write one voxel result.

This can reduce atomic contention and is intended for workloads where a voxel-centric query is more efficient. It requires neighbor-radius/KD-tree support.

### 4.4 SPH updates

- SPH dispatch is available in host and device code.
- Kernel implementations were revised for CUDA compatibility.
- Dense radius and fill calculations were corrected.
- CUB-compatible kernel dispatch was added.
- Density normalization no longer depends on the removed `WITH_NO_DATA_TEMP` compile option.
- The secondary dense buffer is allocated at runtime only when the selected normalization requires it.

### 4.5 Fixed radius behavior

The former `particle_fix_size` concept is split into:

- `particle_radius_multiplier`, a scale applied to effective particle radius;
- `particle_radius_const`, a fixed radius selected with `--radius-const`.

This distinguishes scaling native/derived radii from replacing them.

### 4.6 Min/max corrections

Min/max calculation was revised across:

- CPU dense grids;
- GPU dense grids;
- OpenVDB managers;
- NanoVDB managers;
- CUB sparse data;
- remote result metadata.

The corrections ensure that final reduced values represent the actual converted output rather than stale local or pre-normalization ranges.

## 5. Data Model and Output Formats

### 5.1 Particle output selection

`SpaceData` adds `ExtractedParticleType`:

```text
0 None
1 RAW
2 VDB
3 VTP
```

The `--raw-particles` option now takes an integer output-format argument rather than acting only as a Boolean mode switch.

### 5.2 Dense output selection

`SpaceData` adds `ExtractedDenseType`:

```text
0 None
1 RAW
2 VTI
```

`--dense-file X` selects a dense file representation and replaces the previous `--dense2file` Boolean option.

### 5.3 VTK particle output

The new VTP path writes raw particles as VTK PolyData:

- positions become VTK points;
- scalar and vector particle attributes become point-data arrays;
- component counts are preserved;
- output is suitable for ParaView and other VTK consumers.

### 5.4 VTK dense output

The new VTI path writes dense results as VTK ImageData:

- dimensions match the dense grid;
- origin and spacing derive from the conversion bbox/transform;
- density/value data is written as image scalars;
- normalized output is finalized before export.

Both VTP and VTI require VTK support in the build.

### 5.5 Raw particle serialization

Raw-particle structures were moved into `raw_common.h`. Serialization and merge behavior remain field-oriented, while ownership is separated from OpenVDB/NanoVDB declarations.

## 6. New Input Formats

### 6.1 iPIC3D HDF5

New data type:

```text
IPIC3D_HDF5
```

New CLI arguments:

```text
--hdf5-file FILE
--num-files COUNT
--settings-file FILE
```

The reader supports up to four particle species:

```text
Species_0
Species_1
Species_2
Species_3
```

Supported particle and grid blocks include:

```text
Position
Velocity
Charge
ID
Electric field
Magnetic field
Charge density
Pxx
Pxy
Pxz
Pyy
Pyz
Pzz
```

Key behavior:

- loads all available species;
- discovers the latest cycle;
- supports multiple restart files;
- loads grid geometry from a companion `settings.hdf`, defaulting to the file beside the restart input;
- maps grid fields and moments to the common converter interface;
- publishes only actually available species/block combinations;
- distributes loaded data across MPI ranks;
- supports both particle and grid-derived fields through the common extraction workflow.

Build requirements are controlled by `WITH_IPIC3D`; HDF5 discovery is part of this module.

### 6.2 PLUTO VTK

New data type:

```text
PLUTO_VTK
```

New CLI arguments:

```text
--vtk-file FILE
--scalar-names NAME...
```

The reader:

- parses PLUTO legacy VTK rectilinear-grid output;
- reads coordinate arrays and cell-centered scalar fields;
- exposes one logical particle/cell type named `PLUTO`;
- dynamically exposes loaded scalar names as blocks;
- treats cell centers as common-interface positions;
- supports optional scalar-name filtering;
- partitions cells for MPI conversion.

The module is controlled by `WITH_PLUTO` and `WITH_VTK`.

## 7. Input Module Reorganization

### 7.1 Unified HACC module

GenericIO and HACC binary readers moved into `src/hacc`:

```text
src/hacc/hacc_genericio_*
src/hacc/haccbin_*
src/hacc/utility/
```

Data-type names changed:

| Previous | Development branch |
|---|---|
| `GENERICIO` | `HACC_GENERICIO` |
| `HACCBIN` | `HACC_BIN` |

The build option is now `WITH_HACC`. GenericIO, SZ, and related bundled sources remain under the new HACC module.

### 7.2 Removed generic readers

The following generic inputs were removed:

- `CSV`
- generic `HDF5`

Their old source directories and CMake targets are no longer built. HDF5 is now used specifically by the iPIC3D adapter.

Applications or scripts using `--data-type CSV`, `--csv-file`, or the old generic `HDF5` mode must migrate to a supported scientific format adapter.

### 7.3 GADGET readers

The GADGET and simplified GADGET readers received:

- corrected block reading;
- fixes to local particle access;
- improved reader cleanup;
- updated radius and cache integration;
- OpenMP fixes;
- new configuration definitions in the bundled OpenGadget3 layer;
- support for the shared namespace and converter interfaces.

`GADGET_MAX_HSML` is no longer a user-facing root build option and is disabled in updated platform configurations.

### 7.4 ChaNGa readers

Tipsy and NChilada adapters were migrated to the shared namespace and manager interfaces. XDR size handling was corrected, and sparse GPU compatibility fixes were applied without changing the high-level input names.

## 8. Command-Line Changes

### 8.1 Added options

```text
--cub
--dense-file X
--radius-const X
--dense-loop-over-voxels
--gpu
--skip-cache-manager
--sort-by-radius
--sort-by-non-overlap
--num-files X
--settings-file X
--vtk-file X
--scalar-names [names...]
```

### 8.2 Changed options

| Previous behavior | New behavior |
|---|---|
| `--raw-particles` flag | `--raw-particles X` selects RAW/VDB/VTP |
| `--dense2file` flag | `--dense-file X` selects RAW/VTI |
| `--export-data TYPE BLOCK` | `--export` is also accepted as an alias |
| `--particle-fix-size` semantics through remote state | Radius multiplier plus optional fixed radius |

### 8.3 Removed or inactive options

- `--rawpart2vdb` is no longer an active parser option; use the particle output selector.
- `--multires` and the associated `WITH_MULTIRES` feature were removed from the active build.
- `WITH_NO_DATA_TEMP` was removed; temporary normalization allocation is runtime-controlled.
- The CSV-specific and generic-HDF5 CLI groups were removed.

### 8.4 Configuration diagnostics

`FromCL` and `SpaceData` now provide detailed diagnostic printing of selected modes, radii, output types, bbox values, GPU flags, sorting flags, and animation state.

## 9. MPI and Communication

### 9.1 Large-transfer behavior

Chunked MPI send, receive, and float reduction remain in place, with clearer contracts and diagnostics. Dense and sparse reduction now dispatch through the selected manager and memory location.

### 9.2 GPU reductions

Dense GPU output can be reduced:

- directly from device memory with CUDA-aware MPI; or
- through synchronized host buffers when CUDA-aware MPI is disabled.

Sparse GPU data is serialized or compacted through the CUB manager before tree merging.

### 9.3 Remote request field change

The request field formerly named `particle_fix_size` is now `particle_radius_multiplier`. It remains a 32-bit float in the same position in the remote request stream, preserving wire layout while changing semantics.

### 9.4 New transport file type

The response file-type enum adds:

```text
5 = CUB
```

Clients that validate file types must be updated before requesting `--cub`. Existing identifiers 0 through 4 remain unchanged.

### 9.5 Improved remote behavior

The development work includes:

- corrected remote min/max values;
- corrected result type dispatch;
- clearer rank-0 ownership;
- additional request-field broadcasts;
- improved connection failure handling;
- reduced noisy message logging;
- updated raw-particle handling in BSpace.

## 10. BSpace Add-on

### 10.1 Version requirements

| Property | Main | Development |
|---|---:|---:|
| BSpace version | 0.3.0 | 1.1.0 |
| Minimum Blender version | 4.2.0 | 5.0.0 |

Copyright headers were updated through 2026.

### 10.2 Point-cloud visualization mode

Particle extraction now offers:

```text
Geometry Nodes
Point Cloud
```

The Point Cloud path:

- creates a native Blender point-cloud datablock;
- resizes it to the received particle count;
- writes the built-in `position` attribute;
- creates scalar `FLOAT` and vector `FLOAT_VECTOR` point attributes;
- links the point cloud into the active collection;
- creates a material driven by a selected particle attribute.

This avoids mesh-vertex overhead and provides a direct Blender-native representation for large point sets.

### 10.3 Geometry Nodes update

The Geometry Nodes particle path now:

- builds a UV sphere inside the node graph;
- instances it on received points;
- uses the `radius` named attribute for instance scale when available;
- assigns the generated particle material;
- realizes instances before output;
- names the modifier and node group after the object.

The workflow no longer depends on a separately created instance object.

### 10.4 Particle material

Particle objects receive a generated node material with:

- a named attribute input;
- Map Range normalization;
- a blue/green/red color ramp;
- Principled BSDF base color and emission;
- alpha and emission-strength control.

The post-import update path applies reduced min/max values to the Map Range node and the UI density value to the material multiplier where those nodes exist.

### 10.5 Radius control

The UI property is renamed to `particle_radius_multiplier` and labeled “Radius Multiplier.” The value is sent in the existing float slot and stored on imported objects as `PARTICLE_RADIUS_MULTIPLIER`.

### 10.6 Active collection behavior

Generated particle objects and helper geometry are linked to the active collection instead of always creating a top-level `SPACE` collection.

### 10.7 Raw-particle handling

Remote raw-particle responses now route through the selected Geometry Nodes or Point Cloud path. Imported objects are consistently named and receive the same metadata/update handling as volume results where applicable.

### 10.8 Interactive bbox object

When extraction produces no imported volume payload, BSpace can create or update an edge-only `BSPACE_INTERACTIVE_BBOX` mesh:

- vertices correspond to the configured bbox size;
- only cube edges are created;
- the object is reused when present;
- vertex positions are updated on subsequent operations;
- the object is linked to the active collection.

### 10.9 Error handling

Extraction operators now catch and print exceptions rather than allowing an unhandled Blender operator failure. Volume-only operations are guarded so they do not run against raw particle objects.

## 11. Logging and Instrumentation

A common logging layer was added:

```text
LOG_Init
LOG_MeasureStart
LOG_MeasureStop
```

Measured stages include:

- total execution;
- reader initialization;
- position caching;
- radius calculation;
- voxel KD-tree construction;
- bbox calculation;
- grid creation;
- conversion;
- min/max reduction;
- MPI reduction/merge;
- finalization;
- VDB/raw/VTI/VTP output.

The logging layer integrates with MERIC when enabled and provides a consistent no-op/timing interface otherwise.

## 12. CMake and Packaging

### 12.1 Installable package

The build now exports an install target set:

```text
space_converter_targets
```

It also installs:

- `space_converterConfig.cmake`;
- `space_converterConfigVersion.cmake`;
- exported target definitions;
- public common headers;
- libraries and executables with explicit archive/library/runtime destinations.

Consumers can use:

```cmake
find_package(space_converter CONFIG REQUIRED)
```

The package configuration discovers required transitive dependencies according to the compiled features.

### 12.2 Build and install interfaces

Target include paths now distinguish build-tree and install-tree use through generator expressions. Optional utility targets are included in the export set.

### 12.3 Debug/release dependency selection

OpenVDB and TBB link selection now uses configuration-aware generator expressions, allowing distinct debug and release library variables.

### 12.4 Dependency changes

- Boost link references were removed from active target definitions.
- Embree support was removed from the active root configuration.
- VTK dependency discovery was added.
- HDF5 discovery moved to iPIC3D support.
- HACC controls GenericIO/Blosc integration.
- OpenVDB was advanced from submodule commit `f564d35` to `84fc1e6`, bringing newer OpenVDB/NanoVDB APIs and fixes.
- `submodules/braas-hpc-renderengine` was added at commit `260a7f0`.

## 13. Platform Script Updates

Aurora, Barbora, CS P06, Karolina, Leonardo, LUMI, macOS, and Polaris profiles were updated to match the new options and link structure.

Common changes include:

- replacing `WITH_GENERICIO` with `WITH_HACC`;
- removing `WITH_NO_DATA_TEMP`;
- disabling/removing `GADGET_MAX_HSML`;
- adding GPU conversion flags on supported systems;
- updating OpenVDB/TBB library paths;
- updating runtime examples to new data-type names;
- minor MPI/run-command corrections.

The Leonardo profile received the largest platform-specific update, including dependency and configuration-path changes.

## 14. Utility Changes

### 14.1 OpenVDB/NanoVDB utilities

The existing utilities remain, with compatibility fixes for the updated OpenVDB/NanoVDB revision:

- `nano2vdb`
- `vdb2nano`
- `vdb2histo`
- `vdb2png`
- `vdb_merge`
- `nano2float3`

Updates include revised NanoVDB conversion API use, filtering/min-max fixes, configuration-aware linking, namespace fixes, and warning cleanup.

### 14.2 KD-tree utility

`cudakdtree_tool` now provides:

- revised GPU and host KNN paths;
- updated MPI cycling;
- additional voxel-query tree support;
- shared CPU/GPU dispatch;
- improved radius/density handling;
- corrected buffer and tree lifecycle management.

The nanoflann interface was updated to match the reorganized namespace and cache model.

## 15. Removed Features and Files

The development branch removes:

- the generic CSV adapter and build target;
- the generic HDF5 adapter and build target;
- the old standalone TCP implementation under `src/common`;
- the monolithic embedded grid representations;
- active Embree integration;
- active multiresolution conversion code;
- the `WITH_NO_DATA_TEMP` compile-time path;
- the old separate `genericio` and `haccbin` directory layout.

Most GenericIO/SZ sources are moves into `src/hacc`, not deletions.

## 16. Compatibility and Migration

### 16.1 Required build migration

Replace:

```cmake
-DWITH_GENERICIO=ON
```

with:

```cmake
-DWITH_HACC=ON
```

Enable the new format adapters explicitly:

```cmake
-DWITH_IPIC3D=ON
-DWITH_PLUTO=ON
```

`WITH_IPIC3D` automatically enables the required HDF5 dependency, and
`WITH_PLUTO` automatically enables the required VTK dependency. `WITH_HDF5`
and `WITH_VTK` remain available as lower-level dependency switches.

For complete GPU conversion:

```cmake
-DWITH_GPU_CUDA=ON
```

or, for AMD GPUs (e.g. MI250X / gfx90a on LUMI-G):

```cmake
-DWITH_GPU_HIP=ON -DCMAKE_HIP_ARCHITECTURES=gfx90a
```

Add `WITH_CUDA_AWARE_MPI` / `WITH_HIP_AWARE_MPI` only when the MPI installation explicitly
supports device buffers. On LUMI that means Cray MPICH with the GTL library (loaded by
`partition/G`) and `MPICH_GPU_SUPPORT_ENABLED=1` in the job environment.

Ready-made LUMI scripts live in `scripts/space_converter/lumi/`:
`build_lumi.sh` (LUMI/25.09, partition/G, PrgEnv-gnu, rocm/6.4.4), `run_lumi.sh`, and
`gpu.sh`, which pins one MPI rank per GCD via `ROCR_VISIBLE_DEVICES`.

### 16.2 Required CLI migration

```text
GENERICIO  -> HACC_GENERICIO
HACCBIN    -> HACC_BIN
--dense2file -> --dense-file <format-id>
--raw-particles -> --raw-particles <format-id>
```

Remove uses of:

```text
CSV
HDF5
--csv-file
--rawpart2vdb
--multires
```

Use `IPIC3D_HDF5` for iPIC3D restart data.

### 16.3 Blender migration

- Upgrade to Blender 5.0 or newer.
- Reinstall or reload BSpace 1.1.0.
- Reconfigure add-on executable paths if the installation prefix changed.
- Choose Geometry Nodes or Point Cloud for particle extraction.
- Revisit particle-size values because the control now acts as a radius multiplier.

### 16.4 C++ API migration

- Update namespaces to `space_converter::common`.
- Replace direct `DenseParticles`/embedded grid access with `VoxelDenseManager`.
- Replace direct NanoVDB/OpenVDB members on `VDBParticles` with `sparse_grid`.
- Access cached particle arrays through `ConvertVDBBase::cache_manager`.
- Implement new virtual interfaces when maintaining an external adapter.

### 16.5 Protocol compatibility

Existing message IDs and file type IDs 0–4 are retained. The particle-size float remains in the same request position but now represents a multiplier.

Clients must be updated to understand:

- CUB file type 5;
- the changed raw-particle format-selection UI/CLI;
- any new output path extensions such as VTP and VTI.

## 17. Known Operational Considerations

- CUDA-aware MPI is environment-dependent and must not be enabled on MPI stacks without GPU-buffer support.
- CUB output requires a client or postprocessor that understands file type 5 and the serialized CUB layout.
- The BSpace release requires Blender 5.0, so it is not a drop-in update for Blender 4.2 installations.
- CSV and generic HDF5 users need a format-specific migration path.
- PLUTO VTK support expects the legacy rectilinear-grid organization implemented by the reader.
- iPIC3D loading depends on correct restart/settings file organization and HDF5 availability.
- GPU results should be validated against the CPU path within floating-point tolerance because reduction order differs.
- No dedicated repository-level automated test suite was added in this branch; validation currently depends on builds, sample runs, and CPU/GPU comparison workflows.

## 18. Recommended Release Validation

Before promoting the development branch:

1. Build a CPU-only OpenVDB configuration.
2. Build a NanoVDB-only configuration.
3. Build CUDA conversion with and without CUDA-aware MPI.
4. Run dense and sparse GADGET conversions on CPU and GPU and compare voxel counts, min/max, and aggregate values.
5. Exercise CUB serialization and reduction with more than one MPI rank.
6. Validate both radius and non-overlap sorting.
7. Test the voxel-centric dense loop with KD-tree radius calculation.
8. Load representative HACC GenericIO and HACC binary files under their new data-type names.
9. Validate all iPIC3D species and grid-field discovery.
10. Validate PLUTO scalar discovery and coordinate mapping.
11. Export RAW, VDB, and VTP particles.
12. Export RAW and VTI dense volumes.
13. Connect BSpace 1.1.0 from Blender 5.0 and test volume, Geometry Nodes, and Point Cloud extraction.
14. Test remote min/max updates and the interactive bbox object.
15. Install the CMake package and compile a minimal external `find_package` consumer.
16. Run all updated HPC build and launch profiles on their target platforms.

## 19. Commit Themes

The 86 development commits fall into these broad groups:

- initial cleanup and removal of obsolete CSV/generic-HDF5 paths;
- CMake, Zlib, OpenVDB, and package-link fixes;
- cache-manager introduction;
- common conversion-kernel extraction;
- GPU dense and sparse manager implementation;
- MPI and sparse-reduction fixes;
- radius and spatial sorting;
- raw-particle and Blender point-cloud work;
- iPIC3D and PLUTO adapters;
- OpenVDB/NanoVDB manager refinements;
- voxel-centric dense conversion;
- logging and profiling;
- CUB representation and serialization;
- cross-platform OpenMP, macOS, Linux, and HPC fixes;
- final remote, min/max, transform, and reader corrections.

This release should be treated as a major development milestone rather than a small maintenance update because it changes build options, internal APIs, supported input names, output selection, Blender requirements, and execution architecture.
