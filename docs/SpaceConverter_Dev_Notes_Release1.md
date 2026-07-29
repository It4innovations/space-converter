# SpaceConverter Development Specification (Release1)

## 1. Product Definition

SpaceConverter is a high-performance scientific data conversion system for astronomy and cosmology simulations. It reads distributed particle snapshots, selects particle families and data fields, optionally reconstructs continuous scalar fields with Smoothed Particle Hydrodynamics (SPH) kernels, and produces OpenVDB, NanoVDB, dense raw-volume, or raw-particle results suitable for Blender and Cycles.

The product consists of:

1. A C++17 command-line application named `space_converter`.
2. Format-specific reader and adapter libraries.
3. MPI and OpenMP processing infrastructure, with optional CUDA acceleration.
4. A binary TCP service used by a Blender add-on.
5. The Blender add-on named BSpace.
6. Standalone VDB/NanoVDB inspection and conversion utilities.
7. CMake build definitions and platform-specific HPC build/run scripts.

The implementation shall preserve support for local batch conversion, interactive remote extraction, MPI aggregation, animation-frame distribution, CPU and GPU nearest-neighbor radius calculation, multiple particle file formats, and Blender-side volume or particle visualization.

## 2. Scope and Quality Goals

### 2.1 Functional goals

The system shall:

- Read GADGET/OpenGadget3, simplified GADGET, ChaNGa Tipsy, ChaNGa NChilada, CSV, GenericIO, HDF5, and HACC binary data when their corresponding build options are enabled.
- Discover available particle types and data blocks at runtime.
- Select a particle type and block by integer identifier.
- Compute the global data-space bounding box with MPI.
- Map requested world-space subregions into simulation coordinates.
- Export sparse voxels, dense SPH reconstructions, or raw particles.
- Support scalar and three-component particle fields.
- Support OpenVDB and NanoVDB output.
- Merge per-rank grids without losing overlapping contributions.
- Transfer results to Blender or save them to disk.
- Support a sequence of snapshots by distributing frames across MPI ranks.
- Provide Blender material, volume, particle-instancing, histogram, slice, and bounding-box workflows.

### 2.2 Non-functional goals

- Use C++17 for the converter and Python with Blender's `bpy` API for BSpace.
- Scale particle loading and conversion across MPI ranks.
- Use OpenMP for local parallel loops and timing.
- Chunk very large MPI transfers to avoid count-size limitations.
- Keep format readers behind one polymorphic converter interface.
- Use native-endian binary TCP fields because both existing endpoints use native `struct` packing.
- Treat rank 0 as the TCP endpoint and final-result owner.
- Remain buildable on Linux HPC systems, macOS, and Windows where dependencies are available.
- Preserve GPL-3.0-or-later licensing notices.

## 3. Repository Layout

The required source tree is:

```text
/
├── CMakeLists.txt
├── README.md
├── LICENSE
├── .gitmodules
├── blender_addons/
│   └── bspace/
│       ├── __init__.py
│       ├── bspace_pref.py
│       ├── bspace_remote.py
│       └── bspace_panel.py
├── docs/
│   └── bspace_menu.png
├── scripts/
│   ├── util/draw_raw.py
│   └── space_converter/
│       ├── aurora/
│       ├── barbora/
│       ├── cs/
│       ├── karolina/
│       ├── leonardo/
│       ├── lumi/
│       ├── macos/
│       └── polaris/
├── src/
│   ├── main.cpp
│   ├── args_processing.{h,cpp}
│   ├── data_processing.{h,cpp}
│   ├── data_communication.{h,cpp}
│   ├── common/
│   ├── utility/
│   ├── gadget/
│   ├── changa/
│   ├── csv/
│   ├── genericio/
│   ├── hdf5/
│   └── haccbin/
└── submodules/
    ├── openvdb/
    ├── nanoflann/
    ├── cudaKDTree/
    └── c-blosc/
```

Vendored reader support under `src/changa/utility`, `src/gadget/utility`, and `src/genericio/utility` shall be retained. These directories provide ChaNGa structures/XDR support, OpenGadget3 I/O, GenericIO, and SZ compression support respectively.

## 4. Build System

### 4.1 Root CMake configuration

- Require CMake 3.20 or newer.
- Define project `space_converter`.
- Require C++17.
- Default `WITH_MPI`, `WITH_OPENMP`, and `WITH_NANOVDB` to enabled.
- Expose the following options:

| Option | Default | Meaning |
|---|---:|---|
| `WITH_OPENVDB` | ON | Enable OpenVDB grids and tools |
| `WITH_HDF5` | OFF | Enable the HDF5 adapter |
| `WITH_GENERICIO` | OFF | Enable GenericIO and bundled Blosc |
| `WITH_EMBREE` | OFF | Enable Embree-based spatial support |
| `GADGET_MAX_HSML` | OFF | Limit loaded GADGET smoothing-length values |
| `GADGET_READ_ID` | OFF | Read 32-bit GADGET particle IDs |
| `GADGET_READ_ID64` | OFF | Read 64-bit GADGET particle IDs |
| `WITH_NANOFLANN` | OFF | Enable CPU KD-tree neighbor search |
| `WITH_CUDAKDTREE` | OFF | Enable CUDA/host cudaKDTree search |
| `WITH_MULTIRES` | OFF | Enable experimental multiresolution VDB generation |
| `WITH_NO_DATA_TEMP` | OFF | Remove the secondary dense accumulation array |
| `WITH_MERIC` | OFF | Enable MERIC instrumentation |

- Find OpenMP and append its compile and linker flags. On Windows use `/openmp:llvm` and `/bigobj`.
- Accept dependency paths through CMake cache variables such as `OPENVDB_INCLUDE_DIRS`, `OPENVDB_LIBRARIES`, `OPENVDB_LIBRARIES_DEBUG`, `OPENVDB_VERSION`, `NANOVDB_INCLUDE_DIRS`, `TBB_INCLUDE_DIRS`, `TBB_LIBRARIES`, `MPI_CUSTOM_LIBRARIES`, `MERIC_INCLUDE_DIRS`, and `MERIC_LIBRARIES`.
- If OpenVDB is disabled, use the NanoVDB headers from the OpenVDB submodule and default `OPENVDB_VERSION` to 12.
- If CUDA KD-tree support is enabled, enable the CUDA language and select either the user-supplied architecture list, `70;80` for older CMake, or `all-major`.
- When GenericIO is enabled, build `c-blosc` statically with tests and benchmarks disabled.

### 4.2 Targets

Create these libraries:

- `space_common`
- `gadget_converter`
- `gadget_simple_converter`
- `changa_tipsy_converter`
- `changa_nchilada_converter`
- `csv_converter`
- `haccbin_converter`
- `genericio_converter` when enabled
- `hdf5_converter` when enabled
- `nanoflann_tool` when enabled
- `cudakdtree_tool` when enabled

Create these executables:

- `space_converter`
- `vdb_merge`
- `vdb2histo`
- `vdb2png`
- `vdb2nano`
- `nano2vdb`
- `nano2float3`

Install the main executable and utility executables to `bin`; install libraries to `lib`.

Propagate feature macros (`WITH_OPENMP`, `WITH_MPI`, `WITH_OPENVDB`, `WITH_NANOVDB`, `WITH_TBB`, `WITH_EMBREE`, `WITH_GENERICIO`, `WITH_NANOFLANN`, `WITH_CUDAKDTREE`, `WITH_MULTIRES`, `WITH_NO_DATA_TEMP`, `WITH_MERIC`) to the targets that require them. Define `OPENVDB_VERSION` wherever OpenVDB or NanoVDB APIs are compiled.

### 4.3 External dependencies

- MPI
- OpenMP
- OpenVDB and its TBB dependency
- NanoVDB headers
- Optional CUDA toolkit
- Optional nanoflann
- Optional Embree
- Optional HDF5
- GenericIO and Blosc for GenericIO builds
- Blender 4.2 or newer with Python, NumPy, `bpy`, `gpu`, and `mathutils`

The `.gitmodules` file shall point to the Academy Software Foundation OpenVDB repository, jlblancoc/nanoflann, jar091/cudaKDTree, and Blosc/c-blosc.

## 5. Core Data Model

### 5.1 `FromCL`

Define `space_converter::FromCL` with:

- `data_type = "GADGET"`
- `output_path`
- `server = "localhost"`
- `port = 7000`
- `info = false`
- `remote = true`
- `use_nanovdb`: false when OpenVDB is compiled, true otherwise
- `use_save_mpirank = false`
- `use_rawpart2vdb = false`
- `use_dense2file = false`
- optional `use_cudakdtree`, `use_cudakdtree_cpu`, and `use_nanoflann`
- `use_multires = false`
- `world_rank = 0`
- `world_size = 1`

### 5.2 `SpaceData`

Define `common::SpaceData` as the complete mutable request and result state.

Enums and stable integer values:

```text
MessageType:   Exit=-1, Empty=0, Info=1, Data=2, BBOX=3
DenseType:     None=0, Cubic=1, Quintic=2, WendlandC2=3,
               WendlandC4=4, WendlandC6=5, WendlandC8=6
DenseNorm:     None=0, Count=1, SPHInterpolation=2
ExtractedType: Sparse=0, Dense=1, Particle=2
AnimType:      None=0, AllPath=1, AllMerge=2, FrameExtract=3
```

Required fields and defaults:

| Group | Fields |
|---|---|
| Request identity | `message_type=Empty`, `particle_type=0`, `block_name_id=0` |
| Grid | `grid_transform=1`, `bbox_dim=100`, `object_size=1000` |
| Modes | `extracted_type=Sparse`, `dense_type=None`, `dense_norm=None` |
| Values | `min_value=0`, `max_value=1`, `min_rho=0`, `max_rho=1`, reduced min/max |
| Filters | `particle_fix_size=0`, `filter_min=-FLT_MAX`, `filter_max=FLT_MAX` |
| Counts | global/local particle counts and voxel count |
| Display bbox | `bbox_min={0,0,0}`, `bbox_max={1000,1000,1000}` |
| Integer bbox | global and local integer minima/maxima, three elements each |
| Transform | `bbox_size_orig`, `transform_scale` |
| Animation | frame, type, full file path, task counter, start/end/step |
| Radius search | neighbor count, default density kernel Wendland C6, cache filename |
| Selection | spherical bbox flag/center/radius |
| Behavior | simple-density flag, normalize-value flag, position offset |

### 5.3 Grid and particle containers

`ParticleData` contains a string `name`, integer `num_comp`, and `std::vector<float> values`.

`RawParticles` contains `std::vector<ParticleData> data`. It shall:

- Serialize to a file or byte vector.
- Encode `size_t` counts in native representation.
- For every field encode name length, name bytes, component count, value count, and contiguous float values.
- Deserialize the same format.
- Merge another container by field name, appending values for matching fields.

`DenseParticles` contains:

- `data_density`
- optional `data_temp`
- `dims[3]`
- `offset[3]`

It shall provide `create`, `clear`, `x`, `y`, `z`, `size`, `memsize`, and x-major indexing:

```text
index = x + y * dimX + z * dimX * dimY
```

`VDBParticles` is a tagged container. Stable types are:

```text
Dense, Vector, NanoVDB, OpenVDB, RawParticles
```

It owns the dense grid, scalar and vector NanoVDB build grids, optional OpenVDB float grid, serialized bytes, and raw particles.

## 6. Converter Abstraction

Implement `common::vdb::ConvertVDBBase` as the sole interface consumed by the processing pipeline.

Every concrete adapter shall implement:

- lifecycle: `init_lib`, `finish_lib`
- counts: local/global particle count, number of types, number of blocks
- discovery: types-and-blocks matrix and printable names
- particle access: type, position, field component count, raw field values, normalized field value
- physical access: smoothing length, mass, density, density block identifier
- diagnostics: CPU step reporting

The base class shall provide:

- Bounding-box and min/max traversal.
- Type-aware wrappers that reject or normalize unsuitable accesses.
- Particle-to-grid conversion.
- Grid merging.
- Dense-to-OpenVDB and dense-to-NanoVDB conversion.
- OpenVDB stream serialization/deserialization.
- Neighbor-derived radius storage per particle type.
- Optional CUDA KD-tree and nanoflann radius computation.
- Radius cache read/write.
- Animation filename expansion by replacing every `{}` token.

The format metadata availability vector is block-major: element `numTypes * block + type` indicates whether that type/block combination is present.

## 7. Command-Line Interface

The executable syntax is:

```text
space_converter --data-type TYPE [common options] [format options]
```

Supported `TYPE` values:

```text
GADGET
GADGET_SIMPLE
CHANGA_TIPSY
CHANGA_NCHILADA
CSV
GENERICIO
HDF5
HACCBIN
```

Common options:

| Option | Arguments | Effect |
|---|---|---|
| `--data-type` | type | Select adapter |
| `--grid-dim` | integer | Cubic output resolution |
| `-o`, `--output-path` | path | Output directory |
| `--server` | host | TCP host |
| `--port` | integer | TCP port |
| `--info` | none | Print metadata and disable remote mode |
| `--nanovdb` | none | Prefer NanoVDB |
| `--dense2file` | none | Save dense raw arrays |
| `--anim` | start end step | Enable sequence processing |
| `--anim-merge` | none | Merge sequence ranks into one result |
| `--raw-particles` | none | Select particle extraction |
| `--save-mpi-rank` | none | Save each rank before reduction |
| `--rawpart2vdb` | none | Convert extracted particles to VDB points |
| `--export-data` | type block | Select data and disable remote mode |
| `--dense-type` | integer | Select kernel; 0 means sparse |
| `--dense-norm` | integer | Select normalization |
| `--bbox-sphere` | x y z radius | Enable spherical selection |
| `--bbox` | minX minY minZ maxX maxY maxZ | Select display-space bbox |
| `--simple-density` | none | Dense conversion using unit particle values |
| `--no-norm-value` | none | Use field components instead of magnitude |
| `--offset-position` | x y z | Translate particle positions |
| `--calc-radius-neigh` | count | Derive support radius from neighbors |
| `--calc-radius-neigh-rho-kernel` | integer | Kernel for neighbor-derived density |
| `--calc-radius-neigh-file` | path | Read/write radius cache |
| `--cudakdtree` | none | Use CUDA KD-tree |
| `--cudakdtree-cpu` | none | Use cudaKDTree host implementation |
| `--nanoflann` | none | Use nanoflann |
| `--multires` | none | Enable multiresolution output |
| `-h`, `--help` | none | Print usage and exit |

Format-specific options:

- GADGET: `--param-file`, `--gadget-file`, `--max-mem-size`, `--buffer-size`, `--part-alloc-factor`, `--bh-count`.
- GADGET_SIMPLE: `--gadget-file`.
- CHANGA_TIPSY: `--tipsy-file`, optional `--filter-in`.
- CHANGA_NCHILADA: `--nc-dir`.
- CSV: `--csv-file`.
- GENERICIO: `--genericio-file`, `--pos-names x y z`, `--vel-names vx vy vz`, `--mass-name`, `--rho-name`, `--hsml-name`.
- HDF5: `--hdf5-file`.
- HACCBIN: `--haccbin-file`.

Missing values and malformed numeric values should produce a clear error rather than undefined access in a hardened implementation.

## 8. Main Processing Lifecycle

The executable shall perform these steps in order:

1. Initialize MPI and populate rank/size.
2. Optionally read `CONVERTER_LARGE_TEST_FILE` on every rank as an I/O stress hook.
3. Initialize OpenVDB when compiled.
4. Parse common arguments into `FromCL` and `SpaceData`.
5. Instantiate the selected adapter and load the local particle partition.
6. Gather and print the global type/block availability matrix.
7. Enter the connection lifecycle:
   1. Rank 0 opens the TCP connection when remote.
   2. Compute the initial global bbox.
   3. Repeatedly wait for a message.
   4. On `Info`, send animation metadata and available fields.
   5. On `BBOX`, receive type/block, compute that type's bbox, and return it in display coordinates.
   6. On `Data`, receive conversion parameters and execute extraction.
   7. On `Exit` or socket error, leave the request loop.
   8. Close TCP; remote mode permits a later reconnect.
8. Destroy the adapter.
9. Uninitialize OpenVDB.
10. Finalize MPI and MERIC.

For a data request:

1. Receive and broadcast request fields.
2. Create the appropriate local grid/container.
3. Convert the local particle partition.
4. MPI-reduce min, max, and particle count.
5. Optionally save each rank.
6. Reduce or merge grid data.
7. Finalize and serialize the result.
8. Reduce final min/max across animation ranks where required.
9. Send bytes directly for non-path results, or save and return a path.

## 9. Coordinate and Bounding-Box Rules

- Reader positions are simulation-space doubles.
- `iolib_find_bbox` shall scan only the requested particle type, applying `offset_position`.
- Local floating extrema shall be converted to integer extrema using floor/ceil semantics, then reduced with MPI min/max.
- `bbox_size_orig` is the largest of the three global axis lengths.
- Display-space requests use `[0, object_size]` coordinates relative to the global simulation bbox.
- The requested integer sub-bbox shall be derived independently for each axis.
- `transform_scale` maps grid indices into the requested object-space size and includes `grid_transform`.
- A per-grid integer offset shall preserve the requested subregion's placement.
- Optional spherical filtering evaluates distance from `bbox_sphere_pos` and rejects particles outside `bbox_sphere_r`.
- Position offsets apply before bbox and voxel calculations.

When returning a type-specific bbox to Blender, map it from simulation coordinates into the full dataset's display coordinate system:

```text
displayMin[axis] =
  (typeMin[axis] - globalMin[axis]) * objectSize / globalExtent[axis]

displayMax[axis] =
  (typeMax[axis] - globalMin[axis]) * objectSize / globalExtent[axis]
```

## 10. Particle-to-Grid Conversion

### 10.1 Common filtering

For each local particle:

- Require `get_particle_type(id) == requested particle_type`.
- Read position and apply offsets.
- Reject positions outside the requested cuboid.
- Apply optional spherical rejection.
- Read the requested field:
  - normalized mode uses scalar value or vector magnitude;
  - non-normalized mode retains components, allowing vector sparse/particle output.
- Reject values outside `[filter_min, filter_max]`.
- Track local min/max and accepted count.

### 10.2 Sparse mode

- Produce a scalar NanoVDB/OpenVDB grid for one-component or normalized fields.
- Produce a `Vec3f` NanoVDB grid when a three-component non-normalized value is requested.
- Map particle positions into integer voxel coordinates.
- Write or accumulate values at active voxels.
- Name scalar volume grids `density`.
- Set OpenVDB output to fog-volume class.

### 10.3 Raw-particle mode

- Store position as one three-component `ParticleData`.
- Store the selected block with its native component count.
- Preserve compatible named attributes for merging across ranks.
- Serialize the `RawParticles` container for TCP transfer.
- When `--rawpart2vdb` is active and OpenVDB is available, create an OpenVDB points representation suitable for Blender import.

### 10.4 Dense mode

Allocate `bbox_dim³` float cells. For every accepted particle:

- Obtain field value, mass, density, smoothing length, and effective support radius.
- Use `particle_fix_size` when nonzero; otherwise use reader smoothing length or the neighbor-derived radius.
- Convert physical position and radius to grid units.
- Restrict voxel iteration to the particle's support bounds and clip to the local dense grid.
- Evaluate the selected SPH kernel.
- Add the weighted contribution to `data_density`.
- When the secondary array exists, add the normalization denominator to `data_temp`.

Normalization behavior:

- `None`: emit accumulated weighted values.
- `Count`: normalize by accumulated contribution/count where defined.
- `SPHInterpolation`: use the SPH numerator/denominator relationship and physical particle factors.

Zero, NaN, and non-finite results shall not create active VDB voxels.

### 10.5 SPH kernels

Implement separate kernel classes for:

- Cubic spline
- Quintic spline
- Wendland C2
- Wendland C4
- Wendland C6
- Wendland C8

Each kernel shall expose:

- compact support
- normalization as a function of inverse smoothing length
- value evaluation
- optional neighbor-count bias correction

The shared dispatch function selects by `SpaceData::DenseType`. The Wendland kernels shall preserve their kernel-specific polynomial, three-dimensional normalization constant, and bias correction. Neighbor-derived density is:

```text
rho_i = Σ_j mass_i * W(distance(i,j) / h_i, 1 / h_i)
```

followed by the selected kernel's finite-neighbor bias correction.

### 10.6 Dense-to-VDB conversion

- Create a scalar float grid with background zero.
- Apply the linear transform from `transform_scale`.
- Apply dense-grid offset translation for OpenVDB.
- Iterate all dense cells, compute the selected normalized density, and activate finite nonzero values.
- Name the grid `density` and mark it as a fog volume.
- Support both OpenVDB and NanoVDB build APIs, including the namespace/API difference around OpenVDB version 11.

## 11. Radius Calculation

The base converter stores `radius_particles_per_ptype`, `rho_particles_per_ptype`, and per-type offsets.

Radius sources, in precedence order, are:

1. A valid cached radius file.
2. CUDA KD-tree when selected and compiled.
3. cudaKDTree host mode when selected.
4. nanoflann when selected and compiled.
5. Reader-provided smoothing length.
6. Fixed particle size when explicitly requested.

The K-nearest-neighbor workflows shall:

- Separate points by particle type.
- Preserve local/global type offsets.
- Use `K + 1` candidates so the query particle can be included.
- Support cyclic exchange of rank-local point clouds for distributed datasets.
- Transfer data in bounded MPI chunks.
- Keep the K smallest squared distances for each query.
- Set radius to the square root of the Kth selected distance.
- Optionally compute density from the selected neighbor set.
- Support a maximum-radius constraint in CUDA mode.

Radius cache files shall be split by particle type as:

```text
<base>.<particle-type>.bin
```

and contain enough contiguous float data to restore the per-type radius arrays.

## 12. MPI Semantics

- MPI is initialized unconditionally by the current executable design.
- Rank 0 prints startup and final timing information.
- Large point-to-point transfers are split into chunks no larger than 2,000,000,000 units.
- Dense grids are summed with chunked `MPI_Reduce` to rank 0.
- Local and global extrema use `MPI_Allreduce`.
- Availability matrices use collective aggregation so rank 0 reports all fields present in all loaded partitions.
- Sparse grids and raw particles use a logarithmic binary-tree reduction:
  - at each power-of-two step, receiver ranks read a serialized payload and merge it;
  - sender ranks serialize, transmit to `rank - step`, and leave the merge loop.
- Animation modes other than merged animation keep each rank's frame independent.
- When an animation pattern is active, each rank loads frame `anim_start + anim_step * rank`, then initializes the reader as a single-rank reader to prevent partitioning one frame across the animation ranks.

## 13. Grid Merge and Finalization

`merge_grid` shall accept a destination in native form and a received serialized vector.

- Raw particles: deserialize and merge attributes by name.
- OpenVDB: deserialize and combine scalar grids with component-wise sum.
- NanoVDB: read the grid handle and merge active values into the destination build grid.
- Vector grids: preserve vector component values and active topology.

Finalization shall:

- Convert a dense sum into the selected VDB representation.
- Serialize OpenVDB with `openvdb::io::Stream`.
- Serialize NanoVDB with a `GridHandle`.
- Serialize raw particles with `RawParticles::serialize`.
- Store the bytes in `vector_grid`.
- Compute final active-voxel min/max.
- Optionally invoke multiresolution generation.
- Preserve rank ownership: only rank 0 owns a merged non-animation result.

## 14. File Output

Use the selected field name, particle type, frame/rank context, and a collision-resistant time component to construct output names.

Supported output forms:

- `.vdb` for OpenVDB
- `.nvdb` for NanoVDB
- binary dense arrays when `--dense2file` is active
- serialized raw particle data
- per-rank output when `--save-mpi-rank` is active

For path-based animation, store the full saved path in `SpaceData::full_filepath` and return it to Blender. Path output allows Blender to construct a volume sequence without transferring every frame through one response body.

## 15. TCP Transport

### 15.1 Connection

- Use one IPv4 TCP data socket.
- Rank 0 is the converter-side endpoint.
- Default converter host/port are `localhost:7000`; BSpace defaults to `localhost:5000`.
- Support Windows Winsock and POSIX sockets.
- Disable Nagle where configured and enlarge socket buffers for bulk volume transfer.
- Implement exact-length send/receive loops.
- Treat zero-length reads, socket failures, and timeout failures as connection errors.
- Keep an internal error flag visible through `is_error()`.
- Support environment overrides `SOCKET_SERVER_PORT_DATA`, `SOCKET_SERVER_NAME_DATA`, and the corresponding camera-channel variables retained by the TCP class.

### 15.2 Primitive encoding

The established protocol uses native little-endian values on supported deployments:

- `int`: 4 bytes
- `float`: 4 bytes
- `double`: 8 bytes
- `size_t` payload length: 8 bytes on 64-bit builds
- float3: 12 contiguous bytes

All production endpoints shall be 64-bit and use the same endianness. A future protocol revision may add a version preamble and fixed-width network byte order, but compatibility mode shall preserve the existing layout.

### 15.3 Client-to-converter messages

Every request begins with a 4-byte `MessageType`.

`Info` has no request payload.

`BBOX` payload:

```text
int particle_type
int block_name_id
```

`Data` payload, in exact order:

```text
float bbox_min[3]
float bbox_max[3]
int   bbox_dim
float grid_transform
int   particle_type
int   block_name_id
int   extracted_type
int   dense_type
int   dense_norm
float object_size
float particle_fix_size
float filter_min
float filter_max
int   frame
int   anim_type
int   anim_task_counter
```

Rank 0 shall broadcast every received request field to all MPI ranks. The fixed particle size must also be broadcast; implementations shall include it even if an older code path omitted that broadcast.

### 15.4 Converter-to-client responses

`Info` response:

```text
int  anim_type
int  anim_start
int  anim_end
int  text_size
char text[text_size]
```

The text contains one line per available type/block pair:

```text
particle-name;particle-id;block-name;block-id\n
```

The client replies with a 4-byte acknowledgment.

`BBOX` response:

```text
float bbox_min[3]
float bbox_max[3]
```

followed by a 4-byte client acknowledgment.

`Data` response:

```text
int    file_type
size_t file_size
byte   file_data[file_size]
float  min_value
float  max_value
float  min_value_reduced
float  max_value_reduced
int    frames
```

Stable file types:

```text
0 None
1 OpenVDB
2 NanoVDB
3 Path
4 Raw particles
```

The client replies with a 4-byte acknowledgment. `frames` is one for a single or merged result, and equals MPI world size for independent animation paths/frames.

## 16. Input Adapters

### 16.1 General adapter rules

Each reader shall:

- Partition particles approximately evenly across MPI ranks unless one rank owns one animation frame.
- Keep local data in reader-specific structures.
- Return positions as double precision.
- Return field values as floats with explicit component count.
- Return vector magnitude for normalized access.
- Track whether every particle type/block pair exists.
- Provide stable human-readable type and block names.
- Return zero or an explicit fallback for absent mass, density, or smoothing-length fields.
- Release all loaded memory in `finish_lib`.

### 16.2 GADGET/OpenGadget3

- Support six standard GADGET particle types.
- Load via the bundled OpenGadget3 I/O library.
- Configure IC format 2, snapshot format 2, parallel file count, memory size, buffer size, allocation factor, and black-hole count.
- Initialize reader memory and physical units before reading.
- Expose the reader's native block enumeration and names.
- Default tuning: max memory 1000, buffer 100, allocation factor 1.2, black-hole count 1.
- Optional compile flags control smoothing-length limiting and particle ID width.

### 16.3 Simplified GADGET

- Read one-file or multi-file GADGET snapshots without the full OpenGadget3 application initialization.
- Discover file count and labeled blocks.
- Parse the GADGET header and distribute particles.
- Preserve block labels, component counts, and particle-type membership.
- Expose snapshot redshift and Hubble parameter to the base converter for physical SPH scaling.

### 16.4 ChaNGa Tipsy

Particle types:

```text
0 Gas
1 Dark
2 Stars
```

Canonical blocks:

```text
Position, Mass, Velocity, Softening, Potential,
Smoothing Length, Density, Temperature, Metals, Formation Time
```

- Use the bundled Tipsy structures/readers.
- Support partial particle-range reading for MPI partitioning.
- Accept an optional input filter.
- Map Tipsy gas, dark, and star structures to the canonical fields.

### 16.5 ChaNGa NChilada

- Use the bundled NChilada reader and XDR implementation.
- Discover attributes from the dataset directory.
- Support gas, dark, and star particle families.
- Preserve per-attribute component count and type.
- Use the same canonical semantic accessors as Tipsy.

### 16.6 CSV

- Parse a header row followed by numeric rows.
- Detect delimiter and field layout consistently.
- Store data in field-oriented arrays.
- Divide rows across MPI ranks.
- Provide canonical gas/dark/star and standard block identifiers where present.
- Parse scalar and delimited vector columns.
- Use field names for dataset display names.

### 16.7 GenericIO

Particle type:

```text
0 HACC
```

- Use bundled GenericIO.
- Enumerate available variables dynamically, including numeric type and component semantics.
- Allow explicit position and velocity field names.
- Allow explicit mass, density, and smoothing-length field names.
- Supply sensible position defaults compatible with HACC (`x`, `y`, `z`) when names are not supplied.
- Register fields with GenericIO, read the local section, and expose each variable as a block.
- Link Blosc and SZ support as required by the bundled reader.

### 16.8 HDF5

- Compile only under `WITH_HDF5`.
- Support gas, dark, and star types.
- Expose the canonical ten blocks used by CSV/Tipsy.
- Read configured HDF5 snapshots and distribute particle ranges.
- Preserve scalar/vector values and physical accessors.

### 16.9 HACC binary

Particle types and stable IDs:

```text
0 DarkMatter
1 Baryon
2 BaryonStar
3 BaryonWind
4 BaryonGas
5 DarkMatterAGN
```

Blocks:

```text
0 Pos   (Mpc/h)
1 Vel   (km/s)
2 Mass  (internal units)
3 UU    ((km/s)^2)
4 HH    (Mpc/h)
5 MU    (dimensionless)
6 Rho   (h^2 solar-mass / Mpc^3)
7 Phi   (internal units)
8 Id
9 Mask
```

- Read the HACC binary header and particle records.
- Use file particle count to partition records across MPI ranks.
- Derive particle type from the mask.
- Use `HH` as smoothing length and `Rho` as density.

## 17. Animation

- Snapshot filename patterns may contain one or more `{}` placeholders.
- Replace every placeholder with the requested frame number.
- `AllPath`: every MPI rank loads one frame and writes a separate file; rank 0 returns the sequence information.
- `AllMerge`: ranks load frames and merge them into one result.
- `FrameExtract`: preserve the frame-selection request semantics and merge only when the active flow requires it.
- `anim_start`, `anim_end`, and `anim_step` originate on the command line and are sent to Blender in the information response.
- BSpace may register a `frame_change_post` callback that sets the requested frame and automatically triggers extraction.

## 18. Multiresolution Mode

When compiled and requested:

- Convert active voxels to four-dimensional points `(x,y,z,value)`.
- Determine spatial/value precision and data bounds.
- Recursively subdivide the domain.
- Estimate local fractal dimension using box counts and log-log slope.
- Select a VDB resolution level based on subdivision depth and a fractal-dimension threshold.
- Create one OpenVDB float grid per resolution level.
- Name grids `density0`, `density1`, and so on.
- Give each grid a transform whose voxel size is the data range divided by the level's voxel count and whose translation is the data minimum.
- Store nonempty levels in one VDB file.
- Retain optional reconstruction statistics including min/max, error metrics, and SSIM calculations.

This feature is OpenVDB-dependent and is not required in a NanoVDB-only build.

## 19. BSpace Blender Add-on

### 19.1 Add-on metadata

```text
Name: BSpace
Version: 0.3.0
Blender: 4.2.0+
Location: View3D > Sidebar > BSpace
Category: 3D View
```

Register preferences first and panel classes second. Unregister in reverse-compatible fashion and tolerate Blender reload errors.

### 19.2 Preferences

Expose:

- server host, default `localhost`
- TCP port, default `5000`
- local temporary directory
- NanoVDB-to-VDB converter executable
- VDB merger executable
- VDB histogram executable
- VDB-to-PNG executable

### 19.3 Scene settings

Expose:

- grid dimension, default 100
- bbox size, default 1000
- extraction type
- dense kernel and normalization
- fixed particle size
- min/max value filter
- density display scale
- animation type, frame, start, end, and task counter
- automatic frame-change extraction toggle
- slice axis and number
- path replacement enable/original/new strings
- dynamic collections for available and extracted data
- merge-selection state

### 19.4 Operators

Implement operators for:

- connect/disconnect
- discover particle types and fields
- extract data
- load existing volume data
- create/select/find a dataset bbox
- derive bbox from viewport/object
- move/zoom navigation helpers
- configure world environment
- select extracted records
- merge VDB files
- compute VDB histogram
- export a VDB slice to PNG
- clear merge selection
- update extracted data
- create particle geometry/instances

### 19.5 BBox interaction

- Represent the extraction bbox as a Blender cube with a stable identifiable name.
- Bbox size changes rescale the scene representation.
- Compute min/max from object location and dimensions.
- Allow a bbox to be created from current viewport bounds.
- The `Find BBOX` action requests the selected particle type's true bounds from the converter.

### 19.6 Volume import

For OpenVDB:

- Write received bytes to the configured temporary directory.
- Load with `bpy.ops.object.volume_import`.

For NanoVDB:

- Save `.nvdb`.
- Invoke the configured conversion executable to produce `.vdb`.
- Import the resulting VDB.

For a returned path:

- Decode the path.
- Optionally replace its prefix using preference/settings mapping.
- Import the existing VDB file or sequence.

Every imported object shall receive custom metadata including the particle/block name and IDs, original and reduced min/max values, frame information, and a BSpace marker.

### 19.7 Volume shader

Build a node-based material containing:

- Principled Volume
- Material Output
- Volume Info attributes
- Map Range nodes for density and temperature/value mapping
- Color Ramp

Connect density to Principled Volume density and use the selected field to drive color. Configure blackbody/temperature behavior where appropriate. Use object metadata min/max to set mapping ranges. The viewport legend shall sample the material's color ramp and label the reduced min/max and field name.

### 19.8 Raw-particle import

- Deserialize the C++ `RawParticles` binary representation using native 64-bit sizes.
- Create a collection and mesh with one vertex per particle.
- Populate position and selected attributes.
- Support object instancing and Geometry Nodes instancing.
- Create or reuse an instance object and configure point instancing through named attributes as needed.

### 19.9 Connection behavior

- Use a blocking IPv4 socket.
- `recvall` must continue until the exact expected byte count arrives or fail explicitly.
- Pack/unpack with `struct` using native 32-bit ints, floats, 64-bit sizes, doubles, and float triples.
- Send acknowledgments after info, bbox, and data responses.
- Disconnect by sending `MessageType::Exit`, then closing and clearing the socket.

### 19.10 Auxiliary Blender tools

- VDB merge: execute the configured merger with selected files and output path.
- Histogram: execute `vdb2histo`, then load/show the PNG.
- Slice: sample the active material color ramp, encode it as comma-separated byte RGB values, invoke `vdb2png` with axis, slice, ramp, and log options.
- World setup: create a dark world background suitable for volume visualization.

## 20. Standalone Utilities

### 20.1 `vdb_merge`

```text
vdb_merge file1.vdb file2.vdb ... -o merged.vdb
```

Open every input, append all grids, disambiguate names using a filename-derived prefix when needed, and write one output file.

### 20.2 `vdb2histo`

```text
vdb2histo file.vdb -g grid-name -o histogram.png
```

- Select a float grid by name.
- Compute active-value min/max.
- Create a 256-bin histogram in parallel.
- Draw a 1024×512 RGBA histogram image.
- Save PNG through the embedded `stb_image_write`.

### 20.3 `vdb2png`

```text
vdb2png file.vdb -g grid-name -o slice.png
            --axis x|y|z --slice index
            [--ramp r,g,b,...] [--use_log]
```

- Determine slice dimensions from the active voxel bbox.
- Sample the selected index plane.
- Normalize scalar values globally within the sampled plane.
- Optionally use logarithmic mapping.
- Interpolate the supplied RGB ramp or use the default mapping.
- Write an RGB PNG.

### 20.4 `vdb2nano`

```text
vdb2nano file.vdb -g grid-name -o out.nvdb
vdb2nano file.vdb -f min max -g grid-name -o out.nvdb
vdb2nano file.vdb --print
```

Print grid names and extrema, or convert the selected float grid. Optional filtering keeps only values in the requested interval. Preserve transform, active coordinates, and values in the NanoVDB output.

### 20.5 `nano2vdb`

```text
nano2vdb file.nvdb file.vdb
```

Read a NanoVDB grid, convert through NanoVDB's OpenVDB bridge, and write it as OpenVDB.

### 20.6 `nano2float3`

```text
nano2float3 file.nvdb file.bin [--mode mesh|zip] [--grid-index N]
nano2float3 --channel-files x.nvdb y.nvdb z.nvdb file.bin
           [--mode mesh|zip] [--grid-index N]
```

- Accept one Vec3f grid or three scalar channel files.
- Require matching index bboxes for channel files.
- `mesh` writes all bbox voxels.
- `zip` writes only active voxels, or the union of active channel coordinates.
- Write each output entry as three little-endian IEEE-754 floats.
- Validate grid count, grid type, index, bbox compatibility, output stream, and count.

### 20.7 `draw_raw.py`

Provide NumPy/Matplotlib development diagnostics for:

- loading dense raw float arrays
- array min/max and difference checks
- 2D slices and color ramps
- voxel contribution visualization
- Gaussian and SPH kernel normalization checks

## 21. Platform and HPC Scripts

Retain build and run profiles for Aurora, Barbora, Karolina, Leonardo, LUMI, macOS, Polaris, and the CS P06 environment.

Each build script shall:

- establish a repository root
- load platform compiler, MPI, CUDA, and dependency modules
- specify OpenVDB/TBB/NanoVDB include and library locations
- choose optional features for that machine
- configure a machine-specific build directory and install prefix
- use `RelWithDebInfo` or Release
- build and install

Each run script shall:

- load runtime modules and library paths
- create output directories
- demonstrate representative input format arguments
- launch with the platform scheduler or MPI runner
- configure rank/thread/GPU binding where needed
- expose the TCP port for remote BSpace use

The scripts are deployment profiles and may contain site-specific project paths, allocation placeholders, and scheduler directives.

## 22. Error Handling and Diagnostics

- Reject unknown data types with a message listing supported values.
- Report file-open, grid-type, grid-name, and conversion failures.
- Validate output streams before writing.
- Convert socket failures into an exit message on all ranks.
- Use MPI barriers around major timed stages.
- Print reader load steps, bbox, conversion time, merge time, particle count, and total time.
- Keep rank-heavy diagnostics limited where practical; summary output belongs to rank 0.
- Avoid silently accepting a field with an incompatible component count.
- Guard zero bbox extents, zero smoothing lengths, empty grids, and invalid normalization denominators.
- Use RAII for files, vectors, grids, sockets, and reader allocations in maintained implementations.

## 23. Compatibility Requirements

- Preserve enum integer values and TCP field order.
- Preserve the semicolon-delimited metadata format.
- Preserve `density` as the primary volume grid name.
- Preserve `{}` frame-pattern expansion.
- Preserve native raw-particle serialization for BSpace compatibility.
- Preserve Blender object metadata used by the legend, shader, and data lists.
- Support both OpenVDB 11 build-grid APIs and later `nanovdb::tools::build` APIs.
- Keep NanoVDB-only builds possible when OpenVDB is disabled, excluding utilities and features that inherently require OpenVDB.

## 24. Validation Strategy

### 24.1 Unit tests

Add deterministic tests for:

- CLI parsing and missing-argument rejection
- filename placeholder replacement
- vector magnitude calculation
- `RawParticles` round trips and merge behavior
- dense indexing and allocation
- every SPH kernel at zero, support boundary, and outside support
- kernel normalization on sampled 3D domains
- bbox coordinate mapping
- TCP primitive encoding
- VDB stream round trips
- NanoVDB conversion and filtering
- HACC mask-to-type mapping

### 24.2 MPI tests

Run with one, two, three, and power-of-two rank counts:

- global bbox and min/max
- dense sum reduction
- sparse binary-tree merge
- raw-particle merge
- animation frame distribution
- cyclic nearest-neighbor point exchange
- transfers larger than one configured test chunk

### 24.3 Adapter tests

For a small fixture per format, verify:

- particle counts
- type/block availability
- positions and component counts
- normalized vector values
- mass, density, and smoothing length
- local partitions concatenate to the original dataset

### 24.4 Blender integration tests

Against Blender 4.2+:

- add-on registration/reload
- connection and metadata discovery
- bbox request
- OpenVDB and NanoVDB extraction
- path-based animation import
- raw-particle mesh and instancing
- shader node construction
- legend min/max display
- frame-change extraction registration and removal
- merger, histogram, and slice subprocess invocation

### 24.5 Acceptance scenarios

The implementation is complete when:

1. A single-rank local GADGET conversion produces a readable fog-volume VDB.
2. A multi-rank dense conversion matches the single-rank voxel values within floating-point tolerance.
3. Sparse and raw-particle multi-rank merges preserve all accepted data.
4. BSpace discovers available particle fields and extracts a selected sub-bbox.
5. OpenVDB, NanoVDB, path, and raw-particle response types are all handled.
6. At least one ChaNGa and one HACC-family input produce correctly named Blender data.
7. Animation filename expansion loads the expected frame on each rank.
8. Every standalone utility accepts its documented command and produces a valid output.
9. Feature-disabled builds omit optional adapters cleanly without unresolved symbols.

## 25. Implementation Order

Implement in this dependency order:

1. Core enums, state, grid containers, and serialization.
2. CMake feature structure and dependency discovery.
3. Converter base interface and common math.
4. One simple adapter and local sparse conversion.
5. Dense SPH kernels and dense-to-VDB conversion.
6. MPI partitioning and reductions.
7. Remaining adapters.
8. NanoVDB support.
9. TCP server protocol.
10. BSpace connection, discovery, bbox, and volume import.
11. Raw particles, animation, and nearest-neighbor radius search.
12. Standalone utilities.
13. Multiresolution and MERIC options.
14. Platform scripts and the complete validation suite.

All modules shall remain independently understandable, feature-gated code shall compile in every supported option combination, and protocol or enum changes shall be treated as versioned compatibility changes.
