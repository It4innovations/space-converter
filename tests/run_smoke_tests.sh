#!/bin/bash
# Smoke tests for space-converter (see docs/SpaceConverter_Code_Analysis_2026-08.md §11.4).
#
# Runs ON A COMPUTE NODE (never on a login node). Typical invocation through an
# existing Slurm job:
#
#   srun --jobid=<JOBID> --overlap -N1 -n1 -c16 bash tests/run_smoke_tests.sh
#
# Environment (all optional):
#   SC_BIN   path to the space_converter binary
#            [default: <blender root>/install/space_converter_barng_cpu/bin/space_converter]
#   SC_OUT   scratch directory for test outputs [default: <repo>/temp/tests]
#   SC_NP    max MPI ranks to test with (2 or 3 recommended) [default: 3]
#   SC_NODE  run the MPI ranks on this node of the job (srun -w) [any]
#   SC_LUMI_PARTITION  LUMI only: G (HIP build, default) or C (LUMI-C CPU builds)
#
# MPI ranks are started with `mpirun -n N` inside this (single-task) step when
# mpirun is available, else `srun --overlap -n N` against the surrounding job.
# Exit code: number of failed tests.

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$HERE")"
BLENDER_ROOT="$(cd "$REPO/../.." && pwd)"

SC_OUT="${SC_OUT:-$REPO/temp/tests}"
SC_NP="${SC_NP:-3}"

# ── Site detection: LUMI (Cray) vs Barbora-style (EasyBuild modules) ─────────
if [ -d /appl/lumi ]; then
    # LUMI: the HIP binary doubles as the CPU binary (run without --gpu)
    SC_BIN="${SC_BIN:-$BLENDER_ROOT/install/space_converter_lumi_hip/bin/space_converter}"
    # SC_LUMI_PARTITION=C for the LUMI-C CPU builds (e.g. space_converter_lumic_fil)
    if command -v ml >/dev/null 2>&1; then
        ml LUMI/25.09 partition/${SC_LUMI_PARTITION:-G} 2>/dev/null || true
        ml PrgEnv-gnu 2>/dev/null || true
        [ "${SC_LUMI_PARTITION:-G}" = G ] && { ml rocm/6.4.4 2>/dev/null || true; }
        # h5cc for the iPIC3D / FIL test data generators
        ml cray-hdf5 2>/dev/null || true
    fi
    LIB_DIR="$BLENDER_ROOT/install/lib-linux_x64"
    export LD_LIBRARY_PATH="$LIB_DIR/openvdb/lib:$LIB_DIR/tbb/lib:${ROCM_PATH:-/opt/rocm}/lib:${LD_LIBRARY_PATH:-}"
else
    SC_BIN="${SC_BIN:-$BLENDER_ROOT/install/space_converter_barng_cpu/bin/space_converter}"
    # Runtime environment (mirror run_spaceconverter_barng_cpu.sh)
    if command -v ml >/dev/null 2>&1; then
        ml HDF5/1.14.6-gompi-2025b 2>/dev/null || true
        ml GCC/14.3.0 2>/dev/null || true
    fi
    CYCLES_LIB="$BLENDER_ROOT/src/cyclesphi/lib/linux_x64"
    export LD_LIBRARY_PATH="$CYCLES_LIB/openvdb/lib:$CYCLES_LIB/tbb/lib:${LD_LIBRARY_PATH:-}"
    # VTK runtime for PLUTO-enabled builds (VTK bundled in ParaView, see
    # build_spaceconverter_barng_cpu_full.sh); harmless when absent
    PARAVIEW_LIB="/apps/all/ParaView/6.0.1-foss-2025b/lib64"
    [ -d "$PARAVIEW_LIB" ] && export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PARAVIEW_LIB"
fi
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

# Fail fast (with a clear message) if the binary cannot even start.
# (`-h` exits non-zero when argc < 3, so check for the usage banner instead.)
if ! "$SC_BIN" -h 2>&1 | grep -q "space_converter --data-type"; then
    echo "space_converter cannot start (missing runtime libs/modules?):"
    "$SC_BIN" -h 2>&1 | head -3
    exit 98
fi

# Inside a Slurm job, launch MPI ranks as overlapping job steps. When this
# driver itself runs inside a step, the inherited SLURM_* step variables would
# make the nested srun hang — strip them for the child (keeping the job id).
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SC_JOBID="$SLURM_JOB_ID"
    SLURM_UNSET_FLAGS=$(compgen -e | grep "^SLURM_\|^SRUN_" | sed "s/^/-u /")
    if [ -d /appl/lumi ]; then
        launch() {
            local n=$1; shift
            # shellcheck disable=SC2086
            env $SLURM_UNSET_FLAGS srun --jobid="$SC_JOBID" --overlap ${SC_NODE:+-w "$SC_NODE"} -n "$n" -N1 "$@"
        }
    else
    launch() {
        local n=$1; shift
        # shellcheck disable=SC2086
        env $SLURM_UNSET_FLAGS srun --jobid="$SC_JOBID" --overlap ${SC_NODE:+-w "$SC_NODE"} --mpi=pmix -n "$n" -N1 "$@"
    }
    fi
elif command -v mpirun >/dev/null 2>&1; then
    launch() { local n=$1; shift; mpirun -n "$n" "$@"; }
else
    launch() { local n=$1; shift; [ "$n" = 1 ] || { echo "no MPI launcher for n=$n"; return 97; }; "$@"; }
fi

mkdir -p "$SC_OUT"
SNAP="$SC_OUT/snap_test"
PASS=0; FAIL=0
note() { echo "== $*"; }
ok()   { echo "   PASS: $*"; PASS=$((PASS+1)); }
bad()  { echo "   FAIL: $*"; FAIL=$((FAIL+1)); }

[ -x "$SC_BIN" ] || { echo "space_converter binary not found: $SC_BIN"; exit 99; }
python3 "$HERE/gen_gadget_simple.py" "$SNAP" || exit 99

BASE_ARGS=(--data-type GADGET_SIMPLE --gadget-file "$SNAP" --output-path "$SC_OUT" --grid-dim 64)

# ── T1: argument validation ──────────────────────────────────────────────────
note "T1: argument validation"
if launch 1 "$SC_BIN" --data-type GADGET_SIMPLE --no-such-flag >/dev/null 2>&1; then
    bad "unknown flag was accepted"
else
    ok "unknown flag rejected"
fi
if launch 1 "$SC_BIN" "${BASE_ARGS[@]}" --bbox 1 2 3 >/dev/null 2>&1; then
    bad "truncated --bbox was accepted"
else
    ok "truncated --bbox rejected"
fi
if launch 1 "$SC_BIN" "${BASE_ARGS[@]}" --grid-dim 0 --export-data 0 0 >/dev/null 2>&1; then
    bad "--grid-dim 0 was accepted"
else
    ok "--grid-dim 0 rejected"
fi

# The count reported after conversion (the config dump earlier in the log also
# contains a particles_count line — match the post-conversion one specifically)
count_of() { grep -oE "find minmax mpi - particles_count: [0-9]+" "$1" | head -1 | grep -oE "[0-9]+$"; }

run_case() { # <log> <nranks> <extra args...>
    local log=$1 n=$2; shift 2
    rm -f "$SC_OUT"/*.vdb "$SC_OUT"/*_float.raw 2>/dev/null
    launch "$n" "$SC_BIN" "${BASE_ARGS[@]}" "$@" >"$log" 2>&1
}

# ── T2: single-rank sparse extraction, type 0 (gas) ──────────────────────────
note "T2: sparse extraction, type 0"
if run_case "$SC_OUT/t2.log" 1 --export-data 0 0 && [ "$(count_of "$SC_OUT/t2.log")" = "1000" ] && ls "$SC_OUT"/*.vdb >/dev/null 2>&1; then
    ok "1000 gas particles voxelized, .vdb written"
else
    bad "type-0 extraction (count=$(count_of "$SC_OUT/t2.log" || echo '?')) — see t2.log"
fi

# ── T3: particle type > 0 (broken before the id-convention fix) ──────────────
note "T3: sparse extraction, type 1"
if run_case "$SC_OUT/t3.log" 1 --export-data 1 0 && [ "$(count_of "$SC_OUT/t3.log")" = "500" ] && ls "$SC_OUT"/*.vdb >/dev/null 2>&1; then
    ok "500 halo particles voxelized"
else
    bad "type-1 extraction (count=$(count_of "$SC_OUT/t3.log" || echo '?')) — see t3.log"
fi

# ── T4: rank invariance (dense grid, RAW dump, tolerant compare) ─────────────
# Bit-identity is not expected: float accumulation order differs with rank
# count and OpenMP scheduling, so the grids are compared with a tolerance.
note "T4: rank invariance of the dense grid (1..$SC_NP ranks)"
T4_OK=1
REF_RAW="$SC_OUT/t4_ref.raw"
REF_CNT=""
for n in $(seq 1 "$SC_NP"); do
    if ! run_case "$SC_OUT/t4_$n.log" "$n" --export-data 0 0 --dense-type 6 --dense-norm 0 --radius-const 2 --dense-file 1; then
        bad "T4: run with $n ranks failed — see t4_$n.log"; T4_OK=0; break
    fi
    raw=$(ls "$SC_OUT"/*_float.raw 2>/dev/null | head -1)
    if [ -z "$raw" ]; then bad "T4: no RAW dump with $n ranks"; T4_OK=0; break; fi
    cnt=$(count_of "$SC_OUT/t4_$n.log")
    if [ -z "$REF_CNT" ]; then
        cp "$raw" "$REF_RAW"; REF_CNT="$cnt"
    else
        if [ "$cnt" != "$REF_CNT" ]; then
            bad "T4: particle count with $n ranks: $cnt (expected $REF_CNT)"; T4_OK=0; break
        fi
        if ! python3 "$HERE/compare_raw.py" "$REF_RAW" "$raw" > "$SC_OUT/t4_cmp_$n.log" 2>&1; then
            bad "T4: grid with $n ranks differs beyond tolerance — see t4_cmp_$n.log"; T4_OK=0; break
        fi
    fi
done
[ "$T4_OK" = 1 ] && ok "dense grid equal within tolerance across 1..$SC_NP ranks"

# ── T5: k-NN radius sanity — requires a cudakdtree/nanoflann build ───────────
note "T5: k-NN radius (skipped unless the build has cudakdtree/nanoflann)"
HELP_OUT=$("$SC_BIN" -h 2>&1)
if echo "$HELP_OUT" | grep -q -- "--nanoflann"; then
    KNN_FLAG="--nanoflann"
elif echo "$HELP_OUT" | grep -q -- "--cudakdtree-cpu"; then
    KNN_FLAG="--cudakdtree-cpu"
else
    KNN_FLAG=""
fi
if [ -n "$KNN_FLAG" ]; then
    if run_case "$SC_OUT/t5a.log" 1 --export-data 0 0 $KNN_FLAG --calc-radius-neigh 6 \
        && ! grep -qiE "nan|inf" <(grep "minI" "$SC_OUT/t5a.log"); then
        ok "1-rank k-NN ($KNN_FLAG) produced finite values"
    else
        bad "1-rank k-NN — see t5a.log"
    fi
    if run_case "$SC_OUT/t5b.log" 2 --export-data 0 0 $KNN_FLAG --calc-radius-neigh 6 \
        && ! grep -qiE "nan|inf" <(grep "minI" "$SC_OUT/t5b.log"); then
        ok "2-rank k-NN cycling ($KNN_FLAG) produced finite values"
    else
        bad "2-rank k-NN cycling — see t5b.log (radius=0/NaN was the pre-fix symptom)"
    fi
else
    echo "   SKIP: build has no k-NN support"
fi

# ── T6: streaming (--skip-cache-manager) equivalence ─────────────────────────
note "T6: cached vs streaming equivalence (dense RAW md5)"
if run_case "$SC_OUT/t6a.log" 1 --export-data 1 0 --dense-type 6 --radius-const 2 --dense-file 1; then
    raw_a=$(ls "$SC_OUT"/*_float.raw | head -1); cp "$raw_a" "$SC_OUT/t6_ref.raw"
    if run_case "$SC_OUT/t6b.log" 1 --export-data 1 0 --dense-type 6 --radius-const 2 --dense-file 1 --skip-cache-manager; then
        raw_b=$(ls "$SC_OUT"/*_float.raw | head -1)
        if python3 "$HERE/compare_raw.py" "$SC_OUT/t6_ref.raw" "$raw_b" > "$SC_OUT/t6_cmp.log" 2>&1; then
            ok "cached and streaming grids equal within tolerance (type 1)"
        else
            bad "cached vs streaming grids differ — see t6_cmp.log"
        fi
    else
        bad "streaming run failed — see t6b.log"
    fi
else
    bad "cached run failed — see t6a.log"
fi

# ── T8: TCP protocol round-trip ──────────────────────────────────────────────
note "T8: protocol round-trip (server + python client)"
PORT=$(( 20000 + RANDOM % 20000 ))
"$SC_BIN" --data-type GADGET_SIMPLE --gadget-file "$SNAP" --output-path "$SC_OUT" --port "$PORT" \
    >"$SC_OUT/t8_server.log" 2>&1 &
SRV_PID=$!
CLIENT_RC=1
for i in $(seq 1 30); do
    sleep 1
    if python3 "$HERE/protocol_client.py" localhost "$PORT" 0 0 >"$SC_OUT/t8_client.log" 2>&1; then
        CLIENT_RC=0; break
    fi
    kill -0 "$SRV_PID" 2>/dev/null || break   # server died
done
kill "$SRV_PID" 2>/dev/null; wait "$SRV_PID" 2>/dev/null
if [ "$CLIENT_RC" = 0 ]; then
    ok "info + extraction round-trip over TCP"
else
    bad "protocol round-trip — see t8_server.log / t8_client.log"
fi

# ── T9: non-cubic box symmetrization ─────────────────────────────────────────
note "T9: non-cubic dataset bbox is cube-symmetrized"
BOX_SIZE=$(grep -oE "box_size: [0-9]+" "$SC_OUT/t2.log" | head -1 | grep -oE "[0-9]+$")
if [ -n "$BOX_SIZE" ] && [ "$BOX_SIZE" -ge 370 ] && [ "$BOX_SIZE" -le 420 ]; then
    ok "bbox size follows the longest data axis (~400): $BOX_SIZE"
else
    bad "unexpected bbox size: '${BOX_SIZE:-none}' (expected ~380-410 for a 100x200x400 box)"
fi

# ── T11: weight-normalised dense grid on several ranks ───────────────────────
# --dense-norm 1/2 also reduce the per-voxel weight grid over MPI. That second
# reduction was once skipped on the non-root ranks (it was conditioned on the
# rank-0-only sum grid), so every multi-rank run hung. Each rank runs under
# `timeout`, so a hang fails the test instead of blocking the suite.
note "T11: --dense-norm 2 across 1..2 ranks (MPI reduction of the weight grid)"
T11_OK=1
T11_REF="$SC_OUT/t11_ref.raw"
for n in 1 2; do
    rm -f "$SC_OUT"/*.vdb "$SC_OUT"/*_float.raw 2>/dev/null
    if ! launch "$n" timeout 300 "$SC_BIN" "${BASE_ARGS[@]}" --export-data 0 0 --dense-type 6 --dense-norm 2 \
        --radius-const 2 --dense-file 1 >"$SC_OUT/t11_$n.log" 2>&1; then
        bad "T11: run with $n ranks failed or hung (killed after 300 s) — see t11_$n.log"; T11_OK=0; break
    fi
    raw=$(ls "$SC_OUT"/*_float.raw 2>/dev/null | head -1)
    if [ -z "$raw" ]; then bad "T11: no RAW dump with $n ranks"; T11_OK=0; break; fi
    if [ "$n" = 1 ]; then
        cp "$raw" "$T11_REF"
    elif ! python3 "$HERE/compare_raw.py" "$T11_REF" "$raw" > "$SC_OUT/t11_cmp.log" 2>&1; then
        bad "T11: normalised grid with $n ranks differs beyond tolerance — see t11_cmp.log"; T11_OK=0
    fi
done
[ "$T11_OK" = 1 ] && ok "normalised dense grid completes and matches across 1..2 ranks"

# ── T10: per-reader smoke tests (each runs only when the build supports it) ──
# reader_case <label> <expected_count> <ptype> <block> <converter args...>
reader_case() {
    local label=$1 expect=$2 ptype=$3 block=$4; shift 4
    local log="$SC_OUT/t10_${label}.log"
    rm -f "$SC_OUT"/*.vdb 2>/dev/null
    launch 1 "$SC_BIN" "$@" --output-path "$SC_OUT" --grid-dim 32 --export-data "$ptype" "$block" >"$log" 2>&1
    if grep -q "Unknown data type" "$log"; then
        echo "   SKIP: $label (reader not in this build — use build_spaceconverter_barng_cpu_full.sh)"
        return
    fi
    local cnt; cnt=$(count_of "$log")
    if [ "$cnt" = "$expect" ] && ls "$SC_OUT"/*.vdb >/dev/null 2>&1; then
        ok "$label: $cnt particles voxelized"
    else
        bad "$label: count=${cnt:-?} (expected $expect) — see t10_${label}.log"
    fi
}

note "T10: per-reader smoke tests"
# GADGET (OpenGadget3 CodeBase reader) reads the same format-2 snapshot as
# GADGET_SIMPLE; no parameter file is required (defaults suffice).
reader_case gadget_gas 1000 0 0 --data-type GADGET --gadget-file "$SNAP"
reader_case gadget_dark_mass 500 1 3 --data-type GADGET --gadget-file "$SNAP"
# Block 3 = Masses; every synthetic halo particle has mass 2.5, so min == max
# == 2.5 — a direct regression check for the "non-gas mass is 0" fix (M10)
if grep -q "minI: 2.500000e+00, maxI: 2.500000e+00" "$SC_OUT/t10_gadget_dark_mass.log"; then
    ok "gadget_dark_mass: halo mass 2.5 read correctly (M10 regression)"
else
    bad "gadget_dark_mass: expected min=max=2.5 — $(grep -oE 'minI: [^,]+, maxI: [^,]+' "$SC_OUT/t10_gadget_dark_mass.log" | head -1)"
fi

python3 "$HERE/gen_tipsy.py" "$SC_OUT/tipsy_test" >/dev/null
reader_case tipsy_gas  700 0 0 --data-type CHANGA_TIPSY --tipsy-file "$SC_OUT/tipsy_test"
reader_case tipsy_dark 350 1 0 --data-type CHANGA_TIPSY --tipsy-file "$SC_OUT/tipsy_test"

python3 "$HERE/gen_haccbin.py" "$SC_OUT/haccbin_test" >/dev/null
reader_case haccbin 800 0 0 --data-type HACC_BIN --haccbin-file "$SC_OUT/haccbin_test"

python3 "$HERE/gen_pluto_vtk.py" "$SC_OUT/pluto_test.vtk" >/dev/null
reader_case pluto 1536 0 0 --data-type PLUTO_VTK --vtk-file "$SC_OUT/pluto_test.vtk"

# RAMSES: 2 CPU files, 15 AMR leaf cells over two levels, 40 DM + 20 stars
# (10 of them born before t = 0, selected with --ramses-filter), 5 sinks
# (sink_00001.csv, type 6; 3 of them with birth_time >= 0.015)
python3 "$HERE/gen_ramses.py" "$SC_OUT/ramses_test" >/dev/null
reader_case ramses_gas  15 0 2 --data-type RAMSES --ramses-output "$SC_OUT/ramses_test/output_00001"
reader_case ramses_dm   40 1 1 --data-type RAMSES --ramses-output "$SC_OUT/ramses_test/output_00001" --no-gas
reader_case ramses_old  10 2 1 --data-type RAMSES --ramses-output "$SC_OUT/ramses_test/output_00001" --ramses-filter birth_time -1e30 0
reader_case ramses_sink  5 6 1 --data-type RAMSES --ramses-output "$SC_OUT/ramses_test/output_00001" --no-gas
reader_case ramses_young 3 6 1 --data-type RAMSES --ramses-output "$SC_OUT/ramses_test/output_00001" --no-gas --ramses-filter birth_time 0.015 1e30

python3 "$HERE/gen_nchilada.py" "$SC_OUT/nchilada_test" >/dev/null
reader_case nchilada_gas  600 0 0 --data-type CHANGA_NCHILADA --nc-dir "$SC_OUT/nchilada_test"
reader_case nchilada_dark 300 1 1 --data-type CHANGA_NCHILADA --nc-dir "$SC_OUT/nchilada_test"

# iPIC3D needs an HDF5-generated dataset: compile the C generator with h5cc
# (from the HDF5 module). Counts include the synthetic grid points the reader
# adds per species: the 10^3 tiles carry a 1-cell ghost layer that is not
# deposited, so 8^3 = 512 points each (400+512 and 200+512).
if command -v h5cc >/dev/null 2>&1; then
    if [ ! -x "$SC_OUT/gen_ipic3d" ] || [ "$HERE/gen_ipic3d.c" -nt "$SC_OUT/gen_ipic3d" ]; then
        h5cc -o "$SC_OUT/gen_ipic3d" "$HERE/gen_ipic3d.c" || echo "   SKIP: ipic3d (h5cc compile failed)"
    fi
    if [ -x "$SC_OUT/gen_ipic3d" ]; then
        "$SC_OUT/gen_ipic3d" "$SC_OUT/ipic3d_test" >/dev/null
        reader_case ipic3d_sp0 912 0 0 --data-type IPIC3D_HDF5 --hdf5-file "$SC_OUT/ipic3d_test/restart{}.hdf" --num-files 1 --settings-file "$SC_OUT/ipic3d_test/settings.hdf"
        reader_case ipic3d_sp1 712 1 0 --data-type IPIC3D_HDF5 --hdf5-file "$SC_OUT/ipic3d_test/restart{}.hdf" --num-files 1 --settings-file "$SC_OUT/ipic3d_test/settings.hdf"
    fi
else
    echo "   SKIP: ipic3d (h5cc not available — load the HDF5 module)"
fi

# FIL (Carpet HDF5): C generator compiled with h5cc, like iPIC3D. Two
# refinement levels (level 0 split into two components with ghost zones),
# per-process, per-group and per-variable files plus a 2D slice file that must
# be ignored. Finest level wins: 702 + 729 = 1431 points; levels 0..0: 729;
# ghost zones kept: 1728. Blocks: 0 Pos, 1 Mass, 2 Rho, 3 Level, 4 alp, 5 rho,
# 6..8 vel[0..2], 9 vel (vector, |vel| = 0.5).
if command -v h5cc >/dev/null 2>&1; then
    if [ ! -x "$SC_OUT/gen_fil" ] || [ "$HERE/gen_fil.c" -nt "$SC_OUT/gen_fil" ]; then
        h5cc -o "$SC_OUT/gen_fil" "$HERE/gen_fil.c" || echo "   SKIP: fil (h5cc compile failed)"
    fi
    if [ -x "$SC_OUT/gen_fil" ]; then
        FIL_DIR="$SC_OUT/fil_test"
        rm -rf "$FIL_DIR"
        "$SC_OUT/gen_fil" "$FIL_DIR" >/dev/null
        reader_case fil_rho     1431 0 2 --data-type FIL --fil-dir "$FIL_DIR"
        reader_case fil_coarse   729 0 2 --data-type FIL --fil-dir "$FIL_DIR" --fil-levels 0 0
        reader_case fil_ghosts  1728 0 2 --data-type FIL --fil-dir "$FIL_DIR" --fil-ghosts
        # explicit files, the 2D slice file included (its rank-2 datasets must be skipped)
        reader_case fil_files   1431 0 5 --data-type FIL --fil-vars rho \
            --fil-file "$FIL_DIR/hydrobase-rho.xy.h5" --fil-file "$FIL_DIR/hydrobase-rho.file_0.h5" \
            --fil-file "$FIL_DIR/hydrobase-rho.file_1.h5"
        reader_case fil_vel     1431 0 9 --data-type FIL --fil-dir "$FIL_DIR"
        if grep -q "minI: 5.0000[0-9]*e-01, maxI: 5.0000[0-9]*e-01" "$SC_OUT/t10_fil_vel.log"; then
            ok "fil_vel: vector group vel[0..2] has magnitude 0.5"
        else
            bad "fil_vel: expected |vel| = 0.5 — $(grep -oE 'minI: [^,]+, maxI: [^,]+' "$SC_OUT/t10_fil_vel.log" | head -1)"
        fi
        if grep -q "iteration 4, t 1," "$SC_OUT/t10_fil_rho.log"; then
            ok "fil: the latest iteration (4, t = 1) is the default"
        else
            bad "fil: default iteration is not 4 — see t10_fil_rho.log"
        fi
        reader_case fil_it0     1431 0 3 --data-type FIL --fil-dir "$FIL_DIR" --fil-iteration 0
        if grep -q "minI: 0.000000e+00, maxI: 1.000000e+00" "$SC_OUT/t10_fil_it0.log" \
            && grep -q "iteration 0, t 0," "$SC_OUT/t10_fil_it0.log"; then
            ok "fil_it0: iteration 0 selected, levels 0 and 1 present"
        else
            bad "fil_it0: iteration/level range wrong — see t10_fil_it0.log"
        fi
        # two ranks share the patches: the global count must not change
        rm -f "$SC_OUT"/*.vdb 2>/dev/null
        launch 2 "$SC_BIN" --data-type FIL --fil-dir "$FIL_DIR" --output-path "$SC_OUT" --grid-dim 32 \
            --export-data 0 2 >"$SC_OUT/t10_fil_np2.log" 2>&1
        if [ "$(count_of "$SC_OUT/t10_fil_np2.log")" = "1431" ]; then
            ok "fil_np2: 1431 points over 2 ranks"
        else
            bad "fil_np2: count=$(count_of "$SC_OUT/t10_fil_np2.log" || echo '?') (expected 1431) — see t10_fil_np2.log"
        fi
    fi
else
    echo "   SKIP: fil (h5cc not available — load the HDF5 module)"
fi

# HACC_GENERICIO: the dataset is written by the gen_genericio tool (built and
# installed next to space_converter in WITH_HACC builds — the GenericIO format
# is self-describing with CRCs, so it is synthesized via the bundled library)
GEN_GIO="$(dirname "$SC_BIN")/gen_genericio"
if [ -x "$GEN_GIO" ]; then
    launch 1 "$GEN_GIO" "$SC_OUT/genericio_test" >/dev/null 2>&1
    reader_case genericio 900 0 0 --data-type HACC_GENERICIO --genericio-file "$SC_OUT/genericio_test" \
        --pos-names x y z --vel-names vx vy vz --mass-name mass --rho-name rho --hsml-name hh
else
    echo "   SKIP: genericio (gen_genericio tool not in this build — use build_spaceconverter_barng_cpu_full.sh)"
fi

# Not covered here: GADGET (the OpenGadget3 CodeBase reader needs a full
# parameter file + real snapshot). Use a real dataset for it.

echo
echo "==== smoke tests: $PASS passed, $FAIL failed ===="
exit "$FAIL"
