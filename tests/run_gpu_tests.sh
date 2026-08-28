#!/bin/bash
# GPU tests G1..G12 from docs/GPU_Test_Plan.md.
#
# Runs ON A GPU COMPUTE NODE (never on a login node), e.g.:
#   srun --jobid=<GPU_JOBID> --overlap -N1 -n1 -c8 bash tests/run_gpu_tests.sh
#
# Environment (optional):
#   SC_BIN   space_converter binary with WITH_GPU_CUDA=ON
#            [default: <blender root>/install/space_converter_bar_gpu/bin/space_converter]
#   SC_OUT   scratch directory [default: <repo>/temp/gpu_tests]
#
# Every GPU result is compared against the SAME binary run without --gpu
# (CPU reference). Float accumulation order differs between CPU and GPU, so
# dense grids are compared with tests/compare_raw.py tolerances.
# Exit code = number of failed tests.

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$HERE")"
BLENDER_ROOT="$(cd "$REPO/../.." && pwd)"

SC_BIN="${SC_BIN:-$BLENDER_ROOT/install/space_converter_bar_gpu/bin/space_converter}"
SC_OUT="${SC_OUT:-$REPO/temp/gpu_tests}"

# Runtime environment: toolchain of build_spaceconverter_bar_gpu.sh, plus a
# newer libstdc++ for the prebuilt OpenVDB (built with a newer GCC)
if command -v ml >/dev/null 2>&1; then
    ml GCC/13.3.0 2>/dev/null || true
    ml OpenMPI/5.0.3-GCC-13.3.0 2>/dev/null || true
    ml CUDA/12.6.0 2>/dev/null || true
fi
CYCLES_LIB="$BLENDER_ROOT/src/cyclesphi/lib/linux_x64"
export LD_LIBRARY_PATH="$CYCLES_LIB/openvdb/lib:$CYCLES_LIB/tbb/lib:${LD_LIBRARY_PATH:-}"
GCC14_LIB="/apps/all/GCCcore/14.3.0/lib64"
[ -d "$GCC14_LIB" ] && export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$GCC14_LIB"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

# Nested-srun launcher (same pattern as run_smoke_tests.sh)
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SC_JOBID="$SLURM_JOB_ID"
    SLURM_UNSET_FLAGS=$(compgen -e | grep "^SLURM_\|^SRUN_" | sed "s/^/-u /")
    launch() {
        local n=$1; shift
        # shellcheck disable=SC2086
        env $SLURM_UNSET_FLAGS srun --jobid="$SC_JOBID" --overlap --mpi=pmix -n "$n" -N1 "$@"
    }
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

if ! "$SC_BIN" -h 2>&1 | grep -q "space_converter --data-type"; then
    echo "space_converter cannot start (missing runtime libs/modules?):"
    "$SC_BIN" -h 2>&1 | head -3
    exit 98
fi
if ! "$SC_BIN" -h 2>&1 | grep -q -- "--gpu"; then
    echo "SC_BIN has no --gpu support (need a WITH_GPU_CUDA build): $SC_BIN"
    exit 98
fi
if ! nvidia-smi -L >/dev/null 2>&1; then
    echo "no visible GPU on this node"
    exit 98
fi

python3 "$HERE/gen_gadget_simple.py" "$SNAP" || exit 99

BASE=(--data-type GADGET_SIMPLE --gadget-file "$SNAP" --output-path "$SC_OUT" --grid-dim 64)
DENSE=(--dense-type 6 --radius-const 2 --dense-file 1)

count_of()  { grep -oE "find minmax mpi - particles_count: [0-9]+" "$1" | head -1 | grep -oE "[0-9]+$"; }
voxels_of() { grep -oE "Active voxels in input grid: [0-9]+" "$1" | head -1 | grep -oE "[0-9]+$"; }
minmax_of() { grep -oE "minI: [^,]+, maxI: [^,]+, reduced" "$1" | head -1; }

# NOTE on the shared filesystem: files written on the compute node can take
# many seconds to become visible here (attribute caching), and every dense run
# writes the SAME RAW filename. run_case therefore drains stale RAW files
# before launching, and keep_raw waits for the new one to appear.
run_case() { # <log> <nranks> <args...>
    local log=$1 n=$2; shift 2
    rm -f "$SC_OUT"/*.vdb "$SC_OUT"/*.nvdb "$SC_OUT"/*.cub 2>/dev/null
    for i in $(seq 1 15); do
        rm -f "$SC_OUT"/*_float.raw 2>/dev/null
        ls "$SC_OUT"/*_float.raw >/dev/null 2>&1 || break
        sleep 1
    done
    launch "$n" "$SC_BIN" "${BASE[@]}" "$@" >"$log" 2>&1
}
keep_raw() { # <dest>  (wait for and move the produced RAW dump aside)
    local dest=$1 raw
    for i in $(seq 1 30); do
        raw=$(ls "$SC_OUT"/*_float.raw 2>/dev/null | head -1)
        if [ -n "$raw" ]; then mv "$raw" "$dest"; return 0; fi
        sleep 1
    done
    echo "   (keep_raw: no RAW dump appeared for $dest)"
    return 1
}

# ── G1: GPU sparse extraction equals CPU ─────────────────────────────────────
note "G1: sparse CPU vs GPU (type 0)"
run_case "$SC_OUT/g1_cpu.log" 1 --export-data 0 0
run_case "$SC_OUT/g1_gpu.log" 1 --export-data 0 0 --gpu
if [ "$(count_of "$SC_OUT/g1_cpu.log")" = "$(count_of "$SC_OUT/g1_gpu.log")" ] \
   && [ "$(voxels_of "$SC_OUT/g1_cpu.log")" = "$(voxels_of "$SC_OUT/g1_gpu.log")" ] \
   && [ -n "$(count_of "$SC_OUT/g1_gpu.log")" ]; then
    ok "count $(count_of "$SC_OUT/g1_gpu.log"), voxels $(voxels_of "$SC_OUT/g1_gpu.log") match"
else
    bad "G1: cpu=$(count_of "$SC_OUT/g1_cpu.log")/$(voxels_of "$SC_OUT/g1_cpu.log") gpu=$(count_of "$SC_OUT/g1_gpu.log")/$(voxels_of "$SC_OUT/g1_gpu.log")"
fi
if [ "$(minmax_of "$SC_OUT/g1_cpu.log")" = "$(minmax_of "$SC_OUT/g1_gpu.log")" ] && [ -n "$(minmax_of "$SC_OUT/g1_gpu.log")" ]; then
    ok "per-particle min/max semantics match CPU: $(minmax_of "$SC_OUT/g1_gpu.log" | cut -d, -f1-2)"
else
    bad "G1 min/max: CPU [$(minmax_of "$SC_OUT/g1_cpu.log")] vs GPU [$(minmax_of "$SC_OUT/g1_gpu.log")]"
fi

# ── G2: GPU dense splat equals CPU ───────────────────────────────────────────
note "G2: dense CPU vs GPU (WendlandC6, radius-const)"
run_case "$SC_OUT/g2_cpu.log" 1 --export-data 0 0 "${DENSE[@]}"; keep_raw "$SC_OUT/g2_cpu.raw"
run_case "$SC_OUT/g2_gpu.log" 1 --export-data 0 0 "${DENSE[@]}" --gpu; keep_raw "$SC_OUT/g2_gpu.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g2_cpu.raw" "$SC_OUT/g2_gpu.raw" 1e-3 1e-5 > "$SC_OUT/g2_cmp.log" 2>&1; then
    ok "dense grids equal within tolerance"
else
    bad "G2: dense grids differ — see g2_cmp.log"
fi
if [ "$(minmax_of "$SC_OUT/g2_cpu.log")" = "$(minmax_of "$SC_OUT/g2_gpu.log")" ] && [ -n "$(minmax_of "$SC_OUT/g2_gpu.log")" ]; then
    ok "dense min/max semantics match CPU"
else
    bad "G2 min/max: CPU [$(minmax_of "$SC_OUT/g2_cpu.log")] vs GPU [$(minmax_of "$SC_OUT/g2_gpu.log")]"
fi

# ── G2b: dense GPU with a zoom sub-box ───────────────────────────────────────
note "G2b: dense CPU vs GPU with --bbox sub-region"
SUBBOX=(--bbox 100 100 100 600 600 600)
run_case "$SC_OUT/g2b_cpu.log" 1 --export-data 0 0 "${DENSE[@]}" "${SUBBOX[@]}"; keep_raw "$SC_OUT/g2b_cpu.raw"
run_case "$SC_OUT/g2b_gpu.log" 1 --export-data 0 0 "${DENSE[@]}" "${SUBBOX[@]}" --gpu; keep_raw "$SC_OUT/g2b_gpu.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g2b_cpu.raw" "$SC_OUT/g2b_gpu.raw" 1e-3 1e-5 > "$SC_OUT/g2b_cmp.log" 2>&1; then
    ok "sub-box dense grids equal within tolerance"
else
    bad "G2b: sub-box grids differ — see g2b_cmp.log"
fi

# ── G3: particle type > 0 on GPU (compact-id convention on device) ───────────
note "G3: type 1 on GPU (sparse + dense)"
run_case "$SC_OUT/g3s_gpu.log" 1 --export-data 1 0 --gpu
if [ "$(count_of "$SC_OUT/g3s_gpu.log")" = "500" ]; then
    ok "sparse type 1: 500 particles"
else
    bad "G3 sparse type 1: count=$(count_of "$SC_OUT/g3s_gpu.log") (expected 500)"
fi
run_case "$SC_OUT/g3d_cpu.log" 1 --export-data 1 0 "${DENSE[@]}"; keep_raw "$SC_OUT/g3d_cpu.raw"
run_case "$SC_OUT/g3d_gpu.log" 1 --export-data 1 0 "${DENSE[@]}" --gpu; keep_raw "$SC_OUT/g3d_gpu.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g3d_cpu.raw" "$SC_OUT/g3d_gpu.raw" 1e-3 1e-5 > "$SC_OUT/g3d_cmp.log" 2>&1; then
    ok "dense type 1 equal within tolerance"
else
    bad "G3 dense type 1 differs — see g3d_cmp.log"
fi

# ── G4: GPU sorts must not change the result ─────────────────────────────────
note "G4: GPU sorts (radius / Morton / both)"
for variant in "--sort-by-radius" "--sort-by-non-overlap" "--sort-by-radius --sort-by-non-overlap"; do
    # shellcheck disable=SC2086
    run_case "$SC_OUT/g4.log" 1 --export-data 0 0 "${DENSE[@]}" --gpu $variant; keep_raw "$SC_OUT/g4.raw"
    if python3 "$HERE/compare_raw.py" "$SC_OUT/g2_gpu.raw" "$SC_OUT/g4.raw" 1e-3 1e-5 > "$SC_OUT/g4_cmp.log" 2>&1; then
        ok "sorted run ($variant) matches unsorted"
    else
        bad "G4 ($variant) differs — see g4_cmp.log"
    fi
done

# ── G5: cudaKDTree k-NN, single rank, backend comparison ─────────────────────
note "G5: k-NN radius backends (GPU tree vs host tree vs nanoflann)"
KNN=(--calc-radius-neigh 6 --dense-type 6 --dense-file 1)
run_case "$SC_OUT/g5_gpu.log" 1 --export-data 0 0 "${KNN[@]}" --gpu --cudakdtree; keep_raw "$SC_OUT/g5_gpu.raw"
if [ -s "$SC_OUT/g5_gpu.raw" ] && ! grep -qiE "nan|inf" <(grep "minI" "$SC_OUT/g5_gpu.log"); then
    ok "GPU k-NN produced finite values"
else
    bad "G5 GPU k-NN — see g5_gpu.log"
fi
run_case "$SC_OUT/g5_cpu.log" 1 --export-data 0 0 "${KNN[@]}" --cudakdtree-cpu; keep_raw "$SC_OUT/g5_cpu.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g5_cpu.raw" "$SC_OUT/g5_gpu.raw" 5e-2 1e-4 5e-3 > "$SC_OUT/g5_cmp.log" 2>&1; then
    ok "GPU tree matches host tree within 5%"
else
    bad "G5 GPU vs host tree differ — see g5_cmp.log"
fi
run_case "$SC_OUT/g5_nf.log" 1 --export-data 0 0 "${KNN[@]}" --nanoflann; keep_raw "$SC_OUT/g5_nf.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g5_nf.raw" "$SC_OUT/g5_gpu.raw" 5e-2 1e-4 5e-3 > "$SC_OUT/g5_nf_cmp.log" 2>&1; then
    ok "GPU tree matches nanoflann within 5%"
else
    bad "G5 GPU vs nanoflann differ — see g5_nf_cmp.log"
fi

# ── G6: k-NN cycling across MPI ranks (THE big regression) ───────────────────
note "G6: GPU k-NN cycling at 2 and 4 ranks (+ batched path)"
for n in 2 4; do
    run_case "$SC_OUT/g6_$n.log" "$n" --export-data 0 0 "${KNN[@]}" --gpu --cudakdtree; keep_raw "$SC_OUT/g6_$n.raw"
    if grep -qiE "nan|inf" <(grep "minI" "$SC_OUT/g6_$n.log"); then
        bad "G6: $n-rank cycling produced NaN/inf (pre-fix symptom)"
    elif python3 "$HERE/compare_raw.py" "$SC_OUT/g5_gpu.raw" "$SC_OUT/g6_$n.raw" 5e-2 1e-4 5e-3 > "$SC_OUT/g6_cmp_$n.log" 2>&1; then
        ok "$n-rank cycling matches 1-rank within 5%"
    else
        bad "G6: $n-rank cycling differs from 1 rank — see g6_cmp_$n.log"
    fi
done
CUDAKDTREE_NUM_SPLITS=4 run_case "$SC_OUT/g6_split.log" 2 --export-data 0 0 "${KNN[@]}" --gpu --cudakdtree; keep_raw "$SC_OUT/g6_split.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g5_gpu.raw" "$SC_OUT/g6_split.raw" 5e-2 1e-4 5e-3 > "$SC_OUT/g6_split_cmp.log" 2>&1; then
    ok "2-rank cycling with NUM_SPLITS=4 matches 1-rank"
else
    bad "G6 split: differs — see g6_split_cmp.log"
fi

# ── G7 is covered by G5 (GPU k-NN radii reach the device after the upload-order fix)

# ── G8: GPU multi-rank reduction ─────────────────────────────────────────────
note "G8: GPU rank invariance (dense + sparse, 1 vs 2 ranks)"
run_case "$SC_OUT/g8d_2.log" 2 --export-data 0 0 "${DENSE[@]}" --gpu; keep_raw "$SC_OUT/g8d_2.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g2_gpu.raw" "$SC_OUT/g8d_2.raw" 1e-3 1e-5 > "$SC_OUT/g8_cmp.log" 2>&1; then
    ok "dense 2-rank GPU matches 1-rank"
else
    bad "G8 dense: 2-rank differs — see g8_cmp.log"
fi
run_case "$SC_OUT/g8s_2.log" 2 --export-data 0 0 --gpu
if [ "$(count_of "$SC_OUT/g8s_2.log")" = "1000" ] && [ "$(voxels_of "$SC_OUT/g8s_2.log")" = "$(voxels_of "$SC_OUT/g1_gpu.log")" ]; then
    ok "sparse 2-rank GPU matches 1-rank (count + voxels)"
else
    bad "G8 sparse: count=$(count_of "$SC_OUT/g8s_2.log") voxels=$(voxels_of "$SC_OUT/g8s_2.log") vs 1-rank voxels $(voxels_of "$SC_OUT/g1_gpu.log")"
fi

# ── G9: voxel-centric dense loop (KD-tree per voxel) ─────────────────────────
note "G9: --dense-loop-over-voxels on GPU"
run_case "$SC_OUT/g9.log" 1 --export-data 0 0 --dense-type 6 --dense-file 1 \
    --gpu --cudakdtree --calc-radius-neigh 8 --dense-loop-over-voxels; keep_raw "$SC_OUT/g9.raw"
if [ -s "$SC_OUT/g9.raw" ] && ! grep -qiE "nan|inf" <(grep "minI" "$SC_OUT/g9.log"); then
    ok "voxel-centric GPU loop produced a finite grid"
else
    bad "G9 — see g9.log"
fi

# ── G10: GPU NanoVDB output ──────────────────────────────────────────────────
note "G10: GPU NanoVDB output"
run_case "$SC_OUT/g10.log" 1 --export-data 0 0 --gpu --nanovdb
if ls "$SC_OUT"/*.nvdb >/dev/null 2>&1; then
    ok ".nvdb produced ($(ls "$SC_OUT"/*.nvdb | head -1 | xargs basename))"
else
    bad "G10: no .nvdb produced — see g10.log"
fi

# ── G11: remote protocol with GPU server ─────────────────────────────────────
note "G11: TCP protocol round-trip with --gpu server"
PORT=$(( 20000 + RANDOM % 20000 ))
"$SC_BIN" --data-type GADGET_SIMPLE --gadget-file "$SNAP" --output-path "$SC_OUT" --port "$PORT" --gpu \
    >"$SC_OUT/g11_server.log" 2>&1 &
SRV_PID=$!
CLIENT_RC=1
for i in $(seq 1 30); do
    sleep 1
    if python3 "$HERE/protocol_client.py" localhost "$PORT" 0 0 >"$SC_OUT/g11_client.log" 2>&1; then
        CLIENT_RC=0; break
    fi
    kill -0 "$SRV_PID" 2>/dev/null || break
done
kill "$SRV_PID" 2>/dev/null; wait "$SRV_PID" 2>/dev/null
if [ "$CLIENT_RC" = 0 ]; then
    ok "GPU server round-trip"
else
    bad "G11 — see g11_server.log / g11_client.log"
fi

# ── G12: determinism / robustness extras ─────────────────────────────────────
note "G12: repeatability + empty particle type"
run_case "$SC_OUT/g12a.log" 1 --export-data 0 0 "${DENSE[@]}" --gpu; keep_raw "$SC_OUT/g12a.raw"
if python3 "$HERE/compare_raw.py" "$SC_OUT/g2_gpu.raw" "$SC_OUT/g12a.raw" 1e-3 1e-5 > "$SC_OUT/g12_cmp.log" 2>&1; then
    ok "repeated GPU run within tolerance"
else
    bad "G12 repeatability — see g12_cmp.log"
fi
if run_case "$SC_OUT/g12b.log" 1 --export-data 3 0 --gpu && [ "$(count_of "$SC_OUT/g12b.log")" = "0" ]; then
    ok "empty particle type handled cleanly (count 0, exit 0)"
else
    bad "G12 empty type: rc or count wrong — see g12b.log"
fi

# ── G13: CUB outputs (sparse needs the GPU manager; dense is host-side) ──────
note "G13: CUB output (sparse GPU + dense CPU/GPU)"
cub_count() { python3 -c "
import struct, sys
d = open(sys.argv[1], 'rb').read()
print(struct.unpack('<i', d[72:76])[0])" "$1"; }
run_case "$SC_OUT/g13s.log" 1 --export-data 0 0 --gpu --cub
CUB_S=""
for i in $(seq 1 30); do CUB_S=$(ls "$SC_OUT"/*.cub 2>/dev/null | head -1); [ -n "$CUB_S" ] && break; sleep 1; done
# The CUB finalize branch does not print the active-voxel line; compare the
# header count against the sparse voxel count established in G1
if [ -n "$CUB_S" ] && [ "$(cub_count "$CUB_S")" = "$(voxels_of "$SC_OUT/g1_gpu.log")" ]; then
    ok "sparse GPU .cub written, header count == G1 active voxels ($(cub_count "$CUB_S"))"
else
    bad "G13 sparse: file=$CUB_S count=$(cub_count "$CUB_S" 2>/dev/null) expected=$(voxels_of "$SC_OUT/g1_gpu.log")"
fi
rm -f "$SC_OUT"/*.cub
run_case "$SC_OUT/g13c.log" 1 --export-data 0 0 "${DENSE[@]}" --cub
CUB_C=""; for i in $(seq 1 30); do CUB_C=$(ls "$SC_OUT"/*.cub 2>/dev/null | head -1); [ -n "$CUB_C" ] && break; sleep 1; done
CNT_C=$(cub_count "$CUB_C" 2>/dev/null)
run_case "$SC_OUT/g13g.log" 1 --export-data 0 0 "${DENSE[@]}" --cub --gpu
CUB_G=""; for i in $(seq 1 30); do CUB_G=$(ls "$SC_OUT"/*.cub 2>/dev/null | head -1); [ -n "$CUB_G" ] && break; sleep 1; done
CNT_G=$(cub_count "$CUB_G" 2>/dev/null)
if [ -n "$CNT_C" ] && [ "$CNT_C" = "$CNT_G" ]; then
    ok "dense .cub CPU and GPU voxel counts match ($CNT_C)"
else
    bad "G13 dense: CPU count=$CNT_C GPU count=$CNT_G"
fi

# ── G14: raw-particle export with --gpu (CPU fallback) ───────────────────────
note "G14: raw particles with --gpu"
rm -f "$SC_OUT"/*.part "$SC_OUT"/*.bin
run_case "$SC_OUT/g14.log" 1 --export-data 0 0 --raw-particles 1 --gpu
if ls "$SC_OUT"/*.bin >/dev/null 2>&1 && [ "$(count_of "$SC_OUT/g14.log")" = "1000" ]; then
    ok "raw export completed under --gpu (CPU fallback, 1000 particles)"
else
    bad "G14: count=$(count_of "$SC_OUT/g14.log") — see g14.log"
fi

echo
echo "==== GPU tests: $PASS passed, $FAIL failed ===="
exit "$FAIL"
