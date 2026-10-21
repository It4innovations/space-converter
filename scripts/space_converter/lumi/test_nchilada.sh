#!/bin/bash

#####################################################################################################################
# Copyright(C) 2023-2026 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
#
# This program is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#####################################################################################################################
#
# CHANGA_NCHILADA check on LUMI-G: does the dense pipeline produce the same result
# on the CPU and on the AMD GPUs, and is the gas data actually being read?
#
# Usage: test_nchilada.sh <jobid> [ranks]
#
# The NChilada reader keeps every block in memory, so size the allocation to the
# snapshot: cosmo25.2304g.nc is ~547 GB and needs 2 LUMI-G nodes minimum, 4 to be safe.
#

JOBID=${1:?usage: test_nchilada.sh <jobid> [ranks]}
RANKS=${2:-128}

ml purge
ml LUMI/25.09 partition/G
ml PrgEnv-gnu
ml rocm/6.4.4

ROOT=/pfs/lustref1/flash/project_465002608/jaromila/space
NC=/scratch/project_465002608/jaromila/space/data/cosmo25.2304g.nc
BIN=${ROOT}/install/space_converter_lumi_hip/bin/space_converter
GPUSH=${ROOT}/src/space-converter/scripts/space_converter/lumi/gpu.sh
OUT=/scratch/project_465002608/jaromila/space/nch_test

export LD_LIBRARY_PATH=${ROOT}/install/lib-linux_x64/openvdb/lib:${ROOT}/install/lib-linux_x64/tbb/lib:${ROCM_PATH}/lib:${LD_LIBRARY_PATH}
export MPICH_GPU_SUPPORT_ENABLED=1

# The allocation's inherited CPU mask does not fit an arbitrary task count; let Slurm
# place the ranks freely instead.
SRUN="srun --overlap --jobid=${JOBID} --cpu-bind=none --export=ALL -u"

run_case() {
	local name=$1 ; shift
	local o=${OUT}/${name}
	mkdir -p ${o} ; rm -f ${o}/*
	echo "=========== ${name} ==========="
	${SRUN} --ntasks=${RANKS} -t 90:00 "$@" --output-path ${o} 2>&1 \
		| grep -E "Warning|minI|Active voxels|particles_count:|finished|elapsed time: |error|Error"
	ls -l ${o}
}

COMMON="--data-type CHANGA_NCHILADA --nc-dir ${NC} --grid-dim 100 --dense-type 6 --export-data 0 1"

# CPU reference
run_case dense6_cpu ${BIN} ${COMMON}

# AMD GPU (HIP), one rank per GCD
run_case dense6_gpu ${GPUSH} ${BIN} ${COMMON} --gpu

echo "=========== porovnani ==========="
ls -l ${OUT}/dense6_cpu ${OUT}/dense6_gpu
