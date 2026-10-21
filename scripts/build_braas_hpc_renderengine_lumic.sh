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
# braas-hpc-renderengine for LUMI-C (CPU only, no CUDA / GPUJPEG / epoxy).
# Run on a compute node:
#   srun --jobid=<JOBID> --overlap -w <node> -N1 -n1 -c64 bash build_braas_hpc_renderengine_lumic.sh
cd ${SPACE_ROOT:-${PWD}}   # workspace root with src/, build/, install/ (default: current directory)
ROOT_DIR=${PWD}

ml purge >/dev/null 2>&1
ml LUMI/25.09 partition/C >/dev/null 2>&1
ml PrgEnv-gnu buildtools/25.09 >/dev/null 2>&1

export CC=cc
export CXX=CC

output=${ROOT_DIR}/install/braas-hpc-renderengine_lumic
src=${ROOT_DIR}/src

mkdir -p ${ROOT_DIR}/build/braas-hpc-renderengine_lumic
cd ${ROOT_DIR}/build/braas-hpc-renderengine_lumic

make_d="${src}/braas-hpc-renderengine"
make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"
make_d="${make_d} -DWITH_CLIENT_GPUJPEG=OFF"
make_d="${make_d} -DWITH_CLIENT_EPOXY=OFF"

cmake ${make_d}
make -j 64
make install
