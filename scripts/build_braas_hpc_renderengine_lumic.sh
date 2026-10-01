#!/bin/bash
# braas-hpc-renderengine for LUMI-C (CPU only, no CUDA / GPUJPEG / epoxy).
# Run on a compute node:
#   srun --jobid=<JOBID> --overlap -w <node> -N1 -n1 -c64 bash build_braas_hpc_renderengine_lumic.sh
cd /flash/project_465002608/jaromila/space
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
