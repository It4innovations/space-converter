#!/bin/bash

# Per-rank GPU binding for Polaris (4x NVIDIA A100 per node).
# Based on the gpu_bind_kar.sh pattern; Polaris uses PBS/PALS (mpiexec),
# so prefer PALS_LOCAL_RANKID and fall back to Slurm variables.

LOCAL_RANK=$PALS_LOCAL_RANKID
if [ -z "$LOCAL_RANK" ]; then
  LOCAL_RANK=$SLURM_LOCALID
fi
if [ -z "$LOCAL_RANK" ]; then
  LOCAL_RANK=0
fi

export CUDA_VISIBLE_DEVICES=$(( LOCAL_RANK % 4 ))
echo "Process ${PALS_RANKID:-$SLURM_PROCID}: Local rank $LOCAL_RANK, Using GPU: $CUDA_VISIBLE_DEVICES"
exec "$@"
