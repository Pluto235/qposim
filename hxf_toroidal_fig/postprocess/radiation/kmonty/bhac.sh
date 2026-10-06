#!/bin/bash
set -euo pipefail

make clean
make bhac
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-32}"

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 BHAC_DATA_FILE OUTPUT_INDEX" >&2
  exit 2
fi

mpirun -np 1 ./kmonty model_bhac.in "$1" "$2"
