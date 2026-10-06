#!/bin/bash
set -euo pipefail

make clean
make pluto
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 PLUTO_DATA_FILE OUTPUT_INDEX" >&2
  exit 2
fi

mpirun -np 1 ./kmonty model_pluto.in "$1" "$2"
