#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --job-name=kmonty
#SBATCH --time=5-00:00:00
#SBATCH --mem=40G
#SBATCH --output=kmonty-%j.out
#SBATCH --error=kmonty-%j.err

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch job.slurm PLUTO_DATA_FILE OUTPUT_INDEX" >&2
  exit 2
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
make pluto
srun --ntasks=1 ./kmonty model_pluto.in "$1" "$2"
