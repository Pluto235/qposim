#!/usr/bin/env bash
#SBATCH --job-name=raptorp-lc-pilot
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/luoji/qposim/hxf_toroidal_fig/MHD3/runs/116836/postprocess/raptorp_lc/pilot/logs/%x-%A_%a.out

set -euo pipefail

base=/home/luoji/qposim/hxf_toroidal_fig/MHD3/runs/116836/postprocess/raptorp_lc
deps=/home/luoji/qposim/tools/raptorp-deps
case_file=${CASE_FILE:?CASE_FILE is required}
case_ids_file=${CASE_IDS_FILE:?CASE_IDS_FILE is required}
case_id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$case_ids_file")

export PATH="$deps/bin:$PATH"
export LD_LIBRARY_PATH="$deps/prefix/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

printf 'start=%s job=%s task=%s case=%s\n' "$(date -Is)" "$SLURM_JOB_ID" "$SLURM_ARRAY_TASK_ID" "$case_id"
python3 "$base/scripts/run_case.py" --kind pilot --cases "$case_file" --case-id "$case_id"
printf 'finish=%s job=%s task=%s case=%s\n' "$(date -Is)" "$SLURM_JOB_ID" "$SLURM_ARRAY_TASK_ID" "$case_id"
