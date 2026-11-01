#!/bin/bash
#SBATCH --job-name=precision
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1-8
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --partition=tier2q
#SBATCH --cpus-per-task=1

REP=$1
REP_DIR=$2
TASK_ID=${SLURM_ARRAY_TASK_ID}

module load gcc/12.1.0
module load R/4.5.1
Rscript R/precision.R "${REP}" "${TASK_ID}" "${REP_DIR}" \
  > "${REP_DIR}/logs/precision_${TASK_ID}.out" \
  2> "${REP_DIR}/logs/precision_${TASK_ID}.err"