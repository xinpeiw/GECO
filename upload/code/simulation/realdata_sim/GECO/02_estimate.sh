#!/bin/bash
#SBATCH --job-name=est
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --error=slurm_logs/%x_%A_%a.err
#SBATCH --array=1-255
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=tier1q
#SBATCH --exclude=cri22cn115
#SBATCH --time=240:00:00            

REP=$1
REP_DIR=$2
TASK_ID=${SLURM_ARRAY_TASK_ID}

module load gcc/12.1.0
module load R/4.5.1
Rscript R/estimate.R "${REP}" "${TASK_ID}" "${REP_DIR}" \
  > "${REP_DIR}/logs/estimate_${TASK_ID}.out" \
  2> "${REP_DIR}/logs/estimate_${TASK_ID}.err"