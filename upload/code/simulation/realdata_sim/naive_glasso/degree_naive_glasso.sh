#!/bin/bash
#SBATCH --job-name=degree
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --partition=tier2q
#SBATCH --cpus-per-task=1

REP=$1
REP_DIR=$2

module load gcc/12.1.0
module load R/4.5.1
Rscript R/degree_naive_glasso.R "${REP}" "${REP_DIR}" \
  > "${REP_DIR}/logs/degree.out" \
  2> "${REP_DIR}/logs/degree.err"