#!/bin/bash
#SBATCH --job-name=gen_Y
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=240:00:00
#SBATCH --partition=tier2q  
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1

REP=$1
REP_DIR=$2

module load gcc/12.1.0
module load R/4.5.1
Rscript R/generate_Y.R "${REP}" "${REP_DIR}" \
  > "${REP_DIR}/logs/generate_Y.out" \
  2> "${REP_DIR}/logs/generate_Y.err"