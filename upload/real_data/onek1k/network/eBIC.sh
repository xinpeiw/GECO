#!/bin/bash
#SBATCH --job-name=eBIC
#SBATCH --output=log/eBIC_%A_%a.out
#SBATCH --error=log/eBIC_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=tier2q
#SBATCH --time=240:00:00            
#SBATCH --array=1-8     

module load gcc/12.1.0
module load R/4.5.1

Rscript eBIC.R $SLURM_ARRAY_TASK_ID