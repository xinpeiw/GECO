#!/bin/bash
#SBATCH --job-name=precision
#SBATCH --output=log/precision_%A_%a.out
#SBATCH --error=log/precision_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G  
#SBATCH --partition=tier2q
#SBATCH --time=24:00:00
#SBATCH --array=1-8

module load gcc/12.1.0
module load R/4.5.1

Rscript precision.R $SLURM_ARRAY_TASK_ID