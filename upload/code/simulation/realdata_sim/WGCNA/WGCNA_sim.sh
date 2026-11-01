#!/bin/bash
#SBATCH --job-name=WGCNA
#SBATCH --output=log/WGCNA_%A_%a.out
#SBATCH --error=log/WGCNA_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=tier1q
#SBATCH --time=10:00:00            
#SBATCH --array=1-800

module load gcc/12.1.0
module load R/4.5.1

Rscript WGCNA_sim.R $SLURM_ARRAY_TASK_ID