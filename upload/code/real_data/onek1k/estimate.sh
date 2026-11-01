#!/bin/bash
#SBATCH --job-name=estimate
#SBATCH --output=log/estimate_%A_%a.out
#SBATCH --error=log/estimate_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=tier3q                 
#SBATCH --time=48:00:00
#SBATCH --array=1-255

module load gcc/12.1.0
module load R/4.5.1

Rscript 3_estimate.R $SLURM_ARRAY_TASK_ID