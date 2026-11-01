#!/bin/bash
#SBATCH --job-name=vs12
#SBATCH --output=log/vs12_%A_%a.out
#SBATCH --error=log/vs12_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G  
#SBATCH --partition=tier1q                 
#SBATCH --time=240:00:00            
#SBATCH --array=1-255

module load gcc/12.1.0
module load R/4.5.1

Rscript VS12.R $SLURM_ARRAY_TASK_ID