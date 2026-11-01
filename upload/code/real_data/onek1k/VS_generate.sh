#!/bin/bash
#SBATCH --job-name=VS_generate
#SBATCH --output=log/VS_generate_%A_%a.out
#SBATCH --error=log/VS_generate_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G  
#SBATCH --partition=tier3q
#SBATCH --time=48:00:00            
#SBATCH --array=1-255

module load gcc/12.1.0
module load R/4.5.1

Rscript VS_generate.R $SLURM_ARRAY_TASK_ID