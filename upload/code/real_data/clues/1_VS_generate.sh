#!/bin/bash
#SBATCH --job-name=VS_generate
#SBATCH --output=log/VS_generate_%A_%a.out
#SBATCH --error=log/VS_generate_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G  
#SBATCH --partition=tier2q                 
#SBATCH --time=240:00:00            
#SBATCH --array=1-241

module load gcc/12.1.0
module load R/4.5.1

Rscript 1_VS_generate.R $SLURM_ARRAY_TASK_ID