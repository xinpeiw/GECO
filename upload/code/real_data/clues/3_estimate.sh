#!/bin/bash
#SBATCH --job-name=estimate
#SBATCH --output=log/estimate_%A_%a.out
#SBATCH --error=log/estimate_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G  
#SBATCH --partition=tier2q                 
#SBATCH --time=240:00:00            
#SBATCH --array=1-241

module load gcc/12.1.0
module load R/4.5.1

Rscript 3_estimate.R $SLURM_ARRAY_TASK_ID