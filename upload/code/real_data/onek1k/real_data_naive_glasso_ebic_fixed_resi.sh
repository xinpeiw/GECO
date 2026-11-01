#!/bin/bash
#SBATCH --job-name=ebic
#SBATCH --output=log/ebic_%A_%a.out
#SBATCH --error=log/ebic_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G  
#SBATCH --partition=tier3q                 
#SBATCH --time=240:00:00            
#SBATCH --array=1-8

module load gcc/12.1.0
module load R/4.5.1

Rscript real_data_naive_glasso_ebic_fixed_resi.R $SLURM_ARRAY_TASK_ID