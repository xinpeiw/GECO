#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=240:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

REP=$1
REP_DIR=$2

N=$(ls "${REP_DIR}/estimate" | wc -l)
echo "rep ${REP}: found ${N} estimate files"

awk 'FNR==1 && NR!=1 {next} {print}' "${REP_DIR}"/estimate/cov_*.txt > "${REP_DIR}/merged/cov.txt"
awk 'FNR==1 && NR!=1 {next} {print}' "${REP_DIR}"/estimate/se_*.txt  > "${REP_DIR}/merged/se.txt"