#!/bin/bash
#SBATCH --job-name=reformat
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=tier3q
#SBATCH --time=240:00:00


REP=$1
REP_DIR=$2

module load gcc/12.1.0
module load R/4.5.1
Rscript R/reformat.R "${REP}" "${REP_DIR}" \
  > "${REP_DIR}/logs/reformat.out" \
  2> "${REP_DIR}/logs/reformat.err"

mv "${REP_DIR}"/merged/cov_alpha12.rds "${REP_DIR}"/merged/cov_alpha.rds
mv "${REP_DIR}"/merged/cov_g12_c1.rds "${REP_DIR}"/merged/cov_c1_CD4ET.rds
mv "${REP_DIR}"/merged/cov_g12_c2.rds "${REP_DIR}"/merged/cov_c2_NK.rds
mv "${REP_DIR}"/merged/cov_g12_c3.rds "${REP_DIR}"/merged/cov_c3_CD4NC.rds
mv "${REP_DIR}"/merged/cov_g12_c4.rds "${REP_DIR}"/merged/cov_c4_CD8ET.rds
mv "${REP_DIR}"/merged/cov_g12_c5.rds "${REP_DIR}"/merged/cov_c5_BIN.rds
mv "${REP_DIR}"/merged/cov_g12_c6.rds "${REP_DIR}"/merged/cov_c6_CD8NC.rds
mv "${REP_DIR}"/merged/cov_g12_c7.rds "${REP_DIR}"/merged/cov_c7_BMem.rds

mv "${REP_DIR}"/merged/se_alpha12.rds "${REP_DIR}"/merged/se_alpha.rds
mv "${REP_DIR}"/merged/se_g12_c1.rds "${REP_DIR}"/merged/se_c1_CD4ET.rds
mv "${REP_DIR}"/merged/se_g12_c2.rds "${REP_DIR}"/merged/se_c2_NK.rds
mv "${REP_DIR}"/merged/se_g12_c3.rds "${REP_DIR}"/merged/se_c3_CD4NC.rds
mv "${REP_DIR}"/merged/se_g12_c4.rds "${REP_DIR}"/merged/se_c4_CD8ET.rds
mv "${REP_DIR}"/merged/se_g12_c5.rds "${REP_DIR}"/merged/se_c5_BIN.rds
mv "${REP_DIR}"/merged/se_g12_c6.rds "${REP_DIR}"/merged/se_c6_CD8NC.rds
mv "${REP_DIR}"/merged/se_g12_c7.rds "${REP_DIR}"/merged/se_c7_BMem.rds