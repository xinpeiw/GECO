#!/bin/bash

BASE_DIR=$(pwd)/reps
mkdir -p "${BASE_DIR}"

for REP in $(seq 1 100); do
    REP_DIR=$(printf "%s/rep_%03d" "${BASE_DIR}" "${REP}")
    mkdir -p \
        "${REP_DIR}/data" \
        "${REP_DIR}/estimate" \
        "${REP_DIR}/merged" \
        "${REP_DIR}/transform" \
        "${REP_DIR}/ebic" \
        "${REP_DIR}/precision" \
        "${REP_DIR}/degree" \
        "${REP_DIR}/logs"

    JOB1=$(sbatch --parsable 01_generate_Y.sh "${REP}" "${REP_DIR}")
    JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} 02_estimate.sh "${REP}" "${REP_DIR}")
    JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} 03_merge.sh "${REP}" "${REP_DIR}")
    JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} 04_reformat.sh "${REP}" "${REP_DIR}")
    JOB5=$(sbatch --parsable --dependency=afterok:${JOB4} 05_transform.sh "${REP}" "${REP_DIR}")
    JOB6=$(sbatch --parsable --dependency=afterok:${JOB5} 06_pcal.sh "${REP}" "${REP_DIR}")
    JOB7=$(sbatch --parsable --dependency=afterok:${JOB6} 07_ebic.sh "${REP}" "${REP_DIR}")
    JOB8=$(sbatch --parsable --dependency=afterok:${JOB7} 08_precision.sh "${REP}" "${REP_DIR}")
    JOB9=$(sbatch --parsable --dependency=afterok:${JOB8} 09_degree.sh "${REP}" "${REP_DIR}")

    echo "rep ${REP} submitted: ${JOB1} -> ${JOB2} -> ${JOB3} -> ${JOB4} -> ${JOB5} -> ${JOB6} -> ${JOB7} -> ${JOB8} -> ${JOB9}"
done