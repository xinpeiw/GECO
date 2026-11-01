#!/bin/bash

BASE_DIR=$(pwd)/reps
mkdir -p "${BASE_DIR}"

for REP in $(seq 1 100); do
    REP_DIR=$(printf "%s/rep_%03d" "${BASE_DIR}" "${REP}")

    JOB7=$(sbatch --parsable naive_glasso_ebic.sh "${REP}" "${REP_DIR}")
    JOB8=$(sbatch --parsable --dependency=afterok:${JOB7} naive_glasso_precision.sh "${REP}" "${REP_DIR}")
    JOB9=$(sbatch --parsable --dependency=afterok:${JOB8} degree_naive_glasso.sh "${REP}" "${REP_DIR}")

done