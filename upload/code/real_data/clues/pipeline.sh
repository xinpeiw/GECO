#!/bin/bash

JOB1=$(sbatch --parsable 1_VS_generate.sh)

JOB3=$(sbatch \
    --dependency=afterok:${JOB1} \
    3_estimate.sh | awk '{print $4}')

JOB4=$(sbatch \
    --dependency=afterok:${JOB3} \
    4_merge.sh | awk '{print $4}')

JOB5=$(sbatch \
    --dependency=afterok:${JOB4} \
    5_reformat.sh | awk '{print $4}')

JOB6=$(sbatch \
    --dependency=afterok:${JOB5} \
    6_rename.sh | awk '{print $4}')

JOB7=$(sbatch \
    --dependency=afterok:${JOB6} \
    7_transform.sh | awk '{print $4}')

JOB8=$(sbatch \
    --dependency=afterok:${JOB7} \
    8_pcal.sh | awk '{print $4}')

JOB9=$(sbatch \
    --dependency=afterok:${JOB8} \
    eBIC_argmin.sh | awk '{print $4}')

