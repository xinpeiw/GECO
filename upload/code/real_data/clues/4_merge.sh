#!/bin/bash

awk 'FNR==1 && NR!=1 {next} {print}' ./res_all/cov_*.txt > "./res_all/cov_all.txt"
awk 'FNR==1 && NR!=1 {next} {print}' ./res_all/se_*.txt > "./res_all/se_all.txt"

