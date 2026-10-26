#!/bin/bash
# Third-party LDSC reference, not included in this repository (~1.7 GB download, ~6 GB unpacked).
# This script downloads it into raw/sLDSC/reference/ from
#   https://zenodo.org/records/10515792 (1000G Phase 3 plink files, baseline v1.2 LD scores,
#   HapMap3 no-MHC regression weights, allele frequencies) and https://zenodo.org/records/7773502 (w_hm3.snplist).
set -euo pipefail
REF=raw/sLDSC/reference; mkdir -p "$REF"; cd "$REF"
B=https://zenodo.org/records/10515792/files
for f in 1000G_Phase3_plinkfiles.tgz 1000G_Phase3_baseline_v1.2_ldscores.tgz \
         1000G_Phase3_weights_hm3_no_MHC.tgz 1000G_Phase3_frq.tgz; do
  [ -s "$f" ] || curl -fsSL -O "$B/$f"
done
[ -s w_hm3.snplist ] || { curl -fsSL -o w_hm3.snplist.gz https://zenodo.org/records/7773502/files/w_hm3.snplist.gz; gunzip -f w_hm3.snplist.gz; }
[ -d 1000G_EUR_Phase3_plink ] || tar -xzf 1000G_Phase3_plinkfiles.tgz
[ -d baseline_v1.2 ] || { mkdir -p baseline_v1.2 && tar -xzf 1000G_Phase3_baseline_v1.2_ldscores.tgz -C baseline_v1.2; }
[ -d 1000G_Phase3_weights_hm3_no_MHC ] || tar -xzf 1000G_Phase3_weights_hm3_no_MHC.tgz
[ -d 1000G_Phase3_frq ] || tar -xzf 1000G_Phase3_frq.tgz
awk 'NR>1{print $1}' w_hm3.snplist > hm3_snps.list
echo "reference ready in $REF"
