#!/bin/bash
# Inputs not included in this repository: raw/sLDSC/reference/ (from raw/sLDSC/01_prep_reference.sh) and the 78 munged
# FinnGen R12 summary statistics in raw/sLDSC/sumstats/ (see raw/sLDSC/03_ldscores_and_h2.sh).
set -euo pipefail
REF=raw/sLDSC/reference; ANNOT=$REF/CLUES_annot_mega
RES=raw/CLUES_sLDSC/results; SUMSTATS=raw/sLDSC/sumstats
LDSC="python3 raw/sLDSC/ldsc_fork/ldsc.py"
mkdir -p "$RES"
for s in CLUES_Shared CLUES_Specific CLUES_Universe; do
  for c in $(seq 1 22); do
    o=$ANNOT/$s.$c; [ -s "$o.l2.ldscore.gz" ] && continue
    $LDSC --l2 --bfile "$REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$c" \
      --ld-wind-cm 1 --annot "$ANNOT/$s.$c.annot.gz" \
      --print-snps "$REF/hm3_snps.list" --out "$o" > "$o.l2.log" 2>&1
  done
done
for f in "$SUMSTATS"/*.sumstats.gz; do
  sid=$(basename "$f" .sumstats.gz)
  for hub in Shared Specific; do
    out=$RES/$sid.${hub}_primary; [ -s "$out.results" ] && continue
    $LDSC --h2 "$f" \
      --ref-ld-chr "$REF/baseline_v1.2/baseline_v1.2/baseline.,$ANNOT/CLUES_${hub}.,$ANNOT/CLUES_Universe." \
      --w-ld-chr "$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." \
      --frqfile-chr "$REF/1000G_Phase3_frq/1000G.EUR.QC." \
      --overlap-annot --print-coefficients --out "$out" > /dev/null 2>&1
  done
done
echo "results in $RES"
