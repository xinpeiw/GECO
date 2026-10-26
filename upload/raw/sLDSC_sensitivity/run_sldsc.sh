#!/bin/bash
# Inputs not included in this repository: the LDSC reference (raw/sLDSC/01_prep_reference.sh, raw/sLDSC_sensitivity/prep.sh)
# and the 117 munged FinnGen R12 summary statistics (raw/sLDSC_sensitivity/prep.sh).
set -u
LDSC="python3 raw/sLDSC/ldsc_fork/ldsc.py"
REF=raw/sLDSC/reference
ANNOT=$REF/GECO_annot_sensitivity
PLINK=$REF/1000G_EUR_Phase3_plink/1000G.EUR.QC
FRQ=$REF/1000G_Phase3_frq/1000G.EUR.QC.; HM3=$REF/hm3_snps.list
W_NOMHC=$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
W_MHC=$REF/weights_withMHC/weights.withMHC.
B_V12=$REF/baseline_v1.2/baseline_v1.2/baseline.; B_LD=$REF/baselineLD_v2.2/baselineLD.
MUNGED=raw/sLDSC_sensitivity/sumstats/munged
RES=raw/sLDSC_sensitivity/results; mkdir -p "$RES"
CTS="BIN BMem CD4ET CD4NC CD8ET CD8NC NK"

ldscore () {
  local c o
  for c in {1..22}; do o=$1.$c; [ -s "$o.l2.ldscore.gz" ] && continue
    $LDSC --l2 --bfile "$PLINK.$c" --ld-wind-cm 1 ${2:+--annot "$2.$c.annot.gz"} \
      --print-snps "$HM3" --out "$o" >"$o.l2.log" 2>&1; done; }

h2 () {
  local out=$RES/$(basename "$1" .sumstats.gz).$6; [ -s "$out.results" ] && return
  $LDSC --h2 "$1" --ref-ld-chr "$3,$2,$4" --w-ld-chr "$5" --frqfile-chr "$FRQ" \
    --overlap-annot --print-coefficients --out "$out" >"$out.log" 2>&1; }

Rscript raw/sLDSC_sensitivity/build_annot.R

for s in GECO_Shared GECO_Specific GECO_Universe GECO_Shared_500kb GECO_Specific_500kb \
         GECO_Universe_500kb GECO_Shared_top210 GECO_Specific_top210 GECO_Shared_top420 GECO_Specific_top420 \
         $(printf 'GECO_specific_%s ' $CTS); do ldscore "$ANNOT/$s" "$ANNOT/$s"; done
mkdir -p "$REF/weights_withMHC"; ldscore "$REF/weights_withMHC/weights.withMHC"

NJOBS="${NJOBS:-1}"
U=$ANNOT/GECO_Universe.
runpass () { local fn=$1 i=0
  for s in "$MUNGED"/*.sumstats.gz; do "$fn" "$s" & i=$((i + 1)); (( i % NJOBS == 0 )) && wait; done; wait; }
p_primary () { local s=$1; for hub in Shared Specific; do h2 "$s" "$ANNOT/GECO_${hub}."       "$B_V12" "$U"                          "$W_NOMHC" "${hub}_primary"; done; }
p_S1      () { local s=$1; for hub in Shared Specific; do h2 "$s" "$ANNOT/GECO_${hub}."       "$B_LD"  "$U"                          "$W_NOMHC" "${hub}_S1"; done; }
p_S3      () { local s=$1; for hub in Shared Specific; do h2 "$s" "$ANNOT/GECO_${hub}_500kb." "$B_V12" "$ANNOT/GECO_Universe_500kb." "$W_NOMHC" "${hub}_S3"; done; }
p_S4      () { local s=$1; for hub in Shared Specific; do h2 "$s" "$ANNOT/GECO_${hub}."       "$B_V12" "$U"                          "$W_MHC"   "${hub}_S4"; done; }
p_ct      () { local s=$1; for ct in $CTS; do h2 "$s" "$ANNOT/GECO_specific_$ct." "$B_V12" "$U" "$W_NOMHC" "ct_$ct"; done; }
p_nhub    () { local s=$1; for hub in Shared Specific; do for n in top210 top420; do h2 "$s" "$ANNOT/GECO_${hub}_${n}." "$B_V12" "$U" "$W_NOMHC" "${hub}_${n}"; done; done; }
for pass in p_primary p_S1 p_S3 p_S4 p_ct p_nhub; do runpass "$pass"; done
