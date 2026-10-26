#!/bin/bash
# Third-party inputs not included in this repository; this script downloads them:
#   LDSC reference into raw/sLDSC/reference/ from https://zenodo.org/records/10515792 (1000G Phase 3 plink files,
#     baseline-LD v2.2 LD scores, HapMap3 no-MHC weights, allele frequencies) and https://zenodo.org/records/7773502
#     (w_hm3.snplist); the baseline v1.2 LD scores come from raw/sLDSC/01_prep_reference.sh.
#   FinnGen R12 summary statistics for the 117 GWAS in raw/sLDSC_sensitivity/finngen_neff_117.tsv, from
#     https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_<endpoint>.gz
#     (the files used here were downloaded on 2026-06-22), munged to LDSC format with effective
#     N = 4*Ncase*Ncontrol/(Ncase+Ncontrol); case/control counts are from the Open Targets Platform 25.03 study index.
set -u
REF=raw/sLDSC/reference; HM3=$REF/w_hm3.snplist
RAW=raw/sLDSC_sensitivity/sumstats/raw
MUN=raw/sLDSC_sensitivity/sumstats/munged
LOGD=raw/sLDSC_sensitivity/sumstats/munge_logs; mkdir -p "$REF" "$RAW" "$MUN" "$LOGD"
NEFF=raw/sLDSC_sensitivity/finngen_neff_117.tsv
MUNGE="python3 raw/sLDSC/ldsc_fork/munge_sumstats.py"

B=https://zenodo.org/records/10515792/files
fetch () { [ -e "$REF/$1" ] || curl -fsSL -o "$REF/$1" "$B/$1"; }
[ -f "$HM3" ] || { curl -fsSL -o "$REF/w_hm3.snplist.gz" https://zenodo.org/records/7773502/files/w_hm3.snplist.gz; gunzip "$REF/w_hm3.snplist.gz"; }
awk 'NR>1{print $1}' "$REF/w_hm3.snplist" > "$REF/hm3_snps.list"
fetch 1000G_Phase3_plinkfiles.tgz;               [ -d "$REF/1000G_EUR_Phase3_plink" ]         || tar -xzf "$REF/1000G_Phase3_plinkfiles.tgz" -C "$REF"
fetch 1000G_Phase3_baselineLD_v2.2_ldscores.tgz; [ -d "$REF/baselineLD_v2.2" ]                || { mkdir -p "$REF/baselineLD_v2.2"; tar -xzf "$REF/1000G_Phase3_baselineLD_v2.2_ldscores.tgz" -C "$REF/baselineLD_v2.2"; }
fetch 1000G_Phase3_weights_hm3_no_MHC.tgz;       [ -d "$REF/1000G_Phase3_weights_hm3_no_MHC" ]|| tar -xzf "$REF/1000G_Phase3_weights_hm3_no_MHC.tgz" -C "$REF"
fetch 1000G_Phase3_frq.tgz;                      [ -d "$REF/1000G_Phase3_frq" ]               || tar -xzf "$REF/1000G_Phase3_frq.tgz" -C "$REF"

BK=https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release
for sid in $(awk -F'\t' 'NR>1{print $1}' "$NEFF"); do ep=${sid#FINNGEN_R12_}
  [ -s "$RAW/$sid.gz" ]           || curl -fsSL -o "$RAW/$sid.gz" "$BK/finngen_R12_$ep.gz" || continue
  [ -s "$MUN/$sid.sumstats.gz" ]  && continue
  neff=$(awk -F'\t' -v s=$sid '$1==s{print $5; exit}' "$NEFF"); [ -n "$neff" ] || continue
  $MUNGE --sumstats "$RAW/$sid.gz" --out "$MUN/$sid" --snp rsids --a1 alt --a2 ref \
    --p pval --signed-sumstats beta,0 --N "$neff" --merge-alleles "$HM3" >"$LOGD/$sid.log" 2>&1
done
