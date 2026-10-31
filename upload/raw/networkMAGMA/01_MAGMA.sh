#!/bin/bash
# Third-party inputs not included in this repository (place under raw/networkMAGMA/):
#   magma/magma                      MAGMA executable, https://cncr.nl/research/magma/
#   magma/g1000_eur.{bed,bim,fam}    European 1000 Genomes reference panel distributed with MAGMA (same site)
#   magma/snp_2_gene_10k.genes.annot SNP-to-gene annotation from `magma --annotate` with a 10 kb gene window
#   raw/sLDSC/sumstats/              the 78 munged FinnGen R12 summary statistics (raw/sLDSC/03_ldscores_and_h2.sh header)
# Output: raw/networkMAGMA/magma/results/<trait>.genes.out (MAGMA gene-level statistics, one file per GWAS)
set -euo pipefail
M=raw/networkMAGMA/magma; IN=$M/input; OUT=$M/results; mkdir -p "$IN" "$OUT"
for gwas in raw/sLDSC/sumstats/*.sumstats.gz; do
  trait=$(basename "$gwas" .sumstats.gz); infile="$IN/$trait.magma.txt"
  [ -s "$OUT/$trait.genes.out" ] && continue
  zcat "$gwas" | Rscript -e '
    dat <- read.table(file("stdin"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
    keep <- !is.na(dat$SNP) & dat$SNP != "" & dat$SNP != "." & !is.na(dat$Z) & is.finite(dat$Z) & !is.na(dat$N) & dat$N > 0
    dat <- dat[keep, ]
    write.table(data.frame(SNP = dat$SNP, P = 2 * pnorm(abs(dat$Z), lower.tail = FALSE), N = dat$N),
                file = stdout(), sep = "\t", quote = FALSE, row.names = FALSE)' > "$infile"
  "$M/magma" --bfile "$M/g1000_eur" --pval "$infile" use=SNP,P ncol=N \
    --gene-annot "$M/snp_2_gene_10k.genes.annot" --out "$OUT/$trait"
done
