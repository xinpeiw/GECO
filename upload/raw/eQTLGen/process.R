# Inputs not included in this repository; download into raw/eQTLGen/ (eQTLGen Consortium, Vosa et al. 2021):
#   2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz  (322,775,879 bytes)
#     https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
#   2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz  (1,839,918 bytes)
#     https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz
#   cis_tested_genes.txt and trans_tested_genes.txt are included; they are regenerated from the full result files if deleted.
# media-2.xlsx is the supplementary table file of eQTLGen Phase 2 (Warmerdam et al. 2026, medRxiv 2026.02.04.26343575);
#   Table S6 = all colocalizing cis-eQTL/trans-eQTL credible-set pairs (PP4 >= 0.8, lead-variant LD R2 >= 0.8), Table S3 = per-gene summary.

rm(list=ls())

library(data.table)
library(readxl)
cis_file   <- "raw/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
trans_file <- "raw/eQTLGen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz"

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

cis_out    <- 'data/5a-cis_molGenes/cis_eGene/eQTLGen.txt'
trans_out  <- 'data/6a-trans_molGenes/trans_eGene/eQTLGen.txt'

extract_egenes <- function(path) {
  dt <- data.table::fread( cmd    = paste("gzip -dc", shQuote(path)),
    select = c("Gene"), sep    = "\t", header = TRUE, showProgress = FALSE )
  unique(dt$Gene)
}

cis = extract_egenes(cis_file)
cis_geco <- intersect(cis, geco)
write.table(cis_geco, cis_out, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

cis_clues = intersect(cis, clues)
write.table(cis_clues, 'data/others/molQTL_clues/cis/cis_eQTLGen.txt',
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

trans = extract_egenes(trans_file)
trans_geco <- intersect(trans, geco)
write.table(trans_geco, trans_out, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

cis_full_url   <- "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
trans_full_url <- "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz"

tested_genes <- function(url, cache) {
  if (!file.exists(cache)) {
    system(sprintf(
      "curl -sL %s | gzip -dc | awk -F'\\t' 'NR>1{a[$8]} END{for(k in a) print k}' | sed 's/[.].*//' | sort -u > %s",
      shQuote(url), shQuote(cache)))
  }
  readLines(cache)
}

cis_universe_out   <- 'data/5a-cis_molGenes/eQTLGen_universe_cis.txt'
trans_universe_out <- 'data/6a-trans_molGenes/eQTLGen_universe_trans.txt'

cis_universe   <- tested_genes(cis_full_url,   "raw/eQTLGen/cis_tested_genes.txt")
trans_universe <- tested_genes(trans_full_url, "raw/eQTLGen/trans_tested_genes.txt")

cis_universe_geco = intersect(cis_universe, geco)
trans_universe_geco = intersect(trans_universe, geco)
write.table(cis_universe_geco,   cis_universe_out,   sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(trans_universe_geco, trans_universe_out, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

cis_universe_clues = intersect(cis_universe, clues)
write.table(cis_universe_clues, 'data/others/molQTL_clues/eQTLGen_universe_cis.txt',
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

genes <- as.data.table(read_excel('raw/eQTLGen/media-2.xlsx', sheet = "Table S3", skip = 2))$'ENSEMBL gene ID'
pairs = as.data.table(read_excel('raw/eQTLGen/media-2.xlsx', sheet = "Table S6", skip = 2))
pairs = unique(pairs[,c('cis eGene','trans eGene')])
colnames(pairs) = c('cis','trans')

inter = intersect(genes, geco)
pairs_geco = pairs[which(pairs$cis %in% inter & pairs$trans %in% inter),]
for (o in c('data/6b-edge/cis_trans_pair/eQTLGen.txt', 'data/6c-hub_position/cis_trans_pair/eQTLGen.txt')) write.table(pairs_geco, o,
            row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
for (o in c('data/6b-edge/cis_trans_pair/eQTLGen_universe_pair.txt', 'data/6c-hub_position/cis_trans_pair/eQTLGen_universe_pair.txt')) write.table(inter, o,
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
