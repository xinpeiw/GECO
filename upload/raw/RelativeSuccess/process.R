# Inputs not included in this repository; download into raw/RelativeSuccess/ from Minikel et al. 2024
# (Nature 629:624-629; github.com/ericminikel/genetic_support, CC-BY-4.0; Citeline Pharmaprojects subset released with permission),
# https://raw.githubusercontent.com/ericminikel/genetic_support/master/data/<file> (md5, bytes):
#   merge2.tsv.gz   a67fce1d9db3973fdc6d74ab9bf95d79  21,012,587   target-indication programmes x genetic associations
#   assoc.tsv.gz    1cf266cb122f7355ae24050ed89703e8  16,188,145   genetic associations
#   indic.tsv       b66877769f6f0fea80b25b43877adbcd      63,260   indications (therapy areas, genetic insight)
#   meta_hcat.tsv   910100518bf2633bc541fcc7deb8ad8c          66   phase category orderings
#   meta_acat.tsv   29876d8edd8a9af95a17497d7ca693c5          79
#   meta_ccat.tsv   910100518bf2633bc541fcc7deb8ad8c          66
# minikel_fns.R holds pipeline_best() and abs_or() copied verbatim from their src/gensup_analysis.R.
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
options(dplyr.summarise.inform = FALSE)

DD <- 'raw/RelativeSuccess'
merge2 <- read_tsv(file.path(DD, 'merge2.tsv.gz'), col_types = cols(), progress = FALSE)
assoc  <- read_tsv(file.path(DD, 'assoc.tsv.gz'),  col_types = cols(), progress = FALSE)
indic  <- read_tsv(file.path(DD, 'indic.tsv'), col_types = cols())
meta_hcat <- read_tsv(file.path(DD, 'meta_hcat.tsv'), col_types = cols())
meta_acat <- read_tsv(file.path(DD, 'meta_acat.tsv'), col_types = cols())
meta_ccat <- read_tsv(file.path(DD, 'meta_ccat.tsv'), col_types = cols())
source(file.path(DD, 'minikel_fns.R'))

U  <- unique(trimws(readLines('data/onek1k/GECO_universe.txt')))
Uc <- unique(trimws(readLines('data/clues/allgenes_CLUES.txt')))
hg <- fread('raw/hgnc_complete_set.txt', select = c('symbol', 'ensembl_gene_id'))[ensembl_gene_id != '']
s2e <- setNames(hg$ensembl_gene_id, hg$symbol)

best <- as.data.table(pipeline_best(merge2, phase = 'combined', basis = 'ti', verbose = FALSE))
b <- best[target_status %in% c('genetically supported target', 'unsupported target')]
b[, `:=`(genetic_support = target_status == 'genetically supported target', ensg = s2e[gene])]
imm <- b[grepl('immune', areas)]

pairs <- function(x) { out <- x[, .(ti_uid, gene, ensg, indication_mesh_id, indication_mesh_term, areas, genetic_support, similarity,
                                    assoc_source, ccat, succ_p_1, succ_1_2, succ_2_3, succ_3_a)]; setorder(out, ti_uid); out }
fwrite(pairs(imm[ensg %in% U]),  'data/5f-RS/immune_TI_pairs.tsv', sep = '\t')
fwrite(pairs(imm[ensg %in% Uc]), 'data/others/immune_TI_pairs_CLUES.tsv', sep = '\t')
