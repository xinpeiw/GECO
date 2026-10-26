# Inputs not included in this repository: raw/GWASgene/study_index.tsv (raw/GWASgene/build_ot_inputs.py) and
# raw/CLUES_GWASgene/gene_locus_L2G_allgenes.tsv.gz (raw/CLUES_GWASgene/build_l2g_table.py), both from Open Targets Platform 26.03.
suppressMessages(library(data.table))

U     <- readLines('data/clues/allgenes_CLUES.txt')
panel <- fread('data/5b-ldsc/sensitivity/FinnGen_phenotypes_final.csv')
si    <- fread('raw/GWASgene/study_index.tsv', select = c('studyId', 'efo_id'))
panel <- merge(panel, si, by.x = 'FinnGen ID', by.y = 'studyId')
GROUPS <- split(panel$efo_id, panel$Category)

glf <- fread('raw/CLUES_GWASgene/gene_locus_L2G_allgenes.tsv.gz',
             select = c('geneId', 'efo_id', 'studyId', 'L2G_score'))
glf <- glf[geneId %in% U]

L2G_THR <- 0.8

for (grp in names(GROUPS)) {
  sub   <- glf[efo_id %in% GROUPS[[grp]]]
  genes <- unique(sub[L2G_score > L2G_THR, geneId])
  writeLines(sort(genes), file.path('data/7b-GWASgene', sprintf('genelist_%s.txt', gsub('[^A-Za-z]', '', grp))))
}
