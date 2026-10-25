# Input not included in this repository; download into raw/ImmuNexUT/ (ImmuNexUT, Ota et al. 2021, DDBJ E-GEAD-398):
#   https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-398/E-GEAD-398.processed.zip
#     (564,882,709 bytes, md5 2fc0a61db661d179d1083bb39ec605de). Unzip, extract conditional_eQTL_FDR0.05.tar, and place
#     the 28 *_conditional_eQTL_FDR0.05.txt files in raw/ImmuNexUT/conditional_eQTL_FDR0.05/.

rm(list=ls())

library(data.table)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

in_dir <- "raw/ImmuNexUT/conditional_eQTL_FDR0.05"
files  <- list.files(in_dir, pattern = "_conditional_eQTL_FDR0.05\\.txt$", full.names = TRUE)

genes <- unique(unlist(lapply(files, function(f) {
  unique(fread(f, select = "Gene_id", showProgress = FALSE)$Gene_id)
})))

genes = unique(sub('\\..*', '', genes))

genes_geco = intersect(genes, geco)
write.table(genes_geco, 'data/5a-cis_molGenes/cis_eGene/ImmuNexUT.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
genes_clues = intersect(genes, clues)
write.table(genes_clues, 'data/others/molQTL_clues/cis/ImmuNexUT.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
