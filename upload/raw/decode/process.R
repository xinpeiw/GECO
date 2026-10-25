rm(list=ls())

library(data.table)
library(readxl)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

xlsx_file <- "raw/decode/41588_2021_978_MOESM4_ESM.xlsx"

st1 = as.data.table(read_excel(xlsx_file, sheet = "ST01", skip = 2))
st1 = st1[, c(1, 6, 10)]
setnames(st1, c("SeqId", "included", "ensembl"))
st1 = st1[included == 'Yes',]
seq2ensg = setNames(st1$ensembl, st1$SeqId)

universe_geco = intersect(unique(st1[grepl('^ENSG[0-9]+$', ensembl)]$ensembl), geco)
out_uni = c('data/5a-cis_molGenes/universe_deCODE.txt', 'data/6a-trans_molGenes/universe_deCODE.txt')
for (o in out_uni) write.table(universe_geco, o, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

universe_clues = intersect(unique(st1[grepl('^ENSG[0-9]+$', ensembl)]$ensembl), clues)
out_uni = 'data/others/molQTL_clues/universe_deCODE.txt'
for (o in out_uni) write.table(universe_clues, o, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

st2 = as.data.table(read_excel(xlsx_file, sheet = "ST02", skip = 2))
st2 = st2[, c(3, 10, 11, 22, 28)]
setnames(st2, c("gene", "SeqId", "variant", "cistrans", "mlog10p"))
st2[, ensg := unname(seq2ensg[SeqId])]

st2_geco = st2[grepl('^ENSG[0-9]+$', ensg) & ensg %in% geco]
cis_geco   = sort(unique(st2_geco[cistrans == 'cis']$ensg))
trans_geco = sort(unique(st2_geco[cistrans == 'trans']$ensg))
out_cis   = 'data/5a-cis_molGenes/cis_pGene/deCODE.txt'
out_trans = 'data/6a-trans_molGenes/trans_pGene/deCODE.txt'
write.table(cis_geco,   out_cis,   row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
write.table(trans_geco, out_trans, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

st2_clues = st2[grepl('^ENSG[0-9]+$', ensg) & ensg %in% clues]
cis_clues   = sort(unique(st2_clues[cistrans == 'cis']$ensg))
out_cis   = 'data/others/molQTL_clues/cis/deCODE.txt'
write.table(cis_clues,   out_cis,   row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
