rm(list=ls())

library(data.table)
library(readxl)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

xlsx_file <- "raw/UKB-PPP/41586_2023_6592_MOESM3_ESM.xlsx"
st9 <- as.data.table(read_excel(xlsx_file, sheet = "ST9", skip = 4))

cis_genes   = unique(st9[`cis/trans` == 'cis']$`Assay Target`)
trans_genes = unique(st9[`cis/trans` == 'trans']$`Assay Target`)

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

cis_ensg   = intersect(ensg(unique(cis_genes)),   geco)
trans_ensg = intersect(ensg(unique(trans_genes)), geco)
out_cis   = 'data/5a-cis_molGenes/cis_pGene/UKB-PPP.txt'
out_trans = 'data/6a-trans_molGenes/trans_pGene/UKB-PPP.txt'
write.table(cis_ensg,   out_cis,   row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
write.table(trans_ensg, out_trans, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

cis_clues   = intersect(ensg(unique(cis_genes)), clues)
out_cis   = 'data/others/molQTL_clues/cis/UKB-PPP.txt'
write.table(cis_clues,   out_cis,   row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

st3 <- as.data.table(read_excel(xlsx_file, sheet = "ST3", skip = 2))
universe_genes = unique(st3$`Gene symbol`)
universe_geco  = intersect(ensg(unique(universe_genes)), geco)
universe_clues  = intersect(ensg(unique(universe_genes)), clues)

out_uni = c('data/5a-cis_molGenes/universe_UKB-PPP.txt', 'data/6a-trans_molGenes/universe_UKB-PPP.txt')
for (o in out_uni) write.table(universe_geco, o, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
out_uni = 'data/others/molQTL_clues/universe_UKB-PPP.txt'
for (o in out_uni) write.table(universe_clues, o, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
