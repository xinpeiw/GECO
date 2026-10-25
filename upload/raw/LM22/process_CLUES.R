rm(list=ls())

library(data.table)
library(readxl)

geco = readLines('data/clues/allgenes_CLUES.txt')

df <- read_excel(
  path = 'raw/LM22/41592_2015_BFnmeth3337_MOESM207_ESM.xls',
  sheet = 1,
  range = cell_limits(c(14, 1), c(NA, 3))
)
names(df) <- c("Entrez_ID", "Max_probe_set", "Gene_symbol")

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }
LM22 = ensg(unique(df$Gene_symbol))

LM22 = intersect(LM22, geco)
length(LM22)
write.table(LM22, 'data/others/geneset_clues/LM22.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

