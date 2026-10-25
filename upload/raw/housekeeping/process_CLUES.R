rm(list=ls())

library(data.table)

df = fread('raw/housekeeping/Housekeeping_GenesHuman.csv')
geco = readLines('data/clues/allgenes_CLUES.txt')

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

hk = ensg(unique(df$Gene.name))
length(unique(hk)); length(unique(df$Gene.name))

hk = intersect(hk, geco)
write.table(hk, 'data/others/geneset_clues/housekeeping.txt',  row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
