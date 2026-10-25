rm(list=ls())

library(RegulationIndex)
library(data.table)

geco = readLines('data/clues/allgenes_CLUES.txt')

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

data('CD34', package='RegulationIndex')
counts = t(CD34$counts)
counts = counts[!duplicated(rownames(counts)) & !duplicated(rownames(counts), fromLast=TRUE), ]
sc     = as.data.table(Kprop_inflation_test(counts, grouping=NULL, min.depth.group=500,
                                            cell.depth.ranges=c(500,10000), genemean.cutoff=0.5, filter.gene=FALSE)$Zscores, keep.rownames='symbol')
burst  = intersect(ensg(unique(sc[qvalue<0.05]$symbol)), geco)
burst_universe = intersect(ensg(unique(sc$symbol)), geco)
write.table(burst, 'data/others/geneset_clues/inflated_asone.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
write.table(burst_universe, 'data/others/geneset_clues/inflated_universe.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

