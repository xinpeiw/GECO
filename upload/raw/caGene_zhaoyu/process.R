rm(list=ls())

library(data.table)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

tbl_dir = 'raw/caGene_zhaoyu'

st5 = fread(file.path(tbl_dir, "Supplementary Table 5 Cell type specific caQTLs.csv"))
ca_peaks = unique(st5$QTL_caPeak)

st2 = fread(file.path(tbl_dir, "Supplementary Table 2 Peak gene co-activity.csv"))
genes = unique(st2[peak %in% ca_peaks]$linkedGene)

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

caGenes_geco = intersect(ensg(unique(genes)), geco)
out = 'data/5a-cis_molGenes/caGene/zhaoyu.txt'
write.table(caGenes_geco, out, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

caGenes_clues = intersect(ensg(unique(genes)), clues)
out = 'data/others/molQTL_clues/cis/ca_intestine.txt'
write.table(caGenes_clues, out, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
