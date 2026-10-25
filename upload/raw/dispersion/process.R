rm(list=ls())

library(data.table)
library(readxl)

geco = fread('data/onek1k/GECO_universe.txt', header = FALSE)$V1

disp  = as.data.table(read_excel('raw/dispersion/media-1.xlsx', sheet='dispersion'))
disp_universe = unique(disp$gene)
disp_cols = setdiff(names(disp), "gene")
cells = unique(sub("^NA[0-9]+_", "", disp_cols))
ind   = unique(sub("_.*$", "", disp_cols))
med   = sapply(cells, function(ct) apply(disp[, paste(ind,ct,sep='_'), with=FALSE], 1, median))
high_always = disp$gene[rowSums(apply(med, 2, function(v) v >= quantile(v, .8, na.rm=TRUE)), na.rm=TRUE) == length(cells)]
low_always = disp$gene[rowSums(apply(med, 2, function(v) v <= quantile(v, .2, na.rm=TRUE)), na.rm=TRUE) == length(cells)]

car   = as.data.table(read_excel('raw/dispersion/media-3.xlsx', sheet='cardiac_HDC_results'))
dd = car$gene[rowSums(as.matrix(car[, grep('^dispersion_adj\\.P\\.Val_', names(car), value=TRUE), with=FALSE]) < 0.05, na.rm=TRUE) >= 1]

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

high_always_ensg = ensg(unique(high_always))
low_always_ensg = ensg(unique(low_always))
dd_ensg = ensg(unique(dd))
disp_universe_ensg = ensg(unique(disp_universe))

high_always_ensg = intersect(high_always_ensg, geco)
low_always_ensg = intersect(low_always_ensg, geco)
dd_ensg = intersect(dd_ensg, geco)
universe = intersect(disp_universe_ensg, geco)

write.table(high_always_ensg, 'data/3a-benchmark_genesets/high_dispersion.txt',row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
write.table(low_always_ensg, 'data/3a-benchmark_genesets/low_dispersion.txt',row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
write.table(dd_ensg, 'data/3a-benchmark_genesets/dd.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
write.table(universe, 'data/3a-benchmark_genesets/disp_universe.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
