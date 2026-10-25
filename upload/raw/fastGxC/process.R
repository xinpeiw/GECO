rm(list=ls())

library(data.table)
library(readxl)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

xlsx = "raw/fastGxC/1-s2.0-S2666979X26001126-mmc2.xlsx"
gtex = as.data.table(read_excel(xlsx, sheet = "GTEx"))

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

get_eGenes = function(d, type) intersect(ensg(unique(d[eGene_type == type]$eGene)), geco)

gtex_shared   = get_eGenes(gtex, "FastGxC sh-eGene")
gtex_specific = get_eGenes(gtex, "FastGxC sp-eGene")
gtex_universe = union(gtex_shared, gtex_specific)

outdir = 'data/5a-cis_molGenes/cis_eGene'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
wr = function(x, f) write.table(x, file.path(outdir, f), row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
wr(gtex_shared,   'fastGxC_GTEx_shared.txt')
wr(gtex_specific, 'fastGxC_GTEx_specific.txt')
write.table(gtex_universe, 'data/5a-cis_molGenes/fastGxC_GTEx_universe.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

get_eGenes = function(d, type) intersect(ensg(unique(d[eGene_type == type]$eGene)), clues)

gtex_shared   = get_eGenes(gtex, "FastGxC sh-eGene")
gtex_specific = get_eGenes(gtex, "FastGxC sp-eGene")
gtex_universe = union(gtex_shared, gtex_specific)

outdir = 'data/others/molQTL_clues/cis'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
wr = function(x, f) write.table(x, file.path(outdir, f), row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
wr(gtex_shared,   'fastGxC_GTEx_shared.txt')
wr(gtex_specific, 'fastGxC_GTEx_specific.txt')
write.table(gtex_universe, 'data/others/molQTL_clues/fastGxC_GTEx_universe.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
