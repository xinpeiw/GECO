rm(list=ls())

library(data.table)
library(readxl)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

xlsx_file <- "raw/MemoryT/41586_2022_4713_MOESM3_ESM.xlsx"
st20 <- as.data.table(read_excel(xlsx_file, sheet = "20 PME CV1-7 interact",
                                 col_types = c(rep("guess", 18), "text")))
st20[, `Model LRT p value` := as.numeric(`Model LRT p value`)]
universe_study = unique(as.data.table(read_excel(xlsx_file, sheet = "1 Pseudobulk memory T eQTLs"))$Gene)
st20[, lrt_q := qvalue::qvalue(`Model LRT p value`)$qvalues]
st20 = st20[lrt_q < 0.05]
genes = unique(st20$Gene)

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

genes_ensg = intersect(ensg(unique(genes)), geco)
write.table(genes_ensg, 'data/5a-cis_molGenes/cis_eGene/MemoryT.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
universe = intersect(ensg(unique(universe_study)), geco)
write.table(universe, 'data/5a-cis_molGenes/MemoryT_universe.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

genes_ensg_clues = intersect(ensg(unique(genes)), clues)
write.table(genes_ensg_clues, 'data/others/molQTL_clues/cis/MemoryT.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
universe_clues = intersect(ensg(unique(universe_study)), clues)
write.table(universe_clues, 'data/others/molQTL_clues/MemoryT_universe.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
