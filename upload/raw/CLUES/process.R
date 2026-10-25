rm(list=ls())

library(data.table)

s7_file <- "raw/CLUES/science.abf1970_table_s7.gz"
geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

s7 <- data.table::fread( cmd = paste("gzip -dc", shQuote(s7_file)),
                         sep = "\t", header = FALSE, skip = 1, showProgress = FALSE )
setnames(s7, c("celltype", "gene", "pvalue", "thresh1", "thresh2"))
clues_universe = unique(s7$gene)

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

universe = intersect(ensg(unique(clues_universe)), geco)
write.table(universe, 'data/5a-cis_molGenes/CLUES_universe.txt',
            row.names = FALSE, col.names = FALSE, sep =' ', quote = FALSE)

universe_clues = intersect(ensg(unique(clues_universe)), clues)
write.table(universe_clues, 'data/others/molQTL_clues/clues_universe.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

s9_file <- "raw/CLUES/science.abf1970_table_s9.gz"
s9 <- data.table::fread( cmd = paste("gzip -dc", shQuote(s9_file)),
                         sep = "\t", header = FALSE, skip = 1, showProgress = FALSE )
setnames(s9, c("celltype", "gene", "pvalue", "thresh1", "thresh2"))
s9 = s9[which(s9$thresh2 < 0.05),]
genes_shared = unique(s9[which(s9$celltype=='AverageTissue'),]$gene)
genes_specific = setdiff(unique(s9[which(s9$celltype!='AverageTissue'),]$gene), genes_shared)

genes_shared_ensg = intersect(ensg(unique(genes_shared)), geco)
genes_specific_ensg = intersect(ensg(unique(genes_specific)), geco)
write.table(genes_shared_ensg, 'data/5a-cis_molGenes/cis_eGene/CLUES_shared.txt',
            row.names = FALSE, col.names = FALSE, sep =' ', quote = FALSE)
write.table(genes_specific_ensg, 'data/5a-cis_molGenes/cis_eGene/CLUES_specific.txt',
            row.names = FALSE, col.names = FALSE, sep =' ', quote = FALSE)

genes_shared_clues = intersect(ensg(unique(genes_shared)), clues)
genes_specific_clues = intersect(ensg(unique(genes_specific)), clues)
write.table(genes_shared_clues, 'data/others/molQTL_clues/cis/CLUES_shared.txt',
            row.names = FALSE, col.names = FALSE, sep =' ', quote = FALSE)
write.table(genes_specific_clues, 'data/others/molQTL_clues/cis/CLUES_specific.txt',
            row.names = FALSE, col.names = FALSE, sep =' ', quote = FALSE)
