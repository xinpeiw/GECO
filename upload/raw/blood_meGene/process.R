rm(list=ls())

library(data.table)
library(readxl)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

xlsx_file <- "raw/blood_meGene/41588_2021_923_MOESM5_ESM.xlsx"
me <- as.data.table(read_excel(xlsx_file, sheet = 1, skip = 1))
me <- me[, .(SNP, CpG = `DNAm site`, P = as.numeric(Discovery_P), is_cis = toupper(as.character(Cis)) == 'TRUE')]

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann = as.data.table(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19),
                    keep.rownames = "CpG")[, .(CpG, gene = UCSC_RefGene_Name)]
cpg2gene = ann[nzchar(gene), .(gene = unique(unlist(strsplit(gene, ';')))), by = CpG]

me_g = merge(me, cpg2gene, by = "CpG", allow.cartesian = TRUE)

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }

me_g[, ensg := to_ensg(gene)]
me_geco = me_g[grepl('^ENSG[0-9]+$', ensg) & ensg %in% geco]
cis_meGenes   = sort(unique(me_geco[is_cis == TRUE]$ensg))
trans_meGenes = sort(unique(me_geco[is_cis == FALSE]$ensg))
out_cis   = 'data/5a-cis_molGenes/cis_meGene/GoDMC.txt'
out_trans = 'data/6a-trans_molGenes/trans_meGene/GoDMC.txt'
write.table(cis_meGenes,   out_cis,   row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
write.table(trans_meGenes, out_trans, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

me_clues = me_g[grepl('^ENSG[0-9]+$', ensg) & ensg %in% clues]
cis_meGenes   = sort(unique(me_clues[is_cis == TRUE]$ensg))
out_cis   = 'data/others/molQTL_clues/cis/GoDMC.txt'
write.table(cis_meGenes,   out_cis,   row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
