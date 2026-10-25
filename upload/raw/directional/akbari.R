rm(list=ls())

library(data.table)
library(GenomicRanges)
library(rtracklayer)

akbari = fread('raw/directional/akbari2026_479_loci_grch37.tsv')
noanno = akbari[which(akbari$gene_caption=='nan'),]
anno = akbari[which(akbari$gene_caption!='nan'),]

genes = anno$gene_caption
hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg = function(x){ e = to_ensg(x); e = e[!is.na(e) & nzchar(e)]; sort(unique(e)) }

genes_ensg = ensg(genes)
geco = fread('data/onek1k/GECO_universe.txt', header = FALSE)$V1
length(genes_ensg); length(genes)
genes_ensg = intersect(genes_ensg, geco)
write.table(genes_ensg,'raw/directional/akbari_anno_gene.txt',row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

snps <- noanno[,c('chrom','pos_grch37','rsid')]
colnames(snps) = c('chr','pos37','rsid')
snps[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
gr37 <- GRanges(
  seqnames = snps$chr,
  ranges = IRanges(start = snps$pos37, end = snps$pos37),
  rsid = snps$rsid
)

chain <- import.chain('raw/directional/hg19ToHg38.over.chain')
gr38_list <- liftOver(gr37, chain)

table(lengths(gr38_list))
keep_unique <- lengths(gr38_list) == 1
gr38 <- unlist(gr38_list[keep_unique])
snps38 <- snps[keep_unique, ]
snps38$chr38 <- as.character(seqnames(gr38))
snps38$pos38 <- start(gr38)
region = snps38[,c('chr38','pos38')]
colnames(region) = c('chrom','start')
region$end = region$start
write.table(region,'raw/directional/akbari_noanno_region.txt',row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
