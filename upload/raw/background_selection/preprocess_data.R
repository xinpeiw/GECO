# Inputs not included in this repository; download into raw/background_selection/:
#   bkgd/chr*.bkgd  B-value map, McVicker et al. 2009 (NCBI36/hg18)
#     http://www.phrap.org/software_dir/mcvicker_dir/bkgd.tar.gz  (13,756,836 bytes; extract the chr*.bkgd files into raw/background_selection/bkgd/)
#   hg38ToHg18.over.chain  UCSC liftOver chain
#     https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg18.over.chain.gz  (507,477 bytes; gunzip)

rm(list = ls())

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)

bkgd_dir   <- 'raw/background_selection/bkgd'
chain_file <- 'raw/background_selection/hg38ToHg18.over.chain'

chroms <- sub('\\.bkgd$', '', list.files(bkgd_dir, pattern = '\\.bkgd$'))

bmap_list <- lapply(chroms, function(ch) {
  d <- fread(file.path(bkgd_dir, paste0(ch, '.bkgd')),
             header = FALSE, col.names = c('b1000', 'len'))
  end   <- cumsum(d$len)
  start <- end - d$len + 1L
  GRanges(seqnames = ch,
          ranges   = IRanges(start = start, end = end),
          B        = d$b1000 / 1000)
})
bmap <- do.call(c, bmap_list)
rm(bmap_list); gc()
message('B map: ', length(bmap), ' segments over ',
        length(unique(as.character(seqnames(bmap)))), ' chromosomes')

genes_hg38 <- genes(EnsDb.Hsapiens.v86)
keep <- as.character(seqnames(genes_hg38)) %in% c(as.character(1:22), 'X', 'Y')
genes_hg38 <- genes_hg38[keep]
gene_ids <- names(genes_hg38)
genes_hg38 <- GRanges(
  seqnames = paste0('chr', as.character(seqnames(genes_hg38))),
  ranges   = ranges(genes_hg38),
  strand   = strand(genes_hg38)
)
names(genes_hg38) <- gene_ids

tss_pos <- ifelse(as.character(strand(genes_hg38)) == '-',
                  end(genes_hg38), start(genes_hg38))
tss_hg38 <- GRanges(seqnames = seqnames(genes_hg38),
                    ranges   = IRanges(tss_pos, tss_pos),
                    strand   = strand(genes_hg38))
names(tss_hg38) <- names(genes_hg38)

chain <- import.chain(chain_file)
lift_body <- liftOver(genes_hg38, chain)
lift_tss  <- liftOver(tss_hg38,  chain)
message('Gene bodies lifted: ',
        sum(elementNROWS(lift_body) > 0), ' / ', length(genes_hg38))

flat <- unlist(lift_body)
strand(flat) <- '*'
hits  <- findOverlaps(flat, bmap, ignore.strand = TRUE)
inter <- pintersect(flat[queryHits(hits)], bmap[subjectHits(hits)],
                    ignore.strand = TRUE)
w  <- width(inter)
gn <- names(flat)[queryHits(hits)]
bv <- bmap$B[subjectHits(hits)]
num <- tapply(w * bv, gn, sum)
den <- tapply(w,      gn, sum)
body_B <- num / den

ftss <- unlist(lift_tss)
strand(ftss) <- '*'
th <- findOverlaps(ftss, bmap, ignore.strand = TRUE)
th <- th[!duplicated(queryHits(th))]
tss_B <- setNames(bmap$B[subjectHits(th)], names(ftss)[queryHits(th)])

out <- data.table(ensg = names(genes_hg38))
out[, bvalue     := body_B[ensg]]
out[, bvalue_tss := tss_B[ensg]]
out <- out[!is.na(bvalue)]

universe = readLines('data/onek1k/GECO_universe.txt')

out_geco = out[which(out$ensg %in% universe),]
write.table(out_geco, 'data/3cd-evolution/B.txt',
            row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
