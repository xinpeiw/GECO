# Inputs not included in this repository; download into raw/caGene_pbmc/ (Mu et al., Zenodo record 14652966):
#   https://zenodo.org/records/14652966  files (md5):
#     full_summstats_B_250kb.txt.gz        c332d0e07a125887ecb3aa73d13902e3
#     full_summstats_CD4.T_250kb.txt.gz    f551b6e542cb9944f772727e2b5eaf8d
#     full_summstats_CD8.T_250kb.txt.gz    475180bb7cfa776e7d9581135e984bfb
#     full_summstats_DC_250kb.txt.gz       499d3934bcfe127ab5cd8adbd596ad65
#     full_summstats_other.T_250kb.txt.gz  b8c841c67652a01582341ea7da5fa06b
#     full_summstats_Mono_250kb.txt.gz     3c6211186cd577decafbf7abdb75e3f1
#     full_summstats_NK_250kb.txt.gz       c50f394da7d54b72325d3f32a7f78de5

rm(list=ls())

library(data.table)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

files = list.files('raw/caGene_pbmc', pattern = '^full_summstats_.*_250kb\\.txt\\.gz$', full.names = TRUE)

sig_peaks = unique(unlist(lapply(files, function(f) {
  dt   = fread(cmd = paste('gzip -dc', shQuote(f)),
               select = c('peak', 'PVALUE_RE2'), showProgress = FALSE)
  dt   = dt[!is.na(PVALUE_RE2)]
  lead = dt[, .(p_adj = pmin(min(PVALUE_RE2) * .N, 1)), by = peak]
  lead = lead[is.finite(p_adj)]
  lead[, q := qvalue::qvalue(p_adj)$qvalues]
  lead[q < 0.1, peak]
})))

pk = data.table(peak = sig_peaks)
pk[, c('chrom', 'rng') := tstrsplit(peak, '_', fixed = TRUE)]
pk[, c('start', 'end') := tstrsplit(rng, '-', fixed = TRUE)]
pk[, `:=`(start = as.integer(start), end = as.integer(end))]
pk = pk[!is.na(start) & !is.na(end)]
regs = GRanges(pk$chrom, IRanges(pk$start, pk$end))

genes_gr = genes(EnsDb.Hsapiens.v86)
genes_gr = genes_gr[as.character(seqnames(genes_gr)) %in% c(as.character(1:22), 'X', 'Y')]
gid = names(genes_gr)
genes_gr = GRanges(paste0('chr', as.character(seqnames(genes_gr))), ranges(genes_gr), strand(genes_gr))
names(genes_gr) = gid
tss = ifelse(as.character(strand(genes_gr)) == '-', end(genes_gr), start(genes_gr))
tss_gr = GRanges(seqnames(genes_gr), IRanges(tss, tss), gene_id = names(genes_gr))

d = distanceToNearest(regs, tss_gr, ignore.strand = TRUE)

caGenes_geco = intersect(sort(unique(tss_gr$gene_id[subjectHits(d)])), geco)
out = 'data/5a-cis_molGenes/caGene/pbmc.txt'
write.table(caGenes_geco, out, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

caGenes_clues = intersect(sort(unique(tss_gr$gene_id[subjectHits(d)])), clues)
out = 'data/others/molQTL_clues/cis/ca_pbmc.txt'
write.table(caGenes_clues, out, row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
