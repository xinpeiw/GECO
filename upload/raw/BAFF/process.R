# Gene sets from MSigDB 2026.1.Hs via the msigdbr R package (26.1.0): hallmark M5911 HALLMARK_INTERFERON_ALPHA_RESPONSE,
#   hallmark M5913 HALLMARK_INTERFERON_GAMMA_RESPONSE, C5 GO:BP GOBP_RESPONSE_TO_INTERFERON_BETA (GO:0035456), and all 50 hallmark sets.
# BAFF (TNFSF13B, ENSG00000102524) neighbours = genes with a non-zero entry in the CLUES cM precision matrix (data/others/precision_clues).
# data/others/BAFF_volcano/WebGestalt_GO_BP_ORA.txt is the result table of a WebGestalt ORA run (hsapiens, geneontology_Biological_Process_noRedundant,
#   ID type ensembl_gene_id, reference list = data/clues/allgenes_CLUES.txt, category size 5-2000, BH FDR, Top 1000) on the 55 neighbours.
suppressMessages({ library(data.table); library(msigdbr) })

U  <- readLines('data/clues/allgenes_CLUES.txt')
om <- readRDS('data/others/precision_clues/precision_cM.rds'); diag(om) <- 0
v  <- om['ENSG00000102524', ]; NB <- names(v)[v != 0]
write.table(NB, 'data/7f_BAFF/BAFF_neighbors.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

H  <- as.data.table(msigdbr(species = 'Homo sapiens', collection = 'H'))
G5 <- as.data.table(msigdbr(species = 'Homo sapiens', collection = 'C5', subcollection = 'GO:BP'))
SETS <- c('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'GOBP_RESPONSE_TO_INTERFERON_BETA', 'HALLMARK_INTERFERON_GAMMA_RESPONSE')

geneset <- function(nm) {
  d <- unique(rbind(H, G5)[gs_name == nm, .(gene_symbol, ensembl_gene)])
  d[, `:=`(in_CLUES_universe = ensembl_gene %in% U, in_BAFF_neighbors = ensembl_gene %in% NB)]
  setorder(d, -in_BAFF_neighbors, -in_CLUES_universe, gene_symbol)
  d
}
fisher_row <- function(nm, d) {
  g <- d[in_CLUES_universe == TRUE, ensembl_gene]; a <- sum(NB %in% g); b <- length(NB) - a; c <- length(g) - a; e <- length(U) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, e), 2), alternative = 'greater')
  data.table(gene_set = nm, set_size_total = nrow(d), set_size_in_universe = length(g), overlap = a,
             expected = length(NB) * length(g) / length(U), odds_ratio = unname(ft$estimate), p_value = ft$p.value,
             overlap_genes = paste(sort(d[in_BAFF_neighbors == TRUE & in_CLUES_universe == TRUE, gene_symbol]), collapse = ';'))
}

gs <- lapply(SETS, geneset); names(gs) <- SETS
for (nm in SETS) fwrite(gs[[nm]], file.path('data/7f_BAFF', paste0('geneset_', nm, '.tsv')), sep = '\t')

res <- rbindlist(lapply(SETS, function(nm) fisher_row(nm, gs[[nm]])))
res[, `:=`(source = 'MSigDB 2026.1.Hs', geneset_size = length(NB), universe_size = length(U),
           fold_enrichment = round(overlap / expected, 1), expected = round(expected, 2), odds_ratio = round(odds_ratio, 2),
           p_adj_BH = signif(p.adjust(p_value, 'BH'), 3), p_value = signif(p_value, 3))]
fwrite(res[, .(gene_set, source, geneset_size, universe_size, set_size_total, set_size_in_universe, overlap, expected,
               fold_enrichment, odds_ratio, p_value, overlap_genes, p_adj_BH)], 'data/7f_BAFF/enrichment_results.tsv', sep = '\t')

scan <- rbindlist(lapply(sort(unique(H$gs_name)), function(nm) fisher_row(nm, geneset(nm))))
scan[, fdr := p.adjust(p_value, 'BH')]
setorder(scan, p_value, gene_set)
dir.create('data/others/BAFF_volcano', showWarnings = FALSE)
fwrite(scan[, .(gene_set, set_in_univ = set_size_in_universe, overlap, expected = round(expected, 2), OR = round(odds_ratio, 2),
                p = p_value, genes = overlap_genes, fdr)], 'data/others/BAFF_volcano/hallmark_scan_all50.tsv', sep = '\t')
