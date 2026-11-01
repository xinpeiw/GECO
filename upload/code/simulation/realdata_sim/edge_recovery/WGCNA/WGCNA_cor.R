args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

setwd(file.path('reps', sprintf("rep_%03d", task_id), 'data'))

library(WGCNA)

cts = c('OP','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')

gene_order = readLines('./gene_order.txt')

for (i in 1:length(cts)) {

	datExpr = readRDS(paste0('expression_', cts[i], '.rds'))
	expr_reordered <- datExpr[, gene_order, drop = FALSE]

	cor_mat <- bicor(expr_reordered, use = "pairwise.complete.obs")
	cor_v = cor_mat[upper.tri(cor_mat, diag = FALSE)]


	write.table(cor_v, paste0('cor_',cts[i], '.txt'), 
		row.names = FALSE, col.names = FALSE, sep = ' ', quote = TRUE)
}