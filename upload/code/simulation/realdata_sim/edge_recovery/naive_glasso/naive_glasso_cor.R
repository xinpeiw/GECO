args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
task_id <- as.integer(args[2])
rep_dir <- args[3]

seed <- 1e6 + rep_id * 1000 + task_id
set.seed(seed)

i = task_id

cts = c('OP','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
ct = cts[i]

gene_order = readLines('./gene_order.txt')
gene_order_exp = readLines('./gene_order_exp.txt')

omega = readRDS(file.path(rep_dir, 'naive_glasso', paste0('precision_', ct,'.rds')))
rownames(omega) = colnames(omega) = gene_order_exp

omega_reordered <- omega[gene_order, gene_order, drop = FALSE]

cor_matrix = -cov2cor(omega_reordered)
cor_v = cor_matrix[upper.tri(cor_matrix, diag = FALSE)]

write.table(cor_v, file.path(rep_dir, 'naive_glasso', paste0('cor_', ct,'.txt')),
	row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
