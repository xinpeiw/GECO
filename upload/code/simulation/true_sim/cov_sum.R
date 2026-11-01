cts = c('alpha','c3_CD4NC','c1_CD4ET','c6_CD8NC','c4_CD8ET','c2_NK','c5_BIN','c7_BMem')

for (c in 1:length(cts)) {
	true = readRDS(paste0('cov_matrix_scale_', cts[c], '.rds'))
	true_cov = true[upper.tri(true, diag = FALSE)]
	true_var = diag(true)
	print(cts[c])
	print(summary(true_cov))
	print(summary(true_var))
}