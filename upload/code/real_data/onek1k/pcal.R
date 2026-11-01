
cts = c('alpha','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')

for (i in 1:length(cts)) {

	cov = readRDS(paste0('cov_', cts[i], '.rds'))
	se = readRDS(paste0('se_', cts[i], '.rds'))

	z = cov/se
	p = pnorm(-abs(z))*2

	saveRDS(p, paste0('p_', cts[i], '.rds'))
}