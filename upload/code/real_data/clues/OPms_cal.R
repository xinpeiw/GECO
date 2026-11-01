rm(list=ls())

library(rhdf5)
library(data.table)
library(fst)

group = 'all' # 'SLE' or 'healthy'

opms = function(g) {

	dt <- as.data.table( read_fst(paste0('data_allcell_', group, '.fst'), columns = c("ind_cov", g)) )
	individuals = unique(dt$ind_cov)
	n_ind = length(individuals)
	ind_codes <- setNames(seq_along(individuals), individuals)
	dt[, ind_cov := ind_codes[as.character(ind_cov)] ]

	results <- dt[, .(mu = mean(get(g))), by=.(ind_cov)]
	slice <- matrix(NA_real_, nrow = n_ind)
	slice[cbind(results$ind_cov)] <- results$mu

	m = mean(slice); s = sd(slice)

	return(c(m,s))

}

genes = fread('genelist11087.txt', header = FALSE)$V1

ms_allgene = data.frame(matrix(NA, nrow = length(genes), ncol = 3))
for (i in 1:length(genes)) {
	g = genes[i]
	ms = opms(g)
	ms_allgene[i,1] = g
	ms_allgene[i,2:3] = ms
}

colnames(ms_allgene) = c('g','OPmean','OPsd')
write.table(ms_allgene, paste0('OPms_allgene_', group, '.txt'), row.names = FALSE, sep = ' ', quote = FALSE)



