library(rhdf5)
library(data.table)
library(fst)

opms = function(g) {

	dt <- as.data.table( read_fst("data_allcell.fst", columns = c("individual", g)) )
	individuals = unique(dt$individual)
	n_ind = length(individuals)
	ind_codes <- setNames(seq_along(individuals), individuals)
	dt[, individual := ind_codes[as.character(individual)] ]

	results <- dt[, .(mu = mean(get(g))), by=.(individual)]
	slice <- matrix(NA_real_, nrow = n_ind)
	slice[cbind(results$individual)] <- results$mu

	m = mean(slice); s = sd(slice)

	return(c(m,s))

}


ctmm = fread("he_wald_free_V_p.txt")
genes = ctmm$gene

ms_allgene = data.frame(matrix(NA, nrow = 11526, ncol = 3))
for (i in 1:11526) {
	g = genes[i]
	ms = opms(g)
	ms_allgene[i,1] = g
	ms_allgene[i,2:3] = ms
}

colnames(ms_allgene) = c('g','OPmean','OPsd')
write.table(ms_allgene, 'OPms_allgene.txt', row.names = FALSE, sep = ' ', quote = FALSE)


