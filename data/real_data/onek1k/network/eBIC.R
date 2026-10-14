args = commandArgs(trailingOnly=TRUE)
task_id = as.numeric(args[1])
i = task_id

eBIC = function(Theta, S, E, gamma) {

	log_det <- as.numeric(determinant(Theta, logarithm = TRUE)$modulus)

	l = 928 * ( log_det - sum(S*Theta) )
	BIC = -l + log(928) * E
	eBIC = BIC + 4 * E * log(11526) * gamma 
	return(list('ll' = l/2, 'BIC' = BIC, 'eBIC' = eBIC))
}

degree = function(Omega) {

	Omega_offdiag <- Omega
	diag(Omega_offdiag) <- 0

	threshold <- 1e-10
	nonzero_mat <- abs(Omega_offdiag) > threshold

	num_nonzero_edges <- sum(nonzero_mat) / 2

	degree_vec <- rowSums(nonzero_mat)
	degree_df <- data.frame(Gene = rownames(Omega), Degree = degree_vec)
	degree_df <- degree_df[order(degree_df$Degree, decreasing = TRUE), ]

	return(list('num_nonzero_edges'=num_nonzero_edges, 'degree_df'=degree_df))
}


cts = c('alpha','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
coe_list = seq(0.1,1,0.1)

res = data.frame('ll' = NA, 'BIC' = NA, 'eBIC' = NA, 'num_edges' = NA, 'lambda' = NA)

for (j in 1:length(coe_list)) {
	theta = readRDS(paste0(cts[i],'_',coe_list[j],'.rds'))
	cov = readRDS(paste0('cov_matrix_scale_eig_', cts[i], '.rds'))
	E = degree(theta)$num_nonzero_edges
	a = eBIC(Theta = theta, S = cov, E = E, gamma = 0.1)

	res[j,'ll'] = a$ll
	res[j,'BIC'] = a$BIC
	res[j,'eBIC'] = a$eBIC
	res[j,'num_edges'] = E
	res[j,'lambda'] = coe_list[j]
	print(res)
}

write.csv(res, paste0('eBIC_',cts[i],'.csv'), row.names = FALSE, quote = TRUE)






