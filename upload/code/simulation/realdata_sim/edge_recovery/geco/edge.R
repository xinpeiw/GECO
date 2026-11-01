
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

cts = c('alpha','CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK', 'BIN', 'BMem')

df = matrix(NA, nrow = choose(11526,2), ncol = 101)
true = readRDS(paste0('../../data/onek1k_precision/', cts[i], '.rds'))
cor_true = -cov2cor(true)
df[,1] = cor_true[upper.tri(cor_true, diag = FALSE)]

for (rep in 1:100) {

	file = file.path('reps',
		sprintf("rep_%03d", rep), 'precision', paste0('precision_',cts[i],'.rds'))

	sim = readRDS(file)
	cor_sim = -cov2cor(sim)

	df[,(rep+1)] = cor_sim[upper.tri(cor_sim, diag = FALSE)]
}

df = as.data.frame(df)
colnames(df) = c('true',sprintf("rep_%03d", 1:100))

saveRDS(df, paste0('edge/geco/edge_cor_',
	cts[i], '.rds'))

