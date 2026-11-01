args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

library(rhdf5)
library(data.table)

load('gene_pairs.Rdata')

n = 260465
start = (task_id - 1) * n + 1
end = task_id * n 

opms = as.data.table(fread('OPms_allgene.txt'))
colnames(opms)[1] = 'gene'
setkey(opms, gene)

for (gp in start:end) {
	g1 = gene_pairs[1,gp]
	g2 = gene_pairs[2,gp]

	pair_local_idx = ((gp - 1) %% n) + 1
  cov12 <- h5read(sprintf("VS_data/VS_cov_%03d.h5", task_id), "VS/cov", index = list(NULL, NULL, pair_local_idx))
  dim(cov12) <- dim(cov12)[1:2]
	# scale
	ms1 = opms[g1]; s1 = ms1$OPsd
	ms2 = opms[g2]; s2 = ms2$OPsd
	cov12 = cov12 / s1 / s2
	a = apply(cov12, 2, mean)

	if (pair_local_idx == 1) { df = a } else { df = rbind(df, a) }
}
colnames(df) = c('CD4ET','NK','CD4NC','CD8ET','BIN','CD8NC','BMem')
ord = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
df = df[, ord]

write.table(df, paste0('VS12_', task_id, '.txt'), row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)