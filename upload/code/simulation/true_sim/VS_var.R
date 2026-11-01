library(rhdf5)
library(data.table)


opms = as.data.table(fread('OPms_allgene.txt'))
colnames(opms)[1] = 'gene'
setkey(opms, gene)

genes_all <- opms$gene
gene_to_idx <- setNames(seq_along(genes_all), genes_all)


for (g in genes_all) {
	idx <- unname(gene_to_idx[[g]])
	var <- h5read("VS_data/VS_diag.h5", "VS/var", index = list(NULL, NULL, idx))
	dim(var) <- dim(var)[1:2]

	# scale
	ms = opms[g]; s = ms$OPsd
	var = var / s / s

	a = apply(var, 2, mean)
	if (idx == 1) { df = a } else (df = rbind(df, a))
}
colnames(df) = c('CD4ET','NK','CD4NC','CD8ET','BIN','CD8NC','BMem')
ord = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
df = df[, ord]
write.table(df, 'VS_var_txt', row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)