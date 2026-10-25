rm(list=ls())

library(data.table)

geco_universe = readLines('data/onek1k/GECO_universe.txt')
TF_set = list('Low specificity TF' = intersect(fread('raw/TF_HPA/low_specificity.tsv')$Ensembl, geco_universe),
              'T cell TF' = intersect(fread('raw/TF_HPA/elevated_T-cells_Lineage.tsv')$Ensembl, geco_universe),
              'NK cell TF' = intersect(fread('raw/TF_HPA/elevated_NK-cells_Lineage.tsv')$Ensembl, geco_universe),
              'B cell TF' = intersect(fread('raw/TF_HPA/elevated_B-cells_Lineage.tsv')$Ensembl, geco_universe))
alltfs = unique(unlist(TF_set))

cts = c("alpha", "CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")
for (i in 1:8) {
	M = readRDS(file.path('data/onek1k_precision', paste0(cts[i], '.rds')))
	d = data.frame(gene = rownames(M), degree = rowSums(M != 0) - (diag(M) != 0), ct = cts[i])
	if (i==1) {df = d} else {df = rbind(df, d)}
}
data = df[which(df$gene %in% alltfs),]

g_t = c()
g_b = c()
g_nk = c()
g_shared = c()
for (i in 1:length(alltfs)) {
	g = alltfs[i]
	d = data[which(data$gene == g),]
	m = max(d$degree)
	if (m==0) {
		next
	} else {
		c = d[which(d$degree == m), ]$ct
		if ('alpha' %in% c) {
			g_shared = c(g_shared, g)
		}
		if ('CD4NC' %in% c | 'CD4ET' %in% c | 'CD8NC' %in% c | 'CD8ET' %in% c) {
			g_t = c(g_t, g)
		}
		if ('BIN' %in% c | 'BMem' %in% c) {
			g_b = c(g_b, g)
		}
		if ('NK' %in% c) {
			g_nk = c(g_nk, g)
		}
	}
}

geco_predicted = list('alltfs' = alltfs, 'predicted_Shared' = g_shared, 'predicted_T' = g_t,
                      'predicted_NK' = g_nk, 'predicted_B' = g_b)
for (n in names(geco_predicted)) {
  write.table(geco_predicted[[n]], file.path('data/4a-TF_HPA/GECOpredicted', paste0(n, '.txt')),
              row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
}
