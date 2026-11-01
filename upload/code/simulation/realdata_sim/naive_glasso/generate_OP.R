args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

setwd(file.path('reps', sprintf("rep_%03d", task_id), 'data'))


P = readRDS('P.rds')
pre_data = data.table::as.data.table(fst::read_fst("data_allcell.fst", columns = c("individual")))
individuals = unique(pre_data$individual)
P = P[which(P$ind %in% individuals),c('CD4 NC','CD4 ET','CD8 NC','CD8 ET','NK','B IN','B Mem')]
p = apply(P, 2, median)
p = p / sum(p)
names(p) = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')


cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
for (c in 1:length(cts)) {
	datExpr = readRDS(paste0('expression_', cts[c], '.rds'))
	d = datExpr*p[cts[c]]
	if (c == 1) { df = d } else { df = df+d }
}

saveRDS(df, 'expression_OP.rds')
