args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
rep_dir <- args[2]

cts = c('OP','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')


d_f = function(x) { sum(x != 0) }
r_f = function(x) { sum(abs(x)) }

df_list = list()
genes = readLines('./gene_order_exp.txt')
for (i in 1:length(cts)) {

	v = readRDS(file.path(rep_dir, 'naive_glasso', paste0('precision_', cts[i],'.rds')))
	r = -cov2cor(v)

	diag(v) = 0
	diag(r) = 0
	
	d = apply(v, 2, d_f)
	cor = apply(r, 2, r_f)

	df = cbind(d, cor)
	df = as.data.frame(df)
	df$gene = genes

	df = df[,c('gene','d','cor')]
	colnames(df) = c('gene', 'degree', 'strength')

	df_list[[i]] = df

}

names(df_list) = cts

library(openxlsx)
write.xlsx(df_list, file = file.path(rep_dir, 'naive_glasso', 'degree.xlsx'))



