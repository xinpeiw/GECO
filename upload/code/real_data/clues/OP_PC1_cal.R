rm(list=ls())

library(fst)

group = 'all'


	data = read_fst(paste0('data_allcell_', group, '.fst')) 

	# generate OP data (N*G)
	ind = as.character(unique(data$ind_cov))
	genes = colnames(data)[14:ncol(data)]
	col = c('ind', genes)
	df = data.frame(matrix(ncol = length(col), nrow = length(ind)))
	colnames(df) = col

	for (i in 1:length(ind)) {
		if (i %% 100 == 0) { print(i) }
		d = data[which(data$ind_cov == ind[i]),]
		df[i, 2:ncol(df)] = apply(d[,14:ncol(data)], 2, mean)  
		df[i,1] = ind[i]
	}
	write.table(df, paste0('OP_Y_', group, '.txt'), row.names = FALSE, sep = '\t', quote = FALSE) # no scale

	df_scale = cbind(df$ind, apply(df[,2:ncol(df)],2,scale))
	df_scale = as.data.frame(df_scale)
	names(df_scale)[1] = 'ind'
	write.table(df_scale, paste0('OP_Y_scale_', group, '.txt'), row.names = FALSE, sep = '\t', quote = FALSE) # scale

	pca_result = prcomp(df[, 2:ncol(df)], center = TRUE, scale. = TRUE)
	PC1 = pca_result$x[, 1] 
	head(PC1)
	summary(pca_result)
	a = cbind(df[,1], PC1)
	a = as.data.frame(a)
	names(a) = c('ind','PC1')
	write.table(a, paste0('OP_PC1_', group, '.txt'), row.names = FALSE, sep = '\t', quote = FALSE)





