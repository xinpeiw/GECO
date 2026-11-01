rm(list=ls())

library(fst)

data = read_fst('data_allcell.fst') 

# generate OP data (N*G)
ind = as.character(unique(data$individual))
genes = colnames(data)[14:ncol(data)]
col = c('individual', genes)
df = data.frame(matrix(ncol = length(col), nrow = length(ind)))
colnames(df) = col
for (i in 1:length(ind)) {
	if (i %% 100 == 0) { print(i) }
	d = data[which(data$individual == ind[i]),]
	df[i, 2:ncol(df)] = apply(d[,14:ncol(data)], 2, mean)  
	df[i,1] = ind[i]
}
write.table(df, 'OP_Y.txt', row.names = FALSE, sep = '\t', quote = FALSE) # no scale

df_scale = cbind(df$individual, apply(df[,2:ncol(df)],2,scale))
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
write.table(a, 'OP_PC1.txt', row.names = FALSE, sep = '\t', quote = FALSE)





