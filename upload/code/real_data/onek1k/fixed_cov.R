library(fst)
library(data.table)

data_cov = read_fst("data_allcell.fst", columns = c("individual", "sex", "age")) 
PC1 = fread('OP_PC1.txt')
data_cov = data_cov[!duplicated(data_cov$individual),]

colnames(data_cov)[1] = 'ind'
merged_df = merge(data_cov, PC1, by = "ind", all.x = TRUE)
merged_df = merged_df[match(data_cov$ind, merged_df$ind), ]
merged_df$age_group = cut(merged_df$age, breaks = c(-Inf, 25, seq(30, 90, by = 5), Inf), right = FALSE, labels = rep(1:15))

fixed_cov = merged_df[, c('age_group','sex','PC1')]
fixed_cov$age_group = as.numeric(fixed_cov$age_group)
fixed_cov = as.matrix(fixed_cov)
