library(fst)
library(data.table)


data_cov = read_fst("data_allcell_all.fst", columns = c("ind_cov", "Sex", "Age", "pop_cov", "SLE_status")) 
data_cov = data_cov[!duplicated(data_cov$ind_cov),]
colnames(data_cov)[1] = 'ind'

PC1 = fread('OP_PC1_all.txt')

merged_df = merge(data_cov, PC1, by = "ind", all.x = TRUE)
merged_df = merged_df[match(data_cov$ind, merged_df$ind), ]

merged_df$age_group = cut(as.numeric(merged_df$Age), breaks = c(-Inf, 25, seq(30, 90, by = 5), Inf), right = FALSE, labels = rep(1:15))
merged_df$sex = ifelse(merged_df$Sex == 'Female', 0, 1)
merged_df$pop = ifelse(merged_df$pop_cov == 'European', 0, 1)
merged_df$cases = ifelse(merged_df$SLE_status == 'Healthy', 0, 1)


fixed_cov = merged_df[, c('age_group','sex','pop','PC1','cases')]

fixed_cov$age_group = as.numeric(fixed_cov$age_group)
fixed_cov$sex = as.numeric(fixed_cov$sex)
fixed_cov$PC1 = as.numeric(fixed_cov$PC1)
fixed_cov$pop = as.numeric(fixed_cov$pop)
fixed_cov$cases = as.numeric(fixed_cov$cases)
