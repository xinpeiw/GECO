args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

i = task_id

seed <- 1000 * task_id
set.seed(seed)

library(glassoFast)

cts = c('Shared','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
ct = cts[i]

ebic = read.csv(paste0('adj_exp/eBIC_', ct,'.csv'))
idx = which.min(ebic$y)
lambda = ebic$lambda[idx]

cov = readRDS(paste0('adj_exp/naive_cov_mat_', cts[i],'.rds'))
n_gene = nrow(cov)
rho_ori = matrix(1, nrow = n_gene, ncol = n_gene)
diag(rho_ori) = 0

res = glassoFast(S = cov, rho = rho_ori*lambda, thr = 1e-04, maxIt = 10000)
precision = res$wi

saveRDS(precision, paste0('adj_exp/precision_', ct,'.rds'))
