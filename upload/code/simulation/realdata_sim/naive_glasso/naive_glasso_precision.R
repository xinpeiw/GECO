args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
task_id <- as.integer(args[2])
rep_dir <- args[3]

seed <- 1e6 + rep_id * 1000 + task_id
set.seed(seed)

i = task_id

library(glassoFast)

cts = c('OP','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
ct = cts[i]

ebic = read.csv(file.path(rep_dir, 'naive_glasso', paste0('eBIC_', ct,'.csv')))
idx = which.min(ebic$y)
lambda = ebic$lambda[idx]

cov = readRDS(file.path(rep_dir, 'naive_glasso', paste0('cov_eig_', ct,'.rds')))
n_gene = nrow(cov)
rho_ori = matrix(1, nrow = n_gene, ncol = n_gene)
diag(rho_ori) = 0

res = glassoFast(S = cov, rho = rho_ori*lambda, thr = 1e-04, maxIt = 10000)
precision = res$wi

saveRDS(precision, file.path(rep_dir, 'naive_glasso', paste0('precision_', ct,'.rds')))

