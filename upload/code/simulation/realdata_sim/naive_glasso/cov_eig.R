args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
task_id <- as.integer(args[2])
rep_dir <- args[3]



i = task_id

cts = c('OP','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
exp = readRDS(file.path(rep_dir, 'data', paste0('expression_', cts[i], '.rds')))

cov_matrix = cov(exp)
eig = eigen(cov_matrix, symmetric = TRUE)
val = eig$values
vec = eig$vectors
set_value = 0
thr = 0
val[val < thr] = set_value
cov_new = vec %*% diag(val) %*% t(vec)

outpath = file.path(rep_dir, 'naive_glasso', paste0('cov_eig_', cts[i],'.rds'))
saveRDS(cov_new, file = outpath)
