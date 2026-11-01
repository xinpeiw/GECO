args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

ct_id = task_id

cts = c("alpha", "T4", "T8", "NK", "B", "cM", "ncM", "cDC")

group = 'all'
thr = 0

rep_dir = paste0('res_',group)

input = file.path(rep_dir, paste0('cov_', cts[ct_id],'.rds'))
cov = readRDS(input)
eig = eigen(cov, symmetric = TRUE)
val = eig$values
vec = eig$vectors
set_value = 0
val[val < thr] = set_value
cov_new = vec %*% diag(val) %*% t(vec)


outpath = file.path(rep_dir, paste0('cov_eig_', cts[ct_id],'.rds'))
saveRDS(cov_new, file = outpath)