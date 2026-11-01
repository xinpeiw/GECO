args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
task_id <- as.integer(args[2])
rep_dir <- args[3]

ct_id = task_id

cts = c('alpha','c1_CD4ET','c2_NK','c3_CD4NC','c4_CD8ET','c5_BIN','c6_CD8NC','c7_BMem')
thr = 0

input = file.path(rep_dir, 'merged', paste0('cov_', cts[ct_id],'.rds'))
cov = readRDS(input)
eig = eigen(cov, symmetric = TRUE)
val = eig$values
vec = eig$vectors
set_value = 0
val[val < thr] = set_value
cov_new = vec %*% diag(val) %*% t(vec)


if (ct_id == 1) {
  ct = strsplit(cts[ct_id],'_')[[1]]
} else {
  ct = strsplit(cts[ct_id],'_')[[1]][2]
}

outpath = file.path(rep_dir, 'transform', paste0('cov_eig_', ct,'.rds'))
saveRDS(cov_new, file = outpath)