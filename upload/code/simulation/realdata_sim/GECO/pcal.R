args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
task_id <- as.integer(args[2])
rep_dir <- args[3]

ct_id = task_id

cts = c('alpha','c3_CD4NC','c1_CD4ET','c6_CD8NC','c4_CD8ET','c2_NK','c5_BIN','c7_BMem')


cov = readRDS(file.path(rep_dir, 'merged', paste0('cov_', cts[ct_id],'.rds')))
se = readRDS(file.path(rep_dir, 'merged', paste0('se_', cts[ct_id],'.rds')))
z = cov/se
p = pnorm(-abs(z)) * 2
saveRDS(p, file.path(rep_dir, 'merged', paste0('p_', cts[ct_id],'.rds')))