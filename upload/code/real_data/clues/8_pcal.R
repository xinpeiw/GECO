args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])


ct_id = task_id

cts = c("alpha", "T4", "T8", "NK", "B", "cM", "ncM", "cDC")
group = 'all'

rep_dir = paste0('res_',group)


cov = readRDS(file.path(rep_dir, paste0('cov_', cts[ct_id],'.rds')))
se = readRDS(file.path(rep_dir, paste0('se_', cts[ct_id],'.rds')))
z = cov/se
p = pnorm(-abs(z)) * 2
saveRDS(p, file.path(rep_dir, paste0('p_', cts[ct_id],'.rds')))
