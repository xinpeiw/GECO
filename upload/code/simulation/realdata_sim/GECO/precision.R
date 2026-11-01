args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
task_id <- as.integer(args[2])
rep_dir <- args[3]

seed <- 1e6 + rep_id * 1000 + task_id
set.seed(seed)

i = task_id

library(glassoFast)

cts = c('alpha','c3_CD4NC','c1_CD4ET','c6_CD8NC','c4_CD8ET','c2_NK','c5_BIN','c7_BMem')

bh_q_matrix <- function(P) {
  P <- as.matrix(P)
  p <- P
  diag(p) <- NA
  p[upper.tri(p)] <- pmax(p[upper.tri(p)], .Machine$double.xmin) # protect against 0
  qvec <- p.adjust(p[upper.tri(p)], method = "BH") # BH on upper triangle only
  Q <- matrix(NA_real_, nrow(p), ncol(p), dimnames = dimnames(p))
  Q[upper.tri(Q)] <- qvec
  Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
  diag(Q) <- 0
  Q
}

ct = ifelse(i==1, strsplit(cts[i],'_')[[1]], strsplit(cts[i],'_')[[1]][2])

ebic = read.csv(file.path(rep_dir, 'ebic', paste0('eBIC_', ct,'.csv')))
idx = which.min(ebic$y)
lambda = ebic$lambda[idx]

cov = readRDS(file.path(rep_dir, 'transform', paste0('cov_eig_', ct,'.rds')))
p = readRDS(file.path(rep_dir, 'merged', paste0('p_', cts[i],'.rds')))
q = bh_q_matrix(p)
z = abs(qnorm(q/2))
c = 0.1
rho_ori = 1/(z+c)

res = glassoFast(S = cov, rho = rho_ori*lambda, thr = 1e-04, maxIt = 10000)
precision = res$wi
genes = rownames(p)
rownames(precision) = genes
colnames(precision) = genes

saveRDS(precision, file.path(rep_dir, 'precision', paste0('precision_', ct,'.rds')))

