args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

cts = c('alpha','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
lambdas = c(1, 1.2, 4.1, 2.2, 2.1, 3.4, 6.1, 8)


library(glassoFast)

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

rep_dir = './'

i = task_id
ct = cts[i]
lambda = lambdas[i]

cov = readRDS(file.path(rep_dir, paste0('cov_eig_', ct,'.rds')))
p = readRDS(file.path(rep_dir, paste0('p_', ct,'.rds')))
q = bh_q_matrix(p)
z = abs(qnorm(q/2))
c = 0.1
rho_ori = 1/(z+c)

res = glassoFast(S = cov, rho = rho_ori*lambda, thr = 1e-04, maxIt = 10000)
precision = res$wi
genes = rownames(p)
rownames(precision) = genes
colnames(precision) = genes

saveRDS(precision, file.path(rep_dir, paste0('precision_', ct,'.rds')))

