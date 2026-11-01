args <- commandArgs(trailingOnly = TRUE)
rep_id  <- as.integer(args[1])
task_id <- as.integer(args[2])
rep_dir <- args[3]

seed <- 1e6 + rep_id * 1000 + task_id
set.seed(seed)

i = task_id

library(glassoFast)

cts = c('alpha','c3_CD4NC','c1_CD4ET','c6_CD8NC','c4_CD8ET','c2_NK','c5_BIN','c7_BMem')
ct = ifelse(i==1, strsplit(cts[i],'_')[[1]], strsplit(cts[i],'_')[[1]][2])

n_ind = 928

# helper
eBIC = function(Theta, S, E, gamma) {
  log_det <- as.numeric(determinant(Theta, logarithm = TRUE)$modulus)
  l = n_ind * ( log_det - sum(S*Theta) )
  BIC = -l + log(n_ind) * E
  eBIC = BIC + 4 * E * log(n_gene) * gamma 
  return(list('ll' = l/2, 'BIC' = BIC, 'eBIC' = eBIC))
}

edge = function(Omega) {

  diag(Omega) <- 0
  nonzero_mat = abs(Omega) > 0

  num_nonzero_edges <- sum(nonzero_mat) / 2
  num_nonzero_edges

}

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


cov = readRDS(file.path(rep_dir, 'transform', paste0('cov_eig_', ct,'.rds')))
p = readRDS(file.path(rep_dir, 'merged', paste0('p_', cts[i],'.rds')))
n_gene = nrow(p)
q = bh_q_matrix(p)
z = abs(qnorm(q/2))
c = 0.1
rho_ori = 1/(z+c)

f = function(coe) {

  res = glassoFast(S = cov, rho = rho_ori*coe, thr = 1e-04, maxIt = 10000)

  precision = res$wi
  E = edge(precision)
  a = eBIC(Theta = precision, S = cov, E = E, gamma = 0.1)

  return(a$eBIC)

}

search_lambda <- function(f, start = 1, fine_step = 0.1, coarse_step = 1,lower = 0.1, upper = 50, max_iter = 1000) {
  # ---------- cache ----------
  cache <- new.env(parent = emptyenv())
  eval_lambda <- numeric(0); eval_y <- numeric(0)

  eval_f <- function(lambda) {
    if (lambda < lower) return(Inf)
    key <- sprintf("%.10f", lambda)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      y <- f(lambda); assign(key, y, envir = cache)
      eval_lambda <<- c(eval_lambda, lambda)
      eval_y <<- c(eval_y, y)
    }
    get(key, envir = cache, inherits = FALSE)
  }

  # ---------- initial direction check ----------
  start <- round(start / fine_step) * fine_step

  y0 <- eval_f(start)
  yR <- eval_f(start + fine_step)

  if (yR < y0) { dir <- 1L } else {
    yL <- eval_f(start - fine_step)
    if (yL < y0) {
      dir <- -1L
    } else {
      # start itself is the minimum grid point
      return(list(lambda = start,y = y0,n_eval = length(eval_lambda),history = data.frame(lambda = eval_lambda, y = eval_y)))
    }
  }

  # ---------- coarse search on integer scale ----------
  if (dir == 1L) {
    prev_x <- start
    prev_y <- y0
    cur_x <- start + coarse_step

    for (iter in seq_len(max_iter)) {
      cur_y <- eval_f(cur_x)

      if (cur_y >= prev_y) {
        # minimum is in (prev_x - coarse_step, cur_x)
        left <- prev_x - coarse_step
        right <- cur_x
        break
      }

      prev_x <- cur_x
      prev_y <- cur_y
      cur_x <- cur_x + coarse_step

      if (cur_x > upper) {
        # descended all the way to the upper cap
        left <- prev_x - coarse_step
        right <- upper
        break
      }
    }

  } else {
    prev_x <- start
    prev_y <- y0
    cur_x <- start - coarse_step

    for (iter in seq_len(max_iter)) {
      cur_y <- eval_f(cur_x)

      if (cur_y >= prev_y) {
        # minimum is in (cur_x, prev_x + coarse_step)
        left <- cur_x
        right <- prev_x + coarse_step
        break
      }

      prev_x <- cur_x
      prev_y <- cur_y
      cur_x <- cur_x - coarse_step

      if (cur_x < lower) {
        # descended all the way to the lower cap
        left <- lower
        right <- prev_x + coarse_step
        break
      }
    }
  }

  if (!exists("left") || !exists("right")) {
    stop("Failed to bracket the minimum within max_iter.")
  }

  # ---------- discrete binary search on 0.1 grid inside [left, right] ----------
  k_left  <- as.integer(round(left  / fine_step))
  k_right <- as.integer(round(right / fine_step))

  g <- function(k) eval_f(k * fine_step)

  while (k_left < k_right) {
    mid <- floor((k_left + k_right) / 2)

    y_mid  <- g(mid)
    y_mid1 <- g(mid + 1L)

    if (y_mid <= y_mid1) {
      k_right <- mid
    } else {
      k_left <- mid + 1L
    }
  }

  lambda_opt <- k_left * fine_step
  y_opt <- eval_f(lambda_opt)

  list(lambda = lambda_opt,y = y_opt,bracket = c(left = left, right = right),
    n_eval = length(eval_lambda),history = data.frame(lambda = eval_lambda, y = eval_y))
}

res <- search_lambda(f, start = 1)
df = res$history
write.csv(df, file.path(rep_dir, 'ebic', paste0('eBIC_', ct,'.csv')), row.names = FALSE, quote = TRUE)
