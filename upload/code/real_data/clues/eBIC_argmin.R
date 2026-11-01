args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

ct_id = task_id

cts = c("alpha", "T4", "T8", "NK", "B", "cM", "ncM", "cDC")

group = 'all'
n_ind = 243

rep_dir = paste0('res_',group)

library(glassoFast)

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

ct = cts[ct_id]

cov = readRDS(file.path(rep_dir, paste0('cov_eig_', ct,'.rds')))
p = readRDS(file.path(rep_dir, paste0('p_', ct,'.rds')))
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



find_argmin_lambda <- function(f, start_lambda = 2, step_size = 0.1) {
  
  # Initialize storage for all evaluated points
  all_lambdas <- c()
  all_f_values <- c()
  
  # Helper function to record evaluations
  record_evaluation <- function(lambda_val) {
    f_val <- f(lambda_val)
    all_lambdas <<- c(all_lambdas, lambda_val)
    all_f_values <<- c(all_f_values, f_val)
    return(f_val)
  }
  
  # Start at lambda = 2
  current_lambda <- start_lambda
  current_f_value <- record_evaluation(current_lambda)
  
  # Test both directions to find decreasing direction
  left_lambda <- current_lambda - step_size  # 1.9
  right_lambda <- current_lambda + step_size # 2.1
  
  left_f_value <- record_evaluation(left_lambda)
  right_f_value <- record_evaluation(right_lambda)
  
  # Determine search direction
  if (left_f_value < current_f_value && left_f_value <= right_f_value) {
    # Search in decreasing direction (left)
    search_direction <- -step_size
    best_lambda <- left_lambda
    best_f_value <- left_f_value
  } else if (right_f_value < current_f_value) {
    # Search in increasing direction (right)
    search_direction <- step_size
    best_lambda <- right_lambda
    best_f_value <- right_f_value
  } else {
    # Current point is already optimal among the three tested
    # Create results data frame
    results_df <- data.frame(lambda = all_lambdas, f_value = all_f_values)
    return(list(
      argmin = current_lambda, 
      min_value = current_f_value,
      all_evaluations = results_df,
      search_path = results_df
    ))
  }
  
  # Continue searching in the chosen direction
  while (TRUE) {
    next_lambda <- best_lambda + search_direction
    next_f_value <- record_evaluation(next_lambda)
    
    # If function value increases, stop and return current best
    if (next_f_value >= best_f_value) {
      break
    }
    
    # Update best values and continue
    best_lambda <- next_lambda
    best_f_value <- next_f_value
  }
  
  # Create results data frame
  results_df <- data.frame(lambda = all_lambdas, f_value = all_f_values)
  
  return(list(
    argmin = best_lambda, 
    min_value = best_f_value,
    all_evaluations = results_df,
    search_path = results_df
  ))
}

result <- find_argmin_lambda(f, start_lambda = 2)
res = result$all_evaluations
write.csv(res, file.path(rep_dir, paste0('eBIC_argmin_', ct,'.csv')), row.names = FALSE, quote = TRUE)



