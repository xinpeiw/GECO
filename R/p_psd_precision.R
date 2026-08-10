#' Two-sided p-value matrix from a covariance and SE matrix
#'
#' @param cov Symmetric covariance matrix.
#' @param se  Matching SE matrix.
#' @return A matrix of two-sided p-values with the same dimnames as
#'   `cov`.
#' @export
compute_p_matrix <- function(cov, se) {
  z <- cov / se
  stats::pnorm(-abs(z)) * 2
}

#' Clip eigenvalues to enforce positive (semi-)definiteness
#'
#' Eigendecomposes a symmetric matrix and replaces eigenvalues below
#' `thr` with `set_value` (default 0), reconstructing the matrix.
#'
#' @param mat Symmetric matrix.
#' @param thr Threshold below which eigenvalues are replaced. Default 0.
#' @param set_value Replacement value (default 0).
#'
#' @return A symmetric matrix of the same dimensions.
#' @export
transform_psd <- function(mat, thr = 0, set_value = 0) {
  eig <- eigen(mat, symmetric = TRUE)
  val <- eig$values
  vec <- eig$vectors
  val[val < thr] <- set_value
  out <- vec %*% diag(val) %*% t(vec)
  dimnames(out) <- dimnames(mat)
  out
}

#' Benjamini–Hochberg q-value matrix from a p-value matrix
#'
#' Adjusts the upper triangle of `P` by Benjamini–Hochberg, mirrors
#' the result back into the lower triangle, and sets the diagonal to 0.
#'
#' @param P p-value matrix.
#' @return q-value matrix of the same dimensions.
#' @export
bh_q_matrix <- function(P) {
  P <- as.matrix(P)
  p <- P
  diag(p) <- NA
  p[upper.tri(p)] <- pmax(p[upper.tri(p)], .Machine$double.xmin)
  qvec <- stats::p.adjust(p[upper.tri(p)], method = "BH")
  Q <- matrix(NA_real_, nrow(p), ncol(p), dimnames = dimnames(p))
  Q[upper.tri(Q)] <- qvec
  Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
  diag(Q) <- 0
  Q
}

#' Count off-diagonal non-zero edges of a precision matrix
#'
#' @param Omega A square precision matrix.
#' @return Number of non-zero off-diagonal entries divided by 2.
#' @export
edge_count <- function(Omega) {
  diag(Omega) <- 0
  sum(abs(Omega) > 0) / 2
}

## Internal: two-phase lambda search (coarse walk + discrete binary refine)
##
## Called by fit_precision(); ported from ebic_fast.R. Given a scalar
## eBIC evaluator `f(lambda) -> numeric`, returns the lambda that
## minimizes `f` on the `fine_step` grid inside `[lower, upper]`,
## along with a history of every lambda actually evaluated. Results
## are cached so no lambda is evaluated twice.
.search_lambda <- function(f, start = 1,
                           fine_step   = 0.1,
                           coarse_step = 1,
                           lower = 0.1, upper = 50,
                           max_iter = 1000) {

  cache <- new.env(parent = emptyenv())
  eval_lambda <- numeric(0); eval_y <- numeric(0)

  eval_f <- function(lambda) {
    if (lambda < lower) return(Inf)
    key <- sprintf("%.10f", lambda)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      y <- f(lambda); assign(key, y, envir = cache)
      eval_lambda <<- c(eval_lambda, lambda)
      eval_y      <<- c(eval_y, y)
    }
    get(key, envir = cache, inherits = FALSE)
  }

  ## Initial direction check on the fine grid.
  start <- round(start / fine_step) * fine_step
  y0 <- eval_f(start)
  yR <- eval_f(start + fine_step)

  if (yR < y0) {
    dir <- 1L
  } else {
    yL <- eval_f(start - fine_step)
    if (yL < y0) {
      dir <- -1L
    } else {
      return(list(lambda = start, y = y0,
                  n_eval = length(eval_lambda),
                  history = data.frame(lambda = eval_lambda, y = eval_y)))
    }
  }

  ## Phase 1: coarse walk to bracket the minimum.
  if (dir == 1L) {
    prev_x <- start; prev_y <- y0; cur_x <- start + coarse_step
    for (iter in seq_len(max_iter)) {
      cur_y <- eval_f(cur_x)
      if (cur_y >= prev_y) { left <- prev_x - coarse_step; right <- cur_x; break }
      prev_x <- cur_x; prev_y <- cur_y; cur_x <- cur_x + coarse_step
      if (cur_x > upper) { left <- prev_x - coarse_step; right <- upper; break }
    }
  } else {
    prev_x <- start; prev_y <- y0; cur_x <- start - coarse_step
    for (iter in seq_len(max_iter)) {
      cur_y <- eval_f(cur_x)
      if (cur_y >= prev_y) { left <- cur_x; right <- prev_x + coarse_step; break }
      prev_x <- cur_x; prev_y <- cur_y; cur_x <- cur_x - coarse_step
      if (cur_x < lower) { left <- lower; right <- prev_x + coarse_step; break }
    }
  }
  if (!exists("left") || !exists("right")) {
    stop("Failed to bracket the minimum within max_iter=", max_iter, ".")
  }

  ## Phase 2: discrete binary search on the fine-step grid inside [left, right].
  k_left  <- as.integer(round(left  / fine_step))
  k_right <- as.integer(round(right / fine_step))
  g <- function(k) eval_f(k * fine_step)

  while (k_left < k_right) {
    mid    <- floor((k_left + k_right) / 2)
    y_mid  <- g(mid)
    y_mid1 <- g(mid + 1L)
    if (y_mid <= y_mid1) k_right <- mid else k_left <- mid + 1L
  }

  lambda_opt <- k_left * fine_step
  y_opt <- eval_f(lambda_opt)

  list(lambda = lambda_opt, y = y_opt,
       bracket = c(left = left, right = right),
       n_eval  = length(eval_lambda),
       history = data.frame(lambda = eval_lambda, y = eval_y))
}


#' Fit a sparse precision matrix with glassoFast + eBIC selection
#'
#' Tunes a single scalar multiplier `coe` of a baseline penalty
#' `rho_ori = 1/(z + c)` (where `z = |qnorm(q/2)|` from BH-adjusted
#' p-values) by minimizing eBIC. The search proceeds in two phases:
#' first a **coarse walk** in steps of `coarse_step` brackets the
#' minimum, then a **discrete binary search** on the `fine_step` grid
#' refines it inside the bracket. Every evaluated lambda is cached so
#' `glassoFast()` never runs twice at the same point.
#'
#' @param cov Covariance matrix (should already be PSD).
#' @param p   Matching matrix of two-sided p-values.
#' @param n_ind Number of individuals (sample size used in eBIC).
#' @param gamma eBIC gamma parameter (default 0.1).
#' @param c     Numerator offset for rho_ori (default 0.1).
#' @param start_lambda Initial multiplier for the search (default 1).
#' @param fine_step   Grid spacing for the discrete refinement
#'   (default 0.1). Also determines the resolution of the returned
#'   optimum.
#' @param coarse_step Step size for the coarse bracketing walk
#'   (default 1). Larger values reach the bracket faster; smaller
#'   values give a tighter starting bracket.
#' @param lower,upper Hard bounds on the lambda multiplier. Default
#'   `0.1` and `50`.
#' @param max_iter Maximum coarse-walk iterations before giving up
#'   (default 1000).
#' @param thr   glasso convergence threshold (default 1e-4).
#' @param maxIt glasso max iterations per call (default 10000).
#'
#' @return A list with:
#' \describe{
#'   \item{precision}{Fitted precision matrix at the chosen lambda.}
#'   \item{lambda}{Chosen multiplier.}
#'   \item{rho}{Final penalty matrix `rho_ori * lambda`.}
#'   \item{eBIC_path}{`data.frame(lambda, y)` of every evaluated lambda.}
#'   \item{n_eval}{Number of distinct lambda values `glassoFast()` was
#'     evaluated at.}
#' }
#' @export
fit_precision <- function(cov, p, n_ind,
                          gamma = 0.1, c = 0.1,
                          start_lambda = 1,
                          fine_step    = 0.1,
                          coarse_step  = 1,
                          lower = 0.1, upper = 50,
                          max_iter = 1000,
                          thr = 1e-4, maxIt = 10000) {

  if (!requireNamespace("glassoFast", quietly = TRUE)) {
    stop("Package 'glassoFast' is required. install.packages('glassoFast').")
  }
  n_gene <- nrow(p)

  q   <- bh_q_matrix(p)
  z   <- abs(stats::qnorm(q / 2))
  rho_ori <- 1 / (z + c)

  eBIC <- function(Theta, S, E) {
    log_det <- as.numeric(determinant(Theta, logarithm = TRUE)$modulus)
    l   <- n_ind * (log_det - sum(S * Theta))
    BIC <- -l + log(n_ind) * E
    BIC + 4 * E * log(n_gene) * gamma
  }

  ## Track the best precision matrix alongside the eBIC values so we
  ## don't have to re-run glassoFast at the winning lambda.
  best <- new.env(parent = emptyenv())
  best$lambda    <- NA_real_
  best$eBIC      <- Inf
  best$precision <- NULL

  f <- function(coe) {
    fit <- glassoFast::glassoFast(S = cov, rho = rho_ori * coe,
                                  thr = thr, maxIt = maxIt)
    E   <- edge_count(fit$wi)
    val <- eBIC(fit$wi, cov, E)
    if (is.finite(val) && val < best$eBIC) {
      best$eBIC      <- val
      best$lambda    <- coe
      best$precision <- fit$wi
    }
    val
  }

  res <- .search_lambda(f,
                        start       = start_lambda,
                        fine_step   = fine_step,
                        coarse_step = coarse_step,
                        lower = lower, upper = upper,
                        max_iter = max_iter)

  chosen_lambda <- res$lambda

  ## The two-phase search may have landed on a lambda whose glasso fit
  ## we already ran (best$precision is that one). If not, run once more.
  if (is.null(best$precision) || !identical(best$lambda, chosen_lambda)) {
    fit <- glassoFast::glassoFast(S = cov, rho = rho_ori * chosen_lambda,
                                  thr = thr, maxIt = maxIt)
    best$precision <- fit$wi
    best$lambda    <- chosen_lambda
  }

  list(precision = best$precision,
       lambda    = chosen_lambda,
       rho       = rho_ori * chosen_lambda,
       eBIC_path = res$history,
       n_eval    = res$n_eval)
}

#' Compute node degrees from a precision matrix
#'
#' @param Omega Square precision matrix (rownames = gene IDs).
#' @return Integer vector of off-diagonal non-zero degrees, named by
#'   the rownames of `Omega`.
#' @export
compute_degree <- function(Omega) {
  M <- abs(Omega) > 0
  diag(M) <- FALSE
  deg <- rowSums(M)
  names(deg) <- rownames(Omega)
  deg
}
