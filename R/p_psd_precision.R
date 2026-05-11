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

#' Fit a sparse precision matrix with glassoFast + eBIC selection
#'
#' Tunes a single scalar multiplier `coe` of a baseline penalty
#' `rho_ori = 1/(z + c)` (where `z = |qnorm(q/2)|` from BH-adjusted
#' p-values) by stepping along a 1-D grid and minimizing eBIC.
#'
#' @param cov Covariance matrix (will be made PSD if not already).
#' @param p   Matching matrix of two-sided p-values.
#' @param n_ind Number of individuals (sample size used in eBIC).
#' @param gamma eBIC gamma parameter (default 0.1).
#' @param c     Numerator offset for rho_ori (default 0.1).
#' @param start_lambda Initial multiplier (default 2).
#' @param step_size    Step size in lambda (default 0.1).
#' @param thr   glasso convergence threshold (default 1e-4).
#' @param maxIt glasso max iterations (default 10000).
#'
#' @return A list with:
#' \describe{
#'   \item{precision}{Fitted precision matrix at the chosen lambda.}
#'   \item{lambda}{Chosen multiplier.}
#'   \item{rho}{Final penalty matrix `rho_ori * lambda`.}
#'   \item{eBIC_path}{`data.frame(lambda, f_value)` of evaluated points.}
#' }
#' @export
fit_precision <- function(cov, p, n_ind,
                          gamma = 0.1, c = 0.1,
                          start_lambda = 2, step_size = 0.1,
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

  f <- function(coe) {
    fit <- glassoFast::glassoFast(S = cov, rho = rho_ori * coe,
                                  thr = thr, maxIt = maxIt)
    E   <- edge_count(fit$wi)
    list(eBIC = eBIC(fit$wi, cov, E), precision = fit$wi)
  }

  evals <- list()
  record <- function(lam) {
    res <- f(lam)
    evals[[length(evals) + 1L]] <<- list(lambda = lam,
                                         f_value = res$eBIC,
                                         precision = res$precision)
    res$eBIC
  }

  current <- start_lambda
  current_val <- record(current)
  left  <- current - step_size
  right <- current + step_size
  left_val  <- record(left)
  right_val <- record(right)

  if (left_val < current_val && left_val <= right_val) {
    dir <- -step_size; best_lam <- left;  best_val <- left_val
  } else if (right_val < current_val) {
    dir <-  step_size; best_lam <- right; best_val <- right_val
  } else {
    best_lam <- current; best_val <- current_val
    path <- do.call(rbind,
              lapply(evals, function(e) data.frame(lambda = e$lambda,
                                                   f_value = e$f_value)))
    chosen <- evals[[which.min(sapply(evals, function(e) e$f_value))]]
    return(list(precision = chosen$precision,
                lambda    = chosen$lambda,
                rho       = rho_ori * chosen$lambda,
                eBIC_path = path))
  }

  repeat {
    nxt <- best_lam + dir
    nxt_val <- record(nxt)
    if (nxt_val >= best_val) break
    best_lam <- nxt; best_val <- nxt_val
  }

  path <- do.call(rbind,
            lapply(evals, function(e) data.frame(lambda = e$lambda,
                                                 f_value = e$f_value)))
  chosen <- evals[[which.min(sapply(evals, function(e) e$f_value))]]
  list(precision = chosen$precision,
       lambda    = chosen$lambda,
       rho       = rho_ori * chosen$lambda,
       eBIC_path = path)
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
