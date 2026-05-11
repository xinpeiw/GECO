#' Build a positive-definite matrix from its Cholesky lower-triangular entries
#'
#' @param par numeric vector of length n*(n+1)/2 containing the lower-triangular
#'   entries of L (column-major).
#' @return symmetric positive (semi-) definite matrix V = L %*% t(L).
#' @export
pd_matrix <- function(par) {
  n <- (sqrt(1 + 8 * length(par)) - 1) / 2
  L <- matrix(0, nrow = n, ncol = n)
  L[lower.tri(L, diag = TRUE)] <- par
  V <- L %*% t(L)
  V
}

#' Build a symmetric matrix from its lower-triangular entries
#'
#' @param par numeric vector of length n*(n+1)/2.
#' @return symmetric matrix.
#' @export
sym_matrix <- function(par) {
  n <- (sqrt(1 + 8 * length(par)) - 1) / 2
  V <- matrix(0, n, n)
  V[lower.tri(V, diag = TRUE)] <- par
  V <- V + t(V) - diag(diag(V))
  V
}

#' Vectorize a 3-D tensor with permutation (C, G, N) before flattening
#'
#' @param tensor a 3-D array with dims (N, G, C).
#' @return a numeric vector of length N*G*C.
#' @export
vectorized <- function(tensor) {
  tensor_permuted <- aperm(tensor, perm = c(3, 2, 1))
  as.vector(tensor_permuted)
}

#' Reconstruct a 3-D tensor from a vector produced by `vectorized`
#'
#' @param vec the vector to reshape.
#' @param N,G,C dimensions of the resulting tensor.
#' @return an array with dimensions (N, G, C).
#' @export
tensor_reconstructed <- function(vec, N, G, C) {
  tensor_permuted <- array(vec, dim = c(C, G, N))
  aperm(tensor_permuted, perm = c(3, 2, 1))
}

#' A d x d zero matrix with a single 1 at position (i, i)
#'
#' @param d matrix dimension.
#' @param i index of the diagonal entry to set to 1.
#' @return d x d matrix.
#' @export
only_one <- function(d, i) {
  v <- matrix(0, d, d)
  v[i, i] <- 1
  v
}
