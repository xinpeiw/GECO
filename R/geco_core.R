#' GeCo: per-gene-pair covariance and SE decomposition
#'
#' Given the means tensor `Y` (of shape `[N, G, C]`, with `G = 2` for a
#' single pair) and the within-cell variance tensor `VS`
#' (of shape `[N, C, G, G]`), returns the across-individual covariance
#' (`V_shared`, `G x G`), per-cell-type covariance (`V_GC`, `G*C x G*C`),
#' and their standard errors.
#'
#' @param Y Array `[N, G, C]`: per-(individual, gene, cell-type) means.
#' @param VS Array `[N, C, G, G]`: per-(individual, cell-type) covariance
#'   of gene pair, scaled by 1/n_cells.
#' @param fixed_cov Optional individual-level covariate matrix
#'   (`N x p`) to regress out. `NA` (default) means no regression.
#'
#' @section Errors instead of NA:
#' The estimator needs a complete, finite `Y` and `VS` and at least 2
#' cell types; the covariate regression needs a finite `fixed_cov`
#' with one row per individual. Each of these is checked up front and
#' violated conditions stop with a message naming the offending
#' entries, instead of letting `qr.resid()` fail or NA propagate into
#' the estimates.
#'
#' @return A list with elements:
#' \describe{
#'   \item{V_shared}{`G x G` matrix of across-individual variances/covariances.}
#'   \item{V_GC}{`(G*C) x (G*C)` matrix of within-cell-type variance.}
#'   \item{SE_shared}{Standard errors of `V_shared`.}
#'   \item{SE_GC}{Standard errors of `V_GC`.}
#' }
#' @export
GeCo <- function(Y, VS, fixed_cov = NA) {

  if (length(dim(Y)) != 3L) {
    stop("GeCo(): `Y` must be a 3-D array [N, G, C].", call. = FALSE)
  }
  N <- dim(Y)[1]; G <- dim(Y)[2]; C <- dim(Y)[3]
  if (C < 2L) {
    stop("GeCo(): only ", C, " cell type in `Y`; the shared/specific ",
         "decomposition needs at least 2 cell types (it divides by C*(C-1)).",
         call. = FALSE)
  }
  if (N < 2L) {
    stop("GeCo(): only ", N, " individual in `Y`; need at least 2 to ",
         "estimate across-individual covariance and its SE.", call. = FALSE)
  }
  if (!identical(as.integer(dim(VS)), as.integer(c(N, C, G, G)))) {
    stop("GeCo(): `VS` has dims [", paste(dim(VS), collapse = ", "),
         "] but must be [N, C, G, G] = [", paste(c(N, C, G, G), collapse = ", "),
         "] to match `Y`.", call. = FALSE)
  }
  .check_finite(Y,  "Y",  "GeCo",
                hint = "Some (individual, cell type) mean is missing; GECO needs a complete individual x cell-type grid. See compute_Y().")
  .check_finite(VS, "VS", "GeCo",
                hint = "Some (individual, cell type) within-group variance is missing (group with < 2 cells?). See compute_VS_diag().")

  no_cov <- .no_fixed_cov(fixed_cov)
  if (!no_cov) {
    fixed_cov <- as.matrix(fixed_cov)
    if (nrow(fixed_cov) != N) {
      stop("GeCo(): `fixed_cov` has ", nrow(fixed_cov), " rows but `Y` has ",
           N, " individuals.", call. = FALSE)
    }
    .check_finite(fixed_cov, "fixed_cov", "GeCo",
                  hint = "The covariate regression cannot handle missing covariates. See build_fixed_cov().")
  }

  y <- vectorized(Y)
  if (no_cov) {
    X <- kronecker(matrix(1, ncol = 1, nrow = N), diag(G * C))
  } else {
    X <- cbind(
      kronecker(kronecker(as.matrix(fixed_cov), diag(G)),
                matrix(1, nrow = C, ncol = 1)),
      kronecker(matrix(1, nrow = N, ncol = 1), diag(G * C))
    )
  }
  qrX <- qr(X); y <- qr.resid(qrX, y)
  Y <- aperm(array(y, dim = c(C, G, N)), perm = c(3, 2, 1))

  S <- apply(Y, c(1, 2), sum)

  Vshared  <- matrix(NA_real_, G,    G)
  SEshared <- matrix(NA_real_, G,    G)
  Vgc     <- matrix(NA_real_, G * C, G * C)
  SEgc    <- matrix(NA_real_, G * C, G * C)

  pairs <- which(upper.tri(matrix(TRUE, G, G), diag = TRUE), arr.ind = TRUE)

  for (k in seq_len(nrow(pairs))) {
    j  <- pairs[k, 1]
    jj <- pairs[k, 2]

    a <- (S[, j] * S[, jj] - rowSums(Y[, j, ] * Y[, jj, ])) / (C * (C - 1))
    mean_a <- mean(a)
    se_a   <- sqrt(stats::var(a) / N)

    Vshared[j, jj]  <- Vshared[jj, j]  <- mean_a
    SEshared[j, jj] <- SEshared[jj, j] <- se_a

    E  <- Y[, j, ] * Y[, jj, ] - VS[, , j, jj]
    Gm <- sweep(E, 1, a, `-`)

    m_g <- colMeans(Gm)
    s_g <- sqrt(apply(Gm, 2, stats::var) / N)

    rows_j  <- (0:(C - 1)) * G + j
    cols_jj <- (0:(C - 1)) * G + jj

    Vgc [cbind(rows_j, cols_jj)] <- m_g
    Vgc [cbind(cols_jj, rows_j)] <- m_g
    SEgc[cbind(rows_j, cols_jj)] <- s_g
    SEgc[cbind(cols_jj, rows_j)] <- s_g
  }

  list(V_shared = Vshared, V_GC = Vgc,
       SE_shared = SEshared, SE_GC = SEgc)
}
