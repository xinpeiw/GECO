#' Compute the within-(individual, cell-type) variance/covariance tensor
#' for a single gene pair
#'
#' Given the raw cell-level expression data, returns a 4-D tensor
#' `VS[ind, ct, k, l]` where the inner 2x2 matrix at every (ind, ct)
#' equals \eqn{\widehat{Var}/n} of (gene1, gene2) within that
#' (individual, cell-type) combination.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param gene1,gene2 Gene column names.
#' @param ind_col Individual ID column. Default: "ind_cov".
#' @param ct_col  Cell-type column. Default: "cg_cov".
#' @param individuals Optional explicit ordering of individuals.
#' @param cell_types  Optional explicit ordering of cell types.
#'
#' @return Array with dims `c(n_ind, n_ct, 2, 2)`. NAs fill empty cells.
#' @export
compute_VS_pair <- function(data, gene1, gene2,
                            ind_col = "ind_cov", ct_col = "cg_cov",
                            individuals = NULL, cell_types = NULL) {

  dt <- data.table::as.data.table(data)
  if (is.null(individuals)) individuals <- as.character(unique(dt[[ind_col]]))
  if (is.null(cell_types))  cell_types  <- as.character(unique(dt[[ct_col]]))
  n_ind <- length(individuals)
  n_ct  <- length(cell_types)

  ind_codes <- stats::setNames(seq_along(individuals), individuals)
  ct_codes  <- stats::setNames(seq_along(cell_types),  cell_types)

  ii <- as.integer(ind_codes[as.character(dt[[ind_col]])])
  cc <- as.integer(ct_codes[as.character(dt[[ct_col]])])
  x  <- dt[[gene1]]; y <- dt[[gene2]]

  dt2 <- data.table::data.table(i = ii, c = cc, x = x, y = y)

  res <- dt2[, list(
              v11 = if (.N > 1) stats::var(x) / .N else NA_real_,
              v22 = if (.N > 1) stats::var(y) / .N else NA_real_,
              v12 = if (.N > 1) stats::cov(x, y) / .N else NA_real_),
              by = list(i, c)]

  out <- array(NA_real_, dim = c(n_ind, n_ct, 2, 2))
  out[cbind(res$i, res$c, 1, 1)] <- res$v11
  out[cbind(res$i, res$c, 2, 2)] <- res$v22
  out[cbind(res$i, res$c, 1, 2)] <- res$v12
  out[cbind(res$i, res$c, 2, 1)] <- res$v12
  out
}
