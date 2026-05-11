#' Compute the within-(individual, cell-type) variance tensor for a
#' list of genes
#'
#' For every gene in `genes`, returns the per-(individual, cell-type)
#' sample variance divided by the number of cells in that group
#' (i.e. the diagonal entries `VS[ind, ct, k, k]` for the pair
#' decomposition). Doing this once for all genes is much cheaper than
#' redoing each gene's per-group variance inside every gene-pair loop.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param genes Character vector of gene column names.
#' @param ind_col Individual ID column. Default: "ind_cov".
#' @param ct_col  Cell-type column. Default: "cg_cov".
#' @param individuals Optional explicit ordering of individuals.
#' @param cell_types  Optional explicit ordering of cell types.
#'
#' @return A 3-D array `[n_ind, n_ct, n_gene]` of `var(gene) / .N`,
#'   with NAs for empty (ind, ct) groups.
#' @export
compute_VS_diag <- function(data, genes,
                            ind_col = "ind_cov", ct_col = "cg_cov",
                            individuals = NULL, cell_types = NULL) {

  dt <- data.table::as.data.table(data)
  if (is.null(individuals)) individuals <- as.character(unique(dt[[ind_col]]))
  if (is.null(cell_types))  cell_types  <- as.character(unique(dt[[ct_col]]))
  n_ind <- length(individuals); n_ct <- length(cell_types); G <- length(genes)

  ind_codes <- stats::setNames(seq_along(individuals), individuals)
  ct_codes  <- stats::setNames(seq_along(cell_types),  cell_types)

  ii <- as.integer(ind_codes[as.character(dt[[ind_col]])])
  cc <- as.integer(ct_codes [as.character(dt[[ct_col]])])

  out <- array(NA_real_, dim = c(n_ind, n_ct, G))

  ## Process gene-by-gene to keep memory bounded.
  for (j in seq_along(genes)) {
    g <- genes[j]
    dt2 <- data.table::data.table(i = ii, c = cc, x = dt[[g]])
    res <- dt2[, list(v = if (.N > 1L)
                            sum((x - mean(x))^2) / ((.N - 1L) * .N)
                          else NA_real_),
               by = list(i, c)]
    slice <- matrix(NA_real_, n_ind, n_ct)
    slice[cbind(res$i, res$c)] <- res$v
    out[, , j] <- slice
  }
  out
}


#' Compute the within-(individual, cell-type) covariance for a gene
#' pair (off-diagonal only)
#'
#' Used inside the chunk loop: each pair only needs the cov-12
#' matrix because the var-11 and var-22 diagonals come from a
#' precomputed [compute_VS_diag()] tensor.
#'
#' @param data,gene1,gene2,ind_col,ct_col,individuals,cell_types See
#'   [compute_VS_pair()].
#'
#' @return Matrix `[n_ind, n_ct]` of `cov(gene1, gene2) / .N`,
#'   with NAs for empty groups.
#' @export
compute_VS_cov_pair <- function(data, gene1, gene2,
                                ind_col = "ind_cov", ct_col = "cg_cov",
                                individuals = NULL, cell_types = NULL) {

  dt <- data.table::as.data.table(data)
  if (is.null(individuals)) individuals <- as.character(unique(dt[[ind_col]]))
  if (is.null(cell_types))  cell_types  <- as.character(unique(dt[[ct_col]]))
  n_ind <- length(individuals); n_ct <- length(cell_types)

  ind_codes <- stats::setNames(seq_along(individuals), individuals)
  ct_codes  <- stats::setNames(seq_along(cell_types),  cell_types)

  ii <- as.integer(ind_codes[as.character(dt[[ind_col]])])
  cc <- as.integer(ct_codes [as.character(dt[[ct_col]])])
  dt2 <- data.table::data.table(i = ii, c = cc,
                                x = dt[[gene1]], y = dt[[gene2]])

  res <- dt2[, list(v = if (.N > 1L)
                          sum((x - mean(x)) * (y - mean(y))) /
                            ((.N - 1L) * .N)
                        else NA_real_),
             by = list(i, c)]
  out <- matrix(NA_real_, n_ind, n_ct)
  out[cbind(res$i, res$c)] <- res$v
  out
}


#' Compute the within-(individual, cell-type) variance/covariance
#' 4-D tensor for a single gene pair
#'
#' Convenience wrapper that combines [compute_VS_diag()] (for the
#' two diagonal entries) and [compute_VS_cov_pair()] (for the
#' off-diagonal) into the 4-D tensor that [GeCo()] expects.
#' Most pipelines should precompute `VS_diag` once with
#' [compute_VS_diag()] and call [compute_VS_cov_pair()] inside the
#' per-pair loop instead — see [estimate_geco()].
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param gene1,gene2 Gene column names.
#' @param ind_col Individual ID column. Default: "ind_cov".
#' @param ct_col  Cell-type column. Default: "cg_cov".
#' @param individuals Optional explicit ordering of individuals.
#' @param cell_types  Optional explicit ordering of cell types.
#'
#' @return Array with dims `c(n_ind, n_ct, 2, 2)`.
#' @export
compute_VS_pair <- function(data, gene1, gene2,
                            ind_col = "ind_cov", ct_col = "cg_cov",
                            individuals = NULL, cell_types = NULL) {

  v1 <- compute_VS_diag(data, gene1,
                        ind_col = ind_col, ct_col = ct_col,
                        individuals = individuals, cell_types = cell_types)
  v2 <- compute_VS_diag(data, gene2,
                        ind_col = ind_col, ct_col = ct_col,
                        individuals = individuals, cell_types = cell_types)
  c12 <- compute_VS_cov_pair(data, gene1, gene2,
                             ind_col = ind_col, ct_col = ct_col,
                             individuals = individuals,
                             cell_types  = cell_types)

  n_ind <- dim(v1)[1]; n_ct <- dim(v1)[2]
  out <- array(NA_real_, dim = c(n_ind, n_ct, 2, 2))
  out[, , 1, 1] <- v1[, , 1]
  out[, , 2, 2] <- v2[, , 1]
  out[, , 1, 2] <- c12
  out[, , 2, 1] <- c12
  out
}
