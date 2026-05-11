#' Build the Y tensor of per-cell-type-per-individual gene means
#'
#' For each gene, computes the average expression in every
#' (individual, cell-type) combination, returning a 3-D array
#' indexed as `[n_ind, n_gene, n_ct]`.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param genes Character vector of gene column names.
#' @param ind_col Name of the individual ID column. Default: "ind_cov".
#' @param ct_col Name of the cell-type column. Default: "cg_cov".
#' @param individuals Optional explicit ordering of individuals. If `NULL`,
#'   uses the unique values in the data in their natural order.
#' @param cell_types Optional explicit ordering of cell types. If `NULL`,
#'   uses the unique values in the data.
#'
#' @return A list with:
#' \describe{
#'   \item{Y}{Array `[n_ind, n_gene, n_ct]`, with NAs for missing combos.}
#'   \item{individuals}{Character vector of individual IDs (length n_ind).}
#'   \item{cell_types}{Character vector of cell-type IDs (length n_ct).}
#' }
#' @export
compute_Y <- function(data, genes,
                      ind_col = "ind_cov", ct_col = "cg_cov",
                      individuals = NULL, cell_types = NULL) {

  dt <- data.table::as.data.table(data)

  if (is.null(individuals)) individuals <- as.character(unique(dt[[ind_col]]))
  if (is.null(cell_types))  cell_types  <- as.character(unique(dt[[ct_col]]))
  n_ind <- length(individuals)
  n_ct  <- length(cell_types)
  G     <- length(genes)

  ind_codes <- stats::setNames(seq_along(individuals), individuals)
  ct_codes  <- stats::setNames(seq_along(cell_types),  cell_types)

  # integer indices for fast aggregation
  ii <- ind_codes[as.character(dt[[ind_col]])]
  cc <- ct_codes[as.character(dt[[ct_col]])]
  dt2 <- data.table::data.table(.__i = as.integer(ii),
                                .__c = as.integer(cc))
  for (g in genes) dt2[[g]] <- dt[[g]]

  Y <- array(NA_real_, dim = c(n_ind, G, n_ct))

  for (j in seq_along(genes)) {
    g <- genes[j]
    res <- dt2[, list(mu = mean(get(g))), by = c(".__i", ".__c")]
    slice <- matrix(NA_real_, nrow = n_ind, ncol = n_ct)
    slice[cbind(res[[".__i"]], res[[".__c"]])] <- res$mu
    Y[, j, ] <- slice
  }

  list(Y = Y, individuals = individuals, cell_types = cell_types)
}
