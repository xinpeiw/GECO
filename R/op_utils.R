#' Compute per-gene mean and sd of the individual-mean expression vector
#'
#' For each gene, average expression across cells within each individual
#' (the "OP" representation), then return the mean and sd of those
#' individual-level averages across individuals.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param genes Character vector of gene column names to process.
#' @param ind_col Name of the individual ID column. Default: "ind_cov".
#'
#' @return A `data.frame` with columns g, OPmean, OPsd (one row per gene).
#' @export
compute_OPms <- function(data, genes, ind_col = "ind_cov") {
  if (!(ind_col %in% colnames(data))) {
    stop("ind_col '", ind_col, "' not in data.")
  }
  missing_genes <- setdiff(genes, colnames(data))
  if (length(missing_genes) > 0L) {
    stop(length(missing_genes),
         " gene(s) not in data (first few: ",
         paste(utils::head(missing_genes, 5), collapse = ", "), ").")
  }
  dt <- data.table::as.data.table(data)

  ## One grouped aggregation for all genes at once. `lapply(.SD, mean)`
  ## is dispatched to data.table's GForce fast path (C-level group
  ## aggregation), and grouping by ind_col happens once instead of once
  ## per gene. ~60x faster than a per-gene get(g) loop on realistic
  ## cell counts.
  per_ind <- dt[, lapply(.SD, mean), by = ind_col, .SDcols = genes]
  M <- as.matrix(per_ind[, -1L, with = FALSE])   # n_ind x n_gene

  data.frame(
    g      = genes,
    OPmean = colMeans(M, na.rm = TRUE),
    OPsd   = apply(M, 2, stats::sd, na.rm = TRUE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Compute PC1 of the individual-mean gene-expression matrix
#'
#' Builds a `length(individuals) x length(genes)` matrix of per-individual
#' mean expression, then returns the first PC scores.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param genes Character vector of gene column names to use.
#' @param ind_col Name of the individual ID column. Default: "ind_cov".
#'
#' @return A `data.frame` with columns `ind`, `PC1`.
#' @export
compute_PC1 <- function(data, genes, ind_col = "ind_cov") {
  dt <- data.table::as.data.table(data)
  inds <- as.character(unique(dt[[ind_col]]))

  M <- matrix(NA_real_, nrow = length(inds), ncol = length(genes),
              dimnames = list(inds, genes))
  for (j in seq_along(genes)) {
    sub <- dt[, list(mu = mean(get(genes[j]))), by = ind_col]
    rn  <- as.character(sub[[ind_col]])
    M[rn, j] <- sub$mu
  }

  pca <- stats::prcomp(M, center = TRUE, scale. = TRUE)
  data.frame(ind = inds, PC1 = pca$x[, 1], stringsAsFactors = FALSE)
}
