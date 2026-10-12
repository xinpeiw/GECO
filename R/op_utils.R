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
#' @section Errors instead of NA:
#' `OPsd` is used to standardize every gene before estimation, so it
#' must be finite and non-zero. This function stops and names the
#' offending genes if any gene has NA expression, or if its
#' individual-mean expression is constant across individuals
#' (`OPsd == 0`, e.g. a gene that is zero in every cell) or undefined
#' (fewer than 2 individuals).
#'
#' @return A `data.frame` with columns g, OPmean, OPsd (one row per gene).
#' @export
compute_OPms <- function(data, genes, ind_col = "ind_cov") {
  dt <- data.table::as.data.table(data)
  .check_group_col(dt, ind_col, "individual", "compute_OPms")
  missing_genes <- setdiff(genes, colnames(dt))
  if (length(missing_genes) > 0L) {
    stop("compute_OPms(): ", length(missing_genes),
         " gene(s) not in data: ", .preview(missing_genes), ".",
         call. = FALSE)
  }
  na_genes <- genes[vapply(genes, function(g) anyNA(dt[[g]]), logical(1))]
  if (length(na_genes) > 0L) {
    stop("compute_OPms(): ", length(na_genes), " gene(s) have NA expression ",
         "values: ", .preview(na_genes), ". Expression must be complete.",
         call. = FALSE)
  }
  if (length(unique(dt[[ind_col]])) < 2L) {
    stop("compute_OPms(): only ", length(unique(dt[[ind_col]])),
         " individual in column '", ind_col, "'; the SD of the ",
         "individual-mean expression needs at least 2.", call. = FALSE)
  }

  ## One grouped aggregation for all genes at once. `lapply(.SD, mean)`
  ## is dispatched to data.table's GForce fast path (C-level group
  ## aggregation), and grouping by ind_col happens once instead of once
  ## per gene. ~60x faster than a per-gene get(g) loop on realistic
  ## cell counts.
  per_ind <- dt[, lapply(.SD, mean), by = ind_col, .SDcols = genes]
  M <- as.matrix(per_ind[, -1L, with = FALSE])   # n_ind x n_gene

  out <- data.frame(
    g      = genes,
    OPmean = colMeans(M),
    OPsd   = apply(M, 2, stats::sd),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  bad <- !is.finite(out$OPsd) | out$OPsd == 0
  if (any(bad)) {
    stop("compute_OPms(): ", sum(bad), " gene(s) have zero (or undefined) ",
         "SD of individual-mean expression and cannot be standardized: ",
         .preview(out$g[bad]), ". These genes are constant across ",
         "individuals (typically all-zero). Drop them from `genes`.",
         call. = FALSE)
  }
  out
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
#' @return A `data.frame` with columns `ind`, `PC1`. Stops (naming the
#'   genes) if any gene has NA expression or is constant across
#'   individuals, since `prcomp(scale. = TRUE)` cannot rescale it.
#' @export
compute_PC1 <- function(data, genes, ind_col = "ind_cov") {
  dt <- data.table::as.data.table(data)
  .check_group_col(dt, ind_col, "individual", "compute_PC1")
  missing_genes <- setdiff(genes, colnames(dt))
  if (length(missing_genes) > 0L) {
    stop("compute_PC1(): ", length(missing_genes), " gene(s) not in data: ",
         .preview(missing_genes), ".", call. = FALSE)
  }
  inds <- as.character(unique(dt[[ind_col]]))

  ## Same GForce aggregation as compute_OPms(): one group-by for all genes.
  per_ind <- dt[, lapply(.SD, mean), by = ind_col, .SDcols = genes]
  M <- as.matrix(per_ind[, -1L, with = FALSE])
  rownames(M) <- as.character(per_ind[[ind_col]])
  M <- M[inds, , drop = FALSE]

  .check_finite(M, "individual-mean expression matrix", "compute_PC1",
                hint = "Some gene has NA expression values.")
  sds <- apply(M, 2, stats::sd)
  bad <- !is.finite(sds) | sds == 0
  if (any(bad)) {
    stop("compute_PC1(): ", sum(bad), " gene(s) are constant across ",
         "individuals and cannot be scaled for PCA: ", .preview(genes[bad]),
         ". Drop them from `genes`.", call. = FALSE)
  }

  pca <- stats::prcomp(M, center = TRUE, scale. = TRUE)
  data.frame(ind = inds, PC1 = pca$x[, 1], stringsAsFactors = FALSE)
}
