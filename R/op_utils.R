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
  dt <- data.table::as.data.table(data)
  inds <- unique(dt[[ind_col]])
  out <- data.frame(g = character(length(genes)),
                    OPmean = numeric(length(genes)),
                    OPsd   = numeric(length(genes)),
                    stringsAsFactors = FALSE)

  for (i in seq_along(genes)) {
    g <- genes[i]
    if (!(g %in% colnames(dt))) {
      stop("Gene '", g, "' not in data.")
    }
    # mean per individual
    sub <- dt[, list(mu = mean(get(g))), by = ind_col]
    out$g[i]      <- g
    out$OPmean[i] <- mean(sub$mu, na.rm = TRUE)
    out$OPsd[i]   <- stats::sd(sub$mu, na.rm = TRUE)
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
