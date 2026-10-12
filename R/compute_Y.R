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
#' @section Errors instead of NA:
#' GECO's model needs a complete individual x cell-type grid. If any
#' (individual, cell type) combination has no cells, this function
#' stops and lists the affected groups and individuals rather than
#' returning NA (which would make [GeCo()] fail for every gene pair).
#' It also stops if any requested gene column contains NA.
#'
#' @return A list with:
#' \describe{
#'   \item{Y}{Array `[n_ind, n_gene, n_ct]` (no NA; see above).}
#'   \item{individuals}{Character vector of individual IDs (length n_ind).}
#'   \item{cell_types}{Character vector of cell-type IDs (length n_ct).}
#' }
#' @export
compute_Y <- function(data, genes,
                      ind_col = "ind_cov", ct_col = "cg_cov",
                      individuals = NULL, cell_types = NULL) {

  dt <- data.table::as.data.table(data)
  .check_group_col(dt, ind_col, "individual", "compute_Y")
  .check_group_col(dt, ct_col,  "cell-type",  "compute_Y")

  if (is.null(individuals)) individuals <- as.character(unique(dt[[ind_col]]))
  if (is.null(cell_types))  cell_types  <- as.character(unique(dt[[ct_col]]))
  n_ind <- length(individuals)
  n_ct  <- length(cell_types)
  G     <- length(genes)

  missing_genes <- setdiff(genes, colnames(dt))
  if (length(missing_genes) > 0L) {
    stop("compute_Y(): ", length(missing_genes), " gene(s) not in data: ",
         .preview(missing_genes), ".", call. = FALSE)
  }

  # integer indices for fast aggregation (stops on unmapped labels)
  codes  <- .group_codes(dt, ind_col, ct_col, individuals, cell_types,
                         where = "compute_Y")
  counts <- .group_counts(codes, n_ind, n_ct)
  .check_group_sizes(counts, individuals, cell_types,
                     min_cells = 1L, where = "compute_Y")

  dt2 <- data.table::data.table(.__i = codes$i, .__c = codes$c)
  for (g in genes) {
    x <- dt[[g]]
    if (anyNA(x)) {
      stop("compute_Y(): gene '", g, "' has NA expression in ",
           sum(is.na(x)), " cell(s). Expression values must be complete; ",
           "check the .h5ad assay or drop the gene.", call. = FALSE)
    }
    dt2[[g]] <- x
  }

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
