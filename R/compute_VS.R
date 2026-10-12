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
#' @section Errors instead of NA:
#' The within-group variance is undefined for groups with fewer than 2
#' cells. Rather than returning NA (which would make the whole
#' covariance matrix of that cell type NA downstream), this function
#' stops and lists every (individual, cell type) group with fewer than
#' 2 cells, plus the affected individuals.
#'
#' @return A 3-D array `[n_ind, n_ct, n_gene]` of `var(gene) / .N`
#'   (no NA; see above).
#' @export
compute_VS_diag <- function(data, genes,
                            ind_col = "ind_cov", ct_col = "cg_cov",
                            individuals = NULL, cell_types = NULL) {

  dt <- data.table::as.data.table(data)
  .check_group_col(dt, ind_col, "individual", "compute_VS_diag")
  .check_group_col(dt, ct_col,  "cell-type",  "compute_VS_diag")
  if (is.null(individuals)) individuals <- as.character(unique(dt[[ind_col]]))
  if (is.null(cell_types))  cell_types  <- as.character(unique(dt[[ct_col]]))
  n_ind <- length(individuals); n_ct <- length(cell_types); G <- length(genes)

  missing_genes <- setdiff(genes, colnames(dt))
  if (length(missing_genes) > 0L) {
    stop("compute_VS_diag(): ", length(missing_genes), " gene(s) not in data: ",
         .preview(missing_genes), ".", call. = FALSE)
  }

  codes  <- .group_codes(dt, ind_col, ct_col, individuals, cell_types,
                         where = "compute_VS_diag")
  counts <- .group_counts(codes, n_ind, n_ct)
  .check_group_sizes(counts, individuals, cell_types,
                     min_cells = 2L, where = "compute_VS_diag")
  ii <- codes$i; cc <- codes$c

  out <- array(NA_real_, dim = c(n_ind, n_ct, G))

  ## Process gene-by-gene to keep memory bounded.
  for (j in seq_along(genes)) {
    g <- genes[j]
    x_all <- dt[[g]]
    if (anyNA(x_all)) {
      stop("compute_VS_diag(): gene '", g, "' has NA expression in ",
           sum(is.na(x_all)), " cell(s). Expression values must be complete.",
           call. = FALSE)
    }
    dt2 <- data.table::data.table(i = ii, c = cc, x = x_all)
    res <- dt2[, list(v = sum((x - mean(x))^2) / ((.N - 1L) * .N)),
               by = list(i, c)]
    slice <- matrix(NA_real_, n_ind, n_ct)
    slice[cbind(res$i, res$c)] <- res$v
    out[, , j] <- slice
  }
  out
}


## ---------------------------------------------------------------------
##  Fast per-pair within-group covariance.
##
##  The per-pair loop (estimate_geco(), geco_chunk()) is dominated by
##  computing cov(g1, g2) within every (individual, cell type) group.
##  Doing that with a data.table group-by per pair costs ~100 ms at
##  700k cells, mostly in re-deriving the grouping and in the grouped
##  mean-centering. Instead we:
##    * derive the integer group code of every cell ONCE per chunk
##      (.group_setup), and
##    * per pair do a single C-level rowsum() of the product of the two
##      (globally centred) gene vectors, then correct with the per-group
##      means, which are already available as the Y tensor.
##  Same numbers (differences ~1e-17), ~10x faster.
## ---------------------------------------------------------------------

## Internal: one-time grouping for a data.table of cells.
## Returns list(grp, ng, n_ind, n_ct) where grp[cell] in 1..n_ind*n_ct
## with individual index fastest (so matrix(., n_ind, n_ct) works) and
## ng = cells per group. Stops if any group has < 2 cells.
.group_setup <- function(dt, ind_col, ct_col, individuals, cell_types, where) {
  n_ind <- length(individuals); n_ct <- length(cell_types)
  codes  <- .group_codes(dt, ind_col, ct_col, individuals, cell_types, where)
  counts <- .group_counts(codes, n_ind, n_ct)
  .check_group_sizes(counts, individuals, cell_types,
                     min_cells = 2L, where = where)
  list(grp   = codes$i + n_ind * (codes$c - 1L),
       ng    = as.vector(counts),
       n_ind = n_ind, n_ct = n_ct)
}

## Internal: per-group sum of a vector, as a length n_ind*n_ct vector.
## All groups are guaranteed non-empty by .group_setup(), so rowsum()'s
## sorted output lines up with group codes 1..n_ind*n_ct.
.group_sum <- function(x, setup) {
  as.vector(rowsum(x, setup$grp, reorder = TRUE))
}

## Internal: within-group cov(x, y) / n for one gene pair, given the
## per-group means xbar, ybar ([n_ind, n_ct] matrices, e.g. slices of
## the Y tensor computed from the SAME cells). x and y are centred by
## their global means first so the Sxy - n*xbar*ybar correction does
## not suffer cancellation when expression means are large.
.group_cov_pair <- function(x, y, xbar, ybar, setup) {
  a <- mean(x); b <- mean(y)
  Sxy <- .group_sum((x - a) * (y - b), setup)
  ng  <- setup$ng
  v   <- (Sxy - ng * (as.vector(xbar) - a) * (as.vector(ybar) - b)) /
         ((ng - 1) * ng)
  matrix(v, setup$n_ind, setup$n_ct)
}


#' Compute the within-(individual, cell-type) covariance for a gene
#' pair (off-diagonal only)
#'
#' Standalone version of the per-pair step used inside
#' [estimate_geco()] / [geco_chunk()]: each pair only needs the
#' cov-12 matrix because the var-11 and var-22 diagonals come from a
#' precomputed [compute_VS_diag()] tensor. Implemented with a single
#' `rowsum()` over the cells rather than a grouped data.table
#' aggregation (about 10x faster at 1e6 cells, identical numbers).
#'
#' @param data,gene1,gene2,ind_col,ct_col,individuals,cell_types See
#'   [compute_VS_pair()].
#'
#' @return Matrix `[n_ind, n_ct]` of `cov(gene1, gene2) / .N`. Stops
#'   (instead of returning NA) if any (individual, cell type) group
#'   has fewer than 2 cells or either gene has NA expression.
#' @export
compute_VS_cov_pair <- function(data, gene1, gene2,
                                ind_col = "ind_cov", ct_col = "cg_cov",
                                individuals = NULL, cell_types = NULL) {

  dt <- data.table::as.data.table(data)
  if (is.null(individuals)) individuals <- as.character(unique(dt[[ind_col]]))
  if (is.null(cell_types))  cell_types  <- as.character(unique(dt[[ct_col]]))

  for (g in c(gene1, gene2)) {
    if (!(g %in% colnames(dt))) {
      stop("compute_VS_cov_pair(): gene '", g, "' not in data.", call. = FALSE)
    }
    if (anyNA(dt[[g]])) {
      stop("compute_VS_cov_pair(): gene '", g, "' has NA expression in ",
           sum(is.na(dt[[g]])), " cell(s).", call. = FALSE)
    }
  }

  setup <- .group_setup(dt, ind_col, ct_col, individuals, cell_types,
                        where = "compute_VS_cov_pair")
  x <- as.numeric(dt[[gene1]]); y <- as.numeric(dt[[gene2]])
  xbar <- .group_sum(x, setup) / setup$ng
  ybar <- .group_sum(y, setup) / setup$ng
  .group_cov_pair(x, y, xbar, ybar, setup)
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
