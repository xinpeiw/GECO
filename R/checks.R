## ---------------------------------------------------------------------
##  Internal validation helpers.
##
##  Package policy: any step that would otherwise produce NA stops
##  immediately with an error that names the cause (which individuals,
##  cell types, genes, or gene pairs are affected) and what to do about
##  it. Downstream steps (GeCo(), glasso) cannot handle NA, and a
##  silent NA would otherwise surface only as an all-NA output or a
##  crash inside glassoFast's C code.
## ---------------------------------------------------------------------

## Format a vector for an error message: first `n` entries + a count.
.preview <- function(x, n = 10L) {
  x <- as.character(x)
  if (length(x) <= n) return(paste(x, collapse = ", "))
  paste0(paste(x[seq_len(n)], collapse = ", "),
         ", ... (", length(x) - n, " more)")
}

## Stop if a grouping column is missing or contains NA.
.check_group_col <- function(dt, col, what, where) {
  if (!(col %in% colnames(dt))) {
    stop(where, "(): ", what, " column '", col, "' is not in the data. ",
         "Available columns: ", .preview(colnames(dt), 15), ".",
         call. = FALSE)
  }
  v <- dt[[col]]
  n_na <- sum(is.na(v))
  if (n_na > 0L) {
    stop(where, "(): ", what, " column '", col, "' is NA for ", n_na,
         " cell(s). Every cell must have a non-missing ", what,
         " label; drop these cells (or fix the metadata) before running GECO.",
         call. = FALSE)
  }
  invisible(TRUE)
}

## Map the data's individual / cell-type labels onto the requested
## orderings. Stops if the data contains labels that are not in
## `individuals` / `cell_types` (the index would be NA).
.group_codes <- function(dt, ind_col, ct_col, individuals, cell_types, where) {
  .check_group_col(dt, ind_col, "individual", where)
  .check_group_col(dt, ct_col,  "cell-type",  where)

  ind_codes <- stats::setNames(seq_along(individuals), individuals)
  ct_codes  <- stats::setNames(seq_along(cell_types),  cell_types)
  ii <- as.integer(ind_codes[as.character(dt[[ind_col]])])
  cc <- as.integer(ct_codes [as.character(dt[[ct_col]])])

  if (anyNA(ii)) {
    bad <- unique(as.character(dt[[ind_col]])[is.na(ii)])
    stop(where, "(): ", length(bad), " individual label(s) in column '",
         ind_col, "' are not in the supplied `individuals` ordering: ",
         .preview(bad), ". Pass `individuals = NULL` to use all ",
         "individuals in the data, or subset the data first.",
         call. = FALSE)
  }
  if (anyNA(cc)) {
    bad <- unique(as.character(dt[[ct_col]])[is.na(cc)])
    stop(where, "(): ", length(bad), " cell-type label(s) in column '",
         ct_col, "' are not in the supplied `cell_types` ordering: ",
         .preview(bad), ". Pass `cell_types = NULL` to use all cell ",
         "types in the data, or subset the data first.",
         call. = FALSE)
  }
  list(i = ii, c = cc)
}

## Per-(individual, cell type) cell counts as an n_ind x n_ct matrix.
.group_counts <- function(codes, n_ind, n_ct) {
  tab <- tabulate(codes$i + n_ind * (codes$c - 1L), nbins = n_ind * n_ct)
  matrix(tab, n_ind, n_ct)
}

## Stop if any (individual, cell type) group has fewer than `min_cells`
## cells. `min_cells = 1` catches missing groups (means would be NA);
## `min_cells = 2` additionally catches singleton groups (within-group
## variance would be NA).
.check_group_sizes <- function(counts, individuals, cell_types,
                               min_cells, where) {
  bad <- which(counts < min_cells, arr.ind = TRUE)
  if (nrow(bad) == 0L) return(invisible(TRUE))

  n_cells  <- counts[bad]
  labels <- sprintf("%s/%s (%s)",
                    individuals[bad[, 1]], cell_types[bad[, 2]],
                    ifelse(n_cells == 0L, "no cells", "1 cell"))
  bad_inds <- unique(individuals[bad[, 1]])
  n_missing <- sum(n_cells == 0L)
  n_single  <- sum(n_cells == 1L)

  what <- c(
    if (n_missing > 0L) paste0(n_missing, " (individual, cell type) group(s) have no cells"),
    if (n_single  > 0L) paste0(n_single,  " (individual, cell type) group(s) have only 1 cell")
  )
  consequence <- c(
    if (n_missing > 0L) "a missing group makes the individual x cell-type grid incomplete, so the per-group mean is NA and GeCo() cannot be fit for ANY gene pair",
    if (n_single  > 0L) "a 1-cell group has no within-group variance (the measurement-noise term), which makes the whole covariance matrix of that cell type NA"
  )
  stop(where, "(): ", paste(what, collapse = " and "), ": ",
       paste(consequence, collapse = "; "),
       ".\n  Affected groups: ", .preview(labels),
       ".\n  Affected individuals (", length(bad_inds), "): ",
       .preview(bad_inds),
       ".\n  GECO requires every individual to have >= 2 cells in every ",
       "cell type. Drop these individuals (or these cell types) before ",
       "running GECO, e.g. keep only individuals present with >= 2 cells ",
       "in all of: ", .preview(cell_types), ".",
       call. = FALSE)
}

## Stop if a numeric array/matrix/vector contains non-finite values.
.check_finite <- function(x, name, where, hint = NULL) {
  bad <- !is.finite(x)
  if (!any(bad)) return(invisible(TRUE))
  n_bad <- sum(bad)
  idx <- which(bad, arr.ind = TRUE)
  idx <- utils::head(idx, 5)
  pos <- if (is.matrix(idx)) {
    dn <- dimnames(x)
    apply(idx, 1, function(r) {
      lab <- vapply(seq_along(r), function(d) {
        if (!is.null(dn) && !is.null(dn[[d]])) dn[[d]][r[d]] else as.character(r[d])
      }, character(1))
      paste0("[", paste(lab, collapse = ", "), "]")
    })
  } else {
    paste0("[", idx, "]")
  }
  stop(where, "(): `", name, "` contains ", n_bad, " NA/NaN/Inf value(s) ",
       "(first few: ", paste(pos, collapse = ", "), ").",
       if (!is.null(hint)) paste0("\n  ", hint) else "",
       call. = FALSE)
}

## Is `fixed_cov` the "no covariates" sentinel (NULL or a single NA)?
.no_fixed_cov <- function(fixed_cov) {
  is.null(fixed_cov) ||
    (length(fixed_cov) == 1L && !is.list(fixed_cov) && is.na(fixed_cov))
}
