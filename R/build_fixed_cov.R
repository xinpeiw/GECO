## Internal: validate a user-supplied fixed-covariate matrix.
##
## `fixed_covariates` is either NULL / NA (no covariate regression) or
## a numeric matrix / data.frame with one row per individual. Returns
## the numeric matrix with rows in the order of unique(data[[ind_col]])
## (the order used by compute_Y() and the Y tensor), or NA.
##
## Errors instead of NA: the matrix must be numeric, contain no
## NA/NaN/Inf, and cover every individual in the data. If it has row
## names they are matched to the individual IDs and reordered; without
## row names the rows are taken positionally and must equal the number
## of individuals. Any violation stops with the offending individuals /
## columns named, because a single bad entry would make the covariate
## regression in GeCo() fail for every gene pair.
.resolve_fixed_cov <- function(data, fixed_covariates = NULL,
                               ind_col = "ind_cov") {
  if (.no_fixed_cov(fixed_covariates)) return(NA)

  if (!(is.matrix(fixed_covariates) || is.data.frame(fixed_covariates))) {
    stop("`fixed_covariates` must be NULL (no regression) or a numeric ",
         "matrix / data.frame with one row per individual; got an object of ",
         "class '", paste(class(fixed_covariates), collapse = "/"), "'.",
         call. = FALSE)
  }
  m <- as.matrix(fixed_covariates)
  if (!is.numeric(m)) {
    stop("`fixed_covariates` is not numeric (storage mode '",
         storage.mode(m), "'). Encode factors as 0/1 columns first.",
         call. = FALSE)
  }
  if (ncol(m) == 0L) {
    stop("`fixed_covariates` has no columns. Pass NULL to skip the ",
         "covariate regression.", call. = FALSE)
  }

  inds <- as.character(unique(data[[ind_col]]))
  if (!is.null(rownames(m))) {
    miss <- setdiff(inds, rownames(m))
    if (length(miss) > 0L) {
      stop(length(miss), " individual(s) in the data have no row in ",
           "`fixed_covariates`: ", .preview(miss), ".", call. = FALSE)
    }
    if (anyDuplicated(rownames(m))) {
      stop("`fixed_covariates` has duplicated row names: ",
           .preview(unique(rownames(m)[duplicated(rownames(m))])), ".",
           call. = FALSE)
    }
    m <- m[inds, , drop = FALSE]
  } else if (nrow(m) != length(inds)) {
    stop("`fixed_covariates` has ", nrow(m), " rows but the data has ",
         length(inds), " individuals. Give it row names = individual IDs, ",
         "or one row per individual in the order of unique(data[['",
         ind_col, "']]).", call. = FALSE)
  } else {
    message("[GECO] `fixed_covariates` has no row names; rows are taken in ",
            "the order of unique(data[['", ind_col, "']]).")
    rownames(m) <- inds
  }

  bad <- which(!is.finite(m), arr.ind = TRUE)
  if (nrow(bad) > 0L) {
    cols <- colnames(m)
    if (is.null(cols)) cols <- paste0("col", seq_len(ncol(m)))
    stop("`fixed_covariates` has ", nrow(bad), " NA/NaN/Inf value(s), e.g. ",
         "individual ", .preview(unique(rownames(m)[bad[, 1]]), 5),
         " in column ", .preview(unique(cols[bad[, 2]]), 5),
         ". The covariate regression cannot handle missing values; impute ",
         "them or drop these individuals before running GECO.",
         call. = FALSE)
  }

  const <- apply(m, 2, function(v) length(unique(v)) == 1L)
  if (any(const)) {
    cols <- colnames(m)
    if (is.null(cols)) cols <- paste0("col", seq_len(ncol(m)))
    stop("`fixed_covariates` column(s) ", .preview(cols[const]), " are ",
         "constant across individuals (an intercept is already included ",
         "per gene and cell type; a constant column makes the design ",
         "rank-deficient). Remove them.", call. = FALSE)
  }
  m
}
