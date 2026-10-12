#' Resolve a `fixed_covariates` argument into a numeric matrix (or NA)
#'
#' Helper used by [geco_prepare()]. Accepts any of the forms documented
#' under `geco_prepare(fixed_covariates = ...)` and returns either:
#'   - the resolved per-individual `fixed_cov` matrix, or
#'   - `NA` if `fixed_covariates` is `NULL`/`NA`/`character(0)`, which
#'     tells [GeCo()] to skip the covariate-regression step.
#'
#' @param data,genes,ind_col Passed through from the caller.
#' @param fixed_covariates See [geco_prepare()].
#'
#' @section Errors instead of NA:
#' A precomputed matrix / `data.frame` is validated before use: it
#' must be numeric, contain no NA, and have one row per individual.
#' If it has row names they are matched to the individual IDs in
#' `data` (and reordered to the data's order); without row names the
#' rows are taken positionally and must equal the number of
#' individuals. Covariates built from the catalog are checked for NA
#' by [build_fixed_cov()]. Any failure stops with the offending
#' individuals / columns named.
#'
#' @return A numeric matrix, or `NA`.
#' @export
resolve_fixed_covariates <- function(data, genes,
                                     fixed_covariates = NULL,
                                     ind_col = "ind_cov") {
  if (.no_fixed_cov(fixed_covariates) ||
      (is.character(fixed_covariates) && length(fixed_covariates) == 0L)) {
    return(NA)
  }
  if (is.matrix(fixed_covariates) || is.data.frame(fixed_covariates)) {
    m <- as.matrix(fixed_covariates)
    if (!is.numeric(m)) {
      stop("resolve_fixed_covariates(): the supplied covariate matrix is not ",
           "numeric (storage mode '", storage.mode(m), "'). Encode factors ",
           "as 0/1 columns first.", call. = FALSE)
    }
    inds <- as.character(unique(data[[ind_col]]))
    if (!is.null(rownames(m))) {
      miss <- setdiff(inds, rownames(m))
      if (length(miss) > 0L) {
        stop("resolve_fixed_covariates(): ", length(miss), " individual(s) in ",
             "the data have no row in the supplied covariate matrix: ",
             .preview(miss), ".", call. = FALSE)
      }
      m <- m[inds, , drop = FALSE]
    } else if (nrow(m) != length(inds)) {
      stop("resolve_fixed_covariates(): the supplied covariate matrix has ",
           nrow(m), " rows but the data has ", length(inds), " individuals. ",
           "Give it row names = individual IDs, or one row per individual ",
           "in the order of unique(data[['", ind_col, "']]).", call. = FALSE)
    } else {
      message("[GECO] Covariate matrix has no row names; rows are taken in ",
              "the order of unique(data[['", ind_col, "']]).")
    }
    bad <- which(!is.finite(m), arr.ind = TRUE)
    if (nrow(bad) > 0L) {
      cols <- colnames(m); if (is.null(cols)) cols <- paste0("col", seq_len(ncol(m)))
      stop("resolve_fixed_covariates(): the supplied covariate matrix has ",
           nrow(bad), " NA/NaN/Inf value(s), e.g. individual ",
           .preview(unique(rownames(m)[bad[, 1]]), 5), " in column ",
           .preview(unique(cols[bad[, 2]]), 5), ". The regression step cannot ",
           "handle missing covariates; impute them or drop these individuals.",
           call. = FALSE)
    }
    return(m)
  }

  ## character vector or named list -> hand off to build_fixed_cov.
  needs_pc1 <- FALSE
  if (is.character(fixed_covariates)) {
    needs_pc1 <- "PC1" %in% fixed_covariates
  } else if (is.list(fixed_covariates)) {
    needs_pc1 <- any(vapply(fixed_covariates,
                            function(s) identical(s$type, "pc1"),
                            logical(1)))
  }

  pc1 <- NULL
  if (needs_pc1) {
    message("[GECO] Computing PC1 for fixed_covariates...")
    pc1 <- compute_PC1(data, genes, ind_col = ind_col)
  }

  message("[GECO] Building fixed covariate matrix...")
  build_fixed_cov(data, covariates = fixed_covariates,
                  ind_col = ind_col, pc1 = pc1)
}


#' Built-in catalog of covariate "recipes"
#'
#' Each entry in the catalog tells [build_fixed_cov()] how to turn a
#' raw metadata column into a numeric covariate column. Pass any
#' subset of these names to `build_fixed_cov(covariates = ...)`,
#' [geco_prepare()] via `fixed_covariates`. Use
#' [available_covariates()] to list what's available.
#'
#' @format A named `list`. Each element is a spec list with at least
#'   a `type` field (one of `"age_bins"`, `"binary"`, `"numeric"`,
#'   `"pc1"`).
#'
#' @details
#' Default entries:
#' \describe{
#'   \item{`age_group`}{Age binned into 1..15 using breaks
#'     `c(-Inf, 25, seq(30, 90, 5), Inf)` — source column `"Age"`.}
#'   \item{`sex`}{1 if `"Sex" == "Male"`, else 0 — source column `"Sex"`.}
#'   \item{`PC1`}{The first principal component of the per-individual
#'     mean expression matrix (see [compute_PC1()]). Requires the
#'     `pc1` argument to [build_fixed_cov()].}
#' }
#'
#' Study-specific binary covariates (population, case/control, etc.)
#' are easy to add via a custom spec passed as a named list to
#' [build_fixed_cov()] — see its examples.
#'
#' @export
GECO_COVARIATE_CATALOG <- list(
  age_group = list(col = "Age", type = "age_bins"),
  sex       = list(col = "Sex", type = "binary", positive = "Male"),
  PC1       = list(type = "pc1")
)


#' Names of covariates available in the built-in catalog
#'
#' @return Character vector of catalog keys.
#' @export
available_covariates <- function() names(GECO_COVARIATE_CATALOG)


#' Build a per-individual fixed-covariate matrix
#'
#' Users typically pick a subset of the built-in
#' [GECO_COVARIATE_CATALOG] by name, e.g.
#' `covariates = c("age_group", "sex", "PC1")`. Advanced users can
#' instead pass a named list of full custom specs.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param covariates Either:
#'   - a character vector of names from [GECO_COVARIATE_CATALOG]
#'     (e.g. `c("age_group", "sex", "PC1")`), or
#'   - a named list of fully-specified covariate "recipes" — each
#'     entry should have `col` (the source metadata column) and
#'     `type` (`"age_bins"`, `"binary"`, `"numeric"`, or `"pc1"`).
#'     For `type = "binary"` you also need either `positive` (1 if
#'     equal to that value) or `positive_not` (1 if NOT equal to that
#'     value). For `type = "age_bins"` you can optionally pass
#'     `breaks` (default `c(-Inf, 25, seq(30, 90, 5), Inf)`).
#'   The character-vector form is the recommended way to pick a
#'   subset of the defaults.
#' @param ind_col Individual ID column. Default: "ind_cov".
#' @param pc1 Optional `data.frame(ind, PC1)` from [compute_PC1()].
#'   Required only if `"PC1"` (or any entry of type `"pc1"`) is in
#'   `covariates`.
#'
#' @return A numeric matrix with one row per individual, one column
#'   per requested covariate. Row names are individual IDs.
#'
#' @section Errors instead of NA:
#' Every covariate column is checked after construction. If it is NA
#' for any individual (missing `Age`, an age outside `breaks`, a
#' non-numeric value, an individual absent from `pc1`, ...), the
#' function stops and names the covariate and the individuals, because
#' a single NA would make the covariate regression fail for every gene
#' pair.
#'
#' @examples
#' \dontrun{
#' # Use just age and sex:
#' fcov <- build_fixed_cov(data, covariates = c("age_group", "sex"))
#'
#' # Use age, sex, and PC1:
#' pc1  <- compute_PC1(data, genes)
#' fcov <- build_fixed_cov(data,
#'                         covariates = c("age_group", "sex", "PC1"),
#'                         pc1 = pc1)
#'
#' # Custom covariate (smoker indicator from a "Smoker" column):
#' fcov <- build_fixed_cov(data, covariates = list(
#'   age_group = GECO_COVARIATE_CATALOG$age_group,
#'   smoker = list(col = "Smoker", type = "binary", positive = "Yes")
#' ))
#' }
#' @export
build_fixed_cov <- function(data,
                            covariates = c("age_group", "sex", "PC1"),
                            ind_col = "ind_cov",
                            pc1 = NULL) {

  dt <- data.table::as.data.table(data)

  ## Resolve covariates argument into a named list of specs.
  if (is.character(covariates)) {
    unknown <- setdiff(covariates, names(GECO_COVARIATE_CATALOG))
    if (length(unknown)) {
      stop("Unknown covariate name(s): ", paste(unknown, collapse = ", "),
           ".\n  Available: ", paste(names(GECO_COVARIATE_CATALOG),
                                     collapse = ", "),
           "\n  See ?GECO_COVARIATE_CATALOG to add your own.")
    }
    specs <- GECO_COVARIATE_CATALOG[covariates]
  } else if (is.list(covariates) && !is.null(names(covariates))) {
    specs <- covariates
  } else if (length(covariates) == 0L) {
    stop("`covariates` is empty. Pass at least one covariate, or set ",
         "`fixed_covariates = NULL` upstream to skip the regression.")
  } else {
    stop("`covariates` must be a character vector of catalog names ",
         "or a named list of specs.")
  }

  ## Unique individuals (preserve first-appearance order).
  ord  <- !duplicated(dt[[ind_col]])
  uniq <- dt[ord, ]
  inds <- as.character(uniq[[ind_col]])

  out  <- data.frame(ind = inds, stringsAsFactors = FALSE)

  for (name in names(specs)) {
    spec <- specs[[name]]
    if (identical(spec$type, "pc1")) {
      if (is.null(pc1)) {
        stop("Covariate '", name, "' has type='pc1' but no `pc1` ",
             "argument supplied. Compute it with compute_PC1() and ",
             "pass it as pc1 = <data.frame(ind, PC1)>.")
      }
      m <- match(out$ind, pc1$ind)
      val <- as.numeric(pc1$PC1[m])
    } else {
      col <- spec$col
      if (is.null(col) || !(col %in% colnames(uniq))) {
        stop("Covariate '", name, "' requires column '", col,
             "', which is not in `data`.")
      }
      raw <- uniq[[col]]
      val <- switch(spec$type,
        "age_bins" = {
          breaks <- if (!is.null(spec$breaks)) spec$breaks
                    else c(-Inf, 25, seq(30, 90, by = 5), Inf)
          as.numeric(cut(as.numeric(raw), breaks = breaks, right = FALSE,
                         labels = seq_len(length(breaks) - 1)))
        },
        "binary" = {
          if (!is.null(spec$positive)) {
            as.numeric(raw == spec$positive)
          } else if (!is.null(spec$positive_not)) {
            as.numeric(raw != spec$positive_not)
          } else {
            stop("Binary covariate '", name,
                 "' needs `positive` or `positive_not`.")
          }
        },
        "numeric" = as.numeric(raw),
        stop("Unknown covariate type for '", name, "': ", spec$type)
      )
    }

    bad <- !is.finite(val)
    if (any(bad)) {
      src <- if (identical(spec$type, "pc1")) {
        "individual not found in `pc1$ind`"
      } else {
        paste0("source column '", spec$col, "' is NA, non-numeric",
               if (identical(spec$type, "age_bins")) " or outside `breaks`" else "")
      }
      stop("build_fixed_cov(): covariate '", name, "' is NA for ", sum(bad),
           " of ", length(val), " individual(s): ", .preview(inds[bad]),
           " (", src, "). The covariate regression cannot handle missing ",
           "values; fill them in or drop these individuals before running GECO.",
           call. = FALSE)
    }
    out[[name]] <- val
  }

  rownames(out) <- out$ind
  out$ind <- NULL
  as.matrix(out)
}
