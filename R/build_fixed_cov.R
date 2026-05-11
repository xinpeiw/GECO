#' Build a per-individual fixed-covariate matrix
#'
#' Generic version of the original `fixed_cov_all.R`. Produces one row
#' per individual from the cell-level metadata and assembles a numeric
#' matrix suitable to pass to [GeCo()].
#'
#' By default behaves like the legacy script: bins age, dummy-codes sex
#' (Female = 0), population (European = 0), and case status (Healthy = 0),
#' and optionally adds a PC1 column.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param ind_col Individual ID column. Default: "ind_cov".
#' @param covariates A named list describing how to build each column.
#'   The default mimics the legacy script:
#'   `list(age_group = list(col = "Age", type = "age_bins"),
#'         sex = list(col = "Sex", type = "binary", positive = "Male"),
#'         pop = list(col = "pop_cov", type = "binary",
#'                   positive_not = "European"),
#'         cases = list(col = "SLE_status", type = "binary",
#'                     positive_not = "Healthy"))`.
#'   Each entry can also be `type = "numeric"` (use the column as-is) or
#'   `type = "binary"` with either `positive = <value>` (1 if equal) or
#'   `positive_not = <value>` (1 if not equal).
#' @param pc1 Optional data.frame with columns `ind, PC1`
#'   (e.g. from [compute_PC1()]). If supplied, a PC1 column is appended.
#'
#' @return A numeric matrix with `length(unique(data[[ind_col]]))` rows.
#' @export
build_fixed_cov <- function(data,
                            ind_col = "ind_cov",
                            covariates = NULL,
                            pc1 = NULL) {

  dt <- data.table::as.data.table(data)

  if (is.null(covariates)) {
    covariates <- list(
      age_group = list(col = "Age",        type = "age_bins"),
      sex       = list(col = "Sex",        type = "binary", positive = "Male"),
      pop       = list(col = "pop_cov",    type = "binary", positive_not = "European"),
      cases     = list(col = "SLE_status", type = "binary", positive_not = "Healthy")
    )
  }

  # Unique individuals (preserve order of first appearance)
  ord  <- !duplicated(dt[[ind_col]])
  uniq <- dt[ord, ]
  inds <- as.character(uniq[[ind_col]])

  out <- data.frame(ind = inds, stringsAsFactors = FALSE)

  for (name in names(covariates)) {
    spec <- covariates[[name]]
    col  <- spec$col
    if (!is.null(col) && !(col %in% colnames(uniq))) {
      stop("Covariate column '", col, "' not in data.")
    }
    raw <- if (is.null(col)) NA else uniq[[col]]

    val <- switch(spec$type,
      "age_bins" = {
        breaks <- if (!is.null(spec$breaks)) spec$breaks else c(-Inf, 25, seq(30, 90, by = 5), Inf)
        as.numeric(cut(as.numeric(raw), breaks = breaks, right = FALSE,
                       labels = seq_len(length(breaks) - 1)))
      },
      "binary" = {
        if (!is.null(spec$positive)) {
          as.numeric(raw == spec$positive)
        } else if (!is.null(spec$positive_not)) {
          as.numeric(raw != spec$positive_not)
        } else {
          stop("Binary covariate '", name, "' needs `positive` or `positive_not`.")
        }
      },
      "numeric" = as.numeric(raw),
      stop("Unknown covariate type: ", spec$type)
    )
    out[[name]] <- val
  }

  if (!is.null(pc1)) {
    m <- match(out$ind, pc1$ind)
    out$PC1 <- as.numeric(pc1$PC1[m])
  }

  rownames(out) <- out$ind
  out$ind <- NULL
  as.matrix(out)
}
