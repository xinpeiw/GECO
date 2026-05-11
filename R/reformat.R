#' Reformat the long pair-wise GeCo results into per-cell-type matrices
#'
#' Turns the long-format `cov` / `se` `data.frame`s returned by
#' [estimate_geco()] into named `n_gene x n_gene` matrices, one for
#' the across-individual component (`shared`) and one for each cell
#' type. Diagonals are filled from the shared11/22 (and g11/g22)
#' columns.
#'
#' @param est_result Output of [estimate_geco()].
#' @param genes Character vector of gene IDs (length `n_gene`).
#' @param cell_types Character names for the cell types (length
#'   `n_ct`). The first matrix is always returned under the name
#'   `shared`; subsequent ones are named after `cell_types`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cov}{Named list of `n_gene x n_gene` covariance matrices.}
#'   \item{se}{Named list of `n_gene x n_gene` SE matrices.}
#' }
#' @export
reformat_to_matrices <- function(est_result, genes, cell_types = NULL) {
  cov_df <- est_result$cov
  se_df  <- est_result$se
  if (is.null(cell_types)) cell_types <- est_result$cell_types
  if (is.null(cell_types)) {
    n_ct <- (ncol(cov_df) - 5) / 3
    cell_types <- paste0("c", seq_len(n_ct))
  }
  n_ct <- length(cell_types)
  n    <- length(genes)

  build_one <- function(df, off_col, diag_cols) {
    M <- matrix(NA_real_, n, n, dimnames = list(genes, genes))
    idx1 <- match(df$g1, genes)
    idx2 <- match(df$g2, genes)
    keep <- !is.na(idx1) & !is.na(idx2)
    M[cbind(idx1[keep], idx2[keep])] <- df[[off_col]][keep]
    M[cbind(idx2[keep], idx1[keep])] <- df[[off_col]][keep]
    # Diagonal: from g1 (col1) and g2 (col2)
    a1 <- df[, c("g1", diag_cols[1])]; names(a1) <- c("g", "var")
    a2 <- df[, c("g2", diag_cols[2])]; names(a2) <- c("g", "var")
    a1 <- a1[!duplicated(a1$g), ]
    a2 <- a2[!duplicated(a2$g), ]
    d  <- rbind(a1, a2)
    d  <- d[!duplicated(d$g), ]
    m  <- match(genes, d$g)
    diag(M) <- as.numeric(d$var[m])
    M
  }

  cov_out <- list()
  se_out  <- list()

  cov_out[["shared"]] <- build_one(cov_df, "shared12", c("shared11", "shared22"))
  se_out [["shared"]] <- build_one(se_df,  "shared12", c("shared11", "shared22"))

  for (c_idx in seq_len(n_ct)) {
    nm <- cell_types[c_idx]
    cov_out[[nm]] <- build_one(cov_df,
      paste0("g12_c", c_idx),
      c(paste0("g11_c", c_idx), paste0("g22_c", c_idx))
    )
    se_out [[nm]] <- build_one(se_df,
      paste0("g12_c", c_idx),
      c(paste0("g11_c", c_idx), paste0("g22_c", c_idx))
    )
  }

  list(cov = cov_out, se = se_out)
}
