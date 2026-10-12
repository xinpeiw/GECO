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
#' @section Errors instead of NA:
#' Every output matrix must be complete: each of the
#' `choose(n_gene, 2)` gene pairs must appear exactly once in
#' `est_result` with finite values, and every gene must get a
#' diagonal entry. If not, the function stops and names the missing /
#' non-finite gene pairs. In cluster mode a missing pair almost always
#' means a `chunk_*.rds` file failed or was not merged.
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

  ## ---- completeness of the long table ---------------------------------
  ## (rows whose genes are outside `genes` are ignored, so a sub-matrix
  ## can be extracted by passing a subset of genes)
  pair_rows <- !is.na(cov_df$g1) & !is.na(cov_df$g2) &
               cov_df$g1 %in% genes & cov_df$g2 %in% genes
  n_expected <- choose(n, 2)
  if (sum(pair_rows) != n_expected ||
      anyDuplicated(paste(cov_df$g1[pair_rows], cov_df$g2[pair_rows]))) {
    key_have <- paste(cov_df$g1[pair_rows], cov_df$g2[pair_rows], sep = "\r")
    all_pairs <- utils::combn(genes, 2)
    key_all <- paste(all_pairs[1, ], all_pairs[2, ], sep = "\r")
    key_all_rev <- paste(all_pairs[2, ], all_pairs[1, ], sep = "\r")
    missing <- which(!(key_all %in% key_have) & !(key_all_rev %in% key_have))
    dup <- unique(key_have[duplicated(key_have)])
    stop("reformat_to_matrices(): expected ", n_expected, " gene pairs for ",
         n, " genes but found ", sum(pair_rows), " rows with a gene pair",
         if (length(missing)) paste0("; ", length(missing), " pair(s) missing: ",
           .preview(paste0(all_pairs[1, missing], "-", all_pairs[2, missing]))) else "",
         if (length(dup)) paste0("; ", length(dup), " pair(s) duplicated: ",
           .preview(gsub("\r", "-", dup))) else "",
         ".\n  In cluster mode this means a chunk_*.rds file failed, was ",
         "lost, or was merged twice. Re-run the affected chunks.",
         call. = FALSE)
  }
  val_cols <- setdiff(colnames(cov_df), c("g1", "g2"))
  for (df_name in c("cov", "se")) {
    df <- if (df_name == "cov") cov_df else se_df
    bad <- !is.finite(as.matrix(df[pair_rows, val_cols, drop = FALSE]))
    if (any(bad)) {
      idx <- which(bad, arr.ind = TRUE)
      rows <- which(pair_rows)[unique(idx[, 1])]
      stop("reformat_to_matrices(): the `", df_name, "` table has ", sum(bad),
           " NA/NaN/Inf value(s) in ", length(rows), " gene pair(s), e.g. ",
           .preview(paste0(df$g1[rows], "-", df$g2[rows]), 5),
           " (columns: ", .preview(unique(val_cols[idx[, 2]]), 5), ").",
           call. = FALSE)
    }
  }

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
    .check_finite(M, paste0("matrix '", off_col, "'"), "reformat_to_matrices",
                  hint = "Some gene pair or gene diagonal is missing from est_result.")
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
