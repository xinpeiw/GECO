# Input not included in this repository (934 MB): raw/perturb_seq/primary_geco_intersection.h5 (md5 d1cc65777e32762af26bf971f2e5c536).
#   It holds the log_fc / zscore / adj_p_value layers of GWCD4i.DE_stats.h5ad (genome-wide CD4+ T cell CRISPRi Perturb-seq,
#   Zhu et al. 2025) for the 17,525 primary contrasts whose target is in the GECO universe x the 6,603 intersection genes;
#   row/column metadata are the included primary_geco_intersection_{rows,cols}.tsv.
suppressPackageStartupMessages({
  library(hdf5r)
})

DIR   <- "raw/perturb_seq"
H5IN  <- file.path(DIR, "primary_geco_intersection.h5")
H5OUT <- file.path(DIR, "primary_geco_intersection_normed_effect.h5")
ROWS  <- file.path(DIR, "primary_geco_intersection_rows.tsv")
COLS  <- file.path(DIR, "primary_geco_intersection_cols.tsv")

main <- function() {
  rows <- read.table(ROWS, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cols <- read.table(COLS, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  n_contr <- nrow(rows); n_var <- nrow(cols)
  cat(sprintf("contrasts: %s   intersection genes: %s\n",
              formatC(n_contr, big.mark=","), formatC(n_var, big.mark=",")))

  abs_eff <- abs(as.numeric(rows$self_log_fc))
  invalid_eff <- !is.finite(abs_eff) | abs_eff == 0
  if (any(invalid_eff))
    stop(sprintf("invalid self_log_fc found in %d rows; expected finite, non-zero values",
                 sum(invalid_eff)))
  cat(sprintf("|self_log_fc|  min=%.4f  median=%.4f  max=%.4f\n",
              min(abs_eff, na.rm = TRUE),
              median(abs_eff, na.rm = TRUE),
              max(abs_eff, na.rm = TRUE)))

  col_idx <- match(rows$target_ensg, cols$gene_ensg)
  self_mask_present <- !is.na(col_idx)
  cat(sprintf("self-loop cells: %s (expected = %s)\n",
              formatC(sum(self_mask_present), big.mark=","),
              formatC(n_contr, big.mark=",")))
  stopifnot(sum(self_mask_present) == n_contr)

  cat("loading log_fc ...\n")
  h <- H5File$new(H5IN, mode = "r")
  on.exit(h$close_all(), add = TRUE)
  log_fc <- h[["log_fc"]][seq_len(n_var), seq_len(n_contr)]
  stopifnot(dim(log_fc) == c(n_var, n_contr))

  cat("computing normed_effect ...\n")
  normed <- abs(log_fc) / matrix(abs_eff, nrow = n_var, ncol = n_contr, byrow = TRUE)

  i_self <- col_idx[self_mask_present]
  j_self <- which(self_mask_present)
  normed[cbind(i_self, j_self)] <- NA_real_

  finite <- is.finite(normed)
  cat(sprintf("normed_effect R-dim %d x %d  finite cells: %s / %s (%.2f%%)\n",
              nrow(normed), ncol(normed),
              formatC(sum(finite), big.mark=","), formatC(length(normed), big.mark=","),
              100 * sum(finite) / length(normed)))
  qs <- quantile(normed[finite], c(0.01, 0.25, 0.50, 0.75, 0.90, 0.99, 0.999, 1.0))
  cat("  finite-value quantiles (1, 25, 50, 75, 90, 99, 99.9, max %):\n")
  print(round(qs, 4))

  if (file.exists(H5OUT)) file.remove(H5OUT)
  cat(sprintf("\nwriting %s ...\n", H5OUT))
  fout <- H5File$new(H5OUT, mode = "w")
  f32 <- h5types$H5T_IEEE_F32LE
  space <- H5S$new("simple", dims = c(n_var, n_contr))
  pl <- H5P_DATASET_CREATE$new()
  chunk_contr <- min(256L, n_contr)
  pl$set_chunk(c(n_var, chunk_contr))$set_deflate(4L)$set_shuffle()
  ds <- fout$create_dataset("normed_effect", dtype = f32, space = space,
                             chunk_dims = NULL, gzip_level = NULL,
                             dataset_create_pl = pl)
  ds[seq_len(n_var), seq_len(n_contr)] <- normed
  ds$close(); space$close(); pl$close()

  h5attr(fout, "source_h5")            <- H5IN
  h5attr(fout, "n_contrasts")          <- as.integer(n_contr)
  h5attr(fout, "n_intersection_genes") <- as.integer(n_var)
  h5attr(fout, "definition")           <- "|log_fc(P, m)| / |self_log_fc(P)|"
  h5attr(fout, "self_loops")           <- "self-loop cells (target ENSG == measured ENSG) set to NaN"
  fout$close_all()

  cat(sprintf("wrote %s  (%.1f MB)\n", H5OUT, file.info(H5OUT)$size / 1e6))
}

main()
