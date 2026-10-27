# Input not included in this repository (934 MB): raw/perturb_seq/primary_geco_intersection.h5 (md5 d1cc65777e32762af26bf971f2e5c536).
#   It holds the log_fc / zscore / adj_p_value layers of GWCD4i.DE_stats.h5ad (genome-wide CD4+ T cell CRISPRi Perturb-seq,
#   Zhu et al. 2025) for the 17,525 primary contrasts whose target is in the GECO universe x the 6,603 intersection genes;
#   row/column metadata are the included primary_geco_intersection_{rows,cols}.tsv.
suppressPackageStartupMessages({
  library(hdf5r)
})

DIR     <- "raw/perturb_seq"
H5_DAT  <- file.path(DIR, "primary_geco_intersection.h5")
H5_NORM <- file.path(DIR, "primary_geco_intersection_normed_effect.h5")
ROWS    <- file.path(DIR, "primary_geco_intersection_rows.tsv")
COLS    <- file.path(DIR, "primary_geco_intersection_cols.tsv")
UNION   <- "data/others/pertseq_union_network_matrix.rds"
OUT_RDS <- "data/6b-edge/perturb_seq/pertseq_normed_effect_network_matrix.rds"
THRESH  <- 0.10

main <- function() {
  rows <- read.table(ROWS, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cols <- read.table(COLS, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  n_contr <- nrow(rows); n_var <- nrow(cols)
  cat(sprintf("contrasts: %s   intersection genes: %s\n",
              formatC(n_contr, big.mark = ","), formatC(n_var, big.mark = ",")))

  target_col_idx <- match(rows$target_ensg, cols$gene_ensg)
  if (anyNA(target_col_idx))
    stop("at least one contrast target is missing from the intersection columns")

  cat("loading layers ...\n")
  hd <- H5File$new(H5_DAT, mode = "r")
  hn <- H5File$new(H5_NORM, mode = "r")
  on.exit({ hd$close_all(); hn$close_all() }, add = TRUE)
  adj_p  <- hd[["adj_p_value"]][seq_len(n_var), seq_len(n_contr)]
  normed <- hn[["normed_effect"]][seq_len(n_var), seq_len(n_contr)]
  stopifnot(dim(adj_p) == c(n_var, n_contr))

  self_cells <- cbind(target_col_idx, seq_len(n_contr))
  adj_p[self_cells] <- NA_real_
  sig <- is.finite(adj_p) & (adj_p < THRESH)
  rm(adj_p); invisible(gc(verbose = FALSE))

  cat("accumulating directed significant pair means ...\n")
  P_sum <- matrix(0, n_var, n_var)
  P_n   <- matrix(0, n_var, n_var)
  rows_by_target <- split(seq_len(n_contr), target_col_idx)
  for (target_name in names(rows_by_target)) {
    focal <- as.integer(target_name)
    jj <- rows_by_target[[target_name]]
    n_block <- normed[, jj, drop = FALSE]
    s_block <- sig[, jj, drop = FALSE]
    ns <- n_block; ns[!s_block] <- NA_real_
    P_sum[, focal] <- rowSums(ns, na.rm = TRUE)
    P_n[, focal]   <- rowSums(s_block & is.finite(n_block))
  }
  P <- P_sum / P_n                      # NaN where a pair has no sig cell
  rm(P_sum, normed, sig); invisible(gc(verbose = FALSE))

  has_val <- is.finite(P)
  Pz <- P; Pz[!has_val] <- 0
  rm(P); invisible(gc(verbose = FALSE))

  cnt <- (has_val + 0) + t(has_val + 0)
  TP <- t(Pz)
  ties <- (cnt == 2) & (abs(Pz) == abs(TP))
  n_ties <- sum(ties[upper.tri(ties)])
  if (n_ties > 0) {
    cat(sprintf("TIE CASES (abs(a) == abs(b) in a bidirectional pair): %d\n", n_ties))
    tie_idx <- which(ties & upper.tri(ties), arr.ind = TRUE)
    for (k in seq_len(min(20L, nrow(tie_idx)))) {
      i <- tie_idx[k, 1]; j <- tie_idx[k, 2]
      cat(sprintf("  %s (%s) -- %s (%s): a = %.6g, b = %.6g\n",
                  cols$gene_ensg[i], cols$gene_symbol[i],
                  cols$gene_ensg[j], cols$gene_symbol[j],
                  Pz[i, j], TP[i, j]))
    }
  } else {
    cat("no tie cases (abs(a) == abs(b)) among bidirectional pairs\n")
  }
  W <- Pz
  take <- abs(TP) > abs(W)
  W[take] <- TP[take]
  rm(TP, take); invisible(gc(verbose = FALSE))
  diag(W) <- 0
  dimnames(W) <- list(cols$gene_ensg, cols$gene_ensg)
  stopifnot(isSymmetric(W))

  n_edges  <- sum(W[upper.tri(W)] != 0)
  n_bidir  <- sum((cnt == 2)[upper.tri(cnt)])
  vals <- W[upper.tri(W)]; vals <- vals[vals != 0]
  cat(sprintf("weighted edges: %s  (of which bidirectional, max-|value| taken: %s)\n",
              formatC(n_edges, big.mark = ","), formatC(n_bidir, big.mark = ",")))
  cat("edge-weight quantiles (1/25/50/75/99%):\n")
  print(round(quantile(vals, c(.01, .25, .50, .75, .99)), 4))

  saveRDS(W, OUT_RDS)
  cat(sprintf("wrote %s  (%d x %d matrix)\n", OUT_RDS, nrow(W), ncol(W)))

  if (file.exists(UNION)) {
    adj <- readRDS(UNION)
    stopifnot(identical(rownames(adj), rownames(W)))
    both  <- sum((W != 0) & (adj != 0)) / 2
    w_only <- sum((W != 0) & (adj == 0)) / 2
    a_only <- sum((W == 0) & (adj != 0)) / 2
    cat(sprintf("union-network comparison: shared edges %s, weight-only %s, union-only %s\n",
                formatC(both, big.mark = ","), formatC(w_only, big.mark = ","),
                formatC(a_only, big.mark = ",")))
    if (w_only > 0)
      cat("  WARNING: edges present here but absent from the union network\n")
    if (a_only > 0)
      cat(sprintf("  note: %s union edges have no finite normed_effect at their significant cells (weight 0 here)\n",
                  formatC(a_only, big.mark = ",")))
  }
}

main()
