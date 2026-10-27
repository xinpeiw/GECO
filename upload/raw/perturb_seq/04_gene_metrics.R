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
OUT_TSV <- "data/6c-hub_position/perturb_seq/aggregate_per_intersection_gene_metrics.tsv"
THRESH  <- 0.10
CONDS   <- c("Rest", "Stim8hr", "Stim48hr")

safe_ratio <- function(numerator, denominator) {
  out <- rep(NA_real_, length(numerator))
  ok <- denominator > 0
  out[ok] <- numerator[ok] / denominator[ok]
  out
}

nan_to_na <- function(x) { x[is.nan(x)] <- NA_real_; x }

pair_dir_means <- function(P_sum, P_n) {
  P <- P_sum / P_n
  P[!is.finite(P)] <- NA_real_
  fin <- is.finite(P)
  list(
    incoming = safe_ratio(rowSums(P, na.rm = TRUE), rowSums(fin)),
    outgoing = safe_ratio(colSums(P, na.rm = TRUE), colSums(fin))
  )
}

main <- function() {
  rows <- read.table(ROWS, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cols <- read.table(COLS, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  n_contr <- nrow(rows); n_var <- nrow(cols)
  cat(sprintf("contrasts: %s   intersection genes: %s\n",
              formatC(n_contr, big.mark = ","), formatC(n_var, big.mark = ",")))
  stopifnot(all(rows$condition %in% CONDS))

  target_col_idx <- match(rows$target_ensg, cols$gene_ensg)
  if (anyNA(target_col_idx))
    stop("at least one contrast target is missing from the intersection columns")

  cat("loading layers ...\n")
  hd <- H5File$new(H5_DAT, mode = "r")
  hn <- H5File$new(H5_NORM, mode = "r")
  on.exit({ hd$close_all(); hn$close_all() }, add = TRUE)
  zscore <- hd[["zscore"]][seq_len(n_var), seq_len(n_contr)]
  adj_p  <- hd[["adj_p_value"]][seq_len(n_var), seq_len(n_contr)]
  normed <- hn[["normed_effect"]][seq_len(n_var), seq_len(n_contr)]
  stopifnot(dim(zscore) == c(n_var, n_contr))

  self_cells <- cbind(target_col_idx, seq_len(n_contr))
  zscore[self_cells] <- NA_real_
  adj_p [self_cells] <- NA_real_
  abs_z <- abs(zscore)
  rm(zscore); invisible(gc(verbose = FALSE))

  sig <- is.finite(adj_p) & (adj_p < THRESH)
  cond_idx_of_contrast <- match(rows$condition, CONDS)

  cat("computing incoming (regulator-side) sums ...\n")
  in_n_z       <- rowSums(is.finite(abs_z))
  in_n_n       <- rowSums(is.finite(normed))
  in_sum_z     <- rowSums(abs_z,  na.rm = TRUE)
  in_sum_n     <- rowSums(normed, na.rm = TRUE)
  in_n_sig     <- rowSums(sig)
  in_n_z_sig   <- rowSums(sig & is.finite(abs_z))
  in_n_n_sig   <- rowSums(sig & is.finite(normed))
  in_sum_z_sig <- rowSums(ifelse(sig, abs_z, 0),  na.rm = TRUE)
  in_sum_n_sig <- rowSums(ifelse(sig, normed, 0), na.rm = TRUE)

  cat("computing outgoing (downstream-side) sums ...\n")
  out_n_z       <- numeric(n_var)
  out_n_n       <- numeric(n_var)
  out_sum_z     <- numeric(n_var)
  out_sum_n     <- numeric(n_var)
  out_n_sig     <- numeric(n_var)
  out_n_z_sig   <- numeric(n_var)
  out_n_n_sig   <- numeric(n_var)
  out_sum_z_sig <- numeric(n_var)
  out_sum_n_sig <- numeric(n_var)

  D <- matrix(FALSE, nrow = n_var, ncol = n_var)

  P_sum_z     <- matrix(0, n_var, n_var)
  P_n_z       <- matrix(0, n_var, n_var)
  P_sum_n     <- matrix(0, n_var, n_var)
  P_n_n       <- matrix(0, n_var, n_var)
  P_sum_z_sig <- matrix(0, n_var, n_var)
  P_n_z_sig   <- matrix(0, n_var, n_var)
  P_sum_n_sig <- matrix(0, n_var, n_var)
  P_n_n_sig   <- matrix(0, n_var, n_var)

  rows_by_target <- split(seq_len(n_contr), target_col_idx)
  for (target_name in names(rows_by_target)) {
    focal <- as.integer(target_name)
    jj <- rows_by_target[[target_name]]

    z_block <- abs_z[, jj, drop = FALSE]
    n_block <- normed[, jj, drop = FALSE]
    s_block <- sig[, jj, drop = FALSE]

    out_n_z[focal]       <- sum(is.finite(z_block))
    out_n_n[focal]       <- sum(is.finite(n_block))
    out_sum_z[focal]     <- sum(z_block, na.rm = TRUE)
    out_sum_n[focal]     <- sum(n_block, na.rm = TRUE)
    out_n_sig[focal]     <- sum(s_block)
    out_n_z_sig[focal]   <- sum(s_block & is.finite(z_block))
    out_n_n_sig[focal]   <- sum(s_block & is.finite(n_block))
    out_sum_z_sig[focal] <- sum(z_block[s_block], na.rm = TRUE)
    out_sum_n_sig[focal] <- sum(n_block[s_block], na.rm = TRUE)

    D[focal, ] <- rowSums(s_block) > 0L

    P_sum_z[, focal] <- rowSums(z_block, na.rm = TRUE)
    P_n_z[, focal]   <- rowSums(is.finite(z_block))
    P_sum_n[, focal] <- rowSums(n_block, na.rm = TRUE)
    P_n_n[, focal]   <- rowSums(is.finite(n_block))
    zb_sig <- z_block; zb_sig[!s_block] <- NA_real_
    nb_sig <- n_block; nb_sig[!s_block] <- NA_real_
    P_sum_z_sig[, focal] <- rowSums(zb_sig, na.rm = TRUE)
    P_n_z_sig[, focal]   <- rowSums(s_block & is.finite(z_block))
    P_sum_n_sig[, focal] <- rowSums(nb_sig, na.rm = TRUE)
    P_n_n_sig[, focal]   <- rowSums(s_block & is.finite(n_block))
  }
  diag(D) <- FALSE

  n_outgoing_neighbors <- rowSums(D)             # KD(g) -> partner significant
  n_incoming_neighbors <- colSums(D)             # KD(partner) -> g significant
  n_neighbors          <- rowSums(D | t(D))      # union of both directions
  n_neighbors_bidir    <- rowSums(D & t(D))      # significant both ways

  pooled_z_in      <- safe_ratio(in_sum_z, in_n_z)
  pooled_n_in      <- safe_ratio(in_sum_n, in_n_n)
  pooled_z_sig_in  <- safe_ratio(in_sum_z_sig, in_n_z_sig)
  pooled_n_sig_in  <- safe_ratio(in_sum_n_sig, in_n_n_sig)
  pooled_z_out     <- safe_ratio(out_sum_z, out_n_z)
  pooled_n_out     <- safe_ratio(out_sum_n, out_n_n)
  pooled_z_sig_out <- safe_ratio(out_sum_z_sig, out_n_z_sig)
  pooled_n_sig_out <- safe_ratio(out_sum_n_sig, out_n_n_sig)

  cat("computing per-condition means for condition averaging ...\n")
  cm_z_in      <- matrix(NA_real_, n_var, length(CONDS))
  cm_n_in      <- matrix(NA_real_, n_var, length(CONDS))
  cm_z_sig_in  <- matrix(NA_real_, n_var, length(CONDS))
  cm_n_sig_in  <- matrix(NA_real_, n_var, length(CONDS))
  cm_z_out     <- matrix(NA_real_, n_var, length(CONDS))
  cm_n_out     <- matrix(NA_real_, n_var, length(CONDS))
  cm_z_sig_out <- matrix(NA_real_, n_var, length(CONDS))
  cm_n_sig_out <- matrix(NA_real_, n_var, length(CONDS))

  for (ci in seq_along(CONDS)) {
    idx <- which(cond_idx_of_contrast == ci)
    if (length(idx) == 0L)
      stop(sprintf("no contrasts found for condition %s", CONDS[ci]))
    z_c <- abs_z[, idx, drop = FALSE]
    n_c <- normed[, idx, drop = FALSE]
    s_c <- sig[, idx, drop = FALSE]

    cm_z_in[, ci]     <- safe_ratio(rowSums(z_c, na.rm = TRUE),
                                    rowSums(is.finite(z_c)))
    cm_n_in[, ci]     <- safe_ratio(rowSums(n_c, na.rm = TRUE),
                                    rowSums(is.finite(n_c)))
    cm_z_sig_in[, ci] <- safe_ratio(rowSums(ifelse(s_c, z_c, 0), na.rm = TRUE),
                                    rowSums(s_c & is.finite(z_c)))
    cm_n_sig_in[, ci] <- safe_ratio(rowSums(ifelse(s_c, n_c, 0), na.rm = TRUE),
                                    rowSums(s_c & is.finite(n_c)))

    o_n_z <- numeric(n_var); o_s_z <- numeric(n_var)
    o_n_n <- numeric(n_var); o_s_n <- numeric(n_var)
    o_n_z_sig <- numeric(n_var); o_s_z_sig <- numeric(n_var)
    o_n_n_sig <- numeric(n_var); o_s_n_sig <- numeric(n_var)
    by_target_c <- split(seq_along(idx), target_col_idx[idx])
    for (target_name in names(by_target_c)) {
      focal <- as.integer(target_name)
      jj <- by_target_c[[target_name]]
      zb <- z_c[, jj, drop = FALSE]
      nb <- n_c[, jj, drop = FALSE]
      sb <- s_c[, jj, drop = FALSE]
      o_n_z[focal]     <- sum(is.finite(zb))
      o_s_z[focal]     <- sum(zb, na.rm = TRUE)
      o_n_n[focal]     <- sum(is.finite(nb))
      o_s_n[focal]     <- sum(nb, na.rm = TRUE)
      o_n_z_sig[focal] <- sum(sb & is.finite(zb))
      o_s_z_sig[focal] <- sum(zb[sb], na.rm = TRUE)
      o_n_n_sig[focal] <- sum(sb & is.finite(nb))
      o_s_n_sig[focal] <- sum(nb[sb], na.rm = TRUE)
    }
    cm_z_out[, ci]     <- safe_ratio(o_s_z, o_n_z)
    cm_n_out[, ci]     <- safe_ratio(o_s_n, o_n_n)
    cm_z_sig_out[, ci] <- safe_ratio(o_s_z_sig, o_n_z_sig)
    cm_n_sig_out[, ci] <- safe_ratio(o_s_n_sig, o_n_n_sig)

    rm(z_c, n_c, s_c)
    invisible(gc(verbose = FALSE))
  }

  cond_z_in      <- nan_to_na(rowMeans(cm_z_in,      na.rm = TRUE))
  cond_n_in      <- nan_to_na(rowMeans(cm_n_in,      na.rm = TRUE))
  cond_z_sig_in  <- nan_to_na(rowMeans(cm_z_sig_in,  na.rm = TRUE))
  cond_n_sig_in  <- nan_to_na(rowMeans(cm_n_sig_in,  na.rm = TRUE))
  cond_z_out     <- nan_to_na(rowMeans(cm_z_out,     na.rm = TRUE))
  cond_n_out     <- nan_to_na(rowMeans(cm_n_out,     na.rm = TRUE))
  cond_z_sig_out <- nan_to_na(rowMeans(cm_z_sig_out, na.rm = TRUE))
  cond_n_sig_out <- nan_to_na(rowMeans(cm_n_sig_out, na.rm = TRUE))

  cat("computing pair-averaged means ...\n")
  pr_z <- pair_dir_means(P_sum_z, P_n_z)
  rm(P_sum_z, P_n_z); invisible(gc(verbose = FALSE))
  pr_n <- pair_dir_means(P_sum_n, P_n_n)
  rm(P_sum_n, P_n_n); invisible(gc(verbose = FALSE))
  pr_z_sig <- pair_dir_means(P_sum_z_sig, P_n_z_sig)
  rm(P_sum_z_sig, P_n_z_sig); invisible(gc(verbose = FALSE))
  pr_n_sig <- pair_dir_means(P_sum_n_sig, P_n_n_sig)
  rm(P_sum_n_sig, P_n_n_sig); invisible(gc(verbose = FALSE))

  out <- data.frame(
    gene_ensg             = cols$gene_ensg,
    gene_symbol           = cols$gene_symbol,
    n_sig_incoming_cells  = as.integer(in_n_sig),
    n_sig_outgoing_cells  = as.integer(out_n_sig),
    n_incoming_neighbors  = as.integer(n_incoming_neighbors),
    n_outgoing_neighbors  = as.integer(n_outgoing_neighbors),
    n_neighbors           = as.integer(n_neighbors),
    n_neighbors_bidir     = as.integer(n_neighbors_bidir),

    mean_abs_zscore_in           = pooled_z_in,
    mean_normed_effect_in        = pooled_n_in,
    mean_abs_zscore_sig_in       = pooled_z_sig_in,
    mean_normed_effect_sig_in    = pooled_n_sig_in,
    mean_abs_zscore_cond_in       = cond_z_in,
    mean_normed_effect_cond_in    = cond_n_in,
    mean_abs_zscore_sig_cond_in   = cond_z_sig_in,
    mean_normed_effect_sig_cond_in = cond_n_sig_in,
    mean_abs_zscore_pair_in      = pr_z$incoming,
    mean_normed_effect_pair_in   = pr_n$incoming,
    mean_abs_zscore_sig_pair_in  = pr_z_sig$incoming,
    mean_normed_effect_sig_pair_in = pr_n_sig$incoming,

    mean_abs_zscore_out          = pooled_z_out,
    mean_normed_effect_out       = pooled_n_out,
    mean_abs_zscore_sig_out      = pooled_z_sig_out,
    mean_normed_effect_sig_out   = pooled_n_sig_out,
    mean_abs_zscore_cond_out      = cond_z_out,
    mean_normed_effect_cond_out   = cond_n_out,
    mean_abs_zscore_sig_cond_out  = cond_z_sig_out,
    mean_normed_effect_sig_cond_out = cond_n_sig_out,
    mean_abs_zscore_pair_out     = pr_z$outgoing,
    mean_normed_effect_pair_out  = pr_n$outgoing,
    mean_abs_zscore_sig_pair_out = pr_z_sig$outgoing,
    mean_normed_effect_sig_pair_out = pr_n_sig$outgoing,

    stringsAsFactors      = FALSE
  )

  no_sig_out <- out$n_sig_outgoing_cells == 0L
  out$mean_normed_effect_sig_out[no_sig_out]      <- 0
  out$mean_normed_effect_sig_cond_out[no_sig_out] <- 0
  out$mean_normed_effect_sig_pair_out[no_sig_out] <- 0
  cat(sprintf("genes with no significant outgoing cell (normed sig_out set to 0): %s\n",
              formatC(sum(no_sig_out), big.mark = ",")))

  write.table(out, OUT_TSV, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("\nwrote %s (%s rows x %d cols)\n",
              OUT_TSV, formatC(nrow(out), big.mark = ","), ncol(out)))
}

main()
