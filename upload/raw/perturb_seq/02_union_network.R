# Input not included in this repository (934 MB): raw/perturb_seq/primary_geco_intersection.h5 (md5 d1cc65777e32762af26bf971f2e5c536).
#   It holds the log_fc / zscore / adj_p_value layers of GWCD4i.DE_stats.h5ad (genome-wide CD4+ T cell CRISPRi Perturb-seq,
#   Zhu et al. 2025) for the 17,525 primary contrasts whose target is in the GECO universe x the 6,603 intersection genes;
#   row/column metadata are the included primary_geco_intersection_{rows,cols}.tsv.
suppressPackageStartupMessages({
  library(hdf5r)
  library(data.table)
})

DIR         <- "raw/perturb_seq"
H5_DAT      <- file.path(DIR, "primary_geco_intersection.h5")
ROWS        <- file.path(DIR, "primary_geco_intersection_rows.tsv")
COLS        <- file.path(DIR, "primary_geco_intersection_cols.tsv")
OUT_MATRIX  <- "data/others/pertseq_union_network_matrix.rds"
OUT_EDGES   <- "data/6c-hub_position/perturb_seq/pertseq_union_network_4col.tsv"
THRESH      <- 0.10
BLOCK_SIZE  <- 256L
CONDS       <- c("Rest", "Stim8hr", "Stim48hr")

main <- function() {
  rows <- fread(ROWS)
  cols <- fread(COLS)
  n_contr <- nrow(rows)
  n_var <- nrow(cols)
  cat(sprintf("contrasts: %s   intersection genes: %s\n",
              formatC(n_contr, big.mark = ","),
              formatC(n_var, big.mark = ",")))

  if (anyDuplicated(cols$gene_ensg))
    stop("duplicated gene_ensg values found in intersection columns")
  if (!all(rows$condition %in% CONDS))
    stop("unexpected culture condition found in row metadata")

  target_col_idx <- match(rows$target_ensg, cols$gene_ensg)
  if (anyNA(target_col_idx))
    stop("at least one contrast target is missing from intersection columns")

  symbol_idx <- setNames(cols$gene_symbol, cols$gene_ensg)

  h <- H5File$new(H5_DAT, mode = "r")
  on.exit(h$close_all(), add = TRUE)
  stopifnot(
    identical(h[["zscore"]]$dims, c(n_var, n_contr)),
    identical(h[["adj_p_value"]]$dims, c(n_var, n_contr))
  )

  cat("extracting significant non-self cells ...\n")
  cell_blocks <- vector("list", ceiling(n_contr / BLOCK_SIZE))
  block_i <- 0L
  total_sig <- 0L

  for (start in seq.int(1L, n_contr, by = BLOCK_SIZE)) {
    end <- min(n_contr, start + BLOCK_SIZE - 1L)
    contrast_idx <- start:end

    adj_block <- h[["adj_p_value"]][seq_len(n_var), contrast_idx]
    z_block <- h[["zscore"]][seq_len(n_var), contrast_idx]

    self_cells <- cbind(target_col_idx[contrast_idx],
                        seq_along(contrast_idx))
    adj_block[self_cells] <- NA_real_
    z_block[self_cells] <- NA_real_

    sig_pos <- which(is.finite(adj_block) & adj_block < THRESH,
                     arr.ind = TRUE)
    block_i <- block_i + 1L

    if (nrow(sig_pos) == 0L) {
      cell_blocks[[block_i]] <- NULL
      next
    }

    measured_idx <- sig_pos[, 1L]
    global_contrast_idx <- contrast_idx[sig_pos[, 2L]]
    target_ensg <- rows$target_ensg[global_contrast_idx]
    measured_ensg <- cols$gene_ensg[measured_idx]
    target_is_a <- target_ensg < measured_ensg

    cell_blocks[[block_i]] <- data.table(
      gene_a_ensg = fifelse(target_is_a, target_ensg, measured_ensg),
      gene_b_ensg = fifelse(target_is_a, measured_ensg, target_ensg),
      condition = rows$condition[global_contrast_idx],
      dir_a_to_b = target_is_a,
      adj_p = adj_block[sig_pos],
      abs_z = abs(z_block[sig_pos])
    )
    total_sig <- total_sig + nrow(sig_pos)
  }

  cells <- rbindlist(cell_blocks, use.names = TRUE)
  rm(cell_blocks)
  invisible(gc(verbose = FALSE))

  cat(sprintf("significant cells (self-loop excluded): %s\n",
              formatC(total_sig, big.mark = ",")))
  if (nrow(cells) == 0L)
    stop("no significant non-self cells found")

  cat(sprintf("aggregating %s significant cells into undirected edges ...\n",
              formatC(nrow(cells), big.mark = ",")))

  edge_df <- cells[, .(
    in_Rest = any(condition == "Rest"),
    in_Stim8hr = any(condition == "Stim8hr"),
    in_Stim48hr = any(condition == "Stim48hr"),
    has_dir_a_to_b = any(dir_a_to_b),
    has_dir_b_to_a = any(!dir_a_to_b),
    n_sig_cells = .N,
    min_adj_p = min(adj_p),
    max_abs_zscore = max(abs_z)
  ), by = .(gene_a_ensg, gene_b_ensg)]

  edge_df[, `:=`(
    gene_a_symbol = unname(symbol_idx[gene_a_ensg]),
    gene_b_symbol = unname(symbol_idx[gene_b_ensg]),
    n_conditions = as.integer(in_Rest) + as.integer(in_Stim8hr) +
                   as.integer(in_Stim48hr),
    n_directions = as.integer(has_dir_a_to_b) +
                   as.integer(has_dir_b_to_a)
  )]

  setcolorder(edge_df, c(
    "gene_a_ensg", "gene_b_ensg",
    "gene_a_symbol", "gene_b_symbol",
    "in_Rest", "in_Stim8hr", "in_Stim48hr", "n_conditions",
    "has_dir_a_to_b", "has_dir_b_to_a", "n_directions",
    "n_sig_cells", "min_adj_p", "max_abs_zscore"
  ))
  setorder(edge_df, gene_a_ensg, gene_b_ensg)

  cat(sprintf("unique undirected edges: %s\n",
              formatC(nrow(edge_df), big.mark = ",")))

  gene_idx <- setNames(seq_len(n_var), cols$gene_ensg)
  edge_i <- unname(gene_idx[edge_df$gene_a_ensg])
  edge_j <- unname(gene_idx[edge_df$gene_b_ensg])
  if (anyNA(edge_i) || anyNA(edge_j))
    stop("network edge contains a gene absent from intersection columns")

  write_adjacency <- function(name, keep_edge, rds_path) {
    adjacency <- matrix(0L, nrow = n_var, ncol = n_var,
                        dimnames = list(cols$gene_ensg, cols$gene_ensg))
    i <- edge_i[keep_edge]
    j <- edge_j[keep_edge]
    adjacency[cbind(i, j)] <- 1L
    adjacency[cbind(j, i)] <- 1L
    diag(adjacency) <- 0L
    saveRDS(adjacency, rds_path)
    n_edges <- sum(keep_edge)
    rm(adjacency)
    invisible(gc(verbose = FALSE))
    cat(sprintf("  wrote %-8s %s (%d x %d matrix; %s undirected edges)\n",
                name, rds_path, n_var, n_var,
                formatC(n_edges, big.mark = ",")))
  }

  cat("writing adjacency matrices as RDS ...\n")
  write_adjacency("union", rep(TRUE, nrow(edge_df)), OUT_MATRIX)
  fwrite(edge_df[, .(gene_a_ensg, gene_b_ensg, has_dir_a_to_b, has_dir_b_to_a)], OUT_EDGES, sep = "\t")
}

main()
