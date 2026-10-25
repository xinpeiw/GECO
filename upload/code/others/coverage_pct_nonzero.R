rm(list=ls())

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

degree_coverage_from_adj_nonisolated <- function(A) {
  A <- as.matrix(A)
  B <- (A != 0)
  diag(B) <- FALSE
  deg_all <- rowSums(B, na.rm = TRUE)
  keep <- deg_all > 0
  B2 <- B[keep, keep, drop = FALSE]
  node_names <- rownames(A)[keep]
  n_noniso <- nrow(B2)
  score <- rowSums(B2, na.rm = TRUE)
  ord <- order(score, decreasing = TRUE)
  B_ord <- B2[ord, ord, drop = FALSE]
  covered <- rep(FALSE, n_noniso)
  covered_count <- integer(n_noniso)
  for (k in seq_len(n_noniso)) {
    covered <- covered | B_ord[k, ]
    covered_count[k] <- sum(covered)
  }
  data.frame(top_pct = seq_len(n_noniso) / n_noniso,
             covered_pct = covered_count / n_noniso,
             k = seq_len(n_noniso), score = score[ord], node = node_names[ord])
}

cts <- c('alpha','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
res <- NULL
for (ct in cts) {
  A <- readRDS(file.path("data/onek1k_precision", paste0(ct, ".rds")))
  df <- degree_coverage_from_adj_nonisolated(A)
  res <- rbind(res, cbind(network = ifelse(ct == 'alpha', 'Shared', ct), df))
  cat(sprintf("%-6s non-isolated genes %5d | top 10%% of genes cover %.1f%%\n",
              ct, nrow(df), 100 * df$covered_pct[which.min(abs(df$top_pct - 0.1))]))
}
write.table(res, file.path(outdir, "coverage_pct_nonzero.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
