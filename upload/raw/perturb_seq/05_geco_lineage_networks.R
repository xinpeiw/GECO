suppressPackageStartupMessages({ library(data.table) })

PREC <- "data/onek1k_precision"
cols <- read.table("raw/perturb_seq/primary_geco_intersection_cols.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
inter <- cols$gene_ensg
sym_of_ensg <- setNames(cols$gene_symbol, cols$gene_ensg)
n_var <- length(inter)

pc_edges <- function(net) {
  om <- readRDS(file.path(PREC, paste0(net, ".rds")))[inter, inter]
  d <- diag(om); diag(om) <- 0
  idx <- which(om != 0 & upper.tri(om), arr.ind = TRUE)
  a <- inter[idx[, 1]]; b <- inter[idx[, 2]]
  data.table(gene_a_ensg = pmin(a, b), gene_b_ensg = pmax(a, b), net = net,
             pc = -om[idx] / sqrt(d[idx[, 1]] * d[idx[, 2]]))
}

build <- function(networks, output_fields, out_rds, out_tsv) {
  e <- rbindlist(lapply(networks, pc_edges))
  edges <- unique(e[, .(gene_a_ensg, gene_b_ensg)])
  setorder(edges, gene_a_ensg, gene_b_ensg)
  df <- data.frame(gene_a_ensg = edges$gene_a_ensg, gene_b_ensg = edges$gene_b_ensg,
                   gene_a_symbol = unname(sym_of_ensg[edges$gene_a_ensg]),
                   gene_b_symbol = unname(sym_of_ensg[edges$gene_b_ensg]), stringsAsFactors = FALSE)
  key <- paste(edges$gene_a_ensg, edges$gene_b_ensg)
  for (nt in output_fields) df[[paste0("in_", nt)]] <- key %in% e[net == nt, paste(gene_a_ensg, gene_b_ensg)]
  for (nt in output_fields) {
    sub <- e[net == nt]
    df[[paste0("pc_", nt)]] <- sub$pc[match(key, paste(sub$gene_a_ensg, sub$gene_b_ensg))]
  }
  df$n_networks <- rowSums(df[paste0("in_", output_fields)])
  adjacency <- matrix(0L, n_var, n_var, dimnames = list(inter, inter))
  gi <- setNames(seq_len(n_var), inter)
  i <- unname(gi[df$gene_a_ensg]); j <- unname(gi[df$gene_b_ensg])
  adjacency[cbind(i, j)] <- 1L; adjacency[cbind(j, i)] <- 1L
  saveRDS(adjacency, out_rds)
  write.table(df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
}

build(c("CD4NC", "CD4ET"), c("CD4NC", "CD4ET", "alpha"),
      "data/others/geco_cd4_union_network_matrix.rds", "data/others/geco_cd4_union_network.tsv")
build(c("CD8NC", "CD8ET"), c("CD8NC", "CD8ET"),
      "data/others/geco_cd8_union_network_matrix.rds", "data/others/geco_cd8_union_network.tsv")
build("NK", "NK", "data/others/geco_nk_network_matrix.rds", "data/others/geco_nk_network.tsv")
