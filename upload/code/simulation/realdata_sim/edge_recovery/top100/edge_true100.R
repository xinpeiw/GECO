cts = c('alpha','CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK', 'BIN', 'BMem')

maxK = 100

# Zeros are not edges. order(.)[seq_len(maxK)] would still return maxK indices
# when fewer than maxK entries are non-zero, padding the list with zero-valued
# pairs. The list is capped at the number of non-zero entries and the remaining
# slots are left NA, so the table keeps its maxK x 8 shape while never passing
# a zero off as an edge. Anything reading this file must drop the NAs.
top_nz <- function(v, K) {
  v <- abs(v)
  nz <- sum(v > 0, na.rm = TRUE)
  out <- rep(NA_integer_, K)
  k <- min(K, nz)
  if (k > 0) out[seq_len(k)] <- order(v, decreasing = TRUE, na.last = NA)[seq_len(k)]
  out
}

edge_true100 = matrix(NA_integer_, nrow = maxK, ncol = length(cts))
for (i in 1:length(cts)) {

	print(i)

	edge <- readRDS(paste0('geco/edge_cor_', cts[i], '.rds'))

	edge_true100[,i] = top_nz(edge[[1]], maxK)

	cat(sprintf('  %-6s non-zero true edges kept: %d / %d\n',
	            cts[i], sum(!is.na(edge_true100[,i])), maxK))
}

edge_true100 = as.data.frame(edge_true100)
colnames(edge_true100) = cts

write.table(edge_true100, 'edge_true100.txt', row.names = FALSE, col.names = TRUE,
	sep = ' ', quote = FALSE)
