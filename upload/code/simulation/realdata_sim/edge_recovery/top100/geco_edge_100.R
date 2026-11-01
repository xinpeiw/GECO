args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

cts = c('alpha','CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK', 'BIN', 'BMem')

edge <- readRDS(paste0('geco/edge_cor_', cts[i], '.rds'))

maxK <- 100

est = edge[,2:101]; rm(edge)

# Zeros are not edges: cap at the number of non-zero entries and leave the rest
# NA rather than padding the top-K list with zero-valued pairs. Anything
# reading this file must drop the NAs.
f <- function(x) {
  x <- abs(x)
  nz <- sum(x > 0, na.rm = TRUE)
  out <- rep(NA_integer_, maxK)
  k <- min(maxK, nz)
  if (k > 0) out[seq_len(k)] <- order(x, decreasing = TRUE, na.last = NA)[seq_len(k)]
  out
}

res <- vapply(seq_along(est), function(j) f(est[[j]]), integer(maxK))
colnames(res) <- colnames(est)

cat(sprintf('%s: non-zero edges kept per replicate, min %d / median %d / max %d\n',
            cts[i], min(colSums(!is.na(res))), median(colSums(!is.na(res))),
            max(colSums(!is.na(res)))))

write.table(res, paste0('geco/edge_100_',cts[i],'.txt'),
	row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
