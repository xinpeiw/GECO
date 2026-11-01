args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

cts_other <- c('OP', 'CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK', 'BIN', 'BMem')

maxK <- 100

edge <- data.table::fread(
  paste0('./WGCNA/edge_cor_WGCNA_', cts_other[i], '.txt')
)

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

res <- vapply(seq_along(edge), function(j) f(edge[[j]]), integer(maxK))
colnames(res) <- sprintf("rep_%03d", seq_len(ncol(res)))

cat(sprintf('%s: non-zero edges kept per replicate, min %d / median %d / max %d\n',
            cts_other[i], min(colSums(!is.na(res))), median(colSums(!is.na(res))),
            max(colSums(!is.na(res)))))

write.table(
  res,
  paste0('./WGCNA/edge_100_', cts_other[i], '.txt'),
  row.names = FALSE,
  col.names = TRUE,
  sep = ' ',
  quote = FALSE
)
