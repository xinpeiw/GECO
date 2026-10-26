PREC = 'data/onek1k_precision'; TBL = 'raw/sLDSC_sensitivity'
CTS  = c('BIN','BMem','CD4ET','CD4NC','CD8ET','CD8NC','NK')
deg  = function(x) { m = readRDS(file.path(PREC, paste0(x, '.rds'))); diag(m) = 0; rowSums(m != 0) }
shared = readLines('data/onek1k/GECO_Shared.txt')
spec   = readLines('data/onek1k/GECO_Specific.txt')
pct    = sapply(CTS, function(ct) { d = deg(ct); (rank(d) / length(d))[spec] })
score  = list(Shared = deg('alpha')[shared], Specific = apply(pct, 1, max, na.rm = TRUE))
for (nm in names(score)) {
  g = names(score[[nm]])[order(-score[[nm]])]
  for (N in c(210, 420)) writeLines(g[seq_len(N)], file.path(TBL, sprintf('GECO_%s_top%d.txt', nm, N)))
}
