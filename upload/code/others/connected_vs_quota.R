rm(list=ls())

library(readxl)

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

wgcna_dir = 'data/2b_hubs_edges/WGCNA'
glasso_file = 'data/2b_hubs_edges/naive_glasso/degree_naive_glasso.xlsx'
geco_file = 'data/2b_hubs_edges/geco/degree_geco_merged_100reps.xlsx'

cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
nets = c('OP', cts)
Ktrue = c(OP=482, CD4NC=224, CD4ET=78, CD8NC=210, CD8ET=272, NK=297, BIN=122, BMem=96)

df = NULL
for (net in nets) {
  M = as.data.frame(read_excel(file.path(wgcna_dir, paste0(net, '_Module_merged_100reps.xlsx'))))
  n = sapply(sprintf('Module_rep%03d', 1:100),
             function(cl) sum(!is.na(M[[cl]]) & M[[cl]] != 'grey'))
  df = rbind(df, data.frame(method = 'WGCNA', net = net, rep = 1:100, n = as.numeric(n)))
}
for (meth in c('standard glasso','GECO')) {
  f = if (meth == 'GECO') geco_file else glasso_file
  for (net in nets) {
    sh = if (meth == 'GECO' && net == 'OP') 'alpha' else net
    d = as.data.frame(read_excel(f, sheet = sh))
    n = sapply(grep('^degree_rep', colnames(d), value = TRUE),
               function(cl) sum(d[[cl]] > 0))
    df = rbind(df, data.frame(method = meth, net = net, rep = 1:100, n = as.numeric(n)))
  }
}
df$K = Ktrue[df$net]
df$ratio = df$n / df$K
df$below = df$n < df$K
write.table(df, file.path(outdir, 'connected_vs_quota.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
