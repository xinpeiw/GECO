rm(list=ls())

library(readxl)

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

R = 'data/others'
cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')

res = NULL

for (net in c('Shared', cts)) {
  p = file.path(R, 'WGCNA', paste0(net, '.xlsx'))
  d = do.call(rbind, lapply(excel_sheets(p),
                            function(s) as.data.frame(read_excel(p, sheet = s))))
  d = d[!duplicated(d$Gene), ]
  res = rbind(res, data.frame(method = 'WGCNA',
                              net = ifelse(net == 'Shared', 'Bulk/Shared', net),
                              connected = sum(d$Module != 'grey'), total = nrow(d)))
}

for (i in 1:2) {
  meth = c('Standard glasso','GECO')[i]
  f = c(file.path(R, 'real_data_naive_glasso_degree.xlsx'),
        file.path(R, 'wglasso_degree.xlsx'))[i]
  gl = c('Shared','alpha')[i]
  for (net in c(gl, cts)) {
    d = as.data.frame(read_excel(f, sheet = net))
    res = rbind(res, data.frame(method = meth,
                                net = ifelse(net == gl, 'Bulk/Shared', net),
                                connected = sum(d$degree > 0), total = nrow(d)))
  }
}
res$pct = res$connected / res$total
write.table(res, file.path(outdir, 'connected_genes_3methods.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
