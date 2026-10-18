rm(list=ls())

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

geco_dir = 'data/onek1k'

univ = readLines(file.path(geco_dir, 'GECO_universe.txt'))
cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
rd = function(p) { x = readLines(p); length(intersect(x[nzchar(x)], univ)) }

spec = list('WGCNA'           = list(d = 'data/others/WGCNA_hubs',
                                     p = '%s.txt', gfile = 'Shared.txt'),
            'Standard glasso' = list(d = 'data/others/glasso_hubs',
                                     p = 'hubs_%s.txt', gfile = 'hubs_Shared.txt'),
            'GECO'            = list(d = geco_dir,
                                     p = 'hubs_%s.txt', gfile = 'GECO_Shared.txt'))

res = NULL
for (m in names(spec)) {
  s = spec[[m]]
  for (net in c('Bulk/Shared', cts)) {
    res = rbind(res, data.frame(method = m,
                                net = net,
                                n = rd(file.path(s$d, if (net == 'Bulk/Shared') s$gfile else sprintf(s$p, net)))))
  }
}
write.table(res, file.path(outdir, 'hub_number_3methods.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
