rm(list=ls())

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

geco_dir = 'data/onek1k'

univ = readLines(file.path(geco_dir, 'GECO_universe.txt'))
N = length(univ)
cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
rd = function(p) { x = readLines(p); intersect(x[nzchar(x)], univ) }

REF = rd(file.path(geco_dir, 'GECO_Shared.txt'))

spec = list(
  'WGCNA'           = list(dir = 'data/others/WGCNA_hubs',
                           pat = '%s.txt', gfile = 'Shared.txt'),
  'Standard glasso' = list(dir = 'data/others/glasso_hubs',
                           pat = 'hubs_%s.txt', gfile = 'hubs_Shared.txt'),
  'GECO'            = list(dir = geco_dir,
                           pat = 'hubs_%s.txt', gfile = 'GECO_Shared.txt'))

fe_test = function(h) {
  a = length(intersect(h, REF))
  p = fisher.test(matrix(c(a, length(h) - a, length(REF) - a,
                           N - length(union(h, REF))), 2, byrow = TRUE),
                  alternative = 'greater')$p.value
  data.frame(n = length(h), n_ref = length(REF), overlap = a,
             fe = (a / length(h)) / (length(REF) / N), p = p)
}

res = NULL
for (m in names(spec)) {
  s = spec[[m]]
  for (net in c('Bulk/Shared', cts))
    res = rbind(res, cbind(data.frame(method = m, net = net),
                           fe_test(rd(file.path(s$dir, if (net == 'Bulk/Shared') s$gfile else sprintf(s$pat, net))))))
}
write.table(res, file.path(outdir, 'hub_fe_vs_geco_shared.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
