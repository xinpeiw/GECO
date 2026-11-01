rm(list=ls())

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

geco_dir = 'data/onek1k'

univ = readLines(file.path(geco_dir, 'GECO_universe.txt'))
N = length(univ)
cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')

rd = function(p) { x = readLines(p); intersect(x[nzchar(x)], univ) }

spec = list(
  'WGCNA'           = list(dir = 'data/others/WGCNA_hubs',
                           pat = '%s.txt', gfile = 'Shared.txt'),
  'Standard glasso' = list(dir = 'data/others/glasso_hubs',
                           pat = 'hubs_%s.txt', gfile = 'hubs_Shared.txt'),
  'GECO'            = list(dir = geco_dir,
                           pat = 'hubs_%s.txt', gfile = 'GECO_Shared.txt'))

res = NULL
for (m in names(spec)) {
  s = spec[[m]]
  gl = rd(file.path(s$dir, s$gfile))
  for (ct in cts) {
    h = rd(file.path(s$dir, sprintf(s$pat, ct)))
    a = length(intersect(h, gl))
    fe = (a / length(h)) / (length(gl) / N)
    p = fisher.test(matrix(c(a, length(h) - a, length(gl) - a,
                             N - length(union(h, gl))), 2, byrow = TRUE),
                    alternative = 'greater')$p.value
    res = rbind(res, data.frame(method = m, ct = ct, n = length(h),
                                n_global = length(gl), overlap = a, fe = fe, p = p))
  }
}
write.table(res, file.path(outdir, 'hub_fe_3methods.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
