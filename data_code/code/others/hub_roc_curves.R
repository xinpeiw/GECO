rm(list=ls())

library(readxl)

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

sim = 'data/2b_hubs_edges'
wg_dir  = file.path(sim, 'WGCNA')
gl_file = file.path(sim, 'naive_glasso/degree_naive_glasso.xlsx')
gc_file = file.path(sim, 'geco/degree_geco_merged_100reps.xlsx')
geco_dir = 'data/onek1k'

cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
shared_hub = readLines(file.path(geco_dir, 'GECO_Shared.txt'))
hubs_of = function(net) {
  x = if (net == 'OP') shared_hub else
      readLines(file.path(geco_dir, paste0('hubs_', net, '.txt')))
  x[nzchar(x)]
}
grid = seq(0, 1, length.out = 201)

roc_on_grid = function(score, lab) {
  o = order(score, decreasing = TRUE)
  s = score[o]; l = lab[o]

  last = c(s[-length(s)] != s[-1], TRUE)
  tpr = c(0, cumsum(l)[last] / sum(l))
  fpr = c(0, cumsum(!l)[last] / sum(!l))
  approx(fpr, tpr, xout = grid, ties = 'ordered', rule = 2)$y
}

auroc = function(score, lab) {
  r = rank(score); n1 = sum(lab); n0 = sum(!lab)
  (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

curves = NULL; au = NULL
add = function(method, net, scores, genes) {
  lab = genes %in% hubs_of(net)
  M = vapply(scores, roc_on_grid, numeric(length(grid)), lab = lab)
  curves <<- rbind(curves, data.frame(
    method = method, net = net, fpr = grid,
    tpr = rowMeans(M), sd = apply(M, 1, sd)))
  au <<- rbind(au, data.frame(method = method, net = net,
                              auroc = mean(vapply(scores, auroc, numeric(1), lab = lab))))
}

wgcna_score = function(mod, kme) {
  s = !is.na(mod) & mod != 'grey' & !is.na(kme)
  out = numeric(length(kme))
  if (any(s)) {
    k = kme[s]
    out[s] = unsplit(lapply(split(k, mod[s]), function(x) rank(x) / length(x)),
                     mod[s])
  }
  out
}
for (net in c('OP', cts)) {
  M  = as.data.frame(read_excel(file.path(wg_dir, paste0(net, '_Module_merged_100reps.xlsx'))))
  Km = as.data.frame(read_excel(file.path(wg_dir, paste0(net, '_kME_merged_100reps.xlsx'))))
  stopifnot(identical(M$Gene, Km$Gene))
  sc = lapply(sprintf('rep%03d', 1:100), function(r)
    wgcna_score(M[[paste0('Module_', r)]], Km[[paste0('kME_', r)]]))
  add('WGCNA', net, sc, M$Gene)
}
for (meth in c('Standard glasso','GECO')) {
  f  = if (meth == 'GECO') gc_file else gl_file
  gl = if (meth == 'GECO') 'alpha' else 'OP'
  for (net in c(gl, cts)) {
    d = as.data.frame(read_excel(f, sheet = net))
    reps = sub('^degree_', '', grep('^degree_rep', colnames(d), value = TRUE))
    sc = lapply(reps, function(r)
      rank(d[[paste0('degree_', r)]]) + rank(d[[paste0('strength_', r)]]) / (2 * nrow(d)))
    add(meth, ifelse(net == gl, 'OP', net), sc, d$gene)
  }
}
curves$net = ifelse(curves$net == 'OP', 'Bulk/Shared', curves$net)
au$net = ifelse(au$net == 'OP', 'Bulk/Shared', au$net)
write.table(curves, file.path(outdir, 'hub_roc_curves.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

write.table(au, file.path(outdir, 'hub_roc_curves_auroc.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
print(au, row.names = FALSE)
