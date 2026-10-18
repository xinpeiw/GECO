rm(list=ls())

library(readxl)

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

wgcna_dir = 'data/2b_hubs_edges/WGCNA'
glasso_file = 'data/2b_hubs_edges/naive_glasso/degree_naive_glasso.xlsx'
geco_file = 'data/2b_hubs_edges/geco/degree_geco_merged_100reps.xlsx'
geco_dir = 'data/onek1k'

shared = readLines(file.path(geco_dir, 'GECO_Shared.txt'))
specific_all = readLines(file.path(geco_dir, 'GECO_Specific.txt'))
cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
Ktrue = c(OP=482, CD4NC=224, CD4ET=78, CD8NC=210, CD8ET=272, NK=297, BIN=122, BMem=96)
nm = c('sharedOnly','both','specificOnly','neither')

comp = function(s, sp) {
  n = length(s)
  both = sum(s %in% shared & s %in% sp)
  shO = sum(s %in% shared) - both
  spO = sum(s %in% sp) - both
  out = c(shO, both, spO, n - shO - both - spO) / n
  names(out) = nm
  out
}
sp_of = function(ct) {
  if (ct == 'OP') specific_all else readLines(file.path(geco_dir, paste0('hubs_', ct, '.txt')))
}
mkrow = function(method, net, ss) {
  x = as.data.frame(as.list(ss$m)); x$method = method; x$net = net
  x$valid = ss$valid; x$mean_size = ss$meansz
  x[, c('method','net',nm,'valid','mean_size')]
}
summarise_sets = function(sets, sp) {
  fr = NULL; sz = c()
  for (s in sets) {
    if (length(s) == 0) next
    fr = rbind(fr, comp(s, sp))
    sz = c(sz, length(s))
  }
  list(m = colMeans(fr), valid = nrow(fr), meansz = mean(sz))
}

pick_modprop = function(genes, mod, kme, K) {
  keep = !is.na(mod) & mod != 'grey' & !is.na(kme)
  genes = genes[keep]; mod = mod[keep]; kme = kme[keep]
  n = length(genes)
  if (n == 0) return(character(0))
  if (n <= K) return(genes)
  ms = table(mod)
  quota = K * as.numeric(ms) / n
  n_m = floor(quota); rem = quota - n_m
  short = K - sum(n_m)
  if (short > 0) n_m[order(-rem)[seq_len(short)]] = n_m[order(-rem)[seq_len(short)]] + 1
  names(n_m) = names(ms)
  unlist(lapply(names(ms), function(m) {
    idx = which(mod == m); idx = idx[order(-kme[idx])]
    genes[idx[seq_len(min(n_m[m], length(idx)))]]
  }))
}

res = NULL
for (ct in c('OP', cts)) {
  M = as.data.frame(read_excel(file.path(wgcna_dir, paste0(ct, '_Module_merged_100reps.xlsx'))))
  Km = as.data.frame(read_excel(file.path(wgcna_dir, paste0(ct, '_kME_merged_100reps.xlsx'))))
  sets = lapply(sprintf('rep%03d', 1:100), function(r)
    pick_modprop(M$Gene, M[[paste0('Module_', r)]], Km[[paste0('kME_', r)]], Ktrue[ct]))
  res = rbind(res, mkrow('WGCNA', ct, summarise_sets(sets, sp_of(ct))))
}
for (meth in c('glasso','GECO')) {
  f = if (meth == 'GECO') geco_file else glasso_file
  gl = if (meth == 'GECO') 'alpha' else 'OP'
  for (ct in c(gl, cts)) {
    d = as.data.frame(read_excel(f, sheet = ct))
    K = Ktrue[ifelse(ct == gl, 'OP', ct)]
    sp = sp_of(ifelse(ct == gl, 'OP', ct))
    reps = sub('^degree_', '', grep('^degree_rep', colnames(d), value = TRUE))
    sets = lapply(reps, function(r) {
      deg = d[[paste0('degree_', r)]]; str = d[[paste0('strength_', r)]]
      nz = which(deg > 0)
      ord = nz[order(-deg[nz], -str[nz])]
      d$gene[ord[seq_len(min(K, length(nz)))]]
    })
    res = rbind(res, mkrow(meth, ifelse(ct == gl, 'OP', ct), summarise_sets(sets, sp)))
  }
}
write.table(res, file.path(outdir, 'hub_composition_perct.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
