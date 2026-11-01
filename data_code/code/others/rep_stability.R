rm(list=ls())

library(readxl)

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

wgcna_dir   = 'data/2b_hubs_edges/WGCNA'
glasso_file = 'data/2b_hubs_edges/naive_glasso/degree_naive_glasso.xlsx'
geco_file   = 'data/2b_hubs_edges/geco/degree_geco_merged_100reps.xlsx'

cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
Ktrue = c(OP=482, CD4NC=224, CD4ET=78, CD8NC=210, CD8ET=272, NK=297, BIN=122, BMem=96)

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

jac = function(sets) {
  g = sort(unique(unlist(sets)))
  M = sapply(sets, function(s) as.integer(g %in% s))
  if (is.null(dim(M))) M = matrix(M, nrow = length(g))
  I = crossprod(M); sz = diag(I); U = outer(sz, sz, '+') - I
  J = ifelse(U > 0, I / U, 0)
  up = upper.tri(J)
  ne = sz > 0
  keep = outer(ne, ne, '&') & up
  data.frame(n_rep = length(sets), n_nonempty = sum(ne),
             mean_size = mean(sz[ne]),
             jaccard = mean(J[up]),
             jaccard_nonempty = if (any(keep)) mean(J[keep]) else NA_real_)
}

res = NULL

for (ct in c('OP', cts)) {
  M  = as.data.frame(read_excel(file.path(wgcna_dir, paste0(ct, '_Module_merged_100reps.xlsx'))))
  Km = as.data.frame(read_excel(file.path(wgcna_dir, paste0(ct, '_kME_merged_100reps.xlsx'))))
  sets = lapply(sprintf('rep%03d', 1:100), function(r)
    pick_modprop(M$Gene, M[[paste0('Module_', r)]], Km[[paste0('kME_', r)]], Ktrue[ct]))
  res = rbind(res, cbind(data.frame(method = 'WGCNA', net = ct, K = Ktrue[ct]), jac(sets)))
}

for (meth in c('Standard glasso','GECO')) {
  f  = if (meth == 'GECO') geco_file else glasso_file
  gl = if (meth == 'GECO') 'alpha' else 'OP'
  for (ct in c(gl, cts)) {
    d = as.data.frame(read_excel(f, sheet = ct))
    key = ifelse(ct == gl, 'OP', ct); K = Ktrue[key]
    reps = sub('^degree_', '', grep('^degree_rep', colnames(d), value = TRUE))
    sets = lapply(reps, function(r) {
      deg = d[[paste0('degree_', r)]]; str = d[[paste0('strength_', r)]]
      nz = which(deg > 0)
      ord = nz[order(-deg[nz], -str[nz])]
      d$gene[ord[seq_len(min(K, length(nz)))]]
    })
    res = rbind(res, cbind(data.frame(method = meth, net = key, K = K), jac(sets)))
  }
}
rownames(res) = NULL
res$net = ifelse(res$net == 'OP', 'Bulk/Shared', res$net)
write.table(res, file.path(outdir, 'rep_stability.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
