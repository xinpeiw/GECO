rm(list=ls())

library(readxl)

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

wgcna_dir   = 'data/2b_hubs_edges/WGCNA'
glasso_file = 'data/2b_hubs_edges/naive_glasso/degree_naive_glasso.xlsx'
geco_file   = 'data/2b_hubs_edges/geco/degree_geco_merged_100reps.xlsx'

cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
TOPS = 100

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

jac = function(sets, TOP) {
  g = sort(unique(unlist(sets)))
  M = sapply(sets, function(s) as.integer(g %in% s))
  if (is.null(dim(M))) M = matrix(M, nrow = length(g))
  I = crossprod(M); sz = diag(I); U = outer(sz, sz, '+') - I
  J = ifelse(U > 0, I / U, 0)
  up = upper.tri(J); ne = sz > 0
  keep = outer(ne, ne, '&') & up
  data.frame(n_rep = length(sets), n_nonempty = sum(ne), n_full = sum(sz == TOP),
             mean_size = mean(sz[ne]),
             jaccard = mean(J[up]),
             jaccard_nonempty = if (any(keep)) mean(J[keep]) else NA_real_)
}

wg = lapply(c('OP', cts), function(ct) list(
  M  = as.data.frame(read_excel(file.path(wgcna_dir, paste0(ct, '_Module_merged_100reps.xlsx')))),
  Km = as.data.frame(read_excel(file.path(wgcna_dir, paste0(ct, '_kME_merged_100reps.xlsx'))))))
names(wg) = c('OP', cts)
dg = list('Standard glasso' = lapply(c('OP', cts), function(s)
            as.data.frame(read_excel(glasso_file, sheet = s))),
          'GECO' = lapply(c('alpha', cts), function(s)
            as.data.frame(read_excel(geco_file, sheet = s))))
for (m in names(dg)) names(dg[[m]]) = c('OP', cts)

for (TOP in TOPS) {
  res = NULL
  for (ct in c('OP', cts)) {
    x = wg[[ct]]
    sets = lapply(sprintf('rep%03d', 1:100), function(r)
      pick_modprop(x$M$Gene, x$M[[paste0('Module_', r)]],
                   x$Km[[paste0('kME_', r)]], TOP))
    res = rbind(res, cbind(data.frame(method = 'WGCNA', net = ct), jac(sets, TOP)))
  }
  for (meth in c('Standard glasso','GECO')) {
    for (ct in c('OP', cts)) {
      d = dg[[meth]][[ct]]
      reps = sub('^degree_', '', grep('^degree_rep', colnames(d), value = TRUE))
      sets = lapply(reps, function(r) {
        deg = d[[paste0('degree_', r)]]; str = d[[paste0('strength_', r)]]
        nz = which(deg > 0)
        ord = nz[order(-deg[nz], -str[nz])]
        d$gene[ord[seq_len(min(TOP, length(nz)))]]
      })
      res = rbind(res, cbind(data.frame(method = meth, net = ct), jac(sets, TOP)))
    }
  }
  rownames(res) = NULL
  res$net = ifelse(res$net == 'OP', 'Bulk/Shared', res$net)
  write.table(res, file.path(outdir, sprintf('rep_stability_top%d.txt', TOP)),
              row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
}
