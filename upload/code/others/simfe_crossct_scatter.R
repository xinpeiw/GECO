rm(list=ls())

library(readxl)

outdir = 'results/others'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

wgcna_dir   = 'data/2b_hubs_edges/WGCNA'
glasso_file = 'data/2b_hubs_edges/naive_glasso/degree_naive_glasso.xlsx'
geco_file   = 'data/2b_hubs_edges/geco/degree_geco_merged_100reps.xlsx'
geco_dir    = 'data/onek1k'

N = 11526
cts  = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
nets = c('OP', cts)
Ktrue = c(OP=482, CD4NC=224, CD4ET=78, CD8NC=210, CD8ET=272, NK=297, BIN=122, BMem=96)
lab_of = function(net) ifelse(net == 'OP', 'Bulk/Shared', net)

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

load_sets = function(meth, net) {
  if (meth == 'WGCNA') {
    M  = as.data.frame(read_excel(file.path(wgcna_dir, paste0(net, '_Module_merged_100reps.xlsx'))))
    Km = as.data.frame(read_excel(file.path(wgcna_dir, paste0(net, '_kME_merged_100reps.xlsx'))))
    lapply(sprintf('rep%03d', 1:100), function(r)
      pick_modprop(M$Gene, M[[paste0('Module_', r)]], Km[[paste0('kME_', r)]], Ktrue[net]))
  } else {
    f  = if (meth == 'GECO') geco_file else glasso_file
    sh = if (meth == 'GECO' && net == 'OP') 'alpha' else net
    d = as.data.frame(read_excel(f, sheet = sh))
    reps = sub('^degree_', '', grep('^degree_rep', colnames(d), value = TRUE))
    lapply(reps, function(r) {
      deg = d[[paste0('degree_', r)]]; str = d[[paste0('strength_', r)]]
      nz = which(deg > 0)
      ord = nz[order(-deg[nz], -str[nz])]
      d$gene[ord[seq_len(min(Ktrue[net], length(nz)))]]
    })
  }
}

rd = function(p) { x = readLines(p); x[nzchar(x)] }
truth_sets = c(list(OP = rd(file.path(geco_dir, 'GECO_Shared.txt'))),
               setNames(lapply(cts, function(ct)
                 rd(file.path(geco_dir, paste0('hubs_', ct, '.txt')))), cts))

pairs = t(combn(nets, 2))
methods = c('WGCNA','Standard glasso','GECO')

res = NULL
for (k in seq_len(nrow(pairs))) {
  a = truth_sets[[pairs[k, 1]]]; b = truth_sets[[pairs[k, 2]]]
  ov = length(intersect(a, b))
  res = rbind(res, data.frame(net1 = pairs[k, 1], net2 = pairs[k, 2],
                              series = 'Truth',
                              jac = ov / length(union(a, b)),
                              fe = ov * N / (length(a) * length(b)),
                              valid_reps = NA_integer_))
}
for (meth in methods) {
  sets = setNames(lapply(nets, function(net) load_sets(meth, net)), nets)
  for (k in seq_len(nrow(pairs))) {
    A = sets[[pairs[k, 1]]]; B = sets[[pairs[k, 2]]]
    ok = vapply(seq_along(A), function(r)
      length(A[[r]]) > 0 && length(B[[r]]) > 0, logical(1))
    js = vapply(seq_along(A), function(r) {
      if (!ok[r]) return(NA_real_)
      length(intersect(A[[r]], B[[r]])) / length(union(A[[r]], B[[r]]))
    }, numeric(1))
    fes = vapply(seq_along(A), function(r) {
      if (!ok[r]) return(NA_real_)
      length(intersect(A[[r]], B[[r]])) * N / (length(A[[r]]) * length(B[[r]]))
    }, numeric(1))
    res = rbind(res, data.frame(net1 = pairs[k, 1], net2 = pairs[k, 2],
                                series = meth,
                                jac = if (all(is.na(js)))  NA_real_ else mean(js,  na.rm = TRUE),
                                fe  = if (all(is.na(fes))) NA_real_ else mean(fes, na.rm = TRUE),
                                valid_reps = sum(ok)))
  }
}
res$pair = paste(lab_of(res$net1), lab_of(res$net2), sep = '-')
write.table(res[, c('series','net1','net2','pair','jac','fe','valid_reps')],
            file.path(outdir, 'simfe_crossct_scatter.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

sc = merge(res[res$series != 'Truth', c('pair','series','jac','valid_reps')],
           setNames(res[res$series == 'Truth', c('pair','jac')], c('pair','truth')), by = 'pair')
sc = sc[is.finite(sc$jac) & sc$jac > 0 & is.finite(sc$truth), ]
ann = do.call(rbind, lapply(methods, function(m) {
  z = sc[sc$series == m, ]
  data.frame(series = m, n_pairs = nrow(z),
             spearman_rho = cor(z$jac, z$truth, method = 'spearman'),
             median_est_over_truth = median(z$jac / z$truth))
}))
write.table(ann, file.path(outdir, 'simfe_crossct_scatter_stats.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
print(ann, row.names = FALSE)
