rm(list=ls())


edge_dir = 'data/2b_hubs_edges/edge_top100'
outdir = 'results/2b_hubs_edges'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

true = read.table(file.path(edge_dir, 'edge_true100.txt'), header = TRUE)
true_shared = true$alpha[!is.na(true$alpha)]
cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
true_spec_union = unique(na.omit(unlist(true[, cts])))

methods = list(
  WGCNA = list(dir = 'WGCNA', global = 'OP', global_label = 'Bulk'),
  geco = list(dir = 'geco', global = 'alpha', global_label = 'Shared'),
  naive = list(dir = 'naive_glasso', global = 'OP', global_label = 'Bulk'))

for (m in names(methods)) {
  info = methods[[m]]
  res = NULL
  for (net in c(info$global, cts)) {
    sp = if (net == info$global) true_spec_union else true[[net]][!is.na(true[[net]])]
    d = read.table(file.path(edge_dir, info$dir, paste0('edge_100_', net, '.txt')),
                   header = TRUE)
    fr = NULL; sz = c()
    for (j in 1:ncol(d)) {
      s = d[[j]]
      s = s[!is.na(s)]
      if (length(s) == 0) next
      both = sum(s %in% true_shared & s %in% sp)
      shO = sum(s %in% true_shared) - both
      spO = sum(s %in% sp) - both
      nei = length(s) - both - shO - spO
      fr = rbind(fr, c(shO, both, spO, nei) / length(s))
      sz = c(sz, length(s))
    }
    mm = colMeans(fr)
    res = rbind(res, data.frame(net = ifelse(net == info$global, info$global_label, net),
                                valid = nrow(fr), mean_size = mean(sz),
                                sharedOnly = mm[1], both = mm[2],
                                specificOnly = mm[3], neither = mm[4]))
  }
  write.table(res, file.path(outdir, paste0('edge_composition_', m, '.txt')),
              row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
}
