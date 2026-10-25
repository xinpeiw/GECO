rm(list=ls())

outdir = 'results/2b_hubs_edges'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

cats = c('sharedOnly','both','specificOnly','neither')

summarise_method = function(file, global_label, ct_label) {
  d = read.table(file, header = TRUE)
  out = as.data.frame(rbind(as.numeric(d[d$net %in% c('Bulk','Shared'), cats]),
                            colMeans(d[!d$net %in% c('Bulk','Shared'), cats])))
  colnames(out) = cats
  out$row = c(global_label, ct_label)
  out
}

res = rbind(
  cbind(method = 'WGCNA',
        summarise_method(file.path(outdir, 'edge_composition_WGCNA.txt'), 'Bulk', 'Pseudobulk')),
  cbind(method = 'standard glasso',
        summarise_method(file.path(outdir, 'edge_composition_naive.txt'), 'Bulk', 'Pseudobulk')),
  cbind(method = 'GECO',
        summarise_method(file.path(outdir, 'edge_composition_geco.txt'), 'Shared', 'Specific')))

write.table(res, file.path(outdir, 'edge_composition_summary.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
