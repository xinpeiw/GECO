rm(list=ls())
suppressMessages(library(data.table))

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

bv <- fread("data/onek1k/beta_var.txt")
bv <- bv[order(bv$beta_var, decreasing = TRUE), ]

cts  <- c('Shared', 'CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK', 'BIN', 'BMem')
hubf <- c(Shared = "GECO_Shared.txt", setNames(paste0("hubs_", cts[-1], ".txt"), cts[-1]))

res <- rbindlist(lapply(cts, function(ct) {
  hubs <- fread(file.path("data/onek1k", hubf[[ct]]), header = FALSE)$V1
  topv <- bv$gene[seq_len(length(hubs))]
  data.table(network = ct, n_hubs = length(hubs), n_top_var = length(topv),
             hub_only = length(setdiff(hubs, topv)), overlap = length(intersect(hubs, topv)),
             top_var_only = length(setdiff(topv, hubs)), union = length(union(hubs, topv)))
}))
res[, `:=`(prop_hub_only = hub_only / union, prop_overlap = overlap / union, prop_top_var_only = top_var_only / union)]
write.table(res, file.path(outdir, "barplot_varbeta.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
print(res, row.names = FALSE)
