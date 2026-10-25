rm(list=ls())
suppressMessages(library(igraph))

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

pert <- readRDS("data/others/pertseq_union_network_matrix.rds")
geco <- readRDS("data/others/geco_nk_network_matrix.rds")
stopifnot(identical(rownames(pert), rownames(geco)))
s <- read.delim("data/others/geco_nk_network.tsv")
id2sym <- unique(rbind(setNames(s[, c(1,3)], c("e","s")), setNames(s[, c(2,4)], c("e","s"))))
map <- setNames(id2sym$s, id2sym$e)

ov  <- geco * pert
g   <- rownames(ov)
idx <- which(ov == 1 & upper.tri(ov), arr.ind=TRUE)
gr  <- graph_from_edgelist(cbind(g[idx[,1]], g[idx[,2]]), directed=FALSE)
comp <- components(gr)
gcc  <- induced_subgraph(gr, V(gr)[comp$membership == which.max(comp$csize)])
set.seed(1); cl <- cluster_louvain(gcc); mem <- membership(cl)
sym <- map[V(gcc)$name]; sym[is.na(sym)] <- V(gcc)$name[is.na(sym)]
deg <- degree(gcc)

nodes <- data.frame(gene = V(gcc)$name, symbol = unname(sym), module = as.integer(mem), degree = unname(deg))
nodes <- nodes[order(nodes$module, -nodes$degree), ]
el <- as_edgelist(gcc, names=TRUE)
edges <- data.frame(gene_a = el[,1], gene_b = el[,2], symbol_a = unname(map[el[,1]]), symbol_b = unname(map[el[,2]]))
write.table(nodes, file.path(outdir, "nk_overlap_module_map_nodes.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(edges, file.path(outdir, "nk_overlap_module_map_edges.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

cat(sprintf("overlap network: %d genes %d edges | largest component: %d genes %d edges | %d modules, modularity %.2f\n",
            gorder(gr), gsize(gr), gorder(gcc), gsize(gcc), length(cl), modularity(cl)))
sz <- sizes(cl)
for (mm in order(sz, decreasing=TRUE)) {
  vs <- which(mem == mm); topg <- vs[order(deg[vs], decreasing=TRUE)][1:min(4, length(vs))]
  cat(sprintf("m%d n=%d: %s\n", mm, sz[mm], paste(sym[topg], collapse=",")))
}
