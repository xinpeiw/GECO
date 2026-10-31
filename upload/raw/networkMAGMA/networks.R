CT_GECO  <- c("CD4NC", "CD4ET", "CD8NC", "CD8ET", "NK", "BIN", "BMem")

VARIANTS <- list(
  geco          = c("alpha", "CD4NC", "CD4ET", "CD8NC", "CD8ET", "NK", "BIN", "BMem"),
  clues         = c("alpha", "T4", "T8", "NK", "B", "cM", "ncM", "cDC"),
  perturbseq    = c("geco_shared", "geco_specific_union", "pertseq_specific_overlap", "pertseq_shared_overlap", "pertseq_union"),
  perturbseq_ct = c("pertseq_cd4", "pertseq_cd8", "pertseq_nk", "pertseq_b"),
  glasso01      = c("glasso_Bulk", "glasso_Specific", "GECO_Shared", "GECO_Specific"))

partial_cor <- function(omega) { r <- -omega / sqrt(outer(diag(omega), diag(omega))); diag(r) <- 0; r }
binary      <- function(m) { diag(m) <- 0; (m != 0) + 0L }
union_bin   <- function(mats) { a <- binary(mats[[1]]); for (m in mats[-1]) a <- a | binary(m); a + 0L }

geco_precision   <- function(ct) readRDS(file.path("data/onek1k_precision", paste0(ct, ".rds")))
clues_precision  <- function(ct) readRDS(file.path("data/others/precision_clues", paste0("precision_", ct, ".rds")))
glasso_precision <- function(ct) {
  m <- readRDS(file.path("data/others/glasso_precision", paste0("precision_", ct, ".rds")))
  g <- readLines("data/others/glasso_precision/gene_order_exp.txt"); dimnames(m) <- list(g, g); m }

build_network <- function(variant, name) {
  if (variant == "geco")  return(partial_cor(geco_precision(name)))
  if (variant == "clues") return(partial_cor(clues_precision(name)))
  if (variant == "glasso01") {
    return(switch(name,
      glasso_Bulk   = binary(glasso_precision("Bulk")),
      glasso_Specific = union_bin(lapply(CT_GECO, glasso_precision)),
      GECO_Shared     = binary(geco_precision("alpha")),
      GECO_Specific   = union_bin(lapply(CT_GECO, geco_precision))))
  }
  pert <- readRDS("data/others/pertseq_union_network_matrix.rds"); g <- rownames(pert)
  sub  <- function(ct) geco_precision(ct)[g, g]
  switch(name,
    pertseq_union            = pert,
    geco_shared              = binary(sub("alpha")),
    geco_specific_union      = union_bin(lapply(CT_GECO, sub)),
    pertseq_shared_overlap   = binary(sub("alpha")) * pert,
    pertseq_specific_overlap = union_bin(lapply(CT_GECO, sub)) * pert,
    pertseq_cd4              = union_bin(lapply(c("CD4NC", "CD4ET"), sub)) * pert,
    pertseq_cd8              = union_bin(lapply(c("CD8NC", "CD8ET"), sub)) * pert,
    pertseq_nk               = binary(sub("NK")) * pert,
    pertseq_b                = union_bin(lapply(c("BIN", "BMem"), sub)) * pert)
}
