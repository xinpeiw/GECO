rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

evo    <- "data/3cd-evolution"
hubdir <- "data/onek1k"
outdir <- "results/others"
CT  <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
rd  <- function(f) unique(trimws(readLines(f)))

U <- rd(file.path(hubdir, "GECO_universe.txt"))
methods <- list(
  GECO   = list(Bulk = rd(file.path(hubdir, "GECO_Shared.txt")),
                Pseudobulk = rd(file.path(hubdir, "GECO_Specific.txt"))),
  WGCNA  = list(Bulk = rd("data/others/WGCNA_hubs/Shared.txt"),
                Pseudobulk = unique(unlist(lapply(CT, function(ct)
                  rd(file.path("data/others/WGCNA_hubs", paste0(ct, ".txt"))))))),
  glasso = list(Bulk = rd("data/others/glasso_hubs/hubs_Shared.txt"),
                Pseudobulk = unique(unlist(lapply(CT, function(ct)
                  rd(file.path("data/others/glasso_hubs",
                               paste0("hubs_", ct, ".txt"))))))))
MTH  <- names(methods)
MLAB <- c(GECO="GECO", WGCNA="WGCNA", glasso="Standard glasso")

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
rc  <- fread(file.path(evo, "GECO_gene_recombination_rate.1000G.hg38.tsv"))[, .(gene = gene_id, recomb = recomb_rate_cM_Mb)]
ex  <- fread(file.path(evo, "exon.txt"))[, .(gene = ensg, nexon = n_exons_unique)]
bv  <- fread("data/onek1k/beta_var.txt")[, .(gene, beta_var)]
covdt <- cov[, .(gene, logE, logL, gc)]
covdt[rc, on = "gene", lrec := log1p(i.recomb)]
covdt[ex, on = "gene", lexn := log(i.nexon)]
covdt[bv, on = "gene", lbv  := log(i.beta_var)]
covdt <- covdt[is.finite(logE) & is.finite(logL) & is.finite(gc) &
               is.finite(lrec) & is.finite(lexn) & is.finite(lbv)]

vars <- list(list(file = "loeuf.txt",  col = "loeuf",  lab = "LOEUF"),
             list(file = "phylop.txt", col = "phylop", lab = "phyloP"),
             list(file = "B.txt",      col = "bvalue", lab = "B statistic"))
VLAB <- vapply(vars, function(v) v$lab, character(1))

MD <- lapply(vars, function(v) {
  m <- fread(file.path(evo, v$file)); setnames(m, 1, "gene")
  m[gene %in% U, .(gene, val = as.numeric(get(v$col)))][is.finite(val)]
})
names(MD) <- VLAB

grid <- CJ(vi = VLAB, mth = MTH, gp = c("Bulk","Pseudobulk"),
           adj = c(FALSE, TRUE), sorted = FALSE)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  d <- copy(MD[[grid$vi[i]]])
  d[, x := as.integer(gene %in% methods[[grid$mth[i]]][[grid$gp[i]]])]
  if (grid$adj[i]) {
    d <- merge(d, covdt, by = "gene")
    mm <- lm(val ~ x + logE + logL + gc + lexn + lrec + lbv, d)
  } else mm <- lm(val ~ x, d)
  co <- coef(summary(mm))["x", ]
  data.table(variable = grid$vi[i], method = grid$mth[i], group = grid$gp[i],
             mode = ifelse(grid$adj[i], "Adjusted", "Baseline"),
             n = nrow(d), n_hub = sum(d$x), beta = co[1], se = co[2], p = co[4])
}))
res[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
fwrite(res[order(variable, method, group, mode)],
       file.path(outdir, "evolution_continuous_lm_3methods.csv"))

cat("hub set sizes in the GECO universe (each metric further restricts to genes carrying it):\n")
