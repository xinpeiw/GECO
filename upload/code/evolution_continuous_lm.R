rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

evo    <- "data/3cd-evolution"
hubdir <- "data/onek1k"
outdir <- "results/3cd-evolution"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) unique(trimws(readLines(f)))
U        <- rd(file.path(hubdir, "GECO_universe.txt"))
SHARED   <- intersect(rd(file.path(hubdir, "GECO_Shared.txt")),   U)
SPECIFIC <- intersect(rd(file.path(hubdir, "GECO_Specific.txt")), U)

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
rc  <- fread(file.path(evo, "GECO_gene_recombination_rate.1000G.hg38.tsv"))[, .(gene = gene_id, recomb = recomb_rate_cM_Mb)]
ex  <- fread(file.path(evo, "exon.txt"))[, .(gene = ensg, nexon = n_exons_unique)]
bv  <- fread("data/onek1k/beta_var.txt")[, .(gene, beta_var)]
covdt <- cov[, .(gene, logE, logL, gc)]
covdt[rc, on = "gene", lrec := log1p(i.recomb)]
covdt[ex, on = "gene", lexn := log(i.nexon)]
covdt[bv, on = "gene", lbv := log(i.beta_var)]
covdt <- covdt[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(lrec) & is.finite(lexn) & is.finite(lbv)]

vars <- list(list(file = "loeuf.txt",  col = "loeuf",  lab = "LOEUF"),
             list(file = "phylop.txt", col = "phylop", lab = "phyloP"),
             list(file = "B.txt",      col = "bvalue", lab = "B statistic"))

res <- rbindlist(lapply(vars, function(v) {
  m <- fread(file.path(evo, v$file)); setnames(m, 1, "gene")
  md <- m[gene %in% U, .(gene, val = as.numeric(get(v$col)))][is.finite(val)]
  rbindlist(lapply(list(list("GECO Shared", SHARED), list("GECO Specific", SPECIFIC)), function(g)
    rbindlist(lapply(c(FALSE, TRUE), function(adj) {
      d <- copy(md); d[, x := as.integer(gene %in% g[[2]])]
      if (adj) { d <- merge(d, covdt, by = "gene"); mm <- lm(val ~ x + logE + logL + gc + lexn + lrec + lbv, d) }
      else     { mm <- lm(val ~ x, d) }
      co <- coef(summary(mm))["x", ]
      data.table(variable = v$lab, group = g[[1]], mode = ifelse(adj, "Adjusted", "Baseline"),
                 n = nrow(d), n_hub = sum(d$x), beta = co[1], se = co[2], p = co[4])
    }))))
}))
res[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
fwrite(res[order(variable, group, mode)], file.path(outdir, "evolution_continuous_lm.csv"))
cat("\n===== continuous evolutionary features (linear regression coefficient, hub vs non-hub) =====\n")
print(res[order(variable, group, mode), .(variable, group, mode, n, n_hub, beta = round(beta,3), se = round(se,3), p = signif(p,2), stars)])
