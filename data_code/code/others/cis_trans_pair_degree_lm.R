rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

moldir <- "data/6c-hub_position/cis_trans_pair"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(p) unique(trimws(readLines(p)))

GU <- rd("data/onek1k/GECO_universe.txt")
SH <- rd("data/onek1k/GECO_Shared.txt"); SP <- rd("data/onek1k/GECO_Specific.txt")
cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
cov <- merge(cov, fread("data/onek1k/beta_var.txt")[, .(gene, lbv=log(beta_var))], by="gene")

M_TOT <- "# colocalized genes"
M_IN  <- "# colocalized cis-eGenes (as trans-eGene)"
M_OUT <- "# colocalized trans-eGenes (as cis-eGene)"
mets  <- c(M_TOT, M_IN, M_OUT)

datasets <- list(
  list(f="onek1k.txt",  uni="onek1k_universe.txt",       lab="OneK1K",  stub="OneK1K"),
  list(f="eQTLGen.txt", uni="eQTLGen_universe_pair.txt", lab="eQTLGen", stub="eQTLGen"))

tab <- function(v, U) { t <- table(v); o <- as.integer(t[U]); o[is.na(o)] <- 0L; o }
fit <- function(deg, metric, hub, adj) {
  x <- deg[, .(gene, val=log1p(get(metric)), h=as.integer(gene %in% hub))]
  if (adj) { x <- merge(x, cov, by="gene"); x <- x[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]
             m <- lm(val ~ h + logE + logL + gc + lbv, x) } else m <- lm(val ~ h, x)
  co <- coef(summary(m))
  list(b=co["h","Estimate"], se=co["h","Std. Error"], p=co["h",4], n=nrow(x), nh=sum(x$h))
}

res <- rbindlist(lapply(datasets, function(s) {
  U <- intersect(rd(file.path(moldir, s$uni)), GU)
  d <- fread(file.path(moldir, s$f), header=TRUE); setnames(d, 1:2, c("cis","trans"))
  d <- unique(d[cis %in% U & trans %in% U & cis != trans, .(cis, trans)])
  und <- unique(data.table(a=pmin(d$cis,d$trans), b=pmax(d$cis,d$trans)))
  deg <- data.table(gene=U); deg[[M_TOT]] <- tab(c(und$a, und$b), U)
  deg[[M_IN]] <- tab(d$trans, U); deg[[M_OUT]] <- tab(d$cis, U)
  cat(sprintf("\n== %s | universe=%d | Shared=%d Specific=%d | directed pairs=%d ==\n",
              s$lab, length(U), length(intersect(SH,U)), length(intersect(SP,U)), nrow(d)))
  for (m in mets) cat(sprintf("   %-42s nonzero=%5d max=%4d | mean Shared=%.2f Specific=%.2f rest=%.2f\n", m,
      sum(deg[[m]]>0), max(deg[[m]]), mean(deg[[m]][deg$gene %in% SH]), mean(deg[[m]][deg$gene %in% SP]),
      mean(deg[[m]][!deg$gene %in% c(SH,SP)])))
  rbindlist(lapply(mets, function(mt) rbindlist(lapply(list(list("GECO Shared",SH),list("GECO Specific",SP)), function(g)
    rbindlist(lapply(c(FALSE,TRUE), function(a) { r <- fit(deg, mt, g[[2]], a)
      data.table(dataset=s$lab, metric=mt, group=g[[1]], mode=ifelse(a,"Adjusted","Baseline"),
                 n=r$n, n_hub=r$nh, coef=r$b, se=r$se, p=r$p) })))))) }))
res[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
fwrite(res[, .(dataset, metric, group, mode, n, n_hub, coef=round(coef,4), se=round(se,4), p, stars)],
       file.path(outdir, "cis_trans_pair_degree_lm.csv"))
cat("\n== lm(log1p(degree) ~ in_hub [+cov]) : coefficient (hub vs rest) ==\n")
print(res[, .(dataset, metric, group, mode, coef=round(coef,3), se=round(se,3), p=signif(p,2), stars)], nrows=Inf)
