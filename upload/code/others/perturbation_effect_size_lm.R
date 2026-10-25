rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

OFF <- 1

f      <- "data/6c-hub_position/perturb_seq/aggregate_per_intersection_gene_metrics.tsv"
hubdir <- "data/onek1k"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(p) unique(trimws(readLines(p)))

SHARED_all <- rd(file.path(hubdir,"GECO_Shared.txt"))
SPEC_all   <- rd(file.path(hubdir,"GECO_Specific.txt"))
CD4_all    <- union(rd(file.path(hubdir,"hubs_CD4NC.txt")), rd(file.path(hubdir,"hubs_CD4ET.txt")))
get_grps <- function(U) {
  SP <- intersect(SPEC_all, U); CD4 <- intersect(intersect(CD4_all, SP), U)
  list("Shared"=intersect(SHARED_all, U), "CD4-specific"=CD4, "non-CD4 specific"=setdiff(SP, CD4))
}
cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
cov[fread("data/onek1k/beta_var.txt"), on="gene", lbv:=log(i.beta_var)]
cov <- cov[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]
fit <- function(md, grp, adj) {
  dd <- copy(md); dd[, x := as.integer(gene %in% grp)]
  m  <- if (adj) { dd <- merge(dd, cov, by="gene"); lm(val ~ x + logE + logL + gc + lbv, dd) } else lm(val ~ x, dd)
  co <- coef(summary(m))["x", ]; list(b=co[1], se=co[2], p=co[4], n=nrow(dd), nh=sum(dd$x))
}

metrics <- list(
  list(col="mean_normed_effect_sig_pair_in",  lab="Mean incoming effect (received)"),
  list(col="mean_normed_effect_sig_pair_out", lab="Mean outgoing effect (exerted)"))
dat <- fread(f, select=c("gene_ensg", vapply(metrics, function(m) m$col, character(1))))

res <- rbindlist(lapply(metrics, function(mt) {

  md <- dat[, .(gene=gene_ensg, val=log(get(mt$col) + OFF))][is.finite(val)]
  grps <- get_grps(md$gene)
  cat(sprintf("%-32s U=%d | Shared=%d CD4=%d nonCD4=%d\n", mt$col, nrow(md),
              length(grps$Shared), length(grps[["CD4-specific"]]), length(grps[["non-CD4 specific"]])))
  rbindlist(lapply(names(grps), function(gn) rbindlist(lapply(c(FALSE,TRUE), function(a) {
    r <- fit(md, grps[[gn]], a)
    data.table(metric=mt$lab, col=mt$col, group=gn, mode=ifelse(a,"Adjusted","Baseline"), n=r$n, n_hub=r$nh, coef=r$b, se=r$se, p=r$p)
  }))))
}))
res[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
fwrite(res[, .(metric, col, group, mode, n, n_hub, coef=round(coef,4), se=round(se,4), p, stars)],
       file.path(outdir, "perturbation_effect_size_lm.csv"))
cat("\n== lm(log(effect, sig_pair) ~ in_group [+cov]) : coefficient vs rest ==\n")
print(res[order(metric,group,mode), .(metric, group, mode, coef=round(coef,3), se=round(se,3), p=signif(p,2), stars)])
