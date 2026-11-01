rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

f      <- "data/6c-hub_position/perturb_seq/aggregate_per_intersection_gene_metrics.tsv"
hubdir <- "data/onek1k"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(p) unique(trimws(readLines(p)))

d <- fread(f)[, .(gene=gene_ensg, n_incoming_neighbors, n_outgoing_neighbors, n_neighbors)]
U <- d$gene
SHARED <- intersect(rd(file.path(hubdir,"GECO_Shared.txt")),   U)
SPEC   <- intersect(rd(file.path(hubdir,"GECO_Specific.txt")), U)
CD4    <- intersect(intersect(union(rd(file.path(hubdir,"hubs_CD4NC.txt")), rd(file.path(hubdir,"hubs_CD4ET.txt"))), SPEC), U)
NONCD4 <- setdiff(SPEC, CD4)
grps   <- list("Shared"=SHARED, "CD4-specific"=CD4, "non-CD4 specific"=NONCD4)
cat(sprintf("universe=%d | Shared=%d | CD4-specific=%d | non-CD4 specific=%d (Specific=%d)\n",
            length(U), length(SHARED), length(CD4), length(NONCD4), length(SPEC)))

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
cov[fread("data/onek1k/beta_var.txt"), on="gene", lbv:=log(i.beta_var)]
cov <- cov[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]

fit <- function(md, grp, adj) {
  dd <- copy(md); dd[, x := as.integer(gene %in% grp)]
  m  <- if (adj) { dd <- merge(dd, cov, by="gene"); lm(val ~ x + logE + logL + gc + lbv, dd) } else lm(val ~ x, dd)
  co <- coef(summary(m))["x", ]; list(b=co[1], se=co[2], p=co[4], n=nrow(dd), nh=sum(dd$x))
}

metrics <- list(
  list(col="n_incoming_neighbors", lab="# regulators", stub="n_incoming_neighbors"),
  list(col="n_outgoing_neighbors", lab="# responders", stub="n_outgoing_neighbors"),
  list(col="n_neighbors",          lab="# neighbors",  stub="n_neighbors"))

all <- list()
for (mt in metrics) {
  md  <- d[, .(gene, val=log1p(get(mt$col)))]
  all[[mt$col]] <- rbindlist(lapply(names(grps), function(gn) rbindlist(lapply(c(FALSE,TRUE), function(a) {
    r <- fit(md, grps[[gn]], a)
    data.table(metric=mt$lab, col=mt$col, group=gn, mode=ifelse(a,"Adjusted","Baseline"),
               n=r$n, n_hub=r$nh, coef=r$b, se=r$se, p=r$p)
  }))))
}
res_all <- rbindlist(all)
res_all[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
fwrite(res_all[, .(metric, col, group, mode, n, n_hub, coef=round(coef,4), se=round(se,4), p, stars)],
       file.path(outdir, "perturbation_neighbors_lm.csv"))
cat("\n== lm(log1p(metric) ~ in_group [+cov]) : coefficient vs rest ==\n")
print(res_all[, .(metric, group, mode, coef=round(coef,3), se=round(se,3), p=signif(p,2), stars)])
