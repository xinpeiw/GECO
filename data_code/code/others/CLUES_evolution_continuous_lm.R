rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

COND   <- "mega"
COVTAG <- "all"
DLAB   <- "all"

evo    <- "data/others/evolution_clues"
covdir <- "data/clues"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(trimws(readLines(f)))

U   <- rd("data/clues/allgenes_CLUES.txt")
cov <- fread(file.path(covdir, paste0("covariates_", COVTAG, ".txt")))[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
bv  <- fread(file.path(covdir, paste0("beta_var_",   COVTAG, ".txt")))[, .(gene, logBV=log(beta_var))]
cov <- merge(cov, bv, by="gene")
rc  <- fread(file.path(evo, "CLUES_gene_recombination_rate.1000G.hg38.tsv"))[, .(gene=gene_id, recomb=recomb_rate_cM_Mb)]
ex  <- fread(file.path(evo, "exon.txt"))[, .(gene=ensg, nexon=n_exons_unique)]
cov[rc, on="gene", lrec := log1p(i.recomb)]
cov[ex, on="gene", lexn := log(i.nexon)]

vars <- list(list(file="loeuf.txt",  col="loeuf",  lab="LOEUF"),
             list(file="phylop.txt", col="phylop", lab="phyloP"),
             list(file="B.txt",      col="bvalue", lab="B statistic"))

SH_all <- intersect(rd("data/clues/hubs_alpha.txt"), U)
SP_all <- intersect(rd("data/clues/hubs_specific.txt"), U)

fit1 <- function(hubG, md, adj) {
  d <- copy(md); d[, x := as.integer(gene %in% hubG)]
  if (adj) { d <- merge(d, cov, by="gene"); d <- d[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lexn)&is.finite(lrec)&is.finite(logBV)]
             mm <- lm(val ~ x + logE + logL + gc + lexn + lrec + logBV, d) }
  else     { mm <- lm(val ~ x, d) }
  co <- coef(summary(mm)); if (!("x" %in% rownames(co))) return(list(beta=NA, se=NA, p=NA, n=nrow(d), n_hub=sum(d$x)))
  list(beta=co["x","Estimate"], se=co["x","Std. Error"], p=co["x","Pr(>|t|)"], n=nrow(d), n_hub=sum(d$x))
}

res <- rbindlist(lapply(vars, function(v) {
  m <- fread(file.path(evo, v$file)); setnames(m, 1, "gene")
  md <- m[gene %in% U, .(gene, val=as.numeric(get(v$col)))][is.finite(val)]
  rbindlist(lapply(list(list("CLUES Shared", SH_all), list("CLUES Specific", SP_all)), function(g)
    rbindlist(lapply(c(FALSE, TRUE), function(adj) {
      r <- fit1(g[[2]], md, adj)
      data.table(dataset=DLAB, variable=v$lab, group=g[[1]], mode=ifelse(adj,"Adjusted","Baseline"),
                 n=r$n, n_hub=r$n_hub, beta=r$beta, se=r$se, p=r$p)
    }))))
}))
res[, stars := fifelse(is.na(p),"",fifelse(p<.001,"***",fifelse(p<.01,"**",fifelse(p<.05,"*",""))))]
fwrite(res[order(variable,group,mode)], file.path(outdir, "CLUES_evolution_continuous_lm.csv"))
cat(sprintf("\n===== CLUES continuous evolution (%s / %s) : lm coef, hub vs rest =====\n", DLAB, COND))
print(res[order(variable,group,mode), .(variable, group, mode, n_hub, beta=round(beta,3), se=round(se,3), p=signif(p,2), stars)])
