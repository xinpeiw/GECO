rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

COND   <- "mega"
COVTAG <- "all"
DLAB   <- "all"

padir  <- "data/7c-Mendelian_gene"
covdir <- "data/clues"
outdir <- "results/7c-Mendelian_gene"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(trimws(readLines(f)))

U   <- rd("data/clues/allgenes_CLUES.txt")
cov <- fread(file.path(covdir, paste0("covariates_", COVTAG, ".txt")))[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
bv  <- fread(file.path(covdir, paste0("beta_var_",   COVTAG, ".txt")))[, .(gene, logBV=log(beta_var))]
cov <- merge(cov, bv, by="gene")

SH <- intersect(rd("data/clues/hubs_alpha.txt"),    U)
SP <- intersect(rd("data/clues/hubs_specific.txt"), U)

sets <- list(
  Full = list(imm = rd(file.path(padir, "mendelian_immune_full.txt")),
              non = rd(file.path(padir, "mendelian_nonimmune_full.txt"))))

fit1 <- function(hubG, setG, adj) {
  d <- data.table(gene=U, in_hub=as.integer(U %in% hubG), in_set=as.integer(U %in% setG))
  ov <- sum(d$in_hub & d$in_set)
  if (adj) {
    d <- merge(d, cov, by="gene"); d <- d[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(logBV)]
    a2 <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set)
    naret <- list(l2=NA, se=NA, p=NA, ov=a2)
    if (a2==0 || b2==0) return(naret)
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + logBV, d, family=binomial)), error=function(e) NULL)
    if (is.null(m) || !m$converged) return(naret)
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(naret)
    se <- co["in_hub","Std. Error"]; if (!is.finite(se) || se > 10) return(naret)
    return(list(l2=co["in_hub","Estimate"]/log(2), se=se/log(2), p=co["in_hub","Pr(>|z|)"], ov=a2))
  } else {
    a<-sum(d$in_hub&d$in_set); b<-sum(d$in_hub&!d$in_set); c<-sum(!d$in_hub&d$in_set); e<-sum(!d$in_hub&!d$in_set)
    if (a==0||b==0||c==0||e==0) { a<-a+.5;b<-b+.5;c<-c+.5;e<-e+.5 }
    l2 <- log((a*e)/(b*c))/log(2); se <- sqrt(1/a+1/b+1/c+1/e)/log(2)
    return(list(l2=l2, se=se, p=2*pnorm(-abs(l2/se)), ov=ov))
  }
}

res <- rbindlist(lapply(names(sets), function(v) {
  datasets <- list(list(G = sets[[v]]$imm, lab = "Immune Mendelian genes"),
                   list(G = sets[[v]]$non, lab = "Non-immune Mendelian genes"))
  rbindlist(lapply(datasets, function(s) {
    G <- intersect(s$G, U)
    rbindlist(lapply(list(list("CLUES Shared", SH), list("CLUES Specific", SP)), function(g)
      rbindlist(lapply(c(FALSE, TRUE), function(adj) {
        r <- fit1(g[[2]], G, adj)
        data.table(variant=v, dataset=DLAB, geneset=s$lab, nU=length(U), nSet=length(G), nOv=r$ov,
                   group=g[[1]], mode=ifelse(adj,"Adjusted","Baseline"), l2=r$l2, se=r$se, p=r$p)
      }))))
  }))
}))
res[, OR := 2^l2]; res[, stars := fifelse(is.na(p),"",fifelse(p<.001,"***",fifelse(p<.01,"**",fifelse(p<.05,"*",""))))]
fwrite(res[order(variant,geneset,group,mode)], file.path(outdir, "CLUES_mendelian_panelapp_enrichment.csv"))
cat(sprintf("\nuniverse = CLUES (%d) | Shared=%d Specific=%d\n", length(U), length(SH), length(SP)))
cat(sprintf("===== CLUES PanelApp Mendelian gene enrichment (%s / %s) : baseline + adjusted log2 OR =====\n", DLAB, COND))
print(res[order(variant,geneset,group,mode), .(variant, geneset, group, mode, nSet, nOv, l2=round(l2,2), OR=round(OR,2), p=signif(p,2), stars)])
