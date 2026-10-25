rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

COND   <- "mega"
COVTAG <- "all"
DLAB   <- "all"

ddir   <- "data/7d_drug_target"
covdir <- "data/clues"
outdir <- "results/7d_drug_target"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(grep("^ENSG", trimws(readLines(f)), value=TRUE))

U  <- rd("data/clues/allgenes_CLUES.txt")
SH <- intersect(rd("data/clues/hubs_alpha.txt"),    U)
SP <- intersect(rd("data/clues/hubs_specific.txt"), U)
covall <- fread(file.path(covdir, paste0("covariates_", COVTAG, ".txt")))[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
covall <- merge(covall, fread(file.path(covdir, paste0("beta_var_", COVTAG, ".txt")))[, .(gene, lbv=log(beta_var))], by="gene")

grps  <- list(list(k="Immune",    lab="Immune"), list(k="NonImmune", lab="Non-immune"))
tiers <- list(list(k="approved", lab="Approved"))

fit1 <- function(hubG, setG, adj) {
  d <- data.table(gene=U, in_hub=as.integer(U %in% hubG), in_set=as.integer(U %in% setG))
  if (adj) {
    d <- merge(d, covall, by="gene"); d <- d[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]
    bad <- list(l2=NA_real_, se=NA_real_, p=NA_real_, nov=sum(d$in_hub & d$in_set))
    a2 <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set)
    if (a2 == 0 || b2 == 0) return(bad)
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + lbv, d, family=binomial)), error=function(e) NULL)
    if (is.null(m) || !m$converged) return(bad)
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(bad)
    se <- co["in_hub","Std. Error"]; if (!is.finite(se) || se > 10) return(bad)
    b <- co["in_hub","Estimate"]; p <- co["in_hub","Pr(>|z|)"]
  } else {
    a <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set)
    c2 <- sum(!d$in_hub & d$in_set); e2 <- sum(!d$in_hub & !d$in_set)
    if (a==0 || b2==0 || c2==0 || e2==0) { aa<-a+.5; bb<-b2+.5; cc<-c2+.5; dd<-e2+.5
      b <- log((aa*dd)/(bb*cc)); se <- sqrt(1/aa+1/bb+1/cc+1/dd); p <- 2*pnorm(-abs(b/se))
    } else { m <- glm(in_set ~ in_hub, d, family=binomial); co <- coef(summary(m))
      b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"] }
  }
  list(l2=b/log(2), se=se/log(2), p=p, nov=sum(d$in_hub & d$in_set))
}

out <- list(); k <- 0L
for (gr in grps) for (ti in tiers) {
  G <- intersect(rd(file.path(ddir, sprintf("targets_%s_%s.txt", gr$k, ti$k))), U)
  cat(sprintf("  %-11s %-9s n=%3d\n", gr$lab, ti$lab, length(G)))
  for (g in list(list("CLUES Shared",SH), list("CLUES Specific",SP))) for (adj in c(FALSE,TRUE)) {
    r <- fit1(g[[2]], G, adj); k <- k + 1L
    out[[k]] <- data.table(group=gr$lab, tier=ti$lab, cond=DLAB, nU=length(U), nSet=length(G), nOv=r$nov,
                           hub=g[[1]], mode=ifelse(adj,"Adjusted","Baseline"),
                           l2=r$l2, se=r$se, p=r$p)
  }
}
res <- rbindlist(out)
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p),"", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[order(group, tier, hub, mode)], file.path(outdir, "CLUES_DrugTargets_enrichment.csv"))
cat(sprintf("\nuniverse = CLUES %s (%d) | Shared=%d Specific=%d\n", DLAB, length(U), length(SH), length(SP)))
cat("\n===== CLUES drug-target enrichment (log2 OR, hub vs rest) =====\n")
print(res[order(group, tier, hub, mode), .(group, tier, nSet, hub, mode, nOv,
      l2=round(l2,2), OR=round(OR,2), se=round(se,3), p=signif(p,2), stars)], nrows=Inf)
