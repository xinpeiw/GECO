rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

gdir   <- "data/5c-GWAS_gene"
hubdir <- "data/onek1k"
outdir <- "results/5c-GWAS_gene"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(grep("^ENSG", trimws(readLines(f)), value=TRUE))

U  <- rd(file.path(hubdir, "GECO_universe.txt"))
SH <- intersect(rd(file.path(hubdir, "GECO_Shared.txt")),   U)
SP <- intersect(rd(file.path(hubdir, "GECO_Specific.txt")), U)
covall <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
covall <- merge(covall, fread("data/onek1k/beta_var.txt")[, .(gene, lbv=log(beta_var))], by="gene")

datasets <- list(list(f="genelist_Immune.txt",    lab="Immune GWAS genes"),
                 list(f="genelist_NonImmune.txt", lab="Non-immune GWAS genes"))

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

res <- rbindlist(lapply(datasets, function(s) {
  Graw <- rd(file.path(gdir, s$f)); G <- intersect(Graw, U)
  cat(sprintf("  %-22s file=%4d | in GECO universe=%4d (dropped %d)\n", s$lab, length(Graw), length(G), length(Graw)-length(G)))
  rbindlist(lapply(list(list("GECO Shared",SH), list("GECO Specific",SP)), function(g)
    rbindlist(lapply(c(FALSE,TRUE), function(adj) { r <- fit1(g[[2]], G, adj)
      data.table(dataset=s$lab, nU=length(U), nSet=length(G), nOv=r$nov, group=g[[1]],
                 mode=ifelse(adj,"Adjusted","Baseline"), l2=r$l2, se=r$se, p=r$p) }))))
}))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p),"", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[order(dataset, group, mode)], file.path(outdir, "GWASgene_enrichment.csv"))
cat(sprintf("\nuniverse = GECO (%d) | Shared=%d Specific=%d\n", length(U), length(SH), length(SP)))
cat("\n===== GWAS gene enrichment (log2 OR, hub vs rest) =====\n")
print(res[order(dataset, group, mode), .(dataset, nSet, group, mode, nOv,
      l2=round(l2,2), OR=round(OR,2), se=round(se,3), p=signif(p,2), stars)], nrows=Inf)
