rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

gsdir  <- "data/3a-benchmark_genesets"
hubdir <- "data/onek1k"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
CT  <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
rd  <- function(f) unique(trimws(readLines(f)))

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

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
cov <- merge(cov, fread("data/onek1k/beta_var.txt")[, .(gene, logBV=log(beta_var))], by="gene")

uni <- list(
  GECO     = rd(file.path(hubdir, "GECO_universe.txt")),
  disp     = rd(file.path(gsdir,  "disp_universe.txt")),
  inflated = rd(file.path(gsdir,  "inflated_universe.txt")))

sets <- list(
  list(f="housekeeping.txt",    u="GECO",     lab="Housekeeping genes"),
  list(f="LM22.txt",            u="GECO",     lab="Immune cell marker genes"),
  list(f="low_dispersion.txt",  u="disp",     lab="Lowly dispersed genes"),
  list(f="high_dispersion.txt", u="disp",     lab="Highly dispersed genes"),
  list(f="dd.txt",              u="disp",     lab="Differentially dispersed genes"),
  list(f="inflated_asone.txt",  u="inflated", lab="Bursty genes"))

fit1 <- function(hubG, setG, U, adj) {
  d <- data.table(gene=U, in_hub=as.integer(U %in% hubG), in_set=as.integer(U %in% setG))
  if (adj) {
    d <- merge(d, cov, by="gene")
    d <- d[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(logBV)]
    bad <- list(l2=NA_real_, se=NA_real_, p=NA_real_, n=nrow(d),
                ov=sum(d$in_hub & d$in_set), nhub=sum(d$in_hub))
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + logBV,
                                       d, family=binomial)), error=function(e) NULL)
    if (is.null(m) || !m$converged) return(bad)
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(bad)
    b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
  } else {
    a <- sum(d$in_hub & d$in_set);  b2 <- sum(d$in_hub & !d$in_set)
    c2 <- sum(!d$in_hub & d$in_set); e2 <- sum(!d$in_hub & !d$in_set)
    if (a==0 || b2==0 || c2==0 || e2==0) {
      aa<-a+.5; bb<-b2+.5; cc<-c2+.5; dd<-e2+.5
      b <- log((aa*dd)/(bb*cc)); se <- sqrt(1/aa+1/bb+1/cc+1/dd); p <- 2*pnorm(-abs(b/se))
    } else {
      co <- coef(summary(glm(in_set ~ in_hub, d, family=binomial)))
      b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
    }
  }
  list(l2=b/log(2), se=se/log(2), p=p, n=nrow(d), ov=sum(d$in_hub & d$in_set), nhub=sum(d$in_hub))
}

grid <- CJ(si=seq_along(sets), mth=names(methods), gp=c("Bulk","Pseudobulk"),
           adj=c(FALSE,TRUE), sorted=FALSE)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  s <- sets[[grid$si[i]]]; U <- uni[[s$u]]
  G <- intersect(rd(file.path(gsdir, s$f)), U)
  H <- intersect(methods[[grid$mth[i]]][[grid$gp[i]]], U)
  r <- fit1(H, G, U, grid$adj[i])
  data.table(set=s$lab, universe=s$u, method=grid$mth[i], group=grid$gp[i],
             mode=ifelse(grid$adj[i], "Adjusted", "Baseline"),
             nU=length(U), nSet=length(G), nHub=r$nhub, overlap=r$ov, n=r$n,
             log2OR=r$l2, se=r$se, p=r$p)
}))
res[, OR := 2^log2OR]
res[, stars := fifelse(is.na(p), "", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
SETORD <- vapply(sets, function(s) s$lab, character(1))
setorder(res, set, method, group, mode)
fwrite(res, file.path(outdir, "geneset_enrichment_3methods.csv"))

cat("hub set sizes (in the GECO universe; each set is re-intersected per gene set's universe):\n")
