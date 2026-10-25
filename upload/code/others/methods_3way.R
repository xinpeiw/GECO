rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
CT   <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
rd   <- function(f) unique(trimws(readLines(f)))

U      <- rd("data/onek1k/GECO_universe.txt")
covall <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
covall <- merge(covall, fread("data/onek1k/beta_var.txt")[, .(gene, lbv=log(beta_var))], by="gene")

methods <- list(
  GECO   = list(Bulk = rd("data/onek1k/GECO_Shared.txt"),
                Pseudobulk = rd("data/onek1k/GECO_Specific.txt")),
  WGCNA  = list(Bulk = rd("data/others/WGCNA_hubs/Shared.txt"),
                Pseudobulk = unique(unlist(lapply(CT, function(ct)
                  rd(file.path("data/others/WGCNA_hubs", paste0(ct, ".txt"))))))),
  glasso = list(Bulk = rd("data/others/glasso_hubs/hubs_Shared.txt"),
                Pseudobulk = unique(unlist(lapply(CT, function(ct)
                  rd(file.path("data/others/glasso_hubs", paste0("hubs_", ct, ".txt"))))))))

padir <- "data/5d-Mendelian_gene"
genesets <- list(
  list(an="GWAS genes",     cls="Immune",     f="data/5c-GWAS_gene/genelist_Immune.txt"),
  list(an="GWAS genes",     cls="Non-immune", f="data/5c-GWAS_gene/genelist_NonImmune.txt"),
  list(an="Drug targets",   cls="Immune",     f="data/5e-drug_target/targets_Immune_approved.txt"),
  list(an="Drug targets",   cls="Non-immune", f="data/5e-drug_target/targets_NonImmune_approved.txt"),
  list(an="Mendelian genes",cls="Immune",     f=file.path(padir, "mendelian_immune_full_in_GECO_universe.txt"),    variant="Full"),
  list(an="Mendelian genes",cls="Non-immune", f=file.path(padir, "mendelian_nonimmune_full_in_GECO_universe.txt"), variant="Full"))

fit1 <- function(hubG, setG, adj) {
  d <- data.table(gene=U, in_hub=as.integer(U %in% hubG), in_set=as.integer(U %in% setG))
  ov <- sum(d$in_hub & d$in_set)
  if (adj) {
    d <- merge(d, covall, by="gene")
    d <- d[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(lbv)]
    a2 <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set)
    bad <- list(l2=NA_real_, se=NA_real_, p=NA_real_, ov=a2)
    if (a2 == 0 || b2 == 0) return(bad)
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + lbv, d, family=binomial)),
                  error=function(e) NULL)
    if (is.null(m) || !m$converged) return(bad)
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(bad)
    se <- co["in_hub","Std. Error"]; if (!is.finite(se) || se > 10) return(bad)
    list(l2=co["in_hub","Estimate"]/log(2), se=se/log(2), p=co["in_hub","Pr(>|z|)"], ov=a2)
  } else {
    a <- sum(d$in_hub & d$in_set);  b <- sum(d$in_hub & !d$in_set)
    c <- sum(!d$in_hub & d$in_set); e <- sum(!d$in_hub & !d$in_set)
    if (a==0 || b==0 || c==0 || e==0) { a<-a+.5; b<-b+.5; c<-c+.5; e<-e+.5 }
    l2 <- log((a*e)/(b*c))/log(2); se <- sqrt(1/a+1/b+1/c+1/e)/log(2)
    list(l2=l2, se=se, p=2*pnorm(-abs(l2/se)), ov=ov)
  }
}

grid <- CJ(si = seq_along(genesets), mth = names(methods),
           gp = c("Bulk","Pseudobulk"), adj = c(FALSE, TRUE), sorted = FALSE)
gl <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  s <- genesets[[grid$si[i]]]
  S <- intersect(rd(s$f), U)
  H <- intersect(methods[[grid$mth[i]]][[grid$gp[i]]], U)
  r <- fit1(H, S, grid$adj[i])
  data.table(analysis=s$an, class=s$cls,
             variant=if (is.null(s$variant)) "-" else s$variant,
             method=grid$mth[i], group=grid$gp[i],
             mode=ifelse(grid$adj[i], "Adjusted", "Baseline"),
             nU=length(U), nHub=length(H), nSet=length(S), nOv=r$ov,
             est=r$l2, se=r$se, p=r$p)
}))

sl <- rbindlist(lapply(list(
    list(mth="GECO",   f="data/5b-ldsc/sldsc_meta.tsv"),
    list(mth="WGCNA",  f="data/others/sldsc_meta_WGCNA.tsv"),
    list(mth="glasso", f="data/others/sldsc_meta_glasso.tsv")),
  function(s) {
    d <- fread(s$f)
    d[, group := fcase(hub %in% c("Shared","Bulk"), "Bulk",
                       hub %in% c("Specific","Pseudobulk"), "Pseudobulk")]
    stopifnot(!anyNA(d$group), setequal(d$Category, c("Immune","Non-Immune")))
    data.table(analysis="S-LDSC", class=fifelse(d$Category=="Immune","Immune","Non-immune"),
               variant="-", method=s$mth, group=d$group, mode="Baseline",
               nU=NA_integer_, nHub=NA_integer_, nSet=d$n, nOv=NA_integer_,
               est=d$tau_star, se=d$se, p=d$p)
  }))

res <- rbind(gl, sl)
res[, stars := fifelse(is.na(p), "", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
setorder(res, analysis, class, method, group, mode)
fwrite(res, file.path(outdir, "methods_3way.csv"))

cat("hub set sizes within the GECO universe:\n")
