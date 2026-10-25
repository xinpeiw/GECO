rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

moldir <- "data/6a-trans_molGenes"
hubdir <- "data/onek1k"
outdir <- "results/others"
CT  <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
rd  <- function(f) unique(trimws(readLines(f)))

GECOU <- rd(file.path(hubdir, "GECO_universe.txt"))
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
MTH <- names(methods)

covall <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
covall <- merge(covall, fread("data/onek1k/beta_var.txt")[, .(gene, lbv = log(beta_var))], by = "gene")

datasets <- list(
  list(f = "trans_eGene/trans_onek1k.txt", uni = "onek1k_universe.txt",         lab = "OneK1K",  type = "trans-eGene"),
  list(f = "trans_eGene/eQTLGen.txt",      uni = "eQTLGen_universe_trans.txt",  lab = "eQTLGen", type = "trans-eGene"),
  list(f = "trans_pGene/UKB-PPP.txt",      uni = "universe_UKB-PPP.txt",        lab = "UKB-PPP", type = "trans-pGene"),
  list(f = "trans_pGene/deCODE.txt",       uni = "universe_deCODE.txt",         lab = "deCODE",  type = "trans-pGene"),
  list(f = "trans_meGene/GoDMC.txt",       uni = "GECO",                        lab = "GoDMC",   type = "trans-meGene"))

fit1 <- function(hubG, setG, U, adj) {
  d <- data.table(gene = U, in_hub = as.integer(U %in% hubG), in_set = as.integer(U %in% setG))
  if (adj) {
    d <- merge(d, covall, by = "gene"); d <- d[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(lbv)]
    bad <- list(l2 = NA_real_, se = NA_real_, p = NA_real_, ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub))
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + lbv, d, family = binomial)), error = function(e) NULL)
    if (is.null(m) || !m$converged) return(bad)
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(bad)
    se <- co["in_hub","Std. Error"]; if (!is.finite(se) || se > 10) return(bad)
    b <- co["in_hub","Estimate"]; p <- co["in_hub","Pr(>|z|)"]
  } else {
    a <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set); c2 <- sum(!d$in_hub & d$in_set); e2 <- sum(!d$in_hub & !d$in_set)
    if (a == 0 || b2 == 0 || c2 == 0 || e2 == 0) {
      aa<-a+.5; bb<-b2+.5; cc<-c2+.5; dd<-e2+.5
      or <- (aa*dd)/(bb*cc); b <- log(or); se <- sqrt(1/aa+1/bb+1/cc+1/dd); p <- 2*pnorm(-abs(b/se))
    } else {
      m <- glm(in_set ~ in_hub, d, family = binomial); co <- coef(summary(m))
      b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
    }
  }
  list(l2 = b/log(2), se = se/log(2), p = p, ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub))
}

UG <- lapply(datasets, function(s) {
  U <- if (s$uni == "GECO") GECOU else intersect(rd(file.path(moldir, s$uni)), GECOU)
  list(lab = s$lab, type = s$type, uni = s$uni, U = U,
       G = intersect(rd(file.path(moldir, s$f)), U))
})
grid <- CJ(i = seq_along(UG), mth = MTH, gp = c("Bulk","Pseudobulk"),
           adj = c(FALSE, TRUE), sorted = FALSE)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(k) {
  x <- UG[[grid$i[k]]]
  r <- fit1(intersect(methods[[grid$mth[k]]][[grid$gp[k]]], x$U), x$G, x$U, grid$adj[k])
  data.table(dataset = x$lab, type = x$type,
             universe = ifelse(x$uni == "GECO", "GECO", basename(x$uni)),
             method = grid$mth[k], group = grid$gp[k],
             mode = ifelse(grid$adj[k], "Adjusted", "Baseline"),
             nU = length(x$U), nSet = length(x$G), nHub = r$nhub, nOverlap = r$ov,
             l2 = r$l2, se = r$se, p = r$p)
}))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p), "", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[order(dataset, method, group, mode)],
       file.path(outdir, "trans_QTL_enrichment_3methods.csv"))

cat("hub set sizes (in the GECO universe; re-intersected per dataset universe):\n")
