rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

moldir <- "data/6a-trans_molGenes"
hubdir <- "data/onek1k"
outdir <- "results/6a-trans_molGenes"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) unique(trimws(readLines(f)))
GECOU    <- rd(file.path(hubdir, "GECO_universe.txt"))
SHARED   <- rd(file.path(hubdir, "GECO_Shared.txt"))
SPECIFIC <- rd(file.path(hubdir, "GECO_Specific.txt"))
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
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + lbv, d, family = binomial)), error = function(e) NULL)
    if (is.null(m) || !m$converged) return(list(l2 = NA, se = NA, p = NA))
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(list(l2 = NA, se = NA, p = NA))
    se <- co["in_hub","Std. Error"]; if (!is.finite(se) || se > 10) return(list(l2 = NA, se = NA, p = NA))
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
  list(l2 = b/log(2), se = se/log(2), p = p)
}

res <- rbindlist(lapply(datasets, function(s) {
  U  <- if (s$uni == "GECO") GECOU else intersect(rd(file.path(moldir, s$uni)), GECOU)
  G  <- intersect(rd(file.path(moldir, s$f)), U)
  SH <- intersect(SHARED, U); SP <- intersect(SPECIFIC, U)
  rbindlist(lapply(list(list("GECO Shared", SH), list("GECO Specific", SP)), function(g)
    rbindlist(lapply(c(FALSE, TRUE), function(adj) {
      r <- fit1(g[[2]], G, U, adj)
      data.table(dataset = s$lab, type = s$type, universe = ifelse(s$uni=="GECO","GECO",basename(s$uni)),
                 nU = length(U), nSet = length(G), group = g[[1]], mode = ifelse(adj, "Adjusted", "Baseline"),
                 l2 = r$l2, se = r$se, p = r$p)
    }))))
}))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p), "", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[order(dataset, group, mode)], file.path(outdir, "trans_QTL_enrichment.csv"))
cat("\n===== trans-QTL enrichment (log2 OR, hub vs non-hub) =====\n")
print(res[, .(dataset, type, universe, nU, nSet, group, mode, l2 = round(l2,2), OR = round(OR,2), p = signif(p,2), stars)])
