rm(list=ls())

suppressPackageStartupMessages({ library(data.table) })

COND   <- "mega"
COVTAG <- "all"
DLAB   <- "all"

cisdir <- "data/others/molQTL_clues/cis"
covdir <- "data/clues"
outdir <- "results/others"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) unique(trimws(readLines(f)))
CLUESU   <- rd("data/clues/allgenes_CLUES.txt")
SHARED   <- rd("data/clues/hubs_alpha.txt")
SPECIFIC <- rd("data/clues/hubs_specific.txt")
covall <- fread(file.path(covdir, paste0("covariates_", COVTAG, ".txt")))[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
covall <- merge(covall, fread(file.path(covdir, paste0("beta_var_", COVTAG, ".txt")))[, .(gene, lbv = log(beta_var))], by = "gene")

qtls <- list(
  list(qtl = "cis-eGene", datasets = list(

    list(f = "cis_eQTLGen.txt",           uni = "../eQTLGen_universe_cis.txt",  lab = "eQTLGen"),
    list(f = "cis_onek1k.txt",            uni = "../onek1k_universe.txt",       lab = "OneK1K"),
    list(f = "cis_DICE.txt",              uni = "../DICE_universe.txt",         lab = "DICE"),
    list(f = "ImmuNexUT.txt",             uni = "CLUES",                        lab = "ImmuNexUT"),
    list(f = "CLUES_shared.txt",          uni = "../clues_universe.txt",        lab = "CLUES shared"),
    list(f = "CLUES_specific.txt",        uni = "../clues_universe.txt",        lab = "CLUES specific"),
    list(f = "fastGxC_GTEx_shared.txt",   uni = "../fastGxC_GTEx_universe.txt", lab = "GTEx shared"),
    list(f = "fastGxC_GTEx_specific.txt", uni = "../fastGxC_GTEx_universe.txt", lab = "GTEx specific"),
    list(f = "MemoryT.txt",               uni = "../MemoryT_universe.txt",      lab = "Memory T-cell\nstate-specific"))),
  list(qtl = "cis-pGene", datasets = list(
    list(f = "UKB-PPP.txt", uni = "../universe_UKB-PPP.txt", lab = "UKB-PPP"),
    list(f = "deCODE.txt",  uni = "../universe_deCODE.txt",  lab = "deCODE"))),
  list(qtl = "cis-meGene", datasets = list(
    list(f = "GoDMC.txt", uni = "CLUES", lab = "GoDMC"))),
  list(qtl = "cis-caGene", datasets = list(
    list(f = "ca_pbmc.txt",      uni = "CLUES", lab = "PBMC"),
    list(f = "ca_intestine.txt", uni = "CLUES", lab = "Intestine"))))

fit1 <- function(hubG, setG, U, adj) {
  d <- data.table(gene = U, in_hub = as.integer(U %in% hubG), in_set = as.integer(U %in% setG))
  if (adj) {
    d <- merge(d, covall, by = "gene"); d <- d[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(lbv)]
    bad <- list(l2 = NA_real_, se = NA_real_, p = NA_real_, n = nrow(d), ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub))
    a2 <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set)
    if (a2 == 0 || b2 == 0) return(bad)
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + lbv, d, family = binomial)), error = function(e) NULL)
    if (is.null(m) || !m$converged) return(bad)
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(bad)
    se <- co["in_hub","Std. Error"]; if (!is.finite(se) || se > 10) return(bad)
    b <- co["in_hub","Estimate"]; p <- co["in_hub","Pr(>|z|)"]; n <- nrow(d)
  } else {
    a <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set); c2 <- sum(!d$in_hub & d$in_set); e2 <- sum(!d$in_hub & !d$in_set)
    if (a == 0 || b2 == 0 || c2 == 0 || e2 == 0) {
      aa<-a+.5; bb<-b2+.5; cc<-c2+.5; dd<-e2+.5
      or <- (aa*dd)/(bb*cc); b <- log(or); se <- sqrt(1/aa+1/bb+1/cc+1/dd); p <- 2*pnorm(-abs(b/se))
    } else {
      m <- glm(in_set ~ in_hub, d, family = binomial); co <- coef(summary(m))
      b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
    }
    n <- length(U)
  }
  list(l2 = b/log(2), se = se/log(2), p = p, n = n, ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub))
}

res <- rbindlist(lapply(qtls, function(q) rbindlist(lapply(q$datasets, function(s) {
  U  <- if (s$uni == "CLUES") CLUESU else intersect(rd(file.path(cisdir, s$uni)), CLUESU)
  G  <- intersect(rd(file.path(cisdir, s$f)), U)
  SH <- intersect(SHARED, U); SP <- intersect(SPECIFIC, U)
  rbindlist(lapply(list(list("CLUES Shared", SH), list("CLUES Specific", SP)), function(g)
    rbindlist(lapply(c(FALSE, TRUE), function(adj) {
      r <- fit1(g[[2]], G, U, adj)
      data.table(dataset_grp = DLAB, qtl = q$qtl, dataset = s$lab,
                 universe = ifelse(s$uni == "CLUES", "allgenes_CLUES", basename(s$uni)),
                 nU = length(U), nSet = length(G), nHub = r$nhub, nOverlap = r$ov,
                 group = g[[1]], mode = ifelse(adj, "Adjusted", "Baseline"), l2 = r$l2, se = r$se, p = r$p)
    }))))
}))))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p), "", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))) ]
fwrite(res[order(qtl, dataset, group, mode)], file.path(outdir, "CLUES_cis_QTL_enrichment.csv"))
cat(sprintf("\n===== CLUES cis-molQTL enrichment (%s / %s) : log2 OR, hub vs rest =====\n", DLAB, COND))
cat(sprintf("universe allgenes_CLUES=%d | Shared=%d Specific=%d | Shared n Specific=%d\n",
            length(CLUESU), length(intersect(SHARED,CLUESU)), length(intersect(SPECIFIC,CLUESU)),
            length(intersect(intersect(SHARED,CLUESU), SPECIFIC))))
print(res[order(qtl, dataset, group, mode), .(qtl, dataset, nU, nSet, nHub, nOverlap, group, mode,
                                              l2 = round(l2,2), OR = round(OR,2), p = signif(p,2), stars)], nrows = Inf)
