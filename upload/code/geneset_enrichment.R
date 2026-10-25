rm(list=ls())
suppressPackageStartupMessages({
  library(data.table);
})

gsdir  <- "data/3a-benchmark_genesets"
hubdir <- "data/onek1k"
outdir <- "results/3a-benchmark_genesets"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) unique(trimws(readLines(f)))
SH_all <- rd(file.path(hubdir, "GECO_Shared.txt"))
SP_all <- rd(file.path(hubdir, "GECO_Specific.txt"))

cov <- fread("data/onek1k/covariates.txt")
cov <- cov[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
bv  <- fread("data/onek1k/beta_var.txt")[, .(gene, logBV = log(beta_var))]
cov <- merge(cov, bv, by = "gene")

uni <- list(
  GECO     = rd(file.path(hubdir, "GECO_universe.txt")),
  disp     = rd(file.path(gsdir,  "disp_universe.txt")),
  inflated = rd(file.path(gsdir,  "inflated_universe.txt"))
)

sets <- list(
  list(f = "housekeeping.txt",    u = "GECO",     lab = "Housekeeping genes"),
  list(f = "LM22.txt",            u = "GECO",     lab = "Immune cell marker genes"),
  list(f = "low_dispersion.txt",  u = "disp",     lab = "Lowly dispersed genes"),
  list(f = "high_dispersion.txt", u = "disp",     lab = "Highly dispersed genes"),
  list(f = "dd.txt",              u = "disp",     lab = "Differentially dispersed genes"),
  list(f = "inflated_asone.txt",  u = "inflated", lab = "Bursty genes")

)

fit1 <- function(hubG, setG, U, adj) {
  d <- data.table(gene = U,
                  in_hub = as.integer(U %in% hubG),
                  in_set = as.integer(U %in% setG))
  if (adj) {
    d <- merge(d, cov, by = "gene")
    d <- d[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(logBV)]
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + logBV, d, family = binomial)),
                  error = function(e) NULL)
    if (is.null(m) || !m$converged) return(list(l2 = NA, se = NA, p = NA, n = nrow(d), ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub)))
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) return(list(l2 = NA, se = NA, p = NA, n = nrow(d), ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub)))
    b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
  } else {
    a <- sum(d$in_hub & d$in_set); b2 <- sum(d$in_hub & !d$in_set)
    c2 <- sum(!d$in_hub & d$in_set); e2 <- sum(!d$in_hub & !d$in_set)
    if (a == 0 || b2 == 0 || c2 == 0 || e2 == 0) {
      aa <- a+.5; bb <- b2+.5; cc <- c2+.5; dd <- e2+.5
      or <- (aa*dd)/(bb*cc); b <- log(or); se <- sqrt(1/aa+1/bb+1/cc+1/dd); p <- 2*pnorm(-abs(b/se))
    } else {
      m <- glm(in_set ~ in_hub, d, family = binomial); co <- coef(summary(m))
      b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
    }
  }
  list(l2 = b/log(2), se = se/log(2), p = p, n = nrow(d), ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub))
}

res <- rbindlist(lapply(sets, function(s) {
  U  <- uni[[s$u]]; G <- intersect(rd(file.path(gsdir, s$f)), U)
  SH <- intersect(SH_all, U); SP <- intersect(SP_all, U)
  rbindlist(lapply(list(list("GECO Shared", SH), list("GECO Specific", SP)), function(g)
    rbindlist(lapply(c(FALSE, TRUE), function(adj) {
      r <- fit1(g[[2]], G, U, adj)
      data.table(set = s$lab, universe = s$u, nU = length(U), nSet = length(G),
                 group = g[[1]], nHub = r$nhub,
                 mode = ifelse(adj, "Adjusted", "Baseline"),
                 log2OR = r$l2, se = r$se, p = r$p, n = r$n, overlap = r$ov)
    }))))
}))
res[, OR := 2^log2OR]
res[, stars := fifelse(p < .001, "***", fifelse(p < .01, "**", fifelse(p < .05, "*", "")))]

fwrite(res[order(set, group, mode), .(set, universe, group, mode, nSet, nHub, overlap, n,
        log2OR = round(log2OR,3), OR = round(OR,3), se = round(se,3), p, stars)],
       file.path(outdir, "geneset_enrichment.csv"))
cat("\n===== gene-set enrichment (log2 OR, hub vs non-hub) =====\n")
print(res[order(set, group, mode), .(set, group, mode, nSet, overlap,
        log2OR = round(log2OR,2), OR = round(OR,2), p = signif(p,2), stars)])
