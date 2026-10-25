rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

gsdir  <- "data/3a-benchmark_genesets"
hubdir <- "data/onek1k"
pmdir  <- "data/onek1k_precision"
outdir <- "results/others"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
rd  <- function(f) unique(trimws(readLines(f)))
CT  <- c("CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")
PCT <- c(0.01, 0.05, 0.10, 0.15, 0.20); PLAB <- sprintf("top %d%%", round(100*PCT))

DEG <- lapply(c("alpha", CT), function(n) {
  m <- readRDS(file.path(pmdir, paste0(n, ".rds")))
  d <- data.table(gene = rownames(m), degree = rowSums(m != 0) - 1L); rm(m); gc(FALSE)
  d[degree > 0]
}); names(DEG) <- c("alpha", CT)
hubs_p <- function(net, p) DEG[[net]][degree > quantile(degree, 1 - p), gene]
HUB <- lapply(PCT, function(p) list(Shared   = hubs_p("alpha", p),
                                    Specific = unique(unlist(lapply(CT, hubs_p, p = p)))))
names(HUB) <- PLAB
stopifnot(setequal(HUB[["top 5%"]]$Shared,   rd(file.path(hubdir, "GECO_Shared.txt"))),
          setequal(HUB[["top 5%"]]$Specific, rd(file.path(hubdir, "GECO_Specific.txt"))))
cat("check: the 5% rule reproduces GECO_Shared.txt and GECO_Specific.txt exactly\n")

hc <- rbindlist(lapply(seq_along(PCT), function(i) data.table(
  threshold = PLAB[i], network = names(DEG), n_connected = sapply(DEG, nrow),
  n_hub = sapply(names(DEG), function(n) length(hubs_p(n, PCT[i]))))))
fwrite(hc, file.path(outdir, "geneset_enrichment_toppct_hubcounts.csv"))
cat("\nhubs per network and threshold:\n")
print(dcast(hc, network ~ factor(threshold, levels = PLAB), value.var = "n_hub")[match(c("alpha", CT), network)], row.names = FALSE)
cat("\nhub set sizes in the GECO universe:\n")
for (h in PLAB) cat(sprintf("  %-8s Shared %4d | Specific %4d | overlap %3d\n", h,
  length(HUB[[h]]$Shared), length(HUB[[h]]$Specific), length(intersect(HUB[[h]]$Shared, HUB[[h]]$Specific))))

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
cov <- merge(cov, fread("data/onek1k/beta_var.txt")[, .(gene, logBV = log(beta_var))], by = "gene")
uni <- list(GECO     = rd(file.path(hubdir, "GECO_universe.txt")),
            disp     = rd(file.path(gsdir,  "disp_universe.txt")),
            inflated = rd(file.path(gsdir,  "inflated_universe.txt")))
sets <- list(
  list(f = "housekeeping.txt",    u = "GECO",     lab = "Housekeeping genes"),
  list(f = "LM22.txt",            u = "GECO",     lab = "Immune cell marker genes"),
  list(f = "low_dispersion.txt",  u = "disp",     lab = "Lowly dispersed genes"),
  list(f = "high_dispersion.txt", u = "disp",     lab = "Highly dispersed genes"),
  list(f = "dd.txt",              u = "disp",     lab = "Differentially dispersed genes"),
  list(f = "inflated_asone.txt",  u = "inflated", lab = "Bursty genes"))
SETLAB <- vapply(sets, function(s) s$lab, character(1))

fit1 <- function(hubG, setG, U, adj) {
  d <- data.table(gene = U, in_hub = as.integer(U %in% hubG), in_set = as.integer(U %in% setG))
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
      b <- log((aa*dd)/(bb*cc)); se <- sqrt(1/aa+1/bb+1/cc+1/dd); p <- 2*pnorm(-abs(b/se))
    } else {
      m <- glm(in_set ~ in_hub, d, family = binomial); co <- coef(summary(m))
      b <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
    }
  }
  list(l2 = b/log(2), se = se/log(2), p = p, n = nrow(d), ov = sum(d$in_hub & d$in_set), nhub = sum(d$in_hub))
}

res <- rbindlist(lapply(PLAB, function(h) rbindlist(lapply(sets, function(s) {
  U <- uni[[s$u]]; G <- intersect(rd(file.path(gsdir, s$f)), U)
  rbindlist(lapply(c("Shared","Specific"), function(g) rbindlist(lapply(c(FALSE, TRUE), function(adj) {
    r <- fit1(intersect(HUB[[h]][[g]], U), G, U, adj)
    data.table(threshold = h, set = s$lab, universe = s$u, group = paste("GECO", g),
               mode = ifelse(adj, "Adjusted", "Baseline"),
               nU = length(U), nSet = length(G), nHub = r$nhub, overlap = r$ov, n = r$n,
               log2OR = r$l2, se = r$se, p = r$p)
  }))))
}))))
res[, OR := 2^log2OR]
res[, stars := fifelse(is.na(p), "", fifelse(p < .001, "***", fifelse(p < .01, "**", fifelse(p < .05, "*", ""))))]
fwrite(res[, .(threshold, set, universe, group, mode, nU, nSet, nHub, overlap, n,
               log2OR = round(log2OR, 3), OR = round(OR, 3), se = round(se, 3), p, stars)],
       file.path(outdir, "geneset_enrichment_toppct.csv"))

pub <- fread("results/3a-benchmark_genesets/geneset_enrichment.csv")
mm  <- merge(res[threshold == "top 5%", .(set, group, mode, new = round(log2OR, 3), new_ov = overlap)],
             pub[, .(set, group, mode, old = log2OR, old_ov = overlap)], by = c("set","group","mode"))
stopifnot(nrow(mm) == 24L, identical(mm$new_ov, mm$old_ov), max(abs(mm$new - mm$old)) < 1e-3)
cat(sprintf("\ncheck: top 5%% rows reproduce geneset_enrichment.csv (%d cells, max |diff| = %.1e)\n", nrow(mm), max(abs(mm$new - mm$old))))
cat("\n===== log2 OR by hub threshold =====\n")
print(dcast(res[, .(threshold, set, group, mode, v = sprintf("%5.2f%s", log2OR, stars))],
            set + group + mode ~ factor(threshold, levels = PLAB), value.var = "v")[order(match(set, SETLAB), group, mode)], nrows = Inf)
