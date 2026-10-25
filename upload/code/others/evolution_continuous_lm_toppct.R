rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

evo    <- "data/3cd-evolution"
hubdir <- "data/onek1k"
pmdir  <- "data/onek1k_precision"
outdir <- "results/others"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
rd  <- function(f) unique(trimws(readLines(f)))
CT  <- c("CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")
PCT <- c(0.01, 0.05, 0.10, 0.15, 0.20); PLAB <- sprintf("top %d%%", round(100*PCT))
U   <- rd(file.path(hubdir, "GECO_universe.txt"))

DEG <- lapply(c("alpha", CT), function(n) {
  m <- readRDS(file.path(pmdir, paste0(n, ".rds")))
  d <- data.table(gene = rownames(m), degree = rowSums(m != 0) - 1L); rm(m); gc(FALSE)
  d[degree > 0]
}); names(DEG) <- c("alpha", CT)
hubs_p <- function(net, p) DEG[[net]][degree > quantile(degree, 1 - p), gene]
HUB <- lapply(PCT, function(p) list(Shared   = intersect(hubs_p("alpha", p), U),
                                    Specific = intersect(unique(unlist(lapply(CT, hubs_p, p = p))), U)))
names(HUB) <- PLAB
stopifnot(setequal(HUB[["top 5%"]]$Shared,   intersect(rd(file.path(hubdir, "GECO_Shared.txt")), U)),
          setequal(HUB[["top 5%"]]$Specific, intersect(rd(file.path(hubdir, "GECO_Specific.txt")), U)))
cat("check: the 5% rule reproduces GECO_Shared.txt and GECO_Specific.txt exactly\n")
for (h in PLAB) cat(sprintf("  %-8s Shared %4d | Specific %4d | overlap %3d\n", h,
  length(HUB[[h]]$Shared), length(HUB[[h]]$Specific), length(intersect(HUB[[h]]$Shared, HUB[[h]]$Specific))))

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
rc  <- fread(file.path(evo, "GECO_gene_recombination_rate.1000G.hg38.tsv"))[, .(gene = gene_id, recomb = recomb_rate_cM_Mb)]
ex  <- fread(file.path(evo, "exon.txt"))[, .(gene = ensg, nexon = n_exons_unique)]
bv  <- fread("data/onek1k/beta_var.txt")[, .(gene, beta_var)]
covdt <- cov[, .(gene, logE, logL, gc)]
covdt[rc, on = "gene", lrec := log1p(i.recomb)]
covdt[ex, on = "gene", lexn := log(i.nexon)]
covdt[bv, on = "gene", lbv := log(i.beta_var)]
covdt <- covdt[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(lrec) & is.finite(lexn) & is.finite(lbv)]

vars <- list(list(file = "loeuf.txt",  col = "loeuf",  lab = "LOEUF"),
             list(file = "phylop.txt", col = "phylop", lab = "phyloP"),
             list(file = "B.txt",      col = "bvalue", lab = "B statistic"))
VLAB <- vapply(vars, function(v) v$lab, character(1))
MD <- lapply(vars, function(v) { m <- fread(file.path(evo, v$file)); setnames(m, 1, "gene")
  m[gene %in% U, .(gene, val = as.numeric(get(v$col)))][is.finite(val)] }); names(MD) <- VLAB

res <- rbindlist(lapply(PLAB, function(h) rbindlist(lapply(VLAB, function(v) {
  rbindlist(lapply(c("Shared","Specific"), function(g) rbindlist(lapply(c(FALSE, TRUE), function(adj) {
    d <- copy(MD[[v]]); d[, x := as.integer(gene %in% HUB[[h]][[g]])]
    if (adj) { d <- merge(d, covdt, by = "gene"); mm <- lm(val ~ x + logE + logL + gc + lexn + lrec + lbv, d) }
    else     { mm <- lm(val ~ x, d) }
    co <- coef(summary(mm))["x", ]
    data.table(threshold = h, variable = v, group = paste("GECO", g), mode = ifelse(adj, "Adjusted", "Baseline"),
               n = nrow(d), n_hub = sum(d$x), beta = co[1], se = co[2], p = co[4])
  }))))
}))))
res[, stars := fifelse(p < .001, "***", fifelse(p < .01, "**", fifelse(p < .05, "*", "")))]
fwrite(res, file.path(outdir, "evolution_continuous_lm_toppct.csv"))

pub <- fread("results/3cd-evolution/evolution_continuous_lm.csv")
mm  <- merge(res[threshold == "top 5%", .(variable, group, mode, new = beta, new_n = n_hub)],
             pub[, .(variable, group, mode, old = beta, old_n = n_hub)], by = c("variable","group","mode"))
stopifnot(nrow(mm) == 12L, identical(mm$new_n, mm$old_n), max(abs(mm$new - mm$old)) < 1e-8)
cat(sprintf("\ncheck: top 5%% rows reproduce evolution_continuous_lm.csv (%d cells, max |diff| = %.1e)\n", nrow(mm), max(abs(mm$new - mm$old))))
cat("\n===== lm coefficient by hub threshold =====\n")
print(dcast(res[, .(threshold, variable, group, mode, v = sprintf("%6.3f%s", beta, stars))],
            variable + group + mode ~ factor(threshold, levels = PLAB), value.var = "v")[order(match(variable, VLAB), group, mode)], nrows = Inf)
