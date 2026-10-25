rm(list=ls())

suppressPackageStartupMessages({ library(data.table) })

hubdir <- "data/onek1k"
outdir <- "results/others"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) unique(trimws(readLines(f)))

GECOU    <- rd(file.path(hubdir, "GECO_universe.txt"))
SHARED   <- rd(file.path(hubdir, "GECO_Shared.txt"))
SPECIFIC <- rd(file.path(hubdir, "GECO_Specific.txt"))
cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
bv  <- fread("data/onek1k/beta_var.txt")[, .(gene, logBV = log(beta_var))]
cov <- merge(cov, bv, by = "gene")[gene %in% GECOU]
cov <- cov[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(logBV)]

shared <- intersect(cov$gene, SHARED); specific <- intersect(cov$gene, SPECIFIC)
nonhub <- setdiff(cov$gene, union(shared, specific))
grp <- function(g, lab) cov[gene %in% g][, group := lab]
D <- rbind(grp(nonhub, "Non-hub"), grp(shared, "Shared hub"), grp(specific, "Specific hub"))
D[, group := factor(group, levels = c("Non-hub","Shared hub","Specific hub"))]

covlabs <- c(logE = "log(expression)", logBV = "log(expression variability)", logL = "log(gene length)", gc = "GC content")
Dm <- melt(D, id.vars = c("gene","group"), measure.vars = names(covlabs), variable.name = "cov", value.name = "value")
Dm[, cov := factor(covlabs[as.character(cov)], levels = covlabs)]

wt <- function(h, cvn) {
  a <- D[group == h, get(cvn)]; b <- D[group == "Non-hub", get(cvn)]
  w <- suppressWarnings(wilcox.test(a, b)); A <- as.numeric(w$statistic)/(length(a)*length(b))
  data.table(covariate = covlabs[cvn], comparison = paste(sub(" hub","",h), "vs Non-hub"),
             n_hub = length(a), med_hub = round(median(a),3), med_nonhub = round(median(b),3),
             W = as.numeric(w$statistic), p = w$p.value, A_prob_sup = round(A,3), r_rb = round(2*A-1,3),
             higher = ifelse(A > .5, sub(" hub","",h), "Non-hub"))
}
res <- rbindlist(lapply(names(covlabs), function(cvn) rbind(wt("Shared hub", cvn), wt("Specific hub", cvn))))
res[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","ns")))]
fwrite(res, file.path(outdir, "covariate_violin.csv"))
cat("\n===== Wilcoxon rank-sum: hub vs Non-hub, per covariate =====\n")
print(res[, .(covariate, comparison, n_hub, med_hub, med_nonhub, p = signif(p,3), A_prob_sup, r_rb, higher, stars)])
