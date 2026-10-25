rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

indir  <- "data/5b-ldsc"
outdir <- "results/5b-ldsc"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

tau <- fread(file.path(indir, "sldsc_tau.tsv"))
setnames(tau, "FinnGen ID", "trait")
stopifnot(all(c("hub","Category","tau_star","tau_star_se","h2_z") %in% names(tau)))
tau <- tau[hub %in% c("Shared","Specific")]
cat(sprintf("read %d rows | %d traits x %d hubs\n", nrow(tau), uniqueN(tau$trait), uniqueN(tau$hub)))

ivw <- function(d) {
  w  <- 1/d$tau_star_se^2
  mu <- sum(w*d$tau_star)/sum(w); se <- 1/sqrt(sum(w))
  list(n = nrow(d), tau_star = mu, se = se, p = 2*pnorm(-abs(mu/se)))
}
d <- tau[, ivw(.SD), by=.(hub, Category)]
d[, group := factor(paste("GECO", hub), levels=c("GECO Shared","GECO Specific"))]
d[, cat   := factor(fifelse(Category=="Non-Immune","Non-immune", Category),
                    levels=rev(c("Immune","Non-immune")))]
d[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
setorder(d, cat, group)

up <- fread(file.path(indir, "sldsc_meta.tsv"))[, .(hub, Category, up_tau=tau_star, up_se=se)]
m  <- merge(d[, .(hub, Category, tau_star, se)], up, by=c("hub","Category"))
cat(sprintf("max |diff| vs sldsc_meta.tsv:  tau %.2e   se %.2e\n",
            max(abs(m$tau_star-m$up_tau)), max(abs(m$se-m$up_se))))
stopifnot(max(abs(m$tau_star-m$up_tau)) < 1e-10, max(abs(m$se-m$up_se)) < 1e-10)

cat("===== sLDSC meta-analysis (tau*, hub annotation vs baseline model) =====\n")
print(d[, .(Category, n, hub, tau_star=round(tau_star,4), se=round(se,4), p=signif(p,2), stars)], row.names=FALSE)
fwrite(d[, .(Category, n, hub, tau_star, se, p, stars)], file.path(outdir, "sLDSC.csv"))
