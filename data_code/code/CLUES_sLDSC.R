rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

COND   <- "mega"
indir  <- "data/7a-ldsc"
outdir <- "results/7a-ldsc"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

tau <- fread(file.path(indir, "sldsc_tau.tsv"))
setnames(tau, "FinnGen ID", "trait")
stopifnot(all(c("hub","Category","tau_star","tau_star_se","h2_z") %in% names(tau)))
cat(sprintf("read %d rows | %d traits x %d hubs\n", nrow(tau), uniqueN(tau$trait), uniqueN(tau$hub)))

ivw <- function(d) {
  w  <- 1/d$tau_star_se^2
  mu <- sum(w*d$tau_star)/sum(w); se <- 1/sqrt(sum(w))
  Q  <- sum(w*(d$tau_star - mu)^2); df <- nrow(d) - 1L
  list(n = nrow(d), tau_star = mu, se = se, p = 2*pnorm(-abs(mu/se)),
       Q = Q, Q_p = if (df > 0) pchisq(Q, df, lower.tail=FALSE) else NA_real_,
       I2 = if (df > 0) max(0, (Q - df)/Q)*100 else NA_real_)
}
d <- tau[, ivw(.SD), by=.(hub, Category)]
d[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
setorder(d, Category, hub)
cat("\n===== IVW meta (all 78 traits) =====\n")
print(d[, .(Category, hub, n, tau_star=round(tau_star,4), se=round(se,4),
            p=signif(p,2), stars, I2=round(I2,1), Q_p=signif(Q_p,2))], row.names=FALSE)
fwrite(d, file.path(outdir, "CLUES_sLDSC.csv"))

up <- file.path(indir, "sldsc_meta.tsv")
if (file.exists(up)) {
  u <- fread(up)[, .(hub, Category, up_tau=tau_star, up_se=se)]
  m <- merge(d[, .(hub, Category, tau_star, se)], u, by=c("hub","Category"))
  cat(sprintf("\nmax |diff| vs upstream sldsc_meta.tsv:  tau %.2e   se %.2e\n",
              max(abs(m$tau_star-m$up_tau)), max(abs(m$se-m$up_se))))
}
