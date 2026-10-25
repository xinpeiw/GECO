rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

outdir <- "results/5b-ldsc"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

tau   <- fread("data/5b-ldsc/sensitivity/sldsc_tau.tsv")
panel <- fread("data/5b-ldsc/sensitivity/FinnGen_phenotypes_final.csv")[
           , .(`FinnGen ID`, Category, primary_h2z_gt4)]
tau   <- merge(tau, panel, by = "FinnGen ID")

tau   <- tau[hub_set %in% c("Shared", "Specific")]

ivw <- function(d) {
  w <- 1/d$tau_star_se^2
  mu <- sum(w*d$tau_star)/sum(w); se <- 1/sqrt(sum(w))
  data.table(n = nrow(d), tau_star = mu, se = se, p = 2*pnorm(-abs(mu/se)))
}

MODELS <- c(primary         = "Primary",
            baselineLD2.2   = "Baseline-LD v2.2",
            window500kb     = "Window +/-500 kb",
            MHCincl         = "MHC included",
            top210          = "Top 210 hubs",
            top420          = "Top 420 hubs")
mod <- tau[ldsc_settings %in% names(MODELS) & primary_h2z_gt4 == TRUE,
           ivw(.SD), by = .(setting = ldsc_settings, hub = hub_set, Category)]
mod[, panel := MODELS[setting]]

thr <- rbindlist(lapply(c(0, 7), function(z)
  tau[ldsc_settings == "primary" & h2_z > z,
      ivw(.SD), by = .(hub = hub_set, Category)][
      , panel := sprintf("Primary, h2 z > %d", z)][]))

d <- rbind(mod[, .(panel, hub, Category, n, tau_star, se, p)],
           thr[, .(panel, hub, Category, n, tau_star, se, p)])
d[, panel := factor(panel, levels = c(unname(MODELS), "Primary, h2 z > 0", "Primary, h2 z > 7"))]
d[, group := factor(hub, levels = c("Shared","Specific"))]
d[, cat   := factor(fifelse(Category == "Non-Immune", "Non-immune", Category),
                    levels = rev(c("Immune","Non-immune")))]
d[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
setorder(d, panel, cat, group)
cat("===== sLDSC sensitivity analyses (tau*, IVW meta) =====\n")
print(d[, .(panel, Category, n, hub, tau_star=round(tau_star,4), se=round(se,4),
            p=signif(p,2), stars)], nrows=Inf)
fwrite(d[, .(panel, Category, n, hub, tau_star, se, p, stars)],
       file.path(outdir, "sLDSC_sensitivity.csv"))
