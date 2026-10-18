rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })
resfile <- "data/6f-networkMAGMA_perturbseq/networkMAGMA_78traits_perturbseq_results.txt"
outdir <- "results/6f-networkMAGMA_perturbseq"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

map <- fread("data/6f-networkMAGMA_perturbseq/networkMAGMA_78traits_trait_immune_group_mapping.csv")

stopifnot(nrow(map) == 78L, !anyDuplicated(map$trait_code))
if (!"immune_binary" %in% names(map))
  stop("mapping file has no `immune_binary` column -- run geco/add_immune_binary.R first
  (make_trait_mapping_78.R overwrites the file and strips it; see that script's header)")

keep <- c(pertseq_union="Perturb-seq", pertseq_shared_overlap="Shared overlap", pertseq_specific_overlap="Specific overlap")
dt <- fread(resfile, colClasses=list(character="trait"))[, .(network, trait, p_value)]
dt <- dt[network %in% names(keep)]
dt <- merge(dt, map, by.x="trait", by.y="trait_code")

stopifnot(nrow(dt) == length(keep) * nrow(map))
dt[, neglogp := -log10(p_value)]
net_lv <- c("Perturb-seq","Shared overlap","Specific overlap")
tier_lv <- c("Immune","Non-immune")

dt[, tier := factor(immune_binary, levels=tier_lv)]
dt[, net  := factor(keep[network], levels=net_lv)]

ref <- "Perturb-seq"; others <- c("Shared overlap","Specific overlap")
wt <- rbindlist(lapply(tier_lv, function(tt) rbindlist(lapply(others, function(nn) {
  m <- merge(dt[tier==tt & net==nn,  .(trait, a=neglogp)],
             dt[tier==tt & net==ref, .(trait, r=neglogp)], by="trait")
  data.table(tier=tt, net=nn, p=wilcox.test(m$a, m$r, paired=TRUE)$p.value, n=nrow(m),
             med=median(m$a), med_ref=median(m$r), med_diff=median(m$a - m$r),
             dir=ifelse(median(m$a - m$r) >= 0, "up", "down"))
}))))
wt[, stars := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","ns")))]
wt[, lab := sprintf("P = %s", formatC(p, format="g", digits=2))]
fwrite(wt, file.path(outdir, "networkMAGMA_pertseq_overlap.csv"))
cat("\n===== Paired Wilcoxon signed-rank: each overlap vs Perturb-seq reference (per tier) =====\n")
print(wt[, .(tier, net, n, med, med_ref, med_diff, dir, p=signif(p,2), stars)])
