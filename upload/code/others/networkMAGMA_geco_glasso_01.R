rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

outd <- "results/others"
dir.create(outd, showWarnings=FALSE, recursive=TRUE)
resfile <- "data/others/networkMAGMA_78traits_geco_glasso_01_results.txt"

NET  <- c("GECO_Shared","GECO_Specific","glasso_Bulk","glasso_Specific")
GRPS <- c("GECO Shared","GECO Specific","glasso Bulk","glasso Pseudobulk")
LAB  <- setNames(GRPS, NET)

map <- fread("data/6d-networkMAGMA/networkMAGMA_78traits_trait_immune_group_mapping.csv")
stopifnot(nrow(map) == 78L, !anyDuplicated(map$trait_code))
if (!"immune_binary" %in% names(map))
  stop("mapping file has no `immune_binary` column -- run geco/add_immune_binary.R first")

dt <- fread(resfile, colClasses = list(character = "trait"))[, .(network, trait, N, S0, p_value)]

stopifnot(nrow(dt) == 312L, setequal(dt$network, NET),
          !anyDuplicated(dt[, .(network, trait)]))

nper  <- dt[, .(k = .N), by=trait]
drop  <- nper[k < 4L, trait]
if (length(drop)) {
  miss <- dt[trait %in% drop, .(have = paste(sort(network), collapse=", ")), by=trait]
  miss[, missing := vapply(trait, function(tr) paste(setdiff(NET, dt[trait == tr, network]),
                                                     collapse=", "), character(1))]
  cat(sprintf("\n*** dropping %d of %d traits that lack a complete set of 4 networks ***\n",
              length(drop), uniqueN(dt$trait)))
  print(miss[, .(trait, missing)], row.names=FALSE)
  fwrite(miss, file.path(outd, "networkMAGMA_geco_glasso_01_dropped_traits.csv"))
  dt <- dt[!trait %in% drop]
}
stopifnot(dt[, .N, by=trait][, all(N == 4L)])

stopifnot(dt[, uniqueN(N), by=trait][, all(V1 == 1L)])
cat(sprintf("rows %d | traits %d | genes per trait %s (identical across the 4 networks)\n",
            nrow(dt), uniqueN(dt$trait),
            paste(sort(unique(dt$N)), collapse="/")))

cat("\nnetwork density (S0 = sum of the 0/1 adjacency = 2 x edges; median over traits):\n")
print(dt[, .(S0 = median(S0), edges = median(S0)/2), by=network][match(NET, network)],
      row.names=FALSE)

stopifnot(all(dt$trait %in% map$trait_code))
ntr <- uniqueN(dt$trait)
dt <- merge(dt, map, by.x="trait", by.y="trait_code")
stopifnot(nrow(dt) == 4L * ntr)
dt[, tier    := factor(immune_binary, levels=c("Immune","Non-immune"))]
dt[, grp     := factor(LAB[network], levels=GRPS)]
dt[, neglogp := -log10(p_value)]
cat(sprintf("\np range: [%.3g, %.3g] | permutation floor 1/(1e5+1) = %.3g -> max -log10(P) = %.2f\n",
            min(dt$p_value), max(dt$p_value), 1/(1e5+1), max(dt$neglogp)))

w <- dcast(dt, trait + tier ~ grp, value.var="neglogp")
stopifnot(!anyNA(w), nrow(w) == ntr)
cat("\ntraits analysed per tier:\n"); print(w[, .N, by=tier], row.names=FALSE)

pairs <- as.data.table(t(combn(GRPS, 2))); setnames(pairs, c("a","b"))
wt <- rbindlist(lapply(levels(dt$tier), function(tt) {
  wi <- w[tier == tt]
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    x <- wi[[pairs$a[i]]]; y <- wi[[pairs$b[i]]]
    data.table(tier=tt, a=pairs$a[i], b=pairs$b[i], n_traits=nrow(wi),
               med_a=median(x), med_b=median(y),
               n_a_higher=sum(x > y), n_b_higher=sum(y > x),
               p_paired=wilcox.test(x, y, paired=TRUE)$p.value)
  }))
}))
fwrite(wt, file.path(outd, "networkMAGMA_geco_glasso_01_paired.csv"))
fwrite(dt[order(tier, trait, grp), .(trait, tier, network, grp, N, S0, p_value, neglogp)],
       file.path(outd, "networkMAGMA_geco_glasso_01_values.csv"))

cat("\n===== median -log10(P), and how many traits clear P < 0.05 =====\n")
print(dt[, .(median = round(median(neglogp), 3),
             n_sig  = sum(p_value < 0.05), n = .N), by=.(tier, grp)][order(tier, grp)],
      row.names=FALSE)
cat("\n===== all six pairwise, PAIRED Wilcoxon signed-rank =====\n")
print(wt[order(tier), .(tier, a, b, n_traits, med_a=round(med_a,3), med_b=round(med_b,3),
                        a_higher=n_a_higher, b_higher=n_b_higher,
                        p=signif(p_paired,3))], row.names=FALSE)
