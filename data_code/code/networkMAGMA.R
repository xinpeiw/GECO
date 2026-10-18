rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })
dir0   <- "data/6d-networkMAGMA"
outd   <- "results/6d-networkMAGMA"
dir.create(outd, showWarnings=FALSE, recursive=TRUE)
resfile <- "data/6d-networkMAGMA/networkMAGMA_78traits_results.txt"

map <- fread(file.path(dir0, "networkMAGMA_78traits_trait_immune_group_mapping.csv"))

stopifnot(nrow(map) == 78L, !anyDuplicated(map$trait_code))
if (!"immune_binary" %in% names(map))
  stop("mapping file has no `immune_binary` column -- run geco/add_immune_binary.R first
  (make_trait_mapping_78.R overwrites the file and strips it; see that script's header)")

dt <- fread(resfile, colClasses=list(character="trait"))[, .(network, trait, p_value)]
cat(sprintf("read %d rows | networks: %s\n", nrow(dt),
            paste(sort(unique(dt$network)), collapse=",")))
dt[, net := fifelse(network=="alpha","GECO-shared", network)]

dt <- dt[trait %in% map$trait_code]
dt <- merge(dt, map, by.x="trait", by.y="trait_code")
stopifnot(nrow(dt) == 8L * nrow(map))
dt[, neglogp := -log10(p_value)]
dt[, grp  := factor(fifelse(net=="GECO-shared","GECO-shared","GECO-specific"), levels=c("GECO-shared","GECO-specific"))]
tier_lv <- c("Immune","Non-immune")

dt[, tier := factor(immune_binary, levels=tier_lv)]

wt <- dt[, .(p    = wilcox.test(neglogp[grp=="GECO-shared"], neglogp[grp=="GECO-specific"])$p.value,
             nS   = sum(grp=="GECO-shared"),  nP = sum(grp=="GECO-specific"),
             medS = median(neglogp[grp=="GECO-shared"]), medP = median(neglogp[grp=="GECO-specific"])), by=tier]
wt[, lab := sprintf("P = %s", formatC(p, format="g", digits=2))]
fwrite(wt[order(tier)], file.path(outd, "networkMAGMA.csv"))
cat("\n===== Shared vs Specific -log10(P) within each immune tier (78 traits) =====\n"); print(wt[order(tier)])
