rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

COND   <- "mega"

outdir <- "results/7e_networkMAGMA"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

MAPFILE <- "data/7e_networkMAGMA/networkMAGMA_78traits_trait_immune_group_mapping.csv"
stopifnot(file.exists(MAPFILE))
map <- fread(MAPFILE)

resfile <- "data/7e_networkMAGMA/networkMAGMA_78traits_CLUES_mega_results.txt"
dt <- fread(resfile, colClasses=list(character="trait"))[, .(network, trait, p_value)]
cat(sprintf("[%s] read %d rows | networks: %s\n", COND, nrow(dt),
            paste(sort(unique(dt$network)), collapse=",")))
dt[, net := fifelse(network=="alpha","Shared", network)]
dt <- merge(dt, map, by.x="trait", by.y="trait_code")
stopifnot(nrow(dt) == 624L)
dt[, neglogp := -log10(p_value)]
dt[, grp  := factor(fifelse(net=="Shared","Shared","Specific"), levels=c("Shared","Specific"))]

if (!"immune_binary" %in% names(map))
  stop("mapping has no `immune_binary` column -- run GECO/analysis_plot/",
       "networkMAGMA_78traits/geco/add_immune_binary.R first")
dt[, tier := factor(immune_binary, levels=c("Immune","Non-immune"))]

wt <- dt[, .(p    = wilcox.test(neglogp[grp=="Shared"], neglogp[grp=="Specific"])$p.value,
             nS   = sum(grp=="Shared"),  nP = sum(grp=="Specific"),
             medS = median(neglogp[grp=="Shared"]), medP = median(neglogp[grp=="Specific"])), by=tier]
wt[, lab := sprintf("P = %s", formatC(p, format="g", digits=2))]
fwrite(wt[order(tier)], file.path(outdir, "networkMAGMA_CLUES.csv"))
cat(sprintf("\n===== CLUES %s: Shared vs Specific -log10(P), Immune / Non-immune (78 traits) =====\n", COND))
print(wt[order(tier)])
