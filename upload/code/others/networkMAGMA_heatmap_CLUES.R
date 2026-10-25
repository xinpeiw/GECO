rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

dir0   <- "data/7e_networkMAGMA"
outd   <- "results/others"
dir.create(outd, showWarnings=FALSE, recursive=TRUE)
resfile <- "data/7e_networkMAGMA/networkMAGMA_78traits_CLUES_mega_results.txt"
GRPS   <- c("Immune","Non-immune")

map <- fread(file.path(dir0, "networkMAGMA_78traits_trait_immune_group_mapping.csv"))

if (!"immune_binary" %in% names(map))
  stop("mapping file has no `immune_binary` column")
map <- map[, .(column_order, trait_code, trait_name, immune_binary)]
stopifnot(nrow(map) == 78L, !anyDuplicated(map$trait_code),
          setequal(map$immune_binary, GRPS))
map[, grp := factor(immune_binary, levels=GRPS)]
setorder(map, grp, column_order)
cat(sprintf("groups: %s (all %d traits kept)\n",
            paste(sprintf("%s %d", GRPS, map[, .N, by=grp][match(GRPS, grp), N]), collapse=" | "),
            nrow(map)))

dt <- fread(resfile, colClasses=list(character="trait"))[, .(network, trait, p_value, Z_perm)]
cat(sprintf("read %d rows | networks: %s\n", nrow(dt),
            paste(sort(unique(dt$network)), collapse=",")))
dt <- dt[trait %in% map$trait_code]
dt <- merge(dt, map, by.x="trait", by.y="trait_code")
stopifnot(nrow(dt) == 8L * nrow(map))
cat(sprintf("kept %d rows (%d traits x 8 networks)\n", nrow(dt), nrow(map)))

dt[, net     := fifelse(network == "alpha", "Shared", network)]
dt[, neglogp := -log10(p_value)]
dt[, q       := p.adjust(p_value, "BH"), by=grp]
dt[, star    := fifelse(q<.01, "***", fifelse(q<.05, "**", fifelse(q<.1, "*", "")))]

net_lv <- c("Shared","T4","T8","NK","B","cM","ncM","cDC")
stopifnot(setequal(dt$net, net_lv))
dt[, net        := factor(net, levels=net_lv)]
dt[, trait_name := factor(trait_name, levels=rev(map$trait_name))]

fwrite(dt[order(grp, column_order, net),
          .(group=grp, trait, trait_name, network=net, Z_perm, p_value, q, star)],
       file.path(outd, "networkMAGMA_heatmap_CLUES.csv"))
cat(sprintf("FDR<0.1 cells: %d of %d | Shared column: %d\n",
            dt[q < .1, .N], nrow(dt), dt[q < .1 & net == "Shared", .N]))
