rm(list=ls())
suppressMessages(library(data.table))

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
CELLORD <- c("Shared","CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")
PCOLS   <- c("p_global","p_strat_10","p_strat_20","p_strat_50","p_strat_100","p_strat_200")
LABS    <- c("Global","Strat 10","Strat 20","Strat 50","Strat 100","Strat 200")

d <- fread("data/others/sim_power_binary.txt")
summ <- d[, c(list(n_rep = .N), lapply(.SD, function(p) mean(p < 0.05))), by = .(network, q), .SDcols = PCOLS]
long <- melt(summ, id.vars = c("network","q","n_rep"), measure.vars = PCOLS, variable.name = "null", value.name = "power")
long[, null    := factor(null, levels = PCOLS, labels = LABS)]
long[, network := factor(network, levels = CELLORD)]
setorder(long, network, null, q)
fwrite(long, file.path(outdir, "sim_power_binary.txt"), sep = "\t")
print(dcast(long, network + null ~ q, value.var = "power"), nrows = Inf)
