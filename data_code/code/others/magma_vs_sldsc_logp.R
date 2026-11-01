rm(list=ls())
suppressMessages(library(data.table))

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
NET_ORDER <- c("Shared","CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")

res <- fread("data/6d-networkMAGMA/networkMAGMA_78traits_results.txt",
             colClasses=list(character="trait"))[, .(network, trait, p_value)]
map <- fread("data/6d-networkMAGMA/networkMAGMA_78traits_trait_immune_group_mapping.csv")
M <- merge(res, map[, .(trait = trait_code, trait_name, group = immune_binary)], by = "trait")
M[, network := fifelse(network == "alpha", "Shared", network)]
stopifnot(nrow(M) == 8L * nrow(map))

S <- fread("data/5b-ldsc/sensitivity/sldsc_tau.tsv")[ldsc_settings == "primary" & hub_set %in% NET_ORDER]
S[, z := tau_star / tau_star_se]
S[, p_ldsc := pnorm(-z)]

d <- merge(M[, .(trait, trait_name, group, network, p_magma = p_value)],
           S[, .(trait = `FinnGen ID`, network = hub_set, p_ldsc, h2_z)],
           by = c("trait","network"))
d[, `:=`(network = factor(network, NET_ORDER),
         group   = factor(group, c("Immune","Non-immune")),
         x = -log10(p_magma), y = -log10(p_ldsc))]
fwrite(d[, .(trait, trait_name, group, network, p_magma, p_ldsc, h2_z)],
       file.path(outdir, "magma_vs_sldsc_logp.tsv"), sep = "\t")

rho <- d[, .(n_traits = .N, spearman_rho = cor(x, y, method = "spearman")), by = network][order(network)]
fwrite(rho, file.path(outdir, "magma_vs_sldsc_logp_spearman.tsv"), sep = "\t")
cat("pairs:", nrow(d), " overall Spearman:", d[, cor(x, y, method = "spearman")], "\n")
print(rho, row.names = FALSE)
