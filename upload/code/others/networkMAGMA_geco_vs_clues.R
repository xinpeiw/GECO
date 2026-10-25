rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

outd <- "results/others"
dir.create(outd, showWarnings=FALSE, recursive=TRUE)

map <- fread("data/6d-networkMAGMA/networkMAGMA_78traits_trait_immune_group_mapping.csv")
stopifnot(nrow(map) == 78L, !anyDuplicated(map$trait_code))
if (!"immune_binary" %in% names(map))
  stop("mapping file has no `immune_binary` column")

load_res <- function(resfile, tag) {
  d <- fread(resfile, colClasses = list(character = "trait"))[, .(network, trait, p_value)]
  d <- d[trait %in% map$trait_code]
  cat(sprintf("  %-8s %4d mapped rows | networks: %s\n",
              tag, nrow(d), paste(sort(unique(d$network)), collapse=",")))
  d[, .(network, trait, neglogp = -log10(p_value))]
}
cat("reading:\n")
G <- load_res("data/6d-networkMAGMA/networkMAGMA_78traits_results.txt",            "OneK1K")
C <- load_res("data/7e_networkMAGMA/networkMAGMA_78traits_CLUES_mega_results.txt", "CLUES")

PAIRS <- data.table(
  onek1k = c("CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem","alpha"),
  clues  = c("T4",   "T4",   "T8",   "T8",   "NK","B",  "B",   "alpha"))
PAIRS[, lab := sprintf("%s  vs  %s", fifelse(onek1k == "alpha", "Shared", onek1k),
                                     fifelse(clues  == "alpha", "Shared", clues))]
stopifnot(all(PAIRS$onek1k %in% G$network), all(PAIRS$clues %in% C$network))
cat(sprintf("\nCLUES networks with no OneK1K counterpart: %s\n",
            paste(setdiff(unique(C$network), PAIRS$clues), collapse=", ")))

d <- rbindlist(lapply(seq_len(nrow(PAIRS)), function(i) {
  g <- G[network == PAIRS$onek1k[i], .(trait, onek1k = neglogp)]
  c <- C[network == PAIRS$clues[i],  .(trait, clues  = neglogp)]
  m <- merge(g, c, by = "trait")
  stopifnot(nrow(m) == nrow(map))
  m[, `:=`(pair = PAIRS$lab[i], onek1k_net = PAIRS$onek1k[i], clues_net = PAIRS$clues[i])]
}))
d <- merge(d, map[, .(trait_code, trait_name, immune_binary)], by.x="trait", by.y="trait_code")
d[, pair := factor(pair, levels = PAIRS$lab)]
d[, tier := factor(immune_binary, levels = c("Immune","Non-immune"))]
stopifnot(nrow(d) == nrow(PAIRS) * nrow(map))

agree <- function(x) x[, .(n = .N,
    r_pearson   = cor(onek1k, clues),
    rho_spear   = cor(onek1k, clues, method = "spearman"),
    med_onek1k  = median(onek1k), med_clues = median(clues),
    n_onek1k_hi = sum(onek1k > clues), n_clues_hi = sum(clues > onek1k),
    p_paired    = wilcox.test(onek1k, clues, paired = TRUE)$p.value)]
st  <- d[, agree(.SD), by = pair]
stt <- d[, agree(.SD), by = .(pair, tier)]
fwrite(st,  file.path(outd, "networkMAGMA_geco_vs_clues.csv"))
fwrite(stt, file.path(outd, "networkMAGMA_geco_vs_clues_bytier.csv"))
fwrite(d[order(pair, trait)], file.path(outd, "networkMAGMA_geco_vs_clues_values.csv"))

cat("\n===== agreement per matched network (all 78 traits) =====\n")
print(st[, .(pair, n, r = round(r_pearson,3), rho = round(rho_spear,3),
             med_OneK1K = round(med_onek1k,2), med_CLUES = round(med_clues,2),
             OneK1K_higher = n_onek1k_hi, CLUES_higher = n_clues_hi,
             p_paired = signif(p_paired,2))], row.names = FALSE)
cat("\n===== the same, split by disease group =====\n")
print(stt[order(pair, tier), .(pair, tier, n, r = round(r_pearson,3), rho = round(rho_spear,3),
             OneK1K_higher = n_onek1k_hi, CLUES_higher = n_clues_hi,
             p_paired = signif(p_paired,2))], row.names = FALSE)
