rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

ddir   <- "data/5f-RS"
hubdir <- "data/others"
outdir <- "results/others"
CT <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(trimws(readLines(f)))

x0 <- fread(file.path(ddir, "immune_TI_pairs.tsv"))
HUBS <- list(
  `WGCNA Pseudobulk`  = unique(unlist(lapply(CT, function(ct) rd(file.path(hubdir, "WGCNA_hubs", paste0(ct, ".txt")))))),
  `WGCNA Bulk`        = rd(file.path(hubdir, "WGCNA_hubs/Shared.txt")),
  `glasso Pseudobulk` = unique(unlist(lapply(CT, function(ct) rd(file.path(hubdir, "glasso_hubs", paste0("hubs_", ct, ".txt")))))),
  `glasso Bulk`       = rd(file.path(hubdir, "glasso_hubs/hubs_Shared.txt")))

rs_row <- function(y, v, ref) { x1 <- sum(y & v); n1 <- sum(v); x0 <- sum(y & ref); n0 <- sum(ref)
  r <- (x1/n1)/(x0/n0); se <- sqrt(1/x1 - 1/n1 + 1/x0 - 1/n0)
  data.table(x_yes=x1, n_yes=n1, x_no=x0, n_no=n0, rs=r, rs_l95=exp(log(r)-1.96*se), rs_u95=exp(log(r)+1.96*se)) }
TR <- c("pre-I"="succ_p_1", "I-II"="succ_1_2", "II-III"="succ_2_3", "III-launch"="succ_3_a", "I-launch"="cumulative")

res <- rbindlist(lapply(names(TR), function(tn) {
  if (tn == "I-launch") { x <- x0[!is.na(succ_1_2) | !is.na(succ_2_3) | !is.na(succ_3_a)]; y <- !is.na(x$succ_3_a) & x$succ_3_a == TRUE
  } else { x <- x0[!is.na(x0[[TR[tn]]])]; y <- x[[TR[tn]]] == TRUE }
  g <- x$genetic_support
  rbindlist(c(lapply(names(HUBS), function(hn) { h <- x$ensg %in% HUBS[[hn]]
    rbind(cbind(data.table(transition=tn, hub_set=hn, comparison="hub vs all other targets"),                 rs_row(y, h, !h)),
          cbind(data.table(transition=tn, hub_set=hn, comparison="genetic & hub vs no genetic support"),     rs_row(y, g & h, !g)),
          cbind(data.table(transition=tn, hub_set=hn, comparison="genetic & hub vs neither"),                rs_row(y, g & h, !g & !h)),
          cbind(data.table(transition=tn, hub_set=hn, comparison="within genetic support: hub vs non-hub"), rs_row(y, g & h, g & !h))) }),
    list(cbind(data.table(transition=tn, hub_set="(reference) genetic support, all", comparison="genetic vs no genetic support"), rs_row(y, g, !g))))) }))
res[, n_pairs_transition := n_yes + n_no]
res[comparison == "within genetic support: hub vs non-hub",
    p_fisher := mapply(function(a, b, c, d) if (b == 0 || d == 0) NA_real_ else fisher.test(matrix(c(a, b - a, c, d - c), 2))$p.value, x_yes, n_yes, x_no, n_no)]
res[, stars := fifelse(is.na(p_fisher), "", fifelse(p_fisher < .001, "***", fifelse(p_fisher < .01, "**", fifelse(p_fisher < .05, "*", ""))))]
fwrite(res, file.path(outdir, "relative_success_transitions_altmethods.csv"))
print(res[, .(transition, hub_set, comparison, lab=sprintf("%d/%d vs %d/%d", x_yes, n_yes, x_no, n_no),
              rs=round(rs,2), l=round(rs_l95,2), u=round(rs_u95,2), p_fisher=signif(p_fisher,2), stars)], nrows=Inf)
