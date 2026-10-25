rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

moldir <- "data/6c-hub_position/cis_trans_pair"
hubdir <- "data/onek1k"
outdir <- "results/others"
CT  <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
rd  <- function(f) unique(trimws(readLines(f)))

GU <- rd(file.path(hubdir, "GECO_universe.txt"))
methods <- list(
  GECO   = list(Bulk = rd(file.path(hubdir, "GECO_Shared.txt")),
                Pseudobulk = rd(file.path(hubdir, "GECO_Specific.txt"))),
  WGCNA  = list(Bulk = rd("data/others/WGCNA_hubs/Shared.txt"),
                Pseudobulk = unique(unlist(lapply(CT, function(ct)
                  rd(file.path("data/others/WGCNA_hubs", paste0(ct, ".txt"))))))),
  glasso = list(Bulk = rd("data/others/glasso_hubs/hubs_Shared.txt"),
                Pseudobulk = unique(unlist(lapply(CT, function(ct)
                  rd(file.path("data/others/glasso_hubs",
                               paste0("hubs_", ct, ".txt"))))))))
MTH  <- names(methods)
MLAB <- c(GECO="GECO", WGCNA="WGCNA", glasso="Standard glasso")

drop_bidir <- function(d) {
  k <- paste(d$source, d$target, sep="|"); rk <- paste(d$target, d$source, sep="|")
  n <- sum(rk %in% k)
  if (n) cat(sprintf("      dropped %d rows in %d bidirectional link(s)\n", n, n/2))
  d[!(rk %in% k)]
}
qtl_loader <- function(f, uni) function() {
  U <- intersect(rd(file.path(moldir, uni)), GU)
  d <- fread(file.path(moldir, f), header=TRUE); setnames(d, 1:2, c("source","target"))
  d <- unique(d[source %in% U & target %in% U & source != target, .(source, target)])
  list(U=U, d=drop_bidir(d))
}
pertseq_loader <- function() {
  pb <- "data/6c-hub_position/perturb_seq"
  U <- intersect(fread(file.path(pb, "aggregate_per_intersection_gene_metrics.tsv"))$gene_ensg, GU)
  e <- fread(file.path(pb, "pertseq_union_network_4col.tsv"),
             select=c("gene_a_ensg","gene_b_ensg","has_dir_a_to_b","has_dir_b_to_a"))
  e <- e[xor(has_dir_a_to_b, has_dir_b_to_a)]
  d <- data.table(source = fifelse(e$has_dir_a_to_b, e$gene_a_ensg, e$gene_b_ensg),
                  target = fifelse(e$has_dir_a_to_b, e$gene_b_ensg, e$gene_a_ensg))
  list(U=U, d=unique(d[source %in% U & target %in% U & source != target]))
}
datasets <- list(
  list(lab="OneK1K",      get=qtl_loader("onek1k.txt",  "onek1k_universe.txt")),
  list(lab="eQTLGen",     get=qtl_loader("eQTLGen.txt", "eQTLGen_universe_pair.txt")),
  list(lab="Perturb-seq", get=pertseq_loader))
DSORD <- vapply(datasets, function(s) s$lab, character(1))

res <- rbindlist(lapply(datasets, function(s) {
  L <- s$get(); U <- L$U; d <- L$d
  cat(sprintf("  %-12s universe %5d | directed links %7d\n", s$lab, length(U), nrow(d)))
  grid <- CJ(mth=MTH, gp=c("Bulk","Pseudobulk"), sorted=FALSE)
  rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    H <- intersect(methods[[grid$mth[i]]][[grid$gp[i]]], U)
    x <- as.integer(d$source %in% H); y <- as.integer(d$target %in% H)
    n11 <- sum(x & y); n10 <- sum(x & !y); n01 <- sum(!x & y); n00 <- sum(!x & !y)
    tt <- t.test(x, y, paired=TRUE)
    bt <- if (n10 + n01 > 0) binom.test(n10, n10 + n01, 0.5) else NULL
    data.table(dataset=s$lab, method=grid$mth[i], group=grid$gp[i],
               nPair=nrow(d), nHub=length(H),
               n_src_only=n10, n_tgt_only=n01, n_both=n11, n_neither=n00,
               n_disc=n10+n01,
               prop_src = if (is.null(bt)) NA_real_ else unname(bt$estimate),
               pci_lo   = if (is.null(bt)) NA_real_ else bt$conf.int[1],
               pci_hi   = if (is.null(bt)) NA_real_ else bt$conf.int[2],
               mean_src=mean(x), mean_tgt=mean(y), diff=mean(x)-mean(y),
               t=unname(tt$statistic), ci_lo=tt$conf.int[1], ci_hi=tt$conf.int[2],
               p_t=tt$p.value, p_mcnemar=if (is.null(bt)) NA_real_ else bt$p.value)
  }))
}))
res[, lean  := fifelse(is.na(diff), NA_character_, fifelse(diff > 0, "regulator", "regulatee"))]
res[, stars := fifelse(is.na(p_mcnemar), "",
                fifelse(p_mcnemar<.001,"***", fifelse(p_mcnemar<.01,"**",
                fifelse(p_mcnemar<.05,"*",""))))]
fwrite(res[order(match(dataset, DSORD), match(method, MTH), group)],
       file.path(outdir, "cis_trans_pair_position_3methods.csv"))
