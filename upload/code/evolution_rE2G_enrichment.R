rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

evo    <- "data/3cd-evolution"
hubdir <- "data/onek1k"
outdir <- "results/3cd-evolution"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) unique(trimws(readLines(f)))
ct <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
U  <- rd(file.path(hubdir, "GECO_universe.txt"))
SHARED <- intersect(rd(file.path(hubdir, "GECO_Shared.txt")), U)
HUBc   <- lapply(ct, function(c) intersect(rd(file.path(hubdir, paste0("hubs_", c, ".txt"))), U)); names(HUBc) <- ct

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]

rc  <- fread(file.path(evo, "GECO_gene_recombination_rate.1000G.hg38.tsv"))[, .(gene = gene_id, recomb = recomb_rate_cM_Mb)]
ex  <- fread(file.path(evo, "exon.txt"))[, .(gene = ensg, nexon = n_exons_unique)]
bv  <- fread("data/onek1k/beta_var.txt")[, .(gene, beta_var)]
Ud  <- data.table(gene = U, logE = NA_real_, logL = NA_real_, gc = NA_real_, lrec = NA_real_, lexn = NA_real_, lbv = NA_real_)
Ud[cov, on = "gene", `:=`(logE = i.logE, logL = i.logL, gc = i.gc)]
Ud[rc,  on = "gene", lrec := log1p(i.recomb)]
Ud[ex,  on = "gene", lexn := log(i.nexon)]
Ud[bv,  on = "gene", lbv := log(i.beta_var)]
cc  <- is.finite(Ud$logE) & is.finite(Ud$logL) & is.finite(Ud$gc) & is.finite(Ud$lrec) & is.finite(Ud$lexn) & is.finite(Ud$lbv)
cat(sprintf("covariate-complete (6-cov): %d of %d universe genes\n", sum(cc), length(U)))

fit_base <- function(inhub, inset) {
  a <- sum(inhub & inset); b <- sum(inhub & !inset); c <- sum(!inhub & inset); d <- sum(!inhub & !inset)
  if (a == 0 || b == 0 || c == 0 || d == 0) { a<-a+.5; b<-b+.5; c<-c+.5; d<-d+.5 }
  list(b = log((a*d)/(b*c)), se = sqrt(1/a+1/b+1/c+1/d), ov = sum(inhub & inset))
}
fit_adj <- function(inhub, inset) {
  d <- data.table(y = as.integer(inset[cc]), x = as.integer(inhub[cc]), logE = Ud$logE[cc], logL = Ud$logL[cc],
                  gc = Ud$gc[cc], lexn = Ud$lexn[cc], lrec = Ud$lrec[cc], lbv = Ud$lbv[cc])
  m <- tryCatch(suppressWarnings(glm(y ~ x + logE + logL + gc + lexn + lrec + lbv, d, family = binomial)), error = function(e) NULL)
  if (is.null(m) || !m$converged) return(list(b = NA, se = NA))
  co <- coef(summary(m)); if (!("x" %in% rownames(co))) return(list(b = NA, se = NA))
  se <- co["x","Std. Error"]; if (!is.finite(se) || se > 10) return(list(b = NA, se = NA))
  list(b = co["x","Estimate"], se = se)
}
meta_fe <- function(bs, ses) {
  ok <- is.finite(bs) & is.finite(ses) & ses > 0; bs <- bs[ok]; ses <- ses[ok]
  if (!length(bs)) return(list(b = NA, se = NA, p = NA, k = 0))
  w <- 1/ses^2; bm <- sum(w*bs)/sum(w); sem <- sqrt(1/sum(w))
  list(b = bm, se = sem, p = 2*pnorm(-abs(bm/sem)), k = length(bs))
}
zp <- function(b, se) if (is.na(b) || is.na(se)) NA else 2*pnorm(-abs(b/se))

anns <- list(list(dir = "zooHAR_positive", lab = "zooHAR"),
             list(dir = "HAQER",           lab = "HAQER"),
             list(dir = "akbari",          lab = "Akbari"))

pc <- list(); res <- list()
for (a in anns) {
  rdir <- file.path(evo, a$dir)
  Gc <- lapply(ct, function(c) intersect(rd(file.path(rdir, paste0(c, "_enhanceronly.txt"))), U)); names(Gc) <- ct
  Gunion <- Reduce(union, Gc)

  bb <- sapply(ct, function(c) { r <- fit_base(U %in% HUBc[[c]], U %in% Gc[[c]]); c(r$b, r$se, r$ov) })
  aa <- sapply(ct, function(c) { r <- fit_adj (U %in% HUBc[[c]], U %in% Gc[[c]]); c(r$b, r$se) })
  for (i in seq_along(ct)) pc[[length(pc)+1]] <- data.table(
    annotation = a$lab, celltype = ct[i], nSet = length(Gc[[ct[i]]]), overlap = bb[3,i],
    base_l2 = bb[1,i]/log(2), base_se = bb[2,i]/log(2), adj_l2 = aa[1,i]/log(2), adj_se = aa[2,i]/log(2))
  mb <- meta_fe(bb[1,], bb[2,]); ma <- meta_fe(aa[1,], aa[2,])

  sb <- fit_base(U %in% SHARED, U %in% Gunion); sa <- fit_adj(U %in% SHARED, U %in% Gunion)

  res[[length(res)+1]] <- rbindlist(list(
    data.table(annotation=a$lab, group="GECO Specific", mode="Baseline", l2=mb$b/log(2), se=mb$se/log(2), p=mb$p,        k=mb$k, nSet=length(Gunion)),
    data.table(annotation=a$lab, group="GECO Specific", mode="Adjusted", l2=ma$b/log(2), se=ma$se/log(2), p=ma$p,        k=ma$k, nSet=length(Gunion)),
    data.table(annotation=a$lab, group="GECO Shared",   mode="Baseline", l2=sb$b/log(2), se=sb$se/log(2), p=zp(sb$b,sb$se), k=1L, nSet=length(Gunion)),
    data.table(annotation=a$lab, group="GECO Shared",   mode="Adjusted", l2=sa$b/log(2), se=sa$se/log(2), p=zp(sa$b,sa$se), k=1L, nSet=length(Gunion))))
}
res <- rbindlist(res); pc <- rbindlist(pc)
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p), "", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]

fwrite(res[order(annotation,group,mode)], file.path(outdir, "evolution_rE2G_enrichment.csv"))
fwrite(pc,  file.path(outdir, "evolution_rE2G_enrichment_percelltype.csv"))
cat("\n===== meta/union enrichment (log2 OR) =====\n")
print(res[order(annotation,group,mode), .(annotation, group, mode, nSet_union=nSet, k, l2=round(l2,2), OR=round(OR,2), p=signif(p,2), stars)])
cat("\n----- per-cell-type (Specific), baseline -----\n")
print(pc[, .(annotation, celltype, nSet, overlap, base_l2=round(base_l2,2), adj_l2=round(adj_l2,2))])
