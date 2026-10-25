rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

evo    <- "data/3cd-evolution"
hubdir <- "data/onek1k"
pmdir  <- "data/onek1k_precision"
outdir <- "results/others"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
rd  <- function(f) unique(trimws(readLines(f)))
ct  <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
PCT <- c(0.01, 0.05, 0.10, 0.15, 0.20); PLAB <- sprintf("top %d%%", round(100*PCT))
U   <- rd(file.path(hubdir, "GECO_universe.txt"))

DEG <- lapply(c("alpha", ct), function(n) {
  m <- readRDS(file.path(pmdir, paste0(n, ".rds")))
  d <- data.table(gene = rownames(m), degree = rowSums(m != 0) - 1L); rm(m); gc(FALSE)
  d[degree > 0]
}); names(DEG) <- c("alpha", ct)
hubs_p <- function(net, p) intersect(DEG[[net]][degree > quantile(degree, 1 - p), gene], U)
HUB <- lapply(PCT, function(p) { l <- lapply(ct, hubs_p, p = p); names(l) <- ct; list(Shared = hubs_p("alpha", p), ct = l) })
names(HUB) <- PLAB
stopifnot(setequal(HUB[["top 5%"]]$Shared, intersect(rd(file.path(hubdir, "GECO_Shared.txt")), U)),
          all(vapply(ct, function(c) setequal(HUB[["top 5%"]]$ct[[c]], intersect(rd(file.path(hubdir, paste0("hubs_", c, ".txt"))), U)), logical(1))))
cat("check: the 5% rule reproduces GECO_Shared.txt and the 7 hubs_<ct>.txt exactly\n")

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
ALAB <- vapply(anns, function(a) a$lab, character(1))

pc <- list(); res <- list()
for (a in anns) {
  rdir <- file.path(evo, a$dir)
  Gc <- lapply(ct, function(c) intersect(rd(file.path(rdir, paste0(c, "_enhanceronly.txt"))), U)); names(Gc) <- ct
  Gunion <- Reduce(union, Gc)
  for (h in PLAB) {
    HUBc <- HUB[[h]]$ct; SHARED <- HUB[[h]]$Shared
    bb <- sapply(ct, function(c) { r <- fit_base(U %in% HUBc[[c]], U %in% Gc[[c]]); c(r$b, r$se, r$ov) })
    aa <- sapply(ct, function(c) { r <- fit_adj (U %in% HUBc[[c]], U %in% Gc[[c]]); c(r$b, r$se) })
    for (i in seq_along(ct)) pc[[length(pc)+1]] <- data.table(
      threshold = h, annotation = a$lab, celltype = ct[i], nHub = length(HUBc[[ct[i]]]), nSet = length(Gc[[ct[i]]]), overlap = bb[3,i],
      base_l2 = bb[1,i]/log(2), base_se = bb[2,i]/log(2), adj_l2 = aa[1,i]/log(2), adj_se = aa[2,i]/log(2))
    mb <- meta_fe(bb[1,], bb[2,]); ma <- meta_fe(aa[1,], aa[2,])
    sb <- fit_base(U %in% SHARED, U %in% Gunion); sa <- fit_adj(U %in% SHARED, U %in% Gunion)
    res[[length(res)+1]] <- rbindlist(list(
      data.table(threshold=h, annotation=a$lab, group="GECO Specific", mode="Baseline", l2=mb$b/log(2), se=mb$se/log(2), p=mb$p,        k=mb$k, nSet=length(Gunion)),
      data.table(threshold=h, annotation=a$lab, group="GECO Specific", mode="Adjusted", l2=ma$b/log(2), se=ma$se/log(2), p=ma$p,        k=ma$k, nSet=length(Gunion)),
      data.table(threshold=h, annotation=a$lab, group="GECO Shared",   mode="Baseline", l2=sb$b/log(2), se=sb$se/log(2), p=zp(sb$b,sb$se), k=1L, nSet=length(Gunion)),
      data.table(threshold=h, annotation=a$lab, group="GECO Shared",   mode="Adjusted", l2=sa$b/log(2), se=sa$se/log(2), p=zp(sa$b,sa$se), k=1L, nSet=length(Gunion))))
  }
}
res <- rbindlist(res); pc <- rbindlist(pc)
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p), "", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[order(match(threshold, PLAB), match(annotation, ALAB), group, mode)], file.path(outdir, "evolution_rE2G_enrichment_toppct.csv"))
fwrite(pc,  file.path(outdir, "evolution_rE2G_enrichment_toppct_percelltype.csv"))

pub <- fread("results/3cd-evolution/evolution_rE2G_enrichment.csv")
mm  <- merge(res[threshold == "top 5%", .(annotation, group, mode, new = l2, new_k = k)],
             pub[, .(annotation, group, mode, old = l2, old_k = k)], by = c("annotation","group","mode"))
stopifnot(nrow(mm) == 12L, all(mm$new_k == mm$old_k), max(abs(mm$new - mm$old), na.rm = TRUE) < 1e-8,
          identical(is.na(mm$new), is.na(mm$old)))
cat(sprintf("\ncheck: top 5%% rows reproduce evolution_rE2G_enrichment.csv (%d cells, max |diff| = %.1e)\n", nrow(mm), max(abs(mm$new - mm$old), na.rm = TRUE)))
cat("\n===== log2 OR by hub threshold (k = estimable cell types in the meta) =====\n")
print(dcast(res[, .(threshold, annotation, group, mode, v = sprintf("%5.2f%s (k%d)", l2, stars, k))],
            annotation + group + mode ~ factor(threshold, levels = PLAB), value.var = "v")[order(match(annotation, ALAB), group, mode)], nrows = Inf)
