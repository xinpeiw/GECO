rm(list=ls())
suppressPackageStartupMessages({ library(data.table); library(readxl) })

dir0 <- "data/4a-TF_HPA"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd   <- function(f) {
  g <- unique(grep("^ENSG", trimws(readLines(f, warn = FALSE)), value = TRUE))
  if (!length(g)) stop(sprintf("%s yielded 0 genes -- unreadable or an unhydrated placeholder", f))
  g
}
CT   <- c("CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")
LIN  <- list(T = c("CD4NC","CD4ET","CD8NC","CD8ET"), NK = "NK", B = c("BIN","BMem"))
U    <- rd(file.path(dir0, "GECOpredicted/alltfs.txt"))
cat(sprintf("alltfs universe: %d\n", length(U)))

predict_sets <- function(xlsx, pooled) {
  sh <- excel_sheets(xlsx)
  stopifnot(all(c(pooled, CT) %in% sh))
  deg <- rbindlist(lapply(c(pooled, CT), function(s) {
    d <- as.data.table(read_excel(xlsx, sheet = s))[, .(gene, degree)]
    d[, net := s][] }))
  deg <- deg[gene %in% U]
  stopifnot(uniqueN(deg$net) == 8L)
  deg[, mx := max(degree), by = gene]
  win <- deg[mx > 0 & degree == mx]
  list(Shared = sort(unique(win[net == pooled, gene])),
       T      = sort(unique(win[net %in% LIN$T,  gene])),
       NK     = sort(unique(win[net %in% LIN$NK, gene])),
       B      = sort(unique(win[net %in% LIN$B,  gene])))
}

gp <- predict_sets("data/others/wglasso_degree.xlsx", "alpha")
for (k in names(gp)) {
  pub <- rd(file.path(dir0, "GECOpredicted", paste0("predicted_", k, ".txt")))
  if (!setequal(gp[[k]], pub))
    stop(sprintf("re-implemented rule does NOT reproduce predicted_%s.txt (%d vs %d genes)",
                 k, length(gp[[k]]), length(pub)))
}
cat("check: re-implemented rule reproduces all four GECOpredicted/predicted_*.txt exactly\n")

preds <- predict_sets("data/others/real_data_naive_glasso_degree.xlsx",
                      "Shared")
cat("\npredicted set sizes (TFs; sets overlap because ties assign to several):\n")
for (k in c("Shared","T","NK","B"))
  cat(sprintf("  %-7s GECO %4d | glasso %4d\n", k, length(gp[[k]]), length(preds[[k]])))
for (k in names(preds))
  writeLines(preds[[k]], file.path(outdir, paste0("TF_HPA_predicted_enrichment_glasso_predicted_", k, ".txt")))

tfs  <- list("B cell TF" = rd(file.path(dir0,"TFset/B cell TF.txt")),
             "NK cell TF" = rd(file.path(dir0,"TFset/NK cell TF.txt")),
             "T cell TF" = rd(file.path(dir0,"TFset/T cell TF.txt")),
             "Low specificity TF" = rd(file.path(dir0,"TFset/Low specificity TF.txt")))
corr  <- c(B="B cell TF", NK="NK cell TF", T="T cell TF", Shared="Low specificity TF")
covdt <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
covdt[fread("data/onek1k/beta_var.txt"), on="gene", lbv:=log(i.beta_var)]
covdt <- covdt[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]

fit_base <- function(pr, tf) {
  ih <- U %in% pr; is <- U %in% tf
  a<-sum(ih&is); b<-sum(ih&!is); c<-sum(!ih&is); d<-sum(!ih&!is); ov<-a
  if (min(a,b,c,d)==0) { a<-a+.5;b<-b+.5;c<-c+.5;d<-d+.5 }
  l2<-log2((a*d)/(b*c)); se<-sqrt(1/a+1/b+1/c+1/d)/log(2)
  list(l2=l2, se=se, p=2*pnorm(-abs(l2/se)), ov=ov)
}

SETS <- list(GECO = gp, glasso = preds)
MLAB <- c(GECO = "GECO", glasso = "Standard glasso")
grid <- CJ(m = names(SETS), pn = c("B","NK","T","Shared"), tn = names(tfs), sorted = FALSE)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  b <- fit_base(SETS[[grid$m[i]]][[grid$pn[i]]], tfs[[grid$tn[i]]])
  data.table(method = grid$m[i], pred = grid$pn[i], tf = grid$tn[i],
             corresponding = (corr[[grid$pn[i]]] == grid$tn[i]),
             nPred = length(SETS[[grid$m[i]]][[grid$pn[i]]]),
             nTF = length(tfs[[grid$tn[i]]]), overlap = b$ov,
             l2 = b$l2, se = b$se, p = b$p)
}))
res[, OR := 2^l2]
res[, star := fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*","")))]
fwrite(res, file.path(outdir, "TF_HPA_predicted_enrichment_glasso.csv"))
