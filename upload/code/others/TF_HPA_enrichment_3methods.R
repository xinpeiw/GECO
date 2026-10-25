rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

dir0 <- "data/4a-TF_HPA"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

rd <- function(f) {
  g <- unique(grep("^ENSG", trimws(readLines(f, warn = FALSE)), value = TRUE))
  if (!length(g))
    stop(sprintf("%s yielded 0 genes -- unreadable or empty", f))
  g
}
U    <- rd("data/onek1k/GECO_universe.txt")

hubsrc <- list(
  GECO   = list(bulk = "data/onek1k/GECO_Shared.txt",
                cell = function(c) file.path("data/onek1k", paste0("hubs_", c, ".txt"))),
  WGCNA  = list(bulk = "data/others/WGCNA_hubs/Shared.txt",
                cell = function(c) file.path("data/others/WGCNA_hubs", paste0(c, ".txt"))),
  glasso = list(bulk = "data/others/glasso_hubs/hubs_Shared.txt",
                cell = function(c) file.path("data/others/glasso_hubs",
                                             paste0("hubs_", c, ".txt"))))
MTH  <- names(hubsrc)
MLAB <- c(GECO = "GECO", WGCNA = "WGCNA", glasso = "Standard glasso")
LIN  <- list(Tcell = c("CD4NC","CD4ET","CD8NC","CD8ET"), NKcell = "NK", Bcell = c("BIN","BMem"))

groups_of <- function(m) {
  h <- hubsrc[[m]]
  g <- lapply(LIN, function(cs) unique(unlist(lapply(cs, function(c) rd(h$cell(c))))))
  c(list(Shared = rd(h$bulk)), g)
}
HG <- lapply(MTH, groups_of); names(HG) <- MTH

for (g in names(HG$GECO))
  stopifnot(setequal(HG$GECO[[g]], rd(file.path(dir0, "GECOhubs_group", paste0(g, ".txt")))))
cat("check: GECO lineage collapse reproduces GECOhubs_group/ exactly\n")

tfs  <- list("B cell TF"          = rd(file.path(dir0, "TFset/B cell TF.txt")),
             "NK cell TF"         = rd(file.path(dir0, "TFset/NK cell TF.txt")),
             "T cell TF"          = rd(file.path(dir0, "TFset/T cell TF.txt")),
             "Low specificity TF" = rd(file.path(dir0, "TFset/Low specificity TF.txt")))
corr <- c(Bcell = "B cell TF", NKcell = "NK cell TF",
          Tcell = "T cell TF", Shared = "Low specificity TF")

covdt <- fread("data/onek1k/covariates.txt")[, .(gene, logE = log(expression), logL = log(gene_length), gc = gc_content)]
covdt[fread("data/onek1k/beta_var.txt"), on = "gene", lbv := log(i.beta_var)]
covdt <- covdt[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(lbv)]

fit_base <- function(hub, tf) {
  ih <- U %in% hub; is <- U %in% tf
  a <- sum(ih & is); b <- sum(ih & !is); c <- sum(!ih & is); d <- sum(!ih & !is); ov <- a
  if (min(a,b,c,d) == 0) { a<-a+.5; b<-b+.5; c<-c+.5; d<-d+.5 }
  l2 <- log2((a*d)/(b*c)); se <- sqrt(1/a+1/b+1/c+1/d)/log(2)
  list(l2 = l2, se = se, p = 2*pnorm(-abs(l2/se)), ov = ov)
}
fit_adj <- function(hub, tf) {
  x <- data.table(gene = U, y = as.integer(U %in% tf), z = as.integer(U %in% hub))
  x <- merge(x, covdt, by = "gene")
  m <- tryCatch(suppressWarnings(glm(y ~ z + logE + logL + gc + lbv, x, family = binomial)),
                error = function(e) NULL)
  if (is.null(m) || !m$converged) return(list(l2 = NA, se = NA, p = NA))
  co <- coef(summary(m)); if (!("z" %in% rownames(co))) return(list(l2 = NA, se = NA, p = NA))
  se <- co["z","Std. Error"]; if (!is.finite(se) || se > 10) return(list(l2 = NA, se = NA, p = NA))
  list(l2 = co["z","Estimate"]/log(2), se = se/log(2), p = co["z","Pr(>|z|)"])
}

grid <- CJ(mth = MTH, hn = names(corr), tn = names(tfs), sorted = FALSE)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  H <- HG[[grid$mth[i]]][[grid$hn[i]]]; G <- tfs[[grid$tn[i]]]
  b <- fit_base(H, G); a <- fit_adj(H, G)
  data.table(method = grid$mth[i], hub = grid$hn[i], tf = grid$tn[i],
             corresponding = (corr[[grid$hn[i]]] == grid$tn[i]),
             nHub = length(intersect(H, U)), nTF = length(G), overlap = b$ov,
             base_l2 = b$l2, base_se = b$se, base_p = b$p,
             adj_l2  = a$l2, adj_se  = a$se, adj_p  = a$p)
}))
res[, `:=`(base_OR = 2^base_l2, adj_OR = 2^adj_l2)]
res[, base_star := fifelse(base_p<.001,"***", fifelse(base_p<.01,"**", fifelse(base_p<.05,"*","")))]
res[, adj_star  := fifelse(is.na(adj_p), "", fifelse(adj_p<.001,"***", fifelse(adj_p<.01,"**", fifelse(adj_p<.05,"*",""))))]
fwrite(res, file.path(outdir, "TF_HPA_enrichment_3methods.csv"))

cat("\nhub group sizes in the GECO universe:\n")
