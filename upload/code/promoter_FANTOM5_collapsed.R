rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

B <- 5000

set.seed(1)

hubdir <- "data/onek1k"
outdir <- "results/4b-promoter_enhancer_TFBS"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

rd <- function(f) unique(trimws(readLines(f)))
U       <- rd(file.path(hubdir, "GECO_universe.txt")); nU <- length(U)
SHARED  <- intersect(rd(file.path(hubdir, "GECO_Shared.txt")), U)
ct      <- c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK")
HUBc    <- lapply(ct, function(c) intersect(rd(file.path(hubdir, paste0("hubs_", c, ".txt"))), U)); names(HUBc) <- ct
ALLSPEC <- Reduce(union, HUBc)

cv  <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
cv  <- merge(cv, fread("data/onek1k/beta_var.txt")[, .(gene, lbv=log(beta_var))], by="gene")
COV <- as.matrix(cv[match(U, gene), .(logE, logL, gc, lbv)])
cov_ok <- is.finite(rowSums(COV))

getv <- function(f) {
  d <- fread(file.path("data/4b-promoter_enhancer_TFBS/promoter_FANTOM5", f)); setnames(d, 1:2, c("g","v"))
  log1p(d$v[match(U, d$g)])
}

assays <- list(promoter = list(
  lab="Promoter", xlab="Promoter enrichment", stub="promoter_FANTOM5",
  units=lapply(list(list("CD4", c("CD4NC","CD4ET"), "CD4.txt"),
                    list("CD8", c("CD8NC","CD8ET"), "CD8.txt"),
                    list("NK",  "NK",               "NK.txt"),
                    list("B",   c("BIN","BMem"),    "B.txt")),
               function(z) list(lab=z[[1]], cells=z[[2]], v=getv(z[[3]])))))

GRP  <- c("Focal specific","Non-focal specific","Shared")
MODE <- c("Baseline","Adjusted")

for (an in names(assays)) for (u in seq_along(assays[[an]]$units)) {
  un <- assays[[an]]$units[[u]]; foc <- Reduce(union, HUBc[un$cells])
  assays[[an]]$units[[u]]$G <- cbind(as.integer(U %in% foc),
                                     as.integer(U %in% setdiff(ALLSPEC, foc)),
                                     as.integer(U %in% SHARED))
}

mean_over_units <- function(A, idx) {
  acc <- array(NA_real_, c(3, 2, length(A$units)))
  for (u in seq_along(A$units)) {
    un <- A$units[[u]]; v <- un$v
    for (m in 1:2) {
      keep <- if (m == 1L) is.finite(v) else is.finite(v) & cov_ok
      sel  <- idx[keep[idx]]; y <- v[sel]
      for (g in 1:3) {
        x <- un$G[sel, g]
        if (sum(x) < 2L || sum(x) > length(x) - 2L) next
        X  <- if (m == 1L) cbind(1, x) else cbind(1, x, COV[sel, , drop=FALSE])
        cf <- tryCatch(.lm.fit(X, y)$coefficients, error=function(e) NULL)
        if (!is.null(cf) && is.finite(cf[2])) acc[g, m, u] <- cf[2]
      }
    }
  }
  m2 <- matrix(NA_real_, 3, 2)
  for (g in 1:3) for (m in 1:2) m2[g, m] <- mean(acc[g, m, ], na.rm=TRUE)
  list(mean = m2, per = acc)
}

allres <- list()
for (an in names(assays)) {
  A <- assays[[an]]; k <- length(A$units)
  o  <- mean_over_units(A, seq_len(nU))

  per <- rbindlist(lapply(seq_along(A$units), function(u) rbindlist(lapply(seq_len(6), function(i) {
    g <- ((i-1) %% 3) + 1; m <- ((i-1) %/% 3) + 1
    data.table(assay=A$lab, unit=A$units[[u]]$lab, group=GRP[g], mode=MODE[m], beta=o$per[g, m, u])
  }))))
  fwrite(per, file.path(outdir, sprintf("%s_collapsed_perunit.csv", A$stub)))

  bt <- array(NA_real_, c(3, 2, B)); t0 <- proc.time()[3]
  for (b in seq_len(B)) bt[, , b] <- mean_over_units(A, sample.int(nU, nU, replace=TRUE))$mean
  cat(sprintf("%-16s k=%d  B=%d in %.0f s\n", A$lab, k, B, proc.time()[3]-t0))

  R <- rbindlist(lapply(seq_len(6), function(i) {
    g <- ((i-1) %% 3) + 1; m <- ((i-1) %/% 3) + 1
    d <- bt[g, m, ]; d <- d[is.finite(d)]; ci <- quantile(d, c(.025,.975))
    data.table(assay=A$lab, group=GRP[g], mode=MODE[m], k=k, beta=o$mean[g, m],
               se_boot=sd(d), lo=ci[[1]], hi=ci[[2]],
               p_perc=min(1, 2*min((1+sum(d<=0))/(length(d)+1), (1+sum(d>=0))/(length(d)+1))),
               p_wald=2*pnorm(-abs(o$mean[g, m]/sd(d))), nrep=length(d))
  }))
  R[, `:=`(lo_se = beta - se_boot, hi_se = beta + se_boot)]
  R[, stars := fifelse(p_wald<.001,"***", fifelse(p_wald<.01,"**", fifelse(p_wald<.05,"*","")))]
  allres[[an]] <- R
  cat(sprintf("  %-19s %-8s beta %9.3f  +/-1SE [%9.3f, %9.3f]  p %9.3g %s\n",
              R$group, R$mode, R$beta, R$lo_se, R$hi_se, R$p_wald, R$stars), sep="")
  fwrite(R, file.path(outdir, sprintf("%s_collapsed.csv", A$stub)))
}
