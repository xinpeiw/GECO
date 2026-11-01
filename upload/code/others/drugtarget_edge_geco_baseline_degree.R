rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(trimws(readLines(f)))
MD <- list(Immune = rd("data/7c-Mendelian_gene/mendelian_immune_full.txt"), `Non-immune` = rd("data/7c-Mendelian_gene/mendelian_nonimmune_full.txt"))
CFG <- list(
  OneK1K = list(pmdir="data/onek1k_precision", prefix="", cts=c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK"),
    dt="data/5e-drug_target", gw="data/5c-GWAS_gene"),
  CLUES  = list(pmdir="data/others/precision_clues", prefix="precision_", cts=c("B","NK","T4","T8","cDC","cM","ncM"),
    dt="data/7d_drug_target", gw="data/7b-GWASgene"))

fit_base <- function(a, nP, nE, N) { a <- as.numeric(a); b <- nP-a; cc <- nE-a; dd <- N-a-b-cc
  if (a == 0) return(list(l2=NA_real_, se=NA_real_, p=NA_real_))
  if (min(a,b,cc,dd) == 0) { a<-a+.5; b<-b+.5; cc<-cc+.5; dd<-dd+.5 }
  l2 <- log2((a*dd)/(b*cc)); se <- sqrt(1/a+1/b+1/cc+1/dd)/log(2); list(l2=l2, se=se, p=2*pnorm(-abs(l2/se))) }

run_cohort <- function(co) {
  cfg <- CFG[[co]]
  rdm <- function(nm) readRDS(file.path(cfg$pmdir, paste0(cfg$prefix, nm, ".rds")))
  genes <- rownames(rdm("alpha")); m <- length(genes); gmap <- setNames(seq_len(m), genes); N <- choose(m, 2)
  pid <- function(lo, hi) (as.numeric(lo) - 1) * m + hi
  eint <- function(nm) { M <- rdm(nm); rn <- rownames(M); stopifnot(setequal(rn, genes))
    idx <- which(M != 0, arr.ind=TRUE); idx <- idx[idx[,1] < idx[,2], , drop=FALSE]
    a <- unname(gmap[rn[idx[,1]]]); b <- unname(gmap[rn[idx[,2]]]); rm(M, idx); gc(FALSE); unique(cbind(pmin(a,b), pmax(a,b))) }
  NETS <- list(Shared = eint("alpha"), Specific = unique(do.call(rbind, lapply(cfg$cts, eint))))
  cat(sprintf("%s: %d nodes | Shared %d edges | Specific (union) %d edges\n", co, m, nrow(NETS$Shared), nrow(NETS$Specific)))
  DT <- list(Immune = rd(file.path(cfg$dt, "targets_Immune_approved.txt")), `Non-immune` = rd(file.path(cfg$dt, "targets_NonImmune_approved.txt")))
  GW <- list(Immune = rd(file.path(cfg$gw, "genelist_Immune.txt")),        `Non-immune` = rd(file.path(cfg$gw, "genelist_NonImmune.txt")))
  DS <- list(Mendelian = MD, GWAS = GW)
  cross_ids <- function(A, B) { A <- unname(gmap[intersect(A, genes)]); B <- unname(gmap[intersect(B, genes)])
    a <- rep(A, each=length(B)); b <- rep(B, times=length(A)); ok <- a != b; unique(pid(pmin(a,b)[ok], pmax(a,b)[ok])) }
  fit_degree <- function(Pid, E, nb=50) {
    d <- tabulate(c(E[,1], E[,2]), nbins=m); pos <- d > 0; br <- unique(quantile(d[pos], probs=seq(0, 1, length.out=nb+1)))
    bin <- integer(m); bin[pos] <- cut(d[pos], breaks=br, include.lowest=TRUE, labels=FALSE); K <- max(bin) + 1L; nk <- tabulate(bin + 1L, K)
    key <- function(b1, b2) pmin(b1, b2) * K + pmax(b1, b2)
    tot <- CJ(lo=0:(K-1), hi=0:(K-1))[lo <= hi]
    tot[, n := fifelse(lo == hi, nk[lo+1]*(nk[lo+1]-1)/2, as.numeric(nk[lo+1])*nk[hi+1])]; tot[, k := lo*K + hi]
    Eid <- pid(E[,1], E[,2]); plo <- ((Pid-1) %/% m) + 1; phi <- ((Pid-1) %% m) + 1; pe <- Pid %in% Eid
    if (!any(pe)) return(list(l2=NA_real_, se=NA_real_, p=NA_real_))
    tot[data.table(k=key(bin[E[,1]], bin[E[,2]]))[, .(nE=.N), by=k], on="k", nE := i.nE]
    tot[data.table(k=key(bin[plo], bin[phi]), e=pe)[, .(nP=.N, nPE=sum(e)), by=k], on="k", `:=`(nP=i.nP, nPE=i.nPE)]
    for (cn in c("nE","nP","nPE")) tot[is.na(get(cn)), (cn) := 0]
    stopifnot(sum(tot$n) == N, sum(tot$nE) == nrow(E), sum(tot$nP) == length(Pid))
    cells <- rbind(tot[, .(lo, hi, inE=1, succ=nPE,      trials=nE)],
                   tot[, .(lo, hi, inE=0, succ=nP - nPE, trials=n - nE)])[trials > 0]
    Xb <- sapply(1:(K-1), function(k) ((cells$lo == k) + (cells$hi == k)) / 2)
    f  <- suppressWarnings(glm(cbind(cells$succ, cells$trials - cells$succ) ~ cells$inE + Xb, family=binomial()))
    co <- summary(f)$coefficients["cells$inE", ]
    b <- unname(co[1]); s <- unname(co[2])
    if (!is.finite(b) || !is.finite(s) || abs(b) > 25 || s > 10) return(list(l2=NA_real_, se=NA_real_, p=NA_real_))
    list(l2=b/log(2), se=s/log(2), p=2*pnorm(-abs(b/s))) }
  rbindlist(lapply(names(DS), function(src) rbindlist(lapply(names(DT), function(cl) {
    Pid <- cross_ids(DT[[cl]], DS[[src]][[cl]])
    rbindlist(lapply(names(NETS), function(nw) { E <- NETS[[nw]]; a <- length(intersect(Pid, pid(E[,1], E[,2])))
      b0 <- fit_base(a, length(Pid), nrow(E), N); b1 <- fit_degree(Pid, E)
      cat(sprintf("  %-9s %-10s %-8s obs %4d | baseline %5.2f | degree-adjusted %5.2f (p %.2g)\n", src, cl, nw, a, 2^b0$l2, 2^b1$l2, b1$p))
      base <- data.table(cohort=co, source=src, class=cl, network=nw, nPair=length(Pid), nEdge=nrow(E), obs=a, expected=as.numeric(length(Pid))*nrow(E)/N)
      rbind(copy(base)[, `:=`(mode="Baseline",        l2=b0$l2, se=b0$se, p=b0$p)],
            copy(base)[, `:=`(mode="Degree-adjusted", l2=b1$l2, se=b1$se, p=b1$p)]) })) }))))
}
res <- rbind(run_cohort("OneK1K"), run_cohort("CLUES"))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p), "", fifelse(p<.001, "***", fifelse(p<.01, "**", fifelse(p<.05, "*", ""))))]

fwrite(res[order(cohort, source, class, network, mode)], file.path(outdir, "drugtarget_edge_geco_baseline_degree.csv"))
print(dcast(res[, .(cohort, source, class, network, mode, v=sprintf("%.2f%s", OR, stars))], cohort + source + class + network ~ mode, value.var="v"))
