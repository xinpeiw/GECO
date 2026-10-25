rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

psf    <- "data/6b-edge/perturb_seq/pertseq_normed_effect_network_matrix.rds"
pmdir  <- "data/onek1k_precision"
outdir <- "results/6b-edge"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

ekey <- function(a, b) { ok <- a!=b & !is.na(a) & !is.na(b); a<-a[ok]; b<-b[ok]; paste(pmin(a,b), pmax(a,b), sep="|") }
edges_of <- function(name) {
  m <- readRDS(file.path(pmdir, paste0(name, ".rds"))); rn <- rownames(m)
  idx <- which(m != 0, arr.ind=TRUE); idx <- idx[idx[,1] < idx[,2], , drop=FALSE]
  k <- unique(ekey(rn[idx[,1]], rn[idx[,2]])); rm(m, idx); gc(FALSE); k
}

M  <- readRDS(psf); rn <- rownames(M)
idx <- which(M != 0, arr.ind=TRUE); idx <- idx[idx[,1] < idx[,2], , drop=FALSE]
P_all <- unique(ekey(rn[idx[,1]], rn[idx[,2]])); rm(M, idx); gc(FALSE)
gnodes <- rownames(readRDS(file.path(pmdir, "alpha.rds")))
ns <- intersect(rn, gnodes); N <- choose(length(ns), 2)
pab <- tstrsplit(P_all, "|", fixed=TRUE); P <- P_all[pab[[1]] %in% ns & pab[[2]] %in% ns]
cat(sprintf("perturb-seq edges: %s | node set: %d | choose(n,2)=%s\n",
            format(length(P), big.mark=","), length(ns), format(N, big.mark=",")))

Gnets <- list(
  Shared = edges_of("alpha"),
  CD4    = union(edges_of("CD4NC"), edges_of("CD4ET")),
  CD8    = union(edges_of("CD8NC"), edges_of("CD8ET")),
  NK     = edges_of("NK"),
  B      = union(edges_of("BIN"), edges_of("BMem")))
restr <- function(keys, ns) { ab <- tstrsplit(keys, "|", fixed=TRUE); keys[ab[[1]] %in% ns & ab[[2]] %in% ns] }

cov <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
cov[fread("data/onek1k/beta_var.txt"), on="gene", lbv:=log(i.beta_var)]
cov <- cov[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]; COVG <- cov$gene

fit_edge <- function(P, G, N) {
  a <- sum(P %in% G); b <- length(P) - a; c <- length(G) - a; d <- N - a - b - c
  if (min(a,b,c,d) == 0) { A<-a+.5;B<-b+.5;C<-c+.5;D<-d+.5; lor<-log((A*D)/(B*C)); se<-sqrt(1/A+1/B+1/C+1/D)
    return(list(nov=a, l2=lor/log(2), se=se/log(2), p=2*pnorm(-abs(lor/se)))) }
  co <- coef(summary(suppressWarnings(glm(in_P ~ in_G, binomial, weights=w,
        data=data.frame(in_P=c(1L,1L,0L,0L), in_G=c(1L,0L,1L,0L), w=c(a,b,c,d))))))["in_G",]
  list(nov=a, l2=co["Estimate"]/log(2), se=co["Std. Error"]/log(2), p=co["Pr(>|z|)"])
}
irls <- function(X, y, w, maxit=60, tol=1e-9) {
  b <- rep(0, ncol(X)); XtWX <- NULL
  for (it in 1:maxit) { eta <- as.vector(X %*% b); mu <- pmin(pmax(1/(1+exp(-eta)),1e-12),1-1e-12)
    Wd <- w*mu*(1-mu); z <- eta+(y-mu)/(mu*(1-mu)); XtWX <- crossprod(X, Wd*X)
    bn <- tryCatch(solve(XtWX, crossprod(X, Wd*z)), error=function(e) NULL)
    if (is.null(bn)) return(list(b=rep(NA_real_,ncol(X)), se=rep(NA_real_,ncol(X))))
    if (max(abs(bn-b)) < tol) { b <- as.vector(bn); break }; b <- as.vector(bn) }
  V <- tryCatch(solve(XtWX), error=function(e) matrix(NA_real_,ncol(X),ncol(X))); list(b=as.vector(b), se=sqrt(diag(V)))
}
to_ids <- function(keys, gmap, m) { if (!length(keys)) return(numeric(0))
  ab <- tstrsplit(keys, "|", fixed=TRUE); ia <- gmap[ab[[1]]]; ib <- gmap[ab[[2]]]
  ok <- !is.na(ia)&!is.na(ib); ia<-ia[ok]; ib<-ib[ok]; lo<-pmin(ia,ib); hi<-pmax(ia,ib); unique((lo-1)*m+hi) }
fit_edge_adj <- function(P, G, ns, seedv, Sbg=3e6, exact_thresh=5e6) {
  nsc <- intersect(ns, COVG); m <- length(nsc); gmap <- setNames(seq_len(m), nsc)
  mi <- match(nsc, cov$gene); cE<-cov$logE[mi]; cL<-cov$logL[mi]; cG<-cov$gc[mi]; cBV<-cov$lbv[mi]
  Pid <- to_ids(P, gmap, m); Gid <- to_ids(G, gmap, m); a <- length(intersect(Pid, Gid))
  if (a == 0) return(list(l2=NA_real_, se=NA_real_, p=NA_real_, nocov=NA_real_, nOv=0L))
  Iid <- union(Pid, Gid); nI <- length(Iid); Nc <- choose(m,2); dcount <- Nc - nI
  Ihi <- ((Iid-1)%%m)+1; Ilo <- ((Iid-1)%/%m)+1; inP_I <- as.numeric(Iid%in%Pid); inG_I <- as.numeric(Iid%in%Gid)
  if (Nc <= exact_thresh) { up <- which(upper.tri(matrix(TRUE,m,m))); r <- ((up-1)%%m)+1; cc <- ((up-1)%/%m)+1
    allid <- (r-1)*m+cc; bgid <- allid[!(allid%in%Iid)]; wbg <- rep(1, length(bgid))
  } else { set.seed(seedv); i1 <- sample.int(m,Sbg,TRUE); i2 <- sample.int(m,Sbg,TRUE); ok <- i1!=i2
    lo <- pmin(i1,i2)[ok]; hi <- pmax(i1,i2)[ok]; cand <- unique((lo-1)*m+hi); bgid <- cand[!(cand%in%Iid)]; wbg <- rep(dcount/length(bgid), length(bgid)) }
  bhi <- ((bgid-1)%%m)+1; blo <- ((bgid-1)%/%m)+1
  y <- c(inP_I, rep(0,length(bgid))); inG <- c(inG_I, rep(0,length(bgid)))
  aE <- c((cE[Ilo]+cE[Ihi])/2,(cE[blo]+cE[bhi])/2); aL <- c((cL[Ilo]+cL[Ihi])/2,(cL[blo]+cL[bhi])/2)
  aG <- c((cG[Ilo]+cG[Ihi])/2,(cG[blo]+cG[bhi])/2); aBV <- c((cBV[Ilo]+cBV[Ihi])/2,(cBV[blo]+cBV[bhi])/2)
  w <- c(rep(1,nI), wbg); X <- cbind(1, inG, aE, aL, aG, aBV)
  fit <- irls(X,y,w); nc <- irls(cbind(1,inG),y,w)
  bg2 <- fit$b[2]; se2 <- fit$se[2]
  if (!is.finite(bg2)||!is.finite(se2)||abs(bg2)>25) return(list(l2=NA_real_, se=NA_real_, p=NA_real_, nocov=nc$b[2]/log(2), nOv=a))
  list(l2=bg2/log(2), se=se2/log(2), p=2*pnorm(-abs(bg2/se2)), nocov=nc$b[2]/log(2), nOv=a)
}

res <- rbindlist(lapply(seq_along(Gnets), function(gi) {
  gn <- names(Gnets)[gi]; G <- restr(Gnets[[gn]], ns)
  b <- fit_edge(P, G, N); aj <- fit_edge_adj(P, G, ns, seedv=2000+gi)
  cat(sprintf("  %-7s nGedge=%s baseline log2OR=%6.2f | adj=%6s (sanity no-cov=%6s)\n", gn,
              format(length(G), big.mark=","), b$l2, ifelse(is.na(aj$l2),"NA",sprintf("%.2f",aj$l2)),
              ifelse(is.na(aj$nocov),"NA",sprintf("%.2f",aj$nocov))))
  rbind(data.table(network=gn, mode="Baseline", nGedge=length(G), nPedge=length(P), nOv=b$nov, l2=b$l2, se=b$se, p=b$p),
        data.table(network=gn, mode="Adjusted", nGedge=length(G), nPedge=length(P), nOv=aj$nOv, l2=aj$l2, se=aj$se, p=aj$p))
}))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p),"", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[, .(network, mode, nGedge, nPedge, nOv, log2OR=l2, se, OR, p, stars)], file.path(outdir, "perturbation_edge_enrichment.csv"))
cat("\n===== perturb-seq edge enrichment of GECO networks (log2 OR) =====\n")
print(res[, .(network, mode, nGedge, nOv, l2=round(l2,2), OR=round(OR,2), p=signif(p,2), stars)])
