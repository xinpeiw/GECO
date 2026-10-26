rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(trimws(readLines(f)))
ekey <- function(a,b){ ok <- a!=b & !is.na(a) & !is.na(b); a<-a[ok]; b<-b[ok]; paste(pmin(a,b),pmax(a,b),sep="|") }

DTG <- list(`with metformin`    = list(Immune = rd("data/5e-drug_target/targets_Immune_approved.txt"),
                                       `Non-immune` = rd("data/5e-drug_target/targets_NonImmune_approved.txt")))
DTC <- list(`WITHOUT metformin` = list(Immune = rd("data/7d_drug_target/targets_Immune_approved.txt"),
                                       `Non-immune` = rd("data/7d_drug_target/targets_NonImmune_approved.txt")))
MD <- list(Immune = rd("data/7c-Mendelian_gene/mendelian_immune_full.txt"),
           `Non-immune` = rd("data/7c-Mendelian_gene/mendelian_nonimmune_full.txt"))

fit <- function(P, G, N) {
  a <- as.numeric(sum(P %in% G)); b <- as.numeric(length(P)) - a
  cc <- as.numeric(length(G)) - a; dd <- N - a - b - cc
  if (a == 0) return(list(a=0, expd=NA_real_, l2=NA_real_, se=NA_real_, p=NA_real_))
  if (min(a,b,cc,dd) == 0) { a<-a+.5; b<-b+.5; cc<-cc+.5; dd<-dd+.5 }
  l2 <- log2((a*dd)/(b*cc)); se <- sqrt(1/a+1/b+1/cc+1/dd)/log(2)
  list(a=sum(P %in% G), expd=as.numeric(length(P))*as.numeric(length(G))/N,
       l2=l2, se=se, p=2*pnorm(-abs(l2/se)))
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
fit_adj <- function(P, G, ns, cov, seedv, Sbg=3e6, exact_thresh=5e6) {
  nsc <- intersect(ns, cov$gene); m <- length(nsc); gmap <- setNames(seq_len(m), nsc)
  mi <- match(nsc, cov$gene); cE<-cov$logE[mi]; cL<-cov$logL[mi]; cG<-cov$gc[mi]; cBV<-cov$lbv[mi]
  Pid <- to_ids(P, gmap, m); Gid <- to_ids(G, gmap, m); a <- length(intersect(Pid, Gid))
  if (a == 0) return(list(l2=NA_real_, se=NA_real_, p=NA_real_, nOv=0L))
  Iid <- union(Pid, Gid); nI <- length(Iid); Nc <- choose(m,2); dcount <- Nc - nI
  Ihi <- ((Iid-1)%%m)+1; Ilo <- ((Iid-1)%/%m)+1; inP_I <- as.numeric(Iid%in%Pid); inG_I <- as.numeric(Iid%in%Gid)
  if (Nc <= exact_thresh) { up <- which(upper.tri(matrix(TRUE,m,m))); r <- ((up-1)%%m)+1; cc <- ((up-1)%/%m)+1
    allid <- (r-1)*m+cc; bgid <- allid[!(allid%in%Iid)]; wbg <- rep(1, length(bgid))
  } else { set.seed(seedv); i1 <- sample.int(m,Sbg,TRUE); i2 <- sample.int(m,Sbg,TRUE); ok <- i1!=i2
    lo <- pmin(i1,i2)[ok]; hi <- pmax(i1,i2)[ok]; cand <- unique((lo-1)*m+hi); bgid <- cand[!(cand%in%Iid)]; wbg <- rep(dcount/length(bgid), length(bgid)) }
  bhi <- ((bgid-1)%%m)+1; blo <- ((bgid-1)%/%m)+1
  y <- c(inP_I, rep(0,length(bgid))); inG <- c(inG_I, rep(0,length(bgid)))
  pe <- function(cv) c((cv[Ilo]+cv[Ihi])/2, (cv[blo]+cv[bhi])/2)
  aE<-pe(cE); aL<-pe(cL); aG<-pe(cG); aBV<-pe(cBV)
  w <- c(rep(1,nI), wbg); X <- cbind(1, inG, aE, aL, aG, aBV)
  f <- irls(X,y,w); bg2 <- f$b[2]; se2 <- f$se[2]
  if (!is.finite(bg2)||!is.finite(se2)||abs(bg2)>25||se2>10) return(list(l2=NA_real_, se=NA_real_, p=NA_real_, nOv=a))
  list(l2=bg2/log(2), se=se2/log(2), p=2*pnorm(-abs(bg2/se2)), nOv=a)
}
load_cov <- function(covf, bvf) {
  cv <- fread(covf)[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
  cv[fread(bvf), on="gene", lbv := log(i.beta_var)]
  cv[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]
}

run_cohort <- function(cohort, pmdir, prefix, cts, DTL, covtab) {
  nodes <- rownames(readRDS(file.path(pmdir, paste0(prefix,"alpha.rds")))); N <- choose(length(nodes),2)
  eo <- function(nm){ m <- readRDS(file.path(pmdir, paste0(prefix,nm,".rds"))); rn <- rownames(m)
    idx <- which(m!=0, arr.ind=TRUE); idx <- idx[idx[,1]<idx[,2],,drop=FALSE]
    k <- unique(ekey(rn[idx[,1]], rn[idx[,2]])); rm(m,idx); gc(FALSE); k }
  NETS <- list(Shared = eo("alpha"), Specific = unique(unlist(lapply(cts, eo))))
  cross <- function(A,B){ A<-intersect(A,nodes); B<-intersect(B,nodes)
    k <- unique(ekey(rep(A,each=length(B)), rep(B,times=length(A))))
    ab <- tstrsplit(k,"|",fixed=TRUE); k[(ab[[1]]%in%A&ab[[2]]%in%B)|(ab[[1]]%in%B&ab[[2]]%in%A)] }
  combos <- expand.grid(mv="Full", mt=names(DTL), dn=names(MD), mn=names(MD), nn=names(NETS),
                        KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  combos <- unique(combos[, c("mv","mt","dn","mn","nn")])
  rbindlist(lapply(seq_len(nrow(combos)), function(i) {
    mt <- combos$mt[i]; dn <- combos$dn[i]; mn <- combos$mn[i]; nn <- combos$nn[i]
    P <- cross(DTL[[mt]][[dn]], MD[[mn]])
    r  <- fit(P, NETS[[nn]], N)
    aj <- fit_adj(P, NETS[[nn]], nodes, covtab, seedv = 8000 + i)
    base <- data.table(cohort=cohort, metformin=mt, DT=dn, MD=mn, network=nn,
                       nDT=length(intersect(DTL[[mt]][[dn]], nodes)), nPair=length(P),
                       nEdge=length(NETS[[nn]]), obs=r$a, expected=r$expd)
    rbind(copy(base)[, `:=`(mode="Baseline", l2=r$l2,  se=r$se,  p=r$p)],
          copy(base)[, `:=`(mode="Adjusted", l2=aj$l2, se=aj$se, p=aj$p)])
  }))
}
res <- rbind(
  run_cohort("OneK1K", "data/onek1k_precision", "",
             c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK"), DTG,
             load_cov("data/onek1k/covariates.txt", "data/onek1k/beta_var.txt")),
  run_cohort("CLUES", "data/others/precision_clues", "precision_",
             c("B","NK","T4","T8","cDC","cM","ncM"), DTC,
             load_cov("data/clues/covariates_all.txt", "data/clues/beta_var_all.txt")))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p),"", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[order(cohort, metformin, DT, MD, network, mode)], file.path(outdir,"drugtarget_mendelian_edge_2x2.csv"))
print(res[order(cohort, network, DT, MD, mode), .(cohort, network, DT, MD, mode, nDT, obs, expected=round(expected,1),
                                                  OR=round(OR,2), p=signif(p,3), stars)], nrows=40)
