rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(trimws(readLines(f)))
ekey <- function(a,b){ ok <- a!=b & !is.na(a) & !is.na(b); a<-a[ok]; b<-b[ok]; paste(pmin(a,b),pmax(a,b),sep="|") }

DTG <- list(
  `WITHOUT metformin` = list(Immune = rd("data/5e-drug_target/targets_Immune_approved.txt"),
                             `Non-immune` = rd("data/others/targets_onek1k_noMetformin_NonImmune_approved.txt")),
  `with metformin`    = list(Immune = rd("data/5e-drug_target/targets_Immune_approved.txt"),
                             `Non-immune` = rd("data/5e-drug_target/targets_NonImmune_approved.txt")))
DTC <- list(
  `WITHOUT metformin` = list(Immune = rd("data/7d_drug_target/targets_Immune_approved.txt"),
                             `Non-immune` = rd("data/7d_drug_target/targets_NonImmune_approved.txt")),
  `with metformin`    = list(Immune = rd("data/7d_drug_target/targets_Immune_approved.txt"),
                             `Non-immune` = rd("data/others/targets_clues_withMetformin_NonImmune_approved.txt")))
MD <- list(Immune = rd("data/others/mendelian_immune_exclusive.txt"),
           `Non-immune` = rd("data/others/mendelian_nonimmune_exclusive.txt"))

fit <- function(P, G, N) {
  a <- as.numeric(sum(P %in% G)); b <- as.numeric(length(P)) - a
  cc <- as.numeric(length(G)) - a; dd <- N - a - b - cc
  if (a == 0) return(list(a=0, expd=NA_real_, l2=NA_real_, se=NA_real_, p=NA_real_))
  if (min(a,b,cc,dd) == 0) { a<-a+.5; b<-b+.5; cc<-cc+.5; dd<-dd+.5 }
  l2 <- log2((a*dd)/(b*cc)); se <- sqrt(1/a+1/b+1/cc+1/dd)/log(2)
  list(a=sum(P %in% G), expd=as.numeric(length(P))*as.numeric(length(G))/N,
       l2=l2, se=se, p=2*pnorm(-abs(l2/se)))
}
run_cohort <- function(cohort, pmdir, prefix, cts, DTL) {
  nodes <- rownames(readRDS(file.path(pmdir, paste0(prefix,"alpha.rds")))); N <- choose(length(nodes),2)
  eo <- function(nm){ m <- readRDS(file.path(pmdir, paste0(prefix,nm,".rds"))); rn <- rownames(m)
    idx <- which(m!=0, arr.ind=TRUE); idx <- idx[idx[,1]<idx[,2],,drop=FALSE]
    k <- unique(ekey(rn[idx[,1]], rn[idx[,2]])); rm(m,idx); gc(FALSE); k }
  NETS <- list(Shared = eo("alpha"), Specific = unique(unlist(lapply(cts, eo))))
  cross <- function(A,B){ A<-intersect(A,nodes); B<-intersect(B,nodes)
    k <- unique(ekey(rep(A,each=length(B)), rep(B,times=length(A))))
    ab <- tstrsplit(k,"|",fixed=TRUE); k[(ab[[1]]%in%A&ab[[2]]%in%B)|(ab[[1]]%in%B&ab[[2]]%in%A)] }
  combos <- expand.grid(mt=names(DTL), dn=names(MD), mn=names(MD), nn=names(NETS),
                        KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  combos <- unique(combos[, c("mt","dn","mn","nn")])
  rbindlist(lapply(seq_len(nrow(combos)), function(i) {
    mt <- combos$mt[i]; dn <- combos$dn[i]; mn <- combos$mn[i]; nn <- combos$nn[i]
    P <- cross(DTL[[mt]][[dn]], MD[[mn]])
    r <- fit(P, NETS[[nn]], N)
    data.table(cohort=cohort, metformin=mt, DT=dn, MD=mn, network=nn,
               nDT=length(intersect(DTL[[mt]][[dn]], nodes)), nPair=length(P),
               nEdge=length(NETS[[nn]]), obs=r$a, expected=r$expd,
               l2=r$l2, se=r$se, p=r$p)
  }))
}
res <- rbind(
  run_cohort("OneK1K", "data/onek1k_precision", "",
             c("BIN","BMem","CD4ET","CD4NC","CD8ET","CD8NC","NK"), DTG),
  run_cohort("CLUES", "data/others/precision_clues", "precision_",
             c("B","NK","T4","T8","cDC","cM","ncM"), DTC))
res[, OR := 2^l2]
res[, stars := fifelse(is.na(p),"", fifelse(p<.001,"***", fifelse(p<.01,"**", fifelse(p<.05,"*",""))))]
fwrite(res[order(cohort, metformin, DT, MD, network)], file.path(outdir,"drugtarget_mendelian_edge_matched.csv"))

cat("\nPRIMARY figure = locked convention (OneK1K all drugs / CLUES no metformin)\n")
print(res[(cohort=="OneK1K" & metformin=="with metformin") | (cohort=="CLUES" & metformin=="WITHOUT metformin")][
      network=="Specific"][order(cohort, DT, MD), .(cohort, DT, MD, nDT, obs, OR=round(OR,2), p=signif(p,3), stars)])
