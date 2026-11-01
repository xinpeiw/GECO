rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })
dir0   <- "data/4a-TF_HPA"
outdir <- "results/4a-TF_HPA"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(grep("^ENSG", trimws(readLines(f)), value = TRUE))
U  <- rd("data/onek1k/GECO_universe.txt")
hubs <- list(Bcell = rd(file.path(dir0,"GECOhubs_group/Bcell.txt")),  NKcell = rd(file.path(dir0,"GECOhubs_group/NKcell.txt")),
             Tcell = rd(file.path(dir0,"GECOhubs_group/Tcell.txt")),  Shared = rd(file.path(dir0,"GECOhubs_group/Shared.txt")))
tfs  <- list("B cell TF" = rd(file.path(dir0,"TFset/B cell TF.txt")), "NK cell TF" = rd(file.path(dir0,"TFset/NK cell TF.txt")),
             "T cell TF" = rd(file.path(dir0,"TFset/T cell TF.txt")), "Low specificity TF" = rd(file.path(dir0,"TFset/Low specificity TF.txt")))
corr <- c(Bcell="B cell TF", NKcell="NK cell TF", Tcell="T cell TF", Shared="Low specificity TF")
covdt <- fread("data/onek1k/covariates.txt")[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
covdt[fread("data/onek1k/beta_var.txt"), on="gene", lbv:=log(i.beta_var)]; covdt <- covdt[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lbv)]

fit_base <- function(hub, tf) {
  ih <- U %in% hub; is <- U %in% tf
  a<-sum(ih&is); b<-sum(ih&!is); c<-sum(!ih&is); d<-sum(!ih&!is); ov<-a
  if (min(a,b,c,d)==0) { a<-a+.5;b<-b+.5;c<-c+.5;d<-d+.5 }
  l2<-log2((a*d)/(b*c)); se<-sqrt(1/a+1/b+1/c+1/d)/log(2); list(l2=l2, se=se, p=2*pnorm(-abs(l2/se)), ov=ov)
}
fit_adj <- function(hub, tf) {
  x<-data.table(gene=U, y=as.integer(U%in%tf), z=as.integer(U%in%hub)); x<-merge(x,covdt,by="gene")
  m<-tryCatch(suppressWarnings(glm(y~z+logE+logL+gc+lbv, x, family=binomial)), error=function(e) NULL)
  if(is.null(m)||!m$converged) return(list(l2=NA,se=NA,p=NA)); co<-coef(summary(m)); if(!("z"%in%rownames(co))) return(list(l2=NA,se=NA,p=NA))
  se<-co["z","Std. Error"]; if(!is.finite(se)||se>10) return(list(l2=NA,se=NA,p=NA)); list(l2=co["z","Estimate"]/log(2), se=se/log(2), p=co["z","Pr(>|z|)"])
}

res <- rbindlist(lapply(names(hubs), function(hn) rbindlist(lapply(names(tfs), function(tn) {
  b <- fit_base(hubs[[hn]], tfs[[tn]]); diag <- (corr[[hn]] == tn)
  a <- fit_adj(hubs[[hn]], tfs[[tn]])
  data.table(hub=hn, tf=tn, corresponding=diag, nTF=length(tfs[[tn]]), overlap=b$ov,
             base_l2=b$l2, base_se=b$se, base_p=b$p, adj_l2=a$l2, adj_se=a$se, adj_p=a$p)
}))))
res[, base_OR := 2^base_l2]; res[, adj_OR := 2^adj_l2]
res[, base_star := fifelse(base_p<.001,"***",fifelse(base_p<.01,"**",fifelse(base_p<.05,"*","")))]
res[, adj_star  := fifelse(is.na(adj_p),"",fifelse(adj_p<.001,"***",fifelse(adj_p<.01,"**",fifelse(adj_p<.05,"*",""))))]
fwrite(res, file.path(outdir, "TF_HPA_enrichment.csv"))
cat("\n===== GECO hub group x HPA TF set enrichment (baseline log2 OR; * on diagonal = corresponding) =====\n")
print(res[, .(hub, tf, corresponding, nTF, overlap, base_log2OR=round(base_l2,2), base_OR=round(base_OR,2), base_p=signif(base_p,2), base_star)])
cat("\n----- corresponding (diagonal): baseline vs 4-cov adjusted -----\n")
print(res[corresponding==TRUE, .(hub, tf, overlap, base_OR=round(base_OR,2), base_p=signif(base_p,2), base_star, adj_OR=round(adj_OR,2), adj_p=signif(adj_p,2), adj_star)])
