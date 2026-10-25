rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

COND   <- "mega"
COVTAG <- "all"
DLAB   <- "all"

evo    <- "data/others/evolution_clues"
covdir <- "data/clues"
outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) unique(trimws(readLines(f)))

U   <- rd("data/clues/allgenes_CLUES.txt")
cov <- fread(file.path(covdir, paste0("covariates_", COVTAG, ".txt")))[, .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
bv  <- fread(file.path(covdir, paste0("beta_var_",   COVTAG, ".txt")))[, .(gene, logBV=log(beta_var))]
cov <- merge(cov, bv, by="gene")
rc  <- fread(file.path(evo, "CLUES_gene_recombination_rate.1000G.hg38.tsv"))[, .(gene=gene_id, recomb=recomb_rate_cM_Mb)]
ex  <- fread(file.path(evo, "exon.txt"))[, .(gene=ensg, nexon=n_exons_unique)]
cov[rc, on="gene", lrec := log1p(i.recomb)]
cov[ex, on="gene", lexn := log(i.nexon)]

ctmap  <- c(B="B", CD4="T4", CD8="T8", NK="NK", cDC="cDC", cM="cM")
enh_ct <- names(ctmap)
anns   <- list(list(dir="zooHAR_positive", lab="zooHAR"), list(dir="HAQER", lab="HAQER"), list(dir="akbari", lab="Akbari"))

fit_base <- function(inhub, inset) {
  a<-sum(inhub&inset); b<-sum(inhub&!inset); c<-sum(!inhub&inset); d<-sum(!inhub&!inset)
  if (a==0||b==0||c==0||d==0){a<-a+.5;b<-b+.5;c<-c+.5;d<-d+.5}
  list(b=log((a*d)/(b*c)), se=sqrt(1/a+1/b+1/c+1/d), ov=sum(inhub&inset))
}
fit_adj <- function(hubG, setG) {
  d <- merge(data.table(gene=U, in_hub=as.integer(U%in%hubG), in_set=as.integer(U%in%setG)), cov, by="gene")
  d <- d[is.finite(logE)&is.finite(logL)&is.finite(gc)&is.finite(lexn)&is.finite(lrec)&is.finite(logBV)]
  a2<-sum(d$in_hub&d$in_set); b2<-sum(d$in_hub&!d$in_set)
  if (a2==0||b2==0) return(list(b=NA,se=NA,ov=a2))
  m <- tryCatch(suppressWarnings(glm(in_set~in_hub+logE+logL+gc+lexn+lrec+logBV, d, family=binomial)), error=function(e) NULL)
  if (is.null(m)||!m$converged) return(list(b=NA,se=NA,ov=a2))
  co <- coef(summary(m)); if(!("in_hub"%in%rownames(co))) return(list(b=NA,se=NA,ov=a2))
  se <- co["in_hub","Std. Error"]; if(!is.finite(se)||se>10) return(list(b=NA,se=NA,ov=a2))
  list(b=co["in_hub","Estimate"], se=se, ov=a2)
}
meta_fe <- function(bs, ses){ ok<-is.finite(bs)&is.finite(ses)&ses>0; bs<-bs[ok]; ses<-ses[ok]
  if(!length(bs)) return(list(b=NA,se=NA,p=NA,k=0)); w<-1/ses^2; bm<-sum(w*bs)/sum(w); sem<-sqrt(1/sum(w))
  list(b=bm, se=sem, p=2*pnorm(-abs(bm/sem)), k=length(bs)) }
zp <- function(b,se) if(is.na(b)||is.na(se)) NA else 2*pnorm(-abs(b/se))

ALPHA <- intersect(rd("data/clues/hubs_alpha.txt"), U)
HUBc  <- lapply(enh_ct, function(ec) intersect(rd(file.path("data/clues", paste0("hubs_", ctmap[ec], ".txt"))), U)); names(HUBc) <- enh_ct

res <- list(); pc <- list()
for (a in anns) {
  Gc <- lapply(enh_ct, function(ec) intersect(rd(file.path(evo, a$dir, paste0(ec,"_enhanceronly.txt"))), U)); names(Gc) <- enh_ct
  Gunion <- Reduce(union, Gc)
  for (adj in c(FALSE, TRUE)) {
    bb <- sapply(enh_ct, function(ec) {
      r <- if (adj) fit_adj(HUBc[[ec]], Gc[[ec]]) else fit_base(U%in%HUBc[[ec]], U%in%Gc[[ec]]); c(r$b, r$se, r$ov)
    })
    if (!adj) for (i in seq_along(enh_ct)) pc[[length(pc)+1]] <- data.table(dataset=DLAB, annotation=a$lab,
      enh_ct=enh_ct[i], hub_ct=ctmap[enh_ct[i]], nSet=length(Gc[[enh_ct[i]]]), overlap=bb[3,i], l2=bb[1,i]/log(2))
    mb <- meta_fe(bb[1,], bb[2,])
    sb <- if (adj) fit_adj(ALPHA, Gunion) else fit_base(U%in%ALPHA, U%in%Gunion)
    md <- ifelse(adj, "Adjusted", "Baseline")
    res[[length(res)+1]] <- rbindlist(list(
      data.table(dataset=DLAB, annotation=a$lab, mode=md, group="CLUES Specific", l2=mb$b/log(2), se=mb$se/log(2), p=mb$p,          k=mb$k, nSet=length(Gunion), overlap=NA_integer_),
      data.table(dataset=DLAB, annotation=a$lab, mode=md, group="CLUES Shared",   l2=sb$b/log(2), se=sb$se/log(2), p=zp(sb$b,sb$se), k=1L,    nSet=length(Gunion), overlap=sb$ov)))
  }
}
res <- rbindlist(res); pc <- rbindlist(pc)
res[, OR := 2^l2]; res[, stars := fifelse(is.na(p),"",fifelse(p<.001,"***",fifelse(p<.01,"**",fifelse(p<.05,"*",""))))]
fwrite(res[order(annotation,mode,group)], file.path(outdir, "CLUES_evolution_rE2G_enrichment.csv"))
fwrite(pc,  file.path(outdir, "CLUES_evolution_rE2G_enrichment_percelltype.csv"))
cat(sprintf("\n===== CLUES rE2G enrichment (%s / %s) : baseline + adjusted log2 OR =====\n", DLAB, COND))
print(res[order(annotation,mode,group), .(annotation, mode, group, k, nSet, overlap, l2=round(l2,2), OR=round(OR,2), p=signif(p,2), stars)])
