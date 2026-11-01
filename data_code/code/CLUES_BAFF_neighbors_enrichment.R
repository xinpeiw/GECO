rm(list=ls())
suppressPackageStartupMessages({ library(data.table) })

src    <- "data/7f_BAFF"
outdir <- "results/7f_BAFF"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

rd <- function(f) {
  g <- unique(trimws(readLines(f, warn=FALSE)))
  g <- g[nzchar(g)]
  if (!length(g)) stop(sprintf("%s yielded 0 genes -- unreadable or empty", f))
  g
}

U  <- rd("data/clues/allgenes_CLUES.txt")
NB <- rd(file.path(src,  "BAFF_neighbors.txt"))
stopifnot(length(U) == 11087L, length(NB) == 55L, all(NB %in% U))

SETS <- list(
  list(f="geneset_HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv",
       lab="Interferon alpha response",   key="HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  list(f="geneset_GOBP_MYELOID_LEUKOCYTE_ACTIVATION.tsv",
       lab="Myeloid leukocyte activation", key="GOBP_MYELOID_LEUKOCYTE_ACTIVATION"),
  list(f="geneset_IMMUNOPROTEASOME_COMPLEX_curated.tsv",
       lab="Immunoproteasome complex",     key="IMMUNOPROTEASOME_COMPLEX_curated"))

gs <- lapply(SETS, function(s) {
  d <- fread(file.path(src, s$f))
  stopifnot(all(c("ensembl_gene","in_CLUES_universe","in_BAFF_neighbors") %in% names(d)))
  g <- unique(d[in_CLUES_universe == TRUE, ensembl_gene])

  stopifnot(setequal(intersect(g, NB), unique(d[in_CLUES_universe & in_BAFF_neighbors, ensembl_gene])))
  stopifnot(all(g %in% U))
  g
})
names(gs) <- vapply(SETS, function(s) s$key, character(1))

up <- fread(file.path(src, "enrichment_results.tsv"))
stopifnot(nrow(up) == 3L, all(names(gs) %in% up$gene_set))
chk <- data.table(gene_set = names(gs),
                  n_set_here = vapply(gs, length, integer(1)),
                  n_ov_here  = vapply(gs, function(g) length(intersect(g, NB)), integer(1)))
chk <- merge(chk, up[, .(gene_set, set_size_in_universe, overlap, odds_ratio, p_value)], by="gene_set")
stopifnot(chk[, all(n_set_here == set_size_in_universe)], chk[, all(n_ov_here == overlap)])
cat("check: set sizes and overlaps reproduce enrichment_results.tsv exactly\n")

xl <- file.path(src, "BAFF_neighbors_genesets.xlsx")
if (requireNamespace("readxl", quietly=TRUE) && file.exists(xl)) {
  x <- as.data.table(readxl::read_excel(xl))
  setnames(x, 1:2, c("symbol","gene"))
  x <- x[gene %in% NB]
  stopifnot(nrow(x) == 55L)

  xk <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE","IMMUNOPROTEASOME_COMPLEX_curated",
          "GOBP_MYELOID_LEUKOCYTE_ACTIVATION")
  for (j in seq_along(xk)) {
    mem <- x[!is.na(x[[2 + j]]) & x[[2 + j]] == "Yes", gene]
    stopifnot(setequal(mem, intersect(gs[[xk[j]]], NB)))
  }
  cat("check: the xlsx membership table agrees with the gene-set tsvs on all 3 sets\n")
} else cat("note : readxl unavailable or xlsx missing -- membership cross-check skipped\n")

cov <- fread("data/clues/covariates_all.txt")[
         , .(gene, logE=log(expression), logL=log(gene_length), gc=gc_content)]
bv  <- fread("data/clues/beta_var_all.txt")[, .(gene, logBV=log(beta_var))]
cov <- merge(cov, bv, by="gene")
cov <- cov[is.finite(logE) & is.finite(logL) & is.finite(gc) & is.finite(logBV)]
cat(sprintf("covariate-complete universe genes: %d of %d | BAFF neighbours kept: %d of %d\n",
            sum(U %in% cov$gene), length(U), sum(NB %in% cov$gene), length(NB)))
stopifnot(all(NB %in% cov$gene))

fit1 <- function(setG, adj) {
  d <- data.table(gene=U, in_hub=as.integer(U %in% NB), in_set=as.integer(U %in% setG))
  if (adj) d <- merge(d, cov, by="gene")
  a <- sum(d$in_hub & d$in_set); b <- sum(d$in_hub & !d$in_set)
  c2<- sum(!d$in_hub & d$in_set); e <- sum(!d$in_hub & !d$in_set)
  out <- list(l2=NA_real_, se=NA_real_, p=NA_real_, n=nrow(d), ov=a, nset=a+c2, sep=FALSE)
  if (adj) {
    if (a == 0 || b == 0) { out$sep <- TRUE; return(out) }
    m <- tryCatch(suppressWarnings(glm(in_set ~ in_hub + logE + logL + gc + logBV,
                                       d, family=binomial)), error=function(e) NULL)
    if (is.null(m) || !m$converged) { out$sep <- TRUE; return(out) }
    co <- coef(summary(m)); if (!("in_hub" %in% rownames(co))) { out$sep <- TRUE; return(out) }
    bb <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
    if (!is.finite(se) || se > 10) { out$sep <- TRUE; return(out) }
  } else {
    if (a==0||b==0||c2==0||e==0) {
      aa<-a+.5; bb2<-b+.5; cc<-c2+.5; dd<-e+.5
      bb <- log((aa*dd)/(bb2*cc)); se <- sqrt(1/aa+1/bb2+1/cc+1/dd); p <- 2*pnorm(-abs(bb/se))
    } else {
      m <- glm(in_set ~ in_hub, d, family=binomial); co <- coef(summary(m))
      bb <- co["in_hub","Estimate"]; se <- co["in_hub","Std. Error"]; p <- co["in_hub","Pr(>|z|)"]
    }
  }
  out$l2 <- bb/log(2); out$se <- se/log(2); out$p <- p
  out
}

res <- rbindlist(lapply(SETS, function(s) rbindlist(lapply(c(FALSE, TRUE), function(adj) {
  r <- fit1(gs[[s$key]], adj)
  data.table(gene_set=s$key, lab=s$lab, mode=ifelse(adj,"Adjusted","Baseline"),
             n=r$n, nSet=r$nset, overlap=r$ov, log2OR=r$l2, se=r$se, p=r$p, unstable=r$sep)
}))))
res[, q     := p.adjust(p, "BH"), by=mode]
res[, stars := fifelse(is.na(q), "", fifelse(q<.001,"***", fifelse(q<.01,"**",
                       fifelse(q<.05,"*",""))))]
res[, OR := 2^log2OR]

res[, npar := fifelse(mode == "Adjusted", 5L, 1L)]
res[, epv  := nSet/npar]

bl <- merge(res[mode=="Baseline", .(gene_set, OR, p_wald=p)],
            up[, .(gene_set, odds_ratio, p_fisher=p_value)], by="gene_set")
bl[, rel := abs(OR - odds_ratio)/odds_ratio]
cat(sprintf("check: baseline OR vs the upstream Fisher OR, max rel diff = %.2f%%\n", 100*max(bl$rel)))
stopifnot(max(bl$rel) < 0.02)
cat("\n*** same 2x2, very different p: Wald is anticonservative in this regime ***\n")
print(bl[, .(gene_set, OR_glm=round(OR,1), OR_fisher=round(odds_ratio,1),
             p_wald=signif(p_wald,2), p_fisher=signif(p_fisher,2))], row.names=FALSE)

fwrite(res[, .(gene_set, mode, n_universe=n, nSet, overlap, log2OR=round(log2OR,3),
               OR=round(OR,2), se=round(se,3), p_wald=p, q_BH=q, stars,
               n_params=npar, events_per_var=round(epv,1), unstable)],
       file.path(outdir, "CLUES_BAFF_neighbors_enrichment.csv"))
cat("\n===== BAFF neighbours (n = 55) vs gene set, project-standard glm =====\n")
print(res[, .(gene_set, mode, n_universe=n, nSet, overlap, log2OR=round(log2OR,2),
              OR=round(OR,1), se=round(se,2), p_wald=signif(p,2), q=signif(q,2), stars,
              epv=round(epv,1))], row.names=FALSE)
