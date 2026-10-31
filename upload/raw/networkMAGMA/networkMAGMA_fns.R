library(Matrix)
library(dplyr)
library(tidyr)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)

.map_entrez_to_ensg <- function(entrez_ids) {
  entrez <- as.character(unique(entrez_ids))
  m <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = entrez,
    columns = c("ENTREZID", "ENSEMBL"),
    keytype = "ENTREZID"
  )
  m %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarise(
      ENSG = paste(unique(na.omit(ENSEMBL)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(GENE = ENTREZID)
}

.compute_gene_stats <- function(magma_path) {
  raw <- readr::read_table(
    magma_path,
    col_types = readr::cols(.default = readr::col_character())
  ) %>%
    dplyr::mutate(
      GENE  = as.character(GENE),
      P     = as.numeric(P),
      CHR   = as.character(CHR),
      START = suppressWarnings(as.numeric(START))
    ) %>%
    dplyr::filter(!is.na(P)) %>%
    dplyr::mutate(
      neg_log10_P = -log10(pmax(P, .Machine$double.xmin))
    )

  emap <- .map_entrez_to_ensg(raw$GENE)

  raw %>%
    dplyr::inner_join(emap, by = "GENE") %>%
    tidyr::separate_rows(ENSG, sep = ";") %>%
    dplyr::filter(nzchar(ENSG), startsWith(ENSG, "ENSG")) %>%
    dplyr::arrange(P) %>%
    dplyr::distinct(ENSG, .keep_all = TRUE)
}

.mhc_ensg_set <- function(gene_stats, mhc_chr, mhc_start, mhc_end, mhc_genes) {
  if (!is.null(mhc_genes)) {
    return(unique(as.character(mhc_genes)))
  }
  if (!all(c("CHR", "START") %in% colnames(gene_stats))) {
    stop("gene_stats lacks CHR/START columns; cannot compute MHC by coordinates. ",
         "Either pass mhc_genes explicitly, or ensure the MAGMA .genes.out file ",
         "includes CHR and START columns.")
  }
  gene_stats %>%
    dplyr::filter(
      as.character(CHR) == as.character(mhc_chr),
      !is.na(START),
      START >= mhc_start,
      START <= mhc_end
    ) %>%
    dplyr::pull(ENSG) %>%
    unique() %>%
    as.character()
}

.general_G_perm <- function(W, x, n_perm, seed) {
  good <- is.finite(x)
  if (!all(good)) { W <- W[good, good]; x <- x[good] }
  n <- length(x)
  if (n < 3) stop("Need at least three finite x values.")

  S0  <- sum(W)
  m1  <- sum(x); m2 <- sum(x^2)
  Den <- m1^2 - m2
  if (Den == 0) stop("Denominator (sum(x)^2 - sum(x^2)) is zero.")
  num   <- as.numeric(crossprod(x, W %*% x))
  G_obs <- num / Den
  E_G   <- S0 / (n * (n - 1))

  set.seed(seed)
  G_perm <- numeric(n_perm)
  for (k in seq_len(n_perm)) {
    xp <- sample(x)
    G_perm[k] <- as.numeric(crossprod(xp, W %*% xp)) / Den
  }
  mu_p <- mean(G_perm); sd_p <- sd(G_perm)
  Z_perm <- if (is.finite(sd_p) && sd_p > 0) (G_obs - mu_p) / sd_p else NA_real_

  p_perm_right <- (sum(G_perm >= G_obs) + 1) / (n_perm + 1)
  p_perm_two   <- (sum(abs(G_perm - mu_p) >= abs(G_obs - mu_p)) + 1) / (n_perm + 1)

  list(N = n, S0 = S0, G = G_obs, E_G = E_G,
       mean_G_perm = mu_p, sd_G_perm = sd_p, Z_perm = Z_perm,
       p_perm_one_sided = p_perm_right, p_perm_two_sided = p_perm_two)
}

.general_G_perm_stratified <- function(W, x, n_perm, seed,
                                        n_bins      = 10,
                                        degree_type = c("unweighted","weighted")) {
  degree_type <- match.arg(degree_type)
  good <- is.finite(x)
  if (!all(good)) { W <- W[good, good]; x <- x[good] }
  n <- length(x)
  if (n < 3) stop("Need at least three finite x values.")

  if (degree_type == "unweighted") {
    deg <- Matrix::rowSums(W != 0)
  } else {
    deg <- as.numeric(Matrix::rowSums(W))
  }

  set.seed(seed)
  ranks <- rank(deg, ties.method = "random")
  bin   <- ceiling(ranks / n * n_bins)
  bin[bin == 0] <- 1L
  bin_idx <- split(seq_len(n), bin)

  S0  <- sum(W)
  m1  <- sum(x); m2 <- sum(x^2)
  Den <- m1^2 - m2
  if (Den == 0) stop("Denominator (sum(x)^2 - sum(x^2)) is zero.")
  num   <- as.numeric(crossprod(x, W %*% x))
  G_obs <- num / Den
  E_G   <- S0 / (n * (n - 1))

  G_perm <- numeric(n_perm)
  for (k in seq_len(n_perm)) {
    xp <- x
    for (idx in bin_idx) {
      xp[idx] <- x[idx][sample.int(length(idx))]
    }
    G_perm[k] <- as.numeric(crossprod(xp, W %*% xp)) / Den
  }
  mu_p <- mean(G_perm); sd_p <- sd(G_perm)
  Z_perm <- if (is.finite(sd_p) && sd_p > 0) (G_obs - mu_p) / sd_p else NA_real_

  p_perm_right <- (sum(G_perm >= G_obs) + 1) / (n_perm + 1)
  p_perm_two   <- (sum(abs(G_perm - mu_p) >= abs(G_obs - mu_p)) + 1) / (n_perm + 1)

  rho_dx <- suppressWarnings(cor(deg, x, method = "spearman"))

  list(N = n, S0 = S0, G = G_obs, E_G = E_G,
       mean_G_perm = mu_p, sd_G_perm = sd_p, Z_perm = Z_perm,
       p_perm_one_sided = p_perm_right, p_perm_two_sided = p_perm_two,
       n_bins = n_bins, degree_type = degree_type,
       spearman_deg_x = rho_dx)
}

networkMAGMA <- function(network, magma_path,
                         n_perm      = 1000,
                         seed        = 111,
                         null        = c("global", "stratified"),
                         n_bins      = 10,
                         degree_type = "unweighted",
                         mhc_filter  = c("keep",
                                         "drop_internal_edges",
                                         "drop_genes"),
                         mhc_chr     = "6",
                         mhc_start   = 25000000,
                         mhc_end     = 35000000,
                         mhc_genes   = NULL,
                         verbose     = TRUE) {

  null       <- match.arg(null)
  mhc_filter <- match.arg(mhc_filter)

  gene_stats <- .compute_gene_stats(magma_path)
  common <- intersect(rownames(network), gene_stats$ENSG)
  if (length(common) < 3) {
    stop("Only ", length(common), " common ENSGs between network and MAGMA result.")
  }

  mhc_set <- character(0)
  if (mhc_filter != "keep") {
    mhc_set <- .mhc_ensg_set(gene_stats, mhc_chr, mhc_start, mhc_end, mhc_genes)
    if (verbose) {
      cat(sprintf("[mhc] %d MHC ENSGs (from %s)\n", length(mhc_set),
                  if (!is.null(mhc_genes)) "user-supplied list"
                  else sprintf("chr%s:%s-%s",
                               mhc_chr, format(mhc_start, big.mark=",", scientific=FALSE),
                               format(mhc_end,   big.mark=",", scientific=FALSE))))
    }
  }

  if (mhc_filter == "drop_genes") {
    keep_common <- setdiff(common, mhc_set)
    n_dropped   <- length(common) - length(keep_common)
    if (verbose) cat(sprintf("[mhc] drop_genes: removed %d MHC genes; %d remain\n",
                             n_dropped, length(keep_common)))
    if (length(keep_common) < 3) {
      stop("After dropping MHC genes only ", length(keep_common),
           " genes remain — too few.")
    }
    common <- keep_common
  }

  W <- abs(network[common, common])
  diag(W) <- 0

  n_mhc_in_common    <- 0L
  n_mhc_edges_zeroed <- 0L
  if (mhc_filter == "drop_internal_edges") {
    is_mhc <- common %in% mhc_set
    n_mhc_in_common <- sum(is_mhc)
    if (n_mhc_in_common >= 2) {
      W_sub_pre <- W[is_mhc, is_mhc]
      n_mhc_edges_zeroed <- as.integer(sum(W_sub_pre != 0) / 2L)
      W[is_mhc, is_mhc] <- 0
      if (inherits(W, "Matrix")) W <- Matrix::drop0(W)
    }
    if (verbose) {
      cat(sprintf("[mhc] drop_internal_edges: %d MHC genes in common, zeroed %d MHC-MHC edges\n",
                  n_mhc_in_common, n_mhc_edges_zeroed))
    }
  }

  x <- gene_stats[match(common, gene_stats$ENSG), 'neg_log10_P', drop = TRUE]

  res <- if (null == "stratified") {
    .general_G_perm_stratified(W, x, n_perm = n_perm, seed = seed,
                                n_bins = n_bins, degree_type = degree_type)
  } else {
    .general_G_perm(W, x, n_perm = n_perm, seed = seed)
  }

  out <- data.frame(
    N                  = res$N,
    observed_G         = res$G,
    Z_perm             = res$Z_perm,
    p_value            = res$p_perm_one_sided,
    p_value_two_sided  = res$p_perm_two_sided,
    mean_G_perm        = res$mean_G_perm,
    sd_G_perm          = res$sd_G_perm,
    E_G                = res$E_G,
    S0                 = res$S0,
    n_perm             = n_perm,
    null               = null,
    n_bins             = if (null == "stratified") n_bins else NA_integer_,
    degree_type        = if (null == "stratified") degree_type else NA_character_,
    spearman_deg_x     = if (null == "stratified") res$spearman_deg_x else NA_real_,
    mhc_filter         = mhc_filter,
    n_mhc_in_common    = n_mhc_in_common,
    n_mhc_edges_zeroed = n_mhc_edges_zeroed,
    stringsAsFactors   = FALSE
  )
  rownames(out) <- NULL
  out
}
