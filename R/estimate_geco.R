#' Estimate GeCo across all gene pairs
#'
#' Loops over every pair of genes, applies the OP rescaling
#' (subtracting `OPmean` and dividing by `OPsd`), assembles the
#' per-pair `Y` and `VS` tensors, and runs [GeCo()]. Returns long
#' `data.frame`s of covariance and SE entries indexed by gene pair.
#'
#' @param data A `data.table` returned by [load_h5ad()].
#' @param genes Character vector of genes to include.
#' @param opms `data.frame` with columns `g, OPmean, OPsd` (from
#'   [compute_OPms()]).
#' @param fixed_cov Optional individual-level covariate matrix
#'   (or `NA` to skip the regression step).
#' @param ind_col,ct_col Names of the individual and cell-type columns.
#' @param cell_types Optional ordering of cell types.
#' @param verbose Print progress every `verbose` pairs (0 = silent).
#'
#' @details The per-(individual, cell type) mean tensor `Y` and the
#' within-group variance tensor `VS_diag` are always computed here from
#' `data` (via [compute_Y()] / [compute_VS_diag()]); they cannot be
#' supplied from outside, because the estimator requires them to come
#' from exactly the same cells as the per-pair covariance.
#'
#' @section Errors instead of NA:
#' Nothing in this loop is allowed to fail silently. Before the loop,
#' the data are checked for a complete individual x cell-type grid
#' with >= 2 cells per group (via [compute_Y()] / [compute_VS_diag()]),
#' and every gene must have a finite, non-zero `OPsd`. Inside the
#' loop, any error from [GeCo()] is re-raised with the gene pair named.
#'
#' @return A list with two `data.frame`s, `cov` and `se`, each having
#'   `2 + 3 + 3*C` columns: `g1, g2, shared11, shared12, shared22,
#'   g11_c<i>, g12_c<i>, g22_c<i>` for each cell-type index `i`.
#' @export
estimate_geco <- function(data, genes, opms,
                          fixed_cov = NA,
                          ind_col = "ind_cov",
                          ct_col  = "cg_cov",
                          cell_types = NULL,
                          verbose = 0) {

  dt <- data.table::as.data.table(data)
  .check_group_col(dt, ind_col, "individual", "estimate_geco")
  .check_group_col(dt, ct_col,  "cell-type",  "estimate_geco")

  if (is.null(cell_types)) cell_types <- as.character(unique(dt[[ct_col]]))
  individuals <- as.character(unique(dt[[ind_col]]))
  n_ind <- length(individuals)
  n_ct  <- length(cell_types)
  G     <- length(genes)

  if (n_ct < 2L) {
    stop("estimate_geco(): only ", n_ct, " cell type; GECO needs at least 2.",
         call. = FALSE)
  }
  if (G < 2L) stop("estimate_geco(): need at least 2 genes.", call. = FALSE)

  # Pre-compute Y for all genes (means per ind x gene x ct)
  Y_obj <- compute_Y(dt, genes,
                     ind_col = ind_col, ct_col = ct_col,
                     individuals = individuals, cell_types = cell_types)
  Y_all <- Y_obj$Y  # [n_ind, n_gene, n_ct]

  # Pre-compute VS diagonal (per-(ind, ct) variance) once for all genes.
  # The off-diagonal cov-12 is computed inside the per-pair loop below.
  message("[estimate_geco] Precomputing VS_diag for ", length(genes), " genes...")
  VS_diag <- compute_VS_diag(dt, genes,
                             ind_col = ind_col, ct_col = ct_col,
                             individuals = individuals,
                             cell_types  = cell_types)

  # Per-gene mean/sd lookup
  opms_dt <- data.table::as.data.table(opms)
  if (!all(c("g", "OPmean", "OPsd") %in% colnames(opms_dt))) {
    stop("estimate_geco(): `opms` must have columns g, OPmean, OPsd ",
         "(see compute_OPms()).", call. = FALSE)
  }
  miss_op <- setdiff(genes, opms_dt$g)
  if (length(miss_op) > 0L) {
    stop("estimate_geco(): ", length(miss_op), " gene(s) have no row in ",
         "`opms`: ", .preview(miss_op), ".", call. = FALSE)
  }
  op_bad <- opms_dt$g[opms_dt$g %in% genes &
                        (!is.finite(opms_dt$OPsd) | opms_dt$OPsd == 0)]
  if (length(op_bad) > 0L) {
    stop("estimate_geco(): ", length(op_bad), " gene(s) have zero or NA ",
         "`OPsd` and cannot be standardized: ", .preview(op_bad),
         ". Drop them from `genes` (see compute_OPms()).", call. = FALSE)
  }
  if (anyDuplicated(opms_dt$g)) {
    stop("estimate_geco(): `opms$g` has duplicated gene IDs: ",
         .preview(unique(opms_dt$g[duplicated(opms_dt$g)])), ".", call. = FALSE)
  }
  data.table::setkey(opms_dt, g)

  n_pairs <- choose(G, 2)

  # Output column names
  shared_names <- c("shared11", "shared12", "shared22")
  ct_block <- function(prefix, c_idx) {
    c(paste0(prefix, "11_c", c_idx),
      paste0(prefix, "12_c", c_idx),
      paste0(prefix, "22_c", c_idx))
  }
  ct_names <- unlist(lapply(seq_len(n_ct), function(c_idx) ct_block("g", c_idx)))
  out_cols <- c("g1", "g2", shared_names, ct_names)

  res_cov <- matrix(NA, nrow = n_pairs, ncol = length(out_cols))
  res_se  <- matrix(NA, nrow = n_pairs, ncol = length(out_cols))
  colnames(res_cov) <- out_cols
  colnames(res_se)  <- out_cols

  lt_shared <- lower.tri(matrix(1, 2, 2), diag = TRUE)
  idxs <- split(seq_len(2 * n_ct), rep(seq_len(n_ct), each = 2))

  ## Re-raise any per-pair error with the gene pair named, so a failure
  ## deep inside GeCo() is attributable.
  with_pair <- function(expr, g1, g2) {
    tryCatch(expr, error = function(e) {
      stop("estimate_geco(): gene pair (", g1, ", ", g2, ") failed: ",
           conditionMessage(e), call. = FALSE)
    })
  }

  ## One-time grouping of the cells for the fast per-pair covariance
  ## (also enforces >= 2 cells per group). Gene columns are pulled out
  ## once as plain numeric vectors.
  setup <- .group_setup(dt, ind_col, ct_col, individuals, cell_types,
                        where = "estimate_geco")
  gene_vec <- function(g) {
    x <- dt[[g]]
    if (anyNA(x)) {
      stop("estimate_geco(): gene '", g, "' has NA expression in ",
           sum(is.na(x)), " cell(s).", call. = FALSE)
    }
    as.numeric(x)
  }

  k <- 0L
  for (i in seq_len(G - 1)) {
    g1 <- genes[i]
    ms1 <- opms_dt[g1]; m1 <- ms1$OPmean; s1 <- ms1$OPsd

    Y1raw <- Y_all[, i, , drop = FALSE]; dim(Y1raw) <- c(n_ind, n_ct)
    Y1 <- (Y1raw - m1) / s1
    x1 <- gene_vec(g1)

    for (j in (i + 1):G) {
      k <- k + 1L
      g2 <- genes[j]

      if (verbose > 0 && (k %% verbose) == 0) {
        message(sprintf("[estimate_geco] pair %d / %d : %s, %s",
                        k, n_pairs, g1, g2))
      }

      ms2 <- opms_dt[g2]; m2 <- ms2$OPmean; s2 <- ms2$OPsd

      Y2raw <- Y_all[, j, , drop = FALSE]; dim(Y2raw) <- c(n_ind, n_ct)
      Y2 <- (Y2raw - m2) / s2

      cov12 <- .group_cov_pair(x1, gene_vec(g2), Y1raw, Y2raw, setup)

      VS <- array(NA_real_, dim = c(n_ind, n_ct, 2, 2))
      VS[, , 1, 1] <- VS_diag[, , i] / s1^2
      VS[, , 2, 2] <- VS_diag[, , j] / s2^2
      VS[, , 1, 2] <- cov12 / s1 / s2
      VS[, , 2, 1] <- VS[, , 1, 2]

      Y_pair <- array(NA_real_, dim = c(n_ind, 2, n_ct))
      Y_pair[, 1, ] <- Y1
      Y_pair[, 2, ] <- Y2

      res <- with_pair(GeCo(Y = Y_pair, VS = VS, fixed_cov = fixed_cov),
                       g1, g2)

      res_cov[k, 1:2] <- c(g1, g2)
      res_se [k, 1:2] <- c(g1, g2)

      Vshared  <- res$V_shared;  SEshared <- res$SE_shared
      Vgc     <- res$V_GC;     SEgc    <- res$SE_GC

      res_cov[k, 3:5] <- as.numeric(Vshared[lt_shared])
      res_se [k, 3:5] <- as.numeric(SEshared[lt_shared])

      lower_cov <- lapply(idxs, function(id) {
        subm <- Vgc[id, id]; subm[lower.tri(subm, diag = TRUE)]
      })
      lower_se <- lapply(idxs, function(id) {
        subm <- SEgc[id, id]; subm[lower.tri(subm, diag = TRUE)]
      })
      res_cov[k, 6:(5 + 3 * n_ct)] <- as.numeric(unlist(lower_cov))
      res_se [k, 6:(5 + 3 * n_ct)] <- as.numeric(unlist(lower_se))
    }
  }

  cov_df <- as.data.frame(res_cov, stringsAsFactors = FALSE)
  se_df  <- as.data.frame(res_se,  stringsAsFactors = FALSE)

  # numeric columns
  for (col in out_cols[-c(1, 2)]) {
    cov_df[[col]] <- as.numeric(cov_df[[col]])
    se_df [[col]] <- as.numeric(se_df [[col]])
  }

  list(cov = cov_df, se = se_df, cell_types = cell_types)
}
