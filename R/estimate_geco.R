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

  if (is.null(cell_types)) cell_types <- as.character(unique(dt[[ct_col]]))
  individuals <- as.character(unique(dt[[ind_col]]))
  n_ind <- length(individuals)
  n_ct  <- length(cell_types)

  # Pre-compute Y for all genes
  Y_obj <- compute_Y(dt, genes,
                     ind_col = ind_col, ct_col = ct_col,
                     individuals = individuals, cell_types = cell_types)
  Y_all <- Y_obj$Y  # [n_ind, n_gene, n_ct]

  # Per-gene mean/sd lookup
  opms_dt <- data.table::as.data.table(opms)
  data.table::setkey(opms_dt, g)

  G <- length(genes)
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

  k <- 0L
  for (i in seq_len(G - 1)) {
    g1 <- genes[i]
    ms1 <- opms_dt[g1]; m1 <- ms1$OPmean; s1 <- ms1$OPsd
    if (is.na(s1) || s1 == 0) next

    Y1 <- Y_all[, i, , drop = FALSE]; dim(Y1) <- c(n_ind, n_ct)
    Y1 <- (Y1 - m1) / s1

    for (j in (i + 1):G) {
      k <- k + 1L
      g2 <- genes[j]

      if (verbose > 0 && (k %% verbose) == 0) {
        message(sprintf("[estimate_geco] pair %d / %d : %s, %s",
                        k, n_pairs, g1, g2))
      }

      ms2 <- opms_dt[g2]; m2 <- ms2$OPmean; s2 <- ms2$OPsd
      if (is.na(s2) || s2 == 0) next

      Y2 <- Y_all[, j, , drop = FALSE]; dim(Y2) <- c(n_ind, n_ct)
      Y2 <- (Y2 - m2) / s2

      VS <- try(compute_VS_pair(dt, g1, g2,
                                ind_col = ind_col, ct_col = ct_col,
                                individuals = individuals,
                                cell_types  = cell_types),
                silent = TRUE)
      if (inherits(VS, "try-error")) next

      VS[, , 1, 1] <- VS[, , 1, 1] / s1^2
      VS[, , 2, 2] <- VS[, , 2, 2] / s2^2
      VS[, , 1, 2] <- VS[, , 1, 2] / s1 / s2
      VS[, , 2, 1] <- VS[, , 2, 1] / s1 / s2

      Y_pair <- array(NA_real_, dim = c(n_ind, 2, n_ct))
      Y_pair[, 1, ] <- Y1
      Y_pair[, 2, ] <- Y2

      res <- try(GeCo(Y = Y_pair, VS = VS, fixed_cov = fixed_cov),
                 silent = TRUE)
      if (inherits(res, "try-error")) next

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
