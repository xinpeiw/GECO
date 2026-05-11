## ---------------------------------------------------------------------
##  Cluster / HPC mode
##
##  The in-memory pipeline (run_geco()) keeps the full cell-level
##  expression matrix in RAM, which only works for at most a few
##  thousand genes. For ~10k genes a node typically cannot hold the
##  full matrix, so the original SLURM-based workflow precomputes
##  per-gene-pair statistics in HDF5 chunks and runs them as a SLURM
##  array. This file ports that style into a small set of functions
##  that read/write to a single `work_dir`:
##
##    geco_prepare()             # one-time setup
##    geco_chunk()               # SLURM array (1..n_jobs), gene-pair estimation
##    geco_finalize_estimates()  # merge chunks, reformat, p, PSD
##    geco_precision_chunk()     # SLURM array (1..n_ct+1), glasso+eBIC
##    geco_finalize_outputs()    # degrees + final outputs (incl. degree.xlsx)
##
##  The cluster mode uses the {fst} package for column-selectable
##  on-disk cell data. {fst} is in Suggests; install it before
##  running cluster mode.
## ---------------------------------------------------------------------


#' Prepare a working directory for the cluster pipeline
#'
#' One-time setup that reads the .h5ad, writes the cell-level data to
#' an `fst` file (so chunks can pull only the columns they need),
#' computes per-gene `OPmean`/`OPsd`, builds the fixed-covariate
#' matrix, and saves the gene list, gene pairs, and per-pair
#' bookkeeping metadata to `work_dir`.
#'
#' @param h5ad_file Path to the .h5ad file.
#' @param work_dir  Directory to use. Created if missing.
#' @param genes Character vector of gene IDs to extract from the
#'   .h5ad (required). See [run_geco()]/[load_h5ad()] for details and
#'   [read_gene_list()] for loading from a file.
#' @param ind_col,ct_col See [run_geco()].
#' @param subset_col,subset_value Optional metadata filter (see
#'   [load_h5ad()]).
#' @param fixed_covariates Which covariates to regress out before
#'   GeCo estimation. See [run_geco()] for accepted forms (a character
#'   vector of catalog names is the common case, e.g.
#'   `c("age_group", "sex", "PC1")`).
#' @param Y_path Optional path to a precomputed Y tensor `.rds`
#'   (shape `[n_ind, n_gene, n_ct]`). If supplied, copied into
#'   `work_dir/Y.rds` instead of recomputing from the .h5ad.
#' @param VS_diag_path Optional path to a precomputed VS_diag tensor
#'   `.rds` (shape `[n_ind, n_ct, n_gene]`). If supplied, copied into
#'   `work_dir/VS_diag.rds` instead of recomputing.
#'
#' @return Invisibly, the prepared `prep` metadata list (also saved to
#'   `work_dir/prep.rds`).
#' @export
geco_prepare <- function(h5ad_file,
                         work_dir,
                         genes,
                         ind_col = "ind_cov",
                         ct_col  = "cg_cov",
                         subset_col = NULL, subset_value = NULL,
                         fixed_covariates = NULL,
                         Y_path = NULL,
                         VS_diag_path = NULL) {

  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Cluster mode requires the 'fst' package. install.packages('fst').")
  }
  if (missing(genes) || is.null(genes)) {
    stop("`genes` is required. Pass a character vector of gene IDs ",
         "(or use read_gene_list() to load one from a .txt file).",
         call. = FALSE)
  }
  if (length(genes) < 2L) stop("Need at least 2 genes.")
  genes <- unique(as.character(genes))

  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(work_dir, "chunks"),    showWarnings = FALSE)
  dir.create(file.path(work_dir, "precision"), showWarnings = FALSE)

  message("[geco_prepare] Loading h5ad: ", h5ad_file)
  data <- load_h5ad(h5ad_file,
                    genes = genes,
                    subset_col = subset_col, subset_value = subset_value)

  message("[geco_prepare] Writing cell-level data to fst (",
          nrow(data), " cells x ", length(genes), " genes)...")
  fst::write_fst(as.data.frame(data),
                 path = file.path(work_dir, "cell_data.fst"),
                 compress = 50)

  message("[geco_prepare] Computing OPms...")
  opms <- compute_OPms(data, genes, ind_col = ind_col)

  fixed_cov <- resolve_fixed_covariates(data, genes,
                                        fixed_covariates = fixed_covariates,
                                        ind_col = ind_col)

  individuals <- as.character(unique(data[[ind_col]]))
  cts_raw     <- as.character(unique(data[[ct_col]]))   # straight from h5ad

  ## Precompute Y (per-(ind, gene, ct) means) and VS_diag (per-(ind, ct,
  ## gene) within-group variance). Doing both once here means every
  ## chunk worker only has to compute the off-diagonal cov12 per pair.
  Y_out_path <- file.path(work_dir, "Y.rds")
  if (!is.null(Y_path)) {
    message("[geco_prepare] Copying user-supplied Y from ", Y_path)
    file.copy(Y_path, Y_out_path, overwrite = TRUE)
  } else {
    message("[geco_prepare] Computing Y tensor [n_ind=", length(individuals),
            ", n_gene=", length(genes), ", n_ct=", length(cts_raw), "]...")
    Y_obj <- compute_Y(data, genes,
                       ind_col = ind_col, ct_col = ct_col,
                       individuals = individuals, cell_types = cts_raw)
    saveRDS(Y_obj$Y, Y_out_path)
  }

  VS_out_path <- file.path(work_dir, "VS_diag.rds")
  if (!is.null(VS_diag_path)) {
    message("[geco_prepare] Copying user-supplied VS_diag from ", VS_diag_path)
    file.copy(VS_diag_path, VS_out_path, overwrite = TRUE)
  } else {
    message("[geco_prepare] Computing VS_diag tensor [n_ind=",
            length(individuals), ", n_ct=", length(cts_raw),
            ", n_gene=", length(genes), "]...")
    VS_diag <- compute_VS_diag(data, genes,
                               ind_col = ind_col, ct_col = ct_col,
                               individuals = individuals,
                               cell_types  = cts_raw)
    saveRDS(VS_diag, VS_out_path)
  }

  gene_pairs <- utils::combn(genes, 2)

  prep <- list(
    h5ad_file       = normalizePath(h5ad_file, mustWork = FALSE),
    work_dir        = normalizePath(work_dir, mustWork = FALSE),
    genes           = genes,
    individuals     = individuals,
    cell_types_raw  = cts_raw,
    cell_type_names = cts_raw,
    ind_col         = ind_col,
    ct_col          = ct_col,
    fixed_cov       = fixed_cov,
    opms            = opms,
    n_pairs         = ncol(gene_pairs)
  )

  saveRDS(prep,       file.path(work_dir, "prep.rds"))
  saveRDS(gene_pairs, file.path(work_dir, "gene_pairs.rds"))
  message("[geco_prepare] Wrote ", work_dir,
          "/{cell_data.fst, prep.rds, gene_pairs.rds, Y.rds, VS_diag.rds}. ",
          ncol(gene_pairs), " gene pairs total.")
  invisible(prep)
}


# Internal: which gene-pair indices belong to a given (task_id, n_jobs)?
.chunk_indices <- function(n_pairs, task_id, n_jobs) {
  if (task_id < 1 || task_id > n_jobs) {
    stop("task_id must be in 1..", n_jobs)
  }
  base  <- n_pairs %/% n_jobs
  extra <- n_pairs %%  n_jobs
  # First `extra` chunks get one extra pair.
  sizes <- rep(base, n_jobs) + c(rep(1L, extra), rep(0L, n_jobs - extra))
  starts <- c(0L, cumsum(sizes))[seq_len(n_jobs)] + 1L
  ends   <- cumsum(sizes)
  list(start = starts[task_id], end = ends[task_id], size = sizes[task_id])
}


#' Run one chunk of the gene-pair estimation (SLURM array element)
#'
#' Reads `prep.rds` + the relevant gene columns from `cell_data.fst`,
#' computes `compute_VS_pair()` and `GeCo()` for every pair in this
#' chunk, and writes `work_dir/chunks/chunk_<task_id>.rds` containing
#' a list of two long-format data.frames (`cov` and `se`).
#'
#' @param work_dir Working directory built by [geco_prepare()].
#' @param task_id  Integer in `1..n_jobs`. Often supplied via
#'   `commandArgs()` / `$SLURM_ARRAY_TASK_ID`.
#' @param n_jobs   Total number of chunks.
#' @param verbose  Print progress every `verbose` pairs (0 = silent).
#'
#' @return Invisibly the chunk's results list.
#' @export
geco_chunk <- function(work_dir, task_id, n_jobs, verbose = 100) {

  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Cluster mode requires the 'fst' package.")
  }

  prep       <- readRDS(file.path(work_dir, "prep.rds"))
  gene_pairs <- readRDS(file.path(work_dir, "gene_pairs.rds"))
  Y_all      <- readRDS(file.path(work_dir, "Y.rds"))        # [n_ind, n_gene, n_ct]
  VS_diag    <- readRDS(file.path(work_dir, "VS_diag.rds"))  # [n_ind, n_ct, n_gene]
  cell_fst   <- file.path(work_dir, "cell_data.fst")

  rng <- .chunk_indices(prep$n_pairs, task_id, n_jobs)
  message(sprintf("[geco_chunk] task %d / %d: pairs %d..%d (%d total)",
                  task_id, n_jobs, rng$start, rng$end, rng$size))

  ind_col <- prep$ind_col;  ct_col <- prep$ct_col
  inds    <- prep$individuals
  cts     <- prep$cell_types_raw
  n_ind   <- length(inds);  n_ct   <- length(cts)

  ## Map gene name -> column index in Y_all / VS_diag (same order as prep$genes)
  gene_to_idx <- stats::setNames(seq_along(prep$genes), prep$genes)

  opms_dt <- data.table::as.data.table(prep$opms)
  data.table::setkey(opms_dt, g)

  shared_names <- c("shared11", "shared12", "shared22")
  ct_names <- unlist(lapply(seq_len(n_ct), function(c_idx)
    c(paste0("g11_c", c_idx), paste0("g12_c", c_idx), paste0("g22_c", c_idx))))
  out_cols <- c("g1", "g2", shared_names, ct_names)

  res_cov <- matrix(NA, nrow = rng$size, ncol = length(out_cols))
  res_se  <- matrix(NA, nrow = rng$size, ncol = length(out_cols))
  colnames(res_cov) <- out_cols
  colnames(res_se)  <- out_cols

  lt_shared <- lower.tri(matrix(1, 2, 2), diag = TRUE)
  idxs <- split(seq_len(2 * n_ct), rep(seq_len(n_ct), each = 2))

  # We need ind_col + ct_col every time; pre-read them.
  meta_df <- fst::read_fst(cell_fst, columns = c(ind_col, ct_col))
  meta_df <- data.table::as.data.table(meta_df)

  for (k in seq_len(rng$size)) {
    j  <- rng$start + k - 1L
    g1 <- gene_pairs[1, j]
    g2 <- gene_pairs[2, j]

    if (verbose > 0 && (k %% verbose) == 0) {
      message(sprintf("[geco_chunk] pair %d / %d : %s, %s",
                      k, rng$size, g1, g2))
    }

    ms1 <- opms_dt[g1]; ms2 <- opms_dt[g2]
    m1 <- ms1$OPmean; s1 <- ms1$OPsd
    m2 <- ms2$OPmean; s2 <- ms2$OPsd
    if (is.na(s1) || s1 == 0 || is.na(s2) || s2 == 0) next

    idx1 <- gene_to_idx[[g1]]
    idx2 <- gene_to_idx[[g2]]
    if (is.null(idx1) || is.null(idx2)) next

    ## Pull pre-computed Y means for these two genes; rescale by OP.
    Y1 <- Y_all[, idx1, , drop = FALSE]; dim(Y1) <- c(n_ind, n_ct)
    Y2 <- Y_all[, idx2, , drop = FALSE]; dim(Y2) <- c(n_ind, n_ct)
    Y1 <- (Y1 - m1) / s1
    Y2 <- (Y2 - m2) / s2
    Y_pair <- array(NA_real_, dim = c(n_ind, 2, n_ct))
    Y_pair[, 1, ] <- Y1; Y_pair[, 2, ] <- Y2

    ## Read just the two gene columns for the off-diagonal cov12.
    dat <- try({
      gene_df <- fst::read_fst(cell_fst, columns = c(g1, g2))
      cbind(meta_df, data.table::as.data.table(gene_df))
    }, silent = TRUE)
    if (inherits(dat, "try-error")) next

    cov12 <- try(compute_VS_cov_pair(dat, g1, g2,
                                     ind_col = ind_col, ct_col = ct_col,
                                     individuals = inds, cell_types = cts),
                 silent = TRUE)
    if (inherits(cov12, "try-error")) next

    VS <- array(NA_real_, dim = c(n_ind, n_ct, 2, 2))
    VS[, , 1, 1] <- VS_diag[, , idx1] / s1^2
    VS[, , 2, 2] <- VS_diag[, , idx2] / s2^2
    VS[, , 1, 2] <- cov12 / s1 / s2
    VS[, , 2, 1] <- VS[, , 1, 2]

    res <- try(GeCo(Y = Y_pair, VS = VS, fixed_cov = prep$fixed_cov),
               silent = TRUE)
    if (inherits(res, "try-error")) next

    res_cov[k, 1:2] <- c(g1, g2)
    res_se [k, 1:2] <- c(g1, g2)

    res_cov[k, 3:5] <- as.numeric(res$V_shared [lt_shared])
    res_se [k, 3:5] <- as.numeric(res$SE_shared[lt_shared])
    lower_cov <- lapply(idxs, function(id) {
      subm <- res$V_GC [id, id]; subm[lower.tri(subm, diag = TRUE)]
    })
    lower_se <- lapply(idxs, function(id) {
      subm <- res$SE_GC[id, id]; subm[lower.tri(subm, diag = TRUE)]
    })
    res_cov[k, 6:(5 + 3 * n_ct)] <- as.numeric(unlist(lower_cov))
    res_se [k, 6:(5 + 3 * n_ct)] <- as.numeric(unlist(lower_se))
  }

  cov_df <- as.data.frame(res_cov, stringsAsFactors = FALSE)
  se_df  <- as.data.frame(res_se,  stringsAsFactors = FALSE)
  for (col in out_cols[-c(1, 2)]) {
    cov_df[[col]] <- as.numeric(cov_df[[col]])
    se_df [[col]] <- as.numeric(se_df [[col]])
  }

  out_path <- file.path(work_dir, "chunks",
                        sprintf("chunk_%05d.rds", task_id))
  saveRDS(list(cov = cov_df, se = se_df, task_id = task_id,
               start = rng$start, end = rng$end),
          out_path)
  message("[geco_chunk] Wrote ", out_path)
  invisible(list(cov = cov_df, se = se_df))
}


#' Merge chunks and reformat into per-cell-type matrices
#'
#' Concatenates all `work_dir/chunks/chunk_*.rds` files, calls
#' [reformat_to_matrices()], computes p-value matrices, applies
#' [transform_psd()], and saves a single `work_dir/mats.rds`
#' containing all the inputs that downstream precision/degree steps
#' need.
#'
#' @param work_dir Working directory built by [geco_prepare()].
#' @param psd_thr Threshold for [transform_psd()] (default 0).
#'
#' @return Invisibly the merged result list.
#' @export
geco_finalize_estimates <- function(work_dir, psd_thr = 0) {

  prep <- readRDS(file.path(work_dir, "prep.rds"))
  chunk_files <- list.files(file.path(work_dir, "chunks"),
                            pattern = "^chunk_.*\\.rds$",
                            full.names = TRUE)
  if (length(chunk_files) == 0L) stop("No chunk files found in ", work_dir)
  message("[geco_finalize_estimates] Merging ", length(chunk_files), " chunks...")

  cov_parts <- vector("list", length(chunk_files))
  se_parts  <- vector("list", length(chunk_files))
  for (i in seq_along(chunk_files)) {
    chunk <- readRDS(chunk_files[i])
    cov_parts[[i]] <- chunk$cov
    se_parts [[i]] <- chunk$se
  }
  cov_df <- do.call(rbind, cov_parts)
  se_df  <- do.call(rbind, se_parts)

  est <- list(cov = cov_df, se = se_df,
              cell_types = prep$cell_types_raw)
  mats <- reformat_to_matrices(est, prep$genes,
                               cell_types = prep$cell_type_names)

  message("[geco_finalize_estimates] Computing p-value + PSD matrices...")
  p_mats   <- mapply(compute_p_matrix, mats$cov, mats$se, SIMPLIFY = FALSE)
  names(p_mats) <- names(mats$cov)
  cov_psd  <- lapply(mats$cov, transform_psd, thr = psd_thr)

  out <- list(cov = mats$cov, se = mats$se,
              p = p_mats, cov_psd = cov_psd,
              cell_type_names = prep$cell_type_names,
              n_ind = length(prep$individuals),
              genes = prep$genes)
  saveRDS(out, file.path(work_dir, "mats.rds"))
  message("[geco_finalize_estimates] Wrote ", file.path(work_dir, "mats.rds"))
  invisible(out)
}


#' Fit the precision matrix for a single cell type (SLURM array element)
#'
#' Loads `work_dir/mats.rds`, picks the cell type at index `ct_idx`
#' (1 = "shared", then cell types in the order saved by
#' [geco_prepare()]), and runs [fit_precision()]. The result is saved
#' to `work_dir/precision/precision_<name>.rds`.
#'
#' @param work_dir Working directory built by [geco_prepare()].
#' @param ct_idx   Integer index into the ordered list
#'   `c("shared", <cell types from h5ad>)`. So valid values are
#'   `1..(n_cell_types + 1)`.
#' @param fit_args Named list of arguments forwarded to
#'   [fit_precision()].
#'
#' @return Invisibly the precision fit result.
#' @export
geco_precision_chunk <- function(work_dir, ct_idx, fit_args = list()) {
  mats   <- readRDS(file.path(work_dir, "mats.rds"))
  all_nm <- names(mats$cov)
  if (ct_idx < 1 || ct_idx > length(all_nm)) {
    stop("ct_idx must be in 1..", length(all_nm))
  }
  nm <- all_nm[ct_idx]
  message("[geco_precision_chunk] Fitting precision for: ", nm)

  fit <- do.call(fit_precision,
    c(list(cov = mats$cov_psd[[nm]],
           p = mats$p[[nm]],
           n_ind = mats$n_ind),
      fit_args))
  dimnames(fit$precision) <- list(mats$genes, mats$genes)

  out_path <- file.path(work_dir, "precision",
                        paste0("precision_", nm, ".rds"))
  saveRDS(fit, out_path)
  message("[geco_precision_chunk] Wrote ", out_path,
          " (lambda = ", fit$lambda, ")")
  invisible(fit)
}


#' Finalize cluster outputs: save matrices, build degree.xlsx
#'
#' Reads `work_dir/mats.rds` and `work_dir/precision/*.rds`, writes
#' per-cell-type `cov/se/p/precision_<ct>.rds` to `output_dir`, and
#' assembles a single `degree.xlsx` workbook with one column per cell
#' type.
#'
#' @param work_dir Working directory built by [geco_prepare()].
#' @param output_dir Final outputs directory (created if missing).
#'
#' @return Invisibly the final result list.
#' @export
geco_finalize_outputs <- function(work_dir, output_dir) {
  mats <- readRDS(file.path(work_dir, "mats.rds"))
  pdir <- file.path(work_dir, "precision")
  precision_files <- list.files(pdir,
                                pattern = "^precision_.*\\.rds$",
                                full.names = TRUE)
  if (length(precision_files) == 0L) {
    stop("No precision files found in ", pdir)
  }

  precisions <- list()
  lambdas    <- list()
  eBIC_paths <- list()
  for (pf in precision_files) {
    nm  <- sub("^precision_(.*)\\.rds$", "\\1", basename(pf))
    fit <- readRDS(pf)
    precisions[[nm]] <- fit$precision
    lambdas[[nm]]    <- fit$lambda
    eBIC_paths[[nm]] <- fit$eBIC_path
  }

  genes <- mats$genes
  deg_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (nm in names(mats$cov)) {
    if (is.null(precisions[[nm]])) {
      warning("Missing precision for cell type '", nm,
              "' - degree column set to NA.")
      deg_df[[nm]] <- NA_integer_
    } else {
      deg <- compute_degree(precisions[[nm]])
      deg_df[[nm]] <- as.integer(deg[genes])
    }
  }

  result <- list(
    cov        = mats$cov,
    se         = mats$se,
    p          = mats$p,
    cov_psd    = mats$cov_psd,
    precision  = precisions,
    lambda     = lambdas,
    eBIC_paths = eBIC_paths,
    degree     = deg_df,
    genes      = genes,
    cell_types = mats$cell_type_names
  )

  save_results(result, output_dir)
  message("[geco_finalize_outputs] Wrote outputs to ", output_dir)
  invisible(result)
}
