## ---------------------------------------------------------------------
##  The GECO pipeline: five SLURM-chained stages
##
##  A genome-wide run (~10k genes, ~50M gene pairs) cannot be held in
##  one R process, so the pipeline is split into stages that read and
##  write a single `work_dir`:
##
##    geco_prepare()             # one-time setup
##    geco_chunk()               # SLURM array (1..n_jobs), gene-pair estimation
##    geco_finalize_estimates()  # merge chunks, reformat, p, PSD
##    geco_precision_chunk()     # SLURM array (1..n_ct+1), glasso+eBIC
##    geco_finalize_outputs()    # degrees + final outputs (incl. degree.xlsx)
##
##  Cell-level data are stored with the {fst} package for
##  column-selectable on-disk access. {fst} is in Suggests; install it
##  before running the pipeline.
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
#' @param genes Character vector of gene IDs (required). Exactly these
#'   columns are selected from the .h5ad's expression matrix and every
#'   other gene is ignored; any requested gene that is not in the
#'   matrix is an error. Load from a text file with [read_gene_list()]
#'   if needed.
#' @param ind_col,ct_col Names of the individual-ID and cell-type
#'   columns in the .h5ad `obs` metadata.
#' @param subset_col,subset_value Optional metadata filter (see
#'   [load_h5ad()]).
#' @param fixed_covariates Per-individual covariates to regress out
#'   before GeCo estimation: either `NULL` (default, no regression) or
#'   a numeric matrix / `data.frame` with one row per individual and
#'   one column per covariate, row names = individual IDs as they
#'   appear in `ind_col`. Factors must already be encoded as numeric
#'   (e.g. 0/1) columns; do not include an intercept column (one is
#'   fitted per gene and cell type). The matrix is validated (numeric,
#'   no NA, every individual present, no constant column) and its rows
#'   are reordered to the data's individual order. See
#'   [compute_PC1()] if you want an expression PC as a covariate.
#' @param backend Which .h5ad reader to use; see [load_h5ad()].
#'
#' @details The `Y.rds` and `VS_diag.rds` tensors are always computed
#' here from the loaded cells (the same cells written to
#' `cell_data.fst`); they cannot be supplied from outside, because the
#' chunk workers rely on them being the group means / variances of
#' exactly those cells.
#'
#' @section Errors instead of NA:
#' Before writing anything, the loaded data are checked for a complete
#' individual x cell-type grid with at least 2 cells per group.
#' Violations stop with the affected individuals / cell types named,
#' so the failure surfaces here rather than in every chunk worker
#' hours later.
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
                         backend = c("auto", "anndata", "zellkonverter")) {

  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Cluster mode requires the 'fst' package. install.packages('fst').")
  }
  backend <- match.arg(backend)
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
                    backend = backend,
                    subset_col = subset_col, subset_value = subset_value)

  ## Basic bookkeeping (cheap; fails fast if ind_col/ct_col are wrong).
  if (!(ind_col %in% colnames(data)))
    stop("ind_col '", ind_col, "' not in the loaded data.")
  if (!(ct_col %in% colnames(data)))
    stop("ct_col '",  ct_col,  "' not in the loaded data.")

  .check_group_col(data, ind_col, "individual", "geco_prepare")
  .check_group_col(data, ct_col,  "cell-type",  "geco_prepare")
  individuals <- as.character(unique(data[[ind_col]]))
  cts_raw     <- as.character(unique(data[[ct_col]]))   # straight from h5ad
  if (length(cts_raw) < 2L) {
    stop("geco_prepare(): only ", length(cts_raw), " cell type in column '",
         ct_col, "'; GECO needs at least 2.", call. = FALSE)
  }
  ## Fail fast (before writing anything) if the individual x cell-type
  ## grid is incomplete or has singleton groups — otherwise every chunk
  ## worker would fail with the same error hours later.
  codes  <- .group_codes(data, ind_col, ct_col, individuals, cts_raw,
                         where = "geco_prepare")
  .check_group_sizes(.group_counts(codes, length(individuals), length(cts_raw)),
                     individuals, cts_raw, min_cells = 2L, where = "geco_prepare")

  ## Write cell_data.fst FIRST — worker jobs need it to compute cov12
  ## per pair, and having it on disk early lets us free the in-memory
  ## copy for the subsequent big computations if needed. We only write
  ## the columns workers actually use (ind_col, ct_col, gene columns);
  ## other metadata columns are dropped.
  worker_cols <- c(ind_col, ct_col, genes)
  keep <- intersect(worker_cols, colnames(data))
  cell_data <- as.data.frame(data)[, keep, drop = FALSE]
  message("[geco_prepare] Writing cell_data.fst (",
          nrow(cell_data), " cells x ", ncol(cell_data),
          " cols = ind_col + ct_col + ", length(genes), " genes)...")
  fst::write_fst(cell_data,
                 path = file.path(work_dir, "cell_data.fst"),
                 compress = 50)
  rm(cell_data); invisible(gc(verbose = FALSE))

  message("[geco_prepare] Computing OPms...")
  opms <- compute_OPms(data, genes, ind_col = ind_col)

  message("[geco_prepare] Validating fixed_covariates...")
  fixed_cov <- .resolve_fixed_cov(data, fixed_covariates, ind_col = ind_col)

  ## Precompute Y (per-(ind, gene, ct) means) and VS_diag (per-(ind, ct,
  ## gene) within-group variance). Doing both once here means every
  ## chunk worker only has to compute the off-diagonal cov12 per pair.
  message("[geco_prepare] Computing Y tensor [n_ind=", length(individuals),
          ", n_gene=", length(genes), ", n_ct=", length(cts_raw), "]...")
  Y_obj <- compute_Y(data, genes,
                     ind_col = ind_col, ct_col = ct_col,
                     individuals = individuals, cell_types = cts_raw)
  saveRDS(Y_obj$Y, file.path(work_dir, "Y.rds"))
  rm(Y_obj)

  message("[geco_prepare] Computing VS_diag tensor [n_ind=",
          length(individuals), ", n_ct=", length(cts_raw),
          ", n_gene=", length(genes), "]...")
  VS_diag <- compute_VS_diag(data, genes,
                             ind_col = ind_col, ct_col = ct_col,
                             individuals = individuals,
                             cell_types  = cts_raw)
  saveRDS(VS_diag, file.path(work_dir, "VS_diag.rds"))
  rm(VS_diag)

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
          "/{cell_data.fst, Y.rds, VS_diag.rds, prep.rds, gene_pairs.rds}. ",
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

  ## Sanity-check the tensors once per worker (cheap) so a corrupt or
  ## mismatched Y.rds / VS_diag.rds fails here with a clear message
  ## rather than inside GeCo() for every pair.
  if (!identical(as.integer(dim(Y_all)),
                 as.integer(c(n_ind, length(prep$genes), n_ct)))) {
    stop("geco_chunk(): Y.rds has dims [", paste(dim(Y_all), collapse = ", "),
         "] but prep.rds implies [", n_ind, ", ", length(prep$genes), ", ",
         n_ct, "].", call. = FALSE)
  }
  if (!identical(as.integer(dim(VS_diag)),
                 as.integer(c(n_ind, n_ct, length(prep$genes))))) {
    stop("geco_chunk(): VS_diag.rds has dims [",
         paste(dim(VS_diag), collapse = ", "), "] but prep.rds implies [",
         n_ind, ", ", n_ct, ", ", length(prep$genes), "].", call. = FALSE)
  }
  .check_finite(Y_all,   "Y.rds",       "geco_chunk")
  .check_finite(VS_diag, "VS_diag.rds", "geco_chunk")
  if (!.no_fixed_cov(prep$fixed_cov)) {
    .check_finite(prep$fixed_cov, "prep$fixed_cov", "geco_chunk")
  }

  ## Map gene name -> column index in Y_all / VS_diag (same order as prep$genes)
  gene_to_idx <- stats::setNames(seq_along(prep$genes), prep$genes)

  opms_dt <- data.table::as.data.table(prep$opms)
  op_bad <- opms_dt$g[!is.finite(opms_dt$OPsd) | opms_dt$OPsd == 0]
  if (length(op_bad) > 0L) {
    stop("geco_chunk(): ", length(op_bad), " gene(s) in prep$opms have zero ",
         "or NA OPsd and cannot be standardized: ", .preview(op_bad), ".",
         call. = FALSE)
  }
  data.table::setkey(opms_dt, g)

  ## Re-raise any per-pair error with the gene pair named.
  with_pair <- function(expr, g1, g2) {
    tryCatch(expr, error = function(e) {
      stop("geco_chunk(): gene pair (", g1, ", ", g2, ") failed: ",
           conditionMessage(e), call. = FALSE)
    })
  }

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

  ## Group every cell ONCE (integer codes + counts; enforces >= 2 cells
  ## per group). Per pair we then only need one rowsum() over the cells.
  meta_df <- fst::read_fst(cell_fst, columns = c(ind_col, ct_col))
  meta_df <- data.table::as.data.table(meta_df)
  setup <- .group_setup(meta_df, ind_col, ct_col, inds, cts,
                        where = "geco_chunk")
  rm(meta_df)

  ## Gene columns are read from the fst file on demand. Pairs in a
  ## chunk are contiguous in combn() order, so g1 stays the same for
  ## long runs: cache its column and only re-read g2 each pair.
  read_gene <- function(g) {
    x <- fst::read_fst(cell_fst, columns = g)[[1L]]
    if (anyNA(x)) {
      stop("geco_chunk(): gene '", g, "' has NA expression in ",
           sum(is.na(x)), " cell(s) in cell_data.fst.", call. = FALSE)
    }
    as.numeric(x)
  }
  g1_cached <- NA_character_; x1 <- NULL

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
    if (is.na(m1) || is.na(m2)) {
      stop("geco_chunk(): gene ", if (is.na(m1)) g1 else g2,
           " has no row in prep$opms.", call. = FALSE)
    }

    idx1 <- gene_to_idx[g1]
    idx2 <- gene_to_idx[g2]
    if (is.na(idx1) || is.na(idx2)) {
      stop("geco_chunk(): gene ", if (is.na(idx1)) g1 else g2,
           " from gene_pairs.rds is not in prep$genes.", call. = FALSE)
    }

    ## Pull pre-computed Y means for these two genes (raw = per-group
    ## means of the fst data, used for the cov12 centering), then
    ## rescale by OP for the estimator.
    Y1raw <- Y_all[, idx1, , drop = FALSE]; dim(Y1raw) <- c(n_ind, n_ct)
    Y2raw <- Y_all[, idx2, , drop = FALSE]; dim(Y2raw) <- c(n_ind, n_ct)
    Y1 <- (Y1raw - m1) / s1
    Y2 <- (Y2raw - m2) / s2
    Y_pair <- array(NA_real_, dim = c(n_ind, 2, n_ct))
    Y_pair[, 1, ] <- Y1; Y_pair[, 2, ] <- Y2

    ## Off-diagonal cov12: read the gene column(s) and do one rowsum().
    if (!identical(g1, g1_cached)) {
      x1 <- with_pair(read_gene(g1), g1, g2)
      g1_cached <- g1
    }
    x2 <- with_pair(read_gene(g2), g1, g2)
    cov12 <- .group_cov_pair(x1, x2, Y1raw, Y2raw, setup)

    VS <- array(NA_real_, dim = c(n_ind, n_ct, 2, 2))
    VS[, , 1, 1] <- VS_diag[, , idx1] / s1^2
    VS[, , 2, 2] <- VS_diag[, , idx2] / s2^2
    VS[, , 1, 2] <- cov12 / s1 / s2
    VS[, , 2, 1] <- VS[, , 1, 2]

    res <- with_pair(GeCo(Y = Y_pair, VS = VS, fixed_cov = prep$fixed_cov),
                     g1, g2)

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


#' Merge chunks and write per-cell-type cov/se/p matrices
#'
#' Concatenates all `work_dir/chunks/chunk_*.rds` files, calls
#' [reformat_to_matrices()], computes p-value matrices, and writes
#' per-cell-type `cov_<ct>.rds`, `se_<ct>.rds`, and `p_<ct>.rds`
#' files directly into `output_dir`. Also saves a small
#' `work_dir/cell_type_order.rds` (a character vector) so step 03
#' knows the canonical ordering of cell types — no consolidated
#' `mats.rds` file is written.
#'
#' @param work_dir Working directory built by [geco_prepare()].
#' @param output_dir Directory to write the final per-cell-type
#'   matrices into. Created if missing.
#' @param cleanup_chunks Logical. If `TRUE` (the default), the
#'   `work_dir/chunks/chunk_*.rds` files are deleted after the merged
#'   `cov_*/se_*/p_*` matrices have been written to `output_dir` —
#'   they're not needed by any later step. Set to `FALSE` to keep them
#'   for inspection / re-runs.
#'
#' @return Invisibly the merged matrices list (`cov`, `se`, `p`).
#' @export
geco_finalize_estimates <- function(work_dir, output_dir,
                                    cleanup_chunks = TRUE) {

  prep <- readRDS(file.path(work_dir, "prep.rds"))
  chunk_files <- list.files(file.path(work_dir, "chunks"),
                            pattern = "^chunk_.*\\.rds$",
                            full.names = TRUE)
  if (length(chunk_files) == 0L) stop("No chunk files found in ", work_dir)
  message("[geco_finalize_estimates] Merging ", length(chunk_files), " chunks...")

  cov_parts <- vector("list", length(chunk_files))
  se_parts  <- vector("list", length(chunk_files))
  task_ids  <- integer(length(chunk_files))
  starts    <- integer(length(chunk_files))
  ends      <- integer(length(chunk_files))
  for (i in seq_along(chunk_files)) {
    chunk <- readRDS(chunk_files[i])
    cov_parts[[i]] <- chunk$cov
    se_parts [[i]] <- chunk$se
    task_ids[i] <- chunk$task_id
    starts[i]   <- chunk$start
    ends[i]     <- chunk$end
  }

  ## ---- every gene pair must be covered exactly once -------------------
  ## A SLURM task that was killed (OOM / time limit) leaves no file, and
  ## the merged table would silently lack its pairs. Fail here instead.
  if (anyDuplicated(task_ids)) {
    stop("geco_finalize_estimates(): duplicated task_id in chunk files: ",
         .preview(unique(task_ids[duplicated(task_ids)])),
         ". Delete the stale chunk_*.rds and re-run.", call. = FALSE)
  }
  covered <- logical(prep$n_pairs)
  for (i in seq_along(starts)) covered[starts[i]:ends[i]] <- TRUE
  if (!all(covered)) {
    miss <- which(!covered)
    ## collapse to contiguous ranges for readability
    brk <- c(0L, which(diff(miss) != 1L), length(miss))
    ranges <- vapply(seq_len(length(brk) - 1L), function(k) {
      a <- miss[brk[k] + 1L]; b <- miss[brk[k + 1L]]
      if (a == b) as.character(a) else paste0(a, "-", b)
    }, character(1))
    n_jobs_guess <- if (length(starts)) max(ends - starts + 1L) else NA
    stop("geco_finalize_estimates(): ", length(miss), " of ", prep$n_pairs,
         " gene pairs are not covered by any chunk_*.rds file (pair index ",
         "ranges: ", .preview(ranges), "). Found ", length(chunk_files),
         " chunk file(s) with task_ids ", .preview(sort(task_ids)),
         if (!is.na(n_jobs_guess)) paste0("; chunk size suggests n_jobs ~ ",
           ceiling(prep$n_pairs / n_jobs_guess)) else "",
         ". Some SLURM array tasks did not finish; check their .out files ",
         "and resubmit the missing task_ids (see README, 'Resuming after a ",
         "partial failure').", call. = FALSE)
  }

  cov_df <- do.call(rbind, cov_parts)
  se_df  <- do.call(rbind, se_parts)

  est <- list(cov = cov_df, se = se_df,
              cell_types = prep$cell_types_raw)
  mats <- reformat_to_matrices(est, prep$genes,
                               cell_types = prep$cell_type_names)

  message("[geco_finalize_estimates] Computing p-value matrices...")
  p_mats <- mapply(compute_p_matrix, mats$cov, mats$se, SIMPLIFY = FALSE)
  names(p_mats) <- names(mats$cov)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(mats$cov)) {
    saveRDS(mats$cov[[nm]], file.path(output_dir, paste0("cov_", nm, ".rds")))
    saveRDS(mats$se [[nm]], file.path(output_dir, paste0("se_",  nm, ".rds")))
    saveRDS(p_mats  [[nm]], file.path(output_dir, paste0("p_",   nm, ".rds")))
  }

  ## Tiny metadata file used by step 03 to size the precision array.
  saveRDS(names(mats$cov), file.path(work_dir, "cell_type_order.rds"))

  message("[geco_finalize_estimates] Wrote per-cell-type ",
          "cov_*.rds / se_*.rds / p_*.rds to ", output_dir, ".")

  ## chunk_*.rds files have served their purpose — delete to save disk
  ## unless the caller explicitly asked to keep them.
  if (isTRUE(cleanup_chunks)) {
    removed <- file.remove(chunk_files)
    n_ok <- sum(removed)
    message("[geco_finalize_estimates] Cleaned up ", n_ok,
            " chunk_*.rds file(s) from ", file.path(work_dir, "chunks"),
            ". Pass cleanup_chunks = FALSE to keep them.")
  }

  invisible(list(cov = mats$cov, se = mats$se, p = p_mats))
}


#' Fit the precision matrix for a single cell type (SLURM array element)
#'
#' Reads `output_dir/cov_<ct>.rds` and `output_dir/p_<ct>.rds`,
#' PSD-clips the covariance on the fly with [transform_psd()], runs
#' [fit_precision()], and writes the result to
#' `output_dir/precision_<ct>.rds`. The full fit (including lambda
#' path) is also saved to `work_dir/precision/fit_<ct>.rds` for
#' diagnostics.
#'
#' @param work_dir Working directory built by [geco_prepare()].
#' @param output_dir Directory with `cov_<ct>.rds` and `p_<ct>.rds`
#'   written by [geco_finalize_estimates()]. The precision matrix is
#'   written here too.
#' @param ct_idx   Integer index into the ordered list
#'   `c("shared", <cell types from h5ad>)`. So valid values are
#'   `1..(n_cell_types + 1)`.
#' @param psd_thr Threshold for [transform_psd()] (default 0).
#' @param fit_args Named list of arguments forwarded to
#'   [fit_precision()].
#'
#' @return Invisibly the precision fit result.
#' @export
geco_precision_chunk <- function(work_dir, output_dir, ct_idx,
                                 psd_thr = 0, fit_args = list()) {
  prep   <- readRDS(file.path(work_dir, "prep.rds"))
  order  <- readRDS(file.path(work_dir, "cell_type_order.rds"))
  if (ct_idx < 1 || ct_idx > length(order)) {
    stop("ct_idx must be in 1..", length(order))
  }
  nm <- order[ct_idx]
  message("[geco_precision_chunk] Fitting precision for: ", nm)

  cov <- readRDS(file.path(output_dir, paste0("cov_", nm, ".rds")))
  p   <- readRDS(file.path(output_dir, paste0("p_",   nm, ".rds")))
  cov_psd <- transform_psd(cov, thr = psd_thr)

  n_ind <- length(prep$individuals)
  fit <- do.call(fit_precision,
    c(list(cov = cov_psd, p = p, n_ind = n_ind), fit_args))
  dimnames(fit$precision) <- list(prep$genes, prep$genes)

  saveRDS(fit$precision,
          file.path(output_dir, paste0("precision_", nm, ".rds")))
  saveRDS(list(lambda = fit$lambda, eBIC_path = fit$eBIC_path),
          file.path(work_dir, "precision", paste0("fit_", nm, ".rds")))

  message("[geco_precision_chunk] Wrote ",
          file.path(output_dir, paste0("precision_", nm, ".rds")),
          " (lambda = ", fit$lambda, ")")
  invisible(fit)
}


#' Finalize cluster outputs: build degree.xlsx
#'
#' Reads `output_dir/precision_*.rds` files (written by
#' [geco_precision_chunk()]), computes node degrees, and writes
#' `output_dir/degree.xlsx` with one column per cell type. All other
#' per-cell-type outputs (`cov_*`, `se_*`, `p_*`, `precision_*`) are
#' already in `output_dir` by the time this runs — this function only
#' adds the degree workbook.
#'
#' @param work_dir Working directory built by [geco_prepare()].
#' @param output_dir Final outputs directory.
#'
#' @return Invisibly the degree `data.frame`.
#' @export
geco_finalize_outputs <- function(work_dir, output_dir) {
  prep  <- readRDS(file.path(work_dir, "prep.rds"))
  order <- readRDS(file.path(work_dir, "cell_type_order.rds"))

  genes <- prep$genes
  pfiles <- file.path(output_dir, paste0("precision_", order, ".rds"))
  missing <- !file.exists(pfiles)
  if (any(missing)) {
    stop("geco_finalize_outputs(): precision file(s) missing for ",
         sum(missing), " of ", length(order), " cell type(s): ",
         .preview(order[missing]), " (expected ",
         .preview(basename(pfiles[missing])), " in ", output_dir, "). ",
         "The corresponding step-03 SLURM tasks (ct_idx = ",
         .preview(which(missing)), ") did not finish; check their .out ",
         "files and resubmit them before running this step.", call. = FALSE)
  }

  deg_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (k in seq_along(order)) {
    nm <- order[k]
    Omega <- readRDS(pfiles[k])
    .check_finite(Omega, paste0("precision_", nm, ".rds"), "geco_finalize_outputs")
    deg   <- compute_degree(Omega)
    if (is.null(names(deg)) || !all(genes %in% names(deg))) {
      stop("geco_finalize_outputs(): precision_", nm, ".rds has no (or ",
           "incomplete) gene row names; cannot assign degrees to genes.",
           call. = FALSE)
    }
    deg_df[[nm]] <- as.integer(deg[genes])
  }

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "degree")
    openxlsx::writeData(wb, "degree", deg_df)
    openxlsx::saveWorkbook(wb,
      file = file.path(output_dir, "degree.xlsx"),
      overwrite = TRUE)
  } else {
    utils::write.csv(deg_df,
      file = file.path(output_dir, "degree.csv"),
      row.names = FALSE)
    message("openxlsx not available; wrote degree.csv instead of .xlsx")
  }

  message("[geco_finalize_outputs] Wrote ",
          file.path(output_dir, "degree.xlsx"))
  invisible(deg_df)
}
