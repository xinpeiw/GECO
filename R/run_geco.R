#' Run the full GECO pipeline end-to-end
#'
#' Convenience wrapper around all of the package's exported functions.
#' Takes a QC'd `.h5ad` file and a vector of gene names, performs the
#' per-pair GeCo estimation, reformats into matrices, computes p-values,
#' fits sparse precision matrices via glassoFast + eBIC for every cell
#' type (plus shared), computes degrees, and optionally writes everything
#' to `output_dir` including a single `degree.xlsx` with one sheet per
#' cell type.
#'
#' @param h5ad_file Path to the `.h5ad` file.
#' @param genes Character vector of gene IDs (required). The package
#'   selects exactly these columns from the .h5ad's expression matrix
#'   and ignores every other column. Load from a text file with
#'   [read_gene_list()] if needed.
#' @param ind_col,ct_col Names of the individual and cell-type columns.
#' @param subset_col,subset_value Optional filter on a metadata column
#'   (see [load_h5ad()]).
#' @param fixed_covariates Which covariates to regress out before
#'   computing GeCo. Accepts:
#'   - `NULL` (default): no regression.
#'   - A character vector of names from [GECO_COVARIATE_CATALOG],
#'     e.g. `c("age_group", "sex", "PC1")`. Call [available_covariates()]
#'     to see all options.
#'   - A named list of fully-specified covariate recipes (see
#'     [build_fixed_cov()] for the spec format).
#'   - A precomputed matrix / `data.frame` with one row per individual,
#'     used directly.
#' @param output_dir Directory to save all results. If `NULL` nothing
#'   is written and everything is returned.
#' @param psd_thr Threshold for [transform_psd()]. Default 0.
#' @param fit_precision_args Named list of extra arguments forwarded to
#'   [fit_precision()].
#' @param verbose Pair-level progress interval for [estimate_geco()].
#'
#' @return A list with elements `cov`, `se`, `p`, `cov_psd`, `precision`
#'   (each a named list keyed by cell type, including `shared`),
#'   `degree` (a `data.frame` of degrees with one column per cell type),
#'   and `genes`.
#' @export
run_geco <- function(h5ad_file,
                     genes,
                     ind_col = "ind_cov",
                     ct_col  = "cg_cov",
                     subset_col = NULL, subset_value = NULL,
                     fixed_covariates = NULL,
                     output_dir = NULL,
                     psd_thr = 0,
                     fit_precision_args = list(),
                     verbose = 0) {

  if (missing(genes) || is.null(genes)) {
    stop("`genes` is required. Pass a character vector of gene IDs ",
         "(or use read_gene_list() to load one from a .txt file).",
         call. = FALSE)
  }
  if (length(genes) < 2L) stop("Need at least 2 genes.")
  genes <- unique(as.character(genes))

  message("[GECO] Loading h5ad: ", h5ad_file)
  data <- load_h5ad(h5ad_file,
                    genes = genes,
                    subset_col = subset_col,
                    subset_value = subset_value)
  message("[GECO] Using ", length(genes), " genes across ",
          length(unique(data[[ind_col]])), " individuals and ",
          length(unique(data[[ct_col]])), " cell types.")

  message("[GECO] Computing per-gene OPmean / OPsd...")
  opms <- compute_OPms(data, genes, ind_col = ind_col)

  fixed_cov <- resolve_fixed_covariates(data, genes,
                                        fixed_covariates = fixed_covariates,
                                        ind_col = ind_col)

  message("[GECO] Estimating GeCo across all gene pairs...")
  est <- estimate_geco(data, genes, opms = opms,
                       fixed_cov = fixed_cov,
                       ind_col = ind_col, ct_col = ct_col,
                       verbose = verbose)

  cell_type_names <- est$cell_types   # taken straight from the h5ad

  message("[GECO] Reformatting to per-cell-type matrices...")
  mats <- reformat_to_matrices(est, genes, cell_types = cell_type_names)

  message("[GECO] Computing p-value matrices...")
  p_mats <- mapply(compute_p_matrix, mats$cov, mats$se, SIMPLIFY = FALSE)
  names(p_mats) <- names(mats$cov)

  message("[GECO] Applying PSD transform to covariance matrices...")
  cov_psd <- lapply(mats$cov, transform_psd, thr = psd_thr)

  message("[GECO] Fitting sparse precision matrices...")
  n_ind <- length(unique(data[[ind_col]]))
  precisions <- list()
  lambdas    <- list()
  ebic_paths <- list()
  for (nm in names(cov_psd)) {
    message("    - ", nm)
    res <- do.call(fit_precision,
      c(list(cov = cov_psd[[nm]], p = p_mats[[nm]], n_ind = n_ind),
        fit_precision_args))
    precisions[[nm]] <- res$precision
    lambdas[[nm]]    <- res$lambda
    ebic_paths[[nm]] <- res$eBIC_path
    # ensure dimnames
    dimnames(precisions[[nm]]) <- list(genes, genes)
  }

  message("[GECO] Computing degrees...")
  deg_list <- lapply(precisions, compute_degree)
  deg_df   <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (nm in names(deg_list)) {
    deg_df[[nm]] <- as.integer(deg_list[[nm]][genes])
  }

  result <- list(
    cov        = mats$cov,
    se         = mats$se,
    p          = p_mats,
    cov_psd    = cov_psd,
    precision  = precisions,
    lambda     = lambdas,
    eBIC_paths = ebic_paths,
    degree     = deg_df,
    genes      = genes,
    cell_types = cell_type_names
  )

  if (!is.null(output_dir)) {
    save_results(result, output_dir)
  }

  result
}


#' Save GECO results to a directory
#'
#' Writes per-cell-type covariance/SE/p/precision matrices as .rds
#' files plus a single `degree.xlsx` workbook.
#'
#' @param result Result list returned by [run_geco()].
#' @param output_dir Output directory (will be created if missing).
#'
#' @return Invisibly returns `output_dir`.
#' @export
save_results <- function(result, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(result$cov)) {
    saveRDS(result$cov      [[nm]], file.path(output_dir, paste0("cov_",       nm, ".rds")))
    saveRDS(result$se       [[nm]], file.path(output_dir, paste0("se_",        nm, ".rds")))
    saveRDS(result$p        [[nm]], file.path(output_dir, paste0("p_",         nm, ".rds")))
    saveRDS(result$cov_psd  [[nm]], file.path(output_dir, paste0("cov_psd_",   nm, ".rds")))
    saveRDS(result$precision[[nm]], file.path(output_dir, paste0("precision_", nm, ".rds")))
  }

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "degree")
    openxlsx::writeData(wb, "degree", result$degree)
    openxlsx::saveWorkbook(wb,
      file = file.path(output_dir, "degree.xlsx"),
      overwrite = TRUE)
  } else {
    utils::write.csv(result$degree,
      file = file.path(output_dir, "degree.csv"),
      row.names = FALSE)
    message("openxlsx not available; wrote degree.csv instead of .xlsx")
  }
  invisible(output_dir)
}
