#' Load a .h5ad file into a long-format expression data frame
#'
#' Reads the supplied .h5ad file (via Bioconductor zellkonverter, or the
#' anndata package if available) and returns a single `data.table` whose
#' columns are the cell-level metadata followed by all gene expression
#' columns. The output is the input expected by the rest of the GECO
#' pipeline.
#'
#' @param h5ad_file Path to a `.h5ad` file. Assumed to have already passed
#'   the user's QC filtering.
#' @param assay Name of the assay layer to extract. Defaults to "X".
#' @param gene_id_col Column in `rowData(adata)` that holds gene IDs to
#'   use as gene column names. If `NULL`, `rownames(expr_matrix)` are
#'   kept as-is.
#' @param keep_meta Optional character vector restricting which metadata
#'   columns are returned. By default, all are kept.
#' @param subset_col,subset_value Optional: filter cells where
#'   `metadata[[subset_col]] == subset_value` (e.g. `"SLE_status"` and
#'   `"Healthy"`). NULL means no filtering.
#'
#' @return A `data.table` with one row per cell. Columns: all (kept)
#'   metadata columns first, then one column per gene (numeric).
#' @export
load_h5ad <- function(h5ad_file,
                      assay = "X",
                      gene_id_col = "gene_ids",
                      keep_meta = NULL,
                      subset_col = NULL,
                      subset_value = NULL) {

  if (!file.exists(h5ad_file)) {
    stop("h5ad file not found: ", h5ad_file)
  }

  expr_matrix <- NULL
  metadata    <- NULL

  # Preferred: zellkonverter (Bioconductor)
  if (requireNamespace("zellkonverter", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {

    adata <- zellkonverter::readH5AD(h5ad_file)
    expr_matrix <- SummarizedExperiment::assay(adata, assay)
    metadata    <- as.data.frame(SummarizedExperiment::colData(adata))

    if (!is.null(gene_id_col)) {
      rd <- SummarizedExperiment::rowData(adata)
      if (gene_id_col %in% colnames(rd)) {
        rownames(expr_matrix) <- rd[[gene_id_col]]
      }
    }

  } else if (requireNamespace("anndata", quietly = TRUE)) {
    # Fallback: anndata reticulate
    adata <- anndata::read_h5ad(h5ad_file)
    expr_matrix <- t(as.matrix(adata$X))   # genes x cells
    rownames(expr_matrix) <- adata$var_names
    colnames(expr_matrix) <- adata$obs_names
    metadata <- adata$obs

  } else {
    stop("Need either 'zellkonverter' (Bioconductor) or 'anndata' (CRAN) ",
         "installed to read .h5ad files. Install with: ",
         'BiocManager::install("zellkonverter") or install.packages("anndata")')
  }

  # genes x cells -> cells x genes
  expr_df <- t(as.matrix(expr_matrix))

  # align cell order to metadata
  if (!is.null(rownames(metadata)) && !is.null(rownames(expr_df))) {
    common <- intersect(rownames(metadata), rownames(expr_df))
    if (length(common) > 0) {
      expr_df  <- expr_df[common, , drop = FALSE]
      metadata <- metadata[common, , drop = FALSE]
    }
  }

  if (!is.null(keep_meta)) {
    keep_meta <- intersect(keep_meta, colnames(metadata))
    metadata  <- metadata[, keep_meta, drop = FALSE]
  }

  data_all <- cbind(metadata, expr_df)
  data_all <- data.table::as.data.table(data_all)

  if (!is.null(subset_col) && !is.null(subset_value)) {
    if (!(subset_col %in% colnames(data_all))) {
      stop("subset_col '", subset_col, "' not in metadata columns.")
    }
    data_all <- data_all[data_all[[subset_col]] == subset_value, ]
  }

  data_all
}
