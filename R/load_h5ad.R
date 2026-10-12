## Internal: message used when no .h5ad backend is available. The
## install itself is handled (with consent) by check_dependencies().
.h5ad_install_msg <- function(which = "auto") {
  paste0(
    switch(which,
      auto          = "No .h5ad reader is installed (need the 'anndata' R package, or zellkonverter + SingleCellExperiment + SummarizedExperiment).",
      anndata       = "backend = 'anndata' requested but the 'anndata' R package isn't installed.",
      zellkonverter = "backend = 'zellkonverter' requested but zellkonverter / SingleCellExperiment / SummarizedExperiment aren't all installed."),
    "\n  Run GECO::check_dependencies(backend = \"", which,
    "\") to install it (you will be asked for confirmation)."
  )
}


#' Load a .h5ad file into a long-format expression data frame
#'
#' Reads the supplied .h5ad file and returns a single `data.table`
#' whose columns are the cell-level metadata followed by the gene
#' expression columns.
#'
#' Two backends are supported (in this order of preference):
#'
#' 1. **anndata** (CRAN). Uses `reticulate` + whatever Python is on your
#'   `PATH` (or configured via `reticulate::use_python()` /
#'   `RETICULATE_PYTHON`). This is the recommended backend on HPC
#'   nodes where you already have a Python environment with the
#'   `anndata` Python package.
#'   Install: `install.packages("anndata")` (R) and
#'   `pip install anndata` (Python).
#' 2. **zellkonverter** (Bioconductor). Bundles a private Python
#'   environment via `basilisk`, which builds Python from source on
#'   first use. This works cleanly on most workstations but often
#'   fails on locked-down HPC nodes missing Tk / OpenSSL / build
#'   headers. Prefer the anndata backend there.
#'
#' @param h5ad_file Path to a `.h5ad` file. Assumed to have already
#'   passed the user's QC filtering.
#' @param genes Character vector of gene IDs to extract from the
#'   expression matrix. If `NULL`, every gene in the matrix is kept.
#'   Either pass an in-memory vector or load one from a file with
#'   [read_gene_list()] first. Any requested gene that's not in the
#'   matrix raises an error.
#' @param backend One of `"auto"` (default: try anndata first, then
#'   zellkonverter), `"anndata"`, or `"zellkonverter"`. Force one
#'   backend if the automatic choice picks the wrong one.
#' @param assay Name of the assay layer to extract. Defaults to "X".
#'   (Only used by the zellkonverter backend.)
#' @param gene_id_col Column in `rowData(adata)` that holds gene IDs to
#'   use as gene column names. If `NULL`, `rownames(expr_matrix)` are
#'   kept as-is. (Only used by the zellkonverter backend.)
#' @param keep_meta Optional character vector restricting which metadata
#'   columns are returned. By default, all are kept.
#' @param subset_col,subset_value Optional: filter cells where
#'   `metadata[[subset_col]] == subset_value`. NULL means no filtering.
#'
#' @return A `data.table` with one row per cell. Columns: all (kept)
#'   metadata columns first, then one column per requested gene
#'   (numeric).
#' @export
load_h5ad <- function(h5ad_file,
                      genes = NULL,
                      backend = c("auto", "anndata", "zellkonverter"),
                      assay = "X",
                      gene_id_col = "gene_ids",
                      keep_meta = NULL,
                      subset_col = NULL,
                      subset_value = NULL) {

  if (!file.exists(h5ad_file)) {
    stop("h5ad file not found: ", h5ad_file)
  }
  backend <- match.arg(backend)

  expr_matrix <- NULL
  metadata    <- NULL

  has_anndata <- requireNamespace("anndata", quietly = TRUE)
  has_zk <- requireNamespace("zellkonverter",       quietly = TRUE) &&
            requireNamespace("SingleCellExperiment", quietly = TRUE) &&
            requireNamespace("SummarizedExperiment", quietly = TRUE)

  ## Decide which backend to actually use.
  use <- switch(backend,
    auto          = if (has_anndata) "anndata"
                    else if (has_zk) "zellkonverter"
                    else stop(.h5ad_install_msg("auto"), call. = FALSE),
    anndata       = if (has_anndata) "anndata"
                    else stop(.h5ad_install_msg("anndata"), call. = FALSE),
    zellkonverter = if (has_zk) "zellkonverter"
                    else stop(.h5ad_install_msg("zellkonverter"), call. = FALSE)
  )
  message("[load_h5ad] Using backend: ", use)

  if (use == "anndata") {
    adata <- anndata::read_h5ad(h5ad_file)

    ## var_names / obs_names come through as pandas.Index proxies:
    ## as.character() on them returns the printed repr, not values.
    ## Use $tolist() to force a real Python list -> R char vector.
    as_names <- function(x) {
      if (inherits(x, "python.builtin.object")) {
        as.character(unlist(x$tolist()))
      } else {
        as.character(x)
      }
    }
    var_names <- as_names(adata$var_names)
    obs_names <- as_names(adata$obs_names)

    ## Densify sparse X if needed. Dense X comes through as an R
    ## matrix already in (n_cells x n_genes) orientation.
    X <- adata$X
    if (inherits(X, "python.builtin.object") &&
        any(grepl("scipy\\.sparse", class(X)))) {
      X <- X$toarray()
    }
    X <- as.matrix(X)
    if (nrow(X) != length(obs_names) || ncol(X) != length(var_names)) {
      stop("Unexpected X shape: got ", nrow(X), "x", ncol(X),
           " but obs/var lengths are ", length(obs_names), " / ",
           length(var_names), ".")
    }
    expr_matrix <- t(X)                  # -> genes x cells
    rownames(expr_matrix) <- var_names
    colnames(expr_matrix) <- obs_names

    ## adata$obs is usually already an R data.frame; coerce defensively.
    metadata <- adata$obs
    if (!is.data.frame(metadata)) metadata <- as.data.frame(metadata)
    rownames(metadata) <- obs_names

  } else {
    adata <- zellkonverter::readH5AD(h5ad_file)
    expr_matrix <- SummarizedExperiment::assay(adata, assay)
    metadata    <- as.data.frame(SummarizedExperiment::colData(adata))

    if (!is.null(gene_id_col)) {
      rd <- SummarizedExperiment::rowData(adata)
      if (gene_id_col %in% colnames(rd)) {
        rownames(expr_matrix) <- rd[[gene_id_col]]
      }
    }
  }

  # Optional gene subsetting (done while still genes x cells, so the
  # transpose is cheap).
  if (!is.null(genes)) {
    avail <- rownames(expr_matrix)
    if (is.null(avail)) {
      stop("Cannot select genes by name: expression matrix has no rownames. ",
           "Check the gene_id_col argument.")
    }
    miss <- setdiff(genes, avail)
    if (length(miss) > 0L) {
      n <- length(miss)
      preview <- paste(utils::head(miss, 5), collapse = ", ")
      stop(n, " requested gene", if (n > 1L) "s are" else " is",
           " not in the h5ad: ", preview,
           if (n > 5L) sprintf(" (... and %d more)", n - 5L) else "",
           call. = FALSE)
    }
    expr_matrix <- expr_matrix[genes, , drop = FALSE]
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


#' Read a gene list from a text file
#'
#' Tiny convenience wrapper: one gene ID per line (no header). Useful
#' when you already have something like `genelist11087.txt` from the
#' legacy pipeline.
#'
#' @param path Path to the file.
#' @param comment.char Character marking comment lines. Default `"#"`.
#'
#' @return Character vector of gene IDs (whitespace trimmed, blanks
#'   dropped, duplicates dropped while preserving order).
#'
#' @examples
#' \dontrun{
#' genes <- read_gene_list("genelist11087.txt")
#' length(genes)
#' }
#' @export
read_gene_list <- function(path, comment.char = "#") {
  if (!file.exists(path)) stop("Gene list file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  lines <- sub(paste0("\\s*", comment.char, ".*$"), "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  unique(lines)
}


#' Read a per-individual covariate matrix from a file
#'
#' Loads the `fixed_covariates` matrix for [geco_prepare()] from disk.
#' Used by the `00_prepare.R` script that [generate_slurm_pipeline()]
#' writes, and handy interactively.
#'
#' @param path Path to the covariate table. The default format is a
#'   plain-text `.txt` file: tab- or whitespace-separated, a header row
#'   with the covariate names, and the individual ID in the first
#'   column (one row per individual), e.g. as written by
#'   `write.table(m, "covariates.txt", sep = "\t", quote = FALSE,
#'   col.names = NA)`. Also accepted, detected from the extension:
#'   `.tsv` (tab), `.csv` (comma), and `.rds` (a saved numeric matrix
#'   or `data.frame` with row names).
#'
#' @return A numeric matrix with row names = individual IDs. The
#'   content is validated later by [geco_prepare()] (numeric, no NA,
#'   every individual present, no constant column), so this function
#'   only checks that the file exists and has row names.
#' @export
read_covariates <- function(path) {
  if (!file.exists(path)) stop("Covariate file not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    m <- readRDS(path)
  } else {
    sep <- switch(ext, csv = ",", tsv = "\t", "")   # txt / anything else: whitespace
    m <- utils::read.table(path, header = TRUE, sep = sep, row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE,
                           quote = "\"", comment.char = "")
  }
  if (!(is.matrix(m) || is.data.frame(m))) {
    stop("read_covariates(): ", path, " does not contain a matrix or ",
         "data.frame (got class '", paste(class(m), collapse = "/"), "').",
         call. = FALSE)
  }
  if (is.null(rownames(m)) ||
      identical(rownames(m), as.character(seq_len(nrow(m))))) {
    stop("read_covariates(): ", path, " has no row names. The first column ",
         "must be the individual IDs (for .rds: the matrix row names).",
         call. = FALSE)
  }
  m <- as.matrix(m)
  if (!is.numeric(m)) {
    stop("read_covariates(): ", path, " has non-numeric columns. Encode ",
         "factors as 0/1 columns before saving.", call. = FALSE)
  }
  m
}
