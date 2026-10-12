## ---------------------------------------------------------------------
##  Pure-R .h5ad reader (rhdf5), no Python.
##
##  An .h5ad file is HDF5. We read only what the pipeline needs:
##    * the requested obs columns (individual, cell type),
##    * the requested genes from X / a layer / raw.X,
##  without ever materialising the full expression matrix.
##
##  On-disk layouts handled (anndata >= 0.7):
##    dataframe  : group with attrs `_index`, `column-order`; one child
##                 per column; or (very old) a compound dataset
##    column     : dataset (numeric / string / bool-enum)
##                 group "categorical"  = codes (+ -1 = NA) & categories
##                 group "nullable-*"   = values & mask
##                 old 0.7 categorical  = int codes dataset whose
##                                        categories live in __categories/<col>
##    matrix     : dense dataset, or group "csr_matrix"/"csc_matrix"
##                 (data, indices, indptr; attr `shape`), or the older
##                 h5sparse_format / h5sparse_shape attrs
## ---------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

## Listing of the file as a data.frame with a full `path` column.
.h5_ls <- function(file) {
  ls <- rhdf5::h5ls(file, all = FALSE)
  ls$path <- paste0(sub("/$", "", ls$group), "/", ls$name)
  ls
}
.h5_entry <- function(ls, path) ls[ls$path == path, , drop = FALSE]
.h5_exists <- function(ls, path) any(ls$path == path)
.h5_is_group <- function(ls, path) {
  e <- .h5_entry(ls, path); nrow(e) == 1L && e$otype == "H5I_GROUP"
}
.h5_attrs <- function(file, path) {
  tryCatch(rhdf5::h5readAttributes(file, path), error = function(e) list())
}
.h5_children <- function(ls, path) ls$name[ls$group == path]

## Read a 1-D dataset as a plain R vector (int64 -> double, bool -> logical).
.h5_vec <- function(file, path, ...) {
  x <- rhdf5::h5read(file, path, bit64conversion = "double", ...)
  if (is.factor(x)) {                          # HDF5 enum (numpy bool)
    lv <- levels(x)
    if (all(lv %in% c("FALSE", "TRUE"))) return(as.character(x) == "TRUE")
    return(as.character(x))
  }
  as.vector(x)
}

## Read one column of an anndata dataframe (obs or var) as a vector.
.h5ad_read_col <- function(file, df_path, col, ls) {
  path <- paste0(df_path, "/", col)
  if (!.h5_exists(ls, path)) {
    stop("load_h5ad(): column '", col, "' not found in ", df_path, ". ",
         "Available: ", .preview(setdiff(.h5_children(ls, df_path),
                                         "__categories"), 20), ".",
         call. = FALSE)
  }
  if (.h5_is_group(ls, path)) {
    enc <- .h5_attrs(file, path)[["encoding-type"]] %||% ""
    kids <- .h5_children(ls, path)
    if (enc == "categorical" || all(c("codes", "categories") %in% kids)) {
      codes <- .h5_vec(file, paste0(path, "/codes"))
      cats  <- .h5_vec(file, paste0(path, "/categories"))
      codes[codes < 0] <- NA
      return(cats[codes + 1L])
    }
    if (grepl("^nullable", enc) || all(c("values", "mask") %in% kids)) {
      v <- .h5_vec(file, paste0(path, "/values"))
      m <- .h5_vec(file, paste0(path, "/mask"))
      v[as.logical(m)] <- NA
      return(v)
    }
    stop("load_h5ad(): column '", col, "' in ", df_path, " has an ",
         "unsupported encoding ('", enc, "'; children: ",
         paste(kids, collapse = ", "), ").", call. = FALSE)
  }
  x <- .h5_vec(file, path)
  ## anndata 0.7 categorical: integer codes + categories under __categories
  cat_path <- paste0(df_path, "/__categories/", col)
  if (!is.null(.h5_attrs(file, path)[["categories"]]) && .h5_exists(ls, cat_path)) {
    cats <- .h5_vec(file, cat_path)
    x[x < 0] <- NA
    x <- cats[x + 1L]
  }
  x
}

## Index (row names) of an anndata dataframe.
.h5ad_read_index <- function(file, df_path, ls) {
  e <- .h5_entry(ls, df_path)
  if (nrow(e) == 0L) stop("load_h5ad(): no '", df_path, "' in the file.", call. = FALSE)
  if (e$otype != "H5I_GROUP") {                 # ancient compound layout
    df <- rhdf5::h5read(file, df_path)
    idx <- .h5_attrs(file, df_path)[["_index"]] %||% colnames(df)[1L]
    return(as.character(df[[idx]]))
  }
  idx <- .h5_attrs(file, df_path)[["_index"]] %||% "_index"
  if (!.h5_exists(ls, paste0(df_path, "/", idx))) {
    stop("load_h5ad(): ", df_path, " has no index dataset '", idx, "'.",
         call. = FALSE)
  }
  as.character(.h5ad_read_col(file, df_path, idx, ls))
}

.h5ad_df_col <- function(file, df_path, col, ls) {
  e <- .h5_entry(ls, df_path)
  if (e$otype != "H5I_GROUP") {
    df <- rhdf5::h5read(file, df_path)
    if (!(col %in% colnames(df))) {
      stop("load_h5ad(): column '", col, "' not found in ", df_path,
           ". Available: ", .preview(colnames(df), 20), ".", call. = FALSE)
    }
    return(as.vector(df[[col]]))
  }
  .h5ad_read_col(file, df_path, col, ls)
}

## What kind of matrix lives at `path`, and its (n_obs, n_var) shape.
.h5ad_matrix_info <- function(file, path, ls) {
  e <- .h5_entry(ls, path)
  if (nrow(e) == 0L) stop("load_h5ad(): no matrix at '", path, "'.", call. = FALSE)
  at <- .h5_attrs(file, path)
  if (e$otype == "H5I_GROUP") {
    kind  <- at[["encoding-type"]] %||% at[["h5sparse_format"]] %||% ""
    kind  <- sub("_matrix$", "", kind)
    shape <- at[["shape"]] %||% at[["h5sparse_shape"]]
    if (!(kind %in% c("csr", "csc")) || is.null(shape)) {
      stop("load_h5ad(): '", path, "' is a group but not a recognised ",
           "sparse matrix (encoding '", kind, "').", call. = FALSE)
    }
    return(list(kind = kind, n_obs = as.integer(shape[1]),
                n_var = as.integer(shape[2])))
  }
  d <- as.integer(strsplit(e$dim, " x ")[[1]])   # rhdf5 reports R (reversed) dims
  if (length(d) != 2L) {
    stop("load_h5ad(): '", path, "' is not a 2-D matrix (dim ", e$dim, ").",
         call. = FALSE)
  }
  list(kind = "dense", n_obs = d[2], n_var = d[1])
}

## Read the columns `gene_idx` (1-based, in var order) for all cells into
## an n_obs x length(gene_idx) dense double matrix.
.h5ad_read_genes <- function(file, path, info, gene_idx, cells_per_block) {
  n_obs <- info$n_obs; k <- length(gene_idx)
  out <- matrix(0, n_obs, k)
  blocks <- split(seq_len(n_obs), ceiling(seq_len(n_obs) / cells_per_block))

  if (info$kind == "dense") {
    for (rows in blocks) {
      ## rhdf5 returns the dataset with dims reversed: [n_var, n_obs]
      blk <- rhdf5::h5read(file, path, index = list(gene_idx, rows),
                           bit64conversion = "double")
      out[rows, ] <- t(matrix(as.numeric(blk), nrow = k))
    }
  } else if (info$kind == "csr") {
    indptr <- .h5_vec(file, paste0(path, "/indptr"))
    for (rows in blocks) {
      a <- rows[1L]; b <- rows[length(rows)]
      s <- indptr[a] + 1; n <- indptr[b + 1L] - indptr[a]
      if (n == 0) next
      dat <- .h5_vec(file, paste0(path, "/data"),    start = s, count = n)
      ind <- .h5_vec(file, paste0(path, "/indices"), start = s, count = n)
      row <- rep(rows, diff(indptr[a:(b + 1L)]))
      jj  <- match(ind + 1, gene_idx)
      sel <- !is.na(jj)
      if (any(sel)) out[cbind(row[sel], jj[sel])] <- dat[sel]
    }
  } else {                                        # csc: columns are contiguous
    indptr <- .h5_vec(file, paste0(path, "/indptr"))
    for (j in seq_len(k)) {
      g <- gene_idx[j]
      s <- indptr[g] + 1; n <- indptr[g + 1L] - indptr[g]
      if (n == 0) next
      dat <- .h5_vec(file, paste0(path, "/data"),    start = s, count = n)
      ind <- .h5_vec(file, paste0(path, "/indices"), start = s, count = n)
      out[ind + 1, j] <- dat
    }
  }
  out
}


#' Load the cells x (metadata + selected genes) table from a .h5ad file
#'
#' Reads the file directly as HDF5 with `rhdf5` (no Python, no
#' conversion step). Only the requested `obs` columns and the
#' requested genes are read; the full expression matrix is never
#' loaded. Dense and sparse (CSR / CSC) matrices are supported, as are
#' `X`, any entry of `layers`, and `raw.X`.
#'
#' @param h5ad_file Path to the `.h5ad` file.
#' @param genes Character vector of gene IDs to extract (required).
#'   Matched against `gene_id_col`. Any gene not found is an error
#'   that lists the missing IDs and shows what the file's IDs look like.
#' @param obs_cols Character vector of `obs` column names to return
#'   (e.g. the individual and cell-type columns). Each must exist; the
#'   error lists the available columns otherwise. Categorical columns
#'   are returned as character.
#' @param layer Which matrix to read: `"X"` (default), the name of an
#'   entry in `layers`, or `"raw"` for `raw.X` (whose gene IDs then
#'   come from `raw.var`).
#' @param gene_id_col Where the gene IDs live: `NULL` (default) uses
#'   the `var` index (`var_names` in Python); otherwise the name of a
#'   column in `var` (e.g. `"gene_ids"` or `"gene_symbols"`).
#' @param cells_per_block Cells read per HDF5 call for dense / CSR
#'   matrices (controls peak memory of the read, not of the result).
#'
#' @return A `data.table` with one row per cell in file order: the
#'   `obs_cols` first, then one numeric column per gene (in the order
#'   of `genes`).
#' @export
load_h5ad <- function(h5ad_file, genes, obs_cols,
                      layer = "X", gene_id_col = NULL,
                      cells_per_block = 20000L) {

  .require_pkg("rhdf5", "load_h5ad", source = "Bioconductor")
  if (!file.exists(h5ad_file)) {
    stop("h5ad file not found: ", h5ad_file, call. = FALSE)
  }
  if (missing(genes) || is.null(genes) || length(genes) == 0L) {
    stop("load_h5ad(): `genes` is required.", call. = FALSE)
  }
  if (missing(obs_cols) || is.null(obs_cols) || length(obs_cols) == 0L) {
    stop("load_h5ad(): `obs_cols` is required (e.g. c(ind_col, ct_col)).",
         call. = FALSE)
  }
  genes <- unique(as.character(genes))
  ls <- .h5_ls(h5ad_file)

  ## ---- which matrix -------------------------------------------------
  if (identical(layer, "X")) {
    mat_path <- "/X"; var_path <- "/var"
  } else if (identical(layer, "raw")) {
    mat_path <- "/raw/X"; var_path <- "/raw/var"
    if (!.h5_exists(ls, mat_path)) {
      stop("load_h5ad(): the file has no raw.X.", call. = FALSE)
    }
  } else {
    mat_path <- paste0("/layers/", layer); var_path <- "/var"
    if (!.h5_exists(ls, mat_path)) {
      stop("load_h5ad(): no layer '", layer, "'. Available layers: ",
           .preview(.h5_children(ls, "/layers"), 20),
           if (.h5_exists(ls, "/raw/X")) "; also \"raw\"" else "", ".",
           call. = FALSE)
    }
  }
  info <- .h5ad_matrix_info(h5ad_file, mat_path, ls)

  ## ---- gene IDs -----------------------------------------------------
  var_ids <- if (is.null(gene_id_col)) {
    .h5ad_read_index(h5ad_file, var_path, ls)
  } else {
    as.character(.h5ad_df_col(h5ad_file, var_path, gene_id_col, ls))
  }
  if (length(var_ids) != info$n_var) {
    stop("load_h5ad(): ", var_path, " has ", length(var_ids), " gene IDs but ",
         mat_path, " has ", info$n_var, " columns.", call. = FALSE)
  }
  id_src <- if (is.null(gene_id_col)) paste0(var_path, " index") else
            paste0(var_path, "/", gene_id_col)
  miss <- setdiff(genes, var_ids)
  if (length(miss) > 0L) {
    var_cols <- setdiff(.h5_children(ls, var_path),
                        c("__categories", .h5_attrs(h5ad_file, var_path)[["_index"]] %||% "_index"))
    stop("load_h5ad(): ", length(miss), " of ", length(genes), " requested ",
         "gene(s) are not in ", id_src, ": ", .preview(miss),
         ".\n  IDs in the file look like: ", .preview(var_ids, 5),
         if (length(var_cols)) paste0("\n  Other var columns you could use ",
                                      "via gene_id_col: ", .preview(var_cols, 10)) else "",
         call. = FALSE)
  }
  gene_idx <- match(genes, var_ids)
  dup <- genes[genes %in% var_ids[duplicated(var_ids)]]
  if (length(dup) > 0L) {
    stop("load_h5ad(): ", length(dup), " requested gene ID(s) occur more than ",
         "once in ", id_src, ": ", .preview(dup), ". Use a unique ID column.",
         call. = FALSE)
  }

  ## ---- obs columns --------------------------------------------------
  obs <- lapply(obs_cols, function(col) {
    v <- .h5ad_df_col(h5ad_file, "/obs", col, ls)
    if (length(v) != info$n_obs) {
      stop("load_h5ad(): obs column '", col, "' has ", length(v),
           " values but ", mat_path, " has ", info$n_obs, " rows.", call. = FALSE)
    }
    v
  })
  names(obs) <- obs_cols

  message(sprintf(paste0("[load_h5ad] %s: %d cells x %d genes; %s is %s; ",
                         "gene IDs from %s (e.g. %s); reading %d gene(s) ",
                         "and obs column(s) %s."),
                  basename(h5ad_file), info$n_obs, info$n_var, mat_path,
                  info$kind, id_src, .preview(var_ids, 3), length(genes),
                  paste(obs_cols, collapse = ", ")))

  ## ---- expression ---------------------------------------------------
  expr <- .h5ad_read_genes(h5ad_file, mat_path, info, gene_idx,
                           cells_per_block = as.integer(cells_per_block))
  cols <- vector("list", length(genes)); names(cols) <- genes
  for (j in seq_along(genes)) cols[[j]] <- expr[, j]
  rm(expr)

  data.table::setDT(c(obs, cols))[]
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
#'   `write.table(m, "covariates.txt", sep = "\\t", quote = FALSE,
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
