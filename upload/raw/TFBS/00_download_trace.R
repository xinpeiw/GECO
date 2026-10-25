suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

META_DIR  <- "raw/TFBS/trace_metadata"
TRACE_DIR <- "raw/TFBS/trace_footprints"
dir.create(META_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(TRACE_DIR, showWarnings = FALSE, recursive = TRUE)

CONFIG <- list(
  CD4NC = c("naive thymus-derived CD4-positive, alpha-beta T cell",
            "central memory CD4-positive, alpha-beta T cell"),
  CD4ET = c("effector memory CD4-positive, alpha-beta T cell"),
  CD8NC = c("CD8-positive naive resting alpha-beta T cell",
            "central memory CD8-positive, alpha-beta T cell"),
  CD8ET = c("effector memory CD8-positive, alpha-beta T cell"),
  NK    = c("natural killer cell"),
  BIN   = c("naive B cell"),
  BMem  = c("memory B cell")
)

FILTER_PARAMS <- c(
  "type=Annotation",
  "software_used.software.name=trace",
  "status=released",
  "assembly=GRCh38",
  "biosample_ontology.classification=primary cell",
  "organism.scientific_name=Homo sapiens"
)

metadata_cache_key <- function(ct) {
  paste(c(FILTER_PARAMS, paste0("biosample_ontology.term_name=", CONFIG[[ct]])),
        collapse = "\n")
}

fetch_metadata <- function(biosample) {
  base <- "https://www.encodeproject.org/metadata/"
  q <- paste0(
    "?type=Annotation",
    "&software_used.software.name=trace",
    "&biosample_ontology.term_name=", URLencode(biosample, reserved = TRUE),
    "&status=released",
    "&assembly=GRCh38",
    "&biosample_ontology.classification=primary%20cell",
    "&organism.scientific_name=Homo%20sapiens"
  )
  url <- paste0(base, q)
  tmp <- tempfile(fileext = ".tsv")
  utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
  txt <- readLines(tmp, warn = FALSE)
  unlink(tmp)
  txt
}

metadata_is_valid <- function(path, cache_key = NULL) {
  if (!file.exists(path)) return(FALSE)
  if (!is.null(cache_key)) {
    key_path <- paste0(path, ".query")
    if (!file.exists(key_path)) return(FALSE)
    old_key <- paste(readLines(key_path, warn = FALSE), collapse = "\n")
    if (!identical(old_key, cache_key)) return(FALSE)
  }
  txt <- readLines(path, warn = FALSE)
  if (length(txt) < 2L) return(FALSE)
  header <- strsplit(txt[1L], "\t", fixed = TRUE)[[1]]
  all(c("File accession", "File download URL", "Targets") %in% header)
}

gzip_is_valid <- function(path) {
  file.exists(path) &&
    file.info(path)$size > 0 &&
    system2("gzip", c("-t", path), stdout = FALSE, stderr = FALSE) == 0L
}

download_bed <- function(url, out_path) {
  tmp <- paste0(out_path, ".part")
  if (file.exists(tmp)) unlink(tmp)
  status <- system2("curl",
                    c("-fL", "--retry", "5", "--retry-delay", "2",
                      "-o", tmp, url),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0L || !gzip_is_valid(tmp)) {
    if (file.exists(tmp)) unlink(tmp)
    return(FALSE)
  }
  if (file.exists(out_path)) unlink(out_path)
  file.rename(tmp, out_path)
}

cat("Fetching TRACE metadata for each biosample...\n\n")
for (ct in names(CONFIG)) {
  out_path <- file.path(META_DIR, paste0(ct, ".tsv"))
  cache_key <- metadata_cache_key(ct)
  if (metadata_is_valid(out_path, cache_key)) {
    cat(sprintf("  %-5s (cached): %s\n", ct, out_path))
    next
  }
  if (file.exists(out_path)) unlink(out_path)
  key_path <- paste0(out_path, ".query")
  if (file.exists(key_path)) unlink(key_path)
  header <- NULL
  body   <- character()
  for (bs in CONFIG[[ct]]) {
    txt <- fetch_metadata(bs)
    if (length(txt) == 0L) next
    if (is.null(header)) header <- txt[1L]
    body <- c(body, txt[-1L])
  }
  if (is.null(header) || length(body) == 0L) {
    stop("no metadata returned for ", ct)
  }
  meta_dt <- fread(text = paste(c(header, body), collapse = "\n"), sep = "\t")
  req_cols <- c("File accession", "File download URL", "Targets")
  missing_cols <- setdiff(req_cols, names(meta_dt))
  if (length(missing_cols) > 0L) {
    stop("metadata for ", ct, " is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  meta_dt <- unique(meta_dt, by = "File accession")
  if (any(is.na(meta_dt$`File download URL`) | !nzchar(meta_dt$`File download URL`))) {
    stop("metadata for ", ct, " contains missing File download URL values")
  }
  fwrite(meta_dt, out_path, sep = "\t", quote = FALSE)
  writeLines(cache_key, key_path)
  if (!metadata_is_valid(out_path, cache_key)) stop("downloaded invalid metadata: ", out_path)
  cat(sprintf("  %-5s (%s): %d annotations -> %s\n",
              ct, paste(CONFIG[[ct]], collapse = "; "), nrow(meta_dt), out_path))
}

cat("\nDownloading per-cell-type BED files...\n")
for (ct in names(CONFIG)) {
  out_dir <- file.path(TRACE_DIR, ct)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  meta <- fread(file.path(META_DIR, paste0(ct, ".tsv")), sep = "\t")
  acc  <- meta$`File accession`
  url  <- meta$`File download URL`
  bed_path <- file.path(out_dir, paste0(acc, ".bed.gz"))
  ok <- vapply(bed_path, gzip_is_valid, logical(1))
  if (any(!ok & file.exists(bed_path))) unlink(bed_path[!ok & file.exists(bed_path)])
  todo <- data.table(url = url[!ok], bed_path = bed_path[!ok])
  cat(sprintf("  %-5s: %d total, %d to download\n",
              ct, length(acc), nrow(todo)))
  if (nrow(todo) == 0L) next
  cores <- min(12L, nrow(todo))
  success <- if (.Platform$OS.type == "unix" && cores > 1L) {
    unlist(mclapply(seq_len(nrow(todo)), function(i) {
      download_bed(todo$url[i], todo$bed_path[i])
    }, mc.cores = cores), use.names = FALSE)
  } else {
    vapply(seq_len(nrow(todo)), function(i) {
      download_bed(todo$url[i], todo$bed_path[i])
    }, logical(1))
  }
  if (!all(success)) {
    stop("failed downloads for ", ct, ":\n",
         paste(todo$bed_path[!success], collapse = "\n"))
  }
}

cat("\nDone. Per-CT counts:\n")
for (ct in names(CONFIG)) {
  files <- list.files(file.path(TRACE_DIR, ct), pattern = "\\.bed\\.gz$", full.names = TRUE)
  n_have <- sum(vapply(files, gzip_is_valid, logical(1)))
  meta_n <- nrow(fread(file.path(META_DIR, paste0(ct, ".tsv")), sep = "\t"))
  cat(sprintf("  %-5s: %d / %d\n", ct, n_have, meta_n))
}
