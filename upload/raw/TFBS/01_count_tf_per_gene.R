rm(list=ls())
library(data.table)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)

CELL_TYPES <- c("CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")
WINDOW <- 100000L
SKIP_GZIP_CHECK <- Sys.getenv("SKIP_GZIP_CHECK", "0") %in% c("1", "TRUE", "true", "yes")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  CT <- "CD4NC"
} else if (length(args) == 1L && args[1] %in% CELL_TYPES) {
  CT <- args[1]
} else {
  stop("usage: Rscript 01_count_tf_per_gene.R <cell_type>  (one of ",
       paste(CELL_TYPES, collapse = ", "), ")")
}

ROOT      <- "raw/TFBS"
META_FILE <- file.path(ROOT, "trace_metadata", paste0(CT, ".tsv"))
TRACE_DIR <- file.path(ROOT, "trace_footprints", CT)
GECO_FILE <- "data/onek1k/GECO_universe.txt"
OUT_DIR   <- "data/4b-promoter_enhancer_TFBS/TF_TRACE"
OUT       <- file.path(OUT_DIR, paste0("per_gene_tf_diversity_", CT, ".tsv"))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

read_tsv_local <- function(path, ...) {
  if (!file.exists(path)) stop("missing file: ", path)
  txt <- readLines(path, warn = FALSE)
  if (length(txt) == 0L) stop("empty or unreadable file: ", path)
  fread(text = paste(txt, collapse = "\n"), sep = "\t", ...)
}

read_gene_list <- function(path) {
  if (!file.exists(path)) stop("missing file: ", path)
  genes <- trimws(readLines(path, warn = FALSE))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0L) stop("empty or unreadable gene list: ", path)
  unique(sub("\\..*$", "", genes))
}

parse_targets <- function(s) {
  if (is.na(s) || !nzchar(trimws(s))) return(character(0))
  tokens <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  sub("-human$", "", tokens[nzchar(tokens)])
}
meta <- read_tsv_local(META_FILE)
req_cols <- c("File accession", "Targets")
missing_cols <- setdiff(req_cols, names(meta))
if (length(missing_cols) > 0L) {
  stop("metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
}
meta[, TF_list := lapply(Targets, parse_targets)]
file2tfs <- setNames(meta$TF_list, meta$`File accession`)
all_tfs  <- sort(unique(unlist(meta$TF_list)))
tf2idx   <- setNames(seq_along(all_tfs), all_tfs)
n_multi  <- sum(lengths(file2tfs) > 1L)
cat(sprintf("[%s] distinct TFs: %d;  files in metadata: %s  (multi-TF: %d)\n",
            CT, length(all_tfs), format(length(file2tfs), big.mark = ","), n_multi))

cat(sprintf("[%s] reading TRACE BEDs from %s ...\n", CT, TRACE_DIR))
if (!dir.exists(TRACE_DIR)) stop("missing TRACE directory: ", TRACE_DIR)
beds <- list.files(TRACE_DIR, pattern = "\\.bed\\.gz$", full.names = TRUE)
if (length(beds) == 0L) stop("no .bed.gz files found in: ", TRACE_DIR)
if (!SKIP_GZIP_CHECK) {
  for (bed in beds) {
    if (system2("gzip", c("-t", bed), stdout = FALSE, stderr = FALSE) != 0L) {
      stop("corrupt or incomplete BED gzip file: ", bed)
    }
  }
} else {
  cat(sprintf("[%s] skipping gzip integrity precheck (SKIP_GZIP_CHECK=1)\n", CT))
}
parsed <- 0L; missing_meta <- 0L
density_list   <- vector("list", length(beds))
diversity_list <- vector("list", length(beds))
for (i in seq_along(beds)) {
  acc <- sub("\\.bed\\.gz$", "", basename(beds[i]))
  tfs <- file2tfs[[acc]]
  if (is.null(tfs) || length(tfs) == 0L) { missing_meta <- missing_meta + 1L; next }
  dt <- fread(cmd = paste("gunzip -c", shQuote(beds[i])),
              header = FALSE, sep = "\t", fill = TRUE,
              col.names = c("chrom","start","end","score","tf_name"),
              colClasses = c("character","integer","integer","character","character"))
  parsed <- parsed + 1L
  if (is.null(dt) || nrow(dt) == 0L) next
  mids    <- (dt$start + dt$end) %/% 2L
  indices <- unname(unlist(tf2idx[tfs]))
  density_list[[i]] <- data.table(chrom = dt$chrom, mid = mids)
  n_tfs <- length(indices); n_rows <- nrow(dt)
  diversity_list[[i]] <- data.table(
    chrom  = rep(dt$chrom, each  = n_tfs),
    mid    = rep(mids,     each  = n_tfs),
    tf_idx = rep(indices,  times = n_rows))
}
density_fp   <- rbindlist(density_list,   use.names = TRUE, fill = FALSE)
diversity_fp <- rbindlist(diversity_list, use.names = TRUE, fill = FALSE)
if (nrow(density_fp) == 0L || nrow(diversity_fp) == 0L) {
  stop("no usable footprints were loaded for ", CT)
}
cat(sprintf("[%s] parsed %d BEDs (missing TF metadata for %d); footprints: %s\n",
            CT, parsed, missing_meta, format(nrow(density_fp), big.mark = ",")))

setorder(density_fp,   chrom, mid)
setorder(diversity_fp, chrom, mid)
density_per_chrom_mid   <- split(density_fp$mid,     density_fp$chrom)
diversity_per_chrom_mid <- split(diversity_fp$mid,    diversity_fp$chrom)
diversity_per_chrom_tf  <- split(diversity_fp$tf_idx, diversity_fp$chrom)

geco_universe <- read_gene_list(GECO_FILE)
cat(sprintf("[%s] GECO universe genes: %s\n", CT, format(length(geco_universe), big.mark = ",")))

cat(sprintf("[%s] loading EnsDb.Hsapiens.v86 genes ...\n", CT))
genes_gr <- genes(EnsDb.Hsapiens.v86)
keep <- as.character(seqnames(genes_gr)) %in% c(as.character(1:22), "X", "Y")
genes_gr <- genes_gr[keep]
gene_ids <- names(genes_gr)
genes_gr <- GRanges(
  seqnames = paste0("chr", as.character(seqnames(genes_gr))),
  ranges = ranges(genes_gr),
  strand = strand(genes_gr),
  gene_name = genes_gr$gene_name,
  gene_type = genes_gr$gene_biotype
)
names(genes_gr) <- gene_ids

g <- data.table(
  chrom = as.character(seqnames(genes_gr)),
  strand = as.character(strand(genes_gr)),
  gene_id = names(genes_gr),
  gene_name = genes_gr$gene_name,
  gene_type = genes_gr$gene_type,
  gs = start(genes_gr) - 1L,
  ge = end(genes_gr)
)
g[is.na(gene_name), gene_name := gene_id]
g[, tss := fifelse(strand == "+", gs, ge - 1L)]
cat(sprintf("[%s] %s genes loaded from EnsDb.Hsapiens.v86\n", CT, format(nrow(g), big.mark = ",")))

universe_dt <- data.table(gene_id = geco_universe, universe_order = seq_along(geco_universe))
g <- merge(universe_dt, g, by = "gene_id", all.x = TRUE, sort = FALSE)
setorder(g, universe_order)
n_missing_tss <- sum(is.na(g$tss))
cat(sprintf("[%s] retained %s GECO genes; missing EnsDb.Hsapiens.v86 TSS for %s\n",
            CT, format(nrow(g), big.mark = ","), format(n_missing_tss, big.mark = ",")))

cat(sprintf("[%s] counting TFs in TSS ±%dkb ...\n", CT, WINDOW %/% 1000L))
n_TF <- rep(NA_integer_, nrow(g))
n_FP <- rep(NA_integer_, nrow(g))
has_tss <- !is.na(g$tss)
n_TF[has_tss] <- 0L
n_FP[has_tss] <- 0L
for (i in seq_len(nrow(g))) {
  if (is.na(g$tss[i])) next
  ch <- g$chrom[i]
  if (is.null(density_per_chrom_mid[[ch]])) next
  win_s <- max(0L, g$tss[i] - WINDOW)
  win_e <-          g$tss[i] + WINDOW
  dmids <- density_per_chrom_mid[[ch]]
  dlo <- findInterval(win_s, dmids, left.open = TRUE) + 1L
  dhi <- findInterval(win_e, dmids) + 1L
  n_FP[i] <- max(0L, dhi - dlo)
  vmids <- diversity_per_chrom_mid[[ch]]
  vtfs  <- diversity_per_chrom_tf[[ch]]
  vlo <- findInterval(win_s, vmids, left.open = TRUE) + 1L
  vhi <- findInterval(win_e, vmids) + 1L
  if (vhi > vlo) n_TF[i] <- length(unique(vtfs[vlo:(vhi - 1L)]))
}

out <- data.table(
  gene_id                 = g$gene_id,
  gene_name               = g$gene_name,
  gene_type               = g$gene_type,
  chrom                   = g$chrom,
  tss                     = g$tss,
  strand                  = g$strand,
  n_TF_in_window          = n_TF,
  n_footprints_in_window  = n_FP
)
fwrite(out, OUT, sep = "\t", quote = FALSE)
cat(sprintf("[%s] wrote %s\n", CT, OUT))

s <- summary(out$n_TF_in_window)
cat(sprintf("  n_TF_in_window summary:  min=%g  Q1=%g  median=%g  mean=%.2f  Q3=%g  max=%g  NA=%d\n",
            s["Min."], s["1st Qu."], s["Median"], mean(out$n_TF_in_window, na.rm = TRUE),
            s["3rd Qu."], s["Max."], sum(is.na(out$n_TF_in_window))))
