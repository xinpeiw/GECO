# Input not included in this repository: ENCODE-rE2G BED files listed in raw/directional/re2g_files_manifest.tsv.
#   Download them first with: bash raw/directional/02_download_re2g.sh  (writes raw/directional/encode_re2g_multi/)

rm(list=ls())

geco_dir    <- "data/onek1k"
re2g_dir    <- "raw/directional/encode_re2g_multi"
manifest_fp <- "raw/directional/re2g_files_manifest.tsv"
SCORE_MIN <- 0.164
CELL_TYPES <- c("CD4NC","CD4ET","CD8NC","CD8ET","NK","BIN","BMem")
universe    <- readLines(file.path(geco_dir, "GECO_universe.txt"))

manifest <- read.table(manifest_fp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
manifest$path <- file.path(re2g_dir, sprintf("%s_%s_%s.bed.gz",
                                             manifest$hub, manifest$annotation, manifest$file_accession))

read_re2g <- function(path) {
  if (!file.exists(path)) stop("Missing rE2G file: ", path)
  gzip_status <- system2("gunzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE)
  if (gzip_status != 0L) stop("Invalid gzip file: ", path)
  cmd <- sprintf(
    "gunzip -c '%s' | awk -F'\\t' 'NR>1 && $5!=\"promoter\" && $9==\"FALSE\" && ($NF+0)>=%f {print $7\"\\t\"$1\":\"$2\"-\"$3}'",
    path, SCORE_MIN
  )
  con <- pipe(cmd, "r")
  d <- read.table(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("ensg","region_id"), quote = "")
  close(con)
  if (nrow(d) == 0) return(d)
  d$ensg <- sub("\\..*$", "", d$ensg)
  d
}

count_enh_per_gene <- function(focal_ct, gene_universe = universe) {
  donors <- manifest$path[manifest$hub == focal_ct]
  cat(sprintf("  [%s] reading %d donors...\n", focal_ct, length(donors)))
  if (length(donors) == 0L) {
    return(data.frame( ensg = gene_universe, n_enh = integer(length(gene_universe)), stringsAsFactors = FALSE))
  }
  all_pairs <- lapply(donors, read_re2g)
  all_pairs <- Filter(function(x) !is.null(x) && nrow(x) > 0L, all_pairs)
  if (length(all_pairs) == 0L) {
    return(data.frame(
      ensg = gene_universe, n_enh = integer(length(gene_universe)), stringsAsFactors = FALSE
    ))
  }
  pairs <- do.call(rbind, all_pairs)
  pairs <- pairs[
    pairs$ensg %in% gene_universe & !is.na(pairs$ensg) & nzchar(pairs$ensg) & !is.na(pairs$region_id) & nzchar(pairs$region_id),
    , drop = FALSE
  ]
  pairs <- unique(pairs[, c("ensg", "region_id")])

  agg <- aggregate( region_id ~ ensg, data = pairs, FUN = function(x) length(unique(x)) )
  names(agg)[names(agg) == "region_id"] <- "n_enh"

  out <- merge(
    data.frame(ensg = gene_universe, stringsAsFactors = FALSE), agg,
    by = "ensg", all.x = TRUE, sort = FALSE )
  out$n_enh[is.na(out$n_enh)] <- 0L
  out$n_enh <- as.integer(out$n_enh)
  out
}

for (ct in CELL_TYPES) {
  tmp <- count_enh_per_gene(ct)
  write.table(tmp, file.path('data/4b-promoter_enhancer_TFBS/enhancer_rE2G', paste0(ct,'.txt')),
              row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
}
