# Input not included in this repository: ENCODE-rE2G BED files listed in raw/CLUES_evolution/re2g_files_manifest.tsv.
#   Download them first with: bash raw/CLUES_evolution/download_re2g.sh  (writes raw/CLUES_evolution/encode_re2g_multi/)

rm(list=ls())

library(GenomicRanges)

re2g_dir <- 'raw/CLUES_evolution/encode_re2g_multi'
region_dir = 'raw/directional'
manifest <- read.table('raw/CLUES_evolution/re2g_files_manifest.tsv', header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
colnames(manifest)[1] = 'hub'
manifest$path <- file.path(re2g_dir, sprintf("%s_%s_%s.bed.gz",
                                             manifest$hub, manifest$annotation,
                                             manifest$file_accession))
universe = readLines('data/clues/allgenes_CLUES.txt')

SCORE_MIN <- 0.164
read_regions <- function(path, is.bed = FALSE) {
  df <- read.table(path, header = TRUE, stringsAsFactors = FALSE,
                   sep = "", fill = TRUE, comment.char = "")
  df <- df[complete.cases(df[, c("chrom","start","end")]), ]
  if (is.bed) {
    GRanges(df$chrom, IRanges(start = df$start + 1, end = df$end))
  } else {
    GRanges(df$chrom, IRanges(start = df$start, end = df$end))
  }
}
regions <- list(
  HAQER           = read_regions(file.path(region_dir,"HAQER_region.txt"), is.bed = TRUE),
  zooHAR_positive = read_regions(file.path(region_dir,"zoohar_positive_region.txt"), is.bed = TRUE),
  akbari = read_regions(file.path(region_dir,"akbari_noanno_region.txt"), is.bed = FALSE)
)
read_re2g <- function(path) {
  if (!file.exists(path)) stop("Missing rE2G file: ", path)
  gzip_status <- system2("gunzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE)
  if (gzip_status != 0L) stop("Invalid gzip file: ", path)
  cmd <- sprintf(
    "gunzip -c '%s' | awk -F'\\t' 'NR>1 && $5!=\"promoter\" && $9==\"FALSE\" && ($NF+0)>=%f {print $1\"\\t\"$2\"\\t\"$3\"\\t\"$7}'",
    path, SCORE_MIN
  )
  con <- pipe(cmd, "r")
  d <- read.table(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("chr","start","end","ensg"))
  close(con)
  GRanges(d$chr, IRanges(start = d$start+1, end = d$end),
          ensg = sub("\\..*$", "", d$ensg))
}

per_donor <- list()
for (i in seq_len(nrow(manifest))) {
  rec <- manifest[i, ]
  gr  <- read_re2g(rec$path)
  for (rn in names(regions)) {
    h <- findOverlaps(gr, regions[[rn]], ignore.strand = TRUE)
    per_donor[[rec$hub]][[rec$annotation]][[rn]] <- intersect(unique(gr$ensg[queryHits(h)]), universe)
  }
}

linked <- list()
hub_names   <- unique(manifest$hub)
for (hub in hub_names) for (rn in names(regions)) {
  donor_sets <- lapply(per_donor[[hub]], function(d) d[[rn]])
  linked[[paste(rn, hub, sep = "::")]] <- Reduce(union, donor_sets)
}

akbari_anno = readLines(file.path(region_dir,'akbari_anno_gene.txt'))

for (nm in names(linked)) {
  parts <- strsplit(nm, "::", fixed = TRUE)[[1]]
  region_type <- parts[1]
  hub <- parts[2]

  genes <- linked[[nm]]
  genes <- genes[!is.na(genes) & nzchar(genes)]
  genes <- unique(genes)

  if (region_type == 'akbari') { genes = union(genes, akbari_anno) }

  path = file.path('data/others/evolution_clues', region_type)
  write.table(genes, file.path(path, paste0(hub, "_enhanceronly.txt")),  row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

}
