# Inputs not included in this repository; download into raw/promoter_FANTOM5/ (FANTOM5 reprocessed hg38, v9):
#   hg38_fair+new_CAGE_peaks_phase1and2_ann.txt.gz  (6,033,029 bytes)
#     https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/extra/CAGE_peaks_expression/hg38_fair%2Bnew_CAGE_peaks_phase1and2_ann.txt.gz
#   hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz  (837,601,892 bytes)
#     https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/extra/CAGE_peaks_expression/hg38_fair%2Bnew_CAGE_peaks_phase1and2_tpm.osc.txt.gz

rm(list=ls())

library(data.table)

root      <- "raw/promoter_FANTOM5"
data_dir  <- root
f5_dir    <- data_dir

ann_path     <- file.path(f5_dir, "hg38_fair+new_CAGE_peaks_phase1and2_ann.txt.gz")
tpm_path     <- file.path(f5_dir, "hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz")
geco_path    <- "data/onek1k/GECO_universe.txt"

CELL_TYPES      <- c("B","CD4","CD8","NK")

F5_MAP <- list(
  B   = c("CD19+ B Cells, donor1.CNhs12343.11544-120B5.hg38.nobarcode",
          "CD19+ B Cells, donor2.CNhs12352.11624-122B4.hg38.nobarcode",
          "CD19+ B Cells, donor3.CNhs12354.11705-123B4.hg38.nobarcode"),
  CD4 = c(
    "CD4+CD25-CD45RA+ naive conventional T cells, donor1.CNhs13223.11784-124B2.hg38.nobarcode",
    "CD4+CD25-CD45RA+ naive conventional T cells, donor2.CNhs13205.11795-124C4.hg38.nobarcode",
    "CD4+CD25-CD45RA+ naive conventional T cells, donor3.CNhs13512.11906-125F7.hg38.nobarcode",
    "CD4+CD25-CD45RA- memory conventional T cells, donor1.CNhs13239.11786-124B4.hg38.nobarcode",
    "CD4+CD25-CD45RA- memory conventional T cells, donor2.CNhs13237.11798-124C7.hg38.nobarcode",
    "CD4+CD25-CD45RA- memory conventional T cells, donor3.CNhs13539.11909-125G1.hg38.nobarcode"),
  CD8 = c("CD8+ T Cells, donor1.CNhs10854.11226-116C2.hg38.nobarcode",
          "CD8+ T Cells, donor2.CNhs11956.11307-117C2.hg38.nobarcode",
          "CD8+ T Cells, donor3.CNhs11999.11383-118B6.hg38.nobarcode"),
  NK  = c("Natural Killer Cells, donor1.CNhs10859.11230-116C6.hg38.nobarcode",
          "Natural Killer Cells, donor2.CNhs11957.11311-117C6.hg38.nobarcode",
          "Natural Killer Cells, donor3.CNhs12001.11387-118C1.hg38.nobarcode"),
  control = c(  "Fibroblast - Dermal, donor1.CNhs12499.11269-116G9.hg38.nobarcode",
                "Fibroblast - Dermal, donor2.CNhs11379.11346-117G5.hg38.nobarcode",
                "Fibroblast - Dermal, donor3.CNhs12028.11418-118F5.hg38.nobarcode",
                "Hepatocyte, donor1.CNhs12340.11523-119I2.hg38.nobarcode",
                "Hepatocyte, donor2.CNhs12349.11603-120I1.hg38.nobarcode",
                "Hepatocyte, donor3.CNhs12626.11684-122I1.hg38.nobarcode",
                "Keratinocyte - epidermal, donor1.CNhs11064.11272-116H3.hg38.nobarcode",
                "Keratinocyte - epidermal, donor2.CNhs11381.11349-117G8.hg38.nobarcode",
                "Keratinocyte - epidermal, donor3.CNhs12031.11421-118F8.hg38.nobarcode",
                "Astrocyte - cerebellum, donor1.CNhs11321.11500-119F6.hg38.nobarcode",
                "Astrocyte - cerebellum, donor2.CNhs12081.11580-120F5.hg38.nobarcode",
                "Astrocyte - cerebellum, donor3.CNhs12117.11661-122F5.hg38.nobarcode",
                "Smooth Muscle Cells - Aortic, donor1.CNhs11085.11283-116I5.hg38.nobarcode",
                "Smooth Muscle Cells - Aortic, donor2.CNhs11305.11360-117I1.hg38.nobarcode",
                "Smooth Muscle Cells - Aortic, donor3.CNhs11309.11432-118H1.hg38.nobarcode")
)

ann <- fread(ann_path, header = TRUE, sep = "\t", skip = "00Annotation")
setnames(ann, "00Annotation", "peak_id")
ann <- ann[, .(peak_id, hgnc_id)]
ann[hgnc_id == "NA", hgnc_id := NA_character_]

ann_hgnc <- ann[
  !is.na(hgnc_id) & nzchar(hgnc_id),
  .( hgnc_id = unlist( regmatches( hgnc_id, gregexpr("HGNC:[0-9]+", hgnc_id, perl = TRUE)))),
  by = peak_id
]
ann_hgnc <- unique(ann_hgnc, by = c("peak_id", "hgnc_id"))

ensg2hgnc <- fread('raw/hgnc_complete_set.txt')
ensg2hgnc <- ensg2hgnc[
  !is.na(ensembl_gene_id) & nzchar(ensembl_gene_id) &
  !is.na(hgnc_id) & nzchar(hgnc_id)
]
ensg2hgnc <- unique(ensg2hgnc[, .(ensembl_gene_id, hgnc_id, symbol)], by = "ensembl_gene_id")
geco <- readLines(geco_path)
mappable_geco <- intersect(geco, ensg2hgnc$ensembl_gene_id)

unique_samples <- unique(unlist(F5_MAP))
hdr <- system(sprintf("zcat < %s | grep -v '^#' | head -1", shQuote(tpm_path)),
              intern = TRUE)
all_cols <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
col_decoded <- vapply(all_cols, function(c) {
  if (startsWith(c, "tpm.")) URLdecode(substring(c, 5)) else c
}, character(1))
sel <- c(1L, match(unique_samples, col_decoded))
if (anyNA(sel)) {
  stop("Missing FANTOM5 TPM columns: ",
       paste(unique_samples[is.na(match(unique_samples, col_decoded))],
             collapse = ", "))
}
tpm_dt <- fread(cmd = sprintf("zcat < %s | grep -v '^#'", shQuote(tpm_path)),
                header = TRUE, select = sel, sep = "\t")
setnames(tpm_dt, c("peak_id", unique_samples))
tpm_dt <- tpm_dt[grepl("::", peak_id, fixed = TRUE)]

peak2gene <- merge(
  ann_hgnc[, .(peak_id, hgnc_id)],
  ensg2hgnc[, .(hgnc_id, ensembl_gene_id)],
  by = "hgnc_id",
  allow.cartesian = TRUE
)[, .(peak_id, ensembl_gene_id)]
peak2gene <- unique(peak2gene)
peak2gene <- peak2gene[ensembl_gene_id %in% geco]
peak2gene <- peak2gene[peak_id %in% tpm_dt$peak_id]

per_gene_n <- list()
for (ct in CELL_TYPES) {
  cols <- F5_MAP[[ct]]
  m <- as.matrix(tpm_dt[, ..cols])
  active <- tpm_dt$peak_id[rowSums(m > 0) > 0]
  ct_links <- peak2gene[peak_id %in% active]
  counts <- ct_links[, .(peak_num = uniqueN(peak_id)), by = ensembl_gene_id]

  pg <- data.table(ensembl_gene_id = geco)
  pg <- counts[pg, on = "ensembl_gene_id"]
  mappable <- pg$ensembl_gene_id %in% mappable_geco
  pg[mappable & is.na(peak_num), peak_num := 0L]

  stopifnot(
    nrow(pg) == length(geco),
    uniqueN(pg$ensembl_gene_id) == length(geco),
    all(is.na(pg$peak_num[!mappable]))
  )
  per_gene_n[[ct]] <- pg
}

for (ct in names(per_gene_n)) {
  dt = per_gene_n[[ct]]
  write.table(dt, file.path('data/4b-promoter_enhancer_TFBS/promoter_FANTOM5', paste0(ct, '.txt')),
              row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
}
