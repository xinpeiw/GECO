# Requires the LDSC reference from raw/sLDSC/01_prep_reference.sh (raw/sLDSC/reference/).
suppressMessages(library(data.table))

REF     <- 'raw/sLDSC/reference'
OUT     <- file.path(REF, 'CLUES_annot_mega')
dir.create(OUT, showWarnings = FALSE)
WINDOW  <- 100e3

sets <- list(CLUES_Shared   = readLines('data/clues/hubs_alpha.txt'),
             CLUES_Specific = readLines('data/clues/hubs_specific.txt'),
             CLUES_Universe = readLines('data/clues/allgenes_CLUES.txt'))

pos <- fread('raw/sLDSC/gene_position_GRCh37.tsv',
             select = c(1, 2, 3, 6), col.names = c('chr', 'start', 'end', 'gene'))
pos[, chr := sub('^chr', '', chr)]
pos <- pos[chr %in% as.character(1:22)]

beds <- lapply(sets, function(g)
  pos[gene %in% g, .(chr, start = pmax(0, start - WINDOW), end = end + WINDOW)])

for (c in 1:22) {
  bim <- fread(sprintf('%s/1000G_EUR_Phase3_plink/1000G.EUR.QC.%d.bim', REF, c),
               col.names = c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'))
  snp <- bim[, .(SNP, start = BP, end = BP)]; setkey(snp, start, end)
  for (s in names(sets)) {
    out <- sprintf('%s/%s.%d.annot.gz', OUT, s, c)
    if (file.exists(out)) next
    iv <- beds[[s]][chr == as.character(c), .(start, end)]
    setkey(iv, start, end)
    hit   <- foverlaps(snp, iv, type = 'any', which = TRUE, nomatch = NA)
    annot <- as.integer(!is.na(hit$yid[match(seq_len(nrow(bim)), hit$xid)]))
    fwrite(bim[, .(CHR, BP, SNP, CM)][, (s) := annot][], out,
           sep = '\t', quote = FALSE, compress = 'gzip')
  }
  cat('chr', c, 'done\n')
}
