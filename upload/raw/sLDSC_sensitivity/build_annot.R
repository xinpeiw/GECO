# Requires the LDSC reference (raw/sLDSC/01_prep_reference.sh and raw/sLDSC_sensitivity/prep.sh).
suppressMessages(library(data.table))
REF   = 'raw/sLDSC/reference'
TBL   = 'raw/sLDSC_sensitivity'
CTS   = c('BIN','BMem','CD4ET','CD4NC','CD8ET','CD8NC','NK')

sets = c(
  list(GECO_Shared          = list('data/onek1k/GECO_Shared.txt',   100), GECO_Specific        = list('data/onek1k/GECO_Specific.txt', 100),
       GECO_Universe        = list('data/onek1k/GECO_universe.txt', 100), GECO_Shared_500kb    = list('data/onek1k/GECO_Shared.txt',   500),
       GECO_Specific_500kb  = list('data/onek1k/GECO_Specific.txt', 500), GECO_Universe_500kb  = list('data/onek1k/GECO_universe.txt', 500),
       GECO_Shared_top210   = list(file.path(TBL, 'GECO_Shared_top210.txt'),   100), GECO_Specific_top210 = list(file.path(TBL, 'GECO_Specific_top210.txt'), 100),
       GECO_Shared_top420   = list(file.path(TBL, 'GECO_Shared_top420.txt'),   100), GECO_Specific_top420 = list(file.path(TBL, 'GECO_Specific_top420.txt'), 100)),
  setNames(lapply(CTS, function(ct) list(paste0('data/onek1k/hubs_', ct, '.txt'), 100)), paste0('GECO_specific_', CTS)))

pos = fread('raw/sLDSC/gene_position_GRCh37.tsv', select = c(1,2,3,6),
            col.names = c('chr','start','end','ensembl_gene_id'))
pos[, chr := sub('^chr', '', chr)]; pos = pos[chr %in% as.character(1:22)]

beds = lapply(sets, function(s) {
  g = readLines(s[[1]])
  pos[ensembl_gene_id %in% g, .(chr, start = pmax(0L, start - s[[2]]*1000L), end = end + s[[2]]*1000L)]
})

dir.create(file.path(REF, 'GECO_annot_sensitivity'), showWarnings = FALSE)
for (chrom_i in 1:22) {
  bim_file = sprintf('%s/1000G.EUR.QC.%d.bim', file.path(REF, '1000G_EUR_Phase3_plink'), chrom_i)
  bim = fread(bim_file, col.names = c('CHR','SNP','CM','BP','A1','A2'))
  snp = bim[, .(SNP, start = BP, end = BP)]; setkey(snp, start, end)
  for (set_name in names(sets)) {
    out = sprintf('%s/%s.%d.annot.gz', file.path(REF, 'GECO_annot_sensitivity'), set_name, chrom_i)
    if (file.exists(out)) next
    iv = beds[[set_name]][chr == as.character(chrom_i), .(start, end)]
    if (nrow(iv)) { setkey(iv, start, end)
      hit = foverlaps(snp, iv, type = 'any', which = TRUE, nomatch = NA)
      annot = as.integer(!is.na(hit$yid[match(seq_len(nrow(bim)), hit$xid)]))
    } else annot = rep(0L, nrow(bim))
    fwrite(bim[, .(CHR, BP, SNP, CM)][, (set_name) := annot][], out, sep = '\t', quote = FALSE, compress = 'gzip')
  }
}
