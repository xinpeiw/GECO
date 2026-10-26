# Input: raw/DPG/mmc2.xlsx = Table S1 of Salamone, Tian et al., Cell 2026 (DANDELION DPGs, 27 blood-trait sheets).
suppressMessages({ library(readxl); library(data.table) })

xlsx <- 'raw/DPG/mmc2.xlsx'
all  <- read.table('data/onek1k/GECO_universe.txt')$V1

dpg_sym <- list()
for (sh in excel_sheets(xlsx)) {
  raw <- suppressMessages(read_excel(xlsx, sheet = sh, col_names = FALSE))
  title <- as.character(raw[[1]][1])
  trait <- sub('_GCST.*$', '', sub('^Table S1[A-Z]+[.] ', '', title))
  hdr <- as.character(unlist(raw[2, ]))
  gi <- which(hdr == 'DPGs')
  genes <- unique(na.omit(trimws(as.character(unlist(raw[-(1:2), gi])))))
  dpg_sym[[trait]] <- genes[nzchar(genes)]
}

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }
ensg    = function(x){ e=to_ensg(x); e=e[!is.na(e)&nzchar(e)]; sort(e[!(duplicated(e)|duplicated(e,fromLast=TRUE))]) }

dpg_univ  <- lapply(dpg_sym, function(s) intersect(ensg(unique(s)), all))
immune_tr <- c('WBC','LYMPH','MONO','NEUT','EO','LYMPH_P','MONO_P','NEUT_P','EO_P','BASO_P')
eryPlt_tr <- setdiff(names(dpg_univ), immune_tr)

writeLines(unique(unlist(dpg_univ)),            'data/6e-DPG/DPG_all_blood27.txt')
writeLines(unique(unlist(dpg_univ[immune_tr])), 'data/6e-DPG/DPG_immune.txt')
writeLines(unique(unlist(dpg_univ[eryPlt_tr])), 'data/6e-DPG/DPG_eryPlt.txt')
