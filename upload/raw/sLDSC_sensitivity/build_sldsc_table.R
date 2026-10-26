# Requires raw/sLDSC_sensitivity/results/ from raw/sLDSC_sensitivity/run_sldsc.sh (not included; the regression
# outputs of the run behind data/5b-ldsc/sensitivity/sldsc_tau.tsv were not retained).
suppressMessages(library(data.table))

M <- sum(vapply(Sys.glob('raw/sLDSC/baseline_M_5_50/baseline.*.l2.M_5_50'),
                function(f) scan(f, quiet = TRUE, nmax = 1), numeric(1)))

SETTINGS <- c(primary = 'primary', S1 = 'baselineLD2.2', S3 = 'window500kb', S4 = 'MHCincl',
              top210 = 'top210', top420 = 'top420')

parse <- function(rfile) {
  lf <- sub('\\.results$', '.log', rfile);                                 if (!file.exists(lf)) stop(lf, ': .log missing')
  focal <- fread(rfile)[Category == 'L2_1'];                               if (nrow(focal) != 1)   stop(rfile, ': expected one L2_1 row, got ', nrow(focal))
  h2line <- grep('Total Observed scale h2', readLines(lf), value = TRUE)[1]; if (is.na(h2line))    stop(lf, ': no "Total Observed scale h2" line')
  v <- as.numeric(strsplit(sub('.*h2: ', '', h2line), '[ ()]+')[[1]])
  scale <- sqrt(focal$`Prop._SNPs` * (1 - focal$`Prop._SNPs`)) * M / v[1]
  tag <- sub('.*\\.', '', sub('\\.results$', '', basename(rfile)))
  data.table(
    `FinnGen ID` = sub(paste0('\\.', tag, '\\.results$'), '', basename(rfile)),
    tag,
    h2_z        = v[1] / v[2],
    tau_star    = focal$Coefficient           * scale,
    tau_star_se = focal$Coefficient_std_error * scale
  )
}

hub_tags <- as.vector(outer(c('Shared', 'Specific'), names(SETTINGS), paste, sep = '_'))
ct_tags  <- paste0('ct_', c('BIN', 'BMem', 'CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK'))
files <- unlist(lapply(file.path('raw/sLDSC_sensitivity/results', paste0('*.', c(hub_tags, ct_tags), '.results')), Sys.glob))

d <- rbindlist(lapply(files, parse))
ct <- startsWith(d$tag, 'ct_')
d[ ct, `:=`(hub_set = sub('^ct_', '', tag),       ldsc_settings = 'primary')]
d[!ct, `:=`(hub_set = sub('_.*', '', tag), ldsc_settings = unname(SETTINGS[sub('.*_', '', tag)]))]
d[, tag := NULL]
setcolorder(d, c('FinnGen ID', 'hub_set', 'ldsc_settings', 'h2_z', 'tau_star', 'tau_star_se'))

fwrite(d, 'data/5b-ldsc/sensitivity/sldsc_tau.tsv', sep = '\t')
