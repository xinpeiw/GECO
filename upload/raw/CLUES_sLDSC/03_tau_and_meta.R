# raw/CLUES_sLDSC/results/ (the 156 S-LDSC regressions) is included, so this step runs without the LDSC reference.
suppressMessages(library(data.table))

RES <- 'raw/CLUES_sLDSC/results'
OUT <- 'data/7a-ldsc'

M <- sum(vapply(Sys.glob('raw/sLDSC/baseline_M_5_50/baseline.*.l2.M_5_50'),
                function(f) scan(f, quiet = TRUE, nmax = 1), numeric(1)))

parse <- function(rfile) {
  focal <- fread(rfile)[Category == 'L2_1']
  lf <- sub('\\.results$', '.log', rfile)
  h2line <- grep('Total Observed scale h2', readLines(lf), value = TRUE)[1]
  v <- as.numeric(strsplit(sub('.*h2: ', '', h2line), '[ ()]+')[[1]])
  base <- sub('\\.results$', '', basename(rfile))
  scale <- sqrt(focal$Prop._SNPs * (1 - focal$Prop._SNPs)) * M / v[1]
  data.table(`FinnGen ID` = sub('\\.(Shared|Specific)_primary$', '', base),
             hub = sub('.*\\.(Shared|Specific)_primary$', '\\1', base),
             h2_z = v[1] / v[2],
             tau_star = focal$Coefficient * scale,
             tau_star_se = focal$Coefficient_std_error * scale)
}

tau <- rbindlist(lapply(Sys.glob(file.path(RES, '*.results')), parse))
panel <- fread('data/5b-ldsc/sensitivity/FinnGen_phenotypes_final.csv')
tau <- merge(tau, panel[, .(`FinnGen ID`, Phenotype, Category)], by = 'FinnGen ID')
fwrite(tau, file.path(OUT, 'sldsc_tau.tsv'), sep = '\t')

ivw <- function(d) {
  w <- 1 / d$tau_star_se^2
  mu <- sum(w * d$tau_star) / sum(w); se <- 1 / sqrt(sum(w))
  data.table(n = nrow(d), tau_star = mu, se = se,
             p = 2 * pnorm(-abs(mu / se)))
}
meta <- tau[, ivw(.SD), by = .(hub, Category)]
fwrite(meta, file.path(OUT, 'sldsc_meta.tsv'), sep = '\t')
print(meta[, .(hub, Category, n, tau_star = round(tau_star, 4),
               se = round(se, 4), p = signif(p, 2))])
