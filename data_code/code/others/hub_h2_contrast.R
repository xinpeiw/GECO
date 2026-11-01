library(data.table)
library(quantreg)

B <- 10000

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
d <- fread("data/others/h2_long.tsv")
stopifnot(!anyNA(d$in_shared), !anyNA(d$in_specific))

boot <- function(h2, in_shared, in_specific) {
  stats <- function(i) {
    x <- h2[i]; s <- in_shared[i]; p <- in_specific[i]
    c(shared = median(x[s]), `non-shared` = median(x[!s]),
      specific = median(x[p]), `non-specific` = median(x[!p]))
  }
  obs <- stats(seq_along(h2))
  for (set in list(in_shared, in_specific)) {
    m <- rq(h2 ~ set, tau = 0.5)
    fit <- ifelse(set, median(h2[set]), median(h2[!set]))
    stopifnot(isTRUE(all.equal(sum(abs(h2 - fit)), sum(abs(resid(m))),
                               tolerance = 1e-8)))
  }
  reps <- replicate(B, stats(sample(seq_along(h2), replace = TRUE)))
  dsh <- reps["shared", ] - reps["non-shared", ]
  dsp <- reps["specific", ] - reps["non-specific", ]
  data.table(class = factor(names(obs), names(obs)),
             n = c(sum(in_shared), sum(!in_shared), sum(in_specific), sum(!in_specific)),
             median = as.vector(obs),
             se = apply(reps, 1, sd, na.rm = TRUE),
             lo = apply(reps, 1, quantile, 0.025, na.rm = TRUE),
             hi = apply(reps, 1, quantile, 0.975, na.rm = TRUE),
             bg_median = c(obs[["non-shared"]], NA, obs[["non-specific"]], NA),
             vs_bg = c(obs[["shared"]] - obs[["non-shared"]], NA,
                       obs[["specific"]] - obs[["non-specific"]], NA),
             vs_bg_se = c(sd(dsh, na.rm = TRUE), NA, sd(dsp, na.rm = TRUE), NA))
}

res <- d[, { set.seed(1); boot(h2, in_shared, in_specific) }, by = .(dataset, partition)]
res[, `:=`(z = vs_bg / vs_bg_se,
           dataset = factor(dataset, unique(dataset)),
           partition = factor(partition, unique(partition)))]
fwrite(res, file.path(outdir, "hub_h2_contrast.tsv"), sep = "\t")
