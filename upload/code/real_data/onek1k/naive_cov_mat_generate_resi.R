library(data.table)

cts = c('Shared', 'CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK', 'BIN', 'BMem')
cov = data.frame(readRDS("./covariates.rds"))

OP = read.table('OPms_allgene.txt', header = TRUE)

for (c in 1:length(cts)) {

  print(c)
  if (c == 1) {
    data = readRDS('expression_op.rds')  # already scaled
  } else {
    data = readRDS(paste0('expression_', cts[c], '.rds'))
    idx <- match(colnames(data), OP$g)
    OP2 <- OP[idx, ]
    stopifnot(identical(colnames(data), OP2$g))

    data_scaled <- sweep(data, 2, OP2$OPmean, "-")
    data_scaled <- sweep(data_scaled, 2, OP2$OPsd, "/")

    data = data_scaled
  }

  expr_resid <- apply(data, 2, function(y) {
    fit <- lm(y ~ cov$age_group + cov$sex + cov$PC1)
    resid(fit)
  })

  data = expr_resid

  cov_m = cov(data)
  saveRDS(cov_m, paste0('adj_exp/naive_cov_mat_', cts[c], '.rds'))
}
