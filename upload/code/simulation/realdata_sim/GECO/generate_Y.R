
args <- commandArgs(trailingOnly = TRUE)
print(args)

rep_id  <- as.integer(args[1])
rep_dir <- args[2]

seed <- 1e6 + rep_id * 1000 
set.seed(seed)

N = 928
G = 11526
C = 7

library(data.table)
library(rhdf5)
library(fst)

cts = c('CD4ET','NK', 'CD4NC','CD8ET','BIN','CD8NC','BMem')
pre_data = as.data.table(read_fst("data_allcell.fst", columns = c("individual", "cell_label")))
individuals = unique(pre_data$individual)

# alpha
R = readRDS('./cov_matrix/chol_alpha.rds')
Z <- matrix(rnorm(N * G), nrow = N, ncol = G)
alpha <- Z %*% R # 928*11526, colnames是基因名
Alpha <- array(rep(alpha, C), dim = c(nrow(alpha), ncol(alpha), C))


# Gamma
Gamma <- array(NA, dim = c(N, G, C))

for (c in 1:7) {

	R = readRDS(paste0('./cov_matrix/chol_',cts[c],'.rds'))
	Z <- matrix(rnorm(N * G), nrow = N, ncol = G)
	gamma_c <- Z %*% R # N*G
	Gamma[,,c] = gamma_c

}

# Beta
beta = fread('./beta.txt')       # 11526*8
beta = as.matrix(beta[,2:8])
Beta <- array(rep(beta, each = N), dim = c(N, nrow(beta), ncol(beta)))

# DELTA
Delta <- array(NA, dim = c(N, G, C))
for (c in 1:7) {

	for (i in 1:N) {
		R = readRDS(paste0('./cholVS/cholVS_',
			paste0('ind',i), paste0('_ct',c), '.rds'))
		Z <- matrix(rnorm(1 * 11526), nrow = 1, ncol = 11526)
		delta_ic <- Z %*% R
		Delta[i,,c] = delta_ic
	}

}

dimnames(Delta) <- list(individuals, colnames(alpha), cts)
dimnames(Gamma) <- list(individuals, colnames(alpha), cts)
dimnames(Beta) <- list(individuals, colnames(alpha), cts)
dimnames(Alpha) <- list(individuals, colnames(alpha), cts)

Y = Beta + Alpha + Gamma + Delta

h5file <- file.path(rep_dir, "data", 'Y_data.h5')
if (file.exists(h5file)) file.remove(h5file)
h5createFile(h5file)
h5createDataset(
  file       = h5file,
  dataset    = "Y",
  dims       = c(N, G, C),
  chunk      = c(N, 1, C),   # one gene per chunk
  level      = 7,                   # compression level
  storage.mode = "double"
)

genes_all = colnames(alpha)

for (j in seq_along(genes_all)) {
  g <- genes_all[j]

  slice = Y[,g,]
  
  h5write(slice, h5file, paste0("Y_", g))
}

h5closeAll()
