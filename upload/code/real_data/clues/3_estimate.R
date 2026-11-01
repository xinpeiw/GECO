args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

library(fst)
library(gtools)
library(MASS)
library(Matrix)
library(data.table)
library(rhdf5)

n_job = 241
chunk_id = task_id
group = 'all'

source("../onek1k/basic_func.R")
source(paste0("fixed_cov_",group,'.R')) 
source("../onek1k/GeCo.R")

load('gene_pairs.RData')

all_genes = fread('genelist11087.txt', header = FALSE)$V1
n_gene <- length(all_genes)

n = choose(n_gene, 2) / n_job
start = (chunk_id - 1) * n + 1
end = chunk_id * n


get_VS <- function(gene1, gene2, pair_local_idx) {

  cov12 <- h5read(sprintf(paste0("VS_data_", group, "/VS_cov_%03d.h5"), chunk_id), "VS/cov", index = list(NULL, NULL, pair_local_idx))
  dim(cov12) <- dim(cov12)[1:2]

  idx1 <- unname(gene_to_idx[[gene1]])
  idx2 <- unname(gene_to_idx[[gene2]])

  var1 <- h5read(paste0("VS_data_", group, "/VS_diag.h5"), "VS/var", index = list(NULL, NULL, idx1))
  var2 <- h5read(paste0("VS_data_", group, "/VS_diag.h5"), "VS/var", index = list(NULL, NULL, idx2))
  dim(var1) <- dim(var1)[1:2]
  dim(var2) <- dim(var2)[1:2]

  out <- array(NA_real_, dim = c(n_ind, n_ct, 2, 2))
  out[, , 1, 1] <- var1
  out[, , 1, 2] <- cov12
  out[, , 2, 1] <- cov12
  out[, , 2, 2] <- var2

  out
}


Ngp = ncol(gene_pairs)
opms = as.data.table(fread(paste0("OPms_allgene_", group, ".txt")))
colnames(opms)[1] = 'gene'
setkey(opms, gene)  

#genes_all <- opms$gene
genes_all = fread('genelist11087.txt', header = FALSE)$V1
gene_to_idx <- setNames(seq_along(genes_all), genes_all)
n_ind = 243
n_ct = 7

lt_alpha <- lower.tri(matrix(1, 2, 2), diag = TRUE)
idxs <- split(seq_len(14),rep(seq_len(7), each = 2))

res_cov = matrix(NA, nrow = n, ncol = 26)
res_se = matrix(NA, nrow = n, ncol = 26)


for (gp in start:end) {

	g1 = gene_pairs[1,gp]
	g2 = gene_pairs[2,gp]

	if (gp %% n != 0) {index = gp %% n} else (index = n)

	pair_local_idx = ((gp - 1) %% n) + 1

	res = try({

		VS <- get_VS(gene1 = g1, gene2 = g2, pair_local_idx = pair_local_idx)
		Y1 <- h5read(paste0('Y_data_',group,'.h5'), paste0("Y_", g1))
		Y2 <- h5read(paste0('Y_data_',group,'.h5'), paste0("Y_", g2))
		# scale
		ms1 = opms[g1]; m1 = ms1$OPmean; s1 = ms1$OPsd
		Y1 = (Y1 - m1) / s1
		ms2 = opms[g2]; m2 = ms2$OPmean; s2 = ms2$OPsd
		Y2 = (Y2 - m2) / s2
		VS[,,1,1] = VS[,,1,1] / s1^2
		VS[,,2,2] = VS[,,2,2] / s2^2
		VS[,,1,2] = VS[,,1,2] / s1 / s2
		VS[,,2,1] = VS[,,2,1] / s1 / s2 

		Y <- array(NA, dim = c(nrow(Y1), 2, ncol(Y1)))
		Y[,1,] <- Y1
		Y[,2,] <- Y2

		geco = GeCo(Y = Y, VS = VS, fixed_cov = fixed_cov)
		}, silent = TRUE)

	if (inherits(res, "try-error")) {
    	next 
  	}

	res_cov[index,1:2] = res_se[index, 1:2] = c(g1,g2)

	Valpha = geco$V_alpha; Vgc = geco$V_GC
	res_cov[index,3:5] = as.numeric(Valpha[lt_alpha])
	lower.list <- lapply(idxs, function(id) {
      subm <- Vgc[id, id]
      subm[lower.tri(subm, diag = TRUE)]
    })
    res_cov[index, 6:26] = as.numeric(unlist(lower.list))


    SEalpha = geco$SE_alpha; SEgc = geco$SE_GC
	res_se[index,3:5] = as.numeric(SEalpha[lt_alpha])
	lower.list <- lapply(idxs, function(id) {
      subm <- SEgc[id, id]
      subm[lower.tri(subm, diag = TRUE)]
    })
    res_se[index, 6:26] = as.numeric(unlist(lower.list))
}

res_cov = as.data.frame(res_cov)
res_se = as.data.frame(res_se)

names(res_cov) = c('g1', 'g2', 
	paste0('alpha',c('11','12','22')),
	paste0(paste0('g', c('11','12','22')),c(rep('_c1',3),rep('_c2',3),rep('_c3',3),rep('_c4',3),rep('_c5',3),rep('_c6',3),rep('_c7',3))))

names(res_se) = c('g1', 'g2', 
	paste0('alpha',c('11','12','22')),
	paste0(paste0('g', c('11','12','22')),c(rep('_c1',3),rep('_c2',3),rep('_c3',3),rep('_c4',3),rep('_c5',3),rep('_c6',3),rep('_c7',3))))


write.table(res_cov, sprintf(paste0('res_',group,'/cov_%03d.txt'), chunk_id),row.names = FALSE, quote = FALSE, sep = ' ')
write.table(res_se, sprintf(paste0('res_',group,'/se_%03d.txt'), chunk_id),row.names = FALSE, quote = FALSE, sep = ' ')


