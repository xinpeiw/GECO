args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

library(data.table)

group = 'all'
rep_dir = paste0('./res_',group)

genes = fread('genelist11087.txt', header = FALSE)$V1

for (name in c('cov','se')) {

	input = file.path(rep_dir, paste0(name, '_', group,'.txt'))

	data = as.data.frame(fread(input))

	dim(data)[1] == choose(length(genes),2)
	length(unique(data$g1)) == length(genes)-1
	length(unique(data$g2)) == length(genes)-1

	# prepare
	b = data[,c('g1', 'g2', 'alpha12', paste0('g12_c',1:7))]
	n     <- length(genes)
	idx1  <- match(b$g1, genes)
	idx2  <- match(b$g2, genes)

	# re-organize data
	col_list = c('alpha12', paste0('g12_c',1:7))
	for (c in 1:length(col_list)) {

		col = col_list[c]

		cov_matrix <- matrix(NA_real_, n, n, dimnames = list(genes, genes))
		cov_matrix[cbind(idx1, idx2)] <- b[[ col ]]
		cov_matrix[cbind(idx2, idx1)] <- b[[ col ]]

		if (c==1) {
			col1 = paste0(strsplit(col,'1')[[1]][1],'11')
			col2 = paste0(strsplit(col,'1')[[1]][1],'22')
		} else {
			col1 = paste0('g11_',strsplit(col,'_')[[1]][2])
			col2 = paste0('g22_',strsplit(col,'_')[[1]][2])
		}

		a1 = data[,c('g1',col1)]; a2 = data[,c('g2',col2)]
		names(a1) = c('g','var'); names(a2) = c('g','var')
		a1 = a1[!duplicated(a1$g),]; a2 = a2[!duplicated(a2$g),]
		d = rbind(a1,a2)
		d = d[!duplicated(d$g),]

		idx <- match(genes, d$g)  
		diag(cov_matrix) <- d$var[idx]

		outpath = file.path(rep_dir, paste0(name, '_', col , '.rds'))
		saveRDS(cov_matrix,  file = outpath)

	}

}



