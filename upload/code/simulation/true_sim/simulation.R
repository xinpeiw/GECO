

source("../../real_data/onek1k/basic_func.R")
source("../../real_data/onek1k/GeCo.R")
source("data_generate.R")

sum = read.table('variable_sum.txt', header = TRUE)

P = data.table::fread('P')
pre_data = data.table::as.data.table(fst::read_fst("data_allcell.fst", columns = c("individual")))
individuals = unique(pre_data$individual)
P = P[which(P$ind %in% individuals),c('CD4 NC','CD4 ET','CD8 NC','CD8 ET','NK','B IN','B Mem')]
p = apply(P, 2, median)
p = p / sum(p)

cts = c('CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')

beta = matrix(c(sum[which(sum$variable == 'beta'),]$median,
	sum[which(sum$variable == 'beta'),]$X3rd
	), nrow = 2, byrow = TRUE)

V_generate = function(rho_alpha, rho_omega, rho_delta) {

	# V_alpha
	alpha11 = sum[which(sum$variable == 'V_alpha' & sum$subgroup == 'var'),]$median
	alpha22 = sum[which(sum$variable == 'V_alpha' & sum$subgroup == 'var'),]$X3rd
	V_alpha = matrix(c(alpha11, rho_alpha*sqrt(alpha11 * alpha22), rho_alpha*sqrt(alpha11 * alpha22), alpha22), ncol = 2)

	# gamma
	V_GC = array(NA, dim = c(C, G, G))
	for (c in 1:C) {
		d = sum[which(sum$variable == paste0('V_', cts[c])),]
		d11 = d[which(d$subgroup == 'var'),]$median
		d22 = d[which(d$subgroup == 'var'),]$X3rd
		V_GC[c,,] = matrix(c(d11,rho_omega[c]*sqrt(d11*d22),rho_omega[c]*sqrt(d11*d22), d22), nrow = 2)
	}

	# delta
	V_delta = array(NA, dim = c(C, G, G))
	for (c in 1:C) {
		d = sum[which(sum$subgroup == cts[c]),]
		d11 = d[which(d$variable == 'VS11'),]$median
		d22 = d[which(d$variable == 'VS11'),]$X3rd
		V_delta[c,,] = matrix(c(d11, rho_delta[c]*sqrt(d11*d22),rho_delta[c]*sqrt(d11*d22), d22), nrow = 2)
	}

	list('V_alpha' = V_alpha, 'V_GC' = V_GC, 'V_delta' = V_delta)
}

C = 7
G = 2
n_rep = 1000
lt_alpha <- lower.tri(matrix(1, 2, 2), diag = TRUE)
idxs <- split(seq_len(2*C),rep(seq_len(C), each = 2))

cols = c(paste0('Shared_',c('11','12','22')), paste0(cts[1],c('11','12','22')), paste0(cts[2],c('11','12','22')),
	paste0(cts[3],c('11','12','22')), paste0(cts[4],c('11','12','22')), paste0(cts[5],c('11','12','22')), 
	paste0(cts[6],c('11','12','22')), paste0(cts[7],c('11','12','22')))

simu = function(N, rho_alpha = NA, rho_delta = NA, rho_omega = NA) {

	if (sum(is.na(rho_delta)) > 0) { rho_delta = runif(C, -1, 1) }
	if (sum(is.na(rho_omega)) > 0) { rho_omega = runif(C, -1, 1) }
	if (is.na(rho_alpha)) { rho_alpha = runif(1, -1, 1) }

	Vs = V_generate(rho_alpha = rho_alpha, rho_omega = rho_omega, rho_delta = rho_delta)
	V_alpha = Vs$V_alpha; V_GC = Vs$V_GC; V_delta = Vs$V_delta

	res_cov = matrix(NA, nrow = n_rep, ncol = 3*(C+1))
	res_se = matrix(NA, nrow = n_rep, ncol = 3*(C+1))
	res_naive = matrix(NA, nrow = n_rep, ncol = 3*(C+1))

	for (i in 1:n_rep) {
		index = i
		Y = data_generate(N = N, G = G, C = C, beta = beta, V_alpha = V_alpha, V_GC = V_GC, V_delta = V_delta)

		VS <- array(rep(V_delta, times = N), dim = c(dim(V_delta), N))
		VS <- aperm(VS, c(4, 1, 2, 3))

		# geco
		geco = GeCo(Y = Y, VS = VS, fixed_cov = NA)
		Valpha = geco$V_alpha; Vgc = geco$V_GC
		res_cov[index,1:3] = as.numeric(Valpha[lt_alpha])
		lower.list <- lapply(idxs, function(id) { subm <- Vgc[id, id]; subm[lower.tri(subm, diag = TRUE)]})
		res_cov[index, 4:(3*(C+1))] = as.numeric(unlist(lower.list))
		SEalpha = geco$SE_alpha; SEgc = geco$SE_GC
		res_se[index,1:3] = as.numeric(SEalpha[lt_alpha])
		lower.list <- lapply(idxs, function(id) { subm <- SEgc[id, id]; subm[lower.tri(subm, diag = TRUE)]})
		res_se[index, 4:(3*(C+1))] = as.numeric(unlist(lower.list))

		# naive
		est_V = array(NA, dim = c(C, G, G))
		for (c in 1:C) { y = Y[,,c]; est_V[c,,] = cov(y)
			res_naive[index, (c*3+1):(c*3+3)] = as.numeric(est_V[c,,][lt_alpha])}

		op = apply(sweep(Y, 3, p, "*"), c(1, 2), sum)
		res_naive[index, 1:3] = cov(op)[lt_alpha]
	}

	colnames(res_cov) = cols
	colnames(res_se) = cols
	colnames(res_naive) = cols
	true = Vs$V_alpha[lt_alpha]
	for (c in 1:C) { true = c(true, Vs$V_GC[c,,][lt_alpha]) }
	names(true) = cols

	list('geco_cov' = res_cov, 'geco_se' = res_se, 'naive' = res_naive, 'true' = true)
}

alpha_seq = seq(-1, 1, 0.1)
geco = list()
naive = list()
true = list()
set.seed(111)
rho_omega = c(1,runif(C-1, -1, 1))
rho_delta = c(0,runif(C-1, -1, 1))
for (i in 1:length(alpha_seq)) {
	rho_alpha = alpha_seq[i]
 
	a = simu(N = 1000, rho_omega = rho_omega, rho_alpha = rho_alpha, rho_delta = rho_delta)
	geco[[i]] = a$geco_cov
	naive[[i]] = a$naive 
	true[[i]] = a$true
}
saveRDS(geco, "geco_CD4NC.rds")
saveRDS(naive, "naive_CD4NC.rds")
saveRDS(true, "true_CD4NC.rds")
