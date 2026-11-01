args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

ct_id <- (task_id - 1) %/% 100 + 1
rep_id <- (task_id - 1) %% 100 + 1

setwd(file.path('reps', sprintf("rep_%03d", rep_id), 'data'))

library(WGCNA)
library(writexl)

cts = c('OP','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')

module_summary = function(module_of_interest) {
	
	ME_name <- paste0("ME", module_of_interest)
	
	#kME
	kME <- cor(datExpr, MEs[, ME_name], use = "pairwise.complete.obs")
	kME_p <- corPvalueStudent(as.vector(kME), nSamples = nrow(datExpr))

	hub_kME <- data.frame(
  		Gene = colnames(datExpr),
  		Module = moduleColors,
  		kME = as.numeric(kME),
  		kME_p = as.numeric(kME_p)
	)

	hub_kME_mod <- hub_kME[hub_kME$Module == module_of_interest, ]
	hub_kME_mod <- hub_kME_mod[order(-abs(hub_kME_mod$kME)), ]
	
	#kWithin
	adjacency_mat <- adjacency(
  		datExpr,
  		power = softPower,
  		type = "signed",
  		corFnc = "bicor",
  		corOptions = "use = 'pairwise.complete.obs'"
	)

	IMConn <- intramodularConnectivity(adjacency_mat, moduleColors)

	hub_conn <- data.frame(
  		Gene = colnames(datExpr),
 	 	Module = moduleColors,
  		kWithin = IMConn$kWithin,
  		kOut = IMConn$kOut,
  		kTotal = IMConn$kTotal
	)

	hub_conn_mod <- hub_conn[hub_conn$Module == module_of_interest, ]
	hub_conn_mod <- hub_conn_mod[order(-hub_conn_mod$kWithin), ]
	
	# summary
	hub_summary <- merge(hub_kME, hub_conn, by = c("Gene", "Module"))

	hub_summary_mod <- hub_summary[hub_summary$Module == module_of_interest, ]
	hub_summary_mod <- hub_summary_mod[order(-abs(hub_summary_mod$kME),
                                         -hub_summary_mod$kWithin), ]
                                         
    hub_summary_mod

}



datExpr = readRDS(paste0('expression_', cts[ct_id], '.rds'))

# soft-thresholding power
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(
  	datExpr,
  	powerVector = powers,
  	networkType = "signed",
  	corFnc = "bicor",
  	corOptions = list(use = "pairwise.complete.obs"),
  	verbose = 5
)

softPower <- sft$powerEstimate
if (is.na(softPower)) { softPower = (sft$fitIndices)$Power[which.max((sft$fitIndices)$SFT.R.sq)]}


net <- blockwiseModules(
  	datExpr,
  	power = softPower,
  	TOMType = "signed",
  	networkType = "signed",
  	corType = "bicor",
  	maxBlockSize = ncol(datExpr),
  	minModuleSize = 30,
  	reassignThreshold = 0,
  	mergeCutHeight = 0.25,
  	numericLabels = FALSE,
  	pamRespectsDendro = FALSE,
  	saveTOMs = FALSE,
  	verbose = 3
)
moduleColors <- net$colors
table(moduleColors)

MEs <- orderMEs(moduleEigengenes(datExpr, colors = moduleColors)$eigengenes)

modules = names(table(moduleColors))
sum_df = list()
for (j in 1:length(modules)) {
	dt = module_summary(modules[j])
	sum_df[[j]] = dt
}
names(sum_df) = paste0('mod_', 1:length(modules))

write_xlsx(sum_df, path = paste0(cts[ct_id], '.xlsx'))


