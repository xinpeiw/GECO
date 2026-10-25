rm(list=ls())

GECO = readRDS("data/2a-variance_decomposition/geco_CD4NC.rds")
naive = readRDS("data/2a-variance_decomposition/naive_CD4NC.rds")
true = readRDS("data/2a-variance_decomposition/true_CD4NC.rds")

df = data.frame('true_alpha' = NA,'true_cd4nc' = NA,
                'geco_alpha_mean' = NA,'geco_alpha_sd'=NA,
                'geco_cd4nc_mean' = NA,'geco_cd4nc_sd'=NA,
                'naive_alpha_mean'=NA,'naive_alpha_sd'=NA,
                'naive_cd4nc_mean'=NA,'naive_cd4nc_sd'=NA)

for (i in 1:21) {
	df[i, 'true_alpha'] <- as.numeric(true[[i]]['Shared_12'])
	df[i, 'true_cd4nc'] <- as.numeric(true[[i]]['CD4NC12'])

	df[i, 'geco_alpha_mean'] <- as.numeric(apply(GECO[[i]], 2, mean)['Shared_12'])
	df[i, 'geco_alpha_sd'] <- as.numeric(apply(GECO[[i]], 2, sd)['Shared_12'])
	df[i, 'geco_cd4nc_mean'] <- as.numeric(apply(GECO[[i]], 2, mean)['CD4NC12'])
	df[i, 'geco_cd4nc_sd'] <- as.numeric(apply(GECO[[i]], 2, sd)['CD4NC12'])

	df[i, 'naive_alpha_mean'] <- as.numeric(apply(naive[[i]], 2, mean)['Shared_12'])
	df[i, 'naive_alpha_sd'] <- as.numeric(apply(naive[[i]], 2, sd)['Shared_12'])
	df[i, 'naive_cd4nc_mean'] <- as.numeric(apply(naive[[i]], 2, mean)['CD4NC12'])
	df[i, 'naive_cd4nc_sd'] <- as.numeric(apply(naive[[i]], 2, sd)['CD4NC12'])
}

library(tidyr)
library(dplyr)

true_cd4nc <- df[1,2]
df <- df[,-2]
colnames(df) <- c("true",
                       "GECOshared_mean", "GECOshared_sd",
                       "GECOspecific_mean", "GECOspecific_sd",
                       "Bulk_mean", "Bulk_sd",
                       "Pseudobulk_mean", "Pseudobulk_sd")

df_long <- df %>%
  mutate(id = 1:n()) %>%
  pivot_longer(
    cols = c(GECOshared_mean, GECOshared_sd,
             GECOspecific_mean, GECOspecific_sd,
             Bulk_mean, Bulk_sd,
             Pseudobulk_mean, Pseudobulk_sd),
    names_to = c("method", ".value"),
    names_pattern = "(GECOshared|GECOspecific|Bulk|Pseudobulk)_(mean|sd)"
  )
df_long$method[df_long$method == "GECOshared"] <- "GECO Shared"
df_long$method[df_long$method == "GECOspecific"] <- "GECO Specific"

df_long$method = factor(df_long$method,
                        levels = c('GECO Shared','GECO Specific',
                                   'Bulk','Pseudobulk'))

outdir = 'results/2a-variance_decomposition'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
write.table(as.data.frame(df_long), file.path(outdir, 'variance_decomposition.txt'),
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
