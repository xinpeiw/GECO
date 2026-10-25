rm(list=ls())

library(data.table)

geco_universe = readLines('data/onek1k/GECO_universe.txt')
TF_set = list('Low specificity TF' = intersect(fread('raw/TF_HPA/low_specificity.tsv')$Ensembl, geco_universe),
              'T cell TF' = intersect(fread('raw/TF_HPA/elevated_T-cells_Lineage.tsv')$Ensembl, geco_universe),
              'NK cell TF' = intersect(fread('raw/TF_HPA/elevated_NK-cells_Lineage.tsv')$Ensembl, geco_universe),
              'B cell TF' = intersect(fread('raw/TF_HPA/elevated_B-cells_Lineage.tsv')$Ensembl, geco_universe))

for (n in names(TF_set)) {
  write.table(TF_set[[n]], file.path('data/4a-TF_HPA/TFset', paste0(n, '.txt')),
              row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
}
