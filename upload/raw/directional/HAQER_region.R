rm(list=ls())

library(readxl)
library(data.table)

HAQER <- read_excel(
  path = "raw/directional/1-s2.0-S0092867422013587-mmc1.xlsx",
  sheet = "HaqerBed"
)

HAQER = as.data.frame(HAQER)
HAQER_region = HAQER[,c('CHROM','START','END')]
colnames(HAQER_region) = c('chrom','start','end')

write.table(HAQER_region,'raw/directional/HAQER_region.txt',row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
