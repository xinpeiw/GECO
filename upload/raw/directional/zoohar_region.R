rm(list=ls())

library(readxl)
library(data.table)

zoohar <- read_excel(
  path = "raw/directional/science.abm1696_table_s1.xlsx",
  sheet = "zooHARs"
)

zoohar = as.data.frame(zoohar)
zoohar_positive = zoohar[which(zoohar$translated_category=='positive'),]

zoohar_region = zoohar[,c('chrom','start','end')]
zoohar_positive_region = zoohar_positive[,c('chrom','start','end')]

write.table(zoohar_region,'raw/directional/zoohar_region.txt',row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
write.table(zoohar_positive_region,'raw/directional/zoohar_positive_region.txt',row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
