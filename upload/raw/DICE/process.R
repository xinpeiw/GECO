rm(list=ls())

library(data.table)
library(readxl)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

xlsx_file <- "raw/DICE/1-s2.0-S009286741831331X-mmc3.xlsx"
s3a <- as.data.table(read_excel(xlsx_file, sheet = "Table S3A", skip = 2))

genes = unique(sub('\\..*', '', s3a$`Ensembl ID`))

genes_ensg_geco = intersect(genes, geco)
write.table(genes_ensg_geco, 'data/5a-cis_molGenes/cis_eGene/DICE.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
genes_ensg_clues = intersect(genes, clues)
write.table(genes_ensg_clues, 'data/others/molQTL_clues/cis/cis_DICE.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

celltype_files <- list.files('raw/DICE/per_celltype', pattern = '\\.txt$', full.names = TRUE)
tested_union <- unique(unlist(lapply(celltype_files, readLines)))
tested_union <- sub('\\..*', '', tested_union)

universe_geco <- intersect(tested_union, geco)
write.table(universe_geco, 'data/5a-cis_molGenes/DICE_universe.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
universe_clues <- intersect(tested_union, clues)
write.table(universe_clues, 'data/others/molQTL_clues/DICE_universe.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
