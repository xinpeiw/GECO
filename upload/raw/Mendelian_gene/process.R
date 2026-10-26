suppressMessages(library(data.table))

pairs <- fread('raw/Mendelian_gene/panelapp_green_gene_disease_pairs.csv')
U     <- readLines('data/onek1k/GECO_universe.txt')

imm <- sort(unique(pairs[immune == 1, ensembl_id]))
non <- sort(unique(pairs[immune == 0, ensembl_id]))

writeLines(imm, 'data/7c-Mendelian_gene/mendelian_immune_full.txt')
writeLines(non, 'data/7c-Mendelian_gene/mendelian_nonimmune_full.txt')
writeLines(intersect(imm, U), 'data/5d-Mendelian_gene/mendelian_immune_full_in_GECO_universe.txt')
writeLines(intersect(non, U), 'data/5d-Mendelian_gene/mendelian_nonimmune_full_in_GECO_universe.txt')
