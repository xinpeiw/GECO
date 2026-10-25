rm(list = ls())

suppressMessages({ library(data.table); library(readxl) })

geco = readLines('data/clues/allgenes_CLUES.txt')

g = fread('raw/stabilizing_selection/gnomad_v41_constraint.tsv')
g = g[canonical == TRUE | mane_select == TRUE]
g[, ensg := sub('\\..*', '', gene_id)]
loeuf = g[!is.na(lof.oe_ci.upper),
          .(loeuf = min(lof.oe_ci.upper, na.rm = TRUE)), by = ensg]
fwrite(loeuf, 'raw/stabilizing_selection/gnomad_loeuf_per_gene.tsv', sep = '\t')
loeuf = as.data.frame(loeuf)
loeuf = loeuf[which(loeuf$ensg %in% geco),]
write.table(loeuf, 'data/others/evolution_clues/loeuf.txt',
            row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)

chr = suppressMessages(read_excel(
  'raw/stabilizing_selection/science.abn3943_data_s1_to_s3.xlsx',
  sheet = 'Data S1'))
phylop = data.table(ensg = chr$geneID, phylop = chr$meanPhylop)[
  !is.na(ensg) & !is.na(phylop)]
phylop = phylop[, .( phylop = max(phylop, na.rm = TRUE)),
                by = ensg]
fwrite(phylop, 'raw/stabilizing_selection/christmas_phylop_per_gene.tsv', sep = '\t')
phylop = as.data.frame(phylop)
phylop = phylop[which(phylop$ensg %in% geco),]
write.table(phylop, 'data/others/evolution_clues/phylop.txt',
            row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)
