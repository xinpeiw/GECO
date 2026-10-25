rm(list = ls())

library(readxl)
library(data.table)

geco = readLines('data/onek1k/GECO_universe.txt')
clues = readLines('data/clues/allgenes_CLUES.txt')

ST10 <- as.data.frame(read_excel(
  path = 'raw/onek1k_egene/science.abf3041_tables_s6_to_s19.xlsx',
  sheet = "Table.S10", skip = 2 ))
cis = unique(ST10$`Gene Ensembl ID`)
cis_geco = cis[which(cis %in% geco)]
write.table(cis_geco, 'data/5a-cis_molGenes/cis_eGene/cis_onek1k.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
cis_clues = cis[which(cis %in% clues)]
write.table(cis_clues, 'data/others/molQTL_clues/cis/cis_onek1k.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

hg   = fread('raw/hgnc_complete_set.txt', na.strings=c('','NA'), showProgress=FALSE,
             select=c('symbol','ensembl_gene_id','prev_symbol'))
main = setNames(hg$ensembl_gene_id, hg$symbol)
prev = { d=hg[!is.na(prev_symbol)&nzchar(prev_symbol)]; r=strsplit(d[['prev_symbol']],'\\|'); setNames(rep(d$ensembl_gene_id,lengths(r)),unlist(r)) }
prev = prev[!duplicated(names(prev)) & !duplicated(names(prev), fromLast=TRUE)]
to_ensg = function(x){ o=ifelse(grepl('^ENSG[0-9]+$',x), x, unname(main[x])); m=is.na(o)|!nzchar(o); o[m]=unname(prev[x[m]]); o }

ST16 <- as.data.frame(read_excel(
  path = 'raw/onek1k_egene/science.abf3041_tables_s6_to_s19.xlsx',
  sheet = "Table.S16", skip = 2 ))
trans = unique(ST16$`trans eGene`)
ST10 <- as.data.frame(read_excel(
  path = 'raw/onek1k_egene/science.abf3041_tables_s6_to_s19.xlsx',
  sheet = "Table.S10", skip = 2 ))
ST10 = unique(ST10[,c('Gene ID','Gene Ensembl ID')])
st10_map = setNames(ST10$`Gene Ensembl ID`, ST10$`Gene ID`)
st10_map = st10_map[!duplicated(names(st10_map)) & !duplicated(names(st10_map), fromLast=TRUE)]
to_ensg_final = function(x) {
  e = to_ensg(x)
  m = is.na(e) | !nzchar(e)
  e[m] = unname(st10_map[x[m]])
  e
}
trans_ensg = to_ensg_final(trans)

trans_geco = intersect(trans_ensg, geco)
write.table(trans_geco, 'data/6a-trans_molGenes/trans_eGene/trans_onek1k.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

pairs = unique(ST16[,c('cis eGene','trans eGene')])
genes = union(pairs$`cis eGene`, pairs$`trans eGene`)
genes_ensg = data.frame(genes = genes, ensg = to_ensg_final(genes), stringsAsFactors = FALSE)
gene_map <- setNames(genes_ensg$ensg, genes_ensg$genes)
pairs_ensg <- pairs
pairs_ensg[["cis eGene"]] <- unname(gene_map[pairs[["cis eGene"]]])
pairs_ensg[["trans eGene"]] <- unname(gene_map[pairs[["trans eGene"]]])

pairs_ensg_geco = pairs_ensg[which(!is.na(pairs_ensg$`cis eGene`) & !is.na(pairs_ensg$`trans eGene`) &
                              pairs_ensg$`cis eGene` %in% geco & pairs_ensg$`trans eGene` %in% geco),]
colnames(pairs_ensg_geco) = c('cis_eGene', 'trans_eGene')
pairs_ensg_geco = unique(pairs_ensg_geco)
for (o in c('data/6b-edge/cis_trans_pair/onek1k.txt', 'data/6c-hub_position/cis_trans_pair/onek1k.txt')) write.table(pairs_ensg_geco, o,
            row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)

ST8 = as.data.frame(read_excel(
  path = 'raw/onek1k_egene/science.abf3041_tables_s6_to_s19.xlsx',
  sheet = "Table.S8", skip = 2 ))
onek1k_universe = unique(trimws(na.omit(unlist(ST8, use.names = FALSE))))
onek1k_universe <- onek1k_universe[onek1k_universe != ""]
clues_universe = intersect(to_ensg_final(onek1k_universe), clues)
write.table(clues_universe, 'data/others/molQTL_clues/onek1k_universe.txt',
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
geco_universe = intersect(to_ensg_final(onek1k_universe), geco)
for (o in c('data/5a-cis_molGenes/onek1k_universe.txt', 'data/6a-trans_molGenes/onek1k_universe.txt', 'data/6b-edge/cis_trans_pair/onek1k_universe.txt', 'data/6c-hub_position/cis_trans_pair/onek1k_universe.txt')) write.table(geco_universe, o,
            row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
