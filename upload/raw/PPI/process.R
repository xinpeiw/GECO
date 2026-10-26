# Inputs not included in this repository (place in raw/PPI/; md5, bytes):
#   HuRI.tsv                                   faa7b1cd2138f75904d74a0863df0296    1,681,536  https://www.interactome-atlas.org/data/HuRI.tsv
#   BioPlex_293T_Network_10K_Dec_2019.tsv      3e4d18a9da218583ceff95a0bd1dada3   12,866,060  https://bioplex.hms.harvard.edu/data/ (BioPlex 3.0)
#   BioPlex_HCT116_Network_5.5K_Dec_2019.tsv   25877d33ce9c551df82c0b87c3c965ab    7,714,364  https://bioplex.hms.harvard.edu/data/ (BioPlex 3.0)
#   BioPlex_3p0_293T_baitList.tsv              513af1363975b1f5a6dc13c15594c00f      209,203  https://bioplex.hms.harvard.edu/data/
#   BioPlex_3p0_HCT116_baitList.tsv            ee55d56be4b5d10905047c1997cd76fc      113,703  https://bioplex.hms.harvard.edu/data/
#   HIPPIE-current.mitab.txt                   32a0d1414acc0697f760af9276acea3a  257,038,149  HIPPIE v2.4, https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php
#   9606.protein.links.v12.0.txt.gz            5b7c91d02926e723ad30e81ab9d009a4   83,164,437  STRING v12.0, https://stringdb-downloads.org/download/protein.links.v12.0/
#   9606.protein.physical.links.v12.0.txt.gz   112d22188a72885c00e7baec9dc7d74d    8,954,065  STRING v12.0, https://stringdb-downloads.org/download/protein.physical.links.v12.0/
#   9606.protein.aliases.v12.0.txt.gz          09a0a186e4f8edab3e762cb0562f221e   19,777,800  STRING v12.0, https://stringdb-downloads.org/download/protein.aliases.v12.0/
# HuRI_Luck2020_SupplementaryTable1.txt is Supplementary Table 1 of Luck et al. 2020 (Nature 580:402).
suppressMessages({ library(data.table); library(EnsDb.Hsapiens.v86); library(AnnotationFilter) })

IN  <- 'raw/PPI'
OUT <- 'data/6b-edge/PPI'
geco <- readLines('data/onek1k/GECO_universe.txt')
wt <- function(x, f) write.table(x, file.path(OUT, f), row.names = FALSE, col.names = TRUE, sep = ' ', quote = FALSE)

pc <- unique(as.data.table(as.data.frame(genes(EnsDb.Hsapiens.v86, filter = GeneBiotypeFilter('protein_coding'))))$gene_id)
universe_pc <- intersect(geco, pc)

hg  <- fread('raw/hgnc_complete_set.txt', na.strings = c('', 'NA'), showProgress = FALSE,
             select = c('entrez_id', 'ensembl_gene_id', 'symbol'))
map <- unique(hg[!is.na(entrez_id) & !is.na(ensembl_gene_id), .(entrez_id, ensembl_gene_id)])
map_1to1 <- map[, if (.N == 1) .SD, by = entrez_id]

links <- fread(file.path(IN, 'HuRI.tsv'), header = FALSE, col.names = c('ensg1', 'ensg2'))
space <- fread(file.path(IN, 'HuRI_Luck2020_SupplementaryTable1.txt'))
universe <- intersect(geco, space[which(space$in_space_3 == 1), ]$Ensembl_gene_id)
links <- links[!is.na(ensg1) & !is.na(ensg2) & ensg1 != ensg2]
links[, c('ensg1', 'ensg2') := list(pmin(ensg1, ensg2), pmax(ensg1, ensg2))]
links <- unique(links, by = c('ensg1', 'ensg2'))
links <- links[ensg1 %in% universe & ensg2 %in% universe]
wt(links, 'HuRI_links.txt')
wt(universe, 'HuRI_universe.txt')

for (name in c('293T', 'HCT116')) {
  net  <- fread(file.path(IN, sprintf('BioPlex_%s_Network_%s_Dec_2019.tsv', name, c(`293T` = '10K', HCT116 = '5.5K')[[name]])))[, .(GeneA, GeneB)]
  bait <- fread(file.path(IN, sprintf('BioPlex_3p0_%s_baitList.tsv', name)))
  net[, GeneA_ensg := map_1to1[.SD, on = .(entrez_id = GeneA), x.ensembl_gene_id]]
  net[, GeneB_ensg := map_1to1[.SD, on = .(entrez_id = GeneB), x.ensembl_gene_id]]
  bait[, ensg := map_1to1[.SD, on = .(entrez_id = GeneID), x.ensembl_gene_id]]
  links <- net[!is.na(GeneA_ensg) & !is.na(GeneB_ensg), .(GeneA = GeneA_ensg, GeneB = GeneB_ensg)]
  universe <- intersect(bait$ensg, geco)
  links <- links[GeneA %in% universe & GeneB %in% universe & GeneA != GeneB]
  links[, `:=`(gene1 = pmin(GeneA, GeneB), gene2 = pmax(GeneA, GeneB))]
  links <- unique(links[, .(gene1, gene2)])
  wt(links, sprintf('BioPlex_%s_links.txt', name))
  wt(universe, sprintf('BioPlex_%s_universe.txt', name))
}

hp <- fread(file.path(IN, 'HIPPIE-current.mitab.txt'), select = c('ID Interactor A', 'ID Interactor B', 'Confidence Value'))
setnames(hp, c('GeneA', 'GeneB', 'Value'))
hp[, `:=`(GeneA = sub('.*:', '', GeneA), GeneB = sub('.*:', '', GeneB))]
hp <- hp[GeneA != '-' & GeneB != '-']
hp[, `:=`(GeneA = as.numeric(GeneA), GeneB = as.numeric(GeneB))]
hp[, GeneA_ensg := map_1to1[.SD, on = .(entrez_id = GeneA), x.ensembl_gene_id]]
hp[, GeneB_ensg := map_1to1[.SD, on = .(entrez_id = GeneB), x.ensembl_gene_id]]
for (lv in list(c('high', 0.73), c('medium', 0.63))) {
  links <- hp[Value > as.numeric(lv[2]) & !is.na(GeneA_ensg) & !is.na(GeneB_ensg), .(GeneA = GeneA_ensg, GeneB = GeneB_ensg)]
  links <- links[GeneA %in% universe_pc & GeneB %in% universe_pc & GeneA != GeneB]
  links[, `:=`(gene1 = pmin(GeneA, GeneB), gene2 = pmax(GeneA, GeneB))]
  wt(unique(links[, .(gene1, gene2)]), sprintf('HIPPIE_%s_links.txt', lv[1]))
}
wt(universe_pc, 'protein_coding_universe.txt')

aliases <- fread(file.path(IN, '9606.protein.aliases.v12.0.txt.gz'))
setnames(aliases, '#string_protein_id', 'string_protein_id')
id_map <- unique(aliases[source == 'Ensembl_gene', .(string_protein_id, ensg = alias)])
STR <- c(STRING_links = '9606.protein.links.v12.0.txt.gz', STRING_physical_links = '9606.protein.physical.links.v12.0.txt.gz')
for (nm in names(STR)) {
  l <- fread(file.path(IN, STR[[nm]]))
  l <- id_map[l, on = .(string_protein_id = protein1)]; setnames(l, 'ensg', 'ensg1')
  l <- id_map[l, on = .(string_protein_id = protein2)]; setnames(l, 'ensg', 'ensg2')
  l <- l[, c('ensg1', 'ensg2', 'combined_score')]
  l <- l[ensg1 != ensg2]
  l[, c('ensg1', 'ensg2') := list(pmin(ensg1, ensg2), pmax(ensg1, ensg2))]
  l <- unique(l, by = c('ensg1', 'ensg2'))
  l <- l[combined_score >= 700 & !is.na(ensg1) & !is.na(ensg2)]
  l <- l[ensg1 %in% universe_pc & ensg2 %in% universe_pc]
  wt(l, paste0(nm, '.txt'))
}
wt(universe_pc, 'STRING_universe.txt')
