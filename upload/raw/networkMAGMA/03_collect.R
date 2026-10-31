suppressMessages(library(data.table))
source("raw/networkMAGMA/networks.R")

OUT <- list(
  geco          = list(file = "data/6d-networkMAGMA/networkMAGMA_78traits_results.txt",                   nets = VARIANTS$geco),
  clues         = list(file = "data/7e_networkMAGMA/networkMAGMA_78traits_CLUES_mega_results.txt",        nets = VARIANTS$clues),
  perturbseq    = list(file = "data/6f-networkMAGMA_perturbseq/networkMAGMA_78traits_perturbseq_results.txt",
                       nets = c("pertseq_specific_overlap", "pertseq_shared_overlap", "pertseq_union")),
  perturbseq_ct = list(file = "data/others/networkMAGMA_78traits_perturbseq_each_ct_results.txt",         nets = VARIANTS$perturbseq_ct),
  glasso01      = list(file = "data/others/networkMAGMA_78traits_geco_glasso_01_results.txt",             nets = VARIANTS$glasso01))

for (v in names(OUT)) {
  fs <- list.files(file.path("raw/networkMAGMA", paste0("results_", v)), pattern = "[.]txt$", full.names = TRUE)
  res <- rbindlist(lapply(fs, fread))[network %in% OUT[[v]]$nets]
  stopifnot(nrow(res) == 78 * length(OUT[[v]]$nets))
  setorder(res, network, trait)
  write.table(res, OUT[[v]]$file, row.names = FALSE, col.names = TRUE, sep = " ", quote = FALSE)
}
