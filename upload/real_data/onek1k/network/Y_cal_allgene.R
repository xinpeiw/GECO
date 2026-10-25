library(rhdf5)
library(data.table)
library(fst)

ctmm = fread("he_wald_free_V_p.txt")
genes_all = ctmm$gene

pre_data = as.data.table(read_fst("data_allcell.fst", columns = c("individual", "cell_label")))
individuals = unique(pre_data$individual)
cell_types = unique(pre_data$cell_label)
n_ind = length(individuals)
n_ct = length(cell_types)

ind_codes <- setNames(seq_along(individuals), individuals)
ct_codes  <- setNames(seq_along(cell_types),  cell_types)

h5file <- "Y_data.h5"
if (file.exists(h5file)) file.remove(h5file)
h5createFile(h5file)
h5createDataset(
  file       = h5file,
  dataset    = "Y",
  dims       = c(n_ind, length(genes_all), n_ct),
  chunk      = c(n_ind, 1, n_ct),   # one gene per chunk
  level      = 7,                   # compression level
  storage.mode = "double"
)
for (j in seq_along(genes_all)) {
  g <- genes_all[j]

  dt <- as.data.table(
    read_fst("data_allcell.fst", columns = c("individual", "cell_label", g))
  )

  dt[, individual := ind_codes[as.character(individual)] ]
  dt[, cell_label :=  ct_codes[as.character(cell_label)] ]

  results <- dt[, .(mu = mean(get(g))), by=.(individual, cell_label)]

  slice <- matrix(NA_real_, nrow = n_ind, ncol = n_ct)
  slice[cbind(results$individual, results$cell_label)] <- results$mu
  
  h5write(slice, "Y_data.h5", paste0("Y_", g))
}

h5closeAll()

