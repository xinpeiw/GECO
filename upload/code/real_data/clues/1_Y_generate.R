library(rhdf5)
library(data.table)
library(fst)

df = fread('ctp.final.gz')
genes_all = colnames(df)[3:ncol(df)]
write.table(genes_all, 'genelist11087.txt', row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)

genes_all = fread('genelist11087.txt', header = FALSE)$V1

gene_pairs = combn(genes_all,2)
save(gene_pairs, file = "gene_pairs.RData")

group = 'all'

pre_data = as.data.table(read_fst(paste0("data_allcell_", group, ".fst"), columns = c("ind_cov", "cg_cov")))
individuals = unique(pre_data$ind_cov)
cell_types = unique(pre_data$cg_cov)
n_ind = length(individuals)
n_ct = length(cell_types)

ind_codes <- setNames(seq_along(individuals), individuals)
ct_codes  <- setNames(seq_along(cell_types),  cell_types)

h5file <- paste0("Y_data_", group, ".h5")
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
      read_fst(paste0("data_allcell_", group, ".fst"), columns = c("ind_cov", "cg_cov", g))
    )

    dt[, ind_cov := ind_codes[as.character(ind_cov)] ]
    dt[, cg_cov :=  ct_codes[as.character(cg_cov)] ]

    results <- dt[, .(mu = mean(get(g))), by=.(ind_cov, cg_cov)]

    slice <- matrix(NA_real_, nrow = n_ind, ncol = n_ct)
    slice[cbind(results$ind_cov, results$cg_cov)] <- results$mu
  
    h5write(slice, paste0("Y_data_", group, ".h5"), paste0("Y_", g))
}

h5closeAll()

