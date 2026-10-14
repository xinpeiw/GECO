args <- commandArgs(trailingOnly=TRUE)
task_id <- as.numeric(args[1])

n_job = 255
chunk_id = task_id

library(data.table)
library(fst)
library(rhdf5)

dt_base <- as.data.table(
  read_fst('data_allcell.fst', columns = c("individual", "cell_label")))

individuals = unique(dt_base$individual)
cell_types = unique(dt_base$cell_label)
n_ind = length(individuals)
n_ct = length(cell_types)

ind_codes <- setNames(seq_len(n_ind), individuals)
ct_codes  <- setNames(seq_len(n_ct),  cell_types)

dt_base[, individual := ind_codes[as.character(individual)]]
dt_base[, cell_label := ct_codes[as.character(cell_label)]]

n_gene <- 11526

n = choose(n_gene, 2) / n_job
start = (chunk_id - 1) * n + 1
end = chunk_id * n

OPms = fread("OPms_allgene.txt")
all_genes = OPms$g

load("gene_pairs.Rdata")   ## gene_pairs: 2 x M

dir.create('VS_data', showWarnings = FALSE)

gene_to_idx <- setNames(seq_along(all_genes), all_genes)

if (chunk_id == 1) {
  diag_file <- file.path('VS_data',"VS_diag.h5")
  if (file.exists(diag_file)) file.remove(diag_file)

  h5createFile(diag_file)
  h5createGroup(diag_file, "VS")

  ## VS/var: [n_ind, n_ct, n_gene]
  h5createDataset(
    file = diag_file,
    dataset = "VS/var",
    dims = c(n_ind, n_ct, n_gene),
    storage.mode = "double",
    chunk = c(min(n_ind, 10), min(n_ct, 10), 1)
  )

  for (g in all_genes) {
    dtg <- as.data.table(
      read_fst('data_allcell.fst', columns = g)
    )
    dt <- cbind(dt_base, dtg)

    res <- dt[
      ,
      .(var = if (.N > 1) var(get(g)) / .N else NA_real_),
      by = .(individual, cell_label)
    ]

    slice <- matrix(NA_real_, nrow = n_ind, ncol = n_ct)
    slice[cbind(res$individual, res$cell_label)] <- res$var

    idx <- unname(gene_to_idx[[g]])

    h5write(
      obj = slice,
      file = diag_file,
      name = "VS/var",
      index = list(NULL, NULL, idx)
    )
  }
}

cov_file <- sprintf(paste0("VS_data", "/VS_cov_%03d.h5"), chunk_id)
if (file.exists(cov_file)) file.remove(cov_file)

h5createFile(cov_file)
h5createGroup(cov_file, "VS")

n_pair_task <- end - start + 1

## VS/cov: [n_ind, n_ct, n_pair_task]
h5createDataset(
  file = cov_file,
  dataset = "VS/cov",
  dims = c(n_ind, n_ct, n_pair_task),
  storage.mode = "double",
  chunk = c(min(n_ind, 10), min(n_ct, 10), 1)
)

pair_names <- character(n_pair_task)

for (j in start:end) {
  pair_local_idx <- j - start + 1

  sel <- as.character(gene_pairs[, j])
  g1 <- sel[1]
  g2 <- sel[2]

  pair_names[pair_local_idx] <- paste0(g1, "_", g2)

  dt2 <- as.data.table(
    read_fst('data_allcell.fst', columns = sel)
  )
  dt <- cbind(dt_base, dt2)

  res <- dt[
    ,
    .(cov12 = if (.N > 1) cov(get(g1), get(g2)) / .N else NA_real_),
    by = .(individual, cell_label)
  ]

  slice <- matrix(NA_real_, nrow = n_ind, ncol = n_ct)
  slice[cbind(res$individual, res$cell_label)] <- res$cov12

  h5write(
    obj = slice,
    file = cov_file,
    name = "VS/cov",
    index = list(NULL, NULL, pair_local_idx)
  )
}
