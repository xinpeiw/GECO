pre_data = data.table::as.data.table(fst::read_fst("data_allcell.fst", columns = c("individual","cell_label")))


library(data.table)

setDT(pre_data)

cell_prop_wide <- pre_data[, .N, by = .(individual, cell_label)]
cell_prop_wide[, prop := N / sum(N), by = individual]

cell_prop_wide <- dcast(
  cell_prop_wide,
  individual ~ cell_label,
  value.var = "prop",
  fill = 0
)

colnames(cell_prop_wide)[1] = 'ind'

saveRDS(cell_prop_wide, 'P.rds')