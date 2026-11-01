library(readxl)
library(dplyr)
library(purrr)
library(openxlsx)

base_dir <- "reps"
rep_ids <- 1:100

key_col   <- "Gene"  
cts = c('OP', 'CD4NC', 'CD4ET', 'CD8NC', 'CD8ET', 'NK', 'BIN','BMem')

read_one_rep <- function(file, rep_id) {

  sheet_names <- excel_sheets(file)
  if (length(sheet_names) == 0) { stop("No sheets found in: ", file) }
  dat <- map_dfr(sheet_names, ~ as.data.frame(read_excel(file, sheet = .x)) )
  dat <- dat %>% select(all_of(key_col), all_of(value_col)) %>% filter(!is.na(.data[[key_col]]))
  names(dat)[2] <- paste0(value_col, "_rep", sprintf("%03d", rep_id))

  dat
}


for (value_col in c("Module", "kME")) for (c in 1:length(cts)) {
  file_paths <- file.path(base_dir, paste0('rep_',sprintf("%03d", rep_ids)), 'data', paste0(cts[c], '.xlsx'))
  output_file <- file.path(base_dir, paste0(cts[c], '_', value_col, '_merged_100reps.xlsx'))

  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) { stop("Missing files:\n", paste(missing_files, collapse = "\n")) }

  rep_tables <- map2(file_paths,rep_ids, read_one_rep)
  merged_result <- reduce(rep_tables, full_join, by = key_col)

  write.xlsx(merged_result, file = output_file, overwrite = TRUE, freezePane = TRUE, autoFilter = TRUE )
  message("Merged file saved to: ", output_file)

}









