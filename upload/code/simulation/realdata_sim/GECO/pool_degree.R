library(readxl)
library(openxlsx)

base_dir <- "reps"
rep_ids <- 1:100

file_paths <- file.path(base_dir, paste0('rep_',sprintf("%03d", rep_ids)), 'degree', "degree.xlsx")
output_file <- file.path(base_dir, "degree_geco_merged_100reps.xlsx")

# check if the file exists
if (any(!file.exists(file_paths))) {stop( "Missing files:\n", paste(file_paths[!file.exists(file_paths)], collapse = "\n") )}

# 以第一个文件的 sheet 名称为准
sheet_names <- excel_sheets(file_paths[1])

merged_sheets <- lapply(sheet_names, function(sheet_name) {

  # 读取第一个 replicate，保留第一列基因
  first_dat <- as.data.frame( read_excel(file_paths[1], sheet = sheet_name))
  if (ncol(first_dat) != 3) { stop("Each sheet should contain exactly three columns.") }

  gene_col <- names(first_dat)[1]
  genes <- first_dat[[1]]

  # 给第一个 rep 的结果列重命名
  result <- first_dat
  names(result)[2:3] <- paste0( names(result)[2:3], "_rep", sprintf("%03d", rep_ids[1]))

  # 依次添加 rep 2-100 的后两列
  for (i in 2:length(rep_ids)) {

    dat <- as.data.frame( read_excel(file_paths[i], sheet = sheet_name) )

    if (ncol(dat) != 3) {
      stop( "Sheet ", sheet_name, " in rep ", rep_ids[i], " does not contain three columns." )
    }

    # 检查基因及顺序是否完全一致
    if (!identical(dat[[1]], genes)) { stop( "Gene names or order differ in sheet ", sheet_name, " for rep ", rep_ids[i] )}

    names(dat)[2:3] <- paste0( names(dat)[2:3], "_rep", sprintf("%03d", rep_ids[i]) )

    result <- cbind( result, dat[, 2:3, drop = FALSE] )
  }

  result
})

names(merged_sheets) <- sheet_names

# 直接将 list 写成多个 sheet
write.xlsx(merged_sheets,file = output_file,overwrite = TRUE)
message("Merged workbook saved to: ", output_file)