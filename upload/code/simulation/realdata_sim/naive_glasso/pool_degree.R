library(readxl)
library(dplyr)
library(purrr)
library(openxlsx)

base_dir <- "reps"

reps <- sprintf("rep_%03d", 1:100)
files <- file.path(base_dir, reps, "naive_glasso/degree.xlsx")

stopifnot(all(file.exists(files)))

sheet_names <- excel_sheets(files[1])

merge_one_sheet <- function(sheet) {
  
  dat_list <- map2(files, reps, function(file, rep_id) {
    
    df <- read_excel(file, sheet = sheet)
    
    stopifnot(ncol(df) == 3)
    
    key_name <- names(df)[1]
    
    names(df)[2:3] <- paste0(names(df)[2:3], "_", rep_id)
    
    df
  })
  
  merged <- reduce(dat_list, full_join, by = names(dat_list[[1]])[1])
  
  merged
}

res_list <- map(sheet_names, merge_one_sheet)
names(res_list) <- sheet_names

out_file <- "./degree_naive_glasso.xlsx"

wb <- createWorkbook()

for (sheet in sheet_names) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, res_list[[sheet]])
}

saveWorkbook(wb, out_file, overwrite = TRUE)