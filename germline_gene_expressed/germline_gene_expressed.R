# Load required packages
library(readxl)
library(dplyr)
library(writexl)

setwd("C:/Moein/projects/Mutability/Hipergator/gene_expressed_germline")

# File path
file <- "germline_gene.xlsx"

# Sheets to process
sheets <- c("N2", "PB306", "mev", "C", "E", "G", "J", "M", "XZ7", "XZ8")

# Read germline reference
germline <- read_excel(file, sheet = "GermlineSageLib")

# Clean germline gene names
germline_names <- germline$`Gene Name` %>%
  tolower() %>%
  trimws()

# Function to process each sheet
process_sheet <- function(sheet_name) {
  df <- read_excel(file, sheet = sheet_name)
  
  df %>%
    mutate(
      gene_clean = tolower(trimws(`Gene annotation`)),
      match_flag = ifelse(gene_clean %in% germline_names, 1, NA)
    ) %>%
    select(-gene_clean)
}

# Apply function to all sheets
results_list <- lapply(sheets, process_sheet)

# Name the list (important for Excel sheet names)
names(results_list) <- sheets

# Write to Excel with multiple sheets
write_xlsx(results_list, "gene_expressed_germline_results.xlsx")



#End
