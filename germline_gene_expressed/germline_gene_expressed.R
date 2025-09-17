library(readxl)

setwd("C:/Moein/projects/Mutability/5_14_2025/gene_expressed_germline")

# Read the two sheets
germline <- read_excel("germline_gene.xls", sheet = "GermlineSageLib")
ma_bs    <- read_excel("germline_gene.xls", sheet = "MA_BS")

# Merge based on matching columns
matched <- merge(ma_bs, germline,
                 by.x = "Gene annotation", 
                 by.y = "Gene Name",         
                 all = FALSE)              

head(matched)

library(readxl)
library(dplyr)

germline <- read_excel("germline_gene.xls", sheet = "GermlineSageLib")
ma_bs    <- read_excel("germline_gene.xls", sheet = "MA_BS")

matched_2 <- ma_bs %>%
  inner_join(germline, by = c("Gene annotation" = "Gene Name"))

head(matched)