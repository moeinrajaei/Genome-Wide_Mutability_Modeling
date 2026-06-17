setwd("C:/Moein/projects/Mutability/Hipergator")

set.seed(123)

# Genome definition
chrom <- c("I", "II", "III", "IV", "V", "X", "MtDNA")
chr_size <- c(
  15072434,
  15279421,
  13783801,
  17493829,
  20924180,
  17718942,
  13794
)

n_random <- 1000000

# Read mutant sites
mut <- read.delim(
  "all_combined_master_file.tsv",
  stringsAsFactors = FALSE
)

# Create unique identifiers for mutant sites
mut_sites <- unique(
  paste(mut$CHROM, mut$POS, sep = ":")
)

genome_size <- sum(chr_size)

random_sites <- data.frame()

while (nrow(random_sites) < n_random) {
  
  n_needed <- n_random - nrow(random_sites)
  
  # Oversample to compensate for filtering
  sampled_chr <- sample(
    chrom,
    size = n_needed * 2,
    replace = TRUE,
    prob = chr_size / genome_size
  )
  
  sampled_pos <- mapply(
    function(chr) sample.int(chr_size[chrom == chr], 1),
    sampled_chr
  )
  
  candidate <- data.frame(
    CHROM = sampled_chr,
    POS = sampled_pos,
    stringsAsFactors = FALSE
  )
  
  candidate$site_id <- paste(
    candidate$CHROM,
    candidate$POS,
    sep = ":"
  )
  
  # Remove mutant sites
  candidate <- candidate[
    !(candidate$site_id %in% mut_sites),
  ]
  
  # Remove duplicate sites
  candidate <- unique(candidate)
  
  random_sites <- unique(
    rbind(random_sites, candidate)
  )
  
  # Keep only the required number of sites
  if (nrow(random_sites) > n_random) {
    random_sites <- random_sites[1:n_random, ]
  }
}

# Keep desired columns
random_sites <- random_sites[, c("CHROM", "POS")]

# Add strain column
random_sites$Strain <- "nonmutant"

# Reorder columns
random_sites <- random_sites[, c("Strain", "CHROM", "POS")]

# Save output
write.table(
  random_sites,
  file = "random_1000000_non_mutant_sites.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
#End


###nonmutated_OR
# Set Working Directory:
setwd("C:/Moein/projects/Mutability/Hipergator/non_mutated")


## 1- Read data from the files into a data frame:
library(readxl)
library(dplyr)
library(data.table)

## 1- Read data
bed <- read_excel("Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx",  sheet = "bed")
nonmutated <- read.delim("random_1000000_non_mutant_sites.txt", header = TRUE,  stringsAsFactors = FALSE)

## 2- Make column names compatible
bed <- bed %>%
  rename(CHROM = chr)

## 3- Convert to data.table
bed_dt <- as.data.table(bed)
nonmutated_dt <- as.data.table(nonmutated)

## 4- Keep needed columns
bed_dt <- bed_dt[, .(
  CHROM,
  start,
  bed_start = start
)]

nonmutated_dt <- nonmutated_dt[, .(
  CHROM,
  POS
)]

## 5- Make numeric and sort
bed_dt[, start := as.numeric(start)]
bed_dt[, bed_start := as.numeric(bed_start)]
nonmutated_dt[, POS := as.numeric(POS)]

setorder(bed_dt, CHROM, start)
setorder(nonmutated_dt, CHROM, POS)

## 6- Find closest upstream bed start for each nonmutated site
result <- bed_dt[
  nonmutated_dt,
  on = .(CHROM, start < POS),
  mult = "last"
]

## 7- Fix output column names
result <- result[, .(
  CHROM,
  POS = start,
  start = bed_start
)]

## 7- Calculate difference
result[, difference := POS - start]

## 8- Rename output column
final_result <- result[, .(
  CHROM,
  POS,
  start,
  Max_Negative_Difference = difference
)]


## 9- Save output
write.table(
  final_result,
  file = "nonmutated_OR.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


###Local_Recombination_Rate:

library(readxl)
library(data.table)

# =========================
# 1. Load recombination map
# =========================
recmap <- read_excel("C:/Moein/projects/Mutability/Hipergator/Recombination/recmaps_Paree_2024_Genetics.xlsx", sheet = "ABCD")
setDT(recmap)

recmap <- recmap[, .(
  chrom = as.character(chrom),
  pos = as.numeric(pos),
  genetic = as.numeric(genetic)
)]

recmap[chrom == "1", chrom := "I"]
recmap[chrom == "2", chrom := "II"]
recmap[chrom == "3", chrom := "III"]
recmap[chrom == "4", chrom := "IV"]
recmap[chrom == "5", chrom := "V"]
recmap[chrom == "6", chrom := "X"]
recmap[chrom == "7", chrom := "MtDNA"]

# =========================
# 2. Chromosome sizes
# =========================
fai <- data.table(
  chrom = c("I", "II", "III", "IV", "V", "X", "MtDNA"),
  chr_size = c(
    15072434,
    15279421,
    13783801,
    17493829,
    20924180,
    17718942,
    13794
  )
)

# =========================
# 3. Add chromosome boundaries
# =========================
setorder(recmap, chrom, pos)

boundary_rows <- recmap[, {
  
  chr_len <- fai[chrom == .BY$chrom, chr_size]
  
  if (length(chr_len) == 0 || is.na(chr_len) || .N < 2) {
    return(data.table(chrom = character(), pos = numeric(), genetic = numeric()))
  }
  
  ord <- order(pos)
  first_two <- ord[1:2]
  last_two <- tail(ord, 2)
  
  start_slope <- diff(genetic[first_two]) / diff(pos[first_two])
  end_slope <- diff(genetic[last_two]) / diff(pos[last_two])
  
  first_pos <- pos[first_two[1]]
  first_genetic <- genetic[first_two[1]]
  
  last_pos <- pos[last_two[2]]
  last_genetic <- genetic[last_two[2]]
  
  genetic_at_start <- first_genetic - start_slope * first_pos
  genetic_at_end <- last_genetic + end_slope * (chr_len - last_pos)
  
  data.table(
    chrom = c(.BY$chrom, .BY$chrom),
    pos = c(0, chr_len),
    genetic = c(genetic_at_start, genetic_at_end)
  )
  
}, by = chrom]

recmap2 <- rbind(recmap, boundary_rows, fill = TRUE)
recmap2 <- unique(recmap2, by = c("chrom", "pos"))
setorder(recmap2, chrom, pos)

# =========================
# 4. Create recombination intervals
# =========================
rec_intervals <- recmap2[, .(
  recmap_start = pos,
  recmap_end = shift(pos, type = "lead"),
  genetic_start = genetic,
  genetic_end = shift(genetic, type = "lead")
), by = chrom]

rec_intervals <- rec_intervals[!is.na(recmap_end)]
rec_intervals <- rec_intervals[recmap_end > recmap_start]

rec_intervals[, recomb_rate_cM_Mb :=
                (genetic_end - genetic_start) /
                ((recmap_end - recmap_start) / 1e6)
]

# =========================
# 5. Function to annotate one file
# =========================
annotate_recomb_rate <- function(input_file) {
  
  message("Processing: ", input_file)
  
  variants <- fread(input_file)
  
  variants[, CHROM := as.character(CHROM)]
  variants[, POS := as.numeric(POS)]
  variants[, chrom := CHROM]
  
  variants[chrom == "1", chrom := "I"]
  variants[chrom == "2", chrom := "II"]
  variants[chrom == "3", chrom := "III"]
  variants[chrom == "4", chrom := "IV"]
  variants[chrom == "5", chrom := "V"]
  variants[chrom == "6", chrom := "X"]
  variants[chrom == "7", chrom := "MtDNA"]
  
  variants[, `:=`(
    recmap_start = NA_real_,
    recmap_end = NA_real_,
    genetic_start = NA_real_,
    genetic_end = NA_real_,
    recomb_rate_cM_Mb = NA_real_
  )]
  
  for (chr in unique(variants$chrom)) {
    
    v_idx <- which(variants$chrom == chr)
    intervals_chr <- rec_intervals[chrom == chr]
    
    if (nrow(intervals_chr) == 0) next
    
    idx <- findInterval(
      variants$POS[v_idx],
      intervals_chr$recmap_start
    )
    
    valid <- idx > 0 &
      idx <= nrow(intervals_chr) &
      variants$POS[v_idx] <= intervals_chr$recmap_end[idx]
    
    variants$recmap_start[v_idx[valid]] <- intervals_chr$recmap_start[idx[valid]]
    variants$recmap_end[v_idx[valid]] <- intervals_chr$recmap_end[idx[valid]]
    variants$genetic_start[v_idx[valid]] <- intervals_chr$genetic_start[idx[valid]]
    variants$genetic_end[v_idx[valid]] <- intervals_chr$genetic_end[idx[valid]]
    variants$recomb_rate_cM_Mb[v_idx[valid]] <- intervals_chr$recomb_rate_cM_Mb[idx[valid]]
  }
  
  output_file <- sub("\\.[^.]+$", "_with_recomb_rate.tsv", input_file)
  
  fwrite(variants, output_file, sep = "\t")
  
  message("Saved: ", output_file)
  message("Variants total: ", nrow(variants))
  message("With recombination rate: ", sum(!is.na(variants$recomb_rate_cM_Mb)))
  message("Missing recombination rate: ", sum(is.na(variants$recomb_rate_cM_Mb)))
  message("")
  
  return(variants)
}

# =========================
# 6. Files to process
# =========================
variant_files <- c(
  "nonmutated_OR.txt")
# =========================
# 7. Run annotation for all files
# =========================
results <- lapply(variant_files, annotate_recomb_rate)

# =========================
# 8. Summary table
# =========================
summary_table <- rbindlist(lapply(seq_along(variant_files), function(i) {
  data.table(
    file = variant_files[i],
    total_variants = nrow(results[[i]]),
    variants_with_rate = sum(!is.na(results[[i]]$recomb_rate_cM_Mb)),
    variants_missing_rate = sum(is.na(results[[i]]$recomb_rate_cM_Mb))
  )
}))

fwrite(summary_table, "recombination_rate_summary.tsv", sep = "\t")




library(data.table)

# -----------------------------
# 1. Read files
# -----------------------------

base <- fread("random_1000000_non_mutant_sites.txt")
wormbase <- fread("random_1000000_non_mutant_sites.wormbase.unique.bed")
seqcontext <- fread("random_1000000_non_mutant_sites.all_merged.with_header.txt")
gc <- fread("random_1000000_non_mutant_sites.final.GC.tsv")

# -----------------------------
# 2. Standardize CHROM / POS columns
# -----------------------------

base[, CHROM := as.character(CHROM)]
base[, POS := as.integer(POS)]

seqcontext[, CHROM := as.character(CHROM)]
seqcontext[, POS := as.integer(POS)]

gc[, CHROM := as.character(CHROM)]
gc[, POS := as.integer(POS)]

# WormBase BED may not have header
# BED format: CHROM start0 end1 ...
if (!all(c("CHROM", "POS") %in% names(wormbase))) {
  setnames(
    wormbase,
    old = names(wormbase)[1:3],
    new = c("CHROM", "start0", "end1")
  )
  wormbase[, POS := as.integer(end1)]
}

wormbase[, CHROM := as.character(CHROM)]
wormbase[, POS := as.integer(POS)]

# -----------------------------
# 3. Remove duplicates
# -----------------------------

base <- unique(base, by = c("CHROM", "POS"))
seqcontext <- unique(seqcontext, by = c("CHROM", "POS"))
gc <- unique(gc, by = c("CHROM", "POS"))
wormbase <- unique(wormbase, by = c("CHROM", "POS"))

# -----------------------------
# 4. Sort by biological chromosome order
# -----------------------------

chr_order <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

sort_by_chr_pos <- function(dt) {
  dt[, chr_rank := match(CHROM, chr_order)]
  setorder(dt, chr_rank, POS)
  dt[, chr_rank := NULL]
  return(dt)
}

base <- sort_by_chr_pos(base)
seqcontext <- sort_by_chr_pos(seqcontext)
gc <- sort_by_chr_pos(gc)
wormbase <- sort_by_chr_pos(wormbase)

# -----------------------------
# 5. Check matching
# -----------------------------

cat("Base sites:", nrow(base), "\n")
cat("Sequence-context sites:", nrow(seqcontext), "\n")
cat("GC sites:", nrow(gc), "\n")
cat("WormBase annotated sites:", nrow(wormbase), "\n\n")

cat("Base sites missing from sequence-context file:",
    nrow(base[!seqcontext, on = .(CHROM, POS)]), "\n")

cat("Base sites missing from GC file:",
    nrow(base[!gc, on = .(CHROM, POS)]), "\n")

cat("Base sites without WormBase annotation:",
    nrow(base[!wormbase, on = .(CHROM, POS)]), "\n\n")

cat("Extra sequence-context sites not in base:",
    nrow(seqcontext[!base, on = .(CHROM, POS)]), "\n")

cat("Extra GC sites not in base:",
    nrow(gc[!base, on = .(CHROM, POS)]), "\n")

cat("Extra WormBase sites not in base:",
    nrow(wormbase[!base, on = .(CHROM, POS)]), "\n\n")

cat("Same row order for sequence-context file after sorting:",
    nrow(base) == nrow(seqcontext) &&
      all(base$CHROM == seqcontext$CHROM & base$POS == seqcontext$POS),
    "\n")

cat("Same row order for GC file after sorting:",
    nrow(base) == nrow(gc) &&
      all(base$CHROM == gc$CHROM & base$POS == gc$POS),
    "\n\n")

# -----------------------------
# 6. Rename overlapping feature columns
# -----------------------------

rename_overlaps <- function(x, existing_names, suffix) {
  overlap <- intersect(
    setdiff(names(x), c("CHROM", "POS")),
    existing_names
  )
  
  if (length(overlap) > 0) {
    setnames(x, overlap, paste0(overlap, suffix))
  }
  
  return(x)
}

seqcontext <- rename_overlaps(
  seqcontext,
  setdiff(names(base), c("CHROM", "POS")),
  "_seqcontext"
)

gc <- rename_overlaps(
  gc,
  union(
    setdiff(names(base), c("CHROM", "POS")),
    setdiff(names(seqcontext), c("CHROM", "POS"))
  ),
  "_GC"
)

wormbase <- rename_overlaps(
  wormbase,
  union(
    union(
      setdiff(names(base), c("CHROM", "POS")),
      setdiff(names(seqcontext), c("CHROM", "POS"))
    ),
    setdiff(names(gc), c("CHROM", "POS"))
  ),
  "_wormbase"
)

# -----------------------------
# 7. Merge all features by CHROM/POS
# -----------------------------

final <- merge(
  base,
  seqcontext,
  by = c("CHROM", "POS"),
  all.x = TRUE,
  sort = FALSE
)

final <- merge(
  final,
  gc,
  by = c("CHROM", "POS"),
  all.x = TRUE,
  sort = FALSE
)

final <- merge(
  final,
  wormbase,
  by = c("CHROM", "POS"),
  all.x = TRUE,
  sort = FALSE
)

# -----------------------------
# 8. Sort final file again
# -----------------------------

final <- sort_by_chr_pos(final)

# -----------------------------
# 9. Reorder columns
# -----------------------------

main_cols <- intersect(c("Strain", "CHROM", "POS"), names(final))

setcolorder(
  final,
  c(main_cols, setdiff(names(final), main_cols))
)

# -----------------------------
# 10. Final checks
# -----------------------------

cat("Final merged rows:", nrow(final), "\n")
cat("Final merged columns:", ncol(final), "\n")

cat("Final rows equal to base rows:",
    nrow(final) == nrow(base), "\n")

cat("Duplicated CHROM/POS in final:",
    anyDuplicated(final[, .(CHROM, POS)]), "\n")

# -----------------------------
# 11. Save final merged file
# -----------------------------

fwrite(
  final,
  "random_1000000_non_mutant_sites.all_features_merged.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Saved: random_1000000_non_mutant_sites.all_features_merged.tsv\n")



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



library(data.table)

snpEff <- fread("nonmutant_fake_mutations.highest_impact_snpEff.with_mono.tsv")

snpEff_keep <- snpEff[, .(
  CHROM,
  POS,
  REF,
  ALT,
  `Coding context` = `ANN[0].EFFECT`,
  `mononuc?` = mono,
  `putative effect` = `ANN[0].IMPACT`,
  `Gene annotation` = `ANN[0].GENE`
)]

fwrite(
  snpEff_keep,
  "nonmutant_fake_mutations.snpEff_selected_columns.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)




### # Load required packages
library(readxl)
library(dplyr)
library(writexl)
library(data.table)

# -----------------------------
# Files
# -----------------------------

germline_file <- "germline_gene.xlsx"

input_file <- "nonmutant_fake_mutations.snpEff_selected_columns.tsv"

# -----------------------------
# Read germline reference
# -----------------------------

germline <- read_excel(germline_file, sheet = "GermlineSageLib")

germline_names <- germline$`Gene Name` %>%
  tolower() %>%
  trimws()

# -----------------------------
# Read snpEff selected output
# -----------------------------

df <- fread(input_file)

# -----------------------------
# Add germline match flag
# -----------------------------

result <- df %>%
  mutate(
    gene_clean = tolower(trimws(`Gene annotation`)),
    germline_match = ifelse(gene_clean %in% germline_names, 1, NA)
  ) %>%
  select(-gene_clean)



# -----------------------------
# Save outputs
# -----------------------------

fwrite(
  result,
  "nonmutant_gene_expressed_germline_results.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)







#End

library(data.table)

# -----------------------------
# 1. Read files
# -----------------------------

features <- fread("random_1000000_non_mutant_sites.all_features_merged.tsv")
snpEff_germline <- fread("nonmutant_gene_expressed_germline_results.tsv")
recomb <- fread("nonmutated_OR_with_recomb_rate.tsv")

# -----------------------------
# 2. Standardize CHROM/POS
# -----------------------------

features[, CHROM := as.character(CHROM)]
features[, POS := as.integer(POS)]

snpEff_germline[, CHROM := as.character(CHROM)]
snpEff_germline[, POS := as.integer(POS)]

recomb[, CHROM := as.character(CHROM)]
recomb[, POS := as.integer(POS)]

# -----------------------------
# 3. Remove duplicates
# -----------------------------

features <- unique(features, by = c("CHROM", "POS"))
snpEff_germline <- unique(snpEff_germline, by = c("CHROM", "POS"))
recomb <- unique(recomb, by = c("CHROM", "POS"))

# -----------------------------
# 4. Check coordinate matching
# -----------------------------

cat("Feature file sites:", nrow(features), "\n")
cat("snpEff/germline sites:", nrow(snpEff_germline), "\n")
cat("Recombination sites:", nrow(recomb), "\n\n")

cat("Feature sites missing from snpEff/germline:",
    nrow(features[!snpEff_germline, on = .(CHROM, POS)]), "\n")

cat("Feature sites missing from recombination:",
    nrow(features[!recomb, on = .(CHROM, POS)]), "\n\n")

cat("Extra snpEff/germline sites not in features:",
    nrow(snpEff_germline[!features, on = .(CHROM, POS)]), "\n")

cat("Extra recombination sites not in features:",
    nrow(recomb[!features, on = .(CHROM, POS)]), "\n\n")

# -----------------------------
# 5. Rename overlapping columns before merge
# -----------------------------

rename_overlaps <- function(x, existing_names, suffix) {
  overlap <- intersect(
    setdiff(names(x), c("CHROM", "POS")),
    existing_names
  )
  
  if (length(overlap) > 0) {
    setnames(x, overlap, paste0(overlap, suffix))
  }
  
  x
}

snpEff_germline <- rename_overlaps(
  snpEff_germline,
  setdiff(names(features), c("CHROM", "POS")),
  "_snpEff"
)

recomb <- rename_overlaps(
  recomb,
  union(
    setdiff(names(features), c("CHROM", "POS")),
    setdiff(names(snpEff_germline), c("CHROM", "POS"))
  ),
  "_recomb"
)

# -----------------------------
# 6. Merge all features
# -----------------------------

final <- merge(
  features,
  snpEff_germline,
  by = c("CHROM", "POS"),
  all.x = TRUE,
  sort = FALSE
)

final <- merge(
  final,
  recomb,
  by = c("CHROM", "POS"),
  all.x = TRUE,
  sort = FALSE
)

# -----------------------------
# 7. Remove unwanted columns
# -----------------------------

cols_to_remove <- c(
  "chrom",
  "recmap_start",
  "recmap_end",
  "genetic_start",
  "genetic_end",
  "start",
  "V5",
  "V6",
  "start0",
  "end1"
)

final <- final[, !names(final) %in% cols_to_remove, with = FALSE]

# -----------------------------
# 8. Rename final columns
# -----------------------------

rename_map <- c(
  "germline_match" = "gene_expressed_germline",
  "recomb_rate_cM_Mb" = "recomb_rate",
  "Max_Negative_Difference" = "distance_to_ORI",
  "V4" = "lead/lag",
  "seq_3bp_first" = "motif_3bp_mutfirst",
  "seq_3bp_last" = "motif_3bp_mutlast",
  "seq_3bp_middle" = "motif_3bp_mutcenter",
  "seq_5bp_middle" = "motif_5bp_center",
  "seq_7bp_middle" = "motif_7bp_center",
  "seq_9bp_middle" = "motif_9bp_center",
  "seq_99bp_middle" = "seq_complexity_99bp_mutcenter",
  "GC_10bp" = "GC_11bp",
  "GC_100bp" = "GC_101bp",
  "GC_1kb" = "GC_1001bp",
  "GC_10kb" = "GC_10001bp",
  "GC_100kb" = "GC_100001bp"
)

existing <- intersect(names(rename_map), names(final))

setnames(
  final,
  old = existing,
  new = rename_map[existing]
)

# -----------------------------
# 9. Sort final file
# -----------------------------

chr_order <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

final[, chr_rank := match(CHROM, chr_order)]
setorder(final, chr_rank, POS)
final[, chr_rank := NULL]

# -----------------------------
# 10. Reorder main columns
# -----------------------------

main_cols <- intersect(
  c("Strain", "CHROM", "POS", "REF", "ALT"),
  names(final)
)

setcolorder(
  final,
  c(main_cols, setdiff(names(final), main_cols))
)

# -----------------------------
# 11. Final checks
# -----------------------------

cat("Final rows:", nrow(final), "\n")
cat("Final columns:", ncol(final), "\n")
cat("Duplicated CHROM/POS:",
    anyDuplicated(final[, .(CHROM, POS)]), "\n")

# -----------------------------
# 12. Save final output
# -----------------------------

fwrite(
  final,
  "random_1000000_non_mutant_sites.final_all_features.cleaned.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Saved: random_1000000_non_mutant_sites.final_all_features.cleaned.tsv\n")

