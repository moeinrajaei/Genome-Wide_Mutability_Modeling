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
