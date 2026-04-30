setwd("C:/Moein/projects/Mutability/Hipergator/Recombination")

library(readxl)
library(data.table)

# =========================
# 1. Load recombination map
# =========================
recmap <- read_excel("recmaps_Paree_2024_Genetics.xlsx", sheet = "ABCD")
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
  
  output_file <- sub(
    "\\.tsv$",
    "_with_recomb_rate.tsv",
    input_file
  )
  
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
  "XZ7.snp.highest_impact_snpEff.tsv",
  "XZ7.indel.highest_impact_snpEff.tsv",
  "XZ8.snp.highest_impact_snpEff.tsv",
  "XZ8.indel.highest_impact_snpEff.tsv",
  "PB306.snp.highest_impact_snpEff.tsv",
  "PB306.indel.highest_impact_snpEff.tsv",
  "N2.snp.highest_impact_snpEff.tsv",
  "N2.indel.highest_impact_snpEff.tsv",
  "mev.snp.highest_impact_snpEff.tsv",
  "mev.indel.highest_impact_snpEff.tsv",
  "J.snp.highest_impact_snpEff.tsv",
  "J.indel.highest_impact_snpEff.tsv",
  "G.snp.highest_impact_snpEff.tsv",
  "G.indel.highest_impact_snpEff.tsv",
  "E.snp.highest_impact_snpEff.tsv",
  "E.indel.highest_impact_snpEff.tsv",
  "C.snp.highest_impact_snpEff.tsv",
  "C.indel.highest_impact_snpEff.tsv",
  "M.snp.highest_impact_snpEff.tsv",
  "M.indel.highest_impact_snpEff.tsv"
)

# Optional: keep only files that exist
variant_files <- variant_files[file.exists(variant_files)]

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

fwrite(summary_table, "recomb_rate_annotation_summary.tsv", sep = "\t")

summary_table