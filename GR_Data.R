setwd("C:/Moein/projects/Mutability/5_14_2025/rawdata")

library(vcfR)

# Your input files
vcf_files <- c("mev.snp.vcf", "N2.snp.vcf", "PB306.snp.vcf", "rare.snp.vcf")

outdir <- "vcf_exports"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

process_vcf_to_long <- function(infile) {
  message("Processing: ", infile)
  strain <- sub("\\.snp\\.vcf(\\.gz)?$", "", basename(infile), perl = TRUE)  # "mev", "N2", ...

  vcf <- read.vcfR(infile, verbose = FALSE)

  # Fixed fields
  fix <- as.data.frame(getFIX(vcf), stringsAsFactors = FALSE)
  # Variant column REF>ALT (uses first ALT if multiple; for strictness, normalize VCF beforehand)
  Variant <- paste0(fix$REF, ">", fix$ALT)

  # Genotypes (variants x samples)
  gt <- as.data.frame(extract.gt(vcf, element = "GT", as.numeric = FALSE), check.names = FALSE)

  # Flatten GT into a vector (column-wise), build matching sample vector
  gt_vec   <- as.vector(as.matrix(gt))
  samples  <- rep(colnames(gt), each = nrow(gt))
  # Replicate fixed per-sample
  CHROM_rep <- rep(fix$CHROM, times = ncol(gt))
  POS_rep   <- rep(as.integer(fix$POS), times = ncol(gt))
  Var_rep   <- rep(Variant, times = ncol(gt))

  # Keep rows where genotype shows any ALT allele (not 0/0 or ./.)
  has_alt <- !is.na(gt_vec) & gt_vec != "./." & grepl("[1-9]", gt_vec)

  out <- data.frame(
    strain  = strain,
    line    = samples[has_alt],
    `#CHROM`= CHROM_rep[has_alt],
    POS     = POS_rep[has_alt],
    Variant = Var_rep[has_alt],
    stringsAsFactors = FALSE
  )

  # Per-file (optional): write the filtered long table
  per_file <- file.path(outdir, paste0(strain, ".long_by_sample.tsv"))
  write.table(out, file = per_file, sep = "\t", quote = FALSE, row.names = FALSE)

  return(out)
}

# Process all and combine
all_long <- do.call(rbind, lapply(vcf_files, process_vcf_to_long))

# Write the single final file you asked for
final_path <- file.path(outdir, "FINAL_strain_line_CHROM_POS_Variant.tsv")
write.table(all_long, file = final_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("Done. Wrote final table to: ", final_path)

#End
