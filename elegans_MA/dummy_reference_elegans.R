# Set current folder of the R script as the working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()
# =============================================================================
# Dummy Reference Generator for C. elegans FN Rate Estimation
# =============================================================================
# Produces:
#   1. elegans_snp_dummy.fa            - reference with 1200 substitutions only
#   2. elegans_indel_dummy.fa          - reference with 12000 indels only
#   3. elegans_ground_truth_snps.vcf   - injected SNPs in ORIGINAL coordinates
#   4. elegans_ground_truth_indels.vcf - injected indels in ORIGINAL coordinates
#   5. elegans_indel_offset_table.csv  - per-indel offset map 
#   6. elegans_simulation_log.csv      - full details of every injected variant
# Spectrum source: C. elegans SR 10X (E+G lines combined, n=8172 SNPs, n=1932 indels)
# SNP spectrum from VCF files (E_bp_GVCFs + G_bp_GVCFs mutation_filtered)
# Indel length probs: empirical bin proportions (1bp|2bp|3-5bp|6-10bp|11-20bp)
#   with equal weight within each bin
# Chromosome names preserved exactly as in reference (I, II, III, IV, V, X)
# MtDNA excluded
# Positions allocated proportional to chromosome length
# =============================================================================

library(Biostrings)
library(dplyr)

set.seed(42)

# =============================================================================
# USER SETTINGS - edit here only
# =============================================================================
REF_FA     <- "c_elegans.PRJNA13758.WS292.genomic.fa"
OUTDIR     <- "."
N_SNP      <- 1200
N_INDEL    <- 12000
MTDNA_PATT <- "MtDNA|chrM|MT|mitochondri"

# --- C. elegans SR 10X SNP spectrum (E+G combined, pyrimidine-referenced) ---
# n = 8172 SNPs
snp_spectrum <- c("C>A" = 0.1327, "C>G" = 0.0852, "C>T" = 0.3064,
                  "T>A" = 0.2393, "T>C" = 0.1421, "T>G" = 0.0945)
snp_spectrum <- snp_spectrum / sum(snp_spectrum)   # renormalize

# --- C. elegans indel parameters---
# Indels >20 bp excluded (146 skipped: mostly mobile elements / repeat expansions)
# Ins/Del ratio: 734/1932 = 38.0% ins, 1198/1932 = 62.0% del
# Note: C. elegans is strongly deletion-biased vs C. briggsae (~52/48)
INS_PROB <- 0.3800
DEL_PROB <- 0.6200

# Length probabilities use empirical bin proportions with equal weight
# within each bin: 1bp | 2bp | 3-5bp | 6-10bp | 11-20bp
# Insertions bin proportions: 53.81% | 12.81% | 18.39% | 9.81% | 5.18%
# Deletions  bin proportions: 51.92% |  9.18% | 16.53% | 12.35% | 10.02%
# Expected counts per length at n=12000 (4560 ins, 7440 del):
#   lengths 11-20: ~24 ins each, ~75 del each (~98 total) - sufficient for FNR
ins_len_probs <- c(0.53814714, 0.12806540,
                   0.06130790, 0.06130790, 0.06130790,
                   0.01961853, 0.01961853, 0.01961853, 0.01961853, 0.01961853,
                   0.00517711, 0.00517711, 0.00517711, 0.00517711, 0.00517711,
                   0.00517711, 0.00517711, 0.00517711, 0.00517711, 0.00517711)

del_len_probs <- c(0.51919866, 0.09181970,
                   0.05509182, 0.05509182, 0.05509182,
                   0.02470785, 0.02470785, 0.02470785, 0.02470785, 0.02470785,
                   0.01001669, 0.01001669, 0.01001669, 0.01001669, 0.01001669,
                   0.01001669, 0.01001669, 0.01001669, 0.01001669, 0.01001669)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

complement_base <- function(b) {
    c(A="T", T="A", C="G", G="C")[b]
}

draw_alt <- function(ref_base) {
    if (ref_base == "C") {
        alts    <- c("A","G","T")
        weights <- snp_spectrum[c("C>A","C>G","C>T")]
        return(sample(alts, 1, prob=weights))
    } else if (ref_base == "T") {
        alts    <- c("A","C","G")
        weights <- snp_spectrum[c("T>A","T>C","T>G")]
        return(sample(alts, 1, prob=weights))
    } else if (ref_base == "G") {
        alts    <- c("T","C","A")
        weights <- snp_spectrum[c("C>A","C>G","C>T")]
        return(sample(alts, 1, prob=weights))
    } else if (ref_base == "A") {
        alts    <- c("T","G","C")
        weights <- snp_spectrum[c("T>A","T>C","T>G")]
        return(sample(alts, 1, prob=weights))
    } else {
        return(NA)
    }
}

proportional_allocation <- function(lengths, total) {
    raw    <- lengths / sum(lengths) * total
    floors <- floor(raw)
    remain <- raw - floors
    deficit <- total - sum(floors)
    top_idx <- order(remain, decreasing=TRUE)[seq_len(deficit)]
    floors[top_idx] <- floors[top_idx] + 1
    return(as.integer(floors))
}

random_seq <- function(n) {
    paste(sample(c("A","T","C","G"), n, replace=TRUE), collapse="")
}

write_vcf <- function(df, filepath, ref_fa_name, genome) {
    contig_lines <- paste0("##contig=<ID=", names(genome),
                           ",length=", width(genome), ">")
    header <- c(
        "##fileformat=VCFv4.2",
        "##source=dummy_reference_elegans.R",
        paste0("##reference=", ref_fa_name),
        contig_lines,
        paste0('##INFO=<ID=SIMTYPE,Number=1,Type=String,',
               'Description="Simulated variant type">'),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    writeLines(header, filepath)
    write.table(df, filepath, sep="\t",
                col.names=FALSE, row.names=FALSE,
                quote=FALSE, append=TRUE)
    cat("  Written:", filepath, "\n")
}

# =============================================================================
# STEP 1: Load reference, identify chromosomes
# =============================================================================
cat("Loading reference genome...\n")
genome_full <- readDNAStringSet(REF_FA)

mtdna_idx <- grep(MTDNA_PATT, names(genome_full), ignore.case=TRUE)
if (length(mtdna_idx) > 0) {
    cat("  Excluding MtDNA contig(s):",
        paste(names(genome_full)[mtdna_idx], collapse=", "), "\n")
    genome <- genome_full[-mtdna_idx]
} else {
    cat("  No MtDNA contig detected.\n")
    genome <- genome_full
}

chr_names   <- names(genome)
chr_lengths <- width(genome)
names(chr_lengths) <- chr_names

cat("  Chromosomes:", paste(chr_names, collapse=", "), "\n")
cat("  Lengths:", paste(chr_lengths, collapse=", "), "\n")

# =============================================================================
# STEP 2: Allocate counts per chromosome proportional to length
# =============================================================================
snp_per_chr   <- proportional_allocation(chr_lengths, N_SNP)
indel_per_chr <- proportional_allocation(chr_lengths, N_INDEL)
names(snp_per_chr)   <- chr_names
names(indel_per_chr) <- chr_names

cat("\nSNP allocation per chromosome:\n");   print(snp_per_chr)
cat("Total SNPs:", sum(snp_per_chr), "\n")
cat("\nIndel allocation per chromosome:\n"); print(indel_per_chr)
cat("Total indels:", sum(indel_per_chr), "\n")

# =============================================================================
# STEP 3: Simulate SNPs
# =============================================================================
cat("\n--- Simulating SNPs ---\n")

snp_vcf_rows <- list()
snp_log_rows <- list()
snp_genome   <- genome

for (chr in chr_names) {
    cat("  Processing", chr, "...\n")
    seq_str <- as.character(genome[[chr]])
    n_sites <- snp_per_chr[chr]
    chr_len <- nchar(seq_str)

    positions <- c()
    attempts  <- 0
    while (length(positions) < n_sites && attempts < 20) {
        candidates <- sample(2:(chr_len - 1), min(n_sites * 5, chr_len - 2),
                             replace=FALSE)
        ref_bases  <- substr(seq_str, candidates, candidates)
        valid      <- candidates[ref_bases %in% c("A","T","C","G")]
        positions  <- unique(c(positions, valid))
        attempts   <- attempts + 1
    }
    positions <- sort(positions[seq_len(min(n_sites, length(positions)))])

    seq_vec <- strsplit(seq_str, "")[[1]]
    for (pos in positions) {
        ref_base <- seq_vec[pos]
        alt_base <- draw_alt(ref_base)
        if (is.na(alt_base)) next

        if (ref_base %in% c("C","T")) {
            mut_type <- paste0(ref_base, ">", alt_base)
        } else {
            mut_type <- paste0(complement_base(ref_base), ">",
                               complement_base(alt_base))
        }

        seq_vec[pos] <- alt_base

        snp_vcf_rows[[length(snp_vcf_rows)+1]] <- data.frame(
            CHROM=chr, POS=pos, ID=".", REF=ref_base, ALT=alt_base,
            QUAL=".", FILTER="PASS",
            INFO=paste0("SIMTYPE=SNP;MUTTYPE=", mut_type),
            stringsAsFactors=FALSE)

        snp_log_rows[[length(snp_log_rows)+1]] <- data.frame(
            CHROM=chr, POS=pos, TYPE="SNP", MUT_TYPE=mut_type,
            REF=ref_base, ALT=alt_base, INDEL_LEN=NA,
            stringsAsFactors=FALSE)
    }

    snp_genome[[chr]] <- DNAString(paste(seq_vec, collapse=""))
    cat("    Injected", length(positions), "SNPs\n")
}

# =============================================================================
# STEP 4: Simulate Indels
# =============================================================================
cat("\n--- Simulating Indels ---\n")

indel_vcf_rows    <- list()
indel_log_rows    <- list()
offset_table_rows <- list()
indel_genome      <- genome

for (chr in chr_names) {
    cat("  Processing", chr, "...\n")
    seq_str <- as.character(genome[[chr]])
    n_sites <- indel_per_chr[chr]
    chr_len <- nchar(seq_str)

    positions <- sort(sample(2:(chr_len - 22), n_sites, replace=FALSE))

    is_ins  <- sample(c(TRUE, FALSE), n_sites, replace=TRUE,
                      prob=c(INS_PROB, DEL_PROB))
    ins_len <- sample(1:20, n_sites, replace=TRUE, prob=ins_len_probs)
    del_len <- sample(1:20, n_sites, replace=TRUE, prob=del_len_probs)

    # Remove overlapping positions
    used_ranges <- c()
    keep        <- rep(TRUE, n_sites)
    for (i in seq_len(n_sites)) {
        pos <- positions[i]
        len <- if (is_ins[i]) 1 else del_len[i]
        end <- pos + len
        conflict <- any(sapply(used_ranges, function(r) {
            !(end < r[1] || pos > r[2])
        }))
        if (conflict) {
            keep[i] <- FALSE
        } else {
            used_ranges <- c(used_ranges, list(c(pos, end)))
        }
    }

    positions <- positions[keep]
    is_ins    <- is_ins[keep]
    ins_len   <- ins_len[keep]
    del_len   <- del_len[keep]
    n_actual  <- sum(keep)
    if (n_actual < n_sites)
        cat("    Note:", n_sites - n_actual,
            "indels removed due to position overlap\n")

    seq_vec    <- strsplit(seq_str, "")[[1]]
    cum_offset <- 0L

    for (i in seq_along(positions)) {
        orig_pos    <- positions[i]
        dummy_pos   <- orig_pos + cum_offset
        anchor_base <- seq_vec[dummy_pos]

        if (is_ins[i]) {
            ins_seq   <- random_seq(ins_len[i])
            ref_field <- anchor_base
            alt_field <- paste0(anchor_base, ins_seq)

            seq_vec <- c(seq_vec[seq_len(dummy_pos)],
                         strsplit(ins_seq, "")[[1]],
                         seq_vec[(dummy_pos+1):length(seq_vec)])

            cum_offset <- cum_offset + ins_len[i]

            offset_table_rows[[length(offset_table_rows)+1]] <- data.frame(
                CHROM=chr, ORIG_POS=orig_pos, DUMMY_POS=dummy_pos,
                INDEL_TYPE="INS", INDEL_LEN=ins_len[i],
                CUM_OFFSET=cum_offset, stringsAsFactors=FALSE)

            indel_vcf_rows[[length(indel_vcf_rows)+1]] <- data.frame(
                CHROM=chr, POS=orig_pos, ID=".", REF=ref_field,
                ALT=alt_field, QUAL=".", FILTER="PASS",
                INFO=paste0("SIMTYPE=INS;INDELSIZE=", ins_len[i]),
                stringsAsFactors=FALSE)

            indel_log_rows[[length(indel_log_rows)+1]] <- data.frame(
                CHROM=chr, POS=orig_pos, TYPE="INS", MUT_TYPE="INS",
                REF=ref_field, ALT=alt_field, INDEL_LEN=ins_len[i],
                stringsAsFactors=FALSE)

        } else {
            dlen    <- del_len[i]
            end_idx <- dummy_pos + dlen
            if (end_idx > length(seq_vec)) next

            deleted_bases <- paste(seq_vec[(dummy_pos+1):end_idx], collapse="")
            ref_field <- paste0(anchor_base, deleted_bases)
            alt_field <- anchor_base

            seq_vec <- c(seq_vec[seq_len(dummy_pos)],
                         seq_vec[(end_idx+1):length(seq_vec)])

            cum_offset <- cum_offset - dlen

            offset_table_rows[[length(offset_table_rows)+1]] <- data.frame(
                CHROM=chr, ORIG_POS=orig_pos, DUMMY_POS=dummy_pos,
                INDEL_TYPE="DEL", INDEL_LEN=dlen,
                CUM_OFFSET=cum_offset, stringsAsFactors=FALSE)

            indel_vcf_rows[[length(indel_vcf_rows)+1]] <- data.frame(
                CHROM=chr, POS=orig_pos, ID=".", REF=ref_field,
                ALT=alt_field, QUAL=".", FILTER="PASS",
                INFO=paste0("SIMTYPE=DEL;INDELSIZE=", dlen),
                stringsAsFactors=FALSE)

            indel_log_rows[[length(indel_log_rows)+1]] <- data.frame(
                CHROM=chr, POS=orig_pos, TYPE="DEL", MUT_TYPE="DEL",
                REF=ref_field, ALT=alt_field, INDEL_LEN=dlen,
                stringsAsFactors=FALSE)
        }
    }

    indel_genome[[chr]] <- DNAString(paste(seq_vec, collapse=""))
    cat("    Injected", length(positions), "indels  |",
        "Ins:", sum(is_ins), " Del:", sum(!is_ins), "\n")
}

# =============================================================================
# STEP 5: Write dummy FASTA files
# =============================================================================
cat("\n--- Writing FASTA files ---\n")

write_dummy_fa <- function(modified_chrs, original_full, mtdna_idx, outpath) {
    out_genome <- original_full
    for (chr in names(modified_chrs)) {
        out_genome[[chr]] <- modified_chrs[[chr]]
    }
    writeXStringSet(out_genome, outpath)
    cat("  Written:", outpath, "\n")
}

write_dummy_fa(snp_genome,   genome_full, mtdna_idx,
               file.path(OUTDIR, "elegans_snp_dummy.fa"))
write_dummy_fa(indel_genome, genome_full, mtdna_idx,
               file.path(OUTDIR, "elegans_indel_dummy.fa"))

# =============================================================================
# STEP 6: Write VCF files
# =============================================================================
cat("\n--- Writing VCF files ---\n")

snp_vcf_df   <- bind_rows(snp_vcf_rows)   %>% arrange(CHROM, POS)
indel_vcf_df <- bind_rows(indel_vcf_rows) %>% arrange(CHROM, POS)

write_vcf(snp_vcf_df,
          file.path(OUTDIR, "elegans_ground_truth_snps.vcf"),   REF_FA, genome_full)
write_vcf(indel_vcf_df,
          file.path(OUTDIR, "elegans_ground_truth_indels.vcf"), REF_FA, genome_full)

# =============================================================================
# STEP 7: Write offset table and simulation log
# =============================================================================
cat("\n--- Writing offset table and simulation log ---\n")

offset_df <- bind_rows(offset_table_rows)
write.csv(offset_df,
          file.path(OUTDIR, "elegans_indel_offset_table.csv"), row.names=FALSE)
cat("  Written: indel_offset_table.csv\n")

log_df <- bind_rows(snp_log_rows, indel_log_rows) %>% arrange(CHROM, POS)
write.csv(log_df,
          file.path(OUTDIR, "elegans_simulation_log.csv"), row.names=FALSE)
cat("  Written: simulation_log.csv\n")

# =============================================================================
# STEP 8: Summary
# =============================================================================
cat("\n")
cat("============================================================\n")
cat(" SIMULATION COMPLETE - C. elegans\n")
cat("============================================================\n")
cat(sprintf("  SNPs injected:   %d\n", nrow(snp_vcf_df)))
cat(sprintf("  Indels injected: %d\n", nrow(indel_vcf_df)))
cat(sprintf("    Insertions:    %d\n", sum(grepl("SIMTYPE=INS", indel_vcf_df$INFO))))
cat(sprintf("    Deletions:     %d\n", sum(grepl("SIMTYPE=DEL", indel_vcf_df$INFO))))
cat("\n  SNPs per chromosome:\n");    print(table(snp_vcf_df$CHROM))
cat("\n  Indels per chromosome:\n");  print(table(indel_vcf_df$CHROM))
cat("\n  Output files:\n")
cat("    elegans_snp_dummy.fa\n")
cat("    elegans_indel_dummy.fa\n")
cat("    elegans_ground_truth_snps.vcf\n")
cat("    elegans_ground_truth_indels.vcf\n")
cat("    elegans_indel_offset_table.csv\n")
cat("    elegans_simulation_log.csv\n")
cat("============================================================\n")
