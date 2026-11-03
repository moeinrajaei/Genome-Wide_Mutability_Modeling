
# ============================================================
# XZ7 vs XZ8: Sequencing-lane comparison for one MA line
# ============================================================

# Install packages if necessary
packages <- c("tidyverse", "patchwork")

new_packages <- packages[
  !packages %in% installed.packages()[, "Package"]
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(tidyverse)
library(patchwork)

set.seed(123)

# ------------------------------------------------------------
# 1. INPUT AND OUTPUT
# ------------------------------------------------------------

input_file <- "XZ7_XZ8_Masterfile.txt"
output_dir <- "XZ7_XZ8_Analysis"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

graph_dir <- file.path(output_dir, "Graphs")
table_dir <- file.path(output_dir, "Tables")

dir.create(graph_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

# Save all graphs as JPG at 300 dpi
save_graph <- function(plot, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(graph_dir, paste0(filename, ".jpg")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

# Color palette
lane_colors <- c(
  "Shared" = "#78A899",
  "XZ7 only" = "#4C78A8",
  "XZ8 only" = "#E28E65"
)

status_levels <- c("Shared", "XZ7 only", "XZ8 only")

# ------------------------------------------------------------
# 2. LOAD AND CLEAN DATA
# ------------------------------------------------------------

df <- read.delim(
  input_file,
  header = TRUE,
  check.names = FALSE,
  na.strings = c("NA", "", "."),
  stringsAsFactors = FALSE
)

df <- df %>%
  filter(Strain %in% c("XZ7", "XZ8")) %>%
  mutate(
    Strain = factor(Strain, levels = c("XZ7", "XZ8")),
    POS = as.integer(POS),
    REF = toupper(trimws(REF)),
    ALT = toupper(trimws(ALT))
  )

# Keep first coding context before "&"
# Example:
# frameshift_variant&splice_donor_variant
# becomes frameshift_variant

df <- df %>%
  mutate(
    `Coding context` = stringr::str_trim(
      stringr::str_extract(
        `Coding context`,
        "^[^&]+"
      )
    )
  )

# Remove duplicate records within each lane.
# Review duplicates if they contain different annotations.

duplicate_records <- df %>%
  count(Strain, CHROM, POS, REF, ALT) %>%
  filter(n > 1)

write.csv(
  duplicate_records,
  file.path(table_dir, "Duplicate_records.csv"),
  row.names = FALSE
)

df <- df %>%
  distinct(
    Strain, CHROM, POS, REF, ALT,
    .keep_all = TRUE
  )

# ------------------------------------------------------------
# 3. CLASSIFY MUTATIONS BY DETECTION STATUS
# ------------------------------------------------------------

variant_status <- df %>%
  group_by(CHROM, POS, REF, ALT) %>%
  summarise(
    in_XZ7 = any(Strain == "XZ7"),
    in_XZ8 = any(Strain == "XZ8"),
    .groups = "drop"
  ) %>%
  mutate(
    Status = case_when(
      in_XZ7 & in_XZ8 ~ "Shared",
      in_XZ7 ~ "XZ7 only",
      in_XZ8 ~ "XZ8 only"
    ),
    Status = factor(Status, levels = status_levels)
  )

annotated <- df %>%
  left_join(
    variant_status,
    by = c("CHROM", "POS", "REF", "ALT")
  )

# One record per unique mutation
unique_mutations <- annotated %>%
  arrange(Strain) %>%
  distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)

lane_specific <- unique_mutations %>%
  filter(Status %in% c("XZ7 only", "XZ8 only")) %>%
  mutate(Status = droplevels(Status))

write.csv(
  unique_mutations,
  file.path(table_dir, "All_unique_mutations.csv"),
  row.names = FALSE
)

write.csv(
  lane_specific,
  file.path(table_dir, "Lane_specific_mutations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. MUTATION OVERLAP
# ------------------------------------------------------------

overlap_summary <- variant_status %>%
  count(Status, .drop = FALSE)

n_shared <- sum(variant_status$Status == "Shared")
n_XZ7 <- sum(variant_status$in_XZ7)
n_XZ8 <- sum(variant_status$in_XZ8)
n_total <- nrow(variant_status)

jaccard <- if (n_total > 0) {
  n_shared / n_total
} else {
  NA_real_
}

overlap_metrics <- tibble(
  Metric = c(
    "XZ7 mutations",
    "XZ8 mutations",
    "Shared mutations",
    "Unique mutations",
    "Jaccard similarity"
  ),
  Value = c(
    n_XZ7, n_XZ8, n_shared, n_total, jaccard
  )
)

write.csv(
  overlap_metrics,
  file.path(table_dir, "Overlap_metrics.csv"),
  row.names = FALSE
)

p_overlap <- ggplot(
  overlap_summary,
  aes(x = Status, y = n, fill = Status)
) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = lane_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Mutation overlap between XZ7 and XZ8",
    subtitle = paste(
      "Jaccard similarity =",
      round(jaccard, 3)
    ),
    x = NULL,
    y = "Number of mutations"
  )

save_graph(p_overlap, "01_Mutation_overlap")

# ------------------------------------------------------------
# 5. NUMERICAL FEATURES
# ------------------------------------------------------------

numeric_cols <- c(
  "GC_11 bp",
  "GC_101 bp",
  "GC_1001 bp",
  "GC_10001 bp",
  "GC_100001 bp",
  "recomb_rate",
  "distance_to_ORI"
)

numeric_cols <- intersect(
  numeric_cols,
  names(unique_mutations)
)

# Convert numerical columns
unique_mutations <- unique_mutations %>%
  mutate(
    across(
      all_of(numeric_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

lane_specific <- unique_mutations %>%
  filter(Status %in% c("XZ7 only", "XZ8 only")) %>%
  mutate(Status = droplevels(Status))

# Wilcoxon rank-sum tests
numeric_results <- purrr::map_dfr(
  numeric_cols,
  function(feature) {
    
    dat <- lane_specific %>%
      transmute(
        Status,
        Value = .data[[feature]]
      ) %>%
      filter(is.finite(Value))
    
    x <- dat$Value[dat$Status == "XZ7 only"]
    y <- dat$Value[dat$Status == "XZ8 only"]
    
    p <- if (length(x) > 0 && length(y) > 0) {
      tryCatch(
        wilcox.test(
          x, y,
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }
    
    tibble(
      Feature = feature,
      Test = "Wilcoxon",
      N_XZ7 = length(x),
      N_XZ8 = length(y),
      Median_XZ7 = if (length(x)) median(x) else NA_real_,
      Median_XZ8 = if (length(y)) median(y) else NA_real_,
      P_value = p
    )
  }
)

# Numerical plots
numeric_long <- unique_mutations %>%
  select(Status, all_of(numeric_cols)) %>%
  pivot_longer(
    cols = -Status,
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  filter(is.finite(Value))

p_numeric <- ggplot(
  numeric_long,
  aes(x = Status, y = Value, fill = Status)
) +
  geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.25
  ) +
  facet_wrap(
    ~Feature,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(values = lane_colors) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Numerical genomic features",
    x = NULL,
    y = "Feature value"
  )

save_graph(
  p_numeric,
  "02_Numerical_features",
  width = 13,
  height = 10
)

# Individual numerical plots
for (feature in numeric_cols) {
  
  plot_data <- unique_mutations %>%
    filter(is.finite(.data[[feature]]))
  
  p <- ggplot(
    plot_data,
    aes(
      x = Status,
      y = .data[[feature]],
      fill = Status
    )
  ) +
    geom_boxplot(
      width = 0.6,
      outlier.alpha = 0.3
    ) +
    scale_fill_manual(values = lane_colors) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none") +
    labs(
      title = feature,
      x = NULL,
      y = feature
    )
  
  filename <- paste0(
    "Numeric_",
    gsub("[^A-Za-z0-9]+", "_", feature)
  )
  
  save_graph(p, filename, 7, 5)
}

# ------------------------------------------------------------
# 6. CATEGORICAL FEATURES
# ------------------------------------------------------------

# REF, ALT, and Mut_type are excluded from graphs.
# Mutation type is retained for statistical testing.

categorical_test_cols <- c(
  "Mut_type (SNV/indel)",
  "CHROM",
  "Coding context",
  "mononuc?",
  "lead/lag",
  "putative effect",
  "gene_expressed_germline",
  "motif_3bp_mutfirst",
  "motif_3bp_mutlast",
  "motif_3bp_mutcenter"
)

categorical_plot_cols <- setdiff(
  categorical_test_cols,
  c("REF", "ALT", "Mut_type (SNV/indel)")
)

categorical_test_cols <- intersect(
  categorical_test_cols,
  names(unique_mutations)
)

categorical_plot_cols <- intersect(
  categorical_plot_cols,
  names(unique_mutations)
)

# Fisher's exact tests
categorical_results <- purrr::map_dfr(
  categorical_test_cols,
  function(feature) {
    
    dat <- lane_specific %>%
      transmute(
        Status,
        Category = as.character(.data[[feature]])
      ) %>%
      filter(!is.na(Category))
    
    tab <- table(dat$Status, dat$Category)
    
    p <- if (nrow(tab) == 2 && ncol(tab) >= 2) {
      tryCatch(
        fisher.test(
          tab,
          simulate.p.value = TRUE,
          B = 99999
        )$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }
    
    tibble(
      Feature = feature,
      Test = "Fisher",
      Categories = ncol(tab),
      P_value = p
    )
  }
)

# Plot each categorical feature separately
for (feature in categorical_plot_cols) {
  
  dat <- unique_mutations %>%
    transmute(
      Status,
      Category = as.character(.data[[feature]])
    ) %>%
    filter(!is.na(Category)) %>%
    count(Status, Category) %>%
    group_by(Status) %>%
    mutate(
      Proportion = n / sum(n)
    ) %>%
    ungroup()
  
  p <- ggplot(
    dat,
    aes(
      x = reorder(Category, -Proportion),
      y = Proportion,
      fill = Status
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    scale_fill_manual(values = lane_colors) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1
      ),
      legend.position = "top"
    ) +
    labs(
      title = feature,
      x = NULL,
      y = "Proportion within detection group",
      fill = "Detection status"
    )
  
  filename <- paste0(
    "Categorical_",
    gsub("[^A-Za-z0-9]+", "_", feature)
  )
  
  save_graph(
    p,
    filename,
    width = 12,
    height = 7
  )
}

# ------------------------------------------------------------
# 7. SEQUENCE COMPLEXITY: SHANNON ENTROPY
# ------------------------------------------------------------

entropy_dna <- function(sequence) {
  
  if (is.na(sequence)) return(NA_real_)
  
  bases <- strsplit(
    toupper(sequence), ""
  )[[1]]
  
  bases <- bases[
    bases %in% c("A", "C", "G", "T")
  ]
  
  if (length(bases) == 0) return(NA_real_)
  
  p <- as.numeric(table(bases)) / length(bases)
  
  -sum(p * log2(p))
}

seq_col <- "seq_complexity (99 bp, mut center)"

if (seq_col %in% names(unique_mutations)) {
  
  unique_mutations <- unique_mutations %>%
    mutate(
      Shannon_entropy = vapply(
        .data[[seq_col]],
        entropy_dna,
        numeric(1)
      )
    )
  
  entropy_specific <- unique_mutations %>%
    filter(Status %in% c("XZ7 only", "XZ8 only"))
  
  x <- entropy_specific$Shannon_entropy[
    entropy_specific$Status == "XZ7 only"
  ]
  y <- entropy_specific$Shannon_entropy[
    entropy_specific$Status == "XZ8 only"
  ]
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  entropy_p <- if (length(x) > 0 && length(y) > 0) {
    wilcox.test(x, y, exact = FALSE)$p.value
  } else {
    NA_real_
  }
  
  entropy_result <- tibble(
    Feature = "Shannon_entropy",
    Test = "Wilcoxon",
    N_XZ7 = length(x),
    N_XZ8 = length(y),
    Median_XZ7 = if (length(x)) median(x) else NA_real_,
    Median_XZ8 = if (length(y)) median(y) else NA_real_,
    P_value = entropy_p
  )
  
  numeric_results <- bind_rows(
    numeric_results,
    entropy_result
  )
  
  p_entropy <- ggplot(
    unique_mutations,
    aes(
      x = Status,
      y = Shannon_entropy,
      fill = Status
    )
  ) +
    geom_boxplot(outlier.alpha = 0.25) +
    scale_fill_manual(values = lane_colors) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none") +
    labs(
      title = "Sequence complexity around mutations",
      x = NULL,
      y = "Shannon entropy (bits)"
    )
  
  save_graph(
    p_entropy,
    "03_Shannon_entropy",
    8, 6
  )
}

# ------------------------------------------------------------
# 8. SIX-CLASS SNV MUTATION SPECTRUM
# ------------------------------------------------------------

complement <- c(
  A = "T",
  T = "A",
  C = "G",
  G = "C"
)

snv <- unique_mutations %>%
  filter(
    `Mut_type (SNV/indel)` == "SNV",
    REF %in% names(complement),
    ALT %in% names(complement),
    REF != ALT
  ) %>%
  mutate(
    REF_norm = if_else(
      REF %in% c("A", "G"),
      unname(complement[REF]),
      REF
    ),
    ALT_norm = if_else(
      REF %in% c("A", "G"),
      unname(complement[ALT]),
      ALT
    ),
    Substitution = paste0(
      REF_norm, ">", ALT_norm
    )
  )

substitution_levels <- c(
  "C>A", "C>G", "C>T",
  "T>A", "T>C", "T>G"
)

spectrum <- snv %>%
  count(Status, Substitution) %>%
  complete(
    Status = factor(
      status_levels,
      levels = status_levels
    ),
    Substitution = substitution_levels,
    fill = list(n = 0)
  ) %>%
  group_by(Status) %>%
  mutate(
    Proportion = if (sum(n) > 0) {
      n / sum(n)
    } else {
      rep(NA_real_, n())
    }
  ) %>%
  ungroup()

write.csv(
  spectrum,
  file.path(table_dir, "SNV_spectrum.csv"),
  row.names = FALSE
)

p_spectrum <- ggplot(
  spectrum,
  aes(
    x = Substitution,
    y = Proportion,
    fill = Status
  )
) +
  geom_col(
    position = "dodge",
    width = 0.75
  ) +
  scale_fill_manual(values = lane_colors) +
  theme_classic(base_size = 13) +
  labs(
    title = "Six-class SNV mutation spectrum",
    x = "Substitution type",
    y = "Proportion",
    fill = "Detection status"
  )

save_graph(
  p_spectrum,
  "04_SNV_spectrum",
  9, 6
)

# Fisher test for lane-specific SNVs
snv_specific <- snv %>%
  filter(Status %in% c("XZ7 only", "XZ8 only")) %>%
  mutate(Status = droplevels(Status))

tab_snv <- table(
  snv_specific$Status,
  snv_specific$Substitution
)

snv_p <- if (
  nrow(tab_snv) == 2 &&
  ncol(tab_snv) >= 2
) {
  fisher.test(
    tab_snv,
    simulate.p.value = TRUE,
    B = 99999
  )$p.value
} else {
  NA_real_
}

snv_result <- tibble(
  Feature = "Six_class_SNV_spectrum",
  Test = "Fisher",
  Categories = ncol(tab_snv),
  P_value = snv_p
)

categorical_results <- bind_rows(
  categorical_results,
  snv_result
)

# ------------------------------------------------------------
# 9. EXTENDED MOTIF FREQUENCY COMPARISONS
# ------------------------------------------------------------

motif_cols <- c(
  "motif_5bp_center",
  "motif_7bp_center",
  "motif_9bp_center"
)

motif_cols <- intersect(
  motif_cols,
  names(unique_mutations)
)

# Only display the 20 most common motifs per feature.
# Restrict to SNVs for sequence-context comparisons.

snv_motifs <- unique_mutations %>%
  filter(`Mut_type (SNV/indel)` == "SNV")

for (feature in motif_cols) {
  
  dat <- snv_motifs %>%
    transmute(
      Status,
      Motif = as.character(.data[[feature]])
    ) %>%
    filter(!is.na(Motif))
  
  top_motifs <- dat %>%
    count(Motif, sort = TRUE) %>%
    slice_head(n = 20) %>%
    pull(Motif)
  
  plot_data <- dat %>%
    count(Status, Motif) %>%
    group_by(Status) %>%
    mutate(
      Proportion = n / sum(n)
    ) %>%
    ungroup() %>%
    filter(Motif %in% top_motifs)
  
  p <- ggplot(
    plot_data,
    aes(
      x = reorder(Motif, Proportion),
      y = Proportion,
      fill = Status
    )
  ) +
    geom_col(
      position = "dodge",
      width = 0.75
    ) +
    coord_flip() +
    scale_fill_manual(values = lane_colors) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Top 20 motifs:", feature),
      x = "Sequence motif",
      y = "Proportion of all SNVs in group",
      fill = "Detection status"
    )
  
  save_graph(
    p,
    paste0("Motif_", feature),
    width = 11,
    height = 9
  )
  
  write.csv(
    plot_data,
    file.path(
      table_dir,
      paste0("Motif_", feature, ".csv")
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 10. FINAL STATISTICAL RESULTS
# ------------------------------------------------------------

all_results <- bind_rows(
  numeric_results %>%
    select(Feature, Test, P_value),
  
  categorical_results %>%
    select(Feature, Test, P_value)
) %>%
  mutate(
    FDR = p.adjust(P_value, method = "BH"),
    Significant = case_when(
      is.na(FDR) ~ "Not tested",
      FDR < 0.05 ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  arrange(FDR)

write.csv(
  numeric_results,
  file.path(table_dir, "Numeric_results.csv"),
  row.names = FALSE
)

write.csv(
  categorical_results,
  file.path(table_dir, "Categorical_results.csv"),
  row.names = FALSE
)

write.csv(
  all_results,
  file.path(table_dir, "All_statistical_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. VISUALIZE FEATURE-LEVEL P-VALUES
# ------------------------------------------------------------

p_values <- all_results %>%
  filter(!is.na(FDR)) %>%
  mutate(
    Feature = forcats::fct_reorder(
      Feature,
      -log10(pmax(FDR, .Machine$double.xmin))
    ),
    LogFDR = -log10(pmax(FDR, .Machine$double.xmin))
  )

p_fdr <- ggplot(
  p_values,
  aes(
    x = Feature,
    y = LogFDR,
    fill = Significant
  )
) +
  geom_col(width = 0.7) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Yes" = "#E28E65",
      "No" = "#4C78A8"
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Exploratory genomic feature comparisons",
    subtitle = "XZ7-only versus XZ8-only mutations",
    x = NULL,
    y = expression(-log[10]("FDR")),
    fill = "FDR < 0.05"
  )

save_graph(
  p_fdr,
  "05_Feature_FDR_comparison",
  width = 11,
  height = 9
)

# ------------------------------------------------------------
# 12. PRINT SUMMARY
# ------------------------------------------------------------

cat("\n====================================\n")
cat("XZ7 vs XZ8 ANALYSIS COMPLETE\n")
cat("====================================\n\n")

cat("XZ7 mutations:", n_XZ7, "\n")
cat("XZ8 mutations:", n_XZ8, "\n")
cat("Shared mutations:", n_shared, "\n")
cat("Unique mutations:", n_total, "\n")
cat("Jaccard similarity:", round(jaccard, 4), "\n\n")

cat("Feature-level results:\n")
print(all_results, n = Inf)

cat("\nGraphs saved in:\n")
cat(normalizePath(graph_dir), "\n")

cat("\nTables saved in:\n")
cat(normalizePath(table_dir), "\n")

cat("\nDone!\n")
