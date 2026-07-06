# Plot total normalized miRNA counts per sample
# Project: TCGA-miRNA-biogenesis-analysis
#
# This script calculates the total normalized miRNA count for each sample
# by summing normalized miRNA counts across all detected miRNAs.
#
# For each cancer type:
# Normal and tumor samples are compared using Welch's t-test.
#
# Input:
# data/processed/miRNA/normalized_counts/
# data/processed/miRNA/sample_annotation/
#
# Outputs:
# data/processed/miRNA/total_normalized_counts/
# results/tables/miRNA/tcga_miRNA_total_normalized_count_statistics.csv
# results/figures/miRNA/tcga_miRNA_total_normalized_counts_boxplot.png
# results/figures/miRNA/tcga_miRNA_total_normalized_counts_boxplot.pdf
# results/figures/miRNA/tcga_miRNA_total_normalized_counts_boxplot.tiff


# Required libraries
#_________________________________________________________________

library(ggplot2)
library(dplyr)
library(magrittr)


# Define TCGA cancer types
#_________________________________________________________________

cancer_types <- c(
  "COAD", "ESCA", "BLCA", "BRCA", "LUAD", "LUSC",
  "KIRC", "LIHC", "THCA", "UCEC", "GBM", "STAD", "PRAD"
)


# Define input directories
#_________________________________________________________________

normalized_counts_dir <- file.path(
  getwd(),
  "data",
  "processed",
  "miRNA",
  "normalized_counts"
)

sample_annotation_dir <- file.path(
  getwd(),
  "data",
  "processed",
  "miRNA",
  "sample_annotation"
)


# Define output directories
#_________________________________________________________________

total_counts_dir <- file.path(
  getwd(),
  "data",
  "processed",
  "miRNA",
  "total_normalized_counts"
)

tables_dir <- file.path(
  getwd(),
  "results",
  "tables",
  "miRNA"
)

figures_dir <- file.path(
  getwd(),
  "results",
  "figures",
  "miRNA"
)

dir.create(total_counts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# Check input directories
#_________________________________________________________________

if (!dir.exists(normalized_counts_dir)) {
  stop("Normalized miRNA count directory was not found.")
}

if (!dir.exists(sample_annotation_dir)) {
  stop("Sample annotation directory was not found.")
}


# Function to convert p-value to significance label
#_________________________________________________________________

get_significance_label <- function(p_value) {
  
  if (is.na(p_value)) {
    return("NA")
  }
  
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}


# Empty objects for combined results
#_________________________________________________________________

all_total_counts <- data.frame()
all_statistics <- data.frame()


# Process each cancer type
#_________________________________________________________________

for (cancer in cancer_types) {
  
  message("------------------------------------------------------------")
  message("Processing total normalized miRNA counts for: ", cancer)
  message("------------------------------------------------------------")
  
  normalized_counts_file <- file.path(
    normalized_counts_dir,
    paste0("normalized_miRNA_counts_log2_", cancer, ".csv")
  )
  
  sample_annotation_file <- file.path(
    sample_annotation_dir,
    paste0("sample_annotation_miRNA_", cancer, ".csv")
  )
  
  if (!file.exists(normalized_counts_file)) {
    warning("Normalized count file was not found for ", cancer)
    next
  }
  
  if (!file.exists(sample_annotation_file)) {
    warning("Sample annotation file was not found for ", cancer)
    next
  }
  
  
  # Read normalized count matrix and sample annotation
  #_________________________________________________________________
  
  normalized_counts <- read.csv(
    normalized_counts_file,
    row.names = 1,
    check.names = FALSE
  )
  
  sample_annotation <- read.csv(
    sample_annotation_file,
    check.names = FALSE
  )
  
  
  # Keep only common samples
  #_________________________________________________________________
  
  common_samples <- intersect(
    colnames(normalized_counts),
    sample_annotation$Sample
  )
  
  if (length(common_samples) == 0) {
    warning("No common samples were found between count matrix and annotation for ", cancer)
    next
  }
  
  normalized_counts <- normalized_counts[
    ,
    common_samples,
    drop = FALSE
  ]
  
  sample_annotation <- sample_annotation[
    match(common_samples, sample_annotation$Sample),
    ,
    drop = FALSE
  ]
  
  
  # Keep only Normal and Tumor samples
  #_________________________________________________________________
  
  keep_samples <- sample_annotation$Group %in% c("Normal", "Tumor")
  
  normalized_counts <- normalized_counts[
    ,
    keep_samples,
    drop = FALSE
  ]
  
  sample_annotation <- sample_annotation[
    keep_samples,
    ,
    drop = FALSE
  ]
  
  sample_annotation$Group <- factor(
    sample_annotation$Group,
    levels = c("Normal", "Tumor")
  )
  
  
  # Calculate total normalized miRNA count per sample
  #_________________________________________________________________
  
  total_counts <- colSums(
    normalized_counts,
    na.rm = TRUE
  )
  
  total_counts_df <- data.frame(
    Cancer = cancer,
    Sample = colnames(normalized_counts),
    Group = sample_annotation$Group,
    Total_normalized_miRNA_count = as.numeric(total_counts),
    stringsAsFactors = FALSE
  )
  
  
  # Save cancer-specific total count table
  #_________________________________________________________________
  
  write.csv(
    total_counts_df,
    file = file.path(
      total_counts_dir,
      paste0("total_normalized_miRNA_counts_", cancer, ".csv")
    ),
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # Statistical comparison: Normal vs Tumor
  #_________________________________________________________________
  
  normal_values <- total_counts_df$Total_normalized_miRNA_count[
    total_counts_df$Group == "Normal"
  ]
  
  tumor_values <- total_counts_df$Total_normalized_miRNA_count[
    total_counts_df$Group == "Tumor"
  ]
  
  if (length(normal_values) >= 2 && length(tumor_values) >= 2) {
    
    test_result <- t.test(
      normal_values,
      tumor_values
    )
    
    p_value <- test_result$p.value
    
  } else {
    
    p_value <- NA_real_
  }
  
  significance_label <- get_significance_label(p_value)
  
  cancer_statistics <- data.frame(
    Cancer = cancer,
    Normal_samples = length(normal_values),
    Tumor_samples = length(tumor_values),
    Mean_Normal = mean(normal_values, na.rm = TRUE),
    Mean_Tumor = mean(tumor_values, na.rm = TRUE),
    Median_Normal = median(normal_values, na.rm = TRUE),
    Median_Tumor = median(tumor_values, na.rm = TRUE),
    P_value = p_value,
    Significance = significance_label,
    Test = "Welch t-test",
    stringsAsFactors = FALSE
  )
  
  all_total_counts <- rbind(
    all_total_counts,
    total_counts_df
  )
  
  all_statistics <- rbind(
    all_statistics,
    cancer_statistics
  )
}


# Save combined total count and statistics tables
#_________________________________________________________________

write.csv(
  all_total_counts,
  file = file.path(
    total_counts_dir,
    "total_normalized_miRNA_counts_all_cancers.csv"
  ),
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  all_statistics,
  file = file.path(
    tables_dir,
    "tcga_miRNA_total_normalized_count_statistics.csv"
  ),
  quote = FALSE,
  row.names = FALSE
)


# Prepare significance annotation positions for plotting
#_________________________________________________________________

plot_statistics <- all_total_counts %>%
  group_by(Cancer) %>%
  summarise(
    y_max = max(Total_normalized_miRNA_count, na.rm = TRUE),
    y_min = min(Total_normalized_miRNA_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    all_statistics,
    by = "Cancer"
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_position = y_max + 0.08 * y_range,
    y_text = y_max + 0.13 * y_range,
    y_tick = y_max + 0.04 * y_range
  )


# Set cancer order
#_________________________________________________________________

all_total_counts$Cancer <- factor(
  all_total_counts$Cancer,
  levels = cancer_types
)

plot_statistics$Cancer <- factor(
  plot_statistics$Cancer,
  levels = cancer_types
)

all_total_counts$Group <- factor(
  all_total_counts$Group,
  levels = c("Normal", "Tumor")
)


# Generate faceted boxplot
#_________________________________________________________________

p <- ggplot(
  all_total_counts,
  aes(
    x = Group,
    y = Total_normalized_miRNA_count,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.95
  ) +
  geom_jitter(
    width = 0.12,
    size = 0.7,
    alpha = 0.45,
    color = "black"
  ) +
  geom_segment(
    data = plot_statistics,
    aes(
      x = 1,
      xend = 2,
      y = y_position,
      yend = y_position
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = plot_statistics,
    aes(
      x = 1,
      xend = 1,
      y = y_tick,
      yend = y_position
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = plot_statistics,
    aes(
      x = 2,
      xend = 2,
      y = y_tick,
      yend = y_position
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_text(
    data = plot_statistics,
    aes(
      x = 1.5,
      y = y_text,
      label = Significance
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  facet_wrap(
    ~ Cancer,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#1f78b4",
      "Tumor" = "#e31a1c"
    ),
    name = "Sample"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  labs(
    x = NULL,
    y = "Total normalized miRNA count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.major = element_line(
      linewidth = 0.25,
      color = "grey90"
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      size = 8
    ),
    axis.text.y = element_text(
      size = 8
    ),
    axis.title.y = element_text(
      size = 10
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 10
    )
  )


# Save figure
#_________________________________________________________________

ggsave(
  filename = file.path(
    figures_dir,
    "tcga_miRNA_total_normalized_counts_boxplot.png"
  ),
  plot = p,
  width = 9,
  height = 11,
  dpi = 300
)

ggsave(
  filename = file.path(
    figures_dir,
    "tcga_miRNA_total_normalized_counts_boxplot.pdf"
  ),
  plot = p,
  width = 9,
  height = 11
)

ggsave(
  filename = file.path(
    figures_dir,
    "tcga_miRNA_total_normalized_counts_boxplot.tiff"
  ),
  plot = p,
  width = 9,
  height = 11,
  dpi = 300,
  compression = "lzw"
)


