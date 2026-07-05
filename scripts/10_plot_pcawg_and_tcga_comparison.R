# ============================================================
# Script 10: Plot PCAWG cBioPortal alteration frequency data
# and compare PCAWG with TCGA
#
# Project: TCGA-miRNA-biogenesis-analysis
#
# Purpose:
#   Generate PCAWG cBioPortal alteration frequency plots and
#   compare overall gene alteration frequency between PCAWG and
#   TCGA using comparable cancer categories.
#
# Inputs:
#   data/processed/pcawg_alteration_frequency_long.csv
#   data/processed/cbioportal_alteration_frequency_long.csv
#
# Outputs:
#   results/figures/pcawg_gene_alteration_counts_stacked_bar.*
#   results/figures/pcawg_alteration_type_pie_chart.*
#   results/figures/pcawg_gene_alteration_bar_and_pie.*
#   results/figures/pcawg_tcga_back_to_back_gene_alteration_frequency.*
#
#   results/tables/tcga_to_pcawg_cancer_mapping_check.csv
#   results/tables/tcga_pcawg_common_cancer_sample_size_check.csv
#   results/tables/tcga_pcawg_common_cancer_gene_frequency_summary.csv
#   results/tables/tcga_pcawg_common_gene_frequency_summary.csv
# ============================================================


# Load packages
#_________________________________________________________________

required_packages <- c("tidyverse")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' is required. Please install it first."))
  }
}

library(tidyverse)

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)


# Define paths
#_________________________________________________________________

pcawg_file <- "data/processed/pcawg_alteration_frequency_long.csv"
tcga_file <- "data/processed/cbioportal_alteration_frequency_long.csv"

figures_dir <- "results/figures"
tables_dir <- "results/tables"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pcawg_file)) {
  stop("PCAWG processed file was not found: ", pcawg_file)
}

if (!file.exists(tcga_file)) {
  stop("TCGA processed file was not found: ", tcga_file)
}


# Define orders
#_________________________________________________________________

gene_order <- c(
  "AGO1", "AGO2", "AGO3", "AGO4",
  "DDX17", "DDX20", "DDX5",
  "DGCR8", "DICER1", "DROSHA",
  "GEMIN4", "PRKRA", "TARBP2",
  "TNRC6A", "XPO1", "XPO5"
)

selected_pcawg_cancers <- c(
  "Bladder Cancer",
  "Breast Cancer",
  "Colorectal Cancer",
  "Non-Small Cell Lung Cancer",
  "Uterine Endometrioid Carcinoma",
  "Hepatobiliary Cancer",
  "Esophagogastric Cancer",
  "Renal Cell Carcinoma",
  "Prostate Cancer",
  "Thyroid Cancer"
)

alteration_type_order <- c(
  "Amplification",
  "Deep deletion",
  "Multiple alterations",
  "Point mutation"
)


# Colors
# Same colors as TCGA cBioPortal plotting script
#_________________________________________________________________

alteration_colors <- c(
  "Amplification" = "#40B3A2",
  "Deep deletion" = "#F2C94C",
  "Multiple alterations" = "#2D9CDB",
  "Point mutation" = "#EB5757"
)

study_colors <- c(
  "PCAWG" = "#40B3A2",
  "TCGA" = "#F2C94C"
)


# Helper function to save ggplot in 3 formats
#_________________________________________________________________

save_ggplot_all_formats <- function(plot_object, file_stem, width, height, dpi = 300) {
  
  png_file <- file.path(figures_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(figures_dir, paste0(file_stem, ".pdf"))
  tiff_file <- file.path(figures_dir, paste0(file_stem, ".tiff"))
  
  ggsave(
    filename = png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    filename = tiff_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )
}


# Helper functions for labels
#_________________________________________________________________

format_percent_label <- function(x) {
  ifelse(
    abs(x - round(x)) < 0.05,
    paste0(round(x), "%"),
    paste0(round(x, 1), "%")
  )
}

format_number_label <- function(x) {
  ifelse(
    abs(x - round(x)) < 0.05,
    as.character(round(x)),
    as.character(round(x, 1))
  )
}


# Helper function to normalize alteration type names
#_________________________________________________________________

normalize_alteration_type <- function(x) {
  
  x_clean <- str_to_lower(str_trim(as.character(x)))
  
  case_when(
    x_clean %in% c("amp", "amplification") ~ "Amplification",
    x_clean %in% c("homdel", "deep deletion", "deep_deletion") ~ "Deep deletion",
    x_clean %in% c("multiple", "multiple alterations") ~ "Multiple alterations",
    x_clean %in% c("mutated", "mutation", "point mutation", "point mutations") ~ "Point mutation",
    TRUE ~ as.character(x)
  )
}


# Helper function to standardize TCGA and PCAWG tables
#_________________________________________________________________

normalize_cbioportal_table <- function(df, dataset_name) {
  
  names(df) <- names(df) %>%
    str_replace_all("\\s+", "_")
  
  cancer_col <- intersect(
    c("Cancer", "Cancer_Type", "Cancer_Study"),
    names(df)
  )
  
  if (length(cancer_col) == 0) {
    stop(
      "No cancer column found in ", dataset_name,
      ". Expected one of: Cancer, Cancer_Type, Cancer_Study."
    )
  }
  
  required_cols <- c(
    "Gene",
    "Alteration_Frequency",
    "Alteration_Type",
    "Alteration_Count"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", dataset_name, " table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  sample_col <- intersect(
    c("Total_Samples", "Total_Samples_Inferred", "Inferred_Total_Samples"),
    names(df)
  )
  
  out <- df %>%
    transmute(
      Dataset = dataset_name,
      Gene = as.character(.data[["Gene"]]),
      Cancer_Label = as.character(.data[[cancer_col[1]]]),
      Alteration_Type = normalize_alteration_type(.data[["Alteration_Type"]]),
      Alteration_Frequency = readr::parse_number(as.character(.data[["Alteration_Frequency"]])),
      Alteration_Count = readr::parse_number(as.character(.data[["Alteration_Count"]]))
    )
  
  if (length(sample_col) > 0) {
    out$Total_Samples_Input <- readr::parse_number(as.character(df[[sample_col[1]]]))
  } else {
    out$Total_Samples_Input <- NA_real_
  }
  
  out <- out %>%
    mutate(
      Row_Level_Total_Samples = if_else(
        Alteration_Frequency > 0 & Alteration_Count > 0,
        Alteration_Count / (Alteration_Frequency / 100),
        NA_real_
      ),
      Total_Samples_Candidate = coalesce(
        Total_Samples_Input,
        Row_Level_Total_Samples
      )
    )
  
  sample_size_lookup <- out %>%
    filter(!is.na(Total_Samples_Candidate)) %>%
    group_by(Dataset, Cancer_Label) %>%
    summarise(
      Total_Samples_Inferred = round(median(Total_Samples_Candidate, na.rm = TRUE)),
      .groups = "drop"
    )
  
  out <- out %>%
    left_join(
      sample_size_lookup,
      by = c("Dataset", "Cancer_Label")
    ) %>%
    select(
      Dataset,
      Gene,
      Cancer_Label,
      Alteration_Type,
      Alteration_Frequency,
      Alteration_Count,
      Total_Samples_Inferred
    )
  
  out
}


# Helper function to map TCGA cancer codes to PCAWG-comparable cancer groups
#_________________________________________________________________

map_tcga_to_common_cancer <- function(x) {
  
  x_clean <- str_to_lower(as.character(x))
  
  case_when(
    str_detect(x_clean, "bladder|blca") ~ "Bladder Cancer",
    str_detect(x_clean, "breast|brca") ~ "Breast Cancer",
    str_detect(x_clean, "colon|colorectal|coad|read") ~ "Colorectal Cancer",
    str_detect(x_clean, "lung adenocarcinoma|lung squamous|non-small cell lung|luad|lusc") ~ "Non-Small Cell Lung Cancer",
    str_detect(x_clean, "uterine|endometrioid|endometrial|ucec") ~ "Uterine Endometrioid Carcinoma",
    str_detect(x_clean, "liver|hepatocellular|hepatobiliary|lihc") ~ "Hepatobiliary Cancer",
    str_detect(x_clean, "esophageal|esophagus|gastric|stomach|esophagogastric|esca|stad") ~ "Esophagogastric Cancer",
    str_detect(x_clean, "kidney|renal|kirc") ~ "Renal Cell Carcinoma",
    str_detect(x_clean, "prostate|prad") ~ "Prostate Cancer",
    str_detect(x_clean, "thyroid|thca") ~ "Thyroid Cancer",
    TRUE ~ NA_character_
  )
}


# Read processed data
#_________________________________________________________________

pcawg_raw <- readr::read_csv(pcawg_file, show_col_types = FALSE)
tcga_raw <- readr::read_csv(tcga_file, show_col_types = FALSE)

pcawg_data <- normalize_cbioportal_table(pcawg_raw, "PCAWG")
tcga_data <- normalize_cbioportal_table(tcga_raw, "TCGA")


# Check alteration type names
#_________________________________________________________________

unexpected_pcawg_types <- setdiff(
  unique(pcawg_data$Alteration_Type),
  alteration_type_order
)

unexpected_tcga_types <- setdiff(
  unique(tcga_data$Alteration_Type),
  alteration_type_order
)

if (length(unexpected_pcawg_types) > 0) {
  stop(
    "Unexpected PCAWG alteration types found: ",
    paste(unexpected_pcawg_types, collapse = ", ")
  )
}

if (length(unexpected_tcga_types) > 0) {
  stop(
    "Unexpected TCGA alteration types found: ",
    paste(unexpected_tcga_types, collapse = ", ")
  )
}


# Keep selected PCAWG cancer types
#_________________________________________________________________

pcawg_common <- pcawg_data %>%
  filter(
    Gene %in% gene_order,
    Cancer_Label %in% selected_pcawg_cancers
  ) %>%
  mutate(
    Common_Cancer = Cancer_Label,
    Gene = factor(Gene, levels = gene_order),
    Common_Cancer = factor(Common_Cancer, levels = selected_pcawg_cancers),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )


# Map TCGA cancer types to PCAWG-comparable groups
#_________________________________________________________________

tcga_common <- tcga_data %>%
  filter(Gene %in% gene_order) %>%
  mutate(
    Common_Cancer = map_tcga_to_common_cancer(Cancer_Label)
  ) %>%
  filter(!is.na(Common_Cancer)) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Common_Cancer = factor(Common_Cancer, levels = selected_pcawg_cancers),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )


# Save TCGA to PCAWG cancer mapping check table
#_________________________________________________________________

tcga_mapping_check <- tcga_data %>%
  distinct(Cancer_Label) %>%
  mutate(
    Mapped_Common_Cancer = map_tcga_to_common_cancer(Cancer_Label),
    Included_in_Comparison = !is.na(Mapped_Common_Cancer)
  ) %>%
  arrange(desc(Included_in_Comparison), Cancer_Label)

readr::write_csv(
  tcga_mapping_check,
  file.path(tables_dir, "tcga_to_pcawg_cancer_mapping_check.csv")
)


# Combine common data
#_________________________________________________________________

combined_common <- bind_rows(pcawg_common, tcga_common) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("PCAWG", "TCGA")),
    Gene = factor(Gene, levels = gene_order),
    Common_Cancer = factor(Common_Cancer, levels = selected_pcawg_cancers),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )


# Save sample size check table
#_________________________________________________________________

sample_size_check <- combined_common %>%
  distinct(
    Dataset,
    Cancer_Label,
    Common_Cancer,
    Total_Samples_Inferred
  ) %>%
  arrange(Dataset, Common_Cancer, Cancer_Label)

readr::write_csv(
  sample_size_check,
  file.path(tables_dir, "tcga_pcawg_common_cancer_sample_size_check.csv")
)

missing_sample_sizes <- sample_size_check %>%
  filter(is.na(Total_Samples_Inferred))

if (nrow(missing_sample_sizes) > 0) {
  warning(
    "Some cancer types have missing inferred sample sizes. ",
    "Check results/tables/tcga_pcawg_common_cancer_sample_size_check.csv"
  )
}


# ============================================================
# Plot 1:
# PCAWG stacked bar plot of alteration counts by gene
# Similar to previous TCGA stacked bar plot
# ============================================================

total_pcawg_samples <- pcawg_common %>%
  distinct(Cancer_Label, Total_Samples_Inferred) %>%
  summarise(
    Total_PCAWG_Samples = sum(Total_Samples_Inferred, na.rm = TRUE)
  ) %>%
  pull(Total_PCAWG_Samples)

if (is.na(total_pcawg_samples) || total_pcawg_samples == 0) {
  stop("Total PCAWG sample size is missing or zero. Please check PCAWG sample size inference.")
}

pcawg_gene_count_data <- pcawg_common %>%
  group_by(Gene, Alteration_Type) %>%
  summarise(
    Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )

pcawg_gene_label_data <- pcawg_common %>%
  group_by(Gene) %>%
  summarise(
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    Overall_Percent = 100 * Total_Alteration_Count / total_pcawg_samples,
    .groups = "drop"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Percent_Label = format_percent_label(Overall_Percent)
  )

p_pcawg_stacked_bar <- ggplot(
  pcawg_gene_count_data,
  aes(
    x = Gene,
    y = Alteration_Count,
    fill = Alteration_Type
  )
) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.15
  ) +
  geom_text(
    data = pcawg_gene_label_data,
    aes(
      x = Gene,
      y = Total_Alteration_Count,
      label = Percent_Label
    ),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = alteration_colors,
    drop = FALSE,
    name = "Condition"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "PCAWG",
    x = NULL,
    y = "Number of alterations"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

save_ggplot_all_formats(
  plot_object = p_pcawg_stacked_bar,
  file_stem = "pcawg_gene_alteration_counts_stacked_bar",
  width = 8,
  height = 5
)


# ============================================================
# Plot 2:
# PCAWG pie chart of alteration type proportions
# Similar to previous TCGA pie chart
# ============================================================

pcawg_pie_data <- pcawg_common %>%
  group_by(Alteration_Type) %>%
  summarise(
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  ) %>%
  arrange(Alteration_Type) %>%
  filter(Total_Alteration_Count > 0)

if (sum(pcawg_pie_data$Total_Alteration_Count, na.rm = TRUE) == 0) {
  stop("No PCAWG alterations were found for the selected cancer types.")
}

pcawg_pie_data <- pcawg_pie_data %>%
  mutate(
    Percent = 100 * Total_Alteration_Count / sum(Total_Alteration_Count),
    Label = format_percent_label(Percent),
    fraction = Total_Alteration_Count / sum(Total_Alteration_Count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    ymid = (ymin + ymax) / 2
  )

p_pcawg_pie <- ggplot(pcawg_pie_data) +
  geom_rect(
    aes(
      ymin = ymin,
      ymax = ymax,
      xmin = 0,
      xmax = 1,
      fill = Alteration_Type
    ),
    color = "black",
    linewidth = 0.3
  ) +
  coord_polar(theta = "y") +
  geom_label(
    aes(
      x = 1.18,
      y = ymid,
      label = Label,
      fill = Alteration_Type
    ),
    color = "black",
    fontface = "bold",
    size = 3.5,
    label.size = 0.25,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = alteration_colors,
    drop = FALSE
  ) +
  xlim(0, 1.35) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

save_ggplot_all_formats(
  plot_object = p_pcawg_pie,
  file_stem = "pcawg_alteration_type_pie_chart",
  width = 5,
  height = 5
)


# ============================================================
# Plot 3:
# Combined PCAWG stacked bar plot and pie chart
# ============================================================

if (has_patchwork) {
  
  p_pcawg_combined <- patchwork::wrap_plots(
    p_pcawg_stacked_bar,
    p_pcawg_pie,
    nrow = 1,
    widths = c(1.6, 1)
  )
  
  save_ggplot_all_formats(
    plot_object = p_pcawg_combined,
    file_stem = "pcawg_gene_alteration_bar_and_pie",
    width = 13,
    height = 5
  )
  
} else {
  message(
    "Package 'patchwork' is not installed. ",
    "The separate PCAWG bar and pie plots were saved, ",
    "but the combined plot was skipped."
  )
}


# ============================================================
# PCAWG vs TCGA comparison
# Overall gene alteration frequency across common cancer groups
# ============================================================

# ============================================================
# PCAWG vs TCGA comparison
# Overall gene alteration frequency across common cancer groups
# Exact tumor sample sizes are manually entered from cBioPortal
# Cancer Type Detailed summary files.
# No inferred sample sizes are used.
# ============================================================

# ============================================================
# Plot 4:
# Back-to-back PCAWG vs TCGA overall gene alteration frequency
# Single ggplot version using ggh4x
# This fixes:
#   1. duplicated/incomplete legends
#   2. overlapping zero labels
#   3. x-axis title position
# ============================================================

if (!requireNamespace("ggh4x", quietly = TRUE)) {
  stop(
    "Package 'ggh4x' is required for this back-to-back plot. ",
    "Please install it using: install.packages('ggh4x')"
  )
}

comparison_plot_data <- study_gene_summary %>%
  mutate(
    Dataset = factor(Dataset, levels = c("PCAWG", "TCGA")),
    Gene = factor(Gene, levels = gene_order),
    Frequency_Label = format_number_label(Overall_Alteration_Frequency)
  )

max_frequency <- max(
  comparison_plot_data$Overall_Alteration_Frequency,
  na.rm = TRUE
)

if (is.na(max_frequency) || max_frequency == 0) {
  stop(
    "Maximum comparison frequency is missing or zero. ",
    "Please check TCGA/PCAWG summary data."
  )
}

axis_limit <- ceiling((max_frequency * 1.15) / 5) * 5
label_offset <- axis_limit * 0.025

comparison_plot_data <- comparison_plot_data %>%
  mutate(
    Label_X = Overall_Alteration_Frequency + label_offset,
    Label_Hjust = if_else(Dataset == "PCAWG", 1, 0)
  )

p_comparison <- ggplot(
  comparison_plot_data,
  aes(
    x = Overall_Alteration_Frequency,
    y = forcats::fct_rev(Gene),
    fill = Dataset
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.15
  ) +
  geom_text(
    aes(
      x = Label_X,
      label = Frequency_Label,
      hjust = Label_Hjust
    ),
    size = 3.5,
    color = "black",
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Dataset,
    scales = "free_x",
    nrow = 1
  ) +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_reverse(
        limits = c(axis_limit, 0),
        labels = function(x) format_number_label(abs(x)),
        expand = expansion(mult = c(0, 0))
      ),
      scale_x_continuous(
        limits = c(0, axis_limit),
        labels = function(x) format_number_label(abs(x)),
        expand = expansion(mult = c(0, 0))
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      "PCAWG" = "#40B3A2",
      "TCGA" = "#F2C94C"
    ),
    breaks = c("PCAWG", "TCGA"),
    drop = FALSE,
    name = "Study"
  ) +
  labs(
    x = "Alteration frequency (%)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text.x = element_text(size = 12, color = "black", face = "bold"),
    
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black", face = "bold"),
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_blank(),
    
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    
    panel.grid.minor = element_blank(),
    
    # Increase this value if the two zero labels are still too close
    panel.spacing.x = grid::unit(1.2, "lines"),
    
    plot.margin = margin(5, 10, 5, 5)
  )

save_ggplot_all_formats(
  plot_object = p_comparison,
  file_stem = "pcawg_tcga_back_to_back_gene_alteration_frequency",
  width = 11.5,
  height = 6
)

