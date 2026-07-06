# ============================================================
# Script 11: TCGA expression-alteration integration plot
#
# Project: TCGA-miRNA-biogenesis-analysis
#
# Purpose:
#   Integrate TCGA cBioPortal alteration frequency data with
#   TCGA DESeq2 expression log2FC results.
#
#   This script generates one final multi-panel figure containing
#   all selected TCGA cancer types.
#
#   For each cancer type:
#     Top panel: stacked alteration frequency by gene
#     Bottom panel: DESeq2 log2FC heatmap tile by gene
#
# Important:
#   This is a cancer-level integration summary.
#   It is not a patient-level matched mutation-expression analysis.
#
# Inputs:
#   data/processed/cbioportal_alteration_frequency_long.csv
#   results/combined_results/combined_DESeq2_results_selected_16_genes_all_cancers.csv
#
# Outputs:
#   data/processed/tcga_expression_alteration_integration_long.csv
#   data/processed/tcga_expression_alteration_integration_summary.csv
#   results/tables/tcga_expression_alteration_missing_log2fc_check.csv
#   results/figures/integration/tcga_expression_alteration_all_cancers.png
#   results/figures/integration/tcga_expression_alteration_all_cancers.pdf
#   results/figures/integration/tcga_expression_alteration_all_cancers.tiff
# ============================================================


# Load packages
#_________________________________________________________________

required_packages <- c(
  "tidyverse",
  "patchwork",
  "cowplot"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' is required. Please install it first."))
  }
}

library(tidyverse)
library(patchwork)
library(cowplot)

# Define paths
#_________________________________________________________________

cbioportal_file <- "data/processed/cbioportal_alteration_frequency_long.csv"

expression_file <- "results/combined_results/combined_DESeq2_results_selected_16_genes_all_cancers.csv"

processed_dir <- "data/processed"
figures_dir <- "results/figures/integration"
tables_dir <- "results/tables"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cbioportal_file)) {
  stop("cBioPortal processed file was not found: ", cbioportal_file)
}

if (!file.exists(expression_file)) {
  stop("Expression log2FC file was not found: ", expression_file)
}


# Define cancer and gene orders
#_________________________________________________________________

cancer_order <- c(
  "COAD", "ESCA", "BLCA", "BRCA",
  "LUAD", "LUSC", "KIRC", "LIHC",
  "THCA", "UCEC", "GBM", "STAD", "PRAD"
)

gene_order_all <- c(
  "DROSHA", "DGCR8", "DICER1", "TARBP2",
  "XPO1", "XPO5",
  "AGO1", "AGO2", "AGO3", "AGO4",
  "TNRC6A", "PRKRA",
  "GEMIN4", "DDX5", "DDX17", "DDX20"
)

# The example integration figure uses these six genes.
# To include all 16 genes, replace this vector with gene_order_all.
selected_genes_for_integration <- c(
  "AGO2", "DICER1", "DROSHA", "TARBP2", "XPO1", "XPO5"
)

alteration_type_order <- c(
  "Amplification",
  "Deep deletion",
  "Multiple alterations",
  "Point mutation"
)


# Colors
#_________________________________________________________________

alteration_colors <- c(
  "Amplification" = "#40B3A2",
  "Deep deletion" = "#F2C94C",
  "Multiple alterations" = "#2D9CDB",
  "Point mutation" = "#EB5757"
)


# Helper functions
#_________________________________________________________________

format_number_label <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(
      abs(x - round(x)) < 0.05,
      as.character(round(x)),
      as.character(round(x, 2))
    )
  )
}


standardize_cancer_code <- function(x) {
  
  x_clean <- toupper(as.character(x))
  
  case_when(
    str_detect(x_clean, "COAD|COLON") ~ "COAD",
    str_detect(x_clean, "ESCA|ESOPH") ~ "ESCA",
    str_detect(x_clean, "BLCA|BLADDER") ~ "BLCA",
    str_detect(x_clean, "BRCA|BREAST") ~ "BRCA",
    str_detect(x_clean, "LUAD|LUNG ADENOCARCINOMA") ~ "LUAD",
    str_detect(x_clean, "LUSC|LUNG SQUAMOUS") ~ "LUSC",
    str_detect(x_clean, "KIRC|RENAL CLEAR|KIDNEY") ~ "KIRC",
    str_detect(x_clean, "LIHC|HEPATOCELLULAR|LIVER") ~ "LIHC",
    str_detect(x_clean, "THCA|THYROID") ~ "THCA",
    str_detect(x_clean, "UCEC|UTERINE ENDOMETRIOID|ENDOMETRIAL") ~ "UCEC",
    str_detect(x_clean, "GBM|GLIOBLASTOMA") ~ "GBM",
    str_detect(x_clean, "STAD|STOMACH|GASTRIC") ~ "STAD",
    str_detect(x_clean, "PRAD|PROSTATE") ~ "PRAD",
    TRUE ~ NA_character_
  )
}


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


# ============================================================
# Read and standardize cBioPortal alteration data
# ============================================================

cbio_data_raw <- readr::read_csv(
  cbioportal_file,
  show_col_types = FALSE
)

names(cbio_data_raw) <- names(cbio_data_raw) %>%
  str_replace_all("\\s+", "_")

cancer_col_cbio <- intersect(
  c("Cancer", "Cancer_Type", "Cancer_Study"),
  names(cbio_data_raw)
)

if (length(cancer_col_cbio) == 0) {
  stop("No cancer column found in cBioPortal table.")
}

required_cbio_cols <- c(
  "Gene",
  "Alteration_Frequency",
  "Alteration_Type",
  "Alteration_Count"
)

missing_cbio_cols <- setdiff(required_cbio_cols, names(cbio_data_raw))

if (length(missing_cbio_cols) > 0) {
  stop(
    "Missing required columns in cBioPortal table: ",
    paste(missing_cbio_cols, collapse = ", ")
  )
}

cbio_data <- cbio_data_raw %>%
  transmute(
    Cancer = standardize_cancer_code(.data[[cancer_col_cbio[1]]]),
    Gene = as.character(Gene),
    Alteration_Type = normalize_alteration_type(Alteration_Type),
    Alteration_Frequency = readr::parse_number(as.character(Alteration_Frequency)),
    Alteration_Count = readr::parse_number(as.character(Alteration_Count))
  ) %>%
  filter(
    !is.na(Cancer),
    Cancer %in% cancer_order,
    Gene %in% selected_genes_for_integration,
    Alteration_Type %in% alteration_type_order
  )


# Add missing zero alteration rows
#_________________________________________________________________

all_cbio_combinations <- expand_grid(
  Cancer = cancer_order,
  Gene = selected_genes_for_integration,
  Alteration_Type = alteration_type_order
)

cbio_data_complete <- all_cbio_combinations %>%
  left_join(
    cbio_data,
    by = c("Cancer", "Gene", "Alteration_Type")
  ) %>%
  mutate(
    Alteration_Frequency = replace_na(Alteration_Frequency, 0),
    Alteration_Count = replace_na(Alteration_Count, 0),
    Cancer = factor(Cancer, levels = cancer_order),
    Gene = factor(Gene, levels = selected_genes_for_integration),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )


# ============================================================
# Read and standardize DESeq2 expression data
# ============================================================

expression_raw <- readr::read_csv(
  expression_file,
  show_col_types = FALSE
)

names(expression_raw) <- names(expression_raw) %>%
  str_replace_all("\\s+", "_")

required_expression_cols <- c(
  "Cancer",
  "Gene",
  "log2FoldChange"
)

missing_expression_cols <- setdiff(
  required_expression_cols,
  names(expression_raw)
)

if (length(missing_expression_cols) > 0) {
  stop(
    "The expression file is missing required columns: ",
    paste(missing_expression_cols, collapse = ", "),
    "\nExpected file: results/combined_results/combined_DESeq2_results_selected_16_genes_all_cancers.csv"
  )
}

padj_col_expr <- if ("padj" %in% names(expression_raw)) {
  "padj"
} else {
  NA_character_
}

expression_data <- expression_raw %>%
  transmute(
    Cancer = standardize_cancer_code(Cancer),
    Gene = as.character(Gene),
    log2FC = as.numeric(log2FoldChange),
    padj = if (!is.na(padj_col_expr)) {
      as.numeric(.data[[padj_col_expr]])
    } else {
      NA_real_
    }
  ) %>%
  filter(
    !is.na(Cancer),
    Cancer %in% cancer_order,
    Gene %in% selected_genes_for_integration
  ) %>%
  group_by(Cancer, Gene) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    padj = suppressWarnings(min(padj, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    padj = if_else(is.infinite(padj), NA_real_, padj),
    Cancer = factor(Cancer, levels = cancer_order),
    Gene = factor(Gene, levels = selected_genes_for_integration)
  )


# ============================================================
# Create integration tables
# ============================================================

integration_long <- cbio_data_complete %>%
  left_join(
    expression_data,
    by = c("Cancer", "Gene")
  ) %>%
  mutate(
    Dataset = "TCGA",
    Cancer = factor(Cancer, levels = cancer_order),
    Gene = factor(Gene, levels = selected_genes_for_integration),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  ) %>%
  select(
    Dataset,
    Cancer,
    Gene,
    Alteration_Type,
    Alteration_Frequency,
    Alteration_Count,
    log2FC,
    padj
  ) %>%
  arrange(Cancer, Gene, Alteration_Type)

integration_summary <- integration_long %>%
  group_by(
    Dataset,
    Cancer,
    Gene
  ) %>%
  summarise(
    Total_Alteration_Frequency = sum(Alteration_Frequency, na.rm = TRUE),
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    log2FC = first(log2FC),
    padj = first(padj),
    .groups = "drop"
  ) %>%
  arrange(Cancer, Gene)

readr::write_csv(
  integration_long,
  file.path(processed_dir, "tcga_expression_alteration_integration_long.csv")
)

readr::write_csv(
  integration_summary,
  file.path(processed_dir, "tcga_expression_alteration_integration_summary.csv")
)


# Check missing expression values
#_________________________________________________________________

missing_expression_check <- integration_summary %>%
  filter(is.na(log2FC)) %>%
  arrange(Cancer, Gene)

readr::write_csv(
  missing_expression_check,
  file.path(tables_dir, "tcga_expression_alteration_missing_log2fc_check.csv")
)

if (nrow(missing_expression_check) > 0) {
  warning(
    "Some cancer-gene combinations have missing log2FC values. ",
    "Check results/tables/tcga_expression_alteration_missing_log2fc_check.csv"
  )
}


# ============================================================
# Function to create one cancer-specific panel
# ============================================================

make_integration_panel <- function(cancer_code) {
  
  bar_data <- integration_long %>%
    filter(Cancer == cancer_code) %>%
    mutate(
      Gene = factor(Gene, levels = selected_genes_for_integration),
      Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
    )
  
  heatmap_data <- integration_summary %>%
    filter(Cancer == cancer_code) %>%
    mutate(
      Gene = factor(Gene, levels = selected_genes_for_integration),
      log2FC_Label = format_number_label(log2FC)
    )
  
  gene_total_frequency <- bar_data %>%
    group_by(Gene) %>%
    summarise(
      Total_Alteration_Frequency = sum(Alteration_Frequency, na.rm = TRUE),
      .groups = "drop"
    )
  
  max_alt_freq <- max(
    gene_total_frequency$Total_Alteration_Frequency,
    na.rm = TRUE
  )
  
  if (is.na(max_alt_freq) || max_alt_freq == 0) {
    max_alt_freq <- 1
  }
  
  y_limit <- max_alt_freq * 1.25
  
  max_abs_log2fc <- max(
    abs(integration_summary$log2FC),
    na.rm = TRUE
  )
  
  if (is.na(max_abs_log2fc) || max_abs_log2fc == 0) {
    max_abs_log2fc <- 1
  }
  
  p_bar <- ggplot(
    bar_data,
    aes(
      x = Gene,
      y = Alteration_Frequency,
      fill = Alteration_Type
    )
  ) +
    geom_col(
      width = 0.8,
      color = "black",
      linewidth = 0.15
    ) +
    scale_fill_manual(
      values = alteration_colors,
      drop = FALSE,
      name = "Alteration type"
    ) +
    scale_y_continuous(
      limits = c(0, y_limit),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = NULL,
      y = "Alteration frequency (%)"
    ) +
    annotate(
      "text",
      x = length(selected_genes_for_integration),
      y = y_limit * 0.92,
      label = cancer_code,
      fontface = "bold",
      size = 4,
      hjust = 1
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 0, 5)
    )
  
  p_heatmap <- ggplot(
    heatmap_data,
    aes(
      x = Gene,
      y = "log2FC",
      fill = log2FC
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.4,
      height = 0.9
    ) +
    geom_text(
      aes(label = log2FC_Label),
      size = 3.5,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-max_abs_log2fc, max_abs_log2fc),
      name = "log2FC"
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 5, 5, 5)
    )
  
  p_combined <- p_bar / p_heatmap +
    patchwork::plot_layout(
      heights = c(4, 0.9)
    )
  
  p_combined
}


# ============================================================
# Generate only one final combined figure
# ============================================================

# ============================================================
# Generate only one final combined figure
# with compact manual color codes for alteration type and log2FC
# ============================================================

integration_panels <- map(
  cancer_order,
  make_integration_panel
)

names(integration_panels) <- cancer_order

panel_grid <- patchwork::wrap_plots(
  integration_panels,
  ncol = 3
)


# -----------------------------
# Compact square color code for alteration types
# -----------------------------

alteration_color_code_data <- tibble(
  Alteration_Type = factor(
    alteration_type_order,
    levels = alteration_type_order
  ),
  y = rev(seq(
    from = 1,
    by = 0.42,
    length.out = length(alteration_type_order)
  )),
  x = 1
)

p_alteration_color_code <- ggplot(
  alteration_color_code_data,
  aes(
    x = x,
    y = y,
    fill = Alteration_Type
  )
) +
  geom_point(
    shape = 22,
    size = 5.5,
    color = "black",
    stroke = 0.25
  ) +
  geom_text(
    aes(
      x = 1.18,
      label = Alteration_Type
    ),
    hjust = 0,
    size = 3.2,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = alteration_colors,
    drop = FALSE
  ) +
  annotate(
    "text",
    x = 0.9,
    y = max(alteration_color_code_data$y) + 0.35,
    label = "Alteration type",
    hjust = 0,
    fontface = "bold",
    size = 3.8
  ) +
  coord_cartesian(
    xlim = c(0.85, 3.6),
    ylim = c(
      min(alteration_color_code_data$y) - 0.25,
      max(alteration_color_code_data$y) + 0.55
    ),
    clip = "off"
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )


# -----------------------------
# Short compact color code for log2FC
# -----------------------------

global_max_abs_log2fc <- max(
  abs(integration_summary$log2FC),
  na.rm = TRUE
)

if (is.na(global_max_abs_log2fc) || global_max_abs_log2fc == 0) {
  global_max_abs_log2fc <- 1
}

log2fc_color_code_data <- tibble(
  x = 1,
  y = seq(
    -global_max_abs_log2fc,
    global_max_abs_log2fc,
    length.out = 120
  ),
  log2FC = seq(
    -global_max_abs_log2fc,
    global_max_abs_log2fc,
    length.out = 120
  )
)

p_log2fc_color_code <- ggplot(
  log2fc_color_code_data,
  aes(
    x = x,
    y = y,
    fill = log2FC
  )
) +
  geom_tile(
    width = 0.22,
    height = (2 * global_max_abs_log2fc) / 120
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-global_max_abs_log2fc, global_max_abs_log2fc),
    guide = "none"
  ) +
  annotate(
    "text",
    x = 0.75,
    y = global_max_abs_log2fc,
    label = format_number_label(global_max_abs_log2fc),
    hjust = 1,
    size = 3
  ) +
  annotate(
    "text",
    x = 0.75,
    y = 0,
    label = "0",
    hjust = 1,
    size = 3
  ) +
  annotate(
    "text",
    x = 0.75,
    y = -global_max_abs_log2fc,
    label = format_number_label(-global_max_abs_log2fc),
    hjust = 1,
    size = 3
  ) +
  annotate(
    "text",
    x = 0.72,
    y = global_max_abs_log2fc * 1.25,
    label = "log2FC",
    hjust = 0,
    fontface = "bold",
    size = 3.8
  ) +
  coord_cartesian(
    xlim = c(0.45, 1.6),
    ylim = c(
      -global_max_abs_log2fc,
      global_max_abs_log2fc * 1.35
    ),
    clip = "off"
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )


# -----------------------------
# Combine manual color codes
# The spacer keeps the log2FC bar short instead of stretching vertically
# -----------------------------

manual_color_codes <- p_alteration_color_code /
  p_log2fc_color_code /
  patchwork::plot_spacer() +
  patchwork::plot_layout(
    heights = c(0.55, 0.75, 3.2)
  )


# -----------------------------
# Combine main panel grid and compact color codes
# -----------------------------

combined_all_cancers <- panel_grid | manual_color_codes

combined_all_cancers <- combined_all_cancers +
  patchwork::plot_layout(
    widths = c(1, 0.17)
  )


# -----------------------------
# Save final combined figure
# -----------------------------

ggsave(
  filename = file.path(figures_dir, "tcga_expression_alteration_all_cancers.pdf"),
  plot = combined_all_cancers,
  width = 15.5,
  height = 16,
  bg = "white"
)

ggsave(
  filename = file.path(figures_dir, "tcga_expression_alteration_all_cancers.png"),
  plot = combined_all_cancers,
  width = 15.5,
  height = 16,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(figures_dir, "tcga_expression_alteration_all_cancers.tiff"),
  plot = combined_all_cancers,
  width = 15.5,
  height = 16,
  dpi = 300,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)
