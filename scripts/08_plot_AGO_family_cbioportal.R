# ============================================================
# 08_plot_AGO_family_cbioportal.R
#
# Project: TCGA-miRNA-biogenesis-analysis
#
# Purpose:
#   Generate cBioPortal-based alteration frequency plots
#   specifically for the AGO gene family:
#   AGO1, AGO2, AGO3, AGO4
#
# Input:
#   data/processed/cbioportal_alteration_frequency_long.csv
#
# Outputs:
#   results/figures/cbioportal_AGO_family_faceted_alteration_frequency.*
#   results/figures/cbioportal_AGO_family_gene_frequency_distance_heatmap.*
#   results/figures/cbioportal_AGO_family_combined_panel.*
#   results/tables/cbioportal_AGO_family_gene_frequency_distance_matrix.csv
# ============================================================



# Load packages
#_________________________________________________________________

required_packages <- c("tidyverse", "pheatmap")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' is required. Please install it first."))
  }
}

library(tidyverse)
library(pheatmap)


# Define paths
#_________________________________________________________________

processed_file <- "data/processed/cbioportal_alteration_frequency_long.csv"

figures_dir <- "results/figures"
tables_dir <- "results/tables"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)


# Read processed data
#_________________________________________________________________

cbio_data <- readr::read_csv(processed_file, show_col_types = FALSE)


# Define AGO family genes and plot order
#_________________________________________________________________

ago_genes <- c("AGO1", "AGO2", "AGO3", "AGO4")

cancer_order <- c(
  "BLCA", "BRCA", "COAD", "ESCA", "GBM",
  "KIRC", "LIHC", "LUAD", "LUSC", "PRAD",
  "STAD", "THCA", "UCEC"
)

alteration_type_order <- c(
  "Amplification",
  "Deep deletion",
  "Multiple alterations",
  "Point mutation"
)


# Prepare data
#_________________________________________________________________


cbio_data <- cbio_data %>%
  mutate(
    Gene = as.character(Gene),
    Cancer = factor(Cancer, levels = cancer_order),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )



# Colors for AGO family
#_________________________________________________________________

ago_gene_colors <- c(
  "AGO1" = "#40B3A2",
  "AGO2" = "#56C7D9",
  "AGO3" = "#2D75B6",
  "AGO4" = "#F05A3D"
)


# Helper function to save ggplot in 3 formats
#_________________________________________________________________

save_ggplot_all_formats <- function(plot_object, file_stem, width, height, dpi = 300) {
  
  ggsave(
    filename = file.path(figures_dir, paste0(file_stem, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(figures_dir, paste0(file_stem, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(figures_dir, paste0(file_stem, ".tiff")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )
}


# ============================================================
# Plot A:
# Faceted alteration frequency plot for AGO family
# ============================================================

ago_faceted_data <- cbio_data %>%
  filter(Gene %in% ago_genes) %>%
  mutate(
    Gene = factor(Gene, levels = rev(ago_genes)),
    Cancer = factor(Cancer, levels = cancer_order),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )

p_ago_faceted <- ggplot(
  ago_faceted_data,
  aes(
    x = Alteration_Frequency,
    y = Gene,
    fill = Gene
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.15
  ) +
  facet_grid(
    rows = vars(Cancer),
    cols = vars(Alteration_Type),
    scales = "free_x",
    switch = "y"
  ) +
  scale_fill_manual(
    values = ago_gene_colors,
    drop = FALSE
  ) +
  labs(
    x = "Alteration frequency",
    y = NULL,
    fill = "Gene"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text = element_text(face = "bold", size = 8),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 7),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "grey92"),
    panel.spacing = unit(0.08, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold.italic"),
    plot.margin = margin(5, 5, 5, 5)
  )

save_ggplot_all_formats(
  plot_object = p_ago_faceted,
  file_stem = "cbioportal_AGO_family_faceted_alteration_frequency",
  width = 8,
  height = 10
)


# ============================================================
# Plot B:
# Distance-based heatmap for AGO family
# ============================================================
# This heatmap compares AGO1, AGO2, AGO3, and AGO4 based on
# their pan-cancer alteration frequency profiles.
#
# Lower distance = more similar alteration pattern
# Higher distance = more different alteration pattern
# ============================================================

ago_distance_input <- cbio_data %>%
  filter(Gene %in% ago_genes) %>%
  group_by(Cancer, Gene) %>%
  summarise(
    Total_Alteration_Frequency = sum(Alteration_Frequency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = ago_genes),
    Cancer = factor(Cancer, levels = cancer_order)
  ) %>%
  arrange(Cancer) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Total_Alteration_Frequency,
    values_fill = 0
  )

frequency_matrix <- ago_distance_input %>%
  select(all_of(ago_genes)) %>%
  as.matrix()

rownames(frequency_matrix) <- ago_distance_input$Cancer

# Rows = genes
# Columns = cancer types

ago_gene_frequency_matrix <- t(frequency_matrix)

# Euclidean distance between AGO genes

ago_distance_matrix <- as.matrix(
  dist(
    ago_gene_frequency_matrix,
    method = "euclidean"
  )
)

# Save distance matrix

ago_distance_table <- as.data.frame(ago_distance_matrix) %>%
  tibble::rownames_to_column("Gene")

readr::write_csv(
  ago_distance_table,
  file.path(tables_dir, "cbioportal_AGO_family_gene_frequency_distance_matrix.csv"),
  na = ""
)


# ----------------------------
# Save AGO heatmap in png, pdf, and tiff
# ----------------------------

save_ago_distance_heatmap_all_formats <- function(matrix_object,
                                                  file_stem,
                                                  width = 4.5,
                                                  height = 4.5,
                                                  dpi = 300) {
  
  heatmap_colors <- colorRampPalette(
    c("#2166AC", "white", "#F46D43")
  )(100)
  
  max_value <- max(matrix_object, na.rm = TRUE)
  
  if (max_value <= 0) {
    max_value <- 1
  }
  
  heatmap_breaks <- seq(
    0,
    max_value,
    length.out = 101
  )
  
  png(
    filename = file.path(figures_dir, paste0(file_stem, ".png")),
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  
  pheatmap::pheatmap(
    matrix_object,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = as.dist(matrix_object),
    clustering_distance_cols = as.dist(matrix_object),
    border_color = NA,
    fontsize = 10,
    main = "AGO family alteration frequency distance"
  )
  
  dev.off()
  
  pdf(
    file = file.path(figures_dir, paste0(file_stem, ".pdf")),
    width = width,
    height = height
  )
  
  pheatmap::pheatmap(
    matrix_object,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = as.dist(matrix_object),
    clustering_distance_cols = as.dist(matrix_object),
    border_color = NA,
    fontsize = 10,
    main = "AGO family alteration frequency distance"
  )
  
  dev.off()
  
  tiff(
    filename = file.path(figures_dir, paste0(file_stem, ".tiff")),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  
  pheatmap::pheatmap(
    matrix_object,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = as.dist(matrix_object),
    clustering_distance_cols = as.dist(matrix_object),
    border_color = NA,
    fontsize = 10,
    main = "AGO family alteration frequency distance"
  )
  
  dev.off()
}

save_ago_distance_heatmap_all_formats(
  matrix_object = ago_distance_matrix,
  file_stem = "cbioportal_AGO_family_gene_frequency_distance_heatmap",
  width = 4.5,
  height = 4.5
)

