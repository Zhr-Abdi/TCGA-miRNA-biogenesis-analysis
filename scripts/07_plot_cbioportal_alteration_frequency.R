# ============================================================
# 07_plot_cbioportal_alteration_frequency.R
#
# Project: TCGA-miRNA-biogenesis-analysis
#
# Purpose:
#   Generate cBioPortal-based alteration frequency plots from the
#   processed long alteration table.
#
# Input:
#   data/processed/cbioportal_alteration_frequency_long.csv
#
# Output:
#   results/figures/cbioportal_faceted_alteration_frequency.*
#   results/figures/cbioportal_gene_alteration_counts_stacked_bar.*
#   results/figures/cbioportal_alteration_type_pie_chart.*
#   results/figures/cbioportal_gene_alteration_bar_and_pie.*
#   results/figures/cbioportal_gene_frequency_correlation_heatmap.*
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


# Define order
#_________________________________________________________________

gene_order <- c(
  "AGO1", "AGO2", "AGO3", "AGO4",
  "DDX17", "DDX20", "DDX5",
  "DGCR8", "DICER1", "DROSHA",
  "GEMIN4", "PRKRA", "TARBP2",
  "TNRC6A", "XPO1", "XPO5"
)

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


cbio_data <- cbio_data %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Cancer = factor(Cancer, levels = cancer_order),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )



# Colors
#_________________________________________________________________

alteration_colors <- c(
  "Amplification" = "#40B3A2",
  "Deep deletion" = "#F2C94C",
  "Multiple alterations" = "#2D9CDB",
  "Point mutation" = "#EB5757"
)

selected_gene_colors <- c(
  "AGO2" = "#40B3A2",
  "DICER1" = "#56C7D9",
  "DROSHA" = "#2D75B6",
  "TARBP2" = "#F05A3D",
  "XPO1" = "#F2D84B",
  "XPO5" = "#E8A0D8"
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


# Helper function for percent labels
#_________________________________________________________________

format_percent_label <- function(x) {
  ifelse(
    abs(x - round(x)) < 0.05,
    paste0(round(x), "%"),
    paste0(round(x, 1), "%")
  )
}


# ============================================================
# Plot 1:
# Faceted alteration frequency plot
# Similar to panel A in your example
# ============================================================

selected_genes_for_faceted_plot <- c(
  "AGO2", "DICER1", "DROSHA", "TARBP2", "XPO1", "XPO5"
)

faceted_data <- cbio_data %>%
  filter(Gene %in% selected_genes_for_faceted_plot) %>%
  mutate(
    Gene = factor(Gene, levels = rev(selected_genes_for_faceted_plot)),
    Cancer = factor(Cancer, levels = cancer_order),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )

p_faceted <- ggplot(
  faceted_data,
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
    values = selected_gene_colors,
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
    legend.text = element_text(face = "bold.italic")
  )

save_ggplot_all_formats(
  plot_object = p_faceted,
  file_stem = "cbioportal_faceted_alteration_frequency",
  width = 8,
  height = 10
)


# ============================================================
# Plot 2:
# Stacked bar plot of alteration counts by gene
# Similar to your TCGA stacked bar plot
# ============================================================

total_tcga_samples <- cbio_data %>%
  distinct(Cancer, Total_Samples) %>%
  summarise(Total_TCGA_Samples = sum(Total_Samples, na.rm = TRUE)) %>%
  pull(Total_TCGA_Samples)

gene_count_data <- cbio_data %>%
  group_by(Gene, Alteration_Type) %>%
  summarise(
    Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  )

gene_label_data <- cbio_data %>%
  group_by(Gene) %>%
  summarise(
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    Overall_Percent = 100 * Total_Alteration_Count / total_tcga_samples,
    .groups = "drop"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Percent_Label = format_percent_label(Overall_Percent)
  )

p_stacked_bar <- ggplot(
  gene_count_data,
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
    data = gene_label_data,
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
    title = "TCGA",
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
  plot_object = p_stacked_bar,
  file_stem = "cbioportal_gene_alteration_counts_stacked_bar",
  width = 8,
  height = 5
)

# ============================================================
# Plot 3:
# Pie chart of alteration type proportions
# Corrected label positions
# ============================================================

pie_data <- cbio_data %>%
  group_by(Alteration_Type) %>%
  summarise(
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_order)
  ) %>%
  arrange(Alteration_Type) %>%
  filter(Total_Alteration_Count > 0) %>%
  mutate(
    Percent = 100 * Total_Alteration_Count / sum(Total_Alteration_Count),
    Label = format_percent_label(Percent),
    fraction = Total_Alteration_Count / sum(Total_Alteration_Count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    ymid = (ymin + ymax) / 2
  )

p_pie <- ggplot(pie_data) +
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
  plot_object = p_pie,
  file_stem = "cbioportal_alteration_type_pie_chart",
  width = 5,
  height = 5
)
# ============================================================
# Plot 4:
# Gene alteration frequency distance heatmap
# Distance-based heatmap instead of correlation heatmap
# ============================================================
# Important:
# This heatmap shows the distance between genes based on their
# total alteration frequency profiles across cancer types.
#
# Lower distance = more similar alteration frequency pattern
# Higher distance = more different alteration frequency pattern
#
# This is not patient-level co-occurrence analysis.
# Patient-level cBioPortal data would be needed for true
# co-occurrence or mutual exclusivity analysis.
# ============================================================

selected_genes_for_heatmap <- c(
  "DROSHA", "AGO2", "XPO5", "TARBP2", "XPO1", "DICER1"
)

distance_input <- cbio_data %>%
  filter(Gene %in% selected_genes_for_heatmap) %>%
  group_by(Cancer, Gene) %>%
  summarise(
    Total_Alteration_Frequency = sum(Alteration_Frequency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = selected_genes_for_heatmap),
    Cancer = factor(Cancer, levels = cancer_order)
  ) %>%
  arrange(Cancer) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Total_Alteration_Frequency,
    values_fill = 0
  )

# Convert to matrix:
# Rows = cancer types
# Columns = genes

frequency_matrix <- distance_input %>%
  select(all_of(selected_genes_for_heatmap)) %>%
  as.matrix()

rownames(frequency_matrix) <- distance_input$Cancer

# Transpose matrix:
# Rows = genes
# Columns = cancer types

gene_frequency_matrix <- t(frequency_matrix)

# Calculate Euclidean distance between genes

distance_matrix <- as.matrix(
  dist(
    gene_frequency_matrix,
    method = "euclidean"
  )
)

# Save distance matrix as table

distance_table <- as.data.frame(distance_matrix) %>%
  tibble::rownames_to_column("Gene")

readr::write_csv(
  distance_table,
  file.path(tables_dir, "cbioportal_gene_frequency_distance_matrix.csv"),
  na = ""
)

# Save distance heatmap in 3 formats
#_________________________________________________________________

save_distance_heatmap_all_formats <- function(matrix_object,
                                              file_stem,
                                              width = 5,
                                              height = 5,
                                              dpi = 300) {
  
  heatmap_colors <- colorRampPalette(
    c("#2166AC", "white", "#F46D43")
  )(100)
  
  max_value <- max(matrix_object, na.rm = TRUE)
  
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
    main = "Gene alteration frequency distance"
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
    main = "Gene alteration frequency distance"
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
    main = "Gene alteration frequency distance"
  )
  
  dev.off()
}

save_distance_heatmap_all_formats(
  matrix_object = distance_matrix,
  file_stem = "cbioportal_gene_frequency_distance_heatmap",
  width = 5,
  height = 5
)
