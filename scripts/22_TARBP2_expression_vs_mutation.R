# ============================================================
# Expression dysregulation vs genomic alteration
# miRNA biogenesis genes across TCGA cancers
#
# X-axis:
#   Somatic alteration frequency (%)
#
# Y-axis:
#   Median absolute log2 fold change across 13 cancers
#
# Bubble size:
#   Number of cancer types with significant DE
#   Significance criterion: raw p-value < 0.05
#
# Bubble color:
#   Predominant direction of significant expression changes
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)


# ============================================================
# 2. Project directory
# ============================================================

# Run this script from the root directory of the project

project_dir <- getwd()


# ============================================================
# 3. Input directories
# ============================================================

expression_results_dir <- file.path(
  project_dir,
  "01_TCGA_expression",
  "results_by_cancer"
)


mutation_file <- file.path(
  project_dir,
  "03_TCGA_mutation",
  "combined_results",
  "cbioportal_gene_summary.csv"
)


# ============================================================
# 4. Output directories
# ============================================================

analysis_dir <- file.path(
  project_dir,
  "02_TARBP2_DICER1_analysis"
)

dir.create(
  analysis_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


figures_dir <- file.path(
  project_dir,
  "Figures",
  "07_TARBP2_DICER1"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 5. TCGA cancer projects
# ============================================================

cancer_projects <- c(
  COAD = "TCGA-COAD",
  ESCA = "TCGA-ESCA",
  BLCA = "TCGA-BLCA",
  BRCA = "TCGA-BRCA",
  LUAD = "TCGA-LUAD",
  LUSC = "TCGA-LUSC",
  KIRC = "TCGA-KIRC",
  LIHC = "TCGA-LIHC",
  THCA = "TCGA-THCA",
  UCEC = "TCGA-UCEC",
  GBM  = "TCGA-GBM",
  STAD = "TCGA-STAD",
  PRAD = "TCGA-PRAD"
)


# ============================================================
# 6. miRNA biogenesis genes
# ============================================================

genes_16 <- c(
  "AGO1",
  "AGO2",
  "AGO3",
  "AGO4",
  "DDX5",
  "DDX17",
  "DDX20",
  "DGCR8",
  "DICER1",
  "DROSHA",
  "GEMIN4",
  "PRKRA",
  "TARBP2",
  "TNRC6A",
  "XPO1",
  "XPO5"
)


# ============================================================
# 7. Significance threshold
#
# IMPORTANT:
# Significant DE is defined ONLY by raw p-value < 0.05
# No log2FC cutoff is applied.
# ============================================================

pvalue_cutoff <- 0.05


# ============================================================
# 8. Read DESeq2 results from all 13 cancer types
# ============================================================

expression_list <- list()


for (cancer_code in names(cancer_projects)) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  
  file_path <- file.path(
    expression_results_dir,
    project_id,
    paste0(
      "DESeq2_results_selected_16_genes_",
      cancer_code,
      ".csv"
    )
  )
  
  
  # ----------------------------------------------------------
  # Check file
  # ----------------------------------------------------------
  
  if (!file.exists(file_path)) {
    
    stop(
      paste(
        "Expression file not found:",
        file_path
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Read file
  # ----------------------------------------------------------
  
  df <- read.csv(
    file_path,
    check.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Check required columns
  # ----------------------------------------------------------
  
  required_columns <- c(
    "Gene",
    "log2FoldChange",
    "pvalue"
  )
  
  
  if (!all(required_columns %in% colnames(df))) {
    
    stop(
      paste(
        "Required columns missing in:",
        file_path
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Keep the 16 genes
  # ----------------------------------------------------------
  
  df_selected <- df %>%
    
    filter(
      Gene %in% genes_16
    ) %>%
    
    transmute(
      Cancer = cancer_code,
      Gene,
      log2FoldChange,
      pvalue
    )
  
  
  expression_list[[cancer_code]] <- df_selected
}


# ============================================================
# 9. Combine expression data
# ============================================================

expression_all <- bind_rows(
  expression_list
)

print(expression_all)


# ============================================================
# 10. Classify differential expression
#
# p-value < 0.05 + positive log2FC = Upregulated
# p-value < 0.05 + negative log2FC = Downregulated
# otherwise = Not significant
# ============================================================

expression_all <- expression_all %>%
  
  mutate(
    
    DE_status = case_when(
      
      !is.na(pvalue) &
        pvalue < pvalue_cutoff &
        log2FoldChange > 0 ~
        "Upregulated",
      
      !is.na(pvalue) &
        pvalue < pvalue_cutoff &
        log2FoldChange < 0 ~
        "Downregulated",
      
      TRUE ~
        "Not significant"
    )
  )


# ============================================================
# 11. Summarize expression behavior for each gene
# ============================================================

expression_summary <- expression_all %>%
  
  group_by(Gene) %>%
  
  summarise(
    
    # Total number of cancers analyzed
    N_cancers = n(),
    
    # Number of significant cancers
    N_significant_DE = sum(
      DE_status != "Not significant",
      na.rm = TRUE
    ),
    
    # Number significantly upregulated
    N_up = sum(
      DE_status == "Upregulated",
      na.rm = TRUE
    ),
    
    # Number significantly downregulated
    N_down = sum(
      DE_status == "Downregulated",
      na.rm = TRUE
    ),
    
    # Expression magnitude used on Y-axis
    Median_abs_log2FC = median(
      abs(log2FoldChange),
      na.rm = TRUE
    ),
    
    # Additional values for reference
    Median_log2FC = median(
      log2FoldChange,
      na.rm = TRUE
    ),
    
    Mean_abs_log2FC = mean(
      abs(log2FoldChange),
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    # --------------------------------------------------------
    # Predominant direction among significant cancers
    # --------------------------------------------------------
    
    Predominant_direction = case_when(
      
      N_up > N_down ~
        "Upregulated",
      
      N_down > N_up ~
        "Downregulated",
      
      TRUE ~
        "Mixed / no clear direction"
    )
  )


print(expression_summary)


# ============================================================
# 12. Read genomic alteration data
# ============================================================

if (!file.exists(mutation_file)) {
  
  stop(
    paste(
      "Mutation summary file not found:",
      mutation_file
    )
  )
}


mutation_summary <- read.csv(
  mutation_file,
  check.names = FALSE
)


# ------------------------------------------------------------
# Check required mutation columns
# ------------------------------------------------------------

mutation_required_columns <- c(
  "Gene",
  "Overall_Alteration_Frequency_Percent"
)


if (!all(
  mutation_required_columns %in%
  colnames(mutation_summary)
)) {
  
  stop(
    paste(
      "Required mutation columns are missing."
    )
  )
}


# ------------------------------------------------------------
# Keep needed columns
# ------------------------------------------------------------

mutation_summary <- mutation_summary %>%
  
  filter(
    Gene %in% genes_16
  ) %>%
  
  select(
    Gene,
    Overall_Alteration_Frequency_Percent
  ) %>%
  
  rename(
    AlterationFrequency =
      Overall_Alteration_Frequency_Percent
  )


print(mutation_summary)


# ============================================================
# 13. Merge expression and alteration data
# ============================================================

plot_data <- expression_summary %>%
  
  left_join(
    mutation_summary,
    by = "Gene"
  )


# ============================================================
# 14. Check merged data
# ============================================================

print(
  plot_data %>%
    
    select(
      Gene,
      AlterationFrequency,
      Median_abs_log2FC,
      N_significant_DE,
      N_up,
      N_down,
      Predominant_direction
    )
)


# ------------------------------------------------------------
# Check missing mutation values
# ------------------------------------------------------------

if (any(
  is.na(plot_data$AlterationFrequency)
)) {
  
  warning(
    paste(
      "Missing alteration frequency for:",
      paste(
        plot_data$Gene[
          is.na(
            plot_data$AlterationFrequency
          )
        ],
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 15. Specifically check TARBP2 and DICER1
# ============================================================

print(
  plot_data %>%
    
    filter(
      Gene %in% c(
        "TARBP2",
        "DICER1"
      )
    ) %>%
    
    select(
      Gene,
      AlterationFrequency,
      Median_abs_log2FC,
      N_significant_DE,
      N_up,
      N_down,
      Predominant_direction
    )
)


# ============================================================
# 16. Save analysis table
# ============================================================

write.csv(
  plot_data,
  file.path(
    analysis_dir,
    "miRNA_genes_expression_vs_TCGA_alteration_pvalue005.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. Bubble plot
#
# X = alteration frequency
# Y = median |log2FC|
# Size = number of cancers with p < 0.05
# Color = predominant direction
# ============================================================
dicer1_color <- "#3B67B1"
tarbp2_color <- "#D94A24"

p_bubble <- ggplot(
  plot_data,
  aes(
    x = AlterationFrequency,
    y = Median_abs_log2FC
  )
) +
  

# ------------------------------------------------------------
# Main bubbles - all grey
# ------------------------------------------------------------

# Main bubbles - all genes in grey
geom_point(
  aes(
    size = N_significant_DE
  ),
  shape = 21,
  fill = "#BDBDBD",
  color = "white",
  stroke = 0.8,
  alpha = 0.9
) +
# ----------------------------------------------------------
# Highlight TARBP2
# ----------------------------------------------------------

# Highlight TARBP2
geom_point(
  data = plot_data %>%
    filter(Gene == "TARBP2"),
  aes(
    x = AlterationFrequency,
    y = Median_abs_log2FC,
    size = N_significant_DE
  ),
  shape = 21,
  fill = "#D94A24",
  color = "white",
  stroke = 1.0,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# Highlight DICER1
# ----------------------------------------------------------
# Highlight DICER1
geom_point(
  data = plot_data %>%
    filter(Gene == "DICER1"),
  aes(
    x = AlterationFrequency,
    y = Median_abs_log2FC,
    size = N_significant_DE
  ),
  shape = 21,
  fill = "#3B67B1",
  color = "white",
  stroke = 1.0,
  show.legend = FALSE
) +
  
  
  # ----------------------------------------------------------
# Labels for other genes
# ----------------------------------------------------------

geom_text_repel(
  data = plot_data %>%
    
    filter(
      !Gene %in% c(
        "TARBP2",
        "DICER1"
      )
    ),
  
  aes(
    label = Gene
  ),
  
  size = 3.6,
  
  color = "grey25",
  
  box.padding = 0.35,
  
  point.padding = 0.25,
  
  max.overlaps = Inf,
  
  seed = 123
) +
  
  # ----------------------------------------------------------
# TARBP2 label
# ----------------------------------------------------------

geom_text_repel(
  data = plot_data %>%
    filter(Gene == "TARBP2"),
  aes(label = Gene),
  size = 4.5,
  fontface = "bold.italic",
  color = "#8B1A1A",
  box.padding = 0.5,
  point.padding = 0.4,
  seed = 123,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# DICER1 label
# ----------------------------------------------------------

geom_text_repel(
  data = plot_data %>%
    filter(Gene == "DICER1"),
  aes(label = Gene),
  size = 4.5,
  fontface = "bold.italic",
  color = "#253494",
  box.padding = 0.5,
  point.padding = 0.4,
  seed = 123,
  show.legend = FALSE
) +  

  
# ==========================================================
# Bubble sizes
# ==========================================================

scale_size_continuous(
  
  range = c(
    3,
    12
  ),
  
  breaks = c(
    0,
    2,
    4,
    6,
    8,
    10,
    12
  ),
  
  limits = c(
    0,
    13
  ),
  
  name =
    "Cancer types with\nsignificant DE (p < 0.05)"
) +
  
  # ==========================================================
# Axes
# ==========================================================

scale_x_continuous(
  
  expand = expansion(
    mult = c(
      0.03,
      0.08
    )
  )
) +
  
  scale_y_continuous(
    
    expand = expansion(
      mult = c(
        0.05,
        0.10
      )
    )
  ) +
  
  labs(
    
    x =
      "Somatic alteration frequency (%)",
    
    y = expression(
      "Median " *
        "|" *
        log[2] *
        "FC" *
        "|" *
        " across 13 TCGA cancer types"
    )
  ) +
  
  # ==========================================================
# Theme
# ==========================================================

theme_classic(
  base_size = 13
) +
  
  theme(
    
    panel.grid.major =
      element_line(
        color = "grey92",
        linewidth = 0.45
      ),
    
    panel.grid.minor =
      element_blank(),
    
    # Full border
    panel.border =
      element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
    
    axis.line =
      element_blank(),
    
    axis.title.x =
      element_text(
        size = 14,
        face = "bold"
      ),
    
    axis.title.y =
      element_text(
        size = 14,
        face = "bold"
      ),
    
    axis.text =
      element_text(
        size = 11,
        color = "black"
      ),
    
    legend.title =
      element_text(
        size = 11,
        face = "bold"
      ),
    
    legend.text =
      element_text(
        size = 10
      ),
    
    legend.position =
      "right",
    
    plot.margin =
      margin(
        10,
        15,
        10,
        10
      )
  )


# ============================================================
# 18. Display plot
# ============================================================

print(p_bubble)


# ============================================================
# 19. Save plot
# ============================================================

figure_name <-
  "TARBP2_Scatter"


# ------------------------------------------------------------
# TIFF
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figures_dir,
    paste0(
      figure_name,
      ".tiff"
    )
  ),
  
  plot = p_bubble,
  
  width = 10,
  
  height = 7.5,
  
  units = "in",
  
  dpi = 600,
  
  compression = "lzw",
  
  bg = "white"
)


# ------------------------------------------------------------
# PNG
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figures_dir,
    paste0(
      figure_name,
      ".png"
    )
  ),
  
  plot = p_bubble,
  
  width = 10,
  
  height = 7.5,
  
  units = "in",
  
  dpi = 600,
  
  bg = "white"
)


# ============================================================
# 20. Final messages
# ============================================================

message(
  "Analysis table saved in: ",
  normalizePath(
    analysis_dir
  )
)

message(
  "Figures saved in: ",
  normalizePath(
    figures_dir
  )
)
# ============================================================
# Bubble plot
# Only TARBP2 highlighted in red
# All other genes in grey
# ============================================================

p_bubble_tarbp2_only <- ggplot(
  plot_data,
  aes(
    x = AlterationFrequency,
    y = Median_abs_log2FC
  )
) +
  
  # ----------------------------------------------------------
# All genes in grey
# ----------------------------------------------------------
geom_point(
  aes(size = N_significant_DE),
  shape = 21,
  fill = "#BDBDBD",
  color = "white",
  stroke = 0.8,
  alpha = 0.9
) +
  
  # ----------------------------------------------------------
# Highlight TARBP2 in red
# ----------------------------------------------------------
geom_point(
  data = plot_data %>%
    filter(Gene == "TARBP2"),
  aes(
    x = AlterationFrequency,
    y = Median_abs_log2FC,
    size = N_significant_DE
  ),
  shape = 21,
  fill = "#D94A24",
  color = "white",
  stroke = 1.0,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# Labels for all genes except TARBP2
# ----------------------------------------------------------
geom_text_repel(
  data = plot_data %>%
    filter(Gene != "TARBP2"),
  aes(label = Gene),
  size = 3.6,
  color = "grey25",
  box.padding = 0.35,
  point.padding = 0.25,
  max.overlaps = Inf,
  seed = 123
) +
  
  # ----------------------------------------------------------
# TARBP2 label
# ----------------------------------------------------------
geom_text_repel(
  data = plot_data %>%
    filter(Gene == "TARBP2"),
  aes(label = Gene),
  size = 4.7,
  fontface = "bold.italic",
  color = "#8B1A1A",
  box.padding = 0.5,
  point.padding = 0.4,
  seed = 123,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# Bubble sizes
# ----------------------------------------------------------
scale_size_continuous(
  range = c(3, 12),
  breaks = c(0, 2, 4, 6, 8, 10, 12),
  limits = c(0, 13),
  name = "Cancer types with\nsignificant DE (p < 0.05)"
) +
  
  # ----------------------------------------------------------
# Axes
# ----------------------------------------------------------
scale_x_continuous(
  expand = expansion(mult = c(0.03, 0.08))
) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  
  labs(
    x = "Somatic alteration frequency (%)",
    y = expression(
      "Median " * "|" * log[2] * "FC" * "|" *
        " across 13 TCGA cancer types"
    )
  ) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------
theme_classic(base_size = 13) +
  
  theme(
    panel.grid.major = element_line(
      color = "grey92",
      linewidth = 0.45
    ),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.line = element_blank(),
    
    axis.title.x = element_text(
      size = 14,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 14,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

print(p_bubble_tarbp2_only)
figure_name <- "miRNA_genes_expression_vs_TCGA_alteration_TARBP2_only"

ggsave(
  filename = file.path(
    figures_dir,
    paste0(figure_name, ".tiff")
  ),
  plot = p_bubble_tarbp2_only,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(
    figures_dir,
    paste0(figure_name, ".png")
  ),
  plot = p_bubble_tarbp2_only,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)
