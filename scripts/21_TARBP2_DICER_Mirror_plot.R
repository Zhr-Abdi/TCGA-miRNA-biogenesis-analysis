# ============================================================
# TARBP2 - DICER1 analysis
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# Cancer types
# ------------------------------------------------------------

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

results_dir <- file.path(
  getwd(),
  "01_TCGA_expression",
  "results_by_cancer"
)

output_dir <- file.path(
  getwd(),
  "02_TARBP2_DICER1_analysis"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# PART 1
# Extract TARBP2 and DICER1 DESeq2 results
# ============================================================

gene_results <- list()

for (cancer_code in names(cancer_projects)) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  file_path <- file.path(
    results_dir,
    project_id,
    paste0(
      "DESeq2_results_selected_16_genes_",
      cancer_code,
      ".csv"
    )
  )
  
  df <- read.csv(
    file_path,
    check.names = FALSE
  )
  
  df_2genes <- df %>%
    dplyr::filter(Gene %in% c("TARBP2", "DICER1")) %>%
    dplyr::select(
      Cancer,
      Gene,
      log2FoldChange,
      padj
    )
  
  gene_results[[cancer_code]] <- df_2genes
}

gene_results <- bind_rows(gene_results)

print(gene_results)

# ============================================================
# TARBP2 vs DICER1
# Horizontal grouped bar plot
# ============================================================


library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Cancer order
# ------------------------------------------------------------

cancer_order <- c(
  "BLCA",
  "BRCA",
  "COAD",
  "ESCA",
  "GBM",
  "KIRC",
  "LIHC",
  "LUAD",
  "LUSC",
  "PRAD",
  "STAD",
  "THCA",
  "UCEC"
)

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

plot_data <- gene_results %>%
  mutate(
    Cancer = factor(
      Cancer,
      levels = cancer_order
    ),
    Gene = factor(
      Gene,
      levels = c("DICER1", "TARBP2")
    )
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = log2FoldChange,
    y = Cancer,
    fill = Gene
  )
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey35"
  ) +
  scale_fill_manual(
    values = c(
      "DICER1" = "#3B67B1",
      "TARBP2" = "#D94A24"
    )
  ) +
  scale_x_continuous(
    limits = c(-0.75, 1.25),
    breaks = seq(-0.75, 1.25, by = 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = expression(log[2] ~ "fold change (Tumor vs Normal)"),
    y = "Cancer type",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.5,
      linetype = "dotted"
    ),
    panel.grid.major.y = element_line(
      color = "grey92",
      linewidth = 0.5
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
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "right",
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 20, 10, 10)
  )

print(p)


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

project_dir <- getwd()
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

# TIFF
ggsave(
  filename = file.path(figures_dir, "TARBP2_DICER1_barplot.tiff"),
  plot = p,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# PNG
ggsave(
  filename = file.path(figures_dir, "TARBP2_DICER1_barplot.png"),
  plot = p,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)



