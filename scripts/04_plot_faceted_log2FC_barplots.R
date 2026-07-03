# Required libraries
#_________________________________________________________________

library(dplyr)
library(ggplot2)


# Define TCGA cancer types
#_________________________________________________________________

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


# Define genes of interest
#_________________________________________________________________

genes_of_interest <- c(
  "AGO1", "AGO2", "AGO3", "AGO4",
  "DDX5", "DDX17", "DDX20",
  "DGCR8", "DROSHA", "DICER1",
  "GEMIN4", "PRKRA", "TNRC6A",
  "TARBP2", "XPO1", "XPO5"
)


# Define directories
#_________________________________________________________________

results_dir <- file.path(getwd(), "results")
combined_dir <- file.path(results_dir, "combined_results")
figures_dir <- file.path(results_dir, "figures")

dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# Read selected 16-gene DESeq2 results from all cancer types
#_________________________________________________________________

all_selected_results <- list()

for (cancer_code in names(cancer_projects)) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  input_file <- file.path(
    results_dir,
    project_id,
    paste0("DESeq2_results_selected_16_genes_", cancer_code, ".csv")
  )
  
  if (!file.exists(input_file)) {
    stop(paste("File not found:", input_file))
  }
  
  cancer_result <- read.csv(input_file, stringsAsFactors = FALSE)
  
  all_selected_results[[cancer_code]] <- cancer_result
}

all_selected_results <- bind_rows(all_selected_results)


# Prepare plotting data
#_________________________________________________________________

plot_data <- all_selected_results %>%
  transmute(
    Genes = Gene,
    CancerType = Cancer,
    log2FC = log2FoldChange,
    padj = padj
  ) %>%
  mutate(
    Genes = factor(Genes, levels = genes_of_interest),
    CancerType = factor(CancerType, levels = names(cancer_projects)),
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj <= 0.05 ~ "*",
      TRUE ~ ""
    ),
    label_position = ifelse(
      log2FC >= 0,
      log2FC + 0.10,
      log2FC - 0.17
    )
  )


# Save plotting input table
#_________________________________________________________________

write.csv(
  plot_data,
  file = file.path(combined_dir, "faceted_barplot_input_log2FC_padj.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Define colors
#_________________________________________________________________

gene_colors <- c(
  "#25015F", "#45419A", "#4D90BD", "#51B3B3",
  "#16A087", "#3CBA94", "#6AE19D", "#E9C820",
  "#E98E2D", "#DC4B08", "#762D5A", "#AE155C",
  "#D63882", "#DF7AA9", "#F5AB90", "#B1ADA9"
)

cancer_colors <- c(
  "#25015F", "#45419A", "#4D90BD", "#51B3B3",
  "#16A087", "#3CBA94", "#6AE19D", "#E9C820",
  "#E98E2D", "#DC4B08", "#762D5A", "#AE155C",
  "#D63882"
)


# Plot 1: cancer-wise faceted bar plot
#_________________________________________________________________
# Each panel is one cancer type.
# X-axis shows genes.
# Y-axis shows log2FC.

plot_by_cancer <- ggplot(
  plot_data,
  aes(
    x = Genes,
    y = log2FC,
    fill = CancerType
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    show.legend = TRUE
  ) +
  geom_text(
    aes(label = sig_label, y = label_position),
    size = 4
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ CancerType,
    nrow = 5,
    ncol = 3
  ) +
  scale_fill_manual(values = cancer_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 6,
      angle = 45,
      hjust = 1,
      colour = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black",
      face = "bold"
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 15,
      colour = "black",
      face = "bold"
    ),
    panel.border = element_rect(
      color = "black",
      linewidth = 1
    ),
    strip.background = element_rect(
      color = "black",
      fill = "lightgray",
      linewidth = 1,
      linetype = "solid"
    ),
    strip.text.x = element_text(
      size = 12,
      color = "black",
      face = "bold"
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(
    y = "Log2FC",
    fill = "Cancer type"
  )


# Save cancer-wise plot
#_________________________________________________________________

ggsave(
  filename = file.path(figures_dir, "faceted_barplot_log2FC_by_cancer.pdf"),
  plot = plot_by_cancer,
  width = 11,
  height = 10
)

ggsave(
  filename = file.path(figures_dir, "faceted_barplot_log2FC_by_cancer.png"),
  plot = plot_by_cancer,
  width = 11,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(figures_dir, "faceted_barplot_log2FC_by_cancer.tiff"),
  plot = plot_by_cancer,
  width = 11,
  height = 10,
  dpi = 300,
  compression = "lzw"
)


# Plot 2: gene-wise faceted bar plot
#_________________________________________________________________
# Each panel is one gene.
# X-axis shows cancer types.
# Y-axis shows log2FC.

plot_by_gene <- ggplot(
  plot_data,
  aes(
    x = CancerType,
    y = log2FC,
    fill = Genes
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    show.legend = TRUE
  ) +
  geom_text(
    aes(label = sig_label, y = label_position),
    size = 4
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ Genes,
    nrow = 6,
    ncol = 3
  ) +
  scale_fill_manual(values = gene_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 7,
      angle = 45,
      hjust = 1,
      colour = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black",
      face = "bold"
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 15,
      colour = "black",
      face = "bold"
    ),
    panel.border = element_rect(
      color = "black",
      linewidth = 1
    ),
    strip.background = element_rect(
      color = "black",
      fill = "lightgray",
      linewidth = 1,
      linetype = "solid"
    ),
    strip.text.x = element_text(
      size = 12,
      color = "black",
      face = "bold.italic"
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(
    y = "Log2FC",
    fill = "Gene"
  )


# Save gene-wise plot
#_________________________________________________________________

ggsave(
  filename = file.path(figures_dir, "faceted_barplot_log2FC_by_gene.pdf"),
  plot = plot_by_gene,
  width = 11,
  height = 13
)

ggsave(
  filename = file.path(figures_dir, "faceted_barplot_log2FC_by_gene.png"),
  plot = plot_by_gene,
  width = 11,
  height = 13,
  dpi = 300
)

ggsave(
  filename = file.path(figures_dir, "faceted_barplot_log2FC_by_gene.tiff"),
  plot = plot_by_gene,
  width = 11,
  height = 13,
  dpi = 300,
  compression = "lzw"
)
