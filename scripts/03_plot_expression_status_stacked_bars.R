# Required libraries
#_________________________________________________________________

library(dplyr)
library(tidyr)
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
  "DROSHA", "DGCR8", "DICER1", "TARBP2",
  "XPO1", "XPO5",
  "AGO1", "AGO2", "AGO3", "AGO4",
  "TNRC6A", "PRKRA",
  "GEMIN4", "DDX5", "DDX17", "DDX20"
)


# Define directories
#_________________________________________________________________

results_dir <- file.path(getwd(), "results")
combined_dir <- file.path(results_dir, "combined_results")
figures_dir <- file.path(results_dir, "figures")

dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# Define thresholds
#_________________________________________________________________
# Up   = significant positive log2FC
# Down = significant negative log2FC
# ns   = not significant

padj_cutoff <- 0.05
log2fc_cutoff <- 0


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


# Classify expression status
#_________________________________________________________________

expression_status_data <- all_selected_results %>%
  mutate(
    Cancer = factor(Cancer, levels = names(cancer_projects)),
    Gene = factor(Gene, levels = genes_of_interest),
    Condition = case_when(
      padj < padj_cutoff & log2FoldChange > log2fc_cutoff  ~ "Up",
      padj < padj_cutoff & log2FoldChange < -log2fc_cutoff ~ "Down",
      TRUE ~ "ns"
    ),
    Condition = factor(Condition, levels = c("Down", "Up", "ns"))
  )


# Save expression status table
#_________________________________________________________________

write.csv(
  expression_status_data,
  file = file.path(combined_dir, "expression_status_16_genes_13_cancers.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Count expression status per cancer type
#_________________________________________________________________

stacked_cancers <- expression_status_data %>%
  count(Cancer, Condition, name = "Number_of_genes") %>%
  complete(
    Cancer,
    Condition = factor(c("Down", "Up", "ns"), levels = c("Down", "Up", "ns")),
    fill = list(Number_of_genes = 0)
  )

write.csv(
  stacked_cancers,
  file = file.path(combined_dir, "stacked_bar_input_by_cancer.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Count expression status per gene
#_________________________________________________________________

stacked_genes <- expression_status_data %>%
  count(Gene, Condition, name = "Number_of_cancer_types") %>%
  complete(
    Gene,
    Condition = factor(c("Down", "Up", "ns"), levels = c("Down", "Up", "ns")),
    fill = list(Number_of_cancer_types = 0)
  )

write.csv(
  stacked_genes,
  file = file.path(combined_dir, "stacked_bar_input_by_gene.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Define colors
#_________________________________________________________________

condition_colors <- c(
  "Down" = "#20639B",
  "Up"   = "#ED553B",
  "ns"   = "#FDB93B"
)


# Plot A: cancer-wise stacked bar plot
#_________________________________________________________________

plot_cancer_wise <- ggplot(
  stacked_cancers,
  aes(
    x = Cancer,
    y = Number_of_genes,
    fill = Condition,
    label = Number_of_genes
  )
) +
  geom_bar(
    position = "stack",
    stat = "identity",
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    position = position_stack(vjust = 0.5),
    size = 4.5,
    color = "black"
  ) +
  scale_fill_manual(values = condition_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 12,
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
    axis.title.y = element_text(
      size = 14,
      colour = "black",
      face = "bold"
    ),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(
    y = "Number of genes",
    fill = "Expression status"
  )


# Save cancer-wise plot
#_________________________________________________________________

ggsave(
  filename = file.path(figures_dir, "stacked_bar_expression_status_by_cancer.pdf"),
  plot = plot_cancer_wise,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(figures_dir, "stacked_bar_expression_status_by_cancer.png"),
  plot = plot_cancer_wise,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figures_dir, "stacked_bar_expression_status_by_cancer.tiff"),
  plot = plot_cancer_wise,
  width = 8,
  height = 5,
  dpi = 300,
  compression = "lzw"
)


# Plot B: gene-wise stacked bar plot
#_________________________________________________________________

plot_gene_wise <- ggplot(
  stacked_genes,
  aes(
    x = Gene,
    y = Number_of_cancer_types,
    fill = Condition,
    label = Number_of_cancer_types
  )
) +
  geom_bar(
    position = "stack",
    stat = "identity",
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    position = position_stack(vjust = 0.5),
    size = 4.5,
    color = "black"
  ) +
  scale_fill_manual(values = condition_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 12,
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
    axis.title.y = element_text(
      size = 14,
      colour = "black",
      face = "bold"
    ),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(
    y = "Number of cancer types",
    fill = "Expression status"
  )


# Save gene-wise plot
#_________________________________________________________________

ggsave(
  filename = file.path(figures_dir, "stacked_bar_expression_status_by_gene.pdf"),
  plot = plot_gene_wise,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(figures_dir, "stacked_bar_expression_status_by_gene.png"),
  plot = plot_gene_wise,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figures_dir, "stacked_bar_expression_status_by_gene.tiff"),
  plot = plot_gene_wise,
  width = 8,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

