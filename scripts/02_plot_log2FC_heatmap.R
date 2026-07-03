# Required libraries
#_________________________________________________________________

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)


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

# Save combined table
#_________________________________________________________________

write.csv(
  all_selected_results,
  file = file.path(combined_dir, "combined_DESeq2_results_selected_16_genes_all_cancers.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Prepare matrix for classical heatmap
#_________________________________________________________________
# Rows: genes
# Columns: cancer types
# Values: log2FoldChange

heatmap_matrix <- all_selected_results %>%
  select(Cancer, Gene, log2FoldChange) %>%
  mutate(
    Gene = factor(Gene, levels = genes_of_interest),
    Cancer = factor(Cancer, levels = names(cancer_projects))
  ) %>%
  arrange(Gene, Cancer) %>%
  pivot_wider(
    names_from = Cancer,
    values_from = log2FoldChange
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()


# Save log2FC matrix
#_________________________________________________________________

heatmap_matrix_output <- data.frame(
  Gene = rownames(heatmap_matrix),
  heatmap_matrix,
  check.names = FALSE
)

write.csv(
  heatmap_matrix_output,
  file = file.path(combined_dir, "log2FC_matrix_16_genes_13_cancers.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Define heatmap colors centered around zero
#_________________________________________________________________

max_abs_log2fc <- max(abs(heatmap_matrix), na.rm = TRUE)

heatmap_breaks <- seq(
  -max_abs_log2fc,
  max_abs_log2fc,
  length.out = 101
)

heatmap_colors <- colorRampPalette(
  c("#08306B", "#F3E1A0", "#D7191C")
)(100)

# Draw clustered heatmap as PDF
#_________________________________________________________________

pheatmap(
  heatmap_matrix,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  border_color = "grey70",
  main = "log2FC",
  filename = file.path(figures_dir, "clustered_heatmap_log2FC_16_genes_13_cancers.pdf"),
  width = 10,
  height = 7
)


# Draw clustered heatmap as PNG
#_________________________________________________________________

pheatmap(
  heatmap_matrix,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  border_color = "grey70",
  main = "log2FC",
  filename = file.path(figures_dir, "clustered_heatmap_log2FC_16_genes_13_cancers.png"),
  width = 10,
  height = 7
)

# Save clustered heatmap as TIFF
#_________________________________________________________________

tiff(
  filename = file.path(figures_dir, "clustered_heatmap_log2FC_16_genes_13_cancers.tiff"),
  width = 10,
  height = 7,
  units = "in",
  res = 300,
  compression = "lzw"
)

pheatmap(
  heatmap_matrix,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = heatmap_colors,
  breaks = heatmap_breaks,
  border_color = "grey70",
  main = "log2FC"
)

dev.off()
