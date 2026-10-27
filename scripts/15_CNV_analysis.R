# ==============================================================================
# Sample-level CNV-expression integration in TCGA-BRCA
#
# Aim:
# To investigate whether gene-level copy number variation is associated
# with RNA expression of major miRNA biogenesis genes in individual
# TCGA-BRCA tumor samples.
#
# Analysis:
# Spearman correlation between absolute gene-level copy number
# and log2-normalized RNA expression.
# ==============================================================================


# Load required packages
#_________________________________________________________________

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readr)


# Define project and genes of interest
#_________________________________________________________________

project <- "TCGA-BRCA"

genes_of_interest <- c(
  "DROSHA",
  "DICER1",
  "XPO1",
  "XPO5",
  "TARBP2",
  "AGO2"
)


# Define directories
#_________________________________________________________________

# Previously downloaded TCGA RNA-seq data
tcga_data_dir <- file.path(
  getwd(),
  "Data",
  "TCGA_raw_data"
)

# Project-specific expression folder
expression_project_dir <- file.path(
  tcga_data_dir,
  project
)

# Directory for CNV data
cnv_data_dir <- file.path(
  getwd(),
  "Data",
  "TCGA_CNV_data"
)

# Project-specific CNV folder
cnv_project_dir <- file.path(
  cnv_data_dir,
  project
)

# Results directory
results_dir <- file.path(
  getwd(),
  "02_CNV_expression",
  "results_by_cancer",
  project
)


dir.create(
  cnv_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Check downloaded expression data directory
#_________________________________________________________________

if (!dir.exists(tcga_data_dir)) {
  stop(
    "The TCGA raw data directory does not exist. Please run 00_download_TCGA_data.R first."
  )
}


if (!dir.exists(expression_project_dir)) {
  stop(
    paste0(
      "Expression data for ",
      project,
      " were not found in:\n",
      expression_project_dir
    )
  )
}


# Helper function to identify the gene symbol column
#_________________________________________________________________

get_gene_symbols <- function(se_object) {
  
  annotation <- as.data.frame(
    rowData(se_object)
  )
  
  possible_columns <- c(
    "gene_name",
    "gene_symbol",
    "external_gene_name",
    "symbol",
    "hugo_gene_symbol",
    "Hugo_Symbol"
  )
  
  symbol_column <- possible_columns[
    possible_columns %in% colnames(annotation)
  ]
  
  
  if (length(symbol_column) == 0) {
    
    stop(
      paste0(
        "Gene symbol column could not be identified.\n",
        "Available columns are:\n",
        paste(
          colnames(annotation),
          collapse = ", "
        )
      )
    )
  }
  
  
  annotation[[symbol_column[1]]]
}


# Query RNA-seq metadata
#_________________________________________________________________

query_expression <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)


# Prepare previously downloaded RNA-seq data
#_________________________________________________________________

expression_se <- GDCprepare(
  query = query_expression,
  directory = tcga_data_dir,
  summarizedExperiment = TRUE
)


# Check available expression assays
#_________________________________________________________________

print(
  assayNames(expression_se)
)


if (!"unstranded" %in% assayNames(expression_se)) {
  
  stop(
    paste0(
      "'unstranded' STAR counts were not found.\n",
      "Available assays: ",
      paste(
        assayNames(expression_se),
        collapse = ", "
      )
    )
  )
}


# Extract raw counts
#_________________________________________________________________

raw_counts <- assay(
  expression_se,
  "unstranded"
)


# Extract gene symbols
#_________________________________________________________________

expression_gene_symbols <- get_gene_symbols(
  expression_se
)


# Normalize RNA-seq counts
#_________________________________________________________________

# Size factors are estimated using all genes,
# before extracting the genes of interest.

size_factors <- estimateSizeFactorsForMatrix(
  raw_counts
)


normalized_counts <- sweep(
  raw_counts,
  MARGIN = 2,
  STATS = size_factors,
  FUN = "/"
)


log2_expression <- log2(
  normalized_counts + 1
)


# Extract genes of interest
#_________________________________________________________________

expression_gene_index <- match(
  genes_of_interest,
  expression_gene_symbols
)


if (anyNA(expression_gene_index)) {
  
  missing_genes <- genes_of_interest[
    is.na(expression_gene_index)
  ]
  
  stop(
    paste0(
      "The following genes were not found in the expression data: ",
      paste(
        missing_genes,
        collapse = ", "
      )
    )
  )
}


expression_subset <- log2_expression[
  expression_gene_index,
  ,
  drop = FALSE
]


rownames(expression_subset) <- genes_of_interest


# Convert expression data to long format
#_________________________________________________________________

expression_long <- as.data.frame(
  t(expression_subset)
) %>%
  
  rownames_to_column(
    "sample_barcode"
  ) %>%
  
  mutate(
    patient_id = substr(
      sample_barcode,
      1,
      12
    )
  ) %>%
  
  pivot_longer(
    cols = all_of(
      genes_of_interest
    ),
    names_to = "gene",
    values_to = "expression"
  )


# Average duplicate aliquots from the same patient, if present
#_________________________________________________________________

expression_long <- expression_long %>%
  
  group_by(
    patient_id,
    gene
  ) %>%
  
  summarise(
    expression = mean(
      expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


cat(
  "\nNumber of patients with expression data:",
  n_distinct(expression_long$patient_id),
  "\n"
)


# Query gene-level CNV data
#_________________________________________________________________

query_cnv <- GDCquery(
  project = project,
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  workflow.type = "ABSOLUTE LiftOver",
  access = "open",
  sample.type = "Primary Tumor"
)
# Download CNV data
#_________________________________________________________________

# Download only if the project folder does not already exist

if (!dir.exists(cnv_project_dir)) {
  
  GDCdownload(
    query = query_cnv,
    method = "api",
    directory = cnv_data_dir,
    files.per.chunk = 20
  )
  
} else {
  
  message(
    "CNV data directory already exists. Download skipped."
  )
}


# Prepare CNV data
#_________________________________________________________________

cnv_se <- GDCprepare(
  query = query_cnv,
  directory = cnv_data_dir,
  summarizedExperiment = TRUE
)


# Check available CNV assays
#_________________________________________________________________

print(
  assayNames(cnv_se)
)


if (!"copy_number" %in% assayNames(cnv_se)) {
  
  stop(
    paste0(
      "'copy_number' assay was not found.\n",
      "Available assays: ",
      paste(
        assayNames(cnv_se),
        collapse = ", "
      )
    )
  )
}


# Extract gene-level copy number matrix
#_________________________________________________________________

copy_number_matrix <- assay(
  cnv_se,
  "copy_number"
)


# Extract CNV gene symbols
#_________________________________________________________________

cnv_gene_symbols <- get_gene_symbols(
  cnv_se
)


# Extract genes of interest
#_________________________________________________________________

cnv_gene_index <- match(
  genes_of_interest,
  cnv_gene_symbols
)


if (anyNA(cnv_gene_index)) {
  
  missing_genes <- genes_of_interest[
    is.na(cnv_gene_index)
  ]
  
  stop(
    paste0(
      "The following genes were not found in the CNV data: ",
      paste(
        missing_genes,
        collapse = ", "
      )
    )
  )
}


cnv_subset <- copy_number_matrix[
  cnv_gene_index,
  ,
  drop = FALSE
]


rownames(cnv_subset) <- genes_of_interest


# Convert CNV data to long format
#_________________________________________________________________

cnv_long <- as.data.frame(
  t(cnv_subset)
) %>%
  
  rownames_to_column(
    "sample_barcode"
  ) %>%
  
  mutate(
    patient_id = substr(
      sample_barcode,
      1,
      12
    )
  ) %>%
  
  pivot_longer(
    cols = all_of(
      genes_of_interest
    ),
    names_to = "gene",
    values_to = "copy_number"
  ) %>%
  
  mutate(
    copy_number = as.numeric(
      copy_number
    )
  )


# Average duplicate aliquots from the same patient, if present
#_________________________________________________________________

cnv_long <- cnv_long %>%
  
  group_by(
    patient_id,
    gene
  ) %>%
  
  summarise(
    copy_number = mean(
      copy_number,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


cat(
  "\nNumber of patients with CNV data:",
  n_distinct(cnv_long$patient_id),
  "\n"
)


# Integrate expression and CNV data
#_________________________________________________________________

integrated_data <- inner_join(
  expression_long,
  cnv_long,
  by = c(
    "patient_id",
    "gene"
  )
)

integrated_data <- dplyr::filter(
  integrated_data,
  is.finite(expression) & is.finite(copy_number)
)

cat(
  "\nNumber of patients with matched expression and CNV data:",
  n_distinct(integrated_data$patient_id),
  "\n"
)


# Check matched sample numbers for each gene
#_________________________________________________________________

matched_samples <- integrated_data %>%
  
  count(
    gene,
    name = "n_matched"
  )


print(
  matched_samples
)


# Spearman correlation analysis
#_________________________________________________________________

correlation_results <- integrated_data %>%
  
  group_by(
    gene
  ) %>%
  
  summarise(
    
    n = sum(
      complete.cases(
        expression,
        copy_number
      )
    ),
    
    rho = cor(
      expression,
      copy_number,
      method = "spearman",
      use = "complete.obs"
    ),
    
    p_value = suppressWarnings(
      cor.test(
        expression,
        copy_number,
        method = "spearman",
        exact = FALSE
      )$p.value
    ),
    
    .groups = "drop"
  )


# Multiple-testing correction
#_________________________________________________________________

correlation_results <- correlation_results %>%
  
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  
  arrange(
    FDR
  )


print(
  correlation_results
)


# Save integrated data
#_________________________________________________________________

write_csv(
  integrated_data,
  file.path(
    results_dir,
    "TCGA_BRCA_CNV_expression_matched_data.csv"
  )
)


# Save correlation results
#_________________________________________________________________

write_csv(
  correlation_results,
  file.path(
    results_dir,
    "TCGA_BRCA_CNV_expression_correlations.csv"
  )
)


# Plot CNV-expression relationships
#_________________________________________________________________

cnv_expression_plot <- ggplot(
  integrated_data,
  aes(
    x = copy_number,
    y = expression
  )
) +
  
  geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.4,
    size = 1.4
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 3
  ) +
  
  labs(
    x = "Absolute gene-level copy number",
    y = "log2(normalized RNA expression + 1)"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    # Complete rectangular border around each panel
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    # Remove the default axis lines because panel.border replaces them
    axis.line = element_blank(),
    
    # Keep gene-name boxes clean
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.7
    ),
    
    strip.text = element_text(
      size = 11
    )
  )


print(cnv_expression_plot)

# Save figure
#_________________________________________________________________

ggsave(
  filename = file.path(
    results_dir,
    "TCGA_BRCA_CNV_expression_relationships.pdf"
  ),
  plot = cnv_expression_plot,
  width = 9,
  height = 6
)


ggsave(
  filename = file.path(
    results_dir,
    "TCGA_BRCA_CNV_expression_relationships.png"
  ),
  plot = cnv_expression_plot,
  width = 9,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(
    results_dir,
    "TCGA_BRCA_CNV_expression_relationships.tiff"
  ),
  plot = cnv_expression_plot,
  device = "tiff",
  width = 9,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

