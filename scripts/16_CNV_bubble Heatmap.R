# ==============================================================================
# Pan-cancer sample-level CNV-expression integration
#
# Aim:
# To investigate whether gene-level copy number variation is associated
# with RNA expression of major miRNA biogenesis genes across 13 TCGA cancers.
#
# Analysis:
# Spearman correlation between absolute gene-level copy number
# and log2-normalized RNA expression in matched individual tumors.
#
# Cancer types:
# BLCA, BRCA, COAD, ESCA, GBM, KIRC, LIHC, LUAD, LUSC,
# PRAD, STAD, THCA, UCEC
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


# Define TCGA projects
#_________________________________________________________________

projects <- c(
  "TCGA-BLCA",
  "TCGA-BRCA",
  "TCGA-COAD",
  "TCGA-ESCA",
  "TCGA-GBM",
  "TCGA-KIRC",
  "TCGA-LIHC",
  "TCGA-LUAD",
  "TCGA-LUSC",
  "TCGA-PRAD",
  "TCGA-STAD",
  "TCGA-THCA",
  "TCGA-UCEC"
)


# Define genes of interest
#_________________________________________________________________

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

# Gene-level CNV data
cnv_data_dir <- file.path(
  getwd(),
  "Data",
  "TCGA_CNV_data"
)

# Main analysis directory
analysis_dir <- file.path(
  getwd(),
  "02_CNV_expression"
)

# Cancer-specific results
results_root_dir <- file.path(
  analysis_dir,
  "results_by_cancer"
)

# Combined pan-cancer results
combined_results_dir <- file.path(
  analysis_dir,
  "combined_results"
)

# Define figure directory

figures_dir <- file.path(
  getwd(),
  "Figures",
  "06_TCGA_CNV"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  cnv_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_root_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  combined_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Check input expression data directory
#_________________________________________________________________

if (!dir.exists(tcga_data_dir)) {
  stop(
    paste0(
      "The TCGA raw data directory does not exist:\n",
      tcga_data_dir,
      "\nPlease run the TCGA download script first."
    )
  )
}


# Helper function: identify gene-symbol column
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


# Function: prepare RNA-seq expression data
#_________________________________________________________________

prepare_expression_data <- function(project) {
  
  message(
    "\nPreparing expression data for ",
    project,
    "..."
  )
  
  
  expression_project_dir <- file.path(
    tcga_data_dir,
    project
  )
  
  
  # Query RNA-seq metadata
  #_________________________________________________________________
  
  query_expression <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  
  
  # Ensure that all files required by the current GDC query are available
  #_________________________________________________________________
  
  message(
    "Checking local RNA-seq files for ",
    project,
    "..."
  )
  
  GDCdownload(
    query = query_expression,
    method = "api",
    directory = tcga_data_dir,
    files.per.chunk = 20
  )
  
  
  # Prepare as a data.frame rather than SummarizedExperiment.
  #
  # This avoids unnecessary retrieval of clinical metadata,
  # which is not required for the present analysis.
  #_________________________________________________________________
  
  expression_df <- GDCprepare(
    query = query_expression,
    directory = tcga_data_dir,
    summarizedExperiment = FALSE
  )
  
  expression_df <- as.data.frame(
    expression_df
  )
  
  # Keep only Ensembl genes.
  #
  # STAR count files also contain technical summary rows such as
  # N_unmapped and N_multimapping, which should not be included
  # in normalization.
  #_________________________________________________________________
  
  expression_df <- expression_df[
    grepl(
      "^ENSG",
      expression_df$gene_id
    ),
  ]
  
  
  # Identify unstranded raw-count columns
  #_________________________________________________________________
  
  count_columns <- grep(
    "^unstranded_",
    colnames(expression_df),
    value = TRUE
  )
  
  
  if (length(count_columns) == 0) {
    
    stop(
      paste0(
        "No unstranded STAR-count columns were found for ",
        project,
        "."
      )
    )
  }
  
  
  # Build count matrix
  #_________________________________________________________________
  
  raw_counts <- as.matrix(
    expression_df[
      ,
      count_columns
    ]
  )
  
  storage.mode(raw_counts) <- "numeric"
  
  
  # Sample barcodes
  #_________________________________________________________________
  
  sample_barcodes <- sub(
    "^unstranded_",
    "",
    count_columns
  )
  
  colnames(raw_counts) <- sample_barcodes
  
  
  # Gene symbols
  #_________________________________________________________________
  
  expression_gene_symbols <- expression_df$gene_name
  
  
  # Normalize RNA-seq counts
  #_________________________________________________________________
  
  size_factors <- DESeq2::estimateSizeFactorsForMatrix(
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
        "The following genes were not found in expression data for ",
        project,
        ": ",
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
  
  
  # Convert to patient-level long format
  #_________________________________________________________________
  
  expression_long <- as.data.frame(
    t(expression_subset)
  ) %>%
    
    tibble::rownames_to_column(
      "sample_barcode"
    ) %>%
    
    dplyr::mutate(
      patient_id = substr(
        sample_barcode,
        1,
        12
      )
    ) %>%
    
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        genes_of_interest
      ),
      names_to = "gene",
      values_to = "expression"
    ) %>%
    
    dplyr::group_by(
      patient_id,
      gene
    ) %>%
    
    dplyr::summarise(
      expression = mean(
        expression,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  message(
    "Expression patients: ",
    dplyr::n_distinct(
      expression_long$patient_id
    )
  )
  
  
  return(
    expression_long
  )
}

# Function: prepare gene-level CNV data
#_________________________________________________________________

prepare_cnv_data <- function(project) {
  
  message(
    "\nPreparing CNV data for ",
    project,
    "..."
  )
  
  
  cnv_project_dir <- file.path(
    cnv_data_dir,
    project
  )
  
  
  # Query gene-level CNV
  #_________________________________________________________________
  
  query_cnv <- GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number",
    workflow.type = "ABSOLUTE LiftOver",
    access = "open",
    sample.type = "Primary Tumor"
  )
  
  
  # Ensure that the required CNV files are available
  #_________________________________________________________________
  
  GDCdownload(
    query = query_cnv,
    method = "api",
    directory = cnv_data_dir,
    files.per.chunk = 20
  )
  
  
  # Prepare as data.frame.
  #
  # Clinical metadata are not needed for this analysis.
  #_________________________________________________________________
  
  cnv_df <- GDCprepare(
    query = query_cnv,
    directory = cnv_data_dir,
    summarizedExperiment = FALSE
  )
  
  cnv_df <- as.data.frame(
    cnv_df
  )
  
  # Check required gene annotation
  #_________________________________________________________________
  
  if (!"gene_name" %in% colnames(cnv_df)) {
    
    stop(
      paste0(
        "'gene_name' column was not found in CNV data for ",
        project,
        "."
      )
    )
  }
  
  
  # Identify absolute copy-number columns.
  #
  # Exclude min_copy_number and max_copy_number.
  #_________________________________________________________________
  
  copy_number_columns <- grep(
    "_copy_number$",
    colnames(cnv_df),
    value = TRUE
  )
  
  
  copy_number_columns <- copy_number_columns[
    !grepl(
      "_min_copy_number$|_max_copy_number$",
      copy_number_columns
    )
  ]
  
  
  if (length(copy_number_columns) == 0) {
    
    stop(
      paste0(
        "No absolute copy-number columns were found for ",
        project,
        "."
      )
    )
  }
  
  
  # Extract genes of interest
  #_________________________________________________________________
  
  cnv_gene_index <- match(
    genes_of_interest,
    cnv_df$gene_name
  )
  
  
  if (anyNA(cnv_gene_index)) {
    
    missing_genes <- genes_of_interest[
      is.na(cnv_gene_index)
    ]
    
    stop(
      paste0(
        "The following genes were not found in CNV data for ",
        project,
        ": ",
        paste(
          missing_genes,
          collapse = ", "
        )
      )
    )
  }
  
  
  cnv_subset <- cnv_df[
    cnv_gene_index,
    c(
      "gene_name",
      copy_number_columns
    )
  ]
  
  
  # Convert to matrix
  #_________________________________________________________________
  
  copy_number_matrix <- as.matrix(
    cnv_subset[
      ,
      copy_number_columns
    ]
  )
  
  storage.mode(
    copy_number_matrix
  ) <- "numeric"
  
  
  rownames(
    copy_number_matrix
  ) <- cnv_subset$gene_name
  
  
  # Extract sample barcodes from column names
  #_________________________________________________________________
  
  sample_barcodes <- sub(
    "_copy_number$",
    "",
    copy_number_columns
  )
  
  colnames(
    copy_number_matrix
  ) <- sample_barcodes
  
  
  # Convert to long format
  #_________________________________________________________________
  
  cnv_long <- as.data.frame(
    t(copy_number_matrix)
  ) %>%
    
    tibble::rownames_to_column(
      "sample_barcode"
    ) %>%
    
    dplyr::mutate(
      patient_id = substr(
        sample_barcode,
        1,
        12
      )
    ) %>%
    
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        genes_of_interest
      ),
      names_to = "gene",
      values_to = "copy_number"
    ) %>%
    
    dplyr::mutate(
      copy_number = as.numeric(
        copy_number
      )
    ) %>%
    
    dplyr::group_by(
      patient_id,
      gene
    ) %>%
    
    dplyr::summarise(
      copy_number = mean(
        copy_number,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  message(
    "CNV patients: ",
    dplyr::n_distinct(
      cnv_long$patient_id
    )
  )
  
  
  return(
    cnv_long
  )
}

# Function: Spearman correlation
#_________________________________________________________________

calculate_correlations <- function(
    integrated_data
) {
  
  correlation_results <- integrated_data %>%
    
    dplyr::group_by(
      gene
    ) %>%
    
    dplyr::summarise(
      
      n = dplyr::n(),
      
      rho = stats::cor(
        expression,
        copy_number,
        method = "spearman",
        use = "complete.obs"
      ),
      
      p_value = suppressWarnings(
        stats::cor.test(
          expression,
          copy_number,
          method = "spearman",
          exact = FALSE
        )$p.value
      ),
      
      .groups = "drop"
    )
  
  
  return(
    correlation_results
  )
}


# Function: create and save cancer-specific figure
#_________________________________________________________________

save_cancer_plot <- function(
    integrated_data,
    project,
    results_dir
) {
  
  cancer_code <- sub(
    "^TCGA-",
    "",
    project
  )
  
  
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
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.7
      ),
      
      axis.line = element_blank(),
      
      strip.background = element_rect(
        colour = "black",
        fill = "white",
        linewidth = 0.7
      ),
      
      strip.text = element_text(
        size = 11
      )
    )
  
  
  # Save PDF
  
  ggsave(
    filename = file.path(
      results_dir,
      paste0(
        cancer_code,
        "_CNV_expression_relationships.pdf"
      )
    ),
    plot = cnv_expression_plot,
    width = 9,
    height = 6
  )
  
  
  # Save PNG
  
  ggsave(
    filename = file.path(
      results_dir,
      paste0(
        cancer_code,
        "_CNV_expression_relationships.png"
      )
    ),
    plot = cnv_expression_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  
  # Save TIFF
  
  ggsave(
    filename = file.path(
      results_dir,
      paste0(
        cancer_code,
        "_CNV_expression_relationships.tiff"
      )
    ),
    plot = cnv_expression_plot,
    device = "tiff",
    width = 9,
    height = 6,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  
  return(
    cnv_expression_plot
  )
}


# Function: run complete analysis for one cancer
#_________________________________________________________________

run_cnv_expression_analysis <- function(
    project
) {
  
  message(
    "\n\n",
    "============================================================\n",
    "Starting analysis: ",
    project,
    "\n",
    "============================================================"
  )
  
  
  cancer_code <- sub(
    "^TCGA-",
    "",
    project
  )
  
  
  results_dir <- file.path(
    results_root_dir,
    project
  )
  
  
  dir.create(
    results_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # Prepare expression data
  
  expression_long <- prepare_expression_data(
    project
  )
  
  
  # Prepare CNV data
  
  cnv_long <- prepare_cnv_data(
    project
  )
  
  
  # Integrate expression and CNV
  
  integrated_data <- dplyr::inner_join(
    expression_long,
    cnv_long,
    by = c(
      "patient_id",
      "gene"
    )
  )
  
  
  # Remove missing or non-finite values
  
  integrated_data <- dplyr::filter(
    integrated_data,
    is.finite(expression) &
      is.finite(copy_number)
  )
  
  
  # Add cancer information
  
  integrated_data <- integrated_data %>%
    
    dplyr::mutate(
      project = project,
      cancer = cancer_code
    )
  
  
  message(
    "Matched patients: ",
    dplyr::n_distinct(
      integrated_data$patient_id
    )
  )
  
  
  # Matched sample numbers for each gene
  
  matched_samples <- integrated_data %>%
    
    dplyr::count(
      gene,
      name = "n_matched"
    )
  
  
  print(
    matched_samples
  )
  
  
  # Calculate correlations
  
  correlation_results <- calculate_correlations(
    integrated_data
  ) %>%
    
    dplyr::mutate(
      project = project,
      cancer = cancer_code
    )
  
  
  print(
    correlation_results
  )
  
  
  # Save integrated data
  
  readr::write_csv(
    integrated_data,
    file.path(
      results_dir,
      paste0(
        cancer_code,
        "_CNV_expression_matched_data.csv"
      )
    )
  )
  
  
  # Save preliminary correlation results
  # Global FDR will be added after all cancers are analyzed.
  
  readr::write_csv(
    correlation_results,
    file.path(
      results_dir,
      paste0(
        cancer_code,
        "_CNV_expression_correlations_raw.csv"
      )
    )
  )
  
  
  # Save matched sample summary
  
  readr::write_csv(
    matched_samples,
    file.path(
      results_dir,
      paste0(
        cancer_code,
        "_matched_sample_numbers.csv"
      )
    )
  )
  
  
  # Create and save plot
  
  cancer_plot <- save_cancer_plot(
    integrated_data = integrated_data,
    project = project,
    results_dir = results_dir
  )
  
  
  message(
    "Completed: ",
    project
  )
  
  
  return(
    list(
      integrated_data = integrated_data,
      correlation_results = correlation_results,
      matched_samples = matched_samples
    )
  )
}


# Run analysis across all 13 cancer types
#_________________________________________________________________

analysis_results <- list()

failed_projects <- character()


for (project in projects) {
  
  result <- tryCatch(
    
    {
      
      run_cnv_expression_analysis(
        project
      )
      
    },
    
    error = function(e) {
      
      message(
        "\nERROR in ",
        project,
        ":\n",
        conditionMessage(e)
      )
      
      failed_projects <<- c(
        failed_projects,
        project
      )
      
      return(NULL)
    }
  )
  
  analysis_results[[project]] <- result
}

 # Remove failed analyses
#_________________________________________________________________

successful_results <- analysis_results[
  !vapply(
    analysis_results,
    is.null,
    logical(1)
  )
]


# Combine all cancer-specific matched data
#_________________________________________________________________

all_integrated_data <- dplyr::bind_rows(
  lapply(
    successful_results,
    function(x) {
      x$integrated_data
    }
  )
)


# Combine all correlation results
#_________________________________________________________________

all_correlation_results <- dplyr::bind_rows(
  lapply(
    successful_results,
    function(x) {
      x$correlation_results
    }
  )
)


# Global multiple-testing correction
#_________________________________________________________________

# 13 cancers x 6 genes = up to 78 tests.
# BH correction is applied across all cancer-gene combinations.

all_correlation_results <- all_correlation_results %>%
  
  dplyr::mutate(
    
    FDR_global = stats::p.adjust(
      p_value,
      method = "BH"
    ),
    
    significance = dplyr::case_when(
      FDR_global < 0.001 ~ "***",
      FDR_global < 0.01  ~ "**",
      FDR_global < 0.05  ~ "*",
      TRUE               ~ ""
    )
  ) %>%
  
  dplyr::arrange(
    cancer,
    gene
  )


# Print combined results
#_________________________________________________________________

print(
  all_correlation_results,
  n = Inf
)


# Save combined pan-cancer datasets
#_________________________________________________________________

readr::write_csv(
  all_integrated_data,
  file.path(
    combined_results_dir,
    "PanCancer_CNV_expression_matched_data.csv"
  )
)


readr::write_csv(
  all_correlation_results,
  file.path(
    combined_results_dir,
    "PanCancer_CNV_expression_correlations.csv"
  )
)


# Save failed-project information
#_________________________________________________________________

if (length(failed_projects) > 0) {
  
  failed_projects_table <- tibble(
    project = failed_projects
  )
  
  
  readr::write_csv(
    failed_projects_table,
    file.path(
      combined_results_dir,
      "failed_projects.csv"
    )
  )
  
  
  message(
    "\nThe following projects failed:\n",
    paste(
      failed_projects,
      collapse = ", "
    )
  )
  
} else {
  
  message(
    "\nAll 13 cancer types were analyzed successfully."
  )
}


# Prepare pan-cancer heatmap data
#_________________________________________________________________

heatmap_data <- all_correlation_results %>%
  
  dplyr::mutate(
    
    cancer = factor(
      cancer,
      levels = c(
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
    ),
    
    gene = factor(
      gene,
      levels = c(
        "AGO2",
        "DROSHA",
        "XPO1",
        "XPO5",
        "TARBP2",
        "DICER1"
      )
    ),
    
    heatmap_label = paste0(
      sprintf(
        "%.2f",
        rho
      ),
      significance
    )
  )


# Pan-cancer correlation heatmap
#_________________________________________________________________

pan_cancer_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = cancer,
    y = gene,
    fill = rho
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label = heatmap_label
    ),
    size = 3.2
  ) +
  
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Spearman\nrho"
  ) +
  
  labs(
    x = "Cancer type",
    y = "Gene"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    axis.line = element_blank()
  )


print(
  pan_cancer_heatmap
)


# Save pan-cancer heatmap as PDF
#_________________________________________________________________

ggsave(
  filename = file.path(
    combined_results_dir,
    "PanCancer_CNV_expression_heatmap.pdf"
  ),
  plot = pan_cancer_heatmap,
  width = 10,
  height = 5
)


# Save pan-cancer heatmap as PNG
#_________________________________________________________________

ggsave(
  filename = file.path(
    combined_results_dir,
    "PanCancer_CNV_expression_heatmap.png"
  ),
  plot = pan_cancer_heatmap,
  width = 10,
  height = 5,
  dpi = 300
)


# Save pan-cancer heatmap as TIFF
#_________________________________________________________________

ggsave(
  filename = file.path(
    combined_results_dir,
    "PanCancer_CNV_expression_heatmap.tiff"
  ),
  plot = pan_cancer_heatmap,
  device = "tiff",
  width = 10,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# Save R session information
#_________________________________________________________________

capture.output(
  sessionInfo(),
  file = file.path(
    combined_results_dir,
    "sessionInfo.txt"
  )
)


message(
  "\n============================================================\n",
  "Pan-cancer CNV-expression analysis completed.\n",
  "Combined results saved in:\n",
  combined_results_dir,
  "\n============================================================"
)



# Read all cancer-specific correlation results from disk
#_________________________________________________________________

correlation_files <- lapply(
  projects,
  function(project) {
    
    cancer_code <- sub(
      "^TCGA-",
      "",
      project
    )
    
    file.path(
      results_root_dir,
      project,
      paste0(
        cancer_code,
        "_CNV_expression_correlations_raw.csv"
      )
    )
  }
)


all_correlation_results <- dplyr::bind_rows(
  lapply(
    correlation_files,
    readr::read_csv,
    show_col_types = FALSE
  )
)


print(
  all_correlation_results,
  n = Inf
)


# Global multiple-testing correction
#_________________________________________________________________

all_correlation_results <- all_correlation_results %>%
  
  dplyr::mutate(
    
    FDR_global = stats::p.adjust(
      p_value,
      method = "BH"
    ),
    
    significance = dplyr::case_when(
      FDR_global < 0.001 ~ "***",
      FDR_global < 0.01  ~ "**",
      FDR_global < 0.05  ~ "*",
      TRUE               ~ ""
    )
    
  ) %>%
  
  dplyr::arrange(
    cancer,
    gene
  )


print(
  all_correlation_results,
  n = Inf
)
# Define combined result directory
#_________________________________________________________________

combined_results_dir <- file.path(
  getwd(),
  "02_CNV_expression",
  "combined_results"
)


dir.create(
  combined_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save final pan-cancer correlation table
#_________________________________________________________________

readr::write_csv(
  all_correlation_results,
  file.path(
    combined_results_dir,
    "PanCancer_CNV_expression_correlations.csv"
  )
)

# Summarize CNV-expression coupling across cancers
#_________________________________________________________________

gene_summary <- all_correlation_results %>%
  
  dplyr::group_by(
    gene
  ) %>%
  
  dplyr::summarise(
    
    n_cancers = dplyr::n(),
    
    median_rho = median(
      rho,
      na.rm = TRUE
    ),
    
    min_rho = min(
      rho,
      na.rm = TRUE
    ),
    
    max_rho = max(
      rho,
      na.rm = TRUE
    ),
    
    n_positive = sum(
      rho > 0,
      na.rm = TRUE
    ),
    
    n_significant = sum(
      FDR_global < 0.05,
      na.rm = TRUE
    ),
    
    n_positive_significant = sum(
      rho > 0 &
        FDR_global < 0.05,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      median_rho
    )
  )


print(
  gene_summary
)
readr::write_csv(
  gene_summary,
  file.path(
    combined_results_dir,
    "PanCancer_CNV_expression_gene_summary.csv"
  )
)
# Define figure directory
#_________________________________________________________________

figures_dir <- file.path(
  getwd(),
  "Figures",
  "06_TCGA_CNV"
)


dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# Prepare bubble heatmap data
#_________________________________________________________________

bubble_data <- all_correlation_results %>%
  
  dplyr::mutate(
    
    cancer = factor(
      cancer,
      levels = c(
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
    ),
    
    gene = factor(
      gene,
      levels = c(
        "AGO2",
        "DROSHA",
        "XPO1",
        "XPO5",
        "TARBP2",
        "DICER1"
      )
    ),
    
    neg_log10_FDR = -log10(FDR_global)
  )


# Define observed rho range for stronger colour contrast
#_________________________________________________________________

rho_min <- min(bubble_data$rho, na.rm = TRUE)
rho_mid <- median(bubble_data$rho, na.rm = TRUE)
rho_max <- max(bubble_data$rho, na.rm = TRUE)

# Pan-cancer CNV-expression bubble heatmap
#_________________________________________________________________

pan_cancer_bubble_heatmap <- ggplot(
  bubble_data,
  aes(
    x = cancer,
    y = gene
  )
) +
  
  geom_point(
    aes(
      size = neg_log10_FDR,
      fill = rho
    ),
    shape = 21,
    colour = "grey30",
    stroke = 0.45,
    alpha = 0.95
  ) +
  
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(3.5, 14),
    breaks = c(
      20,
      40,
      60,
      80
    )
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = rho_mid,
    limits = c(rho_min, rho_max),
    breaks = round(
      c(rho_min, rho_mid, rho_max),
      2
    ),
    name = "Spearman rho"
  ) +
  
  labs(
    x = "Cancer type",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.45
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.line = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    axis.title.x = element_text(
      size = 11
    ),
    
    legend.position = "top",
    legend.box = "horizontal"
  )


print(
  pan_cancer_bubble_heatmap
)

# Define figure directory
#_________________________________________________________________

figures_dir <- file.path(
  getwd(),
  "Figures",
  "06_TCGA_CNV"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Save as PDF
#_________________________________________________________________

ggsave(
  filename = file.path(
    figures_dir,
    "PanCancer_CNV_expression_bubble_heatmap.pdf"
  ),
  plot = pan_cancer_bubble_heatmap,
  width = 10,
  height = 5.5,
  units = "in"
)


# Save as PNG
#_________________________________________________________________

ggsave(
  filename = file.path(
    figures_dir,
    "PanCancer_CNV_expression_bubble_heatmap.png"
  ),
  plot = pan_cancer_bubble_heatmap,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 300
)


# Save as TIFF
#_________________________________________________________________

ggsave(
  filename = file.path(
    figures_dir,
    "PanCancer_CNV_expression_bubble_heatmap.tiff"
  ),
  plot = pan_cancer_bubble_heatmap,
  device = "tiff",
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)