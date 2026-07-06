# Normalize TCGA miRNA expression count data
# Project: TCGA-miRNA-biogenesis-analysis
#
# This script prepares TCGA mature miRNA Expression Quantification data
# downloaded from GDC and generates DESeq2-normalized miRNA count matrices.
#
# Important:
# This script only performs normalization.
# It does not perform differential expression analysis.
#
# Input:
# data/raw/TCGA_miRNA/
#
# Outputs:
# data/processed/miRNA/raw_counts/
# data/processed/miRNA/normalized_counts/
# data/processed/miRNA/sample_annotation/
# results/tables/miRNA/tcga_miRNA_normalization_log.csv


# Required libraries
#_________________________________________________________________

library(TCGAbiolinks)
library(DESeq2)


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


# Choose cancer types for normalization
#_________________________________________________________________

# To process all 13 cancer types, use:
selected_cancers <- names(cancer_projects)

# To process only one cancer type, use for example:
# selected_cancers <- c("COAD")

# To process selected cancer types, use for example:
# selected_cancers <- c("COAD", "BRCA", "LIHC")


# Define input and output directories
#_________________________________________________________________

tcga_miRNA_data_dir <- file.path(getwd(), "data", "raw", "TCGA_miRNA")

processed_miRNA_dir <- file.path(getwd(), "data", "processed", "miRNA")

raw_counts_dir <- file.path(processed_miRNA_dir, "raw_counts")
normalized_counts_dir <- file.path(processed_miRNA_dir, "normalized_counts")
sample_annotation_dir <- file.path(processed_miRNA_dir, "sample_annotation")

miRNA_table_dir <- file.path(getwd(), "results", "tables", "miRNA")

dir.create(
  raw_counts_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  normalized_counts_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  sample_annotation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  miRNA_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Check downloaded data directory
#_________________________________________________________________

if (!dir.exists(tcga_miRNA_data_dir)) {
  stop("The TCGA miRNA raw data directory does not exist. Please run 12_download_TCGA_miRNA_data.R first.")
}


# Check selected cancer types
#_________________________________________________________________

if (length(selected_cancers) == 0) {
  stop("No cancer type was selected. Please check selected_cancers.")
}

if (!all(selected_cancers %in% names(cancer_projects))) {
  stop("One or more selected cancer codes are not valid. Use codes such as COAD, BRCA, PRAD.")
}


# Helper function: extract raw miRNA count matrix
#_________________________________________________________________

extract_miRNA_raw_counts <- function(prepared_data) {
  
  if (!is.data.frame(prepared_data)) {
    stop(
      paste0(
        "For TCGA miRNA data, GDCprepare was expected to return a data.frame, but it returned: ",
        paste(class(prepared_data), collapse = ", ")
      )
    )
  }
  
  message("GDCprepare returned a data.frame. Extracting read_count columns.")
  
  if (!"miRNA_ID" %in% colnames(prepared_data)) {
    stop("The column miRNA_ID was not found in the prepared miRNA data.")
  }
  
  read_count_columns <- grep(
    pattern = "^read_count_",
    x = colnames(prepared_data),
    value = TRUE
  )
  
  if (length(read_count_columns) == 0) {
    stop("No read_count columns were found in the prepared miRNA data.")
  }
  
  raw_count_matrix <- prepared_data[
    ,
    read_count_columns,
    drop = FALSE
  ]
  
  rownames(raw_count_matrix) <- make.unique(
    as.character(prepared_data$miRNA_ID)
  )
  
  colnames(raw_count_matrix) <- sub(
    pattern = "^read_count_",
    replacement = "",
    x = colnames(raw_count_matrix)
  )
  
  raw_count_matrix <- as.matrix(raw_count_matrix)
  
  raw_count_matrix <- apply(
    raw_count_matrix,
    2,
    as.numeric
  )
  
  rownames(raw_count_matrix) <- make.unique(
    as.character(prepared_data$miRNA_ID)
  )
  
  raw_count_matrix <- round(raw_count_matrix)
  
  storage.mode(raw_count_matrix) <- "integer"
  
  return(raw_count_matrix)
}


# Helper function: make sample annotation table
#_________________________________________________________________

make_sample_annotation <- function(sample_barcodes) {
  
  sample_type_code <- substr(sample_barcodes, 14, 15)
  
  sample_group <- ifelse(
    sample_type_code == "01",
    "Tumor",
    ifelse(
      sample_type_code == "11",
      "Normal",
      "Other"
    )
  )
  
  sample_annotation <- data.frame(
    Sample = sample_barcodes,
    Sample_type_code = sample_type_code,
    Group = sample_group,
    stringsAsFactors = FALSE
  )
  
  rownames(sample_annotation) <- sample_barcodes
  
  return(sample_annotation)
}


# Function for miRNA count normalization
#_________________________________________________________________

normalize_miRNA_for_cancer <- function(cancer_code, project_id, data_dir) {
  
  message("------------------------------------------------------------")
  message("Starting miRNA normalization for: ", project_id)
  message("------------------------------------------------------------")
  
  
  # Query downloaded TCGA mature miRNA expression quantification data
  #_________________________________________________________________
  
  query <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    workflow.type = "BCGSC miRNA Profiling",
    experimental.strategy = "miRNA-Seq"
  )
  
  
  # Prepare and merge downloaded files
  #_________________________________________________________________
  
  query_miRNA <- GDCprepare(
    query,
    directory = data_dir
  )
  
  
  # Extract raw miRNA count matrix
  #_________________________________________________________________
  
  raw_count_matrix <- extract_miRNA_raw_counts(query_miRNA)
  
  
  # Create sample annotation
  #_________________________________________________________________
  
  sample_annotation <- make_sample_annotation(
    sample_barcodes = colnames(raw_count_matrix)
  )
  
  
  # Keep only primary tumor and solid tissue normal samples
  #_________________________________________________________________
  
  selected_samples <- sample_annotation$Sample[
    sample_annotation$Group %in% c("Normal", "Tumor")
  ]
  
  selected_samples <- intersect(
    selected_samples,
    colnames(raw_count_matrix)
  )
  
  if (length(selected_samples) == 0) {
    stop(paste0("No primary tumor or solid tissue normal samples were found for ", project_id))
  }
  
  raw_count_matrix <- raw_count_matrix[
    ,
    selected_samples,
    drop = FALSE
  ]
  
  sample_annotation <- sample_annotation[
    selected_samples,
    ,
    drop = FALSE
  ]
  
  
  # Sort samples: Normal first, Tumor second
  #_________________________________________________________________
  
  normal_samples <- rownames(sample_annotation)[
    sample_annotation$Group == "Normal"
  ]
  
  tumor_samples <- rownames(sample_annotation)[
    sample_annotation$Group == "Tumor"
  ]
  
  sorted_samples <- c(normal_samples, tumor_samples)
  
  raw_count_matrix <- raw_count_matrix[
    ,
    sorted_samples,
    drop = FALSE
  ]
  
  sample_annotation <- sample_annotation[
    sorted_samples,
    ,
    drop = FALSE
  ]
  
  message("Number of normal samples for ", project_id, ": ", length(normal_samples))
  message("Number of tumor samples for ", project_id, ": ", length(tumor_samples))
  
  
  # Save raw miRNA count matrix
  #_________________________________________________________________
  
  raw_counts_file <- file.path(
    raw_counts_dir,
    paste0("raw_miRNA_counts_", cancer_code, ".csv")
  )
  
  write.csv(
    raw_count_matrix,
    file = raw_counts_file,
    quote = FALSE
  )
  
  
  # Save sample annotation
  #_________________________________________________________________
  
  sample_annotation_file <- file.path(
    sample_annotation_dir,
    paste0("sample_annotation_miRNA_", cancer_code, ".csv")
  )
  
  write.csv(
    sample_annotation,
    file = sample_annotation_file,
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # Create DESeq2 object for normalization only
  #_________________________________________________________________
  
  dds <- DESeqDataSetFromMatrix(
    countData = raw_count_matrix,
    colData = sample_annotation,
    design = ~ 1
  )
  
  
  # Filter very low-count miRNAs
  #_________________________________________________________________
  
  keep <- rowSums(counts(dds)) >= 10
  
  dds <- dds[keep, ]
  
  
  # Estimate size factors for DESeq2 normalization
  #_________________________________________________________________
  
  # type = "poscounts" is used because miRNA count matrices can be sparse
  # and may contain many zero values.
  
  dds <- estimateSizeFactors(
    dds,
    type = "poscounts"
  )
  
  
  # Extract DESeq2-normalized counts
  #_________________________________________________________________
  
  normalized_counts <- counts(
    dds,
    normalized = TRUE
  )
  
  
  # Log2-transform normalized counts
  #_________________________________________________________________
  
  normalized_counts_log2 <- log2(
    normalized_counts + 1
  )
  
  
  # Save normalized miRNA count matrix
  #_________________________________________________________________
  
  normalized_counts_file <- file.path(
    normalized_counts_dir,
    paste0("normalized_miRNA_counts_log2_", cancer_code, ".csv")
  )
  
  write.csv(
    normalized_counts_log2,
    file = normalized_counts_file,
    quote = FALSE
  )
  
  
  # Create output summary for this cancer
  #_________________________________________________________________
  
  cancer_summary <- data.frame(
    Cancer = cancer_code,
    Project = project_id,
    Total_samples = ncol(raw_count_matrix),
    Normal_samples = length(normal_samples),
    Tumor_samples = length(tumor_samples),
    Total_miRNAs_before_filtering = nrow(raw_count_matrix),
    Total_miRNAs_after_filtering = nrow(normalized_counts_log2),
    Raw_counts_file = raw_counts_file,
    Normalized_counts_file = normalized_counts_file,
    Sample_annotation_file = sample_annotation_file,
    Status = "Completed",
    Message = "miRNA normalization completed successfully.",
    stringsAsFactors = FALSE
  )
  
  message("Finished miRNA normalization for: ", project_id)
  
  return(cancer_summary)
}


# Run miRNA normalization
#_________________________________________________________________

normalization_log <- data.frame(
  Cancer = character(),
  Project = character(),
  Total_samples = integer(),
  Normal_samples = integer(),
  Tumor_samples = integer(),
  Total_miRNAs_before_filtering = integer(),
  Total_miRNAs_after_filtering = integer(),
  Raw_counts_file = character(),
  Normalized_counts_file = character(),
  Sample_annotation_file = character(),
  Status = character(),
  Message = character(),
  stringsAsFactors = FALSE
)


for (cancer_code in selected_cancers) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  result <- tryCatch({
    
    normalize_miRNA_for_cancer(
      cancer_code = cancer_code,
      project_id = project_id,
      data_dir = tcga_miRNA_data_dir
    )
    
  }, error = function(e) {
    
    message("miRNA normalization failed for: ", project_id)
    message("Error message: ", e$message)
    
    data.frame(
      Cancer = cancer_code,
      Project = project_id,
      Total_samples = NA_integer_,
      Normal_samples = NA_integer_,
      Tumor_samples = NA_integer_,
      Total_miRNAs_before_filtering = NA_integer_,
      Total_miRNAs_after_filtering = NA_integer_,
      Raw_counts_file = NA_character_,
      Normalized_counts_file = NA_character_,
      Sample_annotation_file = NA_character_,
      Status = "Failed",
      Message = as.character(e$message),
      stringsAsFactors = FALSE
    )
  })
  
  normalization_log <- rbind(
    normalization_log,
    result
  )
}


message("------------------------------------------------------------")
