# Required libraries
#_________________________________________________________________

library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(AnnotationDbi)
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


# Choose cancer types for analysis
#_________________________________________________________________

# To download all 13 cancer types, use:

selected_cancers <- names(cancer_projects)

# To download only one cancer type, use for example:
# selected_cancers <- c("PRAD")
# To download selected cancer types, use for example:
# selected_cancers <- c("COAD", "BRCA", "LIHC")

# Define directories
#_________________________________________________________________

tcga_data_dir <- file.path(getwd(), "data", "TCGA_raw_data")
results_dir <- file.path(getwd(), "results")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Check downloaded data directory
#_________________________________________________________________

if (!dir.exists(tcga_data_dir)) {
  stop("The TCGA raw data directory does not exist. Please run 00_download_TCGA_data.R first.")
}


# Check selected cancer types
#_________________________________________________________________

if (length(selected_cancers) == 0) {
  stop("No cancer type was selected. Please check selected_cancers.")
}

if (!all(selected_cancers %in% names(cancer_projects))) {
  stop("One or more selected cancer codes are not valid. Use codes such as PRAD, COAD, BRCA.")
}


# Define genes of interest
#_________________________________________________________________

genes_of_interest <- c(
  "DROSHA", "DGCR8", "DICER1", "TARBP2",
  "XPO1", "XPO5",
  "AGO1", "AGO2", "AGO3", "AGO4",
  "TNRC6A", "PRKRA",
  "GEMIN4", "DDX5", "DDX17", "DDX20"
)


# Function for DESeq2 analysis
#_________________________________________________________________

run_deseq2_for_cancer <- function(cancer_code, project_id, data_dir, output_dir, genes_of_interest) {
  
  message("Starting DESeq2 analysis for: ", project_id)
  
  cancer_output_dir <- file.path(output_dir, project_id)
  
  dir.create(
    cancer_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  

  # Query TCGA RNA-seq raw count data
  #_________________________________________________________________
  
  query <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    experimental.strategy = "RNA-Seq"
  )
  
  
  # Prepare and merge downloaded files
  #_________________________________________________________________
  
  query_counts <- GDCprepare(
    query,
    directory = data_dir
  )
  
  

  # Extract raw count matrix
  #_________________________________________________________________
  
  
  pri_matrix <- as.data.frame(
    SummarizedExperiment::assay(query_counts)
  )
  
  
  # Convert Ensembl IDs to gene symbols
  #_________________________________________________________________
  
  ens <- rownames(pri_matrix)
  
  ens_clean <- sub("\\..*", "", ens)
  
  ens_to_symbol <- mapIds(
    org.Hs.eg.db,
    keys = ens_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  count_matrix_annotated <- pri_matrix
  
  count_matrix_annotated$Gene.Symbol <- as.character(ens_to_symbol)
  
  
  # Remove genes without valid gene symbols
  #_________________________________________________________________
  
  count_matrix_annotated <- count_matrix_annotated[
    !is.na(count_matrix_annotated$Gene.Symbol) &
      count_matrix_annotated$Gene.Symbol != "",
  ]
  
  

  # Remove duplicated gene symbols
  #_________________________________________________________________
  
  
  count_matrix_annotated <- count_matrix_annotated[
    !duplicated(count_matrix_annotated$Gene.Symbol),
  ]
  
  

  # Create final expression matrix
  #_________________________________________________________________
  
  exp_matrix <- count_matrix_annotated[
    , setdiff(colnames(count_matrix_annotated), "Gene.Symbol")
  ]
  
  rownames(exp_matrix) <- count_matrix_annotated$Gene.Symbol
  
  exp_matrix <- round(as.matrix(exp_matrix))
  
  

  # Identify primary tumor and normal samples
  #_________________________________________________________________
  
  
  dataSmTP <- TCGAquery_SampleTypes(
    getResults(query, cols = "cases"),
    "TP"
  )
  
  dataSmNT <- TCGAquery_SampleTypes(
    getResults(query, cols = "cases"),
    "NT"
  )
  
  dataSmTP <- intersect(dataSmTP, colnames(exp_matrix))
  dataSmNT <- intersect(dataSmNT, colnames(exp_matrix))
  
  message("Number of primary tumor samples for ", project_id, ": ", length(dataSmTP))
  message("Number of normal samples for ", project_id, ": ", length(dataSmNT))
  
  if (length(dataSmTP) == 0 | length(dataSmNT) == 0) {
    stop(paste0("Tumor or normal samples were not found for ", project_id))
  }
  
  

  # Subset expression matrix into normal and tumor samples
  #_________________________________________________________________
  
  col_normal <- exp_matrix[, dataSmNT, drop = FALSE]
  col_tumor  <- exp_matrix[, dataSmTP, drop = FALSE]
  
  sorted_matrix <- cbind(col_normal, col_tumor)
  
  

  # Create sample metadata for DESeq2
  #_________________________________________________________________
  
  group <- factor(
    c(
      rep("Normal", length(dataSmNT)),
      rep("Tumor", length(dataSmTP))
    ),
    levels = c("Normal", "Tumor")
  )
  
  coldata <- data.frame(group = group)
  
  rownames(coldata) <- colnames(sorted_matrix)
  
  

  # Create DESeq2 dataset
  #_________________________________________________________________
  
  dds <- DESeqDataSetFromMatrix(
    countData = sorted_matrix,
    colData = coldata,
    design = ~ group
  )
  
  
  # Filter low-count genes
  #_________________________________________________________________
  
  
  keep <- rowSums(counts(dds)) >= 10
  
  dds <- dds[keep, ]
  
  

  # Run DESeq2 normalization and differential expression analysis
  #_________________________________________________________________
  
  dds <- DESeq(dds)
  
  
  # Save normalized count matrix
  #_________________________________________________________________
  
  normalized_counts <- log2(
    1 + counts(dds, normalized = TRUE)
  )
  
  write.csv(
    normalized_counts,
    file = file.path(cancer_output_dir, paste0("NormalizedCount_", cancer_code, ".csv")),
    quote = FALSE
  )
  
  
  # Differential expression analysis
  #_________________________________________________________________

  # The p-value and adjusted p-value are calculated by DESeq2.
  # The log2 fold change is then shrunk for more stable estimates.
  
  res_raw <- results(
    dds,
    contrast = c("group", "Tumor", "Normal")
  )
  
  res_lfc <- data.frame(
    lfcShrink(
      dds,
      coef = "group_Tumor_vs_Normal",
      res = res_raw,
      type = "normal"
    )
  )
  
  

  # Sort differential expression results by adjusted p-value
  #_________________________________________________________________
  
  
  res_lfc <- res_lfc[
    order(res_lfc$padj, na.last = TRUE),
    ,
    drop = FALSE
  ]
  
  

  # Save differential expression results for all genes
  #_________________________________________________________________
  
  
  res_lfc_output <- cbind(
    Gene = rownames(res_lfc),
    res_lfc
  )
  
  write.csv(
    res_lfc_output,
    file = file.path(cancer_output_dir, paste0("DESeq2_results_all_genes_", cancer_code, ".csv")),
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # Extract the 16 genes of interest
  #_________________________________________________________________
  
  available_genes <- genes_of_interest[
    genes_of_interest %in% rownames(res_lfc)
  ]
  
  missing_genes <- genes_of_interest[
    !genes_of_interest %in% rownames(res_lfc)
  ]
  
  if (length(missing_genes) > 0) {
    warning(
      paste(
        "The following genes were not found in",
        project_id,
        ":",
        paste(missing_genes, collapse = ", ")
      )
    )
  }
  
  my_dif_list <- res_lfc[
    available_genes,
    ,
    drop = FALSE
  ]
  
  my_dif_list_output <- cbind(
    Cancer = cancer_code,
    Gene = rownames(my_dif_list),
    my_dif_list
  )
  
  
  # Save differential expression results for the 16 genes
  #_________________________________________________________________
 
   # Important:
  # The padj column is kept from the DESeq2 result across all genes.
  # It is not recalculated only among the selected 16 genes.
  
  write.csv(
    my_dif_list_output,
    file = file.path(cancer_output_dir, paste0("DESeq2_results_selected_16_genes_", cancer_code, ".csv")),
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Finished DESeq2 analysis for: ", project_id)
}


# Run DESeq2 analysis
#_________________________________________________________________

for (cancer_code in selected_cancers) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  run_deseq2_for_cancer(
    cancer_code = cancer_code,
    project_id = project_id,
    data_dir = tcga_data_dir,
    output_dir = results_dir,
    genes_of_interest = genes_of_interest
  )
}
