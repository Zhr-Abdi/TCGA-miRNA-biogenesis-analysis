# Download TCGA miRNA expression quantification data
# Project: TCGA-miRNA-biogenesis-analysis
#
# This script downloads TCGA miRNA Expression Quantification data
# for the selected TCGA cancer types using TCGAbiolinks.
#
# Data source:
# GDC / TCGA
#
# Data type:
# miRNA Expression Quantification
#
# Workflow:
# BCGSC miRNA Profiling
#
# Output:
# data/raw/TCGA_miRNA/
# results/tables/miRNA/tcga_miRNA_download_log.csv


# Required libraries
#_________________________________________________________________

library(TCGAbiolinks)
library(SummarizedExperiment)


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


# Choose cancer types to download
#_________________________________________________________________

# To download all 13 cancer types, use:
selected_cancers <- names(cancer_projects)

# To download only one cancer type, use for example:
# selected_cancers <- c("PRAD")

# To download selected cancer types, use for example:
# selected_cancers <- c("COAD", "BRCA", "LIHC")


# Define output directories
#_________________________________________________________________

miRNA_raw_directory <- file.path(getwd(), "data", "raw", "TCGA_miRNA")

miRNA_log_directory <- file.path(getwd(), "results", "tables", "miRNA")

dir.create(
  miRNA_raw_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  miRNA_log_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# Check selected cancer types
#_________________________________________________________________

if (length(selected_cancers) == 0) {
  stop("No cancer type was selected. Please check selected_cancers.")
}

if (!all(selected_cancers %in% names(cancer_projects))) {
  stop("One or more selected cancer codes are not valid.")
}


# Create empty download log
#_________________________________________________________________

download_log <- data.frame(
  Cancer = character(),
  Project = character(),
  Number_of_files = integer(),
  Download_status = character(),
  Message = character(),
  stringsAsFactors = FALSE
)


# Download TCGA miRNA data
#_________________________________________________________________

for (cancer_code in selected_cancers) {
  
  project_id <- cancer_projects[cancer_code]
  
  message("------------------------------------------------------------")
  message("Starting miRNA download for: ", project_id)
  message("------------------------------------------------------------")
  
  tryCatch({
    
    query <- GDCquery(
      project = project_id,
      data.category = "Transcriptome Profiling",
      data.type = "miRNA Expression Quantification",
      workflow.type = "BCGSC miRNA Profiling",
      experimental.strategy = "miRNA-Seq"
    )
    
    query_results <- getResults(query)
    number_of_files <- nrow(query_results)
    
    message("Number of files found for ", project_id, ": ", number_of_files)
    
    if (number_of_files == 0) {
      
      download_log <- rbind(
        download_log,
        data.frame(
          Cancer = cancer_code,
          Project = project_id,
          Number_of_files = number_of_files,
          Download_status = "No files found",
          Message = "No miRNA Expression Quantification files were found.",
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    GDCdownload(
      query,
      method = "api",
      directory = miRNA_raw_directory
    )
    
    download_log <- rbind(
      download_log,
      data.frame(
        Cancer = cancer_code,
        Project = project_id,
        Number_of_files = number_of_files,
        Download_status = "Downloaded",
        Message = "Download completed successfully.",
        stringsAsFactors = FALSE
      )
    )
    
    message("Finished miRNA download for: ", project_id)
    
  }, error = function(e) {
    
    download_log <- rbind(
      download_log,
      data.frame(
        Cancer = cancer_code,
        Project = project_id,
        Number_of_files = NA_integer_,
        Download_status = "Failed",
        Message = as.character(e$message),
        stringsAsFactors = FALSE
      )
    )
    
    message("Download failed for: ", project_id)
    message("Error message: ", e$message)
  })
}
