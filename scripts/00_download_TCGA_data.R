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


# Define output directory
#_________________________________________________________________

CancerType_directory <- file.path(getwd(), "data", "TCGA_raw_data")

dir.create(
  CancerType_directory,
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


# Download TCGA data
#_________________________________________________________________

for (cancer_code in selected_cancers) {
  
  project_id <- cancer_projects[cancer_code]
  
  message("Starting download for: ", project_id)
  
  query <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    experimental.strategy = "RNA-Seq"
  )
  
  GDCdownload(
    query,
    method = "api",
    directory = CancerType_directory
  )
  
  message("Finished download for: ", project_id)
}
