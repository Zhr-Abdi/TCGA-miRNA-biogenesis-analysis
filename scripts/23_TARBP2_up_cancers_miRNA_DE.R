# ============================================================
# Mature miRNA differential-expression analysis
# in cancer types with TARBP2 upregulation
#
# TCGA mature miRNA:
# Tumor vs Normal
#
# Input:
#   - Existing TCGA BCGSC miRNA Profiling files
#   - Existing mRNA DESeq2 results
#
# Output:
#   - TARBP2-up cancer list
#   - miRNA DESeq2 results for each selected cancer
#   - Combined miRNA DE table
#   - Summary of up/downregulated miRNAs
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)


# ============================================================
# 2. Project directory
# ============================================================

project_dir <- getwd()


# ============================================================
# 3. Input directories
# ============================================================

# mRNA DESeq2 results
mrna_results_dir <- file.path(
  project_dir,
  "01_TCGA_expression",
  "results_by_cancer"
)


# Previously downloaded TCGA miRNA raw data
mirna_raw_dir <- file.path(
  project_dir,
  "Data",
  "TCGA_miR_raw_data"
)


# ============================================================
# 4. Output directory
# ============================================================

output_dir <- file.path(
  project_dir,
  "08_TARBP2_miRNA_analysis"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 5. TCGA cancer projects
# ============================================================

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


# ============================================================
# 6. Thresholds
# ============================================================

# ------------------------------------------------------------
# miRNA differential-expression criterion
# ------------------------------------------------------------

mirna_padj_cutoff <- 0.05
mirna_log2fc_cutoff <- 0.7


# ------------------------------------------------------------
# miRNA filtering
#
# Retain miRNAs with at least 30 reads across all samples
# ------------------------------------------------------------

mirna_min_total_counts <- 30


# ============================================================
# PART 1
# Identify cancer types with TARBP2 upregulation
# ============================================================

tarbp2_results_list <- list()


for (cancer_code in names(cancer_projects)) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  file_path <- file.path(
    mrna_results_dir,
    project_id,
    paste0(
      "DESeq2_results_selected_16_genes_",
      cancer_code,
      ".csv"
    )
  )
  
  
  if (!file.exists(file_path)) {
    
    stop(
      paste(
        "mRNA DESeq2 file not found:",
        file_path
      )
    )
  }
  
  
  df <- read.csv(
    file_path,
    check.names = FALSE
  )
  
  
  tarbp2_row <- df %>%
    
    filter(
      Gene == "TARBP2"
    ) %>%
    
    transmute(
      Cancer = cancer_code,
      log2FoldChange,
      pvalue,
      padj
    )
  
  
  tarbp2_results_list[[cancer_code]] <- tarbp2_row
}


tarbp2_results <- bind_rows(
  tarbp2_results_list
)


# ------------------------------------------------------------
# Classify TARBP2 expression
# ------------------------------------------------------------

tarbp2_results <- tarbp2_results %>%
  
  mutate(
    
    TARBP2_status = case_when(
      
      !is.na(padj) &
        padj < tarbp2_padj_cutoff &
        log2FoldChange >= tarbp2_log2fc_cutoff ~
        "Upregulated",
      
      !is.na(padj) &
        padj < tarbp2_padj_cutoff &
        log2FoldChange <= -tarbp2_log2fc_cutoff ~
        "Downregulated",
      
      TRUE ~
        "Not significantly changed"
    )
  )


print(tarbp2_results)



# Check current columns
colnames(tarbp2_results)

# Create TARBP2 status
tarbp2_results <- tarbp2_results %>%
  dplyr::mutate(
    TARBP2_status = dplyr::case_when(
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange >= 0.7 ~ "Upregulated",
      
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange <= -0.7 ~ "Downregulated",
      
      TRUE ~ "Not significantly changed"
    )
  )

print(tarbp2_results)

# ============================================================
# 7. Select TARBP2-up cancers
# ============================================================

tarbp2_up_cancers <- tarbp2_results %>%
  
  filter(
    TARBP2_status == "Upregulated"
  ) %>%
  
  pull(
    Cancer
  )


cat(
  "\nCancer types with TARBP2 upregulation:\n"
)

print(
  tarbp2_up_cancers
)


cat(
  "\nNumber of selected cancer types:",
  length(tarbp2_up_cancers),
  "\n"
)


# Save selection table
write_csv(
  tarbp2_results,
  file.path(
    output_dir,
    "TARBP2_expression_status_across_TCGA_cancers.csv"
  )
)


# ============================================================
# PART 2
# Function for mature miRNA DE analysis
# ============================================================

run_mirna_de <- function(
    cancer_code,
    project_id
) {
  
  message(
    "\n============================================"
  )
  
  message(
    "Processing ",
    cancer_code,
    " (",
    project_id,
    ")"
  )
  
  message(
    "============================================"
  )
  
  
  # ==========================================================
  # 8. Query GDC metadata
  #
  # No download is performed.
  # Existing local files will be used by GDCprepare().
  # ==========================================================
  
  query <- GDCquery(
    
    project = project_id,
    
    data.category =
      "Transcriptome Profiling",
    
    data.type =
      "Isoform Expression Quantification",
    
    workflow.type =
      "BCGSC miRNA Profiling",
    
    experimental.strategy =
      "miRNA-Seq",
    
    sample.type = c(
      "Primary Tumor",
      "Solid Tissue Normal"
    )
  )
  
  
  # ==========================================================
  # 9. Prepare previously downloaded data
  # ==========================================================
  
  prepared_data <- GDCprepare(
    
    query,
    
    directory = mirna_raw_dir
  )
  
  
  mirna_isoform_data <- as.data.frame(
    prepared_data
  )
  
  
  # ==========================================================
  # 10. Check required columns
  # ==========================================================
  
  required_columns <- c(
    "barcode",
    "miRNA_region",
    "read_count"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(mirna_isoform_data)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "Missing required columns for",
        cancer_code,
        ":",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
  
  
  # ==========================================================
  # 11. Keep mature miRNAs
  # ==========================================================
  
  mature_counts <- mirna_isoform_data %>%
    
    filter(
      grepl(
        "^mature,",
        miRNA_region
      )
    ) %>%
    
    transmute(
      
      Mature_miRNA_ID = sub(
        "^mature,",
        "",
        miRNA_region
      ),
      
      Barcode = barcode,
      
      Read_Count =
        as.numeric(
          read_count
        )
    ) %>%
    
    filter(
      !is.na(Mature_miRNA_ID),
      Mature_miRNA_ID != "",
      !is.na(Barcode)
    ) %>%
    
    group_by(
      Mature_miRNA_ID,
      Barcode
    ) %>%
    
    summarise(
      
      Read_Count = sum(
        Read_Count,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  # ==========================================================
  # 12. Create count matrix
  # ==========================================================
  
  count_table <- mature_counts %>%
    
    pivot_wider(
      
      names_from =
        Barcode,
      
      values_from =
        Read_Count,
      
      values_fill = 0
    )
  
  
  count_matrix <- as.matrix(
    
    count_table[
      ,
      setdiff(
        colnames(count_table),
        "Mature_miRNA_ID"
      ),
      drop = FALSE
    ]
  )
  
  
  rownames(count_matrix) <-
    count_table$Mature_miRNA_ID
  
  
  storage.mode(
    count_matrix
  ) <- "numeric"
  
  
  count_matrix <- round(
    count_matrix
  )
  
  
  # ==========================================================
  # 13. Sample metadata
  # ==========================================================
  
  sample_metadata <- data.frame(
    
    Barcode =
      colnames(
        count_matrix
      ),
    
    Sample_Type_Code =
      substr(
        colnames(
          count_matrix
        ),
        14,
        15
      ),
    
    stringsAsFactors = FALSE
  ) %>%
    
    mutate(
      
      Sample_Type = case_when(
        
        Sample_Type_Code == "01" ~
          "Tumor",
        
        Sample_Type_Code == "11" ~
          "Normal",
        
        TRUE ~
          NA_character_
      )
    ) %>%
    
    filter(
      !is.na(
        Sample_Type
      )
    )
  
  
  # ----------------------------------------------------------
  # Keep only tumor and normal samples
  # ----------------------------------------------------------
  
  count_matrix <- count_matrix[
    ,
    sample_metadata$Barcode,
    drop = FALSE
  ]
  
  
  sample_metadata <-
    sample_metadata %>%
    
    mutate(
      
      Sample_Type = factor(
        Sample_Type,
        levels = c(
          "Normal",
          "Tumor"
        )
      )
    ) %>%
    
    as.data.frame()
  
  
  rownames(
    sample_metadata
  ) <- sample_metadata$Barcode
  
  
  # ==========================================================
  # 14. Check sample numbers
  # ==========================================================
  
  group_counts <- table(
    sample_metadata$Sample_Type
  )
  
  
  print(
    group_counts
  )
  
  
  if (
    
    !all(
      c(
        "Normal",
        "Tumor"
      ) %in%
      names(
        group_counts
      )
    )
    
  ) {
    
    stop(
      paste(
        cancer_code,
        "does not contain both Tumor and Normal samples."
      )
    )
  }
  
  
  if (
    
    any(
      group_counts[
        c(
          "Normal",
          "Tumor"
        )
      ] < 2
    )
    
  ) {
    
    stop(
      paste(
        cancer_code,
        "does not contain at least two samples in each group."
      )
    )
  }
  
  
  # ==========================================================
  # 15. Create DESeq2 object
  # ==========================================================
  
  dds <- DESeqDataSetFromMatrix(
    
    countData =
      count_matrix,
    
    colData =
      sample_metadata,
    
    design =
      ~ Sample_Type
  )
  
  
  # ==========================================================
  # 16. Filter low-count miRNAs
  # ==========================================================
  
  keep_mirnas <- rowSums(
    counts(
      dds
    )
  ) >= mirna_min_total_counts
  
  
  dds <- dds[
    keep_mirnas,
  ]
  
  
  message(
    "Mature miRNAs retained after filtering: ",
    nrow(
      dds
    )
  )
  
  
  # ==========================================================
  # 17. Run DESeq2
  # ==========================================================
  
  dds <- DESeq(
    dds,
    quiet = TRUE
  )
  
  
  # Tumor relative to Normal
  de_result <- results(
    
    dds,
    
    contrast = c(
      "Sample_Type",
      "Tumor",
      "Normal"
    ),
    
    alpha =
      mirna_padj_cutoff
  )
  
  
  # ==========================================================
  # 18. Convert results to data frame
  # ==========================================================
  
  de_table <- as.data.frame(
    de_result
  ) %>%
    
    rownames_to_column(
      "Mature_miRNA_ID"
    ) %>%
    
    mutate(
      
      Cancer =
        cancer_code,
      
      TCGA_Project =
        project_id,
      
      DE_status = case_when(
        
        !is.na(padj) &
          padj < mirna_padj_cutoff &
          log2FoldChange >= mirna_log2fc_cutoff ~
          "Upregulated",
        
        !is.na(padj) &
          padj < mirna_padj_cutoff &
          log2FoldChange <= -mirna_log2fc_cutoff ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      )
    ) %>%
    
    select(
      Cancer,
      TCGA_Project,
      Mature_miRNA_ID,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      DE_status
    ) %>%
    
    arrange(
      padj
    )
  
  
  # ==========================================================
  # 19. Save cancer-specific results
  # ==========================================================
  
  cancer_output_dir <- file.path(
    output_dir,
    cancer_code
  )
  
  
  dir.create(
    cancer_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  write_csv(
    
    de_table,
    
    file.path(
      cancer_output_dir,
      paste0(
        cancer_code,
        "_mature_miRNA_Tumor_vs_Normal_DESeq2.csv"
      )
    )
  )
  
  
  # ==========================================================
  # 20. Save normalized counts for later analyses
  # ==========================================================
  
  normalized_counts <- counts(
    dds,
    normalized = TRUE
  )
  
  
  saveRDS(
    
    normalized_counts,
    
    file.path(
      cancer_output_dir,
      paste0(
        cancer_code,
        "_mature_miRNA_normalized_counts.rds"
      )
    )
  )
  
  
  # ==========================================================
  # 21. Summary
  # ==========================================================
  
  summary_table <- data.frame(
    
    Cancer =
      cancer_code,
    
    Normal_n =
      as.integer(
        group_counts[
          "Normal"
        ]
      ),
    
    Tumor_n =
      as.integer(
        group_counts[
          "Tumor"
        ]
      ),
    
    miRNAs_tested =
      nrow(
        de_table
      ),
    
    miRNAs_up =
      sum(
        de_table$DE_status ==
          "Upregulated",
        na.rm = TRUE
      ),
    
    miRNAs_down =
      sum(
        de_table$DE_status ==
          "Downregulated",
        na.rm = TRUE
      ),
    
    miRNAs_significant =
      sum(
        de_table$DE_status !=
          "Not significant",
        na.rm = TRUE
      )
  )
  
  
  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------
  
  return(
    list(
      DE_table =
        de_table,
      
      Summary =
        summary_table
    )
  )
}


# ============================================================
# PART 3
# Run miRNA DE only for TARBP2-up cancers
# ============================================================

mirna_results <- list()


for (cancer_code in tarbp2_up_cancers) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  mirna_results[[cancer_code]] <- run_mirna_de(
    cancer_code = cancer_code,
    project_id = project_id
  )
}

# ============================================================
# PART 4
# Combine all miRNA DE results
# ============================================================

combined_mirna_de <- bind_rows(
  
  lapply(
    mirna_results,
    function(x) {
      x$DE_table
    }
  )
)

# ============================================================
# Convert MIMAT accessions to mature miRNA names
# and generate readable DE tables
# ============================================================


# ============================================================
# 1. Install/load miRBaseConverter
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("miRBaseConverter", quietly = TRUE)) {
  BiocManager::install("miRBaseConverter")
}

library(miRBaseConverter)
library(dplyr)
library(tidyr)
library(readr)


# ============================================================
# 2. Get unique mature miRNA accessions
# ============================================================

mimat_ids <- unique(
  combined_mirna_de$Mature_miRNA_ID
)

length(mimat_ids)

head(mimat_ids)


# ============================================================
# 3. Convert MIMAT accession -> mature miRNA name
#
# targetVersion = v22
# ============================================================

mirbase_mapping <- miRNA_AccessionToName(
  Accessions = mimat_ids,
  targetVersion = "v22"
)

head(mirbase_mapping)

colnames(mirbase_mapping)


# ============================================================
# 4. Prepare clean mapping table
# ============================================================

mirbase_mapping_clean <- mirbase_mapping %>%
  
  dplyr::transmute(
    Mature_miRNA_ID = Accession,
    Mature_miRNA_Name = TargetName
  ) %>%
  
  dplyr::distinct()


print(
  head(
    mirbase_mapping_clean,
    20
  )
)


# ============================================================
# 5. Add miRNA names to DESeq2 results
# ============================================================

combined_mirna_de_named <- combined_mirna_de %>%
  
  dplyr::left_join(
    mirbase_mapping_clean,
    by = "Mature_miRNA_ID"
  ) %>%
  
  dplyr::relocate(
    Mature_miRNA_Name,
    .after = Mature_miRNA_ID
  )


# Check
print(
  head(
    combined_mirna_de_named,
    20
  )
)


# ============================================================
# 6. Check unmapped MIMAT accessions
# ============================================================

unmapped_mirnas <- combined_mirna_de_named %>%
  
  dplyr::filter(
    is.na(Mature_miRNA_Name) |
      Mature_miRNA_Name == ""
  ) %>%
  
  dplyr::distinct(
    Mature_miRNA_ID
  )


cat(
  "\nNumber of unmapped miRNAs:",
  nrow(unmapped_mirnas),
  "\n"
)

print(unmapped_mirnas)


# ============================================================
# 7. Create clean significant-miRNA direction table
# ============================================================

mirna_direction_table_named <- combined_mirna_de_named %>%
  
  dplyr::filter(
    DE_status %in% c(
      "Upregulated",
      "Downregulated"
    )
  ) %>%
  
  dplyr::select(
    Cancer,
    Mature_miRNA_ID,
    Mature_miRNA_Name,
    log2FoldChange,
    pvalue,
    padj,
    DE_status
  ) %>%
  
  dplyr::arrange(
    Cancer,
    DE_status,
    padj
  )


print(
  mirna_direction_table_named
)


# ============================================================
# 8. Save readable direction table
# ============================================================

write_csv(
  mirna_direction_table_named,
  file.path(
    output_dir,
    "TARBP2_up_cancers_significant_miRNAs_with_names.csv"
  )
)


# ============================================================
# 9. Create direction matrix across cancers
#
# Rows    = mature miRNAs
# Columns = TARBP2-up cancer types
# ============================================================

mirna_direction_matrix_named <- combined_mirna_de_named %>%
  
  dplyr::select(
    Mature_miRNA_ID,
    Mature_miRNA_Name,
    Cancer,
    DE_status
  ) %>%
  
  tidyr::pivot_wider(
    names_from = Cancer,
    values_from = DE_status,
    values_fill = "Not significant"
  )


print(
  mirna_direction_matrix_named
)


# ============================================================
# 10. Save direction matrix
# ============================================================

write_csv(
  mirna_direction_matrix_named,
  file.path(
    output_dir,
    "TARBP2_up_cancers_miRNA_direction_matrix_with_names.csv"
  )
)


# ============================================================
# 11. Optional: log2FC matrix
#
# Useful later for heatmaps
# ============================================================

mirna_log2FC_matrix <- combined_mirna_de_named %>%
  
  dplyr::select(
    Mature_miRNA_ID,
    Mature_miRNA_Name,
    Cancer,
    log2FoldChange
  ) %>%
  
  tidyr::pivot_wider(
    names_from = Cancer,
    values_from = log2FoldChange
  )


write_csv(
  mirna_log2FC_matrix,
  file.path(
    output_dir,
    "TARBP2_up_cancers_miRNA_log2FC_matrix_with_names.csv"
  )
)


# ============================================================
# 12. Save accession-name mapping
# ============================================================

write_csv(
  mirbase_mapping_clean,
  file.path(
    output_dir,
    "MIMAT_to_mature_miRNA_name_mapping_miRBase_v22.csv"
  )
)

# ============================================================
# Export significant miRNAs to one Excel file
# with one sheet per cancer
# ============================================================

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(openxlsx)
library(dplyr)


# Create workbook
wb <- createWorkbook()


# Cancer order
cancer_list <- c(
  "COAD",
  "BLCA",
  "BRCA",
  "LUSC",
  "UCEC",
  "GBM"
)


# Create one sheet for each cancer
for (cancer_code in cancer_list) {
  
  cancer_data <- mirna_direction_table_named %>%
    filter(Cancer == cancer_code) %>%
    select(
      Mature_miRNA_Name,
      Mature_miRNA_ID,
      log2FoldChange,
      pvalue,
      padj,
      DE_status
    ) %>%
    arrange(
      DE_status,
      desc(abs(log2FoldChange))
    )
  
  
  # Add sheet
  addWorksheet(
    wb,
    sheetName = cancer_code
  )
  
  
  # Write data
  writeData(
    wb,
    sheet = cancer_code,
    x = cancer_data
  )
  
  
  # Adjust column width
  setColWidths(
    wb,
    sheet = cancer_code,
    cols = 1:ncol(cancer_data),
    widths = "auto"
  )
  
  
  # Freeze first row
  freezePane(
    wb,
    sheet = cancer_code,
    firstRow = TRUE
  )
}


# Save Excel file
saveWorkbook(
  wb,
  file = file.path(
    output_dir,
    "TARBP2_up_cancers_significant_miRNAs_by_cancer.xlsx"
  ),
  overwrite = TRUE
)

