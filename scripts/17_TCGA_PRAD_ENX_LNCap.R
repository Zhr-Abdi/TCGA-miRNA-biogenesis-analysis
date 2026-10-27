# ============================================================
# TCGA-PRAD
# Mature miRNA differential expression
# Tumor vs Normal
# ============================================================


# ============================================================
# 1. Libraries
# ============================================================

library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(tidyr)
library(readr)


# ============================================================
# 2. Directories
# ============================================================

project_dir <- getwd()

raw_data_directory <- file.path(
  project_dir,
  "Data",
  "TCGA_miR_raw_data"
)

output_directory <- file.path(
  project_dir,
  "02_TCGA_miRNA",
  "differential_expression",
  "PRAD"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. TCGA project
# ============================================================

project_id <- "TCGA-PRAD"


# ============================================================
# 4. Query GDC metadata
#
# IMPORTANT:
# This does NOT download the files again.
# ============================================================

query <- GDCquery(
  project = project_id,
  data.category = "Transcriptome Profiling",
  data.type = "Isoform Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  experimental.strategy = "miRNA-Seq",
  sample.type = c(
    "Primary Tumor",
    "Solid Tissue Normal"
  )
)


# ============================================================
# 5. Prepare the files that are already downloaded
# ============================================================

prepared_data <- GDCprepare(
  query,
  directory = raw_data_directory
)

mirna_isoform_data <- as.data.frame(
  prepared_data
)


# ============================================================
# 6. Keep mature miRNAs only
#
# Example:
# mature,MIMAT0000062
#
# becomes:
# MIMAT0000062
# ============================================================

mature_counts <- mirna_isoform_data %>%
  
  dplyr::filter(
    grepl(
      "^mature,",
      miRNA_region
    )
  ) %>%
  
  dplyr::transmute(
    
    Mature_miRNA_ID = sub(
      "^mature,",
      "",
      miRNA_region
    ),
    
    Barcode = barcode,
    
    Read_Count = as.numeric(
      read_count
    )
  ) %>%
  
  dplyr::filter(
    !is.na(Mature_miRNA_ID),
    Mature_miRNA_ID != "",
    !is.na(Barcode)
  ) %>%
  
  dplyr::group_by(
    Mature_miRNA_ID,
    Barcode
  ) %>%
  
  dplyr::summarise(
    
    Read_Count = sum(
      Read_Count,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# ============================================================
# 7. Convert to count matrix
#
# Rows    = mature miRNAs
# Columns = samples
# ============================================================

count_table <- mature_counts %>%
  
  tidyr::pivot_wider(
    names_from = Barcode,
    values_from = Read_Count,
    values_fill = 0
  )


count_matrix <- as.matrix(
  count_table[
    ,
    setdiff(
      names(count_table),
      "Mature_miRNA_ID"
    ),
    drop = FALSE
  ]
)


rownames(count_matrix) <- count_table$Mature_miRNA_ID


# DESeq2 requires integer raw counts

count_matrix <- round(
  count_matrix
)

storage.mode(count_matrix) <- "integer"


# ============================================================
# 8. Identify Tumor and Normal samples
#
# 01 = Primary Tumor
# 11 = Solid Tissue Normal
# ============================================================

sample_metadata <- data.frame(
  
  Barcode = colnames(
    count_matrix
  ),
  
  Sample_Type_Code = substr(
    colnames(count_matrix),
    14,
    15
  ),
  
  stringsAsFactors = FALSE
) %>%
  
  dplyr::mutate(
    
    Sample_Type = dplyr::case_when(
      
      Sample_Type_Code == "01" ~ "Tumor",
      
      Sample_Type_Code == "11" ~ "Normal",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  dplyr::filter(
    !is.na(Sample_Type)
  )


# Keep only Tumor and Normal samples

count_matrix <- count_matrix[
  ,
  sample_metadata$Barcode,
  drop = FALSE
]


# Normal is the reference group

sample_metadata$Sample_Type <- factor(
  sample_metadata$Sample_Type,
  levels = c(
    "Normal",
    "Tumor"
  )
)


rownames(sample_metadata) <- sample_metadata$Barcode


# ============================================================
# 9. Check sample numbers
# ============================================================

cat(
  "\nSample numbers:\n"
)

print(
  table(
    sample_metadata$Sample_Type
  )
)


# Check that metadata and count matrix match

stopifnot(
  all(
    rownames(sample_metadata) ==
      colnames(count_matrix)
  )
)


# ============================================================
# 10. Create DESeq2 object
# ============================================================

dds <- DESeqDataSetFromMatrix(
  
  countData = count_matrix,
  
  colData = sample_metadata,
  
  design = ~ Sample_Type
)


# ============================================================
# 11. Remove very low-count miRNAs
#
# Keep a mature miRNA if total raw counts
# across samples are >= 10
#
# This is consistent with your previous pipeline.
# ============================================================

keep_mirnas <- rowSums(
  counts(dds)
) >= 10


dds <- dds[
  keep_mirnas,
]


cat(
  "\nMature miRNAs retained after filtering:",
  nrow(dds),
  "\n"
)


# ============================================================
# 12. Run DESeq2
# ============================================================

dds <- DESeq(
  dds,
  quiet = TRUE
)


# ============================================================
# 13. Tumor vs Normal
#
# Positive log2FC = higher in Tumor
# Negative log2FC = lower in Tumor
# ============================================================

res <- results(
  
  dds,
  
  contrast = c(
    "Sample_Type",
    "Tumor",
    "Normal"
  ),
  
  alpha = 0.05
)


# ============================================================
# 14. Shrink log2 fold changes
# ============================================================

res_shrunk <- lfcShrink(
  
  dds,
  
  contrast = c(
    "Sample_Type",
    "Tumor",
    "Normal"
  ),
  
  res = res,
  
  type = "normal"
)


# ============================================================
# 15. Build final results table
# ============================================================

result_df <- as.data.frame(
  res_shrunk
)


result_df$Mature_miRNA_ID <- rownames(
  result_df
)


rownames(result_df) <- NULL


result_df <- result_df %>%
  
  dplyr::mutate(
    
    Cancer = "PRAD",
    
    Direction = dplyr::case_when(
      
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange > 0 ~
        "Up",
      
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange < 0 ~
        "Down",
      
      TRUE ~
        "Not significant"
    ),
    
    Strong_DE = dplyr::case_when(
      
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange >= 0.5 ~
        "Up",
      
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange <= -0.5 ~
        "Down",
      
      TRUE ~
        "No"
    )
  ) %>%
  
  dplyr::select(
    
    Cancer,
    
    Mature_miRNA_ID,
    
    baseMean,
    
    log2FoldChange,
    
    lfcSE,
    
    stat,
    
    pvalue,
    
    padj,
    
    Direction,
    
    Strong_DE
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ============================================================
# 16. Significant mature miRNAs
#
# FDR < 0.05
# ============================================================

significant_df <- result_df %>%
  
  dplyr::filter(
    !is.na(padj),
    padj < 0.05
  )


# ============================================================
# 17. Save normalized expression
#
# We will need this later for comparing with
# the Enoxacin miRNA signature.
# ============================================================

normalized_counts <- counts(
  dds,
  normalized = TRUE
)


log2_normalized_counts <- log2(
  normalized_counts + 1
)


normalized_df <- data.frame(
  
  Mature_miRNA_ID =
    rownames(log2_normalized_counts),
  
  log2_normalized_counts,
  
  check.names = FALSE
)


# ============================================================
# 18. Save output files
# ============================================================

write_csv(
  
  result_df,
  
  file.path(
    output_directory,
    "PRAD_mature_miRNA_DE_all.csv"
  )
)


write_csv(
  
  significant_df,
  
  file.path(
    output_directory,
    "PRAD_mature_miRNA_DE_significant.csv"
  )
)


write_csv(
  
  normalized_df,
  
  file.path(
    output_directory,
    "PRAD_mature_miRNA_log2_normalized_expression.csv"
  )
)


# ============================================================
# 19. Print summary
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "PRAD miRNA differential expression\n"
)

cat(
  "====================================\n"
)


cat(
  "Total mature miRNAs tested:",
  nrow(result_df),
  "\n"
)


cat(
  "Significant miRNAs (FDR < 0.05):",
  nrow(significant_df),
  "\n"
)


cat(
  "Upregulated:",
  sum(
    significant_df$log2FoldChange > 0
  ),
  "\n"
)


cat(
  "Downregulated:",
  sum(
    significant_df$log2FoldChange < 0
  ),
  "\n"
)


cat(
  "Strongly upregulated (log2FC >= 0.5):",
  sum(
    significant_df$log2FoldChange >= 0.5
  ),
  "\n"
)


cat(
  "Strongly downregulated (log2FC <= -0.5):",
  sum(
    significant_df$log2FoldChange <= -0.5
  ),
  "\n"
)


# ============================================================
# 20. Show top results
# ============================================================

print(
  head(
    result_df,
    20
  )
)

library(readr)
library(dplyr)


# Current working directory
project_dir <- getwd()

# Input file
input_file <- file.path(
  project_dir,
  "02_TCGA_miRNA",
  "differential_expression",
  "LNCaP_ENX_SuppTable1_legacy.csv"
)

# Read data
lncap_enx <- read_csv(
  input_file,
  show_col_types = FALSE
)

# ============================================================
# Quality control of LNCaP enoxacin miRNA dataset
# ============================================================

# Check dimensions
dim(lncap_enx)

# Check column names
colnames(lncap_enx)

# Preview the dataset
head(lncap_enx)

# Number of ENX-up and ENX-down miRNAs
table(lncap_enx$ENX_direction)

# Check duplicated miRNA names
sum(duplicated(lncap_enx$legacy_name_reported))

# Check missing miRNA names
sum(is.na(lncap_enx$legacy_name_reported))

# Check missing fold-change values
sum(is.na(lncap_enx$ENX_fold_change))

# Check missing p-values
sum(is.na(lncap_enx$ENX_pvalue))

# ============================================================
# Step 2: Map legacy miRNA names to miRBase accessions
# ============================================================

# Separate miRBase-compatible miRNAs from miRPlus probes

lncap_mirbase <- lncap_enx %>%
  dplyr::filter(
    !grepl("miRPlus", legacy_name_reported)
  )

lncap_mirplus <- lncap_enx %>%
  dplyr::filter(
    grepl("miRPlus", legacy_name_reported)
  )

library(miRBaseConverter)

mapping_v11 <- miRBaseConverter::miRNA_NameToAccession(
  lncap_mirbase$legacy_name_reported,
  version = "v11"
)
head(mapping_v11)

dim(mapping_v11)

colnames(mapping_v11)

# ============================================================
# Step 2.1: Standardize miRBase mapping table
# ============================================================

mapping_v11 <- mapping_v11 %>%
  dplyr::rename(
    legacy_name_reported = miRNAName_v11,
    MIMAT = Accession
  )

head(mapping_v11)

# ============================================================
# Step 2.2: Check MIMAT mapping
# ============================================================

# Number of miRNAs without a MIMAT accession
sum(is.na(mapping_v11$MIMAT))

# Show unresolved miRNAs, if any
unresolved_miRNAs <- mapping_v11 %>%
  dplyr::filter(is.na(MIMAT))

unresolved_miRNAs

# ============================================================
# Step 2.3: Manually resolve miRNAs not mapped by miRBase v11
# ============================================================

manual_mapping <- tibble::tribble(
  ~legacy_name_reported, ~MIMAT_manual,    ~mapping_note,
  "hsa-miR-675b",        "MIMAT0006790",   "Legacy alias of hsa-miR-675-3p",
  "hsa-miR-1912",        "MIMAT0007887",   "Legacy unsuffixed name; current hsa-miR-1912-3p",
  "hsa-miR-1539",        "MIMAT0007401",   "Mature miRNA accession"
)

mapping_final <- mapping_v11 %>%
  dplyr::left_join(
    manual_mapping,
    by = "legacy_name_reported"
  ) %>%
  dplyr::mutate(
    mapping_method = dplyr::if_else(
      is.na(MIMAT),
      "manual_rescue",
      "miRBaseConverter_v11"
    ),
    MIMAT = dplyr::coalesce(
      MIMAT,
      MIMAT_manual
    )
  ) %>%
  dplyr::select(
    legacy_name_reported,
    MIMAT,
    mapping_method,
    mapping_note
  )
# Check unresolved miRNAs after manual rescue

sum(is.na(mapping_final$MIMAT))

# ============================================================
# Step 3: Match LNCaP ENX miRNAs with TCGA-PRAD
# ============================================================

# Add final MIMAT mapping to the LNCaP ENX dataset

lncap_mapped <- lncap_mirbase %>%
  dplyr::left_join(
    mapping_final,
    by = "legacy_name_reported"
  )

# ============================================================
# Step 3: Match LNCaP ENX miRNAs with TCGA-PRAD
# ============================================================

prad_file <- file.path(
  project_dir,
  "02_TCGA_miRNA",
  "differential_expression",
  "PRAD",
  "PRAD_mature_miRNA_DE_all.csv"
)

prad_de <- readr::read_csv(
  prad_file,
  show_col_types = FALSE
)

# Check duplicated mature miRNA IDs
sum(duplicated(prad_de$Mature_miRNA_ID))

# Check missing mature miRNA IDs
sum(is.na(prad_de$Mature_miRNA_ID))

enx_prad <- lncap_mapped %>%
  dplyr::left_join(
    prad_de,
    by = c("MIMAT" = "Mature_miRNA_ID")
  )

sum(!is.na(enx_prad$log2FoldChange))

sum(is.na(enx_prad$log2FoldChange))

unmatched_prad <- enx_prad %>%
  dplyr::filter(is.na(log2FoldChange)) %>%
  dplyr::select(
    legacy_name_reported,
    MIMAT,
    ENX_fold_change,
    ENX_direction
  )

unmatched_prad



# ============================================================
# Correct let-7d precursor accession to mature MIMAT accession
# ============================================================

mapping_final <- mapping_final %>%
  dplyr::mutate(
    MIMAT = dplyr::if_else(
      legacy_name_reported == "hsa-let-7d",
      "MIMAT0000065",
      MIMAT
    ),
    mapping_method = dplyr::if_else(
      legacy_name_reported == "hsa-let-7d",
      "manual_rescue",
      mapping_method
    ),
    mapping_note = dplyr::if_else(
      legacy_name_reported == "hsa-let-7d",
      "Corrected precursor MI0000065 to mature hsa-let-7d-5p MIMAT0000065",
      mapping_note
    )
  )

lncap_mapped <- lncap_mirbase %>%
  dplyr::left_join(
    mapping_final,
    by = "legacy_name_reported"
  )
enx_prad <- lncap_mapped %>%
  dplyr::left_join(
    prad_de,
    by = c("MIMAT" = "Mature_miRNA_ID")
  )
unmatched_prad <- enx_prad %>%
  dplyr::filter(is.na(log2FoldChange)) %>%
  dplyr::select(
    legacy_name_reported,
    MIMAT,
    ENX_fold_change,
    ENX_direction
  )

unmatched_prad

# ============================================================
# Step 6: Directional heatmap of ENX vs PRAD miRNA profiles
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# ------------------------------------------------------------
# 1. Prepare data
# ------------------------------------------------------------

heatmap_data <- matched_only %>%
  dplyr::arrange(dplyr::desc(ENX_fold_change)) %>%
  dplyr::select(
    legacy_name_reported,
    ENX_fold_change,
    ENX_direction,
    PRAD_direction
  )


# Keep miRNAs ordered by ENX response
miRNA_order <- heatmap_data$legacy_name_reported


# Convert to long format
heatmap_long <- heatmap_data %>%
  dplyr::select(
    legacy_name_reported,
    ENX_direction,
    PRAD_direction
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      ENX_direction,
      PRAD_direction
    ),
    names_to = "Dataset",
    values_to = "Direction"
  )


# Rename columns for the figure
heatmap_long <- heatmap_long %>%
  dplyr::mutate(
    
    Dataset = dplyr::recode(
      Dataset,
      ENX_direction = "LNCaP + ENX",
      PRAD_direction = "PRAD Tumor vs Normal"
    ),
    
    legacy_name_reported = factor(
      legacy_name_reported,
      levels = rev(miRNA_order)
    )
  )


# ------------------------------------------------------------
# 2. Plot directional heatmap
# ------------------------------------------------------------

p_heatmap <- ggplot(
  heatmap_long,
  aes(
    x = Dataset,
    y = legacy_name_reported,
    fill = Direction
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.2
  ) +
  
  scale_fill_manual(
    values = c(
      "Up" = "#D73027",
      "Down" = "#4575B4"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = "Direction",
    title = "Directional comparison of ENX-induced and PRAD-associated miRNA profiles"
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 5
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    
    legend.position = "right"
  )


p_heatmap

# ------------------------------------------------------------
# 3. Save figure
# ------------------------------------------------------------

figure_dir <- file.path(
  project_dir,
  "02_TCGA_miRNA",
  "differential_expression",
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  filename = file.path(
    figure_dir,
    "LNCaP_ENX_vs_PRAD_directional_heatmap.pdf"
  ),
  plot = p_heatmap,
  width = 5,
  height = 12
)


ggsave(
  filename = file.path(
    figure_dir,
    "LNCaP_ENX_vs_PRAD_directional_heatmap.png"
  ),
  plot = p_heatmap,
  width = 5,
  height = 12,
  dpi = 600
)


# ------------------------------------------------------------
# Save figure as TIFF
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    "LNCaP_ENX_vs_PRAD_directional_heatmap.tiff"
  ),
  plot = p_heatmap,
  width = 5,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# ============================================================
# Step 4: Compare ENX and PRAD miRNA directions
# ============================================================

# Define miRNA direction in PRAD
# Positive log2FC = higher in Tumor than Normal
# Negative log2FC = lower in Tumor than Normal

enx_prad <- enx_prad %>%
  dplyr::mutate(
    
    PRAD_direction = dplyr::case_when(
      is.na(log2FoldChange) ~ NA_character_,
      log2FoldChange > 0 ~ "Up",
      log2FoldChange < 0 ~ "Down",
      TRUE ~ "No_change"
    )
    
  )

table(
  enx_prad$PRAD_direction,
  useNA = "ifany"
)
matched_only <- enx_prad %>%
  dplyr::filter(
    !is.na(PRAD_direction),
    PRAD_direction %in% c("Up", "Down")
  )
matched_only <- enx_prad %>%
  dplyr::filter(
    !is.na(PRAD_direction),
    PRAD_direction %in% c("Up", "Down")
  )

# Determine whether ENX and PRAD change in the same direction

matched_only <- matched_only %>%
  dplyr::mutate(
    
    direction_match = dplyr::case_when(
      
      ENX_direction == PRAD_direction ~ "Concordant",
      
      ENX_direction != PRAD_direction ~ "Discordant"
      
    )
    
  )



