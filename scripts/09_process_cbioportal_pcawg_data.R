# ============================================================
# Script 09: Process PCAWG cBioPortal alteration frequency data
# Project: TCGA-miRNA-biogenesis-analysis
#
# This script processes PCAWG cBioPortal alteration frequency files
# for miRNA biogenesis-related genes.
#
# Input:
#   data/raw/cbioportal/pcawg/alteration_frequency_by_gene/*.txt
#
# Output:
#   data/processed/pcawg_alteration_frequency_long.csv
#   data/processed/pcawg_alteration_frequency_wide_frequency.csv
#   data/processed/pcawg_alteration_frequency_wide_count.csv
#   data/processed/pcawg_gene_summary.csv
#   data/processed/pcawg_alteration_type_summary.csv
#   data/processed/pcawg_alteration_frequency_tables.xlsx
#   results/tables/pcawg_inferred_sample_sizes.csv
#   results/tables/pcawg_excluded_cancer_types.csv
# ============================================================


# -----------------------------
# 1. Load required packages
# -----------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "openxlsx"
)

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are required but not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(openxlsx)



# 2. Define input/output paths
#_________________________________________________________________


input_dir <- "data/raw/cbioportal/pcawg/alteration_frequency_by_gene"

processed_dir <- "data/processed"
tables_dir <- "results/tables"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)



# 3. Define genes and PCAWG cancer types
#_________________________________________________________________


genes_of_interest <- c(
  "DROSHA", "DGCR8", "DICER1", "TARBP2",
  "XPO1", "XPO5",
  "AGO1", "AGO2", "AGO3", "AGO4",
  "TNRC6A", "PRKRA", "GEMIN4",
  "DDX5", "DDX17", "DDX20"
)

selected_pcawg_cancers <- c(
  "Bladder Cancer",
  "Breast Cancer",
  "Colorectal Cancer",
  "Non-Small Cell Lung Cancer",
  "Uterine Endometrioid Carcinoma",
  "Hepatobiliary Cancer",
  "Esophagogastric Cancer",
  "Renal Cell Carcinoma",
  "Prostate Cancer",
  "Thyroid Cancer"
)

pcawg_tcga_mapping <- tibble(
  Cancer_Type = selected_pcawg_cancers,
  TCGA_Comparable_Cancer = c(
    "BLCA",
    "BRCA",
    "COAD",
    "LUAD/LUSC",
    "UCEC",
    "LIHC",
    "ESCA/STAD",
    "KIRC",
    "PRAD",
    "THCA"
  )
)

alteration_type_map <- c(
  "amp" = "Amplification",
  "homdel" = "Deep deletion",
  "multiple" = "Multiple alterations",
  "mutated" = "Point mutation"
)

alteration_type_levels <- c(
  "Amplification",
  "Deep deletion",
  "Multiple alterations",
  "Point mutation"
)


# 4. Check input files
#_________________________________________________________________

input_files <- list.files(
  input_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

if (length(input_files) == 0) {
  stop("No .txt files were found in: ", input_dir)
}

input_gene_names <- input_files %>%
  basename() %>%
  str_remove("\\.txt$")

missing_gene_files <- setdiff(genes_of_interest, input_gene_names)
unexpected_gene_files <- setdiff(input_gene_names, genes_of_interest)

if (length(missing_gene_files) > 0) {
  stop(
    "The following expected gene files are missing:\n",
    paste(missing_gene_files, collapse = ", "),
    "\n\nExpected file format: GENE.txt, for example AGO2.txt"
  )
}

if (length(unexpected_gene_files) > 0) {
  stop(
    "The following unexpected file names were found:\n",
    paste(unexpected_gene_files, collapse = ", "),
    "\n\nPlease rename files cleanly as GENE.txt, for example AGO2.txt, not AGO2(1).txt."
  )
}



# 5. Read and process one file
#_________________________________________________________________

read_pcawg_file <- function(file_path) {
  
  gene_name <- file_path %>%
    basename() %>%
    str_remove("\\.txt$")
  
  df <- read_tsv(
    file_path,
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "Cancer Type",
    "Alteration Frequency",
    "Alteration Type",
    "Alteration Count"
  )
  
  missing_columns <- setdiff(required_columns, colnames(df))
  
  if (length(missing_columns) > 0) {
    stop(
      "File ", basename(file_path), " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  df %>%
    transmute(
      Gene = gene_name,
      Cancer_Type = `Cancer Type`,
      Alteration_Frequency = parse_number(`Alteration Frequency`),
      Alteration_Type_Raw = str_to_lower(str_trim(`Alteration Type`)),
      Alteration_Count = parse_number(`Alteration Count`)
    )
}



# 6. Combine all PCAWG files
#_________________________________________________________________


pcawg_raw <- map_dfr(input_files, read_pcawg_file)

all_pcawg_cancer_types <- pcawg_raw %>%
  distinct(Cancer_Type) %>%
  arrange(Cancer_Type)

excluded_cancer_types <- pcawg_raw %>%
  filter(!Cancer_Type %in% selected_pcawg_cancers) %>%
  distinct(Gene, Cancer_Type) %>%
  arrange(Cancer_Type, Gene)

excluded_cancer_type_summary <- excluded_cancer_types %>%
  count(Cancer_Type, name = "Number_of_Gene_Files") %>%
  arrange(Cancer_Type)

write_csv(
  excluded_cancer_type_summary,
  file.path(tables_dir, "pcawg_excluded_cancer_types.csv")
)



# 7. Keep selected PCAWG cancer types
#_________________________________________________________________


pcawg_selected <- pcawg_raw %>%
  filter(Cancer_Type %in% selected_pcawg_cancers)

unexpected_alteration_types <- setdiff(
  unique(pcawg_selected$Alteration_Type_Raw),
  names(alteration_type_map)
)

if (length(unexpected_alteration_types) > 0) {
  stop(
    "Unexpected alteration types were found:\n",
    paste(unexpected_alteration_types, collapse = ", "),
    "\n\nPlease check the cBioPortal files before continuing."
  )
}

pcawg_selected <- pcawg_selected %>%
  mutate(
    Alteration_Type = alteration_type_map[Alteration_Type_Raw],
    Total_Samples_Inferred_Raw = if_else(
      Alteration_Frequency > 0 & Alteration_Count > 0,
      Alteration_Count / (Alteration_Frequency / 100),
      NA_real_
    )
  ) %>%
  select(
    Gene,
    Cancer_Type,
    Alteration_Type,
    Alteration_Frequency,
    Alteration_Count,
    Total_Samples_Inferred_Raw
  )


# 8. Infer sample size per PCAWG cancer type
#_________________________________________________________________

pcawg_inferred_sample_sizes <- pcawg_selected %>%
  filter(!is.na(Total_Samples_Inferred_Raw)) %>%
  group_by(Cancer_Type) %>%
  summarise(
    Inferred_Total_Samples = round(median(Total_Samples_Inferred_Raw, na.rm = TRUE)),
    Minimum_Row_Level_Estimate = round(min(Total_Samples_Inferred_Raw, na.rm = TRUE)),
    Maximum_Row_Level_Estimate = round(max(Total_Samples_Inferred_Raw, na.rm = TRUE)),
    Number_of_Rows_Used = n(),
    .groups = "drop"
  ) %>%
  left_join(pcawg_tcga_mapping, by = "Cancer_Type") %>%
  arrange(match(Cancer_Type, selected_pcawg_cancers))

write_csv(
  pcawg_inferred_sample_sizes,
  file.path(tables_dir, "pcawg_inferred_sample_sizes.csv")
)


# 9. Add missing zero rows
#_________________________________________________________________

all_expected_combinations <- expand_grid(
  Gene = genes_of_interest,
  Cancer_Type = selected_pcawg_cancers,
  Alteration_Type = alteration_type_levels
)

pcawg_long <- all_expected_combinations %>%
  left_join(
    pcawg_selected,
    by = c("Gene", "Cancer_Type", "Alteration_Type")
  ) %>%
  mutate(
    Alteration_Frequency = replace_na(Alteration_Frequency, 0),
    Alteration_Count = replace_na(Alteration_Count, 0)
  ) %>%
  left_join(
    pcawg_inferred_sample_sizes %>%
      select(Cancer_Type, TCGA_Comparable_Cancer, Inferred_Total_Samples),
    by = "Cancer_Type"
  ) %>%
  mutate(
    Total_Samples_Inferred = Inferred_Total_Samples,
    Dataset = "PCAWG",
    Gene = factor(Gene, levels = genes_of_interest),
    Cancer_Type = factor(Cancer_Type, levels = selected_pcawg_cancers),
    Alteration_Type = factor(Alteration_Type, levels = alteration_type_levels)
  ) %>%
  arrange(Gene, Cancer_Type, Alteration_Type) %>%
  select(
    Dataset,
    Gene,
    Cancer_Type,
    TCGA_Comparable_Cancer,
    Alteration_Type,
    Alteration_Frequency,
    Alteration_Count,
    Total_Samples_Inferred
  )


# 10. Create wide tables
#_________________________________________________________________

pcawg_wide_frequency <- pcawg_long %>%
  select(
    Gene,
    Cancer_Type,
    TCGA_Comparable_Cancer,
    Alteration_Type,
    Alteration_Frequency
  ) %>%
  pivot_wider(
    names_from = Alteration_Type,
    values_from = Alteration_Frequency,
    values_fill = 0
  ) %>%
  arrange(Gene, Cancer_Type)

pcawg_wide_count <- pcawg_long %>%
  select(
    Gene,
    Cancer_Type,
    TCGA_Comparable_Cancer,
    Alteration_Type,
    Alteration_Count
  ) %>%
  pivot_wider(
    names_from = Alteration_Type,
    values_from = Alteration_Count,
    values_fill = 0
  ) %>%
  arrange(Gene, Cancer_Type)



# 11. Create summary tables
# -----------------------------

pcawg_gene_by_cancer_summary <- pcawg_long %>%
  group_by(Gene, Cancer_Type, TCGA_Comparable_Cancer) %>%
  summarise(
    Total_Alteration_Frequency = sum(Alteration_Frequency, na.rm = TRUE),
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    Total_Samples_Inferred = first(Total_Samples_Inferred),
    .groups = "drop"
  )

pcawg_gene_summary <- pcawg_gene_by_cancer_summary %>%
  group_by(Gene) %>%
  summarise(
    Mean_Alteration_Frequency = mean(Total_Alteration_Frequency, na.rm = TRUE),
    Median_Alteration_Frequency = median(Total_Alteration_Frequency, na.rm = TRUE),
    Maximum_Alteration_Frequency = max(Total_Alteration_Frequency, na.rm = TRUE),
    Total_Alteration_Count = sum(Total_Alteration_Count, na.rm = TRUE),
    Number_of_Cancer_Types = n_distinct(Cancer_Type),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Alteration_Frequency))

pcawg_alteration_type_summary <- pcawg_long %>%
  group_by(Alteration_Type) %>%
  summarise(
    Mean_Alteration_Frequency = mean(Alteration_Frequency, na.rm = TRUE),
    Median_Alteration_Frequency = median(Alteration_Frequency, na.rm = TRUE),
    Maximum_Alteration_Frequency = max(Alteration_Frequency, na.rm = TRUE),
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Alteration_Count))


# 12. Save CSV outputs
#_________________________________________________________________

write_csv(
  pcawg_long,
  file.path(processed_dir, "pcawg_alteration_frequency_long.csv")
)

write_csv(
  pcawg_wide_frequency,
  file.path(processed_dir, "pcawg_alteration_frequency_wide_frequency.csv")
)

write_csv(
  pcawg_wide_count,
  file.path(processed_dir, "pcawg_alteration_frequency_wide_count.csv")
)

write_csv(
  pcawg_gene_summary,
  file.path(processed_dir, "pcawg_gene_summary.csv")
)

write_csv(
  pcawg_alteration_type_summary,
  file.path(processed_dir, "pcawg_alteration_type_summary.csv")
)

