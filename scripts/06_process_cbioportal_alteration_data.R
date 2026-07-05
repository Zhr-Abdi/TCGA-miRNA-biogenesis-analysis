
# 06_process_cbioportal_alteration_data.R
#
# Project: TCGA-miRNA-biogenesis-analysis
#
# Purpose:
#   Read raw cBioPortal alteration frequency files downloaded by gene,
#   standardize cancer names and alteration types,
#   add missing zero-frequency rows,
#   and generate processed tables for plotting and reproducibility.
#
# Input:
#   data/raw/cbioportal/alteration_frequency_by_gene/*.txt
#
# Output:
#   data/processed/cbioportal_alteration_frequency_long.csv
#   data/processed/cbioportal_alteration_frequency_wide_frequency.csv
#   data/processed/cbioportal_alteration_frequency_wide_count.csv
#   data/processed/cbioportal_gene_summary.csv
#   data/processed/cbioportal_alteration_type_summary.csv
#   data/processed/cbioportal_alteration_frequency_tables.xlsx
#_________________________________________________________________
#_________________________________________________________________




# Load packages
#_________________________________________________________________

required_packages <- c("tidyverse")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' is required. Please install it first."))
  }
}

library(tidyverse)



# Define paths
#_________________________________________________________________


raw_dir <- "data/raw/cbioportal/alteration_frequency_by_gene"
processed_dir <- "data/processed"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)



# Define genes, cancers, and alteration types
#_________________________________________________________________

genes <- c(
  "DROSHA", "DGCR8", "DICER1", "TARBP2",
  "XPO1", "XPO5",
  "AGO1", "AGO2", "AGO3", "AGO4",
  "TNRC6A", "PRKRA", "GEMIN4",
  "DDX5", "DDX17", "DDX20"
)

cancer_order <- c(
  "BLCA", "BRCA", "COAD", "ESCA", "GBM",
  "KIRC", "LIHC", "LUAD", "LUSC", "PRAD",
  "STAD", "THCA", "UCEC"
)

alteration_type_order <- c(
  "Amplification",
  "Deep deletion",
  "Multiple alterations",
  "Point mutation"
)



# Map cBioPortal cancer study names to short cancer codes
#_________________________________________________________________

# Note:
# cBioPortal uses "Colorectal Adenocarcinoma (TCGA, Firehose Legacy)".
# Here we label it as COAD to match the manuscript figure label.

cancer_map <- tibble::tribble(
  ~Cancer, ~Cancer_Study,
  "BLCA", "Bladder Urothelial Carcinoma (TCGA, Firehose Legacy)",
  "BRCA", "Breast Invasive Carcinoma (TCGA, Firehose Legacy)",
  "COAD", "Colorectal Adenocarcinoma (TCGA, Firehose Legacy)",
  "ESCA", "Esophageal Carcinoma (TCGA, Firehose Legacy)",
  "GBM",  "Glioblastoma Multiforme (TCGA, Firehose Legacy)",
  "KIRC", "Kidney Renal Clear Cell Carcinoma (TCGA, Firehose Legacy)",
  "LIHC", "Liver Hepatocellular Carcinoma (TCGA, Firehose Legacy)",
  "LUAD", "Lung Adenocarcinoma (TCGA, Firehose Legacy)",
  "LUSC", "Lung Squamous Cell Carcinoma (TCGA, Firehose Legacy)",
  "PRAD", "Prostate Adenocarcinoma (TCGA, Firehose Legacy)",
  "STAD", "Stomach Adenocarcinoma (TCGA, Firehose Legacy)",
  "THCA", "Thyroid Carcinoma (TCGA, Firehose Legacy)",
  "UCEC", "Uterine Corpus Endometrial Carcinoma (TCGA, Firehose Legacy)"
)


# Map cBioPortal alteration type labels
#_________________________________________________________________


alteration_type_map <- tibble::tribble(
  ~Alteration_Type_Raw, ~Alteration_Type,
  "amp",      "Amplification",
  "homdel",   "Deep deletion",
  "multiple", "Multiple alterations",
  "mutated",  "Point mutation"
)



# Helper function: mode for inferred sample size
#_________________________________________________________________

get_mode_integer <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_integer_)
  }
  
  as.integer(names(sort(table(x), decreasing = TRUE))[1])
}


# Read one raw cBioPortal file
#_________________________________________________________________

read_cbioportal_file <- function(file_path) {
  
  gene_name <- tools::file_path_sans_ext(basename(file_path))
  gene_name <- toupper(gene_name)
  
  file_ext <- tolower(tools::file_ext(file_path))
  delimiter <- ifelse(file_ext == "csv", ",", "\t")
  
  df <- readr::read_delim(
    file = file_path,
    delim = delimiter,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  expected_cols <- c(
    "Cancer Study",
    "Alteration Frequency",
    "Alteration Type",
    "Alteration Count"
  )
  
  missing_cols <- setdiff(expected_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The file ", basename(file_path),
        " is missing these required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  df %>%
    transmute(
      Gene = gene_name,
      Cancer_Study = `Cancer Study`,
      Alteration_Frequency = as.numeric(`Alteration Frequency`),
      Alteration_Type_Raw = tolower(as.character(`Alteration Type`)),
      Alteration_Count = as.numeric(`Alteration Count`),
      Source_File = basename(file_path)
    )
}


# Read all raw files
#_________________________________________________________________


raw_files <- list.files(
  path = raw_dir,
  pattern = "\\.(txt|csv)$",
  full.names = TRUE,
  ignore.case = TRUE
)

raw_files <- raw_files[!grepl("^README", basename(raw_files), ignore.case = TRUE)]

if (length(raw_files) == 0) {
  stop(paste0("No .txt or .csv files were found in: ", raw_dir))
}

raw_data <- purrr::map_dfr(raw_files, read_cbioportal_file)


# Check expected gene files
#_________________________________________________________________

genes_found <- sort(unique(raw_data$Gene))

missing_genes <- setdiff(genes, genes_found)
extra_genes <- setdiff(genes_found, genes)

if (length(missing_genes) > 0) {
  warning(
    paste0(
      "These expected genes were not found as raw files: ",
      paste(missing_genes, collapse = ", ")
    )
  )
}

if (length(extra_genes) > 0) {
  warning(
    paste0(
      "These extra gene files were found but are not in the predefined gene list: ",
      paste(extra_genes, collapse = ", ")
    )
  )
}


# Standardize cancer names and alteration types
#_________________________________________________________________

raw_standardized <- raw_data %>%
  left_join(cancer_map, by = "Cancer_Study") %>%
  left_join(alteration_type_map, by = "Alteration_Type_Raw") %>%
  filter(Gene %in% genes) %>%
  filter(Cancer %in% cancer_order)


# Check unmapped values
#_________________________________________________________________

unmapped_studies <- raw_standardized %>%
  filter(is.na(Cancer)) %>%
  distinct(Cancer_Study)

if (nrow(unmapped_studies) > 0) {
  print(unmapped_studies)
  stop("Some cancer studies were not mapped. Please update cancer_map.")
}

unmapped_alterations <- raw_standardized %>%
  filter(is.na(Alteration_Type)) %>%
  distinct(Alteration_Type_Raw)

if (nrow(unmapped_alterations) > 0) {
  print(unmapped_alterations)
  stop("Some alteration types were not mapped. Please update alteration_type_map.")
}


# Infer total sample size for each cancer study
#_________________________________________________________________
# Formula:
#   Alteration Frequency = Alteration Count / Total Samples * 100
#
# Therefore:
#   Total Samples = Alteration Count / (Alteration Frequency / 100)

sample_sizes <- raw_standardized %>%
  mutate(
    Total_Samples_Inferred = if_else(
      Alteration_Frequency > 0,
      round(Alteration_Count / (Alteration_Frequency / 100)),
      NA_real_
    )
  ) %>%
  group_by(Cancer, Cancer_Study) %>%
  summarise(
    Total_Samples = get_mode_integer(Total_Samples_Inferred),
    .groups = "drop"
  )

missing_sample_sizes <- sample_sizes %>%
  filter(is.na(Total_Samples))

if (nrow(missing_sample_sizes) > 0) {
  warning("Some cancer sample sizes could not be inferred.")
  print(missing_sample_sizes)
}


# Create complete grid
#_________________________________________________________________
# This step adds missing zero rows.
# For example, if multiple alterations are absent from a raw cBioPortal file,
# they are added here as zero.

complete_grid <- tidyr::expand_grid(
  Gene = genes,
  Cancer = cancer_order,
  Alteration_Type = alteration_type_order
) %>%
  left_join(cancer_map, by = "Cancer")


processed_long <- complete_grid %>%
  left_join(
    raw_standardized %>%
      select(
        Gene,
        Cancer,
        Cancer_Study,
        Alteration_Type,
        Alteration_Frequency,
        Alteration_Count,
        Source_File
      ),
    by = c("Gene", "Cancer", "Cancer_Study", "Alteration_Type")
  ) %>%
  left_join(sample_sizes, by = c("Cancer", "Cancer_Study")) %>%
  mutate(
    Alteration_Frequency = replace_na(Alteration_Frequency, 0),
    Alteration_Count = as.integer(replace_na(Alteration_Count, 0)),
    Source_File = replace_na(Source_File, "Added as zero during processing"),
    Row_Status = if_else(
      Source_File == "Added as zero during processing",
      "Added zero row because absent from raw cBioPortal file",
      "Observed in raw cBioPortal file"
    ),
    Source = "cBioPortal"
  ) %>%
  arrange(
    match(Gene, genes),
    match(Cancer, cancer_order),
    match(Alteration_Type, alteration_type_order)
  )



# Create wide tables
#_________________________________________________________________


processed_wide_frequency <- processed_long %>%
  select(Gene, Cancer, Cancer_Study, Alteration_Type, Alteration_Frequency) %>%
  pivot_wider(
    names_from = Alteration_Type,
    values_from = Alteration_Frequency
  )

processed_wide_count <- processed_long %>%
  select(Gene, Cancer, Cancer_Study, Alteration_Type, Alteration_Count) %>%
  pivot_wider(
    names_from = Alteration_Type,
    values_from = Alteration_Count
  )



# Create summary tables
#_________________________________________________________________

total_tcga_samples <- sample_sizes %>%
  summarise(Total_TCGA_Samples = sum(Total_Samples, na.rm = TRUE)) %>%
  pull(Total_TCGA_Samples)

gene_summary <- processed_long %>%
  group_by(Gene) %>%
  summarise(
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    Overall_Alteration_Frequency_Percent =
      100 * Total_Alteration_Count / total_tcga_samples,
    .groups = "drop"
  ) %>%
  arrange(match(Gene, genes))

alteration_type_summary <- processed_long %>%
  group_by(Alteration_Type) %>%
  summarise(
    Total_Alteration_Count = sum(Alteration_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Percent_of_All_Alterations =
      100 * Total_Alteration_Count / sum(Total_Alteration_Count)
  ) %>%
  arrange(match(Alteration_Type, alteration_type_order))


# Save CSV files
#_________________________________________________________________

readr::write_csv(
  processed_long,
  file.path(processed_dir, "cbioportal_alteration_frequency_long.csv"),
  na = ""
)

readr::write_csv(
  processed_wide_frequency,
  file.path(processed_dir, "cbioportal_alteration_frequency_wide_frequency.csv"),
  na = ""
)

readr::write_csv(
  processed_wide_count,
  file.path(processed_dir, "cbioportal_alteration_frequency_wide_count.csv"),
  na = ""
)

readr::write_csv(
  gene_summary,
  file.path(processed_dir, "cbioportal_gene_summary.csv"),
  na = ""
)

readr::write_csv(
  alteration_type_summary,
  file.path(processed_dir, "cbioportal_alteration_type_summary.csv"),
  na = ""
)
