# ============================================================
# Stage-based expression analysis across TCGA cancer types
# Project: TCGA-miRNA-biogenesis-analysis
#
# Input:
#   data/processed/normalized_counts_log2_<CANCER>.csv
#
# Main outputs:
#   results/tables/stage_expression_all_cancers.csv
#   results/tables/stage_ANOVA_all_cancers.csv
#   results/tables/stage_patient_counts_all_cancers.csv
#   results/tables/stage_analysis_log.csv
#
#   results/figures/stage_expression_pan_cancer.png
#   results/figures/stage_expression_pan_cancer.pdf
#   results/figures/stage_expression_pan_cancer.tiff
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

# -----------------------------
# 1. Project paths
# -----------------------------

project_dir <- getwd()

processed_dir <- file.path(project_dir, "data", "processed")
clinical_download_dir <- file.path(project_dir, "GDCdata_clinical")

table_dir <- file.path(project_dir, "results", "tables")
figure_dir <- file.path(project_dir, "results", "figures")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_download_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Cancer types and genes
# -----------------------------

cancers <- c(
  "COAD", "ESCA", "BLCA", "BRCA",
  "LUAD", "LUSC", "KIRC", "LIHC",
  "THCA", "UCEC", "GBM", "STAD", "PRAD"
)

selected_genes <- c(
  "DROSHA", "DGCR8", "DICER1", "TARBP2",
  "XPO1", "XPO5",
  "AGO1", "AGO2", "AGO3", "AGO4",
  "TNRC6A", "PRKRA", "GEMIN4",
  "DDX5", "DDX17", "DDX20"
)

plot_genes <- c("AGO2", "DICER1", "DROSHA", "TARBP2", "XPO1", "XPO5")

plot_cancer_order <- c(
  "STAD", "LIHC", "COAD",
  "UCEC", "THCA", "BLCA",
  "BRCA", "ESCA", "KIRC",
  "LUSC", "LUAD", "PRAD",
  "GBM"
)

# -----------------------------
# 3. Check input files
# -----------------------------

expected_files <- paste0("normalized_counts_log2_", cancers, ".csv")

file_check <- tibble(
  Cancer = cancers,
  file = expected_files,
  path = file.path(processed_dir, expected_files),
  exists = file.exists(file.path(processed_dir, expected_files))
)

cat("\nChecking normalized count files:\n")
print(file_check)

if (any(!file_check$exists)) {
  stop("Some normalized count files are missing. Check the data/processed folder.")
}

# -----------------------------
# 4. Stage cleaning function
# -----------------------------

clean_clinical_stage <- function(stage_vector) {
  
  stage_upper <- toupper(trimws(as.character(stage_vector)))
  
  stage_upper <- gsub("\\[|\\]", "", stage_upper)
  stage_upper <- gsub("_", " ", stage_upper)
  stage_upper <- gsub("\\s+", " ", stage_upper)
  
  stage_upper[
    stage_upper %in% c(
      "",
      "NA",
      "N/A",
      "NOT REPORTED",
      "NOT AVAILABLE",
      "NOT APPLICABLE",
      "UNKNOWN",
      "NONE",
      "DISCREPANCY"
    )
  ] <- NA
  
  stage_group <- rep(NA_character_, length(stage_upper))
  
  stage_group[
    grepl("^STAGE\\s+IV[A-C0-9]*$|^IV[A-C0-9]*$", stage_upper)
  ] <- "IV"
  
  stage_group[
    grepl("^STAGE\\s+III[A-C0-9]*$|^III[A-C0-9]*$", stage_upper)
  ] <- "III"
  
  stage_group[
    grepl("^STAGE\\s+II[A-C0-9]*$|^II[A-C0-9]*$", stage_upper)
  ] <- "II"
  
  stage_group[
    grepl("^STAGE\\s+I[A-C0-9]*$|^I[A-C0-9]*$", stage_upper)
  ] <- "I"
  
  return(stage_group)
}

# -----------------------------
# 5. P-value to significance stars
# -----------------------------

p_to_star <- function(p) {
  
  dplyr::case_when(
    is.na(p) ~ "",
    p <= 0.0001 ~ "****",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ ""
  )
}

# -----------------------------
# 6. Choose best stage column
# -----------------------------

choose_best_stage_column <- function(clinical_patient, cancer) {
  
  known_stage_cols <- c(
    "clinical_stage",
    "pathologic_stage",
    "ajcc_pathologic_tumor_stage",
    "ajcc_clinical_tumor_stage",
    "stage_event_clinical_stage",
    "stage_event_pathologic_stage",
    "stage_other"
  )
  
  all_stage_like_cols <- grep(
    "stage",
    colnames(clinical_patient),
    value = TRUE,
    ignore.case = TRUE
  )
  
  candidate_stage_cols <- unique(c(known_stage_cols, all_stage_like_cols))
  candidate_stage_cols <- candidate_stage_cols[
    candidate_stage_cols %in% colnames(clinical_patient)
  ]
  
  if (length(candidate_stage_cols) == 0) {
    message("No candidate stage columns found for ", cancer)
    return(NULL)
  }
  
  stage_column_summary <- lapply(candidate_stage_cols, function(col_name) {
    
    stage_group <- clean_clinical_stage(clinical_patient[[col_name]])
    
    tibble(
      Cancer = cancer,
      stage_column = col_name,
      n_valid_stage = sum(!is.na(stage_group)),
      n_stage_I = sum(stage_group == "I", na.rm = TRUE),
      n_stage_II = sum(stage_group == "II", na.rm = TRUE),
      n_stage_III = sum(stage_group == "III", na.rm = TRUE),
      n_stage_IV = sum(stage_group == "IV", na.rm = TRUE)
    )
  }) %>%
    bind_rows()
  
  cat("\nCandidate stage columns for", cancer, ":\n")
  print(stage_column_summary)
  
  best_stage_column <- stage_column_summary %>%
    arrange(desc(n_valid_stage)) %>%
    slice(1)
  
  if (best_stage_column$n_valid_stage == 0) {
    message("No valid Stage I-IV values found for ", cancer)
    return(NULL)
  }
  
  message(
    "Selected stage column for ",
    cancer,
    ": ",
    best_stage_column$stage_column,
    " with ",
    best_stage_column$n_valid_stage,
    " valid patients"
  )
  
  return(best_stage_column$stage_column)
}

# -----------------------------
# 7. Download and prepare clinical stage table
# -----------------------------

get_clinical_stage_table <- function(cancer) {
  
  message("\n==============================")
  message("Clinical data: ", cancer)
  message("==============================")
  
  project <- paste0("TCGA-", cancer)
  
  query.clinical <- tryCatch(
    {
      GDCquery(
        project = project,
        data.category = "Clinical",
        data.type = "Clinical Supplement",
        data.format = "BCR Biotab"
      )
    },
    error = function(e) {
      message("GDCquery failed for ", cancer, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(query.clinical)) {
    return(NULL)
  }
  
  clinical_manifest <- tryCatch(
    {
      getResults(query.clinical)
    },
    error = function(e) {
      message("Could not get clinical manifest for ", cancer, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(clinical_manifest) || nrow(clinical_manifest) == 0) {
    message("No BCR Biotab clinical files found for ", cancer)
    return(NULL)
  }
  
  message("Number of BCR Biotab files found: ", nrow(clinical_manifest))
  
  download_status <- tryCatch(
    {
      GDCdownload(
        query.clinical,
        method = "api",
        directory = clinical_download_dir
      )
      TRUE
    },
    error = function(e) {
      message("GDCdownload failed for ", cancer, ": ", e$message)
      FALSE
    }
  )
  
  if (!download_status) {
    return(NULL)
  }
  
  clinical_tab_all <- tryCatch(
    {
      GDCprepare(
        query = query.clinical,
        directory = clinical_download_dir
      )
    },
    error = function(e) {
      message("GDCprepare failed for ", cancer, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(clinical_tab_all)) {
    return(NULL)
  }
  
  patient_table_name <- grep(
    "clinical_patient",
    names(clinical_tab_all),
    value = TRUE
  )[1]
  
  if (is.na(patient_table_name)) {
    message("No clinical_patient table found for ", cancer)
    return(NULL)
  }
  
  clinical_patient <- clinical_tab_all[[patient_table_name]]
  
  if (!"bcr_patient_barcode" %in% colnames(clinical_patient)) {
    message("bcr_patient_barcode column not found for ", cancer)
    return(NULL)
  }
  
  selected_stage_col <- choose_best_stage_column(
    clinical_patient = clinical_patient,
    cancer = cancer
  )
  
  if (is.null(selected_stage_col)) {
    return(NULL)
  }
  
  clinical_stage_table <- clinical_patient %>%
    filter(grepl("^TCGA-", bcr_patient_barcode)) %>%
    transmute(
      Cancer = cancer,
      bcr_patient_barcode = toupper(bcr_patient_barcode),
      stage_source_column = selected_stage_col,
      stage_raw = .data[[selected_stage_col]],
      stage = clean_clinical_stage(.data[[selected_stage_col]])
    ) %>%
    filter(!is.na(stage)) %>%
    distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
    mutate(
      stage = factor(stage, levels = c("I", "II", "III", "IV"))
    )
  
  if (nrow(clinical_stage_table) == 0) {
    message("No valid Stage I-IV patients found for ", cancer)
    return(NULL)
  }
  
  message("Final clinical stage distribution for ", cancer, ":")
  print(table(clinical_stage_table$stage, useNA = "ifany"))
  
  return(clinical_stage_table)
}

# -----------------------------
# 8. Read expression and merge with clinical stage
# -----------------------------

create_stage_expression_table <- function(cancer, clinical_stage_table) {
  
  message("\nExpression data: ", cancer)
  
  normalized_file <- file.path(
    processed_dir,
    paste0("normalized_counts_log2_", cancer, ".csv")
  )
  
  normalized_counts <- read.csv(
    normalized_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!"Gene" %in% colnames(normalized_counts)) {
    stop(
      "The normalized count file for ",
      cancer,
      " must contain a column named 'Gene'."
    )
  }
  
  sample_cols <- grep(
    "^TCGA[.-]|^XTCGA[.-]",
    colnames(normalized_counts),
    value = TRUE
  )
  
  if (length(sample_cols) == 0) {
    message("No TCGA sample columns found in ", normalized_file)
    return(NULL)
  }
  
  counts_long <- normalized_counts %>%
    filter(Gene %in% selected_genes) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "NormalizedCount"
    ) %>%
    mutate(
      Sample_clean = sub("^X(?=TCGA)", "", Sample, perl = TRUE),
      Sample_clean = gsub("\\.", "-", Sample_clean),
      bcr_patient_barcode = substr(Sample_clean, 1, 12),
      sample_type_code = ifelse(
        nchar(Sample_clean) >= 15,
        substr(Sample_clean, 14, 15),
        NA
      ),
      NormalizedCount = as.numeric(NormalizedCount)
    ) %>%
    filter(is.na(sample_type_code) | sample_type_code == "01")
  
  if (nrow(counts_long) == 0) {
    message("No selected gene expression rows found for ", cancer)
    return(NULL)
  }
  
  nested_stage <- counts_long %>%
    inner_join(
      clinical_stage_table,
      by = "bcr_patient_barcode"
    ) %>%
    group_by(Cancer, bcr_patient_barcode, Gene, stage, stage_source_column) %>%
    summarise(
      NormalizedCount = mean(NormalizedCount, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(
      Cancer,
      bcr_patient_barcode,
      Gene,
      stage,
      NormalizedCount,
      stage_source_column
    )
  
  if (nrow(nested_stage) == 0) {
    message("No matching patients between clinical stage and expression for ", cancer)
    return(NULL)
  }
  
  message("Final matched patients per stage for ", cancer, ":")
  print(
    nested_stage %>%
      distinct(bcr_patient_barcode, stage) %>%
      count(stage)
  )
  
  return(nested_stage)
}

# -----------------------------
# 9. Loop over all cancer types
# -----------------------------

all_stage_expression <- list()
analysis_log <- list()

for (cancer in cancers) {
  
  clinical_stage_table <- get_clinical_stage_table(cancer)
  
  if (is.null(clinical_stage_table)) {
    
    analysis_log[[cancer]] <- tibble(
      Cancer = cancer,
      Status = "Skipped: no valid stage table",
      Stage_source_column = NA_character_,
      Clinical_patients = NA_integer_,
      Matched_patients = NA_integer_
    )
    
    next
  }
  
  stage_expression <- create_stage_expression_table(
    cancer = cancer,
    clinical_stage_table = clinical_stage_table
  )
  
  if (is.null(stage_expression)) {
    
    analysis_log[[cancer]] <- tibble(
      Cancer = cancer,
      Status = "Skipped: no matching expression data",
      Stage_source_column = unique(clinical_stage_table$stage_source_column)[1],
      Clinical_patients = n_distinct(clinical_stage_table$bcr_patient_barcode),
      Matched_patients = NA_integer_
    )
    
    next
  }
  
  all_stage_expression[[cancer]] <- stage_expression
  
  analysis_log[[cancer]] <- tibble(
    Cancer = cancer,
    Status = "Included",
    Stage_source_column = unique(stage_expression$stage_source_column)[1],
    Clinical_patients = n_distinct(clinical_stage_table$bcr_patient_barcode),
    Matched_patients = n_distinct(stage_expression$bcr_patient_barcode)
  )
}

stage_expression_all <- bind_rows(all_stage_expression)
stage_analysis_log <- bind_rows(analysis_log)

if (nrow(stage_expression_all) == 0) {
  stop("No stage-expression data were generated. Check clinical and expression input files.")
}

# -----------------------------
# 10. Save combined expression table
# -----------------------------

stage_expression_all <- stage_expression_all %>%
  mutate(
    Cancer = factor(Cancer, levels = plot_cancer_order),
    Gene = factor(Gene, levels = selected_genes),
    stage = factor(stage, levels = c("I", "II", "III", "IV"))
  ) %>%
  arrange(Cancer, Gene, stage, bcr_patient_barcode)

write.csv(
  stage_expression_all,
  file = file.path(table_dir, "stage_expression_all_cancers.csv"),
  row.names = FALSE
)

write.csv(
  stage_analysis_log,
  file = file.path(table_dir, "stage_analysis_log.csv"),
  row.names = FALSE
)

cat("\nStage analysis log:\n")
print(stage_analysis_log)

# -----------------------------
# 11. Save patient counts
# -----------------------------

stage_patient_counts <- stage_expression_all %>%
  distinct(Cancer, bcr_patient_barcode, stage) %>%
  count(Cancer, stage, name = "n_patients") %>%
  complete(
    Cancer,
    stage = factor(c("I", "II", "III", "IV"), levels = c("I", "II", "III", "IV")),
    fill = list(n_patients = 0)
  ) %>%
  arrange(Cancer, stage)

write.csv(
  stage_patient_counts,
  file = file.path(table_dir, "stage_patient_counts_all_cancers.csv"),
  row.names = FALSE
)

cat("\nPatient counts by cancer and stage:\n")
print(stage_patient_counts)

# -----------------------------
# 12. ANOVA for each Cancer + Gene
# -----------------------------

get_anova_p <- function(df) {
  
  df <- df %>%
    filter(!is.na(NormalizedCount), !is.na(stage))
  
  if (length(unique(df$stage)) < 2) {
    return(NA_real_)
  }
  
  p_value <- tryCatch(
    {
      fit <- aov(NormalizedCount ~ stage, data = df)
      summary(fit)[[1]][["Pr(>F)"]][1]
    },
    error = function(e) NA_real_
  )
  
  return(p_value)
}

anova_results <- stage_expression_all %>%
  group_by(Cancer, Gene) %>%
  group_modify(~ tibble(p_value = get_anova_p(.x))) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    p_signif = p_to_star(p_value)
  )

write.csv(
  anova_results,
  file = file.path(table_dir, "stage_ANOVA_all_cancers.csv"),
  row.names = FALSE
)

cat("\nANOVA results saved:\n")
cat(file.path(table_dir, "stage_ANOVA_all_cancers.csv"), "\n")

# -----------------------------
# 13. Pan-cancer stage-based expression boxplot
# -----------------------------

plot_data <- stage_expression_all %>%
  filter(Gene %in% plot_genes) %>%
  mutate(
    Cancer = factor(as.character(Cancer), levels = plot_cancer_order),
    Gene = factor(as.character(Gene), levels = plot_genes),
    stage = factor(stage, levels = c("I", "II", "III", "IV"))
  ) %>%
  filter(!is.na(Cancer))

if (nrow(plot_data) == 0) {
  stop("plot_data is empty. Check plot_genes or cancer order.")
}

y_positions <- plot_data %>%
  group_by(Cancer, Gene) %>%
  summarise(
    y_pos = max(NormalizedCount, na.rm = TRUE) + 0.35,
    .groups = "drop"
  )

anova_labels <- anova_results %>%
  filter(Gene %in% plot_genes) %>%
  left_join(y_positions, by = c("Cancer", "Gene")) %>%
  filter(p_signif != "") %>%
  mutate(
    Cancer = factor(as.character(Cancer), levels = plot_cancer_order),
    Gene = factor(as.character(Gene), levels = plot_genes)
  ) %>%
  filter(!is.na(Cancer))

p <- ggplot(
  plot_data,
  aes(x = Gene, y = NormalizedCount, fill = stage)
) +
  geom_boxplot(
    outlier.size = 0.35,
    linewidth = 0.25
  ) +
  geom_text(
    data = anova_labels,
    aes(x = Gene, y = y_pos, label = p_signif),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ Cancer,
    ncol = 3,
    scales = "free_y"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(
      size = 16,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      hjust = 1,
      colour = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 16,
      colour = "black",
      face = "bold"
    ),
    axis.title.x = element_blank(),
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8
    ),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    legend.text = element_text(
      size = 14,
      face = "bold"
    )
  ) +
  scale_fill_manual(
    values = c(
      "I" = "#20639B",
      "II" = "#ED553B",
      "III" = "#F6D55C",
      "IV" = "#3CAEA3"
    ),
    drop = FALSE
  ) +
  labs(
    y = "log2 normalized count",
    fill = "Stage"
  )

print(p)

# -----------------------------
# 14. Save final figure
# -----------------------------

ggsave(
  filename = file.path(figure_dir, "stage_expression_pan_cancer.png"),
  plot = p,
  width = 12,
  height = 15,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "stage_expression_pan_cancer.pdf"),
  plot = p,
  width = 12,
  height = 15
)

ggsave(
  filename = file.path(figure_dir, "stage_expression_pan_cancer.tiff"),
  plot = p,
  width = 12,
  height = 15,
  dpi = 600,
  compression = "lzw"
)

cat("\nDone. Stage-based pan-cancer analysis completed successfully.\n")

cat("\nMain outputs:\n")
cat(file.path(table_dir, "stage_expression_all_cancers.csv"), "\n")
cat(file.path(table_dir, "stage_ANOVA_all_cancers.csv"), "\n")
cat(file.path(table_dir, "stage_patient_counts_all_cancers.csv"), "\n")
cat(file.path(table_dir, "stage_analysis_log.csv"), "\n")
cat(file.path(figure_dir, "stage_expression_pan_cancer.png"), "\n")
cat(file.path(figure_dir, "stage_expression_pan_cancer.pdf"), "\n")
cat(file.path(figure_dir, "stage_expression_pan_cancer.tiff"), "\n")
