# ==============================================================================
# AGO2-MYC co-amplification across 13 TCGA cancer types
#
# Main question:
# When AGO2 is amplified (+2), is MYC also amplified (+2)
# in the same tumor?
#
# Discrete CNA interpretation:
# -2 = Deep deletion
# -1 = Shallow loss
#  0 = Diploid
# +1 = Gain
# +2 = Amplification
#
# Source:
# cBioPortal - TCGA Firehose Legacy
# ==============================================================================


# ==============================================================================
# 1. Load packages
# ==============================================================================

library(cBioPortalData)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(ggplot2)


# ==============================================================================
# 2. Connect to cBioPortal
# ==============================================================================

cbio <- cBioPortal()


# ==============================================================================
# 3. Define the 13 TCGA Firehose Legacy studies
# ==============================================================================

study_lookup <- tibble::tribble(
  
  ~Cancer, ~studyId,
  
  "BLCA", "blca_tcga",
  "BRCA", "brca_tcga",
  "COAD", "coadread_tcga",
  "ESCA", "esca_tcga",
  "GBM",  "gbm_tcga",
  "KIRC", "kirc_tcga",
  "LIHC", "lihc_tcga",
  "LUAD", "luad_tcga",
  "LUSC", "lusc_tcga",
  "PRAD", "prad_tcga",
  "STAD", "stad_tcga",
  "THCA", "thca_tcga",
  "UCEC", "ucec_tcga"
)


print(study_lookup)

# ==============================================================================
# 4. Define genes
# ==============================================================================

# NCBI Entrez Gene IDs
AGO2_ID <- 27161
MYC_ID  <- 4609

gene_ids <- c(
  AGO2_ID,
  MYC_ID
)

# ==============================================================================
# 5. Function to retrieve discrete CNA data for one cancer
# ==============================================================================

get_cna_for_one_cancer <- function(cancer_code, study_id) {
  
  
  message(
    "\n==================================================\n",
    "Processing: ", cancer_code,
    " | ", study_id,
    "\n=================================================="
  )
  
  
  # --------------------------------------------------------------------------
  # A. Find discrete CNA molecular profile
  # --------------------------------------------------------------------------
  
  profiles <- molecularProfiles(
    cbio,
    study_id
  )
  
  
  discrete_profiles <- profiles %>%
    dplyr::filter(
      molecularAlterationType == "COPY_NUMBER_ALTERATION",
      datatype == "DISCRETE"
    )
  
  
  # Prefer a GISTIC profile
  gistic_profiles <- discrete_profiles %>%
    dplyr::filter(
      grepl(
        "gistic",
        molecularProfileId,
        ignore.case = TRUE
      )
    )
  
  
  if (nrow(gistic_profiles) == 1) {
    
    profile_id <- gistic_profiles$molecularProfileId[1]
    
  } else if (nrow(discrete_profiles) == 1) {
    
    profile_id <- discrete_profiles$molecularProfileId[1]
    
  } else {
    
    print(discrete_profiles)
    
    stop(
      paste0(
        "Could not uniquely identify the discrete CNA profile for ",
        cancer_code
      )
    )
  }
  
  
  message(
    "CNA profile: ",
    profile_id
  )
  
  
  # --------------------------------------------------------------------------
  # B. Identify samples that have CNA data
  # --------------------------------------------------------------------------
  
  sample_lists <- sampleLists(
    cbio,
    study_id
  )
  
  
  cna_sample_list <- sample_lists %>%
    dplyr::filter(
      category == "all_cases_with_cna_data"
    )
  
  
  if (nrow(cna_sample_list) != 1) {
    
    print(cna_sample_list)
    
    stop(
      paste0(
        "Could not identify CNA sample list for ",
        cancer_code
      )
    )
  }
  
  
  sample_list_id <- cna_sample_list$sampleListId[1]
  
  
  sample_ids <- samplesInSampleLists(
    cbio,
    sample_list_id
  )[[1]]
  
  
  message(
    "Samples with CNA data: ",
    length(sample_ids)
  )
  
  
  # --------------------------------------------------------------------------
  # C. Retrieve AGO2 and MYC discrete CNA calls
  # --------------------------------------------------------------------------
  
  cna_list <- molecularData(
    
    api = cbio,
    
    molecularProfileIds = profile_id,
    
    entrezGeneIds = gene_ids,
    
    sampleIds = sample_ids
  )
  
  
  cna_raw <- cna_list[[1]]
  
  
  # --------------------------------------------------------------------------
  # D. Assign gene symbols
  # --------------------------------------------------------------------------
  
  cna_long <- cna_raw %>%
    
    dplyr::mutate(
      
      Gene = dplyr::case_when(
        
        entrezGeneId == AGO2_ID ~ "AGO2",
        
        entrezGeneId == MYC_ID ~ "MYC",
        
        TRUE ~ NA_character_
      ),
      
      CNA = as.integer(value)
    ) %>%
    
    dplyr::filter(
      !is.na(Gene)
    ) %>%
    
    dplyr::select(
      sampleId,
      Gene,
      CNA
    )
  
  
  # --------------------------------------------------------------------------
  # E. Check duplicate sample-gene rows
  # --------------------------------------------------------------------------
  
  duplicates <- cna_long %>%
    
    dplyr::count(
      sampleId,
      Gene,
      name = "n"
    ) %>%
    
    dplyr::filter(
      n > 1
    )
  
  
  if (nrow(duplicates) > 0) {
    
    print(duplicates)
    
    stop(
      paste0(
        "Duplicate CNA calls detected in ",
        cancer_code
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # F. Convert to one row per tumor
  # --------------------------------------------------------------------------
  
  cna_wide <- cna_long %>%
    
    tidyr::pivot_wider(
      names_from = Gene,
      values_from = CNA
    ) %>%
    
    dplyr::filter(
      !is.na(AGO2),
      !is.na(MYC)
    ) %>%
    
    dplyr::mutate(
      Cancer = cancer_code
    ) %>%
    
    dplyr::select(
      Cancer,
      sampleId,
      AGO2,
      MYC
    )
  
  
  message(
    "Tumors with both AGO2 and MYC data: ",
    nrow(cna_wide)
  )
  
  
  return(cna_wide)
}

# ==============================================================================
# 6. Retrieve AGO2 and MYC CNA data for all 13 cancers
# ==============================================================================

all_cna_data <- purrr::map2_dfr(
  
  study_lookup$Cancer,
  study_lookup$studyId,
  
  ~ get_cna_for_one_cancer(
    cancer_code = .x,
    study_id = .y
  )
)

sample_counts <- all_cna_data %>%
  
  dplyr::group_by(Cancer) %>%
  
  dplyr::summarise(
    Total_samples = dplyr::n(),
    .groups = "drop"
  )

print(sample_counts)

# ==============================================================================
# 9. Calculate AGO2-MYC co-amplification
# ==============================================================================

amplification_summary <- all_cna_data %>%
  
  dplyr::group_by(Cancer) %>%
  
  dplyr::summarise(
    
    Total_samples = dplyr::n(),
    
    AGO2_amplified = sum(
      AGO2 == 2
    ),
    
    MYC_amplified = sum(
      MYC == 2
    ),
    
    Both_amplified = sum(
      AGO2 == 2 &
        MYC == 2
    ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    # % of all tumors with AGO2 amplification
    AGO2_amp_percent =
      100 *
      AGO2_amplified /
      Total_samples,
    
    
    # % of all tumors with MYC amplification
    MYC_amp_percent =
      100 *
      MYC_amplified /
      Total_samples,
    
    
    # % of all tumors with both amplifications
    Coamplification_percent =
      100 *
      Both_amplified /
      Total_samples,
    
    
    # Main metric:
    # Among AGO2-amplified tumors,
    # what % are also MYC-amplified?
    
    Percent_AGO2_amp_with_MYC_amp =
      dplyr::if_else(
        
        AGO2_amplified > 0,
        
        100 *
          Both_amplified /
          AGO2_amplified,
        
        NA_real_
      ),
    
    
    # Reverse question:
    # Among MYC-amplified tumors,
    # what % are also AGO2-amplified?
    
    Percent_MYC_amp_with_AGO2_amp =
      dplyr::if_else(
        
        MYC_amplified > 0,
        
        100 *
          Both_amplified /
          MYC_amplified,
        
        NA_real_
      )
  )

# ==============================================================================
# 10. Fisher's exact test for each cancer
# ==============================================================================

run_fisher <- function(df) {
  
  
  # Force TRUE/FALSE levels so table is always 2 x 2
  AGO2_status <- factor(
    df$AGO2 == 2,
    levels = c(FALSE, TRUE)
  )
  
  
  MYC_status <- factor(
    df$MYC == 2,
    levels = c(FALSE, TRUE)
  )
  
  
  tab <- table(
    AGO2_status,
    MYC_status
  )
  
  
  # If one gene has no variation in amplification status,
  # Fisher's exact test is not informative
  
  if (
    sum(tab[2, ]) == 0 ||
    sum(tab[, 2]) == 0
  ) {
    
    return(
      tibble::tibble(
        Odds_ratio = NA_real_,
        P_value = NA_real_
      )
    )
  }
  
  
  test <- fisher.test(tab)
  
  
  tibble::tibble(
    
    Odds_ratio =
      unname(test$estimate),
    
    P_value =
      test$p.value
  )
}
fisher_results <- all_cna_data %>%
  
  split(.$Cancer) %>%
  
  purrr::imap_dfr(
    
    ~ run_fisher(.x) %>%
      
      dplyr::mutate(
        Cancer = .y
      )
  )
fisher_results <- fisher_results %>%
  
  dplyr::mutate(
    
    FDR = p.adjust(
      P_value,
      method = "BH"
    )
  )

# ==============================================================================
# 11. Final result table
# ==============================================================================

final_results <- amplification_summary %>%
  
  dplyr::left_join(
    fisher_results,
    by = "Cancer"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Percent_AGO2_amp_with_MYC_amp
    )
  )


print(
  final_results,
  n = 13
)

# ==============================================================================
# 12. Save tables
# ==============================================================================

output_dir <- file.path(
  getwd(),
  "07_AGO2_MYC_CNV"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


readr::write_csv(
  
  all_cna_data,
  
  file.path(
    output_dir,
    "AGO2_MYC_sample_level_discrete_CNA.csv"
  )
)


readr::write_csv(
  
  final_results,
  
  file.path(
    output_dir,
    "AGO2_MYC_coamplification_results.csv"
  )
)

# ==============================================================================
# 13. Prepare data for visualization
# ==============================================================================

plot_data <- final_results %>%
  
  dplyr::filter(
    AGO2_amplified > 0
  ) %>%
  
  dplyr::arrange(
    Percent_AGO2_amp_with_MYC_amp
  ) %>%
  
  dplyr::mutate(
    
    Cancer = factor(
      Cancer,
      levels = Cancer
    ),
    
    n_label = paste0(
      Both_amplified,
      "/",
      AGO2_amplified
    )
  )
# ==============================================================================
# 14. Co-amplification bar plot
# ==============================================================================

p_coamp <- ggplot(
  plot_data,
  aes(
    x = Cancer,
    y = Percent_AGO2_amp_with_MYC_amp
  )
) +
  
  geom_col(
    width = 0.68,
    fill = "#4C6272"
  ) +
  
  # percentage outside bars
  geom_text(
    aes(
      label = paste0(
        round(Percent_AGO2_amp_with_MYC_amp, 1),
        "%"
      )
    ),
    hjust = -0.15,
    size = 3.8,
    color = "black"
  ) +
  
  # counts inside bars
  geom_text(
    aes(
      y = 3,
      label = n_label
    ),
    hjust = 0,
    size = 3.2,
    color = "white",
    fontface = "bold"
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(0, 108),
    breaks = seq(0, 100, 20)
  ) +
  
  labs(
    x = NULL,
    y = "MYC-amplified tumors among AGO2-amplified tumors (%)"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_coamp


output_dir_fig<- file.path( 
  getwd(), 
  "Figures",
  "08_AGO2_MYC_CNV"
  ) 

dir.create( 
  output_dir_fig, 
  recursive = TRUE,
  showWarnings = FALSE
  )


ggsave(
  filename = file.path(
    output_dir_fig,
    "AGO2_MYC_coamplification_across_TCGA.png"
  ),
  plot = p_coamp,
  width = 8,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir_fig,
    "AGO2_MYC_coamplification_across_TCGA.tiff"
  ),
  plot = p_coamp,
  width = 8,
  height = 6.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
