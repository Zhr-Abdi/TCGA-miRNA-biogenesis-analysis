# ==============================================================================
# AGO2 neighborhood co-amplification analysis across 13 TCGA cancers
# ==============================================================================


# ==============================================================================
# 1. Load packages
# ==============================================================================

library(cBioPortalData)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)


# ==============================================================================
# 2. Connect to cBioPortal
# ==============================================================================

cbio <- cBioPortal()


# ==============================================================================
# 3. TCGA Firehose Legacy studies
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



# ==============================================================================
# 4. Define AGO2 neighborhood genes
# ==============================================================================

# Entrez IDs from hg19 annotation

neighbor_genes <- c(
  
  AGO2     = 27161,
  CHRAC1   = 54108,
  PTK2     = 5747,
  TRAPPC9  = 83696,
  DENND3   = 22898,
  SLC45A4  = 57210,
  GPR20    = 2843,
  PTP4A3   = 11156,
  MROH5    = 389690,
  KCNK9    = 51305
  
)

gene_ids <- unname(neighbor_genes)



# ==============================================================================
# 5. Function to retrieve CNA for one cancer
# ==============================================================================

get_cna_for_one_cancer <- function(cancer_code, study_id) {
  
  
  message("Processing: ", cancer_code)
  
  
  
  # -----------------------------
  # Find CNA profile
  # -----------------------------
  
  profiles <- molecularProfiles(
    cbio,
    study_id
  )
  
  
  discrete_profiles <- profiles %>%
    filter(
      molecularAlterationType == "COPY_NUMBER_ALTERATION",
      datatype == "DISCRETE"
    )
  
  
  gistic_profiles <- discrete_profiles %>%
    filter(
      grepl(
        "gistic",
        molecularProfileId,
        ignore.case = TRUE
      )
    )
  
  
  if(nrow(gistic_profiles) == 1){
    
    profile_id <- gistic_profiles$molecularProfileId[1]
    
  } else {
    
    profile_id <- discrete_profiles$molecularProfileId[1]
    
  }
  
  
  
  # -----------------------------
  # Samples with CNA data
  # -----------------------------
  
  sample_lists <- sampleLists(
    cbio,
    study_id
  )
  
  
  sample_list_id <- sample_lists %>%
    filter(
      category == "all_cases_with_cna_data"
    ) %>%
    pull(sampleListId)
  
  
  
  sample_ids <- samplesInSampleLists(
    cbio,
    sample_list_id
  )[[1]]
  
  
  
  # -----------------------------
  # Retrieve CNA
  # -----------------------------
  
  cna_list <- molecularData(
    api = cbio,
    molecularProfileIds = profile_id,
    entrezGeneIds = gene_ids,
    sampleIds = sample_ids
  )
  
  
  cna_raw <- cna_list[[1]]
  
  
  
  # -----------------------------
  # Convert Entrez ID to gene name
  # -----------------------------
  
  cna_long <- cna_raw %>%
    
    mutate(
      
      Gene = names(neighbor_genes)[
        match(entrezGeneId, neighbor_genes)
      ],
      
      CNA = as.integer(value)
    ) %>%
    
    dplyr::select(
      sampleId,
      Gene,
      CNA
    )
  
  
  
  # -----------------------------
  # Make wide format
  # -----------------------------
  
  cna_wide <- cna_long %>%
    
    pivot_wider(
      names_from = Gene,
      values_from = CNA
    ) %>%
    
    filter(
      !is.na(AGO2)
    ) %>%
    
    mutate(
      Cancer = cancer_code
    ) %>%
    
    dplyr::select(
      Cancer,
      sampleId,
      everything()
    )
  
  
  return(cna_wide)
  
}



# ==============================================================================
# 6. Retrieve CNA for all cancers
# ==============================================================================


all_cna_data <- map2_dfr(
  
  study_lookup$Cancer,
  study_lookup$studyId,
  
  ~get_cna_for_one_cancer(
    cancer_code = .x,
    study_id = .y
  )
)



# Check
head(all_cna_data)



# ==============================================================================
# 7. Calculate co-amplification
# ==============================================================================


partner_genes <- names(neighbor_genes)[
  names(neighbor_genes) != "AGO2"
]



coamp_results <- map_dfr(
  
  partner_genes,
  
  function(gene){
    
    all_cna_data %>%
      
      group_by(Cancer) %>%
      
      summarise(
        
        AGO2_amplified =
          sum(AGO2 == 2, na.rm = TRUE),
        
        
        Partner_amplified =
          sum(.data[[gene]] == 2, na.rm = TRUE),
        
        
        Both_amplified =
          sum(
            AGO2 == 2 &
              .data[[gene]] == 2,
            na.rm = TRUE
          ),
        
        
        Percent_AGO2_amp_with_partner_amp =
          
          ifelse(
            AGO2_amplified > 0,
            100 *
              Both_amplified /
              AGO2_amplified,
            NA
          ),
        
        
        .groups = "drop"
        
      ) %>%
      
      mutate(
        Partner_gene = gene
      )
    
  }
)



# ==============================================================================
# 8. Final table
# ==============================================================================

coamp_results <- coamp_results %>%
  
  dplyr::select(
    Cancer,
    Partner_gene,
    AGO2_amplified,
    Partner_amplified,
    Both_amplified,
    Percent_AGO2_amp_with_partner_amp
  )



print(coamp_results)

# Create output directory
output_dir <- file.path(
  getwd(),
  "07_AGO2_CNV"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# Save final co-amplification table

write.csv(
  coamp_results,
  file.path(
    output_dir,
    "AGO2_neighbor_coamplification_TCGA13.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_cna_data,
  file.path(
    output_dir,
    "AGO2_neighbor_CNA_sample_level_TCGA13.csv"
  ),
  row.names = FALSE
)


run_fisher_coamp <- function(df, partner_gene){
  
  
  # create 2x2 table
  
  tab <- table(
    AGO2_amp = df$AGO2 == 2,
    Partner_amp = df[[partner_gene]] == 2
  )
  
  
  # Fisher test
  
  test <- fisher.test(tab)
  
  
  data.frame(
    Partner_gene = partner_gene,
    Odds_ratio = unname(test$estimate),
    P_value = test$p.value
  )
  
}
partner_genes <- c(
  "CHRAC1",
  "PTK2",
  "TRAPPC9",
  "DENND3",
  "SLC45A4",
  "GPR20",
  "PTP4A3",
  "MROH5",
  "KCNK9"
)


run_fisher_coamp <- function(df, partner_gene){
  
  
  tab <- table(
    AGO2_amp = df$AGO2 == 2,
    Partner_amp = df[[partner_gene]] == 2
  )
  
  

  if(
    nrow(tab) < 2 ||
    ncol(tab) < 2
  ){
    
    return(
      data.frame(
        Partner_gene = partner_gene,
        Odds_ratio = NA,
        P_value = NA
      )
    )
  }
  
  
  test <- fisher.test(tab)
  
  
  data.frame(
    Partner_gene = partner_gene,
    Odds_ratio = unname(test$estimate),
    P_value = test$p.value
  )
}

fisher_results <- purrr::map_dfr(
  
  partner_genes,
  
  function(gene){
    
    all_cna_data %>%
      split(.$Cancer) %>%
      purrr::imap_dfr(
        
        ~run_fisher_coamp(
          df = .x,
          partner_gene = gene
        ) %>%
          mutate(
            Cancer = .y
          )
      )
  }
)


fisher_results <- fisher_results %>%
  mutate(
    FDR = p.adjust(
      P_value,
      method = "BH"
    )
  )


print(fisher_results)

library(dplyr)
library(ggplot2)


gene_order <- c(
  "CHRAC1",
  "PTK2",
  "TRAPPC9",
  "DENND3",
  "SLC45A4",
  "GPR20",
  "PTP4A3",
  "MROH5",
  "KCNK9"
)


cancer_order <- c(
  "BLCA", "BRCA", "COAD", "ESCA", "GBM", "KIRC",
  "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC"
)

heatmap_df <- coamp_results %>%
  mutate(
    Partner_gene = factor(Partner_gene, levels = gene_order),
    Cancer = factor(Cancer, levels = cancer_order)
  )

ggplot(heatmap_df, aes(x = Cancer, y = Partner_gene,
                       fill = Percent_AGO2_amp_with_partner_amp)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Percent_AGO2_amp_with_partner_amp, 1)), size = 3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  labs(
    x = "Cancer type",
    y = "Gene",
    fill = "% co-amplification",
) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      face = "italic",
      color = "black"
    ),
    panel.grid = element_blank()
  )
# Create output directory
output_dir <- file.path(
  getwd(),
  "Figures",
  "08_AGO2_MYC_CNV"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Save PNG
ggsave(
  filename = file.path(
    output_dir,
    "AGO2_neighbor_coamplification_heatmap.png"
  ),
  plot = last_plot(),
  width = 8,
  height = 5,
  dpi = 300,
  units = "in"
)


# Save TIFF (publication quality)
ggsave(
  filename = file.path(
    output_dir,
    "AGO2_neighbor_coamplification_heatmap.tiff"
  ),
  plot = last_plot(),
  width = 8,
  height = 5,
  dpi = 600,
  units = "in",
  compression = "lzw"
)

