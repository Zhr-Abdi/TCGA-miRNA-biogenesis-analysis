

comparison_table <- gene_results %>%
  
  tidyr::pivot_wider(
    names_from = Gene,
    values_from = c(
      log2FoldChange,
      padj
    )
  ) %>%
  
  dplyr::rename(
    TARBP2_log2FC = log2FoldChange_TARBP2,
    DICER1_log2FC = log2FoldChange_DICER1,
    TARBP2_padj = padj_TARBP2,
    DICER1_padj = padj_DICER1
  ) %>%
  
  dplyr::mutate(
    Divergence_score =
      TARBP2_log2FC - DICER1_log2FC
  )

print(comparison_table)
comparison_table <- comparison_table %>%
  
  dplyr::mutate(
    
    TARBP2_status = dplyr::case_when(
      TARBP2_padj < 0.05 & TARBP2_log2FC > 0 ~ "Up",
      TARBP2_padj < 0.05 & TARBP2_log2FC < 0 ~ "Down",
      TRUE ~ "Unchanged"
    ),
    
    DICER1_status = dplyr::case_when(
      DICER1_padj < 0.05 & DICER1_log2FC > 0 ~ "Up",
      DICER1_padj < 0.05 & DICER1_log2FC < 0 ~ "Down",
      TRUE ~ "Unchanged"
    ),
    
    Regulation_pattern = paste0(
      "TARBP2 ", TARBP2_status,
      " / DICER1 ", DICER1_status
    )
  )

print(
  comparison_table %>%
    dplyr::select(
      Cancer,
      TARBP2_log2FC,
      TARBP2_padj,
      DICER1_log2FC,
      DICER1_padj,
      Regulation_pattern
    )
)

wilcox.test(
  comparison_table$TARBP2_log2FC,
  comparison_table$DICER1_log2FC,
  paired = TRUE,
  alternative = "greater",
  exact = FALSE
)
median(
  comparison_table$TARBP2_log2FC
)

median(
  comparison_table$DICER1_log2FC
)

comparison_table <- comparison_table %>%
  dplyr::mutate(
    Delta_log2FC = TARBP2_log2FC - DICER1_log2FC
  )

comparison_table %>%
  dplyr::select(
    Cancer,
    TARBP2_log2FC,
    DICER1_log2FC,
    Delta_log2FC
  )
median(
  comparison_table$Delta_log2FC
)
plot_data <- comparison_table %>%
  
  dplyr::select(
    Cancer,
    TARBP2_log2FC,
    DICER1_log2FC
  ) %>%
  
  tidyr::pivot_longer(
    cols = c(
      TARBP2_log2FC,
      DICER1_log2FC
    ),
    names_to = "Gene",
    values_to = "log2FC"
  ) %>%
  
  dplyr::mutate(
    Gene = dplyr::recode(
      Gene,
      TARBP2_log2FC = "TARBP2",
      DICER1_log2FC = "DICER1"
    )
  )
p <- ggplot(
  plot_data,
  aes(
    x = log2FC,
    y = Cancer,
    group = Cancer
  )
) +
  
  # Connecting line between DICER1 and TARBP2
  geom_line(
    color = "grey65",
    linewidth = 0.8
  ) +
  
  # Points
  geom_point(
    aes(
      shape = Gene,
      color = Gene
    ),
    size = 4
  ) +
  
  # Zero line
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.7
  ) +
  
  # Two colors
  scale_color_manual(
    values = c(
      "DICER1" = "#377EB8",
      "TARBP2" = "#E64B35"
    )
  ) +
  
  # Keep circle and triangle
  scale_shape_manual(
    values = c(
      "DICER1" = 16,
      "TARBP2" = 17
    )
  ) +
  
  labs(
    x = "log2 fold change (Tumor vs Normal)",
    y = "Cancer type",
    color = "Gene",
    shape = "Gene"
  ) +
  
  theme_classic() +
  
  theme(
    
    # Complete rectangular border
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Remove classic-theme axis lines to avoid double borders
    axis.line = element_blank(),
    
    # Title
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    # Axis titles
    axis.title = element_text(
      size = 13,
      face = "bold"
    ),
    
    # Axis labels
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    
    # Legend
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 11
    )
  )

print(p)
# ============================================================
# Save figure for publication and GitHub
# ============================================================
figure_dir <- file.path(
  getwd(),
  "Figures",
  "07_TARBP2_DICER1"
)
# Create directory if it does not exist
dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Base file name
figure_name <- "TARBP2_DICER1_divergent_regulation"



# ------------------------------------------------------------
# TIFF - publication quality
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    paste0(figure_name, ".tiff")
  ),
  plot = p,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# PNG - high resolution / GitHub preview
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    paste0(figure_name, ".png")
  ),
  plot = p,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)


message(
  "Figure saved in: ",
  normalizePath(figure_dir)
)

# ============================================================
# PART 2
# Patient-level correlation between TARBP2 and DICER1
# Tumor vs Normal
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

correlation_results <- list()


for (cancer_code in names(cancer_projects)) {
  
  project_id <- cancer_projects[[cancer_code]]
  
  file_path <- file.path(
    results_dir,
    project_id,
    paste0(
      "NormalizedCount_",
      cancer_code,
      ".csv"
    )
  )
  
  # ----------------------------------------------------------
  # Read normalized expression data
  # ----------------------------------------------------------
  
  expr <- read.csv(
    file_path,
    check.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Keep only TARBP2 and DICER1
  # ----------------------------------------------------------
  
  expr_2genes <- expr %>%
    dplyr::filter(
      Gene %in% c("TARBP2", "DICER1")
    )
  
  
  # ----------------------------------------------------------
  # Transpose:
  # rows = samples
  # columns = genes
  # ----------------------------------------------------------
  
  expr_t <- as.data.frame(
    t(expr_2genes[, -1])
  )
  
  colnames(expr_t) <- expr_2genes$Gene
  
  expr_t$Sample <- rownames(expr_t)
  
  rownames(expr_t) <- NULL
  
  
  # ----------------------------------------------------------
  # Identify sample type
  #
  # 01 = Primary Tumor
  # 11 = Solid Tissue Normal
  # ----------------------------------------------------------
  
  expr_t <- expr_t %>%
    
    dplyr::mutate(
      
      SampleTypeCode = substr(
        Sample,
        14,
        15
      ),
      
      Group = dplyr::case_when(
        SampleTypeCode == "01" ~ "Tumor",
        SampleTypeCode == "11" ~ "Normal",
        TRUE ~ NA_character_
      ),
      
      Patient = substr(
        Sample,
        1,
        12
      ),
      
      Cancer = cancer_code
    ) %>%
    
    dplyr::filter(
      !is.na(Group)
    )
  
  
  # ----------------------------------------------------------
  # Collapse duplicate samples from same patient
  # ----------------------------------------------------------
  
  patient_data <- expr_t %>%
    
    dplyr::group_by(
      Cancer,
      Patient,
      Group
    ) %>%
    
    dplyr::summarise(
      TARBP2 = mean(TARBP2, na.rm = TRUE),
      DICER1 = mean(DICER1, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # Tumor correlation
  # ----------------------------------------------------------
  
  tumor_data <- patient_data %>%
    dplyr::filter(
      Group == "Tumor"
    )
  
  if (nrow(tumor_data) >= 5) {
    
    tumor_test <- cor.test(
      tumor_data$TARBP2,
      tumor_data$DICER1,
      method = "spearman",
      exact = FALSE
    )
    
    tumor_rho <- unname(
      tumor_test$estimate
    )
    
    tumor_p <- tumor_test$p.value
    
  } else {
    
    tumor_rho <- NA
    tumor_p <- NA
  }
  
  
  # ----------------------------------------------------------
  # Normal correlation
  # ----------------------------------------------------------
  
  normal_data <- patient_data %>%
    dplyr::filter(
      Group == "Normal"
    )
  
  if (nrow(normal_data) >= 5) {
    
    normal_test <- cor.test(
      normal_data$TARBP2,
      normal_data$DICER1,
      method = "spearman",
      exact = FALSE
    )
    
    normal_rho <- unname(
      normal_test$estimate
    )
    
    normal_p <- normal_test$p.value
    
  } else {
    
    normal_rho <- NA
    normal_p <- NA
  }
  
  
  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------
  
  correlation_results[[cancer_code]] <- data.frame(
    
    Cancer = cancer_code,
    
    Tumor_n = nrow(tumor_data),
    Tumor_rho = tumor_rho,
    Tumor_p = tumor_p,
    
    Normal_n = nrow(normal_data),
    Normal_rho = normal_rho,
    Normal_p = normal_p,
    
    Delta_rho = tumor_rho - normal_rho
  )
}


# ============================================================
# Combine results
# ============================================================

correlation_table <- dplyr::bind_rows(
  correlation_results
)

print(correlation_table)
