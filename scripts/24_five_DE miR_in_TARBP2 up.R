# ============================================================
# Fig. X3
# Recurrently upregulated mature miRNAs across
# TARBP2-upregulated TCGA cancer types
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


# ============================================================
# 2. Project directory
# ============================================================

project_dir <- getwd()


# ============================================================
# 3. Input directory
# ============================================================

input_dir <- file.path(
  project_dir,
  "08_TARBP2_miRNA_analysis"
)


# ============================================================
# 4. Figure output directory
# ============================================================

figure_output_dir <- file.path(
  project_dir,
  "Figures",
  "10_TARBP2_miR"
)

dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 5. Read miRNA direction matrix
# ============================================================

mirna_matrix <- read_csv(
  file.path(
    input_dir,
    "TARBP2_up_cancers_miRNA_direction_matrix_with_names.csv"
  ),
  show_col_types = FALSE
)


# Check input
dim(mirna_matrix)
head(mirna_matrix)
colnames(mirna_matrix)


# ============================================================
# 6. Define TARBP2-upregulated cancer types
# ============================================================

cancer_cols <- c(
  "COAD",
  "BLCA",
  "BRCA",
  "LUSC",
  "UCEC",
  "GBM"
)


# ============================================================
# 7. Count miRNA recurrence across cancers
#
# N_up          = number of cancers in which the miRNA is
#                 significantly upregulated
#
# N_down        = number of cancers in which the miRNA is
#                 significantly downregulated
#
# N_significant = total number of cancers showing significant
#                 differential expression in either direction
# ============================================================

mirna_recurrence <- mirna_matrix %>%
  
  rowwise() %>%
  
  mutate(
    
    N_up = sum(
      c_across(
        all_of(cancer_cols)
      ) == "Upregulated"
    ),
    
    N_down = sum(
      c_across(
        all_of(cancer_cols)
      ) == "Downregulated"
    ),
    
    N_significant = sum(
      c_across(
        all_of(cancer_cols)
      ) != "Not significant"
    )
    
  ) %>%
  
  ungroup()


# ============================================================
# 8. Check recurrence distribution
# ============================================================

print(
  table(
    mirna_recurrence$N_up
  )
)


# ============================================================
# 9. Select recurrently upregulated candidate miRNAs
#
# Candidate definition:
# Upregulated in at least 5 of 6 TARBP2-up cancers
# and not significantly downregulated in any cancer.
# ============================================================

recurrent_up_mirnas <- mirna_recurrence %>%
  
  filter(
    N_up >= 5,
    N_down == 0
  ) %>%
  
  arrange(
    desc(N_up),
    Mature_miRNA_Name
  )


# Number of recurrent candidates
cat(
  "\nNumber of recurrently upregulated miRNAs:",
  nrow(recurrent_up_mirnas),
  "\n"
)


# View selected candidates
print(
  recurrent_up_mirnas %>%
    select(
      Mature_miRNA_Name,
      Mature_miRNA_ID,
      N_up,
      N_down,
      N_significant
    )
)


# ============================================================
# 10. Add cancers in which each miRNA is upregulated
# ============================================================

recurrent_up_mirnas <- recurrent_up_mirnas %>%
  
  rowwise() %>%
  
  mutate(
    
    Upregulated_in = paste(
      
      cancer_cols[
        c_across(
          all_of(cancer_cols)
        ) == "Upregulated"
      ],
      
      collapse = "; "
    )
    
  ) %>%
  
  ungroup()


# ============================================================
# 11. Export candidates for functional/oncomiR annotation
# ============================================================

annotation_candidates <- recurrent_up_mirnas %>%
  
  select(
    Mature_miRNA_Name,
    Mature_miRNA_ID,
    N_up,
    N_down,
    N_significant,
    Upregulated_in,
    all_of(cancer_cols)
  )


write_csv(
  annotation_candidates,
  file.path(
    figure_output_dir,
    "Fig_X3_recurrent_up_miRNAs_for_oncomiR_annotation.csv"
  )
)


# ============================================================
# 12. Identify miRNAs upregulated in ALL 6 cancers
#
# These miRNAs constitute the core recurrent signature
# used for Fig. X3.
# ============================================================

up_in_all_6 <- recurrent_up_mirnas %>%
  
  filter(
    N_up == 6,
    N_down == 0
  ) %>%
  
  select(
    Mature_miRNA_Name,
    Mature_miRNA_ID,
    N_up,
    N_down,
    Upregulated_in
  )


cat(
  "\nmiRNAs upregulated in all six cancers:\n"
)

print(
  up_in_all_6
)


# ============================================================
# 13. Quality-control check
#
# Expected result from the current analysis = 5 miRNAs
# ============================================================

stopifnot(
  nrow(up_in_all_6) == 5
)


# ============================================================
# 14. Export the five-miRNA core signature
# ============================================================

write_csv(
  up_in_all_6,
  file.path(
    figure_output_dir,
    "Fig_X3_core_5_miRNA_signature.csv"
  )
)


# ============================================================
# 15. Read significant miRNA differential-expression results
# ============================================================

mirna_sig <- read_csv(
  file.path(
    input_dir,
    "TARBP2_up_cancers_significant_miRNAs_with_names.csv"
  ),
  show_col_types = FALSE
)


# ============================================================
# 16. Extract log2FC values for the five core miRNAs
# ============================================================

selected_mirs <- up_in_all_6$Mature_miRNA_Name


plot_df <- mirna_sig %>%
  
  filter(
    Mature_miRNA_Name %in% selected_mirs,
    Cancer %in% cancer_cols
  ) %>%
  
  select(
    Cancer,
    Mature_miRNA_ID,
    Mature_miRNA_Name,
    log2FoldChange,
    padj
  ) %>%
  
  mutate(
    
    # Shorter display labels for the figure
    miRNA = sub(
      "^hsa-",
      "",
      Mature_miRNA_Name
    )
    
  )


# ============================================================
# 17. Quality-control checks for heatmap data
#
# 5 miRNAs x 6 cancers = 30 observations
# ============================================================

stopifnot(
  nrow(plot_df) == 30
)

stopifnot(
  all(
    plot_df$log2FoldChange > 0
  )
)

stopifnot(
  all(
    plot_df$padj < 0.05
  )
)


# ============================================================
# 18. Define plotting order
# ============================================================

plot_df <- plot_df %>%
  
  mutate(
    
    Cancer = factor(
      Cancer,
      levels = cancer_cols
    ),
    
    miRNA = factor(
      miRNA,
      levels = rev(
        sub(
          "^hsa-",
          "",
          selected_mirs
        )
      )
    )
    
  )


# ============================================================
# 19. Create Fig. X3 heatmap
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = Cancer,
    y = miRNA,
    fill = log2FoldChange
  )
) +
  
  # ----------------------------------------------------------
# Heatmap cells
# ----------------------------------------------------------

geom_tile(
  color = "white",
  linewidth = 0.8
) +
  
  
  # ----------------------------------------------------------
# Display log2FC values inside each cell
# ----------------------------------------------------------

geom_text(
  aes(
    label = sprintf(
      "%.2f",
      log2FoldChange
    )
  ),
  size = 4,
  color = "black"
) +
  
  
  # ----------------------------------------------------------
# Sequential color scale
#
# All five selected miRNAs are upregulated.
# Darker red indicates greater positive log2FC.
# ----------------------------------------------------------

scale_fill_gradient(
  low = "#F7E5E2",
  high = "#A71930",
  name = expression(log[2] * "FC")
) +
  
  
  # ----------------------------------------------------------
# Remove spacing around heatmap cells
# ----------------------------------------------------------

scale_x_discrete(
  expand = c(0, 0)
) +
  
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  
  
  # ----------------------------------------------------------
# Figure labels
# ----------------------------------------------------------

labs(
  x = NULL,
  y = NULL,
) +
  
  
  # ----------------------------------------------------------
# Base theme
# ----------------------------------------------------------

theme_classic(
  base_size = 12
) +
  
  
  # ----------------------------------------------------------
# Publication-style formatting
# ----------------------------------------------------------

theme(
  
  # Full border around plotting panel
  panel.border = element_rect(
    color = "black",
    fill = NA,
    linewidth = 1.2
  ),
  
  # Remove default partial axis lines
  axis.line = element_blank(),
  
  # Square plotting panel
  aspect.ratio = 1,
  
  # X-axis labels
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1,
    size = 12,
    color = "black"
  ),
  
  # Y-axis labels
  axis.text.y = element_text(
    size = 12,
    color = "black"
  ),
  
  # Remove tick marks
  axis.ticks = element_blank(),
  
  # Figure title
  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 15,
    margin = margin(
      b = 15
    )
  ),
  
  # Legend
  legend.position = "right",
  
  legend.title = element_text(
    size = 11
  ),
  
  legend.text = element_text(
    size = 10
  ),
  
  # White backgrounds
  panel.background = element_rect(
    fill = "white",
    color = NA
  ),
  
  plot.background = element_rect(
    fill = "white",
    color = NA
  ),
  
  legend.background = element_rect(
    fill = "white",
    color = NA
  ),
  
  # Outer margins
  plot.margin = margin(
    t = 15,
    r = 15,
    b = 15,
    l = 15
  )
  
)


# ============================================================
# 20. Display Fig. X3
# ============================================================

print(p)


# ============================================================
# 21. Figure file name
# ============================================================

figure_name <- "recurrent_miRNAs_TARBP2_up_cancers"


# ============================================================
# 22. Save PNG
# ============================================================

ggsave(
  filename = file.path(
    figure_output_dir,
    paste0(
      figure_name,
      ".png"
    )
  ),
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)




# ============================================================
# 23. Save TIFF
# ============================================================

ggsave(
  filename = file.path(
    figure_output_dir,
    paste0(
      figure_name,
      ".tiff"
    )
  ),
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

