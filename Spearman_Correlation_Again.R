#Recreating Spearman Correlation Again...

# ================================================================ #
#     Rebuild DIABLO Model from Saved Omics Matrices (YNG vs OLD)  #
# ================================================================ #
# Author: Matthew D. Bruss
# ================================================================ #

setwd("/Users/mdbruss/Documents/RStudioProjects_2/Rapa_PwR")

suppressPackageStartupMessages({
  library(dplyr)
  library(mixOmics)
})

# ================================================================
# 1. Load previously saved omics data
# ================================================================
logCPM_symbol <- readRDS("logCPM_symbol_yngOLDsed.RDS")   # RNA
metabo_log    <- readRDS("metabo_log_yngOLDsed.RDS")      # Metabolite
lipid_log     <- readRDS("lipid_log_yngOLDsed.RDS")       # Lipid

# ================================================================
# 2. Define sample groups
# ================================================================
yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sedveh_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
all_samples <- c(yng_sed_samples, old_sedveh_samples)

# Check alignment across blocks
stopifnot(identical(colnames(logCPM_symbol), all_samples))
stopifnot(identical(colnames(metabo_log), all_samples))
stopifnot(identical(colnames(lipid_log), all_samples))

# ================================================================
# 3. Build block matrices (samples in rows, features in columns)
# ================================================================
RNA_block  <- t(logCPM_symbol)
Metabo_block <- t(metabo_log)
Lipid_block  <- t(lipid_log)

# Check numeric
stopifnot(is.numeric(RNA_block[1,1]))
stopifnot(is.numeric(Metabo_block[1,1]))
stopifnot(is.numeric(Lipid_block[1,1]))

# ================================================================
# 4. Assemble blocks and outcome variable
# ================================================================
X_blocks <- list(
  RNA = RNA_block,
  Metabolite = Metabo_block,
  Lipid = Lipid_block
)

Y_train_aging <- factor(c(
  rep("YNG", length(yng_sed_samples)),
  rep("OLD", length(old_sedveh_samples))
), levels = c("YNG", "OLD"))
names(Y_train_aging) <- all_samples

# ================================================================
# 5. Design matrix and model parameters
# ================================================================
design_mat <- matrix(1, ncol = 3, nrow = 3,
                     dimnames = list(names(X_blocks), names(X_blocks)))
diag(design_mat) <- 0

# ================================================================
# 6. Run DIABLO model
# ================================================================
set.seed(123)
diablo_model_aging <- block.splsda(
  X_blocks,
  Y_train_aging,
  ncomp = 1,
  keepX = list(RNA = 500, Metabolite = 100, Lipid = 100),  # same as plan
  design = design_mat
)

# ================================================================
# 7. Save and inspect
# ================================================================
saveRDS(diablo_model_aging, "diablo_model_aging.RDS")
summary(diablo_model_aging)
plotDiablo(diablo_model_aging)

# Optional: extract loadings for later use
rna_load <- diablo_model_aging$loadings$RNA[,1]
met_load <- diablo_model_aging$loadings$Metabolite[,1]
lip_load <- diablo_model_aging$loadings$Lipid[,1]

saveRDS(list(rna_load = rna_load, met_load = met_load, lip_load = lip_load),
        "diablo_loadings_aging.RDS")

cat("\n✅ DIABLO model built successfully. Saved as 'diablo_model_aging.RDS'.\n")

diablo_model_aging <- block.splsda(
  X_blocks,
  Y_train_aging,
  ncomp = 2,  # allow plotting Comp1 vs Comp2
  keepX = list(RNA = 500, Metabolite = 100, Lipid = 100),
  design = design_mat
)

plotIndiv(diablo_model_aging,
          comp = c(1,2),
          ind.names = FALSE,
          legend = TRUE,
          title = "DIABLO Aging Axis (Comp1 vs Comp2)",
          ellipse = TRUE)
#=========================================================#
#=========================================================#

#------Spearman Correlation Again------------------------#

#============================================================#
#      Integrated Multi-Omic + Physiology Correlation Map    #
#============================================================#
#============================================================#
#   Integrated Multi-Omic + Physiology Correlation (All Samples)
#============================================================#

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(cluster)
})

#============================================================#
# 1. Load DIABLO model and extract top features
#============================================================#
diablo_model_aging <- readRDS("diablo_model_aging.RDS")

rna_load <- diablo_model_aging$loadings$RNA[, 1, drop = TRUE]
met_load <- diablo_model_aging$loadings$Metabolite[, 1, drop = TRUE]
lip_load <- diablo_model_aging$loadings$Lipid[, 1, drop = TRUE]

top_rna <- names(sort(abs(rna_load), decreasing = TRUE))[1:500]
top_met <- names(sort(abs(met_load), decreasing = TRUE))[1:100]
top_lip <- names(sort(abs(lip_load), decreasing = TRUE))[1:100]

# Clean RNA names (remove possible "RNA_" prefixes or fix casing)
top_rna_clean <- gsub("^RNA_", "", top_rna)
top_rna_clean <- make.unique(top_rna_clean)

#------------------------------------------------------------#
# 2. Load full RNA matrix (logCPM_symbol_all)
#------------------------------------------------------------#
logCPM_symbol_all <- readRDS("logCPM_symbol_final.RDS")

# Match DIABLO-selected genes to available symbols
common_rna <- intersect(top_rna_clean, rownames(logCPM_symbol_all))
missing_rna <- setdiff(top_rna_clean, rownames(logCPM_symbol_all))
cat("Matched", length(common_rna), "of", length(top_rna_clean), "RNA features to full matrix\n")
if (length(missing_rna) > 0) cat("Missing genes (first 10):", paste(head(missing_rna, 10), collapse = ", "), "\n")

# Use only matched genes
RNA_all <- t(logCPM_symbol_all[common_rna, , drop = FALSE])

#------------------------------------------------------------#
# 3. Load and process full metabolomics (Konopka_Muscle_HILIC.csv)
#------------------------------------------------------------#
metabo_raw <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)
metabo_wide <- metabo_raw %>% dplyr::select(-Group)

metabo_matrix <- metabo_wide %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()
rownames(metabo_matrix) <- metabo_matrix$Metabolite
metabo_matrix$Metabolite <- NULL

# Rename samples to T_* codes
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_01","OFR1","T_04","OFR2","T_22","OFR3","T_24","OFR4","T_28","OFR5","T_34","OFR6","T_40","OFR7","T_44","OFR8",
  "T_02","OIR1","T_03","OIR2","T_14","OIR3","T_31","OIR4","T_33","OIR5","T_37","OIR6","T_47","OIR7","T_48","OIR8",
  "T_08","OV1","T_10","OV2","T_11","OV3","T_16","OV4","T_41","OV5","T_43","OV6","T_45","OV7","T_46","OV8",
  "T_06","OS1","T_09","OS2","T_13","OS3","T_19","OS4","T_21","OS5","T_25","OS6","T_29","OS7","T_39","OS8",
  "T_07","YV1","T_12","YV2","T_18","YV3","T_20","YV4","T_27","YV5","T_30","YV6","T_36","YV7","T_42","YV8",
  "T_05","YS1","T_15","YS2","T_17","YS3","T_23","YS4","T_26","YS5","T_32","YS6","T_35","YS7","T_38","YS8"
)
name_map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
colnames(metabo_matrix) <- ifelse(is.na(name_map[colnames(metabo_matrix)]),
                                  colnames(metabo_matrix),
                                  name_map[colnames(metabo_matrix)])

# Impute and log2-transform
metabo_matrix[] <- lapply(metabo_matrix, function(x) as.numeric(as.character(x)))
metabo_imputed <- t(apply(metabo_matrix, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))
metabo_log_all <- as.data.frame(log2(metabo_imputed))

Metabo_all <- t(metabo_log_all[top_met, , drop = FALSE])

#------------------------------------------------------------#
# 4. Load and process full lipidomics (Konopka_muscle_lipids.csv)
#------------------------------------------------------------#
lipid_df <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)
lipid_cols <- setdiff(names(lipid_df), c("Sample","Group"))
lipid_df[ , lipid_cols] <- lapply(lipid_df[ , lipid_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})
lipid_matrix <- t(as.matrix(lipid_df[ , lipid_cols]))
colnames(lipid_matrix) <- lipid_df$Sample
lipid_matrix <- lipid_matrix[ , !grepl("^YV", colnames(lipid_matrix)), drop = FALSE]
map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
colnames(lipid_matrix) <- ifelse(is.na(map[colnames(lipid_matrix)]),
                                 colnames(lipid_matrix),
                                 map[colnames(lipid_matrix)])

Lipid_all <- t(lipid_matrix[top_lip, , drop = FALSE])
dim(lipid_matrix)
head(lipid_matrix)
#------------------------------------------------------------#
# 5. Align sample order across all omics
#------------------------------------------------------------#
common_samples <- Reduce(intersect, list(rownames(RNA_all),
                                         rownames(Metabo_all),
                                         rownames(Lipid_all)))
RNA_all    <- RNA_all[common_samples, , drop = FALSE]
Metabo_all <- Metabo_all[common_samples, , drop = FALSE]
Lipid_all  <- Lipid_all[common_samples, , drop = FALSE]

head(Lipid_all_pwr)
dim(Lipid_all_pwr)

cat("Common samples:", length(common_samples), "\n")

#------------------------------------------------------------#
# 6. Load and reshape physiology data (Triceps_Physio_Data_v2.csv)
#------------------------------------------------------------#

triceps_physio <- read.csv("Triceps_Physio_Data_v2.csv", check.names = FALSE)
# Convert everything except "Sample" to numeric (some columns might be character)
# 2. Replace missing or blank names with synthetic ones
names(physio_long) <- ifelse(
  is.na(names(physio_long)) | names(physio_long) == "",
  paste0("V", seq_along(names(physio_long))),
  names(physio_long)
)

# 3. Now safely convert all but "Sample" to numeric
physio_long <- physio_long %>%
  dplyr::mutate(across(-Sample, as.numeric))

# Check structure again
str(physio_long)

# ------------------------------------------------------------ #
# 0) Libs
# ------------------------------------------------------------ #
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ------------------------------------------------------------ #
# 1) Physio: ensure numeric + rows = samples
# ------------------------------------------------------------ #
# (you already fixed names and numeric conversion on physio_long)
# physio_long has columns V1..V15 and 'Sample' as character IDs (T_05, etc.)

physio_numeric <- physio_long %>%
  tibble::column_to_rownames("Sample") %>%
  as.matrix()           # SAMPLES x VARIABLES (good)

# ------------------------------------------------------------ #
# 2) Align common samples across all blocks
#     RNA_all / Metabo_all / Lipid_all currently have SAMPLES as rows
# ------------------------------------------------------------ #
common_samples <- Reduce(
  intersect,
  list(rownames(RNA_all), rownames(Metabo_all), rownames(Lipid_all), rownames(physio_numeric))
)

stopifnot(length(common_samples) > 1)

RNA_all     <- RNA_all[common_samples, , drop = FALSE]
Metabo_all  <- Metabo_all[common_samples, , drop = FALSE]
Lipid_all   <- Lipid_all[common_samples, , drop = FALSE]
physio_numeric <- physio_numeric[common_samples, , drop = FALSE]

# ------------------------------------------------------------ #
# 3) Convert to FEATURES x SAMPLES for integration
#    (we’ll rbind features later, so rows must be features)
# ------------------------------------------------------------ #
RNA_fx   <- t(as.matrix(RNA_all))        # FEATURES x SAMPLES
Met_fx   <- t(as.matrix(Metabo_all))
Lip_fx   <- t(as.matrix(Lipid_all))
Phys_fx  <- t(as.matrix(physio_numeric))

# ------------------------------------------------------------ #
# 4) Row-wise z-score helper (robust to NAs / constant rows)
# ------------------------------------------------------------ #
row_zscore <- function(m) {
  m <- as.matrix(m)
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}

# ------------------------------------------------------------ #
# 5) Z-score each block (rows = features)
# ------------------------------------------------------------ #
RNA_z  <- row_zscore(RNA_fx)
Met_z  <- row_zscore(Met_fx)
Lip_z  <- row_zscore(Lip_fx)
Phys_z <- row_zscore(Phys_fx)

# Combine all features
integrated_all <- rbind(RNA_z, Met_z, Lip_z, Phys_z)

# ------------------------------------------------------------ #
# 6) Spearman correlation (feature × feature) and clustering
#    cor() computes column-wise correlations, so pass t() so columns=features
# ------------------------------------------------------------ #
cor_matrix <- cor(t(integrated_all), method = "spearman", use = "pairwise.complete.obs")
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

# Force 10 clusters
optimal_k <- 10
cluster_assign <- cutree(hc, k = optimal_k)

# ------------------------------------------------------------ #
# 7) Prepare heatmap annotations
# ------------------------------------------------------------ #
feature_order <- hc$order
cor_matrix_ord <- cor_matrix[feature_order, feature_order]

# upper triangle for plotting (optional aesthetic)
cor_matrix_tri <- cor_matrix_ord
cor_matrix_tri[lower.tri(cor_matrix_tri)] <- NA

# Omics type per feature (by row name membership)
feature_types <- dplyr::case_when(
  rownames(cor_matrix) %in% rownames(RNA_z)  ~ "RNA",
  rownames(cor_matrix) %in% rownames(Met_z)  ~ "Metabolite",
  rownames(cor_matrix) %in% rownames(Lip_z)  ~ "Lipid",
  rownames(cor_matrix) %in% rownames(Phys_z) ~ "Physiology",
  TRUE ~ "Other"
)
feature_types_ord   <- feature_types[feature_order]
cluster_assign_ord  <- cluster_assign[feature_order]

type_colors <- c(
  "RNA" = "#1f78b4",
  "Metabolite" = "#33a02c",
  "Lipid" = "#ff7f00",
  "Physiology" = "#6a3d9a",
  "Other" = "grey70"
)

row_ha <- rowAnnotation(
  Omics = feature_types_ord,
  Cluster = as.factor(cluster_assign_ord),
  col = list(
    Omics = type_colors,
    Cluster = structure(
      circlize::rand_color(optimal_k),
      names = as.character(1:optimal_k)
    )
  ),
  annotation_name_gp = gpar(fontsize = 9)
)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# ------------------------------------------------------------ #
# 8) Draw Heatmap with visible legends for Omics + Clusters
# ------------------------------------------------------------ #

# Build legends
omics_lgd <- Legend(
  title = "Omics Type",
  labels = names(type_colors),
  legend_gp = gpar(fill = type_colors),
  grid_height = unit(4, "mm"),
  grid_width = unit(4, "mm")
)

cluster_colors <- structure(
  circlize::rand_color(optimal_k),
  names = as.character(1:optimal_k)
)

cluster_lgd <- Legend(
  title = "Cluster",
  labels = names(cluster_colors),
  legend_gp = gpar(fill = cluster_colors),
  grid_height = unit(4, "mm"),
  grid_width = unit(4, "mm")
)

# Update right annotation with both bars
row_ha <- rowAnnotation(
  Omics = feature_types_ord,
  Cluster = as.factor(cluster_assign_ord),
  col = list(
    Omics = type_colors,
    Cluster = cluster_colors
  ),
  annotation_name_gp = gpar(fontsize = 9)
)

# Draw heatmap + legends
pdf("Spearman_corr_heatmap_with_physio_FULL.pdf", width = 12, height = 12)
draw(
  Heatmap(
    cor_matrix_tri,
    name = "Spearman\nCorr",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    right_annotation = row_ha,
    column_title = sprintf(
      "Full multi-omic + physiology correlation map (%d clusters, n=%d features)",
      optimal_k, nrow(cor_matrix)
    ),
    na_col = "white"
  ),
  annotation_legend_list = list(omics_lgd, cluster_lgd)
)
dev.off()
# ------------------------------------------------------------ #
# 9) Export correlation and cluster membership tables
# ------------------------------------------------------------ #
cor_df <- as.data.frame(as.table(cor_matrix), stringsAsFactors = FALSE) %>%
  dplyr::filter(Var1 != Var2) %>%
  dplyr::mutate(abs_r = abs(Freq)) %>%
  dplyr::arrange(desc(abs_r)) %>%
  dplyr::rename(Feature1 = Var1, Feature2 = Var2, Corr = Freq)

utils::write.csv(head(cor_df, 100), "top_correlations_physio_FULL.csv", row.names = FALSE)

cluster_table <- data.frame(
  Feature   = rownames(cor_matrix)[feature_order],
  OmicsType = feature_types_ord,
  Cluster   = cluster_assign_ord,
  stringsAsFactors = FALSE
)
utils::write.csv(cluster_table, "feature_cluster_assignments_FULL.csv", row.names = FALSE)

cat("✅ Completed full multi-omic + physiology Spearman correlation map with 10 clusters.\n")
