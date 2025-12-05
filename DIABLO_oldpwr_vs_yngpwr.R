# YngPWR vs OldPWR DIABLO Aging Axis

#==========================================================#
#---Create oldpwr vs yngpwr aging axis---------------------#
#==========================================================#

#=========Construct RNA Maxrix=============================#

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(tidyr)
  library(edgeR);  library(AnnotationDbi); library(org.Mm.eg.db)
})

# ---- Load counts ----
op_vs_ys_counts <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx")
yp_vs_ys_counts <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_YNG_PWR_VEH-YNG_SED_VEH.xlsx")

head(op_vs_ys_counts)
head(yp_vs_ys_counts)

# 1) Rename the first column to "Ensembl"
rename_first_col <- function(df) {
  cn <- colnames(df)
  cn[1] <- "Ensembl"
  colnames(df) <- cn
  df
}
op_vs_ys_counts <- rename_first_col(op_vs_ys_counts)
yp_vs_ys_counts <- rename_first_col(yp_vs_ys_counts)

# ---- Define sample groups ----
yng_pwr_samples <- c("T_07","T_12","T_18","T_20","T_27","T_30","T_36","T_42")
old_pwr_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")


#--- Combine the two count matrices ---
# They both have the same 'Ensembl' column, so we can full_join or inner_join safely.
combined_counts <- full_join(yp_vs_ys_counts, op_vs_ys_counts, by = "Ensembl")

#--- Subset only the columns we want ---
# keep Ensembl + desired sample columns
counts_yngOldPwr <- combined_counts %>%
  dplyr::select(Ensembl, all_of(c(yng_pwr_samples, old_pwr_samples)))

#--- Check output ---
counts_yngOldPwr %>% glimpse()
head(counts_yngOldPwr)


# ---- Load libraries ----
library(dplyr)
library(stringr)
library(org.Mm.eg.db)
library(edgeR)
library(tibble)
library(tidyr)

# ---- Combine and prepare counts ----
# Assumes you already created `counts_yngOldPwr` as before
merged_counts_pwr <- counts_yngOldPwr

# ---- Map Ensembl -> ENTREZ ----
merged_counts_pwr <- merged_counts_pwr %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db,
                           keys = Ensembl_noDec,
                           keytype = "ENSEMBL",
                           column = "ENTREZID",
                           multiVals = "first")) %>%
  tidyr::drop_na(ENTREZID)

# ---- Counts matrix ----
all_pwr_samples <- c(yng_pwr_samples, old_pwr_samples)

counts_matrix_pwr <- merged_counts_pwr %>%
  dplyr::select(dplyr::all_of(all_pwr_samples)) %>%
  replace(is.na(.), 0) %>%
  as.matrix()
rownames(counts_matrix_pwr) <- merged_counts_pwr$ENTREZID

# ---- Define groups ----
group_pwr <- factor(c(rep("YNG_PWR", length(yng_pwr_samples)),
                      rep("OLD_PWR", length(old_pwr_samples))),
                    levels = c("YNG_PWR", "OLD_PWR"))

stopifnot(length(group_pwr) == ncol(counts_matrix_pwr))

# ---- edgeR: TMM normalization, filtering, logCPM ----
dge_pwr <- DGEList(counts = counts_matrix_pwr, group = group_pwr)
dge_pwr <- calcNormFactors(dge_pwr, method = "TMM")
keep_pwr <- filterByExpr(dge_pwr)
dge_pwr <- dge_pwr[keep_pwr, , keep.lib.sizes = FALSE]
logCPM_matrix_pwr <- cpm(dge_pwr, log = TRUE, prior.count = 1)

cat("Genes kept:", nrow(logCPM_matrix_pwr), " Samples:", ncol(logCPM_matrix_pwr), "\n")

# ---- Map ENTREZ -> SYMBOL and collapse duplicates ----
symbols_pwr <- mapIds(org.Mm.eg.db,
                      keys    = rownames(logCPM_matrix_pwr),
                      keytype = "ENTREZID",
                      column  = "SYMBOL",
                      multiVals = "first")

sym_df_pwr <- as.data.frame(logCPM_matrix_pwr) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  mutate(Symbol = symbols_pwr) %>%
  tidyr::drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

logCPM_symbol_pwr <- sym_df_pwr %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()

# ---- Keep exact sample order ----
logCPM_symbol_pwr <- logCPM_symbol_pwr[, all_pwr_samples, drop = FALSE]

# ---- Sanity checks ----
stopifnot(!anyNA(logCPM_symbol_pwr))
stopifnot(identical(colnames(logCPM_symbol_pwr), all_pwr_samples))

# ---- Save normalized matrix ----
saveRDS(logCPM_symbol_pwr, "logCPM_symbol_yngOLDpwr.RDS")

# ---- Quick summary ----
dim(logCPM_symbol_pwr)
head(logCPM_symbol_pwr)
#---------------------------------------------------------#

#============================================================#
#---------Build metab_log_yngOLDpwr -------------------------#
#============================================================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

# ---- Load metabolomics data ----
pwr_muscle_metabo_data <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)

# ---- Drop Group column, keep Sample + metabolites ----
metabo_wide_pwr <- pwr_muscle_metabo_data %>%
  dplyr::select(-Group)

# ---- Convert long → wide: rows = metabolites, cols = samples ----
metabo_matrix_pwr <- metabo_wide_pwr %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()

rownames(metabo_matrix_pwr) <- metabo_matrix_pwr$Metabolite
metabo_matrix_pwr$Metabolite <- NULL

# ---- Map metabolomics names (YV*/OV*) → T_* tube codes ----
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_07","YV1","T_12","YV2","T_18","YV3","T_20","YV4","T_27","YV5","T_30","YV6","T_36","YV7","T_42","YV8",
  "T_08","OV1","T_10","OV2","T_11","OV3","T_16","OV4","T_41","OV5","T_43","OV6","T_45","OV7","T_46","OV8"
)

name_map_pwr <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
new_names_pwr <- name_map_pwr[colnames(metabo_matrix_pwr)]
colnames(metabo_matrix_pwr) <- ifelse(is.na(new_names_pwr), colnames(metabo_matrix_pwr), new_names_pwr)

# ---- Define groups ----
yng_pwr_samples <- c("T_07","T_12","T_18","T_20","T_27","T_30","T_36","T_42")
old_pwr_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
all_pwr_samples <- c(yng_pwr_samples, old_pwr_samples)

# ---- Subset to only PoWeR samples ----
present_metabo_pwr <- intersect(all_pwr_samples, colnames(metabo_matrix_pwr))
missing_metabo_pwr <- setdiff(all_pwr_samples, colnames(metabo_matrix_pwr))
if (length(missing_metabo_pwr)) message("Metabolomics missing: ", paste(missing_metabo_pwr, collapse = ", "))

metabo_matrix_subset_pwr <- metabo_matrix_pwr[, present_metabo_pwr, drop = FALSE]

# ---- Ensure numeric ----
metabo_matrix_subset_pwr[] <- lapply(metabo_matrix_subset_pwr, function(x) as.numeric(as.character(x)))

# ---- Impute zeros/NA with half min positive per metabolite ----
metabo_imputed_pwr <- t(apply(metabo_matrix_subset_pwr, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))

# ---- Log2 transform ----
metabo_log_pwr <- as.data.frame(log2(metabo_imputed_pwr))
colnames(metabo_log_pwr) <- present_metabo_pwr
rownames(metabo_log_pwr) <- rownames(metabo_matrix_subset_pwr)

# ---- Align with RNA samples if available ----
if (exists("logCPM_symbol_pwr")) {
  common_samples_pwr <- intersect(all_pwr_samples, colnames(logCPM_symbol_pwr))
  metabo_log_pwr     <- metabo_log_pwr[, common_samples_pwr, drop = FALSE]
  logCPM_symbol_pwr  <- logCPM_symbol_pwr[, common_samples_pwr, drop = FALSE]
  stopifnot(identical(colnames(metabo_log_pwr), colnames(logCPM_symbol_pwr)))
}

# ---- Save ----
saveRDS(metabo_log_pwr, "metabo_log_yngOLDpwr.RDS")

# ---- Quick peek ----
dim(metabo_log_pwr)
head(metabo_log_pwr)
#--------------------------------------------------------------------------#

#============================================================#
#---------Build lipid_log_yngOLDpwr --------------------------#
#============================================================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble)
})

# ---- Load lipidomics data ----
lipid_df_pwr <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)
# Keep all samples (no filtering yet)

# ---- Identify lipid columns ----
lipid_cols_pwr <- setdiff(names(lipid_df_pwr), c("Sample", "Group"))

# ---- Impute zeros/NA per lipid, then log2 ----
lipid_df_pwr[ , lipid_cols_pwr] <- lapply(lipid_df_pwr[ , lipid_cols_pwr], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})

# ---- Build matrix: rows = lipids, cols = samples ----
lipid_matrix_pwr <- t(as.matrix(lipid_df_pwr[ , lipid_cols_pwr]))
colnames(lipid_matrix_pwr) <- lipid_df_pwr$Sample

# ---- Map lipidomics sample names (YV*/OV*) → T_* tube codes ----
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_07","YV1","T_12","YV2","T_18","YV3","T_20","YV4","T_27","YV5","T_30","YV6","T_36","YV7","T_42","YV8",
  "T_08","OV1","T_10","OV2","T_11","OV3","T_16","OV4","T_41","OV5","T_43","OV6","T_45","OV7","T_46","OV8"
)

lipid_name_map_pwr <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
mapped_names_pwr <- lipid_name_map_pwr[colnames(lipid_matrix_pwr)]
colnames(lipid_matrix_pwr) <- ifelse(is.na(mapped_names_pwr), colnames(lipid_matrix_pwr), mapped_names_pwr)

# ---- Define sample groups ----
yng_pwr_samples <- c("T_07","T_12","T_18","T_20","T_27","T_30","T_36","T_42")
old_pwr_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
all_pwr_samples <- c(yng_pwr_samples, old_pwr_samples)

# ---- Subset to only PoWeR samples ----
present_lipid_pwr <- intersect(all_pwr_samples, colnames(lipid_matrix_pwr))
missing_lipid_pwr <- setdiff(all_pwr_samples, colnames(lipid_matrix_pwr))
if (length(missing_lipid_pwr)) message("Lipidomics missing: ", paste(missing_lipid_pwr, collapse = ", "))

lipid_log_pwr <- lipid_matrix_pwr[ , present_lipid_pwr, drop = FALSE]

# ---- Align with RNA samples if available ----
if (exists("logCPM_symbol_pwr")) {
  common_samples_pwr <- intersect(all_pwr_samples, colnames(logCPM_symbol_pwr))
  lipid_log_pwr      <- lipid_log_pwr[, common_samples_pwr, drop = FALSE]
  logCPM_symbol_pwr  <- logCPM_symbol_pwr[, common_samples_pwr, drop = FALSE]
  stopifnot(identical(colnames(lipid_log_pwr), colnames(logCPM_symbol_pwr)))
}

# ---- Save ----
saveRDS(lipid_log_pwr, "lipid_log_yngOLDpwr.RDS")

# ---- Sanity check ----
dim(lipid_log_pwr)
head(lipid_log_pwr)
#-------------------------------------------------------------#

#============================================================#
#--------- DIABLO: Young vs Old PoWeR (Aging Axis) ----------#
#============================================================#

suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
})

#------------------------------------------------------------#
# 1. Build omics list (samples in rows, features in columns)
#------------------------------------------------------------#

X_pwr <- list(
  RNA        = t(logCPM_symbol_pwr),
  Metabolite = t(metabo_log_pwr),
  Lipid      = t(lipid_log_pwr)
)

Y_pwr <- factor(c(
  rep("YNG", length(yng_pwr_samples)),
  rep("OLD", length(old_pwr_samples))
), levels = c("YNG", "OLD"))

names(Y_pwr) <- rownames(X_pwr$RNA)

# ---- Alignment checks ----
stopifnot(
  identical(rownames(X_pwr$RNA), rownames(X_pwr$Metabolite)),
  identical(rownames(X_pwr$RNA), rownames(X_pwr$Lipid)),
  identical(rownames(X_pwr$RNA), names(Y_pwr))
)

#------------------------------------------------------------#
# 2. Design matrix (full integration)
#------------------------------------------------------------#

design_pwr <- matrix(1, ncol = length(X_pwr), nrow = length(X_pwr),
                     dimnames = list(names(X_pwr), names(X_pwr)))
diag(design_pwr) <- 0
design_pwr

#------------------------------------------------------------#
# 3. Run DIABLO
#------------------------------------------------------------#

diablo_res_pwr <- block.splsda(
  X = X_pwr,
  Y = Y_pwr,
  ncomp = 1,
  keepX = list(RNA = 50, Metabolite = 20, Lipid = 20),
  design = design_pwr
)

#------------------------------------------------------------#
# 4. Extract selected features
#------------------------------------------------------------#

extract_features <- function(sel, block) {
  data.frame(
    block   = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    abs_loading = abs(sel[[block]]$value$value.var),
    stringsAsFactors = FALSE
  )
}

rna_sel_pwr <- selectVar(diablo_res_pwr, block = "RNA", comp = 1)
met_sel_pwr <- selectVar(diablo_res_pwr, block = "Metabolite", comp = 1)
lip_sel_pwr <- selectVar(diablo_res_pwr, block = "Lipid", comp = 1)

selected_df_pwr <- bind_rows(
  extract_features(rna_sel_pwr, "RNA"),
  extract_features(met_sel_pwr, "Metabolite"),
  extract_features(lip_sel_pwr, "Lipid")
) %>%
  arrange(block, desc(abs_loading))

# Save for inspection
write.csv(selected_df_pwr, "DIABLO_selected_features_yngOLDpwr.csv", row.names = FALSE)
head(selected_df_pwr, 15)

#------------------------------------------------------------#
# 5. Visualization
#------------------------------------------------------------#

# Cross-block circos plot
circosPlot(diablo_res_pwr, comp = 1, cutoff = 0.7, line = TRUE)

# Network views
network(diablo_res_pwr, comp = list(RNA = 1, Metabolite = 1), cutoff = 0.6)
network(diablo_res_pwr, comp = list(RNA = 1, Lipid = 1), cutoff = 0.6)
network(diablo_res_pwr, comp = list(Metabolite = 1, Lipid = 1), cutoff = 0.6)
network(diablo_res_pwr, comp = list(1,1), cutoff = 0.6)

# Loadings per block
plotLoadings(diablo_res_pwr, block = "RNA", comp = 1, method = "mean")
plotLoadings(diablo_res_pwr, block = "Metabolite", comp = 1, method = "mean")
plotLoadings(diablo_res_pwr, block = "Lipid", comp = 1, method = "mean")

#------------------------------------------------------------#
# 6. Ranked loadings visualization
#------------------------------------------------------------#

ranked_df_pwr <- selected_df_pwr %>%
  group_by(block) %>%
  arrange(desc(abs_loading), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Histogram of loadings
ggplot(ranked_df_pwr, aes(x = abs_loading)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~block, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of absolute loadings (PoWeR Aging Axis)",
       x = "Absolute loading", y = "Count")

# Ranked loadings plot
ggplot(ranked_df_pwr, aes(x = rank, y = abs_loading, color = block)) +
  geom_point() + geom_line(aes(group = block)) +
  theme_minimal(base_size = 14) +
  labs(title = "Ranked feature loadings per block (PoWeR)",
       x = "Feature rank within block", y = "|Loading|") +
  scale_color_brewer(palette = "Dark2")

# Top drivers (>0.2)
strong_feats_pwr <- ranked_df_pwr %>%
  filter(abs_loading > 0.2) %>%
  arrange(block, desc(abs_loading))

ggplot(strong_feats_pwr, aes(x = reorder(feature, abs_loading),
                             y = abs_loading, fill = block)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = round(abs_loading,3)), hjust = -0.1, size = 3) +
  facet_wrap(~block, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Top DIABLO drivers (|loading| > 0.2)",
       x = "Feature", y = "|Loading|")

#------------------------------------------------------------#
# Save DIABLO object
saveRDS(diablo_res_pwr, "DIABLO_res_yngOLDpwr.RDS")

diablo_res_pwr <- block.splsda(
  X = X_pwr,
  Y = Y_pwr,
  ncomp = 2,  # 👈 two components for plotting
  keepX = list(RNA = c(50, 50),   # same number per comp is fine
               Metabolite = c(20, 20),
               Lipid = c(20, 20)),
  design = design_pwr
)

# Base correlation circle plot (scores correlation between blocks)
plotDiablo(diablo_res_pwr)

# 3×3 sample correlation matrix across blocks
# Each panel shows sample projections colored by group
plotIndiv(
  diablo_res_pwr,
  ind.names = FALSE,
  legend = TRUE,
  ellipse = TRUE,
  title = "Young vs Old PoWeR: Cross-block sample projections",
  col = c("#1F78B4", "#E31A1C"),   # Blue = YNG, Red = OLD
  cex = 1.2,
  style = "graphics"
)
#-------------------------------------------------------------------#

selected_df_pwr
#--------------------------------------------------------------------#

#============================================================#
#---------DIABLO Ingegrated Aging Axis------------------#
#============================================================#

#============================================================#
#----- Combine transcriptomic matrices (all interventions) ---#
#============================================================#

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- Define sample groups ----
yng_pwr_samples      <- c("T_07","T_12","T_18","T_20","T_27","T_30","T_36","T_42")
old_pwr_samples      <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_sed_samples      <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwr_irap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwr_frap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

# ---- Define final desired column order ----
final_order <- c(
  yng_pwr_samples,
  old_pwr_samples,
  old_sed_samples,
  old_pwr_irap_samples,
  old_pwr_frap_samples
)

# ---- Merge the two matrices by rownames (gene symbols) ----
# Make sure rownames are identical type
common_genes <- intersect(rownames(logCPM_symbol_all), rownames(logCPM_symbol_pwr))

logCPM_symbol_all <- logCPM_symbol_all[common_genes, , drop = FALSE]
logCPM_symbol_pwr <- logCPM_symbol_pwr[common_genes, , drop = FALSE]

# Combine columns from both matrices
logCPM_symbol_full <- cbind(
  logCPM_symbol_all,
  logCPM_symbol_pwr[, setdiff(colnames(logCPM_symbol_pwr), colnames(logCPM_symbol_all)), drop = FALSE]
)

# ---- Keep only desired samples in defined order ----
present_samples <- intersect(final_order, colnames(logCPM_symbol_full))
missing_samples <- setdiff(final_order, colnames(logCPM_symbol_full))

if (length(missing_samples) > 0) {
  message("⚠️ Missing samples: ", paste(missing_samples, collapse = ", "))
}

logCPM_symbol_full <- logCPM_symbol_full[, present_samples, drop = FALSE]

# ---- Sanity checks ----
stopifnot(identical(colnames(logCPM_symbol_full), present_samples))
cat("Final matrix dimensions:", dim(logCPM_symbol_full), "\n")

# ---- Save merged matrix ----
saveRDS(logCPM_symbol_full, "logCPM_symbol_full_interventions.RDS")

# Preview
head(logCPM_symbol_full)
#--------------------------------------------------------------------#

#============================================================#
#----- Subset RNA matrix to 50 DIABLO-selected features ------#
#============================================================#

suppressPackageStartupMessages({
  library(dplyr)
})

# 1) Verify DIABLO RNA feature set from the current analysis
rna_features_50_pwr <- rna_sel_pwr$RNA$name

cat("DIABLO-selected RNA features (PoWeR):", length(rna_features_50_pwr), "\n")

# 2) Confirm presence in the full logCPM matrix
common_rna_feats_pwr  <- intersect(rna_features_50_pwr, rownames(logCPM_symbol_full))
missing_rna_feats_pwr <- setdiff(rna_features_50_pwr, rownames(logCPM_symbol_full))

cat("Overlap with full dataset:", length(common_rna_feats_pwr), "\n")
if (length(missing_rna_feats_pwr) > 0) {
  cat("⚠️ Missing RNA features:", paste(missing_rna_feats_pwr, collapse = ", "), "\n")
}

# 3) Subset and transpose: rows = samples, cols = 50 genes
RNA_all_pwr <- t(logCPM_symbol_full[common_rna_feats_pwr, , drop = FALSE])

# 4) Sanity checks
cat("Final RNA_all_pwr dimensions:", dim(RNA_all_pwr), "\n")
cat("Samples:", nrow(RNA_all_pwr), " | Features:", ncol(RNA_all_pwr), "\n")

head(rownames(RNA_all_pwr))   # should be T_* sample IDs
head(colnames(RNA_all_pwr))   # should be gene symbols
stopifnot(is.numeric(RNA_all_pwr[1,1]))

# 5) Save subsetted RNA matrix
saveRDS(RNA_all_pwr, "RNA_all_50features_yngOLDpwr.RDS")
#---------------------------------------------------------#

#============================================================#
#----- Build full metabolomics subset (20 DIABLO features) ---#
#============================================================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

#------------------------------------------------------------#
# 1. Load metabolomics data
#------------------------------------------------------------#
metabo_raw_pwr <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)

# Keep only Sample + metabolite features (drop Group)
metabo_wide_pwr <- metabo_raw_pwr %>% dplyr::select(-Group)

#------------------------------------------------------------#
# 2. Reshape: rows = metabolites, columns = samples
#------------------------------------------------------------#
metabo_matrix_pwr <- metabo_wide_pwr %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()

rownames(metabo_matrix_pwr) <- metabo_matrix_pwr$Metabolite
metabo_matrix_pwr$Metabolite <- NULL

#------------------------------------------------------------#
# 3. Map metabolomics sample names (YV, OV, OS, OIR, OFR) → T_* codes
#------------------------------------------------------------#
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_01", "OFR1", "T_04", "OFR2", "T_22", "OFR3", "T_24", "OFR4",
  "T_28", "OFR5", "T_34", "OFR6", "T_40", "OFR7", "T_44", "OFR8",
  "T_02", "OIR1", "T_03", "OIR2", "T_14", "OIR3", "T_31", "OIR4",
  "T_33", "OIR5", "T_37", "OIR6", "T_47", "OIR7", "T_48", "OIR8",
  "T_08", "OV1", "T_10", "OV2", "T_11", "OV3", "T_16", "OV4",
  "T_41", "OV5", "T_43", "OV6", "T_45", "OV7", "T_46", "OV8",
  "T_06", "OS1", "T_09", "OS2", "T_13", "OS3", "T_19", "OS4",
  "T_21", "OS5", "T_25", "OS6", "T_29", "OS7", "T_39", "OS8",
  "T_07", "YV1", "T_12", "YV2", "T_18", "YV3", "T_20", "YV4",
  "T_27", "YV5", "T_30", "YV6", "T_36", "YV7", "T_42", "YV8"
)

name_map_pwr <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
new_names_pwr <- name_map_pwr[colnames(metabo_matrix_pwr)]
colnames(metabo_matrix_pwr) <- ifelse(is.na(new_names_pwr),
                                      colnames(metabo_matrix_pwr),
                                      new_names_pwr)

#------------------------------------------------------------#
# 4. Ensure numeric, impute zeros/NA, then log2-transform
#------------------------------------------------------------#
metabo_matrix_pwr[] <- lapply(metabo_matrix_pwr, function(x) as.numeric(as.character(x)))

metabo_imputed_pwr <- t(apply(metabo_matrix_pwr, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))

metabo_log_full_pwr <- as.data.frame(log2(metabo_imputed_pwr))

#------------------------------------------------------------#
# 5. Subset to DIABLO-selected metabolites (20 features)
#------------------------------------------------------------#
met_features_20_pwr <- met_sel_pwr$Metabolite$name
common_met_feats_pwr  <- intersect(met_features_20_pwr, rownames(metabo_log_full_pwr))
missing_met_feats_pwr <- setdiff(met_features_20_pwr, rownames(metabo_log_full_pwr))

cat("DIABLO-selected metabolites (PoWeR):", length(met_features_20_pwr), "\n")
cat("Overlap with full dataset:", length(common_met_feats_pwr), "\n")
if (length(missing_met_feats_pwr) > 0) {
  cat("⚠️ Missing metabolite features:", paste(missing_met_feats_pwr, collapse = ", "), "\n")
}

# Subset to the matched 20
metabo_log_sub_pwr <- metabo_log_full_pwr[common_met_feats_pwr, , drop = FALSE]

#------------------------------------------------------------#
# 6. Align with RNA_all_pwr sample order
#------------------------------------------------------------#
Metabo_all_pwr <- t(metabo_log_sub_pwr)

common_samples_pwr <- intersect(rownames(RNA_all_pwr), rownames(Metabo_all_pwr))
missing_samples_pwr <- setdiff(rownames(RNA_all_pwr), rownames(Metabo_all_pwr))
if (length(missing_samples_pwr)) {
  message("⚠️ Missing metabolomics samples: ", paste(missing_samples_pwr, collapse = ", "))
}

# Subset and reorder to match RNA sample order
Metabo_all_pwr <- Metabo_all_pwr[common_samples_pwr, , drop = FALSE]
Metabo_all_pwr <- Metabo_all_pwr[rownames(RNA_all_pwr), , drop = FALSE]

#------------------------------------------------------------#
# 7. Save and sanity check
#------------------------------------------------------------#
saveRDS(Metabo_all_pwr, "Metabo_all_20features_yngOLDpwr.RDS")

cat("Final Metabo_all_pwr dimensions:", dim(Metabo_all_pwr), "\n")
head(rownames(Metabo_all_pwr))
head(colnames(Metabo_all_pwr))
#---------------------------------------------#

#============================================================#
#----- Build full lipidomics subset (20 DIABLO features) -----#
#============================================================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(stringr)
})

#------------------------------------------------------------#
# 1. Load lipidomics data
#------------------------------------------------------------#
lipid_df_pwr <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)

#------------------------------------------------------------#
# 2. Identify lipid columns
#------------------------------------------------------------#
lipid_cols_pwr <- setdiff(names(lipid_df_pwr), c("Sample", "Group"))

# Build matrix: rows = lipids, cols = samples
lipid_matrix_pwr <- t(as.matrix(lipid_df_pwr[, lipid_cols_pwr]))
colnames(lipid_matrix_pwr) <- lipid_df_pwr$Sample
rownames(lipid_matrix_pwr) <- lipid_cols_pwr

#------------------------------------------------------------#
# 3. Map lipidomics sample names (YV, OV, OS, OIR, OFR) → T_* codes
#------------------------------------------------------------#
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_01", "OFR1", "T_04", "OFR2", "T_22", "OFR3", "T_24", "OFR4",
  "T_28", "OFR5", "T_34", "OFR6", "T_40", "OFR7", "T_44", "OFR8",
  "T_02", "OIR1", "T_03", "OIR2", "T_14", "OIR3", "T_31", "OIR4",
  "T_33", "OIR5", "T_37", "OIR6", "T_47", "OIR7", "T_48", "OIR8",
  "T_08", "OV1", "T_10", "OV2", "T_11", "OV3", "T_16", "OV4",
  "T_41", "OV5", "T_43", "OV6", "T_45", "OV7", "T_46", "OV8",
  "T_06", "OS1", "T_09", "OS2", "T_13", "OS3", "T_19", "OS4",
  "T_21", "OS5", "T_25", "OS6", "T_29", "OS7", "T_39", "OS8",
  "T_07", "YV1", "T_12", "YV2", "T_18", "YV3", "T_20", "YV4",
  "T_27", "YV5", "T_30", "YV6", "T_36", "YV7", "T_42", "YV8"
)

lipid_name_map_pwr <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
new_names_pwr <- lipid_name_map_pwr[colnames(lipid_matrix_pwr)]
colnames(lipid_matrix_pwr) <- ifelse(is.na(new_names_pwr),
                                     colnames(lipid_matrix_pwr),
                                     new_names_pwr)

#------------------------------------------------------------#
# 4. Subset to DIABLO-selected lipid features (20 features)
#------------------------------------------------------------#
lip_features_20_pwr <- lip_sel_pwr$Lipid$name
common_lip_feats_pwr  <- intersect(lip_features_20_pwr, rownames(lipid_matrix_pwr))
missing_lip_feats_pwr <- setdiff(lip_features_20_pwr, rownames(lipid_matrix_pwr))

cat("DIABLO-selected lipids (PoWeR):", length(lip_features_20_pwr), "\n")
cat("Overlap with lipidomics dataset:", length(common_lip_feats_pwr), "\n")
if (length(missing_lip_feats_pwr) > 0) {
  cat("⚠️ Missing lipid features:", paste(missing_lip_feats_pwr, collapse = ", "), "\n")
}

# Subset to matched 20 lipids
lipid_log_sub_pwr <- lipid_matrix_pwr[common_lip_feats_pwr, , drop = FALSE]

#------------------------------------------------------------#
# 5. Align with RNA_all_pwr sample order
#------------------------------------------------------------#
Lipid_all_pwr <- t(lipid_log_sub_pwr)

common_samples_pwr <- intersect(rownames(RNA_all_pwr), rownames(Lipid_all_pwr))
missing_samples_pwr <- setdiff(rownames(RNA_all_pwr), rownames(Lipid_all_pwr))
if (length(missing_samples_pwr)) {
  message("⚠️ Missing lipidomics samples: ", paste(missing_samples_pwr, collapse = ", "))
}

# Subset and reorder to match RNA sample order
Lipid_all_pwr <- Lipid_all_pwr[common_samples_pwr, , drop = FALSE]
Lipid_all_pwr <- Lipid_all_pwr[rownames(RNA_all_pwr), , drop = FALSE]

#------------------------------------------------------------#
# 6. Save and sanity check
#------------------------------------------------------------#
saveRDS(Lipid_all_pwr, "Lipid_all_20features_yngOLDpwr.RDS")

cat("Final Lipid_all_pwr dimensions:", dim(Lipid_all_pwr), "\n")
head(rownames(Lipid_all_pwr))
head(colnames(Lipid_all_pwr))
#------------------------------------------------------------------#

#============================================================#
#---- Check structure of full DIABLO PoWeR data set ----------#
#============================================================#

cat("### Checking structure and alignment of PoWeR DIABLO matrices ###\n\n")

# 1) Dimensions
cat("RNA_all_pwr:     ", dim(RNA_all_pwr)[1], "samples x", dim(RNA_all_pwr)[2], "features\n")
cat("Metabo_all_pwr:  ", dim(Metabo_all_pwr)[1], "samples x", dim(Metabo_all_pwr)[2], "features\n")
cat("Lipid_all_pwr:   ", dim(Lipid_all_pwr)[1], "samples x", dim(Lipid_all_pwr)[2], "features\n")

# 2) Samples must match across all omics
stopifnot(identical(rownames(RNA_all_pwr), rownames(Metabo_all_pwr)))
stopifnot(identical(rownames(RNA_all_pwr), rownames(Lipid_all_pwr)))

cat("\n✅ Sample names match perfectly across RNA, Metabolite, and Lipid matrices.\n")

# 3) Quick preview of identifiers
head_samples_pwr <- head(rownames(RNA_all_pwr))
cat("\nFirst few samples:\n")
print(head_samples_pwr)

cat("\nRNA features (first 5):\n"); print(head(colnames(RNA_all_pwr)))
cat("\nMetabolite features (first 5):\n"); print(head(colnames(Metabo_all_pwr)))
cat("\nLipid features (first 5):\n"); print(head(colnames(Lipid_all_pwr)))

# 4) Numeric structure check
cat("\nData type checks:\n")
cat("RNA_all_pwr numeric:", is.numeric(RNA_all_pwr[1, 1]), "\n")
cat("Metabo_all_pwr numeric:", is.numeric(Metabo_all_pwr[1, 1]), "\n")
cat("Lipid_all_pwr numeric:", is.numeric(Lipid_all_pwr[1, 1]), "\n")

cat("\n✅ All omics matrices are numeric and aligned.\n")
#////////////////////////////////////////////////////////////////#

#===========================================================
# DIABLO PoWeR Aging Axis Projection (RNA + Metabolite + Lipid)
#===========================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
  library(broom)
})

#-----------------------------------------------------------
# 0) Define sample groups
#-----------------------------------------------------------
yng_pwr_samples     <- c("T_07","T_12","T_18","T_20","T_27","T_30","T_36","T_42")
old_pwr_samples     <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwr_irap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwr_frap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")
old_sed_samples      <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")

#-----------------------------------------------------------
# 1) Build multi-omics list (already normalized)
#-----------------------------------------------------------
X_all_pwr <- list(
  RNA        = RNA_all_pwr,
  Metabolite = Metabo_all_pwr,
  Lipid      = Lipid_all_pwr
)

#-----------------------------------------------------------
# 2) Define outcome Y for training (YNG vs OLD PoWeR)
#-----------------------------------------------------------
Y_train_pwr <- factor(c(
  rep("YNG", length(yng_pwr_samples)),
  rep("OLD", length(old_pwr_samples))
), levels = c("YNG", "OLD"))

# Restrict training data to YNG vs OLD PoWeR
train_samples_pwr <- c(yng_pwr_samples, old_pwr_samples)
X_train_pwr <- lapply(X_all_pwr, function(m) m[train_samples_pwr, , drop = FALSE])

#-----------------------------------------------------------
# 3) Run DIABLO
#-----------------------------------------------------------
# Full correlation design
design_pwr <- matrix(1, ncol = 3, nrow = 3,
                     dimnames = list(names(X_all_pwr), names(X_all_pwr)))
diag(design_pwr) <- 0

# NOTE: 44 RNA overlapped, so adjust keepX accordingly
n_rna_features <- ncol(RNA_all_pwr)
diablo_res_full <- block.splsda(
  X_train_pwr,
  Y_train_pwr,
  ncomp = 1,
  keepX = list(RNA = n_rna_features, Metabolite = 20, Lipid = 20),
  design = design_pwr
)

#-----------------------------------------------------------
# 4) Project ALL samples into DIABLO space
#-----------------------------------------------------------
proj_pwr <- predict(diablo_res_full, newdata = X_all_pwr)
scores_all_blocks <- lapply(proj_pwr$variates, function(block) block[, 1])

# Combine into a sample x block matrix
scores_mat_pwr <- do.call(cbind, scores_all_blocks)

# Average across blocks to obtain a single consensus component 1 per sample
Comp1_pwr <- rowMeans(scores_mat_pwr)

# Build projection dataframe
proj_df_pwr <- data.frame(
  Sample = rownames(scores_mat_pwr),
  Comp1  = Comp1_pwr
)

#-----------------------------------------------------------
# 5) Annotate experimental groups
#-----------------------------------------------------------
proj_df_pwr <- proj_df_pwr %>%
  mutate(Group = case_when(
    Sample %in% yng_pwr_samples      ~ "Young_PoWeR",
    Sample %in% old_pwr_samples      ~ "Old_PoWeR",
    Sample %in% old_pwr_irap_samples ~ "Old_PoWeR_IRap",
    Sample %in% old_pwr_frap_samples ~ "Old_PoWeR_FRap",
    Sample %in% old_sed_samples      ~ "Old_Sed",
    TRUE                             ~ "Other"
  ))

# Ensure consistent order
proj_df_pwr <- proj_df_pwr %>%
  mutate(Group = factor(Group,
                        levels = c("Young_PoWeR", "Old_PoWeR", "Old_Sed",
                                   "Old_PoWeR_IRap", "Old_PoWeR_FRap")))

#-----------------------------------------------------------
# 6) Plot: Boxplot + jitter of projection scores
#-----------------------------------------------------------
ggplot(proj_df_pwr, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Projection of Samples onto PoWeR DIABLO Aging Axis (Comp1)",
       y = "DIABLO Component 1 Score", x = "") +
  theme(legend.position = "none")



#-----------------------------------------------------------
# 7) Statistics: ANOVA + Tukey post-hoc
#-----------------------------------------------------------
anova_res_pwr <- anova_test(data = proj_df_pwr, dv = Comp1, between = Group)
anova_res_pwr

tukey_res_pwr <- proj_df_pwr %>%
  tukey_hsd(Comp1 ~ Group)
tukey_res_pwr

#-----------------------------------------------------------
# 8) Define comparisons to display
#-----------------------------------------------------------
comparisons_pwr <- list(
  c("Old_PoWeR", "Young_PoWeR"),
  c("Old_PoWeR", "Old_Sed"),
  c("Old_PoWeR", "Old_PoWeR_IRap"),
  c("Old_PoWeR", "Old_PoWeR_FRap")
)

#-----------------------------------------------------------
# 9) Final annotated boxplot
#-----------------------------------------------------------
DIABLO_PoWeR_axis_plot <- ggplot(proj_df_pwr, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Projection of All Samples onto PoWeR DIABLO Aging Axis (Comp1)",
       y = "DIABLO Component 1 Score", x = "") +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = comparisons_pwr,
                     method = "t.test",
                     label = "p.signif")

# Save to PDF
pdf("DIABLO_PoWeR_AgingAxis_Boxplot.pdf", width = 6, height = 4)
print(DIABLO_PoWeR_axis_plot)
dev.off()

#--------------------------------------------------------#

#--------------------------------------------
# DIABLO Aging Axis Direction Check
#--------------------------------------------
train_scores <- diablo_res_full$variates$RNA[, 1]
label_numeric <- as.numeric(Y_train_pwr) - 1

correlation_check <- cor(train_scores, label_numeric)
cat("Correlation between Comp1 and OLD label:", correlation_check, "\n")

# Flip direction if needed
if (correlation_check < 0) {
  message("🔄 Flipping DIABLO axis so that OLD = higher scores")
  for (blk in names(diablo_res_full$variates)) {
    diablo_res_full$variates[[blk]][, 1] <- -diablo_res_full$variates[[blk]][, 1]
    diablo_res_full$loadings[[blk]][, 1]  <- -diablo_res_full$loadings[[blk]][, 1]
  }
  proj_pwr <- predict(diablo_res_full, newdata = X_all_pwr)
  scores_all_blocks <- lapply(proj_pwr$variates, function(block) block[, 1])
  scores_mat_pwr <- do.call(cbind, scores_all_blocks)
  Comp1_pwr <- rowMeans(scores_mat_pwr)
  proj_df_pwr$Comp1 <- Comp1_pwr
}


























#////////////////////////////////////////////////////////////#
head(logCPM_symbol_all)
head(logCPM_symbol_pwr)

# ---- Define groups ----
yng_pwr_samples <- c("T_07","T_12","T_18","T_20","T_27","T_30","T_36","T_42")
old_pwr_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwr_irap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwr_frap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")
head(logCPM_symbol_all)

# Mapping table based on your provided key
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_01", "OFR1",
  "T_04", "OFR2",
  "T_22", "OFR3",
  "T_24", "OFR4",
  "T_28", "OFR5",
  "T_34", "OFR6",
  "T_40", "OFR7",
  "T_44", "OFR8",
  "T_02", "OIR1",
  "T_03", "OIR2",
  "T_14", "OIR3",
  "T_31", "OIR4",
  "T_33", "OIR5",
  "T_37", "OIR6",
  "T_47", "OIR7",
  "T_48", "OIR8",
  "T_08", "OV1",
  "T_10", "OV2",
  "T_11", "OV3",
  "T_16", "OV4",
  "T_41", "OV5",
  "T_43", "OV6",
  "T_45", "OV7",
  "T_46", "OV8",
  "T_06", "OS1",
  "T_09", "OS2",
  "T_13", "OS3",
  "T_19", "OS4",
  "T_21", "OS5",
  "T_25", "OS6",
  "T_29", "OS7",
  "T_39", "OS8",
  "T_07", "YV1",
  "T_12", "YV2",
  "T_18", "YV3",
  "T_20", "YV4",
  "T_27", "YV5",
  "T_30", "YV6",
  "T_36", "YV7",
  "T_42", "YV8"
)



