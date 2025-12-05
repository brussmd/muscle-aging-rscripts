#muliomics integration

#/////////////////////////////////////////////////////////#
#----Transcriptomics--------------------------------------#
#/////////////////////////////////////////////////////////#
oldsedveh_vs_yngsedveh_counts <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
oldpwrveh_vs_yngsedveh_counts <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx")
oldpwrirap_vs_yngsedveh_counts <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_IRAP-YNG_SED_VEH.xlsx")
oldpwrfrap_vs_yngsedveh_counts <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_FRAP-YNG_SED_VEH.xlsx")

yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwr_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwrirap_samples <- c("T_03", "T_14", "T_31", "T_33", "T_37", "T_47", "T_48")
old_pwrfrap_samples <- c("T_01", "T_04", "T_22", "T_24", "T_28", "T_34", "T_40", "T_44")

head(oldsedveh_vs_yngsedveh_counts)
head(oldpwrveh_vs_yngsedveh_counts)
head(oldpwrirap_vs_yngsedveh_counts)
head(oldpwrfrap_vs_yngsedveh_counts)


# rename gene column
colnames(oldsedveh_vs_yngsedveh_counts)[1] <- "Ensembl"
colnames(oldpwrveh_vs_yngsedveh_counts)[1] <- "Ensembl"
colnames(oldpwrirap_vs_yngsedveh_counts)[1] <- "Ensembl"
colnames(oldpwrfrap_vs_yngsedveh_counts)[1] <- "Ensembl"

# merge: keep one set of yng_sed counts, add unique columns from each group
merged_counts <- oldsedveh_vs_yngsedveh_counts %>%
  dplyr::select(Ensembl, all_of(yng_sed_samples), all_of(old_sed_samples)) %>%
  left_join(oldpwrveh_vs_yngsedveh_counts %>% dplyr::select(Ensembl, all_of(setdiff(old_pwr_samples, yng_sed_samples))),
            by="Ensembl") %>%
  left_join(oldpwrirap_vs_yngsedveh_counts %>% dplyr::select(Ensembl, all_of(setdiff(old_pwrirap_samples, yng_sed_samples))),
            by="Ensembl") %>%
  left_join(oldpwrfrap_vs_yngsedveh_counts %>% dplyr::select(Ensembl, all_of(setdiff(old_pwrfrap_samples, yng_sed_samples))),
            by="Ensembl")

head(merged_counts)

merged_counts <- merged_counts %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db,
                           keys = Ensembl_noDec,
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")) %>%
  tidyr::drop_na(ENTREZID)

head(merged_counts)

counts_matrix <- merged_counts %>%
  dplyr::select(-Ensembl, -Ensembl_noDec, -ENTREZID) %>%
  replace(is.na(.), 0) %>%   # <--- convert NA to 0
  as.matrix()
rownames(counts_matrix) <- merged_counts$ENTREZID

# ----------------------------------------------------
# 2. Define groups (update if your sample design changes)
# ----------------------------------------------------
group <- factor(c(
  rep("YNG_SED",      length(yng_sed_samples)),
  rep("OLD_SED",      length(old_sed_samples)),
  rep("OLD_PWR",      length(old_pwr_samples)),
  rep("OLD_PWR_IRAP", length(old_pwrirap_samples)),
  rep("OLD_PWR_FRAP", length(old_pwrfrap_samples))
))

# ----------------------------------------------------
# 3. Build DGEList & normalize
# ----------------------------------------------------
dge <- DGEList(counts = counts_matrix, group = group)
dge <- calcNormFactors(dge, method = "TMM")

# ----------------------------------------------------
# 4. Filter lowly expressed genes
# ----------------------------------------------------
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes kept after filtering:", sum(keep), "\n")

# ----------------------------------------------------
# 5. Calculate logCPM matrix
# ----------------------------------------------------
logCPM_matrix <- cpm(dge, log = TRUE, prior.count = 1)

# ----------------------------------------------------
# 6. Quick checks
# ----------------------------------------------------
cat("Final matrix dimensions: ", dim(logCPM_matrix)[1], " genes × ", dim(logCPM_matrix)[2], " samples\n")
anyNA(logCPM_matrix)   # Should be FALSE
anyDuplicated(colnames(logCPM_matrix)) # Should be 0

head(logCPM_matrix)

# Map ENTREZ IDs to SYMBOL
entrez_ids <- rownames(logCPM_matrix)

symbols <- mapIds(org.Mm.eg.db,
                  keys = entrez_ids,
                  column = "SYMBOL",
                  keytype = "ENTREZID",
                  multiVals = "first")

# Replace row names with SYMBOL, keeping ENTREZID if SYMBOL is missing
symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
rownames(logCPM_matrix) <- symbols

logCPM_df <- as.data.frame(logCPM_matrix, row.names = NULL)
logCPM_df$Symbol <- rownames(logCPM_matrix)

logCPM_symbol <- logCPM_df %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()

rownames(logCPM_symbol) <- logCPM_symbol$Symbol
logCPM_symbol$Symbol <- NULL

head(logCPM_df)

sum(duplicated(rownames(logCPM_matrix)))
anyDuplicated(rownames(logCPM_symbol))  # should return 0
dim(logCPM_symbol) 

saveRDS(logCPM_symbol, file = "logCPM_symbol_final.RDS")

head(logCPM_symbol)
dim(logCPM_symbol)
#--------------------------------------------------#

#/////////////////////////////////////////////////////////#
#----Metabolomics--------------------------------------#
#/////////////////////////////////////////////////////////#

# Read the CSV file (adjust the path as needed)
pwr_rapa_muscle_metabo_data <- read_csv("Konopka_Muscle_HILIC.csv")

# View the first few rows
head(pwr_rapa_muscle_metabo_data)

# 1. Drop Group (keep separately for later)
metabo_wide <- pwr_rapa_muscle_metabo_data %>%
  dplyr::select(-Group)

# 2. Pivot longer so each row is (Sample, Metabolite, Abundance)
metabo_long <- metabo_wide %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance")

# 3. Pivot back to wide so metabolites are rows and samples are columns
metabo_matrix <- metabo_long %>%
  pivot_wider(names_from = Sample, values_from = Abundance)

# 4. Set metabolite names as rownames and convert to matrix
metabo_matrix <- as.data.frame(metabo_matrix)
rownames(metabo_matrix) <- metabo_matrix$Metabolite
metabo_matrix <- metabo_matrix %>% dplyr::select(-Metabolite) %>% as.matrix()

head(metabo_matrix)

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
  "T_42", "YV8",
  "T_05", "YS1",
  "T_15", "YS2",
  "T_17", "YS3",
  "T_23", "YS4",
  "T_26", "YS5",
  "T_32", "YS6",
  "T_35", "YS7",
  "T_38", "YS8"
)

# Convert column names using the mapping table
name_map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)

# Rename columns
colnames(metabo_matrix) <- dplyr::recode(colnames(metabo_matrix), !!!name_map)

head(colnames(metabo_matrix))

colnames(metabo_matrix)
colnames(logCPM_symbol)

# Combine your vectors into one vector of samples in desired order
desired_order <- c(
  yng_sed_samples,
  old_sed_samples,
  old_pwr_samples,
  old_pwrirap_samples,
  old_pwrfrap_samples
)

# Subset metabo_matrix to only include matching columns, in the correct order
metabo_matrix_subset <- metabo_matrix[, desired_order, drop = FALSE]

# Quick check: columns should match exactly
identical(colnames(metabo_matrix_subset), colnames(logCPM_symbol))

head(metabo_matrix_subset)


# 2. Replace NAs and zeros with half minimum non-zero value for each metabolite
metabo_imputed <- t(apply(metabo_matrix_subset, 1, function(x) {
  x[is.na(x)] <- 0  # treat NA as zero
  min_val <- min(x[x > 0], na.rm = TRUE)
  if (is.infinite(min_val)) min_val <- 1e-6 # in case all values are zero
  x[x == 0] <- min_val / 2
  return(x)
}))

# 3. Log transform (log2 recommended, can also use log10)
metabo_log <- log2(metabo_imputed)

# Optional: Convert back to data frame
metabo_log <- as.data.frame(metabo_log)
rownames(metabo_log) <- rownames(metabo_matrix_subset)
colnames(metabo_log) <- desired_order

head(metabo_log)
dim(metabo_log)
#-------------------------------------#

#/////////////////////////////////////////////////////////#
#----Lipidomics--------------------------------------#
#/////////////////////////////////////////////////////////#

# Remove outlier
lipid_data <- read_csv("Konopka_muscle_lipids.csv") %>%
  filter(Sample != "YS6")

# Identify lipid columns
lipid_cols <- setdiff(colnames(lipid_data), c("Sample", "Group"))

# Impute zeros
lipid_imputed <- lipid_data
for (lip in lipid_cols) {
  x <- lipid_imputed[[lip]]
  x[is.na(x)] <- 0
  min_val <- min(x[x > 0], na.rm = TRUE)
  if (is.infinite(min_val)) min_val <- 1e-6
  x[x == 0] <- min_val / 2
  lipid_imputed[[lip]] <- x
}

# Log2 transform
lipid_imputed[, lipid_cols] <- log2(lipid_imputed[, lipid_cols])

# Create matrix: lipids as rows, samples as columns
lipid_matrix <- t(as.matrix(lipid_imputed[, lipid_cols]))
colnames(lipid_matrix) <- lipid_imputed$Sample

head(lipid_matrix)

# --- 1. Remove YV samples ---
lipid_matrix <- lipid_matrix[, !grepl("^YV", colnames(lipid_matrix))]

# --- 2. Rename samples using tube code key ---
# Assuming `sample_key` has columns: TubeCode (e.g., "T_05") and MetabolomicsName (e.g., "YS1")
name_map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
colnames(lipid_matrix) <- dplyr::recode(colnames(lipid_matrix), !!!name_map)

# --- 3. Subset & reorder to match transcriptomics ---
desired_order <- c(
  yng_sed_samples,
  old_sed_samples,
  old_pwr_samples,
  old_pwrirap_samples,
  old_pwrfrap_samples
)

#for lipidomics heatmap
desired_order2 <- c(
  old_sed_samples,
  old_pwr_samples,
  old_pwrirap_samples,
  old_pwrfrap_samples
)

desired_order2_lipids <- desired_order2[desired_order2 != "T_32"]

lipid_matrix <- lipid_matrix[, desired_order2_lipids, drop = FALSE]

# Quick check again
stopifnot(identical(colnames(lipid_matrix), desired_order2_lipids))

# --- Quick checks ---
stopifnot(identical(colnames(lipid_matrix), desired_order)) # should be TRUE
anyNA(lipid_matrix)   # should be FALSE

setdiff(desired_order, colnames(lipid_matrix))
setdiff(colnames(lipid_matrix), desired_order)

head(lipid_matrix)
dim(lipid_matrix)
head(metabo_log)
dim(metabo_log)
head(logCPM_symbol)
dim(logCPM_symbol)
#---------------------------------------------#

#/////////////////////////////////////////////////////////#
#----Combine/Stack all data--------------------------------#
#/////////////////////////////////////////////////////////#

# --- 1. Identify common samples across all omics ---
common_samples <- Reduce(intersect, list(
  colnames(logCPM_symbol),
  colnames(metabo_log),
  colnames(lipid_matrix)
))
cat("Common samples across all omics:", length(common_samples), "\n")

# --- 2. Subset each dataset to these common samples ---
transcriptomics_common <- logCPM_symbol[, common_samples, drop = FALSE]
metabolomics_common    <- metabo_log[, common_samples, drop = FALSE]
lipidomics_common      <- lipid_matrix[, common_samples, drop = FALSE]

# --- 3. Add prefix to rownames to avoid collisions (optional but recommended) ---
rownames(transcriptomics_common) <- paste0("RNA_", rownames(transcriptomics_common))
rownames(metabolomics_common)    <- paste0("Met_", rownames(metabolomics_common))
rownames(lipidomics_common)      <- paste0("Lip_", rownames(lipidomics_common))

# --- 4. Stack by rows ---
integrated_matrix <- rbind(
  transcriptomics_common,
  metabolomics_common,
  lipidomics_common
)

# --- 5. Quick checks ---
cat("Final integrated dimensions: ", dim(integrated_matrix)[1], " features × ", dim(integrated_matrix)[2], " samples\n")
anyNA(integrated_matrix)

saveRDS(integrated_matrix, "integrated_matrix_all_omics_v2.RDS")
head(integrated_matrix)
dim(integrated_matrix)

#--------------------------------------------------------------#

#/////////////////////////////////////////////////////////////#
#-----------GO Core Enrichment Genes (498)--------------------#
#////////////////////////////////////////////////////////////#

library(org.Mm.eg.db)

# Convert ENTREZ IDs to gene symbols
core_immune_symbols <- mapIds(
  org.Mm.eg.db,
  keys = as.character(upreg_old_sed_core_entrez_vec_12jun),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Remove any NA values
core_immune_symbols <- na.omit(core_immune_symbols)

# Check length
length(core_immune_symbols)

# Save to CSV for future use
write.csv(core_immune_symbols, "core_enrichment_symbols.csv", row.names = FALSE)
#----------------------------------------------------------------------#

# Copy transcriptomics data
rna_matrix <- logCPM_symbol[, common_samples]

# Add prefix for integrated version
rownames(rna_matrix) <- paste0("RNA_", rownames(rna_matrix))

# Make a vector of prefixed gene names
core_prefixed <- paste0("RNA_", core_immune_symbols)

# Subset
rna_core <- rna_matrix[rownames(rna_matrix) %in% core_prefixed, ]

# Check dimensions
dim(rna_core)  

integrated_core <- rbind(
  rna_core,
  metabolomics_common,
  lipidomics_common
)
dim(integrated_core)  # should now be (498 + 411 + 248) × 36
#------------------------------------------------------------#

#-------Get Significant Lipids------------------------------#

sig_lipid_features <- aging_lipid_stats_woYS6 %>%
  filter(FDR <0.2 & Log2_FC_OS_vs_YS > 0.4) %>%
  arrange(Lipid)

# Save to CSV for future use
write.csv(sig_lipid_features, "sig_lipid_features.csv", row.names = FALSE)

sig_lipid_features.vec <- sig_lipid_features %>%
  pull(Lipid)

sig_lipid_features.vec

#--------------------------------------------------------------------#

#-----------Get Significant Metabolites------------------------------#

sig_aging_metabolites <-aging_stats %>% 
  filter(p_value <0.05 & log2FC >0.3)

aging_metabolite.vec <-sig_aging_metabolites %>%
  pull(Metabolite)

# Save to CSV for future use
write.csv(aging_metabolite.vec, "aging_metabolite.vec", row.names = FALSE)

aging_metabolite.vec
#---------------------------------------------------------------#

#------Create Integrated Core Omics Set------------------------#

# --- 1. Subset RNA (already immune core selected earlier) ---
rna_core <- rna_matrix[rownames(rna_matrix) %in% paste0("RNA_", core_immune_symbols), ]

# --- 2. Subset Lipids (prefix Lip_) ---
lipid_core <- lipid_matrix[rownames(lipid_matrix) %in% paste0("Lip_", sig_lipid_features.vec), ]

# --- 3. Subset Metabolites (prefix Met_) ---
# Ensure metabolite matrix has same sample columns as RNA/Lipid
metabolite_core <- metabo_log[rownames(metabo_log) %in% paste0("Met_", aging_metabolite.vec), 
                              colnames(rna_core)]

# --- 4. Combine (ensure sample order is identical across blocks) ---
stopifnot(identical(colnames(rna_core), colnames(lipid_core)))
stopifnot(identical(colnames(rna_core), colnames(metabolite_core)))

# No prefixes
rownames(metabo_log) <- gsub("^Met_", "", rownames(metabo_log))
rownames(lipid_matrix) <- gsub("^Lip_", "", rownames(lipid_matrix))

# Subset using raw vectors
metabolite_core <- metabo_log[rownames(metabo_log) %in% aging_metabolite.vec, colnames(rna_core)]
lipid_core <- lipid_matrix[rownames(lipid_matrix) %in% sig_lipid_features.vec, colnames(rna_core)]

integrated_core <- rbind(rna_core, metabolite_core, lipid_core)
cat("Final integrated core dimensions:", dim(integrated_core)[1], "features ×", dim(integrated_core)[2], "samples\n")

# --- 6. Save for future use ---
write.csv(integrated_core, file = "integrated_core_omics.csv")

head(integrated_core)
dim(integrated_core)

#----------------------------------------------------------------#

#-------Inspecting data set-----------------------------------#

# Phenotype vector for PCA coloring
sample_groups <- ifelse(colnames(integrated_core) %in% yng_sed_samples, "YNG_SED",
                        ifelse(colnames(integrated_core) %in% old_sed_samples, "OLD_SED",
                               ifelse(colnames(integrated_core) %in% old_pwr_samples, "OLD_PWR",
                                      ifelse(colnames(integrated_core) %in% old_pwrirap_samples, "OLD_PWR_IRAP",
                                             ifelse(colnames(integrated_core) %in% old_pwrfrap_samples, "OLD_PWR_FRAP", "UNKNOWN")))))


# PCA on samples
pca_res <- prcomp(t(integrated_core), center = TRUE, scale. = TRUE)

# Variance explained for labeling
var_expl <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# Plot with ggplot2
library(ggplot2)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = sample_groups
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       title = "PCA of integrated core omics features") +
  theme_minimal()
#--------------------------------------------------------#

#//////////////////////////////////////////////////////////////#
#-------------DIABLO--------------------------------------------#
#//////////////////////////////////////////////////////////////#

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics")
library(mixOmics)

# Remove T_32 from your group vectors
yng_sed_samples <- setdiff(c("T_05","T_17","T_26","T_32","T_35","T_38"), "T_32")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")

# Keep only these samples for DIABLO
keep_samples <- c(yng_sed_samples, old_sed_samples)

# Create group labels
group_labels <- ifelse(keep_samples %in% yng_sed_samples, "YNG_SED", "OLD_SED")

rna_diablo  <- rna_core[, keep_samples]
met_diablo  <- metabolite_core[, keep_samples]
lipid_diablo <- lipid_core[, keep_samples]

# --- 2. Create sample group factor ---
group_labels <- factor(
  ifelse(keep_samples %in% yng_sed_samples, "YNG_SED",
         ifelse(keep_samples %in% old_sed_samples, "OLD_SED", "OLD_PWR")),
  levels = c("YNG_SED", "OLD_SED", "OLD_PWR")
)

# --- 3. Build block list ---
omics_list <- list(
  RNA = t(rna_diablo),          # samples as rows
  Metabolite = t(met_diablo),
  Lipid = t(lipid_diablo)
)

# --- 4. Create design matrix (correlation strength between blocks) ---
design <- matrix(1, ncol = length(omics_list), nrow = length(omics_list),
                 dimnames = list(names(omics_list), names(omics_list)))
diag(design) <- 0  # no self-correlation


# Run DIABLO (keep all features initially)
diablo_model <- block.splsda(
  X = omics_list,
  Y = group_labels,
  ncomp = 2,          # two components for visualization
  design = design
)

# Visualize samples on the first two components
plotIndiv(diablo_model, legend = TRUE, title = "DIABLO: YNG_SED vs OLD_SED")

# Shows which features from each omics block correlate across datasets
circosPlot(diablo_model, comp = 1, cutoff = 0.7)
network(diablo_model, comp = 1, threshold = 0.3)

# Features selected for each block and component
selected_features <- selectVar(diablo_model, comp = 1)
selected_features$RNA     # top RNA features
selected_features$Metabolite
selected_features$Lipid

# define how many features to keep per block (tune as desired)
list.keepX <- list(
  RNA = length(selected_features$RNA$name),
  Metabolite = length(selected_features$Metabolite$name),
  Lipid = length(selected_features$Lipid$name)
)

diablo_selected <- block.splsda(
  X = omics_list,
  Y = group_labels,
  ncomp = 2,
  design = design,
  keepX = list.keepX
)

# visualize again
plotIndiv(diablo_selected, legend = TRUE, title = "DIABLO Selected Signature")
circosPlot(diablo_selected, comp = 1, cutoff = 0.7)

# PWR sample matrix (subset them from rna_core, metabolite_core, lipid_core)
pwr_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")  # adjust if needed

# PWR data in same structure as training
pwr_data <- list(
  RNA = t(rna_core[, pwr_samples]),
  Metabolite = t(metabolite_core[, pwr_samples]),
  Lipid = t(lipid_core[, pwr_samples])
)

# Predict using DIABLO model
pwr_pred <- predict(diablo_selected, pwr_data)  # no 'newdata='

# Extract scores (variates) for component 1 & 2
pwr_scores <- lapply(pwr_pred$variates, function(x) x[, 1:2])
#-----------------------------------------------------------------#

#-----Good stuff maybe------------------------------#

head(combined_blocks)

# --- 1. Combine all samples (Young, Old, PWR) ---
all_samples <- c(keep_samples, pwr_samples)
all_labels  <- c(as.character(group_labels), rep("PWR", length(pwr_samples)))

combined_blocks <- list(
  RNA = as.matrix(t(rna_core[rna_feats, all_samples, drop = FALSE])),
  Metabolite = as.matrix(t(metabolite_core[met_feats, all_samples, drop = FALSE])),
  Lipid = as.matrix(t(lipid_core[lipid_feats, all_samples, drop = FALSE]))
)

# --- 3. Run DIABLO on all samples (labels still only train supervised separation) ---
diablo_all <- block.splsda(
  X = combined_blocks,
  Y = factor(all_labels, levels = c("YNG_SED", "OLD_SED", "PWR")),
  ncomp = 2,
  design = design,
  keepX = list.keepX
)

# --- 4. Visualize ---
plotIndiv(diablo_all, legend = TRUE, title = "Aging Signature + PWR Projection")











# --- 1. Run DIABLO with feature selection ---
list.keepX <- list(RNA = 30, Metabolite = 10, Lipid = 10)

diablo_signature <- block.splsda(
  X = omics_list,
  Y = group_labels,
  ncomp = 2,
  design = design,
  keepX = list.keepX
)

# --- 2. Visualize sample separation ---
plotIndiv(diablo_signature, legend = TRUE, title = "Aging Signature (DIABLO)")

# --- 3. Cross-omics relationships ---
circosPlot(diablo_signature, comp = 1, cutoff = 0.7)
# plot network for component 1 of all blocks
network(diablo_signature, comp = list(1, 1))

# --- 4. Extract feature list (your “aging signature”) ---
aging_signature <- lapply(names(omics_list), function(block) {
  feats <- selectVar(diablo_signature, comp = 1)[[block]]$name
  data.frame(Block = block, Feature = feats)
})
aging_signature <- do.call(rbind, aging_signature)

# Save to CSV
write.csv(aging_signature, "aging_signature_diablo.csv", row.names = FALSE)


#//////////////////////////////////////////////////////////////#
#-------------MOFA2--------------------------------------------#
#//////////////////////////////////////////////////////////////#

# Install MOFA2 (if not installed)
if (!requireNamespace("MOFA2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("MOFA2")
}
library(MOFA2)
#------------------------#

# --- 1. Create a list of matrices ---
omics_list <- list(
  RNA = as.matrix(transcriptomics_common),
  Metabolite = as.matrix(metabolomics_common),
  Lipid = as.matrix(lipidomics_common)
)

# Replace any NA values with 0
omics_list <- lapply(omics_list, function(x) {
  x[is.na(x)] <- 0
  x
})

# --- 2. Create MOFA object ---
MOFAobject <- create_mofa(omics_list)

# --- 3. Define options ---
data_opts  <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5   # Adjust as needed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"

# --- 4. Prepare the object ---
MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# --- 5. Run MOFA with basilisk-managed Python ---
MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)

# --- 6. Save trained model (optional) ---
saveRDS(MOFAobject, file = "MOFAobject_trained.rds")

# --- 7. Visualize variance explained ---
plot_variance_explained(MOFAobject)

# --- 8. Sample clustering in factor space ---
plot_factors(MOFAobject, factors = 1:2, color_by = "Group")

# --- 9. Feature weights (example for first factor, RNA view) ---
plot_weights(MOFAobject, view = "RNA", factor = 1, nfeatures = 20)

# 1. Pull existing sample metadata from the MOFA object
current_metadata <- samples_metadata(MOFAobject)

# 2. Check what it looks like (just to confirm structure)
head(current_metadata)

# Use your predefined vectors to classify samples
sample_groups <- ifelse(current_metadata$sample %in% yng_sed_samples, "YNG_SED",
                        ifelse(current_metadata$sample %in% old_sed_samples, "OLD_SED",
                               ifelse(current_metadata$sample %in% old_pwr_samples, "OLD_PWR",
                                      ifelse(current_metadata$sample %in% old_pwrirap_samples, "OLD_PWR_IRAP",
                                             ifelse(current_metadata$sample %in% old_pwrfrap_samples, "OLD_PWR_FRAP", "UNKNOWN")))))

# Add this as a new column
current_metadata$Group <- factor(sample_groups,
                                 levels = c("YNG_SED", "OLD_SED", "OLD_PWR",
                                            "OLD_PWR_IRAP", "OLD_PWR_FRAP"))

# Re-assign metadata
samples_metadata(MOFAobject) <- current_metadata

# Plot factors colored by the new Group labels
plot_factors(MOFAobject, factors = 1:2, color_by = "Group")

# Remove Factor 1 from trained MOFA
MOFAobject_filtered <- remove_factors(MOFAobject, factors = 1)

#======================================================================
# 1. Define sample groups
#======================================================================
yng_sed_samples <- c("T_05","T_17","T_26","T_35","T_38")  # T_32 removed
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
pwr_samples     <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")

keep_samples <- c(yng_sed_samples, old_sed_samples)
group_labels <- factor(ifelse(keep_samples %in% yng_sed_samples, "YNG_SED", "OLD_SED"))

#======================================================================
# 2. Build omics blocks (features = columns, samples = rows)
#======================================================================
omics_list <- list(
  RNA        = t(rna_core[, keep_samples, drop = FALSE]),        # samples x features
  Metabolite = t(metabolite_core[, keep_samples, drop = FALSE]), # samples x features
  Lipid      = t(lipid_core[, keep_samples, drop = FALSE])       # samples x features
)

# Check numeric and no NA
stopifnot(all(sapply(omics_list, function(x) is.numeric(x) && !anyNA(x))))

# Design matrix (link all blocks equally)
design <- matrix(1, ncol = length(omics_list), nrow = length(omics_list),
                 dimnames = list(names(omics_list), names(omics_list)))
diag(design) <- 0

#======================================================================
# 3. Train DIABLO (feature selection)
#======================================================================
list.keepX <- list(RNA = 30, Metabolite = 10, Lipid = 10)

diablo_signature <- block.splsda(
  X = omics_list,
  Y = group_labels,
  ncomp = 2,
  design = design,
  keepX = list.keepX
)

# Plots
plotIndiv(diablo_signature, legend = TRUE, title = "Aging Signature (DIABLO)")
circosPlot(diablo_signature, comp = 1, cutoff = 0.7)

# Extract feature list
aging_signature <- do.call(rbind, lapply(names(omics_list), function(block) {
  feats <- selectVar(diablo_signature, comp = 1)[[block]]$name
  data.frame(Block = block, Feature = feats)
}))
write.csv(aging_signature, "aging_signature_diablo.csv", row.names = FALSE)

#======================================================================
# 4. Prepare PWR data (subset to selected features)
#======================================================================
sel_feats <- selectVar(diablo_signature, comp = 1)

pwr_data <- list(
  RNA        = t(rna_core[sel_feats$RNA$name,        pwr_samples, drop = FALSE]),
  Metabolite = t(metabolite_core[sel_feats$Metabolite$name, pwr_samples, drop = FALSE]),
  Lipid      = t(lipid_core[sel_feats$Lipid$name,    pwr_samples, drop = FALSE])
)

# Check sample alignment
stopifnot(identical(rownames(omics_list$RNA), keep_samples))
stopifnot(identical(rownames(pwr_data$RNA), pwr_samples))

#======================================================================
# 5. Combine training + PWR (for visualization)
#======================================================================
all_samples <- c(keep_samples, pwr_samples)
all_labels  <- factor(c(as.character(group_labels), rep("PWR", length(pwr_samples))),
                      levels = c("YNG_SED", "OLD_SED", "PWR"))

combined_blocks <- list(
  RNA        = rbind(omics_list$RNA[, sel_feats$RNA$name, drop = FALSE],
                     pwr_data$RNA[, sel_feats$RNA$name, drop = FALSE]),
  Metabolite = rbind(omics_list$Metabolite[, sel_feats$Metabolite$name, drop = FALSE],
                     pwr_data$Metabolite[, sel_feats$Metabolite$name, drop = FALSE]),
  Lipid      = rbind(omics_list$Lipid[, sel_feats$Lipid$name, drop = FALSE],
                     pwr_data$Lipid[, sel_feats$Lipid$name, drop = FALSE])
)

#======================================================================
# 6. Final DIABLO including PWR projection
#======================================================================
diablo_all <- block.splsda(
  X = combined_blocks,
  Y = all_labels,
  ncomp = 2,
  design = design,
  keepX = list.keepX
)

# Visualization
plotIndiv(diablo_all, legend = TRUE, title = "Aging Signature + PWR Projection")
#---------------------------------------------------#

# --- 1. Extract scores (global integrated space) ---
scores <- diablo_all$variates$RNA[,1:2]   # comp1 & comp2
scores_df <- data.frame(scores, Group = diablo_all$Y)

# --- 2. Compute centroids ---
centroids <- aggregate(scores_df[,1:2], by = list(Group=scores_df$Group), FUN = mean)
rownames(centroids) <- centroids$Group
centroids <- centroids[, -1]  # remove label col

# --- 3. Euclidean distances ---
dist_YNG_OLD <- sqrt(sum((centroids["YNG_SED",] - centroids["OLD_SED",])^2))
dist_YNG_PWR <- sqrt(sum((centroids["YNG_SED",] - centroids["PWR",])^2))
cat("Distance YNG vs OLD:", dist_YNG_OLD, "\n")
cat("Distance YNG vs PWR:", dist_YNG_PWR, "\n")

# --- 4. Permutation test ---
set.seed(123)
nperm <- 10000
perm_diff <- numeric(nperm)

for(i in 1:nperm){
  perm_labels <- sample(scores_df$Group)
  perm_centroids <- aggregate(scores_df[,1:2], by = list(Group=perm_labels), FUN = mean)
  rownames(perm_centroids) <- perm_centroids$Group
  perm_centroids <- perm_centroids[, -1]
  if(all(c("YNG_SED","OLD_SED","PWR") %in% rownames(perm_centroids))){
    perm_dist_YNG_OLD <- sqrt(sum((perm_centroids["YNG_SED",] - perm_centroids["OLD_SED",])^2))
    perm_dist_YNG_PWR <- sqrt(sum((perm_centroids["YNG_SED",] - perm_centroids["PWR",])^2))
    perm_diff[i] <- perm_dist_YNG_OLD - perm_dist_YNG_PWR
  } else {
    perm_diff[i] <- NA
  }
}
perm_diff <- na.omit(perm_diff)
obs_diff <- dist_YNG_OLD - dist_YNG_PWR
pval <- mean(perm_diff <= obs_diff)
cat("Observed distance difference (OLD - PWR):", obs_diff, "\n")
cat("Permutation p-value:", pval, "\n")

# --- 5. Visualization ---
library(ggplot2)
ggplot(scores_df, aes(x = comp1, y = comp2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_point(data = data.frame(centroids, Group=rownames(centroids)),
             aes(x=comp1, y=comp2, color=Group), size=5, shape=4, stroke=2) +
  geom_segment(aes(x=centroids["YNG_SED","comp1"], y=centroids["YNG_SED","comp2"],
                   xend=centroids["OLD_SED","comp1"], yend=centroids["OLD_SED","comp2"]),
               arrow=arrow(length=unit(0.3,"cm")), color="red", size=1) +
  geom_segment(aes(x=centroids["YNG_SED","comp1"], y=centroids["YNG_SED","comp2"],
                   xend=centroids["PWR","comp1"], yend=centroids["PWR","comp2"]),
               arrow=arrow(length=unit(0.3,"cm")), color="blue", size=1) +
  annotate("text", x=mean(c(centroids["YNG_SED","comp1"],centroids["OLD_SED","comp1"])),
           y=mean(c(centroids["YNG_SED","comp2"],centroids["OLD_SED","comp2"])),
           label=paste0("Dist=",round(dist_YNG_OLD,2)), color="red", vjust=-1) +
  annotate("text", x=mean(c(centroids["YNG_SED","comp1"],centroids["PWR","comp1"])),
           y=mean(c(centroids["YNG_SED","comp2"],centroids["PWR","comp2"])),
           label=paste0("Dist=",round(dist_YNG_PWR,2)), color="blue", vjust=1.5) +
  theme_minimal() +
  labs(title="Integrated multi-omics latent space\n(PWR closer to YNG than OLD)",
       x="Component 1", y="Component 2")
#---------------------------------------------------------#

#////////////////////////////////////////////////////////#
#-----Composite Aging Index-----------------------------#
#///////////////////////////////////////////////////////#

old_pwr_irap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwr_frap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

validated_inflamm_geneset.vec
aging_metabolite.vec
sig_lipid_features.vec

# ---------------------------
# 1. Subset integrated dataset
# ---------------------------

# Remove "RNA_" prefix for genes
rownames(integrated_core) <- gsub("^RNA_", "", rownames(integrated_core))

# Combine selected features
selected_features <- c(validated_inflamm_geneset.vec,
                       aging_metabolite.vec,
                       sig_lipid_features.vec)

# Subset integrated_core to these features
aging_subset <- integrated_core[rownames(integrated_core) %in% selected_features, ]

cat("Selected features:", nrow(aging_subset), "\n")  # Expect ~55+17+34
aging_subset
# ---------------------------
# 2. Z-score normalize each feature
# ---------------------------
aging_z <- t(scale(t(aging_subset)))  # z-score by row (feature)
stopifnot(!anyNA(aging_z))

# ---------------------------
# 3. Composite aging index per sample
# ---------------------------
aging_index <- colMeans(aging_z)  # mean z-score for each sample
aging_index_df <- data.frame(Sample = names(aging_index),
                             AgingIndex = aging_index,
                             Group = ifelse(names(aging_index) %in% yng_sed_samples, "YNG_SED",
                                            ifelse(names(aging_index) %in% old_sed_samples, "OLD_SED",
                                                   ifelse(names(aging_index) %in% pwr_samples, "OLD_PWR", "Other"))))
aging_index_df
aging_z
# ---------------------------
# 4. Statistical tests
# ---------------------------
library(dplyr)
library(ggplot2)

# Test OLD vs YNG (baseline aging effect)
wilcox_old_yng <- wilcox.test(AgingIndex ~ Group, data = aging_index_df %>% 
                                filter(Group %in% c("YNG_SED", "OLD_SED")))
cat("OLD vs YNG p-value:", wilcox_old_yng$p.value, "\n")

# Test if OLD_PWR differs from OLD_SED and moves toward YNG
wilcox_pwr_old <- wilcox.test(AgingIndex ~ Group, data = aging_index_df %>% 
                                filter(Group %in% c("OLD_SED", "OLD_PWR")))
cat("PWR vs OLD p-value:", wilcox_pwr_old$p.value, "\n")

# ---------------------------
# 5. Visualization: Aging Index boxplot
# ---------------------------
ggplot(aging_index_df, aes(x = Group, y = AgingIndex, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  theme_minimal() +
  labs(title = "Composite Aging Index",
       subtitle = "Positive values = more aging-like signature",
       y = "Aging Index (mean z-score)", x = "") +
  scale_fill_manual(values = c("YNG_SED" = "#1f78b4", 
                               "OLD_SED" = "#e31a1c",
                               "OLD_PWR" = "#33a02c"))

# ---------------------------
# 6. PCA projection for visualization
# ---------------------------
pca_res <- prcomp(t(aging_z), center = TRUE, scale. = FALSE)
pca_df <- data.frame(pca_res$x[, 1:2],
                     Group = aging_index_df$Group[match(rownames(pca_res$x),
                                                        aging_index_df$Sample)])

var_expl <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95, linetype = 2) +
  theme_minimal() +
  labs(title = "PCA of Aging Signature Features",
       x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"))

anova_model <- aov(AgingIndex ~ Group, data = aging_index_df)
summary(anova_model)  # overall p-value

# Post-hoc pairwise comparisons with adjustment
TukeyHSD(anova_model)


aging_index_df_sub <- subset(aging_index_df, Group %in% c("YNG_SED", "OLD_SED", "OLD_PWR"))

ggboxplot(aging_index_df_sub,
          x = "Group",
          y = "AgingIndex",
          fill = "Group",
          add = "jitter",        # <-- adds individual sample dots
          add.params = list(size = 2, alpha = 0.8),
          palette = c("YNG_SED" = "#1f78b4",
                      "OLD_SED" = "#e31a1c",
                      "OLD_PWR" = "#33a02c")) +
  stat_compare_means(method = "anova", 
                     label.y = max(aging_index_df_sub$AgingIndex) + 0.5) +
  stat_compare_means(comparisons = list(c("YNG_SED","OLD_SED"),
                                        c("OLD_SED","OLD_PWR"),
                                        c("YNG_SED","OLD_PWR")),
                     method = "t.test",
                     label = "p.signif") +
  labs(title = "Composite Aging Index",
       y = "Aging Index (mean z-score)",
       x = "")
#----------------------------------------------------#

#----------weighted----------------------------------#
integrated_core %>%
  filter(rownames(integrated_core) %in% validated_inflamm_geneset.vec)

# --- Subset features by block (already have rownames without RNA_ prefix)
rna_features    <- intersect(validated_inflamm_geneset.vec, rownames(integrated_core))
met_features    <- intersect(aging_metabolite.vec, rownames(integrated_core))
lipid_features  <- intersect(sig_lipid_features.vec, rownames(integrated_core))

# --- Z-score normalize each feature
aging_z <- t(scale(t(integrated_core[c(rna_features, met_features, lipid_features), ])))

# --- Compute per-block means
rna_index   <- colMeans(aging_z[rna_features, , drop = FALSE])
met_index   <- colMeans(aging_z[met_features, , drop = FALSE])
lipid_index <- colMeans(aging_z[lipid_features, , drop = FALSE])

rna_index

# --- Equal block-weight composite score
aging_index <- (rna_index + met_index + lipid_index) / 3

# --- Data frame (keeps old structure, so downstream plots still work)
aging_index_df <- data.frame(
  Sample = names(aging_index),
  AgingIndex = aging_index,
  RNA_Index = rna_index,
  Met_Index = met_index,
  Lipid_Index = lipid_index,
  Group = ifelse(names(aging_index) %in% yng_sed_samples, "YNG_SED",
                 ifelse(names(aging_index) %in% old_sed_samples, "OLD_SED",
                        ifelse(names(aging_index) %in% pwr_samples, "OLD_PWR", "Other")))
)

anova_model <- aov(AgingIndex ~ Group, data = aging_index_df)
summary(anova_model)  # overall p-value

# Post-hoc pairwise comparisons with adjustment
TukeyHSD(anova_model)


aging_index_df_sub <- subset(aging_index_df, Group %in% c("YNG_SED", "OLD_SED", "OLD_PWR"))

ggboxplot(aging_index_df_sub,
          x = "Group",
          y = "AgingIndex",
          fill = "Group",
          add = "jitter",        # <-- adds individual sample dots
          add.params = list(size = 2, alpha = 0.8),
          palette = c("YNG_SED" = "#1f78b4",
                      "OLD_SED" = "#e31a1c",
                      "OLD_PWR" = "#33a02c")) +
  stat_compare_means(method = "anova", 
                     label.y = max(aging_index_df_sub$AgingIndex) + 0.5) +
  stat_compare_means(comparisons = list(c("YNG_SED","OLD_SED"),
                                        c("OLD_SED","OLD_PWR"),
                                        c("YNG_SED","OLD_PWR")),
                     method = "t.test",
                     label = "p.signif") +
  labs(title = "Composite Aging Index",
       y = "Aging Index (mean z-score)",
       x = "")
#----------------------------------------------------------#

# Combine indices into one dataframe (replace your own calculated indices here)
block_index_df <- data.frame(
  Sample = colnames(rna_index_z),
  Group = aging_index_df$Group[match(colnames(rna_index_z), aging_index_df$Sample)],
  RNA_Index = colMeans(rna_index_z),
  Metabolite_Index = colMeans(met_index_z),
  Lipid_Index = colMeans(lipid_index_z),
  Composite_Index = aging_index_df$AgingIndex
)

summary_table <- aging_index_df %>%
  group_by(Group) %>%
  summarise(
    RNA_Index = sprintf("%.2f ± %.2f", mean(RNA_Index), sd(RNA_Index)),
    Met_Index = sprintf("%.2f ± %.2f", mean(Met_Index), sd(Met_Index)),
    Lipid_Index = sprintf("%.2f ± %.2f", mean(Lipid_Index), sd(Lipid_Index)),
    Composite_Index = sprintf("%.2f ± %.2f", mean(AgingIndex), sd(AgingIndex)),
    .groups = "drop"
  )

print(summary_table)

# Convert to long format for plotting (no rename needed)
long_df <- aging_index_df %>%
  dplyr::select(Sample, Group, RNA_Index, Met_Index, Lipid_Index, AgingIndex) %>%
  pivot_longer(cols = c(RNA_Index, Met_Index, Lipid_Index, AgingIndex),
               names_to = "IndexType",
               values_to = "Score")

# Plot all indices side-by-side
ggplot(long_df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  facet_wrap(~ IndexType, scales = "free_y", nrow = 1) +
  theme_minimal() +
  labs(title = "Per-Block and Composite Aging Indices",
       y = "Index Score (z-score mean)",
       x = "") +
  scale_fill_manual(values = c("YNG_SED" = "#1f78b4",
                               "OLD_SED" = "#e31a1c",
                               "OLD_PWR" = "#33a02c",
                               "Other"   = "#6a3d9a"))
#---------------------------------------------------------#

#----Most refined composite index-------------------------#

#----------------------------------------------------------
# 1. Define rapamycin sample IDs
#----------------------------------------------------------
old_pwr_irap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwr_frap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

#----------------------------------------------------------
# 2. Assign all sample group labels including IRAP and FRAP
#----------------------------------------------------------
aging_index_df <- data.frame(
  Sample = names(aging_index),
  AgingIndex = aging_index,
  RNA_Index = rna_index,
  Met_Index = met_index,
  Lipid_Index = lipid_index,
  Group = case_when(
    names(aging_index) %in% yng_sed_samples ~ "YNG_SED",
    names(aging_index) %in% old_sed_samples ~ "OLD_SED",
    names(aging_index) %in% pwr_samples     ~ "OLD_PWR",
    names(aging_index) %in% old_pwr_irap_samples ~ "OLD_PWR_IRAP",
    names(aging_index) %in% old_pwr_frap_samples ~ "OLD_PWR_FRAP",
    TRUE ~ "Other"
  )
)


#----------------------------------------------------------
# 3. ANOVA + Tukey post-hoc for all groups
#----------------------------------------------------------
anova_model <- aov(AgingIndex ~ Group, data = aging_index_df)
summary(anova_model)
TukeyHSD(anova_model)

#----------------------------------------------------------
# 4. Boxplot (all groups visible)
#----------------------------------------------------------
library(ggpubr)
# Define pairwise comparisons (you can trim to only the ones of interest)
comparisons <- list(
  c("YNG_SED", "OLD_SED"),
  c("OLD_SED", "OLD_PWR"),
  c("OLD_SED", "OLD_PWR_IRAP"),
  c("OLD_SED", "OLD_PWR_FRAP")
)

AgingIndex_composite_Boxplot <- ggboxplot(aging_index_df,
          x = "Group",
          y = "AgingIndex",
          fill = "Group",
          add = "jitter",
          add.params = list(size = 2, alpha = 0.8),
          palette = c("YNG_SED" = "#1f78b4",
                      "OLD_SED" = "#e31a1c",
                      "OLD_PWR" = "#33a02c",
                      "OLD_PWR_IRAP" = "#6a3d9a",
                      "OLD_PWR_FRAP" = "#ff7f00")) +
  stat_compare_means(method = "anova",
                     label.y = max(aging_index_df$AgingIndex) + 0.5) +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test",
                     label = "p.signif") +
  labs(title = "Composite Aging Index",
       y = "Aging Index (mean z-score)",
       x = "")

pdf("AgingIndex_composite_Boxplot.pdf", width = 6, height = 5)  # size in inches
print(AgingIndex_composite_Boxplot)
dev.off()

aging_index_df

#----------------------------------------------------------
# 5. Per-block and composite indices (long format)
#----------------------------------------------------------
long_df <- aging_index_df %>%
  dplyr::select(Sample, Group, RNA_Index, Met_Index, Lipid_Index, AgingIndex) %>%
  pivot_longer(cols = c(RNA_Index, Met_Index, Lipid_Index, AgingIndex),
               names_to = "IndexType",
               values_to = "Score")

long_df

ggplot(long_df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  facet_wrap(~ IndexType, scales = "free_y", nrow = 1) +
  theme_minimal() +
  labs(title = "Per-Block and Composite Aging Indices",
       y = "Index Score (z-score mean)",
       x = "") +
  scale_fill_manual(values = c("YNG_SED" = "#1f78b4",
                               "OLD_SED" = "#e31a1c",
                               "OLD_PWR" = "#33a02c",
                               "OLD_PWR_IRAP" = "#6a3d9a",
                               "OLD_PWR_FRAP" = "#ff7f00"))

#-------------------------------------------------------------#
#Test variables
rna_features
aging_index_df
rna_index
integrated_core
#-------------------------------------------------------------#

#=============================================================#
#------Composite Index Heatmaps-------------------------------#
#=============================================================#

aging_subset
aging_index_df
aging_z
rna_features
lipid_features
met_features
rna_index


old_pwr_irap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwr_frap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")
yng_sed_samples <- c("T_05","T_17","T_26","T_35","T_38")  # T_32 removed
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
pwr_samples     <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")


aging_subset %>%
  filter(rownames(aging_subset) %in% rna_features)

library(ComplexHeatmap)
library(circlize)

# 1) Subset & order rows by your gene list
mat <- as.matrix(aging_subset[rna_features, , drop = FALSE])

# 2) Z-score per gene (row) in one line; then guard against NaN/Inf
aging_z <- t(scale(t(mat)))
aging_z[!is.finite(aging_z)] <- 0

# (Optional) cap extremes for nicer colors
aging_z_cap <- pmin(pmax(aging_z, -3), 3)

# 3) Plot heatmap
col_fun <- colorRamp2(c(-3, 0, 3), c("#2c7bb6", "white", "#d7191c"))

# Build the heatmap object
ht <- Heatmap(
  aging_z_cap,
  name = "Z",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "Z-scored expression (genes × samples)"
)

# Render it to the Plots pane
print(ht)

# Generate heatmap
pheatmap(aging_z_cap,
         cluster_rows = FALSE,    # Disable row clustering to preserve order
         cluster_cols = FALSE,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# Pick a tighter symmetric range (tune 1.0–1.5–2.0)
max_abs <- 2.25
cols <- colorRampPalette(c("blue","white","red"))(50)
brks <- seq(-max_abs, max_abs, length.out = length(cols) + 1)

musage_z_map <- pheatmap(
  aging_z_cap,
  scale = "none",                 # already z-scored
  color = cols,
  breaks = brks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE
)

# Save as PDF
pdf("musage_z_map.pdf", width = 4.25, height = 5.5)  # Adjust height
print(musage_z_map)
dev.off()

#-----------------------------------------------#
#-----single row heatmap of RNA column averages----#
# Named vector -> 1-row matrix
mat <- matrix(rna_index, nrow = 1,
              dimnames = list("RNA index", names(rna_index)))

# Palette and breaks (50 colors -> 51 breaks)
max_abs <- 2.25
cols <- colorRampPalette(c("blue","white","red"))(50)
brks <- seq(-max_abs, max_abs, length.out = length(cols) + 1)

musage_col_avg_map <- pheatmap(
  mat,
  scale         = "none",
  color         = cols,
  breaks        = brks,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main          = "Composite RNA Index (z-score)"
)

# Save as PDF
pdf("musage_col_avg_map.pdf", width = 4.25, height = 5.5)  # Adjust height
print(musage_col_avg_map)
dev.off()
#---------------------------------------------#

#--------------------------------------------------#
#----Lipid Aging Signature heatmap for composite---#
#--------------------------------------------------#

aging_subset %>%
  filter(rownames(aging_subset) %in% lipid_features)

library(ComplexHeatmap)
library(circlize)

# 1) Subset & order rows by your gene list
lipid_mat <- as.matrix(aging_subset[lipid_features, , drop = FALSE])

# 2) Z-score per gene (row) in one line; then guard against NaN/Inf
lipid_aging_z <- t(scale(t(lipid_mat)))
lipid_aging_z[!is.finite(lipid_aging_z)] <- 0

# (Optional) cap extremes for nicer colors
lipid_aging_z_cap <- pmin(pmax(lipid_aging_z, -3), 3)

# Pick a tighter symmetric range (tune 1.0–1.5–2.0)
lipid_max_abs <- 2.25
lipid_cols <- colorRampPalette(c("#FF7F00", "black", "#00ffff"))(100)
lipid_brks <- seq(-lipid_max_abs, lipid_max_abs, length.out = length(lipid_cols) + 1)


lipid_z_map <- pheatmap(
  lipid_aging_z_cap,
  scale = "none",                 # already z-scored
  color = lipid_cols,
  breaks = lipid_brks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE
)

# Save as PDF
pdf("lipid_z_map.pdf", width = 4.25, height = 5.5)  # Adjust height
print(lipid_z_map)
dev.off()

#-----------------------------------------------#
#-----single row heatmap of Lipid column averages----#

aging_index_df
# Named vector -> 1-row matrix
lipid_mat_avg <- matrix(lipid_index, nrow = 1,
              dimnames = list("Lipid index", names(lipid_index)))

lipid_mat_avg

# Palette and breaks (50 colors -> 51 breaks)
lipid_max_abs <- 2.25
lipid_cols <- colorRampPalette(c("#FF7F00", "black", "#00ffff"))(100)
lipid_brks <- seq(-lipid_max_abs, lipid_max_abs, length.out = length(lipid_cols) + 1)

lipid_col_avg_map <- pheatmap(
  lipid_mat_avg,
  scale         = "none",
  color         = lipid_cols,
  breaks        = lipid_brks,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main          = "Composite Lipid Index (z-score)"
)

# Save as PDF
pdf("lipid_col_avg_map.pdf", width = 4.25, height = 5.5)  # Adjust height
print(lipid_col_avg_map)
dev.off()
#--------------------------------------------------------------#

#--------------------------------------------------#
#----Metabolite Aging Signature heatmap for composite---#
#--------------------------------------------------#

aging_subset %>%
  filter(rownames(aging_subset) %in% met_features)

library(ComplexHeatmap)
library(circlize)

# 1) Subset & order rows by your gene list
metab_mat <- as.matrix(aging_subset[met_features, , drop = FALSE])

# 2) Z-score per gene (row) in one line; then guard against NaN/Inf
metab_aging_z <- t(scale(t(metab_mat)))
metab_aging_z[!is.finite(metab_aging_z)] <- 0

# (Optional) cap extremes for nicer colors
metab_aging_z_cap <- pmin(pmax(metab_aging_z, -3), 3)

# Pick a tighter symmetric range (tune 1.0–1.5–2.0)
metab_max_abs <- 2.25
metab_brks <- seq(-metab_max_abs, metab_max_abs, length.out = length(metab_cols) + 1)


metab_z_map <- pheatmap(
  metab_aging_z,
  scale = "none",                 # already z-scored
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE
)

# Save as PDF
pdf("metab_z_map.pdf", width = 4.25, height = 5.5)  # Adjust height
print(metab_z_map)
dev.off()
#-------------------------------------------#

aging_mat_avg <- matrix(
  aging_index_df$AgingIndex,
  nrow = 1,
  dimnames = list("Aging index", aging_index_df$Sample)
)

aging_mat_avg


# 1) Set a tight symmetric range so the low group saturates
max_abs <- 1.25  # tune 1.0–1.5; smaller = more contrast among the 4 groups

# 2) Diverging palette with white at 0
aging_cols <- colorRampPalette(c("#0D0887", "#FFFFFF", "#F0F921"))(100)  # blue–white–yellow
aging_cols <- colorRampPalette(c("#003366", "#FFFFFF", "#FF7F0E"))(100)
aging_brks <- seq(-max_abs, max_abs, length.out = length(aging_cols) + 1)

aging_avg_map <- pheatmap(
  aging_mat_avg,
  scale         = "none",
  color         = aging_cols,
  breaks        = aging_brks,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color  = NA,
  main          = "Composite Aging Index (z-score)"
)

pdf("aging_avg_map.pdf", width = 4.25, height = 5.5)
print(aging_avg_map) 
dev.off()
#---------------------------------------#















#=============================================================#
#/////////////////////////////////////////////////////////////#
#--------correlation matrix-----------------------------------#
#/////////////////////////////////////////////////////////////#

if (!requireNamespace("dendextend", quietly = TRUE))
  install.packages("dendextend")

library(cluster)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dplyr)

#--------------------------------------------------------
# 1. Remove global aging effect (regress out aging index)
#--------------------------------------------------------
aging_score
aging_score <- colMeans(aging_z)
aging_residual <- t(apply(aging_z, 1, function(x) {
  lm_fit <- lm(x ~ aging_score)
  residuals(lm_fit)
}))
rownames(aging_residual) <- rownames(aging_z)
colnames(aging_residual) <- colnames(aging_z)

#--------------------------------------------------------
# 2. Compute residual correlation matrix & clustering
#--------------------------------------------------------
resid_cor_matrix <- cor(t(aging_residual), method = "spearman")
hc <- hclust(as.dist(1 - resid_cor_matrix), method = "complete")

# Optimal clusters via silhouette (2–10 clusters)
d <- dist(1 - resid_cor_matrix)
sil_scores <- sapply(2:10, function(k) {
  cluster_assign <- cutree(hc, k)
  sil <- silhouette(cluster_assign, d)
  mean(sil[, 3])
})
optimal_k <- which.max(sil_scores)
cat("Optimal number of clusters based on silhouette:", optimal_k, "\n")
cluster_assign <- cutree(hc, k = optimal_k)

#--------------------------------------------------------
# 3. Order matrix & annotate omics type
#--------------------------------------------------------
feature_order <- hc$order
cor_matrix_ord <- resid_cor_matrix[feature_order, feature_order]
feature_types <- ifelse(rownames(aging_residual) %in% rna_features, "RNA",
                        ifelse(rownames(aging_residual) %in% met_features, "Metabolite", "Lipid"))
feature_types_ord <- feature_types[feature_order]
cluster_assign_ord <- cluster_assign[feature_order]

# Keep only upper triangle
cor_matrix_tri <- cor_matrix_ord
cor_matrix_tri[lower.tri(cor_matrix_tri)] <- NA

# Row annotations
type_colors <- c("RNA"="#1f78b4","Metabolite"="#33a02c","Lipid"="#ff7f00")
row_ha <- rowAnnotation(
  Omics = feature_types_ord,
  Cluster = as.factor(cluster_assign_ord),
  col = list(
    Omics = type_colors,
    Cluster = structure(
      circlize::rand_color(optimal_k), 
      names = as.character(1:optimal_k)
    )
  )
)

# Heatmap colors
col_fun <- colorRamp2(c(-1,0,1), c("blue","white","red"))

#--------------------------------------------------------
# 4. Plot triangular correlation heatmap
#--------------------------------------------------------
Heatmap(
  cor_matrix_tri,
  name = "Residual\nSpearman\nCorr",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  right_annotation = row_ha,
  column_title = paste0("Residual feature-feature correlation (", optimal_k, " clusters)"),
  na_col = "white"
)

#--------------------------------------------------------
# 5. Top feature pairs table (sorted by abs correlation)
#--------------------------------------------------------
cor_df <- as.data.frame(as.table(resid_cor_matrix))
cor_df <- cor_df[cor_df$Var1 < cor_df$Var2, ]  # remove duplicates & self
cor_df <- cor_df %>%
  arrange(desc(abs(Freq)))

head(cor_df, 20)  # top 20 most correlated feature pairs
write.csv(cor_df, "top_residual_feature_pairs.csv", row.names = FALSE)

#--------------------------------------------------------
# 6. Save cluster membership for interpretation
#--------------------------------------------------------
cluster_table <- data.frame(
  Feature = rownames(resid_cor_matrix)[feature_order],
  OmicsType = feature_types_ord,
  Cluster = cluster_assign_ord
)
write.csv(cluster_table, "feature_residual_cluster_assignments.csv", row.names = FALSE)
head(cluster_table)
#-------------------------------------------------------------------------------#

#------Residualized Triangle Heatmap--------------------------------------------#

#--------------------------------------------------------------
# 1. Residualize each feature by Group
#--------------------------------------------------------------
residual_matrix <- t(apply(aging_z, 1, function(x) {
  lm_fit <- lm(x ~ aging_index_df$Group)
  residuals(lm_fit)
}))
rownames(residual_matrix) <- rownames(aging_z)
colnames(residual_matrix) <- colnames(aging_z)

#--------------------------------------------------------------
# 2. Compute Spearman correlation & hierarchical clustering
#--------------------------------------------------------------
cor_resid <- cor(t(residual_matrix), method = "spearman")
hc_resid  <- hclust(as.dist(1 - cor_resid), method = "complete")

#--------------------------------------------------------------
# 3. Optimal number of clusters via silhouette score
#--------------------------------------------------------------
d_resid <- dist(1 - cor_resid)
sil_scores <- sapply(2:10, function(k) {
  cluster_assign <- cutree(hc_resid, k)
  sil <- silhouette(cluster_assign, d_resid)
  mean(sil[, 3])
})
optimal_k_resid <- which.max(sil_scores)
cat("Optimal number of clusters (residualized):", optimal_k_resid, "\n")

# Final cluster assignment
cluster_assign_resid <- cutree(hc_resid, k = optimal_k_resid)

#--------------------------------------------------------------
# 4. Order correlation matrix and create annotations
#--------------------------------------------------------------
feature_order_resid <- hc_resid$order
cor_resid_ord <- cor_resid[feature_order_resid, feature_order_resid]

# Create triangular form
cor_resid_tri <- cor_resid_ord
cor_resid_tri[lower.tri(cor_resid_tri)] <- NA

# Omics type assignment
feature_types_resid <- ifelse(rownames(residual_matrix) %in% rna_features, "RNA",
                              ifelse(rownames(residual_matrix) %in% met_features, "Metabolite", "Lipid"))
feature_types_resid <- factor(feature_types_resid, levels = c("RNA", "Metabolite", "Lipid"))

# Annotation colors
type_colors <- c("RNA" = "#1f78b4", "Metabolite" = "#33a02c", "Lipid" = "#ff7f00")

row_ha_resid <- rowAnnotation(
  Omics = feature_types_resid[feature_order_resid],
  Cluster = as.factor(cluster_assign_resid[feature_order_resid]),
  col = list(
    Omics = type_colors,
    Cluster = structure(circlize::rand_color(optimal_k_resid),
                        names = as.character(1:optimal_k_resid))
  )
)

#--------------------------------------------------------------
# 5. Plot triangular correlation heatmap
#--------------------------------------------------------------
Heatmap(
  cor_resid_tri,
  name = "Spearman Corr (resid)",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  right_annotation = row_ha_resid,
  column_title = paste0("Residualized Feature Correlations (", optimal_k_resid, " clusters)"),
  na_col = "white"
)

#--------------------------------------------------------------
# 6. Export cluster membership table
#--------------------------------------------------------------
cluster_table <- data.frame(
  Feature = rownames(residual_matrix)[feature_order_resid],
  OmicsType = feature_types_resid[feature_order_resid],
  Cluster = cluster_assign_resid[feature_order_resid]
)
write.csv(cluster_table, "feature_cluster_assignments.csv", row.names = FALSE)
cat("Cluster membership table saved to feature_cluster_assignments.csv\n")

#--------------------------------------------------------------
# 7. Identify top cross-omics feature pairs
#--------------------------------------------------------------
all_pairs <- which(upper.tri(cor_resid), arr.ind = TRUE)
pair_data <- data.frame(
  Feature1 = rownames(cor_resid)[all_pairs[, 1]],
  Feature2 = rownames(cor_resid)[all_pairs[, 2]],
  Corr = cor_resid[all_pairs]
)

# Omics type of each feature
pair_data$Type1 <- ifelse(pair_data$Feature1 %in% rna_features, "RNA",
                          ifelse(pair_data$Feature1 %in% met_features, "Metabolite", "Lipid"))
pair_data$Type2 <- ifelse(pair_data$Feature2 %in% rna_features, "RNA",
                          ifelse(pair_data$Feature2 %in% met_features, "Metabolite", "Lipid"))

# Keep only cross-omics pairs
pair_data_cross <- subset(pair_data, Type1 != Type2)

# Sort by absolute correlation strength and get top 20
pair_data_cross <- pair_data_cross[order(-abs(pair_data_cross$Corr)), ]
top_pairs <- head(pair_data_cross, 20)
print(top_pairs)
write.csv(top_pairs, "top_cross_omics_pairs.csv", row.names = FALSE)
cat("Top 20 cross-omics correlated pairs saved to top_cross_omics_pairs.csv\n")
#---------------------------------------------------------------------------------#

#--------------Trying Triangle Heatmap Again--------------------------------------#

# ================================================================
# Correlation Network & Triangle Heatmap for Old Animals
# ================================================================

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(igraph)

#----------------------------------------------------------
# 1. Subset to OLD animals (SED + PWR + IRAP + FRAP)
#----------------------------------------------------------
old_samples <- aging_index_df %>%
  filter(Group %in% c("OLD_SED", "OLD_PWR", "OLD_PWR_IRAP", "OLD_PWR_FRAP")) %>%
  pull(Sample)

aging_old <- aging_z[, old_samples, drop = FALSE]

#----------------------------------------------------------
# 2. Residualize group effect (remove mean differences)
#----------------------------------------------------------
sample_group <- aging_index_df$Group[match(colnames(aging_old), aging_index_df$Sample)]
design <- model.matrix(~ 0 + sample_group)
fit <- lm.fit(design, t(aging_old))
fitted_values <- t(fit$fitted.values)
rownames(fitted_values) <- rownames(aging_old)
residuals_matrix <- aging_old - fitted_values

#----------------------------------------------------------
# 3. Spearman correlation on residualized data
#----------------------------------------------------------
cor_matrix <- cor(t(residuals_matrix), method = "spearman")

#----------------------------------------------------------
# 4. Hierarchical clustering & feature ordering
#----------------------------------------------------------
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")
feature_order <- hc$order
cor_matrix_ord <- cor_matrix[feature_order, feature_order]

#----------------------------------------------------------
# 5. Omics type annotation
#----------------------------------------------------------
feature_types <- ifelse(rownames(aging_z) %in% rna_features, "RNA",
                        ifelse(rownames(aging_z) %in% met_features, "Metabolite", "Lipid"))
feature_types <- factor(feature_types, levels = c("RNA", "Metabolite", "Lipid"))
feature_types_ord <- feature_types[feature_order]
type_colors <- c("RNA" = "#1f78b4", "Metabolite" = "#33a02c", "Lipid" = "#ff7f00")

row_ha <- rowAnnotation(
  Omics = feature_types_ord,
  col = list(Omics = type_colors)
)

#----------------------------------------------------------
# 6. Triangle heatmap (upper triangle only)
#----------------------------------------------------------
cor_matrix_tri <- cor_matrix_ord
cor_matrix_tri[lower.tri(cor_matrix_tri)] <- NA

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

aging_corr_heatmap <- Heatmap(
  cor_matrix_tri,
  name = "Spearman\nCorr",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  right_annotation = row_ha,
  column_title = "Feature-feature correlation (old groups, residualized)",
  na_col = "white",
  row_names_gp = gpar(fontsize = 6)  # <-- adjust font size here
)

# Save as PDF
pdf("aging_corr_heatmap.pdf", width = 8, height = 8)  # Adjust height
print(aging_corr_heatmap)
dev.off()

#----------------------------------------------------------
# 7. Build correlation network (|ρ| >= 0.5)
#----------------------------------------------------------
cor_df <- as.data.frame(as.table(cor_matrix))
cor_df <- cor_df %>% filter(Var1 != Var2, abs(Freq) >= 0.5)

# Ensure each pair appears once (Var1 < Var2)
cor_df <- cor_df %>% 
  rowwise() %>% 
  mutate(pair = paste(sort(c(Var1, Var2)), collapse = "_")) %>% 
  ungroup() %>% 
  distinct(pair, .keep_all = TRUE)

# Create igraph network
edges <- cor_df %>% dplyr::select(Var1, Var2, weight = Freq)
nodes <- data.frame(name = unique(c(edges$Var1, edges$Var2)))
nodes$Omics <- ifelse(nodes$name %in% rna_features, "RNA",
                      ifelse(nodes$name %in% met_features, "Metabolite", "Lipid"))

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# Plot network
plot(g,
     vertex.color = type_colors[nodes$Omics],
     vertex.label = NA,
     vertex.size = 5,
     edge.width = abs(E(g)$weight) * 5,
     main = "Correlation Network (|ρ| ≥ 0.5)")

#----------------------------------------------------------
# 8. Output top cross-omics edges
#----------------------------------------------------------
# Keep edges that connect different omics types
cross_omics_edges <- edges %>%
  mutate(Omics1 = ifelse(Var1 %in% rna_features, "RNA",
                         ifelse(Var1 %in% met_features, "Metabolite", "Lipid")),
         Omics2 = ifelse(Var2 %in% rna_features, "RNA",
                         ifelse(Var2 %in% met_features, "Metabolite", "Lipid"))) %>%
  filter(Omics1 != Omics2) %>%
  arrange(desc(abs(weight)))

top_edges <- head(cross_omics_edges, 50)
write.csv(top_edges, "top_cross_omics_edges.csv", row.names = FALSE)
print(head(top_edges))
#------------------------------------------------#

#----------------------------------------------------------
# 1. Spearman correlation on residualized data (already computed)
#----------------------------------------------------------
# cor_matrix already computed from residuals_matrix

#----------------------------------------------------------
# 2. Convert correlation matrix to long format
#----------------------------------------------------------
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Var1", "Var2", "Corr")
cor_df <- cor_df %>% filter(Var1 != Var2)

# Ensure each pair appears only once (Var1 < Var2)
cor_df <- cor_df %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE)

#----------------------------------------------------------
# 3. Keep only strong correlations (|ρ| >= 0.5)
#----------------------------------------------------------
cor_df <- cor_df %>% filter(abs(Corr) >= 0.5)

#----------------------------------------------------------
# 4. Node table with omics type
#----------------------------------------------------------
nodes <- data.frame(name = unique(c(cor_df$Var1, cor_df$Var2)))
nodes$Omics <- ifelse(nodes$name %in% rna_features, "RNA",
                      ifelse(nodes$name %in% met_features, "Metabolite", "Lipid"))

# Color palette for node types
type_colors <- c("RNA" = "#1f78b4",
                 "Metabolite" = "#33a02c",
                 "Lipid" = "#ff7f00")

#----------------------------------------------------------
# 5. Edge table with sign & type
#----------------------------------------------------------
edges <- cor_df %>%
  mutate(weight = Corr,
         Sign = ifelse(weight > 0, "Positive", "Negative"),
         Omics1 = ifelse(Var1 %in% rna_features, "RNA",
                         ifelse(Var1 %in% met_features, "Metabolite", "Lipid")),
         Omics2 = ifelse(Var2 %in% rna_features, "RNA",
                         ifelse(Var2 %in% met_features, "Metabolite", "Lipid")),
         CrossOmics = Omics1 != Omics2)

# Edge colors based on sign
edge_colors <- ifelse(edges$Sign == "Positive", "red", "blue")

#----------------------------------------------------------
# 6. Create igraph network
#----------------------------------------------------------
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

#----------------------------------------------------------
# 7. Plot network
#----------------------------------------------------------
# Assign absolute weight for layout
E(g)$weight <- abs(E(g)$weight)

# Use Fruchterman-Reingold layout
layout_fr <- layout_with_fr(g, weights = E(g)$weight)

# Plot with improved layout
plot(g,
     layout = layout_fr,
     vertex.color = type_colors[V(g)$Omics],
     vertex.label = NA,
     vertex.size = 5,
     edge.width = abs(edges$weight) * 5,
     edge.color = edge_colors,
     main = "Correlation Network (|ρ| ≥ 0.5, layout improved)")

#----------------------------------------------------------
# 8. Output top cross-omics edges (for supplement)
#----------------------------------------------------------
top_edges <- edges %>%
  filter(CrossOmics) %>%
  arrange(desc(abs(weight))) %>%
  head(50)

write.csv(top_edges, "top_cross_omics_edges.csv", row.names = FALSE)
print(head(top_edges))

library(igraph)
deg <- degree(g, mode = "all")
btw <- betweenness(g, directed = FALSE)
hub_table <- data.frame(Feature = names(deg), Degree = deg, Betweenness = btw,
                        OmicsType = V(g)$Omics) %>% arrange(desc(Degree))
head(hub_table, 10)
#----------------------------------------------------#

# --- 1. Top 5 hubs
top_hubs <- hub_table$Feature[1:5]

# --- 2. Keep edges connected to top hubs
hub_edges <- edges %>%
  filter(Var1 %in% top_hubs | Var2 %in% top_hubs)

# --- 3. Build subgraph
g_sub <- graph_from_data_frame(hub_edges, vertices = nodes, directed = FALSE)

# --- 4. Scale node size by degree
deg_sub <- degree(g_sub)
V(g_sub)$size <- 5 + deg_sub * 2
V(g_sub)$label <- ifelse(names(deg_sub) %in% top_hubs, names(deg_sub), NA)

# --- 5. Edge colors (sign)
edge_colors_sub <- ifelse(E(g_sub)$weight > 0, "red", "blue")

# --- 6. Plot simplified network
layout_sub <- layout_with_fr(g_sub, weights = abs(E(g_sub)$weight))
plot(g_sub,
     layout = layout_sub,
     vertex.color = type_colors[V(g_sub)$Omics],
     vertex.label.color = "black",
     edge.width = abs(E(g_sub)$weight) * 5,
     edge.color = edge_colors_sub,
     main = "Top Hubs Cross-Omics Network (|ρ| ≥ 0.5)")

# --- 7. Save hub table and edge summary
write.csv(hub_table, "top_hub_table.csv", row.names = FALSE)
write.csv(hub_edges, "hub_edge_summary.csv", row.names = FALSE)
#-----------------------------------------------------------------#

#=================================================================#
#-------Simple Feature correlation to ITT results-----------------#
#=================================================================#

### --- 1) Insulin (glucose) named vector ------------------------------------
insulin_tolerance <- c(
  T_29 = 0.9483121927,
  T_06 = 1.048630335,
  T_21 = 1.156550841,
  T_25 = 1.068462587,
  T_09 = 1.038714209,
  T_13 = 1.038714209,
  T_19 = 0.7982481511,
  T_39 = 0.9023674751,
  T_08 = 0.8131223402,
  T_46 = 0.741230426,
  T_41 = 0.590009503,
  T_10 = 0.7660207412,
  T_45 = 0.8169235219,
  T_43 = 0.633475189,
  T_11 = 0.6557864728,
  T_16 = 0.8342767426,
  T_02 = 0.8094864273,
  T_03 = 0.7313142999,
  T_37 = 0.840391687,
  T_47 = 1.028798083,
  T_14 = 0.7698219229,
  T_48 = 0.9197206958,
  T_33 = 0.9122836012,
  T_31 = 0.8279965294,
  T_44 = 0.6817336694,
  T_34 = 1.259348015,
  T_04 = 1.168945998,
  T_01 = 1.054910548,
  T_22 = 1.022683138,
  T_28 = 1.059868611,
  T_40 = 0.6706606619,
  T_24 = 0.892451349
)

## --- 2) Align to matrix: keep overlap only --------------------------------
X_full <- as.matrix(integrated_matrix)   # rows = features, cols = samples
overlap <- intersect(colnames(X_full), names(insulin_tolerance))

# Optional: quick diagnostics
if (length(setdiff(colnames(X_full), names(insulin_tolerance))) > 0)
  message("Samples in matrix but not insulin_tolerance: ",
          paste(setdiff(colnames(X_full), names(insulin_tolerance)), collapse = ", "))
if (length(setdiff(names(insulin_tolerance), colnames(X_full))) > 0)
  message("Samples in insulin_tolerance but not matrix: ",
          paste(setdiff(names(insulin_tolerance), colnames(X_full)), collapse = ", "))

# Subset to common samples and align order
X <- X_full[, overlap, drop = FALSE]
y <- insulin_tolerance[overlap]

## --- 3) Row-wise Spearman correlations (feature vs insulin tolerance) -----
rho <- as.numeric(cor(y, t(X), method = "spearman", use = "pairwise.complete.obs"))

# Per-feature N actually used (handles any NAs in X)
finite_mask <- is.finite(X) & matrix(TRUE, nrow = nrow(X), ncol = ncol(X))
n_obs <- rowSums(finite_mask)  # usually == length(y) unless X has NAs

# Two-sided p-values via t approximation
tval <- rho * sqrt(pmax(n_obs - 2, 1) / pmax(1 - rho^2, 1e-12))
pval <- 2 * pt(-abs(tval), df = pmax(n_obs - 2, 1))
fdr  <- p.adjust(pval, method = "BH")

res <- data.frame(
  Feature   = rownames(X),
  SpearmanR = rho,
  N         = n_obs,
  P         = pval,
  FDR       = fdr,
  stringsAsFactors = FALSE
)

## --- 4) Keep positives only; take top 100 by rho --------------------------
res_pos <- res[res$SpearmanR > 0, ]
res_pos <- res_pos[order(-res_pos$SpearmanR), ]

top100_pos <- head(res_pos, 100)
top100_pos
#----------------------------------------------#

#==============================================#
#----Correlation Triangle All Features (~1,000)--#
#================================================#

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cluster)

# Inspect data
dim(integrated_matrix)

# Step 1: Identify omics types
feature_names <- rownames(integrated_matrix)

head(feature_names)

# Assign omics type
feature_types <- ifelse(grepl("^RNA_", feature_names), "RNA",
                        ifelse(grepl("^Met_", feature_names), "Metabolite",
                               ifelse(grepl("^Lip_", feature_names), "Lipid", "Other")))

head(feature_types)

# Recode "Other" as "Lipid" (based on your clarification)
feature_types[feature_types == "Other"] <- "Lipid"
# Reassign feature types directly from the subset
feature_names <- rownames(integrated_matrix)
feature_types <- case_when(
  grepl("^RNA_", feature_names) ~ "RNA",
  grepl("^Met_", feature_names) ~ "Metabolite",
  grepl("^Lip_", feature_names) ~ "Lipid",
  TRUE ~ "Lipid"
)
feature_types <- factor(feature_types, levels = c("RNA", "Metabolite", "Lipid"))

head(integrated_matrix)

# Step 2: Keep top 1000 most variable features
feature_variances <- apply(integrated_matrix, 1, var, na.rm = TRUE)
top_features <- names(sort(feature_variances, decreasing = TRUE))[1:1000]
integrated_matrix <- integrated_matrix[top_features, ]

# Step 3: Per-omics z-scoring
integrated_matrix_z <- integrated_matrix
omics_levels <- levels(feature_types)

head(omics_levels)

for (otype in omics_levels) {
  idx <- which(feature_types == otype)
  integrated_matrix_z[idx, ] <- t(scale(t(integrated_matrix[idx, ])))
}

# Step 4: Spearman correlation matrix
cor_matrix <- cor(t(integrated_matrix_z), method = "spearman")

# Step 5: Hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")
d <- dist(1 - cor_matrix)

# Determine optimal k (based on silhouette)
sil_scores <- sapply(2:10, function(k) {
  cluster_assign <- cutree(hc, k)
  sil <- silhouette(cluster_assign, d)
  mean(sil[, 3])
})
optimal_k <- which.max(sil_scores)
cat("Optimal number of clusters:", optimal_k, "\n")

cluster_assign <- cutree(hc, k = optimal_k)

# Step 6: Prepare for plotting
feature_order <- hc$order
cor_matrix_ord <- cor_matrix[feature_order, feature_order]
cor_matrix_tri <- cor_matrix_ord
cor_matrix_tri[lower.tri(cor_matrix_tri)] <- NA

feature_types_ord <- feature_types[feature_order]
cluster_assign_ord <- cluster_assign[feature_order]

type_colors <- c("RNA" = "#1f78b4", "Metabolite" = "#33a02c", "Lipid" = "#ff7f00")
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
  annotation_name_gp = gpar(fontsize = 9),
  annotation_legend_param = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 9))
)

# Color scale
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Step 7: Plot triangular correlation heatmap
Heatmap(
  cor_matrix_tri,
  name = "Spearman\nCorr",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  right_annotation = row_ha,
  column_title = paste0("Feature-feature Spearman correlation (", optimal_k, " clusters)"),
  na_col = "white"
)

# Step 8: Export outputs
cor_df <- as.data.frame(as.table(cor_matrix))
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]
cor_df <- cor_df %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(desc(abs(Freq))) %>%
  rename(Corr = Freq)

write.csv(head(cor_df, 50), "top_feature_correlations.csv", row.names = FALSE)

cluster_table <- data.frame(
  Feature = rownames(cor_matrix)[feature_order],
  OmicsType = feature_types_ord,
  Cluster = cluster_assign_ord
)
write.csv(cluster_table, "feature_cluster_assignments.csv", row.names = FALSE)

# Optional: Validate gene cluster assignments
validated_inflamm_geneset_upper.vec <- c(
  "PRKCZ","JCHAIN","ANKRD1","HIP1R","KNG2","CES2C","UVRAG","FOSL2","GAS6","TNFRSF23",
  "MID2","CXCL10","CX3CL1","IL20RB","RAPGEF1","CYP4F14","TRP63","POU4F1","MYC","RAB20",
  "RUNX1","CASP4","LGALS1","DLL1","CCL8","SFRP1","CD55","BAX","POU2F2","CCL7",
  "HSP90AA1","ARID5A","BBC3","MIF","RELB","NINJ1","NR1H3","BID","CDKN1A","SOX4",
  "SBNO2","CCL2","ICAM1","TUBB5","APOD","AKT1","SLC15A4","FAM110A","GPX1","TYK2",
  "PPARD","FOS"
)
validated_inflamm_geneset_rna.vec <- paste0("RNA_", str_to_title(tolower(validated_inflamm_geneset_upper.vec)))

inflamm_cluster_members <- cluster_table %>%
  filter(Feature %in% validated_inflamm_geneset_rna.vec)

print(inflamm_cluster_members)

cluster_table %>%
  filter(Cluster == 3)


getwd()
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")

# Load the integrated matrix from RDS file
integrated_matrix <- readRDS("integrated_matrix_all_omics_v2.RDS")
integrated_matrix
dim(integrated_matrix)
#------------------------------------------------------------------#
#integrated_matrix <- readRDS("integrated_matrix_all_omics_v2.RDS") #don't write over this

#==================================================================#
#------Force include inflamm genes--------------------------------#
#==================================================================#

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cluster)

dim(integrated_matrix)
# Step 1: Identify omics types
feature_names <- rownames(integrated_matrix)
feature_types <- case_when(
  grepl("^RNA_", feature_names) ~ "RNA",
  grepl("^Met_", feature_names) ~ "Metabolite",
  grepl("^Lip_", feature_names) ~ "Lipid",
  TRUE ~ "Lipid"
)
feature_types <- factor(feature_types, levels = c("RNA", "Metabolite", "Lipid"))

# Step 2: Define gene set to force-include
validated_inflamm_geneset_upper.vec <- c(
  "PRKCZ","JCHAIN","ANKRD1","HIP1R","KNG2","CES2C","UVRAG","FOSL2","GAS6","TNFRSF23",
  "MID2","CXCL10","CX3CL1","IL20RB","RAPGEF1","CYP4F14","TRP63","POU4F1","MYC","RAB20",
  "RUNX1","CASP4","LGALS1","DLL1","CCL8","SFRP1","CD55","BAX","POU2F2","CCL7",
  "HSP90AA1","ARID5A","BBC3","MIF","RELB","NINJ1","NR1H3","BID","CDKN1A","SOX4",
  "SBNO2","CCL2","ICAM1","TUBB5","APOD","AKT1","SLC15A4","FAM110A","GPX1","TYK2",
  "PPARD","FOS"
)
validated_inflamm_geneset_rna.vec <- paste0("RNA_", str_to_title(tolower(validated_inflamm_geneset_upper.vec)))

# Step 3: Select features to retain
feature_variances <- apply(integrated_matrix, 1, var, na.rm = TRUE)

# Identify omics types
is_rna <- feature_types == "RNA"
is_met <- feature_types == "Metabolite"
is_lip <- feature_types == "Lipid"

# Get RNA gene names
rna_features <- rownames(integrated_matrix)[is_rna]

# Remove RIKEN genes (typically end in "Rik" or similar)
rna_no_riken <- rna_features[!grepl("Rik$", rna_features, ignore.case = TRUE)]

# Calculate variances only for non-RIKEN RNA genes
top_rna <- names(sort(feature_variances[rna_no_riken], decreasing = TRUE))[1:1000]

# Keep all metabolites and lipids
all_met <- names(feature_variances[is_met])
all_lip <- names(feature_variances[is_lip])

# Force include inflammation gene set
force_include <- validated_inflamm_geneset_rna.vec
force_include <- force_include[force_include %in% rownames(integrated_matrix)]

# Union all selected features
selected_features <- unique(c(top_rna, all_met, all_lip, force_include))
selected_features <- selected_features[selected_features %in% rownames(integrated_matrix)]

# Subset the matrix
integrated_matrix <- integrated_matrix[selected_features, ]

integrated_matrix

# Recalculate feature_types for filtered matrix
feature_names <- rownames(integrated_matrix)
feature_types <- case_when(
  grepl("^RNA_", feature_names) ~ "RNA",
  grepl("^Met_", feature_names) ~ "Metabolite",
  grepl("^Lip_", feature_names) ~ "Lipid",
  TRUE ~ "Lipid"
)
feature_types <- factor(feature_types, levels = c("RNA", "Metabolite", "Lipid"))

# Step 4: Per-omics z-scoring
integrated_matrix_z <- integrated_matrix
omics_levels <- levels(feature_types)
for (otype in omics_levels) {
  idx <- which(feature_types == otype)
  integrated_matrix_z[idx, ] <- t(scale(t(integrated_matrix[idx, ])))
}

# Step 5: Spearman correlation matrix
cor_matrix <- cor(t(integrated_matrix_z), method = "spearman")

# Step 6: Hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")
d <- dist(1 - cor_matrix)

# Optimal cluster number using silhouette
sil_scores <- sapply(2:10, function(k) {
  cluster_assign <- cutree(hc, k)
  sil <- silhouette(cluster_assign, d)
  mean(sil[, 3])
})

# Override if silhouette is flat or max is at 1
if (length(unique(sil_scores)) == 1 || which.max(sil_scores) == 1) {
  warning("Silhouette index flat or favors 1 cluster. Using fallback of k = 10.")
  optimal_k <- 10
} else {
  optimal_k <- which.max(sil_scores)
}
cat("Chosen number of clusters:", optimal_k, "\n")

cluster_assign <- cutree(hc, k = optimal_k)

# Step 7: Prepare for plotting
feature_order <- hc$order
cor_matrix_ord <- cor_matrix[feature_order, feature_order]
cor_matrix_tri <- cor_matrix_ord
cor_matrix_tri[lower.tri(cor_matrix_tri)] <- NA

feature_types_ord <- feature_types[feature_order]
cluster_assign_ord <- cluster_assign[feature_order]

type_colors <- c("RNA" = "#1f78b4", "Metabolite" = "#33a02c", "Lipid" = "#ff7f00")
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
  annotation_name_gp = gpar(fontsize = 9),
  annotation_legend_param = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 9))
)

# Step 8: Plot triangular correlation heatmap
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

png("heatmap_output_v2.png", width = 2000, height = 2000, res = 300)
Heatmap(
  cor_matrix_tri,
  name = "Spearman\nCorr",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  right_annotation = row_ha,
  column_title = paste0("Feature-feature Spearman correlation (", optimal_k, " clusters)"),
  na_col = "white",
  use_raster = TRUE
)
dev.off()

cor_df <- as.data.frame(as.table(cor_matrix))
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]

# Vectorized way to generate unique unordered feature pairs
cor_df <- cor_df %>%
  mutate(
    gene1 = pmin(as.character(Var1), as.character(Var2)),
    gene2 = pmax(as.character(Var1), as.character(Var2)),
    pair = paste(gene1, gene2, sep = "_")
  ) %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(desc(abs(Freq))) %>%
  rename(Corr = Freq)

write.csv(head(cor_df, 50), "top_feature_correlations.csv", row.names = FALSE)

cluster_table <- data.frame(
  Feature = rownames(cor_matrix)[feature_order],
  OmicsType = feature_types_ord,
  Cluster = cluster_assign_ord
)
write.csv(cluster_table, "feature_cluster_assignments.csv", row.names = FALSE)

# Optional: Validate gene cluster assignments
inflamm_cluster_members <- cluster_table %>%
  filter(Feature %in% validated_inflamm_geneset_rna.vec)

print(inflamm_cluster_members)

# Example: list all features in cluster 3
cluster_table
dim(cluster_table %>%
  filter(Cluster == 4))

# List all lipids in cluster 1
cluster_table %>%
  filter(Cluster == 1, OmicsType == "Metabolite") %>%
  arrange(Feature)
#--------------------------------------#

#----Overrepresentation Analysis------#

library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)

# Prepare gene clusters (RNA only)
rna_clusters <- cluster_table %>%
  filter(OmicsType == "RNA") %>%
  mutate(GeneSymbol = sub("^RNA_", "", Feature)) %>%
  group_by(Cluster) %>%
  summarise(GeneSymbols = list(GeneSymbol), .groups = "drop")

rna_clusters

# Create a combined table of GO results for all clusters
go_results_list <- rna_clusters %>%
  mutate(
    ego = map(
      GeneSymbols,
      ~ enrichGO(
        gene         = .x,
        OrgDb        = org.Mm.eg.db,
        keyType      = "SYMBOL",
        ont          = "BP",  # Biological Process
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      )
    )
  )

# Combine all enrichment results into one table
go_enrichment_table <- go_results_list %>%
  mutate(results = map2(Cluster, ego, ~ as_tibble(.y) %>% mutate(Cluster = .x))) %>%
  dplyr::select(results) %>%
  unnest(cols = results)

# Preview result
head(go_enrichment_table)

go_enrichment_table %>%
  filter(Cluster == 2) %>%
    pull(Description)

# Add simplified GO terms per cluster
go_results_simplified <- go_results_list %>%
  mutate(
    ego_simplified = map(
      ego,
      ~ simplify(
        .x,
        cutoff = 0.7,             # Adjust similarity threshold as needed
        by = "p.adjust",
        select_fun = min,
        measure = "Wang"
      )
    )
  )

# Unnest simplified results into one tidy table
go_enrichment_simplified_table <- go_results_simplified %>%
  mutate(results = map2(Cluster, ego_simplified, ~ as_tibble(.y) %>% mutate(Cluster = .x))) %>%
  dplyr::select(results) %>%
  unnest(cols = results)

# Preview
head(go_enrichment_simplified_table)

cluster10_GO <- go_enrichment_simplified_table %>%
  filter(Cluster == 10) %>%
  pull(Description)

cluster10_GO

write.csv(cluster10_GO, "cluster10_GO.csv", row.names = FALSE)

#-----------------------------------------#

cluster8_metabolites <- cluster_table %>%
  filter(Cluster == 8, OmicsType == "Metabolite") %>%
  arrange(Feature) %>%
  pull(Feature)

cluster8_metabolites

write.csv(cluster10_metabolites, "cluster10_metabolites.csv", row.names = FALSE)

metabolite_list <- cluster_table %>%
  filter(Cluster == 1, OmicsType == "Metabolite") %>%
  arrange(Feature) %>%
  pull(Feature) %>%
  str_remove("^Met_")  # Remove "Met_" prefix

# Optional: write to a .txt file for MetaboAnalyst upload
writeLines(metabolite_list, "cluster1_metabolites.txt")
#----------------------------------------------------------#


#=========================================================#
#------Trying to get Lipid and Metabolite heatmaps--------#
#==========================================================#

#---------Lipid Heatmap------------------------------------#
integrated_matrix <- readRDS("integrated_matrix_all_omics_v2.RDS")

integrated_matrix
lipid_matrix
new_lipid_age_signature

# Clean lipid names in signature to match matrix rownames (remove escaped characters)
lipid_signature_clean <- gsub("\\\\", "", new_lipid_age_signature)

# Subset the matrix
lipid_subset_matrix <- lipid_matrix[rownames(lipid_matrix) %in% lipid_signature_clean, , drop = FALSE]

lipid_subset_matrix

library(ComplexHeatmap)
library(circlize)


# Z-score across rows (samples)
lipid_z <- t(scale(t(lipid_subset_matrix)))
lipid_z[!is.finite(lipid_z)] <- 0  # Replace NaN/Inf with 0

# (Optional) cap extremes for nicer colors
lipid_z_cap <- pmin(pmax(lipid_z, -3), 3)

library(pheatmap)

# Define color palette
lipid_max_abs <- 2.25
lipid_cols <- colorRampPalette(c("#FF7F00", "black", "#00ffff"))(100)
lipid_brks <- seq(-lipid_max_abs, lipid_max_abs, length.out = length(lipid_cols) + 1)

# Plot heatmap
pheatmap_obj <- pheatmap(
  lipid_z_cap,
  color = lipid_cols,
  breaks = lipid_brks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 6,
  main = "Z-scored log2 Lipid Values"
)

# Save as PDF
pdf("aged_lipid_signature_heatmap.pdf", width = 5, height = 6)
print(pheatmap_obj)
dev.off()
#----------------------------------------------------#

# Combine your vectors into one vector of samples in desired order
desired_order3 <- c(
  old_sed_samples,
  old_pwr_samples,
  old_pwrirap_samples,
  old_pwrfrap_samples
)

# Subset metabo_matrix to only include matching columns, in the correct order
metabo_matrix_subset <- metabo_matrix[, desired_order3, drop = FALSE]

dim(metabo_matrix_subset)

aging_metabolite_vec <- c(
  "Kynurenic acid",
  "D-Lactose/cellobiose",
  "Isomaltulose",
  "Orotic Acid",
  "Acetyl-L-leucine",
  "D-Maltose",
  "DL-Methionine sulfoxide",
  "Melibiose",
  "N-Acetyl-DL-tryptophan",
  "3-Ureidopropionic acid",
  "Glycerol",
  "1-Aminocyclopropane-1-carboxylic acid",
  "B-Alanine",
  "3R-hydroxy-isobutyric acid",
  "D-Sorbitol",
  "DL-2-Aminoadipic acid",
  "N-Acetyl-L-asparagine"
)

# Subset the matrix
aging_metabo_matrix <- metabo_matrix_subset[rownames(metabo_matrix_subset) %in% aging_metabolite_vec, , drop = FALSE]

aging_metabo_matrix

# Z-score across rows (samples)
metabo_z <- t(scale(t(aging_metabo_matrix)))
metabo_z[!is.finite(metabo_z)] <- 0  # Replace NaN/Inf with 0

# (Optional) cap extremes for nicer colors
metabo_z_cap <- pmin(pmax(metabo_z, -3), 3)

library(pheatmap)

# Define color palette
metabo_max_abs <- 2.25
#metabo_cols <- colorRampPalette(c("#FF7F00", "black", "#00ffff"))(100)
#metabo_brks <- seq(-lipid_max_abs, lipid_max_abs, length.out = length(lipid_cols) + 1)

# Plot heatmap
metabo_intervention_pheatmap <- pheatmap(
  metabo_z_cap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 6,
  main = "Z-scored log2 Lipid Values"
)

# Save as PDF
pdf("metabo_intervention_pheatmap.pdf", width = 5, height = 6)
print(metabo_intervention_pheatmap)
dev.off()
