#logCPM heat map of all old interventions
#use average logCPM for each intervention.

# Intervention Samples
yng_sed_samples <- c("T_05", "T_17", "T_26", "T_32", "T_35", "T_38")
old_sedveh_samples <- c("T_06", "T_09", "T_13", "T_19", "T_21", "T_25", "T_29", "T_39")
old_pwrveh_samples <- c("T_08", "T_10", "T_11", "T_16", "T_41", "T_43", "T_45", "T_46")
old_sedirap_samples <- c("T_49", "T_50", "T_51", "T_52", "T_53")
old_sedfrap_samples <- c("T_54", "T_55", "T_56", "T_57", "T_58")
old_pwrirap_samples <- c("T_03", "T_14", "T_31", "T_33", "T_37", "T_47", "T_48")
old_pwrfrap_samples <- c("T_01", "T_04", "T_22", "T_24", "T_28", "T_34", "T_40", "T_44")

#read in Count files
old_sedveh_counts <-read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
old_pwrveh_counts <-read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx")
old_sedirap_counts <-read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx")
old_sedfrap_counts <-read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx")
old_pwrirap_counts <-read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_IRAP-YNG_SED_VEH.xlsx")
old_pwrfrap_counts <-read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_FRAP-YNG_SED_VEH.xlsx")

head(old_sedveh_counts)
head(old_pwrveh_counts)
head(old_sedirap_counts)
head(old_sedfrap_counts)
head(old_pwrirap_counts)
head(old_pwrfrap_counts)


# Function to clean individual count tables
clean_counts <- function(df) {
  colnames(df)[1] <- "Ensembl"                     # Rename first column
  df <- df %>%
    select(-all_of(yng_sed_samples))               # Remove yng_sed_samples
  return(df)
}

# Read and clean each dataframe
old_sedveh_counts   <- clean_counts(read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx"))
old_pwrveh_counts   <- clean_counts(read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx"))
old_sedirap_counts  <- clean_counts(read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx"))
old_sedfrap_counts  <- clean_counts(read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx"))
old_pwrirap_counts  <- clean_counts(read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_IRAP-YNG_SED_VEH.xlsx"))
old_pwrfrap_counts  <- clean_counts(read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_FRAP-YNG_SED_VEH.xlsx"))

# Join all dataframes by "Ensembl"
combined_counts <- list(
  old_sedveh_counts,
  old_pwrveh_counts,
  old_sedirap_counts,
  old_sedfrap_counts,
  old_pwrirap_counts,
  old_pwrfrap_counts
) %>%
  reduce(full_join, by = "Ensembl")

# Replace NA with 0 (for downstream EdgeR use)
combined_counts[is.na(combined_counts)] <- 0

# View the result
head(combined_counts)
#---------------------------------#

# Step 1: Strip version from Ensembl IDs
all_intervention_counts <- combined_counts %>%
  dplyr::mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"))

# Step 2: Map to ENTREZID
all_intervention_counts$ENTREZID <- mapIds(
  org.Mm.eg.db,
  keys = all_intervention_counts$Ensembl_noDec,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

# Step 3: Remove rows with missing Entrez IDs
all_intervention_counts <- all_intervention_counts %>%
  dplyr::filter(!is.na(ENTREZID))

# Step 4 (alternative): Group by ENTREZID and sum across samples
all_intervention_counts <- all_intervention_counts %>%
  select(-Ensembl, -Ensembl_noDec) %>%
  group_by(ENTREZID) %>%
  summarise(across(everything(), sum)) %>%
  ungroup()

# Step 5: Set ENTREZID as rownames
all_intervention_counts <- all_intervention_counts %>%
  column_to_rownames("ENTREZID")

head(all_intervention_counts)
#-------------------------------------------#

# 1. Combine all sample names in the correct order
all_samples <- c(
  old_sedveh_samples,
  old_pwrveh_samples,
  old_sedirap_samples,
  old_sedfrap_samples,
  old_pwrirap_samples,
  old_pwrfrap_samples
)

# 2. Create a group factor to label each sample
group <- factor(c(rep("OLD_SEDVEH", length(old_sedveh_samples)),
                  rep("OLD_PWRVEH", length(old_pwrveh_samples)),
                  rep("OLD_SEDIRAP", length(old_sedirap_samples)),
                  rep("OLD_SEDFRAP", length(old_sedfrap_samples)),
                  rep("OLD_PWRIRAP", length(old_pwrirap_samples)),
                  rep("OLD_PWRFRAP", length(old_pwrfrap_samples))
))

# 3. Extract counts matrix
counts_matrix <- all_intervention_counts %>%
  dplyr::select(all_of(all_samples)) %>%
  as.matrix()

# 4. Ensure rownames are set to ENTREZID
# This should already be true if you followed the last step correctly
# But just in case:
#rownames(counts_matrix) <- rownames(all_intervention_counts)

# 5. Create DGEList object
dge <- DGEList(counts = counts_matrix, group = group)

# 6. Apply TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# 7. Convert to logCPM (prior.count = 1 is standard)
logCPM_matrix <- cpm(dge, log = TRUE, prior.count = 1)

# ✅ Optional: View part of the result
head(logCPM_matrix)
#------------------------------------------#

inflammaging_gene_vec

inflammaging_gene_vec <- as.character(inflammaging_gene_vec)

# Subset logCPM_matrix to only rows with ENTREZID in inflammaging_gene_vec
inflammaging_logCPM <- logCPM_matrix[rownames(logCPM_matrix) %in% inflammaging_gene_vec, ]

# Optional: Check dimensions and preview
dim(inflammaging_logCPM)
head(inflammaging_logCPM)
#-----------------------------------------------#

# Compute group averages
inflammaging_logCPM_grouped <- data.frame(
  OLD_SEDVEH   = rowMeans(inflammaging_logCPM[, old_sedveh_samples], na.rm = TRUE),
  OLD_PWRVEH   = rowMeans(inflammaging_logCPM[, old_pwrveh_samples], na.rm = TRUE),
  OLD_SEDIRAP  = rowMeans(inflammaging_logCPM[, old_sedirap_samples], na.rm = TRUE),
  OLD_SEDFRAP  = rowMeans(inflammaging_logCPM[, old_sedfrap_samples], na.rm = TRUE),
  OLD_PWRIRAP  = rowMeans(inflammaging_logCPM[, old_pwrirap_samples], na.rm = TRUE),
  OLD_PWRFRAP  = rowMeans(inflammaging_logCPM[, old_pwrfrap_samples], na.rm = TRUE)
)

# Retain ENTREZID as rownames
rownames(inflammaging_logCPM_grouped) <- rownames(inflammaging_logCPM)

# Preview the result
head(inflammaging_logCPM_grouped)

inflammaging_logCPM_grouped
#--------------------------------------------#

# Get gene symbols for rownames (ENTREZ IDs)
entrez_ids <- rownames(inflammaging_logCPM_grouped)
gene_symbols <- mapIds(org.Mm.eg.db, 
                       keys = entrez_ids, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

# Replace rownames with gene symbols (keeping NAs as-is)
rownames(inflammaging_logCPM_grouped) <- gene_symbols

heatmap_matrix <- as.matrix(inflammaging_logCPM_grouped)

# Save heatmap as PDF
pdf("inflammaging_logCPM_heatmap.pdf", width = 8, height = 10)

# Generate heatmap
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",  # standardize each gene
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Inflammaging Gene Expression (logCPM)")

dev.off()




# Generate heatmap with manual column order

pheatmap(inflammaging_logCPM_grouped, 
                                     cluster_rows = TRUE, 
                                     cluster_cols = TRUE, 
                                     scale = "row",  # Standardize by row (gene)
                                     show_rownames = TRUE, 
                                     show_colnames = TRUE, 
                                     color = colorRampPalette(c("blue", "white", "red"))(50))
