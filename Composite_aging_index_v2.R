
#===========================================#
#------Define Omics sets-------------------#
#===========================================#

## Build logCPM_symbol (SYMBOL × samples) across YNG/OLD/OV/IRAP/FRAP

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(tidyr)
  library(edgeR);  library(AnnotationDbi); library(org.Mm.eg.db)
})

# ---- Load counts ----
set01_os_vs_ys  <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
set01_ov_vs_ys  <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx")
set01_oir_vs_ys <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_IRAP-YNG_SED_VEH.xlsx")
set01_ofr_vs_ys <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_FRAP-YNG_SED_VEH.xlsx")

# 1) Rename the first column of each counts tibble to "Ensembl"
rename_first_col <- function(df) {
  cn <- colnames(df)
  cn[1] <- "Ensembl"
  colnames(df) <- cn
  df
}

set01_os_vs_ys  <- rename_first_col(set01_os_vs_ys)
set01_ov_vs_ys  <- rename_first_col(set01_ov_vs_ys)
set01_oir_vs_ys <- rename_first_col(set01_oir_vs_ys)
set01_ofr_vs_ys <- rename_first_col(set01_ofr_vs_ys)

# quick sanity check
stopifnot(all(c("Ensembl") %in% c(names(set01_os_vs_ys)[1],
                                  names(set01_ov_vs_ys)[1],
                                  names(set01_oir_vs_ys)[1],
                                  names(set01_ofr_vs_ys)[1])))

# ---- Sample groups (keep your canonical IDs) ----
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwr_samples     <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

# Helper: keep only columns that actually exist (prevents select() errors)
keep_present <- function(df, cols) intersect(cols, colnames(df))

# Start from OS vs YS (has YNG + OLD SED), then add unique columns from other sets
merged_counts <- set01_os_vs_ys %>%
  dplyr::select(
    Ensembl,
    dplyr::all_of(keep_present(set01_os_vs_ys, c(yng_sed_samples, old_sed_samples)))
  ) %>%
  left_join(
    set01_ov_vs_ys %>%
      dplyr::select(Ensembl, dplyr::all_of(keep_present(set01_ov_vs_ys, old_pwr_samples))),
    by = "Ensembl"
  ) %>%
  left_join(
    set01_oir_vs_ys %>%
      dplyr::select(Ensembl, dplyr::all_of(keep_present(set01_oir_vs_ys, old_pwrirap_samples))),
    by = "Ensembl"
  ) %>%
  left_join(
    set01_ofr_vs_ys %>%
      dplyr::select(Ensembl, dplyr::all_of(keep_present(set01_ofr_vs_ys, old_pwrfrap_samples))),
    by = "Ensembl"
  )

set01_ofr_vs_ys

# Map Ensembl -> ENTREZ
merged_counts <- merged_counts %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db,
                           keys = Ensembl_noDec,
                           keytype = "ENSEMBL",
                           column = "ENTREZID",
                           multiVals = "first")) %>%
  tidyr::drop_na(ENTREZID)

# Counts matrix (NA→0), rows=ENTREZ, columns in exact group order expected by edgeR
all_samples <- c(yng_sed_samples, old_sed_samples, old_pwr_samples,
                 old_pwrirap_samples, old_pwrfrap_samples)
present_samples <- intersect(all_samples, colnames(merged_counts))

counts_matrix <- merged_counts %>%
  dplyr::select(dplyr::all_of(present_samples)) %>%
  replace(is.na(.), 0) %>%
  as.matrix()
rownames(counts_matrix) <- merged_counts$ENTREZID

# Define groups to match column order exactly
group <- factor(c(
  rep("YNG_SED",      sum(present_samples %in% yng_sed_samples)),
  rep("OLD_SED",      sum(present_samples %in% old_sed_samples)),
  rep("OLD_PWR",      sum(present_samples %in% old_pwr_samples)),
  rep("OLD_PWR_IRAP", sum(present_samples %in% old_pwrirap_samples)),
  rep("OLD_PWR_FRAP", sum(present_samples %in% old_pwrfrap_samples))
))

stopifnot(length(group) == ncol(counts_matrix))

# edgeR: TMM, filter, logCPM
dge <- DGEList(counts = counts_matrix, group = group)
dge <- calcNormFactors(dge, method = "TMM")
keep <- filterByExpr(dge)                    # uses the group factor above
dge <- dge[keep, , keep.lib.sizes = FALSE]
logCPM_matrix <- cpm(dge, log = TRUE, prior.count = 1)

cat("Genes kept:", nrow(logCPM_matrix), " Samples:", ncol(logCPM_matrix), "\n")

# Map ENTREZ -> SYMBOL, collapse duplicates by mean, ensure unique SYMBOL rownames
symbols <- mapIds(org.Mm.eg.db,
                  keys    = rownames(logCPM_matrix),
                  keytype = "ENTREZID",
                  column  = "SYMBOL",
                  multiVals = "first")

sym_df <- as.data.frame(logCPM_matrix) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  mutate(Symbol = symbols) %>%
  tidyr::drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

logCPM_symbol <- sym_df %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()

# Optional: keep column order as present_samples
logCPM_symbol <- logCPM_symbol[, present_samples, drop = FALSE]

# Sanity checks
stopifnot(!anyNA(logCPM_symbol))
stopifnot(any(colnames(logCPM_symbol) %in% yng_sed_samples))  # has YNG
stopifnot(any(colnames(logCPM_symbol) %in% old_sed_samples))  # has OLD SED

saveRDS(logCPM_symbol, "logCPM_symbol_final.RDS")
dim(logCPM_symbol)
#-------------------------------------------------#
#--------------------------------------------------#

#----------metabo_log-----------------------------#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

# 0) Load
pwr_rapa_muscle_metabo_data <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)

# 1) Keep Sample + features (stash Group separately if you need it later)
metabo_wide <- pwr_rapa_muscle_metabo_data %>% dplyr::select(-Group)

# 2) Long → Wide (rows = metabolites, cols = samples)
metabo_matrix <- metabo_wide %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()

rownames(metabo_matrix) <- metabo_matrix$Metabolite
metabo_matrix$Metabolite <- NULL

# 3) Map metabolomics sample names (OV1/OIR1/OFR1/OS1/YS1/YV1...) → T_* tube codes
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

# Vectorized rename (keeps original names if not in map)
new_names <- name_map[colnames(metabo_matrix)]
colnames(metabo_matrix) <- ifelse(is.na(new_names), colnames(metabo_matrix), new_names)

# Warn if any columns could not be mapped
unmapped <- setdiff(colnames(metabo_matrix), name_map)
# (Fine to ignore if these are already T_* names)

# 4) Order/align to the RNA-Seq sample order (and intersect, just in case)
desired_order <- c(yng_sed_samples, old_sed_samples, old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples)

present_metabo <- intersect(desired_order, colnames(metabo_matrix))
missing_for_metabo <- setdiff(desired_order, colnames(metabo_matrix))
if (length(missing_for_metabo)) message("Metabolomics missing: ", paste(missing_for_metabo, collapse = ", "))

# If logCPM_symbol exists, also intersect with it to guarantee one-to-one columns
if (exists("logCPM_symbol")) {
  present_rna <- intersect(desired_order, colnames(logCPM_symbol))
  common_samples <- intersect(present_metabo, present_rna)
} else {
  common_samples <- present_metabo
}

stopifnot(length(common_samples) > 0)
metabo_matrix_subset <- metabo_matrix[, common_samples, drop = FALSE]

# Ensure numeric (guard against character columns from CSV)
metabo_matrix_subset[] <- lapply(metabo_matrix_subset, function(x) as.numeric(as.character(x)))

# 5) Impute zeros/NA per metabolite with half the minimum non-zero, then log2
metabo_imputed <- t(apply(metabo_matrix_subset, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))

metabo_log <- as.data.frame(log2(metabo_imputed))
colnames(metabo_log) <- common_samples
rownames(metabo_log) <- rownames(metabo_matrix_subset)

# Optional: line up RNA too (same columns) for downstream composite calculations
if (exists("logCPM_symbol")) {
  logCPM_symbol <- logCPM_symbol[, common_samples, drop = FALSE]
  stopifnot(identical(colnames(metabo_log), colnames(logCPM_symbol)))
}

# Quick peek
dim(metabo_log); head(metabo_log[, 1:4])

# Save for reuse
saveRDS(metabo_log, "metabo_log_final.RDS")
#-------------------------------------------#

#-------lipid_matrix-----------------------#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
})

# 0) Load + drop YS6 outlier
lipid_df <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE) %>%
  dplyr::filter(Sample != "YS6")  # <- remove if you decide to keep it

# 1) Identify lipid columns
lipid_cols <- setdiff(names(lipid_df), c("Sample","Group"))

# 2) Impute zeros/NA per lipid, then log2
lipid_df[ , lipid_cols] <- lapply(lipid_df[ , lipid_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})

# 3) Matrix: rows = lipids, cols = samples
lipid_matrix <- t(as.matrix(lipid_df[ , lipid_cols]))
colnames(lipid_matrix) <- lipid_df$Sample

# 4) Drop YV samples
lipid_matrix <- lipid_matrix[ , !grepl("^YV", colnames(lipid_matrix)), drop = FALSE]

# 5) Rename samples to T_* using a vectorized map (doesn't create NAs)
#    sample_key must have columns: TubeCode, MetabolomicsName
map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
mapped <- map[colnames(lipid_matrix)]
colnames(lipid_matrix) <- ifelse(is.na(mapped), colnames(lipid_matrix), mapped)

# 6) Choose your target order and align by intersection (avoids manual T_32 removal)
desired_order <- c(
  yng_sed_samples,
  old_sed_samples,
  old_pwr_samples,
  old_pwrirap_samples,
  old_pwrfrap_samples
)

present <- intersect(desired_order, colnames(lipid_matrix))
lipid_matrix <- lipid_matrix[ , present, drop = FALSE]

# Optional: for the “lipid heatmap only” view without young groups:
desired_order2 <- c(old_sed_samples, old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples)
present2 <- intersect(desired_order2, colnames(lipid_matrix))
lipid_matrix_heat <- lipid_matrix[ , present2, drop = FALSE]  # use this for Fig 3-style heatmaps

# 7) (Optional but recommended) Align with RNA for composite index
if (exists("logCPM_symbol")) {
  common <- intersect(colnames(lipid_matrix), colnames(logCPM_symbol))
  lipid_matrix   <- lipid_matrix[ , common, drop = FALSE]
  logCPM_symbol  <- logCPM_symbol[ , common, drop = FALSE]
  stopifnot(identical(colnames(lipid_matrix), colnames(logCPM_symbol)))
}

# Quick sanity
dim(lipid_matrix); anyNA(lipid_matrix)




## -----------------------------
## Composite Aging Index (revised — bidirectional)
## -----------------------------

library(dplyr); library(tidyr); library(stringr); library(ggplot2)

# 1) Harmonize: keep only samples present in all three omics
common_samples <- Reduce(intersect, list(
  colnames(logCPM_symbol),
  colnames(metabo_log),
  colnames(lipid_matrix)
))

stopifnot(length(common_samples) >= 6)

rna_mat <- logCPM_symbol[, common_samples, drop = FALSE]
met_mat <- metabo_log[,    common_samples, drop = FALSE]
lip_mat <- lipid_matrix[,  common_samples, drop = FALSE]

# 2) Select features per block using NEW signatures
#    — GENES (93, with direction)
rna_features <- validated_aging_genes %>%
  filter(Symbol %in% rownames(rna_mat))
cat("RNA features matched:", nrow(rna_features), "\n")

#    — LIPIDS (59, direction from Log2_FC_OS_vs_YS)
lip_features <- sig_lipid_features %>%
  filter(Lipid %in% rownames(lip_mat)) %>%
  mutate(aging_direction = ifelse(Log2_FC_OS_vs_YS > 0, "Up", "Down"))
cat("Lipid features matched:", nrow(lip_features), "\n")

#    — METABOLITES (24, direction from Log2_FC_OS_vs_YS)
met_features <- metab_sig_df %>%
  filter(Metabolite %in% rownames(met_mat)) %>%
  mutate(aging_direction = ifelse(Log2_FC_OS_vs_YS > 0, "Up", "Down"))
cat("Metabolite features matched:", nrow(met_features), "\n")

# 3) Z-score by row WITHIN each block
z_by_row <- function(m) t(scale(t(m)))

rna_z <- z_by_row(rna_mat[rna_features$Symbol,    , drop = FALSE])
met_z <- z_by_row(met_mat[met_features$Metabolite, , drop = FALSE])
lip_z <- z_by_row(lip_mat[lip_features$Lipid,      , drop = FALSE])

# 4) Flip sign for features that DECREASE with age
#    so that positive z always = "more aged"
rna_down <- rna_features$Symbol[rna_features$aging_direction == "Down"]
rna_z[rownames(rna_z) %in% rna_down, ] <- -1 * rna_z[rownames(rna_z) %in% rna_down, ]

met_down <- met_features$Metabolite[met_features$aging_direction == "Down"]
met_z[rownames(met_z) %in% met_down, ] <- -1 * met_z[rownames(met_z) %in% met_down, ]

lip_down <- lip_features$Lipid[lip_features$aging_direction == "Down"]
lip_z[rownames(lip_z) %in% lip_down, ] <- -1 * lip_z[rownames(lip_z) %in% lip_down, ]

cat("RNA features flipped (Down):", length(rna_down), "\n")
cat("Lipid features flipped (Down):", length(lip_down), "\n")
cat("Metabolite features flipped (Down):", length(met_down), "\n")

# 5) Per-sample block means, then equal-weight composite
rna_index   <- colMeans(rna_z, na.rm = TRUE)
met_index   <- colMeans(met_z, na.rm = TRUE)
lip_index   <- colMeans(lip_z, na.rm = TRUE)

composite_index <- (rna_index + met_index + lip_index) / 3

# 6) Assemble long table with 5 group labels
group_of <- function(sid) {
  if (sid %in% yng_sed_samples)        return("YNG SED")
  if (sid %in% old_sed_samples)        return("OLD SED")
  if (sid %in% old_pwr_samples)        return("OLD PWR")
  if (sid %in% old_pwrirap_samples)    return("OLD PWR IRAP")
  if (sid %in% old_pwrfrap_samples)    return("OLD PWR FRAP")
  return(NA_character_)
}

df_composite <- data.frame(
  Sample       = common_samples,
  RNA_Index    = rna_index[common_samples],
  Met_Index    = met_index[common_samples],
  Lipid_Index  = lip_index[common_samples],
  AgingIndex   = composite_index[common_samples],
  Group        = vapply(common_samples, group_of, character(1))
) %>% filter(!is.na(Group))

# 7) Boxplot
df_composite$Group <- factor(df_composite$Group,
                             levels = c("YNG SED","OLD SED","OLD PWR","OLD PWR IRAP","OLD PWR FRAP"))

p_box <- ggplot(df_composite, aes(x = Group, y = AgingIndex, fill = Group)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.9) +
  scale_fill_manual(values = c(
    "YNG SED" = "#1f78b4",
    "OLD SED" = "#7f7f7f",
    "OLD PWR" = "#33a02c",
    "OLD PWR IRAP"= "#6a3d9a",
    "OLD PWR FRAP"= "#ff7f00"
  )) +
  labs(title = "Composite Aging Index",
       y = "Mean z-score (equal-weight RNA / Met / Lipid)", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
print(p_box)

# 8) Stats
anova_model <- aov(AgingIndex ~ Group, data = df_composite)
summary(anova_model)
TukeyHSD(anova_model)

# 9) Save
saveRDS(logCPM_symbol, "logCPM_symbol_final.RDS")
saveRDS(metabo_log,    "metabo_log_final.RDS")
saveRDS(lipid_matrix,  "lipid_matrix_final.RDS")
saveRDS(df_composite,  "df_composite_final.RDS")

df_composite

# Long format (one row per sample, easy to paste into Prism)
write.csv(df_composite, "composite_aging_index_for_prism.csv", row.names = FALSE)

yng_mean <- df_composite %>%
  filter(Group == "YNG SED") %>%
  pull(AgingIndex) %>%
  mean()

df_composite <- df_composite %>%
  mutate(AgingIndex_centered = AgingIndex - yng_mean)

df_composite#
