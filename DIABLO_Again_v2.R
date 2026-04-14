###########################################################################
# DIABLO Consolidated Analysis Pipeline
# Combines: DIABLO_Again.R + DIABLO_tuning.R
# Cleaned: removed redundant code, dead ends, and broken sections
###########################################################################

setwd("/Users/mdbruss/Documents/RStudioProjects_2/Rapa_PwR")

#===========================================================================#
# SECTION 0: Libraries
#===========================================================================#

suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(stringr)
  library(tidyr); library(tibble); library(ggplot2); library(purrr)
  library(edgeR); library(AnnotationDbi); library(org.Mm.eg.db)
  library(mixOmics)
})

#===========================================================================#
# SECTION 1: Define Sample Groups (used throughout)
#===========================================================================#

yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwr_samples     <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

train_samples <- c(yng_sed_samples, old_sed_samples)
all_samples   <- c(yng_sed_samples, old_sed_samples, old_pwr_samples,
                   old_pwrirap_samples, old_pwrfrap_samples)

# Shared sample key for metabolomics/lipidomics renaming
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

#===========================================================================#
# SECTION 2: Helper Functions
#===========================================================================#

rename_first_col <- function(df) {
  colnames(df)[1] <- "Ensembl"
  df
}

# Impute zeros/NA with half-min, then log2
impute_log2 <- function(mat) {
  imputed <- t(apply(mat, 1, function(x) {
    x[is.na(x)] <- 0
    min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    if (!is.finite(min_pos)) min_pos <- 1e-6
    x[x == 0] <- min_pos / 2
    x
  }))
  log2(imputed)
}

# Rename columns using sample_key
rename_samples <- function(mat) {
  mapped <- name_map[colnames(mat)]
  colnames(mat) <- ifelse(is.na(mapped), colnames(mat), mapped)
  mat
}

# Assign group label from sample ID
group_of <- function(sid) {
  case_when(
    sid %in% yng_sed_samples     ~ "Young_Sed",
    sid %in% old_sed_samples     ~ "Old_SedVeh",
    sid %in% old_pwr_samples     ~ "Old_PwrVeh",
    sid %in% old_pwrirap_samples ~ "Old_PwrIRap",
    sid %in% old_pwrfrap_samples ~ "Old_PwrFRap",
    TRUE                         ~ "Other"
  )
}

#===========================================================================#
# SECTION 3A: Build logCPM_symbol — TRAINING (YNG + OLD SED only)
#===========================================================================#

set01_os_vs_ys <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx") %>%
  rename_first_col()

counts_train <- set01_os_vs_ys %>%
  dplyr::select(Ensembl, dplyr::all_of(train_samples)) %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db, keys = Ensembl_noDec,
                           keytype = "ENSEMBL", column = "ENTREZID",
                           multiVals = "first")) %>%
  tidyr::drop_na(ENTREZID)

counts_mat_train <- counts_train %>%
  dplyr::select(dplyr::all_of(train_samples)) %>%
  replace(is.na(.), 0) %>% as.matrix()
rownames(counts_mat_train) <- counts_train$ENTREZID

group_train <- factor(c(rep("YNG_SED", length(yng_sed_samples)),
                        rep("OLD_SED", length(old_sed_samples))))

dge <- DGEList(counts = counts_mat_train, group = group_train)
dge <- calcNormFactors(dge, method = "TMM")
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
logCPM_train <- cpm(dge, log = TRUE, prior.count = 1)

symbols <- mapIds(org.Mm.eg.db, keys = rownames(logCPM_train),
                  keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")

logCPM_symbol_train <- as.data.frame(logCPM_train) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  mutate(Symbol = symbols) %>% drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("Symbol") %>% as.matrix()
logCPM_symbol_train <- logCPM_symbol_train[, train_samples, drop = FALSE]

saveRDS(logCPM_symbol_train, "logCPM_symbol_yngOLDsed.RDS")
cat("Training RNA:", nrow(logCPM_symbol_train), "genes x", ncol(logCPM_symbol_train), "samples\n")

#===========================================================================#
# SECTION 3B: Build logCPM_symbol — FULL (all 5 groups)
#===========================================================================#

keep_present <- function(df, cols) intersect(cols, colnames(df))

set01_ov_vs_ys  <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx") %>% rename_first_col()
set01_oir_vs_ys <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_IRAP-YNG_SED_VEH.xlsx") %>% rename_first_col()
set01_ofr_vs_ys <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_FRAP-YNG_SED_VEH.xlsx") %>% rename_first_col()

merged_counts <- set01_os_vs_ys %>%
  dplyr::select(Ensembl, all_of(keep_present(set01_os_vs_ys, c(yng_sed_samples, old_sed_samples)))) %>%
  left_join(set01_ov_vs_ys  %>% dplyr::select(Ensembl, all_of(keep_present(set01_ov_vs_ys,  old_pwr_samples))),     by = "Ensembl") %>%
  left_join(set01_oir_vs_ys %>% dplyr::select(Ensembl, all_of(keep_present(set01_oir_vs_ys, old_pwrirap_samples))), by = "Ensembl") %>%
  left_join(set01_ofr_vs_ys %>% dplyr::select(Ensembl, all_of(keep_present(set01_ofr_vs_ys, old_pwrfrap_samples))), by = "Ensembl") %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db, keys = Ensembl_noDec,
                           keytype = "ENSEMBL", column = "ENTREZID",
                           multiVals = "first")) %>%
  drop_na(ENTREZID)

present_all <- intersect(all_samples, colnames(merged_counts))
counts_mat_all <- merged_counts %>%
  dplyr::select(all_of(present_all)) %>%
  replace(is.na(.), 0) %>% as.matrix()
rownames(counts_mat_all) <- merged_counts$ENTREZID

group_all <- factor(c(
  rep("YNG_SED",      sum(present_all %in% yng_sed_samples)),
  rep("OLD_SED",      sum(present_all %in% old_sed_samples)),
  rep("OLD_PWR",      sum(present_all %in% old_pwr_samples)),
  rep("OLD_PWR_IRAP", sum(present_all %in% old_pwrirap_samples)),
  rep("OLD_PWR_FRAP", sum(present_all %in% old_pwrfrap_samples))
))

dge_all <- DGEList(counts = counts_mat_all, group = group_all)
dge_all <- calcNormFactors(dge_all, method = "TMM")
dge_all <- dge_all[filterByExpr(dge_all), , keep.lib.sizes = FALSE]
logCPM_all <- cpm(dge_all, log = TRUE, prior.count = 1)

symbols_all <- mapIds(org.Mm.eg.db, keys = rownames(logCPM_all),
                      keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")

logCPM_symbol_all <- as.data.frame(logCPM_all) %>%
  rownames_to_column("ENTREZID") %>%
  mutate(Symbol = symbols_all) %>% drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("Symbol") %>% as.matrix()
logCPM_symbol_all <- logCPM_symbol_all[, present_all, drop = FALSE]

saveRDS(logCPM_symbol_all, "logCPM_symbol_final.RDS")
cat("Full RNA:", nrow(logCPM_symbol_all), "genes x", ncol(logCPM_symbol_all), "samples\n")

#===========================================================================#
# SECTION 4: Build Metabolomics Matrices (training + full)
#===========================================================================#

metabo_raw <- read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)
metabo_wide <- metabo_raw %>% dplyr::select(-Group) %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()
rownames(metabo_wide) <- metabo_wide$Metabolite
metabo_wide$Metabolite <- NULL
metabo_wide <- rename_samples(metabo_wide)
metabo_wide[] <- lapply(metabo_wide, function(x) as.numeric(as.character(x)))

# Full (all samples)
metabo_log_all <- as.data.frame(impute_log2(as.matrix(metabo_wide)))
saveRDS(metabo_log_all, "metabo_log_final.RDS")

# Training subset
metabo_log_train <- metabo_log_all[, intersect(train_samples, colnames(metabo_log_all)), drop = FALSE]
saveRDS(metabo_log_train, "metabo_log_yngOLDsed.RDS")

cat("Metabolomics: full =", ncol(metabo_log_all), "samples | train =", ncol(metabo_log_train), "samples\n")

#===========================================================================#
# SECTION 5: Build Lipidomics Matrices (training + full)
#===========================================================================#

lipid_df <- read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)
lipid_cols <- setdiff(names(lipid_df), c("Sample","Group"))

lipid_df[, lipid_cols] <- lapply(lipid_df[, lipid_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})

lipid_matrix <- t(as.matrix(lipid_df[, lipid_cols]))
colnames(lipid_matrix) <- lipid_df$Sample
lipid_matrix <- lipid_matrix[, !grepl("^YV", colnames(lipid_matrix)), drop = FALSE]
lipid_matrix <- rename_samples(lipid_matrix)

# Full (all samples, drop YS6 outlier for consistency)
lipid_matrix <- lipid_matrix[, colnames(lipid_matrix) != "T_32", drop = FALSE]  # T_32 = YS6
present_lip_all <- intersect(all_samples, colnames(lipid_matrix))
lipid_full <- lipid_matrix[, present_lip_all, drop = FALSE]

# Training subset
lipid_train <- lipid_matrix[, intersect(train_samples, colnames(lipid_matrix)), drop = FALSE]
saveRDS(lipid_train, "lipid_log_yngOLDsed.RDS")

cat("Lipidomics: full =", ncol(lipid_full), "samples | train =", ncol(lipid_train), "samples\n")

#===========================================================================#
# SECTION 6: Harmonize Sample Order Across All Omics
#===========================================================================#

# Training
common_train <- Reduce(intersect, list(
  colnames(logCPM_symbol_train),
  colnames(metabo_log_train),
  colnames(lipid_train)
))
cat("Common training samples:", length(common_train), "\n")

RNA_train    <- t(logCPM_symbol_train[, common_train, drop = FALSE])
Metabo_train <- t(as.matrix(metabo_log_train[, common_train, drop = FALSE]))
Lipid_train  <- t(lipid_train[, common_train, drop = FALSE])

# Full
common_full <- Reduce(intersect, list(
  colnames(logCPM_symbol_all),
  colnames(metabo_log_all),
  colnames(lipid_full)
))
cat("Common full samples:", length(common_full), "\n")

RNA_full    <- t(logCPM_symbol_all[, common_full, drop = FALSE])
Metabo_full <- t(as.matrix(metabo_log_all[, common_full, drop = FALSE]))
Lipid_full  <- t(lipid_full[, common_full, drop = FALSE])

# Final alignment checks
stopifnot(identical(rownames(RNA_train), rownames(Metabo_train)))
stopifnot(identical(rownames(RNA_train), rownames(Lipid_train)))
stopifnot(identical(rownames(RNA_full),  rownames(Metabo_full)))
stopifnot(identical(rownames(RNA_full),  rownames(Lipid_full)))

cat("\n--- Data Ready ---\n")
cat("Training:", nrow(RNA_train), "samples |",
    ncol(RNA_train), "RNA |", ncol(Metabo_train), "Met |", ncol(Lipid_train), "Lip\n")
cat("Full:    ", nrow(RNA_full), "samples |",
    ncol(RNA_full), "RNA |", ncol(Metabo_full), "Met |", ncol(Lipid_full), "Lip\n")

#===========================================================================#
# SECTION 7: DIABLO Core Functions
#===========================================================================#

#--- 7A: Run a single DIABLO model, project all samples ---#

run_diablo_once <- function(
    RNA_train, Metabo_train, Lipid_train,
    RNA_full, Metabo_full, Lipid_full,
    yng_sed_samples, old_sed_samples,
    old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
    RNA_keep, Met_keep, Lip_keep,
    design_strength = 1.0
) {
  # Build training objects
  X_train <- list(RNA = RNA_train, Metabolite = Metabo_train, Lipid = Lipid_train)
  Y_train <- factor(c(rep("YNG", sum(rownames(RNA_train) %in% yng_sed_samples)),
                      rep("OLD", sum(rownames(RNA_train) %in% old_sed_samples))),
                    levels = c("YNG","OLD"))
  
  design <- matrix(design_strength, 3, 3,
                   dimnames = list(names(X_train), names(X_train)))
  diag(design) <- 0
  
  # Train
  diablo_fit <- block.splsda(X_train, Y_train, ncomp = 1,
                             keepX = list(RNA = RNA_keep,
                                          Metabolite = Met_keep,
                                          Lipid = Lip_keep),
                             design = design)
  
  # Project all samples
  X_all <- list(RNA = RNA_full, Metabolite = Metabo_full, Lipid = Lipid_full)
  proj <- predict(diablo_fit, newdata = X_all)
  scores <- lapply(proj$variates, function(x) x[, 1])
  scores_mat <- do.call(cbind, scores)
  comp1 <- rowMeans(scores_mat)
  
  # Flip so higher = older
  comp1 <- -comp1
  
  tibble(
    Sample = rownames(scores_mat),
    Comp1  = comp1,
    RNA_keep = RNA_keep,
    Met_keep = Met_keep,
    Lip_keep = Lip_keep,
    DesignStrength = design_strength,
    Group = group_of(rownames(scores_mat))
  )
}

#--- 7B: Run grid of DIABLO models ---#

run_diablo_grid <- function(
    RNA_train, Metabo_train, Lipid_train,
    RNA_full, Metabo_full, Lipid_full,
    yng_sed_samples, old_sed_samples,
    old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
    RNA_vec  = c(30, 50, 100),
    Met_vec  = c(10, 20, 30),
    Lip_vec  = c(10, 20, 30),
    DS_vec   = c(1.0, 0.8)
) {
  grid <- expand.grid(RNA_keep = RNA_vec, Met_keep = Met_vec,
                      Lip_keep = Lip_vec, DS = DS_vec,
                      stringsAsFactors = FALSE)
  
  message("Running DIABLO for ", nrow(grid), " parameter combinations...")
  
  all_scores <- pmap_dfr(grid, function(RNA_keep, Met_keep, Lip_keep, DS) {
    run_diablo_once(
      RNA_train, Metabo_train, Lipid_train,
      RNA_full, Metabo_full, Lipid_full,
      yng_sed_samples, old_sed_samples,
      old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
      RNA_keep, Met_keep, Lip_keep,
      design_strength = DS
    )
  })
  
  group_summary <- all_scores %>%
    group_by(RNA_keep, Met_keep, Lip_keep, DesignStrength, Group) %>%
    summarise(n = n(), mean_Comp1 = mean(Comp1),
              sd_Comp1 = sd(Comp1), se_Comp1 = sd_Comp1/sqrt(n),
              .groups = "drop")
  
  list(sample_scores = all_scores, group_summary = group_summary)
}

#--- 7C: ANOVA + Tukey for every hyperparameter combo ---#

run_tukey_grid <- function(sample_scores) {
  sample_scores %>%
    group_by(RNA_keep, Met_keep, Lip_keep, DesignStrength) %>%
    nest() %>%
    mutate(tukey_results = purrr::map(data, function(df) {
      fit <- aov(Comp1 ~ Group, data = df)
      tuk <- TukeyHSD(fit)
      tuk_df <- as.data.frame(tuk$Group)
      tuk_df$Comparison <- rownames(tuk_df)
      tuk_df %>% dplyr::select(Comparison, p.adj = `p adj`)
    })) %>%
    dplyr::select(-data) %>%
    unnest(tukey_results)
}

#--- 7D: Build metric table with biological quality scores ---#

build_metric_table <- function(group_summary) {
  # Pivot wide
  wide <- group_summary %>%
    dplyr::select(RNA_keep, Met_keep, Lip_keep, DesignStrength, Group, mean_Comp1) %>%
    pivot_wider(names_from = Group, values_from = mean_Comp1) %>%
    arrange(RNA_keep, Met_keep, Lip_keep, DesignStrength)
  
  # Harmonize direction (ensure Old > Young)
  group_cols <- c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap")
  wide <- wide %>%
    rowwise() %>%
    mutate(flip = ifelse(Old_SedVeh < Young_Sed, -1, 1),
           across(all_of(group_cols), ~ .x * flip)) %>%
    ungroup() %>% dplyr::select(-flip)
  
  # Check expected ordering: Young < PwrVeh < PwrIRap < PwrFRap < SedVeh
  wide %>%
    rowwise() %>%
    mutate(
      order_ok = all(diff(c(Young_Sed, Old_PwrVeh, Old_PwrIRap, Old_PwrFRap, Old_SedVeh)) > 0),
      age_separation    = Old_SedVeh - Young_Sed,
      restoration_score = Old_SedVeh - Old_PwrVeh,
      IRap_interference = Old_PwrIRap - Old_PwrVeh,
      FRap_interference = Old_PwrFRap - Old_PwrVeh
    ) %>%
    ungroup() %>%
    mutate(
      score_quality =
        1 * order_ok +
        scales::rescale(age_separation,    to = c(0, 1)) +
        scales::rescale(restoration_score, to = c(0, 1)) +
        scales::rescale(IRap_interference, to = c(0, 1)) +
        scales::rescale(FRap_interference, to = c(0, 1))
    ) %>%
    arrange(desc(score_quality))
}

#===========================================================================#
# SECTION 8: Run the Grid Search
#===========================================================================#

grid_results <- run_diablo_grid(
  RNA_train, Metabo_train, Lipid_train,
  RNA_full,  Metabo_full,  Lipid_full,
  yng_sed_samples, old_sed_samples,
  old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
  RNA_vec = c(25, 50, 100, 150, 200),
  Met_vec = c(10, 20, 30, 40, 50),
  Lip_vec = c(10, 20, 30, 40, 50),
  DS_vec  = c(1.0, 0.8)
)

#===========================================================================#
# SECTION 9: Build Metric + Tukey Tables
#===========================================================================#

metric_table <- build_metric_table(grid_results$group_summary)

tukey_wide <- run_tukey_grid(grid_results$sample_scores) %>%
  pivot_wider(names_from = Comparison, values_from = p.adj, names_prefix = "p_")

metric_table_with_tukey <- metric_table %>%
  left_join(tukey_wide, by = c("RNA_keep","Met_keep","Lip_keep","DesignStrength"))

# Save results
write_csv(metric_table_with_tukey, "metric_table_with_tukey.csv")
write_csv(grid_results$sample_scores, "DIABLO_grid_sample_scores.csv")

# Preview best models
cat("\n--- Top 10 Models by Quality Score ---\n")
print(metric_table_with_tukey %>% slice_max(score_quality, n = 10))

#===========================================================================#
# SECTION 10: Explore a Single Best Model
#===========================================================================#

# << EDIT THESE based on your metric_table results >>
best_params <- list(RNA_keep = 50, Met_keep = 20, Lip_keep = 20, DS = 1.0)

best_df <- grid_results$sample_scores %>%
  filter(RNA_keep == best_params$RNA_keep,
         Met_keep == best_params$Met_keep,
         Lip_keep == best_params$Lip_keep,
         DesignStrength == best_params$DS)

# Harmonize direction
yng_mean <- mean(best_df$Comp1[best_df$Group == "Young_Sed"])
old_mean <- mean(best_df$Comp1[best_df$Group == "Old_SedVeh"])
if (old_mean < yng_mean) best_df$Comp1 <- -best_df$Comp1

best_df$Group <- factor(best_df$Group,
                        levels = c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap"))

# Stats
cat("\n--- ANOVA for Best Model ---\n")
summary(aov(Comp1 ~ Group, data = best_df))
cat("\n--- Tukey HSD ---\n")
print(TukeyHSD(aov(Comp1 ~ Group, data = best_df)))

# Plot
group_colors <- c(Young_Sed = "#1F78B4", Old_SedVeh = "#E31A1C",
                  Old_PwrVeh = "#33A02C", Old_PwrIRap = "#FB9A99",
                  Old_PwrFRap = "#FDBF6F")

p_best <- ggplot(best_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  labs(title = sprintf("DIABLO Aging Axis (RNA=%d, Met=%d, Lip=%d, DS=%.1f)",
                       best_params$RNA_keep, best_params$Met_keep,
                       best_params$Lip_keep, best_params$DS),
       x = "", y = "Comp1 Aging Axis Score") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

print(p_best)
ggsave("DIABLO_best_model_boxplot.pdf", p_best, width = 7, height = 5)

#===========================================================================#
# SECTION 11: Robustness — Distribution Across All Models
#===========================================================================#

model_values_long <- metric_table %>%
  dplyr::select(RNA_keep, Met_keep, Lip_keep, DesignStrength,
                Young_Sed, Old_SedVeh, Old_PwrVeh, Old_PwrIRap, Old_PwrFRap) %>%
  pivot_longer(cols = c(Young_Sed, Old_SedVeh, Old_PwrVeh, Old_PwrIRap, Old_PwrFRap),
               names_to = "Group", values_to = "Comp1") %>%
  mutate(Group = factor(Group,
                        levels = c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap")))

p_robust <- ggplot(model_values_long, aes(x = Group, y = Comp1, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "DIABLO Comp1 Across All Model Configurations",
       subtitle = paste(n_distinct(model_values_long$RNA_keep), "x",
                        "feature combos x design strengths"),
       x = "", y = "Comp1 Score") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

print(p_robust)
ggsave("DIABLO_robustness_violin.pdf", p_robust, width = 8, height = 5)

#===========================================================================#
# SECTION 12: Extract & Save Selected Features from Best Model
#===========================================================================#

# Re-run the best model to get the DIABLO object
X_train_best <- list(RNA = RNA_train, Metabolite = Metabo_train, Lipid = Lipid_train)
Y_train_best <- factor(c(rep("YNG", sum(rownames(RNA_train) %in% yng_sed_samples)),
                         rep("OLD", sum(rownames(RNA_train) %in% old_sed_samples))),
                       levels = c("YNG","OLD"))

design_best <- matrix(best_params$DS, 3, 3,
                      dimnames = list(names(X_train_best), names(X_train_best)))
diag(design_best) <- 0

diablo_best <- block.splsda(X_train_best, Y_train_best, ncomp = 1,
                            keepX = list(RNA = best_params$RNA_keep,
                                         Metabolite = best_params$Met_keep,
                                         Lipid = best_params$Lip_keep),
                            design = design_best)

# Extract features
rna_sel <- selectVar(diablo_best, block = "RNA", comp = 1)
met_sel <- selectVar(diablo_best, block = "Metabolite", comp = 1)
lip_sel <- selectVar(diablo_best, block = "Lipid", comp = 1)

extract_features <- function(sel, block) {
  data.frame(
    block   = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    abs_loading = abs(sel[[block]]$value$value.var),
    stringsAsFactors = FALSE
  )
}

selected_df <- bind_rows(
  extract_features(rna_sel, "RNA"),
  extract_features(met_sel, "Metabolite"),
  extract_features(lip_sel, "Lipid")
) %>% arrange(block, desc(abs_loading))

write_csv(selected_df, "DIABLO_selected_features_best_model.csv")
cat("\nSelected features saved:", nrow(selected_df), "total\n")
print(selected_df, n = 20)

# Loadings plots
plotLoadings(diablo_best, block = "RNA",        comp = 1, method = "mean")
plotLoadings(diablo_best, block = "Metabolite", comp = 1, method = "mean")
plotLoadings(diablo_best, block = "Lipid",      comp = 1, method = "mean")

# Circos plot
circosPlot(diablo_best, comp = 1, cutoff = 0.7, line = TRUE)

cat("\n=== DIABLO Pipeline Complete ===\n")
##############################################################################
##############################################################################
#
#DIABLO with (100 gene, 50 lipid, 20 metabolites)
#-------------------------------------------------#

#===========================================================================#
# Run final DIABLO model (RNA=100, Met=20, Lip=50) and extract scores
#===========================================================================#

# Train
X_train_final <- list(RNA = RNA_train, Metabolite = Metabo_train, Lipid = Lipid_train)
Y_train_final <- factor(c(
  rep("YNG", sum(rownames(RNA_train) %in% yng_sed_samples)),
  rep("OLD", sum(rownames(RNA_train) %in% old_sed_samples))
), levels = c("YNG","OLD"))

design_final <- matrix(1.0, 3, 3,
                       dimnames = list(names(X_train_final), names(X_train_final)))
diag(design_final) <- 0

diablo_final <- block.splsda(
  X_train_final, Y_train_final, ncomp = 1,
  keepX = list(RNA = 100, Metabolite = 20, Lipid = 50),
  design = design_final
)

# Project all samples
X_all_final <- list(RNA = RNA_full, Metabolite = Metabo_full, Lipid = Lipid_full)
proj_final <- predict(diablo_final, newdata = X_all_final)

scores_final <- lapply(proj_final$variates, function(x) x[, 1])
scores_mat_final <- do.call(cbind, scores_final)
comp1_final <- rowMeans(scores_mat_final)

# Flip so higher = older
comp1_final <- -comp1_final

# Build dataframe
diablo_scores_df <- data.frame(
  Sample = names(comp1_final),
  Comp1  = comp1_final,
  Group  = group_of(names(comp1_final))
) %>%
  mutate(Group = factor(Group,
                        levels = c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap")))

# Harmonize direction
yng_m <- mean(diablo_scores_df$Comp1[diablo_scores_df$Group == "Young_Sed"])
old_m <- mean(diablo_scores_df$Comp1[diablo_scores_df$Group == "Old_SedVeh"])
if (old_m < yng_m) diablo_scores_df$Comp1 <- -diablo_scores_df$Comp1

# Center on Young = 0 (like your composite aging index)
yng_mean_diablo <- mean(diablo_scores_df$Comp1[diablo_scores_df$Group == "Young_Sed"])
diablo_scores_df$Comp1_centered <- diablo_scores_df$Comp1 - yng_mean_diablo

# Stats
cat("\n--- ANOVA ---\n")
summary(aov(Comp1_centered ~ Group, data = diablo_scores_df))
cat("\n--- Tukey HSD ---\n")
print(TukeyHSD(aov(Comp1_centered ~ Group, data = diablo_scores_df)))

# Wide format for Prism (columns = groups, rows = replicates)
diablo_prism_wide <- diablo_scores_df %>%
  dplyr::select(Sample, Group, Comp1_centered) %>%
  group_by(Group) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Group, values_from = Comp1_centered) %>%
  dplyr::select(-Sample, -row)

diablo_prism_wide

write.csv(diablo_scores_df, "DIABLO_final_scores_long.csv", row.names = FALSE)
write.csv(diablo_prism_wide, "DIABLO_final_scores_prism.csv", row.names = FALSE)

# Save selected features
rna_sel_final <- selectVar(diablo_final, block = "RNA", comp = 1)
met_sel_final <- selectVar(diablo_final, block = "Metabolite", comp = 1)
lip_sel_final <- selectVar(diablo_final, block = "Lipid", comp = 1)

extract_features <- function(sel, block) {
  data.frame(
    block   = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    stringsAsFactors = FALSE
  )
}

selected_features_final <- bind_rows(
  extract_features(rna_sel_final, "RNA"),
  extract_features(met_sel_final, "Metabolite"),
  extract_features(lip_sel_final, "Lipid")
) %>% arrange(block, desc(abs(loading)))

write.csv(selected_features_final, "DIABLO_final_selected_features.csv", row.names = FALSE)

cat("\nFiles saved:\n")
cat("  - DIABLO_final_scores_long.csv (one row per sample)\n")
cat("  - DIABLO_final_scores_prism.csv (wide format for Prism)\n")
cat("  - DIABLO_final_selected_features.csv (170 features with loadings)\n")







#===========================================================================#
# SECTION 13: Permutation Test — Random Feature Selection
#===========================================================================#

suppressPackageStartupMessages({
  library(mixOmics); library(dplyr); library(purrr); library(ggplot2)
})

# Feature counts to match chosen model
perm_RNA_keep <- 100
perm_Met_keep <- 20
perm_Lip_keep <- 50
perm_DS       <- 1.0
n_perms       <- 1000

set.seed(42)

cat("Running", n_perms, "permutations with random feature selection...\n")

perm_results <- purrr::map_dbl(seq_len(n_perms), function(i) {
  
  if (i %% 100 == 0) cat("  Permutation", i, "/", n_perms, "\n")
  
  # Randomly sample features from each block
  rna_rand <- sample(colnames(RNA_train), perm_RNA_keep)
  met_rand <- sample(colnames(Metabo_train), perm_Met_keep)
  lip_rand <- sample(colnames(Lipid_train), perm_Lip_keep)
  
  # Subset training and full matrices to random features
  X_train_perm <- list(
    RNA        = RNA_train[, rna_rand, drop = FALSE],
    Metabolite = Metabo_train[, met_rand, drop = FALSE],
    Lipid      = Lipid_train[, lip_rand, drop = FALSE]
  )
  
  X_full_perm <- list(
    RNA        = RNA_full[, rna_rand, drop = FALSE],
    Metabolite = Metabo_full[, met_rand, drop = FALSE],
    Lipid      = Lipid_full[, lip_rand, drop = FALSE]
  )
  
  Y_train_perm <- factor(c(
    rep("YNG", sum(rownames(RNA_train) %in% yng_sed_samples)),
    rep("OLD", sum(rownames(RNA_train) %in% old_sed_samples))
  ), levels = c("YNG","OLD"))
  
  design_perm <- matrix(perm_DS, 3, 3,
                        dimnames = list(names(X_train_perm), names(X_train_perm)))
  diag(design_perm) <- 0
  
  # Train on random features (keepX = all, since we already subsetted)
  fit_perm <- tryCatch({
    block.splsda(X_train_perm, Y_train_perm, ncomp = 1,
                 keepX = list(RNA = perm_RNA_keep,
                              Metabolite = perm_Met_keep,
                              Lipid = perm_Lip_keep),
                 design = design_perm)
  }, error = function(e) NULL)
  
  if (is.null(fit_perm)) return(NA_real_)
  
  # Project all samples
  proj_perm <- tryCatch(predict(fit_perm, newdata = X_full_perm), error = function(e) NULL)
  if (is.null(proj_perm)) return(NA_real_)
  
  scores_perm <- lapply(proj_perm$variates, function(x) x[, 1])
  comp1_perm  <- rowMeans(do.call(cbind, scores_perm))
  
  # Compute age separation (|Old_Sed mean - Young_Sed mean|)
  yng_idx <- names(comp1_perm) %in% yng_sed_samples
  old_idx <- names(comp1_perm) %in% old_sed_samples
  
  abs(mean(comp1_perm[old_idx]) - mean(comp1_perm[yng_idx]))
})

# Remove any failed permutations
perm_results <- perm_results[!is.na(perm_results)]
cat("Successful permutations:", length(perm_results), "\n")

#--- Observed age separation from real DIABLO ---#
real_df <- grid_results$sample_scores %>%
  filter(RNA_keep == perm_RNA_keep,
         Met_keep == perm_Met_keep,
         Lip_keep == perm_Lip_keep,
         DesignStrength == perm_DS)

real_age_sep <- abs(
  mean(real_df$Comp1[real_df$Group == "Old_SedVeh"]) -
    mean(real_df$Comp1[real_df$Group == "Young_Sed"])
)

# Empirical p-value
perm_p <- mean(perm_results >= real_age_sep)
cat("\nObserved age separation:", round(real_age_sep, 3), "\n")
cat("Permutation mean:",       round(mean(perm_results), 3), "\n")
cat("Permutation SD:",         round(sd(perm_results), 3), "\n")
cat("Permutation p-value:",    perm_p, "\n")
cat("(proportion of random selections with equal or greater separation)\n")

#--- Plot ---#
p_perm <- ggplot(data.frame(sep = perm_results), aes(x = sep)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(xintercept = real_age_sep, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
           x = real_age_sep, y = Inf, vjust = 2, hjust = -0.1,
           label = paste0("Observed = ", round(real_age_sep, 2),
                          "\np = ", format(perm_p, digits = 3)),
           color = "red", fontface = "bold", size = 4) +
  labs(title = paste0("Permutation Test: Random Feature Selection (n = ", length(perm_results), ")"),
       subtitle = paste0("RNA=", perm_RNA_keep, " / Met=", perm_Met_keep, " / Lip=", perm_Lip_keep),
       x = "Age Separation (|Old_Sed - Young_Sed| on Comp1)",
       y = "Count") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(p_perm)

ggsave("DIABLO_permutation_test.pdf", p_perm, width = 8, height = 5)
saveRDS(perm_results, "DIABLO_perm_results.RDS")
#############################################################################
#############################################################################
#
#predictive power of different omics
#-----------------------------------#

#===========================================================================#
# Omics Contribution Analysis: Do metabolites/lipids add predictive power?
#===========================================================================#

library(dplyr); library(tidyr); library(ggplot2)

# Load metric table (already in R, or reload)
# metric_table_with_tukey <- read.csv("metric_table_with_tukey.csv")

# Add derived columns
metric_table_with_tukey <- metric_table_with_tukey %>%
  mutate(
    FRap_vs_PwrVeh_diff = Old_PwrFRap - Old_PwrVeh,
    PwrVeh_vs_SedVeh_diff = Old_SedVeh - Old_PwrVeh   # restoration effect
  )

# Focus on DS=1.0 for clean interpretation
ds1 <- metric_table_with_tukey %>% filter(DesignStrength == 1.0)

#===========================================================================#
# ANALYSIS 1: FRap vs PwrVeh (rapamycin interference)
#===========================================================================#

#--- 1A: Metabolite effect on FRap vs PwrVeh (hold RNA=100, Lip=20) ---#
met_effect_frap <- ds1 %>%
  filter(RNA_keep == 100, Lip_keep == 20) %>%
  dplyr::select(Met_keep, FRap_vs_PwrVeh_diff, `p_Old_PwrVeh-Old_PwrFRap`) %>%
  dplyr::rename(p_value = `p_Old_PwrVeh-Old_PwrFRap`) %>%
  arrange(Met_keep)

cat("\n=== Metabolite effect on FRap vs PwrVeh (RNA=100, Lip=20) ===\n")
print(met_effect_frap)

#--- 1B: Lipid effect on FRap vs PwrVeh (hold RNA=100, Met=20) ---#
lip_effect_frap <- ds1 %>%
  filter(RNA_keep == 100, Met_keep == 20) %>%
  dplyr::select(Lip_keep, FRap_vs_PwrVeh_diff, `p_Old_PwrVeh-Old_PwrFRap`) %>%
  dplyr::rename(p_value = `p_Old_PwrVeh-Old_PwrFRap`) %>%
  arrange(Lip_keep)

cat("\n=== Lipid effect on FRap vs PwrVeh (RNA=100, Met=20) ===\n")
print(lip_effect_frap)

#--- 1C: Plot metabolite effect ---#
p_met_frap <- ggplot(met_effect_frap, aes(x = Met_keep, y = FRap_vs_PwrVeh_diff)) +
  geom_point(size = 4, color = "#ff7f00") +
  geom_line(color = "#ff7f00", linewidth = 1) +
  geom_text(aes(label = paste0("p=", round(p_value, 3))),
            vjust = -1.2, size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Metabolite Contribution to FRap vs PwrVeh Separation",
       subtitle = "RNA=100, Lipids=20 held constant",
       x = "Number of Metabolite Features",
       y = "Mean Comp1 Difference (FRap - PwrVeh)") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, max(met_effect_frap$FRap_vs_PwrVeh_diff) * 1.3))

#--- 1D: Plot lipid effect ---#
p_lip_frap <- ggplot(lip_effect_frap, aes(x = Lip_keep, y = FRap_vs_PwrVeh_diff)) +
  geom_point(size = 4, color = "#2ca02c") +
  geom_line(color = "#2ca02c", linewidth = 1) +
  geom_text(aes(label = paste0("p=", round(p_value, 3))),
            vjust = -1.2, size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Lipid Contribution to FRap vs PwrVeh Separation",
       subtitle = "RNA=100, Metabolites=20 held constant",
       x = "Number of Lipid Features",
       y = "Mean Comp1 Difference (FRap - PwrVeh)") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, max(lip_effect_frap$FRap_vs_PwrVeh_diff) * 1.3))

print(p_met_frap)
print(p_lip_frap)

#===========================================================================#
# ANALYSIS 2: PwrVeh vs SedVeh (exercise restoration)
#===========================================================================#

#--- 2A: Metabolite effect on PwrVeh vs SedVeh (hold RNA=100, Lip=20) ---#
met_effect_restore <- ds1 %>%
  filter(RNA_keep == 100, Lip_keep == 20) %>%
  dplyr::select(Met_keep, PwrVeh_vs_SedVeh_diff, `p_Old_SedVeh-Old_PwrVeh`) %>%
  dplyr::rename(p_value = `p_Old_SedVeh-Old_PwrVeh`) %>%
  arrange(Met_keep)

cat("\n=== Metabolite effect on PwrVeh vs SedVeh (RNA=100, Lip=20) ===\n")
print(met_effect_restore)

#--- 2B: Lipid effect on PwrVeh vs SedVeh (hold RNA=100, Met=20) ---#
lip_effect_restore <- ds1 %>%
  filter(RNA_keep == 100, Met_keep == 20) %>%
  dplyr::select(Lip_keep, PwrVeh_vs_SedVeh_diff, `p_Old_SedVeh-Old_PwrVeh`) %>%
  dplyr::rename(p_value = `p_Old_SedVeh-Old_PwrVeh`) %>%
  arrange(Lip_keep)

cat("\n=== Lipid effect on PwrVeh vs SedVeh (RNA=100, Met=20) ===\n")
print(lip_effect_restore)

#--- 2C: Plot metabolite effect on restoration ---#
p_met_restore <- ggplot(met_effect_restore, aes(x = Met_keep, y = PwrVeh_vs_SedVeh_diff)) +
  geom_point(size = 4, color = "#ff7f00") +
  geom_line(color = "#ff7f00", linewidth = 1) +
  geom_text(aes(label = paste0("p=", round(p_value, 4))),
            vjust = -1.2, size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Metabolite Contribution to Exercise Restoration (PwrVeh vs SedVeh)",
       subtitle = "RNA=100, Lipids=20 held constant",
       x = "Number of Metabolite Features",
       y = "Mean Comp1 Difference (SedVeh - PwrVeh)") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, max(met_effect_restore$PwrVeh_vs_SedVeh_diff) * 1.3))

#--- 2D: Plot lipid effect on restoration ---#
p_lip_restore <- ggplot(lip_effect_restore, aes(x = Lip_keep, y = PwrVeh_vs_SedVeh_diff)) +
  geom_point(size = 4, color = "#2ca02c") +
  geom_line(color = "#2ca02c", linewidth = 1) +
  geom_text(aes(label = paste0("p=", round(p_value, 4))),
            vjust = -1.2, size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Lipid Contribution to Exercise Restoration (PwrVeh vs SedVeh)",
       subtitle = "RNA=100, Metabolites=20 held constant",
       x = "Number of Lipid Features",
       y = "Mean Comp1 Difference (SedVeh - PwrVeh)") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, max(lip_effect_restore$PwrVeh_vs_SedVeh_diff) * 1.3))

print(p_met_restore)
print(p_lip_restore)

#===========================================================================#
# COMBINED: Both comparisons side by side in one table
#===========================================================================#

# Metabolite contribution summary
cat("\n\n====== METABOLITE CONTRIBUTION SUMMARY (RNA=100, Lip=20, DS=1.0) ======\n")
met_combined <- met_effect_frap %>%
  dplyr::rename(FRap_PwrVeh_diff = FRap_vs_PwrVeh_diff,
                FRap_PwrVeh_p = p_value) %>%
  left_join(
    met_effect_restore %>%
      dplyr::rename(SedVeh_PwrVeh_diff = PwrVeh_vs_SedVeh_diff,
                    SedVeh_PwrVeh_p = p_value),
    by = "Met_keep"
  )
print(met_combined)

# Lipid contribution summary
cat("\n====== LIPID CONTRIBUTION SUMMARY (RNA=100, Met=20, DS=1.0) ======\n")
lip_combined <- lip_effect_frap %>%
  dplyr::rename(FRap_PwrVeh_diff = FRap_vs_PwrVeh_diff,
                FRap_PwrVeh_p = p_value) %>%
  left_join(
    lip_effect_restore %>%
      dplyr::rename(SedVeh_PwrVeh_diff = PwrVeh_vs_SedVeh_diff,
                    SedVeh_PwrVeh_p = p_value),
    by = "Lip_keep"
  )
print(lip_combined)

# Save combined tables
write.csv(met_combined, "metabolite_contribution_analysis.csv", row.names = FALSE)
write.csv(lip_combined, "lipid_contribution_analysis.csv", row.names = FALSE)

# Save plots
ggsave("met_contribution_frap_vs_pwrveh.pdf", p_met_frap, width = 7, height = 5)
ggsave("lip_contribution_frap_vs_pwrveh.pdf", p_lip_frap, width = 7, height = 5)
ggsave("met_contribution_restoration.pdf", p_met_restore, width = 7, height = 5)
ggsave("lip_contribution_restoration.pdf", p_lip_restore, width = 7, height = 5)
###############################################################################
###############################################################################
#
#Each omic's contribution:
#----------------------------#

#===========================================================================#
# Single-block analysis: Does each omics layer independently 
# capture the rapamycin effect?
#===========================================================================#

library(mixOmics); library(dplyr); library(ggplot2); library(tidyr)

run_single_block <- function(train_mat, full_mat, block_name,
                             yng_sed_samples, old_sed_samples,
                             old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
                             keep_n) {
  
  Y_train <- factor(c(
    rep("YNG", sum(rownames(train_mat) %in% yng_sed_samples)),
    rep("OLD", sum(rownames(train_mat) %in% old_sed_samples))
  ), levels = c("YNG","OLD"))
  
  # Single-block sPLS-DA
  fit <- splsda(train_mat, Y_train, ncomp = 1, keepX = keep_n)
  
  # Project all samples using loadings directly
  loadings_vec <- fit$loadings$X[, 1]
  selected_feats <- names(loadings_vec[loadings_vec != 0])
  comp1 <- as.matrix(full_mat[, selected_feats, drop = FALSE]) %*% loadings_vec[selected_feats]
  comp1 <- comp1[, 1]  # matrix to vector
  names(comp1) <- rownames(full_mat)
  
  # Harmonize direction so higher = older
  yng_m <- mean(comp1[names(comp1) %in% yng_sed_samples])
  old_m <- mean(comp1[names(comp1) %in% old_sed_samples])
  if (old_m < yng_m) comp1 <- -comp1
  
  data.frame(
    Sample = names(comp1),
    Comp1  = comp1,
    Block  = block_name,
    keep_n = keep_n,
    Group  = group_of(names(comp1))
  )
}

# Run each block in isolation
rna_only <- run_single_block(RNA_train, RNA_full, "RNA",
                             yng_sed_samples, old_sed_samples,
                             old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
                             keep_n = 100)

met_only <- run_single_block(Metabo_train, Metabo_full, "Metabolite",
                             yng_sed_samples, old_sed_samples,
                             old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
                             keep_n = 20)

lip_only <- run_single_block(Lipid_train, Lipid_full, "Lipid",
                             yng_sed_samples, old_sed_samples,
                             old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples,
                             keep_n = 50)

# Combine
single_block_df <- bind_rows(rna_only, met_only, lip_only) %>%
  mutate(Group = factor(Group,
                        levels = c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap")))

# Group means per block
group_means <- single_block_df %>%
  group_by(Block, Group) %>%
  summarise(mean_Comp1 = mean(Comp1), .groups = "drop")

cat("\n=== Group means by single omics block ===\n")
print(group_means %>% pivot_wider(names_from = Group, values_from = mean_Comp1))

# Key comparisons per block
cat("\n=== FRap vs PwrVeh separation per block ===\n")
frap_pwrveh <- single_block_df %>%
  filter(Group %in% c("Old_PwrVeh","Old_PwrFRap")) %>%
  group_by(Block) %>%
  summarise(
    mean_PwrVeh = mean(Comp1[Group == "Old_PwrVeh"]),
    mean_FRap   = mean(Comp1[Group == "Old_PwrFRap"]),
    diff        = mean_FRap - mean_PwrVeh,
    p_value     = t.test(Comp1 ~ Group)$p.value,
    .groups     = "drop"
  )
print(frap_pwrveh)

cat("\n=== PwrVeh vs SedVeh separation per block ===\n")
restore <- single_block_df %>%
  filter(Group %in% c("Old_PwrVeh","Old_SedVeh")) %>%
  group_by(Block) %>%
  summarise(
    mean_PwrVeh = mean(Comp1[Group == "Old_PwrVeh"]),
    mean_SedVeh = mean(Comp1[Group == "Old_SedVeh"]),
    diff        = mean_SedVeh - mean_PwrVeh,
    p_value     = t.test(Comp1 ~ Group)$p.value,
    .groups     = "drop"
  )
print(restore)

# ANOVA + Tukey per block
cat("\n=== Per-block ANOVA + Tukey ===\n")
for (blk in c("RNA","Metabolite","Lipid")) {
  cat("\n---", blk, "---\n")
  sub <- single_block_df %>% filter(Block == blk)
  fit <- aov(Comp1 ~ Group, data = sub)
  cat("ANOVA p =", summary(fit)[[1]]$`Pr(>F)`[1], "\n")
  tuk <- TukeyHSD(fit)
  # Print just the key comparisons
  key_comps <- c("Old_PwrFRap-Old_PwrVeh", "Old_SedVeh-Old_PwrVeh", "Old_SedVeh-Young_Sed")
  print(tuk$Group[rownames(tuk$Group) %in% key_comps, ])
}

# Plot: faceted boxplot
p_blocks <- ggplot(single_block_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_wrap(~Block, scales = "free_y") +
  scale_fill_manual(values = c(
    "Young_Sed" = "#1F78B4", "Old_SedVeh" = "#E31A1C",
    "Old_PwrVeh" = "#33A02C", "Old_PwrIRap" = "#FB9A99",
    "Old_PwrFRap" = "#FDBF6F"
  )) +
  labs(title = "Single-Block Aging Axis: Each Omics Layer in Isolation",
       x = "", y = "sPLS-DA Comp1 Score") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

print(p_blocks)
ggsave("single_block_aging_axis.pdf", p_blocks, width = 12, height = 5)

# Save for Prism
write.csv(single_block_df, "single_block_scores.csv", row.names = FALSE)

# Wide format per block for Prism
for (blk in c("RNA","Metabolite","Lipid")) {
  wide <- single_block_df %>%
    filter(Block == blk) %>%
    dplyr::select(Sample, Group, Comp1) %>%
    group_by(Group) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = Group, values_from = Comp1) %>%
    dplyr::select(-Sample, -row)
  write.csv(wide, paste0("single_block_", tolower(blk), "_prism.csv"), row.names = FALSE)
}
####################################################################################

oldsedveh_v_yngsedveh_genes_06mar

#===========================================================================#
# logFC Barplots for DIABLO Features (100 RNA, 20 Met, 50 Lip)
#===========================================================================#

library(dplyr); library(ggplot2); library(readr); library(purrr)

# Get selected features from your final DIABLO model (run this if not already done)
rna_sel_final <- selectVar(diablo_final, block = "RNA", comp = 1)
met_sel_final <- selectVar(diablo_final, block = "Metabolite", comp = 1)
lip_sel_final <- selectVar(diablo_final, block = "Lipid", comp = 1)

rna_features_final <- rna_sel_final$RNA$name
met_features_final <- met_sel_final$Metabolite$name
lip_features_final <- lip_sel_final$Lipid$name

#===========================================================================#
# 1) RNA plot data — from your existing DE dataframe
#===========================================================================#

rna_plot_df <- oldsedveh_v_yngsedveh_genes_06mar %>%
  filter(Symbol %in% rna_features_final) %>%
  arrange(logFC) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol),
         neg_logFDR = pmin(-log10(FDR), 5))

#===========================================================================#
# 2) Lipid plot data — compute OS vs YS stats
#===========================================================================#

lip_raw <- read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)
lipid_cols <- setdiff(names(lip_raw), c("Sample","Group"))

# Impute zeros, log2
lip_imp <- lip_raw
for (lip in lipid_cols) {
  x <- lip_imp[[lip]]
  x[x == 0 | is.na(x)] <- suppressWarnings(min(x[x > 0], na.rm = TRUE)) / 2
  lip_imp[[lip]] <- log2(x)
}

lipid_os_ys <- lip_imp %>% filter(Group %in% c("OS","YS"))

lipid_stats <- purrr::map_dfr(lipid_cols, function(lip) {
  vals <- lipid_os_ys[[lip]]
  grp  <- lipid_os_ys$Group
  tt   <- t.test(vals ~ grp)
  tibble(
    Lipid            = lip,
    Log2_FC_OS_vs_YS = mean(vals[grp == "OS"], na.rm = TRUE) - mean(vals[grp == "YS"], na.rm = TRUE),
    P_value          = tt$p.value
  )
}) %>%
  mutate(FDR = p.adjust(P_value, method = "fdr"))

lipid_plot_df <- lipid_stats %>%
  filter(Lipid %in% lip_features_final) %>%
  arrange(Log2_FC_OS_vs_YS) %>%
  mutate(Lipid = factor(Lipid, levels = Lipid),
         neg_logFDR = pmin(-log10(FDR), 5))

#===========================================================================#
# 3) Metabolite plot data — compute OS vs YS stats
#===========================================================================#

metabo_raw <- read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)
metabo_os_ys <- metabo_raw %>% filter(Group %in% c("OS","YS"))
metabo_cols <- setdiff(names(metabo_os_ys), c("Sample","Group"))

metabo_stats <- metabo_os_ys %>%
  pivot_longer(cols = all_of(metabo_cols), names_to = "Metabolite", values_to = "Abundance") %>%
  drop_na(Abundance) %>%
  group_by(Metabolite) %>%
  summarise(
    p_value = tryCatch(t.test(Abundance ~ Group)$p.value, error = function(e) NA_real_),
    mean_OS = mean(Abundance[Group == "OS"], na.rm = TRUE),
    mean_YS = mean(Abundance[Group == "YS"], na.rm = TRUE),
    log2FC  = log2(mean_OS + 1e-6) - log2(mean_YS + 1e-6),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

metabolite_plot_df <- metabo_stats %>%
  filter(Metabolite %in% met_features_final) %>%
  arrange(log2FC) %>%
  mutate(Metabolite = factor(Metabolite, levels = Metabolite),
         neg_logFDR = pmin(-log10(FDR), 5))

#===========================================================================#
# 4) Cap and set shared limits
#===========================================================================#

common_limits <- c(0, 5)

#===========================================================================#
# 5) RNA plot — "mako" (brilliant blue)
#===========================================================================#

diablo_rna_feature_plot <- ggplot(rna_plot_df, aes(x = logFC, y = Symbol, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = "Gene Symbol",
    title = "DIABLO RNA Features: Differential Expression"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#===========================================================================#
# 6) Lipid plot — "inferno" (orange)
#===========================================================================#

diablo_lipid_feature_plot <- ggplot(lipid_plot_df, aes(x = Log2_FC_OS_vs_YS, y = Lipid, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "inferno",
    direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 4)) +
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = "Lipid Species",
    title = "DIABLO Lipid Features: Differential Abundance"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#===========================================================================#
# 7) Metabolite plot — "rocket" (yellow/warm)
#===========================================================================#

diablo_metabolite_feature_plot <- ggplot(metabolite_plot_df, aes(x = log2FC, y = Metabolite, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "rocket",
    direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = "Metabolite",
    title = "DIABLO Metabolite Features: Differential Abundance"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#===========================================================================#
# 8) Print and save
#===========================================================================#

print(diablo_rna_feature_plot)
print(diablo_lipid_feature_plot)
print(diablo_metabolite_feature_plot)

pdf("diablo_rna_feature_plot_v2.pdf", width = 6, height = 10)
print(diablo_rna_feature_plot)
dev.off()

pdf("diablo_lipid_feature_plot_v2.pdf", width = 6, height = 6)
print(diablo_lipid_feature_plot)
dev.off()

pdf("diablo_metabolite_feature_plot_v2.pdf", width = 6, height = 3)
print(diablo_metabolite_feature_plot)
dev.off()
########################################################################

validated_aging_genes
signature59
metab_sig_df
#======================================#

#===========================================================================#
# Find overlaps between DIABLO features and aging signatures
#===========================================================================#

# RNA overlap
rna_overlap <- intersect(rna_features_final, validated_aging_genes$Symbol)
cat("RNA overlap:", length(rna_overlap), "of", length(rna_features_final), "DIABLO genes\n")
print(rna_overlap)

# Lipid overlap
lip_overlap <- intersect(lip_features_final, signature59)
cat("\nLipid overlap:", length(lip_overlap), "of", length(lip_features_final), "DIABLO lipids\n")
print(lip_overlap)

# Metabolite overlap
met_overlap <- intersect(met_features_final, metab_sig_df$Metabolite)
cat("\nMetabolite overlap:", length(met_overlap), "of", length(met_features_final), "DIABLO metabolites\n")
print(met_overlap)

#===========================================================================#
# Add overlap flag to plot dataframes
#===========================================================================#

rna_plot_df <- rna_plot_df %>%
  mutate(in_signature = Symbol %in% rna_overlap)

lipid_plot_df <- lipid_plot_df %>%
  mutate(in_signature = Lipid %in% lip_overlap)

metabolite_plot_df <- metabolite_plot_df %>%
  mutate(in_signature = Metabolite %in% met_overlap)

#===========================================================================#
# Plots with highlighted overlapping features
#===========================================================================#

common_limits <- c(0, 5)

# RNA plot
diablo_rna_feature_plot <- ggplot(rna_plot_df, aes(x = logFC, y = Symbol, fill = neg_logFDR)) +
  geom_col(aes(color = in_signature, linewidth = in_signature)) +
  scale_color_manual(values = c("FALSE" = NA, "TRUE" = "black"), guide = "none") +
  scale_linewidth_manual(values = c("FALSE" = 0, "TRUE" = 0.8), guide = "none") +
  scale_fill_viridis_c(
    option = "mako", direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  theme_classic(base_size = 12) +
  labs(x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
       y = "Gene Symbol",
       title = "DIABLO RNA Features") +
  theme(
    axis.text.y = element_text(
      size = 7,
      face = ifelse(levels(rna_plot_df$Symbol) %in% rna_overlap, "bold", "plain")
    ),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Lipid plot
diablo_lipid_feature_plot <- ggplot(lipid_plot_df, aes(x = Log2_FC_OS_vs_YS, y = Lipid, fill = neg_logFDR)) +
  geom_col(aes(color = in_signature, linewidth = in_signature)) +
  scale_color_manual(values = c("FALSE" = NA, "TRUE" = "black"), guide = "none") +
  scale_linewidth_manual(values = c("FALSE" = 0, "TRUE" = 0.8), guide = "none") +
  scale_fill_viridis_c(
    option = "inferno", direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 4)) +
  theme_classic(base_size = 12) +
  labs(x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
       y = "Lipid Species",
       title = "DIABLO Lipid Features") +
  theme(
    axis.text.y = element_text(
      size = 7,
      face = ifelse(levels(lipid_plot_df$Lipid) %in% lip_overlap, "bold", "plain")
    ),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Metabolite plot
diablo_metabolite_feature_plot <- ggplot(metabolite_plot_df, aes(x = log2FC, y = Metabolite, fill = neg_logFDR)) +
  geom_col(aes(color = in_signature, linewidth = in_signature)) +
  scale_color_manual(values = c("FALSE" = NA, "TRUE" = "black"), guide = "none") +
  scale_linewidth_manual(values = c("FALSE" = 0, "TRUE" = 0.8), guide = "none") +
  scale_fill_viridis_c(
    option = "rocket", direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  theme_classic(base_size = 12) +
  labs(x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
       y = "Metabolite",
       title = "DIABLO Metabolite Features") +
  theme(
    axis.text.y = element_text(
      size = 9,
      face = ifelse(levels(metabolite_plot_df$Metabolite) %in% met_overlap, "bold", "plain")
    ),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(diablo_rna_feature_plot)
print(diablo_lipid_feature_plot)
print(diablo_metabolite_feature_plot)

pdf("diablo_rna_feature_plot_v2.pdf", width = 12, height = 20)
print(diablo_rna_feature_plot)
dev.off()

pdf("diablo_lipid_feature_plot_v2.pdf", width = 12, height = 10)
print(diablo_lipid_feature_plot)
dev.off()

pdf("diablo_metabolite_feature_plot_v2.pdf", width = 12, height = 4)
print(diablo_metabolite_feature_plot)
dev.off()
##############################################################
#################################################################

#===========================================================================#
# Correlate NEW DIABLO Scores (100/20/50) with Physiological Readouts
#===========================================================================#

library(dplyr); library(tidyr); library(ggplot2)

# 1) Load and reshape physiology data
triceps_physio <- read.csv("Triceps_Physio_Data_v2.csv")

physio_long <- as.data.frame(t(triceps_physio[-1]))
colnames(physio_long) <- triceps_physio$Tube.code
physio_long$Sample <- rownames(physio_long)
rownames(physio_long) <- NULL
physio_long <- physio_long %>%
  mutate(across(-Sample, as.numeric))

# 2) Merge with NEW DIABLO scores (diablo_scores_df from 100/20/50 model)
merged <- diablo_scores_df %>%
  dplyr::select(Sample, Comp1_centered, Group) %>%
  left_join(physio_long, by = "Sample") %>%
  filter(Sample %in% physio_long$Sample)

cat("Merged samples:", nrow(merged), "\n")

# 3) Subset to old groups only
old_only <- merged %>% filter(Group != "Young_Sed")

# 4) Spearman correlations — old groups only
physio_vars <- old_only %>%
  dplyr::select(where(is.numeric), -Comp1_centered) %>%
  names()

old_cor_all <- purrr::map_dfr(physio_vars, function(v) {
  x <- old_only$Comp1_centered
  y <- old_only[[v]]
  ok <- complete.cases(x, y)
  if (sum(ok) >= 4) {
    ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
    data.frame(var = v, r = ct$estimate, p = ct$p.value, n = sum(ok))
  } else {
    data.frame(var = v, r = NA, p = NA, n = sum(ok))
  }
}) %>%
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  arrange(FDR)

cat("\n=== Spearman correlations (old groups only) ===\n")
print(old_cor_all, digits = 3)

# 5) Correlation barplot
physio_diablo_corr_barplot <- ggplot(old_cor_all,
                                     aes(x = reorder(var, r), y = r, fill = p < 0.05)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70")) +
  labs(
    x = "",
    y = "Spearman r (DIABLO Comp1 vs Physiology)",
    title = "Physiologic correlates of molecular aging axis (old groups only)",
    fill = "p < 0.05"
  ) +
  geom_text(aes(label = sprintf("n=%d", n)),
            hjust = -0.2, size = 3, color = "black")

print(physio_diablo_corr_barplot)

pdf("physio_diablo_corr_barplot_v2.pdf", width = 8, height = 6)
print(physio_diablo_corr_barplot)
dev.off()

write.csv(old_cor_all, "DIABLO_physio_correlations_v2.csv", row.names = FALSE)

# 6) Add Group3 column for scatter plots (combine rapa groups)
old_plot_df <- old_only %>%
  filter(Sample != "T_25") %>%
  mutate(Group3 = case_when(
    Group %in% c("Old_PwrIRap", "Old_PwrFRap") ~ "Old_PwrRapa",
    TRUE ~ as.character(Group)
  ))

# 7) GTT scatter plot
fit_gtt <- lm(`Post GTT (AOC)` ~ Comp1_centered, data = old_plot_df)
r2_gtt <- round(summary(fit_gtt)$r.squared, 3)
pval_gtt <- signif(coef(summary(fit_gtt))[2, 4], 3)

GTT_diablo_sctrplot <- ggplot(old_plot_df,
                              aes(x = Comp1_centered, y = `Post GTT (AOC)`, color = Group3)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15, color = NA, level = 0.68) +
  geom_smooth(aes(x = Comp1_centered, y = `Post GTT (AOC)`),
              method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 1.2) +
  annotate("text",
           x = min(old_plot_df$Comp1_centered, na.rm = TRUE),
           y = max(old_plot_df$`Post GTT (AOC)`, na.rm = TRUE),
           label = sprintf("R² = %.3f, p = %.3g", r2_gtt, pval_gtt),
           hjust = 0, vjust = 1.2, size = 5, fontface = "italic") +
  scale_color_manual(values = c("Old_SedVeh" = "#2c3e50", "Old_PwrVeh" = "#27ae60", "Old_PwrRapa" = "#e74c3c")) +
  scale_fill_manual(values = c("Old_SedVeh" = "#2c3e50", "Old_PwrVeh" = "#27ae60", "Old_PwrRapa" = "#e74c3c")) +
  labs(x = "DIABLO Comp1 (centered, higher = older)",
       y = "Glucose Tolerance Test AUC",
       title = "Molecular aging axis vs glucose tolerance in old mice") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank())

print(GTT_diablo_sctrplot)
ggsave("GTT_diablo_sctrplot_v2.pdf", GTT_diablo_sctrplot, width = 8, height = 6)

# 8) ITT scatter plot
fit_itt <- lm(`Post Insulin Sensitivity (AUC)` ~ Comp1_centered, data = old_plot_df)
r2_itt <- round(summary(fit_itt)$r.squared, 3)
pval_itt <- signif(coef(summary(fit_itt))[2, 4], 3)

ITT_diablo_sctrplot <- ggplot(old_plot_df,
                              aes(x = Comp1_centered, y = `Post Insulin Sensitivity (AUC)`, color = Group3)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15, color = NA, level = 0.68) +
  geom_smooth(aes(x = Comp1_centered, y = `Post Insulin Sensitivity (AUC)`),
              method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 1.2) +
  annotate("text",
           x = min(old_plot_df$Comp1_centered, na.rm = TRUE),
           y = max(old_plot_df$`Post Insulin Sensitivity (AUC)`, na.rm = TRUE),
           label = sprintf("R² = %.3f, p = %.3g", r2_itt, pval_itt),
           hjust = 0, vjust = 1.2, size = 5, fontface = "italic") +
  scale_color_manual(values = c("Old_SedVeh" = "#2c3e50", "Old_PwrVeh" = "#27ae60", "Old_PwrRapa" = "#e74c3c")) +
  scale_fill_manual(values = c("Old_SedVeh" = "#2c3e50", "Old_PwrVeh" = "#27ae60", "Old_PwrRapa" = "#e74c3c")) +
  labs(x = "DIABLO Comp1 (centered, higher = older)",
       y = "Insulin Sensitivity (AUC)",
       title = "Molecular aging axis vs insulin sensitivity in old mice") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank())

print(ITT_diablo_sctrplot)
ggsave("ITT_diablo_sctrplot_v2.pdf", ITT_diablo_sctrplot, width = 8, height = 6)

# 9) GXT scatter plot
fit_gxt <- lm(`Post GXT time (sec) ` ~ Comp1_centered, data = old_plot_df)
r2_gxt <- round(summary(fit_gxt)$r.squared, 3)
pval_gxt <- signif(coef(summary(fit_gxt))[2, 4], 3)

GXT_diablo_sctrplot <- ggplot(old_plot_df,
                              aes(x = Comp1_centered, y = `Post GXT time (sec) `, color = Group3)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15, color = NA, level = 0.68) +
  geom_smooth(aes(x = Comp1_centered, y = `Post GXT time (sec) `),
              method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 1.2) +
  annotate("text",
           x = min(old_plot_df$Comp1_centered, na.rm = TRUE),
           y = max(old_plot_df$`Post GXT time (sec) `, na.rm = TRUE),
           label = sprintf("R² = %.3f, p = %.3g", r2_gxt, pval_gxt),
           hjust = 0, vjust = 1.2, size = 5, fontface = "italic") +
  scale_color_manual(values = c("Old_SedVeh" = "#2c3e50", "Old_PwrVeh" = "#27ae60", "Old_PwrRapa" = "#e74c3c")) +
  scale_fill_manual(values = c("Old_SedVeh" = "#2c3e50", "Old_PwrVeh" = "#27ae60", "Old_PwrRapa" = "#e74c3c")) +
  labs(x = "DIABLO Comp1 (centered, higher = older)",
       y = "GXT Time (sec)",
       title = "Molecular aging axis vs exercise tolerance in old mice") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank())

print(GXT_diablo_sctrplot)
ggsave("GXT_diablo_sctrplot_v2.pdf", GXT_diablo_sctrplot, width = 8, height = 6)


#

