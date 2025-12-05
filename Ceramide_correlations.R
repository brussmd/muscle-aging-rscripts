# ============================ SETUP ==========================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(FactoMineR)
  library(factoextra)
})

# -------------------- Groups (including YS) ----------------------------------
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")  # YS
old_sedveh_samples  <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwrveh_samples  <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

group_assign <- tibble(
  Sample = c(yng_sed_samples,
             old_sedveh_samples, old_pwrveh_samples,
             old_pwrirap_samples, old_pwrfrap_samples),
  Group = c(rep("Young_Sed", length(yng_sed_samples)),
            rep("Old_SedVeh", length(old_sedveh_samples)),
            rep("Old_PwrVeh", length(old_pwrveh_samples)),
            rep("Old_PwrIRap", length(old_pwrirap_samples)),
            rep("Old_PwrFRap", length(old_pwrfrap_samples)))
)

# -------------------- Sample key --------------------------------------------
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_01","OFR1","T_04","OFR2","T_22","OFR3","T_24","OFR4","T_28","OFR5","T_34","OFR6","T_40","OFR7","T_44","OFR8",
  "T_02","OIR1","T_03","OIR2","T_14","OIR3","T_31","OIR4","T_33","OIR5","T_37","OIR6","T_47","OIR7","T_48","OIR8",
  "T_08","OV1","T_10","OV2","T_11","OV3","T_16","OV4","T_41","OV5","T_43","OV6","T_45","OV7","T_46","OV8",
  "T_06","OS1","T_09","OS2","T_13","OS3","T_19","OS4","T_21","OS5","T_25","OS6","T_29","OS7","T_39","OS8",
  "T_07","YV1","T_12","YV2","T_18","YV3","T_20","YV4","T_27","YV5","T_30","YV6","T_36","YV7","T_42","YV8",
  "T_05","YS1","T_15","YS2","T_17","YS3","T_23","YS4","T_26","YS5","T_32","YS6","T_35","YS7","T_38","YS8"
)

# ====================== LIPID DATA (FILTER OUT YV) ==========================
lipid_df <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE) %>%
  mutate(Sample = str_trim(Sample))

# log2-transform
lipid_cols <- setdiff(names(lipid_df), c("Sample","Group"))
lipid_df[ , lipid_cols] <- lapply(lipid_df[ , lipid_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})

# identify ceramide columns
cer_cols_idx <- grepl("^Cer", names(lipid_df))
stopifnot(sum(cer_cols_idx) > 0)

# filter out YV samples only
lipid_cer_long_all <- lipid_df %>%
  filter(!grepl("^YV", Sample)) %>%
  select(c("Sample", names(lipid_df)[cer_cols_idx])) %>%
  pivot_longer(-Sample, names_to = "Lipid", values_to = "log2_value") %>%
  left_join(sample_key, by = c("Sample" = "MetabolomicsName")) %>%
  mutate(Sample = TubeCode) %>%
  filter(!is.na(Sample)) %>%
  select(Sample, Lipid, log2_value)

# wide numeric matrix
cer_mat_all <- lipid_cer_long_all %>%
  pivot_wider(names_from = Lipid, values_from = log2_value) %>%
  column_to_rownames("Sample") %>%
  as.data.frame()
stopifnot(all(sapply(cer_mat_all, is.numeric)))

# ===================== PHYSIO DATA ==========================================
physio <- read.csv("Triceps_Physio_Data_v2.csv", check.names = FALSE)
names(physio) <- names(physio) %>% gsub("^T_6$", "T_06", .) %>% str_trim()

ins_df <- physio %>%
  filter(`Tube code` == "Post Insulin Sensitivity (AUC)") %>%
  pivot_longer(-`Tube code`, names_to = "Sample", values_to = "Post_Ins_AUC") %>%
  transmute(Sample = str_trim(Sample),
            Post_Ins_AUC = suppressWarnings(as.numeric(Post_Ins_AUC)))

# ===================== PCA (CERAMIDES → PC1) ================================
pca_all <- FactoMineR::PCA(cer_mat_all, graph = FALSE)
pc1_scores_all <- pca_all$ind$coord[, 1]
pc1_df <- tibble(Sample = names(pc1_scores_all), PC1 = as.numeric(pc1_scores_all))

# merge PC1 with physiology + group
dat_all <- ins_df %>%
  inner_join(pc1_df, by = "Sample") %>%
  left_join(group_assign, by = "Sample")

# ===================== CORRELATION ==========================================
ct <- cor.test(dat_all$PC1, dat_all$Post_Ins_AUC, method = "spearman")
rho  <- unname(ct$estimate)
pval <- ct$p.value
lm_fit <- lm(Post_Ins_AUC ~ PC1, data = dat_all)
r2 <- summary(lm_fit)$r.squared

cat(sprintf("Spearman rho = %.3f, p = %.4g | LM R^2 = %.3f\n", rho, pval, r2))

# ===================== PLOT ================================================
pal <- c(
  Young_Sed = "#2E8B57",
  Old_SedVeh = "#1F77B4",
  Old_PwrVeh = "#FF7F0E",
  Old_PwrIRap = "#9467BD",
  Old_PwrFRap = "#D62728"
)

annot_text <- sprintf("Spearman rho = %.2f\nLM R² = %.2f, p = %.3g", rho, r2, pval)

ggplot(dat_all, aes(x = PC1, y = Post_Ins_AUC, color = Group, fill = Group)) +
  stat_ellipse(type = "norm", level = 0.68, geom = "polygon", alpha = 0.15, color = NA, show.legend = FALSE) +
  geom_point(shape = 21, stroke = 2, size = 4, color = "white") +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = pal, na.value = "grey50") +
  scale_fill_manual(values = pal, na.value = "grey80") +
  theme_classic(base_size = 14) +
  labs(
    title = "Ceramide PC1 vs Post Insulin Sensitivity (AUC)",
    subtitle = "YS included, YV excluded",
    x = "Ceramide PC1 (first principal component across ceramides)",
    y = "Post Insulin Sensitivity (AUC)",
    color = "Group", fill = "Group"
  ) +
  annotate("text", x = Inf, y = -Inf, label = annot_text,
           hjust = 1.05, vjust = -0.5, size = 4.2)

# ===================== PC1 LOADINGS (OPTIONAL) ==============================
pc1_loadings_all <- as.data.frame(pca_all$var$coord[, 1, drop = FALSE]) %>%
  tibble::rownames_to_column("Lipid") %>%
  rename(PC1_loading = Dim.1) %>%
  mutate(abs_loading = abs(PC1_loading)) %>%
  arrange(desc(abs_loading))
print(head(pc1_loadings_all, 15))


# ============================================================
# 🔍 CORRELATION WITHIN OLD_PwrFRap GROUP ONLY
# ============================================================
dat_frap <- dat_all %>%
  filter(Group == "Old_PwrFRap") %>%
  filter(is.finite(PC1) & is.finite(Post_Ins_AUC))

# Spearman correlation
ct_frap <- cor.test(dat_frap$PC1, dat_frap$Post_Ins_AUC, method = "spearman")
rho_frap  <- unname(ct_frap$estimate)
pval_frap <- ct_frap$p.value

# Linear model for visualization (R²)
lm_frap <- lm(Post_Ins_AUC ~ PC1, data = dat_frap)
r2_frap <- summary(lm_frap)$r.squared

cat(sprintf("\n--- OLD_PwrFRap ONLY ---\nSpearman rho = %.3f, p = %.4g | LM R² = %.3f\n",
            rho_frap, pval_frap, r2_frap))

# ============================================================
# 🎨 Visualization: OLD_PwrFRap only
# ============================================================
annot_text_frap <- sprintf("Spearman ρ = %.2f\nR² = %.2f, p = %.3g",
                           rho_frap, r2_frap, pval_frap)

cer_ins_plot <- ggplot(dat_frap, aes(x = PC1, y = Post_Ins_AUC)) +
  geom_point(size = 4, color = "#9E1F63", fill = "#9E1F63", alpha = 0.8, shape = 21) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  theme_classic(base_size = 14) +
  labs(
    title = "Ceramide PC1 vs Post Insulin Sensitivity (AUC)",
    subtitle = "Old_PwrFRap samples only",
    x = "Ceramide PC1 (first principal component)",
    y = "Post Insulin Sensitivity (AUC)"
  ) +
  annotate("text", x = Inf, y = -Inf, label = annot_text_frap,
           hjust = 1.05, vjust = -0.5, size = 4.5)

pdf("cer_ins_plot.pdf", width = 6, height = 4)
print(cer_ins_plot)
dev.off()

#====================================================================#

# ============================================================
# 🔍 CORRELATION: Combined Rapamycin Groups (IRap + FRap)
# ============================================================
dat_rapa <- dat_all %>%
  filter(Group %in% c("Old_PwrIRap", "Old_PwrFRap")) %>%
  filter(is.finite(PC1) & is.finite(Post_Ins_AUC))

# Spearman correlation
ct_rapa <- cor.test(dat_rapa$PC1, dat_rapa$Post_Ins_AUC, method = "spearman")
rho_rapa  <- unname(ct_rapa$estimate)
pval_rapa <- ct_rapa$p.value

# Linear model (for visualization and R²)
lm_rapa <- lm(Post_Ins_AUC ~ PC1, data = dat_rapa)
r2_rapa <- summary(lm_rapa)$r.squared

cat(sprintf("\n--- RAPAMYCIN GROUPS COMBINED ---\nSpearman rho = %.3f, p = %.4g | LM R² = %.3f\n",
            rho_rapa, pval_rapa, r2_rapa))

# ============================================================
# 🎨 Visualization: Combined Rapamycin Groups
# ============================================================
pal_rapa <- c(Old_PwrIRap = "#9467BD", Old_PwrFRap = "#D62728")

annot_text_rapa <- sprintf("Spearman ρ = %.2f\nR² = %.2f, p = %.3g",
                           rho_rapa, r2_rapa, pval_rapa)

ggplot(dat_rapa, aes(x = PC1, y = Post_Ins_AUC, color = Group, fill = Group)) +
  stat_ellipse(type = "norm", level = 0.68, geom = "polygon", alpha = 0.18, color = NA) +
  geom_point(shape = 21, stroke = 2, size = 4, color = "white") +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  scale_color_manual(values = pal_rapa) +
  scale_fill_manual(values = pal_rapa) +
  theme_classic(base_size = 14) +
  labs(
    title = "Ceramide PC1 vs Post Insulin Sensitivity (AUC)",
    subtitle = "Old_PwrIRap + Old_PwrFRap combined",
    x = "Ceramide PC1 (first principal component)",
    y = "Post Insulin Sensitivity (AUC)",
    color = "Group", fill = "Group"
  ) +
  annotate("text", x = Inf, y = -Inf, label = annot_text_rapa,
           hjust = 1.05, vjust = -0.5, size = 4.5)
#=========================================================================#

# ====================== CE ANALYSIS =========================================
# 🧭 Switch from ceramides (Cer) → cholesterol esters (CE)
# ============================================================================

# identify CE columns
ce_cols_idx <- grepl("^CE", names(lipid_df))
stopifnot(sum(ce_cols_idx) > 0)

# filter out YV samples only
lipid_ce_long_all <- lipid_df %>%
  filter(!grepl("^YV", Sample)) %>%
  select(c("Sample", names(lipid_df)[ce_cols_idx])) %>%
  pivot_longer(-Sample, names_to = "Lipid", values_to = "log2_value") %>%
  left_join(sample_key, by = c("Sample" = "MetabolomicsName")) %>%
  mutate(Sample = TubeCode) %>%
  filter(!is.na(Sample)) %>%
  select(Sample, Lipid, log2_value)

# wide numeric matrix
ce_mat_all <- lipid_ce_long_all %>%
  pivot_wider(names_from = Lipid, values_from = log2_value) %>%
  column_to_rownames("Sample") %>%
  as.data.frame()
stopifnot(all(sapply(ce_mat_all, is.numeric)))

# ===================== PCA (CE → PC1) ======================================
pca_ce <- FactoMineR::PCA(ce_mat_all, graph = FALSE)
pc1_scores_ce <- pca_ce$ind$coord[, 1]
pc1_ce_df <- tibble(Sample = names(pc1_scores_ce), PC1_CE = as.numeric(pc1_scores_ce))

# merge PC1 with physiology + group
dat_ce_all <- ins_df %>%
  inner_join(pc1_ce_df, by = "Sample") %>%
  left_join(group_assign, by = "Sample")

# ===================== CORRELATION (ALL GROUPS) ============================
ct_ce <- cor.test(dat_ce_all$PC1_CE, dat_ce_all$Post_Ins_AUC, method = "spearman")
rho_ce  <- unname(ct_ce$estimate)
pval_ce <- ct_ce$p.value
lm_ce <- lm(Post_Ins_AUC ~ PC1_CE, data = dat_ce_all)
r2_ce <- summary(lm_ce)$r.squared

cat(sprintf("CE PC1 vs Ins Sensitivity | Spearman rho = %.3f, p = %.4g | LM R² = %.3f\n",
            rho_ce, pval_ce, r2_ce))

# ===================== PLOT (ALL GROUPS) ===================================
pal <- c(
  Young_Sed = "#2E8B57",
  Old_SedVeh = "#1F77B4",
  Old_PwrVeh = "#FF7F0E",
  Old_PwrIRap = "#9467BD",
  Old_PwrFRap = "#D62728"
)

annot_text_ce <- sprintf("Spearman ρ = %.2f\nR² = %.2f, p = %.3g", rho_ce, r2_ce, pval_ce)

ggplot(dat_ce_all, aes(x = PC1_CE, y = Post_Ins_AUC, color = Group, fill = Group)) +
  stat_ellipse(type = "norm", level = 0.68, geom = "polygon", alpha = 0.15, color = NA) +
  geom_point(shape = 21, stroke = 2, size = 4, color = "white") +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = pal, na.value = "grey50") +
  scale_fill_manual(values = pal, na.value = "grey80") +
  theme_classic(base_size = 14) +
  labs(
    title = "Cholesteryl Ester PC1 vs Post Insulin Sensitivity (AUC)",
    subtitle = "YS included, YV excluded",
    x = "CE PC1 (first principal component across cholesteryl esters)",
    y = "Post Insulin Sensitivity (AUC)",
    color = "Group", fill = "Group"
  ) +
  annotate("text", x = Inf, y = -Inf, label = annot_text_ce,
           hjust = 1.05, vjust = -0.5, size = 4.2)

# ===================== PC1 LOADINGS (TOP CE SPECIES) =======================
pc1_loadings_ce <- as.data.frame(pca_ce$var$coord[, 1, drop = FALSE]) %>%
  tibble::rownames_to_column("Lipid") %>%
  rename(PC1_loading = Dim.1) %>%
  mutate(abs_loading = abs(PC1_loading)) %>%
  arrange(desc(abs_loading))
print(head(pc1_loadings_ce, 15))
#-----------------------------------------------------------#

# ============================================================
# 🧠 CE Analysis — Rapamycin Groups (Old_PwrIRap + Old_PwrFRap)
# ============================================================

# 1️⃣ Identify CE columns
ce_cols_idx <- grepl("^CE", names(lipid_df))
stopifnot(sum(ce_cols_idx) > 0)

# 2️⃣ Filter out YV samples, keep all others
lipid_ce_long_all <- lipid_df %>%
  filter(!grepl("^YV", Sample)) %>%
  select(c("Sample", names(lipid_df)[ce_cols_idx])) %>%
  pivot_longer(-Sample, names_to = "Lipid", values_to = "log2_value") %>%
  left_join(sample_key, by = c("Sample" = "MetabolomicsName")) %>%
  mutate(Sample = TubeCode) %>%
  filter(!is.na(Sample)) %>%
  select(Sample, Lipid, log2_value)

# 3️⃣ Restrict to rapamycin groups
rapa_samples <- c(old_pwrirap_samples, old_pwrfrap_samples)
lipid_ce_long_rapa <- lipid_ce_long_all %>%
  filter(Sample %in% rapa_samples)

# 4️⃣ Build wide numeric matrix
ce_mat_rapa <- lipid_ce_long_rapa %>%
  pivot_wider(names_from = Lipid, values_from = log2_value) %>%
  column_to_rownames("Sample") %>%
  as.data.frame()
stopifnot(all(sapply(ce_mat_rapa, is.numeric)))

# 5️⃣ PCA (CE → PC1)
pca_ce_rapa <- FactoMineR::PCA(ce_mat_rapa, graph = FALSE)
pc1_scores_ce_rapa <- pca_ce_rapa$ind$coord[, 1]
pc1_ce_rapa_df <- tibble(Sample = names(pc1_scores_ce_rapa), PC1_CE = as.numeric(pc1_scores_ce_rapa))

# 6️⃣ Merge with insulin sensitivity data + group
dat_ce_rapa <- ins_df %>%
  inner_join(pc1_ce_rapa_df, by = "Sample") %>%
  left_join(group_assign, by = "Sample") %>%
  filter(Group %in% c("Old_PwrIRap", "Old_PwrFRap")) %>%
  filter(is.finite(PC1_CE) & is.finite(Post_Ins_AUC))

# 7️⃣ Correlation
ct_ce_rapa <- cor.test(dat_ce_rapa$PC1_CE, dat_ce_rapa$Post_Ins_AUC, method = "spearman")
rho_ce_rapa  <- unname(ct_ce_rapa$estimate)
pval_ce_rapa <- ct_ce_rapa$p.value

lm_ce_rapa <- lm(Post_Ins_AUC ~ PC1_CE, data = dat_ce_rapa)
r2_ce_rapa <- summary(lm_ce_rapa)$r.squared

cat(sprintf("\n--- CE: RAPAMYCIN GROUPS COMBINED ---\nSpearman rho = %.3f, p = %.4g | LM R² = %.3f\n",
            rho_ce_rapa, pval_ce_rapa, r2_ce_rapa))

# 8️⃣ Visualization
pal_rapa <- c(Old_PwrIRap = "#9467BD", Old_PwrFRap = "#D62728")

annot_text_ce_rapa <- sprintf("Spearman ρ = %.2f\nR² = %.2f, p = %.3g",
                              rho_ce_rapa, r2_ce_rapa, pval_ce_rapa)

ggplot(dat_ce_rapa, aes(x = PC1_CE, y = Post_Ins_AUC, color = Group, fill = Group)) +
  stat_ellipse(type = "norm", level = 0.68, geom = "polygon", alpha = 0.18, color = NA) +
  geom_point(shape = 21, stroke = 2, size = 4, color = "white") +
  geom_point(shape = 21, size = 3, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  scale_color_manual(values = pal_rapa) +
  scale_fill_manual(values = pal_rapa) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cholesteryl Ester PC1 vs Post Insulin Sensitivity (AUC)",
    subtitle = "Old_PwrIRap + Old_PwrFRap combined",
    x = "CE PC1 (first principal component across CE species)",
    y = "Post Insulin Sensitivity (AUC)",
    color = "Group", fill = "Group"
  ) +
  annotate("text", x = Inf, y = -Inf, label = annot_text_ce_rapa,
           hjust = 1.05, vjust = -0.5, size = 4.5)

# 9️⃣ Top CE loadings (to see which species drive PC1)
pc1_loadings_ce_rapa <- as.data.frame(pca_ce_rapa$var$coord[, 1, drop = FALSE]) %>%
  tibble::rownames_to_column("Lipid") %>%
  rename(PC1_loading = Dim.1) %>%
  mutate(abs_loading = abs(PC1_loading)) %>%
  arrange(desc(abs_loading))
print(head(pc1_loadings_ce_rapa, 15))

































































#=============================================================================#
#=============================================================================#
#Ceramide correlation

yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sedveh_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwrveh_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

#-----Get Lipd Matrix------------#

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
#---------------------------------------------------#

triceps_physio <- read.csv("Triceps_Physio_Data_v2.csv", check.names = FALSE)
# Convert everything except "Sample" to numeric (some columns might be character)
# 2. Replace missing or blank names with synthetic ones

triceps_physio
#--------------------------------------------------------#

#=========================================================#

#--------------------------------------------
# 1️⃣ Prepare ceramide subset
#--------------------------------------------
cer_df <- lipid_df %>%
  dplyr::filter(!grepl("^YS", Sample)) %>%  # remove YS samples
  dplyr::select(Sample, starts_with("Cer"))  # keep only ceramides

# Average log2 ceramide signal per sample
cer_df <- cer_df %>%
  tidyr::pivot_longer(-Sample, names_to = "Lipid", values_to = "log2_value") %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarize(mean_ceramide = mean(log2_value, na.rm = TRUE))

#--------------------------------------------
# 2️⃣ Prepare physiologic data
#--------------------------------------------
physio <- read.csv("Triceps_Physio_Data_v2.csv", check.names = FALSE)

# Fix column name T_6 → T_06
names(physio) <- gsub("^T_6$", "T_06", names(physio))

# Extract Post GTT (AOC) row and transpose to long format
gtt_df <- physio %>%
  dplyr::filter(`Tube code` == "Post GTT (AOC)") %>%
  tidyr::pivot_longer(-`Tube code`, names_to = "Sample", values_to = "Post_GTT_AOC") %>%
  dplyr::select(Sample, Post_GTT_AOC)

#--------------------------------------------
# 3️⃣ Merge lipid and physiologic data
#--------------------------------------------
merged_df <- cer_df %>%
  dplyr::inner_join(gtt_df, by = "Sample")

#--------------------------------------------
# 4️⃣ Spearman correlation (aggregate ceramides)
#--------------------------------------------
cor_res <- cor.test(
  merged_df$mean_ceramide,
  merged_df$Post_GTT_AOC,
  method = "spearman",
  exact = FALSE
)

cat("Spearman rho:", cor_res$estimate, "\nP-value:", cor_res$p.value, "\n")

#--------------------------------------------
# 5️⃣ Visualization
#--------------------------------------------
library(ggplot2)

ggplot(merged_df, aes(x = mean_ceramide, y = Post_GTT_AOC)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    x = "Mean log2 Ceramide Abundance",
    y = "Post GTT (AOC)",
    title = "Ceramide Burden vs Glucose Tolerance"
  )
#--------------------------------------------------------------------#


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# --- 0) Basic hygiene on names ----------------------------------------------
lipid_df <- lipid_df %>%
  mutate(Sample = str_trim(Sample))

# How many ceramide columns do we really have?
cer_cols_idx <- grepl("^Cer", names(lipid_df))        # strict
cer_cols_idx2 <- grepl("^Cer[\\(:]", names(lipid_df)) # a bit more permissive

n_cer_strict <- sum(cer_cols_idx)
n_cer_perm   <- sum(cer_cols_idx2)

cat("Ceramide columns (strict ^Cer):", n_cer_strict, "\n")
cat("Ceramide columns (permissive ^Cer[(:]):", n_cer_perm, "\n")

# Choose a selector that actually finds columns
if (n_cer_strict > 0) {
  cer_selector <- cer_cols_idx
} else if (n_cer_perm > 0) {
  cer_selector <- cer_cols_idx2
  message("Using permissive ceramide selector: '^Cer[(:]'")
} else {
  stop("No ceramide columns found. Check lipid column names; e.g., are they 'Cer(d18:1/16:0)'?")
}

cer_cols <- c("Sample", names(lipid_df)[cer_selector])
cat("First few ceramide columns:", paste(head(cer_cols, 6), collapse = ", "), "\n")

# Exclude YS samples
lipid_cer_long <- lipid_df %>%
  filter(!grepl("^YS", Sample)) %>%
  select(all_of(cer_cols)) %>%
  pivot_longer(-Sample, names_to = "Lipid", values_to = "log2_value")

cat("Rows in ceramide-long after filtering YS:", nrow(lipid_cer_long), "\n")
cat("Unique samples in ceramide-long:", length(unique(lipid_cer_long$Sample)), "\n")

if (nrow(lipid_cer_long) == 0) {
  stop("Ceramide-long table has 0 rows. Either no ceramide columns, or all samples were filtered out.")
}

# Average ceramide per-sample (log2 already in your pipeline)
cer_df <- lipid_cer_long %>%
  group_by(Sample) %>%
  summarize(mean_ceramide = mean(log2_value, na.rm = TRUE), .groups = "drop")

# --- 1) Read physio and normalize sample names -------------------------------
physio <- read.csv("Triceps_Physio_Data_v2.csv", check.names = FALSE)

# Normalize column names: fix T_6 -> T_06, trim whitespace
physio_names <- names(physio)
physio_names <- gsub("^T_6$", "T_06", physio_names)
physio_names <- str_trim(physio_names)
names(physio) <- physio_names

# Extract the GTT row
gtt_df <- physio %>%
  filter(`Tube code` == "Post GTT (AOC)") %>%
  pivot_longer(-`Tube code`, names_to = "Sample", values_to = "Post_GTT_AOC") %>%
  select(Sample, Post_GTT_AOC) %>%
  mutate(Sample = str_trim(Sample),
         Post_GTT_AOC = suppressWarnings(as.numeric(Post_GTT_AOC)))

cat("GTT rows:", nrow(gtt_df), "\n")
cat("GTT unique samples:", length(unique(gtt_df$Sample)), "\n")
cat("GTT NAs:", sum(is.na(gtt_df$Post_GTT_AOC)), "\n")
#--------------------------------------------------------------#
#--------------------------------------------------------------#
# Create the sample_key lookup
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_01","OFR1","T_04","OFR2","T_22","OFR3","T_24","OFR4","T_28","OFR5","T_34","OFR6","T_40","OFR7","T_44","OFR8",
  "T_02","OIR1","T_03","OIR2","T_14","OIR3","T_31","OIR4","T_33","OIR5","T_37","OIR6","T_47","OIR7","T_48","OIR8",
  "T_08","OV1","T_10","OV2","T_11","OV3","T_16","OV4","T_41","OV5","T_43","OV6","T_45","OV7","T_46","OV8",
  "T_06","OS1","T_09","OS2","T_13","OS3","T_19","OS4","T_21","OS5","T_25","OS6","T_29","OS7","T_39","OS8",
  "T_07","YV1","T_12","YV2","T_18","YV3","T_20","YV4","T_27","YV5","T_30","YV6","T_36","YV7","T_42","YV8",
  "T_05","YS1","T_15","YS2","T_17","YS3","T_23","YS4","T_26","YS5","T_32","YS6","T_35","YS7","T_38","YS8"
)

# Join to replace metabolomics names (OFR1, etc.) with TubeCode (T_XX)
cer_df <- cer_df %>%
  dplyr::left_join(sample_key, by = c("Sample" = "MetabolomicsName")) %>%
  dplyr::mutate(Sample = TubeCode) %>%
  dplyr::select(Sample, mean_ceramide) %>%
  dplyr::filter(!is.na(Sample))

# Sanity check
cat("After mapping, cer_df sample names:\n")
print(sort(unique(cer_df$Sample)))
#--------------------------------------------#

merged_df <- cer_df %>%
  dplyr::inner_join(gtt_df, by = "Sample") %>%
  dplyr::filter(is.finite(mean_ceramide) & is.finite(Post_GTT_AOC))

cat("Merged rows (complete cases):", nrow(merged_df), "\n")
print(head(merged_df))

cor_res <- cor.test(
  merged_df$mean_ceramide,
  merged_df$Post_GTT_AOC,
  method = "spearman"
)
cor_res

ggplot(merged_df, aes(x = mean_ceramide, y = Post_GTT_AOC)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "steelblue", se = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    x = "Mean log2 Ceramide Abundance",
    y = "Post GTT (AOC)",
    title = "Ceramide Burden vs Glucose Tolerance"
  )
#====================================================#

#------------------------------------------------------------
# 🧠 1️⃣ Prepare ceramide-long data with T_XX IDs
#------------------------------------------------------------
lipid_cer_long_T <- lipid_cer_long %>%
  dplyr::left_join(sample_key, by = c("Sample" = "MetabolomicsName")) %>%
  dplyr::mutate(Sample = TubeCode) %>%
  dplyr::filter(!is.na(Sample)) %>%
  dplyr::select(Sample, Lipid, log2_value)

#------------------------------------------------------------
# 📈 2️⃣ Join with GTT and compute correlations per lipid
#------------------------------------------------------------
cer_corrs <- lipid_cer_long_T %>%
  dplyr::inner_join(gtt_df, by = "Sample") %>%
  dplyr::group_by(Lipid) %>%
  dplyr::summarize(
    n = sum(is.finite(log2_value) & is.finite(Post_GTT_AOC)),
    rho = ifelse(n >= 3, cor(log2_value, Post_GTT_AOC, method = "spearman"), NA_real_),
    pval = ifelse(n >= 3, suppressWarnings(cor.test(log2_value, Post_GTT_AOC, method = "spearman")$p.value), NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  dplyr::arrange(FDR)

# Show top results
print(cer_corrs, n = 10)

#------------------------------------------------------------
# 🎨 3️⃣ Plot correlations
#------------------------------------------------------------
ggplot(cer_corrs, aes(x = reorder(Lipid, rho), y = rho, fill = FDR < 0.05)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 13) +
  scale_fill_manual(values = c("grey70", "steelblue"), name = "FDR < 0.05") +
  labs(
    x = "Ceramide species",
    y = "Spearman correlation (ρ)",
    title = "Correlation of individual ceramides with Post GTT (AOC)"
  )

group_assign <- tibble(
  Sample = c(old_sedveh_samples, old_pwrveh_samples, old_pwrirap_samples, old_pwrfrap_samples),
  Group = c(rep("Old_SedVeh", length(old_sedveh_samples)),
            rep("Old_PwrVeh", length(old_pwrveh_samples)),
            rep("Old_PwrIRap", length(old_pwrirap_samples)),
            rep("Old_PwrFRap", length(old_pwrfrap_samples)))
)

# Example: visualize top correlated ceramide
top_lipid <- cer_corrs$Lipid[1]

lipid_cer_long_T %>%
  dplyr::filter(Lipid == top_lipid) %>%
  dplyr::inner_join(gtt_df, by = "Sample") %>%
  dplyr::inner_join(group_assign, by = "Sample") %>%
  ggplot(aes(x = log2_value, y = Post_GTT_AOC, color = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_classic(base_size = 13) +
  labs(
    x = paste0(top_lipid, " (log2)"),
    y = "Post GTT (AOC)",
    title = paste("Relationship between", top_lipid, "and glucose tolerance")
  )

#===========================================================#
#===========================================================#

#-------------------------------------------------------------------
# 🧩 1️⃣ Define rapamycin-only samples
#-------------------------------------------------------------------
rapa_samples <- c(old_pwrirap_samples, old_pwrfrap_samples)

#-------------------------------------------------------------------
# 🧠 2️⃣ Prepare ceramide-long data with TubeCode mapping
#-------------------------------------------------------------------
lipid_cer_long_T <- lipid_cer_long %>%
  dplyr::left_join(sample_key, by = c("Sample" = "MetabolomicsName")) %>%
  dplyr::mutate(Sample = TubeCode) %>%
  dplyr::filter(!is.na(Sample) & Sample %in% rapa_samples) %>%
  dplyr::select(Sample, Lipid, log2_value)

cat("Unique rapamycin samples in lipid data:", length(unique(lipid_cer_long_T$Sample)), "\n")

#-------------------------------------------------------------------
# 📊 3️⃣ Extract GTT and Insulin Sensitivity data
#-------------------------------------------------------------------
physio <- read.csv("Triceps_Physio_Data_v2.csv", check.names = FALSE)
names(physio) <- gsub("^T_6$", "T_06", names(physio))
names(physio) <- stringr::str_trim(names(physio))

gtt_df <- physio %>%
  dplyr::filter(`Tube code` == "Post GTT (AOC)") %>%
  tidyr::pivot_longer(-`Tube code`, names_to = "Sample", values_to = "Post_GTT_AOC") %>%
  dplyr::mutate(Sample = stringr::str_trim(Sample),
                Post_GTT_AOC = as.numeric(Post_GTT_AOC))

ins_df <- physio %>%
  dplyr::filter(`Tube code` == "Post Insulin Sensitivity (AUC)") %>%
  tidyr::pivot_longer(-`Tube code`, names_to = "Sample", values_to = "Post_Ins_AUC") %>%
  dplyr::mutate(Sample = stringr::str_trim(Sample),
                Post_Ins_AUC = as.numeric(Post_Ins_AUC))

# Merge both into one table for convenience
physio_2traits <- dplyr::inner_join(gtt_df, ins_df, by = "Sample")

cat("Unique samples in physio data:", length(unique(physio_2traits$Sample)), "\n")

#-------------------------------------------------------------------
# 🧮 4️⃣ Compute per-ceramide Spearman correlations for each trait
#-------------------------------------------------------------------
calc_corrs <- function(trait) {
  trait_col <- sym(trait)
  lipid_cer_long_T %>%
    dplyr::inner_join(physio_2traits, by = "Sample") %>%
    dplyr::group_by(Lipid) %>%
    dplyr::summarize(
      n = sum(is.finite(log2_value) & is.finite(!!trait_col)),
      rho = ifelse(n >= 3, cor(log2_value, !!trait_col, method = "spearman"), NA_real_),
      pval = ifelse(n >= 3, suppressWarnings(cor.test(log2_value, !!trait_col, method = "spearman")$p.value), NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Trait = trait,
      FDR = p.adjust(pval, method = "fdr")
    )
}

cer_corrs_gtt  <- calc_corrs("Post_GTT_AOC")
cer_corrs_ins  <- calc_corrs("Post_Ins_AUC")

cer_corrs_all <- dplyr::bind_rows(cer_corrs_gtt, cer_corrs_ins)

#-------------------------------------------------------------------
# 📈 5️⃣ View top results
#-------------------------------------------------------------------
cer_corrs_all %>%
  dplyr::arrange(Trait, FDR) %>%
  dplyr::group_by(Trait) %>%
  dplyr::slice_min(order_by = FDR, n = 10) %>%
  print(n = 20)

print(cer_corrs_all, n = 40)

#-------------------------------------------------------------------
# 🎨 6️⃣ Plot correlations for both traits
#-------------------------------------------------------------------
ggplot(
  cer_corrs_all %>% dplyr::filter(!is.na(rho)),
  aes(x = reorder(Lipid, rho), y = rho, fill = FDR < 0.05)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Trait, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 13) +
  scale_fill_manual(values = c("grey70", "steelblue"), name = "FDR < 0.05") +
  labs(
    x = "Ceramide species",
    y = "Spearman correlation (ρ)",
    title = "Correlations of ceramides with Post GTT and Insulin Sensitivity (AUC)\n(Rapamycin-treated mice only)"
  )
#------------------------------------------------------#
install.packages(c("FactoMineR", "factoextra"))
library(FactoMineR); library(factoextra)

# Subset to rapamycin samples
cer_mat <- lipid_cer_long_T %>%
  filter(Sample %in% rapa_samples) %>%
  pivot_wider(names_from = Lipid, values_from = log2_value) %>%
  column_to_rownames("Sample")

# Option 1: total ceramide mean
cer_mat$total_mean <- rowMeans(cer_mat, na.rm = TRUE)

# Option 2: PCA
pca_res <- PCA(cer_mat[ , grepl("^Cer", colnames(cer_mat))], graph = FALSE)
pc1_scores <- pca_res$ind$coord[,1]

# Correlate PC1 or total_mean with both traits
merged_physio <- inner_join(physio_2traits, tibble(Sample = names(pc1_scores),
                                                   PC1 = pc1_scores,
                                                   meanCer = cer_mat$total_mean),
                            by = "Sample")

cor.test(merged_physio$PC1, merged_physio$Post_Ins_AUC, method="spearman")
ggplot(merged_physio, aes(x = PC1, y = Post_Ins_AUC)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  theme_classic(base_size = 14) +
  labs(
    x = "Ceramide PC1 (overall ceramide burden)",
    y = "Post Insulin Sensitivity (AUC)",
    title = "Ceramide PC1 vs Insulin Sensitivity (rapamycin-treated mice)"
  )

#------------------------------------------------------------
# 📊 1️⃣ Extract PC1 loadings (variable contributions)
#------------------------------------------------------------
pc1_loadings <- as.data.frame(pca_res$var$coord[, 1, drop = FALSE]) %>%
  tibble::rownames_to_column("Lipid") %>%
  dplyr::rename(PC1_loading = Dim.1) %>%
  dplyr::mutate(abs_loading = abs(PC1_loading)) %>%
  dplyr::arrange(desc(abs_loading))

# View top contributors numerically
print(head(pc1_loadings, 30))

#------------------------------------------------------------
# 🎨 2️⃣ Plot top loadings (by absolute value)
#------------------------------------------------------------
library(ggplot2)

top_n <- 15  # choose how many to display
ggplot(pc1_loadings[1:top_n, ], aes(x = reorder(Lipid, PC1_loading), y = PC1_loading)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = "Ceramide species",
    y = "PC1 loading value",
    title = paste0("Top ", top_n, " ceramides driving PC1 (overall ceramide burden)")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")

#============================================================================#
#=============================================================================#

#------------------------------------------------------------
# 🧩 1️⃣ Rebuild cer_mat_all without the TubeCode column
#------------------------------------------------------------
cer_mat_all <- lipid_cer_long %>%
  dplyr::left_join(sample_key, by = c("Sample" = "MetabolomicsName")) %>%
  dplyr::mutate(Sample = TubeCode) %>%
  dplyr::filter(!is.na(Sample) & !grepl("^YS", Sample)) %>%
  tidyr::pivot_wider(names_from = Lipid, values_from = log2_value) %>%
  dplyr::select(-TubeCode) %>%                 # ❗ drop the non-numeric column
  tibble::column_to_rownames("Sample")

# confirm it's now all numeric
str(cer_mat_all[ , 1:3])

#------------------------------------------------------------
# 📊 2️⃣ PCA
#------------------------------------------------------------
pca_all <- FactoMineR::PCA(cer_mat_all[, grepl("^Cer", colnames(cer_mat_all))], graph = FALSE)
pc1_scores_all <- pca_all$ind$coord[, 1]

#------------------------------------------------------------
# 🧪 3️⃣ Merge with physiology
#------------------------------------------------------------
merged_physio_all <- physio_all %>%
  dplyr::inner_join(
    tibble(Sample = names(pc1_scores_all),
           PC1 = pc1_scores_all,
           meanCer = rowMeans(cer_mat_all, na.rm = TRUE)),
    by = "Sample"
  )

# sanity check
print(dim(merged_physio_all))
head(merged_physio_all)

#------------------------------------------------------------
# 📈 4️⃣ Correlations
#------------------------------------------------------------
cor_gtt <- cor.test(merged_physio_all$PC1, merged_physio_all$Post_GTT_AOC, method = "spearman")
cor_ins <- cor.test(merged_physio_all$PC1, merged_physio_all$Post_Ins_AUC, method = "spearman")

cor_gtt
cor_ins

#------------------------------------------------------------
# 🎨 5️⃣  Visualization
#------------------------------------------------------------
library(ggplot2)

ggplot(merged_physio_all, aes(x = PC1, y = Post_Ins_AUC)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  theme_classic(base_size = 14) +
  labs(
    x = "Ceramide PC1 (all samples)",
    y = "Post Insulin Sensitivity (AUC)",
    title = "Ceramide PC1 vs Insulin Sensitivity (all groups)"
  )

ggplot(merged_physio_all, aes(x = PC1, y = Post_GTT_AOC)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  theme_classic(base_size = 14) +
  labs(
    x = "Ceramide PC1 (all samples)",
    y = "Post GTT (AOC)",
    title = "Ceramide PC1 vs Glucose Tolerance (all groups)"
  )

#------------------------------------------------------------
# 🧬 6️⃣  Examine loadings
#------------------------------------------------------------
pc1_loadings_all <- as.data.frame(pca_all$var$coord[, 1, drop = FALSE]) %>%
  tibble::rownames_to_column("Lipid") %>%
  dplyr::rename(PC1_loading = Dim.1) %>%
  dplyr::mutate(abs_loading = abs(PC1_loading)) %>%
  dplyr::arrange(desc(abs_loading))

print(head(pc1_loadings_all, 15))
fviz_contrib(pca_all, choice = "var", axes = 1, top = 20) +
  labs(title = "Top ceramides contributing to PC1 (all samples)")


# Spearman correlation between Post_GTT_AOC and Post_Ins_AUC
cor_physio <- cor.test(
  merged_physio_all$Post_GTT_AOC,
  merged_physio_all$Post_Ins_AUC,
  method = "spearman"
)

cor_physio

# Optional: visualize
library(ggplot2)
ggplot(merged_physio_all, aes(x = Post_Ins_AUC, y = Post_GTT_AOC)) +
  geom_point(size = 3, color = "darkorange") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_classic(base_size = 14) +
  labs(
    x = "Post Insulin Sensitivity (AUC)",
    y = "Post Glucose Tolerance (AOC)",
    title = "Relationship between Insulin Sensitivity and Glucose Tolerance"
  )

#=================================================================#

# Extract both rows of interest from your physio table
ins_sens_df <- triceps_physio %>%
  dplyr::filter(`Tube code` %in% c("Post Insulin Sensitivity (AUC)",
                                   "Delta Insulin Sensitivity (AUC)")) %>%
  tidyr::pivot_longer(-`Tube code`, names_to = "Sample", values_to = "Value") %>%
  tidyr::pivot_wider(names_from = `Tube code`, values_from = Value) %>%
  dplyr::mutate(across(-Sample, as.numeric))

# Check
head(ins_sens_df)

# Compute Spearman correlation between Post and Delta Insulin Sensitivity
cor_ins_check <- cor.test(
  ins_sens_df$`Post Insulin Sensitivity (AUC)`,
  ins_sens_df$`Delta Insulin Sensitivity (AUC)`,
  method = "spearman",
  use = "complete.obs"
)

cor_ins_check

# Optional: visualize
library(ggplot2)
ggplot(ins_sens_df, aes(x = `Post Insulin Sensitivity (AUC)`,
                        y = `Delta Insulin Sensitivity (AUC)`)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_classic(base_size = 14) +
  labs(
    x = "Post Insulin Sensitivity (AUC)",
    y = "Delta Insulin Sensitivity (AUC)",
    title = "Sanity Check: Post vs Delta Insulin Sensitivity"
  )

#=================================================================#

# Extract both rows: Post and Delta GTT (AOC)
gtt_df_check <- triceps_physio %>%
  dplyr::filter(`Tube code` %in% c("Post GTT (AOC)",
                                   "Delta GTT (AOC)")) %>%
  tidyr::pivot_longer(-`Tube code`, names_to = "Sample", values_to = "Value") %>%
  tidyr::pivot_wider(names_from = `Tube code`, values_from = Value) %>%
  dplyr::mutate(across(-Sample, as.numeric))

# Inspect
head(gtt_df_check)

# Spearman correlation between Post and Delta GTT AOC
cor_gtt_check <- cor.test(
  gtt_df_check$`Post GTT (AOC)`,
  gtt_df_check$`Delta GTT (AOC)`,
  method = "spearman",
  use = "complete.obs"
)

cor_gtt_check

# Optional visualization
library(ggplot2)
ggplot(gtt_df_check, aes(x = `Post GTT (AOC)`,
                         y = `Delta GTT (AOC)`)) +
  geom_point(size = 3, color = "firebrick") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_classic(base_size = 14) +
  labs(
    x = "Post Glucose Tolerance (AOC)",
    y = "Delta Glucose Tolerance (AOC)",
    title = "Sanity Check: Post vs Delta GTT (AOC)"
  )

