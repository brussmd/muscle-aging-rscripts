#Sedentary Rapa RNA to ITT/GTT correlations.
setwd("/Users/mdbruss/Documents/RStudioProjects_2/Rapa_PwR")
getwd()


# ---------------------------------------------------------
# 📦 Load required packages
# ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(stringr)
  library(edgeR)
  library(org.Mm.eg.db)
  library(ggplot2)
})

# ---------------------------------------------------------
# 🧱 Helper function: read and clean counts
# ---------------------------------------------------------
read_clean_counts <- function(file, keep_samples) {
  df <- read_xlsx(file)
  cn <- colnames(df)
  cn[1] <- "Ensembl"
  colnames(df) <- cn
  
  df <- df %>%
    dplyr::select(Ensembl, dplyr::all_of(keep_samples)) %>%
    mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
           ENTREZID = mapIds(org.Mm.eg.db,
                             keys = Ensembl_noDec,
                             keytype = "ENSEMBL",
                             column = "ENTREZID",
                             multiVals = "first")) %>%
    tidyr::drop_na(ENTREZID)
  
  counts_mat <- df %>%
    dplyr::select(dplyr::all_of(keep_samples)) %>%
    replace(is.na(.), 0) %>%
    as.matrix()
  rownames(counts_mat) <- df$ENTREZID
  
  return(counts_mat)
}

# ---------------------------------------------------------
# 🧩 Load previous logCPM matrix (YNG + OLD) and features
# ---------------------------------------------------------
logCPM_symbol_yngOLDsed <- readRDS("logCPM_symbol_yngOLDsed.RDS")

# Your existing selected_df is already loaded in the environment
# Filter only RNA features from it
selected_rna <- selected_df %>%
  dplyr::filter(block == "RNA")

selected_rna
head(logCPM_symbol_yngOLDsed)

# -----------------------------------------------------------
# 🧩 1) Inputs
# -----------------------------------------------------------
# selected_rna already in environment
# logCPM_symbol_yngOLDsed already loaded

yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")

# -----------------------------------------------------------
# 🧮 2) Extract RNA features and loadings
# -----------------------------------------------------------
axis_genes  <- selected_rna$feature
axis_loads  <- selected_rna$loading
names(axis_loads) <- axis_genes

# Case-insensitive match between feature names and matrix rownames
common_genes <- intersect(toupper(axis_genes), toupper(rownames(logCPM_symbol_yngOLDsed)))
cat("Matched genes:", length(common_genes), "of", length(axis_genes), "\n")

# -----------------------------------------------------------
# 🧮 Compute sample scores properly
# -----------------------------------------------------------

# Make sure expression rows = genes, columns = samples
expr_mat <- logCPM_symbol_yngOLDsed[match(axis_genes, rownames(logCPM_symbol_yngOLDsed)), , drop = FALSE]

# Align loadings
axis_loads_vec <- axis_loads[match(rownames(expr_mat), axis_genes)]

# Compute weighted sum across genes for each sample
# (each column = a sample)
aging_axis_scores <- apply(expr_mat, 2, function(sample_expr) {
  sum(sample_expr * axis_loads_vec, na.rm = TRUE)
})

# Build dataframe
score_df <- data.frame(
  Sample = names(aging_axis_scores),
  Score  = as.numeric(aging_axis_scores)
) %>%
  mutate(Group = case_when(
    Sample %in% yng_sed_samples ~ "YNG_SED",
    Sample %in% old_sed_samples ~ "OLD_SED",
    TRUE ~ NA_character_
  ))

# Check assignment
table(score_df$Group, useNA = "ifany")

# -----------------------------------------------------------
# Plot
# -----------------------------------------------------------
library(ggplot2)
ggplot(score_df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  theme_classic(base_size = 14) +
  labs(
    title = "RNA Aging Axis (50-gene projection)",
    y = "Weighted Expression (Σ gene × loading)",
    x = ""
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# -----------------------------------------------------------
# Summary
# -----------------------------------------------------------
score_df %>%
  group_by(Group) %>%
  summarise(mean_score = mean(Score), sd_score = sd(Score))

t.test(Score ~ Group, data = score_df)

logCPM_symbol_new

# -----------------------------------------------------------
# 🧩 1) Normalize and prepare combined RNA matrix (already done)
# -----------------------------------------------------------
# You already have logCPM_symbol_new from IRAP + FRAP normalization
# And logCPM_symbol_yngOLDsed from the YNG + OLD normalization
# Merge the two on shared genes:
common_genes <- intersect(rownames(logCPM_symbol_yngOLDsed), rownames(logCPM_symbol_new))
combined_matrix <- cbind(
  logCPM_symbol_yngOLDsed[common_genes, ],
  logCPM_symbol_new[common_genes, ]
)

# -----------------------------------------------------------
# 🧬 2) Prepare RNA axis features and loadings
# -----------------------------------------------------------
axis_genes  <- selected_rna$feature
axis_loads  <- selected_rna$loading
names(axis_loads) <- axis_genes

# Match features (case-insensitive)
common_feats <- intersect(toupper(axis_genes), toupper(rownames(combined_matrix)))

expr_mat <- combined_matrix[match(common_feats, toupper(rownames(combined_matrix))), , drop = FALSE]
rownames(expr_mat) <- common_feats
axis_loads_vec <- axis_loads[match(common_feats, toupper(names(axis_loads)))]

# -----------------------------------------------------------
# 🧮 3) Compute per-sample aging-axis scores
# -----------------------------------------------------------
aging_axis_scores <- apply(expr_mat, 2, function(sample_expr) {
  sum(sample_expr * axis_loads_vec, na.rm = TRUE)
})

# -----------------------------------------------------------
# 🧩 4) Build dataframe with all samples and groups
# -----------------------------------------------------------
score_df <- data.frame(
  Sample = names(aging_axis_scores),
  Score  = as.numeric(aging_axis_scores)
) %>%
  mutate(Group = case_when(
    Sample %in% yng_sed_samples     ~ "YNG_SED",
    Sample %in% old_sed_samples     ~ "OLD_SED",
    Sample %in% old_sedirap_samples ~ "OLD_SED_IRAP",
    Sample %in% old_sedfrap_samples ~ "OLD_SED_FRAP",
    TRUE ~ NA_character_
  ))

# Flip sign if you want “higher = older” consistency
score_df <- score_df %>%
  mutate(Score = -Score)

# -----------------------------------------------------------
# 📊 5) Visualize projection
# -----------------------------------------------------------
ggplot(score_df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic(base_size = 14) +
  labs(
    title = "Projection of IRAP and FRAP Samples onto RNA Aging Axis",
    y = "RNA Aging Axis (higher = older)",
    x = ""
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# -----------------------------------------------------------
# 📈 6) Group means and stats
# -----------------------------------------------------------
score_df %>%
  group_by(Group) %>%
  summarise(mean_score = mean(Score), sd_score = sd(Score))

anova_res <- aov(Score ~ Group, data = score_df)
summary(anova_res)

TukeyHSD(anova_res)


axis_genes  <- selected_rna$feature
axis_loads  <- selected_rna$loading
names(axis_loads) <- axis_genes

# ---------------------------------------------------------
# 🧫 Define sample groups
# ---------------------------------------------------------
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")

# ---------------------------------------------------------
# 🆕 Read in IRAP and FRAP count data
# ---------------------------------------------------------

# Read both Excel files using your existing helper function
counts_irap <- read_clean_counts(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx",
  keep_samples = c(old_sedirap_samples, yng_sed_samples)
)

counts_frap <- read_clean_counts(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx",
  keep_samples = c(old_sedfrap_samples, yng_sed_samples)
)

# ---------------------------------------------------------
# 🧹 Keep only OLD_SED_IRAP and OLD_SED_FRAP samples
# ---------------------------------------------------------
counts_irap <- counts_irap[, old_sedirap_samples, drop = FALSE]
counts_frap <- counts_frap[, old_sedfrap_samples, drop = FALSE]

# ---------------------------------------------------------
# 🔗 Merge the two count matrices (shared genes only)
# ---------------------------------------------------------
common_genes <- intersect(rownames(counts_irap), rownames(counts_frap))
counts_combined <- cbind(
  counts_irap[common_genes, ],
  counts_frap[common_genes, ]
)

# ---------------------------------------------------------
# ⚙️ Normalize (TMM → logCPM)
# ---------------------------------------------------------
dge <- DGEList(counts = counts_combined)
dge <- calcNormFactors(dge, method = "TMM")
logCPM_new <- cpm(dge, log = TRUE, prior.count = 1)

# ---------------------------------------------------------
# 🧬 Map ENTREZ → SYMBOL and collapse duplicates
# ---------------------------------------------------------
symbols_new <- mapIds(org.Mm.eg.db,
                      keys    = rownames(logCPM_new),
                      keytype = "ENTREZID",
                      column  = "SYMBOL",
                      multiVals = "first")

logCPM_symbol_new <- as.data.frame(logCPM_new) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  mutate(Symbol = symbols_new) %>%
  tidyr::drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()

# ---------------------------------------------------------
# 🧠 Quick diagnostic: overlap with RNA aging axis
# ---------------------------------------------------------
common_axis_genes <- intersect(axis_genes, rownames(logCPM_symbol_new))
cat("\n----------------------------------------------\n")
cat("RNA Aging Axis genes found in IRAP+FRAP dataset:",
    length(common_axis_genes), "of", length(axis_genes), "\n")
cat("----------------------------------------------\n\n")

# ---------------------------------------------------------
# 🧩 Merge old + new normalized RNA data
# ---------------------------------------------------------
common_genes <- intersect(rownames(logCPM_symbol_yngOLDsed), rownames(logCPM_symbol_new))

combined_matrix <- cbind(
  logCPM_symbol_yngOLDsed[common_genes, ],
  logCPM_symbol_new[common_genes, ]
)

combined_matrix
colnames(combined_matrix)

selected_df
# -----------------------------------------------------------
# 🧩 1) Extract RNA features and loadings from selected_df
# -----------------------------------------------------------
selected_rna <- selected_df %>%
  dplyr::filter(block == "RNA")

axis_genes  <- selected_rna$feature
axis_loads  <- selected_rna$loading
names(axis_loads) <- axis_genes

axis_genes
expr_mat

# -----------------------------------------------------------
# 🧫 2) Build RNA-only input block for prediction
# -----------------------------------------------------------
# Keep only the 50 genes used in the DIABLO RNA component
common_feats <- intersect(axis_genes, rownames(exp_mat))

X_new <- list(RNA = t(combined_matrix[common_feats, , drop = FALSE]))

# -----------------------------------------------------------
# 🧮 3) Project onto the original DIABLO model
# -----------------------------------------------------------
proj_new <- predict(diablo_res, newdata = X_new)

# Extract RNA block Comp1 scores
scores_new <- proj_new$variates$RNA[, 1, drop = FALSE]
scores_new <- as.data.frame(scores_new)
colnames(scores_new) <- "Comp1"
scores_new$Sample <- rownames(scores_new)

# -----------------------------------------------------------
# 🏷️ 4) Assign experimental group labels
# -----------------------------------------------------------
scores_new <- scores_new %>%
  mutate(Group = case_when(
    Sample %in% yng_sed_samples     ~ "YNG_SED",
    Sample %in% old_sed_samples     ~ "OLD_SED",
    Sample %in% old_sedirap_samples ~ "OLD_SED_IRAP",
    Sample %in% old_sedfrap_samples ~ "OLD_SED_FRAP",
    TRUE ~ NA_character_
  ))

# -----------------------------------------------------------
# 🎨 5) Plot projection onto RNA aging axis
# -----------------------------------------------------------
ggplot(scores_new, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic(base_size = 13) +
  labs(
    y = "DIABLO Comp1 (RNA Aging Axis)",
    x = "",
    title = "Projection of New Samples onto Original DIABLO Axis"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# -----------------------------------------------------------
# 🧪 6) Sanity check: confirm directionality
# -----------------------------------------------------------
scores_new %>%
  filter(Group %in% c("YNG_SED", "OLD_SED")) %>%
  group_by(Group) %>%
  summarise(mean_Comp1 = mean(Comp1))

















#=================================================================#
#===========OLD DIDNT REALLY WORK=================================#
#=================================================================#
# ---------------------------------------------------------
# 🧩 1. Extract RNA Aging Axis features and loadings
# ---------------------------------------------------------
selected_rna <- selected_df %>%
  dplyr::filter(block == "RNA")
axis_genes  <- selected_rna$feature
axis_loads  <- selected_rna$loading
names(axis_loads) <- axis_genes
selected_rna

# ---------------------------------------------------------
# 🧬 2. Subset combined_matrix to RNA axis genes
# ---------------------------------------------------------
expr_mat <- combined_matrix[axis_genes, , drop = FALSE]
expr_mat


# Ensure order consistency
expr_mat <- expr_mat[match(names(axis_loads), rownames(expr_mat)), , drop = FALSE]
expr_mat

# ---------------------------------------------------------
# 🧮 3. Compute RNA aging-axis scores (dot product)
# ---------------------------------------------------------
aging_axis_scores <- as.numeric(crossprod(axis_loads[rownames(expr_mat)], expr_mat))
names(aging_axis_scores) <- colnames(expr_mat)
head(aging_axis_scores)
aging_axis_scores

# ---------------------------------------------------------
# 📊 Build data frame for plotting
# ---------------------------------------------------------
score_df <- data.frame(
  Sample = names(aging_axis_scores),
  Score  = as.numeric(aging_axis_scores),
  stringsAsFactors = FALSE
)

score_df <- score_df %>%
  mutate(Group = case_when(
    Sample %in% yng_sed_samples     ~ "YNG_SED",
    Sample %in% old_sed_samples     ~ "OLD_SED",
    Sample %in% old_sedirap_samples ~ "OLD_SED_IRAP",
    Sample %in% old_sedfrap_samples ~ "OLD_SED_FRAP",
    TRUE ~ NA_character_
  ))

print(table(score_df$Group, useNA = "ifany"))

# ---------------------------------------------------------
# ⚙️ Normalize to YNG/OLD z-score scale
# ---------------------------------------------------------
ref_scores <- score_df %>% filter(Group %in% c("YNG_SED", "OLD_SED"))
mean_ref <- mean(ref_scores$Score)
sd_ref   <- sd(ref_scores$Score)

score_df <- score_df %>%
  mutate(Score_z = (Score - mean_ref) / sd_ref)

# ---------------------------------------------------------
# 🎨 Plot RNA Aging Axis Projection
# ---------------------------------------------------------
library(ggplot2)

ggplot(score_df, aes(x = Group, y = Score_z, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic(base_size = 13) +
  labs(
    y = "RNA Aging Axis (z-score)",
    x = "",
    title = "Projection of Samples onto RNA Aging Axis"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# (Optional) View group means for reporting
score_df %>%
  group_by(Group) %>%
  summarise(mean_z = mean(Score_z), sd_z = sd(Score_z))
#=============================================================#
#=============================================================#
#^Need some of these variables but starting over with new analysis^#
#==============================================================#



suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
})

# Groups
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")

# Quick checks that columns exist
stopifnot(all(yng_sed_samples %in% colnames(logCPM_symbol_yngOLDsed)))
stopifnot(all(old_sed_samples %in% colnames(logCPM_symbol_yngOLDsed)))
stopifnot(all(c(old_sedirap_samples, old_sedfrap_samples) %in% colnames(logCPM_symbol_new)))
#-----------------------------#

expr_train <- logCPM_symbol_yngOLDsed[, c(yng_sed_samples, old_sed_samples), drop = FALSE]
group_train <- factor(c(rep("YNG", length(yng_sed_samples)),
                        rep("OLD", length(old_sed_samples))),
                      levels = c("YNG","OLD"))

cat("Training matrix dim (genes x samples):", paste(dim(expr_train), collapse=" x "), "\n")
cat("Training groups:", table(group_train), "\n")
#---------------------#

set.seed(123)
splsda_rna <- splsda(
  X = t(expr_train),   # mixOmics expects samples as rows
  Y = group_train,
  ncomp = 1,
  keepX = 50
)

cat("Non-zero loadings on Comp1:", sum(splsda_rna$loadings$X[,1] != 0), "\n")

# The 50 selected genes (by name)
sel_tbl <- mixOmics::selectVar(splsda_rna, comp = 1)
sel <- sel_tbl$name
cat("Length of selected gene set:", length(sel), "\n")
cat("First 10 selected genes:\n"); print(head(sel, 10))
#--------------------------------------------------------#

# Build training matrix with only the 50 genes, samples as rows
# (use rows=genes from expr_train; keep drop=FALSE to preserve dims)
expr_train_50 <- t(expr_train[sel, , drop = FALSE])

# Sanity checks
stopifnot(identical(colnames(expr_train_50), sel))           # columns are the selected genes, ordered
stopifnot(identical(rownames(expr_train_50), c(yng_sed_samples, old_sed_samples))) # rows are training samples

# Refit on the 50 genes (keepX = 50 keeps them all)
set.seed(123)
splsda_50 <- splsda(
  X = expr_train_50,
  Y = group_train,
  ncomp = 1,
  keepX = 50
)

# Training scores dataframe
train_df <- data.frame(
  Sample = rownames(splsda_50$variates$X),
  Comp1  = splsda_50$variates$X[, 1],
  Group  = group_train
)

cat("Compact model trained on:", ncol(expr_train_50), "genes and", nrow(expr_train_50), "samples.\n")
#---------------------------------------#

# Column means from training for each of the 50 genes (vector length 50, names=gene)
train_means <- colMeans(expr_train_50)
stopifnot(identical(names(train_means), sel))

# New sample order we want to project
new_samp_order <- c(old_sedirap_samples, old_sedfrap_samples)

# Start an empty matrix (rows = new samples, cols = 50 genes), filled with training means
expr_new_50 <- matrix(rep(train_means, each = length(new_samp_order)),
                      nrow = length(new_samp_order), byrow = TRUE,
                      dimnames = list(new_samp_order, sel))

# Which selected genes are present in the new matrix?
present_genes_new <- intersect(sel, rownames(logCPM_symbol_new))
missing_genes_new <- setdiff(sel, present_genes_new)

cat("New data: present selected genes:", length(present_genes_new), " / 50\n")
if (length(missing_genes_new) > 0) {
  cat("Missing genes in new data (filled with training mean):\n")
  print(missing_genes_new)
}

# Overwrite the present genes with real values from the new matrix
# logCPM_symbol_new is genes x samples, so subset rows=genes, cols=new samples; transpose to samples x genes
if (length(present_genes_new) > 0) {
  expr_new_50[, present_genes_new] <- t(logCPM_symbol_new[present_genes_new, new_samp_order, drop = FALSE])
}

# Final sanity checks
stopifnot(identical(colnames(expr_new_50), sel))
stopifnot(identical(rownames(expr_new_50), new_samp_order))
#------------------------------------------------------------#

proj_new_50 <- predict(splsda_50, newdata = expr_new_50)
str(proj_new_50, max.level = 2)


# Extract Comp1 scores (you already confirmed structure)
scores_new <- proj_new_50$variates[, 1]
samples_new <- rownames(proj_new_50$variates)

# Build the projection dataframe
proj_df <- data.frame(
  Sample = samples_new,
  Comp1  = as.numeric(scores_new),
  Group  = ifelse(samples_new %in% old_sedirap_samples, "OLD_SED_IRAP", "OLD_SED_FRAP")
)

# Quick check
head(proj_df)

# Combine training + projected samples
plot_df <- bind_rows(
  train_df %>% mutate(Source = "Training"),
  proj_df %>% mutate(Source = "Projected")
)

# Plot
ggplot(plot_df, aes(x = Group, y = Comp1, color = Source, shape = Source)) +
  geom_jitter(width = 0.1, size = 3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1, width = 0.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "Projection of OLD_SED_IRAP / OLD_SED_FRAP onto YNG–OLD RNA Axis",
    y = "sPLS-DA Comp1 (RNA Aging Axis)",
    x = "Group"
  ) +
  scale_color_manual(values = c("Training" = "gray40", "Projected" = "steelblue"))

selected_genes <- mixOmics::selectVar(splsda_50, comp = 1)$name

selected_genes

oldsedveh_v_yngsedveh_genes %>%
  filter(Symbol %in% selected_genes) %>%
  arrange(logFC)
#-----------------------------------------------------#
#-----------------------------------------------------#
#=====================================================#

#---------------------------------------------------------------------------#

#=====================================================#
#======200 Features===================================#
#=====================================================#
suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
})
head(logCPM_symbol_all)
head(logCPM_symbol_yngOLDsed)
head(logCPM_symbol_new)
old_pwrveh_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

# -----------------------------------------------------------
# 🧫 Define groups
# -----------------------------------------------------------
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")

# Quick checks that columns exist
stopifnot(all(yng_sed_samples %in% colnames(logCPM_symbol_yngOLDsed)))
stopifnot(all(old_sed_samples %in% colnames(logCPM_symbol_yngOLDsed)))
stopifnot(all(c(old_sedirap_samples, old_sedfrap_samples) %in% colnames(logCPM_symbol_new)))

# -----------------------------------------------------------
# 🧱 Prepare training data (YNG vs OLD)
# -----------------------------------------------------------
expr_train <- logCPM_symbol_yngOLDsed[, c(yng_sed_samples, old_sed_samples), drop = FALSE]
group_train <- factor(c(rep("YNG", length(yng_sed_samples)),
                        rep("OLD", length(old_sed_samples))),
                      levels = c("YNG","OLD"))

cat("Training matrix dim (genes x samples):", paste(dim(expr_train), collapse=" x "), "\n")
cat("Training groups:", table(group_train), "\n")

# -----------------------------------------------------------
# 🧮 Train sparse PLS-DA (sPLS-DA) model with 200 features
# -----------------------------------------------------------
set.seed(123)
splsda_rna <- splsda(
  X = t(expr_train),   # mixOmics expects samples as rows
  Y = group_train,
  ncomp = 1,
  keepX = 200
)

cat("Non-zero loadings on Comp1:", sum(splsda_rna$loadings$X[,1] != 0), "\n")

# The 200 selected genes (by name)
sel_tbl <- mixOmics::selectVar(splsda_rna, comp = 1)
sel <- sel_tbl$name
cat("Length of selected gene set:", length(sel), "\n")
cat("First 10 selected genes:\n"); print(head(sel, 10))

# -----------------------------------------------------------
# 🧬 Build training matrix with only selected genes
# -----------------------------------------------------------
expr_train_200 <- t(expr_train[sel, , drop = FALSE])

# Sanity checks
stopifnot(identical(colnames(expr_train_200), sel))
stopifnot(identical(rownames(expr_train_200), c(yng_sed_samples, old_sed_samples)))

# Refit on the 200 genes (keepX = 200 keeps them all)
set.seed(123)
splsda_200 <- splsda(
  X = expr_train_200,
  Y = group_train,
  ncomp = 1,
  keepX = 200
)

# Training scores dataframe
train_df <- data.frame(
  Sample = rownames(splsda_200$variates$X),
  Comp1  = splsda_200$variates$X[, 1],
  Group  = group_train
)

cat("Compact model trained on:", ncol(expr_train_200), "genes and", nrow(expr_train_200), "samples.\n")

# -----------------------------------------------------------
# ⚙️ Prepare new (IRAP + FRAP) data for projection
# -----------------------------------------------------------
train_means <- colMeans(expr_train_200)
stopifnot(identical(names(train_means), sel))

new_samp_order <- c(old_sedirap_samples, old_sedfrap_samples)
expr_new_200 <- matrix(rep(train_means, each = length(new_samp_order)),
                       nrow = length(new_samp_order), byrow = TRUE,
                       dimnames = list(new_samp_order, sel))

present_genes_new <- intersect(sel, rownames(logCPM_symbol_new))
missing_genes_new <- setdiff(sel, present_genes_new)
cat("New data: present selected genes:", length(present_genes_new), " / 200\n")
if (length(missing_genes_new) > 0) {
  cat("Missing genes in new data (filled with training mean):\n")
  print(missing_genes_new)
}

if (length(present_genes_new) > 0) {
  expr_new_200[, present_genes_new] <- t(logCPM_symbol_new[present_genes_new, new_samp_order, drop = FALSE])
}

stopifnot(identical(colnames(expr_new_200), sel))
stopifnot(identical(rownames(expr_new_200), new_samp_order))

# -----------------------------------------------------------
# 🚀 Project onto training axis
# -----------------------------------------------------------
proj_new_200 <- predict(splsda_200, newdata = expr_new_200)

scores_new <- proj_new_200$variates[, 1]
samples_new <- rownames(proj_new_200$variates)

proj_df <- data.frame(
  Sample = samples_new,
  Comp1  = as.numeric(scores_new),
  Group  = ifelse(samples_new %in% old_sedirap_samples, "OLD_SED_IRAP", "OLD_SED_FRAP")
)

head(proj_df)

# -----------------------------------------------------------
# 🎨 Combine and plot
# -----------------------------------------------------------
plot_df <- bind_rows(
  train_df %>% mutate(Source = "Training"),
  proj_df %>% mutate(Source = "Projected")
)

ggplot(plot_df, aes(x = Group, y = Comp1, color = Source, shape = Source)) +
  geom_jitter(width = 0.1, size = 3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1, width = 0.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "Projection of OLD_SED_IRAP / OLD_SED_FRAP onto YNG–OLD RNA Axis (200 features)",
    y = "sPLS-DA Comp1 (RNA Aging Axis)",
    x = "Group"
  ) +
  scale_color_manual(values = c("Training" = "gray40", "Projected" = "steelblue"))

# -----------------------------------------------------------
# 💾 Export selected genes
# -----------------------------------------------------------
# Extract the table of selected variables with loadings
selected_tbl_200 <- mixOmics::selectVar(splsda_200, comp = 1)$value

# Convert to a tidy data frame with clear column names
selected_genes_200 <- data.frame(
  Gene = rownames(selected_tbl_200),
  Loading = selected_tbl_200[, 1]
) %>%
  arrange(desc(abs(Loading)))

# Preview top contributors
head(selected_genes_200, 10)

# Save to CSV
write.csv(selected_genes_200, "sPLSDA_YNGvsOLD_selected200.csv", row.names = FALSE)
#===================================================================================#


#=====================================================================#
#------Project all 5 groups onto old vs yng sPLS-DA(200) axis---------#
#==================================================================#

suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# -----------------------------------------------------------
# 1. Define sample groups
# -----------------------------------------------------------
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")

old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")

old_pwrveh_samples  <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

# sanity groups list for later looping
proj_groups <- list(
  OLD_SED_IRAP = old_sedirap_samples,
  OLD_SED_FRAP = old_sedfrap_samples,
  OLD_PWR_VEH  = old_pwrveh_samples,
  OLD_PWR_IRAP = old_pwrirap_samples,
  OLD_PWR_FRAP = old_pwrfrap_samples
)

# -----------------------------------------------------------
# 2. Build ONE merged logCPM matrix with all samples we care about
#    Rows = genes, Cols = samples
# -----------------------------------------------------------
# We take only genes present in ALL sources to stay safe
common_genes <- Reduce(
  intersect,
  list(
    rownames(logCPM_symbol_yngOLDsed),
    rownames(logCPM_symbol_new),
    rownames(logCPM_symbol_all)
  )
)

cat("Common genes across all matrices:", length(common_genes), "\n")

# subset each matrix to common genes
m_yngold   <- logCPM_symbol_yngOLDsed[common_genes, , drop = FALSE]
m_new      <- logCPM_symbol_new[common_genes, , drop = FALSE]
m_allpwr   <- logCPM_symbol_all[common_genes, , drop = FALSE]

# now bind columns. Some columns are duplicated between matrices (e.g. YNG/OLD sed are in both m_yngold and m_allpwr)
# we'll make a single combined matrix without duplicate columns:
expr_all <- cbind(
  m_yngold,
  m_new[, setdiff(colnames(m_new), colnames(m_yngold)), drop = FALSE],
  m_allpwr[, setdiff(colnames(m_allpwr), c(colnames(m_yngold), colnames(m_new))), drop = FALSE]
)

cat("Final merged expr_all dim (genes x samples):", paste(dim(expr_all), collapse=" x "), "\n")

# -----------------------------------------------------------
# 3. TRAINING SETUP: YNG_SED_VEH vs OLD_SED_VEH
# -----------------------------------------------------------
train_samples <- c(yng_sed_samples, old_sed_samples)

stopifnot(all(train_samples %in% colnames(expr_all)))

expr_train <- expr_all[, train_samples, drop = FALSE]

group_train <- factor(
  c(rep("YNG", length(yng_sed_samples)),
    rep("OLD", length(old_sed_samples))),
  levels = c("YNG","OLD")
)

cat("Training expr dim:", paste(dim(expr_train), collapse=" x "), "\n")
cat("Group table:\n"); print(table(group_train))

# -----------------------------------------------------------
# 4. Fit sPLS-DA with keepX = 200
# -----------------------------------------------------------
set.seed(123)
splsda_rna <- splsda(
  X = t(expr_train),   # samples as rows
  Y = group_train,
  ncomp = 1,
  keepX = 200
)

cat("Non-zero loadings Comp1:", sum(splsda_rna$loadings$X[,1] != 0), "\n")

# extract the actual 200 selected gene names
sel_tbl <- mixOmics::selectVar(splsda_rna, comp = 1)
sel_genes <- sel_tbl$name
cat("Selected genes length:", length(sel_genes), "\n")

# -----------------------------------------------------------
# 5. Refit model on JUST those 200 genes (locked feature space)
#    This stabilizes projection and matches what we did before
# -----------------------------------------------------------
expr_train_200 <- t(expr_train[sel_genes, , drop = FALSE])

stopifnot(identical(colnames(expr_train_200), sel_genes))
stopifnot(identical(rownames(expr_train_200), train_samples))

set.seed(123)
splsda_200 <- splsda(
  X = expr_train_200,
  Y = group_train,
  ncomp = 1,
  keepX = 200
)

# training scores for plotting
train_df <- data.frame(
  Sample = rownames(splsda_200$variates$X),
  Comp1  = splsda_200$variates$X[, 1],
  Group  = c(rep("YNG_SED_VEH", length(yng_sed_samples)),
             rep("OLD_SED_VEH", length(old_sed_samples))),
  Source = "Training",
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------
# 6. Helper to project ANY new group of samples
#    Uses train_means fallback for genes that might be missing
# -----------------------------------------------------------
train_means <- colMeans(expr_train_200)  # length 200, names(sel_genes)
stopifnot(identical(names(train_means), sel_genes))

project_group <- function(sample_ids, group_label) {
  
  # make sure these samples exist in expr_all
  stopifnot(all(sample_ids %in% colnames(expr_all)))
  
  # start a matrix (rows = samples, cols = sel_genes)
  expr_proj <- matrix(
    rep(train_means, each = length(sample_ids)),
    nrow = length(sample_ids),
    byrow = TRUE,
    dimnames = list(sample_ids, sel_genes)
  )
  
  # overwrite with real expression where we have it
  present_genes <- intersect(sel_genes, rownames(expr_all))
  expr_proj[, present_genes] <- t(expr_all[present_genes, sample_ids, drop = FALSE])
  
  # predict scores on Comp1
  proj_out <- predict(splsda_200, newdata = expr_proj)
  
  scores   <- proj_out$variates[, 1]
  samples  <- rownames(proj_out$variates)
  
  data.frame(
    Sample = samples,
    Comp1  = as.numeric(scores),
    Group  = group_label,
    Source = "Projected",
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------
# 7. Project ALL groups of interest
# -----------------------------------------------------------
proj_list <- list(
  project_group(old_sedirap_samples, "OLD_SED_IRAP"),
  project_group(old_sedfrap_samples, "OLD_SED_FRAP"),
  project_group(old_pwrveh_samples,  "OLD_PWR_VEH"),
  project_group(old_pwrirap_samples, "OLD_PWR_IRAP"),
  project_group(old_pwrfrap_samples, "OLD_PWR_FRAP")
)

proj_df <- bind_rows(proj_list)

# -----------------------------------------------------------
# 8. Combine and plot
# -----------------------------------------------------------
plot_df <- bind_rows(train_df, proj_df)

# nicer x-axis order for interpretation
plot_df$Group <- factor(
  plot_df$Group,
  levels = c("YNG_SED_VEH",
             "OLD_SED_VEH",
             "OLD_SED_IRAP",
             "OLD_SED_FRAP",
             "OLD_PWR_VEH",
             "OLD_PWR_IRAP",
             "OLD_PWR_FRAP")
)

ggplot(plot_df, aes(x = Group, y = Comp1, color = Source, shape = Source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.12, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.8) +
  theme_classic(base_size = 14) +
  labs(
    title = "Projection of treatment groups onto YNG vs OLD Sedentary RNA Aging Axis (200 genes)",
    y = "sPLS-DA Comp1 score (aging axis)",
    x = ""
  ) +
  scale_color_manual(values = c("Training" = "gray40", "Projected" = "steelblue")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

# -----------------------------------------------------------
# 9. Export ranked loadings for the 200 genes
# -----------------------------------------------------------
sel_values <- mixOmics::selectVar(splsda_200, comp = 1)$value
selected_genes_200 <- data.frame(
  Gene    = rownames(sel_values),
  Loading = sel_values[, 1],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(abs(Loading)))

head(selected_genes_200, 10)

write.csv(selected_genes_200,
          "sPLSDA_YNGvsOLDsed_200gene_axis_loadings.csv",
          row.names = FALSE)

# Also return plot_df if you want to look at it in the console
plot_df
#---------------------------------------------------#
#===================================================#

#====================================================#
#------Flipped Axis----------------------------------#
#====================================================#

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
})

# -----------------------------------------------------------
# 1️⃣  Flip axis so YNG < 0 and OLD > 0
# -----------------------------------------------------------
# Compute correlation between Comp1 and numeric group label (YNG=0, OLD=1)
cor_check <- cor(
  train_df$Comp1,
  as.numeric(train_df$Group == "OLD_SED_VEH")
)

if (cor_check < 0) {
  message("🔄 Flipping DIABLO axis so that OLD_SED_VEH = positive scores, YNG_SED_VEH = negative scores")
  
  # flip direction in both training and projected sets
  train_df$Comp1 <- -train_df$Comp1
  proj_df$Comp1  <- -proj_df$Comp1
}

# Combine flipped data
plot_df <- bind_rows(train_df, proj_df)

# -----------------------------------------------------------
# 2️⃣  Set desired group order
# -----------------------------------------------------------
plot_df$Group <- factor(
  plot_df$Group,
  levels = c(
    "YNG_SED_VEH",
    "OLD_SED_VEH",
    "OLD_PWR_VEH",
    "OLD_PWR_IRAP",
    "OLD_PWR_FRAP",
    "OLD_SED_IRAP",
    "OLD_SED_FRAP"
  )
)

# -----------------------------------------------------------
# 3️⃣  Plot updated axis
# -----------------------------------------------------------
ggplot(plot_df, aes(x = Group, y = Comp1, color = Source, shape = Source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.12, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.8) +
  theme_classic(base_size = 14) +
  labs(
    title = "Projection of treatment groups onto YNG–OLD Sedentary RNA Aging Axis (flipped orientation)",
    y = "sPLS-DA Comp1 score (RNA aging axis)",
    x = ""
  ) +
  scale_color_manual(values = c("Training" = "gray40", "Projected" = "steelblue")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))
#---------------------------------------------------------------#

group_colors <- c(
  "YNG_SED_VEH" = "#6baed6",  # baby blue
  "OLD_SED_VEH" = "#969696",  # grey
  "OLD_PWR_VEH" = "#1b9e77",  # strong blue
  "OLD_PWR_IRAP" = "#cc79a7", # violet
  "OLD_PWR_FRAP" = "#9e1f63", # magenta
  "OLD_SED_IRAP" = "#E69F00", # amber/orange
  "OLD_SED_FRAP" = "#D55E00"  # red-orange
)

group_colors <- c(
  "YNG_SED_VEH" = "#4DB8FF",
  "OLD_SED_VEH" = "#636363",
  "OLD_PWR_VEH" = "#0072B2",
  "OLD_PWR_IRAP" = "#A64AC9",
  "OLD_PWR_FRAP" = "#E60073",
  "OLD_SED_IRAP" = "#FF9E00",
  "OLD_SED_FRAP" = "#D73027"
)
p <- ggplot(plot_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, color = "black") +
  geom_jitter(width = 0.15, size = 2.8, alpha = 0.8, color = "black") +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "RNA aging axis (YNG–OLD Sedentary baseline, projected treatments)",
    y = "sPLS-DA Comp1 (flipped, OLD positive)",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

pdf("RNA_AgingAxis_YNG_OLD_projection.pdf.pdf", width = 8, height = 8)
print(p)
dev.off()
# -----------------------------------------------------------
# 4️⃣  Run ANOVA + Tukey post-hoc
# -----------------------------------------------------------
anova_res <- plot_df %>%
  anova_test(Comp1 ~ Group)

tukey_res <- plot_df %>%
  tukey_hsd(Comp1 ~ Group)

# Print results
cat("\n===== ANOVA Results =====\n")
print(anova_res)
cat("\n===== Tukey Post-hoc Results =====\n")
print(tukey_res, n=21)

# Optional compact letter display on plot
pvals <- tukey_res %>%
  multcompLetters4(anova_res, .) %>%
  .$Letters %>%
  as.data.frame() %>%
  rownames_to_column("Group")

# merge letters
plot_df_letters <- plot_df %>%
  group_by(Group) %>%
  summarise(mean = mean(Comp1)) %>%
  left_join(pvals, by = "Group")

ggplot(plot_df, aes(x = Group, y = Comp1, color = Source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.12, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.8) +
  geom_text(
    data = plot_df_letters,
    aes(x = Group, y = mean + 0.5, label = Letters),
    color = "black", size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Group differences on RNA aging axis (ANOVA + Tukey)",
    y = "sPLS-DA Comp1 (flipped, OLD positive)",
    x = ""
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

musage.vec
pls_200 <- selected_genes_200 %>%
  pull(Gene)
intersect(musage.vec, pls_200)
#-----------------------------------------#

# -----------------------------------------------------------
# 📦 Load libraries
# -----------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------------------------------------
# 🧫 Define the gene set (using your selected genes)
# -----------------------------------------------------------
selected_symbols <- selected_genes_200$Gene

# Map to Entrez IDs (required by clusterProfiler)
gene_map <- bitr(selected_symbols,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

# -----------------------------------------------------------
# 🧭 GO Biological Process enrichment
# -----------------------------------------------------------
ego_bp <- enrichGO(
  gene          = gene_map$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",           # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE            # converts back to gene symbols
)

# -----------------------------------------------------------
# 📊 Summarize results
# -----------------------------------------------------------
# View top enriched terms
head(ego_bp, 10)

# Save results to CSV
write.csv(as.data.frame(ego_bp), "GO_BP_overrepresentation_selected200.csv", row.names = FALSE)

# -----------------------------------------------------------
# 🎨 Plot top terms
# -----------------------------------------------------------
dotplot(ego_bp, showCategory = 15, font.size = 12, title = "GO Biological Process Enrichment (200 sPLS-DA genes)")
#------------------------------------------------------------------------------#
#==============================================================================#

#-----------------------------------------------
# Cap the -log10(FDR) values for each dataset
#-----------------------------------------------

# Define overlap between MusAge and PLS-DA 200-gene set
overlap_genes <- intersect(musage.vec, pls_200)

# Shared limits
pls_limits <- c(0, 15)



plot_pls_200 <- oldsedveh_v_yngsed_genes %>%
  filter(Symbol %in% selected_symbols) %>%
  arrange(logFC) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol),
         neg_logFDR = -log10(FDR))  # transform for visualization

ggplot(plot_pls_200, aes(x = logFC, y = Symbol, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "plasma",
    direction = -1,
    limits = pls_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = "Gene Symbol",
    title = "Differential Expression by Gene"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
#---------------------------------------------------------------#

# Define overlap
overlap_genes <- intersect(musage.vec, pls_200)

# Shared fill scale limits
pls_limits <- c(0, 15)

# Build plotting data
plot_pls_200 <- oldsedveh_v_yngsed_genes %>%
  filter(Symbol %in% selected_symbols) %>%   # all 200 selected
  arrange(logFC) %>%
  mutate(
    idx = factor(dplyr::row_number(), levels = dplyr::row_number()), # force discrete y-axis
    neg_logFDR = -log10(FDR),
    is_overlap = Symbol %in% overlap_genes
  )

# Plot
pls_200_dotplot <- ggplot(plot_pls_200, aes(x = logFC, y = idx, fill = neg_logFDR)) +
  geom_col(color = "black", width = 0.8) +
  
  # Add labels only for overlap genes
  geom_text(
    data = plot_pls_200 %>% filter(is_overlap),
    aes(label = Symbol),
    hjust = ifelse(plot_pls_200$logFC[plot_pls_200$is_overlap] > 0, -0.2, 1.2),
    size = 3,
    color = "black"
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    direction = -1,
    limits = pls_limits,
    name = expression(-log[10](FDR))
  ) +
  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = NULL,
    title = "Differential Expression of sPLS-DA 200 Genes\n(MusAge overlap labeled)"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black")
  )

pdf("pls_200_dotplot.pdf", width = 6, height = 4)
print(pls_200_dotplot)
dev.off()
#------------------------------------------------------#

# Plot without black outlines
pls_200_dotplot <- ggplot(plot_pls_200, aes(x = logFC, y = idx, fill = neg_logFDR)) +
  geom_col(width = 0.8, color = NA) +   # remove outlines
  
  # Add labels only for overlap genes
  geom_text(
    data = plot_pls_200 %>% filter(is_overlap),
    aes(label = Symbol),
    hjust = ifelse(plot_pls_200$logFC[plot_pls_200$is_overlap] > 0, -0.2, 1.2),
    size = 3,
    color = "black"
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    direction = -1,
    limits = pls_limits,
    name = expression(-log[10](FDR))
  ) +
  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = NULL,
    title = "Differential Expression of sPLS-DA 200 Genes\n(MusAge overlap labeled)"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black")
  )

# Save high-quality version for Illustrator
ggsave("pls200_MusAge_overlap_dotplot.pdf", pls_200_dotplot,
       width = 8, height = 8, units = "in", dpi = 300)
