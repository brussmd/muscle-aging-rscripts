#Let's try DIABLO again,again....

# Install BiocManager if not already present
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

# Install edgeR (and any other Bioconductor packages you might need)
#BiocManager::install("edgeR")

# While we’re here, also grab org.Mm.eg.db and AnnotationDbi if not already installed
#BiocManager::install(c("org.Mm.eg.db", "AnnotationDbi"))

# Install BiocManager if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

# Install mixOmics
#BiocManager::install("mixOmics")


#============================================================#
#---------Build logCPM_symbol old and young------------------#
#============================================================#
setwd("/Users/mdbruss/Documents/RStudioProjects_2/Rapa_PwR")

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(tidyr)
  library(edgeR);  library(AnnotationDbi); library(org.Mm.eg.db)
})

# ---- Load counts ----
set01_os_vs_ys  <- read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")

# 1) Rename the first column to "Ensembl"
rename_first_col <- function(df) {
  cn <- colnames(df)
  cn[1] <- "Ensembl"
  colnames(df) <- cn
  df
}
set01_os_vs_ys <- rename_first_col(set01_os_vs_ys)

# ---- Define sample groups ----
yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
all_samples     <- c(yng_sed_samples, old_sed_samples)

# ---- Keep only those columns ----
merged_counts <- set01_os_vs_ys %>%
  dplyr::select(Ensembl, dplyr::all_of(all_samples))

# ---- Map Ensembl -> ENTREZ ----
merged_counts <- merged_counts %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db,
                           keys = Ensembl_noDec,
                           keytype = "ENSEMBL",
                           column = "ENTREZID",
                           multiVals = "first")) %>%
  tidyr::drop_na(ENTREZID)

# ---- Counts matrix ----
counts_matrix <- merged_counts %>%
  dplyr::select(dplyr::all_of(all_samples)) %>%
  replace(is.na(.), 0) %>%
  as.matrix()
rownames(counts_matrix) <- merged_counts$ENTREZID

# ---- Define groups ----
group <- factor(c(rep("YNG_SED", length(yng_sed_samples)),
                  rep("OLD_SED", length(old_sed_samples))),
                levels = c("YNG_SED","OLD_SED"))

stopifnot(length(group) == ncol(counts_matrix))

# ---- edgeR: TMM normalization, filtering, logCPM ----
dge <- DGEList(counts = counts_matrix, group = group)
dge <- calcNormFactors(dge, method = "TMM")
keep <- filterByExpr(dge)                   
dge <- dge[keep, , keep.lib.sizes = FALSE]
logCPM_matrix <- cpm(dge, log = TRUE, prior.count = 1)

cat("Genes kept:", nrow(logCPM_matrix), " Samples:", ncol(logCPM_matrix), "\n")

# ---- Map ENTREZ -> SYMBOL and collapse duplicates ----
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

# ---- Keep exact sample order ----
logCPM_symbol <- logCPM_symbol[, all_samples, drop = FALSE]

# ---- Sanity checks ----
stopifnot(!anyNA(logCPM_symbol))
stopifnot(identical(colnames(logCPM_symbol), all_samples))

saveRDS(logCPM_symbol, "logCPM_symbol_yngOLDsed.RDS")
dim(logCPM_symbol)
head(logCPM_symbol)
#----------------------------------------------------------#

#============================================================#
#---------Build metab_log old and young------------------#
#============================================================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

# ---- Load data ----
pwr_rapa_muscle_metabo_data <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)

pwr_rapa_muscle_metabo_data

# ---- Drop Group column, keep Sample + metabolites ----
metabo_wide <- pwr_rapa_muscle_metabo_data %>% dplyr::select(-Group)

# ---- Convert long → wide: rows = metabolites, cols = samples ----
metabo_matrix <- metabo_wide %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()

rownames(metabo_matrix) <- metabo_matrix$Metabolite
metabo_matrix$Metabolite <- NULL

# ---- Map metabolomics sample names (YS1/OS1/…) → T_* tube codes ----
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_06","OS1","T_09","OS2","T_13","OS3","T_19","OS4","T_21","OS5","T_25","OS6","T_29","OS7","T_39","OS8",
  "T_05","YS1","T_15","YS2","T_17","YS3","T_23","YS4","T_26","YS5","T_32","YS6","T_35","YS7","T_38","YS8"
)

name_map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
new_names <- name_map[colnames(metabo_matrix)]
colnames(metabo_matrix) <- ifelse(is.na(new_names), colnames(metabo_matrix), new_names)

# ---- Keep only the YNG + OLD SED samples ----
yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
all_samples     <- c(yng_sed_samples, old_sed_samples)

present_metabo <- intersect(all_samples, colnames(metabo_matrix))
missing_for_metabo <- setdiff(all_samples, colnames(metabo_matrix))
if (length(missing_for_metabo)) message("Metabolomics missing: ", paste(missing_for_metabo, collapse = ", "))

metabo_matrix_subset <- metabo_matrix[, present_metabo, drop = FALSE]

# ---- Ensure numeric ----
metabo_matrix_subset[] <- lapply(metabo_matrix_subset, function(x) as.numeric(as.character(x)))

# ---- Impute zeros/NA with half min positive per metabolite ----
metabo_imputed <- t(apply(metabo_matrix_subset, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))

# ---- Log2 transform ----
metabo_log <- as.data.frame(log2(metabo_imputed))
colnames(metabo_log) <- present_metabo
rownames(metabo_log) <- rownames(metabo_matrix_subset)

# ---- Align with RNA samples if logCPM_symbol exists ----
if (exists("logCPM_symbol")) {
  common_samples <- intersect(all_samples, colnames(logCPM_symbol))
  metabo_log     <- metabo_log[, common_samples, drop = FALSE]
  logCPM_symbol  <- logCPM_symbol[, common_samples, drop = FALSE]
  stopifnot(identical(colnames(metabo_log), colnames(logCPM_symbol)))
}

# ---- Save ----
saveRDS(metabo_log, "metabo_log_yngOLDsed.RDS")

# ---- Quick peek ----
dim(metabo_log)
head(metabo_log[, 1:4])
#------------------------------------------------------#

#============================================================#
#---------Build lipid_matrix old and young------------------#
#============================================================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble)
})

# ---- Load lipidomics ----
lipid_df <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)
# Keep everything, including YS6

# ---- Identify lipid columns ----
lipid_cols <- setdiff(names(lipid_df), c("Sample","Group"))

# ---- Impute zeros/NA per lipid, then log2 ----
lipid_df[ , lipid_cols] <- lapply(lipid_df[ , lipid_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})

# ---- Build matrix: rows = lipids, cols = samples ----
lipid_matrix <- t(as.matrix(lipid_df[ , lipid_cols]))
colnames(lipid_matrix) <- lipid_df$Sample

# ---- Remove YV samples (young vehicle, not needed here) ----
lipid_matrix <- lipid_matrix[ , !grepl("^YV", colnames(lipid_matrix)), drop = FALSE]

# ---- Map sample names (e.g. OS1, YS1) → T_* codes ----
sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "T_06","OS1","T_09","OS2","T_13","OS3","T_19","OS4","T_21","OS5","T_25","OS6","T_29","OS7","T_39","OS8",
  "T_05","YS1","T_15","YS2","T_17","YS3","T_23","YS4","T_26","YS5","T_32","YS6","T_35","YS7","T_38","YS8"
)

map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
mapped <- map[colnames(lipid_matrix)]
colnames(lipid_matrix) <- ifelse(is.na(mapped), colnames(lipid_matrix), mapped)

# ---- Keep only YNG + OLD SED samples ----
yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
all_samples     <- c(yng_sed_samples, old_sed_samples)

present <- intersect(all_samples, colnames(lipid_matrix))
lipid_log <- lipid_matrix[ , present, drop = FALSE]

# ---- Align with RNA samples if logCPM_symbol exists ----
if (exists("logCPM_symbol")) {
  common <- intersect(all_samples, colnames(logCPM_symbol))
  lipid_log     <- lipid_log[, common, drop = FALSE]
  logCPM_symbol <- logCPM_symbol[, common, drop = FALSE]
  stopifnot(identical(colnames(lipid_log), colnames(logCPM_symbol)))
}

# ---- Save ----
saveRDS(lipid_log, "lipid_log_yngOLDsed.RDS")

# ---- Sanity check ----
dim(lipid_log)
head(lipid_log[, 1:4])
#-----------------------------------------------#

library(mixOmics)

# Build omics list: samples in rows, features in columns
X <- list(
  RNA        = t(logCPM_symbol),
  Metabolite = t(metabo_log),
  Lipid      = t(lipid_log)
)

Y <- factor(c(
  rep("YNG", length(yng_sed_samples)),
  rep("OLD", length(old_sed_samples))
), levels = c("YNG","OLD"))

# Assign sample names to match the omics
names(Y) <- rownames(X$RNA)

# Check again
identical(rownames(X$RNA), names(Y))

yng_sed_samples
old_sed_samples
Y
# Check alignment
stopifnot(identical(rownames(X$RNA), rownames(X$Metabolite)))
stopifnot(identical(rownames(X$RNA), rownames(X$Lipid)))
stopifnot(identical(rownames(X$RNA), names(Y)))

#----------------------------------------------#
#----------check names-------------------------#
#----------------------------------------------#

rownames(X$RNA)
rownames(X$Metabolite)
rownames(X$Lipid)
names(Y)
#-----------------------------------#

#---Design matrix-------------------#
# Design matrix (full integration)
design <- matrix(1, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0
design

# Try with 1 component and 30 features per omics block
diablo_res <- block.splsda(X, Y, ncomp = 1,
                           keepX = list(RNA = 30, Metabolite = 20, Lipid = 20),
                           design = design)

plotDiablo(diablo_res)  # quick overview
#-----------------------------------------#

diablo_res$loadings$RNA
diablo_res$loadings$Metabolite
diablo_res$loadings$Lipid


selectVar(diablo_res, block = "RNA", comp = 1)$name
selectVar(diablo_res, block = "Metabolite", comp = 1)$name
selectVar(diablo_res, block = "Lipid", comp = 1)$name

# RNA features selected on component 1
rna_sel <- selectVar(diablo_res, block = "RNA", comp = 1)

# Metabolites selected
met_sel <- selectVar(diablo_res, block = "Metabolite", comp = 1)

# Lipids selected
lip_sel <- selectVar(diablo_res, block = "Lipid", comp = 1)

str(rna_sel)
#----------------------------#

# RNA
rna_features <- rna_sel$RNA$name
rna_loadings <- rna_sel$RNA$value

# Metabolites
met_features <- met_sel$Metabolite$name
met_loadings <- met_sel$Metabolite$value

# Lipids
lip_features <- lip_sel$Lipid$name
lip_loadings <- lip_sel$Lipid$value

library(dplyr)

# helper
extract_features <- function(sel, block) {
  data.frame(
    block = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    stringsAsFactors = FALSE
  )
}

selected_df <- bind_rows(
  extract_features(rna_sel, "RNA"),
  extract_features(met_sel, "Metabolite"),
  extract_features(lip_sel, "Lipid")
)

head(selected_df)
dim(selected_df)

selected_df

circosPlot(diablo_res, cutoff = 0.7)   # cross-block correlation network
network(diablo_res, comp = 1, cutoff = 0.6)  # relevance network
plotVar(diablo_res, comp = 1, var.names = TRUE) # variable plots

#-------------------------------------------------#
#----Run Diablo with X (all features) and Y(sample definitions)
#--------------------------------------------------#
diablo_res <- block.splsda(X, Y, ncomp = 1,
                           keepX = list(RNA = 50,
                                        Metabolite = 20,
                                        Lipid = 20),
                           design = design)

# Extract selected features per block
rna_sel <- selectVar(diablo_res, block = "RNA", comp = 1)
met_sel <- selectVar(diablo_res, block = "Metabolite", comp = 1)
lip_sel <- selectVar(diablo_res, block = "Lipid", comp = 1)

extract_features <- function(sel, block) {
  data.frame(
    block   = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    stringsAsFactors = FALSE
  )
}

selected_df <- bind_rows(
  extract_features(rna_sel, "RNA"),
  extract_features(met_sel, "Metabolite"),
  extract_features(lip_sel, "Lipid")
)

# Sort by block and loading (optional)
selected_df <- selected_df %>%
  arrange(block, desc(abs(loading)))

head(selected_df, 15)  # preview top 15

selected_df

# ---- Extract & save selected features ----
library(dplyr)
library(mixOmics)

# Helper function
extract_features <- function(sel, block) {
  data.frame(
    block   = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    stringsAsFactors = FALSE
  )
}

# 1. Extract selected features for comp 1
rna_sel <- selectVar(diablo_res, block = "RNA", comp = 1)
met_sel <- selectVar(diablo_res, block = "Metabolite", comp = 1)
lip_sel <- selectVar(diablo_res, block = "Lipid", comp = 1)

# 2. Combine into a single dataframe
selected_df <- bind_rows(
  extract_features(rna_sel, "RNA"),
  extract_features(met_sel, "Metabolite"),
  extract_features(lip_sel, "Lipid")
) %>%
  arrange(block, desc(abs(loading)))

# 3. Save to CSV for inspection
write.csv(selected_df, "DIABLO_selected_features_comp1_muscle.csv", row.names = FALSE)

# Preview top features
print(head(selected_df, 20))


# ---- Visualization ----

# Circos plot of cross-block correlations
circosPlot(diablo_res, comp = 1, cutoff = 0.7, line = TRUE)

# Pairwise RNA–Metabolite network
network(diablo_res, comp = list(RNA = 1, Metabolite = 1), cutoff = 0.6)

# Pairwise RNA–Lipid network
network(diablo_res, comp = list(RNA = 1, Lipid = 1), cutoff = 0.6)

# Pairwise Metabolite–Lipid network
network(diablo_res, comp = list(Metabolite = 1, Lipid = 1), cutoff = 0.6)

# Global (all 3 blocks together)
network(diablo_res, comp = list(RNA = 1, Metabolite = 1, Lipid = 1), cutoff = 0.6)

# Global network (all three blocks together)
network(diablo_res, comp = list(1,1), cutoff = 0.6)

# Loadings plots for each block
plotLoadings(diablo_res, block = "RNA", comp = 1, method = "mean")
plotLoadings(diablo_res, block = "Metabolite", comp = 1, method = "mean")
plotLoadings(diablo_res, block = "Lipid", comp = 1, method = "mean")

network(diablo_res, 1, cutoff = 0.6)
#-------------------------------------------------#

# --- Helper: extract loadings with absolute value ---
extract_loadings <- function(sel, block){
  data.frame(
    block   = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    abs_loading = abs(sel[[block]]$value$value.var),
    stringsAsFactors = FALSE
  )
}

# Extract selected features for each block
rna_sel <- selectVar(diablo_res, block = "RNA", comp = 1)
met_sel <- selectVar(diablo_res, block = "Metabolite", comp = 1)
lip_sel <- selectVar(diablo_res, block = "Lipid", comp = 1)

rna_sel

selected_df <- bind_rows(
  extract_loadings(rna_sel, "RNA"),
  extract_loadings(met_sel, "Metabolite"),
  extract_loadings(lip_sel, "Lipid")
)

selected_df

# --- Rank features within each block ---
ranked_df <- selected_df %>%
  group_by(block) %>%
  arrange(desc(abs_loading), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# --- Plot: histogram of loadings per block ---
ggplot(ranked_df, aes(x = abs_loading)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~block, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of absolute loadings",
       x = "Absolute loading",
       y = "Count")

# --- Plot: ranked loadings (top features highlighted) ---
ggplot(ranked_df, aes(x = rank, y = abs_loading, color = block)) +
  geom_point() +
  geom_line(aes(group = block)) +
  theme_minimal(base_size = 14) +
  labs(title = "Ranked feature loadings per block",
       x = "Feature rank within block",
       y = "Absolute loading") +
  scale_color_brewer(palette = "Dark2")

strong_feats <- ranked_df %>%
  filter(abs_loading > 0.2) %>%
  arrange(block, desc(abs_loading))
print(strong_feats, n = 22)

library(ggrepel)

ggplot(strong_feats, aes(x = reorder(feature, abs_loading), 
                         y = abs_loading, fill = block)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = round(abs_loading,3)), 
            hjust = -0.1, size = 3) +
  facet_wrap(~block, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Top drivers (|loading| > 0.2)",
       x = "Feature", y = "Absolute loading")
#------------------------------------------------#

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
})

# Background = all genes expressed in your RNA-seq (SYMBOLs from logCPM_symbol)
bg_genes <- rownames(logCPM_symbol)   

# Map to ENTREZ
bg_entrez <- mapIds(org.Mm.eg.db,
                    keys       = bg_genes,
                    column     = "ENTREZID",
                    keytype    = "SYMBOL",
                    multiVals  = "first") %>%
  unique() %>% na.omit()

cat("Universe size:", length(bg_entrez), "\n")

# Your DIABLO-selected features
rna_features <- rna_sel$RNA$name

entrez_ids <- mapIds(org.Mm.eg.db,
                     keys       = rna_features,
                     column     = "ENTREZID",
                     keytype    = "SYMBOL",
                     multiVals  = "first") %>%
  unique() %>% na.omit()

cat("DIABLO RNA features mapped to ENTREZ:", length(entrez_ids), "\n")

# Very relaxed cutoff
ego_relaxed <- enrichGO(gene          = entrez_ids,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = "ENTREZID",
                        ont           = "BP",
                        universe      = bg_entrez,
                        pvalueCutoff  = 1,    # no p-value restriction
                        qvalueCutoff  = 1,    # no q-value restriction
                        pAdjustMethod = "BH")

head(ego_relaxed)
#-----------------------------------------------------------#

#============================================================#
#---------DIABLO Ingegrated Aging Axis------------------#
#============================================================#

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(mixOmics)
})

# you already have these from earlier work:
# diablo_res, yng_sed_samples, old_sed_samples,
# old_pwr_samples, old_pwrirap_samples, old_pwrfrap_samples

#----------------------------------------------------------#
#----Build All Symbol logCPM-------------------------------#
# Rows = genes (SYMBOL), Cols = samples (T_*); ALL groups included
logCPM_symbol_all <- readRDS("logCPM_symbol_final.RDS")

head(logCPM_symbol_all)

# Training feature set used by DIABLO (RNA block)
train_feats_RNA <- colnames(t(logCPM_symbol))    # features = gene SYMBOLs
# Align columns to training features; drop extras, fill missing with NA if any (shouldn’t happen)
RNA_all <- t(logCPM_symbol_all)
RNA_all <- RNA_all[ , intersect(colnames(RNA_all), train_feats_RNA), drop = FALSE]
RNA_all <- RNA_all[ , train_feats_RNA, drop = FALSE]   # enforce same order

# Find overlap between training features and all features
common_genes <- intersect(train_feats_RNA, colnames(RNA_all))
missing_genes <- setdiff(train_feats_RNA, colnames(RNA_all))

cat("Number of training genes:", length(train_feats_RNA), "\n")
cat("Number of genes in all data:", ncol(RNA_all), "\n")
cat("Number of overlapping genes:", length(common_genes), "\n")
cat("Number of missing genes:", length(missing_genes), "\n")

# Preview a few missing ones
head(missing_genes, 20)

# 50 DIABLO-selected RNA features
rna_features_50 <- rna_sel$RNA$name

# Check against rownames, NOT colnames
common_rna_feats <- intersect(rna_features_50, rownames(logCPM_symbol_all))
missing_rna_feats <- setdiff(rna_features_50, rownames(logCPM_symbol_all))

cat("DIABLO-selected RNA features:", length(rna_features_50), "\n")
cat("Overlap with full dataset:", length(common_rna_feats), "\n")
cat("Missing features:", length(missing_rna_feats), "\n")

# Subset RNA_all correctly: rows = samples, cols = genes
RNA_all <- t(logCPM_symbol_all[common_rna_feats, , drop = FALSE])

logCPM_symbol_all
# Sanity check RNA_all
dim(RNA_all)                # samples x 50 features
head(rownames(RNA_all))     # should be T_* sample IDs
head(colnames(RNA_all))     # should be your 50 gene symbols
is.numeric(RNA_all[1,1])    # should be TRUE

# Number of samples vs features
cat("Samples:", nrow(RNA_all), " | Features:", ncol(RNA_all), "\n")
#---------------------------------------------------------------------#

#-----Rebuild full (all sample) metabo_log---------------------------#

# ---------- Metabolomics log matrix (all samples, all features) ---------- #

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

# 0) Load metabolomics data
metabo_raw <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)

# 1) Keep Sample + features (Group not needed here)
metabo_wide <- metabo_raw %>% dplyr::select(-Group)

# 2) Reshape: rows = metabolites, cols = samples
metabo_matrix <- metabo_wide %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()

rownames(metabo_matrix) <- metabo_matrix$Metabolite
metabo_matrix$Metabolite <- NULL

# 3) Map sample names (YS1/OS1/OV1/OIR1/OFR1/YV1 → T_* codes)
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
new_names <- name_map[colnames(metabo_matrix)]
colnames(metabo_matrix) <- ifelse(is.na(new_names), colnames(metabo_matrix), new_names)

# 4) Ensure numeric and impute zeros/NA per metabolite, then log2-transform
metabo_matrix[] <- lapply(metabo_matrix, function(x) as.numeric(as.character(x)))

metabo_imputed <- t(apply(metabo_matrix, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))

metabo_log_all <- as.data.frame(log2(metabo_imputed))

# 5) Subset to the 20 DIABLO-selected metabolites
met_features_20 <- met_sel$Metabolite$name
common_met_feats <- intersect(met_features_20, rownames(metabo_log_all))
missing_met_feats <- setdiff(met_features_20, rownames(metabo_log_all))

cat("DIABLO-selected metabolites:", length(met_features_20), "\n")
cat("Overlap with full dataset:", length(common_met_feats), "\n")
cat("Missing metabolites:", length(missing_met_feats), "\n")

# Final matrix for DIABLO: rows = samples, cols = 20 metabolites
Metabo_all <- t(metabo_log_all[common_met_feats, , drop = FALSE])

# Sanity check
dim(Metabo_all)            # should be 37 x 20
head(rownames(Metabo_all)) # T_* sample IDs
head(colnames(Metabo_all)) # metabolite names
is.numeric(Metabo_all[1,1])

rownames(Metabo_all)
rownames(RNA_all)

# Ensure Metabo_all has the same samples (rows) as RNA_all
common_samples <- intersect(rownames(RNA_all), rownames(Metabo_all))

# Subset and reorder to match RNA_all
Metabo_all <- Metabo_all[common_samples, , drop = FALSE]
Metabo_all <- Metabo_all[rownames(RNA_all), , drop = FALSE]

# Sanity check
cat("Samples in RNA_all:", nrow(RNA_all), "\n")
cat("Samples in Metabo_all after alignment:", nrow(Metabo_all), "\n")
stopifnot(identical(rownames(RNA_all), rownames(Metabo_all)))
#--------------------------------------------------------------#

#----Rebuild full lipid_matrix-----------------------------#

#------- Build Lipid_all (all samples, keep YS6) -----------------------#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble)
})

# 0) Load lipidomics data (keep all samples, including YS6)
lipid_df <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)

# 1) Identify lipid columns
lipid_cols <- setdiff(names(lipid_df), c("Sample","Group"))

# 2) Impute zeros/NA per lipid, then log2-transform
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

# 4) Drop YV samples (not part of design)
lipid_matrix <- lipid_matrix[ , !grepl("^YV", colnames(lipid_matrix)), drop = FALSE]

# 5) Rename sample IDs (e.g. OS1, YS1 → T_* codes)
map <- setNames(sample_key$TubeCode, sample_key$MetabolomicsName)
mapped <- map[colnames(lipid_matrix)]
colnames(lipid_matrix) <- ifelse(is.na(mapped), colnames(lipid_matrix), mapped)

# 6) Subset to the 20 DIABLO-selected lipid features
lipid_features_20 <- lip_sel$Lipid$name
common_lip_feats  <- intersect(lipid_features_20, rownames(lipid_matrix))
missing_lip_feats <- setdiff(lipid_features_20, rownames(lipid_matrix))

cat("DIABLO-selected lipids:", length(lipid_features_20), "\n")
cat("Overlap with dataset:", length(common_lip_feats), "\n")
cat("Missing:", length(missing_lip_feats), "\n")

# Final Lipid_all matrix: rows = samples, cols = 20 lipids
Lipid_all <- t(lipid_matrix[common_lip_feats, , drop = FALSE])

# 7) Align Lipid_all rows (samples) with RNA_all
common_samples <- intersect(rownames(RNA_all), rownames(Lipid_all))
Lipid_all <- Lipid_all[common_samples, , drop = FALSE]
Lipid_all <- Lipid_all[rownames(RNA_all), , drop = FALSE]  # enforce same order

# Sanity check
cat("Samples in RNA_all:", nrow(RNA_all), "\n")
cat("Samples in Lipid_all after alignment:", nrow(Lipid_all), "\n")
stopifnot(identical(rownames(RNA_all), rownames(Lipid_all)))

dim(Lipid_all)
head(colnames(Lipid_all))
head(Lipid_all)
#------------------------------------------------#

#----Check structure of full DIABLO data set-----#

# ---- Final alignment checks before DIABLO Aging Axis ----

# 1) Dimensions
cat("RNA_all:     ", dim(RNA_all)[1], "samples x", dim(RNA_all)[2], "features\n")
cat("Metabo_all:  ", dim(Metabo_all)[1], "samples x", dim(Metabo_all)[2], "features\n")
cat("Lipid_all:   ", dim(Lipid_all)[1], "samples x", dim(Lipid_all)[2], "features\n")

# 2) Samples must match across all sets
stopifnot(identical(rownames(RNA_all), rownames(Metabo_all)))
stopifnot(identical(rownames(RNA_all), rownames(Lipid_all)))

# 3) Quick preview
head_samples <- head(rownames(RNA_all))
cat("\nFirst few samples:\n")
print(head_samples)

cat("\nRNA features (first 5):\n"); print(head(colnames(RNA_all)))
cat("\nMetabolite features (first 5):\n"); print(head(colnames(Metabo_all)))
cat("\nLipid features (first 5):\n"); print(head(colnames(Lipid_all)))

# 4) Numeric check (all values should be numeric)
cat("\nData type check:\n")
print(is.numeric(RNA_all[1,1]))
print(is.numeric(Metabo_all[1,1]))
print(is.numeric(Lipid_all[1,1]))
#-----------------------------------#

#----DIABLO Aging Axis--------------#

#===========================================================
# DIABLO Aging Axis Projection (RNA + Metabolite + Lipid)
#===========================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
})

#-----------------------------------------------------------
# 0) Define sample groups
#-----------------------------------------------------------
if (!exists("yng_sed_samples", inherits = FALSE))    
  yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")

if (!exists("old_sedveh_samples", inherits = FALSE)) 
  old_sedveh_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")

if (!exists("old_pwrveh_samples", inherits = FALSE)) 
  old_pwrveh_samples <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")

if (!exists("old_pwrirap_samples", inherits = FALSE)) 
  old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")

if (!exists("old_pwrfrap_samples", inherits = FALSE)) 
  old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

#-----------------------------------------------------------
# 1) Build multi-omics list (already preprocessed above)
#    Rows = samples, Cols = features
#-----------------------------------------------------------
X_all <- list(
  RNA        = RNA_all,
  Metabolite = Metabo_all,
  Lipid      = Lipid_all
)

#-----------------------------------------------------------
# 2) Define outcome Y for training (YNG vs OLD)
#    Training only uses young sedentary vs old sedentary vehicle
#-----------------------------------------------------------
Y_train <- factor(c(
  rep("YNG", length(yng_sed_samples)),
  rep("OLD", length(old_sedveh_samples))
), levels = c("YNG","OLD"))

# Restrict X to training samples for the DIABLO fit
train_samples <- c(yng_sed_samples, old_sedveh_samples)
X_train <- lapply(X_all, function(m) m[train_samples, , drop = FALSE])

#-----------------------------------------------------------
# 3) Run DIABLO
#-----------------------------------------------------------
# Full design (all blocks correlated)
design <- matrix(1, ncol = 3, nrow = 3, dimnames = list(names(X_all), names(X_all)))
diag(design) <- 0

diablo_res <- block.splsda(
  X_train,
  Y_train,
  ncomp = 1,
  keepX = list(RNA = 50, Metabolite = 20, Lipid = 20),
  design = design
)

#-----------------------------------------------------------
# 4) Project ALL samples into DIABLO space
#-----------------------------------------------------------
proj <- predict(diablo_res, newdata = X_all)  # scores for all samples
scores <- proj$variates$X   # list by block

#-----------------------------------------------------------
# Extract Comp1 scores for all samples across blocks
#-----------------------------------------------------------
scores_all <- lapply(proj$variates, function(block) block[,1])

# Combine into a matrix (rows = samples, cols = blocks)
scores_mat <- do.call(cbind, scores_all)

# Average across blocks to get one “consensus” axis per sample
comp1 <- rowMeans(scores_mat)

# Build dataframe with Sample + Comp1
proj_df <- data.frame(
  Sample = rownames(scores_mat),
  Comp1  = comp1
)

# Add Group labels
proj_df <- proj_df %>%
  mutate(Group = case_when(
    Sample %in% yng_sed_samples     ~ "Young_Sed",
    Sample %in% old_sedveh_samples  ~ "Old_SedVeh",
    Sample %in% old_pwrveh_samples  ~ "Old_PwrVeh",
    Sample %in% old_pwrirap_samples ~ "Old_PwrIRap",
    Sample %in% old_pwrfrap_samples ~ "Old_PwrFRap",
    TRUE                            ~ "Other"
  ))

ggplot(proj_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Projection of Samples onto DIABLO Aging Axis (Comp1)",
       y = "DIABLO Component 1 Score", x = "") +
  theme(legend.position = "none")
#---------------------------------------------------#

#------With Statistics------------------------------#

#-----------------------------------------------------------
# 5) ANOVA + Tukey to test group differences on Comp1
#-----------------------------------------------------------
anova_res <- aov(Comp1 ~ Group, data = proj_df)
summary(anova_res)

# Tukey HSD post-hoc comparisons
tukey_res <- TukeyHSD(anova_res)
print(tukey_res)

# Optional: tidy results for plotting
if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
library(broom)
tukey_df <- broom::tidy(tukey_res)

#-----------------------------------------------------------
# 6) Plot with stats layered
#-----------------------------------------------------------
ggplot(proj_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Projection of Samples onto DIABLO Aging Axis (Comp1)",
       y = "DIABLO Component 1 Score", x = "") +
  theme(legend.position = "none")
#-------------------------------------------------------------#

#-----Show plot with stats----------------------------#

suppressPackageStartupMessages({
  library(rstatix)
  library(ggpubr)
})

# Ensure factor order
proj_df <- proj_df %>%
  mutate(Group = factor(Group, 
                        levels = c("Young_Sed", "Old_SedVeh", "Old_PwrVeh", "Old_PwrIRap", "Old_PwrFRap")))

# Run one-way ANOVA
anova_res <- anova_test(data = proj_df, dv = Comp1, between = Group)
anova_res

# Tukey posthoc
tukey_res <- proj_df %>%
  tukey_hsd(Comp1 ~ Group)
tukey_res

# Define the pairwise comparisons to display
comparisons <- list(
  c("Old_SedVeh", "Young_Sed"),   # oldsed vs yngsed
  c("Old_SedVeh", "Old_PwrVeh"),  # oldsed vs oldpwr
  c("Old_PwrVeh", "Old_PwrIRap"), # oldpwr vs oldpwrirap
  c("Old_PwrVeh", "Old_PwrFRap")  # oldpwr vs oldpwrfrap
)

# Plot with significance bars
DIABLO_aging_axis <- ggplot(proj_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Projection of Samples onto DIABLO Aging Axis (Comp1)",
       y = "DIABLO Component 1 Score", x = "") +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test", 
                     label = "p.signif")

pdf("DIABLO_AgingAxis_Boxplot.pdf", width = 6, height = 4)
print(DIABLO_aging_axis)
dev.off()
#-----------------------------------------------------#
#======================================================#

# Inspect first few rows
head(proj_df)
#   Sample    Comp1      Group
# 1   T_05  -3.211   Young_Sed
# 2   T_06   2.314   Old_SedVeh
# ...

# Save as CSV for record-keeping or plotting in other tools
write.csv(proj_df, "DIABLO_AgingAxis_SampleScores.csv", row.names = FALSE)




#=======================================================#
#-----Plot DIABLO Features------------------------------#
#=======================================================#

Diablo_latent_features <- read.csv("DIABLO_selected_features_comp1.csv")
Diablo_latent_features_v2 <- read.csv("DIABLO_selected_features_comp1_muscle.csv")

Diablo_latent_features

rna_features <- Diablo_latent_features %>%
  dplyr::filter(block == "RNA") %>%
  dplyr::pull(feature)

rna_features


rna_features_v2 <- Diablo_latent_features_v2 %>%
  dplyr::filter(block == "RNA") %>%
  dplyr::pull(feature)

rna_features_v2


oldsedveh_v_yngsedveh_genes %>%
  filter(Symbol %in% rna_features) %>%
  arrange(logFC)

# Filter and sort your dataframe
plot_df <- oldsedveh_v_yngsedveh_genes %>%
  filter(Symbol %in% rna_features) %>%
  arrange(logFC)

# Reorder factor levels so the bars appear in order
plot_df$Symbol <- factor(plot_df$Symbol, levels = plot_df$Symbol)

plot_df <- oldsedveh_v_yngsedveh_genes %>%
  filter(Symbol %in% rna_features) %>%
  arrange(logFC) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol),
         neg_logFDR = -log10(FDR))  # transform for visualization

diablo_rna_feature_plot <- ggplot(plot_df, aes(x = logFC, y = Symbol, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "plasma",   # alternatives: "viridis", "magma", "inferno", "cividis"
    direction = -1,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
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
#--------------------------------------------------#

#---------Lipid Features---------------------------#

lipid_features_v2 <- Diablo_latent_features_v2 %>%
  dplyr::filter(block == "Lipid") %>%
  dplyr::pull(feature)

lipid_features_v2

lipid_features <- Diablo_latent_features %>%
  dplyr::filter(block == "Lipid") %>%
  dplyr::pull(feature)

lipid_features

# 0) Load
lip_raw <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)

# 1) Identify lipid columns (everything except Sample, Group)
lipid_cols <- setdiff(names(lip_raw), c("Sample","Group"))

# 2) Impute zeros with half of the per-lipid minimum non-zero (robust fallback)
min_nonzero <- sapply(lip_raw[lipid_cols], function(x) {
  m <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (is.infinite(m)) NA_real_ else m
})
# If any lipid had no >0 values, fall back to the global min nonzero across all lipids
global_min <- suppressWarnings(min(unlist(lip_raw[lipid_cols])[unlist(lip_raw[lipid_cols]) > 0], na.rm = TRUE))
min_nonzero[is.na(min_nonzero)] <- global_min

lip_imp <- lip_raw
for (lip in lipid_cols) {
  x <- lip_imp[[lip]]
  x[x == 0] <- min_nonzero[[lip]] / 2
  lip_imp[[lip]] <- x
}

# 3) Log2 transform
lip_log <- lip_imp %>% mutate(across(all_of(lipid_cols), log2))

# 4) Optional: drop YS6 (kept from your workflow; harmless if not present)
lip_log <- lip_log %>% filter(Sample != "YS6")

# 5) Subset to OS vs YS
lipid_os_ys <- lip_log %>% filter(Group %in% c("OS","YS"))

# 6) Per-lipid stats (log2FC = mean(OS) − mean(YS); t-test; FDR)
lipid_os_ys_stats <- purrr::map_dfr(lipid_cols, function(lip) {
  vals  <- lipid_os_ys[[lip]]
  grp   <- lipid_os_ys$Group
  tt    <- t.test(vals ~ grp)
  tibble(
    Lipid              = lip,
    Mean_YS            = mean(vals[grp == "YS"], na.rm = TRUE),
    Mean_OS            = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OS_vs_YS   = Mean_OS - Mean_YS,
    P_value            = tt$p.value
  )
}) %>%
  mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  arrange(FDR)

head(lipid_os_ys_stats)

lipid_os_ys_stats %>%
  filter(Lipid %in% lipid_features) %>%
  arrange(Log2_FC_OS_vs_YS)

# Filter and sort your dataframe
lipid_plot_df <- lipid_os_ys_stats %>%
  filter(Lipid %in% lipid_features) %>%
  arrange(Log2_FC_OS_vs_YS)

# Reorder factor levels so the bars appear in order
lipid_plot_df$Lipid <- factor(lipid_plot_df$Lipid, levels = lipid_plot_df$Lipid)

lipid_plot_df <- lipid_os_ys_stats %>%
  filter(Lipid %in% lipid_features) %>%
  arrange(Log2_FC_OS_vs_YS) %>%
  mutate(Lipid = factor(Lipid, levels = Lipid),
         neg_logFDR = -log10(FDR))  # transform for visualization

diablo_lipid_feature_plot <- ggplot(lipid_plot_df, aes(x = Log2_FC_OS_vs_YS, y = Lipid, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "plasma",   # alternatives: "viridis", "magma", "inferno", "cividis"
    direction = -1,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
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
#-----------------------------------------------------------------#

#---------Metabolite Features---------------------------#

metabolite_features_v2 <- Diablo_latent_features_v2 %>%
  dplyr::filter(block == "Metabolite") %>%
  dplyr::pull(feature)

metabolite_features_v2

metabolite_features <- Diablo_latent_features %>%
  dplyr::filter(block == "Metabolite") %>%
  dplyr::pull(feature)

metabolite_features

# 1) Load metabolomics (HILIC) and keep OS/YS; drop YS6 if present (to mirror lipidomics)
metabo_raw <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)

metabo_df_osys <- metabo_raw %>%
  dplyr::filter(Group %in% c("OS", "YS")) 

# 2) Long format (metabolites = columns other than Sample/Group)
metabo_cols <- setdiff(names(metabo_df_osys), c("Sample","Group"))
metabo_long_osys <- metabo_df_osys %>%
  tidyr::pivot_longer(cols = dplyr::all_of(metabo_cols),
                      names_to = "Metabolite", values_to = "Abundance") %>%
  tidyr::drop_na(Abundance)

# 3) Per-metabolite stats: Welch t-test OS vs YS; log2FC = log2(meanOS) − log2(meanYS)
metabo_stats_osys <- metabo_long_osys %>%
  dplyr::group_by(Metabolite) %>%
  dplyr::summarise(
    p_value = tryCatch(t.test(Abundance ~ Group)$p.value, error = function(e) NA_real_),
    mean_OS = mean(Abundance[Group == "OS"], na.rm = TRUE),
    mean_YS = mean(Abundance[Group == "YS"], na.rm = TRUE),
    log2FC  = log2(mean_OS + 1e-6) - log2(mean_YS + 1e-6),
    .groups = "drop"
  ) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(FDR)

# 4) Define Metabolite Aging Signature (p < 0.05, unadjusted, per spec)
metabo_age_signature <- metabo_stats_osys %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(p_value) %>%
  dplyr::pull(Metabolite)

metabo_age_signature
head(metabo_stats_osys)

metabo_stats_osys %>%
  filter(Metabolite %in% metabolite_features) %>%
  arrange(log2FC)

# Filter and sort your dataframe
metabo_plot_df <- metabo_stats_osys %>%
  filter(Metabolite %in% metabolite_features) %>%
  arrange(log2FC)

# Reorder factor levels so the bars appear in order
metabo_plot_df$Metabolite <- factor(metabo_plot_df$Metabolite, levels = metabo_plot_df$Metabolite)

metabolite_plot_df <- metabo_stats_osys %>%
  filter(Metabolite %in% metabolite_features) %>%
  arrange(log2FC) %>%
  mutate(Metabolite = factor(Metabolite, levels = Metabolite),
         neg_logFDR = -log10(FDR))  # transform for visualization








#=================================================================#
#-------Show plots with common scales-----------------------------#

#-----------------------------------------------
# Cap the -log10(FDR) values for each dataset
#-----------------------------------------------
plot_df <- plot_df %>%
  mutate(neg_logFDR = pmin(neg_logFDR, 5))

lipid_plot_df <- lipid_plot_df %>%
  mutate(neg_logFDR = pmin(neg_logFDR, 5))

metabolite_plot_df <- metabolite_plot_df %>%
  mutate(neg_logFDR = pmin(neg_logFDR, 5))

# Shared limits
common_limits <- c(0, 5)

#-----------------------------------------------
# RNA plot
#-----------------------------------------------
diablo_rna_feature_plot <- ggplot(plot_df, aes(x = logFC, y = Symbol, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "plasma",
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
    title = "Differential Expression by Gene"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#-----------------------------------------------
# Lipid plot
#-----------------------------------------------
diablo_lipid_feature_plot <- ggplot(lipid_plot_df, aes(x = Log2_FC_OS_vs_YS, y = Lipid, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "plasma",
    direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = "Lipid Species",
    title = "Differential Abundance by Lipid"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#-----------------------------------------------
# Metabolite plot
#-----------------------------------------------
diablo_metabolite_feature_plot <- ggplot(metabolite_plot_df, aes(x = log2FC, y = Metabolite, fill = neg_logFDR)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "plasma",
    direction = -1,
    limits = common_limits,
    name = expression(-log[10](FDR))
  ) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  theme_classic(base_size = 12) +
  labs(
    x = "log2 Fold Change (Old Sedentary vs Young Sedentary)",
    y = "Metabolite Species",
    title = "Differential Abundance by Metabolite"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#--------------------------------------------------#
print(diablo_rna_feature_plot)
print(diablo_lipid_feature_plot)
print(diablo_metabolite_feature_plot)

pdf("diablo_rna_feature_plot.pdf", width = 6, height = 5)
print(diablo_rna_feature_plot)
dev.off()

pdf("diablo_lipid_feature_plot.pdf", width = 6, height = 2)
print(diablo_lipid_feature_plot)
dev.off()

pdf("diablo_metabolite_feature_plot.pdf", width = 6, height = 2)
print(diablo_metabolite_feature_plot)
dev.off()










#==============================================================================#
#=======================================================#
#=======================================================#
#=======================================================#
#-----Don't think much of this worked ------------------#

#======================================================#
#----Balance Error Rate--------------------------------#
#======================================================#
#===========================================================#
#---- Step 1: Cross-validation of DIABLO Aging Axis --------#
#===========================================================#

# Train only on YNG vs OLD
diablo_train <- block.splsda(
  X_train, Y_train,
  ncomp = 1,
  keepX = list(RNA = 50, Metabolite = 20, Lipid = 20),
  design = design
)

# Cross-validation
set.seed(123)
perf_res <- perf(
  diablo_train,
  validation = "Mfold",
  folds = 5,
  nrepeat = 50,
  progressBar = TRUE
)

# Extract BER for comp1
ber_comp1 <- perf_res$error.rate$BER[[1]][,"overall"]
cat("Mean BER:", mean(ber_comp1, na.rm = TRUE), "\n")

# Inspect structure
str(perf_res$error.rate$BER)

# Peek at first few entries
perf_res$error.rate$BER[[1]]

names(perf_res)
perf_res$error.rate
perf_res$error.rate.class
perf_res$features$stable
perf_res$choice.ncomp

plot(perf_res)  



# Plot with 3 colors for 3 blocks
plot(perf_res, col = c("darkred","steelblue","darkgreen"), sd = TRUE)

# Projection of all samples (downstream analysis)
proj <- predict(diablo_train, newdata = X_all)











#--------------------------------------------------#
#-----Global Heatmap of 90 DIABLO Features---------#
#--------------------------------------------------#

#--------------------------------------------------#
#-----Global Heatmap of 90 DIABLO Features---------#
#--------------------------------------------------#

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
})

#---------------------------------------------------------
# 1) Extract DIABLO loadings (Comp1) for each block
#---------------------------------------------------------
loadings_rna <- diablo_res$loadings$RNA[,1]
loadings_met <- diablo_res$loadings$Metabolite[,1]
loadings_lip <- diablo_res$loadings$Lipid[,1]

# Combine all loadings into a single named vector
all_loadings <- c(loadings_rna, loadings_met, loadings_lip)

#---------------------------------------------------------
# 2) Subset combined matrix to the 90 selected features
#    (combined_mat = cbind(RNA_all, Metabo_all, Lipid_all))
#---------------------------------------------------------
rna_features <- rna_sel$RNA$name
met_features <- met_sel$Metabolite$name
lip_features <- lip_sel$Lipid$name
selected_feats <- c(rna_features, met_features, lip_features)

common_feats <- intersect(selected_feats, colnames(combined_mat))
heat_mat <- combined_mat[, common_feats, drop = FALSE]

#---------------------------------------------------------
# 3) Z-score per feature (row-wise)
#    Ensures each feature spans blue↔red internally
#---------------------------------------------------------
heat_mat_z <- t(scale(t(heat_mat)))

#---------------------------------------------------------
# 4) Apply DIABLO weights (after z-scoring)
#---------------------------------------------------------
weights <- all_loadings[colnames(heat_mat_z)]
heat_mat_weighted <- sweep(heat_mat_z, 2, weights, FUN = "*")

#---------------------------------------------------------
# 5) Annotation for groups (fixed order, no clustering of cols)
#---------------------------------------------------------
group_order <- c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap")
annot <- data.frame(Group = factor(proj_df$Group, levels = group_order))
rownames(annot) <- proj_df$Sample

#---------------------------------------------------------
# 6) Plot weighted heatmap
#    Rows = features, Cols = samples (ordered by Group)
#---------------------------------------------------------
pheatmap(
  t(heat_mat_weighted),
  annotation_col = annot,
  cluster_cols   = FALSE,   # keep group order fixed
  cluster_rows   = TRUE,    # cluster features
  show_rownames  = FALSE,
  fontsize_col   = 8,
  main = "Weighted Heatmap of 90 DIABLO Aging Axis Features"
)

#-----------------------------------------------------------------#

install.packages("ggraph")

#---------------------------------------------------------
#--------DIABLO Aging Axis (90) Correlation Network-------#
#---------------------------------------------------------
#---------------------------------------------------------
# DIABLO Aging Axis (90) Correlation Network
#---------------------------------------------------------
suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)
})

# 1) Build combined matrix (already available)
combined_mat <- cbind(RNA_all, Metabo_all, Lipid_all)
selected_feats <- c(rna_sel$RNA$name, met_sel$Metabolite$name, lip_sel$Lipid$name)
heat_mat <- combined_mat[, selected_feats, drop = FALSE]

# 2) Correlation matrix
cor_mat <- cor(heat_mat, method = "pearson")
thresh <- 0.7
cor_mat[abs(cor_mat) < thresh] <- 0

# 3) Use absolute correlations (positive weights only)
cor_mat_abs <- abs(cor_mat)

# 4) Convert to igraph object
g <- graph_from_adjacency_matrix(cor_mat_abs, mode = "undirected", weighted = TRUE, diag = FALSE)

# 5) Annotate nodes by data type
node_type <- case_when(
  V(g)$name %in% rna_sel$RNA$name ~ "RNA",
  V(g)$name %in% met_sel$Metabolite$name ~ "Metabolite",
  V(g)$name %in% lip_sel$Lipid$name ~ "Lipid",
  TRUE ~ "Other"
)
V(g)$type <- node_type

# 6) Publication-quality network plot
set.seed(123) # reproducibility
ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), color = "grey50") +
  geom_node_point(aes(color = V(g)$type), size = 5) +
  geom_node_text(aes(label = V(g)$name, color = V(g)$type), repel = TRUE, size = 3) +
  scale_edge_width(range = c(0.2, 2)) +
  scale_edge_alpha(range = c(0.2, 0.8)) +
  scale_color_manual(values = c(RNA = "#1f77b4",
                                Metabolite = "#ff7f0e",
                                Lipid = "#2ca02c")) +
  theme_void(base_size = 14) +
  theme(legend.position = "bottom") +
  guides(edge_width = "none", edge_alpha = "none") +
  labs(title = "Correlation Network of DIABLO Aging Axis Features",
       color = "Data Type")
#-------------------------------------------------------------#

#--------------------------------------------------#
#--- Heatmap with per-feature forced scaling -------#
#--------------------------------------------------#

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
})

# 1) Extract DIABLO loadings (Comp1) as before
loadings_rna <- diablo_res$loadings$RNA[,1]
loadings_met <- diablo_res$loadings$Metabolite[,1]
loadings_lip <- diablo_res$loadings$Lipid[,1]
all_loadings <- c(loadings_rna, loadings_met, loadings_lip)

# 2) Subset to the 90 selected features
rna_features <- rna_sel$RNA$name
met_features <- met_sel$Metabolite$name
lip_features <- lip_sel$Lipid$name
selected_feats <- c(rna_features, met_features, lip_features)
common_feats <- intersect(selected_feats, colnames(combined_mat))
heat_mat <- combined_mat[, common_feats, drop = FALSE]

# 3) Column-wise (feature-wise) scaling to [-1, +1]
col_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0, length(x)))  # flat feature → all zeros
  } else {
    return(2 * (x - rng[1]) / (rng[2] - rng[1]) - 1)
  }
}
heat_mat_scaled <- apply(heat_mat, 2, col_scale)

# 4) Apply DIABLO weights (sign/direction still matters)
weights <- all_loadings[colnames(heat_mat_scaled)]
heat_mat_weighted <- sweep(heat_mat_scaled, 2, weights, FUN = "*")

# 5) Annotation for groups
group_order <- c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap")
annot <- data.frame(Group = factor(proj_df$Group, levels = group_order))
rownames(annot) <- proj_df$Sample

# 6) Plot heatmap
pheatmap(
  t(heat_mat_weighted),        # rows = features, cols = samples
  annotation_col = annot,
  cluster_cols   = FALSE,      # keep group order
  cluster_rows   = TRUE,       # cluster features
  show_rownames  = FALSE,
  fontsize_col   = 8,
  main = "DIABLO Aging Axis Features (forced row min=-1, max=+1)"
)
#----------------------------------------------------------#

#==========================================================#
#------Accentuated Heatmap of DIABLO Aging Axis------------#
#==========================================================#

#--------------------------------------------------#
#--- Publication-Ready Heatmap of DIABLO Features--#
#--------------------------------------------------#

library(colorspace)

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
  library(colorspace)   # for better palettes
})

# 1) Extract DIABLO loadings (Comp1) as before
loadings_rna <- diablo_res$loadings$RNA[,1]
loadings_met <- diablo_res$loadings$Metabolite[,1]
loadings_lip <- diablo_res$loadings$Lipid[,1]
all_loadings <- c(loadings_rna, loadings_met, loadings_lip)

# 2) Subset to the 90 selected features
rna_features <- rna_sel$RNA$name
met_features <- met_sel$Metabolite$name
lip_features <- lip_sel$Lipid$name
selected_feats <- c(rna_features, met_features, lip_features)
common_feats <- intersect(selected_feats, colnames(combined_mat))
heat_mat <- combined_mat[, common_feats, drop = FALSE]

# 3) Row-wise scaling to [-1, +1]
col_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0, length(x)))  # flat feature → all zeros
  } else {
    return(2 * (x - rng[1]) / (rng[2] - rng[1]) - 1)
  }
}
heat_mat_scaled <- apply(heat_mat, 2, col_scale)

# 4) Apply DIABLO weights
weights <- all_loadings[colnames(heat_mat_scaled)]
heat_mat_weighted <- sweep(heat_mat_scaled, 2, weights, FUN = "*")

# 5) Optionally select top-N most variable features (set top_n = 90 for all)
top_n <- 40   # << adjust here (e.g., 30 or 90)
var_order <- order(apply(heat_mat_weighted, 2, var), decreasing = TRUE)
top_feats <- var_order[1:top_n]
heat_mat_final <- heat_mat_weighted[, top_feats, drop = FALSE]

# 6) Annotation for sample groups
group_order <- c("Young_Sed","Old_SedVeh","Old_PwrVeh","Old_PwrIRap","Old_PwrFRap")
annot <- data.frame(Group = factor(proj_df$Group, levels = group_order))
rownames(annot) <- proj_df$Sample

# 7) Row annotation (RNA vs Metabolite vs Lipid)
feat_types <- ifelse(colnames(heat_mat_final) %in% rna_features, "RNA",
                     ifelse(colnames(heat_mat_final) %in% met_features, "Metabolite", "Lipid"))
row_annot <- data.frame(DataType = feat_types)
rownames(row_annot) <- colnames(heat_mat_final)

# 8) Color palette and fixed breaks
# Stronger diverging palette (deep blue → white → deep red)
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)  
breaks <- seq(-.4, .4, length.out = 101)

# 9) Plot heatmap
pheatmap(
  t(heat_mat_final),        # rows = features, cols = samples
  annotation_col = annot,
  annotation_row = row_annot,
  cluster_cols   = FALSE,   # keep group order
  cluster_rows   = TRUE,    # cluster features
  show_rownames  = FALSE,
  fontsize_col   = 8,
  color = colors,
  breaks = breaks,
  main = paste("DIABLO Aging Axis Features (Top", top_n, ")")
)
#----------------------------------------------------------------#

#--------------------------------------------------#
#--- Heatmap of DIABLO Features (Old groups only) --#
#--------------------------------------------------#

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
})

# 1) Extract DIABLO loadings (Comp1) as before
loadings_rna <- diablo_res$loadings$RNA[,1]
loadings_met <- diablo_res$loadings$Metabolite[,1]
loadings_lip <- diablo_res$loadings$Lipid[,1]
all_loadings <- c(loadings_rna, loadings_met, loadings_lip)

# 2) Subset to the 90 selected features
rna_features <- rna_sel$RNA$name
met_features <- met_sel$Metabolite$name
lip_features <- lip_sel$Lipid$name
selected_feats <- c(rna_features, met_features, lip_features)
common_feats <- intersect(selected_feats, colnames(combined_mat))
heat_mat <- combined_mat[, common_feats, drop = FALSE]

# 3) Column-wise scaling to [-1, +1]
col_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0, length(x)))  # flat feature → all zeros
  } else {
    return(2 * (x - rng[1]) / (rng[2] - rng[1]) - 1)
  }
}
heat_mat_scaled <- apply(heat_mat, 2, col_scale)

# 4) Apply DIABLO weights
weights <- all_loadings[colnames(heat_mat_scaled)]
heat_mat_weighted <- sweep(heat_mat_scaled, 2, weights, FUN = "*")

# 5) Optionally select top-N most variable features
top_n <- 40   # adjust as needed
var_order <- order(apply(heat_mat_weighted, 2, var), decreasing = TRUE)
top_feats <- var_order[1:top_n]
heat_mat_final <- heat_mat_weighted[, top_feats, drop = FALSE]

# 6) Subset to OLD groups only (drop Young_Sed)
old_samples <- proj_df %>%
  filter(Group != "Young_Sed") %>%
  pull(Sample)
heat_mat_final <- heat_mat_final[old_samples, , drop = FALSE]

annot <- data.frame(Group = droplevels(proj_df$Group[match(old_samples, proj_df$Sample)]))
rownames(annot) <- old_samples

# 7) Row annotation (RNA vs Metabolite vs Lipid)
feat_types <- ifelse(colnames(heat_mat_final) %in% rna_features, "RNA",
                     ifelse(colnames(heat_mat_final) %in% met_features, "Metabolite", "Lipid"))
row_annot <- data.frame(DataType = feat_types)
rownames(row_annot) <- colnames(heat_mat_final)

# 8) Color palette and fixed breaks
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
breaks <- seq(-.4, .4, length.out = 101)

# 9) Plot heatmap
pheatmap(
  t(heat_mat_final),        # rows = features, cols = samples
  annotation_col = annot,
  annotation_row = row_annot,
  cluster_cols   = FALSE,   # keep groups
  cluster_rows   = TRUE,    # cluster features
  show_rownames  = FALSE,
  fontsize_col   = 8,
  color = colors,
  breaks = breaks,
  main = paste("DIABLO Aging Axis Features (Old groups only, Top", top_n, ")")
)
