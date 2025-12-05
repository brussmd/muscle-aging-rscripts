#Liver DIABLO Analysis

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(tidyr)
  library(edgeR);  library(AnnotationDbi); library(org.Mm.eg.db)
})

#================================================================#
#-------Define Liver Aging Axis----------------------------------#
#================================================================#

#----------------------------------------------------------------#
#=======Set-up old vs yng Liver RNA Seq COUNTS====================#

# ---- Load counts ----
liver_count_os_vs_ys  <- read_xlsx("20250922_M007988_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")

# 1) Rename the first column to "Ensembl"
rename_first_col <- function(df) {
  cn <- colnames(df)
  cn[1] <- "Ensembl"
  colnames(df) <- cn
  df
}
liver_count_os_vs_ys <- rename_first_col(liver_count_os_vs_ys)

head(liver_count_os_vs_ys)

# ---- Define sample groups ----
liver_yng_sed_samples <- c("L05","L17","L26","L32","L35","L38")
liver_old_sed_samples <- c("L06","L09","L13","L19","L21","L25","L29","L39")
liver_all_samples     <- c(liver_yng_sed_samples, liver_old_sed_samples)

# ---- Keep only those columns ----
liver_merged_counts <- liver_count_os_vs_ys %>%
  dplyr::select(Ensembl, dplyr::all_of(liver_all_samples))

# ---- Map Ensembl -> ENTREZ ----
liver_merged_counts <- liver_merged_counts %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db,
                           keys = Ensembl_noDec,
                           keytype = "ENSEMBL",
                           column = "ENTREZID",
                           multiVals = "first")) %>%
  tidyr::drop_na(ENTREZID)
head(liver_merged_counts)

# ---- Counts matrix ----
liver_counts_matrix <- liver_merged_counts %>%
  dplyr::select(dplyr::all_of(liver_all_samples)) %>%
  replace(is.na(.), 0) %>%
  as.matrix()
rownames(liver_counts_matrix) <- liver_merged_counts$ENTREZID

head(liver_counts_matrix)

# ---- Define groups ----
liver_group <- factor(c(rep("YNG_SED", length(liver_yng_sed_samples)),
                  rep("OLD_SED", length(liver_old_sed_samples))),
                levels = c("YNG_SED","OLD_SED"))

stopifnot(length(liver_group) == ncol(liver_counts_matrix))

# ---- edgeR: TMM normalization, filtering, logCPM ----
liver_dge <- DGEList(counts = liver_counts_matrix, group = liver_group)
liver_dge <- calcNormFactors(liver_dge, method = "TMM")
keep <- filterByExpr(liver_dge)                   
liver_dge <- liver_dge[keep, , keep.lib.sizes = FALSE]
liver_logCPM_matrix <- cpm(liver_dge, log = TRUE, prior.count = 1)

cat("Genes kept:", nrow(liver_logCPM_matrix), " Samples:", ncol(liver_logCPM_matrix), "\n")

# ---- Map ENTREZ -> SYMBOL and collapse duplicates ----
liver_symbols <- mapIds(org.Mm.eg.db,
                  keys    = rownames(liver_logCPM_matrix),
                  keytype = "ENTREZID",
                  column  = "SYMBOL",
                  multiVals = "first")

liver_sym_df <- as.data.frame(liver_logCPM_matrix) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  mutate(Symbol = liver_symbols) %>%
  tidyr::drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

liver_logCPM_symbol <- liver_sym_df %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()

# ---- Keep exact sample order ----
liver_logCPM_symbol <- liver_logCPM_symbol[, liver_all_samples, drop = FALSE]

# ---- Sanity checks ----
stopifnot(!anyNA(liver_logCPM_symbol))
stopifnot(identical(colnames(liver_logCPM_symbol), liver_all_samples))

saveRDS(liver_logCPM_symbol, "liver_logCPM_symbol_yngOLDsed.RDS")
dim(liver_logCPM_symbol)
head(liver_logCPM_symbol)
#---------------------------------------------------#


#----------------------------------------------------------------#
#=======Set-up old vs yng Liver Metabo Log ======================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
})

# ---- Load data ----
pwr_rapa_liver_metabo_data <- readr::read_csv("Konopka_Liver_HILIC.csv", show_col_types = FALSE)
head(pwr_rapa_liver_metabo_data)

# ---- Drop Group column, keep Sample + metabolites ----
liver_metabo_wide <- pwr_rapa_liver_metabo_data %>% dplyr::select(-Group)

# ---- Convert long → wide: rows = metabolites, cols = samples ----
liver_metabo_matrix <- liver_metabo_wide %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()

rownames(liver_metabo_matrix) <- liver_metabo_matrix$Metabolite
liver_metabo_matrix$Metabolite <- NULL
head(liver_metabo_matrix)

# ---- Map metabolomics sample names (YS1/OS1/…) → T_* tube codes ----
liver_sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "L06","OS_1","L09","OS_2","L13","OS_3","L19","OS_4","L21","OS_5","L25","OS_6","L29","OS_7","L39","OS_8",
  "L05","YS_1","L15","YS_2","L17","YS_3","L23","YS_4","L26","YS_5","L32","YS_6","L35","YS_7","L38","YS_8"
)

liver_name_map <- setNames(liver_sample_key$TubeCode, liver_sample_key$MetabolomicsName)
liver_new_names <- liver_name_map[colnames(liver_metabo_matrix)]
colnames(liver_metabo_matrix) <- ifelse(is.na(liver_new_names), colnames(liver_metabo_matrix), liver_new_names)

# ---- Define sample groups ----
liver_yng_sed_samples <- c("L05","L17","L26","L32","L35","L38")
liver_old_sed_samples <- c("L06","L09","L13","L19","L21","L25","L29","L39")
liver_all_samples     <- c(liver_yng_sed_samples, liver_old_sed_samples)

liver_present_metabo <- intersect(liver_all_samples, colnames(liver_metabo_matrix))
missing_for_liver_metabo <- setdiff(liver_all_samples, colnames(liver_metabo_matrix))
if (length(missing_for_liver_metabo)) message("Metabolomics missing: ", paste(missing_for_liver_metabo, collapse = ", "))

liver_metabo_matrix_subset <- liver_metabo_matrix[, liver_present_metabo, drop = FALSE]
head(liver_metabo_matrix_subset)

# ---- Ensure numeric ----
liver_metabo_matrix_subset[] <- lapply(liver_metabo_matrix_subset, function(x) as.numeric(as.character(x)))

# ---- Impute zeros/NA with half min positive per metabolite ----
liver_metabo_imputed <- t(apply(liver_metabo_matrix_subset, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))

# ---- Log2 transform ----
liver_metabo_log <- as.data.frame(log2(liver_metabo_imputed))
colnames(liver_metabo_log) <- liver_present_metabo
rownames(liver_metabo_log) <- rownames(liver_metabo_matrix_subset)

# ---- Align with RNA samples if liver_logCPM_symbol exists ----
if (exists("liver_logCPM_symbol")) {
  liver_common_samples <- intersect(liver_all_samples, colnames(liver_logCPM_symbol))
  liver_metabo_log     <- liver_metabo_log[, liver_common_samples, drop = FALSE]
  liver_logCPM_symbol  <- liver_logCPM_symbol[, liver_common_samples, drop = FALSE]
  stopifnot(identical(colnames(liver_metabo_log), colnames(liver_logCPM_symbol)))
}

# ---- Save ----
saveRDS(liver_metabo_log, "liver_metabo_log_yngOLDsed.RDS")

# ---- Quick peek ----
dim(liver_metabo_log)
head(liver_metabo_log)
#-------------------------------------------------#

#============================================================#
#---------Build liver_lipid_matrix old and young------------------#
#============================================================#

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble)
})

# ---- Load lipidomics ----
liver_lipid_df <- readr::read_csv("Konopka_liver_lipids_CORRECTED.csv", show_col_types = FALSE)
head(liver_lipid_df)

# Correct identification of lipid columns
liver_lipid_cols <- setdiff(names(liver_lipid_df), c("Sample", "Group", "Tube code"))

# Re-run the imputation + log2 transform
liver_lipid_df[ , liver_lipid_cols] <- lapply(liver_lipid_df[ , liver_lipid_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})

# Build matrix again
liver_lipid_matrix <- t(as.matrix(liver_lipid_df[ , liver_lipid_cols]))
colnames(liver_lipid_matrix) <- liver_lipid_df$Sample

head(liver_lipid_matrix)

# ---- Map sample names (e.g. OSED1, YSED1) → L* codes ----
liver_lipid_sample_key <- tibble::tribble(
  ~TubeCode, ~LipidName,
  "L06","OSED1","L09","OSED2","L13","OSED3","L19","OSED4","L21","OSED5","L25","OSED6","L29","OSED7","L39","OSED8",
  "L05","YSED1","L15","YSED2","L17","YSED3","L23","YSED4","L26","YSED5","L32","YSED6","L35","YSED7","L38","YSED8"
)

liver_lipid_map <- setNames(liver_lipid_sample_key$TubeCode, liver_lipid_sample_key$LipidName)
liver_lipid_mapped <- liver_lipid_map[colnames(liver_lipid_matrix)]
colnames(liver_lipid_matrix) <- ifelse(is.na(liver_lipid_mapped), colnames(liver_lipid_matrix), liver_lipid_mapped)

# ---- Define sample groups ----
liver_yng_sed_samples <- c("L05","L17","L26","L32","L35","L38")
liver_old_sed_samples <- c("L06","L09","L13","L19","L21","L25","L29","L39")
liver_all_samples     <- c(liver_yng_sed_samples, liver_old_sed_samples)

liver_lipid_present <- intersect(liver_all_samples, colnames(liver_lipid_matrix))
liver_lipid_log <- liver_lipid_matrix[ , liver_lipid_present, drop = FALSE]

# ---- Align with RNA samples if logCPM_symbol exists ----
if (exists("liver_logCPM_symbol")) {
  liver_lipid_common <- intersect(liver_all_samples, colnames(liver_logCPM_symbol))
  liver_lipid_log     <- liver_lipid_log[, liver_lipid_common, drop = FALSE]
  liver_logCPM_symbol <- liver_logCPM_symbol[, liver_lipid_common, drop = FALSE]
  stopifnot(identical(colnames(liver_lipid_log), colnames(liver_logCPM_symbol)))
}

# ---- Save ----
saveRDS(liver_lipid_log, "liver_lipid_log_yngOLDsed.RDS")

# ---- Sanity check ----
dim(liver_lipid_log)
head(liver_lipid_log)
#------------------------------------------------------------#

#============================================================#
#---------Build DIABLO matrix old and young------------------#
#============================================================#
suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
})

#--------------------------------------------------------------#
#---- Align common samples across omics -----------------------#
#--------------------------------------------------------------#
liver_common_samples <- Reduce(intersect, list(
  colnames(liver_logCPM_symbol),
  colnames(liver_metabo_log),
  colnames(liver_lipid_log)
))

liver_logCPM_symbol <- liver_logCPM_symbol[, liver_common_samples, drop = FALSE]
liver_metabo_log    <- liver_metabo_log[, liver_common_samples, drop = FALSE]
liver_lipid_log     <- liver_lipid_log[, liver_common_samples, drop = FALSE]

stopifnot(
  identical(colnames(liver_logCPM_symbol), colnames(liver_metabo_log)),
  identical(colnames(liver_logCPM_symbol), colnames(liver_lipid_log))
)

#--------------------------------------------------------------#
#---- Define sample groups (OLD vs YNG) -----------------------#
#--------------------------------------------------------------#
liver_yng_sed_samples <- c("L05","L17","L26","L32","L35","L38")
liver_old_sed_samples <- c("L06","L09","L13","L19","L21","L25","L29","L39")

Y <- factor(ifelse(colnames(liver_logCPM_symbol) %in% liver_old_sed_samples, "OLD", "YNG"),
            levels = c("YNG","OLD"))
names(Y) <- colnames(liver_logCPM_symbol)

#--------------------------------------------------------------#
#---- Build omics list (samples in rows, features in columns)---#
#--------------------------------------------------------------#
X <- list(
  RNA        = t(liver_logCPM_symbol),
  Metabolite = t(liver_metabo_log),
  Lipid      = t(liver_lipid_log)
)

#--------------------------------------------------------------#
#---- Design matrix (full integration) ------------------------#
#--------------------------------------------------------------#
design <- matrix(1, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

#--------------------------------------------------------------#
#---- Run DIABLO: 1 component, feature selection --------------#
#--------------------------------------------------------------#
diablo_res <- block.splsda(
  X, Y, ncomp = 1,
  keepX = list(
    RNA = 50,
    Metabolite = 20,
    Lipid = 20
  ),
  design = design
)

#--------------------------------------------------------------#
#---- Extract selected features per omic block ----------------#
#--------------------------------------------------------------#
extract_features <- function(sel, block) {
  data.frame(
    block   = block,
    feature = sel[[block]]$name,
    loading = sel[[block]]$value$value.var,
    stringsAsFactors = FALSE
  )
}

rna_sel <- selectVar(diablo_res, block = "RNA", comp = 1)
met_sel <- selectVar(diablo_res, block = "Metabolite", comp = 1)
lip_sel <- selectVar(diablo_res, block = "Lipid", comp = 1)

selected_df <- bind_rows(
  extract_features(rna_sel, "RNA"),
  extract_features(met_sel, "Metabolite"),
  extract_features(lip_sel, "Lipid")
) %>%
  arrange(block, desc(abs(loading)))

#--------------------------------------------------------------#
#---- Save and preview ----------------------------------------#
#--------------------------------------------------------------#
write.csv(selected_df, "Liver_DIABLO_selected_features_comp1.csv", row.names = FALSE)
cat("Top features defining the liver aging axis:\n")
print(head(selected_df, 20))
selected_df

#--------------------------------------------------------------#
#---- Optional diagnostic plots -------------------------------#
#--------------------------------------------------------------#
plotDiablo(diablo_res, legend = TRUE, title = "Liver: DIABLO Aging Axis (Comp 1)")
circosPlot(diablo_res, cutoff = 0.7)
network(diablo_res, comp = 1, cutoff = 0.6)
plotVar(diablo_res, comp = 1, var.names = TRUE)
#-------------------------------------------------------------#

#==============================================================#
#------DIABLO Intervention Analysis----------------------------#
#==============================================================#

#/////////////////////////////////////////////////////////#
#---- Liver Transcriptomics (full integrated dataset) -----#
#/////////////////////////////////////////////////////////#

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(tidyr)
  library(tibble)
})

#----------------------------------------------------------#
# 1. Load liver comparison count files
#----------------------------------------------------------#

liver_oldsedveh_vs_yngsedveh_counts  <- read_xlsx("20250922_M007988_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
liver_oldpwrveh_vs_yngsedveh_counts  <- read_xlsx("20250922_M007988_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx")
liver_oldpwrirap_vs_yngsedveh_counts <- read_xlsx("20250922_M007988_Set01_edgeRglm_Counts_OLD_PWR_IRAP-YNG_SED_VEH.xlsx")
liver_oldpwrfrap_vs_yngsedveh_counts <- read_xlsx("20250922_M007988_Set01_edgeRglm_Counts_OLD_PWR_FRAP-YNG_SED_VEH.xlsx")

#----------------------------------------------------------#
# 2. Rename first column to "Ensembl"
#----------------------------------------------------------#

for (dfname in c("liver_oldsedveh_vs_yngsedveh_counts",
                 "liver_oldpwrveh_vs_yngsedveh_counts",
                 "liver_oldpwrirap_vs_yngsedveh_counts",
                 "liver_oldpwrfrap_vs_yngsedveh_counts")) {
  tmp <- get(dfname)
  colnames(tmp)[1] <- "Ensembl"
  assign(dfname, tmp)
}

#----------------------------------------------------------#
# 3. Define group sample IDs
#----------------------------------------------------------#

liver_yng_sed_samples     <- c("L05","L17","L26","L32","L35","L38")
liver_old_sed_samples     <- c("L06","L09","L13","L19","L21","L25","L29","L39")
liver_old_pwr_samples     <- c("L08","L10","L11","L16","L41","L43","L45","L46")
liver_old_pwrirap_samples <- c("L03","L14","L31","L33","L37","L47","L48")
liver_old_pwrfrap_samples <- c("L01","L04","L22","L24","L28","L34","L40","L44")

#----------------------------------------------------------#
# 4. Merge all count tables, keeping one copy of YNG_SED
#----------------------------------------------------------#

liver_merged_counts <- liver_oldsedveh_vs_yngsedveh_counts %>%
  dplyr::select(Ensembl, all_of(liver_yng_sed_samples), all_of(liver_old_sed_samples)) %>%
  left_join(liver_oldpwrveh_vs_yngsedveh_counts %>%
              dplyr::select(Ensembl, all_of(setdiff(liver_old_pwr_samples, liver_yng_sed_samples))),
            by = "Ensembl") %>%
  left_join(liver_oldpwrirap_vs_yngsedveh_counts %>%
              dplyr::select(Ensembl, all_of(setdiff(liver_old_pwrirap_samples, liver_yng_sed_samples))),
            by = "Ensembl") %>%
  left_join(liver_oldpwrfrap_vs_yngsedveh_counts %>%
              dplyr::select(Ensembl, all_of(setdiff(liver_old_pwrfrap_samples, liver_yng_sed_samples))),
            by = "Ensembl")

#----------------------------------------------------------#
# 5. Map Ensembl → ENTREZ IDs
#----------------------------------------------------------#

liver_merged_counts <- liver_merged_counts %>%
  mutate(
    Ensembl_noDec = str_remove(Ensembl, "\\..+"),
    ENTREZID = mapIds(org.Mm.eg.db,
                      keys = Ensembl_noDec,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")
  ) %>%
  tidyr::drop_na(ENTREZID)

#----------------------------------------------------------#
# 6. Build count matrix and group factors
#----------------------------------------------------------#

liver_counts_matrix <- liver_merged_counts %>%
  dplyr::select(-Ensembl, -Ensembl_noDec, -ENTREZID) %>%
  replace(is.na(.), 0) %>%
  as.matrix()
rownames(liver_counts_matrix) <- liver_merged_counts$ENTREZID

# Define sample groups
liver_group <- factor(c(
  rep("YNG_SED",      length(liver_yng_sed_samples)),
  rep("OLD_SED",      length(liver_old_sed_samples)),
  rep("OLD_PWR",      length(liver_old_pwr_samples)),
  rep("OLD_PWR_IRAP", length(liver_old_pwrirap_samples)),
  rep("OLD_PWR_FRAP", length(liver_old_pwrfrap_samples))
))

stopifnot(ncol(liver_counts_matrix) == length(liver_group))

#----------------------------------------------------------#
# 7. Create DGEList, normalize, and filter
#----------------------------------------------------------#

liver_dge <- DGEList(counts = liver_counts_matrix, group = liver_group)
liver_dge <- calcNormFactors(liver_dge, method = "TMM")

liver_keep <- filterByExpr(liver_dge)
liver_dge <- liver_dge[liver_keep, , keep.lib.sizes = FALSE]
cat("Genes kept after filtering:", sum(liver_keep), "\n")

#----------------------------------------------------------#
# 8. Compute logCPM
#----------------------------------------------------------#

liver_logCPM_matrix <- cpm(liver_dge, log = TRUE, prior.count = 1)

cat("Final matrix dimensions: ",
    nrow(liver_logCPM_matrix), " genes × ",
    ncol(liver_logCPM_matrix), " samples\n")

#----------------------------------------------------------#
# 9. Map ENTREZ → SYMBOL and collapse duplicates
#----------------------------------------------------------#

liver_symbols <- mapIds(org.Mm.eg.db,
                        keys = rownames(liver_logCPM_matrix),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")

liver_logCPM_df <- as.data.frame(liver_logCPM_matrix) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  mutate(Symbol = liver_symbols) %>%
  drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

liver_logCPM_symbol_full <- liver_logCPM_df %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

#----------------------------------------------------------#
# 10. Save and sanity checks
#----------------------------------------------------------#

saveRDS(liver_logCPM_symbol_full, "liver_logCPM_symbol_full.RDS")
head(liver_logCPM_symbol_full)

cat("Genes:", nrow(liver_logCPM_symbol_full), 
    " Samples:", ncol(liver_logCPM_symbol_full), "\n")

anyNA(liver_logCPM_symbol_full)        # should be FALSE
anyDuplicated(rownames(liver_logCPM_symbol_full)) # should be 0
#---------------------------------------------------------------#

#/////////////////////////////////////////////////////////#
#---- Subset liver RNA-seq data to 50 DIABLO features -----#
#/////////////////////////////////////////////////////////#

# 1. Confirm available data
dim(liver_logCPM_symbol_full)
head(rownames(liver_logCPM_symbol_full))[1:10]

# 2. Confirm feature list
length(liver_rna_features)
head(liver_rna_features, 10)

# 3. Identify overlapping genes
common_rna_feats   <- intersect(liver_rna_features, rownames(liver_logCPM_symbol_full))
missing_rna_feats  <- setdiff(liver_rna_features, rownames(liver_logCPM_symbol_full))

cat("Total DIABLO RNA features:", length(liver_rna_features), "\n")
cat("Found in liver dataset:", length(common_rna_feats), "\n")
cat("Missing from liver dataset:", length(missing_rna_feats), "\n")

if (length(missing_rna_feats) > 0) {
  cat("Missing genes:\n")
  print(missing_rna_feats)
}

# 4. Subset to those features (rows = genes, columns = samples)
liver_RNA_subset <- liver_logCPM_symbol_full[common_rna_feats, , drop = FALSE]

# 5. Transpose so that rows = samples, columns = genes
liver_RNA_subset_t <- t(liver_RNA_subset)

# 6. Sanity checks
cat("Samples:", nrow(liver_RNA_subset_t), " | Features:", ncol(liver_RNA_subset_t), "\n")
head(rownames(liver_RNA_subset_t))   # should be sample IDs like L05, L06, etc.
head(colnames(liver_RNA_subset_t))   # should be gene symbols

stopifnot(is.numeric(liver_RNA_subset_t[1,1]))

# 7. Save for downstream DIABLO projection
saveRDS(liver_RNA_subset_t, "liver_RNA_subset_50features.RDS")

cat("✅ Saved 50-feature RNA subset for DIABLO projection.\n")
#--------------------------------------------------------------#

#/////////////////////////////////////////////////////////#
#---- Liver Metabolomics (full integrated dataset) -----#
#/////////////////////////////////////////////////////////#
#/////////////////////////////////////////////////////////#
#---- Liver Metabolomics (HILIC) preprocessing -------------#
#/////////////////////////////////////////////////////////#

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

#----------------------------------------------------------#
# 1) Load raw metabolomics data
#----------------------------------------------------------#
liver_metabo_raw <- read_csv("Konopka_Liver_HILIC.csv", show_col_types = FALSE)

# 2) Keep Sample + feature columns (drop Group if present)
liver_metabo_wide <- liver_metabo_raw %>%
  dplyr::select(-any_of("Group"))

# 3) Reshape to matrix: rows = metabolites, cols = samples
liver_metabo_matrix <- liver_metabo_wide %>%
  pivot_longer(-Sample, names_to = "Metabolite", values_to = "Abundance") %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  as.data.frame()

rownames(liver_metabo_matrix) <- liver_metabo_matrix$Metabolite
liver_metabo_matrix$Metabolite <- NULL

#----------------------------------------------------------#
# 4) Build sample key (L codes and updated group prefixes)
#----------------------------------------------------------#
liver_sample_key <- tibble::tribble(
  ~TubeCode, ~MetabolomicsName,
  "L01","OFR_1","L04","OFR_2","L22","OFR_3","L24","OFR_4","L28","OFR_5","L34","OFR_6","L40","OFR_7","L44","OFR_8",
  "L02","OIR_1","L03","OIR_2","L14","OIR_3","L31","OIR_4","L33","OIR_5","L37","OIR_6","L47","OIR_7","L48","OIR_8",
  "L08","OV_1","L10","OV_2","L11","OV_3","L16","OV_4","L41","OV_5","L43","OV_6","L45","OV_7","L46","OV_8",
  "L06","OS_1","L09","OS_2","L13","OS_3","L19","OS_4","L21","OS_5","L25","OS_6","L29","OS_7","L39","OS_8",
  "L07","YV_1","L12","YV_2","L18","YV_3","L20","YV_4","L27","YV_5","L30","YV_6","L36","YV_7","L42","YV_8",
  "L05","YS_1","L15","YS_2","L17","YS_3","L23","YS_4","L26","YS_5","L32","YS_6","L35","YS_7","L38","YS_8"
)

# Map metabolomics sample names → TubeCode (L##)
liver_name_map <- setNames(liver_sample_key$TubeCode, liver_sample_key$MetabolomicsName)
new_liver_names <- liver_name_map[colnames(liver_metabo_matrix)]
colnames(liver_metabo_matrix) <- ifelse(is.na(new_liver_names),
                                        colnames(liver_metabo_matrix),
                                        new_liver_names)

#----------------------------------------------------------#
# 5) Convert to numeric, impute zeros/NA, log2-transform
#----------------------------------------------------------#
liver_metabo_matrix[] <- lapply(liver_metabo_matrix, function(x) as.numeric(as.character(x)))

liver_metabo_imputed <- t(apply(liver_metabo_matrix, 1, function(x) {
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  x
}))

liver_metabo_log_all <- as.data.frame(log2(liver_metabo_imputed))

#----------------------------------------------------------#
# 6) Subset to DIABLO-selected metabolites
#----------------------------------------------------------#
# Replace this with your actual DIABLO-selected metabolite list:
# e.g., liver_metabo_features <- c("3-Ureidopropionic acid", "Dihydroorotic Acid", ...)
# For now assume you've already loaded it into liver_metabo_features

common_liver_met_feats  <- intersect(liver_metabo_features, rownames(liver_metabo_log_all))
missing_liver_met_feats <- setdiff(liver_metabo_features, rownames(liver_metabo_log_all))

cat("DIABLO-selected metabolites:", length(liver_metabo_features), "\n")
cat("Overlap with full dataset:", length(common_liver_met_feats), "\n")
cat("Missing metabolites:", length(missing_liver_met_feats), "\n")

if (length(missing_liver_met_feats) > 0) {
  cat("Missing metabolites:\n")
  print(missing_liver_met_feats)
}

#----------------------------------------------------------#
# 7) Build final metabolomics matrix for DIABLO
#----------------------------------------------------------#
liver_Metabo_all <- t(liver_metabo_log_all[common_liver_met_feats, , drop = FALSE])

# Sanity check
cat("Samples:", nrow(liver_Metabo_all), " | Metabolite features:", ncol(liver_Metabo_all), "\n")
head(rownames(liver_Metabo_all))  # should be L## sample IDs
head(colnames(liver_Metabo_all))  # metabolite names
stopifnot(is.numeric(liver_Metabo_all[1,1]))

#----------------------------------------------------------#
# 8) Align with RNA samples (optional)
#----------------------------------------------------------#
common_samples <- intersect(rownames(liver_RNA_subset_t), rownames(liver_Metabo_all))
liver_Metabo_all <- liver_Metabo_all[common_samples, , drop = FALSE]
liver_Metabo_all <- liver_Metabo_all[rownames(liver_RNA_subset_t), , drop = FALSE]

stopifnot(identical(rownames(liver_RNA_subset_t), rownames(liver_Metabo_all)))

cat("Samples aligned between RNA and Metabolomics datasets:", nrow(liver_Metabo_all), "\n")

#----------------------------------------------------------#
# 9) Save
#----------------------------------------------------------#
saveRDS(liver_Metabo_all, "liver_Metabo_subset_DIABLO.RDS")
cat("✅ Saved liver metabolomics DIABLO subset.\n")
#-------------------------------------------------------------#

#/////////////////////////////////////////////////////////#
#---- Liver Lipidomics (Corrected sample names) -----------#
#/////////////////////////////////////////////////////////#

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

#----------------------------------------------------------#
# 1) Load lipidomics data
#----------------------------------------------------------#
liver_lipid_df <- read_csv("Konopka_liver_lipids_CORRECTED.csv", show_col_types = FALSE)
liver_lipid_df


# 2) Identify lipid columns (everything except metadata)
liver_lipid_cols <- setdiff(names(liver_lipid_df), c("Sample", "Group", "Tube code"))

# 3) Impute zeros/NA per lipid, then log2-transform
liver_lipid_df[ , liver_lipid_cols] <- lapply(liver_lipid_df[ , liver_lipid_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6
  x[x == 0] <- min_pos / 2
  log2(x)
})

# 4) Convert to matrix: rows = lipids, cols = samples
liver_lipid_matrix <- t(as.matrix(liver_lipid_df[ , liver_lipid_cols]))
colnames(liver_lipid_matrix) <- liver_lipid_df$Sample

#----------------------------------------------------------#
# 5) Build new sample key (with corrected group prefixes)
#----------------------------------------------------------#
liver_lipid_key <- tibble::tribble(
  ~TubeCode, ~LipidomicsName,
  "L01","OFR1","L04","OFR2","L22","OFR3","L24","OFR4","L28","OFR5","L34","OFR6","L40","OFR7","L44","OFR8",
  "L02","OIR1","L03","OIR2","L14","OIR3","L31","OIR4","L33","OIR5","L37","OIR6","L47","OIR7","L48","OIR8",
  "L08","OEX1","L10","OEX2","L11","OEX3","L16","OEX4","L41","OEX5","L43","OEX6","L45","OEX7","L46","OEX8",
  "L06","OSED1","L09","OSED2","L13","OSED3","L19","OSED4","L21","OSED5","L25","OSED6","L29","OSED7","L39","OSED8",
  "L07","YEX1","L12","YEX2","L18","YEX3","L20","YEX4","L27","YEX5","L30","YEX6","L36","YEX7","L42","YEX8",
  "L05","YSED1","L15","YSED2","L17","YSED3","L23","YSED4","L26","YSED5","L32","YSED6","L35","YSED7","L38","YSED8"
)

# Map lipidomics sample names → TubeCode (L##)
liver_lipid_map <- setNames(liver_lipid_key$TubeCode, liver_lipid_key$LipidomicsName)
mapped_lipid_names <- liver_lipid_map[colnames(liver_lipid_matrix)]
colnames(liver_lipid_matrix) <- ifelse(is.na(mapped_lipid_names),
                                       colnames(liver_lipid_matrix),
                                       mapped_lipid_names)

#----------------------------------------------------------#
# 6) Subset to DIABLO-selected lipids
#----------------------------------------------------------#
# e.g. liver_lipid_features <- c("TG(O-52:2)", "LPC(22:4) (a\\b)", ...)
# (replace with your real 20 lipid features)
common_liver_lip_feats  <- intersect(liver_lipid_features, rownames(liver_lipid_matrix))
missing_liver_lip_feats <- setdiff(liver_lipid_features, rownames(liver_lipid_matrix))

cat("DIABLO-selected lipids:", length(liver_lipid_features), "\n")
cat("Overlap with dataset:", length(common_liver_lip_feats), "\n")
cat("Missing:", length(missing_liver_lip_feats), "\n")

if (length(missing_liver_lip_feats) > 0) {
  cat("Missing lipids:\n")
  print(missing_liver_lip_feats)
}

# Final Lipid_all matrix: rows = samples, cols = 20 lipids
liver_Lipid_all <- t(liver_lipid_matrix[common_liver_lip_feats, , drop = FALSE])

#----------------------------------------------------------#
# 7) Align Lipid_all with RNA samples (for DIABLO projection)
#----------------------------------------------------------#
common_samples <- intersect(rownames(liver_RNA_subset_t), rownames(liver_Lipid_all))
liver_Lipid_all <- liver_Lipid_all[common_samples, , drop = FALSE]
liver_Lipid_all <- liver_Lipid_all[rownames(liver_RNA_subset_t), , drop = FALSE]

stopifnot(identical(rownames(liver_RNA_subset_t), rownames(liver_Lipid_all)))

cat("Samples in RNA subset:", nrow(liver_RNA_subset_t), "\n")
cat("Samples in Lipid subset after alignment:", nrow(liver_Lipid_all), "\n")

#----------------------------------------------------------#
# 8) Save
#----------------------------------------------------------#
head(liver_Lipid_all)
head(liver_Metabo_all)
head(liver_RNA_subset_t)

saveRDS(liver_Lipid_all, "liver_Lipid_subset_DIABLO.RDS")
cat("✅ Saved liver lipidomics DIABLO subset.\n")
#-------------------------------------------------------#

#=========================================================#
#---------Intervention Analysis---------------------------#
#=========================================================#

#============================================================#
# Liver DIABLO Aging Axis Projection (RNA + Metabo + Lipid)
#============================================================#

suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
  library(broom)
})

#-----------------------------------------------------------
# 0) Define sample groups (L## identifiers)
#-----------------------------------------------------------
liver_yng_sed_samples <- c("L05","L17","L26","L32","L35","L38")
liver_old_sed_samples <- c("L06","L09","L13","L19","L21","L25","L29","L39")
liver_old_pwr_samples <- c("L08","L10","L11","L16","L41","L43","L45","L46")
liver_old_pwrirap_samples <- c("L03","L14","L31","L33","L37","L47","L48")
liver_old_pwrfrap_samples <- c("L01","L04","L22","L24","L28","L34","L40","L44")

#-----------------------------------------------------------
# 1) Final alignment check of omics matrices
#-----------------------------------------------------------
cat("RNA:", dim(liver_RNA_subset_t)[1], "samples x", dim(liver_RNA_subset_t)[2], "features\n")
cat("Metabolite:", dim(liver_Metabo_all)[1], "samples x", dim(liver_Metabo_all)[2], "features\n")
cat("Lipid:", dim(liver_Lipid_all)[1], "samples x", dim(liver_Lipid_all)[2], "features\n")

stopifnot(identical(rownames(liver_RNA_subset_t), rownames(liver_Metabo_all)))
stopifnot(identical(rownames(liver_RNA_subset_t), rownames(liver_Lipid_all)))

#-----------------------------------------------------------
# 2) Build full omics list (rows = samples, cols = features)
#-----------------------------------------------------------
liver_X_all <- list(
  RNA        = liver_RNA_subset_t,
  Metabolite = liver_Metabo_all,
  Lipid      = liver_Lipid_all
)

#-----------------------------------------------------------
# 3) Define training set (young vs old sedentary vehicle)
#-----------------------------------------------------------
liver_Y_train <- factor(c(
  rep("YNG", length(liver_yng_sed_samples)),
  rep("OLD", length(liver_old_sed_samples))
), levels = c("YNG","OLD"))

liver_train_samples <- c(liver_yng_sed_samples, liver_old_sed_samples)
liver_X_train <- lapply(liver_X_all, function(m) m[liver_train_samples, , drop = FALSE])

#-----------------------------------------------------------
# 4) Run DIABLO on training samples (1 component)
#-----------------------------------------------------------
liver_design <- matrix(1, ncol = 3, nrow = 3,
                       dimnames = list(names(liver_X_all), names(liver_X_all)))
diag(liver_design) <- 0

liver_diablo <- block.splsda(
  liver_X_train,
  liver_Y_train,
  ncomp = 1,
  keepX = list(RNA = 50, Metabolite = 20, Lipid = 20),
  design = liver_design
)

#-----------------------------------------------------------
# 5) Project all samples into the DIABLO space
#-----------------------------------------------------------
liver_proj <- predict(liver_diablo, newdata = liver_X_all)
liver_scores <- lapply(liver_proj$variates, function(b) b[,1])
liver_scores_mat <- do.call(cbind, liver_scores)
liver_comp1 <- rowMeans(liver_scores_mat)

# Build projection dataframe
liver_proj_df <- data.frame(
  Sample = rownames(liver_scores_mat),
  Comp1  = liver_comp1
) %>%
  mutate(Group = case_when(
    Sample %in% liver_yng_sed_samples     ~ "Young_Sed",
    Sample %in% liver_old_sed_samples     ~ "Old_SedVeh",
    Sample %in% liver_old_pwr_samples     ~ "Old_PwrVeh",
    Sample %in% liver_old_pwrirap_samples ~ "Old_PwrIRap",
    Sample %in% liver_old_pwrfrap_samples ~ "Old_PwrFRap",
    TRUE                                 ~ "Other"
  ))

#-----------------------------------------------------------
# 6) Plot Projection
#-----------------------------------------------------------
ggplot(liver_proj_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Liver Samples Projected onto DIABLO Aging Axis (Comp1)",
       y = "DIABLO Component 1 Score", x = "") +
  theme(legend.position = "none")

#-----------------------------------------------------------
# 7) Statistics: ANOVA + Tukey post-hoc
#-----------------------------------------------------------
liver_anova <- aov(Comp1 ~ Group, data = liver_proj_df)
summary(liver_anova)

liver_tukey <- TukeyHSD(liver_anova)
print(liver_tukey)

liver_tukey_df <- broom::tidy(liver_tukey)

# Optional: visualize differences
ggplot(liver_proj_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Liver DIABLO Aging Axis with Statistical Groups",
       y = "Component 1 Score", x = "") +
  theme(legend.position = "none")

#-----------------------------------------------------------
# 8) Save outputs
#-----------------------------------------------------------
saveRDS(liver_proj_df, "Liver_DIABLO_AgingAxis_Projection.rds")
write.csv(liver_tukey_df, "Liver_DIABLO_AgingAxis_Tukey.csv", row.names = FALSE)

cat("✅ Liver DIABLO aging-axis projection completed and saved.\n")


#=============================================================#
# Liver DIABLO Aging Axis Plot with Statistics
#=============================================================#

suppressPackageStartupMessages({
  library(rstatix)
  library(ggpubr)
  library(dplyr)
  library(ggplot2)
})

# Ensure group factor order
liver_proj_df <- liver_proj_df %>%
  mutate(Group = factor(Group, 
                        levels = c("Young_Sed", "Old_SedVeh", 
                                   "Old_PwrVeh", "Old_PwrIRap", "Old_PwrFRap")))

#--------------------------------------------------------------#
# 1) One-way ANOVA and Tukey posthoc
#--------------------------------------------------------------#
liver_anova_res <- anova_test(data = liver_proj_df, dv = Comp1, between = Group)
print(liver_anova_res)

liver_tukey_res <- liver_proj_df %>%
  tukey_hsd(Comp1 ~ Group)
print(liver_tukey_res)

#--------------------------------------------------------------#
# 2) Define pairwise comparisons to display on plot
#--------------------------------------------------------------#
liver_comparisons <- list(
  c("Old_SedVeh", "Young_Sed"),   # oldsed vs youngsed
  c("Old_SedVeh", "Old_PwrVeh"),  # oldsed vs oldpwr
  c("Old_SedVeh", "Old_PwrIRap"), # oldpwr vs oldpwrirap
  c("Old_SedVeh", "Old_PwrFRap")  # oldpwr vs oldpwrfrap
)

#--------------------------------------------------------------#
# 3) Create ggplot with significance annotations
#--------------------------------------------------------------#
liver_DIABLO_plot <- ggplot(liver_proj_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Liver Samples Projected onto DIABLO Aging Axis (Comp1)",
    y = "DIABLO Component 1 Score",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  ) +
  stat_compare_means(comparisons = liver_comparisons,
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = TRUE)

#--------------------------------------------------------------#
# 4) Save to PDF
#--------------------------------------------------------------#
pdf("Liver_DIABLO_AgingAxis_Boxplot.pdf", width = 6, height = 4)
print(liver_DIABLO_plot)
dev.off()

cat("✅ Liver DIABLO Aging-Axis plot with stats saved to 'Liver_DIABLO_AgingAxis_Boxplot.pdf'\n")


head(liver_proj_df)
# Save as CSV for record-keeping or plotting in other tools
write.csv(liver_proj_df, "Liver_DIABLO_AgingAxis_SampleScores.csv", row.names = FALSE)





















#============================================================#
#---------Evaluate Aging Axis Features------------------#
#============================================================#

#==============RNA Evalution================================#

#-------Read in xlsx, Clean Dataframe and add ENTREZID col-----------#
read_clean_xlsx <- function(xlsx_file, logFC_col, FDR_col) {
  # Read in the xlsx file
  df <- read_xlsx(xlsx_file)
  
  # Clean and rename columns
  df <- df %>%
    dplyr::select(Ensembl, Symbol, 
                  logFC = all_of(logFC_col),
                  FDR = all_of(FDR_col)) %>%
    dplyr::mutate(
      Ensembl_noDec = str_remove(Ensembl, "\\..+"),
      ENTREZID = mapIds(org.Mm.eg.db,
                        keys = Ensembl_noDec,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")
    ) %>%
    dplyr::select(-Ensembl) %>%
    tidyr::drop_na(ENTREZID)%>%
    as.data.frame()
  
  return(df)
}
#----------------------------------------------#

#-----------------------------------------------------------------#
#----OldSedVeh vs YngSedVeh ----#
liver_os_v_ys_genes <- read_clean_xlsx(
  xlsx_file = "20250922_M007988_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col = "OLD_SED_VEH-YNG_SED_VEH_FDR")

liver_os_v_ys_genes %>%
  filter(FDR <0.05)

liver_os_v_ys_genes

Diablo_liver_latent_features <- read.csv("Liver_DIABLO_selected_features_comp1.csv")

Diablo_liver_latent_features

liver_rna_features <- Diablo_liver_latent_features %>%
  dplyr::filter(block == "RNA") %>%
  dplyr::pull(feature)

liver_rna_features

liver_os_v_ys_genes %>%
  filter(Symbol %in% liver_rna_features) %>%
  arrange(logFC)
#-------------------------------------------#

#------GO for liver latent RNA features------#
#===========================================================#
#       GO Enrichment of DIABLO RNA Features (Liver)        #
#===========================================================#

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(enrichplot)
})

#-----------------------------------------------------------#
# 1. Prepare gene lists
#-----------------------------------------------------------#
# Background universe = all genes tested
bg_genes <- liver_os_v_ys_genes %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

# DIABLO-selected genes (foreground)
diablo_genes <- liver_os_v_ys_genes %>%
  dplyr::filter(Symbol %in% liver_rna_features) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

#-----------------------------------------------------------#
# 2. Run enrichment (Biological Process)
#-----------------------------------------------------------#
ego <- enrichGO(
  gene          = diablo_genes,
  universe      = bg_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",            # "BP" = Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE             # convert to gene symbols
)

#-----------------------------------------------------------#
# 3. Review and visualize
#-----------------------------------------------------------#
# Top enriched terms
head(ego, 10)

# Dotplot
dotplot(ego, showCategory = 15, title = "GO Enrichment: Liver DIABLO RNA Features")

# Barplot (optional)
barplot(ego, showCategory = 15, title = "GO: Biological Process - DIABLO RNA Features")

#-----------------------------------------------------------#
# 4. Save results
#-----------------------------------------------------------#
ego_df <- as.data.frame(ego)
write.csv(ego_df, "Liver_DIABLO_RNA_GOenrichment.csv", row.names = FALSE)
#-----------------------------------------------------------------------#

## ------------------------------
## Convenience functions from your notes
## ------------------------------
# Read, clean, and create ENTREZID column from an edgeR *_GENE_*.xlsx table
read_clean_xlsx <- function(xlsx_file, logFC_col, FDR_col) {
  df <- read_xlsx(xlsx_file) %>%
    dplyr::select(Ensembl, Symbol,
                  logFC = all_of(logFC_col),
                  FDR   = all_of(FDR_col)) %>%
    dplyr::mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
                  ENTREZID = mapIds(org.Mm.eg.db,
                                    keys = Ensembl_noDec,
                                    column = "ENTREZID",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")) %>%
    dplyr::select(-Ensembl) %>%
    tidyr::drop_na(ENTREZID) %>%
    as.data.frame()
  df
}

# Create a named vector (ENTREZID -> logFC) for GSEA
create_vec_from_df <- function(df, df_name) {
  vec <- df %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(ENTREZID, logFC) %>%
    tidyr::drop_na() %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
    tibble::deframe()
  assign(paste0(df_name, ".vec"), vec, envir = .GlobalEnv)
}


############################################################
# FIGURE 1g – GO Pathway Dotplot (OLD SED VEH vs YNG SED VEH)
############################################################
# Input: 20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx
# Steps: GSEA (GO) -> simplify -> dotplot

# 1) Load edgeR *_GENE_* xlsx and build ranked vector
liver_os_v_ys_genes <- read_clean_xlsx(
  xlsx_file = "20250922_M007988_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col = "OLD_SED_VEH-YNG_SED_VEH_FDR")

# quick sanity check (expected ~701 DEG at FDR<0.05 per your note)
liver_os_v_ys_genes %>% filter(FDR < 0.05) %>% nrow()

create_vec_from_df(liver_os_v_ys_genes, "liver_os_v_ys_genes")

# 2) GSEA (GO Biological Process)
gseGO_liver_os_v_ys_genes.OUTPUT <- gseGO(
  geneList     = liver_os_v_ys_genes.vec,  # or your oldsed_vs_yngsed_genes.vec
  ont          = "BP",
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 300,
  pvalueCutoff = 0.05,      # capture everything, filter later
  eps          = 1e-30,      # better estimation for very small p-values
  verbose      = FALSE
)

# Add Count (genes in leading edge) before plotting
gseGO_liver_os_v_ys.df <- as.data.frame(gseGO_liver_os_v_ys_genes.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )



# Your manuscript dotplot (unchanged, just point it at the simplified df)
gseGO_liver_os_v_ys.df %>%
  dplyr::filter(p.adjust < 0.01) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
#-----------------------------------------------------#

#======================================================#
#-----Evaluate Lipid Features--------------------------#

# Verify column order matches expectations
stopifnot(identical(colnames(liver_lipid_log), liver_all_samples))

#----------------------------------------------------------#
# 2. Create long-format table with Group labels
#----------------------------------------------------------#
liver_lipid_long <- liver_lipid_log %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Lipid") %>%
  tidyr::pivot_longer(-Lipid, names_to = "Sample", values_to = "log2Abund") %>%
  mutate(Group = ifelse(Sample %in% liver_yng_sed_samples, "YNG_SED", "OLD_SED"))

#----------------------------------------------------------#
# 3. Per-lipid t-test, log2FC (OLD − YNG)
#----------------------------------------------------------#
liver_lipid_stats <- liver_lipid_long %>%
  group_by(Lipid) %>%
  summarise(
    Mean_YNG = mean(log2Abund[Group == "YNG_SED"], na.rm = TRUE),
    Mean_OLD = mean(log2Abund[Group == "OLD_SED"], na.rm = TRUE),
    Log2FC_OLD_vs_YNG = Mean_OLD - Mean_YNG,
    P_value = tryCatch(t.test(log2Abund ~ Group)$p.value, error = function(e) NA_real_)
  ) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  arrange(FDR)

#----------------------------------------------------------#
# 4. Save and inspect
#----------------------------------------------------------#
write.csv(liver_lipid_stats, "Liver_Lipid_Stats_OLD_vs_YNG.csv", row.names = FALSE)

head(liver_lipid_stats, 10)

liver_lipid_stats %>%
  filter(FDR <0.05)

liver_lipid_features <- Diablo_liver_latent_features %>%
  dplyr::filter(block == "Lipid") %>%
  dplyr::pull(feature)

liver_lipid_stats %>%
  filter(Lipid %in% liver_lipid_features) %>%
  arrange(Log2FC_OLD_vs_YNG)
#-----------------------------------------------------------#


#======================================================#
#-----Evaluate Metabolite Features--------------------------#
head(liver_metabo_log)

#----------------------------------------------------------#
# 2. Create long-format table with Group labels
#----------------------------------------------------------#
liver_metabo_long <- liver_metabo_log %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Metabolite") %>%
  tidyr::pivot_longer(-Metabolite, names_to = "Sample", values_to = "log2Abund") %>%
  mutate(Group = ifelse(Sample %in% liver_yng_sed_samples, "YNG_SED", "OLD_SED"))

#----------------------------------------------------------#
# 3. Per-metabolite t-test, log2FC (OLD − YNG)
#----------------------------------------------------------#
liver_metabo_stats <- liver_metabo_long %>%
  group_by(Metabolite) %>%
  summarise(
    Mean_YNG = mean(log2Abund[Group == "YNG_SED"], na.rm = TRUE),
    Mean_OLD = mean(log2Abund[Group == "OLD_SED"], na.rm = TRUE),
    Log2FC_OLD_vs_YNG = Mean_OLD - Mean_YNG,
    P_value = tryCatch(t.test(log2Abund ~ Group)$p.value, error = function(e) NA_real_)
  ) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  arrange(FDR)

#----------------------------------------------------------#
# 4. Save and preview
#----------------------------------------------------------#
write.csv(liver_metabo_stats, "Liver_Metabolite_Stats_OLD_vs_YNG.csv", row.names = FALSE)

head(liver_metabo_stats, 10)

liver_metabo_stats %>%
  filter(P_value <0.05)

liver_metabo_features <- Diablo_liver_latent_features %>%
  dplyr::filter(block == "Metabolite") %>%
  dplyr::pull(feature)

liver_metabo_stats %>%
  filter(Metabolite %in% liver_metabo_features) %>%
  arrange(P_value)
