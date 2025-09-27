############################################################
# RapaPWR – Figure 1 consolidated code (minimal edits)
# Goal: Recreate panels 1g–1k using your original code blocks,
#       organized and lightly stitched so it runs top-to-bottom.
# Notes:
#  • This sticks closely to your style; only small add-ons where needed.
#  • Package installs are commented out (journal best practice is to list
#    packages in Methods/Session info rather than auto-install in scripts).
#  • Assumes you’ve set the working directory to the folder that holds the
#    listed files.
#  • Saves simple CSVs of key outputs (core genes, logCPM matrices, etc.).
#  • You must provide `validated_inflamm_geneset.vec` (MusAge 52 genes, SYMBOLs
#    or ENTREZ IDs; see Section 1j / 1k notes).
############################################################
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")

## ------------------------------
## Packages
## ------------------------------
# install.packages(c("readxl","dplyr","tidyverse","readr","stringr","ggplot2","ggrepel","pheatmap","VennDiagram","ggVennDiagram","eulerr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler","org.Mm.eg.db","AnnotationDbi","DOSE","fgsea","ComplexHeatmap","edgeR"))

library(readxl)
library(dplyr)
library(tidyverse)
library(readr)
library(ggVennDiagram)
library(ggplot2)
library(eulerr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DOSE)
library(pheatmap)
library(ComplexHeatmap)
library(fgsea)
library(ggrepel)
library(stringr)
library(edgeR)
library(tibble)

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
oldsedveh_v_yngsedveh_genes_12jun <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col   = "OLD_SED_VEH-YNG_SED_VEH_FDR"
)

# quick sanity check (expected ~701 DEG at FDR<0.05 per your note)
oldsedveh_v_yngsedveh_genes_12jun %>% filter(FDR < 0.05) %>% nrow()

create_vec_from_df(oldsedveh_v_yngsedveh_genes_12jun, "oldsedveh_v_yngsedveh_12jun")

# 2) GSEA (GO Biological Process)
gseGO_oldsed_vs_yngsed.OUTPUT <- gseGO(
  geneList     = oldsedveh_v_yngsedveh_12jun.vec,  # or your oldsed_vs_yngsed_genes.vec
  ont          = "BP",
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 300,
  pvalueCutoff = 0.5,      # capture everything, filter later
  eps          = 1e-30,      # better estimation for very small p-values
  verbose      = FALSE
)

# Add Count (genes in leading edge) before plotting
gseGO_oldsed_vs_yngsed.df <- as.data.frame(gseGO_oldsed_vs_yngsed.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

# Simplify, then compute Count/geneRatio again (simplify changes the table)
oldsed_vs_yngsed_simpl <- simplify(gseGO_oldsed_vs_yngsed.OUTPUT,
                                   cutoff = 0.5, by = "p.adjust", select_fun = min)

oldsed_vs_yngsed_12jun_simpl.df <- as.data.frame(oldsed_vs_yngsed_simpl@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )


# Your manuscript dotplot (unchanged, just point it at the simplified df)
oldsed_vs_yngsed_12jun_simpl.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()

# Save (optional)
# ggsave("Fig1g_GO_dotplot.pdf", oldsed_vs_yngsed_simpl_GO_dotplot, width = 8, height = 8)


############################################################
# FIGURE 1h — Heatmap of Core Enrichment genes (UP GO pathways)
# Matches your original naming/flow + rows ordered by FDR
############################################################

# 1) Core enrichment ENTREZ IDs (from simplified GSEA; upregulated & adj.P<0.05)
upreg_unq_old_sed_core_entrez_vec <- oldsed_vs_yngsed_12jun_simpl.df %>%
  dplyr::filter(p.adjust < 0.05, NES > 0) %>%
  dplyr::pull(core_enrichment) %>%
  strsplit("/") %>% unlist() %>% unique() %>% sort()

# quick check (you had ~498–499)
length(upreg_unq_old_sed_core_entrez_vec)
head(upreg_unq_old_sed_core_entrez_vec)

# 2) Load counts and map to ENTREZID (same as your original)
oldsedveh_vs_yngsedveh_counts <- readxl::read_xlsx(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx"
)
colnames(oldsedveh_vs_yngsedveh_counts)[1] <- "Ensembl"

oldsedveh_vs_yngsedveh_counts <- oldsedveh_vs_yngsedveh_counts %>%
  dplyr::mutate(Ensembl_noDec = stringr::str_remove(Ensembl, "\\..+")) %>%
  dplyr::select(-Ensembl)

oldsedveh_vs_yngsedveh_counts$ENTREZID <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = oldsedveh_vs_yngsedveh_counts$Ensembl_noDec,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

# drop rows without ENTREZID
oldsedveh_vs_yngsedveh_counts <- tidyr::drop_na(oldsedveh_vs_yngsedveh_counts, ENTREZID)

# 3) Your sample order (unchanged)
yng_sed_samples <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
ordered_columns <- c(yng_sed_samples, old_sed_samples)

# counts matrix (rows = ENTREZID, cols = samples)
counts_matrix <- oldsedveh_vs_yngsedveh_counts %>%
  dplyr::select(dplyr::all_of(ordered_columns)) %>%
  as.matrix()
rownames(counts_matrix) <- as.character(oldsedveh_vs_yngsedveh_counts$ENTREZID)  # ensure character IDs

# 4) edgeR normalization → logCPM (same as you had)
group <- factor(c(rep("YNG_SED", length(yng_sed_samples)),
                  rep("OLD_SED", length(old_sed_samples))))
dge <- edgeR::DGEList(counts = counts_matrix, group = group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logCPM_matrix <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# 5) Order core genes by FDR from your genes table, then plot heatmap (blue→white→red)

# Make sure IDs are the same type as the logCPM rownames
core_ids_chr <- as.character(upreg_unq_old_sed_core_entrez_vec)

# Get ENTREZIDs in FDR order from your genes table
ordered_core_ids <- oldsedveh_v_yngsedveh_genes_12jun %>%
  dplyr::filter(ENTREZID %in% core_ids_chr) %>%
  dplyr::arrange(FDR) %>%                             # lowest FDR first
  dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%     # guard against dup rows
  dplyr::pull(ENTREZID) %>%
  as.character()

# Keep only those present in the CPM matrix (just in case)
ordered_core_ids <- ordered_core_ids[ordered_core_ids %in% rownames(logCPM_matrix)]

# Build the matrix in that exact FDR order (columns already defined)
all_oldsed_core_enrich_heatmap <- logCPM_matrix[ordered_core_ids, ordered_columns, drop = FALSE]

# Heatmap: preserve your FDR order (no row clustering)
all_oldsed_core_enrich_heatplot <- pheatmap::pheatmap(
  all_oldsed_core_enrich_heatmap,
  cluster_rows = FALSE,            # ← key change: keep FDR order
  cluster_cols = FALSE,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue","white","red"))(50)
)

# Save as PDF (same filename you used)
pdf("all_oldsed_core_enrich_heatplot.pdf", width = 3.5, height = 6)
print(all_oldsed_core_enrich_heatplot)
dev.off()

############################################################
# FIGURE 1i – Heatmap of 66 core genes elevated in OLD SED (FDR<0.05 & logFC>0)
############################################################
# Inputs: same *_GENE_* xlsx + counts logCPM above
# Steps: subset the 498 core genes to those with FDR<0.05 AND logFC>0 in OLD SED vs YNG SED

# --- compatibility aliases so older blocks run unchanged ---
logCPM <- logCPM_matrix
common_samples <- ordered_columns
upreg_old_sed_core_entrez_vec_12jun <- upreg_unq_old_sed_core_entrez_vec

# 1) Subset gene table to core and significant upregulated
Inflammaging_gene_set_12jun <- oldsedveh_v_yngsedveh_genes_12jun %>%
  filter(ENTREZID %in% upreg_old_sed_core_entrez_vec_12jun) %>%
  filter(FDR < 0.05, logFC > 0)
dim(Inflammaging_gene_set_12jun)

# 2) Heatmap on those (expect ~66 genes; exact number may vary run-to-run)
i66_ids <- Inflammaging_gene_set_12jun$ENTREZID
i66_logCPM <- logCPM[rownames(logCPM) %in% i66_ids, common_samples]

pheatmap(i66_logCPM,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue","white","red"))(100),
         show_rownames = FALSE,
         main = "66 core genes (FDR<0.05 & logFC>0 in OLD SED)")

# Save CSVs
write_csv(Inflammaging_gene_set_12jun %>% select(Symbol, ENTREZID, logFC, FDR), "Fig1i_inflammaging_gene_set.csv")
write_csv(as.data.frame(i66_logCPM) %>% rownames_to_column("ENTREZID"), "Fig1i_i66_logCPM.csv")


############################################################
# FIGURE 1j – Aging mouse re-analysis (GSE226117) – dot plot of MusAge 52
############################################################
# Inputs:
#  • GEO RAW folder: GSE226117_RAW/ (files like GSM7064422_X1385229_raw_counts.txt)
#  • Output RData (optional): all_muscle_age_dfs.RData
#  • Gene set: validated_inflamm_geneset.vec  (== MusAge 52 genes)
# Brief:
#  • Merge all raw counts per muscle
#  • edgeR across ages (12,18,21,24,27 vs 6mo), by sex and muscle
#  • Age-responsive if (≥2 of {21,24,27} significant vs 6) OR Spearman rho increases with age
#  • Dotplot of 52 genes (Female/Male × Gastroc/TA/Soleus); green = age-responsive

## =========================
## Fig. 1j – MusAge dot plot
## =========================

# --- MusAge 52 (symbols) + ENTREZ alias (for anything that needs it later) ---
validated_inflamm_geneset.vec <- readr::read_csv("MusAge_geneset.csv", show_col_types = FALSE) |>
  dplyr::pull(Gene) |>
  as.character() |> stringr::str_trim() |> unique()

Inflammaging_symbols <- validated_inflamm_geneset.vec

inflammaging_gene_vec <- AnnotationDbi::mapIds(
  org.Mm.eg.db, keys = Inflammaging_symbols,
  keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"
) |> unname() |> unique()

# --- Setup (edit if needed) ---
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR/GSE226117_RAW")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(purrr)
  library(tidyr)
})

# --- MusAge 52 (symbols) + ENTREZ alias (for anything that needs it later) ---
validated_inflamm_geneset.vec <- readr::read_csv("../MusAge_geneset.csv", show_col_types = FALSE) |>
  dplyr::pull(Gene) |>
  as.character() |> stringr::str_trim() |> unique()

Inflammaging_symbols <- validated_inflamm_geneset.vec

inflammaging_gene_vec <- AnnotationDbi::mapIds(
  org.Mm.eg.db, keys = Inflammaging_symbols,
  keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"
) |> unname() |> unique()

# ===========================================
# Your helper: edgeR per age vs 6mo (robust)
# ===========================================
run_edger_age_comparison <- function(file_list, group1_ids, group2_ids, muscle_name, older_month) {
  # file_list: character vector of filenames (present in wd)
  # group*_ids: character vector of "X\\d+" IDs that appear in filenames
  # Outputs a data.frame with GeneName, logFC_<month>, FDR_<month>
  
  # Merge all files by GeneName, each column = a sample (Xid)
  get_sample_id <- function(f) stringr::str_extract(basename(f), "X\\d+")
  merge_one <- function(f) {
    dt <- data.table::fread(f)
    # force first column to GeneName & use the 2nd column as counts
    data.table::setnames(dt, old = names(dt)[1], new = "GeneName")
    count_col <- names(dt)[2]
    out <- dt[, .(GeneName, count = get(count_col))]
    # collapse duplicate symbols if any
    out <- out[, .(count = mean(count, na.rm = TRUE)), by = GeneName]
    data.table::setnames(out, "count", get_sample_id(f))
    data.table::setkey(out, GeneName)
    out
  }
  
  # build merged table
  stopifnot(length(file_list) >= 2)
  merged <- merge_one(file_list[1])
  if (length(file_list) > 1) {
    for (f in file_list[-1]) merged <- merged[merge_one(f), on = "GeneName"]
  }
  
  # counts matrix
  counts_mat <- as.matrix(merged[, -1])
  rownames(counts_mat) <- merged$GeneName
  
  all_samples <- c(group1_ids, group2_ids)
  # keep only columns we need (present in counts_mat)
  keep <- intersect(all_samples, colnames(counts_mat))
  if (length(keep) < 2) stop("Not enough matching samples found for this comparison.")
  counts_mat <- counts_mat[, keep, drop = FALSE]
  
  group <- factor(c(rep("mo6", sum(keep %in% group1_ids)),
                    rep(paste0("mo", older_month), sum(keep %in% group2_ids))))
  y <- edgeR::DGEList(counts = counts_mat, group = group)
  y <- y[edgeR::filterByExpr(y), , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  group <- stats::relevel(group, ref = "mo6")
  design <- model.matrix(~ group)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmFit(y, design)
  lrt <- edgeR::glmLRT(fit, coef = 2)
  
  tab <- edgeR::topTags(lrt, n = Inf)$table
  res <- data.frame(
    GeneName = rownames(tab),
    logFC = tab$logFC,
    FDR   = tab$FDR,
    row.names = NULL, check.names = FALSE
  )
  colnames(res)[2:3] <- c(paste0("logFC_", older_month), paste0("FDR_", older_month))
  res
}

# =====================================================
# File lists + sample groups (YOUR provided definitions)
# (leave exactly as you pasted; truncated here for space)
# =====================================================

## ---- FEMALE: Gastroc ----
file_list_gas <- c(
  "GSM7064465_X1385216_raw_counts.txt","GSM7064466_X1385217_raw_counts.txt",
  "GSM7064467_X1385218_raw_counts.txt","GSM7064468_X1385219_raw_counts.txt",
  "GSM7064469_X1385220_raw_counts.txt","GSM7064470_X1385221_raw_counts.txt",
  "GSM7064471_X1385303_raw_counts.txt","GSM7064472_X1385304_raw_counts.txt",
  "GSM7064421_X1385228_raw_counts.txt","GSM7064422_X1385229_raw_counts.txt",
  "GSM7064423_X1385230_raw_counts.txt","GSM7064424_X1385231_raw_counts.txt",
  "GSM7064425_X1385232_raw_counts.txt","GSM7064426_X1385233_raw_counts.txt",
  "GSM7064427_X1385305_raw_counts.txt","GSM7064428_X1385306_raw_counts.txt",
  "GSM7064429_X1385240_raw_counts.txt","GSM7064430_X1385241_raw_counts.txt",
  "GSM7064431_X1385242_raw_counts.txt","GSM7064432_X1385243_raw_counts.txt",
  "GSM7064433_X1385244_raw_counts.txt","GSM7064434_X1385245_raw_counts.txt",
  "GSM7064435_X1385309_raw_counts.txt","GSM7064436_X1385310_raw_counts.txt",
  "GSM7064437_X1385260_raw_counts.txt","GSM7064438_X1385261_raw_counts.txt",
  "GSM7064439_X1385262_raw_counts.txt","GSM7064440_X1385263_raw_counts.txt",
  "GSM7064441_X1385264_raw_counts.txt","GSM7064442_X1385265_raw_counts.txt",
  "GSM7064443_X1385266_raw_counts.txt","GSM7064444_X1385267_raw_counts.txt",
  "GSM7064445_X1385272_raw_counts.txt","GSM7064446_X1385273_raw_counts.txt",
  "GSM7064447_X1385274_raw_counts.txt","GSM7064448_X1385275_raw_counts.txt",
  "GSM7064449_X1385276_raw_counts.txt","GSM7064450_X1385277_raw_counts.txt",
  "GSM7064451_X1385278_raw_counts.txt","GSM7064452_X1385335_raw_counts.txt",
  "GSM7064453_X1385336_raw_counts.txt","GSM7064454_X1385337_raw_counts.txt",
  "GSM7064455_X1385279_raw_counts.txt","GSM7064456_X1385280_raw_counts.txt",
  "GSM7064457_X1385281_raw_counts.txt","GSM7064458_X1385282_raw_counts.txt",
  "GSM7064459_X1385283_raw_counts.txt","GSM7064460_X1385292_raw_counts.txt",
  "GSM7064461_X1385293_raw_counts.txt","GSM7064462_X1385294_raw_counts.txt",
  "GSM7064463_X1385360_raw_counts.txt","GSM7064464_X1385361_raw_counts.txt"
)
sample_groups_gas <- list(
  "6"  = c("X1385216","X1385217","X1385218","X1385219","X1385220","X1385221","X1385303","X1385304"),
  "12" = c("X1385228","X1385229","X1385230","X1385231","X1385232","X1385233","X1385305","X1385306"),
  "18" = c("X1385240","X1385241","X1385242","X1385243","X1385244","X1385245","X1385309","X1385310"),
  "21" = c("X1385260","X1385261","X1385262","X1385263","X1385264","X1385265","X1385266","X1385267"),
  "24" = c("X1385272","X1385273","X1385274","X1385275","X1385276","X1385277","X1385278","X1385335","X1385336","X1385337"),
  "27" = c("X1385279","X1385280","X1385281","X1385282","X1385283","X1385292","X1385293","X1385294","X1385360","X1385361")
)

## ---- FEMALE: Soleus ----
file_list_sol <- c(
  "GSM7064516_X1384990_raw_counts.txt","GSM7064517_X1384991_raw_counts.txt",
  "GSM7064518_X1384992_raw_counts.txt","GSM7064519_X1384993_raw_counts.txt",
  "GSM7064520_X1384994_raw_counts.txt","GSM7064521_X1384995_raw_counts.txt",
  "GSM7064522_X1385286_raw_counts.txt","GSM7064523_X1385287_raw_counts.txt",
  "GSM7064473_X1385002_raw_counts.txt","GSM7064474_X1385003_raw_counts.txt",
  "GSM7064475_X1385004_raw_counts.txt","GSM7064476_X1385005_raw_counts.txt",
  "GSM7064477_X1385006_raw_counts.txt","GSM7064478_X1385007_raw_counts.txt",
  "GSM7064479_X1385288_raw_counts.txt","GSM7064480_X1385289_raw_counts.txt",
  "GSM7064481_X1385014_raw_counts.txt","GSM7064482_X1385015_raw_counts.txt",
  "GSM7064483_X1385016_raw_counts.txt","GSM7064484_X1385017_raw_counts.txt",
  "GSM7064485_X1385018_raw_counts.txt","GSM7064486_X1385019_raw_counts.txt",
  "GSM7064487_X1385313_raw_counts.txt","GSM7064488_X1385314_raw_counts.txt",
  "GSM7064489_X1385034_raw_counts.txt","GSM7064490_X1385035_raw_counts.txt",
  "GSM7064491_X1385036_raw_counts.txt","GSM7064492_X1385037_raw_counts.txt",
  "GSM7064493_X1385038_raw_counts.txt","GSM7064494_X1385039_raw_counts.txt",
  "GSM7064495_X1385040_raw_counts.txt","GSM7064496_X1385041_raw_counts.txt",
  "GSM7064497_X1385046_raw_counts.txt","GSM7064498_X1385048_raw_counts.txt",
  "GSM7064499_X1385049_raw_counts.txt","GSM7064500_X1385050_raw_counts.txt",
  "GSM7064501_X1385051_raw_counts.txt","GSM7064502_X1385052_raw_counts.txt",
  "GSM7064503_X1385344_raw_counts.txt","GSM7064504_X1385345_raw_counts.txt",
  "GSM7064505_X1385346_raw_counts.txt","GSM7064506_X1385053_raw_counts.txt",
  "GSM7064507_X1385054_raw_counts.txt","GSM7064508_X1385055_raw_counts.txt",
  "GSM7064509_X1385056_raw_counts.txt","GSM7064510_X1385057_raw_counts.txt",
  "GSM7064511_X1385058_raw_counts.txt","GSM7064512_X1385059_raw_counts.txt",
  "GSM7064513_X1385060_raw_counts.txt","GSM7064514_X1385366_raw_counts.txt",
  "GSM7064515_X1385367_raw_counts.txt"
)
sample_groups_sol <- list(
  "6"  = c("X1384990","X1384991","X1384992","X1384993","X1384994","X1384995","X1385286","X1385287"),
  "12" = c("X1385002","X1385003","X1385004","X1385005","X1385006","X1385007","X1385288","X1385289"),
  "18" = c("X1385014","X1385015","X1385016","X1385017","X1385018","X1385019","X1385313","X1385314"),
  "21" = c("X1385034","X1385035","X1385036","X1385037","X1385038","X1385039","X1385040","X1385041"),
  "24" = c("X1385046","X1385048","X1385049","X1385050","X1385051","X1385052","X1385344","X1385345","X1385346"),
  "27" = c("X1385053","X1385054","X1385055","X1385056","X1385057","X1385058","X1385059","X1385060","X1385366","X1385367")
)

## ---- FEMALE: TA ----
file_list_ta <- c(
  "GSM7064567_X1385067_raw_counts.txt","GSM7064568_X1385068_raw_counts.txt",
  "GSM7064569_X1385069_raw_counts.txt","GSM7064570_X1385070_raw_counts.txt",
  "GSM7064571_X1385071_raw_counts.txt","GSM7064572_X1385072_raw_counts.txt",
  "GSM7064573_X1385319_raw_counts.txt","GSM7064574_X1385320_raw_counts.txt",
  "GSM7064524_X1385079_raw_counts.txt","GSM7064525_X1385080_raw_counts.txt",
  "GSM7064526_X1385081_raw_counts.txt","GSM7064527_X1385082_raw_counts.txt",
  "GSM7064528_X1385083_raw_counts.txt","GSM7064529_X1385084_raw_counts.txt",
  "GSM7064530_X1385321_raw_counts.txt","GSM7064531_X1385322_raw_counts.txt",
  "GSM7064532_X1385091_raw_counts.txt","GSM7064533_X1385092_raw_counts.txt",
  "GSM7064534_X1385093_raw_counts.txt","GSM7064535_X1385094_raw_counts.txt",
  "GSM7064536_X1385095_raw_counts.txt","GSM7064537_X1385096_raw_counts.txt",
  "GSM7064538_X1385326_raw_counts.txt","GSM7064539_X1385111_raw_counts.txt",
  "GSM7064540_X1385112_raw_counts.txt","GSM7064541_X1385113_raw_counts.txt",
  "GSM7064542_X1385114_raw_counts.txt","GSM7064543_X1385115_raw_counts.txt",
  "GSM7064544_X1385116_raw_counts.txt","GSM7064545_X1385117_raw_counts.txt",
  "GSM7064546_X1385118_raw_counts.txt","GSM7064547_X1385123_raw_counts.txt",
  "GSM7064548_X1385124_raw_counts.txt","GSM7064549_X1385125_raw_counts.txt",
  "GSM7064550_X1385126_raw_counts.txt","GSM7064551_X1385127_raw_counts.txt",
  "GSM7064552_X1385128_raw_counts.txt","GSM7064553_X1385129_raw_counts.txt",
  "GSM7064554_X1385353_raw_counts.txt","GSM7064555_X1385354_raw_counts.txt",
  "GSM7064556_X1385355_raw_counts.txt","GSM7064557_X1385130_raw_counts.txt",
  "GSM7064558_X1385131_raw_counts.txt","GSM7064559_X1385132_raw_counts.txt",
  "GSM7064560_X1385133_raw_counts.txt","GSM7064561_X1385134_raw_counts.txt",
  "GSM7064562_X1385135_raw_counts.txt","GSM7064563_X1385136_raw_counts.txt",
  "GSM7064564_X1385137_raw_counts.txt","GSM7064565_X1385372_raw_counts.txt",
  "GSM7064566_X1385373_raw_counts.txt"
)
sample_groups_ta <- list(
  "6"  = c("X1385067","X1385068","X1385069","X1385070","X1385071","X1385072","X1385319","X1385320"),
  "12" = c("X1385079","X1385080","X1385081","X1385082","X1385083","X1385084","X1385321","X1385322"),
  "18" = c("X1385091","X1385092","X1385093","X1385094","X1385095","X1385096","X1385326"),
  "21" = c("X1385111","X1385112","X1385113","X1385114","X1385115","X1385116","X1385117","X1385118"),
  "24" = c("X1385123","X1385124","X1385125","X1385126","X1385127","X1385128","X1385129","X1385353","X1385354","X1385355"),
  "27" = c("X1385130","X1385131","X1385132","X1385133","X1385134","X1385135","X1385136","X1385137","X1385372","X1385373")
)

## ---- MALE: Gastroc ----
file_list_m_gas <- c(
  "GSM7064667_X1385222_raw_counts.txt","GSM7064668_X1385223_raw_counts.txt",
  "GSM7064669_X1385224_raw_counts.txt","GSM7064670_X1385225_raw_counts.txt",
  "GSM7064671_X1385226_raw_counts.txt","GSM7064672_X1385227_raw_counts.txt",
  "GSM7064673_X1385301_raw_counts.txt","GSM7064674_X1385302_raw_counts.txt",
  "GSM7064623_X1385234_raw_counts.txt","GSM7064624_X1385235_raw_counts.txt",
  "GSM7064625_X1385236_raw_counts.txt","GSM7064626_X1385237_raw_counts.txt",
  "GSM7064627_X1385238_raw_counts.txt","GSM7064628_X1385239_raw_counts.txt",
  "GSM7064629_X1385307_raw_counts.txt","GSM7064630_X1385308_raw_counts.txt",
  "GSM7064631_X1385246_raw_counts.txt","GSM7064632_X1385247_raw_counts.txt",
  "GSM7064633_X1385248_raw_counts.txt","GSM7064634_X1385249_raw_counts.txt",
  "GSM7064635_X1385250_raw_counts.txt","GSM7064636_X1385251_raw_counts.txt",
  "GSM7064637_X1385311_raw_counts.txt","GSM7064638_X1385312_raw_counts.txt",
  "GSM7064639_X1385252_raw_counts.txt","GSM7064640_X1385253_raw_counts.txt",
  "GSM7064641_X1385254_raw_counts.txt","GSM7064642_X1385255_raw_counts.txt",
  "GSM7064643_X1385256_raw_counts.txt","GSM7064644_X1385257_raw_counts.txt",
  "GSM7064645_X1385258_raw_counts.txt","GSM7064646_X1385259_raw_counts.txt",
  "GSM7064647_X1385268_raw_counts.txt","GSM7064648_X1385269_raw_counts.txt",
  "GSM7064649_X1385270_raw_counts.txt","GSM7064650_X1385271_raw_counts.txt",
  "GSM7064651_X1385329_raw_counts.txt","GSM7064652_X1385330_raw_counts.txt",
  "GSM7064653_X1385331_raw_counts.txt","GSM7064654_X1385332_raw_counts.txt",
  "GSM7064655_X1385333_raw_counts.txt","GSM7064656_X1385334_raw_counts.txt",
  "GSM7064657_X1385295_raw_counts.txt","GSM7064658_X1385296_raw_counts.txt",
  "GSM7064659_X1385297_raw_counts.txt","GSM7064660_X1385298_raw_counts.txt",
  "GSM7064661_X1385299_raw_counts.txt","GSM7064662_X1385300_raw_counts.txt",
  "GSM7064663_X1385356_raw_counts.txt","GSM7064664_X1385357_raw_counts.txt",
  "GSM7064665_X1385358_raw_counts.txt","GSM7064666_X1385359_raw_counts.txt"
)
sample_groups_m_gas <- list(
  "6"  = c("X1385222","X1385223","X1385224","X1385225","X1385226","X1385227","X1385301","X1385302"),
  "12" = c("X1385234","X1385235","X1385236","X1385237","X1385238","X1385239","X1385307","X1385308"),
  "18" = c("X1385246","X1385247","X1385248","X1385249","X1385250","X1385251","X1385311","X1385312"),
  "21" = c("X1385252","X1385253","X1385254","X1385255","X1385256","X1385257","X1385258","X1385259"),
  "24" = c("X1385268","X1385269","X1385270","X1385271","X1385329","X1385330","X1385331","X1385332","X1385333","X1385334"),
  "27" = c("X1385295","X1385296","X1385297","X1385298","X1385299","X1385300","X1385356","X1385357","X1385358","X1385359")
)

## ---- MALE: Soleus ----
file_list_m_sol <- c(
  "GSM7064719_X1384996_raw_counts.txt","GSM7064720_X1384997_raw_counts.txt",
  "GSM7064721_X1384998_raw_counts.txt","GSM7064722_X1384999_raw_counts.txt",
  "GSM7064723_X1385000_raw_counts.txt","GSM7064724_X1385001_raw_counts.txt",
  "GSM7064725_X1385284_raw_counts.txt","GSM7064726_X1385285_raw_counts.txt",
  "GSM7064675_X1385008_raw_counts.txt","GSM7064676_X1385009_raw_counts.txt",
  "GSM7064677_X1385010_raw_counts.txt","GSM7064678_X1385011_raw_counts.txt",
  "GSM7064679_X1385012_raw_counts.txt","GSM7064680_X1385013_raw_counts.txt",
  "GSM7064681_X1385290_raw_counts.txt","GSM7064682_X1385291_raw_counts.txt",
  "GSM7064683_X1385020_raw_counts.txt","GSM7064684_X1385021_raw_counts.txt",
  "GSM7064685_X1385022_raw_counts.txt","GSM7064686_X1385023_raw_counts.txt",
  "GSM7064687_X1385024_raw_counts.txt","GSM7064688_X1385025_raw_counts.txt",
  "GSM7064689_X1385315_raw_counts.txt","GSM7064690_X1385316_raw_counts.txt",
  "GSM7064691_X1385026_raw_counts.txt","GSM7064692_X1385027_raw_counts.txt",
  "GSM7064693_X1385028_raw_counts.txt","GSM7064694_X1385029_raw_counts.txt",
  "GSM7064695_X1385030_raw_counts.txt","GSM7064696_X1385031_raw_counts.txt",
  "GSM7064697_X1385032_raw_counts.txt","GSM7064698_X1385033_raw_counts.txt",
  "GSM7064699_X1385042_raw_counts.txt","GSM7064700_X1385043_raw_counts.txt",
  "GSM7064701_X1385044_raw_counts.txt","GSM7064702_X1385045_raw_counts.txt",
  "GSM7064703_X1385338_raw_counts.txt","GSM7064704_X1385339_raw_counts.txt",
  "GSM7064705_X1385340_raw_counts.txt","GSM7064706_X1385341_raw_counts.txt",
  "GSM7064707_X1385342_raw_counts.txt","GSM7064708_X1385343_raw_counts.txt",
  "GSM7064709_X1385061_raw_counts.txt","GSM7064710_X1385062_raw_counts.txt",
  "GSM7064711_X1385063_raw_counts.txt","GSM7064712_X1385064_raw_counts.txt",
  "GSM7064713_X1385065_raw_counts.txt","GSM7064714_X1385066_raw_counts.txt",
  "GSM7064715_X1385362_raw_counts.txt","GSM7064716_X1385363_raw_counts.txt",
  "GSM7064717_X1385364_raw_counts.txt","GSM7064718_X1385365_raw_counts.txt"
)
sample_groups_m_sol <- list(
  "6"  = c("X1384996","X1384997","X1384998","X1384999","X1385000","X1385001","X1385284","X1385285"),
  "12" = c("X1385008","X1385009","X1385010","X1385011","X1385012","X1385013","X1385290","X1385291"),
  "18" = c("X1385020","X1385021","X1385022","X1385023","X1385024","X1385025","X1385315","X1385316"),
  "21" = c("X1385026","X1385027","X1385028","X1385029","X1385030","X1385031","X1385032","X1385033"),
  "24" = c("X1385042","X1385043","X1385044","X1385045","X1385338","X1385339","X1385340","X1385341","X1385342","X1385343"),
  "27" = c("X1385061","X1385062","X1385063","X1385064","X1385065","X1385066","X1385362","X1385363","X1385364","X1385365")
)

## ---- MALE: TA ----
file_list_m_ta <- c(
  "GSM7064770_X1385073_raw_counts.txt","GSM7064771_X1385074_raw_counts.txt",
  "GSM7064772_X1385075_raw_counts.txt","GSM7064773_X1385076_raw_counts.txt",
  "GSM7064774_X1385077_raw_counts.txt","GSM7064775_X1385078_raw_counts.txt",
  "GSM7064776_X1385317_raw_counts.txt","GSM7064777_X1385318_raw_counts.txt",
  "GSM7064727_X1385085_raw_counts.txt","GSM7064728_X1385086_raw_counts.txt",
  "GSM7064729_X1385087_raw_counts.txt","GSM7064730_X1385088_raw_counts.txt",
  "GSM7064731_X1385089_raw_counts.txt","GSM7064732_X1385090_raw_counts.txt",
  "GSM7064733_X1385323_raw_counts.txt","GSM7064734_X1385324_raw_counts.txt",
  "GSM7064735_X1385097_raw_counts.txt","GSM7064736_X1385098_raw_counts.txt",
  "GSM7064737_X1385099_raw_counts.txt","GSM7064738_X1385100_raw_counts.txt",
  "GSM7064739_X1385102_raw_counts.txt","GSM7064740_X1385327_raw_counts.txt",
  "GSM7064741_X1385328_raw_counts.txt","GSM7064742_X1385103_raw_counts.txt",
  "GSM7064743_X1385104_raw_counts.txt","GSM7064744_X1385105_raw_counts.txt",
  "GSM7064745_X1385106_raw_counts.txt","GSM7064746_X1385107_raw_counts.txt",
  "GSM7064747_X1385108_raw_counts.txt","GSM7064748_X1385109_raw_counts.txt",
  "GSM7064749_X1385110_raw_counts.txt","GSM7064750_X1385119_raw_counts.txt",
  "GSM7064751_X1385120_raw_counts.txt","GSM7064752_X1385121_raw_counts.txt",
  "GSM7064753_X1385122_raw_counts.txt","GSM7064754_X1385347_raw_counts.txt",
  "GSM7064755_X1385348_raw_counts.txt","GSM7064756_X1385349_raw_counts.txt",
  "GSM7064757_X1385350_raw_counts.txt","GSM7064758_X1385351_raw_counts.txt",
  "GSM7064759_X1385352_raw_counts.txt","GSM7064760_X1385138_raw_counts.txt",
  "GSM7064761_X1385139_raw_counts.txt","GSM7064762_X1385140_raw_counts.txt",
  "GSM7064763_X1385141_raw_counts.txt","GSM7064764_X1385142_raw_counts.txt",
  "GSM7064765_X1385143_raw_counts.txt","GSM7064766_X1385368_raw_counts.txt",
  "GSM7064767_X1385369_raw_counts.txt","GSM7064768_X1385370_raw_counts.txt",
  "GSM7064769_X1385371_raw_counts.txt"
)
sample_groups_m_ta <- list(
  "6"  = c("X1385073","X1385074","X1385075","X1385076","X1385077","X1385078","X1385317","X1385318"),
  "12" = c("X1385085","X1385086","X1385087","X1385088","X1385089","X1385090","X1385323","X1385324"),
  "18" = c("X1385097","X1385098","X1385099","X1385100","X1385102","X1385327","X1385328"),
  "21" = c("X1385103","X1385104","X1385105","X1385106","X1385107","X1385108","X1385109","X1385110"),
  "24" = c("X1385119","X1385120","X1385121","X1385122","X1385347","X1385348","X1385349","X1385350","X1385351","X1385352"),
  "27" = c("X1385138","X1385139","X1385140","X1385141","X1385142","X1385143","X1385368","X1385369","X1385370","X1385371")
)

# ==========================================================
# Wrapper: run all 5 ages vs 6mo, then classify responsiveness
# ==========================================================
run_stratum <- function(file_list, sample_groups, muscle, sex,
                        ages_old = c(12,18,21,24,27),
                        rho_threshold = 0) {
  # map Xid -> filename for this stratum
  xid <- stringr::str_extract(basename(file_list), "X\\d+")
  file_map <- setNames(file_list, xid)
  
  # run each age vs 6mo
  rlist <- list()
  for (age in ages_old) {
    g6 <- sample_groups[["6"]]
    gA <- sample_groups[[as.character(age)]]
    if (is.null(g6) || is.null(gA)) next
    need_files <- unique(c(file_map[g6], file_map[gA]))
    need_files <- need_files[!is.na(need_files)]
    if (length(need_files) == 0) next
    
    res <- run_edger_age_comparison(
      file_list   = need_files,
      group1_ids  = g6,
      group2_ids  = gA,
      muscle_name = muscle,
      older_month = as.character(age)
    )
    rlist[[as.character(age)]] <- res
  }
  
  # merge results into one table
  stopifnot(length(rlist) >= 1)
  merged <- Reduce(function(a,b) dplyr::full_join(a,b, by = "GeneName"), rlist)
  
  # keep genes present in all available contrasts
  merged <- tidyr::drop_na(merged)
  
  # ---- FIX: create sig_21/24/27 outside mutate() ----
  sig_21 <- if ("FDR_21" %in% names(merged)) merged$FDR_21 < 0.05 else rep(FALSE, nrow(merged))
  sig_24 <- if ("FDR_24" %in% names(merged)) merged$FDR_24 < 0.05 else rep(FALSE, nrow(merged))
  sig_27 <- if ("FDR_27" %in% names(merged)) merged$FDR_27 < 0.05 else rep(FALSE, nrow(merged))
  
  merged <- merged %>%
    dplyr::mutate(
      sig_21 = sig_21,
      sig_24 = sig_24,
      sig_27 = sig_27,
      old_sig_count = sig_21 + sig_24 + sig_27
    )
  
  # Spearman rho across ages using available logFCs
  ages_present <- intersect(c(12,18,21,24,27),
                            as.integer(gsub("logFC_", "", grep("^logFC_", names(merged), value = TRUE))))
  logfc_cols <- paste0("logFC_", ages_present)
  
  merged$SpearmanRho <- apply(as.matrix(merged[, logfc_cols, drop = FALSE]), 1, function(v) {
    suppressWarnings(cor(ages_present, as.numeric(v), method = "spearman", use = "pairwise.complete.obs"))
  }) %>% as.numeric()
  
  merged <- merged %>%
    dplyr::mutate(age_responsive = (old_sig_count >= 2) | (SpearmanRho > rho_threshold))
  
  # subset to MusAge set
  out <- merged %>%
    dplyr::filter(GeneName %in% Inflammaging_symbols) %>%
    dplyr::mutate(muscle = muscle, sex = sex, stratum = paste(sex, muscle, sep = "_"))
  
  list(full = merged, musage = out)
}

# =========================
# Run all 6 strata
# =========================
res_F_Gas <- run_stratum(file_list_gas,    sample_groups_gas,    muscle = "Gastroc", sex = "F")
res_F_Sol <- run_stratum(file_list_sol,    sample_groups_sol,    muscle = "Soleus",  sex = "F")
res_F_TA  <- run_stratum(file_list_ta,     sample_groups_ta,     muscle = "TA",      sex = "F")

res_M_Gas <- run_stratum(file_list_m_gas,  sample_groups_m_gas,  muscle = "Gastroc", sex = "M")
res_M_Sol <- run_stratum(file_list_m_sol,  sample_groups_m_sol,  muscle = "Soleus",  sex = "M")
res_M_TA  <- run_stratum(file_list_m_ta,   sample_groups_m_ta,   muscle = "TA",      sex = "M")

# Save RData bundle (matches your note)
all_mus_age_dfs <- list(
  ta_age_df       = res_F_TA$full,
  gas_age_df      = res_F_Gas$full,
  sol_age_df      = res_F_Sol$full,
  male_ta_age_df  = res_M_TA$full,
  male_gas_age_df = res_M_Gas$full,
  male_sol_age_df = res_M_Sol$full
)
save(all_mus_age_dfs, file = "all_muscle_age_dfs_sep2025.RData")

# =========================
# Dot plot (Fig. 1j)
# =========================
dot_df <- bind_rows(
  res_F_Gas$musage, res_F_TA$musage, res_F_Sol$musage,
  res_M_Gas$musage, res_M_TA$musage, res_M_Sol$musage
) |>
  transmute(
    Gene = GeneName,
    Stratum = factor(stratum, levels = c("F_Gastroc","F_TA","F_Soleus","M_Gastroc","M_TA","M_Soleus")),
    age_responsive
  )

# keep your gene ordering style (top->bottom alphabetical unless you prefer a custom order)
dot_df$Gene <- factor(dot_df$Gene, levels = rev(sort(unique(Inflammaging_symbols))))

Fig1j_musage_dot <- ggplot(dot_df, aes(x = Stratum, y = Gene)) +
  geom_point(aes(fill = ifelse(age_responsive, "yes", "no")), shape = 21, size = 2) +
  scale_fill_manual(values = c(yes = "green3", no = "grey80"), guide = "none") +
  labs(x = NULL, y = NULL, title = "MusAge (52) – age-responsive genes across muscle × sex") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6)
  )

print(Fig1j_musage_dot)
ggsave("Fig1j_musage_dotplot_mouse.pdf", Fig1j_musage_dot, width = 5, height = 9, units = "in")

# CSV summary (TRUE/FALSE by gene×stratum)
wide_summary <- dot_df |>
  tidyr::pivot_wider(names_from = Stratum, values_from = age_responsive) |>
  arrange(Gene)
readr::write_csv(wide_summary, "Fig1j_age_responsive_summary_mouse.csv")


############################################################
# FIGURE 1k – Human aging re-analysis (PMID: 28273480; GSE97084)
############################################################
# 
# • Rows: 52 MusAge genes (mouse symbols) that met age-responsive criteria
# • Col 1: Up in aging (Old > Young)
#      red = FDR<0.05  | gray = FDR≥0.05
# • Col 2: Down in aging (Old < Young)
#      blue = FDR<0.05 | gray = FDR≥0.05
# • GSEA: MusAge gene set enrichment vs ranked Old–Young stats
############################################################
getwd()
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(edgeR)
  library(biomaRt)
  library(clusterProfiler)
  library(enrichplot)
  library(cowplot)
})

#////////////////////////////////////////////////////////////#
#-------Robinson MM 2017 Re-Analysis-------------------------#
#////////////////////////////////////////////////////////////#

# Tab-delimited, gene IDs as rownames
human_aging_counts_robinson_1 <- read.delim("GSE97084_GeneCount_raw.tsv",
                                            header = TRUE,
                                            check.names = FALSE)

# Inspect what you loaded
head(human_aging_counts_robinson_1)
dim(human_aging_counts_robinson_1)

# Tab-delimited, gene IDs as rownames
human_aging_counts_robinson_2 <- read.delim("GSE97084_GeneCount_raw_2.tsv",
                                            header = TRUE,
                                            check.names = FALSE)

# Inspect what you loaded
head(human_aging_counts_robinson_2)
dim(human_aging_counts_robinson_2)

colnames(human_aging_counts_robinson_2)

# Remove columns with empty names
human_aging_counts_robinson_2 <- human_aging_counts_robinson_2[, colnames(human_aging_counts_robinson_2) != ""]

# Now safely select
df1 <- human_aging_counts_robinson_1 %>%
  dplyr::select(-Chr, -Start, -Stop, -CodingLength)

df2 <- human_aging_counts_robinson_2 %>%
  dplyr::select(-Chr, -Start, -Stop, -CodingLength)

# Merge
human_aging_counts_robinson_merged <- full_join(df1, df2, by = c("GeneID","GeneName"))
# ---- 3. Check the result ----
dim(human_aging_counts_robinson_merged)
head(human_aging_counts_robinson_merged)

# ---- 1. Define keys ----
young_ids <- c("2B","3A","6A","10B","11A","13A","14B","15A","16A","17A","18A",
               "19A","23A","24A","25A","26A","27B","29A","30A","31B","32A","33A",
               "35A","36B","38A","39A","40A","41B")

old_ids <- c("1B","4A","5B","7B","8A","9A","20B","21B","22A","28A","34A","37B",
             "42A","43A","44A","45A","46B","47A","48B","49A","50A","51A","52A","53A")

# ---- 2. Extract sample ID from column names ----
# Keep GeneID and GeneName as identifiers
meta_cols <- c("GeneID", "GeneName")

sample_cols <- setdiff(colnames(human_aging_counts_robinson_merged), meta_cols)

# Extract the middle part (like "51A") from sample names
sample_ids <- sub("^s_([^-]+)-.*", "\\1", sample_cols)

# ---- 3. Identify which columns are young vs old ----
young_cols <- sample_cols[sample_ids %in% young_ids]
old_cols   <- sample_cols[sample_ids %in% old_ids]

# ---- 4. Keep only relevant columns ----
counts_filtered <- human_aging_counts_robinson_merged %>%
  dplyr::select(all_of(c(meta_cols, young_cols, old_cols)))

# ---- 5. Optional: create group assignment dataframe ----
sample_group <- data.frame(
  sample = c(young_cols, old_cols),
  group  = c(rep("Young", length(young_cols)), rep("Old", length(old_cols)))
)

counts_filtered

# Get current column names
cn <- colnames(counts_filtered)

# Simplify sample names (but keep GeneID and GeneName)
cn_new <- cn
sample_mask <- !(cn %in% c("GeneID", "GeneName"))
cn_new[sample_mask] <- sub("^s_([^-]+)-.*", "\\1", cn[sample_mask])

# Apply new column names
colnames(counts_filtered) <- cn_new

cn_new
counts_filtered

# Create metadata (Young vs Old)
sample_group <- data.frame(
  sample = cn_new[!(cn_new %in% c("GeneID","GeneName"))],
  group  = ifelse(cn_new[!(cn_new %in% c("GeneID","GeneName"))] %in% young_ids, "Young", "Old")
)

# Build count matrix
count_matrix <- counts_filtered %>%
  dplyr::select(-GeneName) %>%
  column_to_rownames("GeneID")

counts_filtered_unique <- counts_filtered %>%
  group_by(GeneID) %>%
  summarise(across(-GeneName, sum), .groups = "drop")  # sum counts for duplicates

count_matrix <- counts_filtered_unique %>%
  column_to_rownames("GeneID")

count_matrix[is.na(count_matrix)] <- 0

count_matrix

# edgeR pipeline
library(edgeR)
# Ensure group is a factor and relevel it BEFORE design
sample_group$group <- relevel(factor(sample_group$group), ref = "Young")

# DGEList
dge <- DGEList(counts = count_matrix, group = sample_group$group)
dge <- calcNormFactors(dge)

# Design matrix with Young as reference
design <- model.matrix(~ group, data = sample_group)

# Estimate dispersion and fit model
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

# LRT: coef=2 now represents Old vs Young (positive logFC = Old > Young)
lrt <- glmLRT(fit, coef = 2)

topTags(lrt)




# Extract all results as dataframe
deg_results <- topTags(lrt, n = Inf)$table
deg_results <- tibble::rownames_to_column(deg_results, var = "ensembl_gene_id")

head(deg_results)

# Connect to Ensembl archive
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org")

# Map human Ensembl gene IDs to mouse gene symbols
orthologs <- getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = deg_results$ensembl_gene_id,
                    mart = human,
                    attributesL = c("ensembl_gene_id","mgi_symbol"),
                    martL = mouse,
                    uniqueRows = TRUE)

colnames(orthologs) <- c("ensembl_gene_id", "hgnc_symbol",
                         "ensembl_gene_id_mouse", "mgi_symbol")

# Merge with DEG results
deg_mouse_robinson <- merge(deg_results, orthologs, by = "ensembl_gene_id")

head(deg_mouse_robinson)

deg_mouse_robinson %>%
  pull(mgi_symbol)

dim(deg_mouse_robinson)

deg_mouse_robinson %>%
  filter(FDR <0.05 & logFC > 0)
#--------------------------------------------#

#============================================#
#-----Robinson Dotplot-----------------------#
#============================================#
validated_inflamm_geneset.vec

# Step 1: Prepare master gene list
gene_df <- tibble(mgi_symbol = sort(validated_inflamm_geneset.vec))

# Step 2: Reduce deg_mouse_robinson to unique best hits
deg_filtered <- deg_mouse_robinson %>%
  filter(mgi_symbol %in% validated_inflamm_geneset.vec) %>%
  group_by(mgi_symbol) %>%
  slice_min(order_by = FDR, n = 1) %>%
  ungroup()

# Step 3: Join to get all 51 genes with data (or NA if missing)
plot_df <- gene_df %>%
  left_join(deg_filtered, by = "mgi_symbol") %>%
  mutate(
    column = case_when(
      is.na(logFC) ~ NA_character_,
      logFC < 0 ~ "Negative logFC",
      logFC > 0 ~ "Positive logFC",
      TRUE ~ NA_character_
    ),
    color = case_when(
      is.na(logFC) ~ NA_character_,
      logFC < 0 & PValue < 0.05 ~ "blue",
      logFC < 0 & PValue >= 0.05 ~ "gray",
      logFC > 0 & PValue < 0.05 ~ "red",
      logFC > 0 & PValue >= 0.05 ~ "gray",
      TRUE ~ NA_character_
    ),
    size = ifelse(color == "gray", 2, 3)
  )

# Flip the x-axis: Positive on the left (0.8), Negative on the right (1.2)
plot_df_filtered <- plot_df %>%
  mutate(
    x = ifelse(column == "Positive logFC", 0.6, 0.7)
  )

# Updated dotplot
robinson_dotplot <- ggplot(plot_df_filtered, aes(x = x, y = factor(mgi_symbol, levels = rev(sort(mgi_symbol))))) +
  geom_point(aes(color = color, size = size), na.rm = TRUE) +
  scale_color_identity() +
  scale_size_identity() +
  scale_x_continuous(
    breaks = c(0.8, 1.2),
    limits = c(0.6, 1.4)  # tight bounds
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 10)
  )


# =========================================================
# B) GSEA of MusAge set vs ranked Old–Young statistic
# =========================================================
# Rank by signed strength. Two good options:
# 1) -log10(P) * sign(logFC) * |logFC|  (your original)
# 2) sign(logFC) * sqrt(LR)             (often cleaner)
# We'll use (2); switch to (1) by swapping the mutate() line.

geneList <- deg_mouse_robinson %>%
  filter(!is.na(mgi_symbol)) %>%
  mutate(stat = sign(logFC) * sqrt(pmax(LR, 0))) %>%     # <-- ranking metric
  arrange(desc(stat)) %>%
  distinct(mgi_symbol, .keep_all = TRUE) %>%
  { setNames(.$stat, .$mgi_symbol) }

# clusterProfiler GSEA (for NES/Pvals)
gsea_inflamm_human_robinson <- GSEA(
  geneList = geneList,
  TERM2GENE = data.frame(term = "Inflamm", gene = validated_inflamm_geneset.vec),
  pvalueCutoff = 1,
  verbose = FALSE
)

# quick table
print(gsea_inflamm_human_robinson@result %>% dplyr::select(ID, NES, pvalue, p.adjust))

# ----- Custom 3-panel enrichment figure (top: ES; mid: ticks; bottom: stat bar) -----
stats <- sort(geneList, decreasing = TRUE)
pathway <- intersect(names(stats), validated_inflamm_geneset.vec)

N  <- length(stats)
hit <- names(stats) %in% pathway
Nh <- sum(hit)
Nm <- N - Nh
P  <- sum(abs(stats[hit]))   # weighted p = 1

runningES <- cumsum(ifelse(hit, abs(stats) / P, -1 / Nm))

df_top <- data.frame(Rank = seq_along(stats), ES = runningES)
df_hits <- data.frame(Rank = which(hit))
df_bot <- data.frame(Rank = seq_along(stats), Stat = as.numeric(stats))

p_top <- ggplot(df_top, aes(Rank, ES)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(data = df_top[which.max(abs(df_top$ES)), , drop=FALSE],
             aes(Rank, ES), color = "red", size = 2.5) +
  labs(y = "Enrichment Score (ES)", x = NULL,
       title = "GSEA: MusAge set in Old vs Young (human)") +
  theme_minimal(base_size = 11)

p_mid <- ggplot(df_hits, aes(Rank, 1)) +
  geom_segment(aes(xend = Rank, yend = 0), color = "black") +
  labs(y = NULL, x = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = -14, b = -14, l = 5, r = 5))

p_bot <- ggplot(df_bot, aes(Rank, Stat, fill = Stat)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Rank score") +
  labs(x = "Gene Rank", y = "Ranked Statistic") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

Fig1k_gsea <- plot_grid(p_top, p_mid, p_bot, ncol = 1, align = "v",
                        rel_heights = c(3, 0.5, 2))

print(Fig1k_gsea)
# ggsave("Fig1k_human_musage_GSEA.pdf", Fig1k_gsea, width = 7, height = 5, units = "in")

############################################################
# End of Figure 1 consolidated script
############################################################


############################################################
# FIGURE 2a — GO dot plot (ORA) for MusAge 52 (mouse, ENTREZ)
# Goal:
#   • Run enrichGO on the MusAge 52 set against the aging background
#   • Simplify redundant terms (semantic similarity, Wang)
#   • Plot a dot plot:
#       - Rows: significant BP terms (FDR < 0.05), Count > 10
#       - Dot color: -log10(FDR), blue (~4) → red (~8)
#       - Dot size: term hit Count
#       - X-axis: Gene ratio (k / 52)
#
# Inputs expected in the session:
#   • validated_inflamm_geneset.vec  (MusAge 52, mouse SYMBOLs)
#   • oldsedveh_v_yngsedveh_genes_12jun (from Fig. 1g; has ENTREZID column)
#
# Notes:
#   • If simplify() errors about semantic data, install/load GOSemSim.
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr); library(stringr); library(tibble)
  library(ggplot2); library(forcats); library(scales)
  # library(GOSemSim) # <- uncomment if simplify() complains
})

## 1) Ensure MusAge 52 as ENTREZ IDs
if (!exists("validated_inflamm_geneset.vec")) {
  # fallback: read from CSV in project root if missing
  validated_inflamm_geneset.vec <- readr::read_csv("MusAge_geneset.csv", show_col_types = FALSE) |>
    dplyr::pull(Gene) |> as.character() |> stringr::str_trim() |> unique()
}
# SYMBOL -> ENTREZ (reuse prior mapping if present)
if (exists("inflammaging_gene_vec")) {
  validated_inflamm_geneset.entrez <- unname(na.omit(unique(inflammaging_gene_vec)))
} else {
  validated_inflamm_geneset.entrez <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys     = validated_inflamm_geneset.vec,
    keytype  = "SYMBOL",
    column   = "ENTREZID",
    multiVals = "first"
  ) |> unname() |> unique() |> na.omit() |> as.character()
}

## 2) Background universe (aging comparison from Fig. 1g)
if (exists("oldsedveh_v_yngsedveh_genes_12jun")) {
  aging_bckgrnd.entrez <- oldsedveh_v_yngsedveh_genes_12jun %>% dplyr::pull(ENTREZID) %>% unique() %>% na.omit() %>% as.character()
} else if (exists("oldsedveh_v_yngsedveh_genes")) {
  aging_bckgrnd.entrez <- oldsedveh_v_yngsedveh_genes %>% dplyr::pull(ENTREZID) %>% unique() %>% na.omit() %>% as.character()
} else {
  stop("Could not find the aging background object. Make sure Fig. 1g code has run to create the gene table with an ENTREZID column.")
}

## 3) Over-representation analysis (GO BP)
go_enrich <- enrichGO(
  gene          = validated_inflamm_geneset.entrez,
  universe      = aging_bckgrnd.entrez,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

# Early peek (optional)
# head(as.data.frame(go_enrich), 20)

## 4) Reduce redundancy (semantic similarity)
go_enrich_simplified <- simplify(
  go_enrich,
  cutoff     = 0.5,     # 0.4–0.7 typical; 0.5 balances collapse vs. specificity
  by         = "p.adjust",
  select_fun = min,
  measure    = "Wang"
)

go_simpl_df <- as.data.frame(go_enrich_simplified)

## 5) Display subset: Count > 10 and FDR < 0.05 (tune to taste)
inflamm_go_subset <- go_simpl_df %>%
  filter(Count > 10, p.adjust < 0.05) %>%
  arrange(GeneRatio) %>%
  head(20)

# If nothing passes Count>10, relax to >=10:
if (nrow(inflamm_go_subset) == 0) {
  inflamm_go_subset <- go_simpl_df %>%
    filter(Count >= 10, p.adjust < 0.05) %>%
    arrange(GeneRatio) %>% head(20)
}

## 6) Build plotting data
shorten_go <- function(x) x %>%
  str_replace_all("cysteine-type endopeptidase", "caspase") %>%
  str_replace_all("apoptotic process", "apoptosis") %>%
  str_replace_all("\\binvolved in\\b", "in") %>%
  str_replace_all("\\bregulation of\\b", "Reg.") %>%
  str_squish()

df <- inflamm_go_subset %>%
  mutate(
    # GeneRatio in enrichGO is typically "k/|gene|"; convert to numeric fraction
    GeneRatio_num = sapply(strsplit(as.character(GeneRatio), "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2])),
    neglogFDR     = -log10(p.adjust),
    Term          = stringr::str_wrap(shorten_go(Description), width = 36)
  ) %>%
  mutate(Term = forcats::fct_reorder(Term, GeneRatio_num))

## 7) Dot plot
# Map -log10(FDR) with blue (~4) → red (~8); squish outside range
left_pad <- diff(range(df$GeneRatio_num, na.rm = TRUE)) * 0.06

Fig2a_GO_dot <- ggplot(df, aes(x = GeneRatio_num, y = Term)) +
  geom_point(aes(size = Count, color = neglogFDR), alpha = 0.9) +
  scale_size_area(max_size = 10, name = "Hit genes") +
  scale_color_gradientn(
    colours = c("#2B8CBE", "#F46D43"),
    limits  = c(4, 8),
    oob     = scales::squish,
    name    = expression(-log[10]~FDR)
  ) +
  scale_x_continuous(
    name   = "Gene ratio (k / 52)",
    labels = scales::percent_format(accuracy = 0.1),
    limits = c(min(df$GeneRatio_num, na.rm = TRUE) - left_pad, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(y = NULL, title = "GO Biological Process over-representation (MusAge 52)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y  = element_text(size = 9),
    legend.box   = "vertical",
    plot.margin  = margin(8, 16, 8, 20)
  )

print(Fig2a_GO_dot)
#ggsave("Fig2a_GO_dotplot_MusAge52.pdf", Fig2a_GO_dot, width = 6.5, height = 6.5, units = "in")

############################################################
# FIGURE 2b — Dot matrix of cell-type expression (Human Protein Atlas)
# Goal:
#   • 52 MusAge genes (mouse symbols) × 81 HPA cell types
#   • Columns grouped by category (Immune / Muscle / …), colored by category
#   • Rows = genes; per-row (gene) z-score across all 81 cell types
#   • Dot size = positive z-score (z+); zeros hidden
#
# Inputs expected in session:
#   • validated_inflamm_geneset.vec (MusAge 52 mouse symbols; already loaded earlier)
#   • rna_single_cell_type.tsv (HPA single-cell summary; columns Gene.name, Cell.type, nTPM)
############################################################

suppressPackageStartupMessages({
  library(biomaRt); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(pheatmap); library(RColorBrewer); library(scales)
})

# 0) Make sure the MusAge set exists (alias the UPPER vec if that’s what you have)
if (!exists("validated_inflamm_geneset.vec") && exists("validated_inflamm_geneset_upper.vec")) {
  validated_inflamm_geneset.vec <- validated_inflamm_geneset_upper.vec
}

# 1) Load HPA single-cell table
sc_data <- read.delim("rna_single_cell_type.tsv", check.names = FALSE)

# Rename to match the rest of your pipeline
sc_data <- dplyr::rename(
  sc_data,
  Gene.name = `Gene name`,
  Cell.type = `Cell type`
)

# sanity check
stopifnot(all(c("Gene.name","Cell.type","nTPM") %in% names(sc_data)))
head(sc_data)


# 2) Map human symbols -> mouse symbols (Dec 2021 archive for consistency)
human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                          host = "https://dec2021.archive.ensembl.org")
mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                          host = "https://dec2021.archive.ensembl.org")

orthologs <- biomaRt::getLDS(
  attributes = c("hgnc_symbol"),
  filters    = "hgnc_symbol",
  values     = unique(sc_data$Gene.name),
  mart       = human,
  attributesL= c("mgi_symbol"),
  martL      = mouse,
  uniqueRows = TRUE
)
colnames(orthologs) <- c("Gene.name","Mouse.gene")

# 3) Merge mouse symbol and keep MusAge overlap
sc_data_mouse <- sc_data %>% left_join(orthologs, by = "Gene.name")
sc_MusAge <- sc_data_mouse %>% filter(Mouse.gene %in% validated_inflamm_geneset.vec)

# 4) Wide matrix: genes × cell types, log2(nTPM+1)
sc_wide <- sc_MusAge %>%
  group_by(Mouse.gene, Cell.type) %>%
  summarise(nTPM = mean(nTPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cell.type, values_from = nTPM)

sc_matrix <- sc_wide %>% as.data.frame()
rownames(sc_matrix) <- sc_matrix$Mouse.gene
sc_matrix$Mouse.gene <- NULL
sc_matrix <- as.matrix(sc_matrix)
sc_matrix[is.na(sc_matrix)] <- 0
sc_matrix_log <- log2(sc_matrix + 1)

# 5) Category map for the 81 HPA cell types
celltype_categories <- c(
  # Neural
  "Excitatory neurons"="Neural","Inhibitory neurons"="Neural","Astrocytes"="Neural",
  "Microglial cells"="Neural","Oligodendrocytes"="Neural","Oligodendrocyte precursor cells"="Neural",
  "Schwann cells"="Neural","Bipolar cells"="Neural","Horizontal cells"="Neural",
  "Muller glia cells"="Neural","Rod photoreceptor cells"="Neural","Cone photoreceptor cells"="Neural",
  # Muscle
  "Skeletal myocytes"="Muscle","Cardiomyocytes"="Muscle","Smooth muscle cells"="Muscle",
  "Peritubular cells"="Muscle","Breast myoepithelial cells"="Muscle",
  # Immune
  "T-cells"="Immune","B-cells"="Immune","NK-cells"="Immune","Dendritic cells"="Immune",
  "Macrophages"="Immune","Monocytes"="Immune","Granulocytes"="Immune","Plasma cells"="Immune",
  "Langerhans cells"="Immune","Hofbauer cells"="Immune",
  # Epithelial/Secretory
  "Basal keratinocytes"="Epithelial","Suprabasal keratinocytes"="Epithelial",
  "Squamous epithelial cells"="Epithelial","Basal squamous epithelial cells"="Epithelial",
  "Glandular and luminal cells"="Epithelial","Serous glandular cells"="Epithelial",
  "Mucus glandular cells"="Epithelial","Exocrine glandular cells"="Epithelial",
  "Ductal cells"="Epithelial","Salivary duct cells"="Epithelial","Secretory cells"="Epithelial",
  "Ciliated cells"="Epithelial","Club cells"="Epithelial","Ionocytes"="Epithelial",
  "Breast glandular cells"="Epithelial","Prostatic glandular cells"="Epithelial","Basal prostatic cells"="Epithelial",
  # Endothelial/Stromal
  "Endothelial cells"="Endothelial/Stromal","Lymphatic endothelial cells"="Endothelial/Stromal",
  "Mesothelial cells"="Endothelial/Stromal","Fibroblasts"="Endothelial/Stromal",
  "Endometrial stromal cells"="Endothelial/Stromal","Adipocytes"="Endothelial/Stromal",
  # Digestive
  "Pancreatic endocrine cells"="Digestive","Enteroendocrine cells"="Digestive",
  "Gastric mucus-secreting cells"="Digestive","Intestinal goblet cells"="Digestive",
  "Distal enterocytes"="Digestive","Proximal enterocytes"="Digestive",
  "Paneth cells"="Digestive","Hepatocytes"="Digestive",
  # Reproductive/Germline
  "Oocytes"="Reproductive","Granulosa cells"="Reproductive","Ovarian stromal cells"="Reproductive",
  "Sertoli cells"="Reproductive","Leydig cells"="Reproductive","Spermatogonia"="Reproductive",
  "Spermatocytes"="Reproductive","Early spermatids"="Reproductive","Late spermatids"="Reproductive",
  # Other / Specialized
  "Kupffer cells"="Immune",                         # <- minor relabel for grouping
  "Cholangiocytes"="Epithelial",
  "Collecting duct cells"="Epithelial","Distal tubular cells"="Epithelial","Proximal tubular cells"="Epithelial",
  "Cytotrophoblasts"="Reproductive","Extravillous trophoblasts"="Reproductive","Syncytiotrophoblasts"="Reproductive",
  "Undifferentiated cells"="Epithelial"
)

# Keep only cell types that are present; assign "Other" to any leftover columns
valid_cols <- intersect(names(celltype_categories), colnames(sc_matrix_log))
missing_cols <- setdiff(colnames(sc_matrix_log), valid_cols)
if (length(missing_cols)) {
  celltype_categories[missing_cols] <- "Other"
  valid_cols <- colnames(sc_matrix_log)
}

# Order columns by category (alphabetical categories) then within-category by correlation clustering
z_mat <- t(scale(t(sc_matrix_log)))                      # row z-scores
z_mat <- z_mat[, valid_cols, drop = FALSE]

cat_of <- celltype_categories[colnames(z_mat)]
cat_levels <- unique(cat_of[order(cat_of)])              # stable category order

cluster_within <- function(cols) {
  if (length(cols) <= 1) return(cols)
  d <- as.dist(1 - cor(z_mat[, cols, drop = FALSE], method = "spearman", use = "pairwise.complete.obs"))
  h <- hclust(d, method = "average")
  cols[h$order]
}
col_order <- unlist(lapply(cat_levels, function(cat) {
  cluster_within(names(cat_of)[cat_of == cat])
}), use.names = FALSE)

z_mat <- z_mat[, col_order, drop = FALSE]
cat_of <- factor(celltype_categories[col_order], levels = cat_levels)

# Build long df for dot plot (size = positive z; hide negatives)
z_long <- as.data.frame(z_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "z") %>%
  mutate(Category = factor(celltype_categories[CellType], levels = levels(cat_of)),
         z_pos = pmax(z, 0))

# Keep only positive signal (optional but matches “size = z-score” intent)
z_long <- z_long %>% filter(z_pos > 0)

# Category palette
cat_cols <- setNames(brewer.pal(max(3, length(levels(cat_of))), "Set2")[seq_along(levels(cat_of))],
                     levels(cat_of))

# Order axes
z_long$CellType <- factor(z_long$CellType, levels = col_order)
# Order genes by hierarchical clustering on z_mat rows
row_order <- {
  d <- dist(z_mat, method = "euclidean")
  h <- hclust(d, method = "average")
  rownames(z_mat)[h$order]
}
z_long$Gene <- factor(z_long$Gene, levels = row_order)

# Plot
Fig2b_dot <- ggplot(z_long, aes(x = CellType, y = Gene)) +
  geom_point(aes(size = z_pos, fill = Category),
             shape = 21, color = "grey25", stroke = 0.3) +
  scale_size_area(name = "Row z-score (≥0)", max_size = 6) +
  scale_fill_manual(values = cat_cols, name = "Cell type category") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  labs(x = NULL, y = NULL,
       title = "Cell-type expression of MusAge genes (HPA single-cell)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_line(size = 0.2),
    legend.box = "vertical"
  )

print(Fig2b_dot)
#ggsave("Fig2b_dotmatrix_MusAge_HPA.pdf", Fig2b_dot, width = 9.5, height = 8, units = "in")

# (Optional) export top-3 cell types per gene table (from z_mat, not raw nTPM)
top3_by_gene <- as.data.frame(z_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "z") %>%
  group_by(Gene) %>%
  arrange(desc(z), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  mutate(Rank = row_number()) %>%
  pivot_wider(names_from = Rank, values_from = c(CellType, z)) %>%
  ungroup()
readr::write_csv(top3_by_gene, "Fig2b_top3_celltypes_per_gene.csv")

## ============================================================
## Fig. 2c–2f consolidated (uses your helpers & objects)
## ============================================================

# (Make sure these exist from your earlier block)
# - read_clean_xlsx()
# - validated_inflamm_geneset.vec    # 52 SYMBOLs
# - yng_sed_samples, old_sed_samples, old_pwr_samples
# - libraries edgeR, pheatmap, ggplot2, clusterProfiler, cowplot, org.Mm.eg.db, AnnotationDbi

suppressPackageStartupMessages({
  library(cowplot)
})

# -----------------------------
# 0) Load *_GENE_* result tables
# -----------------------------
oldsedveh_v_yngsed_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col   = "OLD_SED_VEH-YNG_SED_VEH_FDR"
)

oldpwrveh_v_oldsed_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_PWR_VEH-OLD_SED_VEH_logFC",
  FDR_col   = "OLD_PWR_VEH-OLD_SED_VEH_FDR"
)

# Shared gene order for all panels (from OLD_PWR vs OLD_SED)
r01_heatmap_order <- oldpwrveh_v_oldsed_genes %>%
  dplyr::filter(Symbol %in% validated_inflamm_geneset.vec) %>%
  dplyr::arrange(logFC) %>%
  dplyr::pull(Symbol)

# --------------------------------------------------
# Helper: Counts -> logCPM (rows SYMBOL, cols samples)
# --------------------------------------------------
counts_to_logCPM_by_symbol <- function(counts_xlsx, sample_order, group_names) {
  # counts_xlsx: path to *_Counts_*.xlsx
  # sample_order: vector of sample IDs in desired order
  # group_names:  vector of group labels (same length as sample_order)
  stopifnot(length(sample_order) == length(group_names))
  
  # read & map to ENTREZID
  df <- readxl::read_xlsx(counts_xlsx)
  colnames(df)[1] <- "Ensembl"
  df <- df %>%
    dplyr::mutate(Ensembl_noDec = stringr::str_remove(Ensembl, "\\..+"),
                  ENTREZID = AnnotationDbi::mapIds(org.Mm.eg.db,
                                                   keys = Ensembl_noDec,
                                                   column = "ENTREZID",
                                                   keytype = "ENSEMBL",
                                                   multiVals = "first")) %>%
    dplyr::select(-Ensembl) %>%
    tidyr::drop_na(ENTREZID)
  
  # tolerate missing samples: intersect with present column names
  present <- intersect(sample_order, colnames(df))
  if (length(present) < 2) stop("Too few matching samples found in counts file.")
  grp_present <- group_names[match(present, sample_order)]
  
  # counts matrix
  cm <- df %>% dplyr::select(all_of(present)) %>% as.matrix()
  rownames(cm) <- df$ENTREZID
  
  # edgeR TMM -> logCPM
  group <- factor(grp_present, levels = unique(grp_present))
  dge <- edgeR::DGEList(counts = cm, group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  
  # map ENTREZID -> SYMBOL and collapse duplicates by mean
  sym <- AnnotationDbi::mapIds(org.Mm.eg.db,
                               keys = rownames(logCPM),
                               column = "SYMBOL",
                               keytype = "ENTREZID",
                               multiVals = "first")
  logCPM_df <- as.data.frame(logCPM) %>%
    tibble::rownames_to_column("ENTREZID") %>%
    dplyr::mutate(Symbol = sym) %>%
    tidyr::drop_na(Symbol) %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), mean), .groups = "drop") %>%
    tibble::column_to_rownames("Symbol")
  logCPM_mat <- as.matrix(logCPM_df)
  
  list(logCPM = logCPM_mat[, present, drop = FALSE], samples = present, groups = grp_present)
}

# ------------------------------------
# 1) FIG. 2c – Heatmap (YNG_SED vs OLD_SED)
# ------------------------------------
group_vec_2c <- c(rep("YNG_SED", length(yng_sed_samples)),
                  rep("OLD_SED", length(old_sed_samples)))
order_vec_2c <- c(yng_sed_samples, old_sed_samples)

res_2c <- counts_to_logCPM_by_symbol(
  counts_xlsx = "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  sample_order = order_vec_2c,
  group_names  = group_vec_2c
)

# subset to the 52 symbols and enforce shared row order
Inflamm_logCPM_matrix_2c <- res_2c$logCPM[rownames(res_2c$logCPM) %in% validated_inflamm_geneset.vec, , drop = FALSE]
row_order_2c <- intersect(r01_heatmap_order, rownames(Inflamm_logCPM_matrix_2c))
Inflamm_logCPM_matrix_2c <- Inflamm_logCPM_matrix_2c[row_order_2c, res_2c$samples, drop = FALSE]

# column annotations (blue vs gray)
ann_col_2c <- data.frame(Group = factor(res_2c$groups, levels = c("YNG_SED","OLD_SED")))
rownames(ann_col_2c) <- res_2c$samples
ann_colors_2c <- list(Group = c(YNG_SED = "#4C78A8", OLD_SED = "#7F7F7F"))

Fig2c_heat <- pheatmap::pheatmap(
  Inflamm_logCPM_matrix_2c,
  cluster_rows = FALSE, cluster_cols = FALSE,
  scale = "row",
  show_rownames = TRUE, show_colnames = TRUE,
  row_names_side = "left",
  annotation_col = ann_col_2c,
  annotation_colors = ann_colors_2c,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

pdf("Fig2c_heatmap_YNGvsOLD.pdf", width = 4.25, height = 5)
print(Fig2c_heat); dev.off()

# Bar plot (OLD_SED − YNG_SED), black bars>0
logfc_bar_data_2c <- oldsedveh_v_yngsed_genes %>%
  dplyr::filter(Symbol %in% row_order_2c) %>%
  dplyr::mutate(Symbol = factor(Symbol, levels = rev(row_order_2c)))

Fig2c_bar <- ggplot2::ggplot(logfc_bar_data_2c, aes(x = logFC, y = Symbol)) +
  ggplot2::geom_col(fill = "black") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "log2 Fold Change (OLD_SED vs YNG_SED)", y = NULL) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                 axis.text.x = ggplot2::element_text(size = 8),
                 panel.grid.major.y = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank())

pdf("Fig2c_barplot_YNGvsOLD.pdf", width = 4.25, height = 5)
print(Fig2c_bar) 
dev.off()

# ------------------------------------
# 2) FIG. 2d – Targeted GSEA (OLD_SED vs YNG_SED)
# ------------------------------------
# ------- Targeted GSEA using logFC as the ranking statistic -------

## 0) Use ONE MusAge ENTREZ vector everywhere
musage_entrez <- unique(na.omit(inflammaging_gene_vec))  # already ENTREZ IDs

## 1) After create_vec_from_df(), force sort + unique names
oldsedveh_v_yngsedveh_12jun.vec <- oldsedveh_v_yngsedveh_12jun.vec[!is.na(oldsedveh_v_yngsedveh_12jun.vec)]
oldsedveh_v_yngsedveh_12jun.vec <- sort(oldsedveh_v_yngsedveh_12jun.vec, decreasing = TRUE)
stopifnot(!any(duplicated(names(oldsedveh_v_yngsedveh_12jun.vec))))

## 2) Sanity: how many MusAge genes are actually ranked?
cat("MusAge hits in rank list:",
    sum(names(oldsedveh_v_yngsedveh_12jun.vec) %in% musage_entrez),
    "of", length(musage_entrez), "\n")

## 3) Run GSEA (same as you did)
set.seed(1)  # reproducible NES/p-values
gsea_inflamm_oldsedveh_v_yngsedveh <- clusterProfiler::GSEA(
  geneList     = oldsedveh_v_yngsedveh_12jun.vec,            # ENTREZ→logFC
  TERM2GENE    = data.frame(term = "Inflamm", gene = musage_entrez),
  pvalueCutoff = 1,
  minGSSize    = 1,  # keep your 52-set even if a few are missing
  maxGSSize    = 10000,
  verbose      = FALSE
)

## Optional: draw the classic curve with enrichplot (matches your “built-in” look)
enrichplot::gseaplot2(gsea_inflamm_oldsedveh_v_yngsedveh, geneSetID = 1,
                      title = "MusAge ~ OLD SED vs YNG SED")


## ======================================================
## Fig. 2e — Heatmap + bar (OLD PWR vs OLD SED)
## ======================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(pheatmap); library(AnnotationDbi); library(org.Mm.eg.db)
  library(cowplot); library(stringr)
})

old_pwr_samples <- c("T_08", "T_10", "T_11", "T_41", "T_43", "T_45", "T_46")

# 1) Load *_GENE_* table and ranked vector for OLD_PWR vs OLD_SED
oldpwrveh_v_oldsed_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_PWR_VEH-OLD_SED_VEH_logFC",
  FDR_col   = "OLD_PWR_VEH-OLD_SED_VEH_FDR"
)

create_vec_from_df(oldpwrveh_v_oldsed_genes, "oldpwrveh_v_oldsed")  # creates oldpwrveh_v_oldsed.vec

# 2) Ensure MusAge 52 ENTREZ vector
if (!exists("validated_inflamm_geneset.vec")) {
  validated_inflamm_geneset.vec <- readr::read_csv("MusAge_geneset.csv", show_col_types = FALSE) |>
    dplyr::pull(Gene) |> as.character() |> stringr::str_trim() |> unique()
}
if (!exists("inflammaging_gene_vec")) {
  inflammaging_gene_vec <- AnnotationDbi::mapIds(
    org.Mm.eg.db, keys = validated_inflamm_geneset.vec,
    keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"
  ) |> unname() |> unique()
}
musage_entrez <- unique(na.omit(inflammaging_gene_vec))

# 3) Counts → logCPM for Set02 (OLD_PWR vs OLD_SED), subset to 52 genes
stopifnot(exists("old_pwr_samples"), exists("old_sed_samples"))
order_vec_2e <- c(old_pwr_samples, old_sed_samples)
group_vec_2e <- c(rep("OLD_PWR", length(old_pwr_samples)),
                  rep("OLD_SED", length(old_sed_samples)))

res_2e <- counts_to_logCPM_by_symbol(
  counts_xlsx  = "20250219_M007853_Set02_edgeRglm_Counts_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
  sample_order = order_vec_2e,
  group_names  = group_vec_2e
)

logCPM_2e <- res_2e$logCPM
Inflamm_logCPM_2e <- logCPM_2e[rownames(logCPM_2e) %in% validated_inflamm_geneset.vec, res_2e$samples, drop = FALSE]

# 4) Shared row order for heatmap = by logFC (OLD_PWR−OLD_SED), low→high
row_order_2e <- oldpwrveh_v_oldsed_genes %>%
  filter(Symbol %in% rownames(Inflamm_logCPM_2e)) %>%
  arrange(logFC) %>% pull(Symbol)

Inflamm_logCPM_2e <- Inflamm_logCPM_2e[row_order_2e, , drop = FALSE]

# 5) Column annotations (green = OLD PWR left; grey = OLD SED right)
ann_col_2e <- data.frame(Group = factor(res_2e$groups, levels = c("OLD_PWR","OLD_SED")))
rownames(ann_col_2e) <- res_2e$samples
ann_colors_2e <- list(Group = c(OLD_PWR = "#2CA02C", OLD_SED = "#7F7F7F"))

# 6) Heatmap (row z-scored)
Fig2e_heat <- pheatmap::pheatmap(
  Inflamm_logCPM_2e,
  cluster_rows = FALSE, cluster_cols = FALSE,
  scale = "row",
  show_rownames = TRUE, show_colnames = TRUE,
  row_names_side = "left",
  annotation_col = ann_col_2e,
  annotation_colors = ann_colors_2e,
  color = colorRampPalette(c("blue","white","red"))(50),
  silent = TRUE
)

# 7) Bar plot of log2FC (OLD_PWR vs OLD_SED); red = logFC<0 (down in OLD_PWR), black = logFC>0
bar_df_2e <- oldpwrveh_v_oldsed_genes %>%
  filter(Symbol %in% row_order_2e) %>%
  mutate(Symbol = factor(Symbol, levels = rev(row_order_2e)),
         dir = ifelse(logFC < 0, "down_in_OLD_PWR", "up_in_OLD_PWR"))

Fig2e_bar <- ggplot(bar_df_2e, aes(x = logFC, y = Symbol, fill = dir)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(down_in_OLD_PWR = "red3", up_in_OLD_PWR = "black"), guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "log2 Fold Change (OLD PWR − OLD SED)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

print(Fig2e_bar)
print(Fig2e_heat)

# 8) Save as separate PDFs and an optional combined panel
pdf("Fig2e_heatmap_OLDPWR_vs_OLDSED.pdf", width = 4.5, height = 5.0); print(Fig2e_heat); dev.off()
ggsave("Fig2e_bar_OLDPWR_vs_OLDSED.pdf", Fig2e_bar, width = 3.2, height = 5.0, units = "in")

# Combine (heatmap left, bars right)
heat_grob <- Fig2e_heat[[4]]  # pheatmap grob
combined_2e <- cowplot::plot_grid(
  ggplotify::as.ggplot(heat_grob),
  Fig2e_bar,
  ncol = 2, rel_widths = c(1.3, 0.9)
)
ggsave("Fig2e_combined_heatmap_bar.pdf", combined_2e, width = 8.2, height = 5.0, units = "in")


## ===============================
## Fig. 2f — Targeted GSEA (MusAge vs OLD PWR − OLD SED)
## ===============================
suppressPackageStartupMessages({
  library(clusterProfiler); library(enrichplot); library(cowplot)
})

# 1) Make sure ranked list is strictly decreasing, unique names (ENTREZ -> logFC)
oldpwrveh_v_oldsed.vec <- oldpwrveh_v_oldsed.vec[!is.na(oldpwrveh_v_oldsed.vec)]
oldpwrveh_v_oldsed.vec <- sort(oldpwrveh_v_oldsed.vec, decreasing = TRUE)
stopifnot(!any(duplicated(names(oldpwrveh_v_oldsed.vec))))

# 2) Run targeted GSEA with MusAge ENTREZ IDs
set.seed(1)
gsea_inflamm_oldpwr_v_oldsed <- clusterProfiler::GSEA(
  geneList     = oldpwrveh_v_oldsed.vec,
  TERM2GENE    = data.frame(term = "MusAge", gene = musage_entrez),
  pvalueCutoff = 1, minGSSize = 1, maxGSSize = 10000,
  verbose      = FALSE
)

# Quick table (NES/p-values)
print(gsea_inflamm_oldpwr_v_oldsed@result %>% dplyr::select(ID, NES, pvalue, p.adjust))

# 3a) Built-in enriched curve (matches your previous “worked fine” plot)
pdf("Fig2f_gseaplot2_MusAge_OLDPWR_vs_OLDSED.pdf", width = 6.5, height = 4.5)
print(enrichplot::gseaplot2(gsea_inflamm_oldpwr_v_oldsed, geneSetID = 1,
                            title = "MusAge ~ OLD PWR vs OLD SED"))
dev.off()

# 3b) Optional 3-panel (top ES, middle ticks, bottom red↔blue stat bars)
stats_vec <- sort(oldpwrveh_v_oldsed.vec, decreasing = TRUE)  # ranked statistic = logFC
pathway   <- intersect(names(stats_vec), musage_entrez)

N   <- length(stats_vec)
hit <- names(stats_vec) %in% pathway
Nh  <- sum(hit); Nm <- N - Nh
P   <- sum(abs(stats_vec[hit]))  # weighted p=1

runningES <- cumsum(ifelse(hit, abs(stats_vec) / P, -1 / Nm))
df_top  <- data.frame(Rank = seq_along(stats_vec), ES = runningES)
df_hits <- data.frame(Rank = which(hit))
df_bot  <- data.frame(Rank = seq_along(stats_vec), Stat = as.numeric(stats_vec))

p_top <- ggplot(df_top, aes(Rank, ES)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(data = df_top[which.max(abs(df_top$ES)), , drop = FALSE],
             aes(Rank, ES), color = "red", size = 2.2) +
  labs(y = "Enrichment Score (ES)", x = NULL,
       title = "GSEA: MusAge set in OLD PWR vs OLD SED") +
  theme_minimal(base_size = 11)

p_mid <- ggplot(df_hits, aes(Rank, 1)) +
  geom_segment(aes(xend = Rank, yend = 0), color = "black") +
  labs(y = NULL, x = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = -14, b = -14, l = 5, r = 5))

p_bot <- ggplot(df_bot, aes(Rank, Stat, fill = Stat)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Ranked statistic (logFC)") +
  labs(x = "Gene Rank", y = "logFC (OLD PWR − OLD SED)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

Fig2f_threepanel <- cowplot::plot_grid(p_top, p_mid, p_bot, ncol = 1, align = "v",
                                       rel_heights = c(3, 0.5, 2))
ggsave("Fig2f_threepanel_MusAge_OLDPWR_vs_OLDSED.pdf", Fig2f_threepanel,
       width = 7.0, height = 5.0, units = "in")
#------------------------------------------------------------------#

## ===============================
## Fig. 2g — MusAge (52) avg logCPM heatmap across OLD interventions
## ===============================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(purrr); library(edgeR); library(pheatmap)
  library(AnnotationDbi); library(org.Mm.eg.db)
})
suppressPackageStartupMessages(library(purrr))

## 0) Sample IDs (define only if not already present)
if (!exists("yng_sed_samples", inherits = FALSE))     yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
if (!exists("old_sedveh_samples", inherits = FALSE))  old_sedveh_samples  <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
if (!exists("old_pwrveh_samples", inherits = FALSE))  old_pwrveh_samples  <- c("T_08","T_10","T_11","T_16","T_41","T_43","T_45","T_46")
if (!exists("old_sedirap_samples", inherits = FALSE)) old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")
if (!exists("old_sedfrap_samples", inherits = FALSE)) old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")
if (!exists("old_pwrirap_samples", inherits = FALSE)) old_pwrirap_samples <- c("T_03","T_14","T_31","T_33","T_37","T_47","T_48")
if (!exists("old_pwrfrap_samples", inherits = FALSE)) old_pwrfrap_samples <- c("T_01","T_04","T_22","T_24","T_28","T_34","T_40","T_44")

## 1) Read + clean counts (drop YNG columns from each Set01 counts file)
clean_counts <- function(path) {
  df <- read_xlsx(path)
  colnames(df)[1] <- "Ensembl"
  df %>% dplyr::select(-all_of(yng_sed_samples))
}

old_sedveh_counts  <- clean_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
old_pwrveh_counts  <- clean_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_VEH-YNG_SED_VEH.xlsx")
old_sedirap_counts <- clean_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx")
old_sedfrap_counts <- clean_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx")
old_pwrirap_counts <- clean_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_IRAP-YNG_SED_VEH.xlsx")
old_pwrfrap_counts <- clean_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_PWR_FRAP-YNG_SED_VEH.xlsx")

combined_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "Ensembl"),
  list(
    old_sedveh_counts, old_pwrveh_counts, old_sedirap_counts,
    old_sedfrap_counts, old_pwrirap_counts, old_pwrfrap_counts
  )
)

combined_counts[is.na(combined_counts)] <- 0

## 2) Map Ensembl -> ENTREZ; collapse duplicates by ENTREZ (sum)
all_intervention_counts <- combined_counts %>%
  mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"),
         ENTREZID = mapIds(org.Mm.eg.db,
                           keys = Ensembl_noDec,
                           keytype = "ENSEMBL",
                           column = "ENTREZID",
                           multiVals = "first")) %>%
  filter(!is.na(ENTREZID)) %>%
  dplyr::select(-Ensembl, -Ensembl_noDec) %>%
  group_by(ENTREZID) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  column_to_rownames("ENTREZID")

## 3) Build counts matrix and TMM → logCPM across all OLD samples
all_samples <- c(old_sedveh_samples, old_pwrveh_samples, old_sedirap_samples,
                 old_sedfrap_samples, old_pwrirap_samples, old_pwrfrap_samples)
present <- intersect(all_samples, colnames(all_intervention_counts))
stopifnot(length(present) >= 2)
counts_matrix <- as.matrix(all_intervention_counts[, present, drop = FALSE])

group <- factor(c(rep("OLD_SEDVEH",  length(intersect(old_sedveh_samples,  present))),
                  rep("OLD_PWRVEH",  length(intersect(old_pwrveh_samples,  present))),
                  rep("OLD_SEDIRAP", length(intersect(old_sedirap_samples, present))),
                  rep("OLD_SEDFRAP", length(intersect(old_sedfrap_samples, present))),
                  rep("OLD_PWRIRAP", length(intersect(old_pwrirap_samples, present))),
                  rep("OLD_PWRFRAP", length(intersect(old_pwrfrap_samples, present)))))

dge <- DGEList(counts = counts_matrix, group = group)
dge <- calcNormFactors(dge, method = "TMM")
logCPM_matrix <- cpm(dge, log = TRUE, prior.count = 1)

## 4) MusAge 52 → ENTREZ; subset to those rows
if (!exists("validated_inflamm_geneset.vec")) {
  validated_inflamm_geneset.vec <- readr::read_csv("MusAge_geneset.csv", show_col_types = FALSE) |>
    pull(Gene) |> as.character() |> stringr::str_trim() |> unique()
}
if (exists("inflammaging_gene_vec")) {
  musAge_52_entrez <- unique(na.omit(inflammaging_gene_vec))
} else {
  musAge_52_entrez <- mapIds(org.Mm.eg.db,
                             keys = validated_inflamm_geneset.vec,
                             keytype = "SYMBOL",
                             column = "ENTREZID",
                             multiVals = "first") |> unname() |> na.omit() |> unique()
}
inflammaging_logCPM <- logCPM_matrix[rownames(logCPM_matrix) %in% musAge_52_entrez, , drop = FALSE]

## 5) Compute group averages (one column per intervention)
avg_logCPM <- data.frame(
  OLD_SedVeh  = rowMeans(inflammaging_logCPM[, intersect(old_sedveh_samples,  colnames(inflammaging_logCPM)),  drop = FALSE], na.rm = TRUE),
  OLD_PwrVeh  = rowMeans(inflammaging_logCPM[, intersect(old_pwrveh_samples,  colnames(inflammaging_logCPM)),  drop = FALSE], na.rm = TRUE),
  OLD_SedIRAP = rowMeans(inflammaging_logCPM[, intersect(old_sedirap_samples, colnames(inflammaging_logCPM)), drop = FALSE], na.rm = TRUE),
  OLD_SedFRAP = rowMeans(inflammaging_logCPM[, intersect(old_sedfrap_samples, colnames(inflammaging_logCPM)), drop = FALSE], na.rm = TRUE),
  OLD_PwrIRAP = rowMeans(inflammaging_logCPM[, intersect(old_pwrirap_samples, colnames(inflammaging_logCPM)), drop = FALSE], na.rm = TRUE),
  OLD_PwrFRAP = rowMeans(inflammaging_logCPM[, intersect(old_pwrfrap_samples, colnames(inflammaging_logCPM)), drop = FALSE], na.rm = TRUE)
)
rownames(avg_logCPM) <- rownames(inflammaging_logCPM)

## 6) Replace ENTREZ rownames with SYMBOLs; collapse duplicate symbols by mean
sym_map <- mapIds(org.Mm.eg.db,
                  keys = rownames(avg_logCPM),
                  keytype = "ENTREZID",
                  column = "SYMBOL",
                  multiVals = "first")
avg_logCPM_sym <- avg_logCPM %>%
  rownames_to_column("ENTREZID") %>%
  mutate(Symbol = sym_map) %>%
  drop_na(Symbol) %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

## 7) Heatmap (rows & cols clustered; row z-scored; blue→white→red)
#pdf("Fig2g_MusAge_avg_logCPM_heatmap.pdf", width = 7.2, height = 8.6)
pheatmap(
  avg_logCPM_sym,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue","white","red"))(60),
  main = "MusAge (52) — Avg logCPM per OLD intervention"
)
#dev.off()
#---------------------------------------------------------#

## =========================
## Fig. 2h — MusAge GSEA NES (interventions vs OLD SED)
## =========================
suppressPackageStartupMessages({ library(purrr); library(readr) })

## 0) Make sure MusAge is ENTREZ
if (exists("inflammaging_gene_vec")) {
  musage_entrez <- unique(na.omit(as.character(inflammaging_gene_vec)))
} else if (exists("validated_inflamm_geneset.vec")) {
  musage_entrez <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys     = validated_inflamm_geneset.vec,
    keytype  = "SYMBOL",
    column   = "ENTREZID",
    multiVals = "first"
  ) |> unname() |> as.character() |> unique() |> na.omit()
} else {
  stop("Provide MusAge 52 genes as SYMBOLs (validated_inflamm_geneset.vec) or ENTREZ (inflammaging_gene_vec).")
}

## 1) Helper: read *_GENE_* XLSX and return ranked ENTREZ->logFC vector
rank_vec_from <- function(xlsx_file, logFC_col, FDR_col, objname) {
  df <- read_clean_xlsx(xlsx_file, logFC_col, FDR_col)
  create_vec_from_df(df, objname)                       # creates <objname>.vec in .GlobalEnv
  v  <- get(paste0(objname, ".vec"), envir = .GlobalEnv)
  v  <- v[!is.na(v)]
  v  <- v[!duplicated(names(v))]
  sort(v, decreasing = TRUE)                            # strictly decreasing
}

## 2) Build ranked lists for each intervention vs OLD SED
ranked <- list(
  OLD_PWR_VEH  = rank_vec_from(
    "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
    "OLD_PWR_VEH-OLD_SED_VEH_logFC", "OLD_PWR_VEH-OLD_SED_VEH_FDR", "oldpwrveh_v_oldsedveh"),
  OLD_SED_FRAP = rank_vec_from(
    "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_FRAP-OLD_SED_VEH.xlsx",
    "OLD_SED_FRAP-OLD_SED_VEH_logFC","OLD_SED_FRAP-OLD_SED_VEH_FDR","oldsedfrap_v_oldsedveh"),
  OLD_SED_IRAP = rank_vec_from(
    "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_IRAP-OLD_SED_VEH.xlsx",
    "OLD_SED_IRAP-OLD_SED_VEH_logFC","OLD_SED_IRAP-OLD_SED_VEH_FDR","oldsedirap_v_oldsedveh"),
  OLD_PWR_IRAP = rank_vec_from(
    "20250530_M007853_Set05_edgeRglm_GENE_OLD_PWR_IRAP-OLD_SED_VEH.xlsx",
    "OLD_PWR_IRAP-OLD_SED_VEH_logFC","OLD_PWR_IRAP-OLD_SED_VEH_FDR","oldpwrirap_v_oldsedveh"),
  OLD_PWR_FRAP = rank_vec_from(
    "20250530_M007853_Set05_edgeRglm_GENE_OLD_PWR_FRAP-OLD_SED_VEH.xlsx",
    "OLD_PWR_FRAP-OLD_SED_VEH_logFC","OLD_PWR_FRAP-OLD_SED_VEH_FDR","oldpwrfrap_v_oldsedveh")
)

## 3) Targeted GSEA (MusAge 52) and extract NES / p-values
set.seed(1)
nes_tbl <- imap_dfr(ranked, function(gl, label) {
  g <- clusterProfiler::GSEA(
    geneList     = gl,
    TERM2GENE    = data.frame(term = "MusAge52", gene = musage_entrez),
    pvalueCutoff = 1, minGSSize = 1, maxGSSize = 10000,
    verbose      = FALSE
  )
  n_in <- sum(names(gl) %in% musage_entrez)
  if (nrow(g@result) == 0) {
    tibble(Intervention = label, NES = NA_real_, pvalue = NA_real_, p.adjust = NA_real_, n_in_set = n_in)
  } else {
    r <- g@result[1, ]
    tibble(Intervention = label, NES = r$NES, pvalue = r$pvalue, p.adjust = r$p.adjust, n_in_set = n_in)
  }
}) |>
  mutate(sig = ifelse(!is.na(pvalue) & pvalue < 0.01, "**", "")) |>
  # order like the figure (tweak if you prefer a different order)
  mutate(Intervention = factor(Intervention,
                               levels = c("OLD_SED_IRAP","OLD_SED_FRAP","OLD_PWR_VEH","OLD_PWR_IRAP","OLD_PWR_FRAP"))) |>
  arrange(Intervention)

print(nes_tbl)

# Save for Prism
readr::write_csv(nes_tbl, "Fig2h_MusAge_GSEA_vs_OLDSED_NES.csv")

############################################################
# End of Figure 2 consolidated script
############################################################

######################################################
#=====================================================#
#Figure 3: Consolidated Code
#=====================================================#
######################################################


## =========================
## Fig. 3a — Lipid dotplot (OLD SED vs YNG SED)
## =========================
library(tidyverse)

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
os_ys <- lip_log %>% filter(Group %in% c("OS","YS"))

# 6) Per-lipid stats (log2FC = mean(OS) − mean(YS); t-test; FDR)
os_ys_stats <- purrr::map_dfr(lipid_cols, function(lip) {
  vals  <- os_ys[[lip]]
  grp   <- os_ys$Group
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

# 7) Lipid class from the part before the first "(" (e.g., "Cer", "PC", "SM", …)
os_ys_stats <- os_ys_stats %>%
  mutate(
    LipidClass  = stringr::str_trim(stringr::str_extract(Lipid, "^[^\\(]+")),
    IsSignature = Log2_FC_OS_vs_YS > 0.4 & FDR < 0.2
  )

# (Optional) export table for Prism
#readr::write_csv(os_ys_stats, "Fig3a_OSvsYS_lipid_stats.csv")

# 8) Choose class order (A) by median logFC descending  — or (B) your fixed order
# (A) data-driven:
class_order <- os_ys_stats %>%
  group_by(LipidClass) %>%
  summarize(med_logFC = median(Log2_FC_OS_vs_YS, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med_logFC)) %>%
  pull(LipidClass)

# (B) fixed order from your notes (uncomment to enforce):
# class_order <- c("CE","Sph","LPE","Cer","FA","TG","PA","DG","LPC","PS","HexCer","SM","PC","PE","PG","PI")

plot_df <- os_ys_stats %>%
  mutate(LipidClass = factor(LipidClass, levels = class_order))

# 9) Dotplot: each point = lipid species; x = log2FC (OS−YS); y = class; red outline = signature
fig3a <- ggplot(plot_df, aes(x = Log2_FC_OS_vs_YS, y = LipidClass)) +
  # base layer: all lipids
  geom_jitter(width = 0, height = 0.22, size = 2, alpha = 0.8, color = "grey35") +
  # highlight: signature (red circles on top)
  geom_point(
    data = subset(plot_df, IsSignature),
    aes(x = Log2_FC_OS_vs_YS, y = LipidClass),
    shape = 21, fill = NA, color = "red3", stroke = 1.1, size = 3.8
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = "log2 Fold Change (OLD SED − YNG SED)",
    y = "Lipid class",
    title = "OLD SED vs YNG SED lipidomics"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

print(fig3a)
#ggsave("Fig3a_OSvsYS_lipid_dotplot.pdf", fig3a, width = 7, height = 5.2, units = "in")

# 10) Also save the signature list for reuse downstream
#sig_lipid_features <- os_ys_stats %>%
#  filter(IsSignature) %>%
#  arrange(Lipid)
#readr::write_csv(sig_lipid_features, "sig_lipid_features.csv")

## =========================
## Fig. 3b — 34-lipid signature heatmap (YS vs OS) + bar plot
## =========================
library(pheatmap)
library(ggplot2)
library(forcats)

# 1) Define the 34-lipid Aging Signature from Fig. 3a stats
sig_lipid_features <- os_ys_stats %>%
  dplyr::filter(Log2_FC_OS_vs_YS > 0.4, FDR < 0.2) %>%
  dplyr::arrange(Lipid)
signature34 <- sig_lipid_features$Lipid

# 2) Row order for all panels: sort by OS−YS log2FC (desc)
lip_order_sig <- os_ys_stats %>%
  dplyr::filter(Lipid %in% signature34) %>%
  dplyr::arrange(dplyr::desc(Log2_FC_OS_vs_YS)) %>%
  dplyr::pull(Lipid)

# 3) Build z-scored matrix (rows = lipids, cols = samples), YS left → OS right
heatmap_data_osys <- lip_log %>%
  dplyr::filter(Group %in% c("YS","OS")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("YS","OS"))) %>%
  dplyr::arrange(Group) %>%
  dplyr::select(Sample, Group, dplyr::all_of(lip_order_sig))

mat_osys <- heatmap_data_osys %>%
  dplyr::select(-Sample, -Group) %>%
  as.matrix()
rownames(mat_osys) <- heatmap_data_osys$Sample
mat_osys_z <- scale(mat_osys)  # z-score per lipid across samples

ann_col_osys <- data.frame(Group = heatmap_data_osys$Group)
rownames(ann_col_osys) <- heatmap_data_osys$Sample
ann_cols <- list(Group = c(YS = "#89CFF0", OS = "#4F4F4F"))

pal <- colorRampPalette(c("#FF7F00","black","#00FFFF"))(100)

hm_3b <- pheatmap(
  t(mat_osys_z),
  annotation_col   = ann_col_osys,
  annotation_colors= ann_cols,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row   = lip_order_sig, row_names_side = "left",
  fontsize = 8, color = pal, breaks = seq(-1, 1, length.out = 101),
  main = "34-lipid Aging Signature: YS (left, blue) vs OS (right, grey)"
)

print(hm_3b)

# 4) Horizontal bar plot (OS − YS); black bars indicate > 0
bar3b_df <- os_ys_stats %>%
  dplyr::filter(Lipid %in% signature34) %>%
  dplyr::mutate(Lipid = factor(Lipid, levels = lip_order_sig))

bar_3b <- ggplot(bar3b_df, aes(x = Log2_FC_OS_vs_YS, y = fct_rev(Lipid))) +
  geom_col(fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "log2 Fold Change (OLD SED − YNG SED)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank())

print(bar_3b)

# ggsave("Fig3b_OSvsYS_signature_heatmap.pdf", width = 4.25, height = 5.5)
# ggsave("Fig3b_OSvsYS_signature_barplot.pdf", bar_3b, width = 4.25, height = 5.5)

## =========================
## Fig. 3e — 34-lipid signature heatmap (OV vs OS) + bar plot
## =========================

# 1) Stats for OV vs OS (log2FC = OV − OS)
ov_os <- lip_log %>% dplyr::filter(Group %in% c("OV","OS"))
ov_os_stats <- purrr::map_dfr(lipid_cols, function(lip) {
  vals <- ov_os[[lip]]; grp <- ov_os$Group; tt <- t.test(vals ~ grp)
  tibble::tibble(
    Lipid = lip,
    Mean_OV = mean(vals[grp == "OV"], na.rm = TRUE),
    Mean_OS = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OV_vs_OS = Mean_OV - Mean_OS,
    P_value = tt$p.value
  )
}) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  dplyr::arrange(FDR)

# 2) Heatmap with SAME lipid order as Fig. 3b (for visual comparability)
heatmap_data_ovos <- lip_log %>%
  dplyr::filter(Group %in% c("OV","OS")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("OV","OS"))) %>%
  dplyr::arrange(Group) %>%
  dplyr::select(Sample, Group, dplyr::all_of(lip_order_sig))

mat_ovos <- heatmap_data_ovos %>%
  dplyr::select(-Sample, -Group) %>%
  as.matrix()
rownames(mat_ovos) <- heatmap_data_ovos$Sample
mat_ovos_z <- scale(mat_ovos)

ann_col_ovos <- data.frame(Group = heatmap_data_ovos$Group)
rownames(ann_col_ovos) <- heatmap_data_ovos$Sample
ann_cols_ov <- list(Group = c(OV = "#98FB98", OS = "#4F4F4F"))

hm_3e <- pheatmap(
  t(mat_ovos_z),
  annotation_col   = ann_col_ovos,
  annotation_colors= ann_cols_ov,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row   = lip_order_sig, row_names_side = "left",
  fontsize = 8, color = pal, breaks = seq(-1, 1, length.out = 101),
  main = "34-lipid Aging Signature: OLD PWR (left, green) vs OS (right, grey)"
)

print(hm_3e)

# 3) Horizontal bar plot (OV − OS); red bars indicate < 0 (down in OV)
bar3e_df <- ov_os_stats %>%
  dplyr::filter(Lipid %in% signature34) %>%
  dplyr::mutate(
    Lipid = factor(Lipid, levels = lip_order_sig),
    Direction = ifelse(Log2_FC_OV_vs_OS < 0, "Down_in_OV", "Up_in_OV")
  )

bar_3e <- ggplot(bar3e_df,
                 aes(x = Log2_FC_OV_vs_OS, y = fct_rev(Lipid), fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = c(Down_in_OV = "red3", Up_in_OV = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "log2 Fold Change (OLD PWR − OLD SED)", y = NULL, fill = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none")

print(bar_3e)

# ggsave("Fig3e_OLDPWRvsOS_signature_heatmap.pdf", width = 4.25, height = 5.5)
# ggsave("Fig3e_OLDPWRvsOS_barplot.pdf", bar_3e, width = 4.25, height = 5.5)
#------------------------------------------------------------------------------#

## =========================
## Fig. 3c–d — Overlap with Riken lipidome (24m vs 2m) + side-by-side bar/points
## =========================
## Files expected:
##  • aging_lipidome.xlsx  (Riken female GF sheet/columns)
##  • OPTIONAL: new_lipid_age_signature_formatted.xlsx (pre-formatted names)
## Uses signature34 from Fig. 3b if the formatted list is absent.

#aging lipidome analysis

# Reads the first sheet by default
aging_lipidome <- read_excel("aging_lipidome.xlsx")

aging_lipidome
aging_lipidome$Metabolite <- sub("\\|.*", "", aging_lipidome$Metabolite)

aging_lipidome$Metabolite <- sub(";.*", "", aging_lipidome$Metabolite)
aging_lipidome

aging_lipidome$Metabolite <- sub("/0:0.*", "", aging_lipidome$Metabolite)
aging_lipidome

formatted_lipid_age_sig <- read_excel("new_lipid_age_signature_formatted.xlsx")
formatted_lipid_age_sig.vec <- formatted_lipid_age_sig %>%
  pull(Lipid)

formatted_lipid_age_sig.vec
aging_lipidome %>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec)

colnames(aging_lipidome)


# Select only columns that end with "F_GF"
female_aging_lipidome <- aging_lipidome %>%
  dplyr::select(Metabolite, contains("F_GF"))


female_aging_lipidome%>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec)

female_aging_lipidome

formatted_lipid_age_sig.vec
#-----------------------------------------------#

library(purrr)
library(tibble)

# Step 1: Define column names for the groups
young_grp <- c("2M_F_GF...9", "2M_F_GF...34", "2M_F_GF...36")
old_grp   <- c("24M_F_GF...21", "24M_F_GF...27", "24M_F_GF...42")

# Step 2: Identify lipid columns (i.e., all group columns)
lipid_cols <- c(young_grp, old_grp)

# Step 3: Impute zero values with half of the minimum non-zero per column
min_nonzero <- sapply(female_aging_lipidome[lipid_cols], function(x) {
  min(x[x > 0], na.rm = TRUE)
})

lipid_imputed <- female_aging_lipidome
for (lip in lipid_cols) {
  x <- lipid_imputed[[lip]]
  x[x == 0] <- min_nonzero[[lip]] / 2
  lipid_imputed[[lip]] <- x
}

# Step 4: Log2-transform the imputed values
lipid_log <- lipid_imputed %>%
  mutate(across(all_of(lipid_cols), log2))

lipid_log

# Step 5: Run t-tests per lipid species
aging_lipidome_stats <- purrr::map_dfr(1:nrow(lipid_log), function(i) {
  lipid_name <- lipid_log$Metabolite[i]
  vals_young <- as.numeric(lipid_log[i, young_grp])
  vals_old   <- as.numeric(lipid_log[i, old_grp])
  
  # Skip rows with zero variance in either group
  if (sd(vals_young) == 0 || sd(vals_old) == 0) {
    return(tibble(
      Metabolite = lipid_name,
      Mean_Young = mean(vals_young),
      Mean_Old   = mean(vals_old),
      Log2_FC_Old_vs_Young = mean(vals_old) - mean(vals_young),
      P_value = NA_real_
    ))
  }
  
  # Otherwise, perform t-test
  ttest <- t.test(vals_old, vals_young)
  
  tibble(
    Metabolite = lipid_name,
    Mean_Young = mean(vals_young),
    Mean_Old   = mean(vals_old),
    Log2_FC_Old_vs_Young = mean(vals_old) - mean(vals_young),
    P_value = ttest$p.value
  )
}) %>%
  mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  arrange(FDR)

# Optional: View top hits
head(aging_lipidome_stats)

aging_lipidome_stats %>%
  filter(FDR <0.05)

print(aging_lipidome_stats %>%
        filter(Metabolite %in% formatted_lipid_age_sig.vec), n=21)

aging_lipidome_stats_unique <- aging_lipidome_stats %>%
  distinct(Metabolite, .keep_all = TRUE)

dim(aging_lipidome_stats_unique)

print(aging_lipidome_stats_unique %>%
        filter(Metabolite %in% formatted_lipid_age_sig.vec), n=21)

ovrlp_aging_sig_lipidome <- aging_lipidome_stats_unique %>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec) %>%
  arrange(desc(Log2_FC_Old_vs_Young))

write.csv(ovrlp_aging_sig_lipidome, "ovrlp_aging_sig_lipidome.csv", row.names = FALSE)

ovrlp_aging_lipid_vec <- aging_lipidome_stats_unique %>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec) %>%
  arrange(desc(Log2_FC_Old_vs_Young)) %>%
  pull(Metabolite)

ovrlp_aging_lipid_vec
formatted_lipid_age_sig.vec
#---------------------------------------------#

#============================================#
#-----lgsea from lipidome--------------------#
#============================================#

# Ensure ranking vector is built from lipids with proper names and values
rank_df <- aging_lipidome_stats_unique %>%
  filter(!is.na(Log2_FC_Old_vs_Young), !is.na(P_value)) %>%
  mutate(
    p_safe = pmax(P_value, 1e-300),
    stat = Log2_FC_Old_vs_Young
  ) %>%
  group_by(Metabolite) %>%
  summarise(stat = mean(stat), .groups = "drop") %>%
  arrange(desc(stat))

# Convert to named vector
lipid_ranks <- setNames(rank_df$stat, rank_df$Metabolite)

# Only use the lipids that are present in the ranked list
lipid_signature_filtered <- formatted_lipid_age_sig.vec[
  formatted_lipid_age_sig.vec %in% names(lipid_ranks)
]

lipid_sets <- list(AgingSignature = lipid_signature_filtered)

fgsea_res <- fgsea(
  pathways = lipid_sets,
  stats = lipid_ranks,
  nperm = 1000
)

print(fgsea_res)

# Plot enrichment
plotEnrichment(lipid_sets$AgingSignature, lipid_ranks) +
  ggtitle("Effect of Age on Aging Lipid Signature")

female_lipidome_signature <- aging_lipidome_stats_unique %>%
  filter(FDR <0.05) %>%
  pull(Metabolite)

female_lipidome_signature

aging_lipidome_stats_unique %>%
  filter (Log2_FC_Old_vs_Young >0.49)

aging_lipidome_stats_unique %>%
  filter(FDR >0.05 & Log2_FC_Old_vs_Young >0)

dim(aging_lipidome_stats_unique)
#--------------------------------------------#

# Subset dataframe
lipid_log_subset <- lipid_log %>%
  dplyr::select(Metabolite, any_of(c(young_grp, old_grp)))

lipid_log_subset %>%
  filter(Metabolite %in% ovrlp_aging_lipid_vec) %>%
  distinct(Metabolite, .keep_all = TRUE)

# Step 1: Compute the young mean
lipid_log_wide_logFC <- lipid_log_subset %>%
  rowwise() %>%
  mutate(young_mean = mean(c_across(all_of(young_grp)), na.rm = TRUE)) %>%
  ungroup()

# Step 2: Add logFC columns for each 24M sample
for (sample in old_grp) {
  new_col <- paste0("logFC_", sample)
  lipid_log_wide_logFC[[new_col]] <- lipid_log_wide_logFC[[sample]] - lipid_log_wide_logFC$young_mean
}

# Optional: Arrange columns for readability
lipid_log_wide_logFC <- lipid_log_wide_logFC %>%
  dplyr::select(Metabolite, all_of(young_grp), young_mean, all_of(old_grp), starts_with("logFC_"))

lipid_log_wide_logFC %>%
  filter(Metabolite %in% ovrlp_aging_lipid_vec) %>%
  distinct(Metabolite, .keep_all = TRUE) %>%
  dplyr::select(Metabolite, starts_with("logFC_"))

lipid_log_wide_logFC %>%
  filter(Metabolite == "PC 38:6") %>%
  dplyr::select(Metabolite, starts_with("logFC_"))

## =========================
## Fig. 3c–d (our side): aggregate to generic names & compute log2FC
## =========================
library(dplyr)
library(tidyr)
library(stringr)

# Helper: normalize our lipid labels to Riken-style generic names (e.g., "PC 38:6", "PE O-38:4", "Sph d18:1")
extract_base <- function(x) {
  x <- gsub("\\|.*|;.*|/0:0.*", "", x)          # drop trailing annotations
  x <- gsub("\\(|\\)|_", " ", x)                # remove parens / underscores
  x <- gsub("\\s+", " ", trimws(x))             # collapse spaces
  # keep Class + optional O-/P- + (optionally 'd')XX:YY
  m <- str_match(x, "^([A-Za-z]+\\s(?:O-|P-)?(?:d)?\\d{2}:\\d{1,2})")[,2]
  ifelse(is.na(m), x, m)
}

# --- overlap list you already computed from Riken side ---
# uses: ovrlp_aging_lipid_vec (the 15 generic names), riken_stats, riken_sample_logfc
stopifnot(exists("ovrlp_aging_lipid_vec"))

# --- (A) Our mean log2FC per generic name (average across species) ---
our_bars <- os_ys_stats %>%
  mutate(Base = extract_base(Lipid)) %>%
  filter(Base %in% ovrlp_aging_lipid_vec) %>%
  group_by(Base) %>%
  summarise(
    MeanLogFC = mean(Log2_FC_OS_vs_YS, na.rm = TRUE),
    n_species = n(),
    .groups = "drop"
  ) %>%
  transmute(Metabolite = Base, Study = "Our study", MeanLogFC)

# --- (B) Our per-sample logFC dots aggregated to generic name ---
os_ys


## =========================
## Fig. 3c–d — overlap vs Riken (robust dplyr-qualified)
## Requires existing objects: os_ys, os_ys_stats,
##   aging_lipidome_stats_unique, lipid_log_wide_logFC,
##   ovrlp_aging_lipid_vec (your 15-lipid generic names)
## =========================
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Helper to collapse species to a generic name (e.g., "PC 38:6", "PE O-38:4", "Sph d18:1")
if (!exists("extract_base")) {
  extract_base <- function(x) {
    x <- gsub("\\|.*|;.*|/0:0.*", "", x)   # drop trailing annotations
    x <- gsub("\\(|\\)|_", " ", x)         # remove parens/underscores
    x <- gsub("\\s+", " ", trimws(x))      # collapse spaces
    m <- stringr::str_match(x, "^([A-Za-z]+\\s(?:O-|P-)?(?:d)?\\d{2}:\\d{1,2})")[,2]
    ifelse(is.na(m), x, m)
  }
}

# If the overlap vector isn't defined, derive a reasonable one from your 34-lipid signature vs Riken
if (!exists("ovrlp_aging_lipid_vec")) {
  riken_stats <- aging_lipidome_stats_unique %>%
    dplyr::transmute(Metabolite, Log2FC_Riken = Log2_FC_Old_vs_Young)
  sig34_generic <- extract_base(signature34)
  ovrlp_aging_lipid_vec <- intersect(sig34_generic, riken_stats$Metabolite)
  ovrlp_aging_lipid_vec <- os_ys_stats %>%
    dplyr::mutate(Base = extract_base(Lipid)) %>%
    dplyr::filter(Base %in% ovrlp_aging_lipid_vec) %>%
    dplyr::arrange(dplyr::desc(Log2_FC_OS_vs_YS)) %>%
    dplyr::distinct(Base) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::pull(Base)
}

# --- Our study: bars (mean of species mapping to each generic name) ---
our_bars <- os_ys_stats %>%
  dplyr::mutate(Base = extract_base(Lipid)) %>%
  dplyr::filter(Base %in% ovrlp_aging_lipid_vec) %>%
  dplyr::group_by(Base) %>%
  dplyr::summarise(MeanLogFC = mean(Log2_FC_OS_vs_YS, na.rm = TRUE), .groups = "drop") %>%
  dplyr::transmute(Metabolite = Base, Study = "Our study", MeanLogFC)

# --- Our study: dots (per-sample logFC averaged within generic name) ---
lipid_cols_osys <- setdiff(names(os_ys), c("Sample", "Group"))

osys_long <- os_ys %>%
  tidyr::pivot_longer(cols = dplyr::all_of(lipid_cols_osys),
                      names_to = "Lipid", values_to = "val") %>%
  dplyr::mutate(Base = extract_base(Lipid)) %>%
  dplyr::filter(Base %in% ovrlp_aging_lipid_vec)

ys_base_mean <- osys_long %>%
  dplyr::filter(Group == "YS") %>%
  dplyr::group_by(Base) %>%
  dplyr::summarise(YS_mean = mean(val, na.rm = TRUE), .groups = "drop")

our_dots <- osys_long %>%
  dplyr::filter(Group == "OS") %>%
  dplyr::left_join(ys_base_mean, by = "Base") %>%
  dplyr::group_by(Base, Sample) %>%
  dplyr::summarise(LogFC = mean(val - YS_mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Study = "Our study") %>%
  dplyr::rename(Metabolite = Base)   # <- explicit dplyr::rename fixes your error

# --- Riken: bars + dots ---
riken_stats <- aging_lipidome_stats_unique %>%
  dplyr::transmute(Metabolite, Log2FC_Riken = Log2_FC_Old_vs_Young)

bars_riken <- riken_stats %>%
  dplyr::filter(Metabolite %in% ovrlp_aging_lipid_vec) %>%
  dplyr::transmute(Metabolite, Study = "Riken", MeanLogFC = Log2FC_Riken)

riken_sample_logfc <- lipid_log_wide_logFC %>%
  dplyr::filter(Metabolite %in% ovrlp_aging_lipid_vec) %>%
  tidyr::pivot_longer(dplyr::starts_with("logFC_"),
                      names_to = "Sample", values_to = "LogFC",
                      names_prefix = "logFC_") %>%
  dplyr::mutate(Study = "Riken") %>%
  dplyr::select(Metabolite, Sample, LogFC, Study)

# --- Combine & order for plotting ---
bars_df <- dplyr::bind_rows(our_bars, bars_riken) %>%
  dplyr::mutate(Study = factor(Study, levels = c("Our study", "Riken")))

dots_df <- dplyr::bind_rows(
  our_dots,
  riken_sample_logfc %>% dplyr::filter(Metabolite %in% ovrlp_aging_lipid_vec)
) %>%
  dplyr::mutate(Study = factor(Study, levels = c("Our study", "Riken")))

lip_order <- our_bars %>% dplyr::arrange(dplyr::desc(MeanLogFC)) %>% dplyr::pull(Metabolite)
bars_df <- bars_df %>% dplyr::mutate(Metabolite = factor(Metabolite, levels = rev(lip_order)))
dots_df <- dots_df %>% dplyr::mutate(Metabolite = factor(Metabolite, levels = rev(lip_order)))

# --- Plot (bars = means; dots = samples) ---
set.seed(1)
fig3d <- ggplot(bars_df, aes(x = MeanLogFC, y = Metabolite, fill = Study)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = NA) +
  geom_point(data = dots_df,
             aes(x = LogFC, y = Metabolite, color = Study),
             position = position_jitterdodge(jitter.height = 0.12,
                                             jitter.width  = 0.03,
                                             dodge.width   = 0.7),
             size = 1.8, alpha = 0.75, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  scale_fill_manual(values = c("Our study" = "#F19CBB", "Riken" = "#98FB98")) +
  scale_color_manual(values = c("Our study" = "#C72267", "Riken" = "#2E8B57")) +
  labs(x = "log2 Fold Change (older − young)", y = NULL,
       title = "Overlap (15 lipids): our study vs Riken (24m vs 2m)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "top")

print(fig3d)
cat("# Fig3d — n bars:", nrow(bars_df), "  n dots:", nrow(dots_df), "\n")

# --- Fisher exact test (Riken: log2FC>0.5 among overlap vs rest) ---
riken_over_05 <- riken_stats %>% dplyr::mutate(hi = Log2FC_Riken > 0.49)
k_overlap <- sum(riken_over_05$hi[riken_over_05$Metabolite %in% ovrlp_aging_lipid_vec], na.rm = TRUE)
n_overlap <- sum(riken_over_05$Metabolite %in% ovrlp_aging_lipid_vec)
k_rest    <- sum(riken_over_05$hi[!(riken_over_05$Metabolite %in% ovrlp_aging_lipid_vec)], na.rm = TRUE)
n_rest    <- sum(!(riken_over_05$Metabolite %in% ovrlp_aging_lipid_vec))
fisher_tab <- matrix(c(k_overlap, n_overlap - k_overlap, k_rest, n_rest - k_rest), nrow = 2, byrow = TRUE)
fisher_out <- fisher.test(fisher_tab, alternative = "greater")
print(list(Fisher_table = fisher_tab, p_value = fisher_out$p.value))
#---------------------------------------------------------------------------#

## =========================
## Fig. 3e — OLD PWR (OV) vs OLD SED (OS): heatmap + bar plot
## (You already had a working version; keeping it here verbatim for continuity.)
## =========================

# Recompute OV vs OS stats (log2FC = OV − OS) just to be explicit here
ov_os <- lip_log %>% dplyr::filter(Group %in% c("OV","OS"))
ov_os_stats <- purrr::map_dfr(lipid_cols, function(lip) {
  vals <- ov_os[[lip]]; grp <- ov_os$Group; tt <- t.test(vals ~ grp)
  tibble::tibble(
    Lipid = lip,
    Mean_OV = mean(vals[grp == "OV"], na.rm = TRUE),
    Mean_OS = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OV_vs_OS = Mean_OV - Mean_OS,
    P_value = tt$p.value
  )
}) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  dplyr::arrange(FDR)

# Heatmap with SAME lipid order as Fig. 3b (lip_order_sig)
heatmap_data_ovos <- lip_log %>%
  dplyr::filter(Group %in% c("OV","OS")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("OV","OS"))) %>%   # OV left, OS right
  dplyr::arrange(Group) %>%
  dplyr::select(Sample, Group, dplyr::all_of(lip_order_sig))

mat_ovos <- heatmap_data_ovos %>% dplyr::select(-Sample, -Group) %>% as.matrix()
rownames(mat_ovos) <- heatmap_data_ovos$Sample
mat_ovos_z <- scale(mat_ovos)

ann_col_ovos <- data.frame(Group = heatmap_data_ovos$Group)
rownames(ann_col_ovos) <- heatmap_data_ovos$Sample
ann_cols_ov <- list(Group = c(OV = "#98FB98", OS = "#4F4F4F"))  # green vs grey

hm_3e <- pheatmap::pheatmap(
  t(mat_ovos_z),
  annotation_col   = ann_col_ovos,
  annotation_colors= ann_cols_ov,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row   = lip_order_sig, row_names_side = "left",
  fontsize = 8, color = pal, breaks = seq(-1, 1, length.out = 101),
  main = "34-lipid Aging Signature: OLD PWR (left, green) vs OLD SED (right, grey)"
)
print(hm_3e)

# Horizontal bar plot (OV − OS); red < 0 (down in OV), black > 0 (up in OV)
bar3e_df <- ov_os_stats %>%
  dplyr::filter(Lipid %in% signature34) %>%
  dplyr::mutate(
    Lipid = factor(Lipid, levels = lip_order_sig),
    Direction = ifelse(Log2_FC_OV_vs_OS < 0, "Down_in_OV", "Up_in_OV")
  )

bar_3e <- ggplot2::ggplot(bar3e_df,
                          ggplot2::aes(x = Log2_FC_OV_vs_OS, y = forcats::fct_rev(Lipid), fill = Direction)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c(Down_in_OV = "red3", Up_in_OV = "black")) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::labs(x = "log2 Fold Change (OLD PWR − OLD SED)", y = NULL, fill = NULL) +
  ggplot2::theme_minimal(base_size = 8) +
  ggplot2::theme(panel.grid.major.y = element_blank(), legend.position = "none")
print(bar_3e)



## =========================
## Fig. 3f — Heatmap across all OLD groups
## Rows: 34-lipid Aging Signature; values = z-scored per lipid
## Orange = low; teal/green = high; units not confirmed
## =========================

library(dplyr)
library(tidyr)
library(pheatmap)
library(stringr)

# helper to sort sample names like OS1..OS8 numerically
.sort_by_suffix_num <- function(x) {
  ord <- suppressWarnings(as.numeric(sub("^\\D+", "", x)))
  x[order(ord, na.last = TRUE)]
}

# desired group display order
grp_order <- c("OS", "OV", "OIR", "OFR")

# color palette: orange -> black -> teal/green
pal_3f <- colorRampPalette(c("#FF7F00", "black", "#00BFAE"))(100)

# group color key
grp_cols <- c(OS = "#4F4F4F", OV = "#98FB98", OIR = "#66CDAA", OFR = "#2E8B57")

# -------- Option A (preferred): build from lip_log you already imputed/log2'ed --------
if (exists("lip_log") && all(c("Sample","Group") %in% names(lip_log))) {
  
  # derive group from Sample prefix if Group doesn't include OIR/OFR explicitly
  sample_groups <- lip_log %>%
    distinct(Sample, Group) %>%
    mutate(Group_derived = case_when(
      str_starts(Sample, "OS")  ~ "OS",
      str_starts(Sample, "OV")  ~ "OV",
      str_starts(Sample, "OIR") ~ "OIR",
      str_starts(Sample, "OFR") ~ "OFR",
      TRUE ~ as.character(Group)
    ))
  
  # use derived when available; keep only OLD groups that exist
  old_groups_found <- intersect(grp_order, unique(sample_groups$Group_derived))
  stopifnot(length(old_groups_found) > 0)
  
  hm_df <- lip_log %>%
    dplyr::select(Sample, Group, all_of(signature34)) %>%
    left_join(sample_groups %>% dplyr::select(Sample, Group_derived), by = "Sample") %>%
    mutate(Group2 = factor(Group_derived, levels = old_groups_found)) %>%
    filter(!is.na(Group2)) %>%
    dplyr::select(Sample, Group = Group2, all_of(signature34))
  
  # order samples: OS -> OV -> OIR -> OFR; numeric within each group
  samp_order <- unlist(lapply(old_groups_found, function(g) {
    .sort_by_suffix_num(hm_df$Sample[hm_df$Group == g])
  }))
  hm_df <- hm_df %>% mutate(Sample = factor(Sample, levels = samp_order)) %>% arrange(Sample)
  
  # use same row order as Fig. 3b for comparability (keep intersection just in case)
  rows_34 <- lip_order_sig[lip_order_sig %in% signature34 & lip_order_sig %in% names(hm_df)]
  
  # matrix: samples x lipids, then z-score per lipid (column)
  mat_allold <- hm_df %>% dplyr::select(-Sample, -Group) %>% as.matrix()
  rownames(mat_allold) <- hm_df$Sample
  mat_allold <- mat_allold[, rows_34, drop = FALSE]
  
  mat_allold_z <- scale(mat_allold)
  mat_allold_z[is.na(mat_allold_z)] <- 0  # handle constant rows defensively
  
  ann_col <- data.frame(Group = hm_df$Group)
  rownames(ann_col) <- hm_df$Sample
  ann_cols_all <- list(Group = grp_cols[names(grp_cols) %in% levels(ann_col$Group)])
  
  hm_3f <- pheatmap::pheatmap(
    t(mat_allold_z),                            # rows = lipids
    annotation_col    = ann_col,
    annotation_colors = ann_cols_all,
    cluster_rows = FALSE, cluster_cols = FALSE,
    labels_row   = rows_34, row_names_side = "left",
    fontsize = 8, color = pal_3f, breaks = seq(-1, 1, length.out = 101),
    main = "34-lipid Aging Signature across OLD groups\n(z-scored per lipid; units of raw abundance not confirmed)"
  )
  print(hm_3f)
  
}
#-------------------------------------------#

## =========================
## Fig. 3g — Three targeted LSEA plots (OV, OIR, OFR vs OS)
## Fig. 3h — NES summary bar graph
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)
library(patchwork)
set.seed(1)

# ---- signature to use ----
lipid_age_sig <- if (exists("new_lipid_age_signature")) {
  new_lipid_age_signature
} else if (exists("signature34")) {
  signature34
} else {
  stop("No lipid aging signature found: define `new_lipid_age_signature` or `signature34`.")
}

# ---- restrict to groups we need & keep only lipid columns + Sample/Group ----
stopifnot(all(c("Sample","Group") %in% names(lip_log)))
lipid_log_filtered <- lip_log %>%
  dplyr::filter(Group %in% c("OS","OV","OIR","OFR")) %>%
  dplyr::select(Sample, Group, dplyr::all_of(lipid_cols))

# ---- helper: compute per-lipid stats and ranks for a two-group contrast ----
.compute_stats_and_ranks <- function(df, g_high, g_low, logfc_name = "LogFC") {
  sub <- df %>% dplyr::filter(Group %in% c(g_high, g_low))
  # per-lipid Welch t-test + logFC (high - low)
  stats <- purrr::map_dfr(lipid_cols, function(lip) {
    vals <- sub[[lip]]
    grp  <- sub$Group
    tt   <- t.test(vals ~ grp)
    tibble::tibble(
      Lipid = lip,
      Mean_high = mean(vals[grp == g_high], na.rm = TRUE),
      Mean_low  = mean(vals[grp == g_low ], na.rm = TRUE),
      !!logfc_name := Mean_high - Mean_low,
      P_value    = tt$p.value
    )
  }) %>%
    dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
    dplyr::arrange(FDR)
  
  # ranking vector: -log10(p) * logFC
  rank_df <- stats %>%
    dplyr::filter(!is.na(.data[[logfc_name]]), !is.na(P_value)) %>%
    dplyr::mutate(
      p_safe = pmax(P_value, 1e-300),
      stat   = -log10(p_safe) * .data[[logfc_name]]
    ) %>%
    dplyr::group_by(Lipid) %>%                  # in case of duplicates
    dplyr::summarise(stat = mean(stat), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(stat))
  
  ranks <- stats::setNames(rank_df$stat, rank_df$Lipid)
  list(stats = stats, ranks = ranks)
}

# ---- helper: enrichment (top) + ranked-stat bars (bottom) for one contrast ----
.make_lsea_panels <- function(ranks, pathway, title_text) {
  # enrichment curve (fgsea provides black ticks for set members)
  p_top <- fgsea::plotEnrichment(pathway, ranks) +
    ggtitle(title_text) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7))
  
  # bottom panel (red to blue across ranked stat; vertical dotted line at sign flip)
  df_bottom <- tibble::tibble(
    Rank  = seq_along(ranks),
    Lipid = names(ranks),
    Stat  = as.numeric(ranks)
  )
  last_pos_rank <- max(which(df_bottom$Stat > 0), na.rm = TRUE)
  if (!is.finite(last_pos_rank)) last_pos_rank <- NA_real_
  flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA_real_
  
  p_bottom <- ggplot(df_bottom, aes(x = Rank, y = Stat, fill = Stat)) +
    geom_col() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Ranked\nstatistic") +
    { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
    labs(x = "Lipid rank", y = "Statistic") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.margin = margin(2, 2, 2, 2))
  
  p_top / p_bottom + plot_layout(heights = c(2, 1))
}

# ---- run the three contrasts (skip gracefully if a group is missing) ----
contrasts <- list(
  `OLD PWR vs OLD SED`      = c("OV",  "OS"),
  `OLD PWR IRAP vs OLD SED` = c("OIR", "OS"),
  `OLD PWR FRAP vs OLD SED` = c("OFR", "OS")
)

results_list <- list()
panels_list  <- list()
nes_rows     <- list()

for (nm in names(contrasts)) {
  g_high <- contrasts[[nm]][1]
  g_low  <- contrasts[[nm]][2]
  if (!all(c(g_high, g_low) %in% lipid_log_filtered$Group)) {
    message("Skipping contrast ", nm, " (missing group in data).")
    next
  }
  
  out <- .compute_stats_and_ranks(lipid_log_filtered, g_high, g_low, logfc_name = "LogFC")
  ranks <- out$ranks
  
  # fgsea on the 34-lipid signature
  lip_sets <- list(AgingSignature = lipid_age_sig)
  fg <- fgsea::fgsea(pathways = lip_sets, stats = ranks, nperm = 1000)
  
  # save NES for Fig. 3h
  nes_rows[[nm]] <- tibble::tibble(
    Contrast = nm,
    NES      = fg$NES[match("AgingSignature", fg$pathway)]
  )
  
  # make the two-panel plot for this contrast
  panels_list[[nm]] <- .make_lsea_panels(
    ranks   = ranks,
    pathway = lip_sets$AgingSignature,
    title_text = nm
  )
  
  results_list[[nm]] <- list(stats = out$stats, ranks = ranks, fgsea = fg)
}

# ---- Fig. 3g: arrange panels left-to-right (only those that ran) ----
panels_present <- panels_list[!vapply(panels_list, is.null, logical(1))]
if (length(panels_present) == 0) stop("No contrasts could be run (check groups in `lip_log`).")

fig3g <- Reduce(`|`, panels_present) + plot_layout(guides = "collect") &
  theme(plot.title = element_text(hjust = 0.5))

print(fig3g)
# ggsave("Fig3g_three_LSEA_panels.pdf", fig3g, width = 12, height = 6.5, device = cairo_pdf)

# ---- Fig. 3h: NES summary (horizontal so 'left = youthful' for negative NES) ----
nes_df <- bind_rows(nes_rows) %>%
  mutate(Direction = ifelse(NES < 0, "Towards youthful", "Away from youthful"))

fig3h <- ggplot(nes_df, aes(x = NES, y = factor(Contrast, levels = names(contrasts)))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = Direction), width = 0.6) +
  scale_fill_manual(values = c("Towards youthful" = "#2E8B57", "Away from youthful" = "#A00000")) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL,
       title = "Summary of NES (center line = OLD SED)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank())

print(fig3h)
# ggsave("Fig3h_NES_summary.pdf", fig3h, width = 6.5, height = 3.8, device = cairo_pdf)

# (optional) Inspect NES values quickly
nes_df
#------------------------------------------------------------------#

## =========================
## Fig. 3i — Volcano plot of metabolome (OLD SED vs YNG SED)
## =========================
library(tidyverse)
library(ggrepel)

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
  dplyr::mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(adj_p)

# 4) Define Metabolite Aging Signature (p < 0.05, unadjusted, per spec)
metabo_age_signature <- metabo_stats_osys %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(p_value) %>%
  dplyr::pull(Metabolite)

# (optional) save for reuse
# readr::write_csv(tibble::tibble(Metabolite = metabo_age_signature),
#                  "Metabolite_Aging_Signature_OSvsYS.csv")

# 5) Volcano plot
signif_metabs <- metabo_stats_osys %>% dplyr::filter(Metabolite %in% metabo_age_signature)

metabo_age_signature

# Symmetric x-limits (robust to outliers)
xlim_val <- max(1, stats::quantile(abs(metabo_stats_osys$log2FC), 0.99, na.rm = TRUE))

fig3i_volcano <- ggplot(metabo_stats_osys, aes(x = log2FC, y = -log10(p_value))) +
  geom_point(color = "grey65", size = 2) +
  # overlay signature in red
  geom_point(data = signif_metabs, color = "red3", size = 2) +
  # label signature (can be noisy; comment out if too dense)
  ggrepel::geom_text_repel(
    data = signif_metabs,
    aes(label = Metabolite),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    seed = 1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  scale_x_continuous(limits = c(-xlim_val, xlim_val)) +
  labs(
    title = "Volcano: age-associated metabolite changes (OLD SED vs YNG SED)",
    x = "log2 Fold Change (OLD SED − YNG SED)",
    y = "−log10(p-value)"
  ) +
  theme_minimal(base_size = 11)

print(fig3i_volcano)
# ggsave("Fig3i_metabolome_volcano_OS_vs_YS.pdf", fig3i_volcano, width = 6.5, height = 6.0)

## =========================
## Fig. 3j,k — 17-metabolite Aging Signature heatmaps + bar plots
## =========================
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(pheatmap)
library(purrr)

# ---------- 0) Ensure imputed_log exists (zeros→half-min per metabolite, then log2) ----------
if (!exists("imputed_log")) {
  metabo_raw <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)
  met_cols <- setdiff(names(metabo_raw), c("Sample","Group"))
  
  # per-metabolite min nonzero (fallback to global)
  min_nonzero <- sapply(metabo_raw[met_cols], function(x) {
    m <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    if (is.infinite(m)) NA_real_ else m
  })
  
  # global min nonzero across ALL metabolites (fallback if a column has no >0)
  flat_vals  <- unlist(metabo_raw[met_cols], use.names = FALSE)
  global_min <- suppressWarnings(min(flat_vals[flat_vals > 0], na.rm = TRUE))
  if (!is.finite(global_min)) global_min <- 1e-6
  
  min_nonzero[is.na(min_nonzero)] <- global_min
  
  # Impute zeros → half the per-metabolite min nonzero (or global if needed)
  imputed_data <- metabo_raw
  for (met in met_cols) {
    x <- imputed_data[[met]]
    x[is.na(x)] <- 0
    if (all(x == 0)) {
      x[x == 0] <- global_min / 2
    } else {
      x[x == 0] <- min_nonzero[[met]] / 2
    }
    imputed_data[[met]] <- x
  }
  
  # Log2 transform
  imputed_log <- imputed_data %>% dplyr::mutate(across(all_of(met_cols), log2))
}

# ---------- 1) Your 17-metabolite signature ----------
if (!exists("metabo_age_signature")) {
  metabo_age_signature <- c(
    "Kynurenic acid","D-Lactose/cellobiose","Isomaltulose","Orotic Acid","Acetyl-L-leucine",
    "D-Maltose","Melibiose","DL-Methionine sulfoxide","N-Acetyl-DL-tryptophan",
    "3-Ureidopropionic acid","Glycerol","1-Aminocyclopropane-1-carboxylic acid",
    "B-Alanine","D-Sorbitol","3R-hydroxy-isobutyric acid","DL-2-Aminoadipic acid","N-Acetyl-L-asparagine"
  )
}

met_cols <- setdiff(names(imputed_log), c("Sample","Group"))
sig_avail   <- intersect(metabo_age_signature, met_cols)
sig_missing <- setdiff(metabo_age_signature, sig_avail)
if (length(sig_missing)) message("Missing (not in data): ", paste(sig_missing, collapse = ", "))
stopifnot(length(sig_avail) >= 1)

# ---------- 2) OS vs YS stats (on log2 data; OS−YS = mean(log2) diff) + row order ----------
osys <- imputed_log %>%
  filter(Group %in% c("YS","OS")) %>%
  mutate(Group = factor(Group, levels = c("YS","OS"))) %>%
  arrange(Group)

osys_stats <- map_dfr(sig_avail, function(met) {
  vals <- osys[[met]]; grp <- osys$Group
  tt <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
  tibble(
    Metabolite = met,
    Mean_YS = mean(vals[grp == "YS"], na.rm = TRUE),
    Mean_OS = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OS_vs_YS = Mean_OS - Mean_YS,
    P_value = if (is.null(tt)) NA_real_ else tt$p.value
  )
}) %>% arrange(desc(Log2_FC_OS_vs_YS))

metab_row_order <- osys_stats$Metabolite  # fixed row order for both panels

# ---------- 3) Fig. 3j — Heatmap YS vs OS ----------
pal_m <- colorRampPalette(c("#FF7F00","black","#00FFFF"))(100)

hm_j <- osys %>%
  dplyr::select(Sample, Group, all_of(metab_row_order))

mat_j <- hm_j %>% dplyr::select(-Sample, -Group) %>% as.matrix()
rownames(mat_j) <- hm_j$Sample
mat_j_z <- scale(mat_j)

ann_col_j <- data.frame(Group = hm_j$Group); rownames(ann_col_j) <- hm_j$Sample
ann_cols_j <- list(Group = c(YS = "#89CFF0", OS = "#4F4F4F"))

hm_3j <- pheatmap::pheatmap(
  t(mat_j_z),
  annotation_col    = ann_col_j,
  annotation_colors = ann_cols_j,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row = metab_row_order, row_names_side = "left",
  fontsize = 8, color = pal_m, breaks = seq(-1, 1, length.out = 101),
  main = "17-Metabolite Aging Signature: YS (left, blue) vs OS (right, grey)"
)
print(hm_3j)

# ---------- 4) Fig. 3j — Bar plot (OS − YS; black = >0) ----------
bar_3j <- osys_stats %>%
  mutate(Metabolite = factor(Metabolite, levels = metab_row_order)) %>%
  ggplot(aes(x = Log2_FC_OS_vs_YS, y = fct_rev(Metabolite))) +
  geom_col(fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "log2 Fold Change (OLD SED − YNG SED)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank())
print(bar_3j)
# ggsave("Fig3j_metabolite_heatmap_YSvsOS.pdf", width = 4.25, height = 5.5)
# ggsave("Fig3j_metabolite_bar_OSminusYS.pdf", bar_3j, width = 4.25, height = 5.5)

# ---------- 5) OV vs OS stats (OV − OS) ----------
ovos <- imputed_log %>%
  filter(Group %in% c("OV","OS")) %>%
  mutate(Group = factor(Group, levels = c("OV","OS"))) %>%
  arrange(Group)

ovos_stats <- map_dfr(sig_avail, function(met) {
  vals <- ovos[[met]]; grp <- ovos$Group
  tt <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
  tibble(
    Metabolite = met,
    Mean_OV = mean(vals[grp == "OV"], na.rm = TRUE),
    Mean_OS = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OV_vs_OS = Mean_OV - Mean_OS,
    P_value = if (is.null(tt)) NA_real_ else tt$p.value
  )
}) %>% arrange(match(Metabolite, metab_row_order))

# ---------- 6) Fig. 3k — Heatmap OV vs OS (same row order as 3j) ----------
hm_k <- ovos %>%
  dplyr::select(Sample, Group, all_of(metab_row_order))

mat_k <- hm_k %>% dplyr::select(-Sample, -Group) %>% as.matrix()
rownames(mat_k) <- hm_k$Sample
mat_k_z <- scale(mat_k)

ann_col_k <- data.frame(Group = hm_k$Group); rownames(ann_col_k) <- hm_k$Sample
ann_cols_k <- list(Group = c(OV = "#98FB98", OS = "#4F4F4F"))

hm_3k <- pheatmap::pheatmap(
  t(mat_k_z),
  annotation_col    = ann_col_k,
  annotation_colors = ann_cols_k,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row = metab_row_order, row_names_side = "left",
  fontsize = 8, color = pal_m, breaks = seq(-1, 1, length.out = 101),
  main = "17-Metabolite Aging Signature: OLD PWR (left, green) vs OLD SED (right, grey)"
)
print(hm_3k)

# ---------- 7) Fig. 3k — Bar plot (OV − OS; red = <0, black = >0) ----------
bar_3k <- ovos_stats %>%
  mutate(
    Metabolite = factor(Metabolite, levels = metab_row_order),
    Direction  = ifelse(Log2_FC_OV_vs_OS < 0, "Down_in_OV", "Up_in_OV")
  ) %>%
  ggplot(aes(x = Log2_FC_OV_vs_OS, y = fct_rev(Metabolite), fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = c(Down_in_OV = "red3", Up_in_OV = "black")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "log2 Fold Change (OLD PWR − OLD SED)", y = NULL, fill = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none")
print(bar_3k)
# ggsave("Fig3k_metabolite_heatmap_OVvsOS.pdf", width = 4.25, height = 5.5)
# ggsave("Fig3k_metabolite_bar_OVminusOS.pdf", bar_3k, width = 4.25, height = 5.5)
#----------------------------------------------------------------------------------#

## =========================
## Fig. 3l — Heatmap across all OLD groups (metabolomics)
## Rows: 17-metabolite Aging Signature; values = z-scored per metabolite
## Blue = low; Red = high; units of raw abundance not confirmed
## =========================

library(dplyr)
library(tidyr)
library(pheatmap)
library(stringr)

# helper: sort sample names like OS1..OS8 numerically by suffix
.sort_by_suffix_num <- function(x) {
  ord <- suppressWarnings(as.numeric(sub("^\\D+", "", x)))
  x[order(ord, na.last = TRUE)]
}

# desired group order and colors (matching your lipidomics)
grp_order <- c("OS", "OV", "OIR", "OFR")
grp_cols  <- c(OS = "#4F4F4F", OV = "#98FB98", OIR = "#66CDAA", OFR = "#2E8B57")

# palette: blue -> white -> red
pal_3l <- colorRampPalette(c("#1F77B4", "white", "#D62728"))(100)

# ---------- 0) Ensure 'imputed_log' exists (zeros→half-min per metabolite, then log2) ----------
if (!exists("imputed_log")) {
  metabo_raw <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)
  met_cols <- setdiff(names(metabo_raw), c("Sample","Group"))
  # per-metabolite min nonzero with global fallback
  min_nonzero <- sapply(metabo_raw[met_cols], function(x) {
    m <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    if (is.infinite(m)) NA_real_ else m
  })
  flat_vals  <- unlist(metabo_raw[met_cols], use.names = FALSE)
  global_min <- suppressWarnings(min(flat_vals[flat_vals > 0], na.rm = TRUE))
  if (!is.finite(global_min)) global_min <- 1e-6
  min_nonzero[is.na(min_nonzero)] <- global_min
  
  imputed_data <- metabo_raw
  for (met in met_cols) {
    x <- imputed_data[[met]]
    x[is.na(x)] <- 0
    if (all(x == 0)) x[x == 0] <- global_min / 2 else x[x == 0] <- min_nonzero[[met]] / 2
    imputed_data[[met]] <- x
  }
  imputed_log <- imputed_data %>% dplyr::mutate(across(all_of(met_cols), log2))
}

# ---------- 1) Ensure we have the 17-metabolite signature ----------
if (!exists("metabo_age_signature")) {
  # build from OS vs YS (unadjusted p < 0.05)
  met_cols <- setdiff(names(imputed_log), c("Sample","Group"))
  osys_df <- imputed_log %>%
    dplyr::filter(Group %in% c("OS","YS")) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(met_cols), names_to = "Metabolite", values_to = "Abundance")
  metabo_stats_osys <- osys_df %>%
    dplyr::group_by(Metabolite) %>%
    dplyr::summarise(
      p_value = tryCatch(t.test(Abundance ~ Group)$p.value, error = function(e) NA_real_),
      mean_OS = mean(Abundance[Group == "OS"], na.rm = TRUE),
      mean_YS = mean(Abundance[Group == "YS"], na.rm = TRUE),
      log2FC  = log2(mean_OS + 1e-6) - log2(mean_YS + 1e-6),
      .groups = "drop"
    ) %>%
    dplyr::arrange(p_value)
  metabo_age_signature <- metabo_stats_osys %>%
    dplyr::filter(p_value < 0.05) %>%
    dplyr::pull(Metabolite) %>%
    unique()
}
# keep only those present in data
met_cols <- setdiff(names(imputed_log), c("Sample","Group"))
sig_avail <- intersect(metabo_age_signature, met_cols)
stopifnot(length(sig_avail) >= 1)

# ---------- 2) Build a dataframe of OLD groups only ----------
# If Group is incomplete, derive from Sample prefix
sample_groups <- imputed_log %>%
  dplyr::distinct(Sample, Group) %>%
  dplyr::mutate(Group_derived = dplyr::case_when(
    stringr::str_starts(Sample, "OS")  ~ "OS",
    stringr::str_starts(Sample, "OV")  ~ "OV",
    stringr::str_starts(Sample, "OIR") ~ "OIR",
    stringr::str_starts(Sample, "OFR") ~ "OFR",
    TRUE ~ as.character(Group)
  ))

old_groups_found <- intersect(grp_order, unique(sample_groups$Group_derived))
stopifnot(length(old_groups_found) > 0)

hm_df <- imputed_log %>%
  dplyr::select(Sample, Group, dplyr::all_of(sig_avail)) %>%
  dplyr::left_join(sample_groups %>% dplyr::select(Sample, Group_derived), by = "Sample") %>%
  dplyr::mutate(Group2 = factor(Group_derived, levels = old_groups_found)) %>%
  dplyr::filter(!is.na(Group2)) %>%
  dplyr::select(Sample, Group = Group2, dplyr::all_of(sig_avail))

# order samples: OS → OV → OIR → OFR; numeric within each group
samp_order <- unlist(lapply(old_groups_found, function(g) {
  .sort_by_suffix_num(hm_df$Sample[hm_df$Group == g])
}))
hm_df <- hm_df %>%
  dplyr::mutate(Sample = factor(Sample, levels = samp_order)) %>%
  dplyr::arrange(Sample)

# ---------- 3) Z-score per metabolite across samples and plot ----------
mat_allold <- hm_df %>% dplyr::select(-Sample, -Group) %>% as.matrix()
rownames(mat_allold) <- hm_df$Sample

# Fix row order (metabolites) for display: keep signature order as provided
row_order <- sig_avail

# z-score per metabolite (columns) across samples
mat_allold <- mat_allold[, row_order, drop = FALSE]
mat_allold_z <- scale(mat_allold)
mat_allold_z[!is.finite(mat_allold_z)] <- 0  # in case of constant metabolites

ann_col <- data.frame(Group = hm_df$Group)
rownames(ann_col) <- hm_df$Sample
ann_cols <- list(Group = grp_cols[names(grp_cols) %in% levels(ann_col$Group)])

hm_3l <- pheatmap::pheatmap(
  t(mat_allold_z),                           # rows = metabolites (17)
  annotation_col    = ann_col,
  annotation_colors = ann_cols,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row   = row_order, row_names_side = "left",
  fontsize = 8, color = pal_3l, breaks = seq(-1, 1, length.out = 101),
  main = "17-metabolite Aging Signature across OLD groups\n(z-scored per metabolite; units of raw abundance not confirmed)"
)
print(hm_3l)

# Optional: save
# pdf("Fig3l_metabolite_signature_all_OLD_groups.pdf", width = 5.0, height = 6.0)
# print(hm_3l)
# dev.off()
#-------------------------------------------------------------#

## ===== Fig. 3m,n — targeted MSEA with p-values in summary table =====
library(dplyr); library(tidyr); library(ggplot2); library(fgsea); library(patchwork)
library(tibble); library(purrr)
set.seed(1)

# ---- Preconditions ----
stopifnot(exists("imputed_log"), all(c("Sample","Group") %in% names(imputed_log)))
met_cols_all <- setdiff(names(imputed_log), c("Sample","Group"))

# If the 17-metabolite signature isn't available, rebuild from OS vs YS (p < 0.05) on log2 data
if (!exists("metabo_age_signature")) {
  osys_long_tmp <- imputed_log %>%
    dplyr::filter(Group %in% c("OS","YS")) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(met_cols_all),
                        names_to = "Metabolite", values_to = "Abundance")
  ms_tmp <- osys_long_tmp %>%
    dplyr::group_by(Metabolite) %>%
    dplyr::summarise(p_value = tryCatch(t.test(Abundance ~ Group)$p.value,
                                        error = function(e) NA_real_), .groups = "drop") %>%
    dplyr::filter(p_value < 0.05)
  metabo_age_signature <- intersect(ms_tmp$Metabolite, met_cols_all)
}
metabo_age_signature <- intersect(metabo_age_signature, met_cols_all)
stopifnot(length(metabo_age_signature) >= 1)

# Restrict to OLD groups
metab_log_filtered <- imputed_log %>%
  dplyr::filter(Group %in% c("OS","OV","OIR","OFR")) %>%
  dplyr::select(Sample, Group, dplyr::all_of(met_cols_all))

# ---- Helper: compute logFC + p and make ranks = -log10(p) * logFC ----
.compute_stats_and_ranks_met <- function(df, g_high, g_low, logfc_name = "LogFC") {
  sub <- df %>% dplyr::filter(Group %in% c(g_high, g_low))
  stats <- purrr::map_dfr(met_cols_all, function(met) {
    vals <- sub[[met]]; grp <- sub$Group
    tt   <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
    tibble::tibble(
      Metabolite = met,
      Mean_high  = mean(vals[grp == g_high], na.rm = TRUE),
      Mean_low   = mean(vals[grp == g_low ], na.rm = TRUE),
      !!logfc_name := Mean_high - Mean_low,                     # on log2 data
      P_value    = if (is.null(tt)) NA_real_ else tt$p.value
    )
  }) %>%
    dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
    dplyr::arrange(FDR)
  
  rank_df <- stats %>%
    dplyr::filter(is.finite(.data[[logfc_name]]), is.finite(P_value)) %>%
    dplyr::mutate(p_safe = pmax(P_value, 1e-300),
                  stat   = -log10(p_safe) * .data[[logfc_name]]) %>%
    dplyr::group_by(Metabolite) %>%
    dplyr::summarise(stat = mean(stat), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(stat))
  
  ranks <- stats::setNames(rank_df$stat, rank_df$Metabolite)
  list(stats = stats, ranks = ranks)
}

# ---- Helper: enrichment (top) + ranked-stat bars (bottom) ----
.make_msea_panels <- function(ranks, pathway, title_text) {
  p_top <- fgsea::plotEnrichment(pathway, ranks) +
    ggtitle(title_text) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7))
  
  df_bottom <- tibble::tibble(Rank = seq_along(ranks),
                              Metabolite = names(ranks),
                              Stat = as.numeric(ranks))
  last_pos_rank <- max(which(df_bottom$Stat > 0), na.rm = TRUE)
  if (!is.finite(last_pos_rank)) last_pos_rank <- NA_real_
  flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA_real_
  
  p_bottom <- ggplot(df_bottom, aes(x = Rank, y = Stat, fill = Stat)) +
    geom_col() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Ranked\nstatistic") +
    { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
    labs(x = "Metabolite rank", y = "Statistic") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))
  
  p_top / p_bottom + plot_layout(heights = c(2, 1))
}

# ---- Run contrasts ----
contrasts_m <- list(
  `OLD PWR vs OLD SED`      = c("OV",  "OS"),
  `OLD PWR IRAP vs OLD SED` = c("OIR", "OS"),
  `OLD PWR FRAP vs OLD SED` = c("OFR", "OS")
)

panels_list_m <- list()
nes_rows_m    <- list()

for (nm in names(contrasts_m)) {
  g_high <- contrasts_m[[nm]][1]
  g_low  <- contrasts_m[[nm]][2]
  if (!all(c(g_high, g_low) %in% metab_log_filtered$Group)) next
  
  out   <- .compute_stats_and_ranks_met(metab_log_filtered, g_high, g_low, logfc_name = "LogFC")
  ranks <- out$ranks
  
  # fgsea on the 17-metabolite signature
  met_sets <- list(AgingSignature = metabo_age_signature)
  # Keep simple permutation to mirror prior behavior; to use fgseaMultilevel, drop nperm.
  fg <- fgsea::fgsea(pathways = met_sets, stats = ranks, nperm = 1000)
  
  # Save NES + p-values
  nes_rows_m[[nm]] <- fg %>%
    dplyr::filter(pathway == "AgingSignature") %>%
    dplyr::transmute(
      Contrast = nm,
      NES      = NES,
      pval     = pval,
      padj     = padj,
      size     = size
    )
  
  # Panels
  panels_list_m[[nm]] <- .make_msea_panels(
    ranks      = ranks,
    pathway    = met_sets$AgingSignature,
    title_text = nm
  )
}

# ---- Fig. 3m ----
panels_present_m <- panels_list_m[!vapply(panels_list_m, is.null, logical(1))]
stopifnot(length(panels_present_m) > 0)
fig3m <- Reduce(`|`, panels_present_m) + plot_layout(guides = "collect") &
  theme(plot.title = element_text(hjust = 0.5))
print(fig3m)

# ---- Fig. 3n + table with p-values ----
nes_df_m <- dplyr::bind_rows(nes_rows_m) %>%
  dplyr::mutate(
    Direction       = ifelse(NES < 0, "Towards youthful", "Away from youthful"),
    pval_overall_BH = p.adjust(pval, method = "BH")  # BH across the three contrasts
  )

fig3n <- ggplot(nes_df_m, aes(x = NES, y = factor(Contrast, levels = names(contrasts_m)))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = Direction), width = 0.6) +
  scale_fill_manual(values = c("Towards youthful" = "#2E8B57", "Away from youthful" = "#A00000")) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL,
       title = "Metabolite Aging Signature — NES summary (center = OLD SED)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())
print(fig3n)

# Show NES + p-values
nes_df_m %>% dplyr::select(Contrast, NES, pval, padj, pval_overall_BH, size, Direction)
## =======================================================================#


###########################################################################
#=====================================================#
#Figure 4: Consolidated Code
#=====================================================#
##########################################################################

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
## Composite Aging Index. Fig. 4a,b (final)
## -----------------------------

library(dplyr); library(tidyr); library(stringr); library(ggplot2)

# 1) Harmonize: keep only samples present in all three omics
common_samples <- Reduce(intersect, list(
  colnames(logCPM_symbol),        # RNA (SYMBOL × T_XX)
  colnames(metabo_log),           # metabolite × T_XX
  colnames(lipid_matrix)          # lipid × T_XX
))

stopifnot(length(common_samples) >= 6)  # expect 6–8 per group overall (T_32 may drop)

rna_mat  <- logCPM_symbol[, common_samples, drop = FALSE]
met_mat  <- metabo_log[,     common_samples, drop = FALSE]
lip_mat  <- lipid_matrix[,   common_samples, drop = FALSE]

# 2) Select features per block
rna_features <- intersect(validated_inflamm_geneset.vec, rownames(rna_mat))
met_features <- intersect(aging_metabolite.vec,           rownames(met_mat))
lip_features <- intersect(sig_lipid_features.vec,         rownames(lip_mat))

stopifnot(length(rna_features) == 52)      # or close (if a symbol maps missing)
stopifnot(length(met_features) >= 15)      # expect ~17
stopifnot(length(lip_features) >= 30)      # expect ~34

# 3) Z-score by row WITHIN each block (feature-wise)
z_by_row <- function(m) t(scale(t(m)))

rna_z <- z_by_row(rna_mat[rna_features, , drop = FALSE])
met_z <- z_by_row(met_mat[met_features, , drop = FALSE])
lip_z <- z_by_row(lip_mat[lip_features, , drop = FALSE])

# 4) Per-sample block means, then equal-weight composite
rna_index   <- colMeans(rna_z, na.rm = TRUE)
met_index   <- colMeans(met_z, na.rm = TRUE)
lip_index   <- colMeans(lip_z, na.rm = TRUE)

composite_index <- (rna_index + met_index + lip_index) / 3

# 5) Assemble long table with 5 group labels
#    (Map IRAP/FRAP onto the “OLD IRAP / OLD FRAP” display names you want)
group_of <- function(sid) {
  if (sid %in% yng_sed_samples)        return("YNG SED")
  if (sid %in% old_sed_samples)        return("OLD SED")
  if (sid %in% old_pwr_samples)        return("OLD PWR")
  if (sid %in% old_pwr_irap_samples)   return("OLD IRAP")
  if (sid %in% old_pwr_frap_samples)   return("OLD FRAP")
  return(NA_character_)
}

df <- data.frame(
  Sample       = common_samples,
  RNA_Index    = rna_index[common_samples],
  Met_Index    = met_index[common_samples],
  Lipid_Index  = lip_index[common_samples],
  AgingIndex   = composite_index[common_samples],
  Group        = vapply(common_samples, group_of, character(1))
) %>% filter(!is.na(Group))

# 6) Boxplot of the Composite Aging Index (your 5 groups)
df$Group <- factor(df$Group, levels = c("YNG SED","OLD SED","OLD PWR","OLD IRAP","OLD FRAP"))

p_box <- ggplot(df, aes(x = Group, y = AgingIndex, fill = Group)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.9) +
  scale_fill_manual(values = c(
    "YNG SED" = "#1f78b4",
    "OLD SED" = "#7f7f7f",
    "OLD PWR" = "#33a02c",
    "OLD IRAP"= "#6a3d9a",
    "OLD FRAP"= "#ff7f00"
  )) +
  labs(title = "Composite Aging Index",
       y = "Mean z-score (equal-weight RNA / Met / Lipid)", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
print(p_box)

# (Optional) one-liner ANOVA + posthoc
anova_model <- aov(AgingIndex ~ Group, data = df)
summary(anova_model); TukeyHSD(anova_model)
#-----------------------------------------------------#

#=====================================================#
#---------Fig. 4e Spearman Correlation Matrix---------#
#=====================================================#


suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(cluster)
})

#------------------------------------------------------------------#
# Step 1: Feature types from rownames of integrated_matrix
#------------------------------------------------------------------#
feature_names <- rownames(integrated_matrix)
feature_types <- case_when(
  grepl("^RNA_", feature_names) ~ "RNA",
  grepl("^Met_", feature_names) ~ "Metabolite",
  grepl("^Lip_", feature_names) ~ "Lipid",
  TRUE ~ "Other"
) %>% factor(levels = c("RNA", "Metabolite", "Lipid"))

#------------------------------------------------------------------#
# Step 2: Force-include MusAge 52 (validated_inflamm_geneset.vec)
#         – already defined in master code as SYMBOLs
#------------------------------------------------------------------#
validated_inflamm_geneset_rna.vec <- paste0("RNA_", validated_inflamm_geneset.vec)
force_include <- intersect(validated_inflamm_geneset_rna.vec, rownames(integrated_matrix))

#------------------------------------------------------------------#
# Step 3: Feature selection
#------------------------------------------------------------------#
feature_variances <- apply(integrated_matrix, 1, var, na.rm = TRUE)

# RNA: top 1,000 variable (excluding RIKENs)
rna_features <- rownames(integrated_matrix)[feature_types == "RNA"]
rna_no_riken <- rna_features[!grepl("Rik$", rna_features, ignore.case = TRUE)]
top_rna <- names(sort(feature_variances[rna_no_riken], decreasing = TRUE))[1:1000]

# Keep all metabolites and lipids
all_met <- rownames(integrated_matrix)[feature_types == "Metabolite"]
all_lip <- rownames(integrated_matrix)[feature_types == "Lipid"]

# Union set
selected_features <- unique(c(top_rna, all_met, all_lip, force_include))
selected_features <- intersect(selected_features, rownames(integrated_matrix))

# Subset matrix
integrated_matrix_sub <- integrated_matrix[selected_features, ]

# Update feature types
feature_names <- rownames(integrated_matrix_sub)
feature_types <- case_when(
  grepl("^RNA_", feature_names) ~ "RNA",
  grepl("^Met_", feature_names) ~ "Metabolite",
  grepl("^Lip_", feature_names) ~ "Lipid",
  TRUE ~ "Other"
) %>% factor(levels = c("RNA", "Metabolite", "Lipid"))

#------------------------------------------------------------------#
# Step 4: Per-omics z-scoring
#------------------------------------------------------------------#
integrated_matrix_z <- integrated_matrix_sub
for (otype in levels(feature_types)) {
  idx <- which(feature_types == otype)
  integrated_matrix_z[idx, ] <- t(scale(t(integrated_matrix_sub[idx, ])))
}

#------------------------------------------------------------------#
# Step 5: Spearman correlation + clustering
#------------------------------------------------------------------#
cor_matrix <- cor(t(integrated_matrix_z), method = "spearman")
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")
d <- dist(1 - cor_matrix)

sil_scores <- sapply(2:10, function(k) {
  cluster_assign <- cutree(hc, k)
  sil <- silhouette(cluster_assign, d)
  mean(sil[, 3])
})

optimal_k <- if (length(unique(sil_scores)) == 1 || which.max(sil_scores) == 1) {
  10
} else {
  which.max(sil_scores)
}
cat("Chosen number of clusters:", optimal_k, "\n")

cluster_assign <- cutree(hc, k = optimal_k)

#------------------------------------------------------------------#
# Step 6: Heatmap
#------------------------------------------------------------------#
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
  annotation_legend_param = list(title_gp = gpar(fontsize = 10),
                                 labels_gp = gpar(fontsize = 9))
)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

png("Spearman_corr_heatmap.png", width = 2000, height = 2000, res = 300)
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

#------------------------------------------------------------------#
# Step 7: Output tables
#------------------------------------------------------------------#
cor_df <- as.data.frame(as.table(cor_matrix)) %>%
  filter(Var1 != Var2) %>%
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
#------------------------------------------------------------------------------#

#=================================================================#
#----- Fig. 4f Ridgeplot of Ceramides (OV, OIR, OFR vs OS mean) --#
#=================================================================#

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggridges)
})

#-------------------------------------------------------------#
# Step 1: Define ceramide set
#-------------------------------------------------------------#
ceramide_vec <- OV_OS_stats %>%
  filter(str_starts(Lipid, "Cer")) %>%
  pull(Lipid)

#-------------------------------------------------------------#
# Step 2: Prepare wide matrix of ceramides
#-------------------------------------------------------------#
lipid_long <- lipid_log_filtered %>%
  pivot_longer(
    cols = -c(Sample, Group),
    names_to = "Lipid",
    values_to = "log_abundance"
  )

lipid_wide <- lipid_long %>%
  filter(Lipid %in% ceramide_vec) %>%
  pivot_wider(
    names_from = Sample,
    values_from = log_abundance
  )

#-------------------------------------------------------------#
# Step 3: Compute OS mean and logFC vs OS
#-------------------------------------------------------------#
os_cols <- grep("^OS", names(lipid_wide), value = TRUE)

cer_df <- lipid_wide %>%
  mutate(OS_mean = rowMeans(across(all_of(os_cols)), na.rm = TRUE))

# Identify comparison groups
logfc_cols <- grep("^(OV|OIR|OFR)", names(cer_df), value = TRUE)

# Subtract OS_mean to get logFC
cer_df_logfc <- cer_df %>%
  mutate(across(all_of(logfc_cols), ~ . - OS_mean, .names = "LogFC_{col}"))

#-------------------------------------------------------------#
# Step 4: Reshape for plotting
#-------------------------------------------------------------#
cer_df_long <- cer_df_logfc %>%
  dplyr::select(Lipid, starts_with("LogFC_")) %>%
  pivot_longer(
    cols = starts_with("LogFC_"),
    names_to = "Sample",
    names_prefix = "LogFC_",
    values_to = "LogFC"
  ) %>%
  mutate(Group = case_when(
    str_detect(Sample, "^OV")  ~ "OV",
    str_detect(Sample, "^OIR") ~ "OIR",
    str_detect(Sample, "^OFR") ~ "OFR"
  ))

# Order ceramides by mean logFC
cer_df_long <- cer_df_long %>%
  group_by(Lipid) %>%
  mutate(Lipid = factor(Lipid, levels = unique(Lipid[order(-mean(LogFC, na.rm = TRUE))])))

# Order groups consistently
cer_df_long <- cer_df_long %>%
  mutate(Group = factor(Group, levels = c("OFR", "OIR", "OV")))

#-------------------------------------------------------------#
# Step 5: Ridgeplot
#-------------------------------------------------------------#
ceramide_ridge_plot <- ggplot(cer_df_long, aes(x = LogFC, y = Lipid, fill = Group)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.2,
    rel_min_height = 0.01,
    color = "white"
  ) +
  scale_fill_manual(
    values = c("OV" = "#1b9e77",
               "OIR" = "#cc79a7",
               "OFR" = "#9e1f63")
  ) +
  labs(
    title = "LogFC of Ceramides Compared to OS",
    x = "LogFC vs OS Mean",
    y = "Ceramide Species"
  ) +
  xlim(-2, 2) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 5)
  )

# Save
pdf("ceramide_ridge_plot.pdf", width = 2.5, height = 6)
print(ceramide_ridge_plot)
dev.off()
#-----------------------------------------------------#

#==============================================================#
#   FIG 4g: Heatmaps (TRI) + logFC barplots for MusAge set     #
#==============================================================#

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(pheatmap); library(ggplot2)
})

# Helper: function to make MusAge heatmap
make_musage_heatmap <- function(tri_wide, geneset, cols_keep, main_title, pdf_file) {
  log_df <- tri_wide %>%
    dplyr::select(Gene, Ensembl_ID, all_of(cols_keep)) %>%
    mutate(across(-c(Gene, Ensembl_ID), ~ log2(.x + 1)))
  
  z_df <- log_df %>%
    column_to_rownames("Gene") %>%
    dplyr::select(-Ensembl_ID) %>%
    as.matrix() %>%
    t() %>% scale(center = TRUE, scale = TRUE) %>% t()
  
  z_mat <- z_df[rownames(z_df) %in% geneset, , drop = FALSE]
  
  custom_palette <- colorRampPalette(c("blue","white","red"))(100)
  max_abs <- 1.25
  custom_breaks <- seq(-max_abs, max_abs, length.out = 101)
  
  hm <- pheatmap(
    z_mat, color = custom_palette, breaks = custom_breaks,
    cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = TRUE,
    border_color = NA, main = main_title
  )
  
  pdf(pdf_file, width = 4.25, height = 5)
  print(hm)
  dev.off()
}

# Helper: function to make MusAge logFC barplot
make_musage_barplot <- function(de_table, geneset, xlab, pdf_file) {
  plot_df <- de_table %>%
    filter(gene_name %in% geneset) %>%
    arrange(gene_name) %>%
    mutate(gene_name = factor(gene_name, levels = rev(gene_name)),
           color = ifelse(logFC < 0, "red", "black"))
  
  p <- ggplot(plot_df, aes(x = logFC, y = gene_name, fill = color)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_identity() +
    theme_minimal() +
    labs(x = xlab, y = NULL) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())
  
  pdf(pdf_file, width = 4.25, height = 5)
  print(p)
  dev.off()
}

#--------------------#
# 10m vs 30m (aging)
#--------------------#
cols_10v30 <- c(paste0("10m_",1:6), paste0("30m_",1:6))
make_musage_heatmap(tri_wide, validated_inflamm_geneset.vec,
                    cols_keep = cols_10v30,
                    main_title = "MusAge: 10m vs 30m TRI",
                    pdf_file = "heatmap_10m_vs_30m.pdf")

make_musage_barplot(old_adlib, validated_inflamm_geneset.vec,
                    xlab = "log2 Fold Change (10m vs 30m)",
                    pdf_file = "barplot_10m_vs_30m.pdf")

#--------------------#
# 30m+CR vs 30m
#--------------------#
cols_30vCR <- c(paste0("30m_",1:6), paste0("30m+CR_",1:6))
make_musage_heatmap(tri_wide, validated_inflamm_geneset.vec,
                    cols_keep = cols_30vCR,
                    main_title = "MusAge: 30m+CR vs 30m TRI",
                    pdf_file = "heatmap_30mCR_vs_30m.pdf")

make_musage_barplot(old_cr, validated_inflamm_geneset.vec,
                    xlab = "log2 Fold Change (30m+CR vs 30m)",
                    pdf_file = "barplot_30mCR_vs_30m.pdf")
#-----------------------------------------------------------#

#==============================================================#
#   FIG 4h,i: Targeted GSEA plots (10m vs 30m; 30m+CR vs 30m)  #
#==============================================================#

suppressPackageStartupMessages({ library(clusterProfiler); library(enrichplot) })

# Rank genes by -log10(pval)*sign(logFC)
make_ranked_list <- function(df) {
  df %>%
    mutate(stat = -log10(pval) * sign(logFC)) %>%
    arrange(desc(stat)) %>%
    { setNames(.$stat, .$gene_name) }
}

# GSEA plot function
make_gsea_plot <- function(geneList, geneset, title, pdf_file) {
  gsea_plot <- plotEnrichment(geneset, geneList) + ggtitle(title)
  pdf(pdf_file, width = 4.25, height = 5.5)
  print(gsea_plot)
  dev.off()
}

#--------------------#
# 10m vs 30m (aging)
#--------------------#
old_adlib_geneList <- make_ranked_list(old_adlib)
make_gsea_plot(old_adlib_geneList,
               validated_inflamm_geneset.vec,
               "Effect of Age (10m vs 30m) on MusAge Gene Set",
               "gsea_10m_vs_30m.pdf")

#--------------------#
# 30m+CR vs 30m
#--------------------#
old_cr_geneList <- make_ranked_list(old_cr)
make_gsea_plot(old_cr_geneList,
               validated_inflamm_geneset.vec,
               "Effect of CR (30m+CR vs 30m) on MusAge Gene Set",
               "gsea_30mCR_vs_30m.pdf")


