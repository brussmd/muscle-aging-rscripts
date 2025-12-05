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

setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")

suppressPackageStartupMessages({
  library(readxl)          # read_xlsx
  library(dplyr)           # data manipulation
  library(tidyr)           # pivot_wider, separate_rows
  library(stringr)         # str_remove, string cleanup
  library(tibble)          # deframe, tibble operations
  library(org.Mm.eg.db)    # annotation mapping (mouse)
  library(AnnotationDbi)   # mapIds
  library(biomaRt)         # mouse → human gene conversion
  library(pheatmap)        # heatmaps
  library(ggplot2)         # plotting expression summary
  library(vroom)           # fast loading of GTEx file
  library(pbapply)         # progress bar with apply
  library(broom)           # tidy stats
  library(logistf)         # Firth logistic regression
  library(glmnet)          # Ridge logistic regression
})
############################################################
# FIGURE 1g – GO Pathway Dotplot (OLD SED VEH vs YNG SED VEH)
############################################################
# Input: 20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx
# Steps: GSEA (GO) -> simplify -> dotplot

# 1) Load edgeR *_GENE_* xlsx and build ranked vector
oldsedveh_v_yngsedveh_genes_16sep <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col   = "OLD_SED_VEH-YNG_SED_VEH_FDR"
)

head(oldsedveh_v_yngsedveh_genes_16sep)

# quick sanity check (expected ~701 DEG at FDR<0.05 per your note)
oldsedveh_v_yngsedveh_genes_16sep %>% 
  filter(FDR < 0.05 & logFC >0.3) %>%
  arrange(logFC) %>%
  pull(Symbol)


# 1) Load edgeR *_GENE_* xlsx and build ranked vector
yngpwrveh_v_yngsedveh_genes_16sep <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_YNG_PWR_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "YNG_PWR_VEH-YNG_SED_VEH_logFC",
  FDR_col   = "YNG_PWR_VEH-YNG_SED_VEH_FDR"
)

yngpwrveh_v_yngsedveh_genes_16sep %>%
  filter(Symbol %in% mouse_targets_1to1)

test_chea3_genes <- yngpwrveh_v_yngsedveh_genes_16sep %>%
  filter(FDR < 0.05 & abs(logFC) >0.4) %>%
  pull(Symbol)

length(test_chea3_genes)
write.csv(test_chea3_genes, "test_chea3_genes.csv",row.names = FALSE)


upreg_pwr_genes <- yngpwrveh_v_yngsedveh_genes_16sep %>%
  filter(FDR < 0.05 & logFC > 0.3) %>%
  pull(Symbol)

length(upreg_pwr_genes)
write.csv(upreg_pwr_genes, "upreg_pwr_genes.csv",row.names = FALSE)

downreg_pwr_genes <- yngpwrveh_v_yngsedveh_genes_16sep %>%
  filter(FDR < 0.05 & logFC < -0.3) %>%
  pull(Symbol)

length(downreg_pwr_genes)
write.csv(downreg_pwr_genes, "downreg_pwr_genes.csv",row.names = FALSE)

yngpwrveh_v_yngsedveh_genes_16sep %>%
  filter(Symbol == "Bcl6b")
#---------------------------------------#

head(yngpwrveh_v_yngsedveh_genes_16sep)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(broom)
})

UP <- function(x) toupper(trimws(as.character(x)))
sc_data <- read.delim("rna_single_cell_type.tsv")
# --- sc_data already loaded/renamed earlier ---
# sc_data must have: Gene.name, Cell.type, nTPM (numeric), one row per gene×celltype
# If not numeric yet:
sc_data$nTPM <- suppressWarnings(as.numeric(sc_data$nTPM))

# De-dup (keep max per gene×cell type), upper-case gene names
sc_data <- sc_data |>
  group_by(Gene.name, Cell.type) |>
  summarise(nTPM = max(nTPM, na.rm = TRUE), .groups = "drop") |>
  mutate(Gene.name = UP(Gene.name))

# Your DE table -> universe
fdr_thr <- 0.05; up_thr <- 0.3; down_thr <- -0.3
de <- yngpwrveh_v_yngsedveh_genes_16sep |>
  transmute(SYMBOL = UP(Symbol),
            logFC, FDR,
            is_up   = FDR < fdr_thr & logFC >  up_thr,
            is_down = FDR < fdr_thr & logFC <  down_thr) |>
  distinct(SYMBOL, .keep_all = TRUE)

# --- Define endothelial with robust stats (mean vs median) ---
endo_types <- unique(sc_data$Cell.type[grepl("endothelial", sc_data$Cell.type, ignore.case = TRUE)])
message("HPA endothelial types: ", paste(endo_types, collapse = ", "))

thr_expr <- 5     # min mean nTPM in endothelium
fc_min   <- 1.5   # fold over non-endothelial median

endo_scores <- sc_data |>
  group_by(Gene.name) |>
  summarise(
    endo_mean = mean(nTPM[Cell.type %in% endo_types], na.rm = TRUE),
    non_median= median(nTPM[!(Cell.type %in% endo_types)], na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    endo_mean  = replace_na(endo_mean,  0),
    non_median = replace_na(non_median, 0),
    SYMBOL     = UP(Gene.name),
    is_endothelial = (endo_mean >= thr_expr) & ((endo_mean + 1) / (non_median + 1) >= fc_min)
  ) |>
  dplyr::select(SYMBOL, is_endothelial)

# --- Merge & answer the two questions ---
df <- de |>
  left_join(endo_scores, by = "SYMBOL") |>
  mutate(is_endothelial = ifelse(is.na(is_endothelial), FALSE, is_endothelial))

# Q1a: % of universe that is endothelial
pct_endo <- mean(df$is_endothelial) * 100
cat(sprintf("Endothelial by rule (endo_mean≥%g & ≥%gx non-endo median): %d of %d (%.1f%%)\n",
            thr_expr, fc_min, sum(df$is_endothelial), nrow(df), pct_endo))

# Q1b: Are UP genes enriched in endothelial?
tab_up <- table(Endothelial = df$is_endothelial, UP = df$is_up)
ft_up  <- fisher.test(tab_up)
cat("UP enrichment: OR =", round(unname(ft_up$estimate), 3),
    "CI95% [", paste(round(unname(ft_up$conf.int),3), collapse=", "), "]",
    "p =", signif(ft_up$p.value, 3), "\n")
cat("% endothelial within UP     =", round(mean(df$is_endothelial[df$is_up])   * 100, 1), "%\n")
cat("% endothelial within not-UP =", round(mean(df$is_endothelial[!df$is_up]) * 100, 1), "%\n")

tab_down <- table(Endothelial = df$is_endothelial, DOWN = df$is_down)
fisher.test(tab_down)

# robustly read your 50-gene overlap list (one gene per line or comma-separated)
bcl6b_vec <- toupper(trimws(unlist(strsplit(paste(readLines("upreg_pwr_bcl6b_ovrlp.csv"), collapse=","), "[,;\\s]+"))))
bcl6b_vec <- unique(bcl6b_vec[nchar(bcl6b_vec) > 0])

df_ec <- df[df$is_endothelial, ]                           # endothelial-only universe
df_ec$in_BCL6B <- df_ec$SYMBOL %in% bcl6b_vec

tab_ec_up <- table(BCL6B = df_ec$in_BCL6B, UP = df_ec$is_up)
fisher.test(tab_ec_up)  # enrichment of BCL6B-overlap among UP within EC-only
#-------------------------------------------------------------#
#-------------------------------------------------------------#

UP <- function(x) toupper(trimws(as.character(x)))

# --- Helper to read your 50-gene BCL6B overlap list robustly (1 per line or comma/space separated) ---
read_gene_list <- function(path){
  v <- tryCatch({
    x <- readLines(path, warn = FALSE)
    x <- paste(x, collapse = ",")
    x <- unlist(strsplit(x, "[,;\\t\\s]+"))
    UP(x)
  }, error = function(e) character(0))
  unique(v[nchar(v) > 0])
}

# --- Build the lookup (SYMBOL -> endothelial flag) from your merged 'df' ---
stopifnot(all(c("SYMBOL","is_endothelial","is_up","is_down") %in% names(df)))
endo_lookup <- setNames(df$is_endothelial, df$SYMBOL)

# --- Helper: % endothelial in an input set, restricted to genes present in df/universe ---
pct_endo_in <- function(genes, label){
  g <- unique(UP(genes))
  g_in <- intersect(g, names(endo_lookup))
  tibble(
    set = label,
    n_input = length(g),
    n_in_universe = length(g_in),
    n_endothelial = sum(endo_lookup[g_in], na.rm = TRUE),
    pct_endothelial = ifelse(length(g_in)>0, round(100*mean(endo_lookup[g_in]),1), NA_real_)
  )
}

# -------------------------------
# 1) % endothelial among ALL genes in your universe
# -------------------------------
res_universe <- tibble(
  set = "Universe (all genes in df)",
  n_input = nrow(df),
  n_in_universe = nrow(df),
  n_endothelial = sum(df$is_endothelial),
  pct_endothelial = round(100*mean(df$is_endothelial), 1)
)

# -------------------------------
# 2) % endothelial among UP-regulated genes
# -------------------------------
res_up <- pct_endo_in(upreg_pwr_genes, "UP-regulated (your list)")

# -------------------------------
# 3) % endothelial among DOWN-regulated genes
# -------------------------------
res_down <- pct_endo_in(downreg_pwr_genes, "DOWN-regulated (your list)")

# -------------------------------
# 4) % endothelial among (UP ∩ BCL6B-overlap)
#     i.e., those ~50 genes returned by ChEA3 from your UP set
# -------------------------------
bcl6b_vec <- read_gene_list("upreg_pwr_bcl6b_ovrlp.csv")

# Make sure we only evaluate overlap WITHIN the UP universe
up_universe <- df %>% filter(is_up)
bcl6b_in_up <- intersect(bcl6b_vec, up_universe$SYMBOL)

res_bcl6b_up <- pct_endo_in(bcl6b_in_up, "UP ∩ BCL6B-overlap")

# Collect and print
res <- bind_rows(res_universe, res_up, res_down, res_bcl6b_up)
print(res, n = nrow(res), width = Inf)

cat("\nPretty counts:\n")
apply(res[, c("set","n_in_universe","n_endothelial","pct_endothelial")], 1, function(r)
  cat(sprintf("%-28s : %4s/%-4s endothelial  (%s%%)\n", r[1], r[3], r[2], r[4])))

# -------------------------------
# Enrichment test you actually want:
# Within UP genes, are BCL6B-overlap genes more likely to be endothelial?
# -------------------------------
up_universe <- up_universe %>%
  mutate(in_BCL6B = SYMBOL %in% bcl6b_in_up)

tab_up_within <- table(BCL6B = up_universe$in_BCL6B,
                       Endothelial = up_universe$is_endothelial)
cat("\n2x2 within UP:\n"); print(tab_up_within)

ft_up_within <- fisher.test(tab_up_within)
cat(sprintf("\nWithin-UP enrichment of endothelial for BCL6B-overlap:\nOR = %.3f, 95%% CI [%0.3f, %0.3f], p = %.3g\n",
            unname(ft_up_within$estimate),
            unname(ft_up_within$conf.int[1]),
            unname(ft_up_within$conf.int[2]),
            ft_up_within$p.value))
#------------------------------------------#
#==========================================#

# Read the TSV you just downloaded from ChEA3
tab <- read.delim("Integrated_meanRank_upreg_pwr_bcl6b_ovrlp.tsv", check.names = FALSE)

# Heuristic to find the right columns across libraries
tf_col <- grep("^(TF|Name)$|^TF$|Regulator", names(tab), value = TRUE)[1]
ov_col <- grep("Overlap|Overlapping.*Genes|Genes", names(tab), value = TRUE)[1]

stopifnot(length(tf_col) == 1, length(ov_col) == 1)

# Split the comma-separated overlap column into vectors
library(tidyr); library(dplyr); library(stringr)
ov_long <- tab %>%
  transmute(TF = .data[[tf_col]],
            overlap_raw = .data[[ov_col]]) %>%
  mutate(overlap_raw = ifelse(is.na(overlap_raw), "", overlap_raw)) %>%
  separate_rows(overlap_raw, sep = "[,;\\s]+") %>%
  filter(overlap_raw != "") %>%
  mutate(GENE = toupper(trimws(overlap_raw))) %>%
  distinct(TF, GENE)

# Now you have TF ↔ overlapping gene pairs for *all* TFs
head(ov_long)
#--------------------------#

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

UP <- function(x) toupper(trimws(as.character(x)))

# -- Ensure basic inputs are present --
stopifnot(all(c("SYMBOL","is_up","is_down","is_endothelial","logFC") %in% names(df)))
stopifnot(all(c("TF","GENE") %in% names(ov_long)))

# Keep only the UP universe
up_universe <- df %>% filter(is_up) %>% mutate(SYMBOL = UP(SYMBOL))
up_syms     <- up_universe$SYMBOL
is_endo     <- setNames(up_universe$is_endothelial, up_universe$SYMBOL)

# Uppercase genes in the ChEA3 overlaps and dedupe
ov_long <- ov_long %>% mutate(GENE = UP(GENE)) %>% distinct(TF, GENE)

# --- Scoring function: per TF, how endothelial is its overlap *within UP*? ---
score_tf <- function(tf) {
  g_tf <- ov_long %>% filter(TF == tf) %>% pull(GENE) %>% unique()
  g_tf <- intersect(g_tf, up_syms)                      # stay within UP
  if (length(g_tf) == 0L) {
    return(data.frame(TF=tf, n_overlap=0, pct_EC_overlap=NA,
                      pct_EC_other=NA, OR=NA, CI_low=NA, CI_high=NA, p_value=NA))
  }
  dat <- up_universe %>% mutate(in_TF = SYMBOL %in% g_tf)
  tab <- table(TF = dat$in_TF, EC = dat$is_endothelial)
  # Guard against degenerate tables
  ft  <- tryCatch(fisher.test(tab), error = function(e) NULL)
  if (is.null(ft)) {
    data.frame(TF=tf, n_overlap=sum(dat$in_TF),
               pct_EC_overlap=mean(dat$is_endothelial[dat$in_TF])*100,
               pct_EC_other  =mean(dat$is_endothelial[!dat$in_TF])*100,
               OR=NA, CI_low=NA, CI_high=NA, p_value=NA)
  } else {
    data.frame(TF=tf, n_overlap=sum(dat$in_TF),
               pct_EC_overlap=mean(dat$is_endothelial[dat$in_TF])*100,
               pct_EC_other  =mean(dat$is_endothelial[!dat$in_TF])*100,
               OR=unname(ft$estimate),
               CI_low=unname(ft$conf.int[1]), CI_high=unname(ft$conf.int[2]),
               p_value=ft$p.value)
  }
}

# --- Score all TFs ---
all_tfs <- sort(unique(ov_long$TF))
res_tf  <- do.call(rbind, lapply(all_tfs, score_tf)) %>%
  mutate(p_adj = p.adjust(p_value, method="BH")) %>%
  arrange(p_adj, desc(pct_EC_overlap))

# --- Show BCL6B quickly + top TFs ---
res_bcl6b <- res_tf %>% filter(TF == "BCL6B")
cat("BCL6B summary:\n"); print(res_bcl6b)

cat("\nTop TFs by EC% in overlap (within UP):\n")
print(res_tf %>% arrange(desc(pct_EC_overlap)) %>% head(15), row.names = FALSE)

cat("\nTop TFs by significance (BH adj):\n")
print(res_tf %>% arrange(p_adj) %>% head(15), row.names = FALSE)

# --- Null: is BCL6B’s EC% extreme vs random UP subsets of same size? ---
n_bcl <- res_bcl6b$n_overlap[1]
obs   <- res_bcl6b$pct_EC_overlap[1]
if (is.finite(n_bcl) && n_bcl > 0) {
  set.seed(1)
  reps <- 5000
  null_pcts <- replicate(reps, mean(is_endo[sample(up_syms, n_bcl)]) ) * 100
  emp_p <- mean(null_pcts >= obs)
  cat(sprintf("\nNull (random %d-from-UP): mean=%.1f%%, sd=%.1f, observed BCL6B=%.1f%%, empirical p=%.4g\n",
              n_bcl, mean(null_pcts), sd(null_pcts), obs, emp_p))
  q <- quantile(null_pcts, c(.5,.9,.95,.99))
  print(q)
} else {
  cat("\nBCL6B had 0 overlap within UP; skipping null.\n")
}

# --- Effect size within EC-only: are BCL6B EC genes more shifted? ---
up_ec <- up_universe %>%
  filter(is_endothelial) %>%
  mutate(in_BCL6B = SYMBOL %in% (ov_long %>% filter(TF=="BCL6B") %>% pull(GENE)))

if (sum(up_ec$in_BCL6B) > 0 && sum(!up_ec$in_BCL6B) > 0) {
  w <- wilcox.test(abs(logFC) ~ in_BCL6B, data = up_ec)
  cat(sprintf("\nEC-only |logFC| difference (BCL6B vs other EC-UP): p=%.3g; n_in= %d vs %d\n",
              w$p.value, sum(up_ec$in_BCL6B), sum(!up_ec$in_BCL6B)))
  cat("Median |logFC| (BCL6B EC)   :", median(abs(up_ec$logFC[up_ec$in_BCL6B])), "\n")
  cat("Median |logFC| (other EC)   :", median(abs(up_ec$logFC[!up_ec$in_BCL6B])), "\n")
} else {
  cat("\nNot enough EC genes in one of the groups for |logFC| comparison.\n")
}
#---------------------------------------------------------#

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(broom); library(stringr)
})

UP <- function(x) toupper(trimws(as.character(x)))

# ---------- Inputs expected ----------
# df: data.frame with SYMBOL, is_up, is_down, is_endothelial, logFC (from your pipeline)
# ov_long: data.frame with TF, GENE (parsed from ChEA3 integrated table)

stopifnot(all(c("SYMBOL","is_up","is_down","is_endothelial","logFC") %in% names(df)))
stopifnot(all(c("TF","GENE") %in% names(ov_long)))

# Keep only UP genes universe for the EC-purity comparisons
up_universe <- df %>% filter(is_up) %>% mutate(SYMBOL = UP(SYMBOL))
up_syms     <- up_universe$SYMBOL
is_endo     <- setNames(up_universe$is_endothelial, up_universe$SYMBOL)

# Make sure ChEA3 overlaps are uppercase and unique
ov_long <- ov_long %>% mutate(TF = UP(TF), GENE = UP(GENE)) %>% distinct(TF, GENE)

# ---------- A) EC-purity & null test per TF (within UP) ----------
# Candidate endothelial TFs to compare head-to-head with BCL6B
cand_tfs <- c("BCL6B","SOX7","SOX18","ERG","EPAS1","TAL1","FOXF1","FOXF2","MEOX1","MEOX2","EBF3","LHX6")
cand_tfs <- intersect(cand_tfs, unique(ov_long$TF))  # keep only those present

# Helper to score one TF
score_tf <- function(tf, reps = 5000, seed = 1L) {
  g_tf <- ov_long %>% filter(TF == tf) %>% pull(GENE) %>% unique()
  g_tf <- intersect(g_tf, up_syms)  # stay within UP universe
  n_overlap <- length(g_tf)
  if (n_overlap == 0L) {
    return(tibble(TF=tf, n_overlap=0, pct_EC_overlap=NA_real_, pct_EC_other=NA_real_,
                  OR=NA_real_, CI_low=NA_real_, CI_high=NA_real_, p_value=NA_real_,
                  null_mean=NA_real_, null_sd=NA_real_, emp_p=NA_real_))
  }
  dat <- up_universe %>% mutate(in_TF = SYMBOL %in% g_tf)
  tab <- table(TF = dat$in_TF, EC = dat$is_endothelial)
  ft  <- fisher.test(tab)
  pct_EC_overlap <- mean(dat$is_endothelial[dat$in_TF]) * 100
  pct_EC_other   <- mean(dat$is_endothelial[!dat$in_TF]) * 100
  
  # Null: random n_overlap genes from UP; empirical p for EC%
  set.seed(seed)
  null_pcts <- replicate(reps, mean(is_endo[sample(up_syms, n_overlap)]) ) * 100
  emp_p <- mean(null_pcts >= pct_EC_overlap)
  
  tibble(
    TF = tf,
    n_overlap = n_overlap,
    pct_EC_overlap = pct_EC_overlap,
    pct_EC_other   = pct_EC_other,
    OR = unname(ft$estimate),
    CI_low = unname(ft$conf.int[1]),
    CI_high= unname(ft$conf.int[2]),
    p_value = ft$p.value,
    null_mean = mean(null_pcts),
    null_sd   = sd(null_pcts),
    emp_p     = emp_p
  )
}

res_A <- bind_rows(lapply(cand_tfs, score_tf)) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, desc(pct_EC_overlap))

cat("\n=== A) EC purity within UP (head-to-head) ===\n")
print(res_A, n = nrow(res_A), width = Inf)

UP <- function(x) toupper(trimws(as.character(x)))

# --- Inputs assumed already in memory ---
# df: data.frame with SYMBOL, is_up (logical), is_down (logical), is_endothelial (logical), logFC
# ov_long: data.frame with TF, GENE (ChEA3 overlaps)
stopifnot(all(c("SYMBOL","is_up","is_down","is_endothelial") %in% names(df)))
stopifnot(all(c("TF","GENE") %in% names(ov_long)))

# Normalize case & dedupe overlaps
ov_long <- ov_long %>% mutate(TF = UP(TF), GENE = UP(GENE)) %>% distinct(TF, GENE)

# Universe for tests
X <- df %>% transmute(SYMBOL = UP(SYMBOL), is_up, is_endothelial)

# A small EC-TF panel for head-to-head comparisons (present in your ChEA3 table)
panel <- c("BCL6B","SOX7","SOX18","ERG","EPAS1","TAL1","FOXF1","FOXF2","MEOX1","MEOX2","EBF3","LHX6")
panel <- intersect(panel, unique(ov_long$TF))

# Helper to build a 2x2 with fixed levels (avoids empty-level issues)
mk2x2 <- function(a, b, rn = c(FALSE, TRUE), cn = c(FALSE, TRUE)) {
  tab <- table(a, b)
  M <- matrix(0L, nrow = 2, ncol = 2, dimnames = list(TF = rn, UP = cn))
  if (length(tab)) M[rownames(tab), colnames(tab)] <- tab
  as.table(M)
}

# ----- CMH (Mantel–Haenszel) for each TF: UP ~ TF_overlap | endothelial -----
mh_one_tf <- function(tf) {
  g <- ov_long %>% filter(TF == tf) %>% pull(GENE) %>% unique()
  ind <- X$SYMBOL %in% g  # TF-overlap indicator over the whole universe
  
  # Build 2x2 tables in each stratum
  mat_ec    <- mk2x2(ind[X$is_endothelial],     X$is_up[X$is_endothelial])
  mat_nonec <- mk2x2(ind[!X$is_endothelial],    X$is_up[!X$is_endothelial])
  
  # 2x2x2 array for CMH
  arr <- array(0L, dim = c(2,2,2),
               dimnames = list(TF = c(FALSE, TRUE), UP = c(FALSE, TRUE), Stratum = c("EC","nonEC")))
  arr[, , "EC"]    <- mat_ec
  arr[, , "nonEC"] <- mat_nonec
  
  mh <- mantelhaen.test(arr, correct = FALSE)
  
  # Also report simple within-EC Fisher (what you found striking for BCL6B)
  ft_ec  <- fisher.test(mat_ec)
  ft_non <- fisher.test(mat_nonec)
  
  data.frame(
    TF = tf,
    n_overlap_universe = sum(ind),
    # counts within EC stratum
    EC_overlap_UP     = mat_ec["TRUE","TRUE"],
    EC_overlap_notUP  = mat_ec["TRUE","FALSE"],
    EC_other_UP       = mat_ec["FALSE","TRUE"],
    EC_other_notUP    = mat_ec["FALSE","FALSE"],
    # CMH (adjusted for EC vs non-EC)
    CMH_common_OR = unname(mh$estimate),
    CMH_p        = mh$p.value,
    # within-EC OR/p (unadjusted, just the EC stratum)
    EC_OR        = unname(ft_ec$estimate),
    EC_p         = ft_ec$p.value,
    # within-nonEC OR/p (usually small for EC TFs)
    nonEC_OR     = unname(ft_non$estimate),
    nonEC_p      = ft_non$p.value,
    row.names = NULL
  )
}

res_cmh <- do.call(rbind, lapply(panel, mh_one_tf)) %>%
  arrange(CMH_p)

cat("\n=== CMH: UP ~ TF_overlap | endothelial (adjusted for EC vs non-EC) ===\n")
print(res_cmh, row.names = FALSE)

# ---------------- OPTIONAL: a small, stable regression ----------------
# Use Firth logistic to avoid separation: is_up ~ is_endothelial + TF_BCL6B + EC_TF_count(others)
# (Install logistf the first time)
if (!requireNamespace("logistf", quietly = TRUE)) {
  message("\nInstalling 'logistf' for bias-reduced logistic (one-time)...")
  install.packages("logistf")
}
# Build indicators
tf_indicator <- function(tf, syms) syms %in% (ov_long %>% filter(TF==tf) %>% pull(GENE) %>% unique())
others <- setdiff(panel, "BCL6B")
X2 <- X %>% mutate(
  TF_BCL6B = tf_indicator("BCL6B", SYMBOL),
  EC_TF_COUNT_OTHERS = if (length(others) == 0) 0L else {
    mats <- sapply(others, function(tf) tf_indicator(tf, SYMBOL))
    rowSums(as.matrix(mats))
  }
)
# Drop constant columns if any
keep <- sapply(X2, function(v) length(unique(v)) > 1)
X2 <- X2[, keep, drop = FALSE]

# Firth logistic
fit_firth <- logistf::logistf(is_up ~ is_endothelial + TF_BCL6B + EC_TF_COUNT_OTHERS, data = X2)
sum_firth <- data.frame(term = rownames(fit_firth$ci),
                        OR = exp(coef(fit_firth)),
                        CI_low = exp(fit_firth$ci[,1]),
                        CI_high= exp(fit_firth$ci[,2]),
                        p_value= fit_firth$prob)
cat("\n=== Firth logistic: is_up ~ is_endothelial + TF_BCL6B + EC_TF_COUNT_OTHERS ===\n")
print(sum_firth, row.names = FALSE)
#-------------------------------------------#

# make arr numeric, not integer (0 not 0L)
mk2x2 <- function(a, b, rn = c(FALSE, TRUE), cn = c(FALSE, TRUE)) {
  tab <- table(a, b)
  M <- matrix(0, nrow = 2, ncol = 2, dimnames = list(TF = rn, UP = cn))
  if (length(tab)) M[rownames(tab), colnames(tab)] <- tab
  as.table(M)
}

mh_one_tf <- function(tf) {
  g   <- ov_long %>% filter(TF == tf) %>% pull(GENE) %>% unique()
  ind <- X$SYMBOL %in% g
  
  mat_ec    <- mk2x2(ind[X$is_endothelial],  X$is_up[X$is_endothelial])
  mat_nonec <- mk2x2(ind[!X$is_endothelial], X$is_up[!X$is_endothelial])
  
  arr <- array(0, dim = c(2,2,2),
               dimnames = list(TF = c(FALSE, TRUE), UP = c(FALSE, TRUE), Stratum = c("EC","nonEC")))
  arr[, , "EC"]    <- mat_ec
  arr[, , "nonEC"] <- mat_nonec
  
  mh    <- mantelhaen.test(arr, correct = FALSE)
  ft_ec <- fisher.test(mat_ec); ft_non <- fisher.test(mat_nonec)
  
  data.frame(
    TF=tf, n_overlap_universe=sum(ind),
    EC_overlap_UP=mat_ec["TRUE","TRUE"], EC_overlap_notUP=mat_ec["TRUE","FALSE"],
    EC_other_UP=mat_ec["FALSE","TRUE"],  EC_other_notUP =mat_ec["FALSE","FALSE"],
    CMH_common_OR=unname(mh$estimate), CMH_p=mh$p.value,
    EC_OR=unname(ft_ec$estimate), EC_p=ft_ec$p.value,
    nonEC_OR=unname(ft_non$estimate), nonEC_p=ft_non$p.value
  )
}

# Count how many *other* EC TFs each gene overlaps (cap at 3+ to keep strata roomy)
others <- setdiff(panel, "BCL6B")
other_count <- function(syms){
  if (length(others)==0) return(integer(length(syms)))
  mats <- sapply(others, function(tf) syms %in% (ov_long %>% filter(TF==tf) %>% pull(GENE) %>% unique()))
  rowSums(as.matrix(mats))
}

ec_only <- X %>% filter(is_endothelial) %>%
  mutate(TF_BCL6B = SYMBOL %in% (ov_long %>% filter(TF=="BCL6B") %>% pull(GENE) %>% unique()),
         OTHER_CT = pmin(other_count(SYMBOL), 3),
         OTHER_BIN = factor(OTHER_CT, levels = 0:3,
                            labels = c("0","1","2","3+")))

# Build 2×2 for each OTHER_BIN stratum: (BCL6B overlap) × (UP) | OTHER_BIN
strata <- levels(ec_only$OTHER_BIN)
arr2 <- array(0, dim=c(2,2,length(strata)),
              dimnames=list(BCL6B=c(FALSE,TRUE), UP=c(FALSE,TRUE), Stratum=strata))

for (s in strata) {
  sub <- ec_only %>% filter(OTHER_BIN==s)
  tab <- table(sub$TF_BCL6B, sub$is_up)
  M <- matrix(0, nrow=2, ncol=2, dimnames=list(BCL6B=c(FALSE,TRUE), UP=c(FALSE,TRUE)))
  if (length(tab)) M[rownames(tab), colnames(tab)] <- tab
  arr2[, , s] <- M
}

mh_ec <- mantelhaen.test(arr2, correct=FALSE)
cat(sprintf("\nEC-only CMH (adjusted for other EC-TF load): common OR = %.2f, p = %.3g\n",
            unname(mh_ec$estimate), mh_ec$p.value))

# (Optional) show the strata tables
apply(arr2, 3, function(M) {print(M); ""})

suppressPackageStartupMessages({ library(glmnet) })

# Build a compact model matrix with EC, BCL6B, and OTHER_CT
Mdat <- ec_only %>% dplyr::select(is_up, TF_BCL6B, OTHER_CT)  # EC-only model
M    <- model.matrix(is_up ~ TF_BCL6B + OTHER_CT, data=Mdat)[, -1]
y    <- as.integer(Mdat$is_up)

set.seed(1)
cv   <- cv.glmnet(M, y, family="binomial", alpha=0)  # ridge
b    <- coef(cv, s="lambda.1se")
b["TF_BCL6BTRUE", , drop=FALSE]
#================================================#

all_genes <- yngpwrveh_v_yngsedveh_genes_16sep %>%
  pull(Symbol)

library(biomaRt)

# Connect to Ensembl archived version if needed
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

# Convert mouse symbols to human symbols
mouse_to_human <- getLDS(attributes = c("mgi_symbol"),
                         filters = "mgi_symbol",
                         values = all_genes,
                         mart = mouse,
                         attributesL = c("hgnc_symbol"),
                         martL = human,
                         uniqueRows = TRUE)

# Extract only human gene symbols and make unique
all_genes_human <- unique(toupper(mouse_to_human$HGNC.symbol))

#Pull out all genes from total Human Protein Atlas data
universe_sc_subset <- sc_data[sc_data$Gene.name %in% all_genes_human, ]

# View result
head(all_genes)
length(all_genes)
dim(universe_sc_subset)
#--------------------------#
#Create Heatmap for expression levels in cell types for each gene

# Convert to wide format: rows = genes, columns = cell types
universe_sc_wide <- universe_sc_subset %>%
  dplyr::select(Gene.name, Cell.type, nTPM) %>%
  pivot_wider(names_from = Cell.type, values_from = nTPM)

# Set Gene.name as rownames
universe_sc_matrix <- as.data.frame(universe_sc_wide)
rownames(universe_sc_matrix) <- universe_sc_matrix$Gene.name
universe_sc_matrix$Gene.name <- NULL

# Turn into numeric matrix
universe_sc_matrix <- as.matrix(universe_sc_matrix)

# log transform
universe_sc_matrix_log <- log2(universe_sc_matrix + 1)

# Plot heatmap
library(pheatmap)
pheatmap(universe_sc_matrix_log, scale = "row", show_rownames = TRUE)
#--------------------------------------------------------------#
#------------Determine nTPM Cutoff values for cell type enrichment analysis------#

expression_summary <- universe_sc_subset %>%
  group_by(Cell.type) %>%
  summarise(
    median_nTPM = median(nTPM, na.rm = TRUE),
    Q1_nTPM = quantile(nTPM, 0.25, na.rm = TRUE),
    Q3_nTPM = quantile(nTPM, 0.75, na.rm = TRUE),
    mean_nTPM = mean(nTPM, na.rm = TRUE),
    max_nTPM = max(nTPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_nTPM))

# View top rows
head(expression_summary)

library(ggplot2)

ggplot(universe_sc_subset, aes(x = reorder(Cell.type, nTPM, FUN = median), y = nTPM)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  labs(x = "Cell Type", y = "nTPM", title = "Distribution of Bcl6b Target Gene Expression Across Cell Types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(ylim = c(0, 1000))

print(expression_summary, n=81)
#----------------------------------------------#

muscle_markers <- c("MYH1", "MYH2", "MYH7", "TNNT3", "TNNC2", "ACTA1", "CKM", "RYR1", "DMD")

muscle_expression <- universe_sc_subset %>%
  filter(Gene.name %in% muscle_markers)

print(muscle_expression)
#========================================================#
#========================================================#

#--------GTEx Analysis----------------------------------#

#--- Load packages
library(vroom)
if (!requireNamespace("pbapply", quietly = TRUE)) install.packages("pbapply")
library(pbapply)

#--- Set working directory and read GTEx TPM (skip header)
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")
gtex_tpm <- vroom("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip = 2)

#--- Rename first 2 columns
colnames(gtex_tpm)[1:2] <- c("EnsemblID", "GeneSymbol")

#--- Filter: valid symbols, not duplicated
gtex_tpm <- gtex_tpm[!duplicated(gtex_tpm$GeneSymbol) & gtex_tpm$GeneSymbol != "", ]

#--- Save BCL6B vector before filtering rest
bcl6b_expr <- as.numeric(gtex_tpm[gtex_tpm$GeneSymbol == "BCL6B", -(1:2)])

#--- Pull matrix of expression values
expr_matrix <- as.matrix(gtex_tpm[, -(1:2)])
rownames(expr_matrix) <- gtex_tpm$GeneSymbol
gene_names <- gtex_tpm$GeneSymbol

#--- Chunked expression filter (TPM > 1 in ≥10% of samples)
chunk_size <- 1000
n_genes <- nrow(expr_matrix)
n_chunks <- ceiling(n_genes / chunk_size)
expressed_filter <- logical(n_genes)

for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, n_genes)
  chunk <- expr_matrix[start:end, , drop = FALSE]
  expressed_filter[start:end] <- pbapply(chunk, 1, function(x) sum(x > 1) >= 0.10 * length(x))
}

#--- Apply filter
expr_matrix <- expr_matrix[expressed_filter, ]
gene_names <- gene_names[expressed_filter]

# Check remaining genes
length(gene_names)



#=======================================================#
#=======================================================#
#--------Original GTEx Analysis--------------------------#
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")
library(vroom)
gtex_tpm <- vroom::vroom("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip = 2)

# Rename columns
colnames(gtex_tpm)[1:2] <- c("EnsemblID", "GeneSymbol")

# Remove duplicates or missing symbols
gtex_tpm <- gtex_tpm[!duplicated(gtex_tpm$GeneSymbol) & !is.na(gtex_tpm$GeneSymbol), ]

bcl6b_expr <- as.numeric(gtex_tpm[gtex_tpm$GeneSymbol == "BCL6B", -(1:2)])

length(bcl6b_expr)

# Create base filter for annotation quality
valid_genes <- gtex_tpm$GeneSymbol != "" & !duplicated(gtex_tpm$GeneSymbol)

# Create expression matrix (exclude first 2 columns)
expr_matrix <- gtex_tpm[valid_genes, -(1:2)]

# Optional: attach gene names for later
gene_names <- gtex_tpm$GeneSymbol[valid_genes]

# Filter: Keep genes expressed at TPM > 1 in at least 10% of samples
# Use pbapply to track progress and avoid freezing
install.packages("pbapply")
library(pbapply)

expressed_filter <- pbapply::pbapply(expr_matrix, 1, function(x) sum(x > 1) >= 0.10 * length(x))

# Apply both filters
expr_matrix <- expr_matrix[expressed_filter, ]
gene_names <- gene_names[expressed_filter]

library(pbapply)
test_mat <- matrix(runif(1e5), nrow = 1000)
test_pb <- pbapply(test_mat, 1, function(x) mean(x) > 0.5)

system.time({
  test_run <- pbapply::pbapply(expr_matrix[1:1000, ], 1, function(x) sum(x > 1) >= 0.10 * length(x))
})
