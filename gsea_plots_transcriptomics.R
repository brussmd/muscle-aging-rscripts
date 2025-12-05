
#-------Read in xlsx, Clean Dataframe and add ENTREZID col-----------#
read_clean_xlsx <- function(xlsx_file, logFC_col, FDR_col, pvalue_col) {
  # Read in the xlsx file
  df <- read_xlsx(xlsx_file)
  
  # Clean and rename columns
  df <- df %>%
    dplyr::select(Ensembl, Symbol, 
                  logFC = all_of(logFC_col),
                  FDR = all_of(FDR_col),
                  pvalue = all_of(pvalue_col)) %>%
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
#---------------------------------------------------------------#

#----OldPwrVeh vs OldSedVeh Inflammaging GeneSet GSEA Analysis----#
oldpwrveh_v_oldsedveh_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_PWR_VEH-OLD_SED_VEH_logFC",
  FDR_col = "OLD_PWR_VEH-OLD_SED_VEH_FDR",
  pvalue_col = "OLD_PWR_VEH-OLD_SED_VEH_PValue")

oldpwrveh_v_oldsedveh_genes %>%
  filter(FDR <0.05)

oldpwrveh_v_oldsedveh_genes %>%
  filter(ENTREZID %in% validated_inflamm_geneset.entrez)

#----------------------------------------------------------------#
# 1) Build the ranking: -log10(p) * logFC
oldpwr_rank_df <- oldpwrveh_v_oldsedveh_genes %>%
  filter(!is.na(logFC), !is.na(pvalue)) %>%
  mutate(
    p_safe = pmax(pvalue, 1e-300),                 # avoid -Inf
    stat   = -log10(p_safe) * logFC
  ) %>%
  group_by(ENTREZID) %>%                                # in case of duplicate names
  summarise(stat = mean(stat), .groups = "drop") %>% # or max(stat), if preferred
  arrange(desc(stat))

# 2) Named numeric vector for fgsea
OP_OS_ranks <- setNames(oldpwr_rank_df$stat, oldpwr_rank_df$ENTREZID)

# 3) Your pathway(s)
gene_sets <- list(MusAge = validated_inflamm_geneset.entrez)


# Run fgsea
OP_fgsea_res <- fgsea(pathways = gene_sets,
                      stats = OP_OS_ranks,
                      nperm = 1000)

print(OP_fgsea_res)

# Plot enrichment

OP_gsea_plot <-plotEnrichment(gene_sets$MusAge, OP_OS_ranks) +
  ggtitle("Effect of Exercise on MusAge Gene Set")

# Save as PDF (adjust size as needed)
pdf("OP_gsea_plot.pdf", width = 4.25, height = 5.5)
print(OP_gsea_plot)
dev.off()

#-----Find flip point-----------------------------------------#

# ---- 1. Create ranked statistic data frame ----
df_bottom_op <- tibble(
  Rank = 1:length(OP_OS_ranks),
  ENTREZID = names(OP_OS_ranks),
  Stat = OP_OS_ranks
)

# Where does the sign switch?
pos_idx <- which(df_bottom_op$Stat > 0)
neg_idx <- which(df_bottom_op$Stat < 0)

last_pos_rank  <- if (length(pos_idx)) max(pos_idx) else NA_integer_
first_neg_rank <- if (length(neg_idx)) min(neg_idx) else NA_integer_

list(
  last_positive_rank  = last_pos_rank,
  last_positive_entrezid = if (!is.na(last_pos_rank)) df_bottom_op$ENTREZID[last_pos_rank] else NA,
  first_negative_rank = first_neg_rank,
  first_negative_entrezid= if (!is.na(first_neg_rank)) df_bottom_op$ENTREZID[first_neg_rank] else NA,
  n_positive = length(pos_idx),
  n_negative = length(neg_idx),
  n_zero     = sum(df_bottom_op$Stat == 0)
)

# ---- 2. Bottom panel (stat values with color scale) ----
# Calculate flip_x
last_pos_rank <- max(which(df_bottom_op$Stat > 0))
flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA

# Define your cap, e.g. 95th percentile or absolute max
stat_cap <- quantile(abs(df_bottom_op$Stat), 0.99)

op_vs_os_p_bottom <- ggplot(df_bottom_op, aes(x = Rank, y = Stat, fill = Stat)) +
  geom_col(width = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       limits = c(-stat_cap, stat_cap),
                       midpoint = 0,
                       name = "Ranked Statistic\n(logFC)") +
  { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
  ylim(-stat_cap, stat_cap) +  # <--- This clips extreme values for plotting
  labs(x = "ENTREZID Rank", y = "Ranked Statistic") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(op_vs_os_p_bottom)

# Save as PDF (adjust size as needed)
pdf("op_vs_os_p_bottom.pdf", width = 4.25, height = 5.5)
print(op_vs_os_p_bottom)
dev.off()
#-----------------------------------------------------------#

#===========================================================#
#------Old Sed vs Yng Sed-----------------------------------#
#===========================================================#

#----OldSedVeh vs YngSedVeh Inflammaging GeneSet GSEA Analysis----#
oldsedveh_v_yngsedveh_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col = "OLD_SED_VEH-YNG_SED_VEH_FDR",
  pvalue_col = "OLD_SED_VEH-YNG_SED_VEH_PValue")

oldsedveh_v_yngsedveh_genes %>%
  filter(FDR <0.05)

oldsedveh_v_yngsedveh_genes %>%
  filter(ENTREZID %in% validated_inflamm_geneset.entrez)

#----------------------------------------------------------------#

#----------------------------------------------------------------#
# 1) Build the ranking: -log10(p) * logFC
oldsed_rank_df <- oldsedveh_v_yngsedveh_genes %>%
  filter(!is.na(logFC), !is.na(pvalue)) %>%
  mutate(
    p_safe = pmax(pvalue, 1e-300),                 # avoid -Inf
    stat   = -log10(p_safe) *logFC
  ) %>%
  group_by(ENTREZID) %>%                                # in case of duplicate names
  summarise(stat = mean(stat), .groups = "drop") %>% # or max(stat), if preferred
  arrange(desc(stat))

# 2) Named numeric vector for fgsea
OS_YS_ranks <- setNames(oldsed_rank_df$stat, oldsed_rank_df$ENTREZID)

# 3) Your pathway(s)
gene_sets <- list(MusAge = validated_inflamm_geneset.entrez)


# Run fgsea
OS_fgsea_res <- fgsea(pathways = gene_sets,
                      stats = OS_YS_ranks)

print(OS_fgsea_res)

# Plot enrichment

OS_gsea_plot <-plotEnrichment(gene_sets$MusAge, OS_YS_ranks) +
  ggtitle("Effect of Age on MusAge Gene Set")

# Save as PDF (adjust size as needed)
pdf("OS_gsea_plot.pdf", width = 4.25, height = 5.5)
print(OS_gsea_plot)
dev.off()
#----------------------------------------------#

#-----Find flip point-----------------------------------------#

# ---- 1. Create ranked statistic data frame ----
df_bottom_os <- tibble(
  Rank = 1:length(OS_YS_ranks),
  ENTREZID = names(OS_YS_ranks),
  Stat = OS_YS_ranks
)

# Where does the sign switch?
pos_idx <- which(df_bottom_os$Stat > 0)
neg_idx <- which(df_bottom_os$Stat < 0)

last_pos_rank  <- if (length(pos_idx)) max(pos_idx) else NA_integer_
first_neg_rank <- if (length(neg_idx)) min(neg_idx) else NA_integer_

list(
  last_positive_rank  = last_pos_rank,
  last_positive_entrezid = if (!is.na(last_pos_rank)) df_bottom_os$ENTREZID[last_pos_rank] else NA,
  first_negative_rank = first_neg_rank,
  first_negative_entrezid= if (!is.na(first_neg_rank)) df_bottom_os$ENTREZID[first_neg_rank] else NA,
  n_positive = length(pos_idx),
  n_negative = length(neg_idx),
  n_zero     = sum(df_bottom_os$Stat == 0)
)

# ---- 2. Bottom panel (stat values with color scale) ----
# Calculate flip_x
last_pos_rank <- max(which(df_bottom_os$Stat > 0))
flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA

# Define your cap, e.g. 95th percentile or absolute max
stat_cap <- quantile(abs(df_bottom_os$Stat), 0.99)

os_vs_ys_p_bottom <- ggplot(df_bottom_os, aes(x = Rank, y = Stat, fill = Stat)) +
  geom_col(width = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       limits = c(-stat_cap, stat_cap),
                       midpoint = 0,
                       name = "Ranked Statistic\n(logFC)") +
  { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
  ylim(-stat_cap, stat_cap) +  # <--- This clips extreme values for plotting
  labs(x = "ENTREZID Rank", y = "Ranked Statistic") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(os_vs_ys_p_bottom)

# Save as PDF (adjust size as needed)
pdf("os_vs_ys_p_bottom.pdf", width = 4.25, height = 5.5)
print(os_vs_ys_p_bottom)
dev.off()
#-----------------------------------------------------------#