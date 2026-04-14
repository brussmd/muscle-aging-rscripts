#Create new Gene Aging Signature
#Core enrichment genes from old vs young GO Pathways
#Keep both increased and decreased with age genes.
#Keep genes with old vs young (FDR <0.05)
#Keep genes with young pwr vs yng sed (FDR >0.05)
#Run gene set enrichment
#Keep aging responsive genes against Glass paper
#re-run gene set enrichment
#------------------------------#
#Code using variables and functions established in GO_comparison_vs_oldsedveh.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


#-----Create old_core_no_pwr_deg ENTREZID gene list----------------#
oldsed_vs_yngsed_GO_core_enrich_genes <- read_csv("oldsed_vs_yngsed_GO_core_enrich_genes.csv")
head(oldsed_vs_yngsed_GO_core_enrich_genes)

oldsed_core_genes <- oldsed_vs_yngsed_GO_core_enrich_genes %>%
  pull(EntrezID) %>%
  unique()
length(oldsed_core_genes) #should be 833 genes

oldrapa_gene_merged_df <- read_csv("oldrapa_gene_merged_df.csv")
head(oldrapa_gene_merged_df)

oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  filter(oldsedveh_FDR <0.05) #119 genes split ~evenly logFC >0 <0
  
yngpwrveh_v_yngsedveh_genes_29Nov <- read_clean_xlsx_simple(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_YNG_PWR_VEH-YNG_SED_VEH.xlsx")

head(yngpwrveh_v_yngsedveh_genes_29Nov)

yngpwrveh_v_yngsedveh_genes_29Nov %>%
  filter(FDR <0.05)

old_core_deg <- oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  filter(oldsedveh_FDR <0.05) %>%
  pull(ENTREZID)

old_core_no_pwr_deg <- yngpwrveh_v_yngsedveh_genes_29Nov %>%
  filter(ENTREZID %in% old_core_deg) %>%
  filter(FDR >0.05) %>%
  pull(ENTREZID)

length(old_core_no_pwr_deg)

oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% old_core_no_pwr_deg) %>%
  filter(oldsedveh_FDR <0.05) %>%
  filter(oldsedveh_logFC <0)

head(old_core_deg)
length(old_core_deg)
#----------------------------------------------#

#==============================================#
#---Targeted GSEA for oldpwrveh vs oldsedveh---#
#==============================================#

head(oldpwrveh_v_oldsedveh_genes_22Nov)

aging_df <- oldsedveh_v_yngsedveh_genes_22Nov

head(aging_df)

#--------------------------------------------------------------#
# 1. AGING DIRECTION: compute using Ensembl_noDec (clean join)
#--------------------------------------------------------------#

aging_direction_df <- aging_df %>%
  mutate(aging_direction = sign(logFC)) %>%     # OLD vs YNG direction
  dplyr::select(Ensembl_noDec, aging_direction)

head(aging_direction_df)
#--------------------------------------------------------------#
# 2. JOIN INTERVENTION DF TO AGING DIRECTION BY ENSEMBL
#--------------------------------------------------------------#

int_df <- oldpwrveh_v_oldsedveh_genes_22Nov %>%
  left_join(aging_direction_df, by = "Ensembl_noDec")

#--------------------------------------------------------------#
# 3. COMPUTE AGING-AWARE RANK METRIC
#--------------------------------------------------------------#

int_df <- int_df %>%
  mutate(rank_metric = aging_direction * logFC *
           ((-log10(PValue + 1e-300))^0.25))

#--------------------------------------------------------------#
# 4. FIRST remove exact duplicate rows (identical rows)
#--------------------------------------------------------------#

int_df_unique_rows <- int_df %>%
  distinct()

#--------------------------------------------------------------#
# 5. THEN deduplicate ENTREZID by *lowest intervention PValue*
#--------------------------------------------------------------#

int_df_dedup <- int_df_unique_rows %>%
  arrange(ENTREZID, PValue) %>%                # lowest PValue first
  distinct(ENTREZID, .keep_all = TRUE)         # keep best ENTREZID row

#--------------------------------------------------------------#
# 6. BUILD RANKED VECTOR FOR GSEA (ENTREZID -> rank_metric)
#--------------------------------------------------------------#

ranked_vec <- int_df_dedup %>%
  dplyr::select(ENTREZID, rank_metric) %>%
  drop_na() %>%
  arrange(desc(rank_metric)) %>%
  deframe()

# Confirm no duplicates
sum(duplicated(names(ranked_vec)))   # should be 0

#--------------------------------------------------------------#
# 7. GENESET PREPARATION
#--------------------------------------------------------------#

AgingCore_geneSet <- list(AgingCore = unique(old_core_deg))

#--------------------------------------------------------------#
# 8. RUN GSEA
#--------------------------------------------------------------#

gsea_res <- GSEA(
  geneList     = ranked_vec,
  TERM2GENE    = data.frame(term = "AgingCore", gene = unique(old_core_deg)),
  pvalueCutoff = 1,
  verbose      = FALSE
)

# View clean NES output
gsea_res@result %>% dplyr::select(ID, NES, pvalue, p.adjust)

head(int_df_dedup)
#--------------------------------------------------------------#

head(oldpwrfrap_v_oldsedveh_genes_22Nov)
#==============================================================#
#//////////////////////////////////////////////////////////////#
#==============================================================#

#=========================================================#
#///Targeted GSEA Function////////////////////////#
#-------------------------------------------------#

run_targeted_aging_gsea <- function(intervention_df,
                                    aging_df,
                                    core_genes,
                                    intervention_name = "Intervention") {
  
  #-------------------------------#
  # 1. Aging direction by Ensembl #
  #-------------------------------#
  aging_direction_df <- aging_df %>%
    mutate(aging_direction = sign(logFC)) %>%
    dplyr::select(Ensembl_noDec, aging_direction)
  
  #-------------------------------#
  # 2. Join intervention by Ensembl
  #-------------------------------#
  int_df <- intervention_df %>%
    left_join(aging_direction_df, by = "Ensembl_noDec")
  
  #-------------------------------#
  # 3. Compute rank metric
  #-------------------------------#
  int_df <- int_df %>%
    mutate(rank_metric = aging_direction * logFC *
             ((-log10(PValue + 1e-300))^0.25))
  
  #-------------------------------#
  # 4. Remove exact duplicated rows
  #-------------------------------#
  int_df_unique_rows <- int_df %>%
    distinct()
  
  #-------------------------------#
  # 5. Deduplicate ENTREZID
  #    (lowest PValue = strongest signal)
  #-------------------------------#
  int_df_dedup <- int_df_unique_rows %>%
    arrange(ENTREZID, PValue) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  #-------------------------------#
  # 6. Build ranked vector
  #-------------------------------#
  ranked_vec <- int_df_dedup %>%
    dplyr::select(ENTREZID, rank_metric) %>%
    drop_na() %>%
    arrange(desc(rank_metric)) %>%
    deframe()
  
  # Safety check
  if (sum(duplicated(names(ranked_vec))) > 0) {
    stop("Duplicate ENTREZIDs remain after deduplication.")
  }
  
  #-------------------------------#
  # 7. Build geneset
  #-------------------------------#
  core_genes <- unique(core_genes)
  
  #-------------------------------#
  # 8. Run GSEA
  #-------------------------------#
  gsea_res <- GSEA(
    geneList     = ranked_vec,
    TERM2GENE    = data.frame(term = intervention_name, gene = core_genes),
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  
  # Return clean result table
  return(gsea_res@result %>% dplyr::select(ID, NES, pvalue, p.adjust))
}

#OldPwrVeh GSEA
res_oldpwrveh <- run_targeted_aging_gsea(
  intervention_df = oldpwrveh_v_oldsedveh_genes_22Nov,
  aging_df = oldsedveh_v_yngsedveh_genes_22Nov,
  core_genes = old_core_no_pwr_deg,
  intervention_name = "PWR_VEH"
)
res_oldpwrveh

#OldSedFrap GSEA
res_oldsedfrap <- run_targeted_aging_gsea(
  intervention_df = oldsedfrap_v_oldsedveh_genes_22Nov,
  aging_df = oldsedveh_v_yngsedveh_genes_22Nov,
  core_genes = old_core_no_pwr_deg,
  intervention_name = "SED_FRAP"
)
res_oldsedfrap

#OldSedIrap GSEA
res_oldsedirap <- run_targeted_aging_gsea(
  intervention_df = oldsedirap_v_oldsedveh_genes_22Nov,
  aging_df = oldsedveh_v_yngsedveh_genes_22Nov,
  core_genes = old_core_no_pwr_deg,
  intervention_name = "SED_IRAP"
)
res_oldsedirap

#OldPwrIrap GSEA
res_oldpwrirap <- run_targeted_aging_gsea(
  intervention_df = oldpwrirap_v_oldsedveh_genes_22Nov,
  aging_df = oldsedveh_v_yngsedveh_genes_22Nov,
  core_genes = old_core_no_pwr_deg,
  intervention_name = "PWR_IRAP"
)
res_oldpwrirap

#OldPwrFrap GSEA
res_oldpwrfrap <- run_targeted_aging_gsea(
  intervention_df = oldpwrfrap_v_oldsedveh_genes_22Nov,
  aging_df = oldsedveh_v_yngsedveh_genes_22Nov,
  core_genes = old_core_no_pwr_deg,
  intervention_name = "PWR_FRAP"
)
res_oldpwrfrap

old_core_no_pwr_deg <- yngpwrveh_v_yngsedveh_genes_29Nov %>%
  filter(ENTREZID %in% old_core_deg) %>%
  filter(FDR >0.05)
head(old_core_no_pwr_deg)
dim(old_core_no_pwr_deg)

head(oldsedveh_v_yngsedveh_genes_22Nov)

down_core_no_pwr <- oldsedveh_v_yngsedveh_genes_22Nov %>%
  filter(ENTREZID %in% old_core_no_pwr_deg) %>%
  filter(logFC < 0) %>%
  pull(ENTREZID) #should be 39 of the 89 genes

down_core_no_pwr

up_core_no_pwr <- oldsedveh_v_yngsedveh_genes_22Nov %>%
  filter(ENTREZID %in% old_core_no_pwr_deg) %>%
  filter(logFC > 0) %>%
  pull(ENTREZID) #should be 50 of the 89 genes

up_core_no_pwr

oldpwrveh_v_oldsedveh_genes_22Nov %>%
  filter(ENTREZID %in% down_core_no_pwr) %>%
  arrange(logFC)

oldpwrveh_v_oldsedveh_genes_22Nov %>%
  filter(ENTREZID %in% up_core_no_pwr) %>%
  arrange(logFC)
#-----------------------------------------------#

head(oldrapa_gene_merged_df)

oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% up_core_no_pwr) %>%
  arrange(oldpwrveh_logFC) 


oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% old_core_no_pwr_deg)

#===========================================================#
#////Vertical bar plot for all old_core_no_pwr_deg genes////#
#===========================================================#

df_plot <- oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% old_core_no_pwr_deg) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%   # safety dedupe
  mutate(
    direction = if_else(oldsedveh_logFC >= 0, "Up", "Down")
  ) %>%
  arrange(oldsedveh_logFC)

# Keep symbol order for nicer plotting
df_plot$Symbol <- factor(df_plot$Symbol, levels = df_plot$Symbol)

core_enrich_gene_plot <- ggplot(df_plot, aes(x = Symbol, y = oldsedveh_logFC, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "black", "Down" = "red")) +
  coord_flip() +
  labs(
    title = "oldsedveh_logFC for Core Enrichment Genes (n = 89)",
    x = "",
    y = "log2 Fold Change (oldsedveh)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title  = element_text(face = "bold", size = 16),
    legend.title = element_blank()
  )

ggsave("core_enrich_gene_plot.pdf", plot = core_enrich_gene_plot, width = 8, height = 8, units = "in", dpi = 600)
#========================================================#

#================================================================#
#////Individual boxplots for upregulated core enrichment genes///#
#================================================================#
oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% up_core_no_pwr)
head(oldrapa_gene_merged_df)

#------------------------------------------
# 1. Filter to genes of interest (deduped)
#------------------------------------------
df_subset <- oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% up_core_no_pwr) %>%
  distinct(ENTREZID, .keep_all = TRUE)

#------------------------------------------
# 2. Pivot both logFC *and* PValue columns
#------------------------------------------
df_long <- df_subset %>%
  dplyr::select(
    ENTREZID, Symbol,
    
    oldpwrveh_logFC,   oldpwrveh_PValue,
    oldsedirap_logFC,  oldsedirap_PValue,
    oldsedfrap_logFC,  oldsedfrap_PValue,
    oldpwrirap_logFC,  oldpwrirap_PValue,
    oldpwrfrap_logFC,  oldpwrfrap_PValue
  ) %>%
  pivot_longer(
    cols = matches("logFC|PValue"),
    names_to = c("intervention", ".value"),
    names_pattern = "(.*)_(logFC|PValue)"
  )

#------------------------------------------
# 3. Recode intervention names & compute stats
#------------------------------------------
df_long$intervention <- factor(
  df_long$intervention,
  levels = c("oldpwrveh", "oldsedirap", "oldsedfrap", "oldpwrirap", "oldpwrfrap"),
  labels = c("OldPwrVeh", "OldSediRapa", "OldSedfRapa", "OldPwriRapa", "OldPwrfRapa")
)

df_long <- df_long %>%
  mutate(
    color_flag = ifelse(logFC < 0, "neg", "pos"),
    size_metric = -log10(PValue + 1e-300)  # avoid log(0)
  )

#------------------------------------------
# 4. Boxplot with jitter sized by -log10(PValue)
#------------------------------------------
up_core_boxplot <- ggplot(df_long, aes(x = intervention, y = logFC)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black") +
  geom_jitter(
    aes(color = color_flag, size = size_metric),
    width = 0.15, alpha = 0.75
  ) +
  scale_color_manual(values = c("pos" = "black", "neg" = "blue")) +
  scale_size_continuous(range = c(0.5, 5)) +   # adjust dot size range
  labs(
    title = "logFC Distribution for Core Genes Across Interventions",
    x = "",
    y = "log2 Fold Change",
    size = "-log10(P-value)"
  ) +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    legend.position = "right"
  )
ggsave("up_core_boxplot.pdf", plot = up_core_boxplot, width = 8, height = 8, units = "in", dpi = 600)
#===================================================================#
#===================================================================#


#================================================================#
#////Individual boxplots for downregulated core enrichment genes///#
#================================================================#
oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% down_core_no_pwr)
head(oldrapa_gene_merged_df)

#------------------------------------------
# 1. Filter to genes of interest (deduped)
#------------------------------------------
df_down_subset <- oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% down_core_no_pwr) %>%
  distinct(ENTREZID, .keep_all = TRUE)

#------------------------------------------
# 2. Pivot both logFC *and* PValue columns
#------------------------------------------
df_down_long <- df_down_subset %>%
  dplyr::select(
    ENTREZID, Symbol,
    
    oldpwrveh_logFC,   oldpwrveh_PValue,
    oldsedirap_logFC,  oldsedirap_PValue,
    oldsedfrap_logFC,  oldsedfrap_PValue,
    oldpwrirap_logFC,  oldpwrirap_PValue,
    oldpwrfrap_logFC,  oldpwrfrap_PValue
  ) %>%
  pivot_longer(
    cols = matches("logFC|PValue"),
    names_to = c("intervention", ".value"),
    names_pattern = "(.*)_(logFC|PValue)"
  )

#------------------------------------------
# 3. Recode intervention names & compute stats
#------------------------------------------
df_down_long$intervention <- factor(
  df_down_long$intervention,
  levels = c("oldpwrveh", "oldsedirap", "oldsedfrap", "oldpwrirap", "oldpwrfrap"),
  labels = c("OldPwrVeh", "OldSediRapa", "OldSedfRapa", "OldPwriRapa", "OldPwrfRapa")
)

df_down_long <- df_down_long %>%
  mutate(
    color_flag = ifelse(logFC < 0, "neg", "pos"),
    size_metric = -log10(PValue + 1e-300)  # avoid log(0)
  )

#------------------------------------------
# 4. Boxplot with jitter sized by -log10(PValue)
#------------------------------------------
down_core_boxplot <- ggplot(df_down_long, aes(x = intervention, y = logFC)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black") +
  geom_jitter(
    aes(color = color_flag, size = size_metric),
    width = 0.15, alpha = 0.75
  ) +
  scale_color_manual(values = c("pos" = "blue", "neg" = "black")) +
  scale_size_continuous(range = c(0.5, 5)) +   # adjust dot size range
  labs(
    title = "logFC Distribution for Core Genes Across Interventions",
    x = "",
    y = "log2 Fold Change",
    size = "-log10(P-value)"
  ) +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    legend.position = "right"
  )
ggsave("down_core_boxplot.pdf", plot = down_core_boxplot, width = 8, height = 8, units = "in", dpi = 600)
#===================================================================#
#===================================================================#


#=============================================================#
#/////Horizontal Bar Graph for Intervention NES//////////////#
#=============================================================#

res_oldpwrveh
res_oldsedirap
res_oldsedfrap
res_oldpwrirap
res_oldpwrfrap

library(dplyr)
library(ggplot2)

#-----------------------------
# Build NES dataframe
#-----------------------------
nes_df <- dplyr::bind_rows(
  OldPwrVeh  = res_oldpwrveh,
  OldSedIRap = res_oldsedirap,
  OldSedFRap = res_oldsedfrap,
  OldPwrIRap = res_oldpwrirap,
  OldPwrFRap = res_oldpwrfrap,
  .id = "Intervention"
) %>%
  dplyr::select(Intervention, NES, p.adjust) %>%
  mutate(
    Intervention = factor(
      Intervention,
      levels = rev(c("OldPwrVeh", "OldSedIRap", "OldSedFRap",
                     "OldPwrIRap", "OldPwrFRap"))  # reversed order
    )
  )

#-----------------------------
# Unique color per intervention
#-----------------------------
intervention_colors <- c(
  "OldPwrVeh"  = "#1B9E77",
  "OldSedIRap" = "#e69f00",
  "OldSedFRap" = "#ff7f00",
  "OldPwrIRap" = "#CC79A7",
  "OldPwrFRap" = "#9E1F63"
)

#-----------------------------#
# Plot
#-----------------------------
NES_intervention_hoz_barplot <- ggplot(nes_df, aes(x = NES, y = Intervention, fill = Intervention)) +
  geom_col(width = 0.65, color = "black") +
  scale_fill_manual(values = intervention_colors) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  labs(
    title = "NES Values Across Interventions",
    x = "Normalized Enrichment Score (NES)",
    y = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )
ggsave("NES_intervention_hoz_barplot.pdf", plot = NES_intervention_hoz_barplot, width = 8, height = 8, units = "in", dpi = 600)



#make horizontal bar chart of all 89 old_core_no_pwr_deg
#then make vertical bars with each interventions logFC