########################################
#//Load Libraries/////#
#######################################
suppressPackageStartupMessages({
  library(tidyverse)      # loads dplyr, readr, ggplot2, tibble, tidyr, stringr
  library(readxl)
  library(ggVennDiagram)
  library(eulerr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(DOSE)
  library(pheatmap)
  library(ComplexHeatmap)
  library(fgsea)
  library(ggrepel)
  library(edgeR)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(httr)
  library(jsonlite)
})

#####################################################
#####################################################
#-------Initial Analysis-----------------------------------#
#Identify old Pwr DEG ----#
#####################################################

setwd("/Users/mdbruss/Documents/RStudioProjects_2/Rapa_PwR")

#Helper function to read in gene files.
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
##################################################################

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
###################################################################
###################################################################

###################################################################
#####///Analysis of Old SED FRAP vs Old SED VEH////////###########
#-----------------------------------------------------------------#
# 1) Load edgeR *_GENE_* xlsx and build ranked vector
oldpwrveh_v_oldsedveh_genes_22feb <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_PWR_VEH-OLD_SED_VEH_logFC",
  FDR_col   = "OLD_PWR_VEH-OLD_SED_VEH_FDR"
)

#remove duplicate ENTREZID and keep most significant
oldpwrveh_v_oldsedveh_genes_22feb <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()

head(oldpwrveh_v_oldsedveh_genes_22feb)

oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(FDR < 0.05) %>%
  arrange(desc(logFC)) %>%
  print(n = 100, width = 120)


oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(FDR < 0.05) %>%
  arrange(logFC) %>%
  print(n = 100, width = 120)

oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(Symbol == "Nmrk2")
############################################
####################################################

# Create significance and direction column (using sig_status, not regulation)
oldpwrveh_v_oldsedveh_genes_22feb <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  mutate(sig_status = case_when(
    FDR < 0.05 & logFC > 0 ~ "Up",
    FDR < 0.05 & logFC < 0 ~ "Down",
    TRUE ~ "NS"
  ))

# Label top genes (5 up, 5 down)
oldpwr_top_up <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(sig_status == "Up") %>%
  arrange(desc(logFC)) %>%
  head(50)

print(oldpwr_top_up, n = 50)

oldpwr_top_down <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(sig_status == "Down") %>%
  arrange(logFC) %>%
  head(50)

print(oldpwr_top_down, n = 50)


################################################
# Curated upregulated labels
up_labels <- c("Ddit4", "Depp1", "Ankrd2", "Csrp3", "Hif3a",
               "Cidea", "Arrdc2", "Myh2", "Mb", "Cebpd")

up_oldpwr_genes <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(Symbol %in% up_labels)

# Curated downregulated labels
down_labels <- c("Fgf21", "Nmrk2", "Apol9a", "Egr1", "Cd28",
                 "Nr4a3", "Ifit3b", "Mybph", "Ppp1r27", "Prkag3")

dwn_oldpwr_genes <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(Symbol %in% down_labels)
#=================================================================#

###########################################################
#//////Volcano Plot for OLD PWR///////////////#
##############################################################
# First, create the top_genes dataframe from leading edge analysis
top_PWR_genes <- bind_rows(up_oldpwr_genes, dwn_oldpwr_genes)

top_PWR_genes

# Create the volcano plot with highlighted top genes
oldpwr_volcano_plot <- ggplot(oldpwrveh_v_oldsedveh_genes_22feb, aes(x = logFC, y = -log10(FDR))) +
  
  # background points
  geom_point(aes(color = sig_status),
             alpha = 0.4, size = 1.2) +
  
  # cutoff lines
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  
  # Highlight upregulated top genes with darker red and black outline
  geom_point(
    data = top_PWR_genes %>% filter(logFC > 0),
    shape = 21,
    size = 3.5,
    stroke = 1.2,
    fill = "#8B0000",      # Dark red
    color = "black"
  ) +
  
  # Highlight downregulated top genes with darker blue and black outline
  geom_point(
    data = top_PWR_genes %>% filter(logFC < 0),
    shape = 21,
    size = 3.5,
    stroke = 1.2,
    fill = "#00008B",      # Dark blue
    color = "black"
  ) +
  
  # label selected genes
  geom_text_repel(
    data = top_PWR_genes,
    aes(label = Symbol),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    nudge_y = 0.5,
    fontface = "bold"
  ) +
  
  scale_color_manual(
    values = c(
      "Up"   = "#D55E00",
      "Down" = "#0072B2",
      "NS"   = "grey75"
    )
  ) +
  
  coord_cartesian(xlim = c(-3, 3), ylim = c(-0.5, 12)) +
  
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Status",
    title = "Volcano plot: Old Pwr Veh vs Old Sed Veh",
    subtitle = "Highlighted: Top leading edge genes from pathway enrichment"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title       = element_text(face = "bold")
  )

print(oldpwr_volcano_plot)


ggsave(
  "FigX_oldpwr_volcano_plot.pdf",
  oldpwr_volcano_plot,
  width = 8,
  height = 8
)

ggsave(
  "FigX_oldpwr_volcano_plot.png",
  oldpwr_volcano_plot,
  width = 8,
  height = 8,
  dpi = 300
)
###############################################################

#########################################################################################
create_vec_from_df(oldpwrveh_v_oldsedveh_genes_22feb, "oldpwrveh_v_oldsedveh_genes_22feb")

# 2) GSEA (GO Biological Process)
gseGO_oldpwrveh_vs_oldsedveh.OUTPUT <- gseGO(
  geneList     = oldpwrveh_v_oldsedveh_genes_22feb.vec, 
  ont          = "BP",
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 300,
  pvalueCutoff = 0.5,      # capture everything, filter later
  eps          = 1e-30,      # better estimation for very small p-values
  verbose      = FALSE
)

head(gseGO_oldpwrveh_vs_oldsedveh.OUTPUT@result)


# Simplify, then compute Count/geneRatio again (simplify changes the table)
oldpwrveh_vs_oldsedveh <- clusterProfiler::simplify(gseGO_oldpwrveh_vs_oldsedveh.OUTPUT,
                                                     cutoff = 0.5, by = "p.adjust", select_fun = min)

oldpwrveh_vs_oldsedveh@result

oldpwrveh_vs_oldsedveh_simpl_22Feb.df <- as.data.frame(oldpwrveh_vs_oldsedveh@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

dim(oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>%
      filter(p.adjust <0.05))

oldpwrveh_vs_oldsedveh_goplot <- oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  scale_size_continuous(range = c(1.5, 9)) +  # Double the size (default is ~1-6)
  theme_bw()


# Display it
print(oldpwrveh_vs_oldsedveh_goplot)

oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>% 
  dplyr::select(ID, Description, NES, Count)

oldpwrveh_vs_oldsedveh_simpl_22Feb.df
###########################################################################

################################################################################
# GO Enrichment Bubble Plot — Exercise (OLD_PWR_VEH vs OLD_SED_VEH)
# Adapted from rapamycin analysis with exercise-relevant pathway categories
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(stringr)
  library(scales)
})

# =============================================================================
# 1. ADD EXERCISE-RELEVANT PATHWAY CATEGORIES
# =============================================================================

oldpwrveh_vs_oldsedveh_simpl_22Feb.df <- oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>%
  mutate(
    Category = case_when(
      # Mitochondria & Oxidative Metabolism
      grepl("electron transport|ATP synthesis|proton.+transport|fatty acid oxidation|ATP metabolic|NADH oxidation|fatty acid derivative|amino acid catabolic|L-amino acid|detoxification|hydrogen peroxide|cold-induced thermogenesis",
            Description, ignore.case = TRUE) ~ "Mitochondria & Metabolism",
      
      # Vascular & Angiogenesis
      grepl("endotheli|vasculogenesis|blood circulation|nitric oxide|angiogen|blood coagulation|coagulation|hemostasis|lipoprotein.+clearance",
            Description, ignore.case = TRUE) ~ "Vascular & Angiogenesis",
      
      # Immune, Inflammation & Interferon
      grepl("immune|interferon|antiviral|viral life cycle|biotic stimulus|T cell|NK T cell|inflammasome|killing.+organism|disruption.+organism|leukocyte|inflammatory|defense response|cGAS|mast cell",
            Description, ignore.case = TRUE) ~ "Immune & Inflammation",
      
      # Translation & Ribosome
      grepl("translat|ribosom|peptidyl-lysine methylation",
            Description, ignore.case = TRUE) ~ "Translation & Ribosome",
      
      # Neuromuscular & Excitation-Contraction
      grepl("membrane depolarization|synaptic transmission|postsynapse|calcium ion concentration|electrical coupling|cyclase activity|regulation of pH",
            Description, ignore.case = TRUE) ~ "Neuromuscular & E-C Coupling",
      
      # Structural & Contractile
      grepl("myofibril|actomyosin|actin filament|mechanical stimulus|muscle hypertrophy|cytoskeleton|lamellipodium|collagen biosynthetic",
            Description, ignore.case = TRUE) ~ "Structural & Contractile",
      
      # Lipid & Nutrient Handling
      grepl("triglyceride|lipid droplet|nutrient storage|nutrient$|steroid hormone|response to nutrient",
            Description, ignore.case = TRUE) ~ "Lipid & Nutrient Handling",
      
      # Everything else
      TRUE ~ "Other"
    )
  )

# Check category assignments
cat("=== Category counts ===\n")
oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>% count(Category) %>% print()

cat("\n=== Category by direction ===\n")
oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>%
  filter(p.adjust < 0.05) %>%
  mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
  count(Category, direction) %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  print()

# =============================================================================
# 2. COMPUTE geneRatio IF NOT ALREADY PRESENT
# =============================================================================
if (!"geneRatio" %in% colnames(oldpwrveh_vs_oldsedveh_simpl_22Feb.df)) {
  oldpwrveh_vs_oldsedveh_simpl_22Feb.df <- oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>%
    mutate(geneRatio = Count / setSize)
}

# =============================================================================
# 3. FILTER TO SIGNIFICANT PATHWAYS AND PREP LABELS
# =============================================================================
df_sig <- oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(
    label = stringr::str_replace_all(Description, "\\s+", " "),
    label = stringr::str_trunc(label, 50)
  )

# =============================================================================
# 4. SPLIT INTO LEFT (negative NES) AND RIGHT (positive NES) PANELS
# =============================================================================
cut <- 1.5
x_span <- 0.75
left_lim  <- c(-(cut + x_span), -cut)
right_lim <- c(cut,  cut + x_span)

df_left  <- df_sig %>% dplyr::filter(NES <= -cut)
df_right <- df_sig %>% dplyr::filter(NES >=  cut)

# =============================================================================
# 5. SELECT LABELS — key pathways from each major category per side
# =============================================================================

# LEFT PANEL: immune/inflammation + translation are the main downregulated themes
lab_left_immune <- df_left %>%
  filter(Category == "Immune & Inflammation") %>%
  slice_max(order_by = abs(NES), n = 4, with_ties = FALSE)

lab_left_translation <- df_left %>%
  filter(Category == "Translation & Ribosome") %>%
  slice_max(order_by = abs(NES), n = 2, with_ties = FALSE)

lab_left <- bind_rows(lab_left_immune, lab_left_translation)

# RIGHT PANEL: mitochondria, vascular, structural are the main upregulated themes
lab_right_mito <- df_right %>%
  filter(Category == "Mitochondria & Metabolism") %>%
  slice_max(order_by = NES, n = 3, with_ties = FALSE)

lab_right_vasc <- df_right %>%
  filter(Category == "Vascular & Angiogenesis") %>%
  slice_max(order_by = NES, n = 3, with_ties = FALSE)

lab_right_struct <- df_right %>%
  filter(Category == "Structural & Contractile") %>%
  slice_max(order_by = NES, n = 2, with_ties = FALSE)

lab_right <- bind_rows(lab_right_mito, lab_right_vasc, lab_right_struct)

# =============================================================================
# 6. SHARED STYLING
# =============================================================================
fill_vals <- c(
  "Mitochondria & Metabolism"      = "#377EB8",   # Blue
  "Vascular & Angiogenesis"        = "#FF7F00",   # Red
  "Immune & Inflammation"          = "#E41A1C",   # Red
  "Translation & Ribosome"         = "#4DAF4A",   # Green
  "Neuromuscular & E-C Coupling"   = "#A65628",   # Brown
  "Structural & Contractile"       = "#984EA3",   # Purple
  "Lipid & Nutrient Handling"      = "#F781BF",   # Pink
  "Other"                          = "#999999"    # Grey
)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

y_scale_shared <- scale_y_continuous(
  limits = c(0.15, 0.80),
  labels = scales::number_format(accuracy = 0.1)
)

#----------------------------------------------------#
dummy_rows <- data.frame(
  NES       = right_lim[1],
  geneRatio = -1,
  Count     = 1,
  Category  = factor(names(fill_vals), levels = names(fill_vals))
)

# =============================================================================
# 7. LEFT PANEL (negative NES — downregulated pathways)
# =============================================================================
p_left <- ggplot(df_left,
                 aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_left,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(15, 30), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category", drop = FALSE) +
  scale_x_continuous(
    limits = left_lim,
    breaks = seq(left_lim[1], left_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(
    x = NULL,
    y = "Gene Ratio (Count / Set Size)"
  ) +
  base_theme +
  theme(legend.position = "none")

# =============================================================================
# 8. RIGHT PANEL (positive NES — upregulated pathways)
# =============================================================================
p_right <- ggplot(df_right,
                  aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_point(data = dummy_rows,
             aes(x = NES, y = geneRatio, fill = Category),
             alpha = 0, size = 0) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1, shape = 21))) +
  geom_text_repel(
    data = lab_right,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(15, 30), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category", drop = FALSE) +
  scale_x_continuous(
    limits = right_lim,
    breaks = seq(right_lim[1], right_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(
    x = NULL,
    y = NULL
    )+
  base_theme +
  theme(legend.position = "none")

# =============================================================================
# 9. COMBINE AND SAVE
# =============================================================================
p_combined <- p_left + p_right + plot_layout(widths = c(1, 1), guides = "collect")

p_combined

ggsave("FigX_oldpwr_go_enrichment.png", p_combined, width = 16, height = 8, dpi = 300)
ggsave("FigX_oldpwr_go_enrichment.pdf", p_combined, width = 15.5, height = 7.32)
#============================================================================#
##############################################################################
###############################################################################
#
#Redrawn Bubble plot (old pwr veh vs old sed veh)
#
#=============================================================#
# GO Enrichment Bubble Plot — Exercise (PWR vs Old Sed)
# Count on y-axis, geneRatio as bubble size (matches aging template)
#=============================================================#

# --- 1) Significant terms, labels ---
df_pwr_sig <- oldpwrveh_vs_oldsedveh_simpl_22Feb.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(
    label = stringr::str_replace_all(Description, "\\s+", " "),
    label = stringr::str_trunc(label, 50)
  )

# --- 2) Split ---
cut_pwr <- 1.5
x_span_pwr <- 0.75
left_lim_pwr  <- c(-(cut_pwr + x_span_pwr), -cut_pwr)
right_lim_pwr <- c(cut_pwr, cut_pwr + x_span_pwr)

df_pwr_left  <- df_pwr_sig %>% dplyr::filter(NES <= -cut_pwr)
df_pwr_right <- df_pwr_sig %>% dplyr::filter(NES >= cut_pwr)

# --- 3) Labels ---
lab_pwr_left <- bind_rows(
  df_pwr_left %>% filter(Category == "Immune & Inflammation"),
  lab_pwr_left_translation
)

lab_pwr_left_translation <- df_pwr_left %>%
  filter(Category == "Translation & Ribosome") %>%
  slice_max(order_by = abs(NES), n = 2, with_ties = FALSE)

lab_pwr_left <- bind_rows(lab_pwr_left_immune, lab_pwr_left_translation)

lab_pwr_right_mito <- df_pwr_right %>%
  filter(Category == "Mitochondria & Metabolism") %>%
  slice_max(order_by = NES, n = 7, with_ties = FALSE)

lab_pwr_right_vasc <- df_pwr_right %>%
  filter(Category == "Vascular & Angiogenesis") %>%
  slice_max(order_by = NES, n = 3, with_ties = FALSE)

lab_pwr_right_struct <- df_pwr_right %>%
  filter(Category == "Structural & Contractile") %>%
  slice_max(order_by = NES, n = 2, with_ties = FALSE)

lab_pwr_right <- bind_rows(lab_pwr_right_mito, lab_pwr_right_vasc, lab_pwr_right_struct)

# --- 4) Fill colors (same as aging plot) ---
fill_vals_pwr <- c(
  "Mitochondria & Metabolism"      = "#377EB8",
  "Vascular & Angiogenesis"        = "#FF7F00",
  "Immune & Inflammation"          = "#E41A1C",
  "Translation & Ribosome"         = "#4DAF4A",
  "Neuromuscular & E-C Coupling"   = "#A65628",
  "Structural & Contractile"       = "#984EA3",
  "Lipid & Nutrient Handling"      = "#F781BF",
  "Other"                          = "#999999"
)

# --- 5) Dummy rows for legend ---
dummy_rows_pwr <- data.frame(
  NES       = right_lim_pwr[1],
  geneRatio = 0.1,
  Count     = 1,
  Category  = factor(names(fill_vals_pwr), levels = names(fill_vals_pwr))
)

# --- 6) Shared y-axis ---
y_count_shared_pwr <- scale_y_continuous(limits = c(0, 110), breaks = seq(0, 110, by = 20))

# --- 7) Left panel ---
p_pwr_left <- ggplot(df_pwr_left,
                     aes(x = NES, y = Count, size = geneRatio, fill = Category)) +
  geom_point(alpha = 1, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_pwr_left,
    aes(label = label),
    size = 3.3, fontface = "bold",
    box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(10, 25), name = "Gene Ratio") +
  scale_fill_manual(values = fill_vals_pwr, name = "Pathway Category", drop = FALSE) +
  scale_x_continuous(limits = left_lim_pwr,
                     breaks = seq(left_lim_pwr[1], left_lim_pwr[2], by = 0.25)) +
  y_count_shared_pwr +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(legend.position = "none")

# --- 8) Right panel ---
p_pwr_right <- ggplot(df_pwr_right,
                      aes(x = NES, y = Count, size = geneRatio, fill = Category)) +
  geom_point(alpha = 1, shape = 21, color = "black", stroke = 0.6) +
  geom_point(data = dummy_rows_pwr,
             aes(x = NES, y = Count, fill = Category),
             alpha = 0, size = 0) +
  geom_text_repel(
    data = lab_pwr_right,
    aes(label = label),
    size = 3.3, fontface = "bold",
    box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(10, 25), name = "Gene Ratio") +
  scale_fill_manual(values = fill_vals_pwr, name = "Pathway Category", drop = FALSE) +
  scale_x_continuous(limits = right_lim_pwr,
                     breaks = seq(right_lim_pwr[1], right_lim_pwr[2], by = 0.25)) +
  y_count_shared_pwr +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(legend.position = "none")

# --- 9) Combine ---
p_combined_pwr <- p_pwr_left + p_pwr_right +
  plot_layout(widths = c(1, 1))

print(p_combined_pwr)

ggsave("FigX_pwr_go_enrichment_alt.pdf", p_combined_pwr,
       width = 10.5, height = 6.75)

#####################################################################################
#///////old pwr chEA3 analysis//////////#
#=======================================#

oldpwrveh_upreg_genes <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(sig_status == "Up") %>%
  pull(Symbol)

oldpwrveh_upreg_genes

write.csv(oldpwrveh_upreg_genes, "oldpwrveh_upreg_genes.csv")  

oldpwrveh_upreg_ENCODE_ChIP_seq <- read_tsv("oldpwrveh_upreg_ENCODE_ChIP_seq.tsv")

oldpwrveh_upreg_ENCODE_ChIP_seq

oldpwrveh_dwnreg_genes <- oldpwrveh_v_oldsedveh_genes_22feb %>%
  filter(sig_status == "Down") %>%
  pull(Symbol)

oldpwrveh_dwnreg_genes

write.csv(oldpwrveh_dwnreg_genes, "oldpwrveh_dwnreg_genes.csv")
###############################################################################
###############################################################################

############################################################
# Heatmap + Aging-Directed Barplot — 93 validated aging genes
# Old PWR VEH vs Old Sed VEH
############################################################
setwd("/Users/mdbruss/Documents/RStudioProjects_2/Rapa_PwR")
# 2) Load and prep all three count files
load_counts <- function(filepath) {
  df <- readxl::read_xlsx(filepath)
  colnames(df)[1] <- "Ensembl"
  df <- df %>%
    dplyr::mutate(Ensembl_noDec = stringr::str_remove(Ensembl, "\\..+")) %>%
    dplyr::select(-Ensembl)
  return(df)
}

# -------------------------------------------------------
# 0) Load and merge counts from two files
# -------------------------------------------------------
counts_veh <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
counts_pwr <- load_counts("20250219_M007853_Set02_edgeRglm_Counts_OLD_PWR_VEH-OLD_SED_VEH.xlsx")

old_sed_samples <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_pwr_samples <- c("T_08","T_10","T_11","T_41","T_43","T_45","T_46")

ordered_columns_pwr <- c(old_sed_samples, old_pwr_samples)

counts_merged_pwr <- counts_veh %>%
  dplyr::select(Ensembl_noDec, dplyr::all_of(old_sed_samples)) %>%
  dplyr::inner_join(
    counts_pwr %>% dplyr::select(Ensembl_noDec, dplyr::all_of(old_pwr_samples)),
    by = "Ensembl_noDec"
  )

# Map Ensembl -> ENTREZID
counts_merged_pwr$ENTREZID <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys      = counts_merged_pwr$Ensembl_noDec,
  keytype   = "ENSEMBL",
  column    = "ENTREZID",
  multiVals = "first"
)
counts_merged_pwr <- tidyr::drop_na(counts_merged_pwr, ENTREZID)

# -------------------------------------------------------
# 1) edgeR normalization -> logCPM
# -------------------------------------------------------
counts_matrix_pwr <- counts_merged_pwr %>%
  dplyr::select(dplyr::all_of(ordered_columns_pwr)) %>%
  as.matrix()
rownames(counts_matrix_pwr) <- as.character(counts_merged_pwr$ENTREZID)

group_pwr <- factor(c(
  rep("OLD_SED_VEH", length(old_sed_samples)),
  rep("OLD_PWR_VEH", length(old_pwr_samples))
))

dge_pwr <- edgeR::DGEList(counts = counts_matrix_pwr, group = group_pwr)
dge_pwr <- edgeR::calcNormFactors(dge_pwr, method = "TMM")
logCPM_matrix_pwr <- edgeR::cpm(dge_pwr, log = TRUE, prior.count = 1)

# -------------------------------------------------------
# 2) Compute aging-directed logFC
# -------------------------------------------------------

# Get symbols for the 93 genes present in the matrix
core_ids_present_pwr <- core_entrez[core_entrez %in% rownames(logCPM_matrix_pwr)]
entrez_to_symbol_pwr <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys    = core_ids_present_pwr,
  column  = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
core_symbols_present_pwr <- unname(na.omit(entrez_to_symbol_pwr))

# Aging direction from old vs young
aging_direction_df_pwr <- oldsedveh_v_yngsedveh_genes_06mar %>%
  dplyr::filter(Symbol %in% core_symbols_present_pwr) %>%
  dplyr::select(Symbol, aging_logFC = logFC)

# PWR VEH vs SED VEH logFC
pwr_logfc_df <- oldpwrveh_v_oldsedveh_genes_06mar %>%
  dplyr::filter(Symbol %in% core_symbols_present_pwr) %>%
  dplyr::select(Symbol, pwr_logFC = logFC, pwr_FDR = FDR)

# Join and compute aging-directed metric
pwr_bar_data <- aging_direction_df_pwr %>%
  dplyr::inner_join(pwr_logfc_df, by = "Symbol") %>%
  dplyr::mutate(
    directed_logFC = sign(aging_logFC) * pwr_logFC,
    bar_fill = ifelse(aging_logFC > 0, "#8B0000", "#00008B")
  )

# Row order: use aging (old vs young) order for direct comparison
row_order_pwr <- row_order_aging[row_order_aging %in% pwr_bar_data$Symbol]

# Now set factor levels for barplot
pwr_bar_data <- pwr_bar_data %>%
  dplyr::mutate(Symbol = factor(Symbol, levels = rev(row_order_pwr))) %>%
  dplyr::filter(!is.na(Symbol))

# -------------------------------------------------------
# 3) Subset heatmap to 93 genes and apply row order
# -------------------------------------------------------
pwr_aging_heatmap_mat <- logCPM_matrix_pwr[core_ids_present_pwr, ordered_columns_pwr, drop = FALSE]
rownames(pwr_aging_heatmap_mat) <- entrez_to_symbol_pwr[rownames(pwr_aging_heatmap_mat)]

# Keep only genes that made it through the bar data join, in the new order
row_order_pwr_heatmap <- intersect(row_order_pwr, rownames(pwr_aging_heatmap_mat))
pwr_aging_heatmap_mat <- pwr_aging_heatmap_mat[row_order_pwr_heatmap, , drop = FALSE]

# Column annotation
ann_col_pwr <- data.frame(
  Group = factor(c(
    rep("Old Sed VEH", length(old_sed_samples)),
    rep("Old PWR VEH", length(old_pwr_samples))
  ))
)
rownames(ann_col_pwr) <- ordered_columns_pwr
ann_colors_pwr <- list(Group = c("Old Sed VEH" = "#999999", "Old PWR VEH" = "#4DAF4A"))

# -------------------------------------------------------
# 4) Heatmap
# -------------------------------------------------------
Fig_pwr_aging_heat <- pheatmap::pheatmap(
  pwr_aging_heatmap_mat,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  scale             = "row",
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  row_names_side    = "left",
  fontsize_row      = 5,
  annotation_col    = ann_col_pwr,
  annotation_colors = ann_colors_pwr,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "93 Aging Axis Genes: Old PWR VEH vs Old Sed VEH"
)

pdf("FigX_heatmap_pwr_aging_axis_93genes.pdf", width = 5, height = 12)
print(Fig_pwr_aging_heat)
dev.off()

# -------------------------------------------------------
# 5) Barplot — aging-directed logFC
#    RIGHT = amplifies aging, LEFT = reverses aging
# -------------------------------------------------------
Fig_pwr_aging_bar <- ggplot(pwr_bar_data, aes(x = directed_logFC, y = Symbol, fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  coord_cartesian(xlim = c(-1.25, 0.75)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.y         = element_text(size = 5),
    axis.text.x         = element_text(size = 8),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank()
  )

print(Fig_pwr_aging_bar)

pdf("FigX_barplot_pwr_aging_directed_93genes.pdf", width = 11.25, height = 6)
print(Fig_pwr_aging_bar)
dev.off()
#########################################################################

###############################################################################
##############################################################################
#
#       #////LIPIDS: OLD PWR VEH vs OLD SED VEh//////////#
#
#
#
##############################################################################

##############################################################################
# Lipid dotplot — Old PWR VEH (OV) vs Old Sed VEH (OS)
# Signature colors preserved from aging (Old Sed vs Young Sed)
##############################################################################

# 1) Stats for OV vs OS
ov_os <- lip_log %>% dplyr::filter(Group %in% c("OV","OS"))

ov_os_stats <- purrr::map_dfr(lipid_cols, function(lip) {
  vals <- ov_os[[lip]]; grp <- ov_os$Group; tt <- t.test(vals ~ grp)
  tibble(
    Lipid            = lip,
    Mean_OV          = mean(vals[grp == "OV"], na.rm = TRUE),
    Mean_OS          = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OV_vs_OS = Mean_OV - Mean_OS,
    P_value          = tt$p.value
  )
}) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  dplyr::arrange(FDR)

# 2) Add lipid class + aging signature info from os_ys_stats
ov_os_stats <- ov_os_stats %>%
  dplyr::left_join(
    os_ys_stats %>% dplyr::select(Lipid, LipidClass, IsSignature, Log2_FC_OS_vs_YS),
    by = "Lipid"
  ) %>%
  dplyr::filter(!is.na(LipidClass)) %>%
  dplyr::mutate(
    neg_log10_p = -log10(P_value),
    # Color based on AGING direction, not OV vs OS direction
    sig_fill = case_when(
      IsSignature & Log2_FC_OS_vs_YS > 0 ~ "age_up",
      IsSignature & Log2_FC_OS_vs_YS < 0 ~ "age_down",
      TRUE ~ "not_sig"
    )
  )

plot_df_ov <- ov_os_stats %>%
  mutate(LipidClass = factor(LipidClass, levels = class_order))

# 3) Dotplot
pwr_lipids_fig <- ggplot(plot_df_ov, aes(x = Log2_FC_OV_vs_OS, y = LipidClass)) +
  # Non-signature: color-coded by class
  geom_jitter(
    data = subset(plot_df_ov, sig_fill == "not_sig"),
    aes(size = neg_log10_p, color = LipidClass),
    width = 0, height = 0.22, alpha = 1, shape = 16
  ) +
  # Signature: up with aging (dark red)
  geom_jitter(
    data = subset(plot_df_ov, sig_fill == "age_up"),
    aes(size = neg_log10_p),
    width = 0, height = 0.22,
    shape = 21, fill = "#8B0000", color = "black", stroke = 0.6, alpha = 0.9
  ) +
  # Signature: down with aging (dark blue)
  geom_jitter(
    data = subset(plot_df_ov, sig_fill == "age_down"),
    aes(size = neg_log10_p),
    width = 0, height = 0.22,
    shape = 21, fill = "#00008B", color = "black", stroke = 0.6, alpha = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_y_discrete(limits = class_order) +
  scale_color_manual(values = class_colors, guide = "none") +
  scale_size_continuous(
    range = c(3, 9),
    name = expression(-log[10](P))
  ) +
  labs(
    x = "log2 Fold Change (Old PWR VEH − Old Sed VEH)",
    y = "Lipid Class",
    title = "Old PWR VEH vs Old Sed VEH Lipidomics",
    subtitle = "Dark red = aging signature (up) | Dark blue = aging signature (down) | Circle size = significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

print(pwr_lipids_fig)

# Save as PDF
ggsave(
  "Fig2d_oldpwrveh_lipid_classes.pdf",
  pwr_lipids_fig,
  width = 13.5,
  height = 6.6
)
#===========================================================================#
#############################################################################

library(ggridges)

# --- Same ordering and clip as OS vs YS ---
x_clip_ov <- 2.0

ridge_df_ov <- ov_os_stats %>%
  filter(!is.na(LipidClass), !LipidClass %in% c("Sph", "PA", "PG", "Ubiquinone")) %>%
  group_by(LipidClass) %>%
  mutate(
    n_class     = n(),
    n_sig_class = sum(IsSignature)
  ) %>%
  ungroup() %>%
  mutate(
    x_display_ov = pmin(pmax(Log2_FC_OV_vs_OS, -x_clip_ov), x_clip_ov),
    # Same dot coloring as aging ridge: red = aging-up, blue = aging-down
    dot_group_ov = case_when(
      IsSignature & Log2_FC_OS_vs_YS > 0 ~ "Sig Up",
      IsSignature & Log2_FC_OS_vs_YS < 0 ~ "Sig Down",
      TRUE ~ "NS"
    ),
    LipidClassLabel_ov = LipidClass
  )

# --- Same factor order as aging ridge ---
label_order_ov <- ridge_df_ov %>%
  distinct(LipidClass, LipidClassLabel_ov) %>%
  mutate(LipidClass = factor(LipidClass, levels = ora_order)) %>%
  arrange(LipidClass) %>%
  pull(LipidClassLabel_ov)

ridge_df_ov <- ridge_df_ov %>%
  mutate(LipidClassLabel_ov = factor(LipidClassLabel_ov, levels = rev(label_order_ov)))

# --- Plot: all ridges grey, dots match aging ridge exactly ---
class_ridge_ov <- ggplot(ridge_df_ov, aes(x = x_display_ov, y = LipidClassLabel_ov)) +
  geom_density_ridges(
    fill = "grey85",
    alpha = 0.7, scale = 0.9, rel_min_height = 0.01,
    color = "black", linewidth = 0.3
  ) +
  geom_jitter(
    data = subset(ridge_df_ov, dot_group_ov == "NS"),
    color = "grey60",
    height = 0.15, width = 0, size = 1.5, alpha = 0.6
  ) +
  geom_jitter(
    data = subset(ridge_df_ov, dot_group_ov == "Sig Up"),
    color = "#8B0000",
    height = 0.15, width = 0, size = 3, alpha = 0.9
  ) +
  geom_jitter(
    data = subset(ridge_df_ov, dot_group_ov == "Sig Down"),
    color = "#00008B",
    height = 0.15, width = 0, size = 3, alpha = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(-2, 1.75), breaks = seq(-2, 1.75, by = 0.5)) +
  labs(
    x = NULL,
    y = "Lipid Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(class_ridge_ov)

ggsave(
  "FigX_oldpwrveh_lipid_class_ridge.pdf",
  class_ridge_ov,
  width = 9.9,
  height = 5.4
)
############################################################
############################################################
#
#===========================================================================#
# FULL PIPELINE: Lipid Class Directed logFC + Colored Ridgeplot
#===========================================================================#

library(dplyr); library(stringr); library(ggplot2); library(ggridges)

# ===== PART 1: Class-level aging-directed analysis =====

# --- 1) Get aging direction per lipid from OS vs YS ---
lip_aging <- os_ys_stats %>%
  dplyr::select(Lipid, LipidClass, aging_logFC = Log2_FC_OS_vs_YS) %>%
  dplyr::mutate(aging_direction = sign(aging_logFC))

# --- 2) Build treatment logFCs ---
add_class <- function(df, fc_col) {
  df %>%
    dplyr::mutate(LipidClass = str_trim(str_extract(Lipid, "^[^\\(]+"))) %>%
    dplyr::select(Lipid, LipidClass, treat_logFC = !!sym(fc_col))
}

treat_dfs <- list(
  "PWR"          = add_class(ov_os_stats,  "Log2_FC_OV_vs_OS"),
  "PWR + iRapa"  = add_class(oir_os_stats, "Log2_FC_vs_OS"),
  "PWR + fRapa"  = add_class(ofr_os_stats, "Log2_FC_vs_OS")
)

# --- 3) Compute directed logFC per lipid per condition ---
lip_directed <- purrr::imap_dfr(treat_dfs, function(df, condition) {
  df %>%
    inner_join(lip_aging %>% dplyr::select(Lipid, aging_direction), by = "Lipid") %>%
    mutate(
      directed_logFC = treat_logFC * aging_direction,
      Condition = condition
    )
})

# --- 4) Drop small classes ---
classes_to_keep <- lip_directed %>%
  distinct(Lipid, LipidClass) %>%
  count(LipidClass) %>%
  filter(n >= 5) %>%
  pull(LipidClass)

lip_directed <- lip_directed %>%
  filter(LipidClass %in% classes_to_keep)

# --- 5) Summarize: one-sample Wilcoxon (uncorrected) ---
class_summary <- lip_directed %>%
  group_by(LipidClass, Condition) %>%
  summarise(
    n = n(),
    mean_directed = mean(directed_logFC, na.rm = TRUE),
    se = sd(directed_logFC, na.rm = TRUE) / sqrt(n()),
    p_value = tryCatch(
      wilcox.test(directed_logFC, mu = 0)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    Condition = factor(Condition, levels = c("PWR + fRapa", "PWR + iRapa", "PWR"))
  )

# --- 6) Order by PWR reversal ---
class_order_plot <- class_summary %>%
  filter(Condition == "PWR") %>%
  arrange(desc(mean_directed)) %>%
  pull(LipidClass) %>%
  as.character()

class_summary <- class_summary %>%
  mutate(LipidClass = factor(as.character(LipidClass), levels = class_order_plot))

# --- Check results ---
cat("\n=== PWR class-level results ===\n")
as.data.frame(class_summary %>%
                filter(Condition == "PWR") %>%
                dplyr::select(LipidClass, mean_directed, p_value, sig_label) %>%
                arrange(p_value))


# ===== PART 2: Ridgeplot with colored ridges =====

# --- Build ridge dataframe (OV vs OS) ---
x_clip_ov <- 2.0

ridge_df_ov <- ov_os_stats %>%
  filter(!is.na(LipidClass), !LipidClass %in% c("Sph", "PA", "PG", "Ubiquinone")) %>%
  group_by(LipidClass) %>%
  mutate(
    n_class     = n(),
    n_sig_class = sum(IsSignature)
  ) %>%
  ungroup() %>%
  mutate(
    x_display_ov = pmin(pmax(Log2_FC_OV_vs_OS, -x_clip_ov), x_clip_ov),
    dot_group_ov = case_when(
      IsSignature & Log2_FC_OS_vs_YS > 0 ~ "Sig Up",
      IsSignature & Log2_FC_OS_vs_YS < 0 ~ "Sig Down",
      TRUE ~ "NS"
    ),
    LipidClassLabel_ov = LipidClass
  )

# --- Factor order (match bar plot order or use ORA order if available) ---
if (exists("ora_order")) {
  label_order_ov <- ridge_df_ov %>%
    distinct(LipidClass, LipidClassLabel_ov) %>%
    mutate(LipidClass = factor(LipidClass, levels = ora_order)) %>%
    arrange(LipidClass) %>%
    pull(LipidClassLabel_ov)
} else {
  label_order_ov <- ridge_df_ov %>%
    distinct(LipidClassLabel_ov) %>%
    pull(LipidClassLabel_ov)
}

ridge_df_ov <- ridge_df_ov %>%
  mutate(LipidClassLabel_ov = factor(LipidClassLabel_ov, levels = rev(label_order_ov)))

# --- Identify significant classes and assign to ridge groups ---
directed_sig_classes <- class_summary %>%
  filter(Condition == "PWR", p_value < 0.05) %>%
  pull(LipidClass) %>%
  as.character()

ridge_df_ov <- ridge_df_ov %>%
  mutate(
    ridge_group_ov = ifelse(LipidClass %in% directed_sig_classes,
                            as.character(LipidClass), "Other")
  )

# --- Red/pink gradient by significance ---
directed_sig_df <- class_summary %>%
  filter(Condition == "PWR", p_value < 0.05) %>%
  arrange(p_value) %>%
  mutate(LipidClass = as.character(LipidClass))

pink_gradient <- colorRampPalette(c("#8B0000", "#F4A0A0"))(nrow(directed_sig_df))
directed_color_map <- setNames(pink_gradient, directed_sig_df$LipidClass)
ridge_fill_vals_ov <- c(directed_color_map, "Other" = "grey85")

# --- Pull significance labels for annotation ---
pwr_class_sig <- class_summary %>%
  filter(Condition == "PWR", sig_label != "") %>%
  dplyr::select(LipidClass, sig_label) %>%
  mutate(LipidClassLabel_ov = factor(
    as.character(LipidClass),
    levels = levels(ridge_df_ov$LipidClassLabel_ov)
  ))

# --- Update fill: white for non-sig ridges ---
ridge_fill_vals_ov <- c(directed_color_map, "Other" = "grey")

# --- Plot ---
class_ridge_ov <- ggplot(ridge_df_ov, aes(x = x_display_ov, y = LipidClassLabel_ov)) +
  geom_density_ridges(
    aes(fill = ridge_group_ov),
    alpha = 0.7, scale = 0.9, rel_min_height = 0.01,
    color = "black", linewidth = 0.3
  ) +
  geom_jitter(
    data = subset(ridge_df_ov, dot_group_ov == "NS"),
    color = "grey40",
    height = 0.15, width = 0, size = 1.5, alpha = 0.6
  ) +
  geom_jitter(
    data = subset(ridge_df_ov, dot_group_ov == "Sig Up"),
    color = "#8B0000",
    height = 0.15, width = 0, size = 3, alpha = 0.9
  ) +
  geom_jitter(
    data = subset(ridge_df_ov, dot_group_ov == "Sig Down"),
    color = "#00008B",
    height = 0.15, width = 0, size = 3, alpha = 0.9
  ) +
  geom_text(
    data = pwr_class_sig,
    aes(x = 1.5, y = LipidClassLabel_ov, label = sig_label),
    inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = ridge_fill_vals_ov, guide = "none") +
  scale_x_continuous(limits = c(-1.3, 0.75), breaks = seq(-2, 1.75, by = 0.5)) +
  labs(x = NULL, y = "Lipid Class") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(class_ridge_ov)

ggsave("FigX_oldpwrveh_lipid_class_ridge_directed.pdf", class_ridge_ov,
       width = 9.9, height = 5.4)
##############################################################
##############################################################


#########################################################################
#############################################################################
## =========================
## Heatmap + Aging-Directed Bar Plot — 59-lipid signature
## Old PWR VEH (OV) vs Old Sed VEH (OS)
## =========================

# 1) Stats for OV vs OS (log2FC = OV − OS)
ov_os <- lip_log %>% dplyr::filter(Group %in% c("OV","OS"))

ov_os_stats <- purrr::map_dfr(lipid_cols, function(lip) {
  vals <- ov_os[[lip]]; grp <- ov_os$Group; tt <- t.test(vals ~ grp)
  tibble::tibble(
    Lipid            = lip,
    Mean_OV          = mean(vals[grp == "OV"], na.rm = TRUE),
    Mean_OS          = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OV_vs_OS = Mean_OV - Mean_OS,
    P_value          = tt$p.value
  )
}) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  dplyr::arrange(FDR)

# 2) Heatmap — same lipid order as old vs young (lip_order_sig)
heatmap_data_ovos <- lip_log %>%
  dplyr::filter(Group %in% c("OS","OV")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("OS","OV"))) %>%
  dplyr::arrange(Group) %>%
  dplyr::select(Sample, Group, dplyr::all_of(lip_order_sig))

mat_ovos <- heatmap_data_ovos %>%
  dplyr::select(-Sample, -Group) %>%
  as.matrix()
rownames(mat_ovos) <- heatmap_data_ovos$Sample

mat_ovos_z <- scale(mat_ovos)

ann_col_ovos <- data.frame(Group = heatmap_data_ovos$Group)
rownames(ann_col_ovos) <- heatmap_data_ovos$Sample
ann_cols_ov <- list(Group = c(OS = "#4F4F4F", OV = "#4DAF4A"))

hm_ovos <- pheatmap(
  t(mat_ovos_z),
  annotation_col    = ann_col_ovos,
  annotation_colors = ann_cols_ov,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row   = lip_order_sig, row_names_side = "left",
  fontsize = 6, color = pal, breaks = seq(-1, 1, length.out = 101),
  main = "59-lipid Aging Signature: Old Sed VEH (left) vs Old PWR VEH (right)"
)
print(hm_ovos)

pdf("FigX_heatmap_lipid_OV_vs_OS_59.pdf", width = 6.75, height = 6.6)
print(hm_ovos)
dev.off()

# 3) Aging-directed bar plot
#    Join aging logFC with OV vs OS logFC, compute directed metric
bar_ovos_df <- os_ys_stats %>%
  dplyr::filter(Lipid %in% signature59) %>%
  dplyr::select(Lipid, aging_logFC = Log2_FC_OS_vs_YS) %>%
  dplyr::inner_join(
    ov_os_stats %>% dplyr::select(Lipid, ov_logFC = Log2_FC_OV_vs_OS),
    by = "Lipid"
  ) %>%
  dplyr::mutate(
    directed_logFC = sign(aging_logFC) * ov_logFC,
    bar_fill = ifelse(aging_logFC > 0, "#8B0000", "#00008B"),
    Lipid = factor(Lipid, levels = lip_order_sig)
  )

bar_ovos <- ggplot(bar_ovos_df, aes(x = directed_logFC, y = fct_rev(Lipid), fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  coord_cartesian(xlim = c(-1.25, 0.75)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(bar_ovos)

pdf("FigX_barplot_lipid_OV_vs_OS_aging_directed_59.pdf", width = 11.25, height = 4.5)
print(bar_ovos)
dev.off()
###############################################################################
##############################################################################

#/////Targeted LSEA///////////////#

## =========================
## Targeted LSEA — 59-lipid aging signature
## Aging-directed ranking (matches RNA logic)
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)
library(patchwork)
set.seed(1)

# ---- Aging direction from OS vs YS ----
aging_lip_direction <- os_ys_stats %>%
  dplyr::select(Lipid, aging_logFC = Log2_FC_OS_vs_YS) %>%
  dplyr::mutate(aging_direction = sign(aging_logFC))

# ---- Signature ----
lipid_age_sig <- signature59

# ---- Restrict to groups we need ----
lipid_log_filtered <- lip_log %>%
  dplyr::filter(Group %in% c("OS","OV","OIR","OFR")) %>%
  dplyr::select(Sample, Group, dplyr::all_of(lipid_cols))

# ---- Helper: compute aging-directed rank metric for a two-group contrast ----
.compute_aging_directed_ranks <- function(df, g_intervention, g_control) {
  sub <- df %>% dplyr::filter(Group %in% c(g_intervention, g_control))
  
  # Per-lipid stats: logFC = intervention - control
  stats <- purrr::map_dfr(lipid_cols, function(lip) {
    vals <- sub[[lip]]; grp <- sub$Group
    tt <- t.test(vals ~ grp)
    # Ensure direction is intervention - control
    mean_int  <- mean(vals[grp == g_intervention], na.rm = TRUE)
    mean_ctrl <- mean(vals[grp == g_control], na.rm = TRUE)
    tibble::tibble(
      Lipid   = lip,
      LogFC   = mean_int - mean_ctrl,
      P_value = tt$p.value
    )
  })
  
  # Join aging direction and compute rank metric (same as RNA)
  stats_ranked <- stats %>%
    dplyr::left_join(aging_lip_direction, by = "Lipid") %>%
    dplyr::mutate(
      p_safe      = pmax(P_value, 1e-300),
      rank_metric = aging_direction * LogFC * ((-log10(p_safe))^0.25)
    ) %>%
    dplyr::filter(!is.na(rank_metric)) %>%
    dplyr::arrange(dplyr::desc(rank_metric))
  
  ranks <- stats::setNames(stats_ranked$rank_metric, stats_ranked$Lipid)
  list(stats = stats, ranks = ranks)
}

# ---- Helper: enrichment (top) + ranked-stat bars (bottom) ----
.make_lsea_panels <- function(ranks, pathway, title_text) {
  p_top <- fgsea::plotEnrichment(pathway, ranks) +
    ggtitle(title_text) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7))
  
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
    labs(x = "Lipid rank", y = "Aging-directed\nstatistic") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.margin = margin(2, 2, 2, 2))
  
  p_top / p_bottom + plot_layout(heights = c(2, 1))
}

# ---- Run the three contrasts ----
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
  
  out <- .compute_aging_directed_ranks(lipid_log_filtered, g_high, g_low)
  ranks <- out$ranks
  
  # fgsea on the 59-lipid signature
  lip_sets <- list(AgingSignature = lipid_age_sig)
  fg <- fgsea::fgsea(pathways = lip_sets, stats = ranks, nperm = 1000)
  
  nes_rows[[nm]] <- tibble::tibble(
    Contrast = nm,
    NES      = fg$NES[match("AgingSignature", fg$pathway)],
    pvalue   = fg$pval[match("AgingSignature", fg$pathway)],
    padj     = fg$padj[match("AgingSignature", fg$pathway)]
  )
  
  panels_list[[nm]] <- .make_lsea_panels(
    ranks      = ranks,
    pathway    = lip_sets$AgingSignature,
    title_text = nm
  )
  
  results_list[[nm]] <- list(stats = out$stats, ranks = ranks, fgsea = fg)
}

# ---- LSEA panels (side by side) ----
panels_present <- panels_list[!vapply(panels_list, is.null, logical(1))]
fig_lsea <- Reduce(`|`, panels_present) + plot_layout(guides = "collect") &
  theme(plot.title = element_text(hjust = 0.5))

print(fig_lsea)

# ---- NES summary bar plot ----
nes_df <- bind_rows(nes_rows) %>%
  mutate(
    Direction = ifelse(NES < 0, "Reverses aging", "Amplifies aging"),
    label = paste0("p=", formatC(pvalue, format = "g", digits = 2))
  )

fig_nes <- ggplot(nes_df, aes(x = NES, y = factor(Contrast, levels = rev(names(contrasts))))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = Direction), width = 0.6) +
  geom_text(aes(label = label), hjust = ifelse(nes_df$NES < 0, 1.1, -0.1), size = 3) +
  scale_fill_manual(values = c("Reverses aging" = "#2166AC", "Amplifies aging" = "#8B0000")) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Lipid Aging Signature: Targeted LSEA Summary",
    subtitle = "Aging-directed ranking (sign × logFC × significance)"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank())

print(fig_nes)

# ---- Print NES table ----
cat("\n=== LSEA Results ===\n")
print(nes_df %>% select(Contrast, NES, pvalue, padj))
#------------------------------------------------------#

fig_lsea_ov <- panels_list[["OLD PWR vs OLD SED"]]
print(fig_lsea_ov)

ggsave("FigX_LSEA_OV_vs_OS.pdf", fig_lsea_ov, width = 6.47, height = 3)
###########################################################################
#############################################################################
#
#
#
#
#//////////Metabolomics Analysis///////////////#
#
#
#
#
##############################################################################
# METABOLOMICS — Old PWR VEH (OV) vs Old Sed VEH (OS)
##############################################################################

##############################################################################
# METABOLOMICS — Old PWR VEH (OV) vs Old Sed VEH (OS)
# Focused on the 17 aging-signature metabolites
##############################################################################

# ---------- 1) Stats: OV vs OS ----------
metab_ovos <- metab_log %>% dplyr::filter(Group %in% c("OV","OS"))

metab_ovos_stats <- purrr::map_dfr(metab_cols, function(met) {
  vals <- metab_ovos[[met]]; grp <- metab_ovos$Group
  tt <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
  tibble(
    Metabolite        = met,
    Mean_OS           = mean(vals[grp == "OS"], na.rm = TRUE),
    Mean_OV           = mean(vals[grp == "OV"], na.rm = TRUE),
    Log2_FC_OV_vs_OS  = mean(vals[grp == "OV"], na.rm = TRUE) - mean(vals[grp == "OS"], na.rm = TRUE),
    P_value           = if (is.null(tt)) NA_real_ else tt$p.value
  )
}) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  dplyr::arrange(FDR)

# Add aging signature info
metab_ovos_stats <- metab_ovos_stats %>%
  dplyr::left_join(
    metab_osys_stats %>% dplyr::select(Metabolite, IsSignature, Log2_FC_OS_vs_YS),
    by = "Metabolite"
  )

# ---------- 2) Volcano: all grey, aging-signature in dark red ----------
metab_ovos_sig_df <- metab_ovos_stats %>% dplyr::filter(IsSignature)

metab_ovos_xlim <- max(1, quantile(abs(metab_ovos_stats$Log2_FC_OV_vs_OS), 0.99, na.rm = TRUE))

metab_ovos_volcano <- ggplot(metab_ovos_stats, aes(x = Log2_FC_OV_vs_OS, y = -log10(P_value))) +
  geom_point(color = "grey65", size = 2) +
  geom_point(
    data = metab_ovos_sig_df,
    shape = 21, fill = "#8B0000", color = "black", stroke = 0.6, size = 5
  ) +
  ggrepel::geom_text_repel(
    data = metab_ovos_sig_df %>% filter(P_value < 0.05),
    aes(label = Metabolite),
    size = 3, max.overlaps = 100,
    box.padding = 0.4, point.padding = 0.3, seed = 1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(-3.5, 2.5), breaks = seq(-3.5, 2.5, by = 1)) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.0, by = 1)) +
  labs(
    x = NULL,
    y = expression(-log[10](P))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

print(metab_ovos_volcano)

ggsave(
  "FigX_oldpwrveh_metabo_volcano_plot.pdf",
  metab_ovos_volcano,
  width = 9.6,
  height = 6.0
)

# ---------- 3) Heatmap: OV vs OS (17 aging-signature metabolites, aging row order) ----------
metab_ovos_hm_data <- metab_log %>%
  dplyr::filter(Group %in% c("OS","OV")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("OS","OV"))) %>%
  dplyr::arrange(Group) %>%
  dplyr::select(Sample, Group, dplyr::all_of(metab_row_order))

metab_ovos_mat <- metab_ovos_hm_data %>% dplyr::select(-Sample, -Group) %>% as.matrix()
rownames(metab_ovos_mat) <- metab_ovos_hm_data$Sample
metab_ovos_mat_z <- scale(metab_ovos_mat)

metab_ovos_ann_col <- data.frame(Group = metab_ovos_hm_data$Group)
rownames(metab_ovos_ann_col) <- metab_ovos_hm_data$Sample
metab_ovos_ann_colors <- list(Group = c(OS = "#4F4F4F", OV = "#4DAF4A"))

metab_ovos_hm_alt <- pheatmap::pheatmap(
  t(metab_ovos_mat_z),
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  legend          = FALSE,
  annotation_col  = NA,
  labels_row      = metab_row_order,
  fontsize_row    = 5,
  color = metab_pal, breaks = seq(-1, 1, length.out = 101)
)
print(metab_ovos_hm_alt)

pdf("FigX_heatmap_metabolite_OV_vs_OS.pdf", width = 4.4, height = 2.8)
print(metab_ovos_hm_alt)
dev.off()

# ---------- 4) Aging-directed barplot (17 signature metabolites, aging row order) ----------
metab_ovos_bar_data <- metab_osys_stats %>%
  dplyr::filter(Metabolite %in% metab_age_signature) %>%
  dplyr::select(Metabolite, aging_logFC = Log2_FC_OS_vs_YS) %>%
  dplyr::inner_join(
    metab_ovos_stats %>% dplyr::select(Metabolite, ov_logFC = Log2_FC_OV_vs_OS),
    by = "Metabolite"
  ) %>%
  dplyr::mutate(
    directed_logFC = sign(aging_logFC) * ov_logFC,
    bar_fill = ifelse(aging_logFC > 0, "#8B0000", "#00008B"),
    Metabolite = factor(Metabolite, levels = metab_row_order)
  )

metab_ovos_bar <- ggplot(metab_ovos_bar_data,
                         aes(x = directed_logFC, y = fct_rev(Metabolite), fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  coord_cartesian(xlim = c(-2.25, 1.25)) +
  labs(
    x = "Aging-Directed log2FC (Old PWR VEH vs Old Sed VEH)\n\u2190 Reverses Aging | Amplifies Aging \u2192",
    y = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(metab_ovos_bar)

pdf("FigX_barplot_metabolite_OV_vs_OS_aging_directed.pdf", width = 7.5, height = 3)
print(metab_ovos_bar)
dev.off()

# ---------- 5) Reversal summary ----------
metab_reversal <- metab_ovos_stats %>%
  filter(IsSignature) %>%
  mutate(
    IsSig_OV     = P_value < 0.05,
    IsReversal   = IsSig_OV & sign(Log2_FC_OV_vs_OS) == -sign(Log2_FC_OS_vs_YS),
    IsConcordant = IsSig_OV & sign(Log2_FC_OV_vs_OS) == sign(Log2_FC_OS_vs_YS)
  )

cat("\n=== Metabolite reversal summary ===\n")
cat("Total aging-signature metabolites:", nrow(metab_reversal), "\n")
cat("Sig in OV vs OS:", sum(metab_reversal$IsSig_OV), "\n")
cat("Reversed:", sum(metab_reversal$IsReversal), "\n")
cat("Concordant:", sum(metab_reversal$IsConcordant), "\n")
###############################################################

#===========================================================================#
# METABOLITE LSEA — 3 contrasts vs Old Sed (mirrors lipid LSEA exactly)
#===========================================================================#

# ---- Aging direction from OS vs YS ----
aging_metab_direction <- metab_osys_stats %>%
  dplyr::select(Metabolite, aging_logFC = Log2_FC_OS_vs_YS) %>%
  dplyr::mutate(aging_direction = sign(aging_logFC))

metab_sig_set <- metab_age_signature  # your 24-feature signature vector

metab_log_filtered <- metab_log %>%
  dplyr::filter(Group %in% c("OS","OV","OIR","OFR")) %>%
  dplyr::select(Sample, Group, dplyr::all_of(metab_cols))

# ---- Helper: aging-directed rank metric for metabolites ----
.compute_metab_aging_directed_ranks <- function(df, g_intervention, g_control) {
  sub <- df %>% dplyr::filter(Group %in% c(g_intervention, g_control))
  
  stats <- purrr::map_dfr(metab_cols, function(met) {
    vals <- sub[[met]]; grp <- sub$Group
    tt <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
    mean_int  <- mean(vals[grp == g_intervention], na.rm = TRUE)
    mean_ctrl <- mean(vals[grp == g_control], na.rm = TRUE)
    tibble::tibble(
      Metabolite = met,
      LogFC      = mean_int - mean_ctrl,
      P_value    = if (is.null(tt)) NA_real_ else tt$p.value
    )
  })
  
  stats_ranked <- stats %>%
    dplyr::left_join(aging_metab_direction, by = "Metabolite") %>%
    dplyr::mutate(
      p_safe      = pmax(P_value, 1e-300),
      rank_metric = aging_direction * LogFC * ((-log10(p_safe))^0.25)
    ) %>%
    dplyr::filter(!is.na(rank_metric)) %>%
    dplyr::arrange(dplyr::desc(rank_metric))
  
  ranks <- stats::setNames(stats_ranked$rank_metric, stats_ranked$Metabolite)
  list(stats = stats, ranks = ranks)
}

# ---- Run 3 contrasts ----
contrasts_metab <- list(
  `OLD PWR vs OLD SED`      = c("OV",  "OS"),
  `OLD PWR IRAP vs OLD SED` = c("OIR", "OS"),
  `OLD PWR FRAP vs OLD SED` = c("OFR", "OS")
)

metab_results_list <- list()
metab_panels_list  <- list()
metab_nes_rows     <- list()

for (nm in names(contrasts_metab)) {
  g_high <- contrasts_metab[[nm]][1]
  g_low  <- contrasts_metab[[nm]][2]
  
  out   <- .compute_metab_aging_directed_ranks(metab_log_filtered, g_high, g_low)
  ranks <- out$ranks
  
  met_sets <- list(AgingSignature = metab_sig_set)
  fg <- fgsea::fgsea(pathways = met_sets, stats = ranks)
  
  metab_nes_rows[[nm]] <- tibble::tibble(
    Contrast = nm,
    NES      = fg$NES[match("AgingSignature", fg$pathway)],
    pvalue   = fg$pval[match("AgingSignature", fg$pathway)],
    padj     = fg$padj[match("AgingSignature", fg$pathway)]
  )
  
  # Reuse your lipid panel helper — works identically for metabolites
  metab_panels_list[[nm]] <- .make_lsea_panels(
    ranks      = ranks,
    pathway    = met_sets$AgingSignature,
    title_text = nm
  )
  
  metab_results_list[[nm]] <- list(stats = out$stats, ranks = ranks, fgsea = fg)
}

metab_nes_df <- bind_rows(metab_nes_rows)
print(metab_nes_df)

#===========================================================================#
# RNA LSEA — 3 contrasts vs Old Sed (testing 93 aging signature genes)
#===========================================================================#

# ---- Aging direction from Old Sed vs Young ----
aging_rna_direction <- oldsedveh_v_yngsedveh_genes_06mar %>%
  dplyr::distinct(Symbol, .keep_all = TRUE) %>%
  dplyr::select(Symbol, aging_logFC = logFC) %>%
  dplyr::mutate(aging_direction = sign(aging_logFC))

rna_sig_set <- validated_aging_genes$Symbol  # your 93-gene signature

# ---- Helper: aging-directed rank metric for RNA ----
.compute_rna_aging_directed_ranks <- function(treatment_df) {
  treatment_df %>%
    dplyr::distinct(Symbol, .keep_all = TRUE) %>%
    dplyr::select(Symbol, treat_logFC = logFC, treat_FDR = FDR) %>%
    dplyr::left_join(aging_rna_direction, by = "Symbol") %>%
    dplyr::mutate(
      p_safe      = pmax(treat_FDR, 1e-300),
      rank_metric = aging_direction * treat_logFC * ((-log10(p_safe))^0.25)
    ) %>%
    dplyr::filter(!is.na(rank_metric)) %>%
    dplyr::arrange(dplyr::desc(rank_metric)) %>%
    { stats::setNames(.$rank_metric, .$Symbol) }
}

rna_dfs <- list(
  `OLD PWR vs OLD SED`      = oldpwrveh_v_oldsedveh_genes_06mar,
  `OLD PWR IRAP vs OLD SED` = oldpwrirap_v_oldsedveh_genes_02mar,
  `OLD PWR FRAP vs OLD SED` = oldpwrfrap_v_oldsedveh_genes_02mar
)

rna_nes_rows   <- list()
rna_panels_list <- list()

for (nm in names(rna_dfs)) {
  ranks <- .compute_rna_aging_directed_ranks(rna_dfs[[nm]])
  
  rna_sets <- list(AgingSignature = rna_sig_set)
  fg <- fgsea::fgsea(pathways = rna_sets, stats = ranks)
  
  rna_nes_rows[[nm]] <- tibble::tibble(
    Contrast = nm,
    NES      = fg$NES[match("AgingSignature", fg$pathway)],
    pvalue   = fg$pval[match("AgingSignature", fg$pathway)],
    padj     = fg$padj[match("AgingSignature", fg$pathway)]
  )
  
  rna_panels_list[[nm]] <- .make_lsea_panels(
    ranks      = ranks,
    pathway    = rna_sets$AgingSignature,
    title_text = nm
  )
}

rna_nes_df <- bind_rows(rna_nes_rows)
print(rna_nes_df)

#===========================================================================#
# COMBINED NES SUMMARY — all 3 layers, all 3 contrasts
#===========================================================================#

all_nes_df <- bind_rows(
  rna_nes_df   %>% mutate(Block = "Transcriptomics"),
  metab_nes_df %>% mutate(Block = "Metabolomics"),
  bind_rows(nes_rows) %>% mutate(Block = "Lipidomics")   # from your existing lipid LSEA
) %>%
  mutate(
    Contrast = factor(Contrast, levels = c(
      "OLD PWR vs OLD SED",
      "OLD PWR IRAP vs OLD SED",
      "OLD PWR FRAP vs OLD SED"
    )),
    Block = factor(Block, levels = c("Transcriptomics","Lipidomics","Metabolomics")),
    sig_label = dplyr::case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  )

p_nes_combined <- ggplot(all_nes_df,
                         aes(x = Contrast, y = NES, fill = Contrast)) +
  geom_col(color = "black", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label,
                vjust = ifelse(NES >= 0, -0.3, 1.3)),
            size = 5) +
  facet_wrap(~Block) +
  scale_fill_manual(values = c(
    "OLD PWR vs OLD SED"      = "#33a02c",
    "OLD PWR IRAP vs OLD SED" = "#6a3d9a",
    "OLD PWR FRAP vs OLD SED" = "#ff7f00"
  )) +
  labs(
    y = "NES (aging-directed)",
    x = "",
    title = "Aging Signature Reversal: FeatureSEA across 3 molecular layers"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 30, hjust = 1),
    strip.text       = element_text(face = "bold")
  )

print(p_nes_combined)
ggsave("Fig5_LSEA_NES_combined.pdf", p_nes_combined, width = 12, height = 5)
###################################################################

# What groups exist in metab_log?
metab_log %>% count(Group)

# What samples are in each group?
metab_log %>% filter(Group == "OFR") %>% pull(Sample)
metab_log %>% filter(Group == "OS") %>% pull(Sample)

# How many samples per contrast?
metab_log_filtered %>% count(Group)

# Look at the ranking for FRAP vs OS
ofr_out <- .compute_metab_aging_directed_ranks(metab_log_filtered, "OFR", "OS")

# Check a few known aging metabolites
ofr_out$stats %>%
  filter(Metabolite %in% metab_age_signature) %>%
  arrange(LogFC) %>%
  print(n = 24)
