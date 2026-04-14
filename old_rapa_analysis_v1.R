
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
#Identify old Rapa DEG ----#
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
####################################################################

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
oldsedfrap_v_oldsedveh_genes_08feb <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_FRAP-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_SED_FRAP-OLD_SED_VEH_logFC",
  FDR_col   = "OLD_SED_FRAP-OLD_SED_VEH_FDR"
)

#remove duplicate ENTREZID and keep most significant
oldsedfrap_v_oldsedveh_genes_08feb <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()

head(oldsedfrap_v_oldsedveh_genes_08feb)

oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(FDR < 0.05)

oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(Symbol == "Cxcl10")
#==========================================#

# Treg signature genes
treg_markers <- c("Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18", "Itgae", "Nrp1", "Icos")

# M2 macrophage markers
m2_mac_markers <- c("Arg1", "Mrc1", "Cd163", "Retnla", "Chil3", "Mgl2", "Ccl22", "Tgfb1", "Il10")

# T-cell exhaustion markers
exhaustion_markers <- c("Pdcd1", "Ctla4", "Havcr2", "Lag3", "Tigit", "Tox", "Entpd1", "Cd244a")

# MDSC-associated markers
mdsc_markers <- c("Itgam", "Ly6g", "Ly6c1", "Arg1", "Nos2", "Il10", "Tgfb1", "Ptgs2")

# Combined ARHGAP4-associated immunosuppressive signature
arhgap4_immune_signature <- unique(c(treg_markers, m2_mac_markers, exhaustion_markers, mdsc_markers))

oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(Symbol %in% arhgap4_immune_signature)

####################################################
############################################

# Create significance and direction column (using sig_status, not regulation)
oldsedfrap_v_oldsedveh_genes_08feb <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  mutate(sig_status = case_when(
    FDR < 0.05 & logFC > 0 ~ "Up",
    FDR < 0.05 & logFC < 0 ~ "Down",
    TRUE ~ "NS"
  ))

# Label top genes (5 up, 5 down)
top_up <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(sig_status == "Up") %>%
  arrange(desc(logFC)) %>%
  head(50)

print(top_up, n=50)

top_down <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(sig_status == "Down") %>%
  arrange(logFC) %>%
  head(50)

print(top_down, n=50)
##############################################################

###############################################################

###########################################################
#//////Volcano Plot for FRAP///////////////#
##############################################################
# First, create the top_genes dataframe from leading edge analysis
top_genes <- bind_rows(top_up, top_down)

top_genes

# Create the volcano plot with highlighted top genes
frap_volcano_plot <- ggplot(oldsedfrap_v_oldsedveh_genes_08feb, aes(x = logFC, y = -log10(FDR))) +
  
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
    data = top_genes %>% filter(logFC > 0),
    shape = 21,
    size = 3.5,
    stroke = 1.2,
    fill = "#8B0000",      # Dark red
    color = "black"
  ) +
  
  # Highlight downregulated top genes with darker blue and black outline
  geom_point(
    data = top_genes %>% filter(logFC < 0),
    shape = 21,
    size = 3.5,
    stroke = 1.2,
    fill = "#00008B",      # Dark blue
    color = "black"
  ) +
  
  # label selected genes
  geom_text_repel(
    data = top_genes,
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
  
  coord_cartesian(xlim = c(-7, 7), ylim = c(-0.5, 18)) +
  
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Status",
    title = "Volcano plot: Old Sed FRAP vs Old Sed Vehicle",
    subtitle = "Highlighted: Top leading edge genes from pathway enrichment"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title       = element_text(face = "bold")
  )

print(frap_volcano_plot)


ggsave(
  "FigX_frap_volcano_plot.pdf",
  frap_volcano_plot,
  width = 13.5,
  height = 6.6
)

ggsave(
  "FigX_frap_volcano_plot.png",
  frap_volcano_plot,
  width = 13.5,
  height = 6.6,
  dpi = 300
)
############################################################################
#################
#Curated volcano plot genes
#---------------------------#

# Curated gene lists
immune_up_genes <- c("Cd19", "Cd79a", "Ms4a1", "Lck", "Ccl5", "Sell", "S100a8", "S100a9")
metabolic_down_genes <- c("Pck1", "Gck", "Scd1", "Acaa1b", "Fgf21", "Cyp2e1", "Aldh1l1")

# Pull from the full dataset
highlight_up <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(Symbol %in% immune_up_genes)

highlight_down <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(Symbol %in% metabolic_down_genes)

highlight_genes <- bind_rows(highlight_up, highlight_down)

# Volcano plot
frap_volcano_plot <- ggplot(oldsedfrap_v_oldsedveh_genes_08feb, aes(x = logFC, y = -log10(FDR))) +
  
  # background points
  geom_point(aes(color = sig_status),
             alpha = 0.4, size = 1.2) +
  
  # cutoff lines
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  
  # Highlight upregulated immune genes
  geom_point(
    data = highlight_up,
    shape = 21, size = 3.5, stroke = 1.2,
    fill = "#8B0000", color = "black"
  ) +
  
  # Highlight downregulated metabolic genes
  geom_point(
    data = highlight_down,
    shape = 21, size = 3.5, stroke = 1.2,
    fill = "#00008B", color = "black"
  ) +
  
  # label highlighted genes
  geom_text_repel(
    data = highlight_genes,
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
  
  coord_cartesian(xlim = c(-6, 6), ylim = c(-0.05, 10)) +
  
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
    legend.position = "none"
  )

print(frap_volcano_plot)

ggsave(
  "FigX_frap_volcano_plot.pdf",
  frap_volcano_plot,
  width = 10.65,
  height = 6.7
)


##################################
#########################################################################################
create_vec_from_df(oldsedfrap_v_oldsedveh_genes_08feb, "oldsedfrap_v_oldsedveh_genes_08feb")

# 2) GSEA (GO Biological Process)
gseGO_oldsedfrap_vs_oldsedveh.OUTPUT <- gseGO(
  geneList     = oldsedfrap_v_oldsedveh_genes_08feb.vec, 
  ont          = "BP",
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 300,
  pvalueCutoff = 0.5,      # capture everything, filter later
  eps          = 1e-30,      # better estimation for very small p-values
  verbose      = FALSE
)

head(gseGO_oldsedfrap_vs_oldsedveh.OUTPUT@result)


# Simplify, then compute Count/geneRatio again (simplify changes the table)
oldsedfrap_vs_oldsedveh <- clusterProfiler::simplify(gseGO_oldsedfrap_vs_oldsedveh.OUTPUT,
                                                    cutoff = 0.5, by = "p.adjust", select_fun = min)

oldsedfrap_vs_oldsedveh@result

oldsedfrap_vs_oldsedveh_simpl_08Feb.df <- as.data.frame(oldsedfrap_vs_oldsedveh@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

dim(oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
      filter(p.adjust <0.05))

oldsedfrap_vs_oldsedveh_goplot <- oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  scale_size_continuous(range = c(2, 12)) +  # Double the size (default is ~1-6)
  theme_bw()


# Display it
print(oldsedfrap_vs_oldsedveh_goplot)

oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>% 
  dplyr::select(ID, Description, NES, Count)
###########################################################################

# Add pathway categories
oldsedfrap_vs_oldsedveh_simpl_08Feb.df <- oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  mutate(
    Category = case_when(
      # Immune & Inflammation (very broad pattern)
      grepl("immune|lymphocyte|leukocyte|antigen|phagocyt|inflammasome|toll-like|inflammatory|defense response|viral", 
            Description, ignore.case = TRUE) ~ "Immune & Inflammation",
      
      # Metabolism & Mitochondria
      grepl("mitochondrial|oxidative|ATP|electron transport|fatty acid|amino acid|gluconeogenesis|metabolic process|beta-oxidation|peroxisomal|cristae", 
            Description, ignore.case = TRUE) ~ "Metabolism & Mitochondria",
      
      # Cell Cycle & Division
      grepl("chromatid|cell division|stem cell|cell cycle", 
            Description, ignore.case = TRUE) ~ "Cell Cycle & Division",
      
      # Cytoskeleton & Cell Movement
      grepl("actin|chemotaxis|cell shape|lamellipodium|cytoskeleton|migration", 
            Description, ignore.case = TRUE) ~ "Cytoskeleton & Movement",
      
      # Hemostasis & Coagulation
      grepl("hemostasis|coagulation|blood|platelet", 
            Description, ignore.case = TRUE) ~ "Hemostasis & Coagulation",
      
      # Development & Differentiation
      grepl("differentiation|development|morphogenesis|embryonic", 
            Description, ignore.case = TRUE) ~ "Development & Differentiation",
      
      # Other/Cellular Processes
      TRUE ~ "Other Cellular Processes"
    )
  )



# Create bubble plot with enrichment-focused axes (only significant pathways)
oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  
  # bubbles with fixed transparency and black outline
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  
  # reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  
  # size scale
  scale_size_continuous(range = c(3, 15), 
                        name = "Gene Count") +
  
  # fill scheme with 7 distinct colors
  scale_fill_manual(
    values = c(
      "Immune & Inflammation" = "#E41A1C",          # Red
      "Metabolism & Mitochondria" = "#377EB8",      # Blue
      "Cell Cycle & Division" = "#4DAF4A",          # Green
      "Cytoskeleton & Movement" = "#984EA3",        # Purple
      "Hemostasis & Coagulation" = "#FF7F00",       # Orange
      "Development & Differentiation" = "#F781BF",  # Pink
      "Other Cellular Processes" = "#999999"        # Grey
    ),
    name = "Pathway Category"
  ) +
  
  # tighter x-axis limits to remove white space
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 0.5)) +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Gene Ratio (Count / Set Size)",
    title = "GO Pathway Enrichment: Old Sed FRAP vs Old Sed Vehicle",
    subtitle = "Bubble size = gene count (p.adjust < 0.05)"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
############################################################################
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(stringr)
  library(scales)
})

# --- 1) Significant terms only ---
df_sig <- oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  dplyr::filter(p.adjust < 0.05)

# Optional: shorten long GO terms for cleaner labels
df_sig <- df_sig %>%
  dplyr::mutate(
    label = stringr::str_replace_all(Description, "\\s+", " "),
    label = stringr::str_trunc(label, 45)
  )

# --- 2) Choose the dead-zone cutoff (same as before) ---
cut <- 1.5
x_span <- 1.0              # width of each side window
left_lim  <- c(-(cut + x_span), -cut)   # -2.5 to -1.5
right_lim <- c(cut,  cut + x_span)      #  1.5 to  2.5


df_left  <- df_sig %>% dplyr::filter(NES <= -cut)
df_right <- df_sig %>% dplyr::filter(NES >=  cut)

# --- 3) Pick top 5 by Count within each category & side ---
lab_mito <- df_left %>%
  dplyr::filter(Category == "Metabolism & Mitochondria") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

lab_imm <- df_right %>%
  dplyr::filter(Category == "Immune & Inflammation") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

# --- 4) Shared styling ---
fill_vals <- c(
  "Immune & Inflammation" = "#E41A1C",
  "Metabolism & Mitochondria" = "#377EB8",
  "Cell Cycle & Division" = "#4DAF4A",
  "Cytoskeleton & Movement" = "#984EA3",
  "Hemostasis & Coagulation" = "#FF7F00",
  "Development & Differentiation" = "#F781BF",
  "Other Cellular Processes" = "#999999"
)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )


y_scale_shared <- scale_y_continuous(
  limits = c(0, 0.9),
  breaks = seq(0, 0.9, by = 0.1),
  labels = scales::number_format(accuracy = 0.1)
)



# --- 5) Left panel (negative NES) ---
p_left <- ggplot(df_left,
                 aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_mito,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = left_lim,
    breaks = seq(left_lim[1], left_lim[2], by = 0.25)
    ) +
  y_scale_shared +
  labs(
    x = "NES (negative)",
    y = "Gene Ratio (Count / Set Size)"
  ) +
  base_theme +
  theme(legend.position = "none")


# --- 6) Right panel (positive NES) ---
p_right <- ggplot(df_right,
                  aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_imm,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = right_lim,
    breaks = seq(right_lim[1], right_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(
    x = "NES (positive)",
    y = NULL,
    title = "GO Pathway Enrichment: Old Sed FRAP vs Old Sed Vehicle",
    subtitle = paste0(
      "p.adjust < 0.05 | Bubble size = gene count | Center region (",
      -cut, " to ", cut, ") omitted"
    )
  ) +
  base_theme

# --- 7) Combine ---#
p_combined <- p_left + p_right + plot_layout(widths = c(1, 1), guides = "collect")

p_combined

# Save as PNG
ggsave(
  "FigX_frap_go_enrichment.png",
  p_combined,
  width = 13.5,
  height = 6.6,
  dpi = 300
)

# Save as PDF
ggsave(
  "FigX_frap_go_enrichment.pdf",
  p_combined,
  width = 13.5,
  height = 6.6
)
###########################################################################
###########################################################################
# Clean GO bubble plot
#----------------------#
# --- 1) Significant terms only ---
df_sig_frap <- oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  dplyr::filter(p.adjust < 0.05)

df_sig_frap <- df_sig_frap %>%
  dplyr::mutate(
    label = stringr::str_replace_all(Description, "\\s+", " "),
    label = stringr::str_trunc(label, 45)
  )

# --- 2) Choose the dead-zone cutoff ---
cut <- 1.5
x_span <- 1.0
left_lim  <- c(-(cut + x_span), -cut)
right_lim <- c(cut,  cut + x_span)

df_left_frap  <- df_sig_frap %>% dplyr::filter(NES <= -cut)
df_right_frap <- df_sig_frap %>% dplyr::filter(NES >=  cut)

# --- 3) Pick top 5 by Count within each category & side ---
lab_mito_frap <- df_left_frap %>%
  dplyr::filter(Category == "Metabolism & Mitochondria") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

lab_imm_frap <- df_right_frap %>%
  dplyr::filter(Category == "Immune & Inflammation") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

# --- 4) Shared styling ---
fill_vals <- c(
  "Immune & Inflammation" = "#E41A1C",
  "Metabolism & Mitochondria" = "#377EB8",
  "Cell Cycle & Division" = "#4DAF4A",
  "Cytoskeleton & Movement" = "#984EA3",
  "Hemostasis & Coagulation" = "#FF7F00",
  "Development & Differentiation" = "#F781BF",
  "Other Cellular Processes" = "#999999"
)

base_theme_frap <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    legend.position = "none"
  )

y_scale_shared <- scale_y_continuous(
  limits = c(0, 0.9),
  breaks = seq(0, 0.9, by = 0.1),
  labels = scales::number_format(accuracy = 0.1)
)

# --- 5) Left panel (negative NES) ---
p_left_frap <- ggplot(df_left_frap,
                      aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_mito_frap,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = left_lim,
    breaks = seq(left_lim[1], left_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(x = NULL, y = NULL) +
  base_theme_frap

# --- 6) Right panel (positive NES) ---
p_right_frap <- ggplot(df_right_frap,
                       aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_imm_frap,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = right_lim,
    breaks = seq(right_lim[1], right_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL) +
  base_theme_frap

# --- 7) Combine ---
p_combined_frap <- p_left_frap + p_right_frap + plot_layout(widths = c(1, 1))
p_combined_frap

# Save as PDF
ggsave(
  "FigX_frap_go_enrichment.pdf",
  p_combined_frap,
  width = 10.65,
  height = 6.7
)

#
#
###########################################################################
#///////Identifying top core enrichment genes from simplified GO//////////#
###########################################################################
# Filter for significant pathways
sig_pathways <- oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  filter(p.adjust < 0.05)

# Separate into upregulated and downregulated pathways
up_pathways <- sig_pathways %>% filter(NES > 0)
down_pathways <- sig_pathways %>% filter(NES < 0)

# Function to count gene frequency in leading edge
count_leading_edge_genes <- function(pathway_df) {
  # Extract all leading edge genes
  all_genes <- pathway_df %>%
    pull(core_enrichment) %>%
    paste(collapse = "/") %>%
    str_split("/") %>%
    unlist()
  
  # Count frequency of each gene
  gene_counts <- table(all_genes) %>%
    as.data.frame() %>%
    dplyr::arrange(desc(Freq)) %>%
    dplyr::rename(ENTREZID = all_genes, pathway_count = Freq)
  
  return(gene_counts)
}

# Get top genes from upregulated pathways
top_up_genes <- count_leading_edge_genes(up_pathways) %>%
  head(20)

# Get top genes from downregulated pathways
top_down_genes <- count_leading_edge_genes(down_pathways) %>%
  head(20)

# Combine and add direction
top_leading_edge_genes <- bind_rows(
  top_up_genes %>% mutate(direction = "Up"),
  top_down_genes %>% mutate(direction = "Down")
)

# Print just the ENTREZID and pathway count
print("Top 10 genes from UPREGULATED pathways (NES > 0):")
print(top_up_genes)

print("\nTop 10 genes from DOWNREGULATED pathways (NES < 0):")
print(top_down_genes)

# Now look them up in the original dataframe
print("\n=== Details from original dataframe ===")
print("\nUpregulated pathway genes:")
oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(ENTREZID %in% top_up_genes$ENTREZID) %>%
  arrange(desc(FDR)) %>%
  print()

print("\nDownregulated pathway genes:")
oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(ENTREZID %in% top_down_genes$ENTREZID) %>%
  arrange(desc(FDR)) %>%
  print()
#==========================================================#
#===========================================================#

# Label top genes (5 up, 5 down)
oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(sig_status == "Up") %>%
  arrange(desc(logFC)) %>%
  head(10)

oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(sig_status == "Down") %>%
  arrange(logFC) %>%
  head(20)

head(oldsedfrap_vs_oldsedveh_simpl_08Feb.df)

oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  filter(p.adjust < 0.05) %>%
  arrange(NES)
##################################################################
####################################################################

#///////chEA3 analysis//////////#
#================================#

oldsedfrap_upreg_genes <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(sig_status == "Up") %>%
  pull(Symbol)

write.csv(oldsedfrap_upreg_genes, "oldsedfrap_upreg_genes.csv")  

oldsedFRAP_upreg_ENCODE_ChIP_seq <- read_tsv("oldsedFRAP_upreg_ENCODE_ChIP_seq.tsv")

oldsedFRAP_upreg_ENCODE_ChIP_seq

################################################################################
# CIRCULAR BARPLOT - Top 15 ENCODE TFs (OldSed FRAP vs VEH; Up genes)
# Order = Odds Ratio
# Bar height = Odds Ratio
# Color = -log10(FET p-value)
################################################################################
######################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(grid)
})

# ---- 1) Prepare data ----
tf_circ <- oldsedFRAP_upreg_ENCODE_ChIP_seq %>%
  dplyr::mutate(
    fet_p = as.numeric(`FET p-value`),
    odds  = as.numeric(`Odds Ratio`),
    fdr   = as.numeric(FDR)
  ) %>%
  dplyr::filter(!is.na(fet_p), !is.na(odds)) %>%
  
  # Step 1: select top 15 by significance
  dplyr::arrange(fet_p) %>%
  dplyr::slice_head(n = 15) %>%
  
  # Step 2: reorder by bar height
  dplyr::arrange(desc(odds)) %>%
  
  # Step 3: finalize plotting fields
  dplyr::mutate(
    id = dplyr::row_number(),
    tf_label = TF,
    sig = -log10(fet_p)
  )

# ---- 2) Label angles ----
tf_circ <- tf_circ %>%
  dplyr::mutate(
    angle = 90 - 360 * (id - 0.5) / dplyr::n(),
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

# ---- 3) Plot ----
pad <- 0.20

p_tf_circular <- ggplot(tf_circ, aes(x = id, y = odds)) +
  geom_col(aes(fill = sig), color = "black", linewidth = 0.5) +
  
  geom_text(
    aes(
      y = odds + pad,
      label = tf_label,
      angle = angle,
      hjust = hjust
    ),
    size = 3,
    fontface = "bold"
  ) +
  
  scale_fill_gradient(
    low = "#4DBBD5",
    high = "#E64B35",
    name = expression(-log[10]("FET p-value"))
  ) +
  
  coord_polar(start = 0) +
  ylim(-0.3, max(tf_circ$odds) + 0.8) +
  
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    legend.position = "right"
  ) +
  
  labs(
    title = "Top ENCODE TF Enrichment (Old Sed FRAP vs VEH; Upregulated Genes)",
    subtitle = "Top 15 by FET p-value | Order = Odds Ratio | Color = significance"
  )

print(p_tf_circular)

ggsave(
  "FigX_circular_tf_plot.pdf",
  p_tf_circular,
  width = 8,
  height = 8
)

ggsave(
  "FigX_circular_tf_plot.png",
  p_tf_circular,
  width = 8,
  height = 8,
  dpi = 300
)

#################################################################################
################################################################################

##############################################################################
####////CIBERSORT Anaylsis FRAP Only//////########################################
#############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(readxl)
})


library(readxl)
read_clean_counts <- function(file, keep_samples) {
  df <- readxl::read_xlsx(file)
  
  cn <- colnames(df)
  cn[1] <- "Ensembl"
  colnames(df) <- cn
  
  df <- df %>%
    dplyr::select(Ensembl, dplyr::all_of(keep_samples)) %>%
    dplyr::mutate(
      Ensembl_noDec = stringr::str_remove(Ensembl, "\\..+"),
      ENTREZID = AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        keys = Ensembl_noDec,
        keytype = "ENSEMBL",
        column = "ENTREZID",
        multiVals = "first"
      )
    ) %>%
    tidyr::drop_na(ENTREZID)
  
  counts_mat <- df %>%
    dplyr::select(dplyr::all_of(keep_samples)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ replace(., is.na(.), 0))) %>%
    as.matrix()
  
  rownames(counts_mat) <- df$ENTREZID
  
  counts_mat
}

# ---------------------------------------------------------
# 🧫 Define sample groups
# ---------------------------------------------------------
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")

# ---------------------------------------------------------
# 🆕 Read in IRAP and FRAP count data
# ---------------------------------------------------------

# Read both Excel files using your existing helper function
counts_oldfrap <- read_clean_counts(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx",
  keep_samples = c(old_sedfrap_samples, yng_sed_samples)
)

head(counts_oldfrap)

counts_oldveh <- read_clean_counts(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  keep_samples = c(old_sed_samples, yng_sed_samples)
)

head(counts_oldveh)



suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
})

# ---------------------------------------------------------
# 🧹 Keep only OLD_SED_FRAP and OLD_SED_VEH samples
# (drops the Young columns you carried along)
# ---------------------------------------------------------
counts_oldfrap_sub <- counts_oldfrap[, old_sedfrap_samples, drop = FALSE]
counts_oldveh_sub  <- counts_oldveh[,  old_sed_samples,     drop = FALSE]

# sanity checks
stopifnot(all(colnames(counts_oldfrap_sub) == old_sedfrap_samples))
stopifnot(all(colnames(counts_oldveh_sub)  == old_sed_samples))

# ---------------------------------------------------------
# 🔗 Merge the two count matrices
# ---------------------------------------------------------
common_genes <- intersect(rownames(counts_oldfrap_sub), rownames(counts_oldveh_sub))
message("Common genes: ", length(common_genes))

counts_combined <- cbind(
  counts_oldfrap_sub[common_genes, , drop = FALSE],
  counts_oldveh_sub[common_genes,  , drop = FALSE]
)

# ensure numeric/integer-ish
storage.mode(counts_combined) <- "numeric"

# ---------------------------------------------------------
# ⚙️ Normalize (TMM → logCPM)
# ---------------------------------------------------------
dge <- edgeR::DGEList(counts = counts_combined)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logCPM_new <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

dim(logCPM_new)
logCPM_new[1:15, 1:13]
#=============================================#

LM22 <- read.delim("LM22.txt",
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

dim(LM22)
head(LM22[, 1:5])
tail(colnames(LM22))
#---------------------------------------------#

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# map Entrez -> mouse SYMBOL
entrez <- rownames(logCPM_new)

mouse_symbol <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = entrez,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

expr_mouse <- as.data.frame(logCPM_new) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  dplyr::mutate(mouse_symbol = unname(mouse_symbol[ENTREZID])) %>%
  dplyr::filter(!is.na(mouse_symbol) & mouse_symbol != "") %>%
  dplyr::select(-ENTREZID)

# sanity
dim(expr_mouse)
head(expr_mouse$mouse_symbol)
#------------------------------------------------------#

expr_mouse  # contains mouse_symbol + expression

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
})

lm_genes <- toupper(LM22[["Gene symbol"]])

expr_symbol_overlap <- expr_mouse %>%
  dplyr::mutate(human_symbol = toupper(mouse_symbol)) %>%   # KEY STEP
  dplyr::filter(human_symbol %in% lm_genes) %>%
  dplyr::select(-mouse_symbol) %>%
  dplyr::group_by(human_symbol) %>%                         # handle any dupes
  dplyr::summarise(dplyr::across(where(is.numeric), mean), .groups = "drop") %>%
  as.data.frame()

rownames(expr_symbol_overlap) <- expr_symbol_overlap$human_symbol
expr_symbol_overlap$human_symbol <- NULL

dim(expr_symbol_overlap)
head(rownames(expr_symbol_overlap))

mix_out <- expr_symbol_overlap %>%
  tibble::rownames_to_column("Gene symbol")

write.table(
  mix_out,
  file = "CIBERSORTx_mixture_logCPM_OLDSED_FRAP_vs_OLDSED_VEH_LM22overlap.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#Create table for CIBERSORT output
readLines("CIBERSORTx_mixture_logCPM_OLDSED_FRAP_vs_OLDSED_VEH_LM22overlap.txt", n = 5)

any(duplicated(mix_out$`Gene symbol`))

oldfrap_sed_cibersort <- read.csv("CIBERSORTx_Job7_Results.csv")

head(oldfrap_sed_cibersort)


library(dplyr)

ciber_collapsed <- oldfrap_sed_cibersort %>%
  mutate(
    Group = ifelse(Mixture %in% old_sedfrap_samples, "FRAP", "VEH"),
    
    total_B_cells =
      B.cells.naive +
      B.cells.memory +
      Plasma.cells,
    
    total_T_cells =
      T.cells.CD8 +
      T.cells.CD4.naive +
      T.cells.CD4.memory.resting +
      T.cells.CD4.memory.activated +
      T.cells.follicular.helper +
      T.cells.regulatory..Tregs. +
      T.cells.gamma.delta,
    
    total_myeloid =
      Monocytes +
      Macrophages.M0 +
      Macrophages.M1 +
      Macrophages.M2 +
      Dendritic.cells.resting +
      Dendritic.cells.activated,
    
    inflammatory_cells =
      Neutrophils +
      Eosinophils +
      Mast.cells.activated +
      NK.cells.activated
  )

library(ggplot2)

ggplot(ciber_collapsed,
       aes(x = Group, y = total_B_cells, color = Group)) +
  geom_jitter(width = 0.1, size = 2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  theme_bw() +
  labs(
    title = "Total B Cell Fraction in Skeletal Muscle",
    subtitle = "CIBERSORTx (LM22) | Collapsed B-cell compartments",
    y = "Estimated Fraction",
    x = ""
  )


library(tidyr)

plot_long <- ciber_collapsed %>%
  dplyr::select(Mixture, Group,
         total_B_cells,
         total_T_cells,
         total_myeloid,
         inflammatory_cells) %>%
  pivot_longer(
    cols = -c(Mixture, Group),
    names_to = "CellClass",
    values_to = "Fraction"
  )

ggplot(plot_long,
       aes(x = Group, y = Fraction, color = Group)) +
  geom_jitter(width = 0.1, size = 1.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  facet_wrap(~ CellClass, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Collapsed Immune Cell Fractions in Skeletal Muscle",
    subtitle = "CIBERSORTx supports immune composition shifts without acute inflammation",
    x = "",
    y = "Estimated Fraction"
  )

wilcox.test(total_B_cells ~ Group, data = ciber_collapsed)
#---------------------------------------------------------------------#

bcell_stats <- ciber_collapsed %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    median = median(total_B_cells),
    mean = mean(total_B_cells),
    sd = sd(total_B_cells),
    iqr = IQR(total_B_cells)
  )

bcell_stats

wilcox_b <- wilcox.test(total_B_cells ~ Group,
                        data = ciber_collapsed,
                        exact = TRUE)

wilcox_b$p.value


p_bcells <- ggplot(ciber_collapsed,
                   aes(x = Group, y = total_B_cells, fill = Group)) +
  
  # Boxplot (median + IQR)
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.5,
    color = "black"
  ) +
  
  # Individual samples
  geom_jitter(
    aes(color = Group),
    width = 0.08,
    size = 2.8,
    alpha = 0.9
  ) +
  
  # Manual colors (subtle, professional)
  scale_fill_manual(values = c("VEH" = "#4DBBD5", "FRAP" = "red")) +
  scale_color_manual(values = c("VEH" = "#4DBBD5", "FRAP" = "red")) +
  
  # Labels
  labs(
    title = "Total B Cell Fraction in Skeletal Muscle",
    subtitle = "CIBERSORTx (LM22) | Collapsed B-cell compartments",
    y = "Estimated B Cell Fraction",
    x = "",
    caption = paste0("Wilcoxon rank-sum test, p = ",
                     signif(wilcox_b$p.value, 2))
  ) +
  
  # Theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

p_bcells

bcell_prism_long <- ciber_collapsed %>%
  dplyr::select(Mixture, Group, total_B_cells) %>%
  arrange(Group)

write.csv(
  bcell_prism_long,
  file = "Total_B_cells_CIBERSORT_FRAP_vs_VEH_longformat.csv",
  row.names = FALSE
)

bcell_prism_long
##############################################################################
###############################################################################

##############################################################################
####////CIBERSORT Anaylsis IRAP Only//////########################################
#############################################################################
##################################################################################


# ---------------------------------------------------------
# 🧫 Define sample groups
# ---------------------------------------------------------
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")

# ---------------------------------------------------------
# 🆕 Read in IRAP and FRAP count data
# ---------------------------------------------------------

# Read both Excel files using your existing helper function
counts_oldirap <- read_clean_counts(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx",
  keep_samples = c(old_sedirap_samples, yng_sed_samples)
)

head(counts_oldirap)

counts_oldveh <- read_clean_counts(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  keep_samples = c(old_sed_samples, yng_sed_samples)
)

head(counts_oldveh)



suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
})

# ---------------------------------------------------------
# 🧹 Keep only OLD_SED_IRAP and OLD_SED_VEH samples
# (drops the Young columns you carried along)
# ---------------------------------------------------------
counts_oldirap_sub <- counts_oldirap[, old_sedirap_samples, drop = FALSE]
counts_oldveh_sub  <- counts_oldveh[,  old_sed_samples,     drop = FALSE]

# sanity checks
stopifnot(all(colnames(counts_oldirap_sub) == old_sedirap_samples))
stopifnot(all(colnames(counts_oldveh_sub)  == old_sed_samples))

# ---------------------------------------------------------
# 🔗 Merge the two count matrices
# ---------------------------------------------------------
common_genes <- intersect(rownames(counts_oldirap_sub), rownames(counts_oldveh_sub))
message("Common genes: ", length(common_genes))

counts_combined_irap <- cbind(
  counts_oldirap_sub[common_genes, , drop = FALSE],
  counts_oldveh_sub[common_genes,  , drop = FALSE]
)

# ensure numeric/integer-ish
storage.mode(counts_combined_irap) <- "numeric"

# ---------------------------------------------------------
# ⚙️ Normalize (TMM → logCPM)
# ---------------------------------------------------------
dge_irap <- edgeR::DGEList(counts = counts_combined_irap)
dge_irap <- edgeR::calcNormFactors(dge_irap, method = "TMM")
logCPM_irap <- edgeR::cpm(dge_irap, log = TRUE, prior.count = 1)

dim(logCPM_irap)
logCPM_irap[1:15, 1:13]
#=============================================#

LM22 <- read.delim("LM22.txt",
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

dim(LM22)
head(LM22[, 1:5])
tail(colnames(LM22))
#---------------------------------------------#

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# map Entrez -> mouse SYMBOL
entrez <- rownames(logCPM_irap)

mouse_symbol <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = entrez,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

expr_mouse_irap <- as.data.frame(logCPM_irap) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  dplyr::mutate(mouse_symbol = unname(mouse_symbol[ENTREZID])) %>%
  dplyr::filter(!is.na(mouse_symbol) & mouse_symbol != "") %>%
  dplyr::select(-ENTREZID)

# sanity
dim(expr_mouse_irap)
head(expr_mouse_irap$mouse_symbol)
#------------------------------------------------------#

expr_mouse_irap  # contains mouse_symbol + expression

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
})

lm_genes <- toupper(LM22[["Gene symbol"]])

expr_symbol_overlap_irap <- expr_mouse_irap %>%
  mutate(GeneSymbol = toupper(mouse_symbol)) %>%
  filter(GeneSymbol %in% lm_genes) %>%
  dplyr::select(-mouse_symbol) %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean), .groups="drop") %>%
  as.data.frame()

rownames(expr_symbol_overlap_irap) <- expr_symbol_overlap_irap$GeneSymbol
expr_symbol_overlap_irap$GeneSymbol <- NULL


dim(expr_symbol_overlap_irap)
head(rownames(expr_symbol_overlap_irap))

mix_out_irap <- expr_symbol_overlap_irap %>%
  tibble::rownames_to_column("Gene symbol")

write.table(
  mix_out_irap,
  file = "CIBERSORTx_mixture_logCPM_OLDSED_IRAP_vs_OLDSED_VEH_LM22overlap_v2.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

head(mix_out_irap$`Gene symbol`)


#Create table for CIBERSORT output
readLines("CIBERSORTx_mixture_logCPM_OLDSED_IRAP_vs_OLDSED_VEH_LM22overlap.txt", n = 5)

any(duplicated(mix_out_irap$`Gene symbol`))

oldirap_sed_cibersort <- read.csv("CIBERSORTx_Job11_Results.csv")

head(oldirap_sed_cibersort)


library(dplyr)

ciber_collapsed_irap <- oldirap_sed_cibersort %>%
  mutate(
    Group = ifelse(Mixture %in% old_sedirap_samples, "IRAP", "VEH"),
    
    total_B_cells =
      B.cells.naive +
      B.cells.memory +
      Plasma.cells,
    
    total_T_cells =
      T.cells.CD8 +
      T.cells.CD4.naive +
      T.cells.CD4.memory.resting +
      T.cells.CD4.memory.activated +
      T.cells.follicular.helper +
      T.cells.regulatory..Tregs. +
      T.cells.gamma.delta,
    
    total_myeloid =
      Monocytes +
      Macrophages.M0 +
      Macrophages.M1 +
      Macrophages.M2 +
      Dendritic.cells.resting +
      Dendritic.cells.activated,
    
    inflammatory_cells =
      Neutrophils +
      Eosinophils +
      Mast.cells.activated +
      NK.cells.activated
  )

library(ggplot2)

ggplot(ciber_collapsed_irap,
       aes(x = Group, y = total_B_cells, color = Group)) +
  geom_jitter(width = 0.1, size = 2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  theme_bw() +
  labs(
    title = "Total B Cell Fraction in Skeletal Muscle",
    subtitle = "CIBERSORTx (LM22) | Collapsed B-cell compartments",
    y = "Estimated Fraction",
    x = ""
  )


library(tidyr)

plot_long_irap <- ciber_collapsed_irap %>%
  dplyr::select(Mixture, Group,
         total_B_cells,
         total_T_cells,
         total_myeloid,
         inflammatory_cells) %>%
  pivot_longer(
    cols = -c(Mixture, Group),
    names_to = "CellClass",
    values_to = "Fraction"
  )

ggplot(plot_long_irap,
       aes(x = Group, y = Fraction, color = Group)) +
  geom_jitter(width = 0.1, size = 1.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  facet_wrap(~ CellClass, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Collapsed Immune Cell Fractions in Skeletal Muscle",
    subtitle = "CIBERSORTx supports immune composition shifts without acute inflammation",
    x = "",
    y = "Estimated Fraction"
  )

wilcox.test(total_B_cells ~ Group, data = ciber_collapsed_irap)
#---------------------------------------------------------------------#

bcell_stats_irap <- ciber_collapsed_irap %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    median = median(total_B_cells),
    mean = mean(total_B_cells),
    sd = sd(total_B_cells),
    iqr = IQR(total_B_cells)
  )

bcell_stats_irap


wilcox_b_irap <- wilcox.test(total_B_cells ~ Group,
                        data = ciber_collapsed_irap,
                        exact = TRUE)

wilcox_b_irap$p.value


p_bcells_irap <- ggplot(ciber_collapsed_irap,
                   aes(x = Group, y = total_B_cells, fill = Group)) +
  
  # Boxplot (median + IQR)
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.5,
    color = "black"
  ) +
  
  # Individual samples
  geom_jitter(
    aes(color = Group),
    width = 0.08,
    size = 2.8,
    alpha = 0.9
  ) +
  
  # Manual colors (subtle, professional)
  scale_fill_manual(values = c("VEH" = "#4DBBD5", "IRAP" = "#E64B35")) +
  scale_color_manual(values = c("VEH" = "#4DBBD5", "IRAP" = "#E64B35")) +
  
  
  # Labels
  labs(
    title = "Total B Cell Fraction in Skeletal Muscle",
    subtitle = "CIBERSORTx (LM22) | Collapsed B-cell compartments",
    y = "Estimated B Cell Fraction",
    x = "",
    caption = paste0("Wilcoxon rank-sum test, p = ",
                     signif(wilcox_b_irap$p.value, 2))
  ) +
  
  # Theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

p_bcells_irap

bcell_prism_long_irap <- ciber_collapsed_irap %>%
  dplyr::select(Mixture, Group, total_B_cells) %>%
  arrange(Group)

write.csv(
  bcell_prism_long_irap,
  file = "Total_B_cells_CIBERSORT_IRAP_vs_VEH_longformat.csv",
  row.names = FALSE
)

bcell_prism_long_irap
##############################################################################

##////////////Both FRAP and IRAP CIBERSORT Plots//////////////////////#

ciber_collapsed$Group <- factor(
  ciber_collapsed$Group,
  levels = c("VEH", "FRAP")
)

ciber_collapsed_irap$Group <- factor(
  ciber_collapsed_irap$Group,
  levels = c("VEH", "IRAP")
)
#------------------------------------------------------------#
p_bcells <- ggplot(
  ciber_collapsed,
  aes(x = Group, y = total_B_cells, fill = Group)
) +
  
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.5,
    color = "black"
  ) +
  
  geom_jitter(
    aes(color = Group),
    width = 0.08,
    size = 2.8,
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = c("VEH" = "#4DBBD5", "FRAP" = "red")) +
  scale_color_manual(values = c("VEH" = "#4DBBD5", "FRAP" = "red")) +
  
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    title = "Total B Cell Fraction in Skeletal Muscle",
    subtitle = "CIBERSORTx (LM22)",
    y = "Estimated B Cell Fraction",
    x = "",
    caption = paste0(
      "Wilcoxon rank-sum test, p = ",
      signif(wilcox_b$p.value, 2)
    )
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

print(p_bcells)
#---------------------------------------------#

p_bcells_irap <- ggplot(
  ciber_collapsed_irap,
  aes(x = Group, y = total_B_cells, fill = Group)
) +
  
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.5,
    color = "black"
  ) +
  
  geom_jitter(
    aes(color = Group),
    width = 0.08,
    size = 2.8,
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = c("VEH" = "#4DBBD5", "IRAP" = "#E64B35")) +
  scale_color_manual(values = c("VEH" = "#4DBBD5", "IRAP" = "#E64B35")) +
  
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    title = "Total B Cell Fraction in Skeletal Muscle",
    subtitle = "CIBERSORTx (LM22)",
    y = "Estimated B Cell Fraction",
    x = "",
    caption = paste0(
      "Wilcoxon rank-sum test, p = ",
      signif(wilcox_b_irap$p.value, 2)
    )
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

print(p_bcells_irap)
#--------------------------------------------------------#

ggsave(
  "FigX_Bcells_FRAP_vs_VEH.pdf",
  p_bcells,
  width = 8,
  height = 8
)

ggsave(
  "FigX_Bcells_FRAP_vs_VEH.png",
  p_bcells,
  width = 8,
  height = 8,
  dpi = 300
)

ggsave(
  "FigX_Bcells_IRAP_vs_VEH.pdf",
  p_bcells_irap,
  width = 8,
  height = 8
)

ggsave(
  "FigX_Bcells_IRAP_vs_VEH.png",
  p_bcells_irap,
  width = 8,
  height = 8,
  dpi = 300
)
#----------------------------------------------------------------------#





#
#
#
################################################################################
###################################################
#/////Network Analysis///////#
#------------------------------------------#

library(enrichplot)
library(clusterProfiler)

# Filter for significant immune pathways
immune_pathways <- oldsedfrap_vs_oldsedveh
immune_pathways@result <- immune_pathways@result %>%
  filter(p.adjust < 0.05) %>%
  filter(grepl("immune|lymphocyte|leukocyte|antigen|phagocyt|inflammasome|toll-like|inflammatory|defense response|viral", 
               Description, ignore.case = TRUE)) %>%
  arrange(p.adjust) %>%
  head(10) %>%
  # Add required columns for cnetplot
  mutate(
    geneID = core_enrichment,  # cnetplot looks for 'geneID'
    Count = stringr::str_count(core_enrichment, "/") + 1  # Add Count column
  )

library(igraph)
library(ggraph)
library(tidygraph)

# Prepare data for custom network
network_data <- immune_pathways@result %>%
  select(pathway = Description, geneID, p.adjust, NES) %>%
  mutate(genes = strsplit(geneID, "/")) %>%
  tidyr::unnest(genes)

# Add gene info
gene_info <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  select(ENTREZID, Symbol, logFC, FDR)

network_data <- network_data %>%
  left_join(gene_info, by = c("genes" = "ENTREZID")) %>%
  filter(!is.na(Symbol))

# Count gene frequency across pathways
gene_counts <- network_data %>%
  count(Symbol, name = "pathway_count") %>%
  filter(pathway_count >= 3)  # Only genes in 3+ pathways

# Filter to multi-pathway genes only
edges <- network_data %>%
  filter(Symbol %in% gene_counts$Symbol) %>%
  select(from = pathway, to = Symbol, logFC) %>%
  distinct()

# Create nodes
pathway_nodes <- data.frame(
  name = unique(edges$from),
  type = "pathway",
  stringsAsFactors = FALSE
)

gene_nodes <- edges %>%
  select(name = to, logFC) %>%
  distinct() %>%
  left_join(gene_counts, by = c("name" = "Symbol")) %>%
  mutate(type = "gene")

all_nodes <- bind_rows(pathway_nodes, gene_nodes)

# Create graph
g <- graph_from_data_frame(d = edges, vertices = all_nodes, directed = FALSE)

# Plot with ggraph
set.seed(123)
ggraph(g, layout = 'fr') +
  geom_edge_link(alpha = 0.3, color = "grey70", width = 0.8) +
  
  # Pathway nodes (red squares)
  geom_node_point(
    aes(filter = type == "pathway"),
    shape = 22,
    size = 10,
    fill = "#E41A1C",
    color = "black",
    stroke = 1.2
  ) +
  
  # Gene nodes (circles colored by logFC, sized by pathway count)
  geom_node_point(
    aes(filter = type == "gene", fill = logFC, size = pathway_count),
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  
  # Gene labels
  geom_node_text(
    aes(filter = type == "gene", label = name),
    size = 3.5,
    repel = TRUE,
    fontface = "bold"
  ) +
  
  # Pathway labels (wrapped for readability)
  geom_node_text(
    aes(filter = type == "pathway", label = str_wrap(name, 25)),
    size = 2.8,
    repel = TRUE
  ) +
  
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    name = "Gene\nlog2FC"
  ) +
  
  scale_size_continuous(range = c(4, 10), name = "# Pathways") +
  
  theme_graph() +
  labs(
    title = "Immune Pathway-Gene Network",
    subtitle = "Genes appearing in 3+ pathways | FRAP vs VEH"
  )
#####################################################################
#######################################################################
#
#
# New Analysis.......................................
#
#
#
#############################################################################
###################################################################
#####///Analysis of Old SED IRAP vs Old SED VEH////////###########
#-----------------------------------------------------------------#
# 1) Load edgeR *_GENE_* xlsx and build ranked vector
oldsedirap_v_oldsedveh_genes_08feb <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_IRAP-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_SED_IRAP-OLD_SED_VEH_logFC",
  FDR_col   = "OLD_SED_IRAP-OLD_SED_VEH_FDR"
)

#remove duplicate ENTREZID and keep most significant
oldsedirap_v_oldsedveh_genes_08feb <- oldsedirap_v_oldsedveh_genes_08feb %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()

head(oldsedirap_v_oldsedveh_genes_08feb)

oldsedirap_v_oldsedveh_genes_08feb %>%
  filter(FDR < 0.05)

oldsedirap_v_oldsedveh_genes_08feb %>%
  filter(Symbol == "Cxcl10")

oldsedfrap_v_oldsedveh_genes_08feb %>%
  filter(FDR < 0.05)

# Significant genes in each contrast
sig_irap <- oldsedirap_v_oldsedveh_genes_08feb %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

sig_frap <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

# Overlapping ENTREZIDs
overlap_ids <- intersect(sig_irap$ENTREZID, sig_frap$ENTREZID)

length(overlap_ids)
overlap_ids[1:10]

oldsedirap_v_oldsedveh_genes_08feb %>%
  filter (ENTREZID %in% overlap_ids)
#############################################################################

###########################################################
#//////Volcano Plot for IRAP///////////////#
##############################################################

# Create significance and direction column (using sig_status, not regulation)
oldsedirap_v_oldsedveh_genes_08feb <- oldsedirap_v_oldsedveh_genes_08feb %>%
  mutate(sig_status = case_when(
    FDR < 0.05 & logFC > 0 ~ "Up",
    FDR < 0.05 & logFC < 0 ~ "Down",
    TRUE ~ "NS"
  ))

# Label top genes (5 up, 5 down)
top_up_irap <- oldsedirap_v_oldsedveh_genes_08feb %>%
  filter(sig_status == "Up") %>%
  arrange(desc(logFC)) %>%
  head(50)

print(top_up_irap, n=50)

top_down_irap <- oldsedirap_v_oldsedveh_genes_08feb %>%
  filter(sig_status == "Down") %>%
  arrange(logFC) %>%
  head(50)

print(top_down_irap, n=50)


# First, create the top_genes dataframe from leading edge analysis
top_irap_genes <- oldsedirap_v_oldsedveh_genes_08feb %>%
  filter(ENTREZID %in% overlap_ids)

top_irap_genes

# Create the volcano plot with highlighted top genes
irap_volcano_plot<- ggplot(oldsedirap_v_oldsedveh_genes_08feb, aes(x = logFC, y = -log10(FDR))) +
  
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
    data = top_irap_genes %>% filter(logFC > 0),
    shape = 21,
    size = 3.5,
    stroke = 1.2,
    fill = "#8B0000",      # Dark red
    color = "black"
  ) +
  
  # Highlight downregulated top genes with darker blue and black outline
  geom_point(
    data = top_irap_genes %>% filter(logFC < 0),
    shape = 21,
    size = 3.5,
    stroke = 1.2,
    fill = "#00008B",      # Dark blue
    color = "black"
  ) +
  
  # label selected genes
  geom_text_repel(
    data = top_irap_genes,
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
  
  coord_cartesian(xlim = c(-7, 7), ylim = c(-0.5, 18)) +   # <-- add here
  
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Status",
    title = "Volcano plot: Old Sed IRAP vs Old Sed Vehicle",
    subtitle = "Highlighted: Top leading edge genes from pathway enrichment"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title       = element_text(face = "bold")
  )


print(irap_volcano_plot)


ggsave(
  "FigX_irap_volcano_plot.pdf",
  irap_volcano_plot,
  width = 13.5,
  height = 6.6
)

ggsave(
  "FigX_irap_volcano_plot.png",
  irap_volcano_plot,
  width = 8,
  height = 8,
  dpi = 300
)
########################################################################
###########################################################################
#Curated IRAP volcano plot
#--------------------------#

# Make sure sig_status exists for IRAP dataset
oldsedirap_v_oldsedveh_genes_08feb <- oldsedirap_v_oldsedveh_genes_08feb %>%
  mutate(sig_status = case_when(
    FDR < 0.05 & logFC > 0 ~ "Up",
    FDR < 0.05 & logFC < 0 ~ "Down",
    TRUE ~ "NS"
  ))

# Curated gene lists for IRAP
immune_up_genes_irap <- c("Ighd", "Fcmr", "Cd22", "Pax5", "Sell", "Ltb", "Saa3")
down_genes_irap <- c("Vegfa", "Retsat")

highlight_up <- oldsedirap_v_oldsedveh_genes_08feb %>%
  filter(Symbol %in% immune_up_genes_irap)
highlight_down <- oldsedirap_v_oldsedveh_genes_08feb %>%
  filter(Symbol %in% down_genes_irap)
highlight_genes <- bind_rows(highlight_up, highlight_down)

# Volcano plot
irap_volcano_plot <- ggplot(oldsedirap_v_oldsedveh_genes_08feb, aes(x = logFC, y = -log10(FDR))) +
  
  geom_point(aes(color = sig_status),
             alpha = 0.4, size = 1.2) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  
  geom_point(
    data = highlight_up,
    shape = 21, size = 3.5, stroke = 1.2,
    fill = "#8B0000", color = "black"
  ) +
  
  geom_point(
    data = highlight_down,
    shape = 21, size = 3.5, stroke = 1.2,
    fill = "#00008B", color = "black"
  ) +
  
  geom_text_repel(
    data = highlight_genes,
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
  
  coord_cartesian(xlim = c(-6, 6), ylim = c(-0.05, 10)) +
  
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
    legend.position = "none"
  )

print(irap_volcano_plot)
ggsave(
  "FigX_irap_volcano_plot.pdf",
  irap_volcano_plot,
  width = 10.65,
  height = 6.7
)





#########################################################################################
create_vec_from_df(oldsedirap_v_oldsedveh_genes_08feb, "oldsedirap_v_oldsedveh_genes_08feb")

# 2) GSEA (GO Biological Process)
gseGO_oldsedirap_vs_oldsedveh.OUTPUT <- gseGO(
  geneList     = oldsedirap_v_oldsedveh_genes_08feb.vec, 
  ont          = "BP",
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 300,
  pvalueCutoff = 0.5,      # capture everything, filter later
  eps          = 1e-30,      # better estimation for very small p-values
  verbose      = FALSE
)

head(gseGO_oldsedirap_vs_oldsedveh.OUTPUT@result)


# Simplify, then compute Count/geneRatio again (simplify changes the table)
oldsedirap_vs_oldsedveh <- clusterProfiler::simplify(gseGO_oldsedirap_vs_oldsedveh.OUTPUT,
                                                     cutoff = 0.5, by = "p.adjust", select_fun = min)

oldsedirap_vs_oldsedveh@result

oldsedirap_vs_oldsedveh_simpl_08Feb.df <- as.data.frame(oldsedirap_vs_oldsedveh@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

dim(oldsedirap_vs_oldsedveh_simpl_08Feb.df %>%
      filter(p.adjust <0.05))

head(oldsedirap_vs_oldsedveh_simpl_08Feb.df)

oldsedirap_vs_oldsedveh_goplot <- oldsedirap_vs_oldsedveh_simpl_08Feb.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  scale_size_continuous(range = c(2, 12)) +  # Double the size (default is ~1-6)
  theme_bw()


# Display it
print(oldsedirap_vs_oldsedveh_goplot)
###############################################################################
###############################################################################


# Add pathway categories
oldsedirap_vs_oldsedveh_simpl_08Feb.df <- oldsedirap_vs_oldsedveh_simpl_08Feb.df %>%
  mutate(
    Category = case_when(
      # Immune & Inflammation (very broad pattern)
      grepl("immune|lymphocyte|leukocyte|antigen|phagocyt|inflammasome|toll-like|inflammatory|defense response|viral", 
            Description, ignore.case = TRUE) ~ "Immune & Inflammation",
      
      # Metabolism & Mitochondria
      grepl("mitochondrial|oxidative|ATP|electron transport|fatty acid|amino acid|gluconeogenesis|metabolic process|beta-oxidation|peroxisomal|cristae", 
            Description, ignore.case = TRUE) ~ "Metabolism & Mitochondria",
      
      # Cell Cycle & Division
      grepl("chromatid|cell division|stem cell|cell cycle", 
            Description, ignore.case = TRUE) ~ "Cell Cycle & Division",
      
      # Cytoskeleton & Cell Movement
      grepl("actin|chemotaxis|cell shape|lamellipodium|cytoskeleton|migration", 
            Description, ignore.case = TRUE) ~ "Cytoskeleton & Movement",
      
      # Hemostasis & Coagulation
      grepl("hemostasis|coagulation|blood|platelet", 
            Description, ignore.case = TRUE) ~ "Hemostasis & Coagulation",
      
      # Development & Differentiation
      grepl("differentiation|development|morphogenesis|embryonic", 
            Description, ignore.case = TRUE) ~ "Development & Differentiation",
      
      # Other/Cellular Processes
      TRUE ~ "Other Cellular Processes"
    )
  )



# Create bubble plot with enrichment-focused axes (only significant pathways)
oldsedirap_vs_oldsedveh_simpl_08Feb.df %>%
  filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  
  # bubbles with fixed transparency and black outline
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  
  # reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  
  # size scale
  scale_size_continuous(range = c(3, 15), 
                        name = "Gene Count") +
  
  # fill scheme with 7 distinct colors
  scale_fill_manual(
    values = c(
      "Immune & Inflammation" = "#E41A1C",          # Red
      "Metabolism & Mitochondria" = "#377EB8",      # Blue
      "Cell Cycle & Division" = "#4DAF4A",          # Green
      "Cytoskeleton & Movement" = "#984EA3",        # Purple
      "Hemostasis & Coagulation" = "#FF7F00",       # Orange
      "Development & Differentiation" = "#F781BF",  # Pink
      "Other Cellular Processes" = "#999999"        # Grey
    ),
    name = "Pathway Category"
  ) +
  
  # tighter x-axis limits to remove white space
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 0.5)) +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Gene Ratio (Count / Set Size)",
    title = "GO Pathway Enrichment: Old Sed IRAP vs Old Sed Vehicle",
    subtitle = "Bubble size = gene count (p.adjust < 0.05)"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
############################################################################
###########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(stringr)
  library(scales)
})

# --- 1) Significant terms only ---
df_sig_irap <- oldsedirap_vs_oldsedveh_simpl_08Feb.df %>%
  dplyr::filter(p.adjust < 0.05)

# Optional: shorten long GO terms for cleaner labels
df_sig_irap <- df_sig_irap %>%
  dplyr::mutate(
    label = stringr::str_replace_all(Description, "\\s+", " "),
    label = stringr::str_trunc(label, 45)
  )

# --- 2) Choose the dead-zone cutoff (same as before) ---
cut <- 1.5
x_span <- 1.0              # width of each side window
left_lim  <- c(-(cut + x_span), -cut)   # -2.5 to -1.5
right_lim <- c(cut,  cut + x_span)      #  1.5 to  2.5


df_left_irap  <- df_sig_irap %>% dplyr::filter(NES <= -cut)
df_right_irap <- df_sig_irap %>% dplyr::filter(NES >=  cut)

# --- 3) Pick top 5 by Count within each category & side ---
lab_mito_irap <- df_left_irap %>%
  dplyr::filter(Category == "Metabolism & Mitochondria") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

lab_imm_irap <- df_right_irap %>%
  dplyr::filter(Category == "Immune & Inflammation") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

# --- 4) Shared styling ---
fill_vals <- c(
  "Immune & Inflammation" = "#E41A1C",
  "Metabolism & Mitochondria" = "#377EB8",
  "Cell Cycle & Division" = "#4DAF4A",
  "Cytoskeleton & Movement" = "#984EA3",
  "Hemostasis & Coagulation" = "#FF7F00",
  "Development & Differentiation" = "#F781BF",
  "Other Cellular Processes" = "#999999"
)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )


y_scale_shared <- scale_y_continuous(
  limits = c(0, 0.9),
  breaks = seq(0, 0.9, by = 0.1),
  labels = scales::number_format(accuracy = 0.1)
)



# --- 5) Left panel (negative NES) ---
# --- 5) Left panel (negative NES) ---
p_left_irap <- ggplot(df_left_irap,
                      aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_mito_irap,   # FIXED
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = left_lim,
    breaks = seq(left_lim[1], left_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(
    x = "NES (negative)",
    y = "Gene Ratio (Count / Set Size)"
  ) +
  base_theme +
  theme(legend.position = "none")


# --- 6) Right panel (positive NES) ---
p_right_irap <- ggplot(df_right_irap,
                       aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_imm_irap,    # FIXED
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = right_lim,
    breaks = seq(right_lim[1], right_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(
    x = "NES (positive)",
    y = NULL,
    title = "GO Pathway Enrichment: Old Sed IRAP vs Old Sed Vehicle",
    subtitle = paste0(
      "p.adjust < 0.05 | Bubble size = gene count | Center region (",
      -cut, " to ", cut, ") omitted"
    )
  ) +
  base_theme


# --- 7) Combine ---
p_combined_irap <- p_left_irap + p_right_irap + plot_layout(widths = c(1, 1), guides = "collect")

print(p_combined_irap)

# Save as PNG
ggsave(
  "FigX_irap_go_enrichment.png",
  p_combined_irap,
  width = 13.5,
  height = 6.6,
  dpi = 300
)

# Save as PDF
ggsave(
  "FigX_irap_go_enrichment.pdf",
  p_combined_irap,
  width = 13.5,
  height = 6.6
)

###############################################################################
###############################################################################
#
# Clean irap GO Bubble Plot
#--------------------------#

# --- 1) Significant terms only ---
df_sig_irap <- oldsedirap_vs_oldsedveh_simpl_08Feb.df %>%
  dplyr::filter(p.adjust < 0.05)

df_sig_irap

df_sig_irap <- df_sig_irap %>%
  dplyr::mutate(
    label = stringr::str_replace_all(Description, "\\s+", " "),
    label = stringr::str_trunc(label, 45)
  )

# --- 2) Choose the dead-zone cutoff ---
cut <- 1.5
x_span <- 1.0
left_lim  <- c(-(cut + x_span), -cut)
right_lim <- c(cut,  cut + x_span)

df_left_irap  <- df_sig_irap %>% dplyr::filter(NES <= -cut)
df_right_irap <- df_sig_irap %>% dplyr::filter(NES >=  cut)

# --- 3) Pick top 5 by Count within each category & side ---
lab_mito_irap <- df_left_irap %>%
  dplyr::filter(Category == "Metabolism & Mitochondria") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

lab_imm_irap <- df_right_irap %>%
  dplyr::filter(Category == "Immune & Inflammation") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

# --- 4) Shared styling ---
fill_vals <- c(
  "Immune & Inflammation" = "#E41A1C",
  "Metabolism & Mitochondria" = "#377EB8",
  "Cell Cycle & Division" = "#4DAF4A",
  "Cytoskeleton & Movement" = "#984EA3",
  "Hemostasis & Coagulation" = "#FF7F00",
  "Development & Differentiation" = "#F781BF",
  "Other Cellular Processes" = "#999999"
)

base_theme_irap <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    legend.position = "none"
  )

y_scale_shared <- scale_y_continuous(
  limits = c(0, 0.9),
  breaks = seq(0, 0.9, by = 0.1),
  labels = scales::number_format(accuracy = 0.1)
)

# --- 5) Left panel (negative NES) ---
p_left_irap <- ggplot(df_left_irap,
                      aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_mito_irap,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = left_lim,
    breaks = seq(left_lim[1], left_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(x = NULL, y = NULL) +
  base_theme_irap

# --- 6) Right panel (positive NES) ---
p_right_irap <- ggplot(df_right_irap,
                       aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_imm_irap,
    aes(label = label),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category") +
  scale_x_continuous(
    limits = right_lim,
    breaks = seq(right_lim[1], right_lim[2], by = 0.25)
  ) +
  y_scale_shared +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL) +
  base_theme_irap

# --- 7) Combine ---
p_combined_irap <- p_left_irap + p_right_irap + plot_layout(widths = c(1, 1))
p_combined_irap

# Save as PDF
ggsave(
  "FigX_irap_go_enrichment.pdf",
  p_combined_irap,
  width = 10.65,
  height = 6.7
)


############################################################################
##############################################################################

####////CIBERSOrT Analysis old only//////###############
#-------------------------------------------------------#

########################################
# Libraries
########################################
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

########################################
# Sample groups
########################################
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")

# ✅ OLD-only samples for normalization + CIBERSORT mixture
old_only_samples <- c(
  old_sed_samples,
  old_sedfrap_samples,
  old_sedirap_samples
)

########################################
# Helper: read counts and map Ensembl->Entrez
########################################
read_clean_counts <- function(file, keep_samples) {
  df <- readxl::read_xlsx(file)
  cn <- colnames(df)
  cn[1] <- "Ensembl"
  colnames(df) <- cn
  
  keep_present <- intersect(keep_samples, colnames(df))
  if (length(keep_present) == 0) {
    stop("None of the requested samples were found in: ", file)
  }
  
  df <- df %>%
    dplyr::select(Ensembl, dplyr::all_of(keep_present)) %>%
    dplyr::mutate(
      Ensembl_noDec = stringr::str_remove(Ensembl, "\\..+"),
      ENTREZID = AnnotationDbi::mapIds(
        org.Mm.eg.db,
        keys = Ensembl_noDec,
        keytype = "ENSEMBL",
        column = "ENTREZID",
        multiVals = "first"
      )
    ) %>%
    tidyr::drop_na(ENTREZID)
  
  counts_mat <- df %>%
    dplyr::select(dplyr::all_of(keep_present)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~replace(., is.na(.), 0))) %>%
    as.matrix()
  
  rownames(counts_mat) <- df$ENTREZID
  storage.mode(counts_mat) <- "numeric"
  counts_mat
}

########################################
# Count files
########################################
count_files <- c(
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx",
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx"
)

########################################
# 2) Read each file; merge by intersecting genes (OLD-only)
########################################
counts_list <- lapply(count_files, read_clean_counts, keep_samples = old_only_samples)

common_genes <- Reduce(intersect, lapply(counts_list, rownames))
message("Common genes across files: ", length(common_genes))

counts_all <- do.call(
  cbind,
  lapply(counts_list, function(m) m[common_genes, , drop = FALSE])
)

# If a sample appears in multiple files, keep the first occurrence
counts_all <- counts_all[, !duplicated(colnames(counts_all)), drop = FALSE]

# ✅ Reorder + restrict to OLD-only desired order
counts_all <- counts_all[, intersect(old_only_samples, colnames(counts_all)), drop = FALSE]

# Sanity checks
stopifnot(all(colnames(counts_all) %in% old_only_samples))
message("Final matrix dims (genes x samples): ", paste(dim(counts_all), collapse=" x "))

########################################
# 3) TMM -> logCPM (OLD-only)
########################################
dge <- edgeR::DGEList(counts = counts_all)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logCPM_old_only <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

########################################
# 4) Load LM22
########################################
LM22 <- read.delim("LM22.txt", header = TRUE, sep = "\t", check.names = FALSE)
lm_genes <- toupper(LM22[["Gene symbol"]])

########################################
# 5) Entrez -> mouse SYMBOL
########################################
entrez <- rownames(logCPM_old_only)

mouse_symbol <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = entrez,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

expr_mouse <- as.data.frame(logCPM_old_only) %>%
  tibble::rownames_to_column("ENTREZID") %>%
  dplyr::mutate(mouse_symbol = unname(mouse_symbol[ENTREZID])) %>%
  dplyr::filter(!is.na(mouse_symbol) & mouse_symbol != "") %>%
  dplyr::select(-ENTREZID)

########################################
# 6) Keep only genes overlapping LM22 (uppercase match)
########################################
expr_symbol_overlap <- expr_mouse %>%
  dplyr::mutate(GeneSymbol = toupper(mouse_symbol)) %>%
  dplyr::filter(GeneSymbol %in% lm_genes) %>%
  dplyr::select(-mouse_symbol) %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), mean), .groups = "drop") %>%
  as.data.frame()

rownames(expr_symbol_overlap) <- expr_symbol_overlap$GeneSymbol
expr_symbol_overlap$GeneSymbol <- NULL

message("LM22-overlap genes x samples: ", paste(dim(expr_symbol_overlap), collapse=" x "))

########################################
# 7) Write mixture file for CIBERSORTx
########################################
mix_out <- expr_symbol_overlap %>%
  tibble::rownames_to_column("Gene symbol")

out_tsv <- "CIBERSORTx_mixture_logCPM_OLDONLY_VEH_FRAP_IRAP_LM22overlap.tsv"
write.table(
  mix_out,
  file = out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("✅ Wrote: ", out_tsv)

# quick preview
readLines(out_tsv, n = 2)
any(duplicated(mix_out$`Gene symbol`))
#--------------------------------------------------#

########################################
# Libraries
########################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

########################################
# Read CIBERSORTx results
########################################
cib <- read.csv("CIBERSORTx_Job9_Results.csv")  # or put exact file name

# Optional sanity: make sure Mixture IDs match your sample IDs
 sort(unique(cib$Mixture))

########################################
# Assign groups (3 groups)
########################################
ciber_collapsed <- cib %>%
  mutate(
    Group = dplyr::case_when(
      Mixture %in% old_sed_samples     ~ "VEH",
      Mixture %in% old_sedfrap_samples ~ "FRAP",
      Mixture %in% old_sedirap_samples ~ "IRAP",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Group = factor(Group, levels = c("VEH", "IRAP", "FRAP")),
    
    # Keep naive B separate (as requested)
    B_naive = B.cells.naive,
    
    # B-lineage (still useful as a second view)
    total_B_lineage = B.cells.naive + B.cells.memory + Plasma.cells,
    
    total_T_cells =
      T.cells.CD8 +
      T.cells.CD4.naive +
      T.cells.CD4.memory.resting +
      T.cells.CD4.memory.activated +
      T.cells.follicular.helper +
      T.cells.regulatory..Tregs. +
      T.cells.gamma.delta,
    
    total_myeloid =
      Monocytes +
      Macrophages.M0 +
      Macrophages.M1 +
      Macrophages.M2 +
      Dendritic.cells.resting +
      Dendritic.cells.activated,
    
    inflammatory_cells =
      Neutrophils +
      Eosinophils +
      Mast.cells.activated +
      NK.cells.activated
  )

########################################
# Summary table
########################################
 total_B_lineage_stats <- ciber_collapsed %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean = mean(total_B_lineage),
    sd = sd(total_B_lineage),
    median = median(total_B_lineage),
    iqr = IQR(total_B_lineage),
    .groups = "drop"
  )

 total_B_lineage_stats

########################################
# Plot: Naive B cells (3 groups)
########################################
group_cols <- c("VEH" = "#4DBBD5", "IRAP" = "#00A087", "FRAP" = "#E64B35")

p_total_B_lineage <- ggplot(ciber_collapsed, aes(x = Group, y = total_B_lineage, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.45, color = "black") +
  geom_jitter(aes(color = Group), width = 0.10, size = 2.8, alpha = 0.9) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "total_B_lineage fraction in skeletal muscle",
    subtitle = "CIBERSORTx (LM22) | 3 groups: VEH vs IRAP vs FRAP",
    y = "Estimated fraction",
    x = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

p_total_B_lineage

########################################
# One-way ANOVA + Tukey (no emmeans needed)
########################################
fit_aov <- aov(total_B_lineage ~ Group, data = ciber_collapsed)
summary(fit_aov)

tuk <- TukeyHSD(fit_aov)
tuk

# If you specifically want "each vs VEH" extracted:
tuk_vs_veh <- as.data.frame(tuk$Group) %>%
  tibble::rownames_to_column("contrast") %>%
  filter(grepl("VEH", contrast))  # keeps VEH-IRAP and VEH-FRAP (direction depends)

tuk_vs_veh

########################################
# Prism tables
########################################

# (A) Long format (best for Prism grouped scatter/box)
total_B_lineage_prism_long <- ciber_collapsed %>%
  select(Mixture, Group, total_B_lineage) %>%
  arrange(Group)

write.csv(
  total_B_lineage_prism_long,
  "total_B_lineage_CIBERSORT_VEH_IRAP_FRAP_long.csv",
  row.names = FALSE
)

# (B) Wide format (one column per group; Prism likes this too)
total_B_lineage_prism_wide <- ciber_collapsed %>%
  select(Mixture, Group, total_B_lineage) %>%
  mutate(row_in_group = ave(total_B_lineage, Group, FUN = seq_along)) %>%
  select(-Mixture) %>%
  pivot_wider(names_from = Group, values_from = total_B_lineage) %>%
  arrange(row_in_group) %>%
  select(-row_in_group)

write.csv(
  total_B_lineage_prism_wide,
  "total_B_lineage_CIBERSORT_VEH_IRAP_FRAP_wide.csv",
  row.names = FALSE
)

total_B_lineage_prism_long
#===========================================================#
###################################################################
#######################################################################

############################################################
# FIGURE — Heatmap of Core Enrichment genes (both directions)
# FRAP vs Old Sed Vehicle | All 3 old sed groups as columns
############################################################

# 1) Core enrichment ENTREZ IDs (both directions, p.adjust < 0.05)
frap_core_entrez_vec <- oldsedfrap_vs_oldsedveh_simpl_08Feb.df %>%
  dplyr::filter(p.adjust < 0.05) %>%          # both NES > 0 and NES < 0
  dplyr::pull(core_enrichment) %>%
  strsplit("/") %>% unlist() %>% unique() %>% sort()

length(frap_core_entrez_vec)
head(frap_core_entrez_vec)

# 2) Load and prep all three count files
load_counts <- function(filepath) {
  df <- readxl::read_xlsx(filepath)
  colnames(df)[1] <- "Ensembl"
  df <- df %>%
    dplyr::mutate(Ensembl_noDec = stringr::str_remove(Ensembl, "\\..+")) %>%
    dplyr::select(-Ensembl)
  return(df)
}

counts_veh  <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
counts_frap <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx")
counts_irap <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx")

# 3) Define samples
old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")
yng_sed_samples     <- c("T_05","T_17","T_26","T_32","T_35","T_38")

ordered_columns <- c(old_sed_samples, old_sedfrap_samples, old_sedirap_samples)

colnames(counts_veh)

# 4) Merge the three count files on Ensembl_noDec
# Each file contributes only its relevant sample columns to avoid duplication
counts_merged <- counts_veh %>%
  dplyr::select(Ensembl_noDec, dplyr::all_of(old_sed_samples)) %>%
  dplyr::inner_join(
    counts_frap %>% dplyr::select(Ensembl_noDec, dplyr::all_of(old_sedfrap_samples)),
    by = "Ensembl_noDec"
  ) %>%
  dplyr::inner_join(
    counts_irap %>% dplyr::select(Ensembl_noDec, dplyr::all_of(old_sedirap_samples)),
    by = "Ensembl_noDec"
  )



# 5) Map Ensembl → ENTREZID
counts_merged$ENTREZID <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys     = counts_merged$Ensembl_noDec,
  keytype  = "ENSEMBL",
  column   = "ENTREZID",
  multiVals = "first"
)

counts_merged <- tidyr::drop_na(counts_merged, ENTREZID)

# 6) Build counts matrix (exclude young samples — already absent, but explicit)
counts_matrix <- counts_merged %>%
  dplyr::select(dplyr::all_of(ordered_columns)) %>%
  as.matrix()
rownames(counts_matrix) <- as.character(counts_merged$ENTREZID)

# 7) edgeR normalization → logCPM
group <- factor(c(
  rep("OLD_SED",      length(old_sed_samples)),
  rep("OLD_SED_FRAP", length(old_sedfrap_samples)),
  rep("OLD_SED_IRAP", length(old_sedirap_samples))
))

dge <- edgeR::DGEList(counts = counts_matrix, group = group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logCPM_matrix <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# 8) Filter logCPM to core enrichment genes present in matrix
core_ids_chr <- as.character(frap_core_entrez_vec)
core_ids_present <- core_ids_chr[core_ids_chr %in% rownames(logCPM_matrix)]

frap_core_heatmap_mat <- logCPM_matrix[core_ids_present, ordered_columns, drop = FALSE]

# 9) Column annotation bar
annotation_col <- data.frame(
  Group = factor(c(
    rep("Old Sed Veh",  length(old_sed_samples)),
    rep("Old Sed FRAP", length(old_sedfrap_samples)),
    rep("Old Sed IRAP", length(old_sedirap_samples))
  ))
)
rownames(annotation_col) <- ordered_columns

annotation_colors <- list(
  Group = c(
    "Old Sed Veh"  = "#999999",
    "Old Sed FRAP" = "#E69F00",
    "Old Sed IRAP" = "#56B4E9"
  )
)

# 10) Plot heatmap (clustered rows, both directions)
frap_core_heatplot <- pheatmap::pheatmap(
  frap_core_heatmap_mat,
  cluster_rows    = TRUE,
  cluster_cols    = FALSE,          # keep group order intact
  scale           = "row",
  show_rownames   = FALSE,
  show_colnames   = TRUE,
  color           = colorRampPalette(c("blue","white","red"))(50),
  annotation_col  = annotation_col,
  annotation_colors = annotation_colors,
  main            = "Core Enrichment Genes: Old Sed FRAP vs Old Sed Vehicle"
)

# 11) Save
ggsave(
  "FigX_frap_core_enrichment_heatmap.png",
  plot   = frap_core_heatplot,
  width  = 10,
  height = 10,
  dpi    = 300
)

ggsave(
  "FigX_frap_core_enrichment_heatmap.pdf",
  plot   = frap_core_heatplot,
  width  = 10,
  height = 10
)

logCPM_matrix
logCPM_matrix["15945", ]

library(ggplot2)
library(dplyr)
library(ggpubr)

# 1) Extract Cxcl10 logCPM and build tidy data frame
cxcl10_df <- data.frame(
  logCPM = logCPM_matrix["15945", ],
  Sample = colnames(logCPM_matrix[, c(old_sed_samples, old_sedfrap_samples, old_sedirap_samples)])
) %>%
  dplyr::mutate(
    Group = dplyr::case_when(
      Sample %in% old_sed_samples     ~ "Old Sed Veh",
      Sample %in% old_sedfrap_samples ~ "Old Sed FRAP",
      Sample %in% old_sedirap_samples ~ "Old Sed IRAP"
    ),
    Group = factor(Group, levels = c("Old Sed Veh", "Old Sed IRAP", "Old Sed FRAP"))
  )

# 2) Summary stats for bar + error bars (mean ± SEM)
cxcl10_summary <- cxcl10_df %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    mean_logCPM = mean(logCPM),
    sem         = sd(logCPM) / sqrt(n()),
    .groups     = "drop"
  )

# 3) One-way ANOVA
cxcl10_anova <- aov(logCPM ~ Group, data = cxcl10_df)
summary(cxcl10_anova)

# Optional: Tukey post-hoc for pairwise comparisons
TukeyHSD(cxcl10_anova)

# 4) Plot
cxcl10_plot <- ggplot(cxcl10_summary, aes(x = Group, y = mean_logCPM, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = mean_logCPM - sem, ymax = mean_logCPM + sem),
    width = 0.2, linewidth = 0.6
  ) +
  geom_jitter(
    data = cxcl10_df,
    aes(y = logCPM),
    width = 0.1, size = 2, alpha = 0.7, color = "black"
  ) +
  # Add ANOVA p-value via ggpubr
  stat_compare_means(
    data        = cxcl10_df,
    aes(x = Group, y = logCPM),
    method      = "anova",
    label       = "p.format",
    label.x     = 1.5,
    label.y     = max(cxcl10_df$logCPM) * 1.15
  ) +
  scale_fill_manual(values = c(
    "Old Sed Veh"  = "#999999",
    "Old Sed IRAP" = "#56B4E9",
    "Old Sed FRAP" = "#C0392B"    # deep crimson, less harsh than pure red
  )) +
  labs(
    x     = NULL,
    y     = "logCPM",
    title = "Cxcl10 Expression",
    subtitle = "Mean ± SEM | One-way ANOVA"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(face = "bold")
  )

print(cxcl10_plot)

# 5) Save
ggsave("Cxcl10_logCPM_barplot.png", cxcl10_plot, width = 8, height = 8, dpi = 300)
ggsave("Cxcl10_logCPM_barplot.pdf", cxcl10_plot, width = 8, height = 8)
#########################################################################
##########################################################################

#
#
#
#
##########################################################################

############################################################
# Heatmap + Aging-Directed Barplot — 93 validated aging genes
# Old Sed FRAP vs Old Sed VEH
############################################################

# -------------------------------------------------------
# 0) Load and merge counts from two files
# -------------------------------------------------------
counts_veh  <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
counts_frap <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_FRAP-YNG_SED_VEH.xlsx")

counts_veh

old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedfrap_samples <- c("T_54","T_55","T_56","T_57","T_58")

ordered_columns_frap <- c(old_sed_samples, old_sedfrap_samples)

counts_merged_frap <- counts_veh %>%
  dplyr::select(Ensembl_noDec, dplyr::all_of(old_sed_samples)) %>%
  dplyr::inner_join(
    counts_frap %>% dplyr::select(Ensembl_noDec, dplyr::all_of(old_sedfrap_samples)),
    by = "Ensembl_noDec"
  )

# Map Ensembl -> ENTREZID
counts_merged_frap$ENTREZID <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys      = counts_merged_frap$Ensembl_noDec,
  keytype   = "ENSEMBL",
  column    = "ENTREZID",
  multiVals = "first"
)
counts_merged_frap <- tidyr::drop_na(counts_merged_frap, ENTREZID)

# -------------------------------------------------------
# 1) edgeR normalization -> logCPM
# -------------------------------------------------------
counts_matrix_frap <- counts_merged_frap %>%
  dplyr::select(dplyr::all_of(ordered_columns_frap)) %>%
  as.matrix()
rownames(counts_matrix_frap) <- as.character(counts_merged_frap$ENTREZID)

group_frap <- factor(c(
  rep("OLD_SED_VEH",  length(old_sed_samples)),
  rep("OLD_SED_FRAP", length(old_sedfrap_samples))
))

dge_frap <- edgeR::DGEList(counts = counts_matrix_frap, group = group_frap)
dge_frap <- edgeR::calcNormFactors(dge_frap, method = "TMM")
logCPM_matrix_frap <- edgeR::cpm(dge_frap, log = TRUE, prior.count = 1)

# -------------------------------------------------------
# 2) Compute aging-directed logFC FIRST (to set row order)
# -------------------------------------------------------

# Get symbols for the 93 genes present in the matrix
core_ids_present_frap <- core_entrez[core_entrez %in% rownames(logCPM_matrix_frap)]
entrez_to_symbol <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys    = core_ids_present_frap,
  column  = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
core_symbols_present <- unname(na.omit(entrez_to_symbol))

# Aging direction from old vs young
aging_direction_df <- oldsedveh_v_yngsedveh_genes_06mar %>%
  dplyr::filter(Symbol %in% core_symbols_present) %>%
  dplyr::select(Symbol, aging_logFC = logFC)

# FRAP vs VEH logFC
frap_logfc_df <- oldsedfrap_v_oldsedveh_genes_08feb %>%
  dplyr::filter(Symbol %in% core_symbols_present) %>%
  dplyr::select(Symbol, frap_logFC = logFC, frap_FDR = FDR)

print(frap_logfc_df, n=93)


# Join and compute aging-directed metric
frap_bar_data <- aging_direction_df %>%
  dplyr::inner_join(frap_logfc_df, by = "Symbol") %>%
  dplyr::mutate(
    directed_logFC = sign(aging_logFC) * frap_logFC,
    bar_fill = ifelse(aging_logFC > 0, "#8B0000", "#00008B")
  )

# Row order: most amplifies aging on top -> most reverses on bottom
row_order_frap <- frap_bar_data %>%
  dplyr::arrange(desc(directed_logFC)) %>%
  dplyr::pull(Symbol) %>%
  as.character()

# Now set factor levels for barplot
frap_bar_data <- frap_bar_data %>%
  dplyr::mutate(Symbol = factor(Symbol, levels = rev(row_order_frap)))

# -------------------------------------------------------
# 3) Subset heatmap to 93 genes and apply row order
# -------------------------------------------------------
frap_aging_heatmap_mat <- logCPM_matrix_frap[core_ids_present_frap, ordered_columns_frap, drop = FALSE]
rownames(frap_aging_heatmap_mat) <- entrez_to_symbol[rownames(frap_aging_heatmap_mat)]

# Keep only genes that made it through the bar data join, in the new order
row_order_frap_heatmap <- intersect(row_order_frap, rownames(frap_aging_heatmap_mat))
frap_aging_heatmap_mat <- frap_aging_heatmap_mat[row_order_frap_heatmap, , drop = FALSE]

# Column annotation
ann_col_frap <- data.frame(
  Group = factor(c(
    rep("Old Sed VEH",  length(old_sed_samples)),
    rep("Old Sed FRAP", length(old_sedfrap_samples))
  ))
)
rownames(ann_col_frap) <- ordered_columns_frap
ann_colors_frap <- list(Group = c("Old Sed VEH" = "#999999", "Old Sed FRAP" = "#E69F00"))

# -------------------------------------------------------
# 4) Heatmap
# -------------------------------------------------------
Fig_frap_aging_heat <- pheatmap::pheatmap(
  frap_aging_heatmap_mat,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  scale             = "row",
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  row_names_side    = "left",
  fontsize_row      = 5,
  annotation_col    = ann_col_frap,
  annotation_colors = ann_colors_frap,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "93 Aging Axis Genes: Old Sed FRAP vs Old Sed VEH"
)

pdf("FigX_heatmap_frap_aging_axis_93genes.pdf", width = 5, height = 12)
print(Fig_frap_aging_heat)
dev.off()

# -------------------------------------------------------
# 5) Barplot — aging-directed logFC
#    RIGHT = amplifies aging, LEFT = reverses aging
# -------------------------------------------------------
Fig_frap_aging_bar <- ggplot(frap_bar_data, aes(x = directed_logFC, y = Symbol, fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  coord_cartesian(xlim = c(-1.2, 1.25)) +
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

pdf("FigX_barplot_frap_aging_directed_93genes.pdf", width = 10.65, height = 6.7)
print(Fig_frap_aging_bar)
dev.off()
###############################################################################
###############################################################################
#
#
#
#
############################################################
# Heatmap + Aging-Directed Barplot — 93 validated aging genes
# Old Sed IRAP vs Old Sed VEH
############################################################

# -------------------------------------------------------
# 0) Load and merge counts from two files
# -------------------------------------------------------
counts_veh  <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
counts_irap <- load_counts("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_IRAP-YNG_SED_VEH.xlsx")

old_sed_samples     <- c("T_06","T_09","T_13","T_19","T_21","T_25","T_29","T_39")
old_sedirap_samples <- c("T_49","T_50","T_51","T_52","T_53")

ordered_columns_irap <- c(old_sed_samples, old_sedirap_samples)

counts_merged_irap <- counts_veh %>%
  dplyr::select(Ensembl_noDec, dplyr::all_of(old_sed_samples)) %>%
  dplyr::inner_join(
    counts_irap %>% dplyr::select(Ensembl_noDec, dplyr::all_of(old_sedirap_samples)),
    by = "Ensembl_noDec"
  )

# Map Ensembl -> ENTREZID
counts_merged_irap$ENTREZID <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys      = counts_merged_irap$Ensembl_noDec,
  keytype   = "ENSEMBL",
  column    = "ENTREZID",
  multiVals = "first"
)
counts_merged_irap <- tidyr::drop_na(counts_merged_irap, ENTREZID)

# -------------------------------------------------------
# 1) edgeR normalization -> logCPM
# -------------------------------------------------------
counts_matrix_irap <- counts_merged_irap %>%
  dplyr::select(dplyr::all_of(ordered_columns_irap)) %>%
  as.matrix()
rownames(counts_matrix_irap) <- as.character(counts_merged_irap$ENTREZID)

group_irap <- factor(c(
  rep("OLD_SED_VEH",  length(old_sed_samples)),
  rep("OLD_SED_IRAP", length(old_sedirap_samples))
))

dge_irap <- edgeR::DGEList(counts = counts_matrix_irap, group = group_irap)
dge_irap <- edgeR::calcNormFactors(dge_irap, method = "TMM")
logCPM_matrix_irap <- edgeR::cpm(dge_irap, log = TRUE, prior.count = 1)

# -------------------------------------------------------
# 2) Compute aging-directed logFC FIRST (to set row order)
# -------------------------------------------------------

# Get symbols for the 93 genes present in the matrix
core_ids_present_irap <- core_entrez[core_entrez %in% rownames(logCPM_matrix_irap)]
entrez_to_symbol_irap <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys    = core_ids_present_irap,
  column  = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
core_symbols_present_irap <- unname(na.omit(entrez_to_symbol_irap))

# Aging direction from old vs young
aging_direction_df_irap <- oldsedveh_v_yngsedveh_genes_06mar %>%
  dplyr::filter(Symbol %in% core_symbols_present_irap) %>%
  dplyr::select(Symbol, aging_logFC = logFC)

# IRAP vs VEH logFC
irap_logfc_df <- oldsedirap_v_oldsedveh_genes_08feb %>%
  dplyr::filter(Symbol %in% core_symbols_present_irap) %>%
  dplyr::select(Symbol, irap_logFC = logFC, irap_FDR = FDR)

# Join and compute aging-directed metric
irap_bar_data <- aging_direction_df_irap %>%
  dplyr::inner_join(irap_logfc_df, by = "Symbol") %>%
  dplyr::mutate(
    directed_logFC = sign(aging_logFC) * irap_logFC,
    bar_fill = ifelse(aging_logFC > 0, "#8B0000", "#00008B")
  )

# Row order: use FRAP order instead of IRAP-specific order
row_order_irap <- row_order_frap

# Now set factor levels for barplot
irap_bar_data <- irap_bar_data %>%
  dplyr::mutate(Symbol = factor(Symbol, levels = rev(row_order_frap)))
irap_bar_data <- irap_bar_data %>%
  dplyr::filter(!is.na(Symbol))
# -------------------------------------------------------
# 3) Subset heatmap to 93 genes and apply row order
# -------------------------------------------------------
irap_aging_heatmap_mat <- logCPM_matrix_irap[core_ids_present_irap, ordered_columns_irap, drop = FALSE]
rownames(irap_aging_heatmap_mat) <- entrez_to_symbol_irap[rownames(irap_aging_heatmap_mat)]

# Keep only genes that made it through the bar data join, in the new order
row_order_irap_heatmap <- intersect(row_order_irap, rownames(irap_aging_heatmap_mat))
irap_aging_heatmap_mat <- irap_aging_heatmap_mat[row_order_irap_heatmap, , drop = FALSE]

# Column annotation
ann_col_irap <- data.frame(
  Group = factor(c(
    rep("Old Sed VEH",  length(old_sed_samples)),
    rep("Old Sed IRAP", length(old_sedirap_samples))
  ))
)
rownames(ann_col_irap) <- ordered_columns_irap
ann_colors_irap <- list(Group = c("Old Sed VEH" = "#999999", "Old Sed IRAP" = "#56B4E9"))

# -------------------------------------------------------
# 4) Heatmap
# -------------------------------------------------------
Fig_irap_aging_heat <- pheatmap::pheatmap(
  irap_aging_heatmap_mat,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  scale             = "row",
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  row_names_side    = "left",
  fontsize_row      = 5,
  annotation_col    = ann_col_irap,
  annotation_colors = ann_colors_irap,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "93 Aging Axis Genes: Old Sed IRAP vs Old Sed VEH"
)

pdf("FigX_heatmap_irap_aging_axis_93genes.pdf", width = 5, height = 12)
print(Fig_irap_aging_heat)
dev.off()

# -------------------------------------------------------
# 5) Barplot — aging-directed logFC
#    RIGHT = amplifies aging, LEFT = reverses aging
# -------------------------------------------------------
Fig_irap_aging_bar <- ggplot(irap_bar_data, aes(x = directed_logFC, y = Symbol, fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  coord_cartesian(xlim = c(-1.25, 1.25)) +
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

pdf("FigX_barplot_irap_aging_directed_93genes.pdf", width = 10.65, height = 6.7)
print(Fig_irap_aging_bar)
dev.off()


