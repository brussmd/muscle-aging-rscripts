#Reanalysis of Old Sed Veh vs Yng Sed Veh.


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





###################################################################

###################################################################
#####///Analysis of Old SED Veh vs Yng SED Veh////////###########
#-----------------------------------------------------------------#
# 1) Load edgeR *_GENE_* xlsx and build ranked vector
oldsedveh_v_yngsedveh_genes_06mar <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col   = "OLD_SED_VEH-YNG_SED_VEH_FDR"
)

#remove duplicate ENTREZID and keep most significant
oldsedveh_v_yngsedveh_genes_06mar <- oldsedveh_v_yngsedveh_genes_06mar %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()

head(oldsedveh_v_yngsedveh_genes_06mar)

oldsedveh_v_yngsedveh_genes_06mar %>%
  filter(FDR < 0.05)

oldsedveh_v_yngsedveh_genes_06mar %>%
  filter(Symbol == "Nmrk2")
#==========================================#

####################################################
############################################
#############################################################
#//////Volcano Plot for Old Sed Veh vs Yng Sed Veh///////////#
##############################################################
#############################################################
#//////Volcano Plot for Old Sed Veh vs Yng Sed Veh///////////#
##############################################################
# --- Prep: get the 93 validated genes from the volcano background data ---
highlight_genes <- oldsedveh_v_yngsedveh_genes_06mar %>%
  filter(Symbol %in% validated_aging_genes$Symbol)
highlight_up   <- highlight_genes %>% filter(logFC > 0)
highlight_down <- highlight_genes %>% filter(logFC < 0)

# --- Select genes to label ---
label_symbols <- c("Jchain", "Ccl7", "Cxcl10", "Ccl2", "Nfkb2", "Cdkn1a",
                   "Aldh1l2", "Hadh", "Sdhb", "Sod2")
label_genes <- highlight_genes %>% filter(Symbol %in% label_symbols)

# --- Volcano plot ---
aging_axis_volcano <- ggplot(oldsedveh_v_yngsedveh_genes_06mar, aes(x = logFC, y = -log10(FDR))) +
  
  geom_point(aes(color = sig_status),
             alpha = 0.4, size = 1.2) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  
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
    data = label_genes,
    aes(label = Symbol),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    nudge_y = 0.5,
    fontface = "bold"
  ) +
  
  scale_color_manual(
    values = c("Up" = "#D55E00", "Down" = "#0072B2", "NS" = "grey75")
  ) +
  
  coord_cartesian(xlim = c(-4.2, 4.2), ylim = c(-0.5, 15)) +
  
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Status",
    title = "Volcano Plot: Old Sed VEH vs Young Sed VEH (Aging)",
    subtitle = "Highlighted: 93 cross-validated aging axis genes"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title       = element_text(face = "bold")
  )

print(aging_axis_volcano)

ggsave(
  "FigX_oldsedveh_volcano_plot_v1.pdf",
  aging_axis_volcano,
  width = 9,
  height = 3.75
)

ggsave(
  "FigX_oldsedveh_volcano_plot.png",
  aging_axis_volcano,
  width = 13.5,
  height = 6.6,
  dpi = 300
)
###################################################################
#/////Mess around with volcano plot//////////#
#===============================================================#

##############################################################
#############################################################
#//////Volcano Plot for Old Sed Veh vs Yng Sed Veh///////////#
##############################################################
# --- Prep: get the 93 validated genes from the volcano background data ---
highlight_genes <- oldsedveh_v_yngsedveh_genes_06mar %>%
  filter(Symbol %in% validated_aging_genes$Symbol)
highlight_up   <- highlight_genes %>% filter(logFC > 0)
highlight_down <- highlight_genes %>% filter(logFC < 0)

# --- Select genes to label ---
label_symbols <- c("Jchain", "Ccl7", "Cxcl10", "Ccl2", "Nfkb2", "Cdkn1a",
                   "Aldh1l2", "Hadh", "Sdhb", "Sod2")
label_genes <- highlight_genes %>% filter(Symbol %in% label_symbols)

# --- Volcano plot ---
aging_axis_volcano_alt <- ggplot(oldsedveh_v_yngsedveh_genes_06mar, aes(x = logFC, y = -log10(FDR))) +
  
  geom_point(aes(color = sig_status),
             alpha = 0.4, size = 1.2) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  
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
    data = label_genes,
    aes(label = Symbol),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    nudge_y = 0.5,
    fontface = "bold"
  ) +
  
  scale_color_manual(
    values = c("Up" = "#D55E00", "Down" = "#0072B2", "NS" = "grey75")
  ) +
  
  coord_cartesian(xlim = c(-4.2, 4.2), ylim = c(-0.5, 15)) +
  
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title       = element_text(face = "bold"),
    legend.position  = "none"
  )

print(aging_axis_volcano_alt)

ggsave(
  "FigX_oldsedveh_volcano_plot_v1.pdf",
  aging_axis_volcano_alt,
  width = 9,
  height = 3.75
)
















#######################################################################
######################################################################
core_entrez






#####################################################################
#////GO Enrichment//////////////#
#########################################################################################
create_vec_from_df(oldsedveh_v_yngsedveh_genes_06mar, "oldsedveh_v_yngsedveh_genes_06mar")

# 2) GSEA (GO Biological Process)
gseGO_oldsedveh_vs_yngsedveh_06mar.OUTPUT <- gseGO(
  geneList     = oldsedveh_v_yngsedveh_genes_06mar.vec, 
  ont          = "BP",
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 300,
  pvalueCutoff = 0.5,      # capture everything, filter later
  eps          = 1e-30,      # better estimation for very small p-values
  verbose      = FALSE
)

head(gseGO_oldsedveh_vs_yngsedveh_06mar.OUTPUT@result)


# Simplify, then compute Count/geneRatio again (simplify changes the table)
oldsedveh_vs_yngsedveh_06mar_simpl <- clusterProfiler::simplify(gseGO_oldsedveh_vs_yngsedveh_06mar.OUTPUT,
                                                     cutoff = 0.5, by = "p.adjust", select_fun = min)

oldsedveh_vs_yngsedveh_06mar_simpl@result

oldsedveh_vs_yngsedveh_06mar_simpl.df <- as.data.frame(oldsedveh_vs_yngsedveh_06mar_simpl@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

dim(oldsedveh_vs_yngsedveh_06mar_simpl.df %>%
      filter(p.adjust <0.05))

oldsedveh_vs_yngsedveh_simpl_06mar_goplot <- oldsedveh_vs_yngsedveh_06mar_simpl.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  scale_size_continuous(range = c(2, 12)) +  # Double the size (default is ~1-6)
  theme_bw()


# Display it
print(oldsedveh_vs_yngsedveh_simpl_06mar_goplot)

oldsedveh_vs_yngsedveh_06mar_simpl.df %>% 
  dplyr::select(ID, Description, NES, Count)
###############################################################################



##########################################################################
#///////GO Enrichment Bubble Plot///////////////////////#
############################################################################
###########################################################################

# Add pathway categories
oldsedveh_vs_yngsedveh_06mar_simpl.df <- oldsedveh_vs_yngsedveh_06mar_simpl.df %>%
  mutate(
    Category = case_when(
      # Immune & Inflammation (very broad pattern)
      grepl("immune|lymphocyte|leukocyte|antigen|phagocyt|inflammasome|toll-like|inflammatory|defense response|viral|organism|killing|burst", 
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
################################################################
################################################################

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
df_oldsed_06mar_sig <- oldsedveh_vs_yngsedveh_06mar_simpl.df %>%
  dplyr::filter(p.adjust < 0.05)

# Optional: shorten long GO terms for cleaner labels
df_oldsed_06mar_sig <- df_oldsed_06mar_sig %>%
  dplyr::mutate(
    label = stringr::str_replace_all(Description, "\\s+", " "),
    label = stringr::str_trunc(label, 45)
  )


# --- 2) Choose the dead-zone cutoff (same as before) ---
cut <- 1.5
x_span <- 1.0              # width of each side window
left_lim  <- c(-(cut + x_span), -cut)   # -2.5 to -1.5
right_lim <- c(cut,  cut + x_span)      #  1.5 to  2.5


df_oldsed_06mar_sig_left  <- df_oldsed_06mar_sig %>% dplyr::filter(NES <= -cut)
df_oldsed_06mar_sig_right <- df_oldsed_06mar_sig %>% dplyr::filter(NES >=  cut)

# --- 3) Pick top 5 by Count within each category & side ---
lab_mito_oldsed_06mar <- df_oldsed_06mar_sig_left %>%
  dplyr::filter(Category == "Metabolism & Mitochondria") %>%
  dplyr::slice_max(order_by = Count, n = 5, with_ties = FALSE)

lab_imm_oldsed_06mar <- df_oldsed_06mar_sig_right %>%
  dplyr::filter(Category == "Immune & Inflammation")

fill_vals <- c(
  "Immune & Inflammation" = "#E41A1C",
  "Metabolism & Mitochondria" = "#377EB8",
  "Cell Cycle & Division" = "#4DAF4A",
  "Cytoskeleton & Movement" = "#984EA3",
  "Hemostasis & Coagulation" = "#FF7F00",
  "Development & Differentiation" = "#F781BF",
  "Other Cellular Processes" = "#999999"
)

# Now convert Category to factor with all levels
df_oldsed_06mar_sig <- df_oldsed_06mar_sig %>%
  mutate(Category = factor(Category, levels = names(fill_vals)))

# Re-split after factoring
df_oldsed_06mar_sig_left  <- df_oldsed_06mar_sig %>% dplyr::filter(NES <= -cut)
df_oldsed_06mar_sig_right <- df_oldsed_06mar_sig %>% dplyr::filter(NES >=  cut)

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

fill_scale <- scale_fill_manual(
  values = fill_vals,
  breaks = names(fill_vals),
  limits = names(fill_vals),
  drop = FALSE,
  name = "Pathway Category"
)

size_scale <- scale_size_continuous(
  range = c(3, 15),
  name = "Gene Count"
)

fill_guide <- guides(
  fill = guide_legend(
    override.aes = list(
      shape = 21,
      size = 5,
      colour = "black",
      alpha = 1,
      stroke = 0.6
    )
  )
)


# --- Create dummy rows for all categories (placed off-screen) ---
dummy_rows <- data.frame(
  NES       = right_lim[1],     # within x limits but...
  geneRatio = -1,               # off-screen on y (below 0)
  Count     = 1,
  Category  = factor(names(fill_vals), levels = names(fill_vals))
)

# --- 5) Left panel (negative NES) ---
p_left_06mar <- ggplot(df_oldsed_06mar_sig_left,
                       aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_mito_oldsed_06mar,
    aes(label = label),
    size = 3.3, fontface = "bold",
    box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category", drop = FALSE) +
  scale_x_continuous(limits = left_lim,
                     breaks = seq(left_lim[1], left_lim[2], by = 0.25)) +
  y_scale_shared +
  labs(x = "NES (negative)", y = "Gene Ratio (Count / Set Size)") +
  base_theme +
  theme(legend.position = "none")

# --- 6) Right panel (positive NES) — with dummy rows ---
p_right_06mar <- ggplot(df_oldsed_06mar_sig_right,
                        aes(x = NES, y = geneRatio, size = Count, fill = Category)) +
  geom_point(alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
  geom_point(data = dummy_rows,
             aes(x = NES, y = geneRatio, fill = Category),
             alpha = 0, size = 0) +
  geom_text_repel(
    data = lab_imm_oldsed_06mar,
    aes(label = label),
    size = 3.3, fontface = "bold",
    box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category", drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1, shape = 21))) +
  scale_x_continuous(limits = right_lim,
                     breaks = seq(right_lim[1], right_lim[2], by = 0.25)) +
  y_scale_shared +
  labs(x = "NES (positive)", y = NULL,
       title = "GO Pathway Enrichment: Old Sed Veh vs Yng Sed Veh",
       subtitle = paste0("p.adjust < 0.05 | Bubble size = gene count | Center region (",
                         -cut, " to ", cut, ") omitted")) +
  base_theme

# --- 7) Combine ---
p_combined_06mar <- p_left_06mar + p_right_06mar +
  plot_layout(widths = c(1, 1))

p_combined_06mar


# Save as PNG
ggsave(
  "FigX_oldsedveh_go_enrichment.png",
  p_combined_06mar,
  width = 13.5,
  height = 6.6,
  dpi = 300
)

# Save as PDF
ggsave(
  "FigX_oldsedveh_go_enrichment.pdf",
  p_combined_06mar,
  width = 13.5,
  height = 6.6
)

ggsave(
  "FigX_oldsedveh_go_enrichment_v1.pdf",
  p_combined_06mar,
  width = 8.3,
  height = 2.4
)
###########################################################################
############################################################################

#=============================================================#
#//////////GO Enrichment Plot Resized/////////////////////////#
#=============================================================#
#########################################################################
#
#Bubbleplot alternate view.
#
#
# Shared y-axis
y_count_shared <- scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 20))

# --- Left panel ---
p_left_06mar_alt <- ggplot(df_oldsed_06mar_sig_left,
                           aes(x = NES, y = Count, size = geneRatio, fill = Category)) +
  geom_point(alpha = 1, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    data = lab_mito_oldsed_06mar,
    aes(label = label),
    size = 3.3, fontface = "bold",
    box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(10, 25), name = "Gene Ratio") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category", drop = FALSE) +
  scale_x_continuous(limits = left_lim,
                     breaks = seq(left_lim[1], left_lim[2], by = 0.25)) +
  y_count_shared +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(legend.position = "none")

# --- Right panel ---
p_right_06mar_alt <- ggplot(df_oldsed_06mar_sig_right,
                            aes(x = NES, y = Count, size = geneRatio, fill = Category)) +
  geom_point(alpha = 1, shape = 21, color = "black", stroke = 0.6) +
  geom_point(data = dummy_rows %>% mutate(geneRatio = 0.1),
             aes(x = NES, y = Count, fill = Category),
             alpha = 0, size = 0) +
  geom_text_repel(
    data = lab_imm_oldsed_06mar,
    aes(label = label),
    size = 3.3, fontface = "bold",
    box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(10, 25), name = "Gene Ratio") +
  scale_fill_manual(values = fill_vals, name = "Pathway Category", drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1, shape = 21))) +
  scale_x_continuous(limits = right_lim,
                     breaks = seq(right_lim[1], right_lim[2], by = 0.25)) +
  y_count_shared +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(legend.position = "none")

# --- Combine ---
p_combined_06mar_alt <- p_left_06mar_alt + p_right_06mar_alt +
  plot_layout(widths = c(1, 1))

print(p_combined_06mar_alt)


# Save as PDF
ggsave(
  "FigX_oldsedveh_go_enrichment_alt.pdf",
  p_combined_06mar_alt,
  width = 13.44,
  height = 4.4292
)

p_for_legend <- ggplot(df_oldsed_06mar_sig_right,
                       aes(x = NES, y = Count, size = geneRatio, fill = Category)) +
  geom_point(shape = 21, color = "black", stroke = 0.6) +
  scale_size_continuous(range = c(10, 25), name = "Gene Ratio") +
  scale_fill_manual(values = fill_vals, guide = "none") +
  guides(size = guide_legend(direction = "vertical")) +
  theme_minimal() +
  theme(legend.position = "right")

ratio_legend <- cowplot::get_legend(p_for_legend)
cowplot::plot_grid(ratio_legend)
ggsave("FigX_gene_ratio_legend.pdf", cowplot::plot_grid(ratio_legend), 
       width = 2, height = 4)
###########################################################################
##########################################################################

old_core_deg_df <- oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  filter(oldsedveh_FDR <0.05)

old_core_deg_df
###########################################################################


# =============================================================================
# Cross-validate aging genes with external MusAge dataset (female Gas + TA)
# =============================================================================

# --- 1) Load the external validation data ---
load("all_muscle_age_dfs.RData")

# Extract female gastroc and TA (the full dataframes, not just MusAge subset)
val_gas <- all_mus_age_dfs$gas_age_df
val_ta  <- all_mus_age_dfs$ta_age_df

# --- 2) Classify direction in validation set using late-age logFC ---
#     Use logFC_27 as the primary direction call; fall back to logFC_24 if needed.
#     Also use SpearmanRho as a secondary confirmation.

classify_direction_val <- function(df, muscle_label) {
  df %>%
    mutate(
      # Direction based on the latest timepoint available
      dir_late = case_when(
        !is.na(logFC_27) & logFC_27 > 0 ~ "Up",
        !is.na(logFC_27) & logFC_27 < 0 ~ "Down",
        !is.na(logFC_24) & logFC_24 > 0 ~ "Up",
        !is.na(logFC_24) & logFC_24 < 0 ~ "Down",
        TRUE ~ "Ambiguous"
      ),
      # Confirm with rho direction
      dir_rho = case_when(
        SpearmanRho >  0.3 ~ "Up",
        SpearmanRho < -0.3 ~ "Down",
        TRUE ~ "Ambiguous"
      ),
      # Consensus direction: both agree
      val_direction = case_when(
        dir_late == "Up"   & dir_rho %in% c("Up", "Ambiguous")   ~ "Up",
        dir_late == "Down" & dir_rho %in% c("Down", "Ambiguous") ~ "Down",
        dir_late == dir_rho ~ dir_late,
        TRUE ~ "Ambiguous"
      ),
      val_muscle = muscle_label
    ) %>%
    # Also flag age-responsive using the corrected abs(rho) criteria
    mutate(
      sig_21 = if ("FDR_21" %in% names(.)) FDR_21 < 0.05 else FALSE,
      sig_24 = if ("FDR_24" %in% names(.)) FDR_24 < 0.05 else FALSE,
      sig_27 = if ("FDR_27" %in% names(.)) FDR_27 < 0.05 else FALSE,
      old_sig_count = sig_21 + sig_24 + sig_27,
      age_responsive = (old_sig_count >= 2) | (abs(SpearmanRho) > 0.3)
    ) %>%
    dplyr::select(GeneName, val_muscle, val_direction, age_responsive,
           SpearmanRho, starts_with("logFC_"), starts_with("FDR_"))
}

val_gas_dir <- classify_direction_val(val_gas, "Gastroc")
val_ta_dir  <- classify_direction_val(val_ta,  "TA")

# Combine and keep only age-responsive genes
val_combined <- bind_rows(val_gas_dir, val_ta_dir) %>%
  filter(age_responsive == TRUE)

# --- 3) Classify direction in your dataset ---
old_core_deg_df <- old_core_deg_df %>%
  mutate(
    my_direction = case_when(
      oldsedveh_logFC > 0 ~ "Up",
      oldsedveh_logFC < 0 ~ "Down",
      TRUE ~ "Ambiguous"
    )
  )

# --- 4) Join: your 119 genes × validation (Gas + TA) ---
cross_val_df <- old_core_deg_df %>%
  inner_join(
    val_combined,
    by = c("Symbol" = "GeneName")
  ) %>%
  mutate(
    concordance = case_when(
      my_direction == val_direction                          ~ "Concordant",
      my_direction != val_direction & val_direction != "Ambiguous" 
      & my_direction != "Ambiguous" ~ "Discordant",
      TRUE                                                   ~ "Ambiguous"
    )
  )

# --- 5) Summary ---
cat("=== Concordance summary ===\n")
cross_val_df %>%
  count(val_muscle, concordance) %>%
  tidyr::pivot_wider(names_from = concordance, values_from = n, values_fill = 0) %>%
  print()

cat("\n=== Concordant genes (same direction in both datasets) ===\n")
concordant_genes <- cross_val_df %>%
  filter(concordance == "Concordant") %>%
  select(Symbol, my_direction, val_direction, val_muscle, 
         oldsedveh_logFC, oldsedveh_FDR, SpearmanRho) %>%
  arrange(val_muscle, my_direction, desc(abs(oldsedveh_logFC)))
print(concordant_genes, n = 50)

cat("\n=== Discordant genes (opposite direction) ===\n")
discordant_genes <- cross_val_df %>%
  filter(concordance == "Discordant") %>%
  select(Symbol, my_direction, val_direction, val_muscle,
         oldsedveh_logFC, oldsedveh_FDR, SpearmanRho) %>%
  arrange(val_muscle, Symbol)
print(discordant_genes, n = 50)

# --- 6) Genes validated in BOTH gastroc and TA ---
cat("\n=== Concordant in both Gastroc AND TA ===\n")
in_both <- concordant_genes %>%
  count(Symbol, my_direction) %>%
  filter(n == 2)
print(in_both, n = 50)

# --- 7) Final validated aging gene set ---
# Concordant in at least one muscle
validated_aging_genes <- concordant_genes %>%
  distinct(Symbol, my_direction) %>%
  rename(aging_direction = my_direction)

cat("\n=== Final validated aging gene count ===\n")
cat("Total:", nrow(validated_aging_genes), "\n")
validated_aging_genes %>% count(aging_direction) %>% print()

validated_aging_genes
write.csv(validated_aging_genes, file = "validated_aging_genes.csv", row.names = FALSE)

in_both %>%
  filter(my_direction == "Up")

MusAge_geneset <- read.csv("MusAge_geneset.csv")      

MusAge_geneset

validated_aging_genes %>%
  filter(aging_direction == "Up")

# Overlap between MusAge (52 genes) and your validated set (93 genes)
overlap <- intersect(MusAge_geneset$Gene, validated_aging_genes$Symbol)
cat("Overlap:", length(overlap), "genes\n")
print(overlap)

# What's in MusAge but NOT validated
musage_only <- setdiff(MusAge_geneset$Gene, validated_aging_genes$Symbol)
cat("\nMusAge only (not validated):", length(musage_only), "\n")
print(musage_only)

# What's in your validated set but NOT in MusAge
validated_only <- setdiff(validated_aging_genes$Symbol, MusAge_geneset$Gene)
cat("\nValidated only (not in MusAge):", length(validated_only), "\n")
print(validated_only)
###############################################################################
###############################################################################
#
#Targeted GSEA Analysis ꜜ
#
#
#############################################################
#///////TARGETED GSEA ANALYSES for each intervention////////#
############################################################

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
#############################################################################
############################################################################

read_clean_pvalue_xlsx <- function(xlsx_file, logFC_col, FDR_col, PValue_col = NULL) {
  cols_to_select <- c("Ensembl", "Symbol", logFC_col, FDR_col)
  if (!is.null(PValue_col)) cols_to_select <- c(cols_to_select, PValue_col)
  
  df <- read_xlsx(xlsx_file) %>%
    dplyr::select(all_of(cols_to_select))
  
  # Rename consistently
  df <- df %>%
    dplyr::rename(logFC = all_of(logFC_col),
                  FDR   = all_of(FDR_col))
  if (!is.null(PValue_col)) df <- df %>% dplyr::rename(PValue = all_of(PValue_col))
  
  df <- df %>%
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
###################################################################

names(read_xlsx("20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx"))

oldpwrveh_v_oldsedveh_genes_06mar <- read_clean_pvalue_xlsx(
  xlsx_file  = "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
  logFC_col  = "OLD_PWR_VEH-OLD_SED_VEH_logFC",
  FDR_col    = "OLD_PWR_VEH-OLD_SED_VEH_FDR",
  PValue_col = "OLD_PWR_VEH-OLD_SED_VEH_PValue"
)

oldpwrveh_v_oldsedveh_genes_06mar <- oldpwrveh_v_oldsedveh_genes_06mar %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()




#===================================================================#
in_both
oldpwrveh_v_oldsedveh_genes_06mar

# --- 1) Get ENTREZIDs for the 54-gene set ---
core_entrez <- old_core_deg_df %>%
  filter(Symbol %in% validated_aging_genes$Symbol) %>%
  pull(ENTREZID) %>%
  as.character() %>%
  unique()

core_entrez

in_both_entrez <- old_core_deg_df %>%
  filter(Symbol %in% in_both$Symbol) %>%
  pull(ENTREZID) %>%
  as.character() %>%
  unique()

in_both_entrez

# --- 3) Run it ---
pwrveh_gsea_result <- run_targeted_aging_gsea(
  intervention_df   = oldpwrveh_v_oldsedveh_genes_06mar,
  aging_df          = oldsedveh_v_yngsedveh_genes_06mar,
  core_genes        = core_entrez,
  intervention_name = "OldPWRVEH_vs_OldSEDVEH"
)

print(pwrveh_gsea_result)

alternate_pwrveh_gsea_result <- run_targeted_aging_gsea(
  intervention_df   = oldpwrveh_v_oldsedveh_genes_06mar,
  aging_df          = oldsedveh_v_yngsedveh_genes_06mar,
  core_genes        = in_both_entrez,
  intervention_name = "OldPWRVEH_vs_OldSEDVEH"
)

print(alternate_pwrveh_gsea_result)
###############################################################################
#############################################################################
#//////old sed FRAP GSEA result/////////#

oldsedfrap_v_oldsedveh_genes_06mar <- read_clean_pvalue_xlsx(
  xlsx_file  = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_FRAP-OLD_SED_VEH.xlsx",
  logFC_col  = "OLD_SED_FRAP-OLD_SED_VEH_logFC",
  FDR_col    = "OLD_SED_FRAP-OLD_SED_VEH_FDR",
  PValue_col = "OLD_SED_FRAP-OLD_SED_VEH_PValue"
)

oldsedfrap_v_oldsedveh_genes_06mar <- oldsedfrap_v_oldsedveh_genes_06mar %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()

# --- 3) Run it ---
sedfrap_gsea_result <- run_targeted_aging_gsea(
  intervention_df   = oldsedfrap_v_oldsedveh_genes_06mar,
  aging_df          = oldsedveh_v_yngsedveh_genes_06mar,
  core_genes        = core_entrez,
  intervention_name = "OldSEDFRAP_vs_OldSEDVEH"
)

print(sedfrap_gsea_result)

# --- 3) Run it ---
alt_sedfrap_gsea_result <- run_targeted_aging_gsea(
  intervention_df   = oldsedfrap_v_oldsedveh_genes_06mar,
  aging_df          = oldsedveh_v_yngsedveh_genes_06mar,
  core_genes        = in_both_entrez,
  intervention_name = "OldSEDFRAP_vs_OldSEDVEH"
)

print(alt_sedfrap_gsea_result)
##############################################################

###############################################################################
#############################################################################
#//////old sed IRAP GSEA result/////////#

oldsedirap_v_oldsedveh_genes_07mar <- read_clean_pvalue_xlsx(
  xlsx_file  = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_IRAP-OLD_SED_VEH.xlsx",
  logFC_col  = "OLD_SED_IRAP-OLD_SED_VEH_logFC",
  FDR_col    = "OLD_SED_IRAP-OLD_SED_VEH_FDR",
  PValue_col = "OLD_SED_IRAP-OLD_SED_VEH_PValue"
)

oldsedirap_v_oldsedveh_genes_07mar <- oldsedirap_v_oldsedveh_genes_07mar %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()

# --- 3) Run it ---
sedirap_gsea_result <- run_targeted_aging_gsea(
  intervention_df   = oldsedirap_v_oldsedveh_genes_07mar,
  aging_df          = oldsedveh_v_yngsedveh_genes_06mar,
  core_genes        = core_entrez,
  intervention_name = "OldSEDIRAP_vs_OldSEDVEH"
)

print(sedirap_gsea_result)
##############################################################



#############################################################
length(core_entrez)
oldsedveh_v_yngsedveh_genes_06mar

validated_aging_genes %>% 
  filter(aging_direction == "Up")

write_csv(validated_aging_genes, file = "Transciptional_aging_sig_93_genes.csv")

print(validated_aging_genes, n=100)

print(oldsedveh_v_yngsedveh_genes_06mar %>%
  filter(ENTREZID %in% core_entrez) %>%
  arrange(desc(abs(logFC))), n= 100)
#############################################################
###############################################################
#
#
#
#
#
#
########################################################
#///////Heatmap + logFC plot///////////////////////////#
########################################################

## ============================================================
## Heatmap + logFC barplot — 93 validated aging axis genes
## Old Sed VEH vs Yng Sed VEH
## ============================================================

suppressPackageStartupMessages({
  library(cowplot)
})

# -------------------------------------------------------
# 0) Map core_entrez to Symbols for row ordering
# -------------------------------------------------------
core_symbols <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys    = core_entrez,
  column  = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
) %>% unname() %>% unique() %>% na.omit()

# Order genes by logFC from the aging comparison
# (most downregulated at top, most upregulated at bottom — or reverse as you prefer)
aging_gene_order <- oldsedveh_v_yngsedveh_genes_06mar %>%
  dplyr::filter(Symbol %in% core_symbols) %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::pull(Symbol)

# -------------------------------------------------------
# 1) Build logCPM heatmap matrix
# -------------------------------------------------------
group_vec_aging <- c(rep("YNG_SED", length(yng_sed_samples)),
                     rep("OLD_SED", length(old_sed_samples)))
order_vec_aging <- c(yng_sed_samples, old_sed_samples)

res_aging <- counts_to_logCPM_by_symbol(
  counts_xlsx  = "20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  sample_order = order_vec_aging,
  group_names  = group_vec_aging
)

# Subset to the 93 validated genes and enforce row order
aging_logCPM_mat <- res_aging$logCPM[rownames(res_aging$logCPM) %in% core_symbols, , drop = FALSE]
row_order_aging <- intersect(aging_gene_order, rownames(aging_logCPM_mat))
aging_logCPM_mat <- aging_logCPM_mat[row_order_aging, res_aging$samples, drop = FALSE]

# Column annotations
ann_col_aging <- data.frame(Group = factor(res_aging$groups, levels = c("YNG_SED", "OLD_SED")))
rownames(ann_col_aging) <- res_aging$samples
ann_colors_aging <- list(Group = c(YNG_SED = "#4C78A8", OLD_SED = "#7F7F7F"))

# -------------------------------------------------------
# 2) Heatmap
# -------------------------------------------------------
Fig_aging_heat <- pheatmap::pheatmap(
  aging_logCPM_mat,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  scale           = "row",
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  row_names_side  = "left",
  fontsize_row    = 5,          # smaller font for 93 genes
  annotation_col  = ann_col_aging,
  annotation_colors = ann_colors_aging,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

pdf("FigX_heatmap_aging_axis_93genes.pdf", width = 4.2, height = 12)
print(Fig_aging_heat)
dev.off()
#============================================================#
#////aging heatmap alternate///////#

Fig_aging_heat_alt <- pheatmap::pheatmap(
  aging_logCPM_mat,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  scale           = "row",
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  legend          = FALSE,
  fontsize_row    = 5,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

pdf("FigX_heatmap_aging_axis_93genes.pdf", width = 4.4, height = 8.4)
print(Fig_aging_heat_alt)
dev.off()



# -------------------------------------------------------
# 3) Barplot — absolute logFC ("aged" direction, all positive)
# -------------------------------------------------------
aging_bar_data <- oldsedveh_v_yngsedveh_genes_06mar %>%
  dplyr::filter(Symbol %in% row_order_aging) %>%
  dplyr::left_join(validated_aging_genes, by = c("Symbol" = "Symbol")) %>%
  dplyr::mutate(
    abs_logFC = abs(logFC),
    Symbol    = factor(Symbol, levels = rev(row_order_aging)),
    bar_fill  = ifelse(aging_direction == "Up", "#8B0000", "#00008B")
  )

Fig_aging_bar <- ggplot(aging_bar_data, aes(x = abs_logFC, y = Symbol, fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(x = "Absolute log2 Fold Change (Old vs Young)",
       y = NULL) +
  theme(
    axis.text.y   = element_blank(),
    axis.text.x   = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor    = element_blank()
  )

pdf("FigX_barplot_aging_axis_93genes.pdf", width = 4.2, height = 8.4)
print(Fig_aging_bar)
dev.off()
#####################################################################
#########################################################################


# =========================================================
# DIAGNOSTIC: What do the 39 extra genes contribute?
# =========================================================

# --- 1) Identify the 39 extra genes ---
only_in_93 <- setdiff(core_entrez, in_both_entrez)
only_in_93_symbols <- AnnotationDbi::mapIds(
  org.Mm.eg.db, keys = only_in_93,
  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
)
cat("Extra genes (in 93 but not 54):", length(only_in_93), "\n")

# --- 2) Get the rank_metric for every gene in the ranked list ---
aging_dir <- oldsedveh_v_yngsedveh_genes_06mar %>%
  mutate(aging_direction = sign(logFC)) %>%
  select(Ensembl_noDec, aging_direction)

ranked_full <- oldpwrveh_v_oldsedveh_genes_06mar %>%
  left_join(aging_dir, by = "Ensembl_noDec") %>%
  mutate(rank_metric = aging_direction * logFC * ((-log10(PValue + 1e-300))^0.25)) %>%
  arrange(ENTREZID, PValue) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  drop_na(rank_metric) %>%
  arrange(desc(rank_metric)) %>%
  mutate(rank_position = row_number())

total_genes <- nrow(ranked_full)

# --- 3) Where do each set's genes sit in the ranked list? ---
ranks_54 <- ranked_full %>% filter(ENTREZID %in% in_both_entrez)
ranks_39 <- ranked_full %>% filter(ENTREZID %in% only_in_93)
ranks_93 <- ranked_full %>% filter(ENTREZID %in% core_entrez)

cat("\n=== Median rank position (out of", total_genes, ") ===\n")
cat("93-gene set:", median(ranks_93$rank_position), "\n")
cat("54-gene set:", median(ranks_54$rank_position), "\n")
cat("39 extra genes:", median(ranks_39$rank_position), "\n")
cat("(Middle of list =", round(total_genes/2), ")\n")

cat("\n=== Mean rank_metric ===\n")
cat("93-gene set:", round(mean(ranks_93$rank_metric), 4), "\n")
cat("54-gene set:", round(mean(ranks_54$rank_metric), 4), "\n")
cat("39 extra genes:", round(mean(ranks_39$rank_metric), 4), "\n")

# --- 4) Show the 39 genes with their rank metrics and positions ---
cat("\n=== 39 extra genes: ranked by rank_metric ===\n")
ranks_39 %>%
  select(Symbol, logFC, PValue, aging_direction, rank_metric, rank_position) %>%
  arrange(rank_metric) %>%
  print(n = 39)

# --- 5) Show the 54 genes for comparison ---
cat("\n=== 54-gene set: ranked by rank_metric ===\n")
ranks_54 %>%
  select(Symbol, logFC, PValue, aging_direction, rank_metric, rank_position) %>%
  arrange(rank_metric) %>%
  print(n = 54)

# --- 6) Visual: where do genes fall in the ranked list? ---
rank_comparison <- bind_rows(
  ranks_54 %>% mutate(set = "54-gene set") %>% select(Symbol, rank_metric, rank_position, set),
  ranks_39 %>% mutate(set = "39 extra genes") %>% select(Symbol, rank_metric, rank_position, set)
)

Fig_rank_dist <- ggplot(rank_comparison, aes(x = rank_metric, fill = set)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("54-gene set" = "#999999", "39 extra genes" = "#E41A1C")) +
  labs(
    x = "Rank Metric (aging_direction × logFC × significance weight)\n← Reverses Aging | Amplifies Aging →",
    y = "Density",
    title = "Where do the gene sets fall in the exercise ranking?"
  ) +
  theme_minimal()

print(Fig_rank_dist)

# --- 7) Wilcoxon test: are the 39 genes more extreme than the 54? ---
cat("\n=== Wilcoxon test: rank_metric 39 vs 54 ===\n")
wilcox.test(ranks_39$rank_metric, ranks_54$rank_metric) %>% print()
###############################################################################
##############################################################################
#
#       #////LIPIDS: OLD SED VEH vs YNG SED VEh//////////#
#
#
#
##############################################################################

##############################################################################
# Fig. 3a — Lipid dotplot (OLD SED vs YNG SED) — REVISED
##############################################################################
library(tidyverse)

# 0) Load
lip_raw <- readr::read_csv("Konopka_muscle_lipids.csv", show_col_types = FALSE)

lip_raw
colnames(lip_raw)

# 1) Identify lipid columns (everything except Sample, Group)
lipid_cols <- setdiff(names(lip_raw), c("Sample","Group"))

# 2) Impute zeros with half of the per-lipid minimum non-zero
min_nonzero <- sapply(lip_raw[lipid_cols], function(x) {
  m <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (is.infinite(m)) NA_real_ else m
})
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

# 4) Drop YS6
lip_log <- lip_log %>% filter(Sample != "YS6")

# 5) Subset to OS vs YS
os_ys <- lip_log %>% filter(Group %in% c("OS","YS"))

# 6) Per-lipid stats
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

# 7) Lipid class extraction with cleanup
os_ys_stats <- os_ys_stats %>%
  mutate(
    LipidClass = stringr::str_trim(stringr::str_extract(Lipid, "^[^\\(]+")),
    # --- Reclassify specific lipids ---
    LipidClass = case_when(
      grepl("^LPC", Lipid)                         ~ "LPC",
      grepl("^HexCer|^Hex2Cer", Lipid)             ~ "HexCer",
      LipidClass == "COH"                           ~ NA_character_,
      LipidClass == "Ubiquinone"                    ~ NA_character_,
      TRUE                                          ~ LipidClass
    ),
    # --- New signature criteria: P_value < 0.05 ---
    IsSignature = P_value < 0.05,
    neg_log10_p = -log10(P_value)
  ) %>%
  # Remove COH
  filter(!is.na(LipidClass))

# 8) Class order: by number of signature lipids (most on top, fewest on bottom)
class_order <- os_ys_stats %>%
  group_by(LipidClass) %>%
  summarize(
    n_sig = sum(IsSignature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n_sig, LipidClass) %>%    # tiebreaker: alphabetical within same count
  pull(LipidClass)

# Verify the order
cat("=== Class order (bottom to top on plot) ===\n")
os_ys_stats %>%
  group_by(LipidClass) %>%
  summarize(n_sig = sum(IsSignature), .groups = "drop") %>%
  arrange(n_sig, LipidClass) %>%
  print(n = 20)

# 9) Category colors
class_colors <- c(
  "CE"      = "#4B0082",
  "Cer"     = "#FF7F00",
  "DG"      = "#984EA3",
  "FA"      = "#A65628",
  "HexCer"  = "#F781BF",
  "LPC"     = "#377EB8",
  "LPE"     = "#4DAF4A",
  "PA"      = "#8C564B",
  "PC"      = "#66C2A5",
  "PE"      = "#FC8D62",
  "PG"      = "#8DA0CB",
  "PI"      = "#E78AC3",
  "PS"      = "#A6D854",
  "SM"      = "#FFD92F",
  "Sph"     = "#E5C494",
  "TG"      = "#B3B3B3"
)
# Add Ubiquinone if still present
if ("Ubiquinone" %in% unique(os_ys_stats$LipidClass)) {
  class_colors["Ubiquinone"] <- "#66BD63"
}

plot_df <- os_ys_stats %>%
  mutate(LipidClass = factor(LipidClass, levels = class_order))

# 10) Dotplot
oldsed_lipids_fig1d <- ggplot(plot_df, aes(x = Log2_FC_OS_vs_YS, y = LipidClass)) +
  geom_jitter(
    data = subset(plot_df, !IsSignature),
    aes(size = neg_log10_p, color = LipidClass),
    width = 0, height = 0.22, alpha = 1, shape = 16
  ) +
  # Signature: increased with age (dark red)
  geom_jitter(
    data = subset(plot_df, IsSignature & Log2_FC_OS_vs_YS > 0),
    aes(size = neg_log10_p),
    width = 0, height = 0.22,
    shape = 21, fill = "#8B0000", color = "black", stroke = 0.6, alpha = 0.9
  ) +
  # Signature: decreased with age (dark blue)
  geom_jitter(
    data = subset(plot_df, IsSignature & Log2_FC_OS_vs_YS < 0),
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
    x = "log2 Fold Change (OLD SED − YNG SED)",
    y = "Lipid Class",
    title = "OLD SED vs YNG SED Lipidomics",
    subtitle = "Dark red = sig up with age | Dark blue = sig down with age | Circle size = significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

print(oldsed_lipids_fig1d)

# Save as PDF
ggsave(
  "Fig1d_oldsedveh_lipid_classes.pdf",
  oldsed_lipids_fig1d,
  width = 13.5,
  height = 6.6
)
#####################################################################
#////Adjusted Axis Lipid plot//////////#
#----------------------------------------#

# Define clip threshold
x_clip <- 1.75

# Flag outliers for annotation
plot_df <- plot_df %>%
  mutate(
    x_display = pmin(pmax(Log2_FC_OS_vs_YS, -x_clip), x_clip),
    is_clipped = abs(Log2_FC_OS_vs_YS) > x_clip,
    clip_label = ifelse(is_clipped, sprintf("%.1f", Log2_FC_OS_vs_YS), NA_character_)
  )

# 10) Dotplot with clipped axis
oldsed_lipids_fig1d <- ggplot(plot_df, aes(x = x_display, y = LipidClass)) +
  geom_jitter(
    data = subset(plot_df, !IsSignature),
    aes(size = neg_log10_p, color = LipidClass),
    width = 0, height = 0.22, alpha = 1, shape = 16
  ) +
  geom_jitter(
    data = subset(plot_df, IsSignature & Log2_FC_OS_vs_YS > 0),
    aes(size = neg_log10_p),
    width = 0, height = 0.22,
    shape = 21, fill = "#8B0000", color = "black", stroke = 0.6, alpha = 0.9
  ) +
  geom_jitter(
    data = subset(plot_df, IsSignature & Log2_FC_OS_vs_YS < 0),
    aes(size = neg_log10_p),
    width = 0, height = 0.22,
    shape = 21, fill = "#00008B", color = "black", stroke = 0.6, alpha = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_x_continuous(
    limits = c(-0.75, x_clip + 0),
    breaks = seq(-0.75, 1.5, by = 0.5)
  ) +
  scale_y_discrete(limits = class_order) +
  scale_color_manual(values = class_colors, guide = "none") +
  scale_size_continuous(
    range = c(3, 9),
    name = expression(-log[10](P))
  ) +
  labs(
    x = "log2 Fold Change (OLD SED − YNG SED)",
    y = "Lipid Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position  = "none"
  )

print(oldsed_lipids_fig1d)

# Save as PDF
ggsave(
  "Fig1d_oldsedveh_lipid_classes_v1.pdf",
  oldsed_lipids_fig1d,
  width = 12,
  height = 8.1
)



###########################################################################
##########################################################################
os_ys_stats %>%
  filter(IsSignature)

## =========================
## Fig. 3b — 59-lipid signature heatmap (YS vs OS) + bar plot
## =========================
library(pheatmap)
library(ggplot2)
library(forcats)

# 1) Define the 59-lipid Aging Signature from updated criteria
sig_lipid_features <- os_ys_stats %>%
  dplyr::filter(IsSignature == TRUE) %>%
  dplyr::arrange(Lipid)

sig_lipid_features

signature59 <- sig_lipid_features$Lipid
signature59

write.csv(sig_lipid_features, "sig_lipid_features.csv", row.names = FALSE)

# 2) Row order: sort by descending Log2_FC (most upregulated on top)
lip_order_sig <- os_ys_stats %>%
  dplyr::filter(Lipid %in% signature59) %>%
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

mat_osys_z <- scale(mat_osys)

ann_col_osys <- data.frame(Group = heatmap_data_osys$Group)
rownames(ann_col_osys) <- heatmap_data_osys$Sample
ann_cols <- list(Group = c(YS = "#89CFF0", OS = "#4F4F4F"))

pal <- colorRampPalette(c("#FF7F00","black","#00FFFF"))(100)

hm_3b <- pheatmap(
  t(mat_osys_z),
  annotation_col    = ann_col_osys,
  annotation_colors = ann_cols,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row   = lip_order_sig, row_names_side = "left",
  fontsize = 6, color = pal, breaks = seq(-1, 1, length.out = 101),
  main = "59-lipid Aging Signature: YS (left, blue) vs OS (right, grey)"
)
print(hm_3b)

pdf("FigX_heatmap_lipid_aging_59.pdf", width = 6.75, height = 6.6)
print(hm_3b)
dev.off()

# 4) Horizontal bar plot — absolute logFC, colored by aging direction
bar3b_df <- os_ys_stats %>%
  dplyr::filter(Lipid %in% signature59) %>%
  dplyr::mutate(
    abs_logFC = abs(Log2_FC_OS_vs_YS),
    bar_fill  = ifelse(Log2_FC_OS_vs_YS > 0, "#8B0000", "#00008B"),
    Lipid     = factor(Lipid, levels = lip_order_sig)
  )

bar_3b <- ggplot(bar3b_df, aes(x = abs_logFC, y = fct_rev(Lipid), fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Absolute log2 Fold Change (Old Sed vs Young Sed)",
    y = NULL
  ) +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(bar_3b)

pdf("FigX_barplot_lipid_aging_59.pdf", width = 6.75, height = 6.6)
print(bar_3b)
dev.off()
while (!is.null(dev.list())) dev.off()
###############################################################################
###############################################################################
#///////Alternate Lipid Heatmap/////////#

hm_3b_alt <- pheatmap(
  t(mat_osys_z),
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  legend          = FALSE,
  annotation_col  = NA,
  labels_row      = lip_order_sig,
  fontsize_row    = 5,
  color = pal, breaks = seq(-1, 1, length.out = 101)
)
print(hm_3b_alt)

pdf("FigX_heatmap_lipid_aging_59_alt.pdf", width = 4.4, height = 5.9)
print(hm_3b_alt)
dev.off()
#====================================================================#

Fig_aging_bar <- ggplot(aging_bar_data, aes(x = abs_logFC, y = Symbol, fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(x = "Absolute log2 Fold Change (Old vs Young)",
       y = NULL) +
  theme(
    axis.text.y   = element_blank(),
    axis.text.x   = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor    = element_blank()
  )


lipid_bar_3b <- ggplot(bar3b_df, aes(x = abs_logFC, y = fct_rev(Lipid), fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Absolute log2 Fold Change (Old Sed vs Young Sed)",
    y = NULL
  ) +
  theme_minimal(base_size = 6) +
  theme(
    axis.text.y   = element_blank(),
    axis.text.x   = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor    = element_blank()
  )


pdf("FigX_barplot_aging_axis_56Lipids.pdf", width = 4.2, height = 5.9)
print(lipid_bar_3b)
dev.off()
###########################################################################
##############################################################################
#/////Lipid Sig Formatting//////////#

##############################################################################
# Export for LipidSig 2.0
##############################################################################

# --- 1) Clean lipid names for rgoslin compatibility ---
# LipidSig uses rgoslin to parse names, which follows shorthand notation
lipid_name_map <- tibble(
  original = lipid_cols
) %>%
  mutate(
    clean = original,
    # Remove (a\b), (a\b\c), [sn1], [SIM], [+OH], trailing (104)
    clean = str_remove_all(clean, "\\s*\\(a\\\\b(\\\\c)?\\)"),
    clean = str_remove_all(clean, "\\s*\\[sn\\d\\]"),
    clean = str_remove_all(clean, "\\s*\\[SIM\\]"),
    clean = str_remove_all(clean, "\\s*\\[\\+OH\\]"),
    clean = str_remove_all(clean, "\\s*\\(\\d+\\)$"),
    clean = str_trim(clean),
    # "LPC 20:3" -> "LPC 20:3" (rgoslin uses space, not parens for some)
    # But most tools want "LPC(20:3)" — keep parens format
    clean = str_replace(clean, "^(LPC|LPE|PC|PE)\\s+(\\d)", "\\1 \\2")
  )

# Check for duplicates after cleaning
dupes <- lipid_name_map %>% count(clean) %>% filter(n > 1)
if (nrow(dupes) > 0) {
  cat("Warning: duplicate names after cleaning:\n")
  print(dupes)
  # Append suffix to duplicates
  lipid_name_map <- lipid_name_map %>%
    group_by(clean) %>%
    mutate(clean = ifelse(n() > 1, paste0(clean, "_", row_number()), clean)) %>%
    ungroup()
}

# Preview the mapping
lipid_name_map %>% filter(original != clean) %>% print(n = 30)

# --- 2) Abundance table (feature + sample columns, RAW values, not log2) ---
# LipidSig wants raw abundance, it handles normalization/transform internally
lip_abundance <- lip_raw %>%
  filter(Sample != "YS6") %>%                    # match your analysis
  filter(Group %in% c("OS", "YS")) %>%           # or include all groups
  select(all_of(lipid_cols)) %>%
  t() %>%
  as.data.frame()

# Set column names to sample names
sample_names <- lip_raw %>%
  filter(Sample != "YS6") %>%
  filter(Group %in% c("OS", "YS")) %>%
  pull(Sample)
colnames(lip_abundance) <- sample_names

# Add feature column with cleaned names
lip_abundance <- lip_abundance %>%
  tibble::rownames_to_column("original") %>%
  left_join(lipid_name_map, by = "original") %>%
  mutate(feature = clean) %>%
  select(feature, all_of(sample_names))

# Remove COH and Ubiquinone
lip_abundance <- lip_abundance %>%
  filter(!feature %in% c("COH", "Ubiquinone"))

lip_abundance <- lip_abundance %>% filter(!feature %in% c("PC(38:6)_1", "PC(38:6)_2"))

# --- 3) Group info table ---
lip_group_info <- lip_raw %>%
  filter(Sample != "YS6") %>%
  filter(Group %in% c("OS", "YS")) %>%
  transmute(
    sample_name = Sample,
    label_name  = Sample,
    group       = Group,
    pair        = NA
  )

write_csv(lip_group_info, "lipidsig_group_info_OS_vs_YS.csv")

# --- 4) Export ---
write_csv(lip_abundance, "lipidsig_abundance_OS_vs_YS.csv")
write_csv(lip_group_info, "lipidsig_group_info_OS_vs_YS.csv")

cat("\n=== Files exported ===\n")
cat("Abundance:", nrow(lip_abundance), "lipids x", ncol(lip_abundance) - 1, "samples\n")
cat("Groups:", nrow(lip_group_info), "samples\n")
cat("\nPreview abundance:\n")
print(lip_abundance[1:5, 1:5])
cat("\nPreview groups:\n")
print(lip_group_info)

library(rgoslin)

test_parse <- rgoslin::parseLipidNames(lip_abundance$feature)

# Check which ones failed
failed <- test_parse %>% filter(Grammar == "NOT_PARSEABLE")
cat("Failed to parse:", nrow(failed), "out of", nrow(test_parse), "\n")
print(failed$Original.Name)
#############################################################################
###############################################################################

##############################################################################
# Lipid Class Enrichment — Fisher's exact test (ORA approach)
# Matches LipidSig methodology: tests over-representation per class
##############################################################################

# Uses os_ys_stats which already has LipidClass and IsSignature defined

# --- 1) Class-level ORA (Fisher's exact, one-sided) ---
class_enrichment <- os_ys_stats %>%
  filter(!is.na(LipidClass)) %>%
  group_by(LipidClass) %>%
  summarize(
    n_total     = n(),
    n_sig       = sum(IsSignature),
    n_sig_up    = sum(IsSignature & Log2_FC_OS_vs_YS > 0),
    n_sig_down  = sum(IsSignature & Log2_FC_OS_vs_YS < 0),
    mean_logFC  = mean(Log2_FC_OS_vs_YS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_lipids = sum(n_total),
    total_sig    = sum(n_sig),
    # Fisher's exact: is this class enriched among signature lipids?
    p_value = purrr::map2_dbl(n_sig, n_total, function(k, n) {
      mat <- matrix(c(
        k,                          # sig in this class
        n - k,                      # not sig in this class
        total_sig[1] - k,           # sig in other classes
        total_lipids[1] - n - total_sig[1] + k  # not sig in other classes
      ), nrow = 2)
      fisher.test(mat, alternative = "greater")$p.value
    }),
    FDR = p.adjust(p_value, method = "fdr"),
    pct_sig = round(100 * n_sig / n_total, 1),
    pct_expected = round(100 * total_sig / total_lipids, 1),
    fold_enrichment = round((n_sig / n_total) / (total_sig / total_lipids), 2)
  ) %>%
  arrange(p_value)

cat("=== Class-level enrichment (Fisher's exact, BH-corrected) ===\n")
print(class_enrichment %>%
        select(LipidClass, n_total, n_sig, n_sig_up, n_sig_down,
               pct_sig, fold_enrichment, p_value, FDR), n = 20)

# --- 2) Class-level t-test (aggregated abundance, like LipidSig) ---
# This tests whether the MEAN abundance of each class differs between OS and YS
class_abundance_test <- os_ys_stats %>%
  filter(!is.na(LipidClass)) %>%
  group_by(LipidClass) %>%
  summarize(
    n_species = n(),
    # For each class, get the mean log2FC across all species
    mean_logFC = mean(Log2_FC_OS_vs_YS, na.rm = TRUE),
    # Aggregate p-value: one-sample t-test on logFC values (H0: mean = 0)
    class_pval = tryCatch(
      t.test(Log2_FC_OS_vs_YS, mu = 0)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    class_FDR = p.adjust(class_pval, method = "fdr"),
    direction = ifelse(mean_logFC > 0, "Up", "Down")
  ) %>%
  arrange(class_pval)

cat("\n=== Class-level abundance test (one-sample t-test on logFC) ===\n")
print(class_abundance_test, n = 20)

# --- 3) Visualization: class enrichment barplot ---
enrich_plot_df <- class_enrichment %>%
  mutate(
    neg_log10_p = -log10(p_value),
    LipidClass = factor(LipidClass, levels = rev(LipidClass))  # ordered by p-value
  )

class_enrichment_plot <- ggplot(enrich_plot_df, aes(x = neg_log10_p, y = LipidClass)) +
  geom_col(aes(fill = fold_enrichment), color = "black", width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_gradient(low = "grey80", high = "#8B0000", name = "Fold\nEnrichment") +
  labs(
    x = expression(-log[10](P)),
    y = "Lipid Class",
    title = "Lipid Class Enrichment in Aging Signature",
    subtitle = "Fisher's exact test (one-sided) | Red line = p = 0.05"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

print(class_enrichment_plot)

# --- 4) Combined summary table ---
cat("\n=== Combined summary ===\n")
combined <- class_enrichment %>%
  select(LipidClass, n_total, n_sig, pct_sig, fold_enrichment,
         Fisher_p = p_value, Fisher_FDR = FDR) %>%
  left_join(
    class_abundance_test %>% select(LipidClass, mean_logFC, class_pval, class_FDR, direction),
    by = "LipidClass"
  ) %>%
  arrange(Fisher_p)

print(combined, n = 20)
##############################################################################

library(tidyverse)

# sample-level class abundance as SUM of raw abundances within class,
# then log2-transform the class total
class_sample_sums <- lip_imp %>%
  filter(Sample != "YS6", Group %in% c("OS", "YS")) %>%
  pivot_longer(cols = all_of(lipid_cols), names_to = "Lipid", values_to = "abund_raw") %>%
  left_join(
    os_ys_stats %>% dplyr::select(Lipid, LipidClass),
    by = "Lipid"
  ) %>%
  filter(!is.na(LipidClass)) %>%
  group_by(Sample, Group, LipidClass) %>%
  summarize(class_sum = sum(abund_raw, na.rm = TRUE), .groups = "drop") %>%
  mutate(log2_class_sum = log2(class_sum))

class_abundance_test2 <- class_sample_sums %>%
  group_by(LipidClass) %>%
  summarize(
    n_OS = sum(Group == "OS"),
    n_YS = sum(Group == "YS"),
    mean_OS = mean(log2_class_sum[Group == "OS"], na.rm = TRUE),
    mean_YS = mean(log2_class_sum[Group == "YS"], na.rm = TRUE),
    class_log2FC = mean_OS - mean_YS,
    class_pval = tryCatch(
      t.test(log2_class_sum ~ Group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    class_FDR = p.adjust(class_pval, method = "fdr"),
    direction = ifelse(class_log2FC > 0, "Up", "Down")
  ) %>%
  arrange(class_pval)

print(class_abundance_test2, n = 20)
##############################################################

# Class-level abundance barplot (no Sph)
class_bar_df <- class_abundance_test2 %>%
  filter(LipidClass != "Sph") %>%
  mutate(
    bar_fill = case_when(
      LipidClass == "LPC" ~ "#8B0000",
      LipidClass == "CE"  ~ "#CC4400",
      TRUE                ~ "grey70"
    ),
    LipidClass = factor(LipidClass, levels = rev(LipidClass))
  )

class_bar_plot <- ggplot(class_bar_df, aes(x = class_log2FC, y = LipidClass, fill = bar_fill)) +
  geom_col(color = "black", width = 0.7, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_identity() +
  labs(
    x = "log2 Fold Change (Old Sed − Young Sed)",
    y = "Lipid Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

print(class_bar_plot)
##############################################################


library(ggridges)

# --- Define ORA-significant classes ---
ora_sig_classes <- class_enrichment %>%
  filter(FDR < 0.05, !LipidClass %in% c("Sph", "PA", "PG")) %>%
  pull(LipidClass)

# --- Assign colors ---
sig_colors <- c("#CC4400", "#8B0000", "#66C2A5", "#984EA3", "#FF7F00")
ora_color_map <- setNames(sig_colors[seq_along(ora_sig_classes)], ora_sig_classes)

# --- ORA-based ordering ---
ora_order <- class_enrichment %>%
  filter(!LipidClass %in% c("Sph", "PA", "PG")) %>%
  arrange(p_value) %>%
  pull(LipidClass)

# --- Clip threshold ---
x_clip <- 2.0

# --- Build ridge dataframe ---
ridge_df <- os_ys_stats %>%
  filter(!is.na(LipidClass), !LipidClass %in% c("Sph", "PA", "PG")) %>%
  group_by(LipidClass) %>%
  mutate(
    n_class     = n(),
    n_sig_class = sum(IsSignature)
  ) %>%
  ungroup() %>%
  mutate(
    x_display = pmin(pmax(Log2_FC_OS_vs_YS, -x_clip), x_clip),
    dot_group = case_when(
      IsSignature & Log2_FC_OS_vs_YS > 0 ~ "Sig Up",
      IsSignature & Log2_FC_OS_vs_YS < 0 ~ "Sig Down",
      TRUE ~ "NS"
    ),
    ridge_group = ifelse(LipidClass %in% ora_sig_classes, LipidClass, "Other"),
    LipidClassLabel = paste0(LipidClass, " (", n_sig_class, "/", n_class, ")")
  )

# --- Build factor levels for y-axis ordering ---
label_order <- ridge_df %>%
  distinct(LipidClass, LipidClassLabel) %>%
  mutate(LipidClass = factor(LipidClass, levels = ora_order)) %>%
  arrange(LipidClass) %>%
  pull(LipidClassLabel)

ridge_df <- ridge_df %>%
  mutate(LipidClassLabel = factor(LipidClassLabel, levels = rev(label_order)))

# --- Build fill palette ---
ridge_fill_vals <- c(ora_color_map, "Other" = "grey85")

# --- Plot ---
class_ridge <- ggplot(ridge_df, aes(x = x_display, y = LipidClassLabel)) +
  geom_density_ridges(
    aes(fill = ridge_group),
    alpha = 0.7,
    scale = 0.9,
    rel_min_height = 0.01,
    color = "black",
    linewidth = 0.3
  ) +
  geom_jitter(
    data = subset(ridge_df, dot_group == "NS"),
    color = "grey60",
    height = 0.15, width = 0, size = 1.5, alpha = 0.6
  ) +
  geom_jitter(
    data = subset(ridge_df, dot_group == "Sig Up"),
    color = "#8B0000",
    height = 0.15, width = 0, size = 3, alpha = 0.9
  ) +
  geom_jitter(
    data = subset(ridge_df, dot_group == "Sig Down"),
    color = "#00008B",
    height = 0.15, width = 0, size = 3, alpha = 0.9
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = ridge_fill_vals, guide = "none") +
  coord_cartesian(xlim = c(-1.2, x_clip + 0.2)) +
  labs(
    x = "log2 Fold Change (Old Sed − Young Sed)",
    y = "Lipid Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(class_ridge)

# Save as PDF
ggsave(
  "Fig1d_oldsedveh_lipid_class_ridge.pdf",
  class_ridge,
  width = 9.9,
  height = 6.9
)

###############################################################

# --- 1) Already have OS vs YS class results in class_abundance_test2 ---

# --- 2) Run the same analysis for OV vs OS ---
class_sample_sums_ov <- lip_imp %>%
  filter(Group %in% c("OV", "OS")) %>%
  pivot_longer(cols = all_of(lipid_cols), names_to = "Lipid", values_to = "abund_raw") %>%
  left_join(
    os_ys_stats %>% dplyr::select(Lipid, LipidClass),
    by = "Lipid"
  ) %>%
  filter(!is.na(LipidClass)) %>%
  group_by(Sample, Group, LipidClass) %>%
  summarize(class_sum = sum(abund_raw, na.rm = TRUE), .groups = "drop") %>%
  mutate(log2_class_sum = log2(class_sum))

class_abundance_ov <- class_sample_sums_ov %>%
  group_by(LipidClass) %>%
  summarize(
    mean_OV = mean(log2_class_sum[Group == "OV"], na.rm = TRUE),
    mean_OS = mean(log2_class_sum[Group == "OS"], na.rm = TRUE),
    class_log2FC = mean_OV - mean_OS,
    class_pval = tryCatch(
      t.test(log2_class_sum ~ Group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(class_FDR = p.adjust(class_pval, method = "fdr"))

# --- 3) Combine into one dataframe ---
combined_class <- bind_rows(
  class_abundance_test2 %>%
    filter(LipidClass != "Sph") %>%
    select(LipidClass, class_log2FC, class_pval, class_FDR) %>%
    mutate(comparison = "Old Sed vs Young Sed"),
  class_abundance_ov %>%
    filter(LipidClass != "Sph") %>%
    select(LipidClass, class_log2FC, class_pval, class_FDR) %>%
    mutate(comparison = "Old PWR vs Old Sed")
)

# Order classes by aging logFC (most changed on top)
class_order_bar <- class_abundance_test2 %>%
  filter(LipidClass != "Sph") %>%
  arrange(class_log2FC) %>%
  pull(LipidClass)

combined_class <- combined_class %>%
  mutate(
    LipidClass = factor(LipidClass, levels = class_order_bar),
    comparison = factor(comparison, levels = c("Old Sed vs Young Sed", "Old PWR vs Old Sed"))
  )

# --- 4) Paired barplot ---
class_paired_bar <- ggplot(combined_class,
                           aes(x = class_log2FC, y = LipidClass, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65,
           color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c("Old Sed vs Young Sed" = "#7F7F7F",
               "Old PWR vs Old Sed"   = "#4DAF4A"),
    name = NULL
  ) +
  labs(
    x = "log2 Fold Change",
    y = "Lipid Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(face = "bold"),
    legend.position    = "top"
  )

print(class_paired_bar)
##########################################################################
#///Lipid species analysis////#

library(tidyverse)

# -----------------------------
# 1) Build structural annotation table from lipid names
# -----------------------------

lipid_annotation <- tibble(Lipid = lipid_cols) %>%
  mutate(
    LipidClass = str_trim(str_extract(Lipid, "^[^\\(]+")),
    LipidClass = case_when(
      grepl("^LPC", Lipid) ~ "LPC",
      grepl("^HexCer|^Hex2Cer", Lipid) ~ "HexCer",
      LipidClass %in% c("COH", "Ubiquinone") ~ NA_character_,
      TRUE ~ LipidClass
    ),
    
    # pull the text inside parentheses if present
    inside_paren = str_match(Lipid, "\\(([^\\)]+)\\)")[,2],
    
    # fallback for names like "LPC 20:3"
    inside_paren = ifelse(
      is.na(inside_paren),
      str_match(Lipid, "\\s(\\d+:\\d+)$")[,2],
      inside_paren
    ),
    
    inside_paren = str_replace_all(inside_paren, "/", "_"),
    
    # extract all n:n patterns
    acyl_matches = stringr::str_extract_all(inside_paren, "\\d+:\\d+"),
    
    n_acyls = purrr::map_int(acyl_matches, length),
    
    total_carbons = purrr::map_dbl(acyl_matches, ~{
      if (length(.x) == 0) return(NA_real_)
      sum(as.numeric(sub(":.*", "", .x)))
    }),
    
    total_double_bonds = purrr::map_dbl(acyl_matches, ~{
      if (length(.x) == 0) return(NA_real_)
      sum(as.numeric(sub(".*:", "", .x)))
    }),
    
    mean_chain_length = ifelse(n_acyls > 0, total_carbons / n_acyls, NA_real_),
    mean_double_bonds = ifelse(n_acyls > 0, total_double_bonds / n_acyls, NA_real_),
    
    # simple unsaturation / PUFA bins
    unsat_bin = case_when(
      is.na(total_double_bonds) ~ NA_character_,
      total_double_bonds == 0 ~ "Saturated",
      total_double_bonds == 1 ~ "Monounsaturated",
      total_double_bonds >= 2 ~ "Polyunsaturated"
    ),
    
    db_bin = case_when(
      is.na(total_double_bonds) ~ NA_character_,
      total_double_bonds <= 1 ~ "0-1 DB",
      total_double_bonds <= 3 ~ "2-3 DB",
      total_double_bonds <= 5 ~ "4-5 DB",
      total_double_bonds >= 6 ~ "6+ DB"
    ),
    
    chain_bin = case_when(
      is.na(total_carbons) ~ NA_character_,
      total_carbons < 34 ~ "<34C",
      total_carbons < 40 ~ "34-39C",
      total_carbons < 46 ~ "40-45C",
      total_carbons >= 46 ~ "46+C"
    ),
    
    # crude PUFA flag
    is_pufa = total_double_bonds >= 2
  ) %>%
  dplyr::select(
    Lipid, LipidClass, n_acyls,
    total_carbons, total_double_bonds,
    mean_chain_length, mean_double_bonds,
    chain_bin, db_bin, unsat_bin, is_pufa
  )

print(lipid_annotation, n = 100)

lipid_annotation %>%
  filter(is.na(total_carbons) | is.na(total_double_bonds)) %>%
  print(n = 50)
#--------------------------------------#

os_ys_struct <- os_ys_stats %>%
  left_join(lipid_annotation, by = c("Lipid", "LipidClass"))

print(os_ys_struct, n = 20)

struct_summary_by_class <- os_ys_struct %>%
  filter(LipidClass %in% c("CE", "LPC", "PC", "PE", "TG")) %>%
  group_by(LipidClass) %>%
  summarize(
    n_species = n(),
    mean_log2FC = mean(Log2_FC_OS_vs_YS, na.rm = TRUE),
    mean_total_carbons = mean(total_carbons, na.rm = TRUE),
    mean_total_double_bonds = mean(total_double_bonds, na.rm = TRUE),
    mean_db_among_up = mean(total_double_bonds[Log2_FC_OS_vs_YS > 0], na.rm = TRUE),
    mean_db_among_down = mean(total_double_bonds[Log2_FC_OS_vs_YS < 0], na.rm = TRUE),
    .groups = "drop"
  )

print(struct_summary_by_class, n = 20)
#-----------------------------------------------#

struct_tests <- os_ys_struct %>%
  filter(LipidClass %in% c("CE", "LPC", "PC", "PE", "TG")) %>%
  mutate(direction = ifelse(Log2_FC_OS_vs_YS > 0, "Up_in_OS", "Down_in_OS")) %>%
  group_by(LipidClass) %>%
  summarize(
    n_up = sum(direction == "Up_in_OS"),
    n_down = sum(direction == "Down_in_OS"),
    
    mean_c_up = mean(total_carbons[direction == "Up_in_OS"], na.rm = TRUE),
    mean_c_down = mean(total_carbons[direction == "Down_in_OS"], na.rm = TRUE),
    p_chain = tryCatch(
      t.test(total_carbons ~ direction)$p.value,
      error = function(e) NA_real_
    ),
    
    mean_db_up = mean(total_double_bonds[direction == "Up_in_OS"], na.rm = TRUE),
    mean_db_down = mean(total_double_bonds[direction == "Down_in_OS"], na.rm = TRUE),
    p_db = tryCatch(
      t.test(total_double_bonds ~ direction)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_chain_fdr = p.adjust(p_chain, method = "fdr"),
    p_db_fdr = p.adjust(p_db, method = "fdr")
  )

print(struct_tests, n = 20)
#-----------------------------------------------#

pufa_enrichment <- os_ys_struct %>%
  filter(LipidClass %in% c("CE", "LPC", "PC", "PE", "TG")) %>%
  mutate(
    sig_up = IsSignature & Log2_FC_OS_vs_YS > 0
  ) %>%
  group_by(LipidClass) %>%
  summarize(
    n_total = n(),
    n_pufa = sum(is_pufa, na.rm = TRUE),
    n_sig_up = sum(sig_up, na.rm = TRUE),
    n_sig_up_pufa = sum(sig_up & is_pufa, na.rm = TRUE),
    fisher_p = tryCatch(
      fisher.test(matrix(c(
        n_sig_up_pufa,
        n_pufa - n_sig_up_pufa,
        n_sig_up - n_sig_up_pufa,
        n_total - n_pufa - n_sig_up + n_sig_up_pufa
      ), nrow = 2), alternative = "greater")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(fisher_p, method = "fdr"))

print(pufa_enrichment, n = 20)
#---------------------------------------------------#

lip_long_struct <- lip_imp %>%
  filter(Sample != "YS6", Group %in% c("OS", "YS")) %>%
  pivot_longer(cols = all_of(lipid_cols), names_to = "Lipid", values_to = "abund_raw") %>%
  left_join(lipid_annotation, by = "Lipid") %>%
  filter(!is.na(LipidClass))

# example: sample-level class x DB-bin totals
sample_struct_sums <- lip_long_struct %>%
  filter(LipidClass %in% c("CE", "LPC", "PC", "PE", "TG"),
         !is.na(db_bin)) %>%
  group_by(Sample, Group, LipidClass, db_bin) %>%
  summarize(bin_sum = sum(abund_raw, na.rm = TRUE), .groups = "drop") %>%
  mutate(log2_bin_sum = log2(bin_sum))

print(sample_struct_sums, n = 20)
#-------------------------------------------------#

struct_bin_tests <- sample_struct_sums %>%
  group_by(LipidClass, db_bin) %>%
  summarize(
    mean_OS = mean(log2_bin_sum[Group == "OS"], na.rm = TRUE),
    mean_YS = mean(log2_bin_sum[Group == "YS"], na.rm = TRUE),
    log2FC = mean_OS - mean_YS,
    p_value = tryCatch(
      t.test(log2_bin_sum ~ Group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)

print(struct_bin_tests, n = 50)
############################################################
##################################################################



###############################################
#
#
#
#Metabolomic Analysis
#
#
#
##############################################################################
###############################################################################
#///////////////METABOLOMIC ANALYSES////////////////#
############################################################################
##############################################################################
# METABOLOMICS — Old Sed VEH vs Yng Sed VEH
# Variable names prefixed with metab_ to avoid lipid collisions
##############################################################################
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(forcats)

# ---------- 0) Load and impute ----------
metab_raw <- readr::read_csv("Konopka_Muscle_HILIC.csv", show_col_types = FALSE)
metab_cols <- setdiff(names(metab_raw), c("Sample","Group"))

# Per-metabolite min nonzero imputation
metab_min_nonzero <- sapply(metab_raw[metab_cols], function(x) {
  m <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (is.infinite(m)) NA_real_ else m
})
metab_global_min <- suppressWarnings(
  min(unlist(metab_raw[metab_cols])[unlist(metab_raw[metab_cols]) > 0], na.rm = TRUE)
)
if (!is.finite(metab_global_min)) metab_global_min <- 1e-6
metab_min_nonzero[is.na(metab_min_nonzero)] <- metab_global_min

metab_imp <- metab_raw
for (met in metab_cols) {
  x <- metab_imp[[met]]
  x[is.na(x)] <- 0
  if (all(x == 0)) {
    x[x == 0] <- metab_global_min / 2
  } else {
    x[x == 0] <- metab_min_nonzero[[met]] / 2
  }
  metab_imp[[met]] <- x
}

# Log2 transform
metab_log <- metab_imp %>% dplyr::mutate(across(all_of(metab_cols), log2))

# ---------- 1) Stats: OS vs YS ----------
metab_osys <- metab_log %>% dplyr::filter(Group %in% c("OS","YS"))

metab_osys_stats <- purrr::map_dfr(metab_cols, function(met) {
  vals <- metab_osys[[met]]; grp <- metab_osys$Group
  tt <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
  tibble(
    Metabolite       = met,
    Mean_YS          = mean(vals[grp == "YS"], na.rm = TRUE),
    Mean_OS          = mean(vals[grp == "OS"], na.rm = TRUE),
    Log2_FC_OS_vs_YS = mean(vals[grp == "OS"], na.rm = TRUE) - mean(vals[grp == "YS"], na.rm = TRUE),
    P_value          = if (is.null(tt)) NA_real_ else tt$p.value
  )
}) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  dplyr::arrange(FDR)

# ---------- 2) Define metabolite aging signature (P < 0.05) ----------
metab_osys_stats <- metab_osys_stats %>%
  dplyr::mutate(IsSignature = P_value < 0.05)

metab_osys_stats
write.csv(metab_osys_stats, "metab_osys_stats.csv", row.names = FALSE)

metab_age_signature <- metab_osys_stats %>%
  dplyr::filter(IsSignature) %>%
  dplyr::pull(Metabolite)

cat("Metabolite aging signature:", length(metab_age_signature), "features\n")
metab_age_signature

# ---------- 3) Volcano plot ----------
metab_sig_df <- metab_osys_stats %>% dplyr::filter(IsSignature)
metab_sig_df

write.csv(metab_sig_df, "metab_sig_df.csv", row.names = FALSE)
metab_xlim <- max(1, quantile(abs(metab_osys_stats$Log2_FC_OS_vs_YS), 0.99, na.rm = TRUE))

metab_volcano <- ggplot(metab_osys_stats, aes(x = Log2_FC_OS_vs_YS, y = -log10(P_value))) +
  geom_point(color = "grey65", size = 2) +
  geom_point(data = metab_sig_df, shape = 21, fill = "#8B0000", color = "black", stroke = 0.6, size = 5) +
  ggrepel::geom_text_repel(
    data = metab_sig_df %>% slice_min(P_value, n = 5),
    aes(label = Metabolite),
    size = 3, max.overlaps = 100,
    box.padding = 0.4, point.padding = 0.3, seed = 1
  )  +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(-metab_xlim, metab_xlim)) +
  labs(
    x = "log2 Fold Change (OLD SED − YNG SED)",
    y = expression(-log[10](P)),
    title = "Metabolome: Old Sed VEH vs Yng Sed VEH"
  ) +
  theme_minimal(base_size = 11)

print(metab_volcano)

ggsave(
  "FigX_oldsedveh_metabo_volcano_plot_v1.pdf",
  metab_volcano,
  width = 9.9,
  height = 6.9
)



# ---------- 4) Row order for heatmap/barplot ----------
metab_row_order <- metab_osys_stats %>%
  dplyr::filter(Metabolite %in% metab_age_signature) %>%
  dplyr::arrange(dplyr::desc(Log2_FC_OS_vs_YS)) %>%
  dplyr::pull(Metabolite)

# ---------- 5) Heatmap: YS vs OS ----------
metab_pal <- colorRampPalette(c("#1F77B4", "white", "#D62728"))(100)

metab_hm_data <- metab_log %>%
  dplyr::filter(Group %in% c("YS","OS")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("YS","OS"))) %>%
  dplyr::arrange(Group) %>%
  dplyr::select(Sample, Group, dplyr::all_of(metab_row_order))

metab_mat <- metab_hm_data %>% dplyr::select(-Sample, -Group) %>% as.matrix()
rownames(metab_mat) <- metab_hm_data$Sample
metab_mat_z <- scale(metab_mat)

metab_ann_col <- data.frame(Group = metab_hm_data$Group)
rownames(metab_ann_col) <- metab_hm_data$Sample
metab_ann_colors <- list(Group = c(YS = "#89CFF0", OS = "#4F4F4F"))

metab_hm <- pheatmap::pheatmap(
  t(metab_mat_z),
  annotation_col    = metab_ann_col,
  annotation_colors = metab_ann_colors,
  cluster_rows = FALSE, cluster_cols = FALSE,
  labels_row   = metab_row_order, row_names_side = "left",
  fontsize = 8, color = metab_pal, breaks = seq(-1, 1, length.out = 101),
  main = "Metabolite Aging Signature: YS (left) vs OS (right)"
)
print(metab_hm)

metab_hm_alt <- pheatmap(
  t(metab_mat_z),
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
print(metab_hm_alt)

pdf("FigX_heatmap_metabolite_aging.pdf", width = 4.4, height = 2.8)
print(metab_hm_alt)
dev.off()

# ---------- 6) Barplot: absolute logFC, colored by aging direction ----------
metab_bar_data <- metab_osys_stats %>%
  dplyr::filter(Metabolite %in% metab_age_signature) %>%
  dplyr::mutate(
    abs_logFC = abs(Log2_FC_OS_vs_YS),
    bar_fill  = ifelse(Log2_FC_OS_vs_YS > 0, "#8B0000", "#00008B"),
    Metabolite = factor(Metabolite, levels = metab_row_order)
  )

metab_bar <- ggplot(metab_bar_data, aes(x = abs_logFC, y = fct_rev(Metabolite), fill = bar_fill)) +
  geom_col() +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Absolute log2 Fold Change (Old Sed vs Young Sed)",
    y = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(metab_bar)

pdf("FigX_barplot_metabolite_aging.pdf", width = 6.75, height = 3)
print(metab_bar)
dev.off()
############################################################
############################################################
#
#Metabolite Analysis
#
#===========================================================================#
# Metabolite Class Directed logFC — parallels lipid class analysis
# Classification based on HMDB chemical taxonomy (Wishart et al., NAR 2022)
#===========================================================================#

library(dplyr); library(stringr); library(ggplot2)

# --- 1) Classify all metabolites by HMDB-style categories ---
# Read the full stats to get all metabolite names
metab_osys_stats <- read.csv("metab_osys_stats.csv")

all_mets <- metab_osys_stats$Metabolite

# Classification function based on HMDB chemical taxonomy
classify_metabolite <- function(met) {
  m <- tolower(met)
  case_when(
    # Amino acids (standard + direct derivatives)
    str_detect(m, "^l-|^dl-|^d-alanine|^b-alanine|^glycine$|^homoserine|^norleucine|sarcosine|^l-norvaline|trans-4-hydroxy-l-proline") &
      str_detect(m, "arginin|glutam|proline|serine|threonine|valine|leucine|isoleucine|histidin|tryptoph|tyrosin|phenylalanin|aspara|aspart|cystath|ornithine|pipecol|alanine|glycine|homoserine|norleucine|norvaline|hydroxy-l-proline|canavanine|allothreonine|hydroxylysine|aminoadipic|sarcosine") ~ "Amino acids",
    str_detect(m, "argininosuccinic|adma|aminobutyric|aminocyclopropane|sulfinoalanine|guanidinosuccinic|citrulline|guanidineacetic|trimethyllysine") ~ "Amino acids",
    # Additional amino acids / derivatives caught from Other
    str_detect(m, "5-hydroxylysine|methionine sulfoxide|l-cystine|cysteic acid|cysteamine|homocysteine|o-phospho-l-serine|nitrotyrosine|4-hydroxy-l-glutamic|ketoisovaleric|prephenic|carnosine|4-methylaminobutyrate|methylguanidine|n-carbamoyl.*aspart|n-carbamyl.*glutam|s-2-aminoethyl-l-cysteine|s-\\(carboxymethyl\\)-l-cysteine|meso-2,6-diaminoheptanedioate|na-acetyl-l-asparagine") ~ "Amino acids",
    
    # Modified amino acids (N-acetyl, methyl, formyl)
    str_detect(m, "n-acetyl.*alanine|n-acetyl.*serine|n-acetyl.*glutam|n-acetyl.*phenyl|n-acetyl.*tryptoph|n-acetyl.*aspara|n-acetylaspart|n-acetylproline|n-formyl|n-methyl.*alanine|n-methyl.*glutam|n-methyl.*aspart|acetyl-l-leucine|acetyl-l-cysteine|methyl.*histidine|methylhistamine|pyroglutamic|acetyl-l-leucine|n-acetylputrescine|aceturic acid") ~ "Modified amino acids",
    
    # Nucleotides
    str_detect(m, "amp |adp |atp|gmp|gdp|gtp|ump|udp|utp|imp|idp|itp|damp|dadp|datp|dgmp|dgdp|dgtp|dcmp|dctp|dtdp|dttp|camp|5-amp|5-phospho-d-ribose|diquafosol|dutp|2,3 cyclic cmp|aicar") ~ "Nucleotides",
    
    # Nucleosides & bases
    str_detect(m, "adenosine|guanosine|inosine|uridine|cytidine|thymidine|2-deoxyadenosine|2-deoxyguanosine|2-deoxyuridine|2-deoxycytidine|5-methyluridine|adenine|guanine|xanthine|hypoxanthine|thymine|uracil|purine|xanthosine|methyladenine|methylcytosine|8-hydroxy-2-deoxyguanosine|uracil-5-carbox|allantoin|orotic acid|lumazine|pterine|dihydrouracil|3-ureidopropionic|b-ureidoisobutyric|dihydroorotic") ~ "Nucleosides & bases",
    
    # Carbohydrates
    str_detect(m, "glucose|fructose|maltose|lactose|cellobiose|galactose|mannose|melibiose|isomaltulose|ribose(?!.*phosphate)|xylose|rhamnose|glucuronic|galacturonic|gluconic|arabinose|tagatose|psicose|glucosamine|galactosamine|sedoheptulose(?!.*phosph)|deoxy-d-ribose|gluconolactone|fucose|n-acetylneuraminic|glucosaminic|methyl-b-d-galactopyranoside") ~ "Carbohydrates",
    
    # Sugar phosphates
    str_detect(m, "glucose.*phosph|fructose.*phosph|mannose.*phosph|ribose.*phosph|ribulose.*phosph|xylulose.*phosph|erythrose.*phosph|glucosamine.*phosph|sedoheptulose.*phosph|deoxy-d-glucose.*phosph|deoxyribose.*phosph|galactose.*phosph|dihydroxyacetone phosph|glycerate.*phosph|arabinose.*phosph|carbamoyl phosphate|dimethylallyl diphosphate") ~ "Sugar phosphates",
    
    # Sugar alcohols / polyols
    str_detect(m, "sorbitol|mannitol|adonitol|dulcitol|arabitol|threitol|inositol|pinitol|glycerol$|glyceraldehyde|l-glyceric") ~ "Sugar alcohols",
    
    # TCA cycle & energy
    str_detect(m, "citric|isocitrate|cis-aconit|trans-aconit|alpha-ketoglutaric|succinic acid|fumaric|oxal(?:ic|oacet)|pyruvic|lactic acid|phosphoenolpyruvic|glyoxylic|citramalic|citraconic|malic acid|ketoglutaric|l-hydroxyglutaric") ~ "TCA cycle intermediates",
    
    # Acyl-CoAs & CoA derivatives
    str_detect(m, "coa$|coa |coenzyme a$|acetoacetyl|acetyl-coa|butyryl|propionyl|malonyl|isobutyryl|isovaleryl|succinyl") ~ "Acyl-CoAs & cofactors",
    
    # Vitamins & redox cofactors
    str_detect(m, "nad$|nadh|nadp|fad |biotin|thiamine|pyridox|folic|folinic|pantothenic|pantolactone|nicotinamide|nicotinic|cobalamin|ascorbic|lipoic|coenzyme q10|lumichrome|mevalonic|dethiobiotin|plp |dihydro-l-biopterin|5\\(d\\)-aminolevulinic") ~ "Vitamins & cofactors",
    
    # Glutathione & redox
    str_detect(m, "glutathione|gsh |ophthalmic acid|s-hexylglutathione") ~ "Glutathione / redox",
    
    # Biogenic amines & neurotransmitters
    str_detect(m, "dopamine|tryptamine|tyramine|histamine|serotonin|norepinephrine|melatonin|octopamine|gaba |acetylcholine|choline$|betaine$|normetanephrine") ~ "Biogenic amines",
    
    # Tryptophan / Phenylalanine metabolism
    str_detect(m, "kynuren|xanthurenic|3-hydroxyanthranilic|indole|5-hydroxy.*tryptoph|tryptophanamide|3-o-methyl-l-dopa|hydroxyphenyl|homovanillic|homogentisic|urocanic|3-\\(3-indolyl\\)") ~ "Trp/Phe metabolism",
    
    # Carnitines
    str_detect(m, "carnitine") ~ "Carnitines",
    
    # Bile acids & steroids
    str_detect(m, "cholic|taurocholic|glycocholic|estradiol|estrone|estriol|progesterone|hydroxyprogesterone|corticosterone|cortisone|hydrocortisone|bilirubin|biliverdin|3,5-diiodo-l-thyronine") ~ "Bile acids & steroids",
    
    # Phospholipid precursors
    str_detect(m, "phosphocholine|phosphorylethanolamine|cdp-ethanolamine|glycerophosphocholine|glycerol.*phosph") ~ "Phospholipid precursors",
    
    # Creatine metabolism
    str_detect(m, "creatine|creatinine") ~ "Creatine metabolism",
    
    # Sulfur metabolism
    str_detect(m, "taurine|hypotaurine|s-\\(5-adenosyl\\)-l-methionine|sah s-adenosyl|5-deoxy-5-\\(methylthio\\)") ~ "Sulfur metabolism",
    
    # Urea cycle / nitrogen
    str_detect(m, "urea$|uric acid") ~ "Urea cycle",
    
    # Organic acids
    str_detect(m, "butyric|heptanoic|hexanoic|adipic|pimelic|mandelic|hydroxybutyric|tartaric|methylglutaric|dimethylsuccinic|methyl glutarate|hydroxyphenylacetic|hydroxybenzoic|hydroxybenzaldehyde|benzoic|cinnamic|cinnamaldehyde|methylacetoacetate|glycolic|malonic|2-isopropylmalic|2-keto|2-oxoadipic|succinic semi|phosphonoacetate|5-valerolactone|3r-hydroxy-isobutyric|\\(r\\)-2,3-dihydroxy|2-hydroxypyridine") ~ "Organic acids",
    
    # Polyphenols
    str_detect(m, "epicatechin") ~ "Polyphenols",
    
    # Catch remaining
    TRUE ~ "Other"
  )
}

metab_osys_stats$MetabClass <- sapply(metab_osys_stats$Metabolite, classify_metabolite)

# Check
as.data.frame(metab_osys_stats %>% count(MetabClass) %>% arrange(desc(n)))

# What's still in Other?
as.data.frame(metab_osys_stats %>%
                filter(MetabClass == "Other") %>%
                dplyr::select(Metabolite))

metab_osys_stats$MetabClass <- sapply(metab_osys_stats$Metabolite, classify_metabolite)

# Check classification
cat("\n=== Metabolite classes ===\n")
metab_osys_stats %>% count(MetabClass) %>% arrange(desc(n)) %>% print(n = 25)

# Check which signature metabolites fall where
cat("\n=== Signature metabolites by class ===\n")
metab_osys_stats %>%
  filter(IsSignature) %>%
  dplyr::select(Metabolite, MetabClass) %>%
  arrange(MetabClass) %>%
  as.data.frame() %>%
  print(n = 24)

# Class counts
as.data.frame(metab_osys_stats %>% count(MetabClass) %>% arrange(desc(n)))

# Signature metabolites by class
as.data.frame(metab_osys_stats %>%
                filter(IsSignature == TRUE) %>%
                dplyr::select(Metabolite, MetabClass) %>%
                arrange(MetabClass))

as.data.frame(metab_osys_stats %>%
                filter(MetabClass == "Other") %>%
                dplyr::select(Metabolite) %>%
                arrange(Metabolite))
#######################################################



