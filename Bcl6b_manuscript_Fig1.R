##Code for Bcl6b Manuscript Figure 1.
#Figure 1. Bcl6b predicted to regulate an exercise-responsive angiogenic transcriptional program.
#This is the analysis of the young sed vs pwr data.
#DEG identify Bcl6b as top TF from ChEA3 analysis. 
#DEG are EC enriched via Human Protein Atlas analysis.

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
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(httr)
library(jsonlite)
library(tidyr)
library(stringr)

#####################################################
#####################################################
#-------Figure 1a-----------------------------------#
#Identify yng PWR DEG and read in ChEA3 results----#
#####################################################

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
#-----------------------------------------------------------------#
# 1) Load edgeR *_GENE_* xlsx and build ranked vector
yngpwrveh_v_yngsedveh_genes_16dec <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_YNG_PWR_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "YNG_PWR_VEH-YNG_SED_VEH_logFC",
  FDR_col   = "YNG_PWR_VEH-YNG_SED_VEH_FDR"
)

#remove duplicate ENTREZID and keep most significant
yngpwrveh_v_yngsedveh_genes_16dec <- yngpwrveh_v_yngsedveh_genes_16dec %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_min(FDR, n = 1) %>%
  dplyr::ungroup()

#----Get all pwr DEG-------------------------#
all_yngpwrdeg_chea3_genes <- yngpwrveh_v_yngsedveh_genes_16dec %>%
  filter(FDR < 0.05 & abs(logFC) >0.3) %>%
  pull(Symbol)

length(all_yngpwrdeg_chea3_genes) #495 Genes/427 ChEA3 valid
write.csv(all_yngpwrdeg_chea3_genes, "all_yngpwrdeg_chea3_genes.csv",row.names = FALSE)
#------------------------------------------------------------------------#

#####################################################
#####################################################
#-------Figure 1b-----------------------------------#
#Identify yng PWR DEG and read in ChEA3 results----#
#####################################################

#---------------------------------------------------------#
#-----Get upregulated pwr DEG-----------------------------#
upreg_yngpwrdeg_chea3_genes <- yngpwrveh_v_yngsedveh_genes_16dec %>%
  filter(FDR < 0.05 & logFC > 0.3) %>%
  pull(Symbol)

length(upreg_yngpwrdeg_chea3_genes) #357 Genes/ 307 ChEA3 valid/ 50 ovrlp w/ Bcl6b
write.csv(upreg_yngpwrdeg_chea3_genes, "upreg_yngpwrdeg_chea3_genes.csv",row.names = FALSE)
#---------------------------------------------------------#

#---------------------------------------------------------#
#-----Get downregulated pwr DEG-----------------------------#
downreg_yngpwrdeg_chea3_genes <- yngpwrveh_v_yngsedveh_genes_16dec %>%
  filter(FDR < 0.05 & logFC < -0.3) %>%
  pull(Symbol)

length(downreg_yngpwrdeg_chea3_genes) #138 Genes/ 120 ChEA3 valid/ 13 ovrlp w/Bcl6b
write.csv(downreg_yngpwrdeg_chea3_genes, "downreg_yngpwrdeg_chea3_genes.csv",row.names = FALSE)
#13 downregulated overlapping genes:
#THY1,FOS,KLF4,APOLD1,CXCL10,CCND1,ADAMTS1,MYC,SLCO2A1,PDK4,SERPINH1,TXNIP,ANGPTL4
#-----------------------------------------------------------#

#-----Create volcano plot w/highlighted overlapping genes---#

head(yngpwrveh_v_yngsedveh_genes_16dec)
dim(yngpwrveh_v_yngsedveh_genes_16dec)
setwd("/Users/mdbruss/Documents/RStudioProjects_2/BCL6B")

upreg_cheas_ovrlp_genes <- read_csv("upreg_yngpwrdeg_chea3_ovrlp_genes.csv")
upreg_cheas_ovrlp_genes <- scan("upreg_yngpwrdeg_chea3_ovrlp_genes.csv", what = character(), sep = ",")
upreg_cheas_ovrlp_genes
upreg_cheas_ovrlp_genes <- upreg_cheas_ovrlp_genes %>%
  tolower() %>%
  str_to_title()

upreg_cheas_ovrlp_genes

downreg_cheas_ovrlp_genes <- c(
  "Thy1",
  "Fos",
  "Klf4",
  "Apold1",
  "Cxcl10",
  "Ccnd1",
  "Adamts1",
  "Myc",
  "Slco2a1",
  "Pdk4",
  "Serpinh1",
  "Txnip",
  "Angptl4"
)
downreg_cheas_ovrlp_genes
#------------------------------------#

library(dplyr)
library(ggplot2)

# 1) Build a volcano dataframe -------------------------------------------
volcano_df <- yngpwrveh_v_yngsedveh_genes_16dec %>%
  dplyr::mutate(
    neg_log10_FDR = -log10(FDR),
    
    sig_status = dplyr::case_when(
      FDR < 0.05 & logFC >  0.3  ~ "Up",
      FDR < 0.05 & logFC < -0.3  ~ "Down",
      TRUE                       ~ "NS"
    ),
    
    highlight = dplyr::case_when(
      Symbol %in% upreg_cheas_ovrlp_genes   ~ "Up_CHEA",
      Symbol %in% downreg_cheas_ovrlp_genes ~ "Down_CHEA",
      TRUE                                  ~ "None"
    )
  )

# 2) Base volcano plot ----------------------------------------------------
library(ggrepel)

genes_to_label <- c("Dll4", "Cdh5", "Flt1", "Efnb2", "Hey1", "Rasip1", "Pecam1")

chea3_ovrlp_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_FDR)) +
  
  # background points
  geom_point(aes(color = sig_status),
             alpha = 0.4, size = 1.2) +
  
  # highlighted CHEA3 overlap genes
  geom_point(
    data = dplyr::filter(volcano_df, highlight != "None"),
    shape = 21,
    size  = 2.0,
    stroke = 0.8,
    fill  = "red",
    color = "black"
  ) +
  
  # cutoff lines
  geom_vline(xintercept = c(-0.3, 0.3),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  
  # label selected genes
  geom_text_repel(
    data = dplyr::filter(volcano_df, Symbol %in% genes_to_label),
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
  
  labs(
    x = "log2 Fold Change (Yng PoWeR Veh vs Yng Sed Veh)",
    y = "-log10(FDR)",
    color = "Status",
    title = "Volcano plot: Yng PoWeR Veh vs Yng Sed Veh",
    subtitle = "Highlighted: ChEA3-overlap genes + labeled angiogenic hits"
  ) +
  
  coord_cartesian(ylim = c(0, 10)) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.title       = element_text(face = "bold")
  )
# Save as PDF
pdf("chea3_ovrlp_volcano.pdf", width = 8, height = 6)  # Adjust height
print(chea3_ovrlp_volcano)
dev.off()
#----------------------------------------------------------------------#
#----------------------------------------------------------------------#

#####################################################
#####################################################
#-------Figure 1c-----------------------------------#
#------Determine Endothelial Cell Enrichment----#
#####################################################
setwd("/Users/mdbruss/Documents/RStudioProjects_2/Rapa_PwR")

# 1) Load HPA single-cell table
sc_data <- read.delim("rna_single_cell_type.tsv", check.names = FALSE)
head(sc_data)

# Rename to match the rest of your pipeline
sc_data <- dplyr::rename(
  sc_data,
  Gene.name = `Gene name`,
  Cell.type = `Cell type`
)
sc_data

# Ranking df. Each gene's top 5 cell express types
top5_celltypes_df <- sc_data %>%
  group_by(Gene.name) %>%
  arrange(desc(nTPM), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  mutate(rank = row_number()) %>%
  pivot_wider(
    names_from = rank,
    values_from = c(Cell.type, nTPM),
    names_glue = "#{rank}_{.value}"
  ) %>%
  ungroup()

n_distinct(sc_data$Gene.name)
nrow(top5_celltypes_df)

head(top5_celltypes_df)

top5_celltypes_df %>%
  filter(Gene.name %in% upreg_cheas_ovrlp_genes)

top5_celltypes_df <- top5_celltypes_df %>%
  rowwise() %>%
  mutate(
    isEC = as.integer(
      any(c_across(ends_with("Cell.type")) == "Endothelial cells")
    )
  ) %>%
  ungroup()

table(top5_celltypes_df$isEC)

top5_celltypes_df %>%
  dplyr::filter(Gene.name %in% upreg_cheas_ovrlp_genes) %>%
  dplyr::count(isEC)

upreg_yngpwrdeg_chea3_genes <- read_csv("upreg_yngpwrdeg_chea3_genes.csv")

upreg_yngpwrdeg_chea3_genes_vec <- upreg_yngpwrdeg_chea3_genes %>%
  dplyr::pull(x)

upper_upreg_yngpwrdeg_chea3_genes_vec <- toupper(upreg_yngpwrdeg_chea3_genes_vec)
length(upper_upreg_yngpwrdeg_chea3_genes_vec)

top5_celltypes_df %>%
  dplyr::filter(Gene.name %in% upper_upreg_yngpwrdeg_chea3_genes_vec) %>%
  dplyr::count(isEC)

top5_celltypes_df %>%
  filter(Gene.name %in% upper_upreg_yngpwrdeg_chea3_genes_vec)


celltype_rank_df <- top5_celltypes_df %>%
  dplyr::filter(Gene.name %in% upper_upreg_yngpwrdeg_chea3_genes_vec) %>%
  tidyr::pivot_longer(
    cols = ends_with("Cell.type"),
    names_to = "rank_position",
    values_to = "Cell.type"
  ) %>%
  dplyr::count(Cell.type, name = "n") %>%
  dplyr::arrange(desc(n))

celltype_rank_df

top_celltypes <- celltype_rank_df %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(Cell.type)

top_celltypes


celltype_cols <- grep("Cell\\.type$", names(top5_celltypes_df), value = TRUE)

for (ct in top_celltypes) {
  new_col <- paste0("is_", gsub(" ", "_", ct))
  
  top5_celltypes_df[[new_col]] <- as.integer(
    apply(top5_celltypes_df[, celltype_cols], 1, function(x) any(x == ct, na.rm = TRUE))
  )
}

# Optional safety: convert any remaining NA flags to 0
top5_celltypes_df <- top5_celltypes_df %>%
  dplyr::mutate(dplyr::across(starts_with("is_"), ~ tidyr::replace_na(.x, 0L)))


top5_celltypes_df %>%
  dplyr::select(starts_with("is_")) %>%
  colSums()

dim(top5_celltypes_df)
#------------------------------------------#

#Subset to overlapping genes
#First just test to see what it would look like.
bcl6b_ovrlp_celltypes_df <- top5_celltypes_df %>%
  filter(Gene.name %in% upper_upreg_yngpwrdeg_chea3_genes_vec) %>%
  filter(Gene.name %in% upreg_cheas_ovrlp_genes)

bcl6b_ovrlp_celltypes_df %>%
  dplyr::select(starts_with("is_")) %>%
  colSums()

upreg_cheas_ovrlp_genes
length(upreg_cheas_ovrlp_genes)

#Subset to up with exercise, but not overlapping genes.

upreg_yngpwrdeg_not_ovrlp_celltypes_df <- top5_celltypes_df %>%
  filter(Gene.name %in% upper_upreg_yngpwrdeg_chea3_genes_vec) %>%
  filter(!Gene.name %in% upreg_cheas_ovrlp_genes)

upreg_yngpwrdeg_not_ovrlp_celltypes_df %>%
  dplyr::select(starts_with("is_")) %>%
  colSums()

upreg_yngpwrdeg_notovrlp_genes <- upreg_yngpwrdeg_not_ovrlp_celltypes_df %>%
  pull(Gene.name)

upreg_yngpwrdeg_notovrlp_genes
length(upreg_yngpwrdeg_notovrlp_genes)

top5_celltypes_upreg_yngpwrdeg_df <- top5_celltypes_df %>%
  filter(Gene.name %in% upper_upreg_yngpwrdeg_chea3_genes_vec)

top5_celltypes_upreg_yngpwrdeg_df %>%
  filter(Gene.name %in% upreg_cheas_ovrlp_genes) %>%
  dplyr::select(starts_with("is_")) %>%
  colSums()

top5_celltypes_upreg_yngpwrdeg_df %>%
  filter(!Gene.name %in% upreg_cheas_ovrlp_genes) %>%
  dplyr::select(starts_with("is_")) %>%
  colSums()
#------------------------------------------------#

#---Create side-by-side chart------#

plot_df <- top5_celltypes_upreg_yngpwrdeg_df %>%
  mutate(Group = ifelse(Gene.name %in% upreg_cheas_ovrlp_genes, "Bcl6b_overlap", "Not_overlap")) %>%
  group_by(Group) %>%
  summarise(
    n_genes = n(),
    across(starts_with("is_"), ~ mean(.x, na.rm = TRUE)),   # fraction 0–1
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("is_"),
    names_to = "CellTypeFlag",
    values_to = "Fraction"
  ) %>%
  mutate(
    Percent = 100 * Fraction,
    Cell.type = CellTypeFlag %>%
      str_remove("^is_") %>%
      str_replace_all("_", " ")
  )

plot_df <- plot_df %>%
  mutate(
    Cell.type = factor(Cell.type, levels = gsub("_", " ", gsub(" ", "_", top_celltypes)))
  )


ggplot(plot_df, aes(x = Cell.type, y = Percent, fill = Group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  coord_flip() +
  labs(
    x = NULL,
    y = "% of genes with cell type in top-5",
    title = "Top-5 cell-type representation: Bcl6b-overlap vs non-overlap"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c(
    "Bcl6b_overlap" = "blue",
    "Not_overlap"  = "red"
  )) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%"))
#-----------------------------------------------#
#-----------------------------------------------#

#============================================================
# 1) Build Gene x Cell.type matrix (nTPM)
#============================================================
mat_df <- sc_data %>%
  dplyr::filter(Gene.name %in% upreg_cheas_ovrlp_genes) %>%
  dplyr::select(Gene.name, Cell.type, nTPM) %>%
  tidyr::pivot_wider(
    names_from  = Cell.type,
    values_from = nTPM,
    values_fill = 0
  )

mat <- mat_df %>%
  tibble::column_to_rownames("Gene.name") %>%
  as.matrix()

#============================================================
# 2) Transform + row Z-score
#============================================================
mat_log <- log2(mat + 1)
mat_z   <- t(scale(t(mat_log)))
mat_z[is.na(mat_z)] <- 0

#============================================================
# 3) Cluster rows and columns
#============================================================
row_hc <- hclust(dist(mat_z), method = "complete")
col_hc <- hclust(dist(t(mat_z)), method = "complete")

row_d <- as.dendrogram(row_hc)
col_d <- as.dendrogram(col_hc)

#============================================================
# 4) Colors (blue -> white -> red)
#============================================================
hm_colors <- colorRampPalette(c("blue", "white", "red"))(101)

#============================================================
# 5) Heatmap WITH dendrograms, but DO NOT reorder dendrogram branches
#============================================================
op <- par(no.readonly = TRUE)
par(mar = c(14, 10, 3, 2))

heatmap(
  mat_z,
  Rowv = row_d,                    # left dendrogram
  Colv = col_d,                    # top dendrogram
  reorderfun = function(d, w) d,   # <-- critical: stop heatmap() from reordering
  scale = "none",
  col = hm_colors,
  labCol = colnames(mat_z),
  cexRow = 0.7,
  cexCol = 0.55,
  margins = c(14, 10)
)

par(op)
#---------------------------------------------------#
#===================================================#


#####################################################
#####################################################
#-------Figure 1d-----------------------------------#
#------Bcl6b EC overlap compared to other TFs----#
#####################################################

all_yngpwrdeg_chea3_ranklist <- read_tsv('all_yngpwrdeg_chea3_ranklist.tsv')

head(all_yngpwrdeg_chea3_ranklist)
top5_celltypes_upreg_yngpwrdeg_df
upper_upreg_yngpwrdeg_chea3_genes_vec

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

#============================================================
# SETTINGS YOU CAN TUNE
#============================================================
OR_CAP <- 50        # <-- (2) tune this: max OR shown in point sizing
TOP_LABEL_N <- 5    # <-- (3) number of TFs to label by %EC

#============================================================
# Helper: parse comma-separated overlap genes from ChEA3
#============================================================
parse_overlap_genes <- function(x) {
  if (is.na(x) || length(x) == 0) return(character(0))
  x %>%
    str_split(",") %>%
    .[[1]] %>%
    str_trim() %>%
    toupper() %>%
    unique()
}

#============================================================
# 1) Define "universe" (genes that went into ChEA3)
#============================================================
universe_genes <- toupper(unique(upper_upreg_yngpwrdeg_chea3_genes_vec))

#============================================================
# 2) Build gene -> isEndothelial lookup within the universe
#    (Using your existing binary column: is_Endothelial_cells)
#============================================================
gene_ec_map <- top5_celltypes_upreg_yngpwrdeg_df %>%
  dplyr::transmute(
    Gene.name = toupper(Gene.name),
    isEC = as.integer(is_Endothelial_cells > 0)
  ) %>%
  dplyr::distinct(Gene.name, .keep_all = TRUE)

# Ensure every universe gene has a value (missing => 0)
gene_ec_map <- tibble(Gene.name = universe_genes) %>%
  dplyr::left_join(gene_ec_map, by = "Gene.name") %>%
  dplyr::mutate(isEC = ifelse(is.na(isEC), 0L, isEC))

isEC_vec <- gene_ec_map$isEC
names(isEC_vec) <- gene_ec_map$Gene.name

#============================================================
# 3) Top 50 TFs from ChEA3 rank list, compute enrichment stats
#============================================================
top50_tf <- all_yngpwrdeg_chea3_ranklist %>%
  dplyr::arrange(Rank) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::select(Rank, TF, Overlapping_Genes)

res_tf50 <- top50_tf %>%
  rowwise() %>%
  mutate(
    overlap_genes = list(parse_overlap_genes(Overlapping_Genes)),
    overlap_genes = list(intersect(overlap_genes, universe_genes)),
    n_overlap = length(overlap_genes),
    
    # counts in overlap
    a_ec    = sum(isEC_vec[overlap_genes] == 1, na.rm = TRUE),
    b_nonec = sum(isEC_vec[overlap_genes] == 0, na.rm = TRUE),
    
    # "other" = universe minus overlap
    other_genes = list(setdiff(universe_genes, overlap_genes)),
    c_ec    = sum(isEC_vec[other_genes] == 1, na.rm = TRUE),
    d_nonec = sum(isEC_vec[other_genes] == 0, na.rm = TRUE),
    
    pct_EC_overlap = ifelse(n_overlap > 0, 100 * a_ec / n_overlap, NA_real_),
    
    fisher = list(fisher.test(matrix(c(a_ec, b_nonec, c_ec, d_nonec), nrow = 2))),
    OR     = unname(fisher$estimate),
    CI_low = fisher$conf.int[1],
    CI_high= fisher$conf.int[2],
    p_value = fisher$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    neglog_padj = -log10(p_adj + 1e-300)   # (1) color by -log10(p_adj)
  )

#============================================================
# 4) Prep for plotting:
#    - cap OR at OR_CAP for size scale
#    - handle Inf OR by setting to OR_CAP
#============================================================
plot_df <- res_tf50 %>%
  mutate(
    OR_plot = dplyr::case_when(
      is.na(OR) ~ NA_real_,
      is.infinite(OR) ~ OR_CAP,
      TRUE ~ as.numeric(OR)
    ),
    OR_plot = pmin(OR_plot, OR_CAP)
  )

# (3) Label top N TFs by %EC in overlap (+ always label BCL6B + selected TFs)
top_to_label <- plot_df %>%
  arrange(desc(pct_EC_overlap), desc(n_overlap)) %>%
  slice_head(n = TOP_LABEL_N) %>%
  pull(TF)

force_label_tfs <- c("BCL6B", "EPAS1", "ERG", "EGR1")

label_set <- union(top_to_label, force_label_tfs)

plot_df <- plot_df %>%
  mutate(label = ifelse(TF %in% label_set, TF, NA_character_))

#============================================================
# 5) Scatter plot
#============================================================
ggplot(plot_df, aes(x = n_overlap, y = pct_EC_overlap)) +
  geom_point(aes(size = OR_plot, color = neglog_padj), alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = label),
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Number of overlapping genes (ChEA3 overlap)",
    y = "% of overlap genes with Endothelial cells in top-5 (HPA scRNA)",
    size  = paste0("Odds ratio (capped at ", OR_CAP, ")"),
    color = expression(-log[10](p[adj])),
    title = "Endothelial enrichment across top 50 ChEA3 TFs",
    subtitle = "Point size = OR of EC enrichment (overlap vs non-overlap within ChEA3 input universe)"
  )

top50_tf
res_tf50 %>% 
  filter(n_overlap >80)

res_tf50 %>% 
  filter(pct_EC_overlap >35)

res_tf50 %>% 
  filter(TF == "MEF2")
