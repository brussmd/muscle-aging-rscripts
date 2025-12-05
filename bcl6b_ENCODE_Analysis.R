# Step 1: Get ENCODE TF target data for BCL6B from Harmonizome (if available)
# Harmonizome has a dataset: "ENCODE Transcription Factor Targets"

# Install packages if needed
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

# Load them
library(httr)
library(jsonlite)
library(tibble)
library(dplyr)


# Define parameters
base_url <- "https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene/"
tf_name <- "BCL6B"
dataset <- "ENCODE+Transcription+Factor+Targets"
url <- paste0(base_url, tf_name, "/dataset/", dataset)

# API call
response <- GET(url)

# Check and parse
if (status_code(response) == 200) {
  tf_targets <- fromJSON(content(response, "text"))
  
  # Check if associations exist
  if (!is.null(tf_targets$associations) && length(tf_targets$associations) > 0) {
    # Print column names for debugging
    print(colnames(tf_targets$associations))
    
    # Try to format output
    bcl6b_targets <- tf_targets$associations %>%
      as_tibble() %>%
      dplyr::select(
        targetGene = `gene.symbol`, 
        score = associationScore
      ) %>%
      arrange(desc(score))
    
    print(bcl6b_targets)
  } else {
    message("No associations found for BCL6B in ENCODE TF dataset.")
  }
} else {
  message("Failed to query Harmonizome API.")
}

#------------------------------------------------#

download.file(
  url = "https://www.encodeproject.org/files/ENCFF165BAG/@@download/ENCFF165BAG.bed",
  destfile = "BCL6B_optimal_IDR_peaks_GRCh38.bed"
)


# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)

# 📦 Load necessary packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
library(dplyr)

# 📥 Import the BCL6B peak file
bcl6b_peaks <- rtracklayer::import("ENCFF165BAG.bed", format = "narrowPeak")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnno <- annotatePeak(
  bcl6b_peaks,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),  # Define promoter region
  annoDb = "org.Hs.eg.db"
)

peak_df <- as.data.frame(peakAnno)

head(peak_df)

# Unique symbols of genes with a nearby BCL6B binding site
bcl6b_target_genes <- unique(peak_df$SYMBOL)

# Optional: inspect first few
head(bcl6b_target_genes)
length(bcl6b_target_genes)

promoter_targets <- peak_df %>%
  filter(grepl("Promoter", annotation)) %>%
  pull(SYMBOL) %>%
  unique()

# Check number
length(promoter_targets)
head(promoter_targets)

write.csv(promoter_targets, "BCL6B_ENCODE_Promoter_Targets.csv",row.names = FALSE)
library(clusterProfiler)
#-----------------------------------------------------#
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

peaks <- rtracklayer::import("ENCFF165BAG.bed", format = "narrowPeak")

# 1) High-confidence peaks only (example thresholds; tune to your file)
m <- mcols(bcl6b_peaks)
peaks_hq <- bcl6b_peaks[
  !is.na(m$qValue) & m$qValue >= 2 &        # q <= 0.01 if ENCODE -log10
    !is.na(m$signalValue) & m$signalValue > 5 # bump as needed
]

# 2) Strict promoter mapping (±1 kb)
anno_prom <- annotatePeak(
  peaks_hq, TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Hs.eg.db"
)
prom_df <- as.data.frame(anno_prom) %>%
  dplyr::filter(grepl("Promoter", annotation)) %>%
  dplyr::filter(!is.na(SYMBOL))

promoter_targets <- unique(prom_df$SYMBOL)
length(promoter_targets); head(promoter_targets)

write.csv(promoter_targets, "Narrowed_BCL6B_ENCODE_Promoter_Targets.csv",row.names = FALSE)
#-----------------------------------------------#

library(biomaRt)
library(dplyr)

# ---- Lock BioMart to the Dec 2021 Ensembl archive (stable mappings) ----
host <- "https://dec2021.archive.ensembl.org"

# ✅ FIX: connect with dataset specified in useMart (no piping into useDataset)
# Human (GRCh38 / Ensembl 105 in Dec 2021)
ensembl_hs <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host    = host
)

# Mouse (GRCm38 / Ensembl 105 in Dec 2021)
ensembl_mm <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host    = host
)

# (Optional) sanity check:
# listMarts(host = host)
# listDatasets(ensembl_hs) %>% dplyr::select(dataset, version) %>% head()

# =========================
# HUMAN: map symbols -> Ensembl + biotype, then make filtered lists
# Input: 'promoter_targets' is your character vector of human SYMBOLs
# =========================
symbols_hs <- unique(promoter_targets)
attrs_hs <- c("ensembl_gene_id", "external_gene_name", "gene_biotype")

map_hs <- getBM(
  attributes = attrs_hs,
  filters    = "external_gene_name",
  values     = symbols_hs,
  mart       = ensembl_hs
) %>%
  dplyr::distinct() %>%
  dplyr::select(
    ENSEMBL     = ensembl_gene_id,
    SYMBOL      = external_gene_name,
    GENEBIOTYPE = gene_biotype
  )

# Protein-coding only (recommended for ChEA3/LISA)
targets_pc <- map_hs %>%
  dplyr::filter(GENEBIOTYPE == "protein_coding") %>%
  dplyr::pull(SYMBOL) %>%
  unique()

# Extended: protein-coding + selected lncRNA classes; drop generic LOC*
targets_extended <- map_hs %>%
  dplyr::filter(
    GENEBIOTYPE == "protein_coding" |
      GENEBIOTYPE %in% c("lncRNA", "antisense", "sense_intronic", "sense_overlapping")
  ) %>%
  dplyr::pull(SYMBOL) %>%
  unique() %>%
  (\(x) x[!grepl("^LOC", x)])()

# Quick peek
length(targets_pc); length(targets_extended); head(targets_pc)

write.csv(targets_pc, "BCL6B_ENCODE_narrowed_protein_coding_Promoter_Targets.csv",row.names = FALSE)

targets_pc

library(biomaRt)
library(dplyr)

host  <- "https://dec2021.archive.ensembl.org"
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)

qs <- sort(unique(na.omit(targets_pc)))

# Step 1: get Ensembl IDs for your human symbols
h_ids <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
               filters    = "hgnc_symbol",
               values     = qs,
               mart       = human) %>%
  dplyr::distinct()

# Step 2: pull mouse homologs + orthology type using the human Ensembl IDs
h2m <- getBM(attributes = c("ensembl_gene_id",
                            "mmusculus_homolog_ensembl_gene",
                            "mmusculus_homolog_associated_gene_name",
                            "mmusculus_homolog_orthology_type"),
             filters    = "ensembl_gene_id",
             values     = h_ids$ensembl_gene_id,
             mart       = human) %>%
  dplyr::distinct()

orth_map <- h_ids %>%
  dplyr::left_join(h2m, by = "ensembl_gene_id") %>%
  dplyr::rename(human_symbol = hgnc_symbol,
                mouse_symbol = mmusculus_homolog_associated_gene_name,
                mouse_ensembl= mmusculus_homolog_ensembl_gene,
                orthology_type = mmusculus_homolog_orthology_type) %>%
  dplyr::filter(mouse_symbol != "") %>%
  dplyr::group_by(human_symbol) %>%
  dplyr::arrange(dplyr::desc(orthology_type == "ortholog_one2one")) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

mouse_targets_1to1 <- unique(orth_map$mouse_symbol)

write.csv(orth_map, "human_to_mouse_orthology_map.csv", row.names = FALSE)
write.csv(mouse_targets_1to1, "mouse_targets_1to1_vector.csv", row.names = FALSE)

length(mouse_targets_1to1); head(mouse_targets_1to1)


mouse_targets_1to1
targets_pc

# 1) uppercase both sets
m_upper <- toupper(mouse_targets_1to1)
h_upper <- toupper(targets_pc)

# 2) exact-name intersection (after uppercasing)
same_upper <- intersect(m_upper, h_upper)

# 3) how many / which ones?
length(same_upper)        # count
sort(same_upper)          # list (uppercased)



#------------------------------------#
ego_promoter <- enrichGO(
  gene = promoter_targets,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

head(ego_promoter, n=30)

head(promoter_targets)

# Load or install biomaRt

library(biomaRt)

# Connect to archived Ensembl (Dec 2021) for stable ortholog mapping
human_mart <- useEnsembl(biomart = "genes",
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://dec2021.archive.ensembl.org")

mouse_mart <- useEnsembl(biomart = "genes",
                         dataset = "mmusculus_gene_ensembl",
                         host = "https://dec2021.archive.ensembl.org")

# Replace this with your vector of human gene symbols
human_genes <- promoter_targets

# Perform ortholog mapping
conversion <- getLDS(attributes = c("hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = human_genes,
                     mart = human_mart,
                     attributesL = c("mgi_symbol"),
                     martL = mouse_mart)

# Rename columns
colnames(conversion) <- c("Human", "Mouse")

# Get unique mouse gene symbols
mouse_genes <- unique(conversion$Mouse)

length(mouse_genes)
head(mouse_genes)

potential_bcl6b_genes.vec <- c(
  "Cdkn1a","Spaar","Flt1","Dock9","Tnfaip2","C1qtnf9","Cdh5","Rgs5","Gja1","Rgs3",
  "Hey1","Ccnd1","Sox18","Adamts1","Vsig2","Myc","Pdk4","Sox7","Pdgfrb","Sema6b",
  "Plekhg1","Egfl7","Rsad2","Ralgapa2","Sh2d3c","Fos","Inmt","Tanc1","Adgrf5","Slco2a1",
  "Pecam1","Txnip","Heca","Angptl4","St6galnac3","Rapgef5","Mocs1","Kank3","Ushbp1","Ccm2l",
  "Lpl","Thy1","Gata2","Fbln2","Apold1","Rasip1","Efnb2","Dll4","Fgd5","Ablim3",
  "Serpinh1","Stom","N4bp3","Gpihbp1","Btnl9","Smad1","Nr1h3","Fmo2","Klf4","Fli1",
  "Bmp6","Efna1","Tjp1","Cxcl10","Smtn","Cxcl12","Tspan18","Pxmp2","Esam"
)

# Overlap genes
overlap_genes <- intersect(mouse_genes, potential_bcl6b_genes.vec)

# How many overlap?
length(overlap_genes)

# Which genes overlap?
overlap_genes

pwr_genes <- yngpwrveh_v_yngsedveh_genes %>%
  filter(FDR < 0.05) %>%
  pull (Symbol)

yngpwrveh_v_yngsedveh_genes
upreg_pwr_genes <- yngpwrveh_v_yngsedveh_genes %>%
  filter(FDR < 0.05 & logFC >0) %>%
  pull (Symbol)

length(pwr_genes)
length(upreg_pwr_genes)

pwr_bcl6b_ovrlp <- intersect(mouse_genes, pwr_genes)
upreg_pwr_bcl6b_ovrlp <- intersect(mouse_genes, upreg_pwr_genes)

length(pwr_bcl6b_ovrlp)
length(upreg_pwr_bcl6b_ovrlp)
pwr_bcl6b_ovrlp

pwr_bcl6b_ovrlp_genes <- as.character(pwr_bcl6b_ovrlp)
ego <- enrichGO(
  gene          = pwr_bcl6b_ovrlp_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",   # because your genes are in gene symbol format
  ont           = "BP",       # Biological Process (can also use "MF" or "CC")
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Quick summary
head(ego, n=20)

dim(yngpwrveh_v_yngsedveh_genes)

tot_yng_genes <- yngpwrveh_v_yngsedveh_genes %>%
  pull(Symbol)

length(tot_yng_genes)
length(mouse_genes)
1872-1571

bcl6b_tot_pool <- intersect(mouse_genes, tot_yng_genes)
length(bcl6b_tot_pool)

valid_pool <- 

degs_tot_pool <- intersect(bcl6b_tot_pool, pwr_genes)
length(degs_tot_pool)


total_genes <- 11318  # or the number of expressed genes in your analysis
bcl6b_targets <- 1571
degs <- 621
overlap <- 70

# Build contingency table
contingency <- matrix(c(overlap,
                        bcl6b_targets - overlap,
                        degs - overlap,
                        total_genes - bcl6b_targets - (degs - overlap)),
                      nrow = 2,
                      byrow = TRUE)

colnames(contingency) <- c("DEG", "Not_DEG")
rownames(contingency) <- c("BCL6B_target", "Not_BCL6B_target")

contingency

# Fisher’s Exact Test
fisher.test(contingency)
#--------------------------------------------#

#============================================#
#----Cell Type Expression of Bcl6b Promoter genes-----#
#=====================================================#

## ---- Setup ----
# Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(stringr)
})

## ---- Inputs you already have ----
# 1) HPA single-cell "cell type" matrix (downloaded from HPA link)
#    Must contain columns: Gene.name, Cell.type, nTPM
sc_data <- read.delim("rna_single_cell_type.tsv", check.names = FALSE)

head(sc_data)
# Rename to safe R-friendly names
colnames(sc_data) <- c("EnsemblID", "GeneName", "CellType", "nTPM")

# 2) Your BCL6B promoter target list (character vector of human symbols)
#    Example: head(promoter_targets) shows "B3GALNT2" "CIPC" ...
target_genes_raw <- promoter_targets

## ---- Clean & harmonize targets to HPA gene naming ----
# HPA Gene.name is uppercase gene symbols; enforce upper-case for matching
target_genes <- unique(na.omit(toupper(target_genes_raw)))

# Keep only targets that exist in HPA table
all_hpa_genes <- unique(toupper(sc_data$GeneName))
target_genes_in_hpa <- intersect(target_genes, all_hpa_genes)

all_hpa_genes

message("Targets provided: ", length(target_genes),
        " | Found in HPA: ", length(target_genes_in_hpa))

# Early exit guard (optional)
if (length(target_genes_in_hpa) < 5) {
  warning("Fewer than 5 targets found in HPA; enrichment may be underpowered.")
}

## ---- Optional: expression heatmap of targets across cell types ----
# Subset to targets and reshape wide
sc_subset <- sc_data %>%
  mutate(GeneName = toupper(GeneName)) %>%
  filter(GeneName %in% target_genes_in_hpa)

# If you want a quick heatmap (log2(nTPM+1)); limit to top-variance genes for readability
if (nrow(sc_subset) > 0) {
  sc_wide <- sc_subset %>%
    select(Gene.name, Cell.type, nTPM) %>%
    pivot_wider(names_from = Cell.type, values_from = nTPM)
  
  sc_mat <- as.data.frame(sc_wide)
  rownames(sc_mat) <- sc_mat$Gene.name
  sc_mat$Gene.name <- NULL
  sc_mat <- as.matrix(sc_mat)
  sc_mat_log <- log2(sc_mat + 1)
  
  # Optionally keep the top 50 most variable genes for plotting
  if (nrow(sc_mat_log) > 50) {
    vars <- apply(sc_mat_log, 1, stats::var, na.rm = TRUE)
    keep <- names(sort(vars, decreasing = TRUE))[1:50]
    sc_mat_log_plot <- sc_mat_log[keep, , drop = FALSE]
  } else {
    sc_mat_log_plot <- sc_mat_log
  }
  
  pheatmap(sc_mat_log_plot, scale = "row", show_rownames = TRUE,
           main = "BCL6B promoter targets (HPA single-cell types)\nlog2(nTPM+1)")
}

## ---- Enrichment via Fisher's Exact Test ----
# Choose an expression threshold; nTPM >= 10 is a common, conservative cut
expr_thresh <- 10

# Background = all genes present in HPA table (recommended)
background_genes <- setdiff(all_hpa_genes, target_genes_in_hpa)  # exclude targets from background to avoid double-counting in label

# Build long table with labels and expression flag
enrich_df <- sc_data %>%
  mutate(Gene.name = toupper(Gene.name)) %>%
  filter(Gene.name %in% c(target_genes_in_hpa, background_genes)) %>%
  mutate(
    GeneSet   = if_else(Gene.name %in% target_genes_in_hpa, "Targets", "Background"),
    Expressed = if_else(nTPM >= expr_thresh, 1L, 0L)
  )

# For each Cell.type, compute 2x2 and run Fisher
fisher_results <- enrich_df %>%
  group_by(Cell.type) %>%
  summarise(
    A = sum(GeneSet == "Targets"    & Expressed == 1L, na.rm = TRUE),  # Targets expressed
    B = sum(GeneSet == "Background" & Expressed == 1L, na.rm = TRUE),  # Background expressed
    C = sum(GeneSet == "Targets"    & Expressed == 0L, na.rm = TRUE),  # Targets not expressed
    D = sum(GeneSet == "Background" & Expressed == 0L, na.rm = TRUE),  # Background not expressed
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    fisher = list(fisher.test(matrix(c(A, B, C, D), nrow = 2))),
    p_value   = fisher$p.value,
    odds_ratio = as.numeric(fisher$estimate),
    conf_low   = as.numeric(fisher$conf.int[1]),
    conf_high  = as.numeric(fisher$conf.int[2])
  ) %>%
  ungroup() %>%
  select(-fisher) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# View the most enriched cell types
print(head(fisher_results, 15))

# Quick visualization of effect sizes (log2 OR) with FDR
fisher_results %>%
  mutate(log2OR = log2(odds_ratio)) %>%
  ggplot(aes(x = reorder(Cell.type, log2OR), y = log2OR, fill = p_adj < 0.05)) +
  geom_col() +
  coord_flip() +
  labs(x = "Cell type", y = "log2(odds ratio)",
       title = paste0("Enrichment of BCL6B promoter targets (nTPM ≥ ", expr_thresh, ")")) +
  scale_fill_manual(values = c("grey70", "black"), name = "FDR < 0.05") +
  theme_minimal(base_size = 12)
