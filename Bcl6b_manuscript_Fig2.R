#Code for Bcl6b Manuscript Figure 2.
#Analysis of the GTEx data set.
#Identify and characterize the Bcl6b gene correlation network.
#Across human tissues Bcl6b is correlated to a distinct angiogenic transcriptional program.

# Core tidyverse and plotting
library(tidyverse)
library(dplyr)
library(ggplot2)
library(vroom)
library(pheatmap)
library(RColorBrewer)

# Bioinformatics
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DOSE)
library(biomaRt)
library(enrichplot)

# Utilities for performance
library(pbapply)
library(matrixStats)
library(data.table)

#--- Set working directory
setwd("/Users/mdbruss/Documents/RStudioProjects_2/BCL6B")

#--- Load packages
library(vroom)
#if (!requireNamespace("pbapply", quietly = TRUE)) install.packages("pbapply")
library(pbapply)

#==============================================#
#/////Create Bcl6b Correlation Network/////////#
#==============================================#
#-----Block #2---------------#

#--- Load GTEx TPM data (skip 2-line header)
gtex_tpm <- vroom("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip = 2)

#--- Rename first 2 columns
colnames(gtex_tpm)[1:2] <- c("EnsemblID", "GeneSymbol")

#--- Filter out unannotated or duplicated gene symbols
gtex_tpm <- gtex_tpm[!duplicated(gtex_tpm$GeneSymbol) & gtex_tpm$GeneSymbol != "", ]

#------Block #3-------------------------#

#--- Save BCL6B vector before filtering rest
bcl6b_expr <- as.numeric(gtex_tpm[gtex_tpm$GeneSymbol == "BCL6B", -(1:2)])

#--- Pull expression matrix and attach gene names
expr_matrix <- as.matrix(gtex_tpm[, -(1:2)])
rownames(expr_matrix) <- gtex_tpm$GeneSymbol
gene_names <- gtex_tpm$GeneSymbol

#----Block #4------------------------#

#--- Filter: TPM > 1 in at least 10% of samples, chunked
chunk_size <- 1000
n_genes <- nrow(expr_matrix)
n_chunks <- ceiling(n_genes / chunk_size)
expressed_filter <- logical(n_genes)

for (i in seq_len(n_chunks)) {
  cat("Filtering chunk", i, "of", n_chunks, "\n")
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, n_genes)
  chunk <- expr_matrix[start:end, , drop = FALSE]
  expressed_filter[start:end] <- pbapply(chunk, 1, function(x) sum(x > 1) >= 0.10 * length(x))
}

#--- Apply filter
expr_matrix <- expr_matrix[expressed_filter, ]
gene_names <- gene_names[expressed_filter]

# Check result
length(gene_names)
#---------------------------------------#

#----Creating correlation matrix--------#
#--- Make sure the dimensions are aligned
stopifnot(length(bcl6b_expr) == ncol(expr_matrix))  # BCL6B vector must match sample count

#--- Compute Pearson correlation for each gene vs BCL6B
cor_vec <- pbapply::pbapply(expr_matrix, 1, function(x) cor(x, bcl6b_expr, method = "pearson"))

#--- Create final coexpression table
bcl6b_net <- tibble(
  GeneSymbol = gene_names,
  BCL6B_cor = cor_vec
) %>%
  arrange(desc(BCL6B_cor))  # most co-expressed at top

#--- Preview top 20
head(bcl6b_net, 20)

bcl6b_net_rev <- tibble(
  GeneSymbol = gene_names,
  BCL6B_cor = cor_vec
) %>%
  arrange(BCL6B_cor)

bcl6b_net_rev

tip_genes <- c("DLL4","KDR","NRP1","FLT4","ESM1","APLN","ANGPT2","CXCR4","PDGFB",
               "MMP14","FSCN1","PFKFB3","ROBO4","PLXND1","UNC5B")
notch_targets <- c("HES1","HEY1","HEY2","NRARP","JAG1")

subset_corr <- bcl6b_net |>
  dplyr::filter(GeneSymbol %in% c(tip_genes, notch_targets)) |>
  dplyr::arrange(dplyr::desc(BCL6B_cor))
subset_corr

dim(bcl6b_net)
head(bcl6b_net)
getwd()
write.csv(bcl6b_net, "GTEx_BCL6B_correlation_network.csv",row.names = FALSE)
bcl6b_net %>%
  filter(BCL6B_cor >0.6 & BCL6B_cor <0.7)
#===================================================#

#==============================================#
#/////Bcl6b Correlation Network histogram///////#
#==============================================#

bcl6b_net %>%
  arrange(BCL6B_cor)

ggplot(
  bcl6b_net %>% filter(BCL6B_cor >= -0.3, BCL6B_cor <= 1.0),
  aes(x = BCL6B_cor)
) +
  geom_histogram(
    binwidth = 0.1,
    boundary = -0.3,
    fill = "steelblue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(-0.3, 1.0, by = 0.1),
    limits = c(-0.3, 1.0),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    x = "BCL6B correlation",
    y = "Number of genes",
    title = "Distribution of gene–BCL6B correlations"
  ) +
  theme_bw(base_size = 13)
#======================================================#

#Zoomed in histogram

library(scales)

ggplot(
  bcl6b_net %>% filter(BCL6B_cor >= 0.5, BCL6B_cor <= 1.0),
  aes(x = BCL6B_cor)
) +
  geom_histogram(
    binwidth = 0.1,
    boundary = 0.5,          # align bins at 0.5, 0.6, ...
    fill = "steelblue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(0.5, 1.0, by = 0.1),
    limits = c(0.5, 1.0),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    x = "BCL6B correlation",
    y = "Number of genes",
    title = "High-correlation tail of gene–BCL6B associations",
    subtitle = "Zoomed view: correlations ≥ 0.5"
  ) +
  theme_bw(base_size = 13)
#======================================================#

head(bcl6b_net)

library(dplyr)
library(ggplot2)
library(ggrepel)

highlight_genes <- c(
  "CDH5", "PECAM1", "ROBO4", "FLT4", "KDR",
  "DLL4", "TEK", "RASIP1", "ESAM", "EGFL7", "VWF"
)

barcode_df <- bcl6b_net %>%
  dplyr::mutate(GeneSymbol = toupper(GeneSymbol)) %>%
  dplyr::arrange(desc(BCL6B_cor)) %>%
  dplyr::mutate(rank = dplyr::row_number())

highlight_df <- barcode_df %>%
  dplyr::filter(GeneSymbol %in% toupper(highlight_genes))

ggplot(barcode_df, aes(x = rank, y = BCL6B_cor)) +
  
  # Background points
  geom_point(size = 0.6, alpha = 0.4, color = "grey40") +
  
  # Highlighted points
  geom_point(
    data = highlight_df,
    size = 2,
    color = "#B2182B"
  ) +
  
  # Labels with leader lines
  ggrepel::geom_text_repel(
    data = highlight_df,
    aes(label = GeneSymbol),
    size = 3.5,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  
  labs(
    x = "Gene rank (by correlation with BCL6B)",
    y = "Pearson correlation (r)",
    title = "Genome-wide BCL6B coexpression landscape (GTEx)",
    subtitle = paste0(
      "All expressed genes ranked by correlation (n = ", nrow(barcode_df), "). ",
      "Highlighted: canonical endothelial / angiogenic genes"
    )
  ) +
  
  theme_classic(base_size = 13)

#==================================================#
#/////GSEA on BCL6B correlation network////////////#

#============================================================
# GSEA (GO Biological Process) on BCL6B coexpression ranking
#   Input: bcl6b_net (columns: GeneSymbol, BCL6B_cor)
#============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
})

#-----------------------------
# 0) Build ranked gene list
#-----------------------------
stopifnot(all(c("GeneSymbol", "BCL6B_cor") %in% colnames(bcl6b_net)))

rank_df <- bcl6b_net %>%
  dplyr::filter(!is.na(GeneSymbol), GeneSymbol != "", !is.na(BCL6B_cor)) %>%
  dplyr::mutate(GeneSymbol = toupper(GeneSymbol)) %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise(BCL6B_cor = max(BCL6B_cor), .groups = "drop") %>%  # just in case
  dplyr::arrange(desc(BCL6B_cor))

rank_vec <- rank_df$BCL6B_cor
names(rank_vec) <- rank_df$GeneSymbol

# (Optional) small jitter to break exact ties (helps fgsea a bit)
# rank_vec <- rank_vec + rnorm(length(rank_vec), mean = 0, sd = 1e-8)

#-----------------------------
# 1) Map SYMBOL -> ENTREZID
#-----------------------------
sym2ent <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = names(rank_vec),
  keytype  = "SYMBOL",
  columns  = c("ENTREZID")
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# Keep only mapped genes; switch names to ENTREZ IDs
rank_vec_entrez <- rank_vec[sym2ent$SYMBOL]
names(rank_vec_entrez) <- sym2ent$ENTREZID

# clusterProfiler expects decreasing order
rank_vec_entrez <- sort(rank_vec_entrez, decreasing = TRUE)

cat("Ranked genes (symbols):", length(rank_vec), "\n")
cat("Mapped to ENTREZ:", length(rank_vec_entrez), "\n")

#-----------------------------
# 2) Run GSEA: GO Biological Process
#-----------------------------
set.seed(1)

gsea_go_bp <- clusterProfiler::gseGO(
  geneList      = rank_vec_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,          # filter on nominal p; use qvalues for interpretation
  pAdjustMethod = "BH",
  eps           = 0,             # keeps very small p-values stable
  verbose       = FALSE
)

# If nothing returns, relax cutoffs:
# gsea_go_bp <- clusterProfiler::gseGO(
#   geneList = rank_vec_entrez, OrgDb = org.Hs.eg.db, keyType="ENTREZID", ont="BP",
#   minGSSize=10, maxGSSize=800, pvalueCutoff=1, pAdjustMethod="BH", eps=0, verbose=FALSE
# )

#-----------------------------
# 3) Inspect results table
#-----------------------------
gsea_res <- as.data.frame(gsea_go_bp) %>%
  dplyr::arrange(p.adjust)

print(head(gsea_res, 15))

#-----------------------------
# 4) Plots (dotplot + top enrichment curve)
#-----------------------------
if (nrow(gsea_res) > 0) {
  
  # Dotplot of top terms (by adjusted p)
  p1 <- enrichplot::dotplot(gsea_go_bp, showCategory = 20, split = ".sign") +
    ggplot2::facet_grid(. ~ .sign) +
    ggplot2::labs(
      title = "GSEA: GO Biological Process (BCL6B coexpression ranking)",
      subtitle = "Positive NES = enriched among genes positively correlated with BCL6B"
    ) +
    ggplot2::theme_classic(base_size = 13)
  
  print(p1)
  
  # Enrichment curve for top term
  top_id <- gsea_res$ID[1]
  p2 <- enrichplot::gseaplot2(
    gsea_go_bp,
    geneSetID = top_id,
    title = gsea_res$Description[1]
  )
  print(p2)
  
} else {
  message("No GO BP terms passed the cutoff. Try relaxing pvalueCutoff or minGSSize.")
}

#-----------------------------
# 5) Save results (optional)
#-----------------------------
# write.csv(gsea_res, "BCL6B_GTEx_GSEA_GOBP.csv", row.names = FALSE)

#============================================================
# Custom GO BP dotplot from GSEA results
#============================================================

library(dplyr)
library(ggplot2)

# Convert gsea result object to dataframe if needed
gsea_df <- as.data.frame(gsea_go_bp)

#============================================================
# Tunable parameters
#============================================================
MIN_SET_SIZE <- 100
MAX_SET_SIZE <- 300
MIN_ABS_NES  <- 1.5
N_TERMS      <- 20

#============================================================
# Dotplot with NES filtering
#============================================================
Fig_BCL6B_GO_dotplot <- gsea_df %>%
  dplyr::filter(
    p.adjust < 0.05,
    setSize >= MIN_SET_SIZE,
    setSize <= MAX_SET_SIZE,
    abs(NES) >= MIN_ABS_NES
  ) %>%
  dplyr::slice_head(n = N_TERMS) %>%
  ggplot(
    aes(
      x     = NES,
      y     = reorder(Description, NES),
      size  = setSize,
      color = p.adjust
    )
  ) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(
    low  = "#2166AC",
    high = "#B2182B",
    name = "Adj. p-value",
    trans = "reverse"
  ) +
  scale_size_continuous(
    name = "Gene set size"
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "GO Biological Process enrichment",
    subtitle = paste0(
      "GSEA on BCL6B coexpression ranking (GTEx)\n",
      "Gene set size: ", MIN_SET_SIZE, "–", MAX_SET_SIZE,
      "; |NES| ≥ ", MIN_ABS_NES
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

Fig_BCL6B_GO_dotplot
#==============================================#
#==============================================#

#/////Create Partial Correlation Network///////#

#============================================================
# Partial correlation network: gene ~ BCL6B | PECAM1
# (Partial Spearman = Pearson on ranks)
#============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(matrixStats)
})

#-----------------------------
# SETTINGS YOU CAN TUNE
#-----------------------------
MARKER_GENE <- "PECAM1"     # try "CDH5" as sensitivity check
CHUNK_SIZE  <- 2000         # increase if you have RAM

#-----------------------------
# Sanity checks
#-----------------------------
stopifnot(exists("expr_matrix"))
stopifnot(exists("bcl6b_expr"))
stopifnot(length(bcl6b_expr) == ncol(expr_matrix))

if (!MARKER_GENE %in% rownames(expr_matrix)) {
  stop(paste0("Marker gene ", MARKER_GENE, " not found in expr_matrix rownames."))
}
if (!"BCL6B" %in% rownames(expr_matrix)) {
  message("Note: BCL6B row is not in expr_matrix rownames (that's OK if bcl6b_expr came from before filtering).")
}

marker <- as.numeric(expr_matrix[MARKER_GENE, ])
y      <- as.numeric(bcl6b_expr)
n      <- ncol(expr_matrix)

# Rank + Z-score covariates (Spearman = Pearson on ranks)
rz <- scale(rank(marker, ties.method = "average"))[, 1]
ry <- scale(rank(y,      ties.method = "average"))[, 1]

# Spearman correlation between BCL6B and PECAM1
r_yz <- as.numeric(cor(y, marker, method = "spearman"))

#-----------------------------
# Chunked Spearman correlations for each gene:
#   r_gm = cor(gene, marker)   (Spearman)
#   r_gb = cor(gene, BCL6B)    (Spearman)
#-----------------------------
ng <- nrow(expr_matrix)
r_gm <- numeric(ng)
r_gb <- numeric(ng)

for (i in seq(1, ng, by = CHUNK_SIZE)) {
  j <- min(i + CHUNK_SIZE - 1, ng)
  X <- expr_matrix[i:j, , drop = FALSE]
  
  # Row-wise ranks for Spearman
  R <- matrixStats::rowRanks(X, ties.method = "average", preserveShape = TRUE)
  
  # Z-score each row of ranks
  mu <- rowMeans(R)
  sd <- matrixStats::rowSds(R)
  sd[sd == 0] <- NA_real_  # constant rows safeguard
  Z <- sweep(R, 1, mu, "-")
  Z <- sweep(Z, 1, sd, "/")
  Z[is.na(Z)] <- 0
  
  # Spearman r = Pearson on ranks = normalized dot product
  r_gm[i:j] <- drop(Z %*% rz) / (n - 1)
  r_gb[i:j] <- drop(Z %*% ry) / (n - 1)
}

#-----------------------------
# Partial Spearman for each gene: r_{g,y·z}
#-----------------------------
denom <- sqrt((1 - r_gm^2) * (1 - r_yz^2))
denom[denom == 0] <- NA_real_

pcor <- (r_gb - r_gm * r_yz) / denom

bcl6b_partial <- tibble(
  GeneSymbol = rownames(expr_matrix),
  r_gene_marker   = r_gm,          # Spearman(gene, PECAM1)
  r_gene_BCL6B    = r_gb,          # Spearman(gene, BCL6B)
  r_partial_BCL6B_given_marker = pcor
) %>%
  arrange(desc(r_partial_BCL6B_given_marker))

#-----------------------------
# Quick checks / examples
#-----------------------------
print(head(bcl6b_partial, 20))

bcl6b_partial %>%
  filter(r_partial_BCL6B_given_marker >0.5)

tip_genes <- c("DLL4","KDR","NRP1","FLT4","ESM1","APLN","ANGPT2","CXCR4","PDGFB",
               "MMP14","FSCN1","PFKFB3","ROBO4","PLXND1","UNC5B")
notch_targets <- c("HES1","HEY1","HEY2","NRARP","JAG1")

bcl6b_partial %>%
  filter(GeneSymbol %in% c(tip_genes, notch_targets, "PECAM1", "CDH5", "VWF", "KDR", "TEK", "ROBO4")) %>%
  arrange(desc(r_partial_BCL6B_given_marker)) %>%
  print(n = Inf)

#-----------------------------
# Save
#-----------------------------
write.csv(bcl6b_partial,
          file = paste0("GTEx_BCL6B_partialSpearman_given_", MARKER_GENE, ".csv"),
          row.names = FALSE)
#====================================================================#

#/////Partial Correlation GSEA Dotplot///////////////#

#============================================================
# GSEA (GO BP) on partial correlation:
# ranking = r_partial_BCL6B_given_marker
#============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

#-----------------------------
# SETTINGS YOU CAN TUNE
#-----------------------------
MIN_SET_SIZE <- 20
MAX_SET_SIZE <- 300
MIN_ABS_NES  <- 1.5
N_TERMS      <- 30

#-----------------------------
# 1) Build ranked vector (partial correlation) - FIXED
#-----------------------------
rank_tbl_partial <- bcl6b_partial %>%
  dplyr::transmute(
    GeneSymbol = toupper(GeneSymbol),
    stat = as.numeric(r_partial_BCL6B_given_marker)
  ) %>%
  dplyr::filter(is.finite(stat)) %>%                 # drops NaN/Inf (e.g. PECAM1)
  dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%  # safety if any duplicates
  dplyr::arrange(desc(stat))

rank_vec_partial <- rank_tbl_partial$stat
names(rank_vec_partial) <- rank_tbl_partial$GeneSymbol

# IMPORTANT: gseGO expects decreasing sort
rank_vec_partial <- sort(rank_vec_partial, decreasing = TRUE)

# Quick sanity check
stopifnot(length(rank_vec_partial) == length(names(rank_vec_partial)))
head(rank_vec_partial)


#-----------------------------
# 2) Run GO BP GSEA
#-----------------------------
gsea_go_bp_partial <- clusterProfiler::gseGO(
  geneList      = rank_vec_partial,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = MIN_SET_SIZE,
  maxGSSize     = MAX_SET_SIZE,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,          # filter later
  verbose       = FALSE
)

# Convert to dataframe for plotting
gsea_df_partial <- as.data.frame(gsea_go_bp_partial)

#-----------------------------
# 3) Dotplot (same style as before)
#-----------------------------
Fig_BCL6B_partial_GO_dotplot <- gsea_df_partial %>%
  dplyr::filter(
    p.adjust < 0.05,
    setSize >= MIN_SET_SIZE,
    setSize <= MAX_SET_SIZE,
    abs(NES) >= MIN_ABS_NES
  ) %>%
  dplyr::slice_head(n = N_TERMS) %>%
  ggplot(
    aes(
      x     = NES,
      y     = reorder(Description, NES),
      size  = setSize,
      color = p.adjust
    )
  ) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(
    low  = "#2166AC",
    high = "#B2182B",
    name = "Adj. p-value",
    trans = "reverse"
  ) +
  scale_size_continuous(name = "Gene set size") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "GO BP enrichment (partial correlation)",
    subtitle = paste0(
      "Ranking = partial Spearman: gene ~ BCL6B | PECAM1\n",
      "Gene set size: ", MIN_SET_SIZE, "–", MAX_SET_SIZE,
      "; |NES| ≥ ", MIN_ABS_NES
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

Fig_BCL6B_partial_GO_dotplot

#====================================================#
#////////////////////////////////////////////////////#
#====================================================#

#/////Partial Correlation Overrepresentation Dotplot/////////#

bcl6b_partial_hicor_gene_vec <- bcl6b_partial %>%
  filter(r_partial_BCL6B_given_marker >0.5) %>%
  pull(GeneSymbol)

bcl6b_partial_hicor_gene_vec

############################################################
# GO BP Over-Representation Analysis (NO term filtering)
# Partial BCL6B correlation gene set (GTEx)
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
})

#-----------------------------
# 1) Inputs
#-----------------------------
genes_of_interest <- toupper(bcl6b_partial_hicor_gene_vec)
background_genes  <- toupper(rownames(expr_matrix))

stopifnot(length(genes_of_interest) >= 5)
stopifnot(length(background_genes) >= 1000)

# Keep only genes that are actually in the background (safety)
genes_of_interest <- intersect(genes_of_interest, background_genes)

cat("Genes of interest (in background):", length(genes_of_interest), "\n")
cat("Background genes:", length(background_genes), "\n")

#-----------------------------
# 2) enrichGO (GO Biological Process)
#-----------------------------
ego_partial <- enrichGO(
  gene          = genes_of_interest,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,     # don't filter here; filter downstream if desired
  qvalueCutoff  = 1,
  readable      = TRUE
)

ego_df <- as.data.frame(ego_partial) %>%
  arrange(p.adjust)

ego_df
#-----------------------------
# 3) View results
#-----------------------------
print(ego_df, n = 50)

#-----------------------------
# USER-TUNABLE FILTERS
#-----------------------------
MIN_COUNT <- 9      # minimum number of genes hitting the GO term
TOP_N     <- 25     # number of terms to show

#-----------------------------
# 4) Bar chart with FoldEnrichment + Count filter
#-----------------------------
if (nrow(ego_df) > 0) {
  
  plot_df <- ego_df %>%
    mutate(
      # Parse ratios like "7/82"
      GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
        as.numeric(sub(".*/", "", GeneRatio)),
      BgRatio_num   = as.numeric(sub("/.*", "", BgRatio)) /
        as.numeric(sub(".*/", "", BgRatio)),
      
      # Continuous enrichment metric
      FoldEnrichment = GeneRatio_num / BgRatio_num,
      
      # For coloring
      neglog_padj = -log10(p.adjust + 1e-300)
    ) %>%
    filter(Count >= MIN_COUNT) %>%
    arrange(desc(FoldEnrichment)) %>%
    slice_head(n = TOP_N) %>%
    mutate(
      Description = factor(Description, levels = rev(Description))  # keep plotted order top->bottom
    )
  
  ggplot(plot_df, aes(x = FoldEnrichment, y = Description, fill = neglog_padj)) +
    geom_col(width = 0.75) +
    scale_fill_gradient(
      low = "#2166AC",
      high = "#B2182B",
      name = expression(-log[10](p[adj]))
    ) +
    labs(
      x = "Fold enrichment (GeneRatio / BgRatio)",
      y = NULL,
      title = "GO BP over-representation analysis",
      subtitle = paste0(
        "Partial BCL6B correlation network (r > 0.5)\n",
        "Min genes per term = ", MIN_COUNT,
        "; showing top ", TOP_N, " terms by fold enrichment"
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
} else {
  message("No GO BP terms returned.")
}


############################################################
# END
############################################################


#====================================================#
#///////Test Code///////////////////////#
#====================================================#

