#Bcl6b Manuscript Figure 4.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

# -------------------------
# 0) File paths (from your snippet)
# -------------------------
raw_dir <- "/Users/mdbruss/Documents/RStudioProjects_2/BCL6B"

sed_prefix <- file.path(raw_dir, "GSM4816919_sample3_sedentary")
run_prefix <- file.path(raw_dir, "GSM4816920_sample1_2weeks")

cache_dir <- file.path(raw_dir, "processed_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cache_rds <- file.path(cache_dir, "PMID34358431_EC_SedRun_processed.rds")

# -------------------------
# 1) Read helper
# -------------------------
read_mtx_triplet <- function(prefix) {
  mtx  <- readMM(paste0(prefix, "_matrix.mtx.gz"))
  gene <- read.delim(paste0(prefix, "_genes.tsv.gz"), header = FALSE)
  bar  <- read.delim(paste0(prefix, "_barcodes.tsv.gz"), header = FALSE)
  
  keep <- !duplicated(gene$V2)
  mtx  <- mtx[keep, , drop = FALSE]
  
  rownames(mtx) <- gene$V2[keep]
  colnames(mtx) <- bar$V1
  
  mtx
}

# -------------------------
# 2) Build or load processed object
# -------------------------
build_or_load <- function(force_rebuild = FALSE) {
  if (file.exists(cache_rds) && !force_rebuild) {
    message("Loading cached object:\n  ", cache_rds)
    return(readRDS(cache_rds))
  }
  
  message("Rebuilding from raw mtx/tsv...")
  
  mat_sed <- read_mtx_triplet(sed_prefix)
  mat_run <- read_mtx_triplet(run_prefix)
  
  sed <- CreateSeuratObject(mat_sed, project = "Sedentary")
  run <- CreateSeuratObject(mat_run, project = "Run")
  
  merged <- merge(sed, y = run, add.cell.ids = c("Sed", "Run"), project = "BCL6B_Exercise")
  merged$Condition <- ifelse(grepl("^Sed_", colnames(merged)), "Sed", "Run")
  
  DefaultAssay(merged) <- "RNA"
  merged <- NormalizeData(merged)
  merged <- FindVariableFeatures(merged)
  merged <- ScaleData(merged)
  merged <- RunPCA(merged)
  merged <- RunUMAP(merged, dims = 1:20)
  merged <- FindNeighbors(merged, dims = 1:20)
  merged <- FindClusters(merged, resolution = 0.4)
  
  saveRDS(merged, cache_rds)
  message("Saved cached object:\n  ", cache_rds)
  
  merged
}

merged <- build_or_load(force_rebuild = FALSE)

# -------------------------
# 3) Make bcl6b_strict_DE with tunable "top_pct"
# -------------------------
make_bcl6b_high_DE <- function(
    merged,
    top_pct = 0.25,                 # 0.25=top 25% of NONZERO; 0.50=top 50%
    min.pct = 0.05,
    logfc.threshold = 0,
    test.use = "wilcox",
    out_csv = file.path(cache_dir, "Run_Bcl6b_high_vs_low_DE.csv")
) {
  stopifnot(top_pct > 0 && top_pct < 1)
  
  run_cells <- subset(merged, subset = Condition == "Run")
  
  bcl <- FetchData(run_cells, "Bcl6b")$Bcl6b
  bcl <- pmax(bcl, 0)
  nonzero <- bcl[bcl > 0]
  if (length(nonzero) < 10) stop("Too few nonzero Bcl6b cells to threshold reliably.")
  
  thr <- as.numeric(quantile(nonzero, probs = 1 - top_pct, names = FALSE))
  run_cells$Bcl6b_high <- bcl > thr
  
  message("Bcl6b-high threshold = ", signif(thr, 4),
          " (top ", top_pct*100, "% of nonzero)")
  
  Idents(run_cells) <- run_cells$Bcl6b_high
  
  de <- FindMarkers(
    run_cells,
    ident.1 = TRUE,
    ident.2 = FALSE,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = test.use
  ) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(avg_log2FC))
  
  write.csv(de, out_csv, row.names = FALSE)
  de
}

# Strict (top 25% of nonzero)
bcl6b_strict_DE <- make_bcl6b_high_DE(merged, top_pct = 0.25,
                                      out_csv = file.path(cache_dir, "Run_Bcl6b_top25pct_DE.csv"))

# Strict_mod (top 33% of nonzero)
bcl6b_midstrict_DE <- make_bcl6b_high_DE(merged, top_pct = 0.33,
                                      out_csv = file.path(cache_dir, "Run_Bcl6b_top33pct_DE.csv"))

# Moderate (top 50% of nonzero)
bcl6b_mod_DE <- make_bcl6b_high_DE(merged, top_pct = 0.50,
                                   out_csv = file.path(cache_dir, "Run_Bcl6b_top50pct_DE.csv"))

head(bcl6b_strict_DE)

bcl6b_midstrict_DE %>%
  filter(p_val_adj <0.05)

bcl6b_mod_DE %>%
  filter(p_val_adj <0.05)

bcl6b_mod_DE %>%
  filter (gene == "Cdh5")

angiogenesis_sets <- list(
  Tip_like = c("Dll4","Kdr","Flt1","Apln","Esm1","Sox17","Sox18","Mcam","Rgcc"),
  Notch_output = c("Hey1","Hey2","Hes1","Nrarp"),
  ECM_remodel = c("Col4a1","Col4a2","Hspg2","Bgn","Eln","Fbln5","Tnc","Mmp2"),
  Permeability = c("Plvap","Kdr","Cdh5","Pecam1","Mfsd2a","Tie1","Tek"),
  Lymphangiogenic = c("Vegfc","Flt4","Lyve1","Prox1")
)

run_fgsea_custom <- function(de_tbl, pathways) {
  ranks <- de_tbl %>%
    mutate(score = avg_log2FC * -log10(p_val_adj + 1e-300)) %>%
    dplyr::select(gene, score) %>%
    distinct() %>%
    tibble::deframe()
  
  # keep only genes present in ranks
  pathways_f <- lapply(pathways, function(g) intersect(g, names(ranks)))
  
  fgsea(pathways = pathways_f, stats = ranks, nperm = 10000) %>%
    arrange(padj)
}

fg_strict_custom <- run_fgsea_custom(bcl6b_strict_DE, angiogenesis_sets)
fg_mod_custom    <- run_fgsea_custom(bcl6b_mod_DE, angiogenesis_sets)
fg_midstrict_custom    <- run_fgsea_custom(bcl6b_midstrict_DE, angiogenesis_sets)

fg_strict_custom
fg_midstrict_custom
fg_mod_custom
#==========================================================================#
############################################################################
#===========================================================================

#===============================================================#
#/////Generate UMAP "cluster maps" for Sed and Run/////////////#
#===============================================================#

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# 1) Define Run-based moderate cutoff (top 50% of nonzero Bcl6b in RUN)
# -----------------------------
run_cells <- subset(merged, subset = Condition == "Run")
bcl_run <- pmax(FetchData(run_cells, "Bcl6b")$Bcl6b, 0)

thr_mod_run <- as.numeric(quantile(bcl_run[bcl_run > 0], probs = 0.50, names = FALSE))
message("Run-based moderate Bcl6b cutoff = ", signif(thr_mod_run, 4))

# Apply cutoff to all cells (for consistent highlighting across panels)
bcl_all <- pmax(FetchData(merged, "Bcl6b")$Bcl6b, 0)
merged$Bcl6b_mod50_runCut <- bcl_all > thr_mod_run

# -----------------------------
# 2) Build plotting dataframe
# -----------------------------
umap_mat <- Embeddings(merged, "umap")
umap_x <- colnames(umap_mat)[1]
umap_y <- colnames(umap_mat)[2]

plot_df <- as.data.frame(umap_mat) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    Condition = merged$Condition[cell],
    cluster   = as.factor(merged$seurat_clusters[cell]),
    bcl_mod   = merged$Bcl6b_mod50_runCut[cell]
  )

# Keep cluster colors consistent across both panels
cluster_levels <- sort(unique(plot_df$cluster))
plot_df$cluster <- factor(plot_df$cluster, levels = cluster_levels)

# -----------------------------
# 3) A small theme helper (optional)
# -----------------------------
base_umap_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

# -----------------------------
# 4) FIG 4a: Sedentary only
# -----------------------------
fig4a_sed <- plot_df %>%
  filter(Condition == "Sed") %>%
  ggplot(aes_string(x = umap_x, y = umap_y)) +
  geom_point(aes(color = cluster), size = 0.40, alpha = 0.35) +
  geom_point(
    data = plot_df %>% filter(Condition == "Sed", bcl_mod),
    color = "red3", size = 1.05, alpha = 0.90
  ) +
  labs(
    title = "Fig. 4a | Sedentary EC clusters",
    subtitle = "Red = Bcl6b-high (Run-defined cutoff; top 50% of nonzero in Run)",
    x = "UMAP 1", y = "UMAP 2", color = "Cluster"
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  base_umap_theme

# -----------------------------
# 5) FIG 4b: Run only
# -----------------------------
fig4b_run <- plot_df %>%
  filter(Condition == "Run") %>%
  ggplot(aes_string(x = umap_x, y = umap_y)) +
  geom_point(aes(color = cluster), size = 0.40, alpha = 0.35) +
  geom_point(
    data = plot_df %>% filter(Condition == "Run", bcl_mod),
    color = "red3", size = 1.05, alpha = 0.90
  ) +
  labs(
    title = "Fig. 4b | Run EC clusters",
    subtitle = "Red = Bcl6b-high (top 50% of nonzero Bcl6b in Run)",
    x = "UMAP 1", y = "UMAP 2", color = "Cluster"
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  base_umap_theme

# Print to screen
fig4a_sed
fig4b_run

ggsave("Fig4a_Sed_UMAP_clusters_Bcl6bHigh.pdf", fig4a_sed, width = 4.5, height = 4.2, units = "in")
ggsave("Fig4b_Run_UMAP_clusters_Bcl6bHigh.pdf", fig4b_run, width = 4.5, height = 4.2, units = "in")

# or PNG if needed
ggsave("Fig4a_Sed_UMAP_clusters_Bcl6bHigh.png", fig4a_sed, width = 4.5, height = 4.2, units = "in", dpi = 400)
ggsave("Fig4b_Run_UMAP_clusters_Bcl6bHigh.png", fig4b_run, width = 4.5, height = 4.2, units = "in", dpi = 400)
#====================================================================================================#
#==================================================================================================#


##########################################################
#//////Fig. 4e Bcl6b+ DEG cluster dotplot////////////////#
#########################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

de_tbl <- bcl6b_mod_DE
padj_cut <- 0.05

merged$seurat_clusters <- factor(merged$seurat_clusters)
cluster_levels <- levels(merged$seurat_clusters)

run_only <- subset(merged, subset = Condition == "Run")
run_only$seurat_clusters <- factor(run_only$seurat_clusters, levels = cluster_levels)

# all genes tested
genes_use <- de_tbl %>%
  filter(!is.na(p_val_adj)) %>%
  pull(gene) %>%
  unique()

# cluster colors from merged UMAP
p_umap_ref <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters")
gb <- ggplot_build(p_umap_ref)
col_df <- gb$data[[1]] %>% dplyr::select(group, colour) %>% distinct()
stopifnot(nrow(col_df) == length(cluster_levels))
cluster_colors <- col_df$colour
names(cluster_colors) <- cluster_levels

# ---- 1) pct expressing per cluster (Run only) ----
pct_mat <- AggregateExpression(
  run_only,
  assays   = "RNA",
  features = intersect(genes_use, rownames(run_only)),
  group.by = "seurat_clusters",
  slot     = "counts",
  return.seurat = FALSE
)$RNA

# pct_mat here is SUM counts (pseudobulk), so we need pct.exp separately:
# Seurat has built-in for DotPlot (it computes pct.exp + avg.exp).
dp <- DotPlot(
  run_only,
  features = intersect(genes_use, rownames(run_only)),
  group.by = "seurat_clusters"
)

# DotPlot returns a ggplot; pull the data it computed:
dp_dat <- dp$data %>%
  transmute(
    gene    = features.plot,
    cluster = as.character(id),
    pct_exp = pct.exp,      # 0-100
    avg_exp = avg.exp.scaled # scaled average used for dotplot; good “specificity-ish”
  ) %>%
  mutate(
    cluster = factor(cluster, levels = cluster_levels)
  ) %>%
  left_join(de_tbl %>% dplyr::select(gene, p_val_adj), by = "gene") %>%
  mutate(
    neglog10_p = -log10(pmax(p_val_adj, 1e-300)),
    sig = p_val_adj < padj_cut
  )

# ---- 2) Top-1 assignment by pct_exp (hard assignment) ----
gene_cluster_top1_pct <- dp_dat %>%
  group_by(gene) %>%
  slice_max(order_by = pct_exp, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---- 2) Top-1 assignment by specificity (DotPlot scaled avg exp) ----
gene_cluster_top1_spec <- dp_dat %>%
  group_by(gene) %>%
  slice_max(order_by = avg_exp, n = 1, with_ties = FALSE) %>%
  ungroup()


# ---- 3) Plot: grey all genes, colored significant genes ----
p_gene_cluster_top1_pct <- ggplot(gene_cluster_top1_pct, aes(x = cluster, y = neglog10_p)) +
  geom_hline(yintercept = -log10(padj_cut),
             linetype = "dotted", linewidth = 0.6, color = "grey40") +
  geom_jitter(
    data = gene_cluster_top1_pct %>% filter(!sig),
    aes(size = pct_exp),
    width = 0.22, height = 0, alpha = 0.25, color = "grey70"
  ) +
  geom_jitter(
    data = gene_cluster_top1_pct %>% filter(sig),
    aes(size = pct_exp, color = cluster),
    width = 0.22, height = 0, alpha = 0.90
  ) +
  scale_color_manual(values = cluster_colors, guide = "none") +
  scale_size_continuous(name = "% cells expressing\n(Run; top-1 cluster)") +
  labs(
    x = "Seurat cluster (top-1 by % cells expressing; Run only)",
    y = expression(-log[10]("adj. p-value (Bcl6b-high vs low, Run)")),
    title = "Top-1 cluster localization of Bcl6b-associated genes",
    subtitle = "Hard assignment by detection rate (% expressing) reduces bias from high-UMI clusters."
  ) +
  coord_cartesian(ylim = c(0, 15)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

p_gene_cluster_top1_pct

library(dplyr)

run_only@meta.data %>%
  dplyr::transmute(
    seurat_clusters = as.factor(seurat_clusters),
    Bcl6b_high = as.logical(Bcl6b_mod50_runCut)
  ) %>%
  dplyr::count(seurat_clusters, Bcl6b_high, name = "n") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Bcl6b_high) %>%
  dplyr::arrange(dplyr::desc(frac))



##############################################################
############################################################
#//////////Different more program associated analysis///////#
#============================================================
# Bcl6b Manuscript – Figure 4
# Cluster-level enrichment of Bcl6b-associated gene program
#============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

#------------------------------------------------------------
# STEP 0: Inputs
#------------------------------------------------------------
padj_cut <- 0.05

# Bcl6b-associated gene set (108 genes)
bcl6b_genes <- bcl6b_mod_DE %>%
  dplyr::filter(p_val_adj < padj_cut) %>%
  dplyr::pull(gene) %>%
  unique()

length(bcl6b_genes)  # sanity check (should be ~108)

#------------------------------------------------------------
# STEP 1: Work in RUN ECs only & define cluster identities
#------------------------------------------------------------
run_only <- subset(merged, subset = Condition == "Run")
Idents(run_only) <- run_only$seurat_clusters

background_genes <- rownames(run_only)

#------------------------------------------------------------
# STEP 2: Identify cluster-defining marker genes
# (These genes DEFINE the cluster identities)
#------------------------------------------------------------
cluster_markers <- FindAllMarkers(
  run_only,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#------------------------------------------------------------
# STEP 3: Enrichment of Bcl6b DEGs in cluster marker programs
#------------------------------------------------------------
cluster_enrichment <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(
    marker_genes = list(unique(gene)),
    n_markers = length(marker_genes[[1]]),
    overlap_genes = list(intersect(marker_genes[[1]], bcl6b_genes)),
    overlap = length(overlap_genes[[1]]),
    .groups = "drop"
  ) %>%
  mutate(
    # Fisher's exact test per cluster
    p_value = map_dbl(marker_genes, function(cl_genes) {
      mat <- matrix(
        c(
          length(intersect(cl_genes, bcl6b_genes)),
          length(setdiff(cl_genes, bcl6b_genes)),
          length(setdiff(bcl6b_genes, cl_genes)),
          length(setdiff(background_genes,
                         union(cl_genes, bcl6b_genes)))
        ),
        nrow = 2
      )
      fisher.test(mat)$p.value
    }),
    p_adj = p.adjust(p_value, method = "BH"),
    neglog10_FDR = -log10(p_adj),
    frac_cluster_markers = overlap / n_markers
  ) %>%
  arrange(desc(neglog10_FDR))

cluster_enrichment

#------------------------------------------------------------
# STEP 4: Plot – Which clusters are defined by Bcl6b DEGs?
#------------------------------------------------------------
ggplot(cluster_enrichment,
       aes(x = factor(cluster),
           y = neglog10_FDR,
           size = overlap,
           color = frac_cluster_markers)) +
  geom_point(alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey40") +
  scale_color_viridis_c(
    name = "Fraction of\ncluster markers"
  ) +
  scale_size_continuous(
    name = "# Bcl6b DEGs"
  ) +
  labs(
    x = "Endothelial cluster",
    y = expression(-log[10]("FDR")),
    title = "Clusters transcriptionally defined by the Bcl6b-associated gene program",
    subtitle = "Enrichment of Bcl6b⁺ DEGs within cluster-defining marker signatures (Run ECs)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
#================================================================#
#=================================================================#

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

#---------------------------------------
# Inputs
#---------------------------------------
top_pct <- 0.50     # match your mod definition; change to 0.25 or 0.33 if desired
min.pct <- 0.05
logfc.threshold <- 0

#---------------------------------------
# Build run_cells and define Bcl6b_high
#---------------------------------------
run_cells <- subset(merged, subset = Condition == "Run")

bcl <- FetchData(run_cells, "Bcl6b")[, 1]
bcl <- pmax(bcl, 0)
nonzero <- bcl[bcl > 0]
stopifnot(length(nonzero) >= 10)

thr <- as.numeric(quantile(nonzero, probs = 1 - top_pct, names = FALSE))
run_cells$Bcl6b_high <- bcl > thr

message("Bcl6b-high threshold = ", signif(thr, 4),
        " (top ", top_pct*100, "% of nonzero)")

# IMPORTANT: set identities so ident.1=TRUE works
Idents(run_cells) <- "Bcl6b_high"
table(Idents(run_cells))  # sanity check you have TRUE/FALSE

#---------------------------------------
# Cluster-adjusted DE (logistic regression)
#---------------------------------------
de_cluster_adj <- FindMarkers(
  run_cells,
  ident.1 = TRUE,
  ident.2 = FALSE,
  test.use = "LR",
  latent.vars = "seurat_clusters",
  min.pct = min.pct,
  logfc.threshold = logfc.threshold
) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj)

de_cluster_adj %>% dplyr::filter(p_val_adj < 0.05) %>% head(30)


# Cluster-adjusted DE (logistic regression)
de_cluster_adj <- FindMarkers(
  run_cells,
  ident.1 = TRUE,
  ident.2 = FALSE,
  test.use = "LR",
  latent.vars = "seurat_clusters",
  min.pct = 0.05,
  logfc.threshold = 0
) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj)

bcl6b_genes <- bcl6b_mod_DE %>% filter(p_val_adj < 0.05) %>% pull(gene) %>% unique()

run_only <- subset(merged, subset = Condition == "Run")

run_only <- AddModuleScore(run_only, features = list(bcl6b_genes), name = "Bcl6bDEG_Score")

FeaturePlot(run_only, features = "Bcl6bDEG_Score1", reduction = "umap")
VlnPlot(run_only, features = "Bcl6bDEG_Score1", group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(run_only, features = "Bcl6b", reduction = "umap")
VlnPlot(run_only, features = "Bcl6b", group.by = "seurat_clusters", pt.size = 0)
#================================================================================#
#================================================================================#

#===============================================================
# Fig 4c–d: Cluster composition changes (paired bars per cluster)
#   Fig 4c = total cells per cluster (Sed vs Run)
#   Fig 4d = Bcl6b+ cells per cluster (Run-defined 50% cutoff)
#   Stats: per-cluster Fisher’s exact (enrichment) + BH correction
#   Plot: paired bars per cluster; cluster colors match UMAP
#===============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(scales)
})

stopifnot("Condition" %in% colnames(merged@meta.data))
stopifnot("seurat_clusters" %in% colnames(merged@meta.data))
stopifnot("Bcl6b_mod50_runCut" %in% colnames(merged@meta.data))

merged$Condition <- factor(merged$Condition, levels = c("Sed","Run"))
merged$seurat_clusters <- factor(merged$seurat_clusters)
cluster_levels <- levels(merged$seurat_clusters)

# ---- Pull cluster colors EXACTLY as in your UMAP (if not already present) ----
if (!exists("cluster_colors")) {
  p_umap_ref <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters")
  gb <- ggplot_build(p_umap_ref)
  col_df <- gb$data[[1]] %>% dplyr::select(group, colour) %>% distinct()
  stopifnot(nrow(col_df) == length(cluster_levels))
  cluster_colors <- col_df$colour
  names(cluster_colors) <- cluster_levels
}

base_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

#===============================================================
# Helper: per-cluster enrichment Fisher test (Sed vs Run)
# For each cluster k:
# [ Run_in_k, Run_not_k
#   Sed_in_k, Sed_not_k ]
#===============================================================
cluster_fisher <- function(wide_df) {
  # wide_df must include: cluster, Sed, Run
  total_run <- sum(wide_df$Run)
  total_sed <- sum(wide_df$Sed)
  
  wide_df %>%
    rowwise() %>%
    mutate(
      fisher_p = {
        mat <- matrix(
          c(Run, total_run - Run,
            Sed, total_sed - Sed),
          nrow = 2, byrow = TRUE
        )
        fisher.test(mat)$p.value
      }
    ) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(fisher_p, method = "BH"),
      sig = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      )
    )
}

#===============================================================
# Make a "safe" metadata tibble with numeric cluster ordering
#===============================================================

# ---- Numeric cluster order (0,1,2,...,10) WITHOUT breaking cluster 0 ----
cluster_levels <- levels(merged$seurat_clusters)
cluster_levels <- as.character(sort(as.integer(cluster_levels)))  # key line

# Apply consistently to Seurat object
merged$seurat_clusters <- factor(
  merged$seurat_clusters,
  levels = cluster_levels
)

# Build safe metadata tibble (plain vectors, no Rle/list columns)
md <- merged@meta.data %>%
  as.data.frame() %>%
  tibble::as_tibble(rownames = "cell") %>%
  mutate(
    Condition = as.character(Condition),
    cluster   = as.character(seurat_clusters),   # keep character first
    bcl_pos   = as.logical(Bcl6b_mod50_runCut)
  ) %>%
  mutate(
    cluster = factor(cluster, levels = cluster_levels)
  )


#===============================================================
# Fig 4c (FIXED): Cluster composition (% of each condition in cluster k)
#   - y-axis: % of ALL cells in Sed (or Run) that fall in cluster k
#   - styling: Sed=blue, Run=red
#   - stats: per-cluster Fisher (in cluster vs not), BH-adjusted
#===============================================================

# ---- totals per condition (compute from md, NOT tot_wide) ----
total_sed <- sum(md$Condition == "Sed", na.rm = TRUE)
total_run <- sum(md$Condition == "Run", na.rm = TRUE)

stopifnot(total_sed > 0, total_run > 0)

# ---- counts per cluster x condition (wide) + Fisher enrichment ----
tot_wide <- md %>%
  dplyr::count(cluster, Condition, name = "n") %>%
  tidyr::complete(
    cluster = cluster_levels,
    Condition = c("Sed","Run"),
    fill = list(n = 0)
  ) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = n, values_fill = 0) %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  arrange(cluster) %>%
  cluster_fisher() %>%
  mutate(
    pct_Sed = 100 * Sed / total_sed,
    pct_Run = 100 * Run / total_run
  )

# sanity check: pct should vary across clusters
stopifnot(dplyr::n_distinct(tot_wide$pct_Sed) > 1, dplyr::n_distinct(tot_wide$pct_Run) > 1)

# ---- long format for paired bars ----
tot_pct_long <- tot_wide %>%
  dplyr::select(cluster, pct_Sed, pct_Run, sig) %>%
  tidyr::pivot_longer(cols = c(pct_Sed, pct_Run),
                      names_to = "Condition", values_to = "pct") %>%
  mutate(
    Condition = dplyr::recode(Condition, pct_Sed = "Sed", pct_Run = "Run"),
    cluster   = factor(cluster, levels = cluster_levels),
    Condition = factor(Condition, levels = c("Sed","Run"))
  )

# ---- sig star positions ----
tot_sig_pos <- tot_wide %>%
  transmute(
    cluster = factor(cluster, levels = cluster_levels),
    y = pmax(pct_Sed, pct_Run, na.rm = TRUE) + 1,
    sig
  )

# ---- plot ----
p_fig4c_pct <- ggplot(tot_pct_long, aes(x = cluster, y = pct, fill = Condition)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.65,
    color = "black",
    linewidth = 0.15
  ) +
  scale_fill_manual(values = c(Sed = "steelblue", Run = "firebrick"),
                    name = "Condition") +
  geom_text(
    data = tot_sig_pos,
    aes(x = cluster, y = y, label = sig),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Fig. 4c | Cluster composition differs between Sed and Run ECs",
    subtitle = "Bars = % of cells in each condition that fall in cluster k; per-cluster Fisher’s exact (in k vs not), BH-adjusted",
    x = "Seurat cluster",
    y = "% of cells (within condition)"
  ) +
  base_theme +
  theme(legend.position = "right")

print(p_fig4c_pct)

# ---- stats report ----
tot_report <- tot_wide %>%
  transmute(
    cluster,
    n_Sed = Sed, pct_Sed,
    n_Run = Run, pct_Run,
    delta_pct_Run_minus_Sed = pct_Run - pct_Sed,
    fisher_p, p_adj, sig
  )

print(tot_report)

#===============================================================
# Fig 4d (REVISED): % Bcl6b+ within each cluster (paired bars)
#   - y-axis: percent Bcl6b+ cells within that cluster
#   - stats: per-cluster Fisher test on 2x2 table (Run vs Sed) of Bcl6b+/-
#===============================================================

# 1) totals per cluster x condition
tot_counts <- md %>%
  dplyr::count(cluster, Condition, name = "n_total") %>%
  tidyr::complete(
    cluster = cluster_levels,
    Condition = c("Sed","Run"),
    fill = list(n_total = 0)
  )

# 2) Bcl6b+ counts per cluster x condition
pos_counts <- md %>%
  dplyr::filter(bcl_pos) %>%
  dplyr::count(cluster, Condition, name = "n_pos") %>%
  tidyr::complete(
    cluster = cluster_levels,
    Condition = c("Sed","Run"),
    fill = list(n_pos = 0)
  )

# 3) combine -> compute percent + run per-cluster Fisher on Bcl6b+/-
bcl_pct_wide <- tot_counts %>%
  left_join(pos_counts, by = c("cluster","Condition")) %>%
  mutate(
    n_pos = ifelse(is.na(n_pos), 0, n_pos),
    n_neg = pmax(n_total - n_pos, 0),
    pct_pos = ifelse(n_total > 0, 100 * n_pos / n_total, NA_real_)
  ) %>%
  tidyr::pivot_wider(
    names_from  = Condition,
    values_from = c(n_total, n_pos, n_neg, pct_pos),
    values_fill = 0
  ) %>%
  rowwise() %>%
  mutate(
    fisher_p = {
      # 2x2 within cluster:
      #          Bcl+   Bcl-
      # Run      n_pos_Run  n_neg_Run
      # Sed      n_pos_Sed  n_neg_Sed
      mat <- matrix(
        c(n_pos_Run, n_neg_Run,
          n_pos_Sed, n_neg_Sed),
        nrow = 2, byrow = TRUE
      )
      # If a cluster has zero cells in a condition, fisher is ill-defined;
      # return NA in that case.
      if (sum(mat[1, ]) == 0 || sum(mat[2, ]) == 0) NA_real_ else fisher.test(mat)$p.value
    }
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(fisher_p, method = "BH"),
    sig = case_when(
      is.na(p_adj)     ~ "",
      p_adj < 0.001    ~ "***",
      p_adj < 0.01     ~ "**",
      p_adj < 0.05     ~ "*",
      TRUE             ~ ""
    )
  ) %>%
  arrange(factor(cluster, levels = cluster_levels))

# 4) long format for plotting paired bars
bcl_pct_long <- bcl_pct_wide %>%
  dplyr::select(cluster, pct_pos_Sed, pct_pos_Run, sig) %>%
  tidyr::pivot_longer(
    cols = c(pct_pos_Sed, pct_pos_Run),
    names_to = "Condition",
    values_to = "pct_pos"
  ) %>%
  mutate(
    Condition = dplyr::recode(Condition,
                              pct_pos_Sed = "Sed",
                              pct_pos_Run = "Run"),
    cluster = factor(cluster, levels = cluster_levels),
    Condition = factor(Condition, levels = c("Sed","Run"))
  )

# where to place sig stars
bcl_pct_sig_pos <- bcl_pct_wide %>%
  transmute(
    cluster = factor(cluster, levels = cluster_levels),
    y = pmax(pct_pos_Sed, pct_pos_Run, na.rm = TRUE) + 2,
    sig
  )

p_fig4d_pct <- ggplot(bcl_pct_long,
                      aes(x = cluster, y = pct_pos, fill = Condition)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.65,
    color = "black",
    linewidth = 0.15
  ) +
  scale_fill_manual(
    values = c(Sed = "steelblue", Run = "firebrick"),
    name = "Condition"
  ) +
  geom_text(
    data = bcl_pct_sig_pos,
    aes(x = cluster, y = y, label = sig),
    inherit.aes = FALSE, size = 5
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Fig. 4d | % Bcl6b+ cells within each cluster",
    subtitle = "Bcl6b+ = Run-defined top 50% of nonzero; per-cluster Fisher’s exact on Bcl6b+/− (Run vs Sed), BH-adjusted",
    x = "Seurat cluster", y = "% Bcl6b+ (within cluster)"
  ) +
  base_theme +
  theme(legend.position = "right")

print(p_fig4d_pct)


# Optional: table to inspect
bcl_pct_report <- bcl_pct_wide %>%
  transmute(
    cluster,
    n_total_Sed, n_pos_Sed, pct_pos_Sed,
    n_total_Run, n_pos_Run, pct_pos_Run,
    fisher_p, p_adj
  )

bcl_pct_report















######################################################################
#//////////////Fig. 4h Test Pseudo-time code/////////////////////////#
#####################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

#------------------------------------------------------------
# USER SETTINGS
#------------------------------------------------------------
padj_cut <- 0.05
target_end_cluster <- "6"   # the cluster you want the lineage to include
use_reduction <- "PCA"      # recommended for slingshot; "UMAP" also possible

#------------------------------------------------------------
# 0) Run-only object
#------------------------------------------------------------
run_only <- subset(merged, subset = Condition == "Run")
run_only$seurat_clusters <- factor(run_only$seurat_clusters)

# sanity: does the Bcl6b_high flag exist?
stopifnot("Bcl6b_mod50_runCut" %in% colnames(run_only@meta.data))

#------------------------------------------------------------
# 1) Build the 108-gene Bcl6b DEG module score (if missing)
#------------------------------------------------------------
bcl6b_genes <- bcl6b_mod_DE %>%
  dplyr::filter(!is.na(p_val_adj), p_val_adj < padj_cut) %>%
  dplyr::pull(gene) %>%
  unique()

message("Bcl6b DEGs used for module score: ", length(bcl6b_genes))

if (!"Bcl6bDEG_Score1" %in% colnames(run_only@meta.data)) {
  feats <- intersect(bcl6b_genes, rownames(run_only))
  run_only <- AddModuleScore(run_only, features = list(feats), name = "Bcl6bDEG_Score")
}

#------------------------------------------------------------
# 2) Choose a start cluster = max fraction of Bcl6b_high cells
#------------------------------------------------------------
start_cluster_tbl <- tibble::tibble(
  seurat_clusters = as.character(run_only$seurat_clusters),
  Bcl6b_high      = as.logical(run_only$Bcl6b_mod50_runCut)
) %>%
  dplyr::count(seurat_clusters, Bcl6b_high, name = "n") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Bcl6b_high) %>%
  dplyr::arrange(dplyr::desc(frac))

start_cluster <- start_cluster_tbl %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::pull(seurat_clusters)

message("Start cluster (max Bcl6b_high fraction): ", start_cluster)
print(start_cluster_tbl)

#------------------------------------------------------------
# 3) Convert to SingleCellExperiment + run slingshot
#------------------------------------------------------------
sce <- as.SingleCellExperiment(run_only)
colData(sce)$cluster <- run_only$seurat_clusters

if (toupper(use_reduction) == "UMAP") {
  emb <- Embeddings(run_only, "umap")
  reducedDims(sce) <- SimpleList(UMAP = emb)
  rd_name <- "UMAP"
} else {
  emb <- Embeddings(run_only, "pca")[, 1:20, drop = FALSE]
  reducedDims(sce) <- SimpleList(PCA = emb)
  rd_name <- "PCA"
}

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim    = rd_name,
  start.clus    = start_cluster
)

#------------------------------------------------------------
# 4) Identify which lineage contains target_end_cluster
#------------------------------------------------------------
sds <- SlingshotDataSet(sce)
lineages <- slingLineages(sds)

message("Number of lineages inferred: ", length(lineages))
print(lineages)

lineage_idx <- which(vapply(lineages, function(x) target_end_cluster %in% x, logical(1)))

if (length(lineage_idx) == 0) {
  stop("No inferred lineage contains target_end_cluster = ", target_end_cluster,
       ". Try using PCA, changing resolution, or pick a different end cluster.")
}
if (length(lineage_idx) > 1) {
  message("Multiple lineages contain cluster ", target_end_cluster,
          ". Using the first one: ", lineage_idx[1])
}
lineage_idx <- lineage_idx[1]
message("Using lineage index: ", lineage_idx,
        " | clusters: ", paste(lineages[[lineage_idx]], collapse = " -> "))

#------------------------------------------------------------
# 5) Pull pseudotime for that lineage into Seurat
#------------------------------------------------------------
pt <- slingPseudotime(sce)
run_only$pseudotime_target <- as.numeric(pt[, lineage_idx])

# QC: which clusters are actually on this lineage (have pseudotime)?
qc_lineage <- tibble::tibble(
  cluster = as.character(run_only$seurat_clusters),
  pt      = run_only$pseudotime_target
) %>%
  dplyr::mutate(has_pt = !is.na(pt)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n = n(),
    n_with_pt = sum(has_pt),
    frac_with_pt = mean(has_pt),
    pt_median = median(pt, na.rm = TRUE),
    pt_q10 = quantile(pt, 0.10, na.rm = TRUE),
    pt_q90 = quantile(pt, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(frac_with_pt), desc(n_with_pt))

print(qc_lineage)

#------------------------------------------------------------
# 6) Build plotting df on UMAP (only cells on this lineage)
#------------------------------------------------------------
um <- Embeddings(run_only, "umap")
df_umap <- as.data.frame(um) %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    cluster = as.character(run_only$seurat_clusters[cell]),
    pt      = run_only$pseudotime_target[cell],
    Bcl6b   = FetchData(run_only, "Bcl6b")[cell, 1],
    mod     = run_only$Bcl6bDEG_Score1[cell]
  ) %>%
  dplyr::filter(!is.na(pt)) %>%
  dplyr::arrange(pt)

umap_x <- colnames(um)[1]
umap_y <- colnames(um)[2]
#------------------------------------------------------------
# 7a) Plots
#------------------------------------------------------------

#======================================================================
# CLEAN UMAP TRACE for lineage 2 (1 -> 2 -> 0 -> 6)
#   Option A: bin-centroid polyline (robust)
#======================================================================
#======================================================================
# EXACT-COLOR lineage trace on UMAP (colors match Seurat DimPlot exactly)
#   - uses per-cell hex colors extracted from DimPlot
#   - highlights clusters 1/2/0/6, others grey
#   - centroid trace + “1 → 2 → 0 → 6” labels at segment medians
#======================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(grid)
})

clusters_keep <- c("1","2","0","6")
n_bins <- 35

#------------------------------------------------------------
# 0) Run-only + trusted DimPlot (color "source of truth")
#------------------------------------------------------------
run_only <- subset(merged, subset = Condition == "Run")
run_only$seurat_clusters <- factor(run_only$seurat_clusters)

p_umap_seurat <- DimPlot(run_only, reduction = "umap", group.by = "seurat_clusters")
gb <- ggplot_build(p_umap_seurat)

# Pull the exact point coordinates & grouping from the DimPlot data
df_dim <- p_umap_seurat$data %>%
  tibble::as_tibble() %>%
  mutate(.row = dplyr::row_number())

# Pull the exact colors used for each point (same row order)
df_col <- gb$data[[1]] %>%
  tibble::as_tibble() %>%
  transmute(.row = dplyr::row_number(), hex = colour)

# Identify which column contains the grouping labels
group_col <- dplyr::case_when(
  "seurat_clusters" %in% colnames(df_dim) ~ "seurat_clusters",
  "ident"           %in% colnames(df_dim) ~ "ident",
  "group"           %in% colnames(df_dim) ~ "group",
  TRUE ~ NA_character_
)

if (is.na(group_col)) {
  stop("Couldn't find cluster labels in DimPlot data. Columns are:\n",
       paste(colnames(df_dim), collapse = ", "))
}

# Identify UMAP coordinate columns robustly (Seurat usually uses UMAP_1/UMAP_2)
umap_cols <- c("UMAP_1", "UMAP_2")
if (!all(umap_cols %in% colnames(df_dim))) {
  # fallback: take the first two numeric columns if UMAP_1/2 aren't present
  num_cols <- names(df_dim)[vapply(df_dim, is.numeric, logical(1))]
  if (length(num_cols) < 2) {
    stop("Couldn't find UMAP columns. Numeric columns are:\n",
         paste(num_cols, collapse = ", "))
  }
  umap_cols <- num_cols[1:2]
}

umap_x <- umap_cols[1]
umap_y <- umap_cols[2]

df_all <- df_dim %>%
  left_join(df_col, by = ".row") %>%
  mutate(
    cluster = as.character(.data[[group_col]])
  )

#------------------------------------------------------------
# 1) df_umap must exist already from your slingshot pipeline:
#    it should contain: cell, cluster, pt, and UMAP coords matching run_only
#------------------------------------------------------------
stopifnot(exists("df_umap"))
stopifnot(all(c("pt","cluster") %in% colnames(df_umap)))

df_path <- df_umap %>%
  tibble::as_tibble() %>%
  mutate(
    x = .data[[umap_x]],
    y = .data[[umap_y]]
  ) %>%
  filter(!is.na(pt), is.finite(pt), is.finite(x), is.finite(y)) %>%
  arrange(pt)

#------------------------------------------------------------
# 2) Bin pseudotime -> centroid trace + segment labels
#------------------------------------------------------------
clusters_keep <- c("1","2","0","6")
n_bins <- 35

df_bins <- df_path %>%
  mutate(pt_bin = cut(pt, breaks = n_bins, include.lowest = TRUE)) %>%
  group_by(pt_bin) %>%
  summarise(
    pt_mid = median(pt),
    x = mean(x),
    y = mean(y),
    .groups = "drop"
  ) %>%
  arrange(pt_mid)

bin_map <- df_path %>%
  mutate(pt_bin = cut(pt, breaks = n_bins, include.lowest = TRUE)) %>%
  group_by(pt_bin) %>%
  summarise(
    pt_mid = median(pt),
    cluster_seg = names(which.max(table(cluster))),
    .groups = "drop"
  )

df_bins2 <- df_bins %>%
  left_join(bin_map, by = c("pt_bin","pt_mid")) %>%
  mutate(cluster_seg = as.character(cluster_seg)) %>%
  filter(cluster_seg %in% clusters_keep)

label_df <- df_bins2 %>%
  group_by(cluster_seg) %>%
  summarise(
    x = median(x),
    y = median(y),
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      cluster_seg == "1" ~ "1",
      cluster_seg == "2" ~ "→ 2",
      cluster_seg == "0" ~ "→ 0",
      cluster_seg == "6" ~ "→ 6",
      TRUE ~ cluster_seg
    )
  )

#------------------------------------------------------------
# 3) Plot: EXACT SAME colors as DimPlot (per-cell hex)
#------------------------------------------------------------
p_exact <- ggplot(df_all, aes(x = .data[[umap_x]], y = .data[[umap_y]])) +
  geom_point(color = "grey85", size = 0.35, alpha = 0.55) +
  geom_point(
    data = df_all %>% filter(cluster %in% clusters_keep),
    aes(color = hex),
    size = 0.50,
    alpha = 0.95
  ) +
  scale_color_identity(guide = "none") +
  geom_path(
    data = df_bins,
    aes(x = x, y = y),
    linewidth = 1.2,
    color = "black",
    lineend = "round",
    arrow = arrow(type = "closed", length = unit(0.12, "in"))
  ) +
  geom_point(
    data = df_bins %>% slice_head(n = 1),
    aes(x = x, y = y),
    size = 2.2, shape = 21, fill = "white", color = "black", stroke = 0.8
  ) +
  geom_point(
    data = df_bins %>% slice_tail(n = 1),
    aes(x = x, y = y),
    size = 2.6, shape = 21, fill = "black", color = "black", stroke = 0.8
  ) +
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4.2,
    linewidth = 0.25,
    label.padding = unit(0.12, "lines"),
    fill = "white",
    color = "black"
  ) +
  labs(
    title = "RUN ECs: lineage 2 traced on UMAP",
    subtitle = "Cluster colors are taken directly from Seurat DimPlot (exact match)",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

print(p_exact)

ggsave("Fig4h_pseudotime_UMAP.pdf", p_exact, width = 4.5, height = 4.2, units = "in")

# or PNG if needed
ggsave("Fig4h_pseudotime_UMAPh.png", p_exact, width = 4.5, height = 4.2, units = "in", dpi = 400)



#//////////other plots for pseudotime///////////////////////#
#------------------------------------------------------------
# 7b) Plots
#------------------------------------------------------------

# A) cells on the lineage, colored by pseudotime
p_lineage_umap <- ggplot(df_umap, aes(x = .data[[umap_x]], y = .data[[umap_y]])) +
  geom_point(aes(color = pt), size = 0.55, alpha = 0.75) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(
    title = paste0("RUN ECs: lineage containing cluster ", target_end_cluster),
    subtitle = paste0("Lineage ", lineage_idx, ": ",
                      paste(lineages[[lineage_idx]], collapse = " -> "),
                      " | start = ", start_cluster),
    x = "UMAP 1", y = "UMAP 2", color = "Pseudotime"
  )

# B) Bcl6b vs pseudotime along this lineage
p_bcl6b_pt <- ggplot(df_umap, aes(x = pt, y = Bcl6b)) +
  geom_point(size = 0.6, alpha = 0.20) +
  geom_smooth(se = TRUE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Bcl6b expression vs pseudotime (cells on selected lineage)",
    x = paste0("Pseudotime (Lineage ", lineage_idx, ")"),
    y = "Bcl6b (log-normalized)"
  )

# C) module vs pseudotime along this lineage
p_mod_pt <- ggplot(df_umap, aes(x = pt, y = mod)) +
  geom_point(size = 0.6, alpha = 0.20) +
  geom_smooth(se = TRUE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Bcl6b-DEG module score vs pseudotime (cells on selected lineage)",
    x = paste0("Pseudotime (Lineage ", lineage_idx, ")"),
    y = "Bcl6bDEG_Score1"
  )

# D) cluster identity along pseudotime (sanity check)
p_cluster_pt <- ggplot(df_umap, aes(x = pt, y = cluster, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.35) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "Cluster membership along pseudotime (selected lineage)",
    x = paste0("Pseudotime (Lineage ", lineage_idx, ")"),
    y = "Seurat cluster"
  )

print(p_lineage_umap)
print(p_bcl6b_pt)
print(p_mod_pt)
print(p_cluster_pt)























###################################################################
###################################################################
#################################################################
#//////Likely do not need this///////////////////////#
###################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(ggplot2)
})

#---------------------------------------
# Inputs
#---------------------------------------
padj_cut <- 0.05
target_end_cluster <- "6"   # set to NULL if you don't want to force an end state
use_reduction <- "PCA"      # "PCA" is recommended; you can also try "UMAP"

#---------------------------------------
# Run-only + module score
#---------------------------------------
run_only <- subset(merged, subset = Condition == "Run")
run_only$seurat_clusters <- factor(run_only$seurat_clusters)

bcl6b_genes <- bcl6b_mod_DE %>%
  dplyr::filter(!is.na(p_val_adj), p_val_adj < padj_cut) %>%
  dplyr::pull(gene) %>%
  unique()

if (!"Bcl6bDEG_Score1" %in% colnames(run_only@meta.data)) {
  feats <- intersect(bcl6b_genes, rownames(run_only))
  run_only <- AddModuleScore(run_only, features = list(feats), name = "Bcl6bDEG_Score")
}

#---------------------------------------
# Start cluster = max fraction Bcl6b_high cells
# (robust: build a plain tibble, avoid meta.data/list weirdness)
#---------------------------------------

# sanity: does the flag exist?
stopifnot("Bcl6b_mod50_runCut" %in% colnames(run_only@meta.data))

start_cluster_tbl <- tibble::tibble(
  seurat_clusters = as.character(run_only$seurat_clusters),
  Bcl6b_high      = as.logical(run_only$Bcl6b_mod50_runCut)
) %>%
  dplyr::count(seurat_clusters, Bcl6b_high, name = "n") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Bcl6b_high) %>%
  dplyr::arrange(dplyr::desc(frac))

print(start_cluster_tbl)

start_cluster <- start_cluster_tbl %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::pull(seurat_clusters)

message("Start cluster: ", start_cluster)

#---------------------------------------
# Build SCE + choose embedding
#---------------------------------------
sce <- as.SingleCellExperiment(run_only)
colData(sce)$cluster <- run_only$seurat_clusters

if (toupper(use_reduction) == "UMAP") {
  emb <- Embeddings(run_only, "umap")
  reducedDims(sce) <- SimpleList(UMAP = emb)
  rd_name <- "UMAP"
} else {
  emb <- Embeddings(run_only, "pca")[, 1:20, drop = FALSE]
  reducedDims(sce) <- SimpleList(PCA = emb)
  rd_name <- "PCA"
}

#---------------------------------------
# Run slingshot (optionally force end cluster)
#---------------------------------------
if (!is.null(target_end_cluster)) {
  sce <- slingshot(
    sce,
    clusterLabels = "cluster",
    reducedDim = rd_name,
    start.clus = start_cluster,
    end.clus   = target_end_cluster
  )
  message("Forced end cluster: ", target_end_cluster)
} else {
  sce <- slingshot(
    sce,
    clusterLabels = "cluster",
    reducedDim = rd_name,
    start.clus = start_cluster
  )
}

pt <- slingPseudotime(sce)
run_only$pseudotime_L1 <- as.numeric(pt[, 1])

suppressPackageStartupMessages({
  library(dplyr)
})

qc_lineage <- tibble::tibble(
  cluster = as.character(run_only$seurat_clusters),
  pt      = run_only$pseudotime_L1
) %>%
  mutate(
    has_pt = !is.na(pt)
  ) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    n_with_pt = sum(has_pt),
    frac_with_pt = mean(has_pt),
    pt_median = median(pt, na.rm = TRUE),
    pt_q10 = quantile(pt, 0.10, na.rm = TRUE),
    pt_q90 = quantile(pt, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(frac_with_pt), desc(n_with_pt))

print(qc_lineage)

# specifically inspect cluster 6
qc_lineage %>% filter(cluster == "6")

# how many lineages?
sds <- SlingshotDataSet(sce)
length(slingLineages(sds))
slingLineages(sds)  # prints clusters in each lineage

# pseudotime matrix: one column per lineage
pt <- slingPseudotime(sce)
dim(pt)

# for each lineage, does cluster 6 have non-NA pseudotime?
for (L in seq_len(ncol(pt))) {
  has6 <- any(!is.na(pt[as.character(run_only$seurat_clusters) == "6", L]))
  message("Lineage ", L, " includes cluster 6? ", has6)
}

