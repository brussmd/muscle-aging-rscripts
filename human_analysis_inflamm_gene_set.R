#Attempt at human data

human_aging_counts <- read.delim("GSE242202_Read_counts.txt", 
                     header = TRUE,       # first row is column names
                     row.names = 1,       # use first column (genes) as row names
                     check.names = FALSE) # keeps original column names

head(human_aging_counts)

library(edgeR)
library(dplyr)
library(biomaRt)

# 1. Extract count matrix
gene_info <- human_aging_counts[,1:2]  # gene_name and gene_biotype
count_matrix <- human_aging_counts[,-c(1,2)]
rownames(count_matrix) <- rownames(human_aging_counts)

# 2. Build group info
sample_groups <- ifelse(grepl("^Y", colnames(count_matrix)), "Young", "Old")
sample_groups <- factor(sample_groups, levels = c("Young", "Old"))
# Ensure group is a factor and relevel it BEFORE design
sample_group$group <- relevel(factor(sample_group$group), ref = "Young")

# 3. Create DGEList
dge <- DGEList(counts = count_matrix, genes = gene_info)
dge <- calcNormFactors(dge)

design <- model.matrix(~ sample_groups)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)   # Old vs Young

# Get results
deg_results <- topTags(lrt, n = Inf)$table
deg_results <- deg_results %>% rownames_to_column("ensembl_gene_id")

# Connect to Ensembl archive (Dec 2021)
human <- useMart("ensembl",
                 dataset = "hsapiens_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org")

mouse <- useMart("ensembl",
                 dataset = "mmusculus_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org")

# Map Ensembl human gene IDs to mouse gene symbols
orthologs <- getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = deg_results$ensembl_gene_id,
                    mart = human,
                    attributesL = c("ensembl_gene_id","mgi_symbol"),
                    martL = mouse,
                    uniqueRows = TRUE)

colnames(orthologs) <- c("ensembl_gene_id", "hgnc_symbol", 
                         "ensembl_gene_id_mouse", "mgi_symbol")

deg_mouse <- merge(deg_results, orthologs, by = "ensembl_gene_id")

# Merge with DEG results
deg_mouse_ferrucci <- merge(deg_results, orthologs, by = "ensembl_gene_id")


deg_mouse

validated_inflamm_geneset_upper.vec <- c(
  "PRKCZ","JCHAIN","ANKRD1","HIP1R","KNG2","CES2C","UVRAG","FOSL2","GAS6","TNFRSF23",
  "MID2","CXCL10","CX3CL1","IL20RB","RAPGEF1","CYP4F14","TRP63","POU4F1","MYC","RAB20",
  "RUNX1","CASP4","LGALS1","DLL1","CCL8","SFRP1","CD55","BAX","POU2F2","CCL7",
  "HSP90AA1","ARID5A","BBC3","MIF","RELB","NINJ1","NR1H3","BID","CDKN1A","SOX4",
  "SBNO2","CCL2","ICAM1","TUBB5","APOD","AKT1","SLC15A4","FAM110A","GPX1","TYK2",
  "PPARD","FOS"
)

validated_inflamm_geneset.vec

deg_mouse_ferrucci %>%
  filter(mgi_symbol %in% validated_inflamm_geneset.vec) %>%
  arrange(logFC)


library(clusterProfiler)
library(dplyr)

# ---- 1. Create ranked gene list ----
geneList <- deg_mouse %>%
  mutate(stat = -log10(PValue) * sign(logFC)) %>%
  arrange(desc(stat)) %>%
  { setNames(.$stat, .$mgi_symbol) }

# ---- 2. Run GSEA with your inflammaging gene set ----
gsea_inflamm_human <- GSEA(geneList = geneList,
                           TERM2GENE = data.frame(term = "Inflamm", 
                                                  gene = validated_inflamm_geneset.vec),
                           pvalueCutoff = 1,
                           verbose = FALSE)

# ---- 3. View results ----
gsea_inflamm_human@result %>%
  dplyr::select(ID, NES, pvalue, p.adjust)
#-----------------------------------------------------------#

# ---- 1. Rank genes by p-value ----
deg_mouse <- deg_mouse %>%
  arrange(PValue) %>%
  mutate(rank = rank(PValue))

# ---- 2. Separate ranks for inflammaging genes vs all others ----
inflammaging_ranks <- deg_mouse$rank[deg_mouse$mgi_symbol %in% validated_inflamm_geneset.vec]
background_ranks   <- deg_mouse$rank[!deg_mouse$mgi_symbol %in% validated_inflamm_geneset.vec]

# ---- 3. Wilcoxon test (one-sided: are inflammaging genes more significant?) ----
wilcox_test <- wilcox.test(inflammaging_ranks, background_ranks, alternative = "less")

# ---- 4. Output results ----
wilcox_test
median(inflammaging_ranks)
median(background_ranks)

rbc <- 1 - (2 * wilcox_test$statistic / (length(inflammaging_ranks) * length(background_ranks)))
rbc
#------------------------------------------------------------#

# ---- 1. Flag significant genes (FDR < 0.05) ----
deg_mouse <- deg_mouse %>%
  mutate(sig = FDR < 0.1,
         inflammaging = mgi_symbol %in% validated_inflamm_geneset.vec)

# ---- 2. Contingency table ----
table_data <- table(deg_mouse$inflammaging, deg_mouse$sig)
colnames(table_data) <- c("Not_Significant","Significant")
rownames(table_data) <- c("Not_Inflammaging","Inflammaging")

# ---- 3. Fisher's exact test (one-sided, enrichment) ----
fisher_result <- fisher.test(table_data, alternative = "greater")

# ---- 4. Percentages ----
pct_inflammaging_sig <- 100 * table_data["Inflammaging","Significant"] /
  sum(table_data["Inflammaging",])
pct_background_sig   <- 100 * table_data["Not_Inflammaging","Significant"] /
  sum(table_data["Not_Inflammaging",])

# ---- 5. Output summary ----
list(
  contingency_table = table_data,
  fisher_test = fisher_result,
  pct_inflammaging_sig = pct_inflammaging_sig,
  pct_background_sig   = pct_background_sig
)

dim(deg_mouse)
#------------------------------------------------------------#

#////////////////////////////////////////////////////////////#
#-------Robinson MM 2017 Re-Analysis-------------------------#
#////////////////////////////////////////////////////////////#

# Tab-delimited, gene IDs as rownames
human_aging_counts_robinson_1 <- read.delim("GSE97084_GeneCount_raw.tsv",
                                          header = TRUE,
                                          check.names = FALSE)

# Inspect what you loaded
head(human_aging_counts_robinson_1)
dim(human_aging_counts_robinson_1)

# Tab-delimited, gene IDs as rownames
human_aging_counts_robinson_2 <- read.delim("GSE97084_GeneCount_raw_2.tsv",
                                            header = TRUE,
                                            check.names = FALSE)

# Inspect what you loaded
head(human_aging_counts_robinson_2)
dim(human_aging_counts_robinson_2)

colnames(human_aging_counts_robinson_2)

# Remove columns with empty names
human_aging_counts_robinson_2 <- human_aging_counts_robinson_2[, colnames(human_aging_counts_robinson_2) != ""]

# Now safely select
df1 <- human_aging_counts_robinson_1 %>%
  dplyr::select(-Chr, -Start, -Stop, -CodingLength)

df2 <- human_aging_counts_robinson_2 %>%
  dplyr::select(-Chr, -Start, -Stop, -CodingLength)

# Merge
human_aging_counts_robinson_merged <- full_join(df1, df2, by = c("GeneID","GeneName"))
# ---- 3. Check the result ----
dim(human_aging_counts_robinson_merged)
head(human_aging_counts_robinson_merged)

# ---- 1. Define keys ----
young_ids <- c("2B","3A","6A","10B","11A","13A","14B","15A","16A","17A","18A",
               "19A","23A","24A","25A","26A","27B","29A","30A","31B","32A","33A",
               "35A","36B","38A","39A","40A","41B")

old_ids <- c("1B","4A","5B","7B","8A","9A","20B","21B","22A","28A","34A","37B",
             "42A","43A","44A","45A","46B","47A","48B","49A","50A","51A","52A","53A")

# ---- 2. Extract sample ID from column names ----
# Keep GeneID and GeneName as identifiers
meta_cols <- c("GeneID", "GeneName")

sample_cols <- setdiff(colnames(human_aging_counts_robinson_merged), meta_cols)

# Extract the middle part (like "51A") from sample names
sample_ids <- sub("^s_([^-]+)-.*", "\\1", sample_cols)

# ---- 3. Identify which columns are young vs old ----
young_cols <- sample_cols[sample_ids %in% young_ids]
old_cols   <- sample_cols[sample_ids %in% old_ids]

# ---- 4. Keep only relevant columns ----
counts_filtered <- human_aging_counts_robinson_merged %>%
  dplyr::select(all_of(c(meta_cols, young_cols, old_cols)))

# ---- 5. Optional: create group assignment dataframe ----
sample_group <- data.frame(
  sample = c(young_cols, old_cols),
  group  = c(rep("Young", length(young_cols)), rep("Old", length(old_cols)))
)

counts_filtered

# Get current column names
cn <- colnames(counts_filtered)

# Simplify sample names (but keep GeneID and GeneName)
cn_new <- cn
sample_mask <- !(cn %in% c("GeneID", "GeneName"))
cn_new[sample_mask] <- sub("^s_([^-]+)-.*", "\\1", cn[sample_mask])

# Apply new column names
colnames(counts_filtered) <- cn_new

cn_new
counts_filtered

# Create metadata (Young vs Old)
sample_group <- data.frame(
  sample = cn_new[!(cn_new %in% c("GeneID","GeneName"))],
  group  = ifelse(cn_new[!(cn_new %in% c("GeneID","GeneName"))] %in% young_ids, "Young", "Old")
)

# Build count matrix
count_matrix <- counts_filtered %>%
  dplyr::select(-GeneName) %>%
  column_to_rownames("GeneID")

counts_filtered_unique <- counts_filtered %>%
  group_by(GeneID) %>%
  summarise(across(-GeneName, sum), .groups = "drop")  # sum counts for duplicates

count_matrix <- counts_filtered_unique %>%
  column_to_rownames("GeneID")

count_matrix[is.na(count_matrix)] <- 0

count_matrix

# edgeR pipeline
library(edgeR)
# Ensure group is a factor and relevel it BEFORE design
sample_group$group <- relevel(factor(sample_group$group), ref = "Young")

# DGEList
dge <- DGEList(counts = count_matrix, group = sample_group$group)
dge <- calcNormFactors(dge)

# Design matrix with Young as reference
design <- model.matrix(~ group, data = sample_group)

# Estimate dispersion and fit model
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

# LRT: coef=2 now represents Old vs Young (positive logFC = Old > Young)
lrt <- glmLRT(fit, coef = 2)

topTags(lrt)




# Extract all results as dataframe
deg_results <- topTags(lrt, n = Inf)$table
deg_results <- tibble::rownames_to_column(deg_results, var = "ensembl_gene_id")

head(deg_results)

# Connect to Ensembl archive
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org")

# Map human Ensembl gene IDs to mouse gene symbols
orthologs <- getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = deg_results$ensembl_gene_id,
                    mart = human,
                    attributesL = c("ensembl_gene_id","mgi_symbol"),
                    martL = mouse,
                    uniqueRows = TRUE)

colnames(orthologs) <- c("ensembl_gene_id", "hgnc_symbol",
                         "ensembl_gene_id_mouse", "mgi_symbol")

# Merge with DEG results
deg_mouse_robinson <- merge(deg_results, orthologs, by = "ensembl_gene_id")

head(deg_mouse_robinson)

deg_mouse_robinson %>%
  pull(mgi_symbol)

dim(deg_mouse_robinson)

deg_mouse_robinson %>%
  filter(FDR <0.05 & logFC > 0)
#--------------------------------------------#

#============================================#
#-----Robinson Dotplot-----------------------#
#============================================#
validated_inflamm_geneset.vec

# Step 1: Prepare master gene list
gene_df <- tibble(mgi_symbol = sort(validated_inflamm_geneset.vec))

# Step 2: Reduce deg_mouse_robinson to unique best hits
deg_filtered <- deg_mouse_robinson %>%
  filter(mgi_symbol %in% validated_inflamm_geneset.vec) %>%
  group_by(mgi_symbol) %>%
  slice_min(order_by = FDR, n = 1) %>%
  ungroup()

# Step 3: Join to get all 51 genes with data (or NA if missing)
plot_df <- gene_df %>%
  left_join(deg_filtered, by = "mgi_symbol") %>%
  mutate(
    column = case_when(
      is.na(logFC) ~ NA_character_,
      logFC < 0 ~ "Negative logFC",
      logFC > 0 ~ "Positive logFC",
      TRUE ~ NA_character_
    ),
    color = case_when(
      is.na(logFC) ~ NA_character_,
      logFC < 0 & PValue < 0.05 ~ "blue",
      logFC < 0 & PValue >= 0.05 ~ "gray",
      logFC > 0 & PValue < 0.05 ~ "red",
      logFC > 0 & PValue >= 0.05 ~ "gray",
      TRUE ~ NA_character_
    ),
    size = ifelse(color == "gray", 2, 3)
  )

# Flip the x-axis: Positive on the left (0.8), Negative on the right (1.2)
plot_df_filtered <- plot_df %>%
  mutate(
    x = ifelse(column == "Positive logFC", 0.6, 0.7)
  )

# Updated dotplot
robinson_dotplot <- ggplot(plot_df_filtered, aes(x = x, y = factor(mgi_symbol, levels = rev(sort(mgi_symbol))))) +
  geom_point(aes(color = color, size = size), na.rm = TRUE) +
  scale_color_identity() +
  scale_size_identity() +
  scale_x_continuous(
    breaks = c(0.8, 1.2),
    limits = c(0.6, 1.4)  # tight bounds
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 10)
  )

# Save as PDF
pdf("robinson_dotplot.pdf", width = 2.5, height = 2.5)
print(robinson_dotplot)
dev.off()
#------------------------------------------#

#-------------------------------------------#
#------GO Enrichment----------------------#
#-------------------------------------------#

# ---- 1. Define gene sets ----
overlap_genes <- c("Nr1h3", "Pou2f2", "Trp63", "Fosl2", "Mid2", 
                   "Sfrp1", "Bbc3", "Tnfrsf23", "Cdkn1a", "Jchain", 
                   "Rab20", "Apod", "Arid5a") %>% unique()

background_genes <- deg_mouse_robinson %>% pull(mgi_symbol) %>% unique()

# ---- 2. Convert to Entrez IDs ----
library(org.Mm.eg.db)
overlap_entrez <- mapIds(org.Mm.eg.db, keys = overlap_genes,
                         column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
background_entrez <- mapIds(org.Mm.eg.db, keys = background_genes,
                            column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Remove NAs
overlap_entrez <- overlap_entrez[!is.na(overlap_entrez)]
background_entrez <- background_entrez[!is.na(background_entrez)]

# ---- 3. Run GO Biological Process enrichment ----
library(clusterProfiler)
ego_bp <- enrichGO(gene         = overlap_entrez,
                   universe     = background_entrez,
                   OrgDb        = org.Mm.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",
                   pAdjustMethod= "BH",
                   qvalueCutoff = 1,
                   readable     = TRUE)

# ---- 4. View results ----
ego_bp@result %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust)

# ---- 5. Optional: Dotplot ----
library(enrichplot)
dotplot(ego_bp, showCategory = 10) + ggtitle("GO Biological Process Enrichment")

#-----------------------------------------#
#----Fisher Exact Test---------------------#
#------------------------------------------#

# 1. Flag significant genes (FDR < 0.05) and membership in inflammaging set
deg_mouse_fisher <- deg_mouse_robinson %>%
  mutate(sig = FDR < 0.05,
         inflammaging = mgi_symbol %in% validated_inflamm_geneset.vec)

# 2. Build contingency table (rows: inflammaging membership, cols: DEG significance)
table_data <- table(deg_mouse_fisher$inflammaging, deg_mouse_fisher$sig)

# 3. Rename for readability
dimnames(table_data) <- list(
  Inflammaging = c("No","Yes"),
  DEG = c("Not_Significant","Significant")
)

# 4. Fisher's exact test (one-sided: enrichment of DEGs in inflammaging genes)
fisher_result <- fisher.test(table_data, alternative = "greater")

# 5. Percentages for context
pct_inflammaging_sig <- 100 * table_data["Yes","Significant"] /
  sum(table_data["Yes",])
pct_background_sig   <- 100 * table_data["No","Significant"] /
  sum(table_data["No",])

# 6. Output summary
list(
  contingency_table = table_data,
  fisher_test = fisher_result,
  pct_inflammaging_sig = pct_inflammaging_sig,
  pct_background_sig   = pct_background_sig
)
#-----------------------------------------------#



#--------Visualize Fisher-----------------------#

library(ggplot2)
library(dplyr)

# ---- Input counts ----
inflammaging_sig    <- 14
inflammaging_total  <- 53
background_sig      <- 3817
background_total    <- 24699

# ---- Calculate percentages ----
pct_inflammaging <- 100 * inflammaging_sig / inflammaging_total
pct_background   <- 100 * background_sig / background_total

# ---- Prepare data for plotting ----
dot_data <- data.frame(
  Category = c("Inflammaging gene set", "Background"),
  Percent  = c(pct_inflammaging, pct_background)
)

# ---- Dot plot ----
ggplot(dot_data, aes(x = Category, y = Percent)) +
  geom_point(size = 5, color = "black") +
  geom_segment(aes(xend = Category, y = 0, yend = Percent), linetype = "dashed", color = "grey50") +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            vjust = -1, size = 5) +
  theme_minimal() +
  labs(y = "% of genes that are significant (FDR < 0.05)",
       x = NULL,
       title = "Enrichment of DEGs in Inflammaging Gene Set vs Background") +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
#-----------------------------------------------#


#-----------------------------------------------#
#--------Directionality Fisher------------------#
#-----------------------------------------------#

# 1. Subset to significant DEGs in inflammaging set
inflam_sig <- deg_mouse_robinson %>%
  filter(FDR < 0.05 & mgi_symbol %in% validated_inflamm_geneset.vec)

# 2. Count up- vs downregulated
table_direction <- table(ifelse(inflam_sig$logFC > 0, "Up", "Down"))

# 3. Compare to background DEGs (optional enrichment check)
background_sig <- deg_mouse_robinson %>%
  filter(FDR < 0.05 & !mgi_symbol %in% validated_inflamm_geneset.vec)

# Compare proportion upregulated
prop_up_inflam <- mean(inflam_sig$logFC > 0)
prop_up_bg     <- mean(background_sig$logFC > 0)

# Fisher's test for directionality bias
direction_table <- matrix(c(
  sum(inflam_sig$logFC > 0), sum(inflam_sig$logFC <= 0),
  sum(background_sig$logFC > 0), sum(background_sig$logFC <= 0)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Inflammaging","Background"), c("Up","Down"))
)
fisher_dir <- fisher.test(direction_table, alternative = "greater")

list(
  table_direction = table_direction,
  prop_up_inflam = prop_up_inflam,
  prop_up_bg = prop_up_bg,
  fisher_dir = fisher_dir
)
#------------------------------------------#


#------------------------------------------#
#-------GSEA-------------------------------#
#------------------------------------------#
# ---- 1. Create ranked gene list ----
geneList <- deg_mouse_robinson %>%
  mutate(stat = -log10(PValue) * logFC) %>%
  arrange(desc(stat)) %>%
  { setNames(.$stat, .$mgi_symbol) }
#------------------------------------------#
#-------alt ranking------------------------#

geneList <- deg_mouse_robinson %>%
  mutate(stat = sign(logFC) * sqrt(LR)) %>%
  arrange(desc(stat)) %>%
  distinct(mgi_symbol, .keep_all = TRUE) %>% # keep unique symbols
  { setNames(.$stat, .$mgi_symbol) }


# ---- 2. Run GSEA with your inflammaging gene set ----
gsea_inflamm_human_robinson <- GSEA(geneList = geneList,
                           TERM2GENE = data.frame(term = "Inflamm", 
                                                  gene = validated_inflamm_geneset.vec),
                           pvalueCutoff = 1,
                           verbose = FALSE)

# ---- 3. View results ----
gsea_inflamm_human_robinson@result %>%
  dplyr::select(ID, NES, pvalue, p.adjust)
#------------------------------------------------#
#------------------------------------------------#

#------Visualize GSEA-----------------------------#

install.packages("cowplot")
library(cowplot)
library(ggplot2)

# ---- 1. Ranked stat values (your precomputed geneList) ----
# geneList should be named vector: names = gene symbols, values = stat
stats <- geneList
pathway <- intersect(names(stats), validated_inflamm_geneset.vec)

# ---- 2. Compute running enrichment score ----
N <- length(stats)
hit <- names(stats) %in% pathway
Nh <- sum(hit)
Nm <- N - Nh
P <- sum(abs(stats[hit]))

runningES <- cumsum(ifelse(hit, abs(stats) / P, -1 / Nm))

# ---- 3. Dataframes for plotting ----
df_top <- data.frame(Rank = seq_along(stats), ES = runningES)
df_hits <- data.frame(Rank = which(hit))

df_bottom <- data.frame(
  Rank = seq_along(stats),
  Stat = stats
)

# ---- 4. Top panel (Enrichment Score curve) ----
p_top <-ggplot(df_top, aes(x = Rank, y = ES)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(data = df_top[which.max(abs(df_top$ES)), ],
             aes(x = Rank, y = ES), color = "red", size = 3) +
  labs(y = "Enrichment Score (ES)", x = NULL,
       title = "GSEA Enrichment Curve") +
  theme_minimal(base_size = 14)

# ---- 5. Middle panel (hit positions) ----
p_mid <-ggplot(df_hits, aes(x = Rank, y = 1)) +
  geom_segment(aes(xend = Rank, yend = 0), color = "black") +
  labs(y = NULL, x = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = -20, b = -20))

# ---- 6. Bottom panel (stat values with color scale) ----
p_bottom <-ggplot(df_bottom, aes(x = Rank, y = Stat, fill = Stat)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, name = "-log10(P)*sign(logFC)") +
  labs(x = "Gene Rank", y = "Ranked Statistic") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# ---- 7. Combine panels (cowplot) ----
plot_combined <- plot_grid(
  p_top,
  p_mid,
  p_bottom,
  ncol = 1,
  align = "v",
  rel_heights = c(3, 0.5, 2)
)

# ---- 8. Display ----
print(plot_combined)

# Save as PDF
pdf("plot_combined.pdf", width = 8, height = 8)  # Adjust height
print(plot_combined)
dev.off()






#---Clustering test------------------------------#

#---- 1. Map Ensembl IDs to human gene symbols ----
human_symbol_map <- orthologs %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol) %>%
  distinct()

count_symbol <- count_matrix %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  left_join(human_symbol_map, by = "ensembl_gene_id") %>%
  filter(!is.na(hgnc_symbol))

# Aggregate duplicate symbols (if any)
count_symbol <- count_symbol %>%
  group_by(hgnc_symbol) %>%
  summarise(across(-ensembl_gene_id, sum), .groups = "drop") %>%
  tibble::column_to_rownames("hgnc_symbol")

#---- 2. Subset to inflammaging genes ----
inflamm_genes <- toupper(validated_inflamm_geneset.vec)
genes_present <- intersect(inflamm_genes, rownames(count_symbol))
expr_subset <- count_symbol[genes_present, ]

#---- 3. Normalize to log2 CPM and z-score ----
dge_tmp <- DGEList(counts = expr_subset)
dge_tmp <- calcNormFactors(dge_tmp)
logCPM <- cpm(dge_tmp, log = TRUE, prior.count = 1)

# Z-score across samples (per gene)
expr_z <- t(scale(t(logCPM)))

#---- 4. Sample annotation ----
ann_col <- data.frame(
  Group = ifelse(colnames(expr_z) %in% young_ids, "Young", "Old")
)
rownames(ann_col) <- colnames(expr_z)

#---- 5. Heatmap ----
pheatmap(expr_z,
         annotation_col = ann_col,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         main = "Inflammaging Gene Set (Robinson Human Dataset)",
         color = colorRampPalette(c("blue", "white", "red"))(50))
#--------------------------------------------------------------------#

install.packages("mcclust")

library(mcclust)

# ---- 1. Hierarchical clustering based on heatmap settings ----
dist_mat <- dist(t(expr_z))                # distance on samples
hc <- hclust(dist_mat, method = "complete")

# Cut into 2 clusters (like heatmap does visually)
cluster_assign <- cutree(hc, k = 2)

# Ensure Young = 1, Old = 2
true_labels <- as.numeric(factor(ann_col$Group, levels = c("Young", "Old")))

# Cluster assignment (k=2)
dist_mat <- dist(t(expr_z))
hc <- hclust(dist_mat, method = "complete")
cluster_assign <- cutree(hc, k = 2)

# Align cluster label orientation
if (mean(true_labels[cluster_assign == 1]) > mean(true_labels[cluster_assign == 2])) {
  cluster_assign <- ifelse(cluster_assign == 1, 2, 1)
}

# Compute ARI
library(mclust)
ari_score <- adjustedRandIndex(cluster_assign, true_labels)
ari_score



#---- 1. Use your existing z-scored data (expr_z) ----
# Transpose: PCA works on samples as rows, genes as columns
expr_z_t <- t(expr_z)

#---- 2. PCA ----
pca_res <- prcomp(expr_z_t, scale. = FALSE)

#---- 3. Prepare dataframe for plotting ----
pca_df <- data.frame(pca_res$x,
                     Sample = rownames(pca_res$x),
                     Group = ann_col$Group[match(rownames(pca_res$x), rownames(ann_col))])

#---- 4. Plot PCA ----
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "PCA of Inflammaging Gene Set (Robinson Human Dataset)",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(text = element_text(size = 14))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 2, size = 1) +  # 95% confidence ellipse
  labs(title = "PCA of Inflammaging Gene Set (Robinson Human Dataset)",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(text = element_text(size = 14))

library(ggforce)  # for geom_mark_ellipse
ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4) +
  geom_mark_ellipse(aes(fill = Group), alpha = 0.2, show.legend = FALSE) +
  labs(title = "PCA of Inflammaging Gene Set (Robinson Human Dataset)") +
  theme_minimal()

table(ann_col$Group)
table(true_labels)
table(cluster_assign, true_labels)

head(colnames(expr_z))
head(rownames(ann_col))

inflamm_score <- colMeans(expr_z)
boxplot(inflamm_score ~ ann_col$Group)
t.test(inflamm_score ~ ann_col$Group)

pca <- prcomp(t(expr_z), scale. = FALSE)
pc1_scores <- pca$x[,1]
t.test(pc1_scores ~ ann_col$Group)

library(pROC)
auc_val <- auc(roc(ann_col$Group, pc1_scores))
auc_val
#-----------------------------------------------#

#------Random testing to check if Inflammaging outperforms random chance-----#

library(edgeR)
library(pROC)

set.seed(123)

# ---- 1. Normalize entire dataset ----
dge_all <- DGEList(counts = count_matrix)  # raw count matrix
dge_all <- calcNormFactors(dge_all)
logCPM_all <- cpm(dge_all, log = TRUE, prior.count = 1)

# ---- 2. Null distribution (random sets of 53 genes) ----
num_genes <- 53
random_abs_auc <- replicate(1000, {
  rand_idx <- sample(nrow(logCPM_all), num_genes)
  expr_rand <- logCPM_all[rand_idx, , drop = FALSE]
  expr_rand_z <- t(scale(t(expr_rand)))
  pca_rand <- prcomp(t(expr_rand_z), scale. = FALSE)
  pc1_rand <- pca_rand$x[, 1]
  roc_rand <- roc(ann_col$Group, pc1_rand, levels = c("Old","Young"))
  auc_rand <- as.numeric(auc(roc_rand))
  max(auc_rand, 1 - auc_rand)  # absolute AUC
})

# ---- 3. Observed value from inflammaging gene set ----
observed_auc <- 0.8534  # previously computed

# ---- 4. Empirical p-value ----
p_val <- mean(random_abs_auc >= observed_auc)

# ---- 5. Plot ----
hist(random_abs_auc, breaks = 30, col = "grey",
     main = "Random 53-gene AUC distribution",
     xlab = "Absolute AUC")
abline(v = observed_auc, col = "red", lwd = 2)
legend("topright", legend = paste("Observed AUC =", round(observed_auc,3),
                                  "\nEmpirical p =", round(p_val,4)),
       bty = "n")

p_val

# ---- 1. Null distribution summary ----
mean_auc <- mean(random_abs_auc)
sd_auc   <- sd(random_abs_auc)

# ---- 2. Z-score of observed result ----
z_score <- (observed_auc - mean_auc) / sd_auc

# ---- 3. Percentile of observed result ----
percentile <- ecdf(random_abs_auc)(observed_auc) * 100

# ---- 4. Output ----
cat("Null mean =", round(mean_auc, 3), 
    "Null sd =", round(sd_auc, 3), "\n")
cat("Observed AUC =", round(observed_auc, 3), "\n")
cat("Z-score =", round(z_score, 2), "\n")
cat("Percentile =", round(percentile, 1), "%\n")
#----------------------------------------------------#

#----------------------------------------------------#
#--------HIIT effect on inflamm gene set-------------#
#----------------------------------------------------#

# Get existing column names (excluding first two annotation columns)
sample_cols <- colnames(human_aging_counts_robinson_merged)[-(1:2)]

# Extract simplified sample IDs (like "47A")
sample_ids <- sub("^s_([^-]+)-.*", "\\1", sample_cols)

# Apply new names
colnames(human_aging_counts_robinson_merged)[-(1:2)] <- sample_ids


old_pre <- c("22A","42A","45A","47A","51A","52A","53A")
old_post <- c("22B","42B","45B","47B","51B","52B","53B")
keep_samples <- c(old_pre, old_post)

# Keep only the first occurrence of each GeneID
count_matrix <- count_matrix[!duplicated(count_matrix$GeneID), ]

# Now set row names
rownames(count_matrix) <- count_matrix$GeneID

# Remove annotation columns for EdgeR
counts_only <- as.matrix(count_matrix[ , -(1:2)])
counts_only[is.na(counts_only)] <- 0


# Build metadata
sample_info <- data.frame(
  SampleID   = keep_samples,
  Individual = gsub("[AB]$", "", keep_samples),  # remove A/B to get individual ID
  Time       = ifelse(grepl("A$", keep_samples), "Pre", "Post")
)
sample_info$Individual <- factor(sample_info$Individual)
sample_info$Time <- factor(sample_info$Time, levels = c("Pre","Post"))

# Build DGEList
dge <- DGEList(counts = counts_only)
dge <- calcNormFactors(dge)

# Paired design: one column for each individual + effect of time
design <- model.matrix(~Individual + Time, data = sample_info)
colnames(design)
# (Intercept), Individual42, Individual45, ..., TimePost

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit model
fit <- glmFit(dge, design)

# Contrast: Post vs Pre (coefficient for TimePost)
lrt <- glmLRT(fit, coef = "TimePost")

# Top differentially expressed genes
topTags(lrt)

# Extract all results as dataframe
deg_results <- topTags(lrt, n = Inf)$table
deg_results <- tibble::rownames_to_column(deg_results, var = "ensembl_gene_id")

head(deg_results)

# Connect to Ensembl archive
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org")

# Map human Ensembl gene IDs to mouse gene symbols
orthologs <- getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = deg_results$ensembl_gene_id,
                    mart = human,
                    attributesL = c("ensembl_gene_id","mgi_symbol"),
                    martL = mouse,
                    uniqueRows = TRUE)

colnames(orthologs) <- c("ensembl_gene_id", "hgnc_symbol",
                         "ensembl_gene_id_mouse", "mgi_symbol")

# Merge with DEG results
deg_mouse_robinson_HIIT <- merge(deg_results, orthologs, by = "ensembl_gene_id")

head(deg_mouse_robinson_HIIT)
dim(deg_mouse_robinson_HIIT)

deg_mouse_robinson_HIIT %>%
  filter(FDR <0.05)

validated_inflamm_geneset.vec

deg_mouse_robinson_HIIT %>%
  filter(mgi_symbol %in% validated_inflamm_geneset.vec) %>%
  arrange(logFC)

#-------GSEA-------------------------------#
#------------------------------------------#
# ---- 1. Create ranked gene list ----
geneList <- deg_mouse_robinson_HIIT %>%
  mutate(stat = -log10(PValue) * sign(logFC)) %>%
  arrange(desc(stat)) %>%
  { setNames(.$stat, .$mgi_symbol) }

# ---- 2. Run GSEA with your inflammaging gene set ----
gsea_inflamm_human_robinson_HIIT <- GSEA(geneList = geneList,
                                    TERM2GENE = data.frame(term = "Inflamm", 
                                                           gene = validated_inflamm_geneset.vec),
                                    pvalueCutoff = 1,
                                    verbose = FALSE)

# ---- 3. View results ----
gsea_inflamm_human_robinson_HIIT@result %>%
  dplyr::select(ID, NES, pvalue, p.adjust)
#------------------------------------------------#
#------------------------------------------------#


#----------------------------------------------------#
#--------Combined ex effect on inflamm gene set-------------#
#----------------------------------------------------#
#----------------------------------------------------#
# 1. Start from the raw merged count matrix
#----------------------------------------------------#
# Keep original counts intact
raw_counts <- human_aging_counts_robinson_merged  
head(human_aging_counts_robinson_merged)
human_aging_counts_robinson_merged %>%
  arrange(desc("1B"))

#----------------------------------------------------#
# 2. Simplify sample column names (keep core ID like "1B")
#----------------------------------------------------#
sample_cols <- colnames(raw_counts)[-(1:2)]  # Exclude GeneID and GeneName
sample_ids <- sub("^s_([^-]+)-.*", "\\1", sample_cols)
colnames(raw_counts)[-(1:2)] <- sample_ids

#----------------------------------------------------#
# 3. Define the paired samples
#----------------------------------------------------#
old_combined_pre  <- c("1B","20B","37B","46B","48B","5B","7B")
old_combined_post <- c("1C","20C","37C","46C","48C","5C","7C")
keep_samples <- c(old_combined_pre, old_combined_post)

#----------------------------------------------------#
# 4. Subset to these samples, remove duplicates, set rownames
#----------------------------------------------------#
count_matrix <- raw_counts %>%
  dplyr::select(GeneID, GeneName, all_of(keep_samples)) %>%
  dplyr::filter(!duplicated(GeneID))

rownames(count_matrix) <- count_matrix$GeneID
counts_only <- as.matrix(count_matrix[, -(1:2)])
counts_only[is.na(counts_only)] <- 0  # replace NAs with 0

#----------------------------------------------------#
# 5. Build sample metadata
#----------------------------------------------------#
sample_info <- data.frame(
  SampleID   = keep_samples,
  Individual = gsub("[BC]$", "", keep_samples),     # "1B" → "1"
  Time       = ifelse(grepl("B$", keep_samples), "Pre", "Post")
)
sample_info$Individual <- factor(sample_info$Individual)
sample_info$Time       <- factor(sample_info$Time, levels = c("Pre","Post"))

# Ensure column order of counts matches sample_info
counts_only <- counts_only[, sample_info$SampleID]

#----------------------------------------------------#
# 6. EdgeR pipeline
#----------------------------------------------------#
library(edgeR)

dge <- DGEList(counts = counts_only)
dge <- calcNormFactors(dge)

# Model matrix with individual blocking + time effect
design <- model.matrix(~Individual + Time, data = sample_info)

dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = "TimePost")

#----------------------------------------------------#
# 7. DEG table
#----------------------------------------------------#
deg_results <- topTags(lrt, n = Inf)$table
deg_results <- tibble::rownames_to_column(deg_results, var = "ensembl_gene_id")

#----------------------------------------------------#
# 8. Optional: Map human → mouse orthologs
#----------------------------------------------------#
library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org")

orthologs <- getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = deg_results$ensembl_gene_id,
                    mart = human,
                    attributesL = c("ensembl_gene_id","mgi_symbol"),
                    martL = mouse,
                    uniqueRows = TRUE)
colnames(orthologs) <- c("ensembl_gene_id", "hgnc_symbol",
                         "ensembl_gene_id_mouse", "mgi_symbol")

deg_mouse_robinson_combined <- merge(deg_results, orthologs, by = "ensembl_gene_id")


head(deg_mouse_robinson_combined)
dim(deg_mouse_robinson_combined)

deg_mouse_robinson_combined %>%
  filter(FDR <0.05) %>%
  filter(abs(logFC) >0.5)

validated_inflamm_geneset.vec

deg_mouse_robinson_combined %>%
  filter(mgi_symbol %in% validated_inflamm_geneset.vec) %>%
  arrange(logFC)

#-------GSEA-------------------------------#
#------------------------------------------#
# ---- 1. Create ranked gene list ----
geneList <- deg_mouse_robinson_HIIT %>%
  mutate(stat = -log10(PValue) * sign(logFC)) %>%
  arrange(desc(stat)) %>%
  { setNames(.$stat, .$mgi_symbol) }

# ---- 2. Run GSEA with your inflammaging gene set ----
gsea_inflamm_human_robinson_HIIT <- GSEA(geneList = geneList,
                                         TERM2GENE = data.frame(term = "Inflamm", 
                                                                gene = validated_inflamm_geneset.vec),
                                         pvalueCutoff = 1,
                                         verbose = FALSE)

# ---- 3. View results ----
gsea_inflamm_human_robinson_HIIT@result %>%
  dplyr::select(ID, NES, pvalue, p.adjust)

deg_mouse_robinson_HIIT
#------------------------------------------------#
#------------------------------------------------#

# Confirm column order alignment
stopifnot(all(colnames(counts_only) == sample_info$SampleID))

# Check factor levels
levels(sample_info$Individual)
table(sample_info$Individual, sample_info$Time)

# Check design matrix
design
qr(design)$rank    # should be (#samples - #unique Individuals) + 1

# Check dispersion
dge <- estimateDisp(dge, design)
summary(dge$common.dispersion)

ENSG00000169245

top_gene <- rownames(topTags(lrt, n = 1)$table)
plot(cpm(dge)[top_gene,], col=sample_info$Time, main=top_gene)

any(counts_only %% 1 != 0)
summary(abs(topTags(lrt, n = 100)$table$logFC))

all(colnames(counts_only) == sample_info$SampleID)
