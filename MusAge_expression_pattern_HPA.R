validated_inflamm_geneset_upper.vec

sc_data <- read.delim("rna_single_cell_type.tsv")

sc_data

# Load necessary library
library(biomaRt)
library(dplyr)

# Connect to Ensembl archive (Dec 2021)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org")

mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org")

# Map Human Gene Symbols to Mouse Gene Symbols
orthologs <- getLDS(attributes = c("hgnc_symbol"),
                    filters = "hgnc_symbol",
                    values = unique(sc_data$Gene.name),
                    mart = human,
                    attributesL = c("mgi_symbol"),
                    martL = mouse,
                    uniqueRows = TRUE)

colnames(orthologs) <- c("Gene.name", "Mouse.gene")

# Merge with your original dataset to get Mouse gene names
sc_data_mouse <- sc_data %>%
  left_join(orthologs, by = "Gene.name")

sc_data_mouse
validated_inflamm_geneset.vec

#Pull out or MusAge overlapping genes from total Human Protein Atlas data
sc_MusAge <- sc_data_mouse[sc_data_mouse$Mouse.gene %in% validated_inflamm_geneset.vec, ]

sc_MusAge %>%
  filter(Mouse.gene == "Sox4")

#Create Heatmap for expression levels in cell types for each gene

# Convert to wide format: rows = genes, columns = cell types
sc_wide <- sc_MusAge %>%
  group_by(Mouse.gene, Cell.type) %>%
  summarise(nTPM = mean(nTPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cell.type, values_from = nTPM)

# Set Gene.name as rownames
sc_matrix <- as.data.frame(sc_wide)
rownames(sc_matrix) <- sc_matrix$Mouse.gene
sc_matrix$Mouse.gene <- NULL

# Turn into numeric matrix
sc_matrix <- as.matrix(sc_matrix)

# log transform
sc_matrix_log <- log2(sc_matrix + 1)

# Plot heatmap
library(pheatmap)
library(RColorBrewer)

# Define custom color palette: blue-white-red
breaks <- seq(-2.5, 2.5, length.out = 101)
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

pheatmap(
  sc_matrix_log,
  scale = "row",
  show_rownames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = my_colors,
  breaks = breaks,
  main = "Z-scored Expression of Genes Across Cell Types"
)

head(sc_matrix_log)
#--------------------------------------------------#

# Cell type -> Category vector
celltype_categories <- c(
  # Neural
  "Excitatory neurons" = "Neural",
  "Inhibitory neurons" = "Neural",
  "Astrocytes" = "Neural",
  "Microglial cells" = "Neural",
  "Oligodendrocytes" = "Neural",
  "Oligodendrocyte precursor cells" = "Neural",
  "Schwann cells" = "Neural",
  "Bipolar cells" = "Neural",
  "Horizontal cells" = "Neural",
  "Muller glia cells" = "Neural",
  "Rod photoreceptor cells" = "Neural",
  "Cone photoreceptor cells" = "Neural",
  
  # Muscle
  "Skeletal myocytes" = "Muscle",
  "Cardiomyocytes" = "Muscle",
  "Smooth muscle cells" = "Muscle",
  "Peritubular cells" = "Muscle",
  "Breast myoepithelial cells" = "Muscle",
  
  # Immune
  "T-cells" = "Immune",
  "B-cells" = "Immune",
  "NK-cells" = "Immune",
  "Dendritic cells" = "Immune",
  "Macrophages" = "Immune",
  "Monocytes" = "Immune",
  "Granulocytes" = "Immune",
  "Plasma cells" = "Immune",
  "Langerhans cells" = "Immune",
  "Hofbauer cells" = "Immune",
  
  # Epithelial/Secretory
  "Basal keratinocytes" = "Epithelial",
  "Suprabasal keratinocytes" = "Epithelial",
  "Squamous epithelial cells" = "Epithelial",
  "Basal squamous epithelial cells" = "Epithelial",
  "Glandular and luminal cells" = "Epithelial",
  "Serous glandular cells" = "Epithelial",
  "Mucus glandular cells" = "Epithelial",
  "Exocrine glandular cells" = "Epithelial",
  "Ductal cells" = "Epithelial",
  "Salivary duct cells" = "Epithelial",
  "Secretory cells" = "Epithelial",
  "Ciliated cells" = "Epithelial",
  "Club cells" = "Epithelial",
  "Ionocytes" = "Epithelial",
  "Breast glandular cells" = "Epithelial",
  "Prostatic glandular cells" = "Epithelial",
  "Basal prostatic cells" = "Epithelial",
  
  # Endothelial/Stromal
  "Endothelial cells" = "Endothelial/Stromal",
  "Lymphatic endothelial cells" = "Endothelial/Stromal",
  "Mesothelial cells" = "Endothelial/Stromal",
  "Fibroblasts" = "Endothelial/Stromal",
  "Endometrial stromal cells" = "Endothelial/Stromal",
  "Adipocytes" = "Endothelial/Stromal",
  
  # Digestive
  "Pancreatic endocrine cells" = "Digestive",
  "Enteroendocrine cells" = "Digestive",
  "Gastric mucus-secreting cells" = "Digestive",
  "Intestinal goblet cells" = "Digestive",
  "Distal enterocytes" = "Digestive",
  "Proximal enterocytes" = "Digestive",
  "Paneth cells" = "Digestive",
  "Hepatocytes" = "Digestive",
  
  # Reproductive/Germline
  "Oocytes" = "Reproductive",
  "Granulosa cells" = "Reproductive",
  "Ovarian stromal cells" = "Reproductive",
  "Sertoli cells" = "Reproductive",
  "Leydig cells" = "Reproductive",
  "Spermatogonia" = "Reproductive",
  "Spermatocytes" = "Reproductive",
  "Early spermatids" = "Reproductive",
  "Late spermatids" = "Reproductive",
  
  # Other / Specialized
  "Kupffer cells" = "Other",
  "Cholangiocytes" = "Other",
  "Collecting duct cells" = "Other",
  "Distal tubular cells" = "Other",
  "Proximal tubular cells" = "Other",
  "Cytotrophoblasts" = "Other",
  "Extravillous trophoblasts" = "Other",
  "Syncytiotrophoblasts" = "Other",
  "Undifferentiated cells" = "Other"
)

# Ensure column names are in the category vector
ordered_cols <- names(sort(celltype_categories[colnames(sc_matrix_log)]))

# Reorder the matrix
sc_matrix_log_ordered <- sc_matrix_log[, ordered_cols]

# Create annotation for heatmap
annotation_col <- data.frame(Category = celltype_categories[ordered_cols])
rownames(annotation_col) <- ordered_cols


# Assign color palette
category_colors <- RColorBrewer::brewer.pal(length(unique(annotation_col$Category)), "Set2")
names(category_colors) <- unique(annotation_col$Category)

ann_colors <- list(Category = category_colors)

# 5. Color palette and breaks
max_abs <- 2.25
cols <- colorRampPalette(c("royalblue4", "white", "firebrick3"))(50)
brks <- seq(-max_abs, max_abs, length.out = length(cols) + 1)


pheatmap(sc_matrix_log_ordered,
         scale = "row",
         color = cols,
         breaks = brks,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,  # Turn off column clustering
         cluster_rows = TRUE,   # Optional
         show_rownames = TRUE,
         fontsize_col = 7,
         border_color = NA,
         main = "Cell-Type Expression Patterns Grouped by Category")
#----------------------------------------------------------------------#

# --- 1) Row-wise z-scores (same as pheatmap(scale="row")) ---
z_mat <- t(scale(t(sc_matrix_log_ordered)))  # rows = genes, cols = cell types

# Optional: drop rows that became all-NA (constant rows)
keep_rows <- apply(z_mat, 1, function(x) all(is.finite(x)))
z_mat <- z_mat[keep_rows, , drop = FALSE]

#column redistribution
# Start from your current mapping (with your Digestive -> Endocrine/Mucosal rename if you applied it)
cat_map <- if (exists("celltype_categories_refined")) {
  celltype_categories_refined
} else {
  celltype_categories
}

# --- Option A (recommended): reassign to existing categories, no new ones ---
relabel_A <- c(
  "Kupffer cells"            = "Immune",
  "Cholangiocytes"           = "Epithelial",  # or "Epithelial" if you prefer
  "Collecting duct cells"    = "Epithelial",
  "Distal tubular cells"     = "Epithelial",
  "Proximal tubular cells"   = "Epithelial",
  "Cytotrophoblasts"         = "Reproductive",
  "Extravillous trophoblasts"= "Reproductive",
  "Syncytiotrophoblasts"     = "Reproductive",
  "Undifferentiated cells"   = "Epithelial"          # avoids a 1-item category
)

celltype_categories_refined <- cat_map
celltype_categories_refined[names(relabel_A)] <- unname(relabel_A)

# --- 2) Long format + category and ordering ---
# --- rebuild Category in z_long using the refined map ---
z_long <- as.data.frame(z_mat) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Cell.type", values_to = "z") %>%
  dplyr::mutate(
    Category = celltype_categories_refined[Cell.type],
    Category = ifelse(is.na(Category), "Epithelial", Category)  # fallback, rarely used
  )

# Category-ordered columns, then your within-category clustering block can run as-is
ordered_cols <- names(sort(celltype_categories_refined[colnames(sc_matrix_log_ordered)]))
z_long$Cell.type <- factor(z_long$Cell.type, levels = ordered_cols)

# (If you’re using the within-category clustering code, keep it but reference celltype_categories_refined)
cat_seq <- unique(celltype_categories_refined[ordered_cols])

cluster_within <- function(cols) {
  cols <- intersect(cols, colnames(z_mat))
  if (length(cols) <= 1) return(cols)
  d  <- as.dist(1 - cor(z_mat[, cols, drop = FALSE],
                        method = "spearman", use = "pairwise.complete.obs"))
  h  <- hclust(d, method = "average")
  cols[h$order]
}

col_order_blocked <- unlist(lapply(cat_seq, function(cat) {
  cols_in_cat <- ordered_cols[celltype_categories_refined[ordered_cols] == cat]
  cluster_within(cols_in_cat)
}), use.names = FALSE)

z_long$Cell.type <- factor(z_long$Cell.type, levels = col_order_blocked)
z_long$Category  <- factor(z_long$Category, levels = unique(celltype_categories_refined[col_order_blocked]))

# --- 3) Threshold + rescale within gene, then optional contrast ---
# knobs to tune:
z_cut <- 0          # absolute z cutoff (e.g., 0 = show only above row mean)
q_cut <- 0.25       # OR quantile cutoff within gene (e.g., keep top 75%)
gamma <- 1.5        # >1 emphasizes high values; 1 = linear; <1 flattens

z_long_thr <- z_long %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(
    # choose a per-gene threshold: stricter of absolute z_cut and quantile q_cut
    thr_z  = max(z_cut, quantile(z, q_cut, na.rm = TRUE)),
    # keep only signal above threshold for sizing
    z_keep = pmax(z - thr_z, 0),
    span   = max(z - thr_z, na.rm = TRUE),
    span   = ifelse(is.finite(span) & span > 0, span, 1e-9),
    # rescale to 0..1 and add optional contrast
    z_rel  = (z_keep / span) ^ gamma
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(z > thr_z)   # hide dots below threshold entirely

# --- 4) Colors for categories (8 categories -> Set2 works nicely) ---
category_colors <- RColorBrewer::brewer.pal(length(unique(z_long$Category)), "Set2")
names(category_colors) <- unique(z_long$Category)

# --- 5) Dotplot ---
p_dot <- ggplot(z_long_thr, aes(x = Cell.type, y = Gene)) +
  geom_point(aes(size = z_rel, fill = Category),
             shape = 21, color = "grey25", stroke = 0.4) +
  scale_size_area(
    name = "Row-scaled (≥ threshold)",
    max_size = 6, limits = c(0, 1),
    breaks = c(0.25, 0.5, 0.75, 1),
    labels = c("low+", "mid", "high", "max")
  ) +
  scale_fill_manual(values = category_colors, name = "Cell type category") +
  scale_x_discrete(drop = FALSE) +  # keep empty columns if any
  scale_y_discrete(drop = FALSE) +  # keep empty rows if any
  labs(x = NULL, y = NULL,
       title = "Dotplot: Row-z expression with thresholded sizing") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(size = 0.2))

p_dot

# Save as PDF
pdf("p_dot.pdf", width = 8, height = 8)
print(p_dot)
dev.off()

#==============================================================#
#-----Top Cell Type Expression Table--------------------------#
#==============================================================#
library(tidyr)

# Convert the matrix to a tidy long format
sc_long <- sc_matrix_log %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Expression")

# For each gene, find the top 3 expressing cell types
top3_by_gene <- sc_long %>%
  group_by(Gene) %>%
  arrange(desc(Expression), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  mutate(Rank = row_number()) %>%
  pivot_wider(names_from = Rank, values_from = c(CellType, Expression)) %>%
  ungroup()

# View the top 3 cell types for each gene
print(top3_by_gene, n=52)

top3_by_gene %>%
  arrange(CellType_1)
#-----------------------------------------#

# Ensure your sc_matrix_log is a data frame with gene names as rownames
sc_df <- as.data.frame(sc_matrix_log)
sc_df


# Cell type -> Category vector
celltype_categories <- c(
  # Neural
  "Excitatory neurons" = "Neural",
  "Inhibitory neurons" = "Neural",
  "Astrocytes" = "Neural",
  "Microglial cells" = "Neural",
  "Oligodendrocytes" = "Neural",
  "Oligodendrocyte precursor cells" = "Neural",
  "Schwann cells" = "Neural",
  "Bipolar cells" = "Neural",
  "Horizontal cells" = "Neural",
  "Muller glia cells" = "Neural",
  "Rod photoreceptor cells" = "Neural",
  "Cone photoreceptor cells" = "Neural",
  
  # Muscle
  "Skeletal myocytes" = "Muscle",
  "Cardiomyocytes" = "Muscle",
  "Smooth muscle cells" = "Muscle",
  "Peritubular cells" = "Muscle",
  "Breast myoepithelial cells" = "Muscle",
  
  # Immune
  "T-cells" = "Immune",
  "B-cells" = "Immune",
  "NK-cells" = "Immune",
  "Dendritic cells" = "Immune",
  "Macrophages" = "Immune",
  "Monocytes" = "Immune",
  "Granulocytes" = "Immune",
  "Plasma cells" = "Immune",
  "Langerhans cells" = "Immune",
  "Hofbauer cells" = "Immune",
  
  # Epithelial/Secretory
  "Basal keratinocytes" = "Epithelial",
  "Suprabasal keratinocytes" = "Epithelial",
  "Squamous epithelial cells" = "Epithelial",
  "Basal squamous epithelial cells" = "Epithelial",
  "Glandular and luminal cells" = "Epithelial",
  "Serous glandular cells" = "Epithelial",
  "Mucus glandular cells" = "Epithelial",
  "Exocrine glandular cells" = "Epithelial",
  "Ductal cells" = "Epithelial",
  "Salivary duct cells" = "Epithelial",
  "Secretory cells" = "Epithelial",
  "Ciliated cells" = "Epithelial",
  "Club cells" = "Epithelial",
  "Ionocytes" = "Epithelial",
  "Breast glandular cells" = "Epithelial",
  "Prostatic glandular cells" = "Epithelial",
  "Basal prostatic cells" = "Epithelial",
  
  # Endothelial/Stromal
  "Endothelial cells" = "Endothelial/Stromal",
  "Lymphatic endothelial cells" = "Endothelial/Stromal",
  "Mesothelial cells" = "Endothelial/Stromal",
  "Fibroblasts" = "Endothelial/Stromal",
  "Endometrial stromal cells" = "Endothelial/Stromal",
  "Adipocytes" = "Endothelial/Stromal",
  
  # Digestive
  "Pancreatic endocrine cells" = "Digestive",
  "Enteroendocrine cells" = "Digestive",
  "Gastric mucus-secreting cells" = "Digestive",
  "Intestinal goblet cells" = "Digestive",
  "Distal enterocytes" = "Digestive",
  "Proximal enterocytes" = "Digestive",
  "Paneth cells" = "Digestive",
  "Hepatocytes" = "Digestive",
  
  # Reproductive/Germline
  "Oocytes" = "Reproductive",
  "Granulosa cells" = "Reproductive",
  "Ovarian stromal cells" = "Reproductive",
  "Sertoli cells" = "Reproductive",
  "Leydig cells" = "Reproductive",
  "Spermatogonia" = "Reproductive",
  "Spermatocytes" = "Reproductive",
  "Early spermatids" = "Reproductive",
  "Late spermatids" = "Reproductive",
  
  # Other / Specialized
  "Kupffer cells" = "Other",
  "Cholangiocytes" = "Other",
  "Collecting duct cells" = "Other",
  "Distal tubular cells" = "Other",
  "Proximal tubular cells" = "Other",
  "Cytotrophoblasts" = "Other",
  "Extravillous trophoblasts" = "Other",
  "Syncytiotrophoblasts" = "Other",
  "Undifferentiated cells" = "Other"
)
# Keep only cell types that exist in the matrix
valid_celltypes <- intersect(names(celltype_categories), colnames(sc_df))

# Step 0: Start fresh with gene expression matrix
sc_df_sub <- sc_df %>%
  dplyr::select(all_of(valid_celltypes)) %>%
  rownames_to_column("Gene")

# Step 1: Convert to long format and map to major categories
sc_long <- sc_df_sub %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Expression") %>%
  mutate(Category = celltype_categories[CellType])

# Step 2: Count how many cell types per category (used for normalization)
celltype_counts <- tibble(
  CellType = names(celltype_categories),
  Category = celltype_categories
) %>%
  distinct() %>%
  count(Category, name = "N_CellTypes")

# Step 3: Sum expression for each Gene × Category
category_scores <- sc_long %>%
  group_by(Gene, Category) %>%
  summarise(Sum = sum(Expression, na.rm = TRUE), .groups = "drop")

# Step 4: Join category sizes and normalize
category_scores <- category_scores %>%
  left_join(celltype_counts, by = "Category") %>%
  mutate(Normalized = Sum / N_CellTypes)

# Step 5: Reshape wide for heatmap annotation
category_wide <- category_scores %>%
  dplyr::select(Gene, Category, Normalized) %>%
  pivot_wider(names_from = Category, values_from = Normalized, values_fill = 0)

# Step 6: Determine top category per gene
gene_category_df <- category_wide %>%
  rowwise() %>%
  mutate(
    TopCategory = names(pick(where(is.numeric)))[which.max(c_across(where(is.numeric)))]
  ) %>%
  ungroup()

# Check result
head(gene_category_df)

# View result
print(gene_category_df, n=52)
print(gene_category_df %>%
  arrange(TopCategory), n=52)

sc_matrix_log
