#Identify Expression pattern of Inflammaging Gene Set

getwd()
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")
getwd()
sc_data <- read.delim("rna_single_cell_type.tsv")
Inflammaging_symbols
Inflammaging_symbols_upper <- toupper(Inflammaging_symbols)
validated_inflamm_geneset_upper.vec

#Pull out or bcl6b overlapping genes from total Human Protein Atlas data
sc_subset <- sc_data[sc_data$Gene.name %in% validated_inflamm_geneset_upper.vec, ]

# View result
head(sc_subset)
dim(sc_subset)
#--------------------------#
#Create Heatmap for expression levels in cell types for each gene

# Convert to wide format: rows = genes, columns = cell types
sc_wide <- sc_subset %>%
  dplyr::select(Gene.name, Cell.type, nTPM) %>%
  pivot_wider(names_from = Cell.type, values_from = nTPM)

# Set Gene.name as rownames
sc_matrix <- as.data.frame(sc_wide)
rownames(sc_matrix) <- sc_matrix$Gene.name
sc_matrix$Gene.name <- NULL

# Turn into numeric matrix
sc_matrix <- as.matrix(sc_matrix)

# log transform
sc_matrix_log <- log2(sc_matrix + 1)

# Plot heatmap
library(pheatmap)
pheatmap(sc_matrix_log, scale = "row", show_rownames = TRUE)

#--------------------------------------------------------------#
#------------Determine nTPM Cutoff values for cell type enrichment analysis------#

expression_summary <- sc_subset %>%
  group_by(Cell.type) %>%
  summarise(
    median_nTPM = median(nTPM, na.rm = TRUE),
    Q1_nTPM = quantile(nTPM, 0.25, na.rm = TRUE),
    Q3_nTPM = quantile(nTPM, 0.75, na.rm = TRUE),
    mean_nTPM = mean(nTPM, na.rm = TRUE),
    max_nTPM = max(nTPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_nTPM))

# View top rows
head(expression_summary)

library(ggplot2)

ggplot(sc_subset, aes(x = reorder(Cell.type, nTPM, FUN = median), y = nTPM)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  labs(x = "Cell Type", y = "nTPM", title = "Distribution of Bcl6b Target Gene Expression Across Cell Types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(ylim = c(0, 1000))

expression_summary
#---------------------------------#

#Use Fisher Exact Test to determine if there is significant cell type
#enrichment in Inflammaging genes#

#Get all genes identified in full data set as background genes
#all_ex_BioCor_data_entrez
oldsed_vs_yngsed_genes

# Extract background gene symbols and convert to uppercase
#background_genes <- toupper(all_ex_BioCor_data_entrez$Symbol.x)
background_genes <- toupper(oldsed_vs_yngsed_genes$Symbol)
head(background_genes)

#Run Fisher Exact Test to determine enrichment
#Set nTPM threshold as 10 for each cell type
#and mark 1 if meets threshold or 0 if does not
fisher_results <- sc_data %>%
  filter(Gene.name %in% c(Inflammaging_symbols_upper, background_genes)) %>%
  mutate(
    GeneSet = ifelse(Gene.name %in% Inflammaging_symbols_upper, "Inflamm", "Background"),
    Expressed = ifelse(nTPM >= 10, 1, 0)
  ) %>%
  group_by(Cell.type) %>%
  summarise(
    A = sum(GeneSet == "Inflamm" & Expressed == 1),
    B = sum(GeneSet == "Background" & Expressed == 1),
    C = sum(GeneSet == "Inflamm" & Expressed == 0),
    D = sum(GeneSet == "Background" & Expressed == 0)
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(A, B, C, D), nrow = 2))$p.value,
    odds_ratio = fisher.test(matrix(c(A, B, C, D), nrow = 2))$estimate
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # FDR correction

fisher_results %>% 
  arrange(desc(odds_ratio))
