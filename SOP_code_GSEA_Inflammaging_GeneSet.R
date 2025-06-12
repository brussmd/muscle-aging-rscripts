#This is the code for creating GSEA plot, identifying core enrichment genes,
#establishing the Inflammaging geneset, then visualizing the heatmaps for
#total core enrichment genes and a separate heatmap for the Inflammaging gene set.
#-------------------------------------------#

#This code is based off of pre-processed RNA Seq data.
#These files have the "..._GENE_..." tag in them.
#Example: 20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx
#----------------------------------------------#

#In this code we are investigtaing the differences between
#Old_Sedentary_Vehicle mice vs Young_Sedentary_Vehicle mice.
#The code does three things.
#1. Create Gene Ontology comaparison with GSEA.
#2. Identifies the Core Enrichment genes from all of the GSEA pathways that
#   are upregulated in old_sed_veh mice vs yng_sed_veh mice.
#3. Finally this code identifies the genes from these Core Enrichment genes
#   that are significantly different between old_sed_veh and yng_sed_veh (FDR <0.05)
#---------------------------------------------------#

#---------SOP--------------------------------#

#Step 1: Upload xlsx GENE files into a folder you will use for coding.

#Step 2: set your working directory in R, then confirm it is set to where
#         your files are stored.
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")
getwd()

#Step 3: Check for required packages that need for this code.
#         Install any necessary packages that you are missing.

#Check packages that you already might have installed
# Define the list of required packages (CRAN and Bioconductor)
required_pkgs <- c(
  "readxl", "dplyr", "tidyverse", "readr", "stringr", "ggplot2",
  "ggrepel", "pheatmap", "VennDiagram", "ggVennDiagram", "eulerr",
  "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "DOSE", "fgsea",
  "ComplexHeatmap"
)

# Function to check installation status
check_installed_pkgs <- function(pkgs) {
  status <- sapply(pkgs, function(pkg) requireNamespace(pkg, quietly = TRUE))
  df <- data.frame(
    Package = names(status),
    Installed = status
  )
  df
}

# Run the check
pkg_status <- check_installed_pkgs(required_pkgs)
print(pkg_status)

# Intall the necessary CRAN packages if you don't already have them installed.
#install.packages(c("readxl", "dplyr", "tidyverse", "readr", "stringr",
#                   "ggplot2", "ggrepel", "pheatmap", "VennDiagram", 
#                   "ggVennDiagram", "eulerr"))

# Bioconductor packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", 
#                       "AnnotationDbi", "DOSE", "fgsea", "ComplexHeatmap"))

#Step 4: Load necessary libraries.
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
library(VennDiagram)
library(pheatmap)
library(ComplexHeatmap)
library(fgsea)
library(ggrepel)
library(stringr)

#Step 5: Load in the Read, Clean and ENTREZID creation function.
# This function 
#1. reads in a processed xlsx file.
#2. Creates columns with relevant names.
#3. Creates an ENTREZID column based off of the ENSEMBL IDs
read_clean_xlsx <- function(xlsx_file, logFC_col, FDR_col) {
  # Read in the xlsx file
  df <- read_xlsx(xlsx_file)
  
  # Clean and rename columns
  df <- df %>%
    dplyr::select(Ensembl, Symbol, 
                  logFC = all_of(logFC_col),
                  FDR = all_of(FDR_col)) %>%
    dplyr::mutate(
      Ensembl_noDec = str_remove(Ensembl, "\\..+"),
      ENTREZID = mapIds(org.Mm.eg.db,
                        keys = Ensembl_noDec,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")
    ) %>%
    dplyr::select(-Ensembl) %>%
    tidyr::drop_na(ENTREZID)%>%
    as.data.frame()
  
  return(df)
}
#------------------------------------#

#Step 6: Load function for creating an ENTREZID vector of the full RNA Seq
#         gene dataframe
#----------Create named vector from dataframe function-----------#
create_vec_from_df <- function(df, df_name) {
  vec <- df %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(ENTREZID, logFC) %>%
    tidyr::drop_na() %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
    tibble::deframe()
  
  assign(paste0(df_name, ".vec"), vec, envir = .GlobalEnv)
}
#-----------------------------------------------------------#

#Step 7: Run both functions to read, clean and create ENTREZID gene vector.
#Function 1.
oldsedveh_v_yngsedveh_genes_12jun <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx",
  logFC_col = "OLD_SED_VEH-YNG_SED_VEH_logFC",
  FDR_col = "OLD_SED_VEH-YNG_SED_VEH_FDR")

#Check to see if function worked properly 
#     (should be 701 DEG for old_sed_veh vs yng_sed_veh)
oldsedveh_v_yngsedveh_genes_12jun %>%
  filter(FDR <0.05)

#Function 2.
create_vec_from_df(oldsedveh_v_yngsedveh_genes_12jun, "oldsedveh_v_yngsedveh_12jun")

#Check to see if created an ordered ENTREZID vector.
#should look like this:
# 108024    56636    16071    66755    16069    14563 
#4.583792 4.375729 4.312056 3.862240 3.612484 3.396991 
head(oldsedveh_v_yngsedveh_12jun.vec)
#----------------------------------------------------#

#Step 8: Run GSEA on Full Gene set vector created in Step 7.
sig_gseGO_oldsed_vs_yngsed_12jun.OUTPUT <- gseGO(geneList = oldsedveh_v_yngsedveh_12jun.vec, 
                                          ont = 'BP',
                                          OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 300,
                                          eps = 1e-30, pvalueCutoff = 0.05)

#check to see if GSEA worked by seeing how many pathways identified.
#should look something like this (but might not be exact)
#[1] 275  11
dim(sig_gseGO_oldsed_vs_yngsed_12jun.OUTPUT)

#Create dataframe of results of GSEA, and add a gene enrichment Count column.
# column called Count
sig_gseGO_oldsed_vs_yngsed_12jun.df <- sig_gseGO_oldsed_vs_yngsed_12jun.OUTPUT@result %>% 
  dplyr::mutate(., Count = str_count(.$core_enrichment, '/')+1)

#check to see if df created properly
head(sig_gseGO_oldsed_vs_yngsed_12jun.df)
#----------------------------------------------#

#Step 9: Run simplified GSEA (this is the parent terms)
oldsed_vs_yngsed_simplified_12jun <- simplify(sig_gseGO_oldsed_vs_yngsed_12jun.OUTPUT, cutoff = 0.5, 
                                        by = "p.adjust", select_fun = min)

#check if simplified ran correctly
dim(oldsed_vs_yngsed_simplified_12jun) # should be: [1] 45 11 (45 significant GO pathways)
head(oldsed_vs_yngsed_simplified_12jun)


#create df of results and add count column
oldsed_vs_yngsed_12jun_simpl.df <- oldsed_vs_yngsed_simplified_12jun@result %>% 
  dplyr::mutate(., Count = str_count(.$core_enrichment, '/')+1)

#Check full simplified GO data set.
oldsed_vs_yngsed_12jun_simpl.df

#//////OPTIONAL///////////////#
#In some of these there is a very long GO Description that makes visualization difficult
#in later steps. So we will shorten it here.
# Modify the description for the specific GO term
oldsed_vs_yngsed_12jun_simpl.df <- oldsed_vs_yngsed_12jun_simpl.df  %>%
  mutate(Description = ifelse(ID == "GO:0043280", "endopeptidase involved in apoptotic process", Description))

# Check if the change was applied correctly
oldsed_vs_yngsed_12jun_simpl.df %>%
  filter(ID == "GO:0043280")
#////////////////////////////////////#
#-------------------------------------------------#

#Step 10: Visualize GO pathways
#   This creates Figure 1 of the working manuscript.
oldsed_vs_yngsed_simplified_12jun_kegg_plot <- oldsed_vs_yngsed_12jun_simpl.df %>%
  ggplot(., aes(x = NES, y=reorder(Description, NES), size = Count, col = p.adjust )) +
  geom_point() +
  theme_bw()  +
  theme(axis.text.y = element_text(size = 7),  # Adjust y-axis label size
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.x = element_text(size = 7))  # Adjust x-axis label size

# Save as PDF
pdf("oldsed_vs_yngsed_simplified_12jun_kegg_plot.pdf", width = 8, height = 6)  # Adjust height
print(oldsed_vs_yngsed_simplified_12jun_kegg_plot)
dev.off()

# Print in RStudio
print(oldsed_vs_yngsed_simplified_12jun_kegg_plot)
#--------------------------------------------------#

#Step 11: Identify all of the core enrichment genes from the upregulated
#         GO pathways (NES >0)

#Get all upregulated GO pathways
upreg_oldsed_simplkegg_12jun <- oldsed_vs_yngsed_12jun_simpl.df %>%
  filter (NES > 0)

#check that pathways match the upregulated pathways from original GSEA figure.
upreg_oldsed_simplkegg_12jun %>%
  arrange(desc(NES))

#------------Get Core Enrichment Genes--------------------------#
# Extract all ENTREZ IDs from the 'core_enrichment' column, split them into a vector
upreg_old_sed_core_entrez_vec_12jun <- upreg_oldsed_simplkegg_12jun$core_enrichment %>%
  str_split(pattern = "/") %>%
  unlist()

# Get unique ENTREZ IDs and sort them
upreg_old_sed_core_entrez_vec_12jun <- unique(upreg_old_sed_core_entrez_vec_12jun) %>%
  as.numeric() %>%  # Convert from character to numeric
  sort()

#Check length and structure of core enrichment gene vector
length(upreg_old_sed_core_entrez_vec_12jun) #Not always the same got 471 this time, 499 last time
head(upreg_old_sed_core_entrez_vec_12jun)
#------------------------------------------#

#Step 12: Create Inflammaging Gene Set.
#       This is the subset of core enrichment genes that are significantly
#       different between old_sed_veh vs yng_sed_veh

Inflammaging_gene_set_12jun<-oldsedveh_v_yngsedveh_genes_12jun %>%
  filter(ENTREZID %in% upreg_old_sed_core_entrez_vec_12jun) %>%
  filter (FDR <0.05)

#Inspect the gene set. Got 63 genes this time. Got 66 genes last time.
#This is because the GSEA function doesn't always return the same exact data.
#So different Core Enrichment genes from subset.
#Can check to see which genes are different than original Inflammaging_symbols vector.
Inflammaging_gene_set_12jun %>%
  pull(Symbol)
#--------------------------------------------#

#Step 13: Create heatmap of full core enrichment genes
#13a. Load in and analyze Raw Count data
#     Upload raw count files into working directory on your computer.

#Load in OLD_SED_VEH-YNG_SED_VEH COUNTS
oldsedveh_vs_yngsedveh_counts_12jun <-read_xlsx("20250219_M007853_Set01_edgeRglm_Counts_OLD_SED_VEH-YNG_SED_VEH.xlsx")
head(oldsedveh_vs_yngsedveh_counts_12jun) # check to see if read in appropriately

# Rename the first column to avoid issues (it may have no name)
colnames(oldsedveh_vs_yngsedveh_counts_12jun)[1] <- "Ensembl"
head(oldsedveh_vs_yngsedveh_counts_12jun) # check to see if rename worked

#create column with decimal removed from Ensembl IDs
oldsedveh_vs_yngsedveh_counts_12jun <- oldsedveh_vs_yngsedveh_counts_12jun %>%
  dplyr::mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+")) #

#Get rid of column with decimal Ensembl ID
oldsedveh_vs_yngsedveh_counts_12jun <- oldsedveh_vs_yngsedveh_counts_12jun %>%
  dplyr::select(-Ensembl)
head(oldsedveh_vs_yngsedveh_counts_12jun) # check to see if rename/remove worked

#Add ENTREZID column mapped from Ensembl. Will returned 1:many mapping between keys and columns warning
oldsedveh_vs_yngsedveh_counts_12jun$ENTREZID = mapIds(org.Mm.eg.db, keys = oldsedveh_vs_yngsedveh_counts_12jun$Ensembl_noDec,
                                                keytype = 'ENSEMBL', column = 'ENTREZID', multiVals = 'first')

#See how many rows/genes there are
dim(oldsedveh_vs_yngsedveh_counts_12jun)

# Remove rows with NA in ENTREZID column
oldsedveh_vs_yngsedveh_counts_12jun <- oldsedveh_vs_yngsedveh_counts_12jun %>% drop_na(ENTREZID)

head(oldsedveh_vs_yngsedveh_counts_12jun)
dim(oldsedveh_vs_yngsedveh_counts_12jun)

#remove duplicate ENTREZID
count_df_clean_12jun <- oldsedveh_vs_yngsedveh_counts_12jun %>%
  group_by(ENTREZID) %>%
  summarise(across(starts_with("T_"), sum), .groups = "drop")

head(count_df_clean_12jun)
dim(count_df_clean_12jun)

# Define expected order
yng_sed_samples <- c("T_05", "T_17", "T_26", "T_32", "T_35", "T_38")
old_sed_samples <- c("T_06", "T_09", "T_13", "T_19", "T_21", "T_25", "T_29", "T_39")

# Combine all sample names
all_samples <- c(yng_sed_samples, old_sed_samples)

# Extract counts matrix (just the sample columns)
counts_matrix <- count_df_clean_12jun %>%
  dplyr::select(all_of(all_samples)) %>%
  as.matrix()

# Set rownames to ENTREZID
rownames(counts_matrix) <- count_df_clean_12jun$ENTREZID

# Create group factor (needed for edgeR object)
group <- factor(c(rep("YNG_SED", length(yng_sed_samples)),
                  rep("OLD_SED", length(old_sed_samples))))

# Create DGEList object
dge <- DGEList(counts = counts_matrix, group = group)

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Convert to logCPM
logCPM_matrix <- cpm(dge, log = TRUE, prior.count = 1)

# View a portion of the matrix
head(logCPM_matrix)
#---------------------------------#

#Step 13b: Create heatmap of oldsedveh v yngsedveh logCPM matrix

#place core enrichment genes in order by FDR (lowest to highest)
ordered_inflam_vec_12jun <-oldsedveh_v_yngsedveh_genes_12jun %>%
  filter(ENTREZID %in% upreg_old_sed_core_entrez_vec_12jun) %>%
  arrange(FDR)%>%
  pull(ENTREZID)


# Define the manually ordered columns (sample groups)
ordered_columns <- c(yng_sed_samples, old_sed_samples)

# Subset and reorder rows (genes) using ordered_inflam_vec
# Note: use intersect to ensure all IDs exist in logCPM_matrix
ordered_entrez_present <- ordered_inflam_vec_12jun[ordered_inflam_vec_12jun %in% rownames(logCPM_matrix)]

# Subset and reorder matrix
heatmap_matrix <- logCPM_matrix[ordered_entrez_present, ordered_columns, drop = FALSE]

# Generate heatmap (preserve row order by setting cluster_rows = FALSE)
heatmap_plot <- pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,              # ❗️Preserves row order
  cluster_cols = FALSE,              # Preserves column order
  scale = "row",                     # Z-score normalize each gene
  show_rownames = FALSE,              # Show gene IDs (ENTREZ)
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Ordered Inflammaging Signature" +
    geom_hline(yintercept = 67.5, color = "black", size = 0.8)
)

# Save to PDF
pdf("ordered_inflammaging_heatmap.pdf", width = 3.5, height = 6)
print(heatmap_plot)
dev.off()

# Print to RStudio plot window
print(heatmap_plot)
#------------------------------------------#

#Step 14: Create Inflammaging Gene Heatmap of oldsedveh vs yngsedveh

inflammaging_entrez_vec_12jun <- Inflammaging_gene_set_12jun$ENTREZID %>% unname()

#14a.----Get Inflammaging gene set and convert to Symbols------#
# Make sure inflammaging_gene_vec is character if rownames are character
inflammaging_gene_vec_12jun <- as.character(inflammaging_entrez_vec_12jun)

# Subset the matrix using rownames
Inflamm_logCPM_matrix <- logCPM_matrix[rownames(logCPM_matrix) %in% inflammaging_entrez_vec_12jun, , drop = FALSE]

# View the result
head(Inflamm_logCPM_matrix)
dim(Inflamm_logCPM_matrix)

#14b------Convert to Symbols-----------------#

# Convert ENTREZIDs (rownames of logCPM_matrix) to gene symbols
entrez_ids <- rownames(Inflamm_logCPM_matrix)

# Map ENTREZIDs to gene symbols
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Remove any NA entries (optional: keeps only matched genes)
valid_idx <- which(!is.na(gene_symbols))
Inflamm_logCPM_matrix <- Inflamm_logCPM_matrix[valid_idx, , drop = FALSE]
rownames(Inflamm_logCPM_matrix) <- gene_symbols[valid_idx]

# View result
head(Inflamm_logCPM_matrix)

#14c.---Create Inflammaging GeneSet Heatmap-----------------------#

# Reorder columns manually
ordered_columns <- c(yng_sed_samples, old_sed_samples)

# Subset logCPM matrix to genes of interest and reorder columns
inflamm_geneset_heatmap <- Inflamm_logCPM_matrix[rownames(Inflamm_logCPM_matrix), ordered_columns, drop = FALSE]

# Generate heatmap with manual column order

inflamm_geneset_heatplot <- pheatmap(inflamm_geneset_heatmap, 
                                     cluster_rows = TRUE, 
                                     cluster_cols = FALSE, 
                                     scale = "row",  # Standardize by row (gene)
                                     show_rownames = TRUE, 
                                     show_colnames = TRUE, 
                                     color = colorRampPalette(c("blue", "white", "red"))(50))

# Save as PDF
pdf("inflamm_geneset_heatplot.pdf", width = 3.5, height = 8)  # Adjust height
print(inflamm_geneset_heatplot)
dev.off()

# Print in RStudio
print(inflamm_geneset_heatplot)
