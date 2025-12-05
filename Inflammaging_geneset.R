#Creating Inflammaging Gene Set

#Load Libraries
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


#-----------------Read in OldSed vs YngSed Genes----------------------------#
oldsed_vs_yngsed_genes <-read_xlsx("20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx")

oldsed_vs_yngsed_genes

oldsed_vs_yngsed_genes <- oldsed_vs_yngsed_genes %>%
  rename(
    old_sed_logFC = `OLD_SED_VEH-YNG_SED_VEH_logFC`,
    old_sed_FDR   = `OLD_SED_VEH-YNG_SED_VEH_FDR`
  )

# now your original lines work:
oldsed_vs_yngsed_genes.vec <- oldsed_vs_yngsed_genes %>%
  arrange(desc(old_sed_logFC)) %>%
  select(ENTREZID, old_sed_logFC) %>%
  tidyr::drop_na() %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  tibble::deframe()

#--------------ENTREZID Vector---------------------#
oldsed_vs_yngsed_genes.vec = oldsed_vs_yngsed_genes %>% 
  dplyr::arrange(., desc(OLD_SED_VEH-YNG_SED_VEH_logFC)) %>% 
  dplyr::select(., ENTREZID, OLD_SED_VEH-YNG_SED_VEH_logFC) %>% 
  tidyr::drop_na() %>% 
  dplyr::distinct(., ENTREZID, .keep_all = T) %>% 
  tibble::deframe()

head(oldsed_vs_yngsed_genes.vec)

#---------------------------------------------------#
#---------Simplified KEGG OldSed vs YngSed----------#
#Simplify
sig_gseGO_oldsed_vs_yngsed.OUTPUT = gseGO(geneList = oldsed_vs_yngsed_genes.vec, 
                                          ont = 'BP',
                                          OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 300,
                                          eps = 1e-30, pvalueCutoff = 0.05)

dim(sig_gseGO_oldsed_vs_yngsed.OUTPUT)

sig_gseGO_oldsed_vs_yngsed.df = sig_gseGO_oldsed_vs_yngsed.OUTPUT@result %>% 
  dplyr::mutate(., Count = str_count(.$core_enrichment, '/')+1)

sig_gseGO_oldsed_vs_yngsed.df


# Run simplify() for old power genes
oldsed_vs_yngsed_simplified <- simplify(sig_gseGO_oldsed_vs_yngsed.OUTPUT, cutoff = 0.5, 
                                        by = "p.adjust", select_fun = min)

dim(oldsed_vs_yngsed_simplified)
oldsed_vs_yngsed_simplified

oldsed_vs_yngsed_simpl.df = oldsed_vs_yngsed_simplified@result %>% 
  dplyr::mutate(., Count = str_count(.$core_enrichment, '/')+1)

# Modify the description for the specific GO term
oldsed_vs_yngsed_simpl.df <- oldsed_vs_yngsed_simpl.df  %>%
  mutate(Description = ifelse(ID == "GO:0043280", "endopeptidase involved in apoptotic process", Description))

# Check if the change was applied correctly
oldsed_vs_yngsed_simpl.df %>%
  filter(ID == "GO:0043280")

upreg_oldsed_simplkegg <- oldsed_vs_yngsed_simpl.df %>%
  filter (NES > 0)

oldsed_vs_yngsed_simpl.df %>%
  filter(NES > 0) %>%
  pull(Description)

#--------------Visualize KEGG Set-------------------#

nested_upreg_oldsed_simplkegg_plot <- upreg_oldsed_simplkegg %>%
  ggplot(., aes(x = NES, y=reorder(Description, NES), size = Count, col = p.adjust )) +
  geom_point() +
  theme_bw()  +
  theme(axis.text.y = element_text(size = 7),  # Adjust y-axis label size
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.x = element_text(size = 7))  # Adjust x-axis label size

# Save as PDF
pdf("nested_upreg_oldsed_simplkegg_plot.pdf", width = 8, height = 6)  # Adjust height
print(nested_upreg_oldsed_simplkegg_plot)
dev.off()

# Print in RStudio
print(nested_upreg_oldsed_simplkegg_plot)

oldsed_vs_yngsed_genes

#------------Get Core Enrichment Genes--------------------------#
# Extract all ENTREZ IDs from the 'core_enrichment' column, split them into a vector
upreg_old_sed_core_entrez_vec <- upreg_oldsed_simplkegg$core_enrichment %>%
  str_split(pattern = "/") %>%
  unlist()

# Get unique ENTREZ IDs and sort them
upreg_unq_old_sed_core_entrez_vec <- unique(upreg_old_sed_core_entrez_vec) %>%
  as.numeric() %>%  # Convert from character to numeric
  sort()

length(upreg_unq_old_sed_core_entrez_vec)
head(upreg_unq_old_sed_core_entrez_vec)

oldsed_vs_yngsed_genes %>%
  filter(ENTREZID %in% upreg_unq_old_sed_core_entrez_vec) %>%
  filter (old_sed_FDR <0.05)

#Get vector for all core enrichment genes from significant upregulated KEGG
#pathways that are in the original dataset and have FDR < 0.05.
inflammaging_gene_vec <- oldsed_vs_yngsed_genes %>%
  filter(ENTREZID %in% upreg_unq_old_sed_core_entrez_vec) %>%
  filter(old_sed_FDR < 0.05) %>%
  pull(ENTREZID)

head(inflammaging_gene_vec)

length(inflammaging_gene_vec)

# Create a gene set list formatted for GSEA
Inflammaging_geneSet <- list(Inflamm = inflammaging_gene_vec)
#--------------------------------------------------------------#

#Create ALL inflammaging gene vector
all_inflammaging_gene_vec <- oldsed_vs_yngsed_genes %>%
  filter(ENTREZID %in% upreg_unq_old_sed_core_entrez_vec) %>%
  pull(ENTREZID)

head(all_inflammaging_gene_vec)

length(all_inflammaging_gene_vec)

#Create ALL inflammaging gene set list for GSEA
All_Inflammaging_geneSet <- list(All_Inflamm = all_inflammaging_gene_vec)

#--------------------------------------------------------------#
#-------------Create upregulated inflammaging gene set---------#
#--------------------------------------------------------------#
#Note...All 500 genes for inflammaging are upregulated in oldsed vs yngsed
#Create ALL inflammaging gene vector
all_upreg_inflammaging_gene_vec <- oldsed_vs_yngsed_genes %>%
  filter(ENTREZID %in% upreg_unq_old_sed_core_entrez_vec) %>%
  filter(old_sed_logFC >0) %>%
  pull(ENTREZID)

head(all_upreg_inflammaging_gene_vec)

length(all_upreg_inflammaging_gene_vec)

#Create ALL inflammaging gene set list for GSEA
All_Inflammaging_geneSet <- list(All_Inflamm = all_inflammaging_gene_vec)

#----------------------------------------#
# Run GSEA using your custom gene set
gsea_inflamm_oldsed <- GSEA(geneList = oldsed_vs_yngsed_genes.vec,
                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                            pvalueCutoff = 1,
                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldsed@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)

oldsed_vs_yngsed_genes %>%
  filter (ENTREZID %in% inflammaging_gene_vec) %>%
  filter (old_sed_FDR <0.05)
#-----------------------------------------------#
#-----------All Inflammaging Gene Set GSEA (500 genes)-----#

# Run GSEA using your custom gene set
gsea_all_inflamm_oldsed <- GSEA(geneList = oldsed_vs_yngsed_genes.vec,
                            TERM2GENE = data.frame(term = "All_Inflamm", gene = all_inflammaging_gene_vec),
                            pvalueCutoff = 1,
                            verbose = FALSE)

# View results: NES and p-value
gsea_all_inflamm_oldsed@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)

#---------------------------------------------#
#--------Test effect of PWR on Inflammaging---#
#---------------------------------------------#

head(oldpwr_vs_yngsed_genes.vec)

gsea_inflamm_oldpwr <- GSEA(geneList = oldpwr_vs_yngsed_genes.vec,
                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                            pvalueCutoff = 1,
                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldpwr@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#----------------------------------------------#
# Run GSEA using your custom gene set
gsea_all_inflamm_oldpwr <- GSEA(geneList = oldpwr_vs_yngsed_genes.vec,
                                TERM2GENE = data.frame(term = "All_Inflamm", gene = all_inflammaging_gene_vec),
                                pvalueCutoff = 1,
                                verbose = FALSE)

# View results: NES and p-value
gsea_all_inflamm_oldpwr@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#-------------------------------------------#
head(oldpwr_vs_yngsed_genes)
oldpwr_vs_yngsed_genes %>%
  filter (ENTREZID %in% inflammaging_gene_vec) %>%
  filter (old_pwr_FDR <0.05)


#---------------------------------------------#
#--------Test effect of IRAP on Inflammaging---#
#---------------------------------------------#

gsea_inflamm_oldsedirap <- GSEA(geneList = oldsedirap_vs_yngsed_genes.vec,
                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                            pvalueCutoff = 1,
                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldsedirap@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)

#----------------------------------------------------------------#
#------------Test OldPwr vs OldSed on Inflammaging gene sets-----#
#----------------------------------------------------------------#

#-----------------Read in OldPwr vs OldSed Genes----------------------------#
oldpwr_vs_oldsed_genes <-read_xlsx("20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx")
oldpwr_vs_oldsed_genes

#---------Clean Dataframe for OldPwr vs OldSed--------------------------------#
oldpwr_vs_oldsed_genes <- oldpwr_vs_oldsed_genes %>%
  dplyr::select(Ensembl, Symbol, "oldpwr_v_oldsed_logFC" = `OLD_PWR_VEH-OLD_SED_VEH_logFC`, 
                "oldpwr_v_oldsed_FDR" = `OLD_PWR_VEH-OLD_SED_VEH_FDR`)
colnames(oldpwr_vs_oldsed_genes)
oldpwr_vs_oldsed_genes

oldpwr_vs_oldsed_genes <- oldpwr_vs_oldsed_genes %>%
  dplyr::mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"))

oldpwr_vs_oldsed_genes <- oldpwr_vs_oldsed_genes %>%
  dplyr::select(-Ensembl)

oldpwr_vs_oldsed_genes$ENTREZID = mapIds(org.Mm.eg.db, keys = oldpwr_vs_oldsed_genes$Ensembl_noDec,
                                             keytype = 'ENSEMBL', column = 'ENTREZID', multiVals = 'first')

oldpwr_vs_oldsed_genes <- as.data.frame(oldpwr_vs_oldsed_genes)

oldpwr_vs_oldsed_genes

# Remove rows with NA in ENTREZID column
oldpwr_vs_oldsed_genes <- oldpwr_vs_oldsed_genes %>% drop_na(ENTREZID)

# Check the updated dataframe
oldpwr_vs_oldsed_genes %>%
  filter (oldpwr_v_oldsed_FDR <0.05)
#-----------------------------#

#--------------ENTREZID Vector---------------------#
oldpwr_vs_oldsed_genes.vec = oldpwr_vs_oldsed_genes %>% 
  dplyr::arrange(., desc(oldpwr_v_oldsed_logFC)) %>% 
  dplyr::select(., ENTREZID, oldpwr_v_oldsed_logFC) %>% 
  tidyr::drop_na() %>% 
  dplyr::distinct(., ENTREZID, .keep_all = T) %>% 
  tibble::deframe()

head(oldpwr_vs_oldsed_genes.vec)
#-----------------------------------#
#--------oldpwr vs oldsed significant inflammaging gene set-------------#
# Run GSEA using your custom gene set
gsea_inflamm_oldpwr_v_oldsed <- GSEA(geneList = oldpwr_vs_oldsed_genes.vec,
                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                            pvalueCutoff = 1,
                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldpwr_v_oldsed@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)

oldpwr_vs_oldsed_genes %>%
  filter (ENTREZID %in% inflammaging_gene_vec)
#------------------------------------------------#

#--------oldpwr vs oldsed All inflammaging gene set-------------#
head(all_inflammaging_gene_vec)
length(all_inflammaging_gene_vec)
oldpwr_vs_oldsed_genes %>%
  filter (ENTREZID %in% all_inflammaging_gene_vec)

# Run GSEA using your custom gene set
gsea_all_inflamm_oldpwr_v_oldsed <- GSEA(geneList = oldpwr_vs_oldsed_genes.vec,
                                     TERM2GENE = data.frame(term = "All_Inflamm", gene = all_inflammaging_gene_vec),
                                     pvalueCutoff = 1,
                                     verbose = FALSE,
                                     nPermSimple = 10000)

# 5. View results: NES and p-value
gsea_all_inflamm_oldpwr_v_oldsed@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)

length(intersect(names(oldpwr_vs_oldsed_genes.vec), all_inflammaging_gene_vec))
#-------------------------------------------------------------------#
#Oldpwr has negative NES vs OldSed on inflammaging, suggesting that these genes are turned down with pwr.
#First determine if this effect holds with SenMayo
#Now need to verify that you don't also see this with IRAP

# Run GSEA using SenMayo on OldPwr vs OldSed
gsea_senmayo_oldpwr_v_oldsed <- GSEA(geneList = oldpwr_vs_oldsed_genes.vec,
                     TERM2GENE = data.frame(term = "SenMayo", gene = SenMayo_entrez),
                     pvalueCutoff = 1,
                     verbose = FALSE)

# 5. View results: NES and p-value
gsea_senmayo_oldpwr_v_oldsed@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#--------------This didn't work Actually higher with power--------------#

#----------------------------------------------------------------#
#-------------OldIrap vs OldSed Control--------------------------#
#-----------------Read in OldIRAP vs OldSed Genes----------------------------#
oldirap_vs_oldsed_genes <-read_xlsx("20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_IRAP-OLD_SED_VEH.xlsx")
oldirap_vs_oldsed_genes

#---------Clean Dataframe for OldPwr vs OldSed--------------------------------#
oldirap_vs_oldsed_genes <- oldirap_vs_oldsed_genes %>%
  dplyr::select(Ensembl, Symbol, "oldirap_v_oldsed_logFC" = `OLD_SED_IRAP-OLD_SED_VEH_logFC`, 
                "oldirap_v_oldsed_FDR" = `OLD_SED_IRAP-OLD_SED_VEH_FDR`)
colnames(oldirap_vs_oldsed_genes)
oldirap_vs_oldsed_genes

oldirap_vs_oldsed_genes <- oldirap_vs_oldsed_genes %>%
  dplyr::mutate(Ensembl_noDec = str_remove(Ensembl, "\\..+"))

oldirap_vs_oldsed_genes <- oldirap_vs_oldsed_genes %>%
  dplyr::select(-Ensembl)

oldirap_vs_oldsed_genes$ENTREZID = mapIds(org.Mm.eg.db, keys = oldirap_vs_oldsed_genes$Ensembl_noDec,
                                         keytype = 'ENSEMBL', column = 'ENTREZID', multiVals = 'first')

oldirap_vs_oldsed_genes <- as.data.frame(oldirap_vs_oldsed_genes)

oldirap_vs_oldsed_genes

# Remove rows with NA in ENTREZID column
oldirap_vs_oldsed_genes <- oldirap_vs_oldsed_genes %>% drop_na(ENTREZID)

# Check the updated dataframe
oldirap_vs_oldsed_genes %>%
  filter (oldirap_v_oldsed_FDR <0.05)

oldirap_vs_oldsed_genes %>%
  filter(ENTREZID %in% inflammaging_gene_vec)
#------------------------------------#

#--------------ENTREZID Vector---------------------#
oldpwr_vs_oldsed_genes.vec = oldpwr_vs_oldsed_genes %>% 
  dplyr::arrange(., desc(oldpwr_v_oldsed_logFC)) %>% 
  dplyr::select(., ENTREZID, oldpwr_v_oldsed_logFC) %>% 
  tidyr::drop_na() %>% 
  dplyr::distinct(., ENTREZID, .keep_all = T) %>% 
  tibble::deframe()

head(oldpwr_vs_oldsed_genes.vec)
str(oldpwr_vs_oldsed_genes.vec)

#-----------------------------------#
#--------oldpwr vs oldsed significant inflammaging gene set-------------#
# Run GSEA using your custom gene set
gsea_inflamm_oldpwr_v_oldsed <- GSEA(geneList = oldpwr_vs_oldsed_genes.vec,
                                     TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                                     pvalueCutoff = 1,
                                     verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldpwr_v_oldsed@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)

oldpwr_vs_oldsed_genes %>%
  filter (ENTREZID %in% inflammaging_gene_vec)




# 4. Run GSEA using your custom gene set
gsea_senmayo_oldsedirap <- GSEA(geneList = oldsedirap_vs_yngsed_genes.vec,
                                TERM2GENE = data.frame(term = "SenMayo", gene = SenMayo_entrez),
                                pvalueCutoff = 1,
                                verbose = FALSE)

# 5. View results: NES and p-value
gsea_senmayo_oldsedirap@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#------------------------------------------#
#---Messing Around-------------------------#

upreg_oldsed_simplkegg$core_enrichment
char_vec <-upreg_unq_old_sed_core_entrez_vec %>%
  as.character()

full_inflamm_geneset <- mapIds(
  org.Mm.eg.db,
  keys = char_vec,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

full_inflamm_geneset %>%
  arrange()

Inflammaging_symbols

#---Notes-------#
#create vector of these gene symbols and run separate gsea
#Use SenMayo as example
#Create logcpm heatmap
#test oldpwr vs oldsed with inflammaging geneset
#oldsed inflammaging genes fdr <0.05 => 66 genes
#oldpwr inflammaging genes fdr <0.05 => 32 genes
#run gsea on 34 genes.
#create function for running gsea on specific set of genes