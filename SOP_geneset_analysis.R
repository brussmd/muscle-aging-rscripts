#This code is for analyzing geneset enrichment in a gene dataset.
#For the Inflammaging analysis, this are datasets comparing oldsedveh vs old interventions.
#The datasets are organized by logFC of old intervenion vs oldsedveh
#Creating the barplot comparison for NES vs oldsedveh was done in Prism
#using the values from #5. of each analysis.

#Step 1: Check for necessary packages if haven't already.
required_pkgs <- c(
  "tidyverse",      # covers dplyr, readr, ggplot2, stringr, etc.
  "readxl",
  "ggrepel",
  "pheatmap",
  "VennDiagram",
  "ggVennDiagram",
  "eulerr",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "DOSE",
  "fgsea",
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
#---------------------------------------#

#Step 2: Load libraries.
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
#--------------------------#

#Step 3: Initialize functions

#-------Read in xlsx, Clean Dataframe and add ENTREZID col-----------#
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
#----------------------------------------------#
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
#-------------------------------------------------------#

#Step 4: Create your geneset as an ENTREZID vector
#An example of this code is shown in SOP_code_Inflammaging_Geneset
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


#---------------------------------------------------------#
#----Run Inflammaging GSEA on Interventions vs OldSedVeh--#
#---------------------------------------------------------#

#-----------------------------------------------------------------#
#----OldPwrVeh vs OldSedVeh Inflammaging GeneSet GSEA Analysis----#
oldpwrveh_v_oldsedveh_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_PWR_VEH-OLD_SED_VEH_logFC",
  FDR_col = "OLD_PWR_VEH-OLD_SED_VEH_FDR")

oldpwrveh_v_oldsedveh_genes %>%
  filter(FDR <0.05)

oldpwrveh_v_oldsedveh_genes %>%
  filter(ENTREZID %in% inflammaging_gene_vec)

create_vec_from_df(oldpwrveh_v_oldsedveh_genes, "oldpwrveh_v_oldsedveh")
head(oldpwrveh_v_oldsedveh.vec)

# Run GSEA using your custom gene set
gsea_inflamm_oldpwrveh_v_oldsedveh <- GSEA(geneList = oldpwrveh_v_oldsedveh.vec,
                                           TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                                           pvalueCutoff = 1,
                                           verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldpwrveh_v_oldsedveh@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#-------------------------------------------------#

#-----------------------------------------------------------------#
#----OldSedFRAP vs OldSedVeh Inflammaging GeneSet GSEA Analysis----#
oldsedfrap_v_oldsedveh_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_FRAP-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_SED_FRAP-OLD_SED_VEH_logFC",
  FDR_col = "OLD_SED_FRAP-OLD_SED_VEH_FDR")

oldsedfrap_v_oldsedveh_genes %>%
  filter(FDR <0.05)

oldsedfrap_v_oldsedveh_genes %>%
  filter(ENTREZID %in% inflammaging_gene_vec)

create_vec_from_df(oldsedfrap_v_oldsedveh_genes, "oldsedfrap_v_oldsedveh")
head(oldsedfrap_v_oldsedveh.vec)

# Run GSEA using your custom gene set
gsea_inflamm_oldsedfrap_v_oldsedveh <- GSEA(geneList = oldsedfrap_v_oldsedveh.vec,
                                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                                            pvalueCutoff = 1,
                                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldsedfrap_v_oldsedveh@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#--------------------------------------------------#

#-----------------------------------------------------------------#
#----OldSedIRAP vs OldSedVeh Inflammaging GeneSet GSEA Analysis----#
oldsedirap_v_oldsedveh_genes <- read_clean_xlsx(
  xlsx_file = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_IRAP-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_SED_IRAP-OLD_SED_VEH_logFC",
  FDR_col = "OLD_SED_IRAP-OLD_SED_VEH_FDR")

oldsedirap_v_oldsedveh_genes %>%
  filter(FDR <0.05)

oldsedirap_v_oldsedveh_genes %>%
  filter(ENTREZID %in% inflammaging_gene_vec)

create_vec_from_df(oldsedirap_v_oldsedveh_genes, "oldsedirap_v_oldsedveh")
head(oldsedirap_v_oldsedveh.vec)

# Run GSEA using your custom gene set
gsea_inflamm_oldsedirap_v_oldsedveh <- GSEA(geneList = oldsedirap_v_oldsedveh.vec,
                                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                                            pvalueCutoff = 1,
                                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldsedirap_v_oldsedveh@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#--------------------------------------------------#

#-----------------------------------------------------------------#
#----OldpwrIRAP vs OldSedVeh Inflammaging GeneSet GSEA Analysis----#
oldpwrirap_v_oldsedveh_genes <- read_clean_xlsx(
  xlsx_file = "20250530_M007853_Set05_edgeRglm_GENE_OLD_PWR_IRAP-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_PWR_IRAP-OLD_SED_VEH_logFC",
  FDR_col = "OLD_PWR_IRAP-OLD_SED_VEH_FDR")

oldpwrirap_v_oldsedveh_genes %>%
  filter(FDR <0.05)

oldpwrirap_v_oldsedveh_genes %>%
  filter(ENTREZID %in% inflammaging_gene_vec)

create_vec_from_df(oldpwrirap_v_oldsedveh_genes, "oldpwrirap_v_oldsedveh")
head(oldsedirap_v_oldsedveh.vec)

# Run GSEA using your custom gene set
gsea_inflamm_oldpwrirap_v_oldsedveh <- GSEA(geneList = oldpwrirap_v_oldsedveh.vec,
                                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                                            pvalueCutoff = 1,
                                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldpwrirap_v_oldsedveh@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#--------------------------------------------------#

#-----------------------------------------------------------------#
#----OldpwrFRAP vs OldSedVeh Inflammaging GeneSet GSEA Analysis----#
oldpwrfrap_v_oldsedveh_genes <- read_clean_xlsx(
  xlsx_file = "20250530_M007853_Set05_edgeRglm_GENE_OLD_PWR_FRAP-OLD_SED_VEH.xlsx",
  logFC_col = "OLD_PWR_FRAP-OLD_SED_VEH_logFC",
  FDR_col = "OLD_PWR_FRAP-OLD_SED_VEH_FDR")

oldpwrfrap_v_oldsedveh_genes %>%
  filter(FDR <0.05)

oldpwrfrap_v_oldsedveh_genes %>%
  filter(ENTREZID %in% inflammaging_gene_vec)

create_vec_from_df(oldpwrfrap_v_oldsedveh_genes, "oldpwrfrap_v_oldsedveh")
head(oldsedfrap_v_oldsedveh.vec)

# Run GSEA using your custom gene set
gsea_inflamm_oldpwrfrap_v_oldsedveh <- GSEA(geneList = oldpwrfrap_v_oldsedveh.vec,
                                            TERM2GENE = data.frame(term = "Inflamm", gene = inflammaging_gene_vec),
                                            pvalueCutoff = 1,
                                            verbose = FALSE)

# 5. View results: NES and p-value
gsea_inflamm_oldpwrfrap_v_oldsedveh@result %>% 
  dplyr::select(ID, NES, pvalue, p.adjust)
#////////////////////////////////////////////#
#////////////////////////////////////////////#

