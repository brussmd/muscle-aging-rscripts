#GO Comparison of all interventions to OldSedVeh


library(readxl)
library(dplyr)
library(stringr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyr)
library(BiocParallel)

#-------------FUNCTIONS----------------------------------------------#
#-------Read in xlsx, Clean Dataframe and add ENTREZID col-----------#
#-------Read in xlsx, Clean Dataframe and add ENTREZID col-----------#
read_clean_xlsx_simple <- function(xlsx_file) {
  # Read in the xlsx file
  df <- read_xlsx(xlsx_file)
  
  # Automatically detect columns
  logFC_col <- names(df)[grepl("logFC", names(df), ignore.case = TRUE)][1]
  FDR_col   <- names(df)[grepl("FDR", names(df), ignore.case = TRUE)][1]
  LR_col<- names(df)[grepl("LR", names(df), ignore.case = TRUE)][1]
  PValue_col <- names(df)[grepl("PValue", names(df), ignore.case = TRUE)][1]
  
  # Check if found
  if (is.na(logFC_col) || is.na(FDR_col)) {
    stop("Could not find logFC or FDR columns automatically.")
  }
  
  # Clean and rename columns
  df <- df %>%
    dplyr::select(
      Ensembl,
      Symbol,
      logFC = all_of(logFC_col),
      FDR = all_of(FDR_col),
      LR = all_of(LR_col),
      PValue = all_of(PValue_col)
    ) %>%
    dplyr::mutate(
      Ensembl_noDec = str_remove(Ensembl, "\\..+"),
      ENTREZID = mapIds(
        org.Mm.eg.db,
        keys = Ensembl_noDec,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    ) %>%
    dplyr::select(-Ensembl) %>%
    tidyr::drop_na(ENTREZID) %>%
    as.data.frame()
  
  return(df)
}
#----------------------------------------------#
#----------Create named vector from dataframe function-----------#
create_vec_from_df <- function(df, df_name) {
  
  vec <- df %>%
    # Ranking metric: logFC × (–log10(PValue))^0.25
    mutate(
      ranking_metric = logFC * ((-log10(PValue + 1e-300))^0.25)
    ) %>%
    
    arrange(desc(ranking_metric)) %>%
    dplyr::select(ENTREZID, ranking_metric) %>%
    tidyr::drop_na() %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    tibble::deframe()
  
  assign(paste0(df_name, ".vec"), vec, envir = .GlobalEnv)
}

#---------------------------------------------------------------#

#///////////////////////////////////////////////////////////////#
#------Get GSEA for each intervention data set------------------#
#///////////////////////////////////////////////////////////////#


#-----oldsedveh vs yngsedveh------------------------------------#
#Function #1
oldsedveh_v_yngsedveh_genes_22Nov <- read_clean_xlsx_simple(
  xlsx_file = "20250219_M007853_Set01_edgeRglm_GENE_OLD_SED_VEH-YNG_SED_VEH.xlsx")

oldsedveh_v_yngsedveh_genes_22Nov %>%
  filter(FDR <0.05)
head(oldsedveh_v_yngsedveh_genes_22Nov)

#Function 2.
create_vec_from_df(oldsedveh_v_yngsedveh_genes_22Nov, "oldsedveh_v_yngsedveh_genes_22Nov")
head(oldsedveh_v_yngsedveh_genes_22Nov.vec)
length(oldsedveh_v_yngsedveh_genes_22Nov.vec)

# 1. Force single-core execution
register(SerialParam())

# 2. Sort the gene list identically each time
oldsedveh_v_yngsedveh_genes_22Nov.vec <- sort(oldsedveh_v_yngsedveh_genes_22Nov.vec, decreasing = TRUE)

# 3. Set seed before GSEA
set.seed(20241122)

# 4. Run GSEA
gseGO_oldsedveh_v_yngsedveh_22Nov.OUTPUT <- gseGO(
  geneList = oldsedveh_v_yngsedveh_genes_22Nov.vec,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  minGSSize = 10,
  maxGSSize = 300,
  eps = 1e-30,
  pvalueCutoff = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)

gseGO_oldsedveh_v_yngsedveh_22Nov.OUTPUT

oldsed_vs_yngsed_22Nov.df <- as.data.frame(gseGO_oldsedveh_v_yngsedveh_22Nov.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

head(oldsed_vs_yngsed_22Nov.df)
write_csv(oldsed_vs_yngsed_22Nov.df, "oldsed_vs_yngsed_22Nov_df.csv")

old_v_yng_all_GO_dotplot <- oldsed_vs_yngsed_22Nov.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
print(old_v_yng_all_GO_dotplot)
# Save (optional)
ggsave("old_v_yng_all_GO_dotplot.pdf", old_v_yng_all_GO_dotplot, width = 8, height = 8)
dim(oldsed_vs_yngsed_22Nov.df)
#================================================================================#
#================================================================================#

#///////////OldPwrVeh vs OldSedVeh//////////////////////////////////#

#Function #1
oldpwrveh_v_oldsedveh_genes_22Nov <- read_clean_xlsx_simple(
  xlsx_file = "20250219_M007853_Set02_edgeRglm_GENE_OLD_PWR_VEH-OLD_SED_VEH.xlsx")

head(oldpwrveh_v_oldsedveh_genes_22Nov)

oldpwrveh_v_oldsedveh_genes_22Nov %>%
  filter(FDR <0.05)

#Function 2.
create_vec_from_df(oldpwrveh_v_oldsedveh_genes_22Nov, "oldpwrveh_v_oldsedveh_genes_22Nov")
head(oldpwrveh_v_oldsedveh_genes_22Nov.vec)

# 1. Force single-core execution
register(SerialParam())

# 2. Sort the gene list identically each time
oldpwrveh_v_oldsedveh_genes_22Nov.vec <- sort(oldpwrveh_v_oldsedveh_genes_22Nov.vec, decreasing = TRUE)

# 3. Set seed before GSEA
set.seed(20241122)

# 4. Run GSEA
gseGO_oldpwrveh_v_oldsedveh_22Nov.OUTPUT <- gseGO(
  geneList = oldpwrveh_v_oldsedveh_genes_22Nov.vec,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  minGSSize = 10,
  maxGSSize = 300,
  eps = 1e-30,
  pvalueCutoff = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)


oldpwrveh_v_oldsedveh_22Nov.df <- as.data.frame(gseGO_oldpwrveh_v_oldsedveh_22Nov.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

head(oldpwrveh_v_oldsedveh_22Nov.df)
dim(oldpwrveh_v_oldsedveh_22Nov.df)
write_csv(oldpwrveh_v_oldsedveh_22Nov.df, "oldpwrveh_v_oldsedveh_22Nov_df.csv")

oldpwrveh_v_oldsedveh_all_GO_dotplot <- oldpwrveh_v_oldsedveh_22Nov.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
print(oldpwrveh_v_oldsedveh_all_GO_dotplot)
# Save (optional)
ggsave("oldpwrveh_v_oldsedveh_GO_dotplot.pdf", oldpwrveh_v_oldsedveh_GO_dotplot, width = 8, height = 8)
#=======================================================================#
#=======================================================================#

#--------Find common GO from age and PWR--------------------------------#

head(oldsed_vs_yngsed_22Nov.df)

sub_GO_oldsed_vs_yngsed_22Nov <- oldsed_vs_yngsed_22Nov.df %>%
  dplyr::select(ID, Description, NES, p.adjust, Count)

head(sub_GO_oldsed_vs_yngsed_22Nov)

sub_GO_oldpwrveh_vs_oldsedveh_22Nov <- oldpwrveh_v_oldsedveh_22Nov.df %>%
  dplyr::select(ID, Description, NES, p.adjust, Count)

head(sub_GO_oldpwrveh_vs_oldsedveh_22Nov)

head(sub_GO_oldsed_vs_yngsed_22Nov)
dim(sub_GO_oldsed_vs_yngsed_22Nov)

# Vector of IDs from each df
ids1 <- sub_GO_oldpwrveh_vs_oldsedveh_22Nov$ID
ids2 <- sub_GO_oldsed_vs_yngsed_22Nov$ID

# Overlapping IDs
overlap_ids <- intersect(ids1, ids2)

# Number of overlapping IDs
length(overlap_ids)
overlap_ids

oldsed_vs_yngsed_22Nov.df %>%
  filter(ID %in% overlap_ids)

#======================================================#
#======================================================#

#///////////OldPwrIrap vs OldSedVeh//////////////////////////////////#

#Function #1
oldpwrirap_v_oldsedveh_genes_22Nov <- read_clean_xlsx_simple(
  xlsx_file = "20250530_M007853_Set05_edgeRglm_GENE_OLD_PWR_IRAP-OLD_SED_VEH.xlsx")

oldpwrirap_v_oldsedveh_genes_22Nov %>%
  filter(FDR <0.05)

#Function 2.
create_vec_from_df(oldpwrirap_v_oldsedveh_genes_22Nov, "oldpwrirap_v_oldsedveh_genes_22Nov")
head(oldpwrirap_v_oldsedveh_genes_22Nov.vec)

# 1. Force single-core execution
register(SerialParam())

# 2. Sort the gene list identically each time
oldpwrirap_v_oldsedveh_genes_22Nov.vec <- sort(oldpwrirap_v_oldsedveh_genes_22Nov.vec, decreasing = TRUE)

# 3. Set seed before GSEA
set.seed(20241122)

# 4. Run GSEA
gseGO_oldpwrirap_v_oldsedveh_22Nov.OUTPUT <- gseGO(
  geneList = oldpwrirap_v_oldsedveh_genes_22Nov.vec,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  minGSSize = 10,
  maxGSSize = 300,
  eps = 1e-30,
  pvalueCutoff = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)


oldpwrirap_v_oldsedveh_22Nov.df <- as.data.frame(gseGO_oldpwrirap_v_oldsedveh_22Nov.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

head(oldpwrirap_v_oldsedveh_22Nov.df)
dim(oldpwrirap_v_oldsedveh_22Nov.df)
write_csv(oldpwrirap_v_oldsedveh_22Nov.df, "oldpwrirap_v_oldsedveh_22Nov_df.csv")

oldpwrirap_v_oldsedveh_all_GO_dotplot <- oldpwrirap_v_oldsedveh_22Nov.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
print(oldpwrirap_v_oldsedveh_all_GO_dotplot)
# Save (optional)
ggsave("oldpwrirap_v_oldsedveh_all_GO_dotplot.pdf", oldpwrirap_v_oldsedveh_all_GO_dotplot, width = 8, height = 8)

# Vector of IDs from each df
pwrirap_ids <- oldpwrirap_v_oldsedveh_22Nov.df$ID
oldsed_ids <- oldsed_vs_yngsed_22Nov.df$ID

# Overlapping IDs
pwrirap_overlap_ids <- intersect(pwrirap_ids, oldsed_ids)

# Number of overlapping IDs
length(pwrirap_overlap_ids)
pwrirap_overlap_ids

oldsed_vs_yngsed_22Nov.df %>%
  filter(ID %in% pwrirap_overlap_ids)
#==================================================#
#==================================================#

#///////////OldPwrFrap vs OldSedVeh//////////////////////////////////#

#Function #1
oldpwrfrap_v_oldsedveh_genes_22Nov <- read_clean_xlsx_simple(
  xlsx_file = "20250530_M007853_Set05_edgeRglm_GENE_OLD_PWR_FRAP-OLD_SED_VEH.xlsx")

oldpwrfrap_v_oldsedveh_genes_22Nov %>%
  filter(FDR <0.05)

#Function 2.
create_vec_from_df(oldpwrfrap_v_oldsedveh_genes_22Nov, "oldpwrfrap_v_oldsedveh_genes_22Nov")
head(oldpwrfrap_v_oldsedveh_genes_22Nov.vec)

# 1. Force single-core execution
register(SerialParam())

# 2. Sort the gene list identically each time
oldpwrfrap_v_oldsedveh_genes_22Nov.vec <- sort(oldpwrfrap_v_oldsedveh_genes_22Nov.vec, decreasing = TRUE)

# 3. Set seed before GSEA
set.seed(20241122)

# 4. Run GSEA
gseGO_oldpwrfrap_v_oldsedveh_22Nov.OUTPUT <- gseGO(
  geneList = oldpwrfrap_v_oldsedveh_genes_22Nov.vec,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  minGSSize = 10,
  maxGSSize = 300,
  eps = 1e-30,
  pvalueCutoff = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)


oldpwrfrap_v_oldsedveh_22Nov.df <- as.data.frame(gseGO_oldpwrfrap_v_oldsedveh_22Nov.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

head(oldpwrfrap_v_oldsedveh_22Nov.df)
dim(oldpwrfrap_v_oldsedveh_22Nov.df)
write_csv(oldpwrfrap_v_oldsedveh_22Nov.df, "oldpwrfrap_v_oldsedveh_22Nov_df.csv")

oldpwrfrap_v_oldsedveh_all_GO_dotplot <- oldpwrfrap_v_oldsedveh_22Nov.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
print(oldpwrfrap_v_oldsedveh_all_GO_dotplot)
# Save (optional)
ggsave("oldpwrfrap_v_oldsedveh_all_GO_dotplot.pdf", oldpwrfrap_v_oldsedveh_all_GO_dotplot, width = 8, height = 8)

# Vector of IDs from each df
pwrfrap_ids <- oldpwrfrap_v_oldsedveh_22Nov.df$ID
oldsed_ids <- oldsed_vs_yngsed_22Nov.df$ID

# Overlapping IDs
pwrfrap_overlap_ids <- intersect(pwrfrap_ids, oldsed_ids)

# Number of overlapping IDs
length(pwrfrap_overlap_ids)
pwrfrap_overlap_ids

oldsed_vs_yngsed_22Nov.df %>%
  filter(ID %in% pwrfrap_overlap_ids)
#==================================================#
#==================================================#

#///////////OldSedIrap vs OldSedVeh//////////////////////////////////#

#Function #1
oldsedirap_v_oldsedveh_genes_22Nov <- read_clean_xlsx_simple(
  xlsx_file = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_IRAP-OLD_SED_VEH.xlsx")

oldsedirap_v_oldsedveh_genes_22Nov %>%
  filter(FDR <0.05)

#Function 2.
create_vec_from_df(oldsedirap_v_oldsedveh_genes_22Nov, "oldsedirap_v_oldsedveh_genes_22Nov")
head(oldsedirap_v_oldsedveh_genes_22Nov.vec)

# 1. Force single-core execution
register(SerialParam())

# 2. Sort the gene list identically each time
oldsedirap_v_oldsedveh_genes_22Nov.vec <- sort(oldsedirap_v_oldsedveh_genes_22Nov.vec, decreasing = TRUE)

# 3. Set seed before GSEA
set.seed(20241122)

# 4. Run GSEA
gseGO_oldsedirap_v_oldsedveh_22Nov.OUTPUT <- gseGO(
  geneList = oldsedirap_v_oldsedveh_genes_22Nov.vec,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  minGSSize = 10,
  maxGSSize = 300,
  eps = 1e-30,
  pvalueCutoff = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)


oldsedirap_v_oldsedveh_22Nov.df <- as.data.frame(gseGO_oldsedirap_v_oldsedveh_22Nov.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

head(oldsedirap_v_oldsedveh_22Nov.df)
dim(oldsedirap_v_oldsedveh_22Nov.df)
write_csv(oldsedirap_v_oldsedveh_22Nov.df, "oldsedirap_v_oldsedveh_22Nov_df.csv")

oldsedirap_v_oldsedveh_all_GO_dotplot <- oldsedirap_v_oldsedveh_22Nov.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
print(oldsedirap_v_oldsedveh_all_GO_dotplot)
# Save (optional)
ggsave("oldsedirap_v_oldsedveh_all_GO_dotplot.pdf", oldsedirap_v_oldsedveh_all_GO_dotplot, width = 8, height = 8)

# Vector of IDs from each df
sedirap_ids <- oldsedirap_v_oldsedveh_22Nov.df$ID
oldsed_ids <- oldsed_vs_yngsed_22Nov.df$ID

# Overlapping IDs
sedirap_overlap_ids <- intersect(sedirap_ids, oldsed_ids)

# Number of overlapping IDs
length(sedirap_overlap_ids)
sedirap_overlap_ids

oldsed_vs_yngsed_22Nov.df %>%
  filter(ID %in% sedirap_overlap_ids)
#==================================================#
#==================================================#

#///////////OldSedFrap vs OldSedVeh//////////////////////////////////#

#Function #1
oldsedfrap_v_oldsedveh_genes_22Nov <- read_clean_xlsx_simple(
  xlsx_file = "20250219_M007853_Set03_edgeRglm_GENE_OLD_SED_FRAP-OLD_SED_VEH.xlsx")

oldsedfrap_v_oldsedveh_genes_22Nov %>%
  filter(FDR <0.05)

#Function 2.
create_vec_from_df(oldsedfrap_v_oldsedveh_genes_22Nov, "oldsedfrap_v_oldsedveh_genes_22Nov")
head(oldsedfrap_v_oldsedveh_genes_22Nov.vec)

# 1. Force single-core execution
register(SerialParam())

# 2. Sort the gene list identically each time
oldsedfrap_v_oldsedveh_genes_22Nov.vec <- sort(oldsedfrap_v_oldsedveh_genes_22Nov.vec, decreasing = TRUE)

# 3. Set seed before GSEA
set.seed(20241122)

# 4. Run GSEA
gseGO_oldsedfrap_v_oldsedveh_22Nov.OUTPUT <- gseGO(
  geneList = oldsedfrap_v_oldsedveh_genes_22Nov.vec,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  minGSSize = 10,
  maxGSSize = 300,
  eps = 1e-30,
  pvalueCutoff = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)


oldsedfrap_v_oldsedveh_22Nov.df <- as.data.frame(gseGO_oldsedfrap_v_oldsedveh_22Nov.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

head(oldsedfrap_v_oldsedveh_22Nov.df)
dim(oldsedfrap_v_oldsedveh_22Nov.df)
write_csv(oldsedfrap_v_oldsedveh_22Nov.df, "oldsedfrap_v_oldsedveh_22Nov_df.csv")

oldsedfrap_v_oldsedveh_all_GO_dotplot <- oldsedfrap_v_oldsedveh_22Nov.df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
print(oldsedfrap_v_oldsedveh_all_GO_dotplot)
# Save (optional)
ggsave("oldsedfrap_v_oldsedveh_all_GO_dotplot.pdf", oldsedfrap_v_oldsedveh_all_GO_dotplot, width = 8, height = 8)

# Vector of IDs from each df
sedfrap_ids <- oldsedfrap_v_oldsedveh_22Nov.df$ID
oldsed_ids <- oldsed_vs_yngsed_22Nov.df$ID

# Overlapping IDs
sedfrap_overlap_ids <- intersect(sedfrap_ids, oldsed_ids)

# Number of overlapping IDs
length(sedfrap_overlap_ids)
sedfrap_overlap_ids

oldsed_vs_yngsed_22Nov.df %>%
  filter(ID %in% sedfrap_overlap_ids)
#==================================================#
#==================================================#

#///////////Join all dataframes//////////////////////////////////#

head(oldsed_vs_yngsed_22Nov.df)
head(oldpwrveh_v_oldsedveh_22Nov.df)
head(oldpwrirap_v_oldsedveh_22Nov.df)
head(oldpwrfrap_v_oldsedveh_22Nov.df)
head(oldsedirap_v_oldsedveh_22Nov.df)
head(oldsedfrap_v_oldsedveh_22Nov.df)

# ============================================================
# 1. Put all GSEA result data frames into a named list
#    (names become prefixes for new columns)
# ============================================================

gsea_list <- list(
  oldpwrveh   = oldpwrveh_v_oldsedveh_22Nov.df,
  oldpwrirap  = oldpwrirap_v_oldsedveh_22Nov.df,
  oldpwrfrap  = oldpwrfrap_v_oldsedveh_22Nov.df,
  oldsedirap  = oldsedirap_v_oldsedveh_22Nov.df,
  oldsedfrap  = oldsedfrap_v_oldsedveh_22Nov.df
)


# ============================================================
# 2. Helper function to rename GSEA columns by prefix
# ============================================================

rename_gsea_cols <- function(df, prefix) {
  df %>%
    dplyr::select(
      ID,
      dplyr::all_of(c("Description", "NES", "p.adjust", "Count"))
    ) %>%
    dplyr::rename(
      !!paste0("NES_", prefix)       := NES,
      !!paste0("p.adjust_", prefix)  := p.adjust,
      !!paste0("Count_", prefix)     := Count
    )
}


# ============================================================
# 3. Merge everything via left_join onto baseline DF
# ============================================================

library(dplyr)
library(purrr)

merged_df <- purrr::reduce(
  .x = names(gsea_list),
  .f = function(out, nm) {
    dplyr::left_join(
      out,
      rename_gsea_cols(gsea_list[[nm]], nm),
      by = "ID"
    )
  },
  .init = oldsed_vs_yngsed_22Nov.df
)


# ============================================================
# 4. Check results
# ============================================================

print(dim(merged_df))
print(head(merged_df))

GO_intervention_merged_df <- merged_df %>%
  dplyr::select(ID, Description.x, NES, p.adjust, Count,
                NES_oldpwrveh, p.adjust_oldpwrveh, Count_oldpwrveh,
                NES_oldpwrirap, p.adjust_oldpwrirap, Count_oldpwrirap,
                NES_oldpwrfrap, p.adjust_oldpwrfrap, Count_oldpwrfrap,
                NES_oldsedirap, p.adjust_oldsedirap, Count_oldsedirap,
                NES_oldsedfrap, p.adjust_oldsedfrap, Count_oldsedfrap)

head(GO_intervention_merged_df)
GO_intervention_merged_df <- GO_intervention_merged_df %>%
  mutate(
    Description.x = ifelse(
      ID == "GO:0002460",
      "adaptive immunity",
      Description.x
    )
  )

head(GO_intervention_merged_df)
readr::write_csv(GO_intervention_merged_df,
                 "GO_intervention_merged_df.csv")
#------------------------------------#

library(dplyr)
library(tidyr)
library(ggplot2)

#--------------------------------------------------------------
# 1. Desired column order
#--------------------------------------------------------------
desired_group_order <- c(
  "OldSed",
  "OldSedIRAP",
  "OldSedFRAP",
  "OldPwrVeh",
  "OldPwrIRAP",
  "OldPwrFRAP"
)

#--------------------------------------------------------------
# 2. Reshape (pivot) to long format
#--------------------------------------------------------------
dot_df <- GO_intervention_merged_df %>%
  dplyr::select(
    Description = Description.x,
    
    NES_OldSed          = NES,
    p.adjust_OldSed     = p.adjust,
    Count_OldSed        = Count,
    
    NES_OldSedIRAP      = NES_oldsedirap,
    p.adjust_OldSedIRAP = p.adjust_oldsedirap,
    Count_OldSedIRAP    = Count_oldsedirap,
    
    NES_OldSedFRAP      = NES_oldsedfrap,
    p.adjust_OldSedFRAP = p.adjust_oldsedfrap,
    Count_OldSedFRAP    = Count_oldsedfrap,
    
    NES_OldPwrVeh       = NES_oldpwrveh,
    p.adjust_OldPwrVeh  = p.adjust_oldpwrveh,
    Count_OldPwrVeh     = Count_oldpwrveh,
    
    NES_OldPwrIRAP      = NES_oldpwrirap,
    p.adjust_OldPwrIRAP = p.adjust_oldpwrirap,
    Count_OldPwrIRAP    = Count_oldpwrirap,
    
    NES_OldPwrFRAP      = NES_oldpwrfrap,
    p.adjust_OldPwrFRAP = p.adjust_oldpwrfrap,
    Count_OldPwrFRAP    = Count_oldpwrfrap
  ) %>%
  pivot_longer(
    cols = -Description,
    names_to = c(".value", "Group"),
    names_pattern = "(NES|p.adjust|Count)_(.*)"
  ) %>%
  
  # Ensure Group ordering
  mutate(Group = factor(Group, levels = desired_group_order))

#--------------------------------------------------------------
# 3. Add NES direction and color intensity
#--------------------------------------------------------------
dot_df <- dot_df %>%
  mutate(
    NES_direction = case_when(
      is.na(NES) ~ "NA",
      NES > 0    ~ "pos",
      NES < 0    ~ "neg"
    ),
    
    # Sig intensity
    color_value = -log10(p.adjust + 1e-300),
    
    # If Count is NA → set to 10 (grey circle)
    Count_plot = ifelse(is.na(Count), 10, Count),
    
    # Override alpha to zero for NA NES (solid grey)
    alpha_plot = ifelse(is.na(NES), 1, scales::rescale(color_value, to = c(0.2,1)))
  )

#--------------------------------------------------------------
# 4. Order rows by original NES (OldSed)
#--------------------------------------------------------------
# Force numeric sort so positive NES go to the top
row_order <- GO_intervention_merged_df %>%
  arrange(desc(NES)) %>%      # highest NES first
  pull(Description.x)

# ggplot has top = last level → reverse here
row_order <- rev(row_order)

# Apply *after* dot_df is fully built
dot_df <- dot_df %>%
  mutate(Description = factor(Description, levels = row_order))

head(dot_df)
#--------------------------------------------------------------
# 5. Plot
#--------------------------------------------------------------
ggplot(dot_df, aes(
  x = Group,
  y = Description,
  size = Count_plot
)) +
  geom_point(aes(
    fill = NES_direction,
    alpha = alpha_plot
  ),
  shape = 21, color = "black") +
  
  scale_fill_manual(
    values = c(
      "pos" = "red",
      "neg" = "blue",
      "NA"  = "grey70"
    ),
    name = "NES direction"
  ) +
  
  scale_size_continuous(
    range = c(1, 8),
    name = "Gene Count"
  ) +
  
  guides(alpha = "none") +  # alpha not shown in legend
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  xlab("Condition") +
  ylab("GO Term")
#===================================#

#---Refine look of ggplot-------#

dot_df <- dot_df %>%
  mutate(
    NES_direction = case_when(
      is.na(NES) ~ "NA",
      NES > 0    ~ "pos",
      TRUE       ~ "neg"
    ),
    
    base_color = case_when(
      NES_direction == "pos" ~ "#FA8072",  # salmon
      NES_direction == "neg" ~ "#2E8BC0",  # sea blue
      TRUE                   ~ "grey80"
    ),
    
    # Larger = more significant = darker
    sig_intensity = ifelse(
      is.na(p.adjust),
      0,  # NA = no darkening
      scales::rescale(-log10(p.adjust + 1e-300), to = c(0, 1))
    ),
    
    # Now apply darkening PROPERLY:
    fill_color = case_when(
      is.na(NES) ~ "grey80",
      TRUE ~ colorspace::darken(base_color, amount = sig_intensity)
    ),
    
    color_color = fill_color,
    
    Count_plot = ifelse(is.na(Count), 10, Count)
  )

GO_intervention_dotplot <- ggplot(dot_df, aes(
  x = Group,
  y = Description,
  size = Count_plot
)) +
  geom_point(
    aes(
      fill  = fill_color,
      color = color_color
    ),
    shape  = 21,
    stroke = 0.3
  ) +
  
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_size_continuous(range = c(2, 5)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  
  xlab("Condition") +
  ylab("GO Term")

ggsave("GO_intervention_dotplot.pdf", GO_intervention_dotplot, width = 8, height = 8)
print(GO_intervention_dotplot)
#=============================================================#
#===========================================================================#

#===========================================================#
#//////Core enrichment genes/////////////////////////////////#
#=========================================================#
library(dplyr)
library(tidyr)
library(stringr)

#-------------------------------------------------------------
# STEP 0 — Select *all* pathways, no NES filtering
#-------------------------------------------------------------
df_all <- oldsed_vs_yngsed_22Nov.df   # keeps all 33 down + all up

#-------------------------------------------------------------
# STEP 1 — Extract + unnest core enrichment Entrez IDs
#-------------------------------------------------------------
gene_long <- df_all %>%
  dplyr::select(ID, core_enrichment) %>%
  mutate(core_enrichment = str_split(core_enrichment, "/")) %>%
  unnest(core_enrichment, keep_empty = FALSE) %>%
  mutate(EntrezID = as.character(core_enrichment)) %>%
  dplyr::select(ID, EntrezID)

#-------------------------------------------------------------
# STEP 2 — Count distinct pathways per Entrez ID
#-------------------------------------------------------------
gene_counts <- gene_long %>%
  distinct(ID, EntrezID) %>%                  # avoid double-counting within a pathway
  dplyr::count(EntrezID, name = "n_pathways") %>%
  arrange(desc(n_pathways))

dim(gene_counts)
head(gene_counts)
dim(gene_counts %>%
      filter(n_pathways > 25))

oldsed_vs_yngsed_GO_core_enrich_genes <- gene_counts
readr::write_csv(oldsed_vs_yngsed_GO_core_enrich_genes,
                 "oldsed_vs_yngsed_GO_core_enrich_genes.csv")

top_GO_oldrapa_genes <- gene_counts %>%
  filter(n_pathways > 20) %>%
  pull(EntrezID)

top_GO_oldrapa_genes
#--------------------------------------------------------#
#========================================================#

#=======================================================================#
#/////////Testing Code with relaxed significance///////////////////#
#----**NOTE** Relaxing significance to 0.25 does not add more pathways to positive NES----#

#------oldpwrveh vs oldsedveh------------#
# 3. Set seed before GSEA
set.seed(20241122)

# 4. Run GSEA
relax_gseGO_oldpwrveh_v_oldsedveh_22Nov.OUTPUT <- gseGO(
  geneList = oldpwrveh_v_oldsedveh_genes_22Nov.vec,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  minGSSize = 10,
  maxGSSize = 300,
  eps = 1e-30,
  pvalueCutoff = 0.25,
  BPPARAM = BiocParallel::SerialParam()
)


relax_oldpwrveh_v_oldsedveh_22Nov.df <- as.data.frame(relax_gseGO_oldpwrveh_v_oldsedveh_22Nov.OUTPUT@result) %>%
  dplyr::mutate(
    Count     = stringr::str_count(core_enrichment, "/") + 1,
    setSize   = as.numeric(setSize),
    geneRatio = Count / setSize
  )

head(relax_oldpwrveh_v_oldsedveh_22Nov.df)
dim(relax_oldpwrveh_v_oldsedveh_22Nov.df)
write_csv(relax_oldpwrveh_v_oldsedveh_22Nov.df, "relax_oldpwrveh_v_oldsedveh_22Nov_df.csv")

relax_oldpwrveh_v_oldsedveh_all_GO_dotplot <- relax_oldpwrveh_v_oldsedveh_22Nov.df %>%
  dplyr::filter(p.adjust < 0.25) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), size = Count, col = p.adjust)) +
  geom_point() +
  theme_bw()
print(relax_oldpwrveh_v_oldsedveh_all_GO_dotplot)

# Vector of IDs from each df
relax_oldpwr_ids <- relax_oldpwrveh_v_oldsedveh_22Nov.df$ID
oldsed_ids <- oldsed_vs_yngsed_22Nov.df$ID

# Overlapping IDs
relax_oldpwr_overlap_ids <- intersect(relax_oldpwr_ids, oldsed_ids)

# Number of overlapping IDs
length(relax_oldpwr_overlap_ids)
relax_oldpwr_overlap_ids

oldsed_vs_yngsed_22Nov.df %>%
  filter(ID %in% relax_oldpwr_overlap_ids)
#------------------------------------------------#
#=======================================================#
#=======================================================#

library(dplyr)
library(stringr)
library(purrr)

#-----------------------------------------------------------
# 0. Helper: deduplicate ENTREZID (keep lowest PValue)
#-----------------------------------------------------------
dedupe_df <- function(df) {
  df %>%
    arrange(PValue) %>%                 # strongest evidence first
    distinct(ENTREZID, .keep_all = TRUE)
}

#-----------------------------------------------------------
# 1. Base comparison (dedupe + rename)
#-----------------------------------------------------------
base_clean <- oldsedveh_v_yngsedveh_genes_22Nov %>%
  dedupe_df() %>%                       # <-- IMPORTANT FIX
  dplyr::select(
    Symbol,
    ENTREZID,
    oldsedveh_logFC = logFC,
    oldsedveh_FDR   = FDR,
    oldsedveh_PValue = PValue
  )

#-----------------------------------------------------------
# 2. List of ALL intervention dfs (deduped automatically)
#-----------------------------------------------------------
dfs <- list(
  oldpwrveh  = oldpwrveh_v_oldsedveh_genes_22Nov,
  oldpwrirap = oldpwrirap_v_oldsedveh_genes_22Nov,
  oldpwrfrap = oldpwrfrap_v_oldsedveh_genes_22Nov,
  oldsedirap = oldsedirap_v_oldsedveh_genes_22Nov,
  oldsedfrap = oldsedfrap_v_oldsedveh_genes_22Nov
) %>%
  # Apply dedupe_df to every element
  purrr::map(dedupe_df)                 # <-- IMPORTANT FIX

#-----------------------------------------------------------
# 3. Helper for renaming each intervention df
#-----------------------------------------------------------
clean_and_prefix <- function(df, prefix) {
  df %>%
    dplyr::select(ENTREZID, logFC, FDR, PValue) %>%
    dplyr::rename(
      !!paste0(prefix, "_logFC")  := logFC,
      !!paste0(prefix, "_FDR")    := FDR,
      !!paste0(prefix, "_PValue") := PValue
    )
}

#-----------------------------------------------------------
# 4. Merge everything onto the base df (one row per ENTREZID)
#-----------------------------------------------------------
oldrapa_gene_merged_df <- base_clean %>%
  purrr::reduce(
    .x = names(dfs),
    .f = function(acc, nm) {
      acc %>%
        left_join(
          clean_and_prefix(dfs[[nm]], nm),
          by = "ENTREZID"
        )
    },
    .init = .
  )

#-----------------------------------------------------------
# 5. Inspect result
#-----------------------------------------------------------
head(oldrapa_gene_merged_df)

#---------------------------------------#

oldrapa_gene_merged_df <- oldrapa_gene_merged_df %>%
  mutate(
    across(
      .cols = everything(),
      .fns = identity
    )
  )

# Identify all prefixes based on existing logFC columns
prefixes <- oldrapa_gene_merged_df %>%
  dplyr::select(ends_with("_logFC")) %>%
  names() %>%
  str_replace("_logFC", "")

# Compute rank scores for each intervention
for (p in prefixes) {
  logfc_col  <- paste0(p, "_logFC")
  pval_col   <- paste0(p, "_PValue")
  rs_col     <- paste0(p, "_rank_score")
  
  oldrapa_gene_merged_df[[rs_col]] <-
    oldrapa_gene_merged_df[[logfc_col]] *
    ((-log10(oldrapa_gene_merged_df[[pval_col]])) ^ 0.25)
}


head(oldrapa_gene_merged_df)
dim(oldrapa_gene_merged_df)

readr::write_csv(oldrapa_gene_merged_df,
                 "oldrapa_gene_merged_df.csv")

head(gene_counts)
dim(gene_counts)

oldsed_core_genes <- gene_counts %>%
  pull(EntrezID) %>%
  unique()

length(oldsed_core_genes)
head(oldsed_core_genes)

oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  filter(oldsedveh_rank_score >0) %>%
  filter(oldsedfrap_rank_score > oldsedveh_rank_score) %>%
  filter(oldsedfrap_FDR < 0.05) %>%
  filter(oldsedfrap_logFC >0) %>%
  filter(oldpwrveh_logFC < 0)


dim(oldrapa_gene_merged_df %>%
      filter(ENTREZID %in% oldsed_core_genes) %>%
      filter(oldsedveh_rank_score <0))

head(oldsedveh_v_yngsedveh_genes_22Nov)

oldsedveh_v_yngsedveh_genes_22Nov %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  arrange(logFC)
#----------------------------------------------#

head(oldsed_vs_yngsed_22Nov.df)
dim(oldsed_vs_yngsed_22Nov.df %>%
      filter(NES < 0))


oldsedveh_v_yngsedveh_genes_22Nov %>%
  filter(Symbol == "Deaf1")
oldpwrveh_v_oldsedveh_genes_22Nov %>%
  filter(Symbol == "Deaf1")
oldpwrveh_v_oldsedveh_genes_22Nov %>%
  filter(Symbol == "Mtor")
oldsedveh_v_yngsedveh_genes_22Nov %>%
  filter(Symbol == "Mtor")

length(oldsed_core_genes)

dim(oldrapa_gene_merged_df %>%
      filter(ENTREZID %in% oldsed_core_genes) %>%
      filter(oldsedveh_rank_score >0))

dim(oldrapa_gene_merged_df %>%
      filter(ENTREZID %in% oldsed_core_genes) %>%
      filter(oldsedveh_rank_score >0) %>%
      filter(oldsedfrap_rank_score > oldsedveh_rank_score))

dim(oldrapa_gene_merged_df %>%
      filter(ENTREZID %in% oldsed_core_genes) %>%
      filter(oldsedveh_rank_score >0) %>%
      filter(oldpwrveh_rank_score < oldsedveh_rank_score))

dwn_w_pwr_genes <- oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  filter(oldsedveh_rank_score >0) %>%
  filter(oldpwrveh_rank_score < oldsedveh_rank_score) %>%
  pull(ENTREZID)

head(oldrapa_gene_merged_df)
head(oldpwrveh_v_oldsedveh_genes_22Nov)
oldpwrveh_v_oldsedveh_genes_22Nov %>%
  filter(ENTREZID %in% dwn_w_pwr_genes) %>%
  arrange(logFC)


oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  filter(oldsedveh_rank_score >0) %>%
  filter(oldsedfrap_rank_score > oldsedveh_rank_score) %>%
  filter(oldsedfrap_FDR < 0.05) %>%
  filter(oldsedfrap_logFC >0) %>%
  pull(ENTREZID)


dim(oldrapa_gene_merged_df %>%
      filter(ENTREZID %in% oldsed_core_genes) %>%
      filter(oldsedveh_rank_score >0) %>%
      filter(oldsedfrap_rank_score < oldsedveh_rank_score) %>%
      filter(oldsedfrap_FDR <0.05))

dim(oldrapa_gene_merged_df %>%
      filter(ENTREZID %in% oldsed_core_genes) %>%
      filter(oldsedveh_rank_score >0) %>%
      filter(oldsedfrap_rank_score > oldsedveh_rank_score) %>%
      filter(oldsedfrap_FDR <0.05))
dim(oldrapa_gene_merged_df)
#--------------------------------------------------#

#==================================================#
#//////Volcano Plots vs OldSedVeh/////////////////#
#==================================================#

#------Genes increased with age-------------------#
pos_age_volcano.df <- oldrapa_gene_merged_df %>%
  filter(ENTREZID %in% oldsed_core_genes) %>%
  filter(oldsedveh_logFC >0)

dim(pos_age_volcano.df)
head(pos_age_volcano.df)

library(ggplot2)
library(dplyr)
#--------------------------------------------#

#/////OldSedFrap Volcano///////////////////////#
#///////Simple volcano plot///////////////////////////////#
oldsedfrap_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldsedfrap_FDR),
    sig = case_when(
      oldsedfrap_FDR < 0.25 & oldsedfrap_logFC > 0  ~ "Up",
      oldsedfrap_FDR < 0.25 & oldsedfrap_logFC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(oldsedfrap_volcano.df, aes(x = oldsedfrap_logFC, y = negLogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: oldsedfrap vs oldsedveh",
    x = "log2 Fold Change (oldsedfrap)",
    y = "-log10(Pvalue)"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6))
#---------------------------------------------------------#
#/////////////////////////////////////////////////////////#

#-------Complex labeled plot------------------------------#
library(ggplot2)
library(ggrepel)
library(dplyr)

highlight_genes <- top_GO_oldrapa_genes   # character vector of ENTREZIDs

oldsedfrap_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldsedfrap_FDR),
    sig = case_when(
      oldsedfrap_FDR < 0.25 & oldsedfrap_logFC > 0  ~ "Up",
      oldsedfrap_FDR < 0.25 & oldsedfrap_logFC < 0  ~ "Down",
      TRUE ~ "NS"
    ),
    highlight = if_else(ENTREZID %in% highlight_genes, "Highlighted", "Other")
  )

ggplot(oldsedfrap_volcano.df, aes(x = oldsedfrap_logFC, y = negLogFDR)) +
  geom_point(
    aes(color = sig),
    alpha = 0.5,
    size = 2
  ) +
  
  # Highlighted genes as larger dark points
  geom_point(
    data = subset(oldsedfrap_volcano.df, highlight == "Highlighted"),
    color = "black",
    size = 3
  ) +
  
  # Optional: label highlighted genes
  geom_text_repel(
    data = subset(oldsedfrap_volcano.df, highlight == "Highlighted"),
    aes(label = Symbol),
    size = 3.5,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey40"
  ) +
  
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  
  labs(
    title = "Volcano Plot: oldsedfrap vs oldsedveh",
    x = "log2 Fold Change (oldsedfrap)",
    y = "-log10(FDR)"
  ) +
  
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6)) +
  theme_minimal(base_size = 14)
#==========================================================#
#==========================================================#

#/////OldSedIrap Volcano///////////////////////#
#///////Simple Volcano Plot////////////////////#
oldsedirap_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldsedirap_FDR),
    sig = case_when(
      oldsedirap_FDR < 0.25 & oldsedirap_logFC > 0  ~ "Up",
      oldsedirap_FDR < 0.25 & oldsedirap_logFC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(oldsedirap_volcano.df, aes(x = oldsedirap_logFC, y = negLogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: oldsedirap vs oldsedveh",
    x = "log2 Fold Change (oldsedirap)",
    y = "-log10(Pvalue)"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6))
#---------------------------------------------------------#
#---------------------------------------------------------#

oldsedirap_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldsedirap_FDR),
    sig = case_when(
      oldsedirap_FDR < 0.25 & oldsedirap_logFC > 0  ~ "Up",
      oldsedirap_FDR < 0.25 & oldsedirap_logFC < 0  ~ "Down",
      TRUE ~ "NS"
    ),
    highlight = if_else(ENTREZID %in% highlight_genes, "Highlighted", "Other")
  )

ggplot(oldsedirap_volcano.df, aes(x = oldsedirap_logFC, y = negLogFDR)) +
  geom_point(
    aes(color = sig),
    alpha = 0.5,
    size = 2
  ) +
  
  # Highlighted genes as larger dark points
  geom_point(
    data = subset(oldsedirap_volcano.df, highlight == "Highlighted"),
    color = "black",
    size = 3
  ) +
  
  # Optional: label highlighted genes
  geom_text_repel(
    data = subset(oldsedirap_volcano.df, highlight == "Highlighted"),
    aes(label = Symbol),
    size = 3.5,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey40"
  ) +
  
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  
  labs(
    title = "Volcano Plot: oldsedirap vs oldsedveh",
    x = "log2 Fold Change (oldsedirap)",
    y = "-log10(FDR)"
  ) +
  
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6)) +
  theme_minimal(base_size = 14)
#============================================================#
#============================================================#


#/////OldPwrVeh Volcano///////////////////////#
#////////Simple Volcano Plot/////////////////#
head(pos_age_volcano.df)
dim(pos_age_volcano.df)
dim(pos_age_volcano.df %>%
      filter(oldpwrveh_rank_score <0))

oldpwrveh_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldpwrveh_FDR),
    sig = case_when(
      oldpwrveh_FDR < 0.25 & oldpwrveh_logFC > 0  ~ "Up",
      oldpwrveh_FDR < 0.25 & oldpwrveh_logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(oldpwrveh_volcano.df, aes(x = oldpwrveh_logFC, y = negLogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: oldpwrveh vs oldsedveh",
    x = "log2 Fold Change (oldpwrveh)",
    y = "-log10(Pvalue)"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6))
#---------------------------------------------------------#
#---------------------------------------------------------#

oldpwrveh_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldpwrveh_FDR),
    sig = case_when(
      oldpwrveh_FDR < 0.25 & oldpwrveh_logFC > 0  ~ "Up",
      oldpwrveh_FDR < 0.25 & oldpwrveh_logFC < 0  ~ "Down",
      TRUE ~ "NS"
    ),
    highlight = if_else(ENTREZID %in% highlight_genes, "Highlighted", "Other")
  )

ggplot(oldpwrveh_volcano.df, aes(x = oldpwrveh_logFC, y = negLogFDR)) +
  geom_point(
    aes(color = sig),
    alpha = 0.5,
    size = 2
  ) +
  
  # Highlighted genes as larger dark points
  geom_point(
    data = subset(oldpwrveh_volcano.df, highlight == "Highlighted"),
    color = "black",
    size = 3
  ) +
  
  # Optional: label highlighted genes
  geom_text_repel(
    data = subset(oldpwrveh_volcano.df, highlight == "Highlighted"),
    aes(label = Symbol),
    size = 3.5,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey40"
  ) +
  
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  
  labs(
    title = "Volcano Plot: oldpwrveh vs oldsedveh",
    x = "log2 Fold Change (oldpwrveh)",
    y = "-log10(FDR)"
  ) +
  
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6)) +
  theme_minimal(base_size = 14)
#=====================================================#
#=====================================================#


#/////OldPwrFrap Volcano///////////////////////#

oldpwrfrap_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldpwrfrap_FDR),
    sig = case_when(
      oldpwrfrap_FDR < 0.25 & oldpwrfrap_logFC > 0  ~ "Up",
      oldpwrfrap_FDR < 0.25 & oldpwrfrap_logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(oldpwrfrap_volcano.df, aes(x = oldpwrfrap_logFC, y = negLogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: oldpwrfrap vs oldsedveh",
    x = "log2 Fold Change (oldpwrfrap)",
    y = "-log10(Pvalue)"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6))
#---------------------------------------------------------#

#/////OldPwrIrap Volcano///////////////////////#

oldpwrirap_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogFDR = -log10(oldpwrirap_FDR),
    sig = case_when(
      oldpwrirap_FDR < 0.25 & oldpwrirap_logFC > 0  ~ "Up",
      oldpwrirap_FDR < 0.25 & oldpwrirap_logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(oldpwrirap_volcano.df, aes(x = oldpwrirap_logFC, y = negLogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: oldpwrirap vs oldsedveh",
    x = "log2 Fold Change (oldpwrirap)",
    y = "-log10(Pvalue)"
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 6))
#----------------------------------------------------#







#/////////Testing effect of Pwr///////////////////////////#

pos_age_volcano.df %>% 
  filter(oldsedveh_FDR < 0.05) %>%
  filter(oldpwrveh_logFC < 0) %>%
  arrange(oldpwrveh_PValue)

head(oldrapa_gene_merged_df)
oldrapa_gene_merged_df %>%
  filter(oldsedveh_logFC > 0) %>%
  filter(oldsedveh_FDR < 0.05) %>%
  filter(oldpwrveh_logFC < 0) %>%
  filter(oldpwrveh_FDR < 0.05)
#---------------------------------------#

#/////OldPwrFRAP Volcano///////////////////////#

oldpwrfrap_volcano.df <- pos_age_volcano.df %>%
  mutate(
    negLogP = -log10(oldpwrfrap_PValue),
    sig = case_when(
      oldpwrfrap_PValue < 0.1 & oldpwrfrap_logFC > 0  ~ "Up",
      oldpwrfrap_PValue < 0.1 & oldpwrfrap_logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(oldpwrfrap_volcano.df, aes(x = oldpwrfrap_logFC, y = negLogP)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "firebrick2", "Down" = "dodgerblue3", "NS" = "grey60")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: oldpwrfrap vs oldsedveh",
    x = "log2 Fold Change (vs oldsedveh)",
    y = "-log10(Pvalue)"
  ) +
  theme_minimal(base_size = 14)
#---------------------------------------------------------#