#aging lipidome analysis

# Reads the first sheet by default
aging_lipidome <- read_excel("aging_lipidome.xlsx")

aging_lipidome
aging_lipidome$Metabolite <- sub("\\|.*", "", aging_lipidome$Metabolite)

aging_lipidome$Metabolite <- sub(";.*", "", aging_lipidome$Metabolite)
aging_lipidome

aging_lipidome$Metabolite <- sub("/0:0.*", "", aging_lipidome$Metabolite)
aging_lipidome

formatted_lipid_age_sig <- read_excel("new_lipid_age_signature_formatted.xlsx")
formatted_lipid_age_sig.vec <- formatted_lipid_age_sig %>%
  pull(Lipid)

formatted_lipid_age_sig.vec
aging_lipidome %>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec)

colnames(aging_lipidome)


# Select only columns that end with "F_GF"
female_aging_lipidome <- aging_lipidome %>%
  dplyr::select(Metabolite, contains("F_GF"))


female_aging_lipidome%>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec)

female_aging_lipidome

formatted_lipid_age_sig.vec
#-----------------------------------------------#

library(purrr)
library(tibble)

# Step 1: Define column names for the groups
young_grp <- c("2M_F_GF...9", "2M_F_GF...34", "2M_F_GF...36")
old_grp   <- c("24M_F_GF...21", "24M_F_GF...27", "24M_F_GF...42")

# Step 2: Identify lipid columns (i.e., all group columns)
lipid_cols <- c(young_grp, old_grp)

# Step 3: Impute zero values with half of the minimum non-zero per column
min_nonzero <- sapply(female_aging_lipidome[lipid_cols], function(x) {
  min(x[x > 0], na.rm = TRUE)
})

lipid_imputed <- female_aging_lipidome
for (lip in lipid_cols) {
  x <- lipid_imputed[[lip]]
  x[x == 0] <- min_nonzero[[lip]] / 2
  lipid_imputed[[lip]] <- x
}

# Step 4: Log2-transform the imputed values
lipid_log <- lipid_imputed %>%
  mutate(across(all_of(lipid_cols), log2))

lipid_log

# Step 5: Run t-tests per lipid species
aging_lipidome_stats <- purrr::map_dfr(1:nrow(lipid_log), function(i) {
  lipid_name <- lipid_log$Metabolite[i]
  vals_young <- as.numeric(lipid_log[i, young_grp])
  vals_old   <- as.numeric(lipid_log[i, old_grp])
  
  # Skip rows with zero variance in either group
  if (sd(vals_young) == 0 || sd(vals_old) == 0) {
    return(tibble(
      Metabolite = lipid_name,
      Mean_Young = mean(vals_young),
      Mean_Old   = mean(vals_old),
      Log2_FC_Old_vs_Young = mean(vals_old) - mean(vals_young),
      P_value = NA_real_
    ))
  }
  
  # Otherwise, perform t-test
  ttest <- t.test(vals_old, vals_young)
  
  tibble(
    Metabolite = lipid_name,
    Mean_Young = mean(vals_young),
    Mean_Old   = mean(vals_old),
    Log2_FC_Old_vs_Young = mean(vals_old) - mean(vals_young),
    P_value = ttest$p.value
  )
}) %>%
  mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  arrange(FDR)

# Optional: View top hits
head(aging_lipidome_stats)

aging_lipidome_stats %>%
  filter(FDR <0.05)

print(aging_lipidome_stats %>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec), n=21)

aging_lipidome_stats_unique <- aging_lipidome_stats %>%
  distinct(Metabolite, .keep_all = TRUE)

dim(aging_lipidome_stats_unique)

print(aging_lipidome_stats_unique %>%
        filter(Metabolite %in% formatted_lipid_age_sig.vec), n=21)

ovrlp_aging_sig_lipidome <- aging_lipidome_stats_unique %>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec) %>%
  arrange(desc(Log2_FC_Old_vs_Young))

write.csv(ovrlp_aging_sig_lipidome, "ovrlp_aging_sig_lipidome.csv", row.names = FALSE)

ovrlp_aging_lipid_vec <- aging_lipidome_stats_unique %>%
  filter(Metabolite %in% formatted_lipid_age_sig.vec) %>%
  arrange(desc(Log2_FC_Old_vs_Young)) %>%
  pull(Metabolite)

ovrlp_aging_lipid_vec
formatted_lipid_age_sig.vec
#---------------------------------------------#

#============================================#
#-----lgsea from lipidome--------------------#
#============================================#

# Ensure ranking vector is built from lipids with proper names and values
rank_df <- aging_lipidome_stats_unique %>%
  filter(!is.na(Log2_FC_Old_vs_Young), !is.na(P_value)) %>%
  mutate(
    p_safe = pmax(P_value, 1e-300),
    stat = Log2_FC_Old_vs_Young
  ) %>%
  group_by(Metabolite) %>%
  summarise(stat = mean(stat), .groups = "drop") %>%
  arrange(desc(stat))

# Convert to named vector
lipid_ranks <- setNames(rank_df$stat, rank_df$Metabolite)

# Only use the lipids that are present in the ranked list
lipid_signature_filtered <- formatted_lipid_age_sig.vec[
  formatted_lipid_age_sig.vec %in% names(lipid_ranks)
]

lipid_sets <- list(AgingSignature = lipid_signature_filtered)

fgsea_res <- fgsea(
  pathways = lipid_sets,
  stats = lipid_ranks,
  nperm = 1000
)

print(fgsea_res)

# Plot enrichment
plotEnrichment(lipid_sets$AgingSignature, lipid_ranks) +
  ggtitle("Effect of Age on Aging Lipid Signature")

female_lipidome_signature <- aging_lipidome_stats_unique %>%
  filter(FDR <0.05) %>%
  pull(Metabolite)

female_lipidome_signature

aging_lipidome_stats_unique %>%
  filter (Log2_FC_Old_vs_Young >0.49)

aging_lipidome_stats_unique %>%
  filter(FDR >0.05 & Log2_FC_Old_vs_Young >0)

dim(aging_lipidome_stats_unique)
#--------------------------------------------#

# Subset dataframe
lipid_log_subset <- lipid_log %>%
  dplyr::select(Metabolite, any_of(c(young_grp, old_grp)))

lipid_log_subset %>%
  filter(Metabolite %in% ovrlp_aging_lipid_vec) %>%
  distinct(Metabolite, .keep_all = TRUE)

# Step 1: Compute the young mean
lipid_log_wide_logFC <- lipid_log_subset %>%
  rowwise() %>%
  mutate(young_mean = mean(c_across(all_of(young_grp)), na.rm = TRUE)) %>%
  ungroup()

# Step 2: Add logFC columns for each 24M sample
for (sample in old_grp) {
  new_col <- paste0("logFC_", sample)
  lipid_log_wide_logFC[[new_col]] <- lipid_log_wide_logFC[[sample]] - lipid_log_wide_logFC$young_mean
}

# Optional: Arrange columns for readability
lipid_log_wide_logFC <- lipid_log_wide_logFC %>%
  dplyr::select(Metabolite, all_of(young_grp), young_mean, all_of(old_grp), starts_with("logFC_"))

lipid_log_wide_logFC %>%
  filter(Metabolite %in% ovrlp_aging_lipid_vec) %>%
  distinct(Metabolite, .keep_all = TRUE) %>%
  dplyr::select(Metabolite, starts_with("logFC_"))

lipid_log_wide_logFC %>%
  filter(Metabolite == "PC 38:6") %>%
  dplyr::select(Metabolite, starts_with("logFC_"))

