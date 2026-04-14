#Tuning DIABLO number of features and correlations.
#This code built off of DIABLO_again.R

#use the following matrices for the oldsedveh vs yngsedveh comparison

#logCPM_symbol_yngOLDsed.RDS
#metabo_log_yngOLDsed.RDS
#lipid_log_yngOLDsed.RDS

run_diablo_once <- function(
    RNA_train, Metabo_train, Lipid_train,   # training-only matrices
    RNA_full, Metabo_full, Lipid_full,      # full matrices (all samples)
    yng_sed_samples, old_sedveh_samples,
    old_pwrveh_samples, old_pwrirap_samples, old_pwrfrap_samples,
    RNA_keep, Met_keep, Lip_keep,
    design_strength = 1.0
){
  library(mixOmics)
  library(dplyr)
  #
  #-------------------------------------------#
  # Step 1: Build training objects
  #-------------------------------------------#
  X_train <- list(
    RNA        = RNA_train,
    Metabolite = Metabo_train,
    Lipid      = Lipid_train
  )
  
  Y_train <- factor(c(
    rep("YNG", length(yng_sed_samples)),
    rep("OLD", length(old_sedveh_samples))
  ), levels = c("YNG","OLD"))
  
  design <- matrix(design_strength, 3, 3,
                   dimnames = list(names(X_train), names(X_train)))
  diag(design) <- 0
  
  #-------------------------------------------#
  # Step 2: Train DIABLO on YNG vs OLD
  #-------------------------------------------#
  diablo_fit <- block.splsda(
    X_train, Y_train, ncomp = 1,
    keepX = list(RNA = RNA_keep,
                 Metabolite = Met_keep,
                 Lipid = Lip_keep),
    design = design
  )
  
  #-------------------------------------------#
  # (Optional) Step 3: Get selected features per block
  # (not needed for prediction, but useful to inspect)
  #-------------------------------------------#
  rna_feats <- selectVar(diablo_fit, block = "RNA",        comp = 1)$RNA$name
  met_feats <- selectVar(diablo_fit, block = "Metabolite", comp = 1)$Metabolite$name
  lip_feats <- selectVar(diablo_fit, block = "Lipid",      comp = 1)$Lipid$name
  # You could save these somewhere if you want, but we don't subset on them.
  
  #-------------------------------------------#
  # Step 4: Use FULL matrices for projection
  #         (features already harmonized to match training)
  #-------------------------------------------#
  X_all <- list(
    RNA        = RNA_full,
    Metabolite = Metabo_full,
    Lipid      = Lipid_full
  )
  
  #-------------------------------------------#
  # Step 5: Project all samples into DIABLO space
  #-------------------------------------------#
  proj <- predict(diablo_fit, newdata = X_all)
  scores <- lapply(proj$variates, function(x) x[, 1])
  scores_mat <- do.call(cbind, scores)
  comp1 <- rowMeans(scores_mat)
  # Flip direction so higher Comp1 = older samples
  comp1 <- -comp1
  
  #-------------------------------------------#
  # Step 6: Build output tibble
  #-------------------------------------------#
  tibble(
    Sample = rownames(scores_mat),
    Comp1  = comp1,
    RNA_keep = RNA_keep,
    Met_keep = Met_keep,
    Lip_keep = Lip_keep,
    DesignStrength = design_strength,
    Group = case_when(
      Sample %in% yng_sed_samples     ~ "Young_Sed",
      Sample %in% old_sedveh_samples  ~ "Old_SedVeh",
      Sample %in% old_pwrveh_samples  ~ "Old_PwrVeh",
      Sample %in% old_pwrirap_samples ~ "Old_PwrIRap",
      Sample %in% old_pwrfrap_samples ~ "Old_PwrFRap",
      TRUE                            ~ "Other"
    )
  )
}
#----------------------------------------------------------------------#

summarize_groups <- function(res_df) {
  res_df %>%
    group_by(RNA_keep, Met_keep, Lip_keep, DesignStrength, Group) %>%
    summarise(
      n = n(),
      mean_Comp1 = mean(Comp1),
      sd_Comp1   = sd(Comp1),
      se_Comp1   = sd_Comp1 / sqrt(n),
      .groups = "drop"
    )
}
#----------------------------------------------------------------------#

run_diablo_grid <- function(
    RNA_train, Metabo_train, Lipid_train,
    RNA_full, Metabo_full, Lipid_full,
    yng_sed_samples, old_sedveh_samples,
    old_pwrveh_samples, old_pwrirap_samples, old_pwrfrap_samples,
    RNA_vec  = c(30, 50, 100),
    Met_vec  = c(10, 20, 30),
    Lip_vec  = c(10, 20, 30),
    DS_vec   = c(1.0, 0.8)
){
  library(purrr)
  library(dplyr)
  
  grid <- expand.grid(
    RNA_keep = RNA_vec,
    Met_keep = Met_vec,
    Lip_keep = Lip_vec,
    DS       = DS_vec,
    stringsAsFactors = FALSE
  )
  
  message("Running DIABLO for ", nrow(grid), " parameter combinations...")
  
  all_scores <- pmap_dfr(grid, function(RNA_keep, Met_keep, Lip_keep, DS){
    run_diablo_once(
      RNA_train, Metabo_train, Lipid_train,
      RNA_full, Metabo_full, Lipid_full,
      yng_sed_samples, old_sedveh_samples,
      old_pwrveh_samples, old_pwrirap_samples, old_pwrfrap_samples,
      RNA_keep, Met_keep, Lip_keep,
      design_strength = DS
    )
  })
  
  group_stats <- summarize_groups(all_scores)
  
  list(
    sample_scores = all_scores,
    group_summary = group_stats
  )
}
#--------------------------------------------------------------------#

grid_results <- run_diablo_grid(
  RNA_train, Metabo_train, Lipid_train,
  RNA_full,  Metabo_full,  Lipid_full,
  yng_sed_samples, old_sedveh_samples,
  old_pwrveh_samples, old_pwrirap_samples, old_pwrfrap_samples,
  RNA_vec = c(25,50,100,150,200),
  Met_vec = c(10,20,30,40,50),
  Lip_vec = c(10,20,30,40,50),
  DS_vec  = c(1.0, 0.8)
)
#--------------------------------------------------------------------#

group_summary_wide <- grid_results$group_summary %>%
  dplyr::select(RNA_keep, Met_keep, Lip_keep, DesignStrength, Group, mean_Comp1) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = mean_Comp1
  ) %>%
  arrange(RNA_keep, Met_keep, Lip_keep, DesignStrength)

harmonize_group_summary <- function(df_wide) {
  
  # Columns containing your group means:
  group_cols <- c(
    "Young_Sed", "Old_SedVeh",
    "Old_PwrVeh", "Old_PwrIRap", "Old_PwrFRap"
  )
  
  df_wide %>%
    rowwise() %>%
    mutate(
      # if Old < Young, flip the sign
      flip = ifelse(Old_SedVeh < Young_Sed, -1, 1),
      across(all_of(group_cols), ~ .x * flip)
    ) %>%
    ungroup() %>%
    dplyr::select(-flip)
}

group_summary_wide_harmonized <-
  harmonize_group_summary(group_summary_wide)

group_summary_wide_harmonized %>%
  filter(RNA_keep == 100)
#-----------------------------------------------------------#

#======================================================================#
#  ADD BIOLOGICAL METRICS + ORDERING CHECK + QUALITY SCORE
#======================================================================#

# Function to check full ordering:
# Young_Sed < Old_PwrVeh < Old_PwrIRap < Old_PwrFRap < Old_SedVeh
check_ordering <- function(df_row) {
  
  vals <- c(
    df_row$Young_Sed,
    df_row$Old_PwrVeh,
    df_row$Old_PwrIRap,
    df_row$Old_PwrFRap,
    df_row$Old_SedVeh
  )
  
  all(diff(vals) > 0)
}

head(grid_results)
head(group_summary_wide_harmonized)

metric_table <- group_summary_wide_harmonized %>%
  
  # 1) Add ordering consistency  
  rowwise() %>%
  mutate(order_ok = check_ordering(cur_data())) %>%
  ungroup() %>%
  
  # 2) Age separation (Old baseline minus Young)
  mutate(age_separation = Old_SedVeh - Young_Sed) %>%
  
  # 3) Restoration score (how much PwrVeh moves toward Young)
  mutate(restoration_score = Old_SedVeh - Old_PwrVeh) %>%
  
  # 4) Interference metrics (how much rapamycin blunts the benefit)
  mutate(
    IRap_interference = Old_PwrIRap - Old_PwrVeh,
    FRap_interference = Old_PwrFRap - Old_PwrVeh
  ) %>%
  
  # 5) Composite quality metric  
  mutate(
    score_quality =
      1*order_ok +
      scales::rescale(age_separation, to = c(0, 1)) +
      scales::rescale(restoration_score, to = c(0, 1)) +
      scales::rescale(IRap_interference, to = c(0, 1)) +
      scales::rescale(FRap_interference, to = c(0, 1))
  ) %>%
  
  arrange(desc(score_quality))  # best parameter sets at the top


# Show the full table
metric_table

# If you want to look at the best rows for sanity:
metric_table %>% slice_max(score_quality, n = 10)

metric_table %>%
  filter(RNA_keep == 30)

readr::write_csv(metric_table, "metric_table.csv")
#--------------------------------------------------------#

#==============================================================
# Extract the sample-level data for the chosen hyperparameters
#==============================================================

# Extract the best model results
best_df <- grid_results$sample_scores %>%
  filter(
    RNA_keep == 200,
    Met_keep == 50,
    Lip_keep == 50,
    DesignStrength == 0.8
  )

# Determine expected direction:
# If Old_SedVeh is LOWER than Young_Sed, flip the sign
young_mean <- mean(best_df$Comp1[best_df$Group == "Young_Sed"])
oldsed_mean <- mean(best_df$Comp1[best_df$Group == "Old_SedVeh"])

flip_factor <- ifelse(oldsed_mean < young_mean, -1, 1)

# Apply sign harmonization
best_df <- best_df %>%
  mutate(Comp1 = Comp1 * flip_factor)

# Inspect again
best_df
table(best_df$Group)


#==============================================================
#               RUN ANOVA
#==============================================================

anova_best <- aov(Comp1 ~ Group, data = best_df)
summary(anova_best)

#==============================================================
#               TUKEY POST-HOC TEST
#==============================================================

tukey_best <- TukeyHSD(anova_best)
tukey_best

#==============================================================
# Optional — Compact Letter Display (which groups differ?)
#==============================================================
if (!require(multcompView)) install.packages("multcompView")
library(multcompView)

tukey_cld <- multcompLetters4(anova_best, tukey_best)
tukey_cld
#------------------------------------------------------#

library(ggplot2)
library(dplyr)

#----------------------------------------------------------
# Order groups biologically
#----------------------------------------------------------
best_df$Group <- factor(
  best_df$Group,
  levels = c("Young_Sed",
             "Old_SedVeh",
             "Old_PwrVeh",
             "Old_PwrIRap",
             "Old_PwrFRap"
             )
)

#----------------------------------------------------------
# A nice custom color palette (optional)
#----------------------------------------------------------
group_colors <- c(
  Young_Sed   = "#1F78B4",  # blue
  Old_SedVeh  = "#E31A1C",  # red
  Old_PwrVeh  = "#33A02C",  # green
  Old_PwrIRap = "#FB9A99",  # light red
  Old_PwrFRap = "#FDBF6F"   # orange
)

#----------------------------------------------------------
# Plot
#----------------------------------------------------------
p <- ggplot(best_df, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = "mean", geom = "point",
               shape = 23, size = 3, fill = "white", color = "black") +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "DIABLO Aging Axis (RNA=100, Met=30, Lip=30, DS=0.8)",
    x = "",
    y = "Comp1 Aging Axis Score"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 18)
  )

p





#-------------------------------------------------------#
#-----summary statistics-------------------------------#

grid_results$sample_scores
grid_results$group_summary
grid_results$group_summary %>%
  group_by(RNA_keep, Met_keep, Lip_keep, DesignStrength) %>%
  arrange(mean_Comp1)


#===========================================================================#
#////Single Test Worked nearly perfectly as it did in the original analysis//#
#============================================================================#
test_out <- run_diablo_once(
  RNA_train, Metabo_train, Lipid_train,
  RNA_full, Metabo_full, Lipid_full,
  yng_sed_samples, old_sedveh_samples,
  old_pwrveh_samples, old_pwrirap_samples, old_pwrfrap_samples,
  RNA_keep = 50,
  Met_keep = 20,
  Lip_keep = 20,
  design_strength = 1.0
)
test_out

group_summary <- test_out %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_Comp1 = mean(Comp1, na.rm = TRUE),
    sd_Comp1   = sd(Comp1, na.rm = TRUE),
    se_Comp1   = sd_Comp1 / sqrt(n)
  )

group_summary

anova_result <- aov(Comp1 ~ Group, data = test_out)
summary(anova_result)

tukey <- TukeyHSD(anova_result)
tukey
#==========================================================#

library(dplyr)
library(purrr)
library(tidyr)
library(broom)

#========================================================================#
# FUNCTION: Run ANOVA + Tukey for one hyperparameter combination
#========================================================================#

run_anova_tukey <- function(df) {
  
  # run ANOVA
  fit <- aov(Comp1 ~ Group, data = df)
  
  # Tukey
  tuk <- TukeyHSD(fit)
  
  # Convert Tukey table into tidy form
  tuk_df <- as.data.frame(tuk$Group)
  tuk_df$Comparison <- rownames(tuk_df)
  rownames(tuk_df) <- NULL
  
  # keep only adj p-values, with readable column name
  tuk_df <- tuk_df %>%
    dplyr::select(Comparison, p.adj = `p adj`)
  
  return(tuk_df)
}

#========================================================================#
# FUNCTION: Build tukey results for every hyperparameter set
#========================================================================#

tukey_table <- grid_results$sample_scores %>%
  
  # group by hyperparameter combination
  group_by(RNA_keep, Met_keep, Lip_keep, DesignStrength) %>%
  
  # nest each subset
  nest() %>%
  
  # run ANOVA + Tukey inside each
  mutate(
    tukey_results = map(data, run_anova_tukey)
  ) %>%
  
  # unnest so each comparison is one row
  dplyr::select(-data) %>%
  unnest(tukey_results)


#========================================================================#
# OPTIONAL: Spread comparisons wide to one row per hyperparameter set
#========================================================================#

tukey_wide <- tukey_table %>%
  pivot_wider(
    names_from = Comparison,
    values_from = p.adj,
    names_prefix = "p_"
  )

#========================================================================#
# MERGE INTO METRIC TABLE
#========================================================================#

metric_table_with_tukey <- metric_table %>%
  left_join(
    tukey_wide,
    by = c("RNA_keep", "Met_keep", "Lip_keep", "DesignStrength")
  )

# Show the enhanced table
metric_table_with_tukey

metric_table_with_tukey %>%
  filter(`p_Old_PwrVeh-Old_PwrFRap` < 0.05)

readr::write_csv(metric_table_with_tukey, "metric_table_with_tukey.csv")

metric_table_with_tukey %>%
  arrange(desc(age_separation))
#========================================================#
#========================================================#

#////////average all feature/strength combinations///////#

#============================================================
# BUILD LONG DATASET WITH 54 VALUES PER GROUP
#============================================================

model_values_long <- metric_table %>%
  dplyr::select(
    RNA_keep, Met_keep, Lip_keep, DesignStrength,
    Young_Sed, Old_SedVeh, Old_PwrVeh, Old_PwrIRap, Old_PwrFRap
  ) %>%
  tidyr::pivot_longer(
    cols = c(Young_Sed, Old_SedVeh, Old_PwrVeh, Old_PwrIRap, Old_PwrFRap),
    names_to = "Group",
    values_to = "Comp1"
  ) %>%
  mutate(Group = factor(
    Group,
    levels = c("Young_Sed", "Old_SedVeh", "Old_PwrVeh", "Old_PwrIRap", "Old_PwrFRap")
  ))

# Confirm 54 points per group:
table(model_values_long$Group)


#============================================================
# 1) BOXPLOT (54 points per intervention)
#============================================================

library(ggplot2)

ggplot(model_values_long, aes(x = Group, y = Comp1, fill = Group)) +
  geom_boxplot(width = 0.7, alpha = 0.8) +
  labs(
    title = "DIABLO Model-Derived Comp1 Scores (54 models per intervention)",
    x = "",
    y = "Comp1 Score"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")


library(ggplot2)

ggplot(model_values_long, aes(x = Group, y = Comp1, fill = Group)) +
  
  # Violin for full distribution
  geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
  
  # Boxplot inside the violin
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, color = "black") +
  
  # Jittered individual points (optional)
  geom_jitter(width = 0.1, alpha = 0.35, size = 1) +
  
  scale_fill_brewer(palette = "Set2") +
  
  labs(
    title = "Distribution of DIABLO Comp1 Scores Across 54 Model Fits",
    subtitle = "Each group shows 54 model-derived values (one per hyperparameter set)",
    x = "",
    y = "Comp1 Score"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

