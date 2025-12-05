#Multi-omic integration...again

triceps_physio <- read.csv("Triceps_Physio_Data_v2.csv")
fraility_data <- read.csv("Fraility_Data.csv")
triceps_physio
proj_df

library(dplyr)
library(tidyr)
library(ggplot2)

#--- Step 1. Reshape the physiology data so samples are rows ---
physio_long <- as.data.frame(t(triceps_physio[-1]))   # transpose, dropping the "Tube.code" label row
colnames(physio_long) <- triceps_physio$Tube.code     # use row 1 as column names
physio_long$Sample <- rownames(physio_long)           # move row names to column
rownames(physio_long) <- NULL

# convert numeric columns to numeric
physio_long <- physio_long %>%
  mutate(across(-Sample, as.numeric))
physio_long

#--- Step 2. Merge with DIABLO Comp1 data ---
merged <- proj_df %>%
  left_join(physio_long, by = "Sample")

# check that sample overlap is correct
cat("Samples in DIABLO but not physio:\n")
print(setdiff(proj_df$Sample, physio_long$Sample))
cat("Samples in physio but not DIABLO:\n")
print(setdiff(physio_long$Sample, proj_df$Sample))

# keep only samples present in both
merged <- merged %>%
  filter(Sample %in% physio_long$Sample)

#--- Step 3. Compute Spearman correlations between Comp1 and all physio measures ---
# identify numeric physio variables (ignore Comp1)
physio_vars <- merged %>%
  dplyr::select(where(is.numeric), -Comp1)

cor_results <- lapply(names(physio_vars), function(v) {
  x <- merged$Comp1
  y <- merged[[v]]
  ok <- complete.cases(x, y)  # remove NAs pairwise
  if (sum(ok) >= 4) {         # require at least 4 samples for correlation
    test <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
    data.frame(
      Variable = v,
      r = test$estimate,
      p = test$p.value,
      n = sum(ok)
    )
  } else {
    data.frame(Variable = v, r = NA, p = NA, n = sum(ok))
  }
}) %>%
  bind_rows() %>%
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  arrange(FDR)

#--- Step 4. Print tidy results ---
print(cor_results, digits = 3)

#--- Step 5. Plot results ---
ggplot(cor_results, aes(x = reorder(Variable, r), y = r,
                        fill = p < 0.05)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70")) +
  labs(x = "", y = "Spearman r (DIABLO Comp1 vs Physiology)",
       title = "Physiologic correlates of molecular aging axis",
       fill = "FDR < 0.05") +
  geom_text(aes(label = sprintf("n=%d", n)), 
            hjust = -0.2, size = 3, color = "black")
#-------------------------------------------------------------#
merged

old_only <- merged %>% dplyr::filter(Group != "Young_Sed")

# ANOVA across old groups
anova_old <- aov(Comp1 ~ Group, data = old_only); summary(anova_old)

# Spearman with key metabolic traits within old only
key_vars <- c("Post Fasting Blood Glucose", "Post GTT (AOC)", "Post Insulin Sensitivity (AUC)")
old_cor <- lapply(key_vars, function(v){
  x <- old_only$Comp1; y <- old_only[[v]]
  ok <- complete.cases(x,y)
  if(sum(ok)>=4){
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    data.frame(var=v, r=ct$estimate, p=ct$p.value, n=sum(ok))
  } else data.frame(var=v, r=NA, p=NA, n=sum(ok))
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(FDR = p.adjust(p,"fdr"))
old_cor

ggplot(old_only, aes(x = Comp1, y = `Post GTT (AOC)`, color = Group)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=FALSE) +
  theme_minimal(base_size=14) +
  labs(x = "DIABLO Comp1 (higher = molecularly older)",
       y = "Glucose Tolerance Test AUC",
       title = "Molecular aging axis correlates with impaired glucose metabolism in old mice")

#---------------------------------------------------------------------------------------------#

group_means <- old_only %>%
  group_by(Group) %>%
  summarize(mean_Comp1 = mean(Comp1),
            se_Comp1 = sd(Comp1)/sqrt(n()),
            mean_GTT = mean(`Post GTT (AOC)`, na.rm=TRUE),
            se_GTT = sd(`Post GTT (AOC)`, na.rm=TRUE)/sqrt(sum(!is.na(`Post GTT (AOC)`))))

ggplot(group_means, aes(x = mean_Comp1, y = mean_GTT, label = Group, color = Group)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = mean_GTT - se_GTT, ymax = mean_GTT + se_GTT), width=0) +
  geom_errorbarh(aes(xmin = mean_Comp1 - se_Comp1, xmax = mean_Comp1 + se_Comp1), height=0) +
  geom_text(nudge_y = 500, show.legend = FALSE) +
  theme_minimal(base_size=14) +
  labs(x = "Molecular aging axis (Comp1)", y = "Glucose tolerance (GTT AUC)",
       title = "Group means: systemic glucose metabolism vs molecular aging signature")
#---------------------------------------------------------------------------------------#

library(broom)
library(dplyr)

# Fit separate linear models by group
group_lms <- old_only %>%
  filter(!is.na(`Post GTT (AOC)`)) %>%
  group_by(Group) %>%
  do({
    fit <- lm(`Post GTT (AOC)` ~ Comp1, data = .)
    glance_fit <- broom::glance(fit)   # R², adj.R², p-value, etc.
    tidy_fit   <- broom::tidy(fit)     # slope/intercept
    data.frame(
      r2 = glance_fit$r.squared,
      p  = tidy_fit$p.value[2],
      slope = tidy_fit$estimate[2]
    )
  })

group_lms
#---------------------------------------#

old_only2 <- old_only %>%
  mutate(Group2 = ifelse(Group %in% c("Old_PwrIRap","Old_PwrFRap"),
                         "Old_PwrRapa", as.character(Group)))

lm_combined <- lm(`Post GTT (AOC)` ~ Comp1 * Group2, data = old_only2)
anova(lm_combined)

# Extract R2 and slope for pooled rapa only
rapa_fit <- lm(`Post GTT (AOC)` ~ Comp1,
               data = filter(old_only2, Group2 == "Old_PwrRapa"))
summary(rapa_fit)

ggplot(filter(old_only2, Group2 == "Old_PwrRapa"),
       aes(x = Comp1, y = `Post GTT (AOC)`)) +
  geom_point(size=3, color="#c0392b") +
  geom_smooth(method="lm", se=TRUE, color="#c0392b", fill="#e74c3c", alpha=0.2) +
  theme_minimal(base_size=14) +
  labs(title = "Rapamycin-treated mice",
       subtitle = sprintf("Slope = %.0f, R² = %.2f, p = %.3f",
                          coef(rapa_fit)[2], summary(rapa_fit)$r.squared,
                          summary(rapa_fit)$coefficients[2,4]),
       x = "Molecular aging score (Comp1)",
       y = "Glucose tolerance (GTT AUC)")
#------------------------------------------------------#

old_only

old_cor_all <- old_only %>%
  dplyr::select(where(is.numeric), -Comp1) %>%
  names() %>%
  lapply(function(v){
    x <- old_only$Comp1; y <- old_only[[v]]
    ok <- complete.cases(x,y)
    if(sum(ok) >= 4){
      ct <- cor.test(x[ok], y[ok], method="spearman", exact=FALSE)
      data.frame(var=v, r=ct$estimate, p=ct$p.value, n=sum(ok))
    } else data.frame(var=v, r=NA, p=NA, n=sum(ok))
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(FDR = p.adjust(p, "fdr")) %>%
  dplyr::arrange(FDR)

print(old_cor_all, digits=3)
#----------------------------------------------------#

# Post Insulin Sensitivity (AUC)
ggplot(old_only, aes(x = Comp1, y = `Post Insulin Sensitivity (AUC)`, color = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Insulin Sensitivity (AUC)",
    title = "Molecular aging axis correlates with impaired insulin sensitivity (Post AUC)"
  )

# Delta Insulin Sensitivity (AUC)
ggplot(old_only, aes(x = Comp1, y = `Delta Insulin Sensitivity (AUC)`, color = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Δ Insulin Sensitivity (AUC)",
    title = "Molecular aging axis correlates with change in insulin sensitivity"
  )
#--------------------------------------------------------------------------------#

#========================================================================#
#---Correlation of Physio traits to aging axis among old groups----------#
#========================================================================#
physio_diablo_corr_barplot <- ggplot(old_cor_all, aes(x = reorder(var, r), y = r,
                        fill = p < 0.05)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70")) +
  labs(
    x = "",
    y = "Spearman r (DIABLO Comp1 vs Physiology)",
    title = "Physiologic correlates of molecular aging axis (old groups only)",
    fill = "p < 0.05"
  ) +
  geom_text(aes(label = sprintf("n=%d", n)),
            hjust = -0.2, size = 3, color = "black")

pdf("physio_diablo_corr_barplot.pdf", width = 8, height = 6)
print(physio_diablo_corr_barplot)
dev.off()
#--------------------------------------------------------------#

# Extract stats from your correlation table
stats_post_auc <- old_cor_all %>% 
  filter(var == "Post Insulin Sensitivity (AUC)") %>%
  dplyr::mutate(r2 = r^2)

ITT_diablo_sctrplot <- ggplot(old_only, aes(x = Comp1, y = `Post Insulin Sensitivity (AUC)`)) +
  geom_point(size = 3, color = "firebrick") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Insulin Sensitivity (AUC)",
    title = "Molecular aging axis vs insulin sensitivity in old mice",
    subtitle = sprintf("Spearman r = %.3f (R² = %.3f, p = %.3g)", 
                       stats_post_auc$r, stats_post_auc$r2, stats_post_auc$p)
  )
pdf("ITT_diablo_sctrplot.pdf", width = 8, height = 6)
print(ITT_diablo_sctrplot)
dev.off()


#----------------------------------------------------------------#

stats_gtt <- old_cor_all %>% 
  filter(var == "Post GTT (AOC)") %>%
  dplyr::mutate(r2 = r^2)

GTT_diablo_sctrplot <- ggplot(old_only, aes(x = Comp1, y = `Post GTT (AOC)`)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Glucose Tolerance (GTT AUC)",
    title = "Molecular aging axis vs glucose tolerance in old mice",
    subtitle = sprintf("Spearman r = %.3f (R² = %.3f, p = %.3g)", 
                       stats_gtt$r, stats_gtt$r2, stats_gtt$p)
  )

pdf("GTT_diablo_sctrplot.pdf", width = 8, height = 6)
print(GTT_diablo_sctrplot)
dev.off()
#-----------------------------------------------------------------#

stats_delta_auc <- old_cor_all %>% 
  filter(var == "Delta Insulin Sensitivity (AUC)") %>%
  dplyr::mutate(r2 = r^2)

ggplot(old_only, aes(x = Comp1, y = `Delta Insulin Sensitivity (AUC)`)) +
  geom_point(size = 3, color = "darkorange3") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Δ Insulin Sensitivity (AUC)",
    title = "Molecular aging axis vs change in insulin sensitivity in old mice",
    subtitle = sprintf("Spearman r = %.3f (R² = %.3f, p = %.3g)", 
                       stats_delta_auc$r, stats_delta_auc$r2, stats_delta_auc$p)
  )
#-------------------------------------------------------------#


library(ggplot2)
library(dplyr)
library(ggpmisc)

# Remove sample T_25 from old_only before everything else
old_plot_df <- old_only %>%
  filter(Sample != "T_25") %>%       # <--- Exclude this sample
  mutate(Group3 = case_when(
    Group %in% c("Old_PwrIRap", "Old_PwrFRap") ~ "Old_PwrRapa",
    TRUE ~ as.character(Group)
  ))

# Recompute linear model for rapamycin-only subset
rapa_fit <- lm(`Post GTT (AOC)` ~ Comp1,
               data = filter(old_plot_df, Group3 == "Old_PwrRapa"))
rapa_summary <- summary(rapa_fit)
r2_rapa <- round(rapa_summary$r.squared, 3)
pval_rapa <- signif(coef(summary(rapa_fit))[2, 4], 3)

# Plot
GTT_diablo_grp_sctrplot <- ggplot(old_plot_df, aes(x = Comp1, y = `Post GTT (AOC)`, color = Group3)) +
  # Scatter points
  geom_point(size = 3, alpha = 0.9) +
  # PCA-style group ellipses
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15, color = NA, level = 0.68) +
  # Regression line for Rapa group only
  geom_smooth(data = filter(old_plot_df, Group3 == "Old_PwrRapa"),
              method = "lm", se = FALSE, color = "#c0392b", size = 1.2, linetype = "solid") +
  # Colors
  scale_color_manual(values = c(
    "Old_SedVeh" = "#2c3e50",  # dark blue-gray
    "Old_PwrVeh" = "#27ae60",  # green
    "Old_PwrRapa" = "#e74c3c"  # red
  )) +
  scale_fill_manual(values = c(
    "Old_SedVeh" = "#2c3e50",
    "Old_PwrVeh" = "#27ae60",
    "Old_PwrRapa" = "#e74c3c"
  )) +
  # Theme and labels
  theme_minimal(base_size = 14) +
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Glucose Tolerance Test AUC",
    title = "Molecular aging axis vs glucose tolerance in old mice",
    subtitle = sprintf("Rapamycin group: R² = %.3f, p = %.3g", r2_rapa, pval_rapa),
    color = "Group",
    fill = "Group"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf("GTT_diablo_grp_sctrplot.pdf", width = 8, height = 6)
print(GTT_diablo_grp_sctrplot)
dev.off()


old_plot_df
#------------------------------------------------#

# Compute the overall regression
fit_all <- lm(`Post GTT (AOC)` ~ Comp1, data = old_plot_df)
summary_all <- summary(fit_all)
r2_all <- round(summary_all$r.squared, 3)
pval_all <- signif(coef(summary_all)[2, 4], 3)

# Plot
GTT_diablo_grp_sctrplot <- ggplot(old_plot_df, aes(x = Comp1, y = `Post GTT (AOC)`, color = Group3)) +
  # Points by group
  geom_point(size = 3, alpha = 0.9) +
  # Ellipses by group
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15, color = NA, level = 0.68) +
  # Single dashed regression line across all points
  geom_smooth(aes(x = Comp1, y = `Post GTT (AOC)`),
              method = "lm", se = FALSE, color = "black",
              linetype = "dashed", size = 1.2) +
  # Annotate R² and p-value
  annotate("text",
           x = min(old_plot_df$Comp1, na.rm = TRUE),
           y = max(old_plot_df$`Post GTT (AOC)`, na.rm = TRUE),
           label = sprintf("R² = %.3f, p = %.3g", r2_all, pval_all),
           hjust = 0, vjust = 1.2, size = 5, fontface = "italic") +
  # Color settings
  scale_color_manual(values = c(
    "Old_SedVeh" = "#2c3e50",  # dark blue-gray
    "Old_PwrVeh" = "#27ae60",  # green
    "Old_PwrRapa" = "#e74c3c"  # red
  )) +
  scale_fill_manual(values = c(
    "Old_SedVeh" = "#2c3e50",
    "Old_PwrVeh" = "#27ae60",
    "Old_PwrRapa" = "#e74c3c"
  )) +
  # Labels and theme
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Glucose Tolerance Test AUC",
    title = "Molecular aging axis vs glucose tolerance in old mice",
    subtitle = "Single regression line across all groups",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save to PDF
pdf("GTT_diablo_grp_sctrplot_singleline_stats.pdf", width = 8, height = 6)
print(GTT_diablo_grp_sctrplot)
dev.off()
#------------------------------------------------#

# Compute the overall regression for Insulin Sensitivity (AUC)
fit_itt_all <- lm(`Post Insulin Sensitivity (AUC)` ~ Comp1, data = old_plot_df)
summary_itt_all <- summary(fit_itt_all)
r2_itt <- round(summary_itt_all$r.squared, 3)
pval_itt <- signif(coef(summary_itt_all)[2, 4], 3)

# Plot: same styling as GTT_diablo_grp_sctrplot
ITT_diablo_grp_sctrplot <- ggplot(old_plot_df, aes(x = Comp1, y = `Post Insulin Sensitivity (AUC)`, color = Group3)) +
  # Scatter points by group
  geom_point(size = 3, alpha = 0.9) +
  # PCA-style group ellipses
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15, color = NA, level = 0.68) +
  # Single dashed regression line across all groups
  geom_smooth(aes(x = Comp1, y = `Post Insulin Sensitivity (AUC)`),
              method = "lm", se = FALSE, color = "black",
              linetype = "dashed", size = 1.2) +
  # Annotate R² and p-value
  annotate("text",
           x = min(old_plot_df$Comp1, na.rm = TRUE),
           y = max(old_plot_df$`Post Insulin Sensitivity (AUC)`, na.rm = TRUE),
           label = sprintf("R² = %.3f, p = %.3g", r2_itt, pval_itt),
           hjust = 0, vjust = 1.2, size = 5, fontface = "italic") +
  # Custom colors (consistent with GTT plot)
  scale_color_manual(values = c(
    "Old_SedVeh" = "#2c3e50",  # dark blue-gray
    "Old_PwrVeh" = "#27ae60",  # green
    "Old_PwrRapa" = "#e74c3c"  # red
  )) +
  scale_fill_manual(values = c(
    "Old_SedVeh" = "#2c3e50",
    "Old_PwrVeh" = "#27ae60",
    "Old_PwrRapa" = "#e74c3c"
  )) +
  # Labels and theme
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Insulin Sensitivity (AUC)",
    title = "Molecular aging axis vs insulin sensitivity in old mice",
    subtitle = "Single regression line across all groups",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save to PDF
pdf("ITT_diablo_grp_sctrplot_singleline_stats.pdf", width = 8, height = 6)
print(ITT_diablo_grp_sctrplot)
dev.off()
#-------------------------------------------------#

names(merged)[grepl("GXT", names(merged))]
# Compute the overall regression for Post GXT time (sec) — note the space at the end
fit_gxt_all <- lm(`Post GXT time (sec) ` ~ Comp1, data = old_plot_df)
summary_gxt_all <- summary(fit_gxt_all)
r2_gxt <- round(summary_gxt_all$r.squared, 3)
pval_gxt <- signif(coef(summary_gxt_all)[2, 4], 3)

# Plot
GXT_diablo_grp_sctrplot <- ggplot(old_plot_df, aes(x = Comp1, y = `Post GXT time (sec) `, color = Group3)) +
  # Scatter points by group
  geom_point(size = 3, alpha = 0.9) +
  # PCA-style ellipses
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15,
               color = NA, level = 0.68) +
  # Single dashed regression line across all groups
  geom_smooth(aes(x = Comp1, y = `Post GXT time (sec) `),
              method = "lm", se = FALSE, color = "black",
              linetype = "dashed", size = 1.2) +
  # Annotate R² and p-value
  annotate("text",
           x = min(old_plot_df$Comp1, na.rm = TRUE),
           y = max(old_plot_df$`Post GXT time (sec) `, na.rm = TRUE),
           label = sprintf("R² = %.3f, p = %.3g", r2_gxt, pval_gxt),
           hjust = 0, vjust = 1.2, size = 5, fontface = "italic") +
  # Custom colors (same as GTT/ITT)
  scale_color_manual(values = c(
    "Old_SedVeh" = "#2c3e50",  # dark blue-gray
    "Old_PwrVeh" = "#27ae60",  # green
    "Old_PwrRapa" = "#e74c3c"  # red
  )) +
  scale_fill_manual(values = c(
    "Old_SedVeh" = "#2c3e50",
    "Old_PwrVeh" = "#27ae60",
    "Old_PwrRapa" = "#e74c3c"
  )) +
  # Labels and theme
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "GXT Time (sec)",
    title = "Molecular aging axis vs exercise tolerance in old mice",
    subtitle = "Single regression line across all groups",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save to PDF
pdf("GXT_diablo_grp_sctrplot_singleline_stats.pdf", width = 8, height = 6)
print(GXT_diablo_grp_sctrplot)
dev.off()
#----------------------------------------#
#===========================================#

#-----------------------------------------------------------#
# Compute regression ONLY for rapamycin-treated groups
#-----------------------------------------------------------#
rapa_data <- old_plot_df %>% 
  filter(Group3 == "Old_PwrRapa")  # combined IRap + FRap group

fit_itt_rapa <- lm(`Post Insulin Sensitivity (AUC)` ~ Comp1, data = rapa_data)
summary_itt_rapa <- summary(fit_itt_rapa)
r2_rapa <- round(summary_itt_rapa$r.squared, 3)
pval_rapa <- signif(coef(summary_itt_rapa)[2, 4], 3)

#-----------------------------------------------------------#
# Plot with regression line for rapamycin group only
#-----------------------------------------------------------#
ITT_diablo_grp_sctrplot_rapa <- ggplot(old_plot_df, 
                                       aes(x = Comp1, 
                                           y = `Post Insulin Sensitivity (AUC)`, 
                                           color = Group3)) +
  # Scatter points by group
  geom_point(size = 3, alpha = 0.9) +
  # Group ellipses
  stat_ellipse(aes(fill = Group3), geom = "polygon", alpha = 0.15, 
               color = NA, level = 0.68) +
  # Regression line for rapamycin group only
  geom_smooth(data = rapa_data,
              aes(x = Comp1, y = `Post Insulin Sensitivity (AUC)`),
              method = "lm", se = FALSE, color = "#e74c3c", 
              size = 1.2, linetype = "solid") +
  # Annotate R² and p-value for rapamycin regression
  annotate("text",
           x = min(rapa_data$Comp1, na.rm = TRUE),
           y = max(rapa_data$`Post Insulin Sensitivity (AUC)`, na.rm = TRUE),
           label = sprintf("Rapamycin only: R² = %.3f, p = %.3g", 
                           r2_rapa, pval_rapa),
           hjust = 0, vjust = 1.2, size = 5, fontface = "italic", 
           color = "#e74c3c") +
  # Colors (consistent with GTT and GXT)
  scale_color_manual(values = c(
    "Old_SedVeh" = "#2c3e50",  # dark blue-gray
    "Old_PwrVeh" = "#27ae60",  # green
    "Old_PwrRapa" = "#e74c3c"  # red
  )) +
  scale_fill_manual(values = c(
    "Old_SedVeh" = "#2c3e50",
    "Old_PwrVeh" = "#27ae60",
    "Old_PwrRapa" = "#e74c3c"
  )) +
  # Labels and theme
  labs(
    x = "DIABLO Comp1 (higher = molecularly older)",
    y = "Insulin Sensitivity (AUC)",
    title = "Molecular aging axis vs insulin sensitivity in old mice",
    subtitle = "Regression line shown for rapamycin-treated mice only",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

#-----------------------------------------------------------#
# Save the plot
#-----------------------------------------------------------#
pdf("ITT_diablo_grp_sctrplot_rapa_only.pdf", width = 8, height = 6)
print(ITT_diablo_grp_sctrplot_rapa)
dev.off()

