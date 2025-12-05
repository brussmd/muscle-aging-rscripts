#Attempt to analyze metabolomic data

setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR")

# Load necessary packages
library(tidyverse)
library(pheatmap)
library(ggplot2)

# Read the CSV file (adjust the path as needed)
pwr_rapa_muscle_metabo_data <- read_csv("Konopka_Muscle_HILIC.csv")

# View the first few rows
head(pwr_rapa_muscle_metabo_data)

# Set metabolite names as rownames and convert to matrix
metabo_matrix <- pwr_rapa_muscle_metabo_data %>%
  column_to_rownames(var = colnames(pwr_rapa_muscle_metabo_data)[1]) %>%
  as.matrix()

head(metabo_matrix)
rownames(metabo_matrix)

# Create a smaller dataframe for analysis
metabo_df <- pwr_rapa_muscle_metabo_data %>%
  dplyr::select(Sample, Group, everything()) %>%
  filter(Group %in% c("OS", "YS"))  # Focus on 2 groups of interest

# Pivot to long format for ANOVA/visualization
metabo_long <- metabo_df %>%
  pivot_longer(
    cols = -c(Sample, Group),
    names_to = "Metabolite",
    values_to = "Abundance"
  )

metabo_long
# Filter to only OS vs YS
aging_df <- metabo_long %>% filter(Group %in% c("OS", "YS"))

# Run t-test per metabolite
aging_stats <- aging_df %>%
  group_by(Metabolite) %>%
  summarise(
    p_value = t.test(Abundance ~ Group)$p.value,
    mean_OS = mean(Abundance[Group == "OS"]),
    mean_YS = mean(Abundance[Group == "YS"]),
    log2FC = log2(mean_OS + 1e-6) - log2(mean_YS + 1e-6)  # add small value to avoid log(0)
  ) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adj_p)


aging_stats

# View top results
aging_stats %>% filter(p_value < 0.05) %>% slice_head(n = 20)
#--------------------------------------------------------------#

#------Volcano Plot for All Metabolites-----------------------#

# Filter significant metabolites
signif_metabs <- aging_stats %>%
  filter(p_value < 0.05)

# Volcano plot
metabo_volcano <- ggplot(aging_stats, aes(x = log2FC, y = -log10(p_value))) +
  geom_point(color = "gray60", size = 2) +
  geom_point(data = signif_metabs, aes(x = log2FC, y = -log10(p_value)),
             color = "red", size = 2) +
  geom_text_repel(data = signif_metabs,
                  aes(label = Metabolite),
                  size = 3.5,
                  max.overlaps = 100,
                  box.padding = 0.4,
                  point.padding = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  scale_x_continuous(limits = c(-3, 3)) +  # Set your desired range here
  theme_minimal() +
  labs(title = "Volcano Plot of Metabolite Changes with Age",
       x = "Log2 Fold Change (Old vs Young)",
       y = "-Log10(p-value)")

# Save as PDF (adjust size as needed)
pdf("metabo_volcano.pdf", width = 6, height = 6)
print(metabo_volcano)
dev.off()


#----------------------------------------------#
# Subset to YS and OS only
ys_os_data <- pwr_rapa_muscle_metabo_data %>%
  filter(Group %in% c("YS", "OS")) %>%
  select(Sample, Group, everything())

ys_os_data

# Save to CSV
write_csv(ys_os_data, "metabo_YS_vs_OS.csv")
#-----------------------------------------------#

# Filter to only OS vs OV
metabo_long
# Filter to only OS vs OV
OV_aging_df <- metabo_long %>% filter(Group %in% c("OS", "OV"))

OV_aging_clean <- OV_aging_df %>%
  group_by(Metabolite) %>%
  filter(n_distinct(Group) == 2) %>%
  ungroup()


# Run t-test per metabolite
OV_aging_stats <- OV_aging_clean %>%
  group_by(Metabolite) %>%
  summarise(
    p_value = t.test(Abundance ~ Group)$p.value,
    mean_OS = mean(Abundance[Group == "OS"], na.rm = TRUE),
    mean_OV = mean(Abundance[Group == "OV"], na.rm = TRUE),
    log2FC = log2(mean_OV + 1e-6) - log2(mean_OS + 1e-6)  # OV vs OS
  ) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adj_p)


OV_aging_stats


















pwr_aging_df <- metabo_long %>% filter(Group %in% c("OS", "OV"))
pwr_aging_df

# 1) Clean & restrict to the two groups of interest
pwr_aging_df2 <- pwr_aging_df %>%
  mutate(Group = trimws(Group)) %>%
  filter(Group %in% c("OS", "OV"))

pwr_aging_df2

# 2) Helper: safe t-test that returns NA on error
safe_t_p <- function(x, g) {
  tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
}

# 3) Compute stats only where both groups are present
pwr_aging_stats <- pwr_aging_df2 %>%
  group_by(Metabolite) %>%
  filter(n_distinct(Group) == 2) %>%                 # require OS & OV present
  summarise(
    p_value = safe_t_p(Abundance, Group),
    mean_OS = mean(Abundance[Group == "OS"], na.rm = TRUE),
    mean_OV = mean(Abundance[Group == "OV"], na.rm = TRUE),
    log2FC  = log2(mean_OV + 1e-6) - log2(mean_OS + 1e-6) # avoid log(0)
  ) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adj_p)

pwr_aging_stats

# View top results
pwr_aging_stats %>% filter(p_value < 0.05) %>% slice_head(n = 20)

inflammaging_metabolites <- c(
  "Kynurenic acid",
  "DL-Methionine sulfoxide",
  "N-Acetyl-DL-tryptophan",
  "DL-2-Aminoadipic acid"
)

pwr_aging_stats %>%
  filter(Metabolite %in% inflammaging_metabolites)
#--------------------------------------------------#

# Assuming you already loaded this tibble:
data <- pwr_rapa_muscle_metabo_data

# View column names to confirm
colnames(data)[1:10]

# 2️⃣ Identify metabolite columns
metabolite_cols <- colnames(data)[3:ncol(data)]

# 3️⃣ Compute per-metabolite minimum non-zero values
min_nonzero <- sapply(data[metabolite_cols], function(x) {
  min(x[x > 0], na.rm = TRUE)
})

# 4️⃣ Replace 0s with half-minimum per metabolite
# Make a copy to modify
imputed_data <- data

for (met in metabolite_cols) {
  x <- imputed_data[[met]]
  # Replace zeros
  x[x == 0] <- min_nonzero[met] / 2
  imputed_data[[met]] <- x
}

# 5️⃣ Inspect a few rows to confirm
head(imputed_data)

# 6️⃣ Save to CSV for MetaboAnalyst
write_csv(imputed_data, "MetaboData_Imputed.csv")
#----------------------------------------------------#

# 1️⃣ Subset to OS and YS groups
compare_data <- imputed_data %>%
  filter(Group %in% c("OS", "YS"))

# 2️⃣ Get the metabolite columns
metabolite_cols <- colnames(compare_data)[3:ncol(compare_data)]

# 3️⃣ Create a function to run t-tests per metabolite
t_test_results <- map_dfr(metabolite_cols, function(met) {
  # Extract values
  vals <- compare_data[[met]]
  group <- compare_data$Group
  
  # Run t-test
  ttest <- t.test(vals ~ group)
  
  # Return as a tibble row
  tibble(
    Metabolite = met,
    Mean_YS = mean(vals[group == "YS"]),
    Mean_OS = mean(vals[group == "OS"]),
    Log2_FC_OS_vs_YS = log2(mean(vals[group == "OS"]) / mean(vals[group == "YS"])),
    P_value = ttest$p.value
  )
})

# 4️⃣ Adjust p-values for multiple testing (Benjamini-Hochberg FDR)
t_test_results <- t_test_results %>%
  mutate(FDR = p.adjust(P_value, method = "fdr"))

# 5️⃣ Arrange by FDR to see the top hits
t_test_results <- t_test_results %>%
  arrange(FDR)

# 6️⃣ Preview the results
print(t_test_results, n=20)

# For example, Kynurenic acid
imputed_data %>%
  filter(Group %in% c("YS", "YV")) %>%
  select(Sample, Group, `3-Methyl-L-Histidine`) %>%
  arrange(Group)


imputed_data %>%
  filter(Group %in% c("YS", "YV")) %>%
  ggplot(aes(x = Group, y = `3-Methyl-L-Histidine`)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "3-Methyl-L-Histidine  levels by group")
#---------------------------------------------------------#

# Identify metabolite columns
metabolite_cols <- colnames(imputed_data)[3:ncol(imputed_data)]

# Create log2-transformed tibble
imputed_log <- imputed_data %>%
  mutate(across(all_of(metabolite_cols), log2))

head(imputed_log)

# Subset to YS and OS
compare_data <- imputed_log %>%
  filter(Group %in% c("YS", "OS"))

table(compare_data$Group)

# Run t-tests per metabolite
t_test_results_log <- map_dfr(metabolite_cols, function(met) {
  vals <- compare_data[[met]]
  group <- compare_data$Group
  ttest <- t.test(vals ~ group)
  tibble(
    Metabolite = met,
    Mean_YS = mean(vals[group == "YS"]),
    Mean_OS = mean(vals[group == "OS"]),
    Log2_FC_OS_vs_YS = mean(vals[group == "OS"]) - mean(vals[group == "YS"]),
    P_value = ttest$p.value
  )
}) %>%
  mutate(FDR = p.adjust(P_value, method = "fdr")) %>%
  arrange(P_value)

print(t_test_results_log, n=20)

aging_metabolite_set_effectsize <- t_test_results_log %>%
  filter(abs(Log2_FC_OS_vs_YS) > 1)

dim(aging_metabolite_set_effectsize %>% filter(Log2_FC_OS_vs_YS>0) %>% filter(P_value < 0.05))

aging_metabolite_set_effectsize %>%
  filter(Log2_FC_OS_vs_YS >0) %>%
  filter(P_value < 0.05)
#-------------------------------------------------#

#----Aging Metabolite Heatmap--------------------#

aging_metabolite_names <- aging_stats %>%
  filter(p_value < 0.05) %>%
  pull(Metabolite) %>%
  unique() %>%
  sort()

length(aging_metabolite_names)

heatmap_data <- imputed_log %>%
  filter(Group %in% c("YS", "OS")) %>%
  dplyr::select(Sample, Group, all_of(aging_metabolite_names))

head(heatmap_data)

# Define desired order of groups
heatmap_data <- heatmap_data %>%
  mutate(Group = factor(Group, levels = c("YS", "OS"))) %>%
  arrange(Group)

# Extract numeric matrix
X <- heatmap_data %>%
  dplyr::select(all_of(aging_metabolite_names)) %>%
  as.matrix()

rownames(X) <- heatmap_data$Sample

# Z-score
X_scaled <- scale(X)

# Transpose
X_scaled_t <- t(X_scaled)

# Annotation for columns (samples)
annotation_col <- data.frame(Group = heatmap_data$Group)
rownames(annotation_col) <- heatmap_data$Sample

# Heatmap
os_metab_heatmap <- pheatmap(
  mat = X_scaled_t,
  annotation_col = annotation_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,                 # hide row dendrogram
  legend = TRUE,                      # keep heatmap legend
  annotation_legend = FALSE,          # hide group legend
  legend_breaks = c(-1, 0, 1),        # fewer ticks = smaller legend
  legend_labels = c("-1", "0", "+1"),
  border_color = NA,                  # cleaner look
  fontsize = 8,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Aging-Associated Upregulated Metabolites (OS vs YS)",
  row_names_side = "left",
  breaks = seq(-1.2, 1.2, length.out = 101)
)

# Save as PDF (adjust size as needed)
pdf("os_metab_heatmap.pdf", width = 4.25, height = 5.5)
print(os_metab_heatmap)
dev.off()
#--------------------------------------------------------------------#

#---------Rank-Based Enrichment Analysis----------------------------#

#--------------- Compare YS--------------------------------#

# 1) Subset to OS and YS samples
compare_data_OS_YS <- imputed_log %>%
  filter(Group %in% c("OS", "YS"))

# 2) Perform t-tests for each metabolite
t_test_OS_YS <- map_dfr(metabolite_cols, function(met) {
  vals <- compare_data_OS_YS[[met]]
  group <- compare_data_OS_YS$Group
  ttest <- t.test(vals ~ group)
  tibble(
    Metabolite = met,
    Mean_OS = mean(vals[group == "OS"]),
    Mean_YS = mean(vals[group == "YS"]),
    Log2_FC_OS_vs_YS = mean(vals[group == "OS"]) - mean(vals[group == "YS"]),
    P_value = ttest$p.value
  )
})

t_test_OS_YS

# 3) Calculate ranked statistic: -log10(p) * logFC
OS_rank_df <- t_test_OS_YS %>%
  filter(!is.na(Log2_FC_OS_vs_YS), !is.na(P_value)) %>%
  mutate(
    p_safe = pmax(P_value, 1e-300), # avoid Inf
    stat   = -log10(p_safe) * Log2_FC_OS_vs_YS
  ) %>%
  group_by(Metabolite) %>%
  summarise(stat = mean(stat), .groups = "drop") %>%
  arrange(desc(stat))

# 4) Named numeric vector for fgsea
OS_ranks <- setNames(OS_rank_df$stat, OS_rank_df$Metabolite)

# 5) Run fgsea
OS_fgsea_res <- fgsea(
  pathways = list(AgingSignature = aging_metabolite_names),
  stats    = OS_ranks,
  nperm    = 1000
)

print(OS_fgsea_res)

# 6) Plot enrichment
OS_lsea_plot <- plotEnrichment(
  pathway = aging_metabolite_names,
  stats   = OS_ranks
) + ggtitle("Enrichment of Aging Metabolite Set in OS vs YS")

# Save enrichment plot
pdf("OS_lsea_plot.pdf", width = 4.25, height = 5.5)
print(OS_lsea_plot)
dev.off()

#-------------------------------------------------------------#
#---- Flip point for annotation -----------------------------#
df_bottom <- tibble(
  Rank = 1:length(OS_ranks),
  Metabolite = names(OS_ranks),
  Stat = OS_ranks
)

last_pos_rank <- max(which(df_bottom$Stat > 0))
flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA

# Bottom panel plot
os_vs_ys_p_bottom <- ggplot(df_bottom, aes(x = Rank, y = Stat, fill = Stat)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       name = "Ranked Statistic\n(logFC)") +
  { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
  labs(x = "Metabolite Rank", y = "Statistic") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(os_vs_ys_p_bottom)

# Save bottom panel plot
pdf("os_vs_ys_p_bottom.pdf", width = 4.25, height = 5.5)
print(os_vs_ys_p_bottom)
dev.off()

#---------------------------------------------------------#
#--------------- Compare OV--------------------------------#

# 1) Subset to OV and OS samples
compare_data_OV_OS <- imputed_log %>%
  filter(Group %in% c("OV", "OS"))

# 2) Perform t-tests for each metabolite
t_test_OV_OS <- map_dfr(metabolite_cols, function(met) {
  vals <- compare_data_OV_OS[[met]]
  group <- compare_data_OV_OS$Group
  ttest <- t.test(vals ~ group)
  tibble(
    Metabolite = met,
    Mean_OV = mean(vals[group == "OV"]),
    Mean_OS = mean(vals[group == "OS"]),
    Log2_FC_OV_vs_OS = mean(vals[group == "OV"]) - mean(vals[group == "OS"]),
    P_value = ttest$p.value
  )
})

t_test_OV_OS

# 3) Calculate ranked statistic: -log10(p) * logFC
OV_rank_df <- t_test_OV_OS %>%
  filter(!is.na(Log2_FC_OV_vs_OS), !is.na(P_value)) %>%
  mutate(
    p_safe = pmax(P_value, 1e-300), # avoid Inf
    stat   = -log10(p_safe) * Log2_FC_OV_vs_OS
  ) %>%
  group_by(Metabolite) %>%
  summarise(stat = mean(stat), .groups = "drop") %>%
  arrange(desc(stat))

# 4) Named numeric vector for fgsea
OV_ranks <- setNames(OV_rank_df$stat, OV_rank_df$Metabolite)

# 5) Run fgsea
OV_fgsea_res <- fgsea(
  pathways = list(AgingSignature = aging_metabolite_names),
  stats    = OV_ranks,
  nperm    = 1000
)

print(OV_fgsea_res)

# 6) Plot enrichment
OV_lsea_plot <- plotEnrichment(
  pathway = aging_metabolite_names,
  stats   = OV_ranks
) + ggtitle("Enrichment of Aging Metabolite Set in OV vs OS")

# Save enrichment plot
pdf("OV_lsea_plot.pdf", width = 4.25, height = 5.5)
print(OV_lsea_plot)
dev.off()
#--------------------------------------------------------------#

#-------------------------------------------------------------#
#---- Flip point for annotation -----------------------------#
df_bottom <- tibble(
  Rank = 1:length(OV_ranks),
  Metabolite = names(OV_ranks),
  Stat = OV_ranks
)

last_pos_rank <- max(which(df_bottom$Stat > 0))
flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA

# Bottom panel plot
ov_vs_os_p_bottom <- ggplot(df_bottom, aes(x = Rank, y = Stat, fill = Stat)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       name = "Ranked Statistic\n(logFC)") +
  { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
  labs(x = "Metabolite Rank", y = "Statistic") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(ov_vs_os_p_bottom)

# Save bottom panel plot
pdf("ov_vs_os_p_bottom.pdf", width = 4.25, height = 5.5)
print(ov_vs_os_p_bottom)
dev.off()

#--------------- Compare OFR--------------------------------#

# 1) Subset to OFR and OS samples
compare_data_OFR_OS <- imputed_log %>%
  filter(Group %in% c("OFR", "OS"))

# 2) Perform t-tests for each metabolite
t_test_OFR_OS <- map_dfr(metabolite_cols, function(met) {
  vals <- compare_data_OFR_OS[[met]]
  group <- compare_data_OFR_OS$Group
  ttest <- t.test(vals ~ group)
  tibble(
    Metabolite = met,
    Mean_OFR = mean(vals[group == "OFR"]),
    Mean_OS = mean(vals[group == "OS"]),
    Log2_FC_OFR_vs_OS = mean(vals[group == "OFR"]) - mean(vals[group == "OS"]),
    P_value = ttest$p.value
  )
})

# 3) Create named vector of log2 fold-changes
OFR_ranks <- t_test_OFR_OS %>%
  arrange(desc(Log2_FC_OFR_vs_OS)) %>%
  dplyr::select(Metabolite, Log2_FC_OFR_vs_OS) %>%
  deframe()

# 4) Remove any NAs
OFR_ranks <- OFR_ranks[!is.na(OFR_ranks)]

# 5) Confirm overlap between ranks and the metabolite set
intersect(names(OFR_ranks), aging_metabolite_names)

# 6) Run fgsea
OFR_fgsea_res <- fgsea(
  pathways = list(AgingSignature = aging_metabolite_names),
  stats = OFR_ranks,
  nperm = 1000
)

# 7) View results
print(OFR_fgsea_res)

# 8) Plot enrichment
OFR_msea_plot <-plotEnrichment(
  pathway = aging_metabolite_names,
  stats = OFR_ranks
) + ggtitle("Enrichment of Aging Metabolite Set in OFR vs OS")

# Save enrichment plot
pdf("OFR_msea_plot.pdf", width = 4.25, height = 5.5)
print(OFR_msea_plot)
dev.off()
#--------------------------------------------------------------#
#---- Flip point for annotation -----------------------------#
OFR_df_bottom <- tibble(
  Rank = 1:length(OFR_ranks),
  Metabolite = names(OFR_ranks),
  Stat = OFR_ranks
)

last_pos_rank <- max(which(OFR_df_bottom$Stat > 0))
flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA

# Bottom panel plot
ofr_vs_os_p_bottom <- ggplot(OFR_df_bottom, aes(x = Rank, y = Stat, fill = Stat)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       name = "Ranked Statistic\n(logFC)") +
  { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
  labs(x = "Metabolite Rank", y = "Statistic") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(ofr_vs_os_p_bottom)

# Save bottom panel plot
pdf("ofr_vs_os_p_bottom.pdf", width = 4.25, height = 5.5)
print(ofr_vs_os_p_bottom)
dev.off()
#-----------------------------------------------------------#

#--------------- Compare OIR--------------------------------#

# 1) Subset to OIR and OS samples
compare_data_OIR_OS <- imputed_log %>%
  filter(Group %in% c("OIR", "OS"))

# 2) Perform t-tests for each metabolite
t_test_OIR_OS <- map_dfr(metabolite_cols, function(met) {
  vals <- compare_data_OIR_OS[[met]]
  group <- compare_data_OIR_OS$Group
  ttest <- t.test(vals ~ group)
  tibble(
    Metabolite = met,
    Mean_OIR = mean(vals[group == "OIR"]),
    Mean_OS = mean(vals[group == "OS"]),
    Log2_FC_OIR_vs_OS = mean(vals[group == "OIR"]) - mean(vals[group == "OS"]),
    P_value = ttest$p.value
  )
})

# 3) Create named vector of log2 fold-changes
OIR_ranks <- t_test_OIR_OS %>%
  arrange(desc(Log2_FC_OIR_vs_OS)) %>%
  dplyr::select(Metabolite, Log2_FC_OIR_vs_OS) %>%
  deframe()

# 4) Remove any NAs
OIR_ranks <- OIR_ranks[!is.na(OIR_ranks)]

# 5) Confirm overlap between ranks and the metabolite set
intersect(names(OIR_ranks), aging_metabolite_names)

# 6) Run fgsea
OIR_fgsea_res <- fgsea(
  pathways = list(AgingSignature = aging_metabolite_names),
  stats = OIR_ranks,
  nperm = 1000
)

# 7) View results
print(OIR_fgsea_res)

# 8) Plot enrichment
OIR_msea_plot <-plotEnrichment(
  pathway = aging_metabolite_names,
  stats = OIR_ranks
) + ggtitle("Enrichment of Aging Metabolite Set in OIR vs OS")

# Save enrichment plot
pdf("OIR_msea_plot.pdf", width = 4.25, height = 5.5)
print(OIR_msea_plot)
dev.off()
#--------------------------------------------------------------#
#---- Flip point for annotation -----------------------------#
OIR_df_bottom <- tibble(
  Rank = 1:length(OIR_ranks),
  Metabolite = names(OIR_ranks),
  Stat = OIR_ranks
)

last_pos_rank <- max(which(OIR_df_bottom$Stat > 0))
flip_x <- if (!is.na(last_pos_rank)) last_pos_rank + 0.5 else NA

# Bottom panel plot
oir_vs_os_p_bottom <- ggplot(OIR_df_bottom, aes(x = Rank, y = Stat, fill = Stat)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       name = "Ranked Statistic\n(logFC)") +
  { if (!is.na(flip_x)) geom_vline(xintercept = flip_x, linetype = "dotted", color = "black") else NULL } +
  labs(x = "Metabolite Rank", y = "Statistic") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(oir_vs_os_p_bottom)

# Save bottom panel plot
pdf("oir_vs_os_p_bottom.pdf", width = 4.25, height = 5.5)
print(oir_vs_os_p_bottom)
dev.off()
#------------------------------------------#








#=========================================================#
#----------PCA OS vs YS-----------------------------------#
#=========================================================#

# 1️⃣ Subset to OS and YS groups
pca_data <- imputed_data %>%
  filter(Group %in% c("OS", "YS"))

# 2️⃣ Extract the numeric matrix of metabolites
metabolite_cols <- colnames(pca_data)[3:ncol(pca_data)]
X <- pca_data %>%
  select(all_of(metabolite_cols)) %>%
  as.matrix()

# 3️⃣ Log-transform and scale
# Note: log2 transform to reduce skew, then center+scale
X_log_scaled <- X %>%
  log2() %>%
  scale(center = TRUE, scale = TRUE)

# 4️⃣ Run PCA
pca <- prcomp(X_log_scaled, center = FALSE, scale. = FALSE)

# 5️⃣ Create a data frame of PCA results for plotting
pca_df <- as.data.frame(pca$x) %>%
  bind_cols(pca_data %>% select(Sample, Group))

# 6️⃣ Inspect variance explained
summary(pca)

# 7️⃣ Visualize PC1 vs PC2
library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of Metabolomics Data (OS vs YS)",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)")
  ) +
  scale_color_manual(values = c("YS" = "steelblue", "OS" = "firebrick"))


aging_stats %>%
  filter(abs(log2FC) >1)


sig_aging_metabolites <-aging_stats %>% 
  filter(p_value <0.05 & log2FC >0.3)

aging_metabolite.vec <-sig_aging_metabolites %>%
  pull(Metabolite)

aging_metabolite_names

aging_metabolite_set_effectsize

t_test_results_log

-log10(.01)

#-----------------------------------------------#

#------old pwr vs old sed heatmap---------------#

#----Aging Metabolite Heatmap--------------------#

aging_metabolite_names <- aging_metabolite_set_effectsize %>%
  filter(Log2_FC_OS_vs_YS > 0, P_value < 0.05) %>%
  pull(Metabolite)

length(aging_metabolite_names)

heatmap_data <- imputed_log %>%
  filter(Group %in% c("OV", "OS")) %>%
  dplyr::select(Sample, Group, all_of(aging_metabolite_names))

head(heatmap_data)

# Define desired order of groups
heatmap_data <- heatmap_data %>%
  mutate(Group = factor(Group, levels = c("OV", "OS"))) %>%
  arrange(Group)

# Extract numeric matrix
X <- heatmap_data %>%
  dplyr::select(all_of(aging_metabolite_names)) %>%
  as.matrix()

rownames(X) <- heatmap_data$Sample

# Z-score
X_scaled <- scale(X)

# Transpose
X_scaled_t <- t(X_scaled)

# Annotation for columns (samples)
annotation_col <- data.frame(Group = heatmap_data$Group)
rownames(annotation_col) <- heatmap_data$Sample

# Heatmap
ov_metab_heatmap <- pheatmap(
  mat = X_scaled_t,
  annotation_col = annotation_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,                 # hide row dendrogram
  legend = TRUE,                      # keep heatmap legend
  annotation_legend = FALSE,          # hide group legend
  legend_breaks = c(-1, 0, 1),        # fewer ticks = smaller legend
  legend_labels = c("-1", "0", "+1"),
  border_color = NA,                  # cleaner look
  fontsize = 8,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Aging-Associated Upregulated Metabolites (OS vs YS)",
  row_names_side = "left",
  breaks = seq(-1.2, 1.2, length.out = 101)
)

# Save as PDF (adjust size as needed)
pdf("ov_metab_heatmap.pdf", width = 4.25, height = 5.5)
print(ov_metab_heatmap)
dev.off()
#--------------------------------------------------#

#-----Horizontal Bar Plot OV vs OS---------------------#

library(ggplot2)
library(forcats)
library(dplyr)

t_test_OV_OS
aging_metabolite_names
# Subset stats for aging lipid signature, alphabetical order
ov_metab_logfc_bar_data <- t_test_OV_OS %>%
  filter(Metabolite %in% aging_metabolite_names) %>%
  arrange(Metabolite) %>%
  mutate(Metabolite = factor(Metabolite, levels = Metabolite))  # preserve alphabetical order

ov_metab_logfc_bar_data

# Create horizontal bar plot
ov_metab_logfc_bar_plot <- ggplot(ov_metab_logfc_bar_data, aes(x = Log2_FC_OV_vs_OS, y = fct_rev(Metabolite))) +
  geom_col(fill = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(x = "log2 Fold Change (OS vs YS)", y = NULL) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save as PDF (adjust size as needed)
pdf("ov_metab_logfc_bar_plot.pdf", width = 4.25, height = 5.5)
print(ov_metab_logfc_bar_plot)
dev.off()
#--------------------------------------------------#

#-----Horizontal Bar Plot OS vs YS---------------------#

library(ggplot2)
library(forcats)
library(dplyr)

aging_stats
aging_metabolite_names
# Subset stats for aging lipid signature, alphabetical order
os_metab_logfc_bar_data <- aging_stats %>%
  filter(Metabolite %in% aging_metabolite_names) %>%
  arrange(Metabolite) %>%
  mutate(Metabolite = factor(Metabolite, levels = Metabolite))  # preserve alphabetical order

os_metab_logfc_bar_data

# Create horizontal bar plot
os_metab_logfc_bar_plot <- ggplot(os_metab_logfc_bar_data, aes(x = log2FC, y = fct_rev(Metabolite))) +
  geom_col(fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(x = "log2 Fold Change (OS vs YS)", y = NULL) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save as PDF (adjust size as needed)
pdf("os_metab_logfc_bar_plot.pdf", width = 4.25, height = 5.5)
print(os_metab_logfc_bar_plot)
dev.off()
print(ov_metab_logfc_bar_plot)
dev.new()   # starts a new RStudio graphics device
print(ov_metab_logfc_bar_plot)
#------------------------------------------#

#==========================================#
#------YV vs YS----------------------------#
#==========================================#

#--------------- Compare YV vs YS--------------------------------#

pwr_stats

# 3) Calculate ranked statistic: -log10(p) * logFC
YV_rank_df <- pwr_stats %>%
  filter(!is.na(log2FC), !is.na(p_value)) %>%
  mutate(
    p_safe = pmax(p_value, 1e-300), # avoid Inf
    stat   = -log10(p_safe) * log2FC
  ) %>%
  group_by(Metabolite) %>%
  summarise(stat = mean(stat), .groups = "drop") %>%
  arrange(desc(stat))

# 4) Named numeric vector for fgsea
YV_ranks <- setNames(YV_rank_df$stat, YV_rank_df$Metabolite)

# 5) Run fgsea
YV_fgsea_res <- fgsea(
  pathways = list(AgingSignature = aging_metabolite_names),
  stats    = YV_ranks,
  nperm    = 1000
)

print(YV_fgsea_res)

# 6) Plot enrichment
YV_lsea_plot <- plotEnrichment(
  pathway = aging_metabolite_names,
  stats   = YV_ranks
) + ggtitle("Enrichment of Aging Metabolite Set in YV vs YS")

# Save enrichment plot
pdf("YV_lsea_plot.pdf", width = 4.25, height = 5.5)
print(YV_lsea_plot)
dev.off()

aging_metabolite_names

write.csv(aging_metabolite_names,
          file = "aging_metabolite_names.csv",
          row.names = FALSE)
