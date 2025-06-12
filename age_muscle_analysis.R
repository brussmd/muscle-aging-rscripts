#Aging muscle analysis
#This is the simple analysis of the aging muscle files.
#The Raw Count files are read in using Aging_RawCountFiles.R
#found here: /Users/brussm/Documents/RStudioProjects/Rapa_PwR
#The functions to perform the analysis are from age_df_function.R
#found here: /Users/brussm/Documents/RStudioProjects/Rapa_PwR


#////////////////////////////#
#---Female Gastrocnemius Analysis---#
#////////////////////////////#

#-----1. Run read and process analysis-------#
gas_age_df <- run_full_age_analysis(
  file_list = file_list_gas,
  sample_groups = sample_groups_gas,
  muscle_name = "gas"
)

head(gas_age_df)

gas_age_df %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  arrange(desc(age_responsive))
#-------------------------------#

#------2. Plot from gas_age_df-----------------#
plot_age_logFC_heatmap(gas_age_df, Inflammaging_symbols)

plot_age_logFC_heatmap(gas_age_df, SenMayo_genes)
#------------------------------------#

#------3. Create age-specific dotplot----------#
plot_gene_set_dotplot(gas_age_df, Inflammaging_symbols)
plot_gene_set_dotplot2(gas_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50))
#------------------------------------#

#------4. Test Fisher Exact Test Function--------------#
test_gene_set_enrichment(gas_age_df, Inflammaging_symbols)


head(gas_age_df)
#-----------------------------------#

#////////////////////////////#
#---Female Soleus Analysis---#
#////////////////////////////#

#-----1. Run read and process analysis-------#
sol_age_df <- run_full_age_analysis(
  file_list = file_list_sol,
  sample_groups = sample_groups_sol,
  muscle_name = "sol"
)

head(sol_age_df)

sol_age_df %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  arrange(desc(age_responsive))
#-------------------------------#

#------2. Plot from sol_age_df-----------------#
plot_age_logFC_heatmap(sol_age_df, Inflammaging_symbols)

plot_age_logFC_heatmap(sol_age_df, SenMayo_genes)
#------------------------------------#

#------3. Create age-specific dotplot----------#
plot_gene_set_dotplot(sol_age_df, Inflammaging_symbols)
plot_gene_set_dotplot2(sol_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50))
#------------------------------------#

#------4. Test Fisher Exact Test Function--------------#
test_gene_set_enrichment(sol_age_df, Inflammaging_symbols)


head(sol_age_df)
#------------------------------#

#////////////////////////////#
#---Female Tibialis Anterior Analysis---#
#////////////////////////////#

#-----1. Run read and process analysis-------#
ta_age_df <- run_full_age_analysis(
  file_list = file_list_ta,
  sample_groups = sample_groups_ta,
  muscle_name = "ta"
)

head(ta_age_df)

ta_age_df %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  arrange(desc(age_responsive))
#-------------------------------#

#------2. Plot from sol_age_df-----------------#
plot_age_logFC_heatmap(ta_age_df, Inflammaging_symbols)

plot_age_logFC_heatmap(ta_age_df, SenMayo_genes)
#------------------------------------#

#------3. Create age-specific dotplot----------#
plot_gene_set_dotplot(ta_age_df, Inflammaging_symbols)
plot_gene_set_dotplot2(ta_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50))
#------------------------------------#

#------4. Test Fisher Exact Test Function--------------#
test_gene_set_enrichment(ta_age_df, Inflammaging_symbols)


head(ta_age_df)
#----------------------------------------#

#////////////////////////////#
#---Male Gastrocnemius Analysis---#
#////////////////////////////#
setwd("/Users/brussm/Documents/RStudioProjects/Rapa_PwR/GSE226117_RAW")

#-----1. Run read and process analysis-------#
male_gas_age_df <- run_full_age_analysis(
  file_list = file_list_m_gas,
  sample_groups = sample_groups_m_gas,
  muscle_name = "gas"
)

head(male_gas_age_df)

male_gas_age_df %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  arrange(desc(age_responsive))
#-------------------------------#

#------2. Plot from gas_age_df-----------------#
plot_age_logFC_heatmap(male_gas_age_df, Inflammaging_symbols)

plot_age_logFC_heatmap(male_gas_age_df, SenMayo_genes)
#------------------------------------#

#------3. Create age-specific dotplot----------#
plot_gene_set_dotplot(male_gas_age_df, Inflammaging_symbols)
plot_gene_set_dotplot2(male_gas_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50))
#------------------------------------#

#------4. Test Fisher Exact Test Function--------------#
test_gene_set_enrichment(male_gas_age_df, Inflammaging_symbols)
test_gene_set_enrichment(male_gas_age_df, SenMayo_genes)


head(male_gas_age_df)
#--------------------------------------------#

#////////////////////////////#
#---Male Soleus Analysis---#
#////////////////////////////#

#-----1. Run read and process analysis-------#
male_sol_age_df <- run_full_age_analysis(
  file_list = file_list_m_sol,
  sample_groups = sample_groups_m_sol,
  muscle_name = "sol"
)

head(male_sol_age_df)

male_sol_age_df %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  arrange(desc(age_responsive))
#-------------------------------#

#------2. Plot from gas_age_df-----------------#
plot_age_logFC_heatmap(male_sol_age_df, Inflammaging_symbols)

plot_age_logFC_heatmap(male_sol_age_df, SenMayo_genes)
#------------------------------------#

#------3. Create age-specific dotplot----------#
plot_gene_set_dotplot(male_sol_age_df, Inflammaging_symbols)
plot_gene_set_dotplot2(male_sol_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50))
#------------------------------------#

#------4. Test Fisher Exact Test Function--------------#
test_gene_set_enrichment(male_sol_age_df, Inflammaging_symbols)


head(male_sol_age_df)
#--------------------------------------------#

#////////////////////////////#
#---Male Tibialis Anterior Analysis---#
#////////////////////////////#

#-----1. Run read and process analysis-------#
male_ta_age_df <- run_full_age_analysis(
  file_list = file_list_m_ta,
  sample_groups = sample_groups_m_ta,
  muscle_name = "ta"
)

head(male_ta_age_df)

male_ta_age_df %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  arrange(desc(age_responsive))
#-------------------------------#

#------2. Plot from gas_age_df-----------------#
plot_age_logFC_heatmap(male_ta_age_df, Inflammaging_symbols)

plot_age_logFC_heatmap(male_ta_age_df, SenMayo_genes)
#------------------------------------#

#------3. Create age-specific dotplot----------#
plot_gene_set_dotplot(male_ta_age_df, Inflammaging_symbols)
plot_gene_set_dotplot2(male_ta_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50))
#------------------------------------#

#------4. Test Fisher Exact Test Function--------------#
test_gene_set_enrichment(male_ta_age_df, Inflammaging_symbols)


head(male_ta_age_df)
#--------------------------------------------#


#--------------------------------------#
#/////////////////////////////////////#
#----Messing around with commonality of Inflammaging Genes-----#
#////////////////////////////////////////////////////#

#---------------------------------#
#----Age Responsive Comparison----#

head(gas_age_df)
head(sol_age_df)
head(ta_age_df)

# Step 1: Select only GeneName and age_responsive, renaming age_responsive
gas_resp <- gas_age_df %>%
  select(GeneName, gas_age_responsive = age_responsive)

sol_resp <- sol_age_df %>%
  select(GeneName, sol_age_responsive = age_responsive)

ta_resp <- ta_age_df %>%
  select(GeneName, ta_age_responsive = age_responsive)

male_gas_resp <- male_gas_age_df %>%
  select(GeneName, m_gas_age_responsive = age_responsive)

head(male_gas_resp)

male_sol_resp <- male_sol_age_df %>%
  select(GeneName, m_sol_age_responsive = age_responsive)
 head(male_sol_resp) 

 male_ta_resp <- male_ta_age_df %>%
   select(GeneName, m_ta_age_responsive = age_responsive)
 head(male_ta_resp) 

 tri_resp <- data.frame(
   GeneName = Inflammaging_symbols,
   tri_age_responsive = TRUE,
   stringsAsFactors = FALSE
 )
 tri_resp

 #-----------------------------------------------#
 # Step 2: Join the three dataframes on GeneName
age_responsive_summary <- gas_resp %>%
  full_join(sol_resp, by = "GeneName") %>%
  full_join(ta_resp, by = "GeneName")

# View the result
head(age_responsive_summary)
#-----------------------------------------------#

# Step 2alternate: Join all dataframes on GeneName
full_age_responsive_summary <- gas_resp %>%
  full_join(sol_resp, by = "GeneName") %>%
  full_join(ta_resp, by = "GeneName") %>%
  full_join(male_gas_resp, by = "GeneName") %>%
  full_join(male_sol_resp, by = "GeneName") %>%
  full_join(male_ta_resp, by = "GeneName") %>%
  full_join(tri_resp, by = "GeneName")

# View the result
head(full_age_responsive_summary)
#------------------------------------#
Inflammaging_symbols


age_responsive_summary %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  filter(gas_age_responsive == TRUE & ta_age_responsive == TRUE)

oldsed_vs_yngsed_genes3 %>%
  filter(Symbol %in% Inflammaging_symbols)

oldpwr_vs_yngsed_genes %>%
  filter(Symbol %in% Inflammaging_symbols)

oldpwrveh_v_oldsedveh_genes %>%
  filter(Symbol %in% Inflammaging_symbols) %>%
  filter(logFC < -0.1)

inflammaging_ex_responsive <- oldpwrveh_v_oldsedveh_genes %>%
  filter(Symbol %in% Inflammaging_symbols) %>%
  rename(GeneName = Symbol) %>%
  mutate(ex_responsive = logFC < -0.1) %>%
  select(GeneName, ex_responsive)
inflammaging_ex_responsive %>%
  arrange(ex_responsive)

library(dplyr)
# View result
head(inflammaging_ex_responsive)
inflammaging_ex_responsive %>%
  filter (ex_responsive == TRUE)

dim(age_responsive_summary)

inflammaging_age_responsive <- age_responsive_summary %>%
  filter(GeneName %in% Inflammaging_symbols)

full_inflammaging_age_responsive <- full_age_responsive_summary %>%
  filter(GeneName %in% Inflammaging_symbols)

inflammaging_age_responsive
full_inflammaging_age_responsive
inflammaging_age_responsive %>%
  summarise(count_true = sum(ta_age_responsive, na.rm = TRUE))


inflammaging_age_ex_summary <- inflammaging_age_responsive %>% 
  full_join(inflammaging_ex_responsive, by = "GeneName")

inflammaging_age_ex_summary %>%
  filter(gas_age_responsive == TRUE &
           ta_age_responsive == TRUE &
           ex_responsive == TRUE)

inflammaging_age_ex_summary
write.csv(inflammaging_age_ex_summary, "inflammaging_age_ex_summary.csv", row.names = FALSE)

inflammaging_age_responsive
#---------------------------------------------#
#--------Print Heatmaps-----------------------#
plot_age_logFC_heatmap(gas_age_df, Inflammaging_symbols)

pdf("gas_age_heatmap.pdf", width = 4, height = 6)
plot_age_logFC_heatmap(gas_age_df, Inflammaging_symbols)
dev.off()

pdf("ta_age_heatmap.pdf", width = 4, height = 6)
print(plot_age_logFC_heatmap(ta_age_df, Inflammaging_symbols))
dev.off()

pdf("sol_age_heatmap.pdf", width = 4, height = 6)
print(plot_age_logFC_heatmap(sol_age_df, Inflammaging_symbols))
dev.off()
#-----------------------------------------------#
#-----Print Dotplots----------------------------#

pdf("gas_age_dotplot.pdf", width = 4, height = 2)
print(plot_gene_set_dotplot2(gas_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50)))
dev.off()

pdf("ta_age_dotplot.pdf", width = 4, height = 2)
print(plot_gene_set_dotplot2(ta_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50)))
dev.off()

pdf("sol_age_dotplot.pdf", width = 4, height = 2)
print(plot_gene_set_dotplot2(sol_age_df, Inflammaging_symbols, logfc_limits = c(0, 1), deg_limits = c(0, 50)))
dev.off()
#--------------------------------------------------#
#----Age Responsiveness Comparison Dot Plot--------#

# Assuming your dataframe is called `inflammaging_age_responsive`
# Step 1: Reshape to long format
dot_data <- inflammaging_age_responsive %>%
  pivot_longer(cols = -GeneName, names_to = "Muscle", values_to = "Responsive") %>%
  mutate(
    Status = case_when(
      Responsive == TRUE ~ "TRUE",
      TRUE ~ "FALSE"
    ),
    DotColor = ifelse(Status == "TRUE", "green3", "gray80"),
    DotSize = 3  # fixed size
  )

# Optional: Order genes by total responsiveness (descending)
gene_order <- dot_data %>%
  group_by(GeneName) %>%
  summarise(Total = sum(Status == "TRUE")) %>%
  arrange(desc(Total)) %>%
  pull(GeneName)

# Ensure factor levels reflect order
dot_data <- dot_data %>%
  mutate(GeneName = factor(GeneName, levels = rev(gene_order))) %>%
  mutate(Muscle = factor(Muscle, levels = c("gas_age_responsive", "ta_age_responsive", "sol_age_responsive")))


#-----------------------------------------------------#

age_response_comparison_dotplot <- ggplot(dot_data, aes(x = Muscle, y = GeneName, color = DotColor)) +
  geom_point(aes(size = DotSize), shape = 16) +
  scale_color_identity() +
  scale_size_identity() +
  scale_x_discrete(expand = expansion(mult = c(0.25,10))) +  # remove left spacing
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_blank(),  # remove vertical grid
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),  # optional: remove horizontal grid too
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "Muscle Type",
    y = "Inflammaging Genes",
    title = "Age-Responsiveness of Inflammaging Genes Across Muscle Types"
  )

print(age_response_comparison_dotplot)

ggsave("age_responsive_dotplot.pdf", plot = age_response_comparison_dotplot,
       width = 4, height = 8, units = "in")

inflammaging_age_responsive
#---------------------------------------------#
#-----Full age responsive summary dotplot-----#

# Reshape to long format
dot_data <- full_inflammaging_age_responsive %>%
  pivot_longer(
    cols = -GeneName,
    names_to = "Muscle",
    values_to = "Responsive",
    values_drop_na = FALSE  # <--- keep NA rows
  ) %>%
  mutate(
    Status = case_when(
      Responsive == TRUE ~ "TRUE",
      is.na(Responsive) ~ "NA",
      TRUE ~ "FALSE"
    )
    ,
    DotColor = case_when(
      Status == "TRUE" ~ "green3",
      Status == "FALSE" ~ "gray80",
      TRUE ~ "white"  # optional: color for NA
    ),
    DotSize = ifelse(Status == "TRUE", 3, 2)
  )

# Optional: Order genes by total responsiveness (excluding NA)
gene_order <- dot_data %>%
  filter(Status != "NA") %>%
  group_by(GeneName) %>%
  summarise(Total = sum(Status == "TRUE")) %>%
  arrange(desc(Total)) %>%
  pull(GeneName)

# Set factors for plotting
dot_data <- dot_data %>%
  mutate(
    GeneName = factor(GeneName, levels = rev(gene_order)),
    Muscle = factor(Muscle, levels = c(
      "tri_age_responsive",
      "gas_age_responsive", "ta_age_responsive", "sol_age_responsive",
      "m_gas_age_responsive", "m_ta_age_responsive", "m_sol_age_responsive"
    ))
  )

#-----------------------------------------------------#

age_response_comparison_dotplot <- ggplot(dot_data, aes(x = Muscle, y = GeneName, color = DotColor)) +
  geom_point(aes(size = DotSize), shape = 16) +
  scale_color_identity() +
  scale_size_identity() +
  scale_x_discrete(expand = expansion(mult = c(0.1,1))) +  # remove left spacing
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_blank(),  # remove vertical grid
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),  # optional: remove horizontal grid too
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "Muscle Type",
    y = "Inflammaging Genes",
    title = "Age-Responsiveness of Inflammaging Genes Across Muscle Types"
  )

print(age_response_comparison_dotplot)

ggsave("age_responsive_dotplot.pdf", plot = age_response_comparison_dotplot,
       width = 4, height = 8, units = "in")


