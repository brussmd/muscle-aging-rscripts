#This is the code for determining if the Inflammaging Gene Set is an accurate
#reflection of the GO pathways from which it was derived.
#This will pull object names from SOP_code_GSEA_Inflammaging_Geneset.R
#Can be found here: /Users/brussm/Documents/RStudioProjects/Rapa_PwR
#or here: https://github.com/brussmd/muscle-aging-rscripts
#so first run the SOP_code_GSEA_Inflammaging_Geneset.R to initialize the objects

#----Revalidate the Inflammaging KEGG enrichment---#

#Check ENTREZID vectors of all genes from RNA Seq run & for Inflammaging gene set
oldsedveh_v_yngsedveh_12jun.vec  #This is actually a named vector
inflammaging_entrez_vec_12jun

#Run GO enrichment on Inflammaging gene set
go_enrich <- enrichGO(gene = as.character(inflammaging_gene_vec),
                      universe = names(oldsedveh_v_yngsedveh.vec),
                      OrgDb = org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP", # biological process
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

head(as.data.frame(go_enrich), 20)

go_enrich_simplified <- simplify(
  go_enrich,
  cutoff = 0.5,         # semantic similarity threshold (adjustable)
  by = "p.adjust",      # keeps the most significant term
  select_fun = min,
  measure = "Wang"      # semantic similarity method
)

head(as.data.frame(go_enrich_simplified), 20)

go_simpl_df <- as.data.frame(go_enrich_simplified)

inflamm_go_subset <- head(go_simpl_df %>%
                            filter(Count >10 & p.adjust <0.05),10)

inflamm_go_subset

# Convert to data frame
go_df <- as.data.frame(go_enrich_simplified)

# Step 1: Split geneID column into gene lists
gene_lists <- strsplit(inflamm_go_subset$geneID, split = "/")
names(gene_lists) <- inflamm_go_subset$Description

# Step 2: Create a long-format presence dataframe (GeneName, GO_Term, Present = 1)
presence_df <- map2_dfr(
  gene_lists,
  names(gene_lists),
  ~ tibble(GeneName = .x, GO_Term = .y, Present = "1")
)

# Step 3: Create full combinations of all Inflammaging genes × GO terms
all_combinations <- expand_grid(
  GeneName = Inflammaging_symbols,
  GO_Term = inflamm_go_subset$Description
)

# Step 4: Merge with presence info and fill missing values with "0"
presence_df <- left_join(all_combinations, presence_df, by = c("GeneName", "GO_Term")) %>%
  mutate(Present = replace_na(Present, "0"))

# Step 5: Order genes by total GO term count (descending)
gene_order <- presence_df %>%
  group_by(GeneName) %>%
  summarise(Total = sum(as.integer(Present))) %>%
  arrange(desc(Total)) %>%
  pull(GeneName)

# Step 6: Order GO terms by total gene count (optional)
go_order <- presence_df %>%
  group_by(GO_Term) %>%
  summarise(Total = sum(as.integer(Present))) %>%
  arrange(desc(Total)) %>%
  pull(GO_Term)

# Step 7: Convert to factors with specified order
presence_df <- presence_df %>%
  mutate(
    GeneName = factor(GeneName, levels = rev(gene_order)),  # top = most represented
    GO_Term = factor(GO_Term, levels = go_order)
  )

# Step 8: Plot using geom_tile
inflamm_kegg_validation_plot <- ggplot(presence_df, aes(x = GO_Term, y = GeneName, fill = Present)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = c("0" = "white", "1" = "green3"), guide = "none") +
  geom_hline(yintercept = seq(0.5, length(gene_order) + 0.5, 1), color = "gray90") +
  geom_vline(xintercept = seq(0.5, length(go_order) + 0.5, 1), color = "gray90") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "GO Term",
    y = "Inflammaging Gene",
    title = "Inflammaging Gene Presence in Top GO Terms"
  )

# Save as PDF
ggsave("inflammaging_go_dotplot.pdf", plot = inflamm_kegg_validation_plot, width = 8, height = 8, units = "in")

print(inflamm_kegg_validation_plot)
