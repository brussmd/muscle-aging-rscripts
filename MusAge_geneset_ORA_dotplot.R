validated_inflamm_geneset.vec
aging_bckgrnd.entrez <- oldsedveh_v_yngsedveh_genes %>%
  pull(ENTREZID)
validated_inflamm_geneset.entrez

go_enrich <- enrichGO(
  gene          = validated_inflamm_geneset.entrez,
  universe      = aging_bckgrnd.entrez,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  minGSSize     = 10,     # adjust if your 52-gene set is small
  maxGSSize     = 500,    # avoid very large, non-specific terms
  readable      = TRUE    # maps result genes back to SYMBOL
)

head(as.data.frame(go_enrich), 20)

## 3) Reduce redundancy (semantic similarity)
go_enrich_simplified <- simplify(
  go_enrich,
  cutoff     = 0.5,      # similarity threshold; 0.4-0.7 common
  by         = "p.adjust",
  select_fun = min,
  measure    = "Wang"
)

go_simpl_df <- as.data.frame(go_enrich_simplified)

# Example display filter: at least 5 genes and FDR < 0.05 (tune to taste)
inflamm_go_subset <- go_simpl_df %>%
  filter(Count >= 9, p.adjust < 0.05) %>%
  arrange(GeneRatio) %>%
  head(20)

inflamm_go_subset
#----------------------------------------------#

#------------DotPlot----------------------#

library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)

# optional: shorten long GO names for plotting
shorten_go <- function(x) x %>%
  str_replace_all("cysteine-type endopeptidase", "caspase") %>%
  str_replace_all("apoptotic process", "apoptosis") %>%
  str_replace_all("\\binvolved in\\b", "in") %>%
  str_replace_all("\\bregulation of\\b", "Reg.") %>%
  str_squish()

df <- inflamm_go_subset %>%
  mutate(
    # turn "13/52" into numeric ratio
    GeneRatio_num = sapply(strsplit(as.character(GeneRatio), "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2])),
    # color by significance as -log10(FDR)
    neglogFDR = -log10(p.adjust),
    # cleaner labels
    Term = str_wrap(shorten_go(Description), width = 36)
  )

# choose what “moves dots to the right”: GeneRatio, RichFactor, or FoldEnrichment
x_metric <- "GeneRatio_num"  # or "RichFactor" or "FoldEnrichment"

df <- df %>% mutate(Term = fct_reorder(Term, !!rlang::sym(x_metric)))

# how much extra room to add on the left (6% of data range)
left_pad <- diff(range(df[[x_metric]], na.rm = TRUE)) * 0.06

p_go <- ggplot(df, aes(x = !!sym(x_metric), y = Term)) +
  geom_point(aes(size = Count, color = neglogFDR), alpha = 0.9) +
  scale_size_area(max_size = 10, name = "Hit genes") +
  # sea blue -> sunset red (with a warm mid)
  scale_color_gradientn(
    colours = c("#2B8CBE", "#73B3D8", "#FEE08B", "#F46D43", "#D7301F"),
    name = expression(-log[10]~FDR)
    # optionally, fix breaks to something nice:
    # breaks = pretty(df$neglogFDR), guide = guide_colorbar(frame.colour = "grey80")
  ) +
  scale_x_continuous(
    name   = if (x_metric == "GeneRatio_num") "Gene ratio (k / 52)" else x_metric,
    labels = if (x_metric == "GeneRatio_num") percent_format(accuracy = 0.1) else label_number(),
    limits = c(min(df[[x_metric]], na.rm = TRUE) - left_pad, NA),  # <- extra room on left
    expand = expansion(mult = c(0, 0.08))                          # no extra mult on left now
  ) +
  coord_cartesian(clip = "off") +  # allow dots to extend into the margin a bit
  labs(y = NULL, title = "GO Biological Process over-representation") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y  = element_text(size = 9),
    legend.box   = "vertical",
    plot.margin  = margin(8, 16, 8, 20)  # a bit more left margin for labels/dots
  )


p_go

# Save as PDF
pdf("p_go.pdf", width = 8, height = 8)
print(p_go)
dev.off()
