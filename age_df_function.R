
#---Create Aging df Function--------------#
#/////////////////////////////////////////#
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(edgeR)

run_full_age_analysis <- function(
    file_list, sample_groups, muscle_name
) {
  # Helper: Extract sample ID
  extract_id <- function(file) str_extract(file, "X\\d+")
  
  # Build count matrix
  all_ids <- unlist(sample_groups)
  first_file <- fread(file_list[1])
  genes <- first_file$GeneName
  count_list <- lapply(file_list, function(f) {
    sample_id <- extract_id(f)
    df <- fread(f)[, .(GeneName, count = get(names(.SD)[2]))]
    df$count[is.na(df$count)] <- 0  # Replace NA with 0
    setNames(df$count, df$GeneName)
  })
  count_mat <- do.call(cbind, count_list)
  colnames(count_mat) <- sapply(file_list, extract_id)
  rownames(count_mat) <- genes
  count_mat <- as.matrix(count_mat)
  
  # Run EdgeR comparisons vs 6mo
  results_list <- list()
  for (month in names(sample_groups)) {
    if (month == "6") next
    g1 <- sample_groups[["6"]]
    g2 <- sample_groups[[month]]
    group <- factor(c(rep("6mo", length(g1)), rep(paste0(month, "mo"), length(g2))), 
                    levels = c("6mo", paste0(month, "mo")))
    dge <- DGEList(counts = count_mat[, c(g1, g2)], group = group)
    dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    design <- model.matrix(~group)
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, coef = 2)
    deg <- topTags(lrt, n = Inf)$table
    df <- data.frame(
      GeneName = rownames(deg),
      logFC = deg$logFC,
      FDR = deg$FDR
    )
    colnames(df)[2:3] <- c(paste0("logFC_", month), paste0("FDR_", month))
    results_list[[month]] <- df
  }
  
  # Merge all results
  merged_df <- results_list[["12"]]
  for (month in c("18", "21", "24", "27")) {
    merged_df <- left_join(merged_df, results_list[[month]], by = "GeneName")
  }
  merged_df <- drop_na(merged_df)
  
  # Age responsiveness (DEG in ≥2 oldest timepoints)
  merged_df <- merged_df %>%
    mutate(
      sig_21 = FDR_21 < 0.05,
      sig_24 = FDR_24 < 0.05,
      sig_27 = FDR_27 < 0.05,
      old_sig_count = sig_21 + sig_24 + sig_27,
      DEG_old = old_sig_count >= 2
    )
  
  # Spearman trend test
  ages <- c(12, 18, 21, 24, 27)
  logFC_mat <- merged_df %>% dplyr::select(starts_with("logFC_"))
  spearman <- t(apply(logFC_mat, 1, function(x) {
    r <- cor.test(x, ages, method = "spearman")
    c(Rho = r$estimate, Pval = r$p.value)
  }))
  spearman_df <- as.data.frame(spearman)
  colnames(spearman_df) <- c("SpearmanRho", "SpearmanPval")
  merged_df <- bind_cols(merged_df, spearman_df)
  
  # Final classification
  merged_df <- merged_df %>%
    mutate(age_responsive = DEG_old | (SpearmanRho > 0.6 & SpearmanPval < 0.1))
  
  # Save to environment
  assign(paste0("full_", muscle_name, "_merged_df"), merged_df, envir = .GlobalEnv)
  assign(paste0(muscle_name, "_age_df"), merged_df, envir = .GlobalEnv)
  return(merged_df)
}
#-----------------------------------------------------#
#////////////////////////////////////////////////////#

#----------LogFC by age heatmap function---------------#
#//////////////////////////////////////////////////////#
plot_age_logFC_heatmap <- function(age_df, gene_set, scale_rows = TRUE, cluster_rows = FALSE, fontsize_row = 6) {
  library(dplyr)
  library(pheatmap)
  library(tibble)
  
  # Title from variable names
  df_name <- deparse(substitute(age_df))
  gene_set_name <- deparse(substitute(gene_set))
  heatmap_title <- paste0(df_name, ": ", gene_set_name)
  
  # Extract logFC columns
  logFC_cols <- grep("^logFC_", names(age_df), value = TRUE)
  
  # Data for matched genes
  plot_df <- age_df %>%
    filter(GeneName %in% gene_set) %>%
    dplyr::select(GeneName, all_of(logFC_cols)) %>%
    column_to_rownames("GeneName")
  
  # Initialize matrix for all input genes, alphabetical order
  all_genes <- sort(gene_set)
  full_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(logFC_cols),
                        dimnames = list(all_genes, logFC_cols))
  
  # Fill with available data
  matched_genes <- intersect(rownames(plot_df), all_genes)
  full_matrix[matched_genes, ] <- as.matrix(plot_df[matched_genes, ])
  
  # Remove rows with all NAs *if* clustering is enabled (to avoid hclust error)
  if (cluster_rows) {
    non_empty_rows <- apply(full_matrix, 1, function(x) !all(is.na(x)))
    full_matrix <- full_matrix[non_empty_rows, ]
  }
  
  # Define palette
  my_palette <- colorRampPalette(c("blue", "white", "red"))(50)
  na_color <- "grey30"
  
  # Plot
  pheatmap(
    full_matrix,
    main = heatmap_title,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    scale = if (scale_rows) "row" else "none",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = fontsize_row,
    cellwidth = 20,
    cellheight = 6,
    color = my_palette,
    na_col = na_color
  )
}
#-----------------------------------------------------#
#////////////////////////////////////////////////////#

#----------Dotplot by age function---------------#
#//////////////////////////////////////////////////////#
plot_gene_set_dotplot <- function(age_df, gene_set) {
  library(dplyr)
  library(ggplot2)
  
  # Extract timepoints from column names
  logFC_cols <- grep("^logFC_", names(age_df), value = TRUE)
  timepoints <- gsub("logFC_", "", logFC_cols)
  
  # Build summary data frame
  summary_df <- lapply(timepoints, function(tp) {
    logfc_col <- paste0("logFC_", tp)
    fdr_col <- paste0("FDR_", tp)
    df <- age_df %>%
      filter(GeneName %in% gene_set) %>%
      mutate(
        logFC = .data[[logfc_col]],
        FDR = .data[[fdr_col]],
        is_sig = FDR < 0.05
      )
    data.frame(
      Age = tp,
      Avg_logFC = mean(df$logFC, na.rm = TRUE),
      SigGeneCount = sum(df$is_sig, na.rm = TRUE)
    )
  }) %>% bind_rows()
  
  # Ensure proper ordering of timepoints
  summary_df$Age <- factor(summary_df$Age, levels = timepoints)
  
  # Create dotplot
  ggplot(summary_df, aes(x = Age, y = 1)) +
    geom_point(aes(size = SigGeneCount, color = Avg_logFC)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size(range = c(3, 12)) +
    scale_y_continuous(expand = c(0.2, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 20, 5, 20)
    ) +
    labs(
      color = "Avg logFC",
      size = "# DEGs (FDR < 0.05)"
    )
}
#---------------------------------#
#//////////////////////////////////////#

#---------------------------------------#
#-----Modified dotplot function---------#
plot_gene_set_dotplot2 <- function(age_df, gene_set,
                                  logfc_limits = c(-1, 1),
                                  deg_limits = c(0, 100)) {
  library(dplyr)
  library(ggplot2)
  
  # Extract timepoints from column names
  logFC_cols <- grep("^logFC_", names(age_df), value = TRUE)
  timepoints <- gsub("logFC_", "", logFC_cols)
  
  # Build summary data frame
  summary_df <- lapply(timepoints, function(tp) {
    logfc_col <- paste0("logFC_", tp)
    fdr_col <- paste0("FDR_", tp)
    df <- age_df %>%
      filter(GeneName %in% gene_set) %>%
      mutate(
        logFC = .data[[logfc_col]],
        FDR = .data[[fdr_col]],
        is_sig = FDR < 0.05
      )
    data.frame(
      Age = tp,
      Avg_logFC = mean(df$logFC, na.rm = TRUE),
      SigGeneCount = sum(df$is_sig, na.rm = TRUE)
    )
  }) %>% bind_rows()
  
  # Ensure consistent ordering of timepoints
  summary_df$Age <- factor(summary_df$Age, levels = timepoints)
  
  # Plot
  ggplot(summary_df, aes(x = Age, y = 1)) +
    geom_point(aes(size = SigGeneCount, color = Avg_logFC)) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = logfc_limits, oob = scales::squish
    ) +
    scale_size(
      range = c(3, 12),
      limits = deg_limits
    )  +
    scale_y_continuous(expand = c(0.2, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 20, 5, 20)
    ) +
    labs(
      color = "Avg logFC",
      size = "# DEGs (FDR < 0.05)"
    )
}


#--Fisher Exact Test Gene Set age_responsive function--#
#//////////////////////////////////////////////////////#
test_gene_set_enrichment <- function(df, gene_set, gene_column = "GeneName", age_resp_column = "age_responsive") {
  library(dplyr)
  library(tidyr)
  
  # Build contingency table
  fisher_df <- df %>%
    mutate(InGeneSet = .data[[gene_column]] %in% gene_set) %>%
    group_by(InGeneSet, .data[[age_resp_column]]) %>%
    summarise(N = n(), .groups = "drop") %>%
    pivot_wider(names_from = {{ age_resp_column }}, values_from = N, values_fill = 0) %>%
    rename(Not_AgeResp = `FALSE`, AgeResp = `TRUE`)
  
  # Build matrix with correct row order: InGeneSet first
  fisher_mat <- fisher_df %>%
    mutate(RowName = ifelse(InGeneSet, "InGeneSet", "Other")) %>%
    arrange(desc(InGeneSet)) %>%  # InGeneSet should be first
    dplyr::select(RowName, AgeResp, Not_AgeResp) %>%
    column_to_rownames("RowName") %>%
    as.matrix()
  
  # Run Fisher's test
  fisher_result <- fisher.test(fisher_mat)
  
  list(
    contingency_table = fisher_mat,
    p_value = fisher_result$p.value,
    odds_ratio = fisher_result$estimate,
    test_result = fisher_result
  )
}
#---------------------------------------#
#///////////////////////////////////////#


#-----------------------------------------------------#
#---Test Function for female gastroc "gas"------------#
#-----------------------------------------------------#

#-----------------------------------#



#-----1. Run read and process analysis-------#
gas_age_df <- run_full_age_analysis(
  file_list = file_list_gas,
  sample_groups = sample_groups,
  muscle_name = "gas"
)

head(gas_age_df)

gas_age_df %>%
  filter(GeneName %in% Inflammaging_symbols) %>%
  arrange(desc(age_responsive))

gas_age_df %>%
  filter(GeneName %in% SenMayo_genes) %>%
  arrange(desc(age_responsive))
#-------------------------------#

#------2. Plot from gas_age_df-----------------#
plot_age_logFC_heatmap(gas_age_df, Inflammaging_symbols)

plot_age_logFC_heatmap(gas_age_df, SenMayo_genes)
#------------------------------------#

#------3. Create age-specific dotplot----------#
plot_gene_set_dotplot(gas_age_df, Inflammaging_symbols)
#------------------------------------#

#------4. Test Fisher Exact Test Function--------------#
test_gene_set_enrichment(gas_age_df, Inflammaging_symbols)


head(gas_age_df)

gas_age_df %>%
  filter(GeneName == "Il6")
