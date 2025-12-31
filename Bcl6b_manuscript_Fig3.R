#Bcl6b Manuscript Figure 3.


#==========================================================#
#//////Get Promoter Targets////////////////////////////////#
#==========================================================#
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(rtracklayer)
  library(dplyr)
  library(readr)
})

txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

# -----------------------------
# Helper: get promoter targets from a BED/narrowPeak
# -----------------------------
get_promoter_targets <- function(bedfile,
                                 tss_up = 3000,
                                 tss_down = 3000) {
  peaks <- rtracklayer::import(bedfile, format = "narrowPeak")
  
  anno <- annotatePeak(
    peaks,
    TxDb = txdb_hg38,
    tssRegion = c(-tss_up, tss_down),
    annoDb = "org.Hs.eg.db"
  )
  
  df <- as.data.frame(anno)
  
  promoter_genes <- df %>%
    filter(!is.na(SYMBOL)) %>%
    filter(grepl("Promoter", annotation, ignore.case = TRUE)) %>%
    pull(SYMBOL) %>%
    unique() %>%
    toupper()
  
  list(
    promoter_genes = promoter_genes,
    anno_df = df
  )
}

# -----------------------------
# Run both files
# -----------------------------
res_hnp <- get_promoter_targets("ENCFF931HNP.bed")
res_bag <- get_promoter_targets("ENCFF165BAG.bed")

hnp_genes <- res_hnp$promoter_genes
bag_genes <- res_bag$promoter_genes

cat("HNP promoter genes:", length(hnp_genes), "\n")
cat("BAG promoter genes:", length(bag_genes), "\n")
cat("Union promoter genes:", length(union(hnp_genes, bag_genes)), "\n")
cat("Intersection (both):", length(intersect(hnp_genes, bag_genes)), "\n")

# -----------------------------
# "Full join" presence table (gene-level)
# -----------------------------
promoters_presence <- full_join(
  tibble(gene = hnp_genes, in_HNP = TRUE),
  tibble(gene = bag_genes, in_BAG = TRUE),
  by = "gene"
) %>%
  mutate(
    in_HNP = ifelse(is.na(in_HNP), FALSE, in_HNP),
    in_BAG = ifelse(is.na(in_BAG), FALSE, in_BAG),
    source = case_when(
      in_HNP & in_BAG ~ "both",
      in_HNP ~ "HNP_only",
      in_BAG ~ "BAG_only",
      TRUE ~ "none"
    )
  ) %>%
  arrange(desc(in_HNP & in_BAG), desc(in_BAG), desc(in_HNP), gene)

# The "all possible promoters" set you asked for:
promoters_union <- promoters_presence %>%
  filter(in_HNP | in_BAG) %>%
  pull(gene) %>%
  unique()

# -----------------------------
# Save outputs
# -----------------------------
write_csv(promoters_presence, "BCL6B_PromoterTargets_HNP_vs_BAG_fulljoin.csv")
write_csv(tibble(gene = promoters_union), "BCL6B_PromoterTargets_UNION.csv")

# Optional: save annotated peak tables too (big-ish)
write_csv(res_hnp$anno_df, "BCL6B_HNP_All_Annotated_Peaks.csv")
write_csv(res_bag$anno_df, "BCL6B_BAG_All_Annotated_Peaks.csv")

# Quick peek
print(head(promoters_presence, 20))
dim(promoters_presence)
head(promoters_presence)
promoters_presence %>%
  filter(gene == "VWF")

promoters_union
#=======================================================================#

###############################################################################
# BCL6B motif PWM from ENCODE peaks (GADEM) +
# Scan UNION promoter targets (HNP ∪ BAG) in chunks +
# Save hit-level + gene-level motif ranking (14–17 “good”)
###############################################################################

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(TFBSTools)
  library(rGADEM)
  library(dplyr)
  library(tidyr)
  library(readr)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
primary_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")

# ----------------------------
# USER SETTINGS
# ----------------------------

# Peak files (narrowPeak)
PEAK_HNP <- "ENCFF931HNP.bed"
PEAK_BAG <- "ENCFF165BAG.bed"

# Learn motif from: "HNP", "BAG", or "BOTH"
MOTIF_SOURCE <- "HNP"   # change to "BAG" or "BOTH"

# Promoter window (same as your annotatePeak)
UPSTREAM <- 3000
DOWNSTREAM <- 3000

# GADEM summit window used to learn motif (±100 = 200 bp)
SUMMIT_HALF_WIN <- 100

# Motif scan threshold for TFBSTools::searchSeq()
MIN_SCORE <- "80%"      # keep as you requested

# Chunking for promoter scanning
CHUNK_SIZE <- 200       # 200 is a good default

# “Good score” bins
GOOD_LO <- 14
GOOD_HI <- 17

# Output prefix
OUT_PREFIX <- paste0(
  "BCL6B_GADEM_", MOTIF_SOURCE,
  "_", gsub("%","",MIN_SCORE), "min",
  "_prom", UPSTREAM, "to", DOWNSTREAM
)

# ----------------------------
# 1) Get UNION promoter gene list
#    - Prefer in-memory object promoters_union if it exists
#    - Else fall back to CSV you wrote earlier
# ----------------------------
if (exists("promoters_union", inherits = FALSE)) {
  gene_vec <- promoters_union %>% as.character() %>% toupper() %>% unique()
} else if (file.exists("BCL6B_PromoterTargets_UNION.csv")) {
  gene_vec <- readr::read_csv("BCL6B_PromoterTargets_UNION.csv", show_col_types = FALSE) %>%
    dplyr::pull(1) %>% as.character() %>% toupper() %>% unique()
} else {
  stop("Couldn't find promoters_union object or BCL6B_PromoterTargets_UNION.csv")
}
message("Genes loaded for promoter scan (UNION): ", length(gene_vec))

# ----------------------------
# 2) Load peaks (HNP/BAG/BOTH) for motif discovery
# ----------------------------
stopifnot(file.exists(PEAK_HNP), file.exists(PEAK_BAG))

peaks_hnp <- rtracklayer::import(PEAK_HNP, format = "narrowPeak")
peaks_bag <- rtracklayer::import(PEAK_BAG, format = "narrowPeak")

peaks_for_motif <- switch(
  toupper(MOTIF_SOURCE),
  "HNP"  = peaks_hnp,
  "BAG"  = peaks_bag,
  "BOTH" = c(peaks_hnp, peaks_bag),
  stop("MOTIF_SOURCE must be 'HNP', 'BAG', or 'BOTH'")
)

message("Peaks used for motif discovery: ", length(peaks_for_motif), " (", MOTIF_SOURCE, ")")

# ----------------------------
# 3) Build GADEM motif from summit-centered windows
# ----------------------------
# narrowPeak 'peak' column = summit offset from start
if (!("peak" %in% colnames(mcols(peaks_for_motif)))) {
  stop("Your imported narrowPeak object doesn't have a 'peak' column in mcols().")
}

summits <- start(peaks_for_motif) + mcols(peaks_for_motif)$peak

peak_win <- GRanges(
  seqnames = seqnames(peaks_for_motif),
  ranges   = IRanges(start = summits - SUMMIT_HALF_WIN,
                     end   = summits + (SUMMIT_HALF_WIN - 1)),
  strand   = "*"
)

# Keep only primary chromosomes (BSgenome-safe)
peak_win <- peak_win[seqnames(peak_win) %in% primary_chrs]

# Extract sequences directly (faster than FASTA round-trip)
seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, peak_win)
names(seqs) <- paste0("peak_", seq_along(seqs))

message("Peak windows used for GADEM: ", length(seqs))

# GADEM (this is typically the slowest step)
gadem_output <- rGADEM::GADEM(seqs, genome = BSgenome.Hsapiens.UCSC.hg38)

motifs <- rGADEM::getPWM(gadem_output)
stopifnot(length(motifs) >= 1)

motif1 <- motifs[[1]]

# Convert PWM safely -> PFMatrix -> PWMatrix
motif1_counts <- round(motif1 * 1000)
motif1_counts[motif1_counts == 0] <- 1

pfm <- PFMatrix(
  ID = paste0("BCL6B_GADEM_", MOTIF_SOURCE),
  name = paste0("BCL6B_GADEM_", MOTIF_SOURCE),
  profileMatrix = motif1_counts
)

motif_pwm <- toPWM(pfm)

# Save PWM so you never rebuild it
pwm_rds <- paste0(OUT_PREFIX, "_motif_pwm.rds")
saveRDS(motif_pwm, pwm_rds)
message("Saved PWM RDS: ", pwm_rds)

# ----------------------------
# 4) Helper: get promoter GRanges for genes (primary chr only)
# ----------------------------
get_promoters_for_genes <- function(gene_symbols, txdb, upstream, downstream) {
  
  gene_symbols <- unique(toupper(gene_symbols))
  
  gene_ids <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_symbols,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  ) %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  if (nrow(gene_ids) == 0) return(NULL)
  
  tx_map <- AnnotationDbi::select(
    txdb,
    keys = gene_ids$ENTREZID,
    keytype = "GENEID",
    columns = "TXID"
  ) %>%
    dplyr::filter(!is.na(TXID))
  
  if (nrow(tx_map) == 0) return(NULL)
  
  promoters_all <- GenomicFeatures::promoters(txdb, upstream = upstream, downstream = downstream)
  promoters_all <- promoters_all[seqnames(promoters_all) %in% primary_chrs]
  
  prom_sel <- promoters_all[mcols(promoters_all)$tx_id %in% tx_map$TXID]
  if (length(prom_sel) == 0) return(NULL)
  
  # map tx_id -> gene symbol
  gene_for_tx <- gene_ids$SYMBOL[match(tx_map$GENEID, gene_ids$ENTREZID)]
  symbol_for_prom <- gene_for_tx[match(mcols(prom_sel)$tx_id, tx_map$TXID)]
  names(prom_sel) <- symbol_for_prom
  
  prom_sel
}

# ----------------------------
# 5) Helper: scan promoter sequences for motif hits
# ----------------------------
scan_promoters_with_motif <- function(prom_gr, pwm, min.score) {
  
  if (is.null(prom_gr) || length(prom_gr) == 0) return(NULL)
  
  prom_seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, prom_gr)
  names(prom_seqs) <- paste0(names(prom_gr), "_", mcols(prom_gr)$tx_name)
  
  hits_list <- lapply(seq_along(prom_seqs), function(i) {
    TFBSTools::searchSeq(
      x = pwm,
      subject = prom_seqs[[i]],
      min.score = min.score,
      strand = "*"
    )
  })
  names(hits_list) <- names(prom_seqs)
  
  hits_df <- do.call(rbind, lapply(names(hits_list), function(nm) {
    x <- hits_list[[nm]]
    if (length(x) == 0) return(NULL)
    
    data.frame(
      gene       = sub("_ENST.*", "", nm),
      transcript = sub(".*_", "", nm),
      start      = start(x),
      end        = end(x),
      width      = end(x) - start(x) + 1,
      strand     = as.character(strand(x)),
      score      = score(x),
      stringsAsFactors = FALSE
    )
  }))
  
  if (!is.null(hits_df) && nrow(hits_df) > 0) {
    hits_df$gene <- toupper(hits_df$gene)
    hits_df <- hits_df[order(-hits_df$score), ]
  }
  
  hits_df
}

# ----------------------------
# 6) Chunked promoter scan across UNION genes
# ----------------------------
idx <- split(seq_along(gene_vec), ceiling(seq_along(gene_vec) / CHUNK_SIZE))
hits_all <- vector("list", length(idx))

for (k in seq_along(idx)) {
  
  genes_k <- gene_vec[idx[[k]]]
  
  prom_k <- get_promoters_for_genes(
    gene_symbols = genes_k,
    txdb = txdb,
    upstream = UPSTREAM,
    downstream = DOWNSTREAM
  )
  
  hits_k <- scan_promoters_with_motif(prom_k, pwm = motif_pwm, min.score = MIN_SCORE)
  hits_all[[k]] <- hits_k
  
  message(
    "Chunk ", k, "/", length(idx),
    " | genes=", length(genes_k),
    " | promoters=", ifelse(is.null(prom_k), 0, length(prom_k)),
    " | hits=", ifelse(is.null(hits_k), 0, nrow(hits_k))
  )
}

motif_hits <- dplyr::bind_rows(hits_all)

message("TOTAL HIT ROWS: ", nrow(motif_hits))
message("TOTAL UNIQUE GENES WITH >=1 HIT: ", dplyr::n_distinct(motif_hits$gene))

# ----------------------------
# 7) Save hit-level results
# ----------------------------
hits_file <- paste0(OUT_PREFIX, "_MotifHits_hitLevel.csv")
readr::write_csv(motif_hits, hits_file)
message("Wrote: ", hits_file)

# ----------------------------
# 8) Gene-level motif ranking df (14–17 = good)
# ----------------------------
motif_rank_df <- motif_hits %>%
  group_by(gene) %>%
  summarise(
    n_hits          = n(),
    max_score       = max(score, na.rm = TRUE),
    mean_score      = mean(score, na.rm = TRUE),
    best_transcript = transcript[which.max(score)],
    best_start      = start[which.max(score)],
    best_end        = end[which.max(score)],
    .groups = "drop"
  ) %>%
  mutate(
    motif_tier = case_when(
      max_score >= GOOD_HI ~ paste0("A (very strong; ≥", GOOD_HI, ")"),
      max_score >= GOOD_LO ~ paste0("B (good; ", GOOD_LO, "–", GOOD_HI - 0.01, ")"),
      max_score >= (GOOD_LO - 2) ~ paste0("C (moderate; ", GOOD_LO - 2, "–", GOOD_LO - 0.01, ")"),
      TRUE ~ paste0("D (weak; <", GOOD_LO - 2, ")")
    ),
    motif_rank_score = max_score + 0.10 * log1p(n_hits)
  ) %>%
  arrange(desc(motif_rank_score), desc(max_score), desc(n_hits))

rank_file <- paste0(OUT_PREFIX, "_MotifRank_geneLevel.csv")
readr::write_csv(motif_rank_df, rank_file)
message("Wrote: ", rank_file)

print(dplyr::slice_head(motif_rank_df, n = 25))
motif_rank_df %>%
  filter(max_score > 17)
head(motif_rank_df)

###############################################################################
#////////////////////////////////////////////////////////////////////////////
###############################################################################

#===================================================================#
#/////Figure 3a. Promoter target motif landscape///////////////////#
#===================================================================#

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- Build anchor df (4 categories) ----
anchor_df <- tibble(GeneSymbol = toupper(promoters_union)) %>%
  left_join(
    motif_rank_df %>%
      transmute(
        GeneSymbol = toupper(gene),
        max_score = max_score
      ),
    by = "GeneSymbol"
  ) %>%
  mutate(
    max_score = ifelse(is.na(max_score), 0, max_score),
    motif_class = case_when(
      max_score >= 17 ~ "≥17 (very strong)",
      max_score >= 14 ~ "14–17 (strong)",
      max_score >= 10 ~ "10–14 (moderate)",
      TRUE            ~ "<10 (weak/none)"
    ),
    # lock legend order
    motif_class = factor(
      motif_class,
      levels = c("≥17 (very strong)", "14–17 (strong)", "10–14 (moderate)", "<10 (weak/none)")
    )
  ) %>%
  arrange(desc(max_score), GeneSymbol) %>%
  mutate(rank = row_number())

# ---- Summary numbers ----
N_total <- nrow(anchor_df)
N_ge10  <- sum(anchor_df$max_score >= 10)
N_ge14  <- sum(anchor_df$max_score >= 14)
N_ge17  <- sum(anchor_df$max_score >= 17)

# ---- Plot: thin bars ----
ggplot(anchor_df, aes(x = rank, y = max_score, fill = motif_class)) +
  geom_col(
    width = 0.9,
    color = NA,
    alpha = 0.95
  ) +
  geom_hline(
    yintercept = c(10, 14, 17),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c(
      "≥17 (very strong)" = "firebrick3",
      "14–17 (strong)"    = "steelblue3",
      "10–14 (moderate)"  = "skyblue2",
      "<10 (weak/none)"   = "grey80"
    ),
    name = "BCL6B motif\nstrength"
  ) +
  annotate(
    "text",
    x = floor(N_total * 0.65),
    y = max(anchor_df$max_score) * 0.98,
    hjust = 0,
    vjust = 1,
    size = 3.8,
    label = paste0(
      "ENCODE promoter targets (±3 kb): ", comma(N_total), "\n",
      "Genes with motif ≥10: ", comma(N_ge10), "\n",
      "Genes with motif ≥14: ", comma(N_ge14), "\n",
      "Genes with motif ≥17: ", comma(N_ge17)
    )
  ) +
  labs(
    x = "ENCODE promoter target genes (ranked by motif max_score)",
    y = "BCL6B motif strength (max PWM score)",
    title = "ENCODE promoter landscape of BCL6B motif strength"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



#=====================================================================#
#/////Figure 3b. Compare ENCODE PWN to UniProbe PWN///////////////////////#
#========================================================================#
suppressPackageStartupMessages({
  library(TFBSTools)
  library(ggseqlogo)
  library(gridExtra)
})

# ----------------------------
# 1) UniPROBE matrix (probabilities)
# ----------------------------
uniprobe_mat <- rbind(
  A = c(0.34655,0.16933,0.08562,0.12876,0.03238,0.01318,0.01840,0.80090,0.14117,0.07347,0.87850,0.87592,0.20568,0.24434,0.20897,0.20161),
  C = c(0.08230,0.59918,0.10582,0.08468,0.00992,0.86865,0.10450,0.01091,0.05540,0.02124,0.01503,0.02701,0.21331,0.22828,0.30444,0.29408),
  G = c(0.20298,0.04651,0.05939,0.07678,0.01648,0.00940,0.45047,0.01598,0.76238,0.81163,0.04496,0.01378,0.04094,0.16823,0.14847,0.33489),
  T = c(0.36817,0.18499,0.74917,0.70978,0.94122,0.10877,0.42663,0.17220,0.04106,0.09366,0.06151,0.08328,0.54007,0.35916,0.33812,0.16942)
)

# Normalize columns (safe)
uniprobe_mat <- apply(uniprobe_mat, 2, function(x) x / sum(x))
uniprobe_mat <- as.matrix(uniprobe_mat)
rownames(uniprobe_mat) <- c("A","C","G","T")

# ----------------------------
# 2) Your GADEM motif PWM -> probability matrix
# ----------------------------
stopifnot(exists("motif_pwm", inherits = FALSE))
my_mat_log  <- motif_pwm@profileMatrix
my_mat_prob <- apply(2^(my_mat_log), 2, function(x) x / sum(x))
my_mat_prob <- as.matrix(my_mat_prob)
my_mat_prob <- my_mat_prob[c("A","C","G","T"), , drop = FALSE]

# ----------------------------
# 3) Reverse-complement helper (for PWM prob matrices)
# ----------------------------
rc_pwm <- function(prob_mat) {
  rev_mat <- prob_mat[, ncol(prob_mat):1, drop = FALSE]     # reverse columns
  # swap A<->T and C<->G
  out <- rev_mat[c("T","G","C","A"), , drop = FALSE]
  rownames(out) <- c("A","C","G","T")
  out
}

# ----------------------------
# 4) Compare function (slide short within long)
# ----------------------------
compare_pwm <- function(short_mat, long_mat) {
  Ls <- ncol(short_mat)
  Ll <- ncol(long_mat)
  if (Ls > Ll) stop("short_mat is longer than long_mat")
  
  best_cor <- -Inf
  best_pos <- NA_integer_
  
  for (i in 1:(Ll - Ls + 1)) {
    window <- long_mat[, i:(i + Ls - 1), drop = FALSE]
    cor_val <- suppressWarnings(cor(as.vector(short_mat), as.vector(window)))
    if (!is.na(cor_val) && cor_val > best_cor) {
      best_cor <- cor_val
      best_pos <- i
    }
  }
  list(best_cor = best_cor, best_pos = best_pos)
}

# Test forward + reverse-complement
res_fwd <- compare_pwm(my_mat_prob, uniprobe_mat)
res_rc  <- compare_pwm(rc_pwm(my_mat_prob), uniprobe_mat)

if (res_fwd$best_cor >= res_rc$best_cor) {
  best_orient <- "forward"
  best_res <- res_fwd
  my_best <- my_mat_prob
} else {
  best_orient <- "reverse-complement"
  best_res <- res_rc
  my_best <- rc_pwm(my_mat_prob)
}

print(list(best_orientation = best_orient,
           best_correlation = best_res$best_cor,
           best_start_pos_in_uniprobe = best_res$best_pos))

# ----------------------------
# 5) Align + pad for stacked logo plot
#    Use NA padding so ggseqlogo doesn't draw letters there
# ----------------------------
align_pwms <- function(short_mat, long_mat, start_pos) {
  Ls <- ncol(short_mat)
  Ll <- ncol(long_mat)
  
  left_pad  <- start_pos - 1
  right_pad <- Ll - (start_pos + Ls - 1)
  
  pad_left  <- matrix(NA_real_, nrow=4, ncol=left_pad,
                      dimnames=list(rownames(short_mat), NULL))
  pad_right <- matrix(NA_real_, nrow=4, ncol=right_pad,
                      dimnames=list(rownames(short_mat), NULL))
  
  padded_short <- cbind(pad_left, short_mat, pad_right)
  list(short = padded_short, long = long_mat)
}

aligned <- align_pwms(my_best, uniprobe_mat, best_res$best_pos)

p1 <- ggseqlogo(aligned$short, method="prob") +
  ggtitle(paste0("GADEM motif (aligned; ", best_orient, ")")) +
  theme_bw()

p2 <- ggseqlogo(aligned$long, method="prob") +
  ggtitle("UniPROBE Bcl6b motif") +
  theme_bw()

grid.arrange(p1, p2, ncol = 1)

# Optional: also show unaligned side-by-side
# grid.arrange(
#   ggseqlogo(my_best, method="prob") + ggtitle("GADEM motif (best orientation)") + theme_bw(),
#   ggseqlogo(uniprobe_mat, method="prob") + ggtitle("UniPROBE motif") + theme_bw(),
#   ncol = 2
# )
#=================================================================#
###################################################################
#=================================================================#


#==============================================================#
#////////GSEA Based on motif max_score////////////////////////#
#==============================================================#
#=================================================#
#--------------------------------------------------#
#==================================================#

#This is good!! Shows cytoskeletal reorganization and GTPase mediated signal transduction.
#this is the Rho/ROCK pathway.
# GO Biological Process (BP) fgsea using max_score ranks (same ranks you just built)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(msigdbr)
  library(fgsea)
  library(readr)
  library(ggplot2)
})

# ---- safety: make sure we have the same ranks from the Hallmark run ----
if (!exists("ranks", inherits = FALSE)) {
  stop("I don't see `ranks` in memory. Re-run the section that builds `rank_df` + `ranks` from max_score.")
}
if (is.null(names(ranks))) stop("`ranks` must be a named numeric vector (names = gene symbols).")

# -----------------------------
# 1) Get GO BP gene sets (msigdbr: C5, subcategory GO:BP)
# -----------------------------
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(gene_symbol = toupper(gene_symbol))

pathways_bp <- split(gobp$gene_symbol, gobp$gs_name)

# Restrict to your universe; keep reasonable sizes
pathways_bp <- lapply(pathways_bp, intersect, y = names(ranks))
pathways_bp <- pathways_bp[lengths(pathways_bp) >= 15]   # bump min a bit for GO
cat("GO:BP pathways retained (>=15 overlapping genes):", length(pathways_bp), "\n")

# -----------------------------
# 2) Run fgsea
# -----------------------------
set.seed(1)
fg_bp <- fgsea(
  pathways = pathways_bp,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500,
  nperm    = 20000
) %>%
  arrange(padj, desc(abs(NES)))

# -----------------------------
# 3) Save + view
# -----------------------------
if (!exists("OUT_PREFIX", inherits = FALSE)) OUT_PREFIX <- "BCL6B"
out_bp <- paste0(OUT_PREFIX, "_fgsea_GOBP_rankByMaxScore.csv")
write_csv(fg_bp, out_bp)
message("Wrote: ", out_bp)

print(fg_bp %>% dplyr::select(pathway, NES, pval, padj, size) %>% head(25))

# -----------------------------
# 4) Quick plot (top 25 by FDR)
# -----------------------------
if (nrow(fg_bp) > 0) {
  top_bp <- fg_bp %>%
    arrange(padj) %>%
    slice_head(n = 25) %>%
    mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  ggplot(top_bp, aes(x = NES, y = pathway)) +
    geom_point() +
    theme_bw() +
    labs(
      title = "GO:BP enrichment of high-scoring BCL6B promoter motif matches",
      x = "NES (fgsea)",
      y = ""
    )
}

# -----------------------------
# Optional: reduce redundancy by collapsing similar GO terms
# (keeps the most significant term in each leading-edge overlap cluster)
# -----------------------------
if (nrow(fg_bp) > 0) {
  fg_bp_sig <- fg_bp %>% filter(padj < 0.05)
  if (nrow(fg_bp_sig) > 1) {
    fg_bp_collapsed <- fgsea::collapsePathways(
      fg_bp_sig,
      pathways = pathways_bp,
      stats = ranks
    )
    # fgsea returns "mainPathways" (keep) and "collapsedPathways" (drop)
    keep <- fg_bp_collapsed$mainPathways
    fg_bp_nonredundant <- fg_bp_sig %>% filter(pathway %in% keep) %>%
      arrange(padj, desc(abs(NES)))
    
    out_bp_nr <- paste0(OUT_PREFIX, "_fgsea_GOBP_rankByMaxScore_nonredundant.csv")
    write_csv(fg_bp_nonredundant, out_bp_nr)
    message("Wrote: ", out_bp_nr)
    
    print(fg_bp_nonredundant %>% dplyr::select(pathway, NES, pval, padj, size) %>% head(25))
  }
}
###################################################################################
###################################################################################

#==============================================================#
#   GO:BP fgsea from max_score ranks + ridge plot of leading-edge
#   (NO redundancy collapsing; ridge fill uses -log10(padj) for contrast)
#==============================================================#

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
  library(ggridges)
  library(stringr)
  library(forcats)
  library(scales)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
TOP_N_TERMS        <- 15      # number of GO terms to display
MIN_PADJ_FOR_PLOT  <- 0.25    # loosen to show more terms if needed (set to 1 to disable)
LABEL_N_PER_TERM   <- 5       # genes per term for Illustrator label table
RIDGE_COLOR_BY     <- "neglog10_padj"  # "neglog10_padj" (recommended) or "NES_rescaled"

# -----------------------------
# 0) Load inputs (if needed)
# -----------------------------
if (!exists("OUT_PREFIX", inherits = FALSE)) OUT_PREFIX <- "BCL6B"

if (!exists("motif_hits", inherits = FALSE)) {
  hits_file_guess <- paste0(OUT_PREFIX, "_MotifHits_hitLevel.csv")
  stopifnot(file.exists(hits_file_guess))
  motif_hits <- readr::read_csv(hits_file_guess, show_col_types = FALSE)
}

if (!exists("promoters_union", inherits = FALSE)) {
  if (file.exists("BCL6B_PromoterTargets_UNION.csv")) {
    promoters_union <- readr::read_csv("BCL6B_PromoterTargets_UNION.csv", show_col_types = FALSE) %>%
      dplyr::pull(1) %>% as.character()
  } else {
    stop("Need promoters_union in memory or BCL6B_PromoterTargets_UNION.csv on disk.")
  }
}

promoters_union <- unique(toupper(promoters_union))
motif_hits <- motif_hits %>% mutate(gene = toupper(gene))

# -----------------------------
# 1) Gene-level max_score
#    (tx_max then gene max)
# -----------------------------
motif_gene_max <- motif_hits %>%
  filter(!is.na(score)) %>%
  group_by(gene, transcript) %>%
  summarise(tx_max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(gene) %>%
  summarise(max_score = max(tx_max_score, na.rm = TRUE), .groups = "drop")

rank_df <- tibble(gene = promoters_union) %>%
  left_join(motif_gene_max, by = "gene") %>%
  mutate(max_score = replace_na(max_score, 0)) %>%
  arrange(desc(max_score), gene)

ranks <- rank_df$max_score
names(ranks) <- rank_df$gene
ranks <- sort(ranks, decreasing = TRUE)

# Break ties deterministically (important because max_score has lots of ties/zeros)
ranks2 <- ranks + (seq_along(ranks) * 1e-10)
ranks2 <- sort(ranks2, decreasing = TRUE)

# -----------------------------
# 2) GO:BP gene sets
# -----------------------------
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(gene_symbol = toupper(gene_symbol))

pathways_bp <- split(gobp$gene_symbol, gobp$gs_name)

# Restrict to your universe; keep sizes reasonable
pathways_bp <- lapply(pathways_bp, intersect, y = names(ranks2))
pathways_bp <- pathways_bp[lengths(pathways_bp) >= 15]
pathways_bp <- pathways_bp[lengths(pathways_bp) <= 500]

# -----------------------------
# 3) Run fgseaMultilevel (returns leadingEdge)
# -----------------------------
set.seed(1)
fg_bp <- fgsea::fgseaMultilevel(
  pathways = pathways_bp,
  stats    = ranks2,
  minSize  = 15,
  maxSize  = 500
) %>%
  arrange(padj, desc(abs(NES)))

# Save full table
out_bp <- paste0(OUT_PREFIX, "_fgsea_GOBP_rankByMaxScore_multilevel.csv")
readr::write_csv(fg_bp, out_bp)
message("Wrote: ", out_bp)

stopifnot("leadingEdge" %in% colnames(fg_bp))

# -----------------------------
# 4) Pick terms for ridge plot (NO collapsing)
# -----------------------------
plot_terms <- fg_bp %>%
  filter(!is.na(padj)) %>%
  filter(padj <= MIN_PADJ_FOR_PLOT) %>%
  arrange(padj) %>%
  slice_head(n = TOP_N_TERMS)

if (nrow(plot_terms) == 0) {
  stop("No pathways passed your plot filter. Try increasing MIN_PADJ_FOR_PLOT (e.g., 1) or lowering minSize.")
}

# -----------------------------
# 5) Build ridge_df:
#    each row = (pathway, gene in leadingEdge, gene_max_score)
# -----------------------------
ridge_df <- plot_terms %>%
  dplyr::select(pathway, NES, padj, leadingEdge) %>%
  tidyr::unnest(leadingEdge) %>%
  dplyr::rename(gene = leadingEdge) %>%
  left_join(motif_gene_max, by = "gene") %>%
  mutate(
    max_score = replace_na(max_score, 0),
    pathway_clean = pathway %>%
      stringr::str_replace("^GOBP_", "") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_title(),
    neglog10_padj = -log10(padj + 1e-300)
  )

# Make sure each pathway has a single fill value repeated across its genes
ridge_df <- ridge_df %>%
  group_by(pathway_clean) %>%
  mutate(
    NES_term = first(NES),
    padj_term = first(padj),
    neglog10_padj_term = first(neglog10_padj),
    NES_color = scales::rescale(NES_term, to = c(0, 1))
  ) %>%
  ungroup()

# Order pathways (most significant on top)
ridge_df <- ridge_df %>%
  mutate(
    pathway_clean = forcats::fct_reorder(
      pathway_clean,
      neglog10_padj_term,
      .desc = FALSE   # <-- flip direction
    )
  )

# -----------------------------
# 6) Ridge plot
#    y = pathway, x = max_score among leading-edge genes
#    fill = -log10(padj) (recommended) OR rescaled NES
# -----------------------------
if (RIDGE_COLOR_BY == "NES_rescaled") {
  p_ridge <- ggplot(ridge_df, aes(x = max_score, y = pathway_clean, fill = NES_color)) +
    ggridges::geom_density_ridges(
      scale = 1.15,
      rel_min_height = 0.01,
      alpha = 0.90,
      color = "grey25",
      size = 0.25
    ) +
    scale_fill_gradient(
      low = "grey90",
      high = "red3",
      name = "NES\n(rescaled)"
    ) +
    labs(
      x = "BCL6B promoter motif strength (gene max PWM score)",
      y = "",
      title = "GO:BP enrichment: distribution of motif strengths in leading-edge genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 10)
    )
} else {
  # default: color by -log10(padj) to avoid “all same color” when NES is narrow
  p_ridge <- ggplot(ridge_df, aes(x = max_score, y = pathway_clean, fill = neglog10_padj_term)) +
    ggridges::geom_density_ridges(
      scale = 1.15,
      rel_min_height = 0.01,
      alpha = 0.90,
      color = "grey25",
      size = 0.25
    ) +
    scale_fill_gradient(
      low = "grey90",
      high = "red3",
      name = expression(-log[10](padj))
    ) +
    labs(
      x = "BCL6B promoter motif strength (gene max PWM score)",
      y = "",
      title = "GO:BP enrichment: distribution of motif strengths in leading-edge genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 10)
    )
}

print(p_ridge)

# Save plot
ggsave(
  filename = paste0(OUT_PREFIX, "_GOBP_ridge_leadingEdge_maxScore_noCollapse.png"),
  plot = p_ridge,
  width = 9, height = 7, dpi = 300
)

# -----------------------------
# 7) Label-candidate table for Illustrator
#    pick top genes by max_score in each pathway (leadingEdge only)
# -----------------------------
label_table <- ridge_df %>%
  group_by(pathway_clean, NES_term, padj_term, neglog10_padj_term) %>%
  arrange(desc(max_score), gene) %>%
  slice_head(n = LABEL_N_PER_TERM) %>%
  ungroup() %>%
  arrange(desc(neglog10_padj_term), desc(NES_term), pathway_clean)

out_labels <- paste0(OUT_PREFIX, "_GOBP_ridge_labelCandidates_topGenes_noCollapse.csv")
readr::write_csv(label_table, out_labels)
message("Wrote: ", out_labels)

label_table %>% print(n = 50)

label_table %>%
  filter(pathway_clean == "Cell Adhesion") %>%
  pull(gene)

# -----------------------------
# 8) (Optional) Save the rank table too
# -----------------------------
out_rank <- paste0(OUT_PREFIX, "_geneRank_maxScore.csv")
readr::write_csv(rank_df, out_rank)
message("Wrote: ", out_rank)

bcl6b_encode_leading_genes <- label_table %>%
  pull(gene)


bcl6b_partial %>%
  filter(GeneSymbol %in% bcl6b_encode_leading_genes) %>%
  arrange(desc(r_gene_BCL6B))
