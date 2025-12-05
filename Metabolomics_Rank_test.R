## Expect: imputed_log and metabo_age_signature are already defined
library(dplyr); library(fgsea); library(tibble); library(purrr)
set.seed(1)

met_cols_all <- setdiff(names(imputed_log), c("Sample","Group"))

# 1) Build OFR vs OS stats once
ofr_os <- imputed_log %>% dplyr::filter(Group %in% c("OFR","OS"))
t_ofr <- purrr::map_dfr(met_cols_all, function(m) {
  vals <- ofr_os[[m]]; grp <- ofr_os$Group
  tt <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
  tibble(
    Metabolite = m,
    LogFC  = mean(vals[grp=="OFR"], na.rm=TRUE) - mean(vals[grp=="OS"], na.rm=TRUE),
    P_value = if (is.null(tt)) NA_real_ else tt$p.value
  )
})

# 2) Make both ranking flavors
ranks_unweighted <- setNames(t_ofr$LogFC, t_ofr$Metabolite)
ranks_unweighted <- ranks_unweighted[is.finite(ranks_unweighted)]

t_ofr2 <- t_ofr %>%
  mutate(p_safe = pmax(P_value, 1e-300),
         stat   = -log10(p_safe) * LogFC)
ranks_weighted <- setNames(t_ofr2$stat, t_ofr2$Metabolite)
ranks_weighted <- ranks_weighted[is.finite(ranks_weighted)]

# 3) Ensure the pathway is the same and present
path <- intersect(metabo_age_signature, names(ranks_unweighted))
path2 <- intersect(metabo_age_signature, names(ranks_weighted))
cat("Set sizes — intended:", length(metabo_age_signature),
    " in unweighted:", length(path),
    " in weighted:", length(path2), "\n")
if (!setequal(path, path2)) cat("WARNING: members missing in one ranking.\n")

# 4) Compute NES both ways
fg_old <- fgsea(pathways = list(Signature = path),  stats = ranks_unweighted, nperm = 1000)
fg_new <- fgsea(pathways = list(Signature = path2), stats = ranks_weighted,  nperm = 1000)

fg_old$NES  # <-- should match your Original sign
fg_new$NES  # <-- should match your New sign

# Optional: inspect leading edges driving each result
fg_old$leadingEdge
fg_new$leadingEdge
