#!/usr/bin/env Rscript
# cliques_step2.R
#
# Classifies gene pairs and orthogroups into co-expression conservation
# categories based on tissue coverage.
#
# Adapted from the original cliques_step2.R / 3_Cliques.Rmd for the
# spruce-pine 2-species, 4-tissue design. The original Angio/Gymno/Cross
# axes (species-clade categories across many species pairs) are replaced
# by stress (cold / drought) and tissue (needle / root) axes.
#
# Input:  results/ComPlEx/RData/weighted_gene_pairs.RData
#         results/ComPlEx/RData/orthogroup_coexpressolog_presence.RData
# Output: results/ComPlEx/RData/  (one RData per category)

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, "../.."))
RESULTS_DIR <- file.path(PROJECT_DIR, "results", "ComPlEx")
RDATA_DIR   <- file.path(RESULTS_DIR, "RData")

cat("Loading weighted_gene_pairs and orthogroup_coexpressolog_presence...\n")
load(file.path(RDATA_DIR, "weighted_gene_pairs.RData"))   # → weighted_gene_pairs
load(file.path(RDATA_DIR, "orthogroup_coexpressolog_presence.RData"))    # → orthogroup_coexpressolog_presence

# ── Parameters (analogous to original cliques_step2.R) ──────────────────────
SELECTION_METHOD  <- "SUM"   # "SUM" = NegLog10CliqueSum  |  "PVAL" = best_pval
TYPE_OF_GENE_SET  <- "c"     # "c" conserved | "cold" cold-specific |
                              # "drought" drought-specific | "needle" needle-specific |
                              # "root" root-specific

# Score thresholds — calibrate against the 4-tissue NegLog10CliqueSum scale.
# Each tissue contributes up to ~300 (the max −log10 for FDR<0.05 pairs near 0).
# Conserved pairs (all 4 tissues) accumulate higher sums than tissue-specific ones.
C_SUM        <- 10    # minimum NegLog10CliqueSum for conserved pairs
C_SUM_STRESS <-  5    # minimum NegLog10CliqueSum for stress-specific (2 tissues max)
P_VAL        <-  0.01 # best_pval < P_VAL threshold when using PVAL method

# ── Orthogroup-level presence summary ────────────────────────────────────────
# Mirrors the orthogroup_coexpressolog_presence usage in the original for filtering OGs.
og_summary <- orthogroup_coexpressolog_presence %>%
  as.data.frame() %>%
  mutate(
    n_cold_tissues    = cold_needle    + cold_root,
    n_drought_tissues = drought_needle + drought_root,
    n_needle_tissues  = cold_needle    + drought_needle,
    n_root_tissues    = cold_root      + drought_root,
    n_total_tissues   = cold_needle + cold_root + drought_needle + drought_root
  )

# ── Category: CONSERVED ───────────────────────────────────────────────────────
# Gene pairs present as co-expressologs in all 4 tissue conditions.
if (TYPE_OF_GENE_SET == "c") {
  cat("\n--- CONSERVED ---\n")

  if (SELECTION_METHOD == "SUM") {
    conserved_genes <- weighted_gene_pairs %>%
      filter(conserved) %>%
      filter(NegLog10CliqueSum >= C_SUM)
  }
  if (SELECTION_METHOD == "PVAL") {
    conserved_genes <- weighted_gene_pairs %>%
      filter(conserved) %>%
      filter(best_pval < P_VAL)
  }

  cat(sprintf("Conserved gene pairs: %d  (%d orthogroups)\n",
              nrow(conserved_genes),
              n_distinct(conserved_genes$OrthoGroup)))
  file_name <- file.path(RDATA_DIR,
    sprintf("conserved_genes_%s_%s.RData", SELECTION_METHOD,
            ifelse(SELECTION_METHOD=="SUM", C_SUM, P_VAL)))
  save(conserved_genes, file = file_name)
  cat(sprintf("Saved: %s\n", file_name))
}

# ── Category: COLD-SPECIFIC ───────────────────────────────────────────────────
# Gene pairs that are co-expressologs only in cold tissues (needle and/or root)
# but NOT in drought needle or drought root.
if (TYPE_OF_GENE_SET == "cold") {
  cat("\n--- COLD-SPECIFIC ---\n")

  if (SELECTION_METHOD == "SUM") {
    cold_genes <- weighted_gene_pairs %>%
      filter(cold_specific) %>%
      filter(ColdSum >= C_SUM_STRESS)
  }
  if (SELECTION_METHOD == "PVAL") {
    cold_genes <- weighted_gene_pairs %>%
      filter(cold_specific) %>%
      filter(best_pval < P_VAL)
  }

  cat(sprintf("Cold-specific gene pairs: %d  (%d orthogroups)\n",
              nrow(cold_genes),
              n_distinct(cold_genes$OrthoGroup)))
  file_name <- file.path(RDATA_DIR,
    sprintf("cold_specific_genes_%s_%s.RData", SELECTION_METHOD,
            ifelse(SELECTION_METHOD=="SUM", C_SUM_STRESS, P_VAL)))
  save(cold_genes, file = file_name)
  cat(sprintf("Saved: %s\n", file_name))
}

# ── Category: DROUGHT-SPECIFIC ───────────────────────────────────────────────
# Gene pairs that are co-expressologs only in drought tissues.
if (TYPE_OF_GENE_SET == "drought") {
  cat("\n--- DROUGHT-SPECIFIC ---\n")

  if (SELECTION_METHOD == "SUM") {
    drought_genes <- weighted_gene_pairs %>%
      filter(drought_specific) %>%
      filter(DroughtSum >= C_SUM_STRESS)
  }
  if (SELECTION_METHOD == "PVAL") {
    drought_genes <- weighted_gene_pairs %>%
      filter(drought_specific) %>%
      filter(best_pval < P_VAL)
  }

  cat(sprintf("Drought-specific gene pairs: %d  (%d orthogroups)\n",
              nrow(drought_genes),
              n_distinct(drought_genes$OrthoGroup)))
  file_name <- file.path(RDATA_DIR,
    sprintf("drought_specific_genes_%s_%s.RData", SELECTION_METHOD,
            ifelse(SELECTION_METHOD=="SUM", C_SUM_STRESS, P_VAL)))
  save(drought_genes, file = file_name)
  cat(sprintf("Saved: %s\n", file_name))
}

# ── Category: NEEDLE-SPECIFIC ────────────────────────────────────────────────
# Gene pairs present in needle tissues (cold_needle and/or drought_needle)
# but NOT in root tissues.
if (TYPE_OF_GENE_SET == "needle") {
  cat("\n--- NEEDLE-SPECIFIC ---\n")

  if (SELECTION_METHOD == "SUM") {
    needle_genes <- weighted_gene_pairs %>%
      filter(needle_specific) %>%
      filter(NeedleSum >= C_SUM_STRESS)
  }
  if (SELECTION_METHOD == "PVAL") {
    needle_genes <- weighted_gene_pairs %>%
      filter(needle_specific) %>%
      filter(best_pval < P_VAL)
  }

  cat(sprintf("Needle-specific gene pairs: %d  (%d orthogroups)\n",
              nrow(needle_genes),
              n_distinct(needle_genes$OrthoGroup)))
  file_name <- file.path(RDATA_DIR,
    sprintf("needle_specific_genes_%s_%s.RData", SELECTION_METHOD,
            ifelse(SELECTION_METHOD=="SUM", C_SUM_STRESS, P_VAL)))
  save(needle_genes, file = file_name)
  cat(sprintf("Saved: %s\n", file_name))
}

# ── Category: ROOT-SPECIFIC ───────────────────────────────────────────────────
# Gene pairs present in root tissues (cold_root and/or drought_root)
# but NOT in needle tissues.
if (TYPE_OF_GENE_SET == "root") {
  cat("\n--- ROOT-SPECIFIC ---\n")

  if (SELECTION_METHOD == "SUM") {
    root_genes <- weighted_gene_pairs %>%
      filter(root_specific) %>%
      filter(RootSum >= C_SUM_STRESS)
  }
  if (SELECTION_METHOD == "PVAL") {
    root_genes <- weighted_gene_pairs %>%
      filter(root_specific) %>%
      filter(best_pval < P_VAL)
  }

  cat(sprintf("Root-specific gene pairs: %d  (%d orthogroups)\n",
              nrow(root_genes),
              n_distinct(root_genes$OrthoGroup)))
  file_name <- file.path(RDATA_DIR,
    sprintf("root_specific_genes_%s_%s.RData", SELECTION_METHOD,
            ifelse(SELECTION_METHOD=="SUM", C_SUM_STRESS, P_VAL)))
  save(root_genes, file = file_name)
  cat(sprintf("Saved: %s\n", file_name))
}

# ── Run all categories at once ────────────────────────────────────────────────
# Convenience: produce an orthogroup-level summary table with category labels.
og_categories <- weighted_gene_pairs %>%
  group_by(OrthoGroup) %>%
  summarise(
    n_pairs            = n(),
    best_pval          = min(best_pval),
    max_NegLog10Sum    = max(NegLog10CliqueSum),
    any_conserved      = any(conserved),
    # conserved AND passing the NegLog10Sum threshold on the SAME pair, matching the
    # Methods definition and the conserved_genes_SUM_10 set (not any_conserved alone)
    any_conserved_sum10 = any(conserved & NegLog10CliqueSum >= C_SUM),
    any_cold_specific  = any(cold_specific),
    any_drought_specific = any(drought_specific),
    any_needle_specific  = any(needle_specific),
    any_root_specific    = any(root_specific),
    .groups = "drop"
  ) %>%
  mutate(
    category = case_when(
      any_conserved_sum10   ~ "conserved",
      any_cold_specific & any_drought_specific ~ "both_stresses",
      any_cold_specific     ~ "cold_specific",
      any_drought_specific  ~ "drought_specific",
      any_needle_specific & any_root_specific  ~ "both_tissues",
      any_needle_specific   ~ "needle_specific",
      any_root_specific     ~ "root_specific",
      TRUE                  ~ "single_tissue"
    )
  )

cat(sprintf("\n--- Orthogroup category summary ---\n"))
print(og_categories %>% count(category, sort = TRUE))

write.table(og_categories,
            file.path(RDATA_DIR, "orthogroup_categories.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSaved orthogroup_categories.tsv to %s/\n", RDATA_DIR))
