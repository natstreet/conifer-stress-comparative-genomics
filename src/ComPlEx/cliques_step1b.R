#!/usr/bin/env Rscript
# cliques_step1b.R
#
# Computes per-gene-pair tissue coverage from the combined co-expressologs
# table produced by cliques_step1.R.
#
# This replaces the igraph max_cliques step in 3_Cliques.Rmd (Steps 1–2).
#
# NOTE on the bipartite clique adaptation:
#   The original pipeline ran igraph::max_cliques() on a graph where nodes
#   were (species, gene) combinations across many species-pair comparisons.
#   Triangles formed because the same gene could appear in nodes from 3+
#   species and those species pairs all produced co-expressologs.
#   In this study there are exactly 2 species (spruce and pine); every edge
#   connects one spruce gene to one pine gene. The resulting graph is strictly
#   bipartite, and bipartite graphs contain no cliques of size ≥ 3. Running
#   max_cliques(net, min=3) on these data returns nothing.
#
#   The equivalent metric for 2 species / N tissue conditions is tissue
#   coverage: for each (spruce_gene, pine_gene) pair within an orthogroup,
#   count how many of the 4 tissue conditions show it as a co-expressolog,
#   and compute NegLog10CliqueSum = sum(−log10(MaxpVal)) across those tissues.
#   This directly parallels the original clique-sum scoring and supports the
#   same conserved / stress-specific / tissue-specific classifications.
#
# Input:  results/ComPlEx/RData/co_expressologs.RData
# Output: results/ComPlEx/RData/weighted_gene_pairs.RData
#           one row per unique (OrthoGroup, Species1, Species2) gene pair
#           columns mirror weighted_max_cliques_filterable from 3_Cliques.Rmd

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, "../.."))
RESULTS_DIR <- file.path(PROJECT_DIR, "results", "ComPlEx")
RDATA_DIR   <- file.path(RESULTS_DIR, "RData")

cat("Loading co_expressologs...\n")
load(file.path(RDATA_DIR, "co_expressologs.RData"))   # → co_expressologs

TISSUES <- c("cold_needle", "cold_root", "drought_needle", "drought_root")
# Stress categories for sum columns
COLD    <- c("cold_needle",    "cold_root")
DROUGHT <- c("drought_needle", "drought_root")
NEEDLE  <- c("cold_needle",    "drought_needle")
ROOT    <- c("cold_root",      "drought_root")

# ── Pivot to one row per gene pair ────────────────────────────────────────────
# For each (OrthoGroup, Species1, Species2) we want:
#   - per-tissue MaxpVal (NA if not a co-expressolog in that tissue)
#   - per-tissue presence indicator (0/1)
#   - aggregate scores
cat("Pivoting to one row per gene pair...\n")

gene_pairs <- co_expressologs %>%
  select(OrthoGroup, Species1, Species2, Tissue, MaxpVal) %>%
  pivot_wider(
    names_from   = Tissue,
    values_from  = MaxpVal,
    names_glue   = "{Tissue}_pval"
  )

# Ensure all 4 tissue p-val columns exist (NA if absent)
for (t in TISSUES) {
  col <- paste0(t, "_pval")
  if (!col %in% colnames(gene_pairs)) gene_pairs[[col]] <- NA_real_
}

# ── Presence indicators (1 = co-expressolog in that tissue) ───────────────────
for (t in TISSUES) {
  gene_pairs[[paste0(t, "_present")]] <- as.integer(
    !is.na(gene_pairs[[paste0(t, "_pval")]])
  )
}

# ── NegLog10 scores ───────────────────────────────────────────────────────────
neg_log10_pval <- function(x) ifelse(is.na(x), 0, -log10(pmax(x, 1e-300)))

for (t in TISSUES) {
  gene_pairs[[paste0(t, "_neglog10")]] <- neg_log10_pval(
    gene_pairs[[paste0(t, "_pval")]]
  )
}

# Per-pair aggregate scores (analogues of AngioSum, GymnoSum, CrossSum)
gene_pairs <- gene_pairs %>%
  mutate(
    ColdSum    = rowSums(across(paste0(COLD,    "_neglog10"))),
    DroughtSum = rowSums(across(paste0(DROUGHT, "_neglog10"))),
    NeedleSum  = rowSums(across(paste0(NEEDLE,  "_neglog10"))),
    RootSum    = rowSums(across(paste0(ROOT,    "_neglog10"))),
    # Total NegLog10CliqueSum across all tissues (analogous to original)
    NegLog10CliqueSum = rowSums(across(paste0(TISSUES, "_neglog10"))),
    n_tissues = rowSums(across(paste0(TISSUES, "_present"))),
    # Best (most significant) MaxpVal across tissues this pair appears in
    best_pval = pmin(
      cold_needle_pval, cold_root_pval,
      drought_needle_pval, drought_root_pval,
      na.rm = TRUE
    ),
    # Tissue pattern label
    tissue_pattern = paste0(
      ifelse(cold_needle_present    == 1, "CN", ""),
      ifelse(cold_root_present      == 1, "CR", ""),
      ifelse(drought_needle_present == 1, "DN", ""),
      ifelse(drought_root_present   == 1, "DR", "")
    )
  )

# ── Per-pair stress/tissue classification (mutually exclusive, hierarchical) ──
# Mirrors the Angio/Gymno/Cross logic from 3_Cliques.Rmd but for tissue axes.
gene_pairs <- gene_pairs %>%
  mutate(
    conserved         = (n_tissues == 4),
    cold_specific     = (cold_needle_present | cold_root_present) &
                        !drought_needle_present & !drought_root_present,
    drought_specific  = !cold_needle_present & !cold_root_present &
                        (drought_needle_present | drought_root_present),
    needle_specific   = (cold_needle_present | drought_needle_present) &
                        !cold_root_present & !drought_root_present,
    root_specific     = !cold_needle_present & !drought_needle_present &
                        (cold_root_present | drought_root_present)
  )

# ── Summary statistics ────────────────────────────────────────────────────────
n_pairs   <- nrow(gene_pairs)
n_ogs     <- n_distinct(gene_pairs$OrthoGroup)

cat(sprintf("\nUnique gene pairs: %d  across  %d orthogroups\n", n_pairs, n_ogs))
cat(sprintf("Tissue coverage:\n"))
for (t in TISSUES) {
  cat(sprintf("  %s: %d pairs (%.1f%%)\n",
              t,
              sum(gene_pairs[[paste0(t, "_present")]]),
              100 * mean(gene_pairs[[paste0(t, "_present")]])))
}
cat(sprintf("\nTissue pattern counts:\n"))
cat(sprintf("  Conserved (all 4 tissues):   %d pairs  (%d OGs)\n",
            sum(gene_pairs$conserved),
            n_distinct(gene_pairs$OrthoGroup[gene_pairs$conserved])))
cat(sprintf("  Cold-specific (cold only):   %d pairs  (%d OGs)\n",
            sum(gene_pairs$cold_specific),
            n_distinct(gene_pairs$OrthoGroup[gene_pairs$cold_specific])))
cat(sprintf("  Drought-specific:            %d pairs  (%d OGs)\n",
            sum(gene_pairs$drought_specific),
            n_distinct(gene_pairs$OrthoGroup[gene_pairs$drought_specific])))
cat(sprintf("  Needle-specific:             %d pairs  (%d OGs)\n",
            sum(gene_pairs$needle_specific),
            n_distinct(gene_pairs$OrthoGroup[gene_pairs$needle_specific])))
cat(sprintf("  Root-specific:               %d pairs  (%d OGs)\n",
            sum(gene_pairs$root_specific),
            n_distinct(gene_pairs$OrthoGroup[gene_pairs$root_specific])))
cat(sprintf("  Mixed (multiple, not all 4): %d pairs  (%d OGs)\n",
            sum(!gene_pairs$conserved & !gene_pairs$cold_specific &
                !gene_pairs$drought_specific & !gene_pairs$needle_specific &
                !gene_pairs$root_specific),
            n_distinct(gene_pairs$OrthoGroup[
              !gene_pairs$conserved & !gene_pairs$cold_specific &
              !gene_pairs$drought_specific & !gene_pairs$needle_specific &
              !gene_pairs$root_specific])))

# ── Best gene pair per orthogroup (lowest best_pval) ─────────────────────────
# Analogous to slice_min(MaxpVal) in the original cliques_step1.R
best_per_og <- gene_pairs %>%
  group_by(OrthoGroup) %>%
  slice_min(order_by = best_pval, n = 1, with_ties = FALSE) %>%
  ungroup()

# ── Save ──────────────────────────────────────────────────────────────────────
weighted_gene_pairs <- gene_pairs
save(weighted_gene_pairs, file = file.path(RDATA_DIR, "weighted_gene_pairs.RData"))
save(best_per_og, file = file.path(RDATA_DIR, "best_gene_pair_per_og.RData"))

# Also write a TSV for easy inspection
write.table(weighted_gene_pairs,
            file.path(RDATA_DIR, "weighted_gene_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nSaved to %s/\n", RDATA_DIR))
cat("  weighted_gene_pairs.RData\n")
cat("  best_gene_pair_per_og.RData\n")
cat("  weighted_gene_pairs.tsv\n")
