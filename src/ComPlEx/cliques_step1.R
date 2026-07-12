#!/usr/bin/env Rscript
# cliques_step1.R
#
# Loads the four per-tissue ComPlEx TSVs produced by complex_py.py and
# combines them into a single co-expressologs table plus an
# orthogroup × tissue presence/absence matrix (orthogroup_coexpressolog_presence).
#
# Replaces the original cliques_step1.R which loaded a single combined
# RData file and filtered by Species1/Species2 tissue codes. In the new
# Python pipeline, each tissue is a separate TSV with gene names (not
# tissue codes) in the Species1 / Species2 columns.
#
# Outputs (saved under results/ComPlEx/RData/):
#   co_expressologs.RData   — combined data frame, one row per gene pair per tissue
#   orthogroup_coexpressolog_presence.RData — orthogroup × tissue matrix (0/1)

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, "../.."))
RESULTS_DIR <- file.path(PROJECT_DIR, "results", "ComPlEx")
RDATA_DIR   <- file.path(RESULTS_DIR, "RData")
dir.create(RDATA_DIR, showWarnings = FALSE, recursive = TRUE)

TISSUES <- c("cold_needle", "cold_root", "drought_needle", "drought_root")

# ── Load 4 tissue TSVs ────────────────────────────────────────────────────────
cat("Loading per-tissue co-expressolog TSVs...\n")

co_expr_list <- lapply(TISSUES, function(tiss) {
  tsv_path <- file.path(RESULTS_DIR, tiss, "RData", "comparison_tables",
                         "comparison_spruce_pine.tsv")
  if (!file.exists(tsv_path)) {
    warning(sprintf("TSV not found: %s", tsv_path))
    return(NULL)
  }
  df <- read.table(tsv_path, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)
  df$Tissue <- tiss
  cat(sprintf("  %s: %d co-expressologs\n", tiss, nrow(df)))
  df
})

co_expressologs <- bind_rows(Filter(Negate(is.null), co_expr_list))

# Rename Max.p.val → MaxpVal to match downstream scripts
co_expressologs <- co_expressologs %>%
  rename(MaxpVal = Max.p.val)

cat(sprintf("\nTotal rows (gene pairs × tissues): %d\n", nrow(co_expressologs)))
cat(sprintf("Unique orthogroups: %d\n", n_distinct(co_expressologs$OrthoGroup)))
cat(sprintf("Unique gene pairs (spruce × pine): %d\n",
            n_distinct(paste(co_expressologs$Species1, co_expressologs$Species2))))

# ── Build orthogroup_coexpressolog_presence (orthogroup × tissue) ───────────────────────────
# 1 = orthogroup has at least one co-expressolog gene pair in that tissue
cat("\nBuilding orthogroup × tissue presence matrix...\n")

orthogroup_coexpressolog_presence <- co_expressologs %>%
  select(OrthoGroup, Tissue) %>%
  distinct() %>%
  mutate(Value = 1L) %>%
  pivot_wider(names_from = Tissue, values_from = Value, values_fill = 0L) %>%
  column_to_rownames("OrthoGroup")

# Ensure all 4 tissue columns present in defined order
for (t in TISSUES) {
  if (!t %in% colnames(orthogroup_coexpressolog_presence)) orthogroup_coexpressolog_presence[[t]] <- 0L
}
orthogroup_coexpressolog_presence <- orthogroup_coexpressolog_presence[, TISSUES]

cat(sprintf("Orthogroups with co-expressologs in each tissue:\n"))
tissue_counts <- colSums(orthogroup_coexpressolog_presence > 0)
for (t in names(tissue_counts)) {
  cat(sprintf("  %s: %d\n", t, tissue_counts[t]))
}

n_tissues_per_og <- rowSums(orthogroup_coexpressolog_presence)
cat(sprintf("\nOrthogroups present in all 4 tissues:    %d\n",
            sum(n_tissues_per_og == 4)))
cat(sprintf("Orthogroups present in >= 2 tissues:     %d\n",
            sum(n_tissues_per_og >= 2)))
cat(sprintf("Orthogroups present in exactly 1 tissue: %d\n",
            sum(n_tissues_per_og == 1)))

# ── Save ──────────────────────────────────────────────────────────────────────
save(co_expressologs,  file = file.path(RDATA_DIR, "co_expressologs.RData"))
save(orthogroup_coexpressolog_presence, file = file.path(RDATA_DIR, "orthogroup_coexpressolog_presence.RData"))

cat(sprintf("\nSaved to %s/\n", RDATA_DIR))
cat("  co_expressologs.RData\n")
cat("  orthogroup_coexpressolog_presence.RData\n")
