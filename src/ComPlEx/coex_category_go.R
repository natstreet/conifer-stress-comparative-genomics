#!/usr/bin/env Rscript
# coex_category_go.R
#
# GO enrichment for the co-expressolog categories (conserved / cold_specific /
# drought_specific / multi_tissue / not_coex) from integration_analysis.R. Three
# comparisons: each category vs the full background; conserved vs not_coex; and
# cold_specific vs drought_specific. Uses the single shared GO method
# (src/lib/go_enrichment.R): topGO weight01/Fisher, nodeSize 5, no BH, BP/MF/CC,
# all terms with a plant_consistent flag. Background = ComPlEx expression universe
# (SC+SD) intersected with the eggNOG+InterPro annotation (data/annotation/gene_annotation.tsv.gz).
#
# Run from AbioticStressConifers/.  Outputs in results/integration/go/.
suppressPackageStartupMessages({ library(data.table); library(here) })
source("src/lib/go_enrichment.R")

INTEG_DIR <- "results/integration"
EXPR_DIR  <- "data/expression"
OUTDIR    <- file.path(INTEG_DIR, "go"); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

gene2go <- build_gene2go(here("data/annotation/gene_annotation.tsv.gz"))
expr_universe <- unique(c(fread(file.path(EXPR_DIR, "SC_expression.txt"), select = 1)[[1]],
                          fread(file.path(EXPR_DIR, "SD_expression.txt"), select = 1)[[1]]))
bg_genes <- intersect(expr_universe, names(gene2go))
cat(sprintf("Background (expression universe ∩ annotated): %d genes\n", length(bg_genes)))

cat_files <- list(
  conserved        = "category_genes_conserved.tsv",
  cold_specific    = "category_genes_cold_specific.tsv",
  drought_specific = "category_genes_drought_specific.tsv",
  multi_tissue     = "category_genes_multi_tissue.tsv",
  not_coex         = "category_genes_not_coex.tsv")
cat_genes <- lapply(cat_files, function(f) intersect(fread(file.path(INTEG_DIR, f))[[1]], bg_genes))
for (nm in names(cat_genes)) cat(sprintf("  %-18s %d genes (in background)\n", nm, length(cat_genes[[nm]])))

# one enrichment group -> data.table with a `category` label (all ontologies + plant flag)
go_group <- function(sig, bg, label) {
  out <- run_go(sig, bg, gene2go, ontologies = c("BP", "MF", "CC"))
  if (nrow(out) == 0) return(NULL)
  out[, category := label][]
}

# 1. each category vs the full annotated background
combined <- rbindlist(lapply(names(cat_genes), function(nm) go_group(cat_genes[[nm]], bg_genes, nm)), fill = TRUE)
fwrite(combined, file.path(OUTDIR, "coex_category_go_vs_bg.tsv"), sep = "\t")
cat(sprintf("Saved coex_category_go_vs_bg.tsv (%d rows)\n", nrow(combined)))

# 2. conserved vs not_coex (shared background = the two sets)
bg1 <- union(cat_genes$conserved, cat_genes$not_coex)
p1  <- rbindlist(list(go_group(cat_genes$conserved, bg1, "conserved_vs_not_coex"),
                      go_group(cat_genes$not_coex,  bg1, "not_coex_vs_conserved")), fill = TRUE)
fwrite(p1, file.path(OUTDIR, "go_conserved_vs_not_coex.tsv"), sep = "\t")
cat(sprintf("Saved go_conserved_vs_not_coex.tsv (%d rows)\n", nrow(p1)))

# 3. cold_specific vs drought_specific
bg2 <- union(cat_genes$cold_specific, cat_genes$drought_specific)
p2  <- rbindlist(list(go_group(cat_genes$cold_specific,    bg2, "cold_specific_vs_drought_bg"),
                      go_group(cat_genes$drought_specific, bg2, "drought_specific_vs_cold_bg")), fill = TRUE)
fwrite(p2, file.path(OUTDIR, "go_cold_vs_drought.tsv"), sep = "\t")
cat(sprintf("Saved go_cold_vs_drought.tsv (%d rows)\n", nrow(p2)))

# console: top plant-consistent BP terms per category
cat("\n-- Top BP (plant-consistent) terms per category vs background --\n")
for (nm in names(cat_genes)) {
  sub <- combined[category == nm & ont == "BP" & plant_consistent][order(p)][1:min(8, .N)]
  if (nrow(sub) == 0 || is.na(sub$GO.ID[1])) next
  cat(sprintf("\n  %s:\n", nm)); print(sub[, .(GO.ID, Term, Significant, p = signif(p, 2))], row.names = FALSE)
}
cat(sprintf("\nOutputs in: %s\n", OUTDIR))
