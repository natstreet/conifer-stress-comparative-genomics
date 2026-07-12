#!/usr/bin/env Rscript
# subset_go.R
#
# GO enrichment for two not-co-expressed DE subsets from not_coex_de_analysis.R:
#   spruce_lineage_specific_de  -- spruce DEGs with no 1:1 pine ortholog
#   diverged_regulation         -- both species DE but not co-expressed
# Uses the single shared GO method (src/lib/go_enrichment.R): topGO weight01/Fisher,
# nodeSize 5, no BH, BP/MF/CC, all terms with a plant_consistent flag. Background =
# ComPlEx expression universe (SC+SD) intersected with the eggNOG+InterPro annotation.
#
# Run from AbioticStressConifers/.  Output: results/integration/go/go_specific_subsets.tsv
suppressPackageStartupMessages({ library(data.table); library(here) })
source("src/lib/go_enrichment.R")

INTEG_DIR <- "results/integration"
EXPR_DIR  <- "data/expression"
OUTDIR    <- file.path(INTEG_DIR, "go"); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

gene2go <- build_gene2go(here("data/annotation/gene_annotation.tsv.gz"))
expr_universe <- unique(c(fread(file.path(EXPR_DIR, "SC_expression.txt"), select = 1)[[1]],
                          fread(file.path(EXPR_DIR, "SD_expression.txt"), select = 1)[[1]]))
bg_genes <- intersect(expr_universe, names(gene2go))

not_coex_de <- fread(file.path(INTEG_DIR, "not_coex_de_genes.tsv"))
lin_spec <- intersect(not_coex_de[mechanism == "no_1to1_ortholog",    pa_gene], bg_genes)
div_reg  <- intersect(not_coex_de[mechanism == "diverged_regulation", pa_gene], bg_genes)
cat(sprintf("background %d ; lineage-specific %d ; diverged-regulation %d (annotated)\n",
            length(bg_genes), length(lin_spec), length(div_reg)))

go_group <- function(sig, label) {
  out <- run_go(sig, bg_genes, gene2go, ontologies = c("BP", "MF", "CC"))
  if (nrow(out) == 0) return(NULL)
  out[, category := label][]
}
combined <- rbindlist(list(go_group(lin_spec, "spruce_lineage_specific_de"),
                           go_group(div_reg,  "diverged_regulation")), fill = TRUE)
fwrite(combined, file.path(OUTDIR, "go_specific_subsets.tsv"), sep = "\t")
cat(sprintf("Saved go_specific_subsets.tsv (%d rows; %d plant-consistent)\n",
            nrow(combined), sum(combined$plant_consistent)))

for (nm in unique(combined$category)) {
  sub <- combined[category == nm & ont == "BP" & plant_consistent][order(p)][1:min(8, .N)]
  if (nrow(sub) == 0 || is.na(sub$GO.ID[1])) next
  cat(sprintf("\n-- %s: top plant-consistent BP --\n", nm))
  print(sub[, .(GO.ID, Term, Significant, p = signif(p, 2))], row.names = FALSE)
}
