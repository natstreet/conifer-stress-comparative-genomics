#!/usr/bin/env Rscript
# table1_2_conserved_dynamics_go.R — main-text Table 1 (drought) + Table 2 (cold)
#
# GO enrichment of the up- and down-regulated conserved co-expressologs
# (from conserved_coexpression_dynamics.R), per stress x direction. Uses the single
# shared GO method in src/lib/go_enrichment.R (topGO weight01/Fisher, nodeSize 5, NO
# multiple-testing correction — weight01 already accounts for GO-graph dependence —
# all terms at p<0.05 across BP/MF/CC, each with a plant_consistent flag).
# The full table is written for transparency; Table 2 (main text) shows the BP,
# plant-consistent terms (printed below).
#
# Run conserved_coexpression_dynamics.R first. Run from AbioticStressConifers/.
# Output: results/integration/fig4/table2_conserved_dynamics_go.tsv
suppressPackageStartupMessages({ library(data.table) })
source("src/lib/go_enrichment.R")

EXPRU <- "data/expression"
F4    <- "results/integration/fig4"

gene2go <- build_gene2go("data/annotation/gene_annotation.tsv.gz")                    # spruce eggNOG + InterPro
univ <- unique(c(fread(file.path(EXPRU, "SC_expression.txt"), select = 1)[[1]],
                 fread(file.path(EXPRU, "SD_expression.txt"), select = 1)[[1]]))
bg <- intersect(univ, names(gene2go))                             # expressed & annotated spruce genes
cat(sprintf("Background: %d spruce genes\n", length(bg)))

sets   <- fread(file.path(F4, "conserved_coexpressolog_sets.tsv"))
combos <- CJ(stress = c("drought", "cold"), direction = c("up", "down"))
res <- rbindlist(lapply(seq_len(nrow(combos)), function(i) {
  st <- combos$stress[i]; dr <- combos$direction[i]
  genes <- sets[stress == st & direction == dr, pa_gene]
  out <- run_go(genes, bg, gene2go, ontologies = c("BP", "MF", "CC"))
  if (nrow(out) == 0) return(NULL)
  out[, `:=`(stress = st, direction = dr, n_genes = length(intersect(genes, bg)))]
  out
}), fill = TRUE)

setcolorder(res, c("stress", "direction", "ont", "GO.ID", "Term", "Annotated", "Significant",
                   "Expected", "p", "plant_consistent", "plant_consistency", "n_genes"))
fwrite(res, file.path(F4, "table2_conserved_dynamics_go.tsv"), sep = "\t")
cat(sprintf("\nSaved table2_conserved_dynamics_go.tsv (%d rows; %d plant-consistent)\n",
            nrow(res), sum(res$plant_consistent)))

# Table 2 (main text) = plant-consistent BP terms per stress x direction
for (st in c("drought", "cold")) for (dr in c("up", "down")) {
  r <- res[stress == st & direction == dr & ont == "BP" & plant_consistent][order(p)]
  cat(sprintf("\n-- %s %s-regulated (%s genes): plant-consistent BP terms --\n",
              st, dr, if (nrow(r)) r$n_genes[1] else "NA"))
  if (nrow(r) == 0) { cat("  (none)\n"); next }
  print(head(r[, .(Term, Significant, p = signif(p, 2))], 8), row.names = FALSE)
}
