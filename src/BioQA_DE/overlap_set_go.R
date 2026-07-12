#!/usr/bin/env Rscript
# overlap_set_go.R
#
# GO enrichment of the Figure 2 DEG-overlap sets, per species (Results-text GO:
# the shared-core / cold-specific / drought-specific terms quoted in the Results,
# e.g. "energy reserve metabolic process" in the Norway spruce shared core).
# Overlap sets are the same gene-level operations as build_deg_overlap_counts.R:
#   cold_specific    = (cold-needle INT cold-root)    \ (drought-needle U drought-root)
#   drought_specific = (drought-needle INT drought-root) \ (cold-needle U cold-root)
#   shared_core      = cold-needle INT cold-root INT drought-needle INT drought-root
# Uses the single shared GO method (src/lib/go_enrichment.R): topGO weight01/Fisher,
# nodeSize 5, no BH, BP/MF/CC, all terms with a plant_consistent flag. gene->GO from
# eggNOG+InterPro: spruce = data/annotation/gene_annotation.tsv.gz, pine = pine_GO_annotation.tsv.gz;
# background = each species' expression universe INT annotated.
#
# Run from AbioticStressConifers/.  Output: results/integration/go/go_overlap_sets_by_species.tsv
suppressPackageStartupMessages({ library(data.table); library(here) })
source("src/lib/go_enrichment.R")

DEGL      <- "data/DEG_lists"
EXPR_DIR  <- "data/expression"
OUTDIR    <- "results/integration/go"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
SPRUCE_ANNOT <- here("data/annotation/gene_annotation.tsv.gz")
PINE_ANNOT   <- "data/annotation/pine_GO_annotation.tsv.gz"

ld <- function(cond) {   # DE gene set for one condition
  e <- new.env(); load(file.path(DEGL, paste0("DE_all_", cond, "_01_2L2FC.RData")), e)
  get(ls(e)[grep("DE_all", ls(e))], e)
}
overlap_sets <- function(CN, CR, DN, DR) list(
  cold_specific    = setdiff(intersect(CN, CR), union(DN, DR)),
  drought_specific = setdiff(intersect(DN, DR), union(CN, CR)),
  shared_core      = Reduce(intersect, list(CN, CR, DN, DR)))

run_species <- function(species, cond_prefix, expr_files, annot_file) {
  cat(sprintf("\n== %s ==\n", species))
  gene2go <- build_gene2go(annot_file)
  expr_universe <- unique(unlist(lapply(expr_files, function(f) fread(file.path(EXPR_DIR, f), select = 1)[[1]])))
  bg <- intersect(expr_universe, names(gene2go))
  cat(sprintf("  background (expression universe INT annotated): %d\n", length(bg)))
  sets <- overlap_sets(ld(paste0(cond_prefix, "CN")), ld(paste0(cond_prefix, "CR")),
                       ld(paste0(cond_prefix, "DN")), ld(paste0(cond_prefix, "DR")))
  rbindlist(lapply(names(sets), function(sn) {
    cat(sprintf("  %-16s: %d genes\n", sn, length(sets[[sn]])))
    out <- run_go(sets[[sn]], bg, gene2go, ontologies = c("BP", "MF", "CC"))
    if (nrow(out) == 0) return(NULL)
    out[, `:=`(overlap_set = sn, species = species)][]
  }), fill = TRUE)
}

spruce <- run_species("Picea_abies",     "S", c("SC_expression.txt", "SD_expression.txt"), SPRUCE_ANNOT)
pine   <- run_species("Pinus_sylvestris", "P", c("PC_expression.txt", "PD_expression.txt"), PINE_ANNOT)
combined <- rbindlist(list(spruce, pine), fill = TRUE)
fwrite(combined, file.path(OUTDIR, "go_overlap_sets_by_species.tsv"), sep = "\t")
cat(sprintf("\nSaved go_overlap_sets_by_species.tsv (%d rows; %d plant-consistent)\n",
            nrow(combined), sum(combined$plant_consistent)))

show <- function(sp, set) {
  cat(sprintf("\n-- %s / %s (top plant-consistent BP) --\n", sp, set))
  r <- combined[species == sp & overlap_set == set & ont == "BP" & plant_consistent][order(p)]
  if (nrow(r) == 0) { cat("  (none)\n"); return() }
  print(head(r[, .(GO.ID, Term, Significant, p = signif(p, 2))], 10), row.names = FALSE)
}
for (s in c("shared_core", "cold_specific", "drought_specific")) { show("Picea_abies", s); show("Pinus_sylvestris", s) }
cat(sprintf("\nOutput in: %s\n", OUTDIR))
