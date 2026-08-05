#!/usr/bin/env Rscript
# not_coex_de_analysis.R
#
# Among spruce genes that are DE in at least one stress tissue but lack a
# conserved cross-species co-expressolog, two non-exclusive explanations exist:
#
#   (a) Diverged regulation: the gene has a 1:1 pine ortholog and that ortholog
#       IS expressed, but the two copies have diverged in how they are regulated
#       under stress (different co-expression neighbourhood).
#
#   (b) Species-specific response: the spruce gene either has no 1:1 pine
#       ortholog (lineage-specific / SD paralog not represented in pine) or the
#       pine ortholog is not DE in the corresponding stress experiment.
#
# Outputs: results/integration/not_coex_de_*.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
PROJ_DIR    <- normalizePath(file.path(SCRIPT_DIR, "../.."))
SD_DIR      <- PROJ_DIR
EXPR_DIR    <- file.path(PROJ_DIR, "data/expression")
DEG_DIR     <- file.path(PROJ_DIR, "data/DEG_lists")
RDATA_DIR   <- file.path(PROJ_DIR, "results/ComPlEx/RData")
OUTDIR      <- file.path(PROJ_DIR, "results/integration")

# ── Load spruce and pine DEG lists ─────────────────────────────────────────────
cat("Loading DEG lists...\n")
load(file.path(DEG_DIR, "DE_all_SCN_01_2L2FC.RData"))
load(file.path(DEG_DIR, "DE_all_SCR_01_2L2FC.RData"))
load(file.path(DEG_DIR, "DE_all_SDN_01_2L2FC.RData"))
load(file.path(DEG_DIR, "DE_all_SDR_01_2L2FC.RData"))
load(file.path(DEG_DIR, "DE_all_PCN_01_2L2FC.RData"))
load(file.path(DEG_DIR, "DE_all_PCR_01_2L2FC.RData"))
load(file.path(DEG_DIR, "DE_all_PDN_01_2L2FC.RData"))
load(file.path(DEG_DIR, "DE_all_PDR_01_2L2FC.RData"))

# Per-tissue DE flags for spruce
spruce_de <- data.table(
  pa_gene = unique(c(DE_all_SCN_01_2L2FC, DE_all_SCR_01_2L2FC,
                     DE_all_SDN_01_2L2FC, DE_all_SDR_01_2L2FC))
)
spruce_de[, cold_needle_DE   := pa_gene %in% DE_all_SCN_01_2L2FC]
spruce_de[, cold_root_DE     := pa_gene %in% DE_all_SCR_01_2L2FC]
spruce_de[, drought_needle_DE := pa_gene %in% DE_all_SDN_01_2L2FC]
spruce_de[, drought_root_DE  := pa_gene %in% DE_all_SDR_01_2L2FC]
spruce_de[, n_de_tissues := cold_needle_DE + cold_root_DE + drought_needle_DE + drought_root_DE]
spruce_de[, stress_type := fifelse(
  (cold_needle_DE | cold_root_DE) & (drought_needle_DE | drought_root_DE), "both",
  fifelse(cold_needle_DE | cold_root_DE, "cold_only", "drought_only"))]
cat(sprintf("  Spruce DE genes (any tissue): %d\n", nrow(spruce_de)))
cat(sprintf("  Pine DE genes   (any tissue): %d\n",
  length(unique(c(DE_all_PCN_01_2L2FC, DE_all_PCR_01_2L2FC,
                  DE_all_PDN_01_2L2FC, DE_all_PDR_01_2L2FC)))))

# Pine per-tissue DE lookup
pine_de_any  <- unique(c(DE_all_PCN_01_2L2FC, DE_all_PCR_01_2L2FC,
                         DE_all_PDN_01_2L2FC, DE_all_PDR_01_2L2FC))
pine_de_cold <- unique(c(DE_all_PCN_01_2L2FC, DE_all_PCR_01_2L2FC))
pine_de_drt  <- unique(c(DE_all_PDN_01_2L2FC, DE_all_PDR_01_2L2FC))

# ── Load expression universe and co-expressolog categories ─────────────────────
cat("Loading co-expressolog data...\n")
sc_g <- fread(file.path(EXPR_DIR, "SC_expression.txt"), select = 1)[[1]]
sd_g <- fread(file.path(EXPR_DIR, "SD_expression.txt"), select = 1)[[1]]
spruce_universe <- unique(c(sc_g, sd_g))
cat(sprintf("  Spruce expression universe: %d genes\n", length(spruce_universe)))

load(file.path(RDATA_DIR, "weighted_gene_pairs.RData"))
wp <- as.data.table(weighted_gene_pairs)
setnames(wp, c("Species1","Species2"), c("pa_gene","ps_gene"))
wp_best <- wp[order(best_pval)][!duplicated(paste(pa_gene, ps_gene))]

wp_gene_cat <- wp_best[order(best_pval)][!duplicated(pa_gene), .(
  pa_gene,
  coex_category = fifelse(conserved, "conserved",
    fifelse(cold_specific, "cold_specific",
    fifelse(drought_specific, "drought_specific", "multi_tissue"))),
  best_ps_gene = ps_gene,
  best_pval
)]

universe_dt <- data.table(pa_gene = spruce_universe)
universe_dt <- merge(universe_dt, wp_gene_cat, by = "pa_gene", all.x = TRUE)
universe_dt[is.na(coex_category), coex_category := "not_coex"]

# ── Load 1:1 ortholog backbone ─────────────────────────────────────────────────
cat("Loading 1:1 ortholog backbone...\n")
backbone <- fread(file.path(SD_DIR, "interspecies_aa_identity.tsv"))
cat(sprintf("  1:1 HOG orthologs: %d\n", nrow(backbone)))

# ── Build analysis table: spruce DE × ortholog status × co-expressolog ─────────
cat("\nBuilding DE × ortholog × co-expressolog table...\n")

# Focus on DE genes in the expression universe
de_tbl <- merge(spruce_de, universe_dt[, .(pa_gene, coex_category, best_ps_gene)],
                by = "pa_gene", all.x = TRUE)
de_tbl[is.na(coex_category), coex_category := "not_coex"]

# 1:1 ortholog presence: does this spruce DE gene have a pine 1:1 ortholog?
de_tbl[, has_1to1_orth := pa_gene %in% backbone$pa_gene]

# If it has a 1:1 ortholog, what is the pine gene and is it DE?
orth_lookup <- backbone[, .(pa_gene, ps_gene_orth = ps_gene)]
de_tbl <- merge(de_tbl, orth_lookup, by = "pa_gene", all.x = TRUE)
de_tbl[, pine_orth_DE_any   := ps_gene_orth %in% pine_de_any]
de_tbl[, pine_orth_DE_cold  := ps_gene_orth %in% pine_de_cold]
de_tbl[, pine_orth_DE_drt   := ps_gene_orth %in% pine_de_drt]

# Is the pine 1:1 ortholog present in the pine expression universe?
pc_g <- fread(file.path(EXPR_DIR, "PC_expression.txt"), select = 1)[[1]]
pd_g <- fread(file.path(EXPR_DIR, "PD_expression.txt"), select = 1)[[1]]
pine_universe <- unique(c(pc_g, pd_g))
de_tbl[, pine_orth_expressed := ps_gene_orth %in% pine_universe]

cat(sprintf("  Spruce DE genes in expression universe: %d\n", nrow(de_tbl)))

# ── Classify not-coex DE genes ─────────────────────────────────────────────────
not_coex_de <- de_tbl[coex_category == "not_coex"]
cat(sprintf("  Not-coex DE spruce genes: %d\n", nrow(not_coex_de)))

# Classify into mechanistic categories
not_coex_de[, mechanism := fifelse(
  !has_1to1_orth,
  "no_1to1_ortholog",        # SD paralog or lineage-specific: no pine counterpart
  fifelse(
    !pine_orth_expressed,
    "pine_orth_not_expressed", # pine ortholog exists but not in expression data
    fifelse(
      !pine_orth_DE_any,
      "pine_not_DE",            # pine ortholog expressed but not stress-responsive
      "diverged_regulation"     # both DE in stress but co-expression diverged
    )
  )
)]

cat("\nMechanism breakdown for not-coex stress-DE genes:\n")
mech_tbl <- not_coex_de[, .(
  n = .N,
  pct = round(100*.N/nrow(not_coex_de), 1)
), by = mechanism][order(-n)]
print(mech_tbl)
fwrite(mech_tbl, file.path(OUTDIR, "not_coex_de_mechanism.tsv"), sep = "\t")

# Full gene table
fwrite(not_coex_de, file.path(OUTDIR, "not_coex_de_genes.tsv"), sep = "\t")
cat(sprintf("  Saved not_coex_de_genes.tsv (%d rows)\n", nrow(not_coex_de)))

# ── Stress-type breakdown per mechanism ────────────────────────────────────────
mech_stress <- not_coex_de[, .N, by = .(mechanism, stress_type)][order(mechanism, stress_type)]
cat("\nMechanism × stress type:\n")
print(mech_stress)
fwrite(mech_stress, file.path(OUTDIR, "not_coex_de_mechanism_stress.tsv"), sep = "\t")

# ── Compare: how are coex categories distributed among DE genes? ───────────────
cat("\nCo-expressolog category distribution among stress-DE genes:\n")
de_cat <- de_tbl[, .(n = .N, pct = round(100*.N/nrow(de_tbl), 1)), by = coex_category][order(-n)]
print(de_cat)
fwrite(de_cat, file.path(OUTDIR, "de_coex_category_distribution.tsv"), sep = "\t")

# Per-gene DE-set co-expressolog category (the assignment underlying the distribution above),
# committed so downstream tables (e.g. Supplementary Table S8) can join the paper's category per gene.
fwrite(de_tbl[, .(pa_gene, coex_category, n_de_tissues, stress_type)][order(pa_gene)],
       file.path(OUTDIR, "de_gene_coex_category.tsv"), sep = "\t")

# Compare DE rate within each co-expressolog category
# (what fraction of each category's genes are stress-DE?)
cat("\nDE rate by co-expressolog category:\n")
de_rate <- merge(
  universe_dt[, .N, by = coex_category],
  de_tbl[, .(n_de = .N), by = coex_category],
  by = "coex_category", all.x = TRUE)
de_rate[is.na(n_de), n_de := 0]
de_rate[, de_rate := round(n_de/N, 4)]
de_rate <- de_rate[order(coex_category)]
print(de_rate)
fwrite(de_rate, file.path(OUTDIR, "de_rate_by_coex_category.tsv"), sep = "\t")

# Fisher: is each co-expressolog category enriched for DE genes?
n_total <- nrow(universe_dt)
n_de_total <- nrow(de_tbl)
cat("\nFisher enrichment: DE rate per co-expressolog category vs universe:\n")
fisher_results <- rbindlist(lapply(unique(de_rate$coex_category), function(cat) {
  n_cat <- universe_dt[coex_category == cat, .N]
  if (n_cat == 0) return(NULL)
  n_cat_de <- de_tbl[coex_category == cat, .N]
  n_bg_de   <- n_de_total - n_cat_de
  n_bg_not  <- n_total - n_de_total - (n_cat - n_cat_de)
  mat <- matrix(c(n_cat_de, n_cat - n_cat_de, n_bg_de, n_bg_not), nrow = 2)
  ft <- fisher.test(mat)
  data.table(category = cat, n_cat = n_cat, n_de = n_cat_de,
             de_rate = round(n_cat_de/n_cat, 4),
             bg_de_rate = round(n_de_total/n_total, 4),
             OR = round(ft$estimate, 3), pvalue = ft$p.value)
}))
fisher_results[, padj := p.adjust(pvalue, method = "BH")]
print(fisher_results[order(pvalue), .(category, n_cat, n_de, de_rate, bg_de_rate, OR, pvalue, padj)])
fwrite(fisher_results, file.path(OUTDIR, "de_enrichment_by_coex_category.tsv"), sep = "\t")

# ── For diverged-regulation genes: compare YN00 dN/dS vs the rest of not_coex ──
# Uses the Yang-Nielsen YN00 cross-species dN/dS (the estimator used for Figure 5a), so this
# comparison is on the same footing as the main dN/dS analysis.
cat("\n── Diverged-regulation subset: YN00 dN/dS comparison ──────────────────\n")
diverg_genes <- not_coex_de[mechanism == "diverged_regulation", pa_gene]
cat(sprintf("  Diverged-regulation genes (both species stress-DE, not co-expressed): %d\n",
            length(diverg_genes)))

backbone_ks <- fread(file.path(OUTDIR, "integration_backbone_1to1.tsv"))
yn00 <- fread(file.path(OUTDIR, "cross_species_dnds_yn00.tsv"))[, .(pa_gene, ps_gene, dNdS)]
backbone_ks <- merge(backbone_ks, yn00, by = c("pa_gene", "ps_gene"), all.x = TRUE)

# Sub-category: diverged_regulation (in our analysis) vs the rest of not_coex.
backbone_ks[, sub_cat := fifelse(
  coex_category == "not_coex" & pa_gene %in% diverg_genes,
  "not_coex_both_DE",
  coex_category)]

dnds_diverg <- backbone_ks[is.finite(dNdS), .(
  n           = .N,
  median_dNdS = round(median(dNdS), 4)
), by = sub_cat][order(sub_cat)]
cat("  YN00 dN/dS by sub-category (not_coex split by whether pine ortholog is also DE):\n")
print(dnds_diverg)
fwrite(dnds_diverg, file.path(OUTDIR, "not_coex_de_dnds.tsv"), sep = "\t")

if (length(diverg_genes) >= 10) {
  g1 <- backbone_ks[sub_cat == "not_coex_both_DE" & is.finite(dNdS), dNdS]
  g2 <- backbone_ks[sub_cat == "not_coex"         & is.finite(dNdS), dNdS]
  if (length(g1) >= 5 && length(g2) >= 5) {
    wt <- wilcox.test(g1, g2)
    cat(sprintf("  Wilcoxon YN00 dN/dS: diverged (n=%d, median=%.4f) vs rest-of-not_coex (n=%d, median=%.4f): W=%.0f p=%.4f\n",
                length(g1), median(g1), length(g2), median(g2), wt$statistic, wt$p.value))
  }
}

cat(sprintf("\nOutputs in: %s\n", OUTDIR))
