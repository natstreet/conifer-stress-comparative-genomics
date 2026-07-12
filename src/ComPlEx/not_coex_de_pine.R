#!/usr/bin/env Rscript
# not_coex_de_pine.R — reciprocal of not_coex_de_analysis.R with SCOTS PINE as the focal
# species. Among pine genes DE in >=1 stress tissue but lacking a conserved cross-species
# co-expressolog, classify each by why (mirrors the spruce logic/thresholds exactly, with the
# focal/counterpart species swapped):
#   no_1to1_ortholog           lineage-specific / SD-expanded: no 1:1 spruce ortholog
#   spruce_orth_not_expressed  spruce ortholog exists but not in the spruce expression universe
#   spruce_not_DE              spruce ortholog expressed but not stress-responsive
#   diverged_regulation        both species stress-DE but the co-expression neighbourhood diverged
# Also: YN00 dN/dS of the diverged set vs the rest of pine not_coex, and a GO enrichment of the
# pine lineage-specific (no-ortholog) stress-responsive set.
# Outputs: results/integration/not_coex_de_{mechanism,genes,dnds}_PINE.tsv and .../go/not_coex_de_go_PINE.tsv
suppressPackageStartupMessages({ library(data.table); library(stats) })

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
PROJ_DIR    <- normalizePath(file.path(SCRIPT_DIR, "../.."))
EXPR <- file.path(PROJ_DIR, "data/expression"); DEG <- file.path(PROJ_DIR, "data/DEG_lists")
RDATA <- file.path(PROJ_DIR, "results/ComPlEx/RData"); INT <- file.path(PROJ_DIR, "results/integration")
dir.create(file.path(INT, "go"), showWarnings = FALSE, recursive = TRUE)

# ── DEG lists ─────────────────────────────────────────────────────────────────
for (f in c("SCN","SCR","SDN","SDR","PCN","PCR","PDN","PDR"))
  load(file.path(DEG, sprintf("DE_all_%s_01_2L2FC.RData", f)))
spruce_de_any <- unique(c(DE_all_SCN_01_2L2FC, DE_all_SCR_01_2L2FC, DE_all_SDN_01_2L2FC, DE_all_SDR_01_2L2FC))
pine_de <- data.table(ps_gene = unique(c(DE_all_PCN_01_2L2FC, DE_all_PCR_01_2L2FC, DE_all_PDN_01_2L2FC, DE_all_PDR_01_2L2FC)))
pine_de[, cold_needle_DE    := ps_gene %in% DE_all_PCN_01_2L2FC]
pine_de[, cold_root_DE      := ps_gene %in% DE_all_PCR_01_2L2FC]
pine_de[, drought_needle_DE := ps_gene %in% DE_all_PDN_01_2L2FC]
pine_de[, drought_root_DE   := ps_gene %in% DE_all_PDR_01_2L2FC]
pine_de[, stress_type := fifelse((cold_needle_DE|cold_root_DE)&(drought_needle_DE|drought_root_DE), "both",
                          fifelse(cold_needle_DE|cold_root_DE, "cold_only", "drought_only"))]
cat(sprintf("Pine DE genes (any tissue): %d ; spruce DE genes: %d\n", nrow(pine_de), length(spruce_de_any)))

# ── expression universes ──────────────────────────────────────────────────────
pine_universe   <- unique(c(fread(file.path(EXPR,"PC_expression.txt"),select=1)[[1]], fread(file.path(EXPR,"PD_expression.txt"),select=1)[[1]]))
spruce_universe <- unique(c(fread(file.path(EXPR,"SC_expression.txt"),select=1)[[1]], fread(file.path(EXPR,"SD_expression.txt"),select=1)[[1]]))

# ── co-expressolog category per PINE gene ─────────────────────────────────────
load(file.path(RDATA, "weighted_gene_pairs.RData"))
wp <- as.data.table(weighted_gene_pairs); setnames(wp, c("Species1","Species2"), c("pa_gene","ps_gene"))
wp_best <- wp[order(best_pval)][!duplicated(paste(pa_gene, ps_gene))]
wp_gene_cat <- wp_best[order(best_pval)][!duplicated(ps_gene), .(ps_gene,
  coex_category = fifelse(conserved,"conserved", fifelse(cold_specific,"cold_specific",
                   fifelse(drought_specific,"drought_specific","multi_tissue"))))]
uni <- merge(data.table(ps_gene = pine_universe), wp_gene_cat, by="ps_gene", all.x=TRUE)
uni[is.na(coex_category), coex_category := "not_coex"]

# ── 1:1 ortholog backbone ─────────────────────────────────────────────────────
backbone <- fread(file.path(PROJ_DIR, "interspecies_aa_identity.tsv"))
de <- merge(pine_de, uni, by="ps_gene", all.x=TRUE); de[is.na(coex_category), coex_category := "not_coex"]
de[, has_1to1_orth := ps_gene %in% backbone$ps_gene]
de <- merge(de, backbone[, .(ps_gene, pa_gene_orth = pa_gene)], by="ps_gene", all.x=TRUE)
de[, spruce_orth_DE_any    := pa_gene_orth %in% spruce_de_any]
de[, spruce_orth_expressed := pa_gene_orth %in% spruce_universe]

# ── classify ──────────────────────────────────────────────────────────────────
nc <- de[coex_category == "not_coex"]
nc[, mechanism := fifelse(!has_1to1_orth, "no_1to1_ortholog",
             fifelse(!spruce_orth_expressed, "spruce_orth_not_expressed",
             fifelse(!spruce_orth_DE_any, "spruce_not_DE", "diverged_regulation")))]
cat(sprintf("\nPine not_coex stress-DE genes: %d\n", nrow(nc)))
mech <- nc[, .(n=.N, pct=round(100*.N/nrow(nc),1)), by=mechanism][order(-n)]
print(mech); fwrite(mech, file.path(INT,"not_coex_de_mechanism_PINE.tsv"), sep="\t")
fwrite(nc, file.path(INT,"not_coex_de_genes_PINE.tsv"), sep="\t")

# ── diverged-set YN00 dN/dS vs rest of pine not_coex ──────────────────────────
yn00 <- fread(file.path(INT,"cross_species_dnds_yn00.tsv"))[, .(pa_gene, ps_gene, dNdS)]
nc_pairs <- merge(nc[has_1to1_orth==TRUE, .(ps_gene, pa_gene=pa_gene_orth, mechanism)], yn00, by=c("pa_gene","ps_gene"), all.x=TRUE)
nc_pairs[, sub_cat := fifelse(mechanism=="diverged_regulation","not_coex_both_DE","not_coex_other")]
dnds <- nc_pairs[is.finite(dNdS), .(n=.N, median_dNdS=round(median(dNdS),4)), by=sub_cat][order(sub_cat)]
cat("\nYN00 dN/dS by sub-category (pine not_coex with a 1:1 ortholog):\n"); print(dnds)
fwrite(dnds, file.path(INT,"not_coex_de_dnds_PINE.tsv"), sep="\t")
g1 <- nc_pairs[sub_cat=="not_coex_both_DE" & is.finite(dNdS), dNdS]
g2 <- nc_pairs[sub_cat=="not_coex_other"   & is.finite(dNdS), dNdS]
if (length(g1)>=5 && length(g2)>=5) { wt <- wilcox.test(g1,g2)
  cat(sprintf("  Wilcoxon: diverged (n=%d, median=%.4f) vs rest (n=%d, median=%.4f): W=%.0f p=%.4f\n",
              length(g1),median(g1),length(g2),median(g2),wt$statistic,wt$p.value)) }

# ── GO of pine lineage-specific (no-ortholog) stress-responsive set ───────────
source(file.path(PROJ_DIR, "src/lib/go_enrichment.R"))
gene2go <- build_gene2go(file.path(PROJ_DIR, "data/annotation/pine_GO_annotation.tsv.gz"))
sig <- nc[mechanism=="no_1to1_ortholog", ps_gene]
bg  <- intersect(pine_universe, names(gene2go))
go  <- run_go(sig, bg, gene2go, ontologies="BP")
go_pc <- go[plant_consistent==TRUE][order(p)]
fwrite(go_pc, file.path(INT,"go","not_coex_de_go_PINE.tsv"), sep="\t")
cat(sprintf("\nGO of lineage-specific set (n=%d): top plant-consistent BP terms:\n", length(sig)))
print(head(go_pc[, .(GO.ID, Term, Significant, p)], 8))
cat("\nDone.\n")
