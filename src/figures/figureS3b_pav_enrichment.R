#!/usr/bin/env Rscript
# figureS3b_pav_enrichment.R — Supplementary Figure S3b
#
# Odds ratio for enrichment of presence-absence-variable (PAV) genes from the
# Picea abies standing population across the five co-expression conservation
# categories (Fisher's exact test, Benjamini-Hochberg correction). Points show
# the odds ratio and 95% confidence interval.
#
# The committed per-category Fisher summary (popgen_category_fisher.tsv) stores
# the odds ratio and adjusted p-value but not the confidence interval; the 2x2
# table is reconstructed here from the category background rate and the total
# non-PAV universe so that the interval can be drawn. The reconstructed OR is
# checked against the committed value (tolerance 0.02).
#
# Input : results/integration/popgen_category_fisher.tsv, popgen_category_enrichment.tsv
# Output: results/integration/figures/fig_pav_categories.pdf
# Run from the AbioticStressConifers directory.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("src/lib/fig_palette.R")

# Run from the AbioticStressConifers/ project root (standalone or sourced by assemble_figures.R);
# getwd() is the project root in both cases. (commandArgs("--file=") is unreliable under source().)
PROJ_DIR    <- normalizePath(getwd())
INTEG       <- file.path(PROJ_DIR, "results/integration")
FIGS        <- file.path(INTEG, "figures")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

CAT_ORD <- c("conserved", "cold_specific", "drought_specific", "multi_tissue", "not_coex")
CAT_LAB <- c(conserved = "Conserved", cold_specific = "Cold-specific",
             drought_specific = "Drought-specific", multi_tissue = "Multi-tissue",
             not_coex = "Not co-expressed")

fish <- fread(file.path(INTEG, "popgen_category_fisher.tsv"))[signal == "pav"]
enr  <- fread(file.path(INTEG, "popgen_category_enrichment.tsv"))
# Significance (asterisks) comes from the ADOPTED 5-test family (PAV x 5 co-expression categories, BH),
# committed in pav_category_bh.tsv — NOT popgen_category_fisher.tsv's 15-test (gwas+pav+selection) padj,
# which over-corrects and hid the not_coex result that the main text reports (p_bh = 0.032, survives).
bh   <- fread(file.path(INTEG, "pav_category_bh.tsv"))
n_pav <- enr[label == "pav", n_group]     # total PAV genes (95)
n_bg  <- enr[label == "pav", n_bg]        # non-PAV universe (42956)

# Reconstruct the 2x2 per category and recompute OR + 95% CI
res <- rbindlist(lapply(CAT_ORD, function(cat) {
  row <- fish[coex_category == cat]
  a <- row$n_sig_in_cat                    # PAV in category
  b <- n_pav - a                           # PAV not in category
  c <- round(row$rate_bg * n_bg)           # non-PAV in category
  d <- n_bg - c                            # non-PAV not in category
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  data.table(category = cat, n_pav_in_cat = a,
             OR = unname(ft$estimate), ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2],
             OR_committed = row$OR, padj = bh[coex_category == cat, p_bh])   # 5-test-family BH p
}))
res[, delta := abs(OR - OR_committed)]
cat("Reconstructed vs committed OR (max delta ",
    sprintf("%.3f", max(res$delta)), "):\n", sep = "")
print(res[, .(category, n_pav_in_cat, OR = round(OR, 3),
              OR_committed, ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi), padj = signif(padj, 2))])
if (max(res$delta) > 0.02)
  warning("Reconstructed OR deviates from committed by > 0.02 — check 2x2 reconstruction")

res[, label := factor(CAT_LAB[category], levels = CAT_LAB[CAT_ORD])]
res[, sig := ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**",
             ifelse(padj < 0.05, "*", "")))]

p <- ggplot(res, aes(x = label, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.5) +
  geom_point(size = 2.6, colour = "#D6604D") +
  geom_text(aes(label = sig, y = ci_hi), vjust = -0.3, size = 4) +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.18))) +   # headroom so the top asterisk isn't clipped
  labs(x = NULL, y = "Odds ratio (PAV enrichment, log scale)",
       title = "PAV enrichment across co-expression categories",
       subtitle = "Fisher's exact test, BH-adjusted (5 categories)\n* padj<0.05  ** <0.01  *** <0.001") +
  theme_paper(base_size = 10, major_y = FALSE) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(FIGS, "fig_pav_categories.pdf"), p, width = 12, height = 9, units = "cm")
cat("Saved fig_pav_categories.pdf\n")
