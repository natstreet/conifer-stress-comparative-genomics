#!/usr/bin/env Rscript
# figureS5_integrative_model.R — Supplementary Figure S5
#
# Three-method comparison of the standardised effects of genomic/evolutionary and
# co-expression-intrinsic features on co-expression conservation breadth, from the
# two complementary integrative models (integrative_conservation_model.R):
#   (a) proportional-odds ordinal regression — odds ratios, 95% CI  [primary]
#   (b) standardised linear regression — coefficients, 95% CI       [robustness]
#   (c) random-forest permutation importance (500 trees, seed 42)   [robustness]
# The text states the ordinal estimates are "complemented by linear regression and
# random-forest importance"; this figure shows all three so that corroboration is
# visible rather than asserted. Panels share the feature ordering; a predictor is
# supported where its ordinal OR / linear coefficient CI excludes the null (1 / 0)
# and it carries appreciable permutation importance.
#
# Run integrative_conservation_model.R first. Run from the AbioticStressConifers dir.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("src/lib/fig_palette.R")
args <- commandArgs(trailingOnly = FALSE)
sp0  <- sub("--file=", "", args[grep("--file=", args)])
PROJ <- normalizePath(file.path(if (length(sp0) > 0) dirname(normalizePath(sp0)) else getwd(), "../.."))
I    <- file.path(PROJ, "results/integration")
OUT  <- file.path(PROJ, "results/final_figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

m   <- fread(file.path(I, "integrative_conservation_model.tsv"))
imp <- fread(file.path(I, "integrative_conservation_model_importance.tsv"))
m   <- merge(m, imp, by = c("model", "term"), all.x = TRUE)

LAB <- c(zdeg = "Network degree", zexpr = "Expression level",
         shared_SD = "Shared segmental duplication", spruce_SD = "Lineage-specific duplication",
         pav = "Presence-absence variation", zdnds = "dN/dS (sequence constraint)",
         zte = "Promoter-TE divergence")
feat_levels <- rev(c("dN/dS (sequence constraint)", "Shared segmental duplication",
                     "Lineage-specific duplication", "Presence-absence variation",
                     "Promoter-TE divergence", "Network degree", "Expression level"))
m[, feature := factor(LAB[term], levels = feat_levels)]
m[, panel := factor(ifelse(grepl("^A", model), "Genome-wide features\n(all co-expression genes)",
                           "Sequence-constraint features\n(1:1 orthologues)"))]
# Fixed structural grouping (not a per-method significance claim): co-expression-intrinsic
# (degree, expression) vs genomic/evolutionary predictor. Significance is read from the CIs.
m[, kind := fifelse(term %in% c("zdeg", "zexpr"), "Co-expression-intrinsic",
                    "Genomic / evolutionary predictor")]
KIND_COL <- c("Co-expression-intrinsic" = "#2166AC", "Genomic / evolutionary predictor" = "#B2182B")

base <- function(p) p +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = KIND_COL, name = NULL) +
  theme_paper(base_size = 10, major_y = FALSE) +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, size = 8, face = "bold"),
        plot.title = element_text(size = 10, face = "bold"))

# (a) ordinal odds ratios ------------------------------------------------------
pa <- base(ggplot(m, aes(odds_ratio, feature, colour = kind)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.22, linewidth = 0.6) +
  geom_point(size = 2.6) +
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  labs(x = "Odds ratio (per SD, 95% CI)", y = NULL, title = "a  Ordinal regression"))

# (b) standardised linear coefficients -----------------------------------------
pb <- base(ggplot(m, aes(lm_std_beta, feature, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lm_ci_lo, xmax = lm_ci_hi), orientation = "y", width = 0.22, linewidth = 0.6) +
  geom_point(size = 2.6) +
  labs(x = "Standardised coefficient (95% CI)", y = NULL, title = "b  Linear regression")) +
  theme(axis.text.y = element_blank(), strip.text.y.left = element_blank())

# (c) random-forest permutation importance (mean over 100 seeds, run-to-run SD) —
# drawn as a point-interval (lollipop) so all three panels share one colour aesthetic and one legend.
pc <- base(ggplot(m, aes(rf_importance_mean, feature, colour = kind)) +
  geom_segment(aes(x = 0, xend = rf_importance_mean, yend = feature), linewidth = 0.6) +
  geom_errorbar(aes(xmin = pmax(0, rf_importance_mean - rf_importance_sd),
                    xmax = rf_importance_mean + rf_importance_sd),
                orientation = "y", width = 0.22, linewidth = 0.5, colour = "grey35") +
  geom_point(size = 2.6) +
  labs(x = "Permutation importance\n(mean +/- run-to-run SD, 100 forests)", y = NULL, title = "c  Random forest") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.10)))) +
  theme(axis.text.y = element_blank(), strip.text.y.left = element_blank(),
        plot.margin = margin(5.5, 9, 5.5, 5.5))

fig <- (pa | pb | pc) +
  plot_layout(guides = "collect", widths = c(1.15, 1, 1)) +
  plot_annotation(
    title = "Independent predictors of co-expression conservation breadth",
    subtitle = "Effects mutually adjusted; three complementary models compared on direction, significance and magnitude",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9))) &
  theme(legend.position = "bottom")

ggsave(file.path(OUT, "FigureS5.pdf"), fig, width = 26, height = 12, units = "cm")
ggsave(file.path(OUT, "FigureS5.png"), fig, width = 26, height = 12, units = "cm", dpi = 300)
cat("Saved FigureS5.pdf and FigureS5.png (3-method comparison)\n")
