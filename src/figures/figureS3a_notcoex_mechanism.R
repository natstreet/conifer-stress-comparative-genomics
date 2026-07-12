#!/usr/bin/env Rscript
# figureS3a_notcoex_mechanism.R — Supplementary Figure S3a
#
# Mechanism classification of not_coex stress-responsive genes: the proportion of
# genes that are stress-responsive in at least one species but classified as
# not co-expressed, broken down by mechanism (no 1:1 ortholog; partner not DE;
# diverged regulation; ortholog not expressed). Percentages and counts shown.
#
# Input : results/integration/not_coex_de_genes.tsv (column `mechanism`)
# Output: results/integration/figures/fig_notcoex_mechanism.pdf
# Run from the AbioticStressConifers directory.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("src/lib/fig_palette.R")

# Run from the AbioticStressConifers/ project root (standalone or sourced by assemble_figures.R);
# getwd() is the project root in both cases. (commandArgs("--file=") is unreliable under source().)
PROJ_DIR    <- normalizePath(getwd())
INTEG       <- file.path(PROJ_DIR, "results/integration")
FIGS        <- file.path(INTEG, "figures")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

d <- fread(file.path(INTEG, "not_coex_de_genes.tsv"))

# Human-readable mechanism labels (order = decreasing frequency in the manuscript)
lab <- c(no_1to1_ortholog        = "No 1:1 ortholog",
         pine_not_DE             = "Partner not DE",
         diverged_regulation     = "Diverged regulation",
         pine_orth_not_expressed = "Ortholog not expressed")
cnt <- d[, .N, by = mechanism]
cnt[, label := factor(lab[mechanism], levels = lab)]
cnt[, pct   := 100 * N / sum(N)]
setorder(cnt, -N)
cat("Mechanism breakdown:\n"); print(cnt[, .(mechanism, N, pct = round(pct, 1))])

p <- ggplot(cnt, aes(x = label, y = pct, fill = label)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, N)),
            vjust = -0.2, size = 3) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, max(cnt$pct) * 1.15), expand = c(0, 0)) +
  labs(x = NULL, y = "% of not_coex stress-responsive genes",
       title = "Mechanism classification of not_coex genes") +
  theme_paper(base_size = 10, major_y = TRUE) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(FIGS, "fig_notcoex_mechanism.pdf"), p, width = 12, height = 9, units = "cm")
cat("Saved fig_notcoex_mechanism.pdf\n")
