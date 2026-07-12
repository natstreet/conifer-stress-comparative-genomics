#!/usr/bin/env Rscript
# figure4_conserved_dynamics.R — Figure 4
#
# Reads conserved_coexpression_dynamics.R outputs and builds Figure 4:
#   (a) overlap of the drought and cold conserved co-expressolog sets;
#   (b) cross-species expression-profile conservation vs a shuffled-partner null;
#   (c,d) heatmaps of median expression (Norway spruce | Scots pine, per tissue),
#         rows split by response direction — the visual impression of conserved timing.
# Run conserved_coexpression_dynamics.R first, then this, from AbioticStressConifers/.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(grid)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
  library(ComplexHeatmap); library(circlize); library(cowplot)
})
source("src/lib/fig_palette.R")

# Run from the AbioticStressConifers/ project root, either standalone (Rscript src/ComPlEx/…)
# or sourced by assemble_figures.R; in both cases getwd() is the project root. (Deriving the
# root from commandArgs("--file=") is unreliable under source(), where --file= is the outer script.)
PROJ <- normalizePath(getwd())
EXPR <- file.path(PROJ, "data/expression/expr_median")
F4   <- file.path(PROJ, "results/integration/fig4")
FIGS <- file.path(PROJ, "results/integration/figures"); dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

sets <- fread(file.path(F4, "conserved_coexpressolog_sets.tsv"))
ov   <- fread(file.path(F4, "stress_set_overlap.tsv"))
STAGES <- list(
  drought = c("FC80","FC60","FC40","FC30","FC30d7","Collapsed","Collapsed2d","Rehydrated"),
  cold    = c("20C_0h","5C_6h","5C_24h","5C_3d","5C_10d","neg5C_6h","neg5C_24h","neg5C_3d","neg5C_10d"))
PREF <- list(drought = c(s = "SD", p = "PD"), cold = c(s = "SC", p = "PC"))
read_med <- function(f) { m <- read.table(f, header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)
  rownames(m) <- m[["Genes"]]; m[, -1] }

# ── (a) overlap bar ──────────────────────────────────────────────────────────
ov_df <- data.frame(part = factor(c("Drought-specific","Shared (both stresses)","Cold-specific"),
                                  levels = c("Drought-specific","Shared (both stresses)","Cold-specific")),
                    n = c(ov$drought_only, ov$shared_both, ov$cold_only))
p_ov <- ggplot(ov_df, aes(x = "", y = n, fill = part)) +
  geom_col(width = 0.5, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.2) +
  scale_fill_manual(values = c("Drought-specific" = PAL$drought, "Shared (both stresses)" = PAL$stress_shared,
                               "Cold-specific" = PAL$cold), name = NULL,
                    guide = guide_legend(reverse = FALSE)) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Conserved co-expressolog orthogroups", title = "a  Stress-set overlap") +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 11))

# ── (b) profile-conservation distributions (observed vs null) ────────────────
long <- rbind(
  sets[, .(stress, r = profile_cor_needle, set = "Observed")],
  sets[, .(stress, r = profile_cor_null_needle, set = "Shuffled null")])
long <- long[!is.na(r)]
long[, stress := factor(stress, c("drought","cold"), c("Drought","Cold"))]
p_pc <- ggplot(long, aes(x = r, fill = set)) +
  geom_density(alpha = 0.55, colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  facet_wrap(~stress, ncol = 1) +
  scale_fill_manual(values = c("Observed" = "#2C7FB8", "Shuffled null" = "grey65"), name = NULL) +
  labs(x = "Spruce-pine expression-profile correlation (needle)", y = "Density",
       title = "b  Cross-species profile conservation") +
  theme_paper(base_size = 10) + theme(legend.position = "bottom")

# ── (c,d) heatmaps: spruce | pine per tissue, rows split by direction ────────
make_hm <- function(stress) {
  cfg_stg <- STAGES[[stress]]; pr <- PREF[[stress]]
  d <- sets[stress == get("stress") & !is.na(direction)]
  d <- d[order(direction, -profile_cor_needle)]
  S <- read_med(file.path(EXPR, paste0(pr["s"], "med_expression.txt")))
  P <- read_med(file.path(EXPR, paste0(pr["p"], "med_expression.txt")))
  gv <- function(mat, pre, g, ti) { c <- paste0(pre, "_", cfg_stg, "_", ti)
    if (!g %in% rownames(mat)) return(rep(NA, length(cfg_stg))); suppressWarnings(as.numeric(mat[g, c])) }
  mat <- t(sapply(seq_len(nrow(d)), function(i) c(
    gv(S, pr["s"], d$pa_gene[i], "N"), gv(P, pr["p"], d$ps_gene[i], "N"),
    gv(S, pr["s"], d$pa_gene[i], "R"), gv(P, pr["p"], d$ps_gene[i], "R"))))
  n <- length(cfg_stg)
  # z-score WITHIN each species x tissue block, per orthogroup, so the panels show the
  # temporal response shape (blue early -> red late for up-regulated) rather than
  # baseline expression-level differences between tissues/species.
  zscore_block <- function(x) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x))); (x - m) / s }
  for (b in list(1:n, (n+1):(2*n), (2*n+1):(3*n), (3*n+1):(4*n)))
    mat[, b] <- t(apply(mat[, b, drop = FALSE], 1, zscore_block))
  mat[!is.finite(mat)] <- 0; mat <- pmin(pmax(mat, -2), 2)
  colsplit <- factor(rep(c("Spruce needle","Pine needle","Spruce root","Pine root"), each = n),
                     levels = c("Spruce needle","Pine needle","Spruce root","Pine root"))
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC","#F7F7F7","#B2182B"))
  Heatmap(mat, name = "Row z", col = col_fun,
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
          row_split = factor(d$direction, c("up","down"), c("Up-regulated","Down-regulated")),
          column_split = colsplit, row_title_gp = gpar(fontsize = 8, fontface = "bold"),
          column_title_gp = gpar(fontsize = 8), column_gap = unit(1.5, "mm"), border = TRUE,
          heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)),
          use_raster = FALSE)
}
hm_d <- grid.grabExpr(draw(make_hm("drought"), column_title = "c  Drought", column_title_gp = gpar(fontsize = 10, fontface = "bold")))
hm_c <- grid.grabExpr(draw(make_hm("cold"),    column_title = "d  Cold",    column_title_gp = gpar(fontsize = 10, fontface = "bold")))

top    <- plot_grid(p_ov, p_pc, ncol = 2, rel_widths = c(1, 1))
bottom <- plot_grid(hm_d, hm_c, ncol = 2)
fig    <- plot_grid(top, bottom, ncol = 1, rel_heights = c(0.9, 1.3))
ggsave(file.path(FIGS, "Figure4_conserved_dynamics.pdf"), fig, width = 24, height = 22, units = "cm")
ggsave(file.path(FIGS, "Figure4_conserved_dynamics.png"), fig, width = 24, height = 22, units = "cm", dpi = 300)
cat("Saved Figure4_conserved_dynamics.pdf/.png to", FIGS, "\n")
