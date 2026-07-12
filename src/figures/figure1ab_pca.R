#!/usr/bin/env Rscript
# Figure 1 a-d: PCA of abiotic-stress expression, one panel per species x stress.
#
# Each panel is a PCA of BOTH tissues (needle + root) together, coloured by
# condition and shaped by tissue (circle = needle, triangle = root). Following the
# DESeq2 plotPCA convention (top-500 most variable genes, centred and UNSCALED),
# PC1 is the needle-vs-root axis and carries ~78-92% of the variance.
#
#   a  P. abies      - cold       b  P. sylvestris - cold
#   c  P. abies      - drought    d  P. sylvestris - drought
#
# Output: results/fig_pca_panels.pdf / .png  (2x2 grid tagged a-d), composed into
# Figure 1 alongside the DEG-bar panel (e) by src/figures/assemble_figures.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(patchwork)
})

DDIR  <- Sys.getenv("EXPR_DIR", "data/expression")  # per-condition VST matrices (ComPlExDataPrep.R)
MDIR  <- "doc"
WDIR  <- "results"
dir.create(WDIR, showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

read_expr <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    column_to_rownames("Genes") %>%
    as.matrix()
}

# PCA on the top-N most variable genes, centred and UNSCALED (DESeq2 plotPCA
# convention). On the combined needle+root matrix, PC1 is the tissue axis.
run_pca <- function(mat, ntop = 500) {
  rv   <- matrixStats::rowVars(mat)
  keep <- order(rv, decreasing = TRUE)[seq_len(min(ntop, nrow(mat)))]
  pca  <- prcomp(t(mat[keep, , drop = FALSE]))
  pct  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  list(scores = as.data.frame(pca$x[, 1:2]) %>% rownames_to_column("sample_id"),
       pct = pct)
}

# ── Colour palettes ───────────────────────────────────────────────────────────
# Cold: control = orange, mild cold = light->mid blue, freezing = purples.
cold_pal <- c(
  "20C_0h"    = "#E6550D",
  "5C_10d"    = "#08519C",
  "5C_24h"    = "#6BAED6",
  "5C_3d"     = "#3182BD",
  "5C_6h"     = "#9ECAE1",
  "neg5C_10d" = "#54278F",
  "neg5C_24h" = "#9E9AC8",
  "neg5C_3d"  = "#756BB1",
  "neg5C_6h"  = "#BCBDDC"
)
# Drought: well-watered = green, progressive drying = yellow->brown, recovery = blue.
drought_pal <- c(
  "30%"         = "#FE9929", "30%7d"     = "#EC7014", "40%"        = "#FEC44F",
  "60%"         = "#ADDD8E", "80%"       = "#2CA25F", "C2d"        = "#8C2D04",
  "Collapse"    = "#CC4C02", "Rehydrate" = "#2166AC",
  "Collapsed"   = "#CC4C02", "Collapsed2d" = "#8C2D04", "FC30"     = "#FE9929",
  "FC30d7"      = "#EC7014", "FC40"      = "#FEC44F", "FC60"       = "#ADDD8E",
  "FC80"        = "#2CA25F", "Rehydrated" = "#2166AC"
)

# ── Load + annotate metadata (maps each sample -> condition + tissue) ──────────

# --- Spruce cold ---
scn_meta <- read_csv(file.path(MDIR, "spruce_cold_needles.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = paste0("SCN", SampleID), tissue = "needle")
scr_meta <- read_csv(file.path(MDIR, "spruce_cold_roots.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = paste0("SCR", SampleID), tissue = "root")

# --- Spruce drought ---
sdn_meta <- read_csv(file.path(MDIR, "spruce_drought_needles.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = paste0("SDN", SciLifeID), Condition = Level, tissue = "needle")
sdr_meta <- read_csv(file.path(MDIR, "spruce_drought_roots.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = paste0("SDR", SciLifeID), Condition = Level, tissue = "root")

# --- Pine cold (single metadata file for both tissues) ---
pine_cold_meta <- read_csv(file.path(MDIR, "pine_cold_stress.csv"), show_col_types = FALSE) %>%
  mutate(tissue = tolower(Tissue),
         sample_id = if_else(tissue == "needle", paste0("PCN", SampleID), paste0("PCR", SampleID)))
pcn_meta <- pine_cold_meta %>% filter(tissue == "needle")
pcr_meta <- pine_cold_meta %>% filter(tissue == "root")

# --- Pine drought (single metadata file, sample names match expression files) ---
pine_drought_meta <- read_csv(file.path(MDIR, "pine_drought_samples.csv"), show_col_types = FALSE) %>%
  mutate(tissue = tolower(Tissue), sample_id = SampleID)
pdn_meta <- pine_drought_meta %>% filter(tissue == "needle")
pdr_meta <- pine_drought_meta %>% filter(tissue == "root")

# ── Read per-condition expression matrices ────────────────────────────────────
cat("Reading expression matrices...\n")
scn_mat <- read_expr(file.path(DDIR, "SCN_expression.txt"))
scr_mat <- read_expr(file.path(DDIR, "SCR_expression.txt"))
sdn_mat <- read_expr(file.path(DDIR, "SDN_expression.txt"))
sdr_mat <- read_expr(file.path(DDIR, "SDR_expression.txt"))
pcn_mat <- read_expr(file.path(DDIR, "PCN_expression.txt"))
pcr_mat <- read_expr(file.path(DDIR, "PCR_expression.txt"))
pdn_mat <- read_expr(file.path(DDIR, "PDN_expression.txt"))
pdr_mat <- read_expr(file.path(DDIR, "PDR_expression.txt"))

# PCR115 is a failed library (99.4% zeros, total counts 788 vs median 32,615)
pcr_mat  <- pcr_mat[, colnames(pcr_mat) != "PCR115"]
pcr_meta <- pcr_meta %>% filter(sample_id != "PCR115")

# T0N3_1 and T0N3_2 are technical replicates of the same biological sample - merge by averaging
pdn_mat  <- cbind(pdn_mat[, !colnames(pdn_mat) %in% c("T0N3_1","T0N3_2")],
                  T0N3 = rowMeans(pdn_mat[, c("T0N3_1","T0N3_2")]))
pdn_meta <- pdn_meta %>%
  filter(!sample_id %in% c("T0N3_1","T0N3_2")) %>%
  bind_rows(tibble(SampleID="T0N3", Tissue="Needle", Condition="FC80", tissue="needle", sample_id="T0N3"))

# ── One panel = combined needle+root PCA, colour = condition, shape = tissue ───
panel_pca <- function(mat_n, mat_r, meta_n, meta_r, pal, species, stress, flip_pc2 = FALSE) {
  shared <- intersect(rownames(mat_n), rownames(mat_r))
  res    <- run_pca(cbind(mat_n[shared, ], mat_r[shared, ]))
  PCA_VAR[[length(PCA_VAR) + 1]] <<- data.frame(                       # record PC1/PC2 variance explained
    dataset = paste(species, stress), PC = c("PC1", "PC2"),
    variance_explained_pct = as.numeric(res$pct[1:2]))
  meta   <- bind_rows(meta_n %>% select(sample_id, Condition, tissue),
                      meta_r %>% select(sample_id, Condition, tissue))
  df <- left_join(res$scores, meta, by = "sample_id")
  # PC-axis sign is arbitrary; flip PC2 so the mild/control end sits low, matching the other panels.
  if (flip_pc2) df$PC2 <- -df$PC2
  df$Condition <- factor(df$Condition, levels = names(pal)[names(pal) %in% df$Condition])
  df$tissue    <- factor(df$tissue, levels = c("needle", "root"))
  ggplot(df, aes(PC1, PC2, colour = Condition, shape = tissue)) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_colour_manual(values = pal, drop = TRUE) +
    scale_shape_manual(values = c(needle = 16, root = 17), name = "Tissue") +
    labs(title = substitute(italic(sp) ~ ds, list(sp = species, ds = paste("-", stress))),
         x = sprintf("PC1 (%.1f%%)", res$pct[1]),
         y = sprintf("PC2 (%.1f%%)", res$pct[2]),
         colour = "Condition") +
    guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 8))
}

cat("Building the four PCA panels...\n")
PCA_VAR <- list()
pa <- panel_pca(scn_mat, scr_mat, scn_meta, scr_meta, cold_pal,    "P. abies",      "cold")
pb <- panel_pca(pcn_mat, pcr_mat, pcn_meta, pcr_meta, cold_pal,    "P. sylvestris", "cold")
pc <- panel_pca(sdn_mat, sdr_mat, sdn_meta, sdr_meta, drought_pal, "P. abies",      "drought")
pd <- panel_pca(pdn_mat, pdr_mat, pdn_meta, pdr_meta, drought_pal, "P. sylvestris", "drought", flip_pc2 = TRUE)

dir.create("results/integration", showWarnings = FALSE, recursive = TRUE)
write.table(do.call(rbind, PCA_VAR), "results/integration/fig1_pca_variance.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote results/integration/fig1_pca_variance.tsv\n")

fig <- (pa | pb) / (pc | pd) + plot_annotation(tag_levels = "a")

ggsave(file.path(WDIR, "fig_pca_panels.pdf"), fig, width = 11, height = 8, device = "pdf")
ggsave(file.path(WDIR, "fig_pca_panels.png"), fig, width = 11, height = 8, dpi = 150)
cat("Saved results/fig_pca_panels.pdf/.png (Figure 1 a-d)\n")
