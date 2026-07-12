#!/usr/bin/env Rscript
## figure2c_orthogroup_overlap.R — Figure 2c: cross-species orthogroup-level
## DEG overlap between Norway spruce and Scots pine at each matched sampling point.
##
## For each sampling point (and overall), an orthogroup is "hit" by a species if any
## of that species' pooled DEGs map into it. Bars show spruce-specific / shared /
## pine-specific orthogroup counts. Shared bars are annotated with the percentage of
## the SMALLER species' DEG-orthogroup set that is shared (NOT Jaccard; Jaccard is
## reported in the supplementary per-timepoint table only). Fold-enrichment over
## chance and a hypergeometric P are written to the log.
##
## Universe = orthogroups containing >=1 spruce AND >=1 pine gene in
## doc/genes_ortholog_categories.tsv (the FULL gene-space universe, N = 11,943 —
## NOT restricted to the co-expressolog backbone).
##
## Run from AbioticStressConifers/:  Rscript src/figures/figure2c_orthogroup_overlap.R
## Inputs : doc/genes_ortholog_categories.tsv
##          data/DEG_lists/DE_all_{SCN,SCR,SDN,SDR,PCN,PCR,PDN,PDR}_01_2L2FC.RData
## Outputs: results/final_figures/Figure2c.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2) })
source("src/lib/fig_palette.R")

OGX  <- "doc/genes_ortholog_categories.tsv"
DEGD <- "data/DEG_lists"
OUT  <- "results/final_figures"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- universe + gene->orthogroup maps (per species) -------------------------
x <- read.table(OGX, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
sp_map <- x[x$species == "Picea_abies",     c("gene", "Ortholog_Group")]
pi_map <- x[x$species == "Pinus_sylvestris", c("gene", "Ortholog_Group")]
universe <- intersect(unique(sp_map$Ortholog_Group), unique(pi_map$Ortholog_Group))
N <- length(universe)
g2og_sp <- split(sp_map$Ortholog_Group, sp_map$gene)
g2og_pi <- split(pi_map$Ortholog_Group, pi_map$gene)

degs <- function(tag) {                       # pooled DEG gene IDs for one dataset
  e <- new.env(); load(file.path(DEGD, sprintf("DE_all_%s_01_2L2FC.RData", tag)), e)
  o <- get(ls(e)[1], e); if (is.data.frame(o)) rownames(o) else as.character(o)
}
ogs_of <- function(genes, map)                # DEG orthogroups within the universe
  intersect(unique(unlist(map[genes[genes %in% names(map)]], use.names = FALSE)), universe)

## ---- per-group overlap ------------------------------------------------------
groups <- list(
  cold_needle    = list(s = "SCN", p = "PCN"),
  cold_root      = list(s = "SCR", p = "PCR"),
  drought_needle = list(s = "SDN", p = "PDN"),
  drought_root   = list(s = "SDR", p = "PDR"),
  overall        = list(s = c("SCN","SCR","SDN","SDR"), p = c("PCN","PCR","PDN","PDR")))

rows <- lapply(names(groups), function(g) {
  so <- ogs_of(unique(unlist(lapply(groups[[g]]$s, degs))), g2og_sp)
  po <- ogs_of(unique(unlist(lapply(groups[[g]]$p, degs))), g2og_pi)
  sh <- length(intersect(so, po)); ss <- length(so); pp <- length(po)
  data.frame(group = g, sp_only = length(setdiff(so, po)), shared = sh,
             pi_only = length(setdiff(po, so)),
             jaccard = round(sh / length(union(so, po)), 3),
             pct_smaller = round(100 * sh / min(ss, pp), 1),
             fold = round(sh / (ss * pp / N), 2),
             P = phyper(sh - 1, pp, N - pp, ss, lower.tail = FALSE),
             # log10 of the same one-sided hypergeometric P, computed in log space so the strongest
             # overlaps do not underflow the double P column to 0 (the pooled "overall" row true P is
             # ~10^-837); this is the committed source for the "P < 10^-300" bound cited in [27].
             log10_P = round(phyper(sh - 1, pp, N - pp, ss, lower.tail = FALSE, log.p = TRUE) / log(10), 1),
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
cat("universe N =", N, "\n"); print(tab, row.names = FALSE)

## ---- write the overlap stats to a committed table (every cited number has a producer) -----
dir.create("results/integration", showWarnings = FALSE, recursive = TRUE)
tab_out <- cbind(universe_N = N, tab)          # N = shared-gene-space universe size (cited in [27])
write.table(tab_out, "results/integration/fig2c_orthogroup_overlap_stats.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote results/integration/fig2c_orthogroup_overlap_stats.tsv\n")

## ---- long form for plotting (explicit x positions for dodge + annotation) ---
glab <- c(cold_needle = "Cold\nneedle", cold_root = "Cold\nroot",
          drought_needle = "Drought\nneedle", drought_root = "Drought\nroot",
          overall = "Overall")
tab$gi <- seq_len(nrow(tab))
cats <- c("spruce-specific", "shared", "pine-specific")
long <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) data.frame(
  gi = tab$gi[i],
  cat = factor(cats, levels = cats),
  x   = tab$gi[i] + c(-0.27, 0, 0.27),
  n   = c(tab$sp_only[i], tab$shared[i], tab$pi_only[i]))))
# place the shared-% label above the tallest bar in each group so it never collides
ann <- data.frame(x = tab$gi, y = pmax(tab$sp_only, tab$shared, tab$pi_only),
                  lab = paste0(tab$pct_smaller, "%"))

p <- ggplot(long, aes(x = x, y = n, fill = cat)) +
  geom_col(width = 0.25) +
  geom_vline(xintercept = 4.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_text(data = ann, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            vjust = -0.4, size = 3, colour = PAL$shared_dark, fontface = "bold") +
  scale_fill_manual(values = PAL_FIG2C, name = NULL) +
  scale_x_continuous(breaks = tab$gi, labels = glab[tab$group], expand = expansion(add = 0.5)) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Orthogroups (DEG-hit, shared universe)",
       title = "Cross-species orthogroup-level DEG overlap") +
  theme_paper(base_size = 11, major_y = TRUE) +
  theme(legend.position = "top")

ggsave(file.path(OUT, "Figure2c.pdf"), p, width = 16/2.54, height = 10/2.54)
ggsave(file.path(OUT, "Figure2c.png"), p, width = 16/2.54, height = 10/2.54, dpi = 300)
cat("\nWrote", file.path(OUT, "Figure2c.pdf"), "and .png\n")
