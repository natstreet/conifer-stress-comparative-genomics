#!/usr/bin/env Rscript
# assemble_figures.R — build the FINAL manuscript figures (Figure1-7, FigureS1-S4)
# for van Zalen et al., "Comparative genomics of abiotic stress response in
# Norway spruce and Scots pine".
#
# Run from the AbioticStressConifers/ directory:  Rscript src/figures/assemble_figures.R
#
# Numbering matches the submitted manuscript:
#   Figure 1  PCA (spruce, pine) + DEG bar plots
#   Figure 2  DEG-orthogroup overlap (UpSet)
#   Figure 3  Co-expressolog conservation UpSet + TF families + directional TF   [direct script]
#   Figure 4  Conserved stress co-expressolog heatmaps (drought, cold)
#   Figure 5  Evolutionary/genomic features: YN00 dN/dS + SD enrichment [direct script]
#   Figure 6  CHS3 segmental-duplicate case study
#   Figure S1 Scots pine drought physiology          [figureS1_physiology.R, run separately]
#   Figure S2 per-timepoint cross-species DEG timing [figureS2_timing_deg_fraction.R, run separately]
#   Figure S3 not_coex mechanism classification + PAV enrichment
#   Figure S4 CHS gene-family ML phylogeny
#
# Figures 3 and 5 are produced directly by their own scripts (sourced below).
# Component panel PDFs come from the plot_*.R scripts in src/ComPlEx/ and src/.

suppressPackageStartupMessages({ library(grid); library(gridExtra) })
for (p in c("magick","pdftools","png"))
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
library(magick); library(png)

FIGS  <- "results"
IFIGS <- "results/integration/figures"
OUT   <- "results/final_figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

img_g <- function(f, dpi=150) {
  if (!file.exists(f)) { warning("Missing component: ", f); return(NULL) }
  tmp <- tempfile(fileext=".png")
  pdftools::pdf_convert(f, format="png", pages=1, filenames=tmp, dpi=dpi, verbose=FALSE)
  rasterGrob(png::readPNG(tmp), interpolate=TRUE)
}
label_g <- function(lbl) textGrob(lbl, x=0.03, y=0.97, just=c("left","top"),
                                  gp=gpar(fontsize=12, fontface="bold"))
panel <- function(f, lbl, dpi=150) {
  g <- img_g(f, dpi); if (is.null(g)) return(NULL)
  if (nchar(lbl)>0) arrangeGrob(label_g(lbl), g, nrow=2, heights=unit(c(0.04,0.96),"npc")) else g
}
save_fig <- function(grobs, out, w, h, ncol=1, wt=NULL, ht=NULL) {
  gl <- Filter(Negate(is.null), grobs)
  if (!length(gl)) { message("No panels — skipping ", basename(out)); return(invisible()) }
  pdf(out, width=w/2.54, height=h/2.54)
  args <- list(grobs=gl, ncol=ncol)
  if (!is.null(wt)) args$widths <- wt
  if (!is.null(ht)) args$heights <- ht
  do.call(grid.arrange, args); dev.off()
  png_out <- sub("\\.pdf$", ".png", out)
  tryCatch({ image_write(image_read_pdf(out, density=200), png_out, format="png") },
           error=function(e) message("  (png copy skipped for ", basename(out), ")"))
  message("Saved: ", basename(out))
}
run_script <- function(path, makes) {
  message("\n>> sourcing ", path, " (produces ", makes, ")")
  tryCatch(sys.source(path, envir=new.env(parent=globalenv())),
           error=function(e) message("  !! ", path, " failed: ", conditionMessage(e)))
}

# ── Figures produced directly by dedicated scripts ───────────────────────────
run_script("src/figures/figure2c_orthogroup_overlap.R", "Figure2c.pdf (Figure 2 panel c)")
run_script("src/figures/figure3_coexpressolog_tf.R",       "Figure3.pdf")
run_script("src/figures/figure5_dnds_evolution.R", "Figure5.pdf")
# The per-timepoint timing supplementary (src/figures/figureS2_timing_deg_fraction.R)
# runs DESeq2 on eight datasets, so — like Figure S1 — it is its own reproduce_paper.sh
# stage (figureS2_timing) rather than assembled here.

# ── Assembled multi-panel figures ────────────────────────────────────────────
# Figure 1 — four PCA panels a-d (fig_pca_panels, tagged a-d inside) above the
# stacked DEG-bar panel e (fig_deg_barplots).
save_fig(list(
  panel(file.path(FIGS,"fig_pca_panels.pdf"),"",dpi=150),      # a-d (tagged inside)
  panel(file.path(FIGS,"fig_deg_barplots.pdf"),"e",dpi=150)),
  file.path(OUT,"Figure1.pdf"), w=20, h=24, ncol=1, ht=c(0.56,0.44))

# Figure 2 — three panels: (a) Picea UpSet + (b) Pinus UpSet (both inside the self-
# labelled fig_deg_overlap.pdf) stacked above (c) the cross-species orthogroup-level
# DEG-overlap bars (Figure2c.pdf).
save_fig(list(
  panel(file.path(FIGS,"fig_deg_overlap.pdf"),"",dpi=150),   # a/b (self-labelled)
  panel(file.path(OUT,"Figure2c.pdf"),"c",dpi=150)),
  file.path(OUT,"Figure2.pdf"), w=18, h=26, ncol=1, ht=c(0.60,0.40))

# Figure 4 is the 4-panel conserved-dynamics figure (a overlap bar, b profile-conservation density,
# c drought heatmap, d cold heatmap), built and self-labelled by figure4_conserved_dynamics.R.
run_script("src/figures/figure4_conserved_dynamics.R", "Figure4_conserved_dynamics.pdf/.png")
for (ext in c("pdf","png")) {
  src4 <- file.path(IFIGS, paste0("Figure4_conserved_dynamics.", ext))
  if (file.exists(src4)) file.copy(src4, file.path(OUT, paste0("Figure4.", ext)), overwrite=TRUE)
}
message("Saved: Figure4.pdf/.png (4-panel conserved dynamics)")

# Figure 6 — CHS3 case study (regenerate component if its script is present)
if (file.exists("src/figures/figure6_chs3_case_study.R"))
  run_script("src/figures/figure6_chs3_case_study.R", "chs3_sd_pair_figure.pdf")
chs6 <- "results/integration/chs3/chs3_sd_pair_figure.pdf"
chs6 <- chs6[file.exists(chs6)][1]
if (!is.na(chs6)) save_fig(list(panel(chs6,"",dpi=120)), file.path(OUT,"Figure6.pdf"), w=18, h=14) else
  message("Figure6 component not found — keeping existing results/final_figures/Figure6.pdf")

# ── Supplementary figures ────────────────────────────────────────────────────
# Regenerate the Figure S3 component panels and the Figure S4 CHS tree from their
# generators (see src/ComPlEx/) so the supplementary figures are reproducible.
run_script("src/figures/figureS3a_notcoex_mechanism.R", "fig_notcoex_mechanism.pdf (Figure S3a)")
run_script("src/figures/figureS3b_pav_enrichment.R",    "fig_pav_categories.pdf (Figure S3b)")
run_script("src/figures/figureS4_chs_phylogeny.R",             "chs3_tree.pdf (Figure S4)")

save_fig(list(
  panel(file.path(IFIGS,"fig_notcoex_mechanism.pdf"),"a",dpi=150),
  panel(file.path(IFIGS,"fig_pav_categories.pdf"),   "b",dpi=150)),
  file.path(OUT,"FigureS3.pdf"), w=24, h=11, ncol=2)

save_fig(list(panel("results/integration/chs3/chs3_tree.pdf","",dpi=150)),
         file.path(OUT,"FigureS4.pdf"), w=18, h=20)

message("\nFigure S1 (Scots pine drought physiology) is built by src/figures/figureS1_physiology.R.")
message("\nAll final figures written to: ", OUT)
