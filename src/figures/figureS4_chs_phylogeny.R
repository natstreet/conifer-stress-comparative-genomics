#!/usr/bin/env Rscript
# figureS4_chs_phylogeny.R — Supplementary Figure S4
#
# Maximum-likelihood phylogeny of the Norway spruce chalcone synthase (CHS) gene
# family: the 20 Picea abies members of the CHS orthogroup (OG0000177). Protein
# sequences were aligned with MAFFT and a maximum-likelihood tree inferred with
# FastTree (LG model); SH-like local support values >= 0.8 are shown at nodes.
# The CHS3 segmental-duplicate pair (PA_chr01_G004115 / _G004116) is highlighted.
#
# Upstream (committed as inputs; regenerate with the commands below if desired):
#   mafft --auto chs3_family.faa > chs3_family.aln
#   FastTree -lg chs3_family.aln > chs3_family.nwk
# This script plots the committed tree so the figure is reproducible without the
# external aligners; if mafft and FastTree are on PATH and REBUILD=TRUE it will
# regenerate the alignment and tree from the sequences first.
#
# Input : results/integration/chs3/chs3_family.{faa,aln,nwk}
# Output: results/integration/chs3/chs3_tree.pdf
# Run from the AbioticStressConifers directory.

suppressPackageStartupMessages({ library(ape) })

# Run from the AbioticStressConifers/ project root (standalone or sourced by assemble_figures.R);
# getwd() is the project root in both cases. (commandArgs("--file=") is unreliable under source().)
PROJ_DIR    <- normalizePath(getwd())
CHS         <- file.path(PROJ_DIR, "results/integration/chs3")   # output dir (regenerated)
CHS_IN      <- file.path(PROJ_DIR, "data/chs3")                  # curated inputs (from the deposit)
FAA <- file.path(CHS_IN, "chs3_family.faa")
ALN <- file.path(CHS_IN, "chs3_family.aln")
NWK <- file.path(CHS_IN, "chs3_family.nwk")
CHS3_PAIR <- c("PA_chr01_G004115", "PA_chr01_G004116")

REBUILD <- toupper(Sys.getenv("REBUILD", "FALSE")) == "TRUE"
have <- function(cmd) nzchar(Sys.which(cmd))
if (REBUILD && have("mafft") && have("FastTree")) {
  cat("REBUILD: regenerating alignment (MAFFT) and tree (FastTree, LG)...\n")
  system(sprintf("mafft --auto %s > %s", shQuote(FAA), shQuote(ALN)))
  system(sprintf("FastTree -lg %s > %s", shQuote(ALN), shQuote(NWK)))
}
if (!file.exists(NWK)) stop("Tree not found: ", NWK, " (run MAFFT+FastTree, or set REBUILD=TRUE)")

tr <- read.tree(NWK)
cat(sprintf("Tree: %d tips\n", length(tr$tip.label)))

# FastTree writes SH-like supports into node labels; show those >= 0.8
supp <- suppressWarnings(as.numeric(tr$node.label))
supp_lab <- ifelse(!is.na(supp) & supp >= 0.8, formatC(supp, format = "f", digits = 2), "")

tip_col <- ifelse(tr$tip.label %in% CHS3_PAIR, "#D6604D", "black")
tip_font <- ifelse(tr$tip.label %in% CHS3_PAIR, 2, 1)

dir.create(CHS, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(CHS, "chs3_tree.pdf"), width = 18/2.54, height = 20/2.54)
par(mar = c(2, 1, 2, 1))
plot(tr, tip.color = tip_col, font = tip_font, cex = 0.8, no.margin = FALSE,
     main = "Norway spruce chalcone synthase gene family (OG0000177)")
nodelabels(supp_lab, frame = "none", adj = c(1.1, -0.3), cex = 0.6, col = "grey30")
add.scale.bar(cex = 0.7)
legend("bottomleft", legend = "CHS3 segmental-duplicate pair",
       text.col = "#D6604D", bty = "n", cex = 0.7)
dev.off()
cat("Saved chs3_tree.pdf\n")
