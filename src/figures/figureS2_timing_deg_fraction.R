#!/usr/bin/env Rscript
## figureS2_timing_deg_fraction.R — Supplementary Figure S2 + Supplementary Table S3: timing
## per-timepoint cross-species DEG orthogroup overlap for spruce vs pine.
##
## For each matched sampling point (cold_needle/cold_root/drought_needle/drought_root)
## and each timepoint, differential expression is computed per contrast from the
## DESeq2 dataset object (DESeq(dds); Wald test per coefficient; independent-filter
## statistic = per-gene median of normalised counts, matching gwas_deg_overlap.R),
## DE = padj < 0.01 & |log2FC| >= 2. A gene's orthogroups are taken within the
## shared universe (orthogroups with >=1 spruce AND >=1 pine gene, N = 11,943);
## an orthogroup is "shared" at a timepoint if it is DE-hit in BOTH species.
##
## FIGURE (FigureS2.{pdf,png}): 4 panels (one per sampling point). Per
## timepoint, one stacked bar per species = fraction of that species' DEGs (gene
## level) lying in shared orthogroups vs species-specific.
## TABLE (Supplementary_Table_S3_timing_overlap.{csv,xlsx}): orthogroup-level counts,
## Jaccard, fold-enrichment, hypergeometric P, and a low-power flag.
##
## Run from AbioticStressConifers/:  Rscript src/figures/figureS2_timing_deg_fraction.R
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats); library(ggplot2) })
source("src/lib/fig_palette.R")

OGX <- "doc/genes_ortholog_categories.tsv"
DDS <- "data/dds"
OUTF <- "results/final_figures"; dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)
OUTT <- "results/integration";   dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)

## ---- universe + gene->orthogroup maps ---------------------------------------
x <- read.table(OGX, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
sp_map <- x[x$species == "Picea_abies",     c("gene", "Ortholog_Group")]
pi_map <- x[x$species == "Pinus_sylvestris", c("gene", "Ortholog_Group")]
universe <- intersect(unique(sp_map$Ortholog_Group), unique(pi_map$Ortholog_Group))
N <- length(universe)
g2og_sp <- split(sp_map$Ortholog_Group, sp_map$gene)
g2og_pi <- split(pi_map$Ortholog_Group, pi_map$gene)

## ---- per-dataset per-timepoint DE genes -------------------------------------
## returns list(labels = non-reference level names in order,
##              de = list of DE gene-ID vectors, one per timepoint)
de_by_timepoint <- function(tag) {
  e <- new.env(); load(file.path(DDS, sprintf("dds_%s.rda", tag)), e)
  dds <- get(ls(e)[1], e)
  dds <- DESeq(dds, quiet = TRUE)
  fac <- all.vars(design(dds)); fac <- fac[length(fac)]
  labels <- levels(colData(dds)[[fac]])[-1]                       # non-reference, in order
  coefs  <- setdiff(resultsNames(dds), "Intercept")              # same order as labels
  filt   <- matrixStats::rowMedians(counts(dds, normalized = TRUE))
  de <- lapply(coefs, function(cf) {
    res <- results(dds, name = cf, filter = filt)
    rownames(res)[which(!is.na(res$padj) & res$padj < 0.01 & abs(res$log2FoldChange) >= 2)]
  })
  list(labels = labels, coefs = coefs, de = de)
}

ogs_of <- function(genes, map)
  intersect(unique(unlist(map[genes[genes %in% names(map)]], use.names = FALSE)), universe)
genes_in_shared <- function(genes, map, shared)                  # DEG genes whose OG is shared
  sum(vapply(genes[genes %in% names(map)],
             function(g) any(map[[g]] %in% shared), logical(1)))

samplings <- list(
  cold_needle    = c(s = "SCN", p = "PCN"),
  cold_root      = c(s = "SCR", p = "PCR"),
  drought_needle = c(s = "SDN", p = "PDN"),
  drought_root   = c(s = "SDR", p = "PDR"))

tab <- list(); frac <- list()
for (sp in names(samplings)) {
  message(">> ", sp)
  S <- de_by_timepoint(samplings[[sp]]["s"])
  P <- de_by_timepoint(samplings[[sp]]["p"])
  n <- min(length(S$de), length(P$de))
  for (i in seq_len(n)) {
    sg <- S$de[[i]]; pg <- P$de[[i]]
    so <- ogs_of(sg, g2og_sp); po <- ogs_of(pg, g2og_pi)
    shared <- intersect(so, po); sh <- length(shared)
    ss <- length(so); pp <- length(po)
    tab[[length(tab) + 1]] <- data.frame(
      sampling_point = sp, timepoint = S$labels[i],
      spruce_DEG_OGs = ss, pine_DEG_OGs = pp, shared = sh,
      jaccard = ifelse(ss + pp - sh > 0, round(sh / (ss + pp - sh), 3), NA),
      fold_enrichment = ifelse(ss > 0 & pp > 0, round(sh / (ss * pp / N), 2), NA),
      P = if (ss > 0 & pp > 0) phyper(sh - 1, pp, N - pp, ss, lower.tail = FALSE) else NA,
      low_power = min(ss, pp) < 15, stringsAsFactors = FALSE)
    # gene-level fractions (denominator = that species' TOTAL DEGs at the timepoint;
    # numerator = DEGs whose orthogroup is shared, which is inherently within-universe)
    s_in <- genes_in_shared(sg, g2og_sp, shared); p_in <- genes_in_shared(pg, g2og_pi, shared)
    frac[[length(frac) + 1]] <- data.frame(
      sampling_point = sp, timepoint = S$labels[i], tp_order = i,
      species = c("spruce", "pine"),
      total = c(length(sg), length(pg)), in_shared = c(s_in, p_in),
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, tab); frac <- do.call(rbind, frac)

## ---- table out --------------------------------------------------------------
tab$P <- signif(tab$P, 3)
write.csv(tab, file.path(OUTT, "Supplementary_Table_S3_timing_overlap.csv"), row.names = FALSE)
dir.create("manuscript/supplementary", showWarnings = FALSE, recursive = TRUE)
# writexl writes a clean minimal xlsx (this openxlsx build emits dangling drawing/printerSettings
# rels that openpyxl / strict Excel cannot open).
writexl::write_xlsx(as.data.frame(tab), "manuscript/supplementary/Supplementary_Table_S3_timing_overlap.xlsx")
cat("universe N =", N, "\n\nSPOT-CHECK rows:\n")
sc <- subset(tab, (sampling_point == "cold_root" & timepoint %in% c("5C_6h", "neg5C_10d")) |
                   (sampling_point == "drought_root" & timepoint == "Collapse"))
print(sc, row.names = FALSE)
cat("\ngene-fraction spot-checks:\n")
gf <- subset(frac, (sampling_point == "cold_needle" & timepoint == "5C_6h") |
                    (sampling_point == "drought_root" & timepoint == "Collapse"))
gf$pct <- round(100 * gf$in_shared / gf$total, 1)
print(gf[, c("sampling_point","timepoint","species","in_shared","total","pct")], row.names = FALSE)

## ---- figure: 4 panels, stacked bars by gene fraction ------------------------
frac$specific <- frac$total - frac$in_shared
long <- rbind(
  data.frame(frac[c("sampling_point","timepoint","tp_order","species")],
             part = ifelse(frac$species == "spruce", "spruce shared", "pine shared"),
             n = frac$in_shared),
  data.frame(frac[c("sampling_point","timepoint","tp_order","species")],
             part = "species-specific", n = frac$specific))
long$part <- factor(long$part, levels = c("spruce shared","pine shared","species-specific"))
long$xlab <- factor(reorder(long$timepoint, long$tp_order))
sp_labels <- c(cold_needle="Cold needle", cold_root="Cold root",
               drought_needle="Drought needle", drought_root="Drought root")
long$sampling_point <- factor(long$sampling_point, levels = names(sp_labels), labels = sp_labels)

p <- ggplot(long, aes(x = interaction(species, xlab), y = n, fill = part)) +
  geom_col(width = 0.85) +
  facet_wrap(~ sampling_point, scales = "free_x", nrow = 2) +
  scale_fill_manual(values = PAL_TIMING, name = NULL) +
  labs(x = "Species x timepoint", y = "DEGs (gene count)",
       title = "Per-timepoint DEG fraction in shared orthogroups") +
  theme_paper(base_size = 10, major_y = TRUE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
        legend.position = "top")

ggsave(file.path(OUTF, "FigureS2.pdf"), p, width = 22/2.54, height = 16/2.54)
ggsave(file.path(OUTF, "FigureS2.png"), p, width = 22/2.54, height = 16/2.54, dpi = 300)
cat("\nWrote FigureS2.{pdf,png} and Supplementary_Table_S3_timing_overlap.{csv,xlsx}\n")
