#!/usr/bin/env Rscript
## no1to1_orphan_vs_expanded.R
## Split the not_coex "no_1to1_ortholog" stress-DE set into true ORPHANS vs EXPANDED_FAMILY, per species,
## using the OrthoFinder orthogroups (a DIFFERENT method from the reciprocal-best 1:1 backbone that set the
## no_1to1 flag — so this is an orthogroup-level view, not a re-derivation of the 1:1 flag).
##   ORPHAN          = gene's orthogroup contains ZERO genes of the OTHER species, OR the gene is in no
##                     orthogroup at all (OrthoFinder-unassigned singleton).
##   EXPANDED_FAMILY = orthogroup contains >=1 gene of the other species (an orthologue exists) but the
##                     relationship is non-1:1 (these genes are, by construction, all outside the 1:1 backbone).
## ID GOTCHA: Orthogroups.tsv carries a ".mRNA.N" transcript suffix; the not_coex files use bare gene IDs.
suppressPackageStartupMessages(library(data.table))
I       <- Sys.getenv("INTEG_DIR", "results/integration")
OG_FILE <- Sys.getenv("ORTHOGROUPS_TSV", "Orthogroups/Orthogroups.tsv")
strip_tx <- function(x) sub("\\.mRNA\\.[0-9]+$", "", x)   # PA_..._G....mRNA.1 -> PA_..._G....

## ---- OrthoFinder membership: OG -> bare gene, per species ---------------------------------------
og <- fread(OG_FILE, sep = "\t", header = TRUE, quote = "")
melt_species <- function(colname) {
  d <- og[get(colname) != "", .(Orthogroup, cell = get(colname))]
  d[, .(gene = strip_tx(trimws(unlist(strsplit(cell, ",", fixed = TRUE))))), by = Orthogroup][gene != ""]
}
sp_long <- unique(melt_species("Picea_abies"),      by = c("Orthogroup", "gene"))  # OG, spruce gene
pi_long <- unique(melt_species("Pinus_sylvestris"), by = c("Orthogroup", "gene"))  # OG, pine gene
og_has_pine   <- unique(pi_long$Orthogroup)   # orthogroups containing >=1 Pinus gene
og_has_spruce <- unique(sp_long$Orthogroup)   # orthogroups containing >=1 Picea gene
g2og_sp <- unique(sp_long, by = "gene")       # spruce gene -> its OG
g2og_pi <- unique(pi_long, by = "gene")       # pine   gene -> its OG

## ---- classify each species' no_1to1 not_coex set ------------------------------------------------
classify <- function(genes, g2og, og_has_other) {
  d <- merge(data.table(gene = genes), g2og[, .(gene, Orthogroup)], by = "gene", all.x = TRUE)
  d[, class := fifelse(is.na(Orthogroup), "orphan",
                fifelse(Orthogroup %in% og_has_other, "expanded_family", "orphan"))]
  d[]
}
sp_nc <- fread(file.path(I, "not_coex_de_genes.tsv"))[mechanism == "no_1to1_ortholog", pa_gene]
pi_nc <- fread(file.path(I, "not_coex_de_genes_PINE.tsv"))[mechanism == "no_1to1_ortholog", ps_gene]
stopifnot(length(sp_nc) == 1859L, length(pi_nc) == 2676L)     # sanity gate

sp_cls <- classify(sp_nc, g2og_sp, og_has_pine)[,   `:=`(species = "Picea_abies")]
pi_cls <- classify(pi_nc, g2og_pi, og_has_spruce)[, `:=`(species = "Pinus_sylvestris")]

## ---- per-gene output (registrable) --------------------------------------------------------------
by_gene <- rbind(sp_cls[, .(species, gene, orthogroup = fcoalesce(Orthogroup, ""), class)],
                 pi_cls[, .(species, gene, orthogroup = fcoalesce(Orthogroup, ""), class)])
fwrite(by_gene, file.path(I, "no1to1_orphan_expanded_by_gene.tsv"), sep = "\t")

## ---- per-species summary (counts + pct) ---------------------------------------------------------
summ <- by_gene[, .(n = .N), by = .(species, class)]
summ[, total := sum(n), by = species]
summ[, pct := round(100 * n / total, 1)]
setorder(summ, species, -n)
fwrite(summ, file.path(I, "no1to1_orphan_expanded.tsv"), sep = "\t")

cat("no1to1 orphan-vs-expanded split (orthogroup method):\n"); print(summ)
cat(sprintf("\nmatched to an orthogroup: spruce %d/%d, pine %d/%d (unmatched -> orphan/no-OG)\n",
            sum(!is.na(sp_cls$Orthogroup)), nrow(sp_cls),
            sum(!is.na(pi_cls$Orthogroup)), nrow(pi_cls)))
