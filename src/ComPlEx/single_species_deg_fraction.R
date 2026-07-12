#!/usr/bin/env Rscript
## single_species_deg_fraction.R — fraction of each species' differentially expressed genes that fall in
## SPECIES-EXCLUSIVE (single-species) orthogroups, i.e. DEGs with no orthologue in the other species.
## Emits the in-text values of Results [89] ("1,027 Norway spruce and 1,461 Scots pine genes, 8.5% and
## 10.9% of the differentially expressed genes of each species"), which were previously computed inline.
##
## Total DEG genes per species = union of that species' four stress-tissue DEG_all sets (the same DEG lists
## the Figure 2a/b UpSet uses). Single-species DEG genes = those DEG genes whose orthogroup (in
## doc/genes_ortholog_categories.tsv, the same map fed to Figure 2c) is NOT shared between the two species
## (present in only one species), plus any DEG gene absent from the orthogroup map entirely.
##
## Output: results/integration/single_species_deg_fraction.tsv (species, n_single_species_genes,
##         n_total_deg_genes, pct)

OGX  <- "doc/genes_ortholog_categories.tsv"
DEGD <- "data/DEG_lists"
INTEG <- "results/integration"; dir.create(INTEG, showWarnings = FALSE, recursive = TRUE)

x <- read.table(OGX, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
sp_og  <- unique(x$Ortholog_Group[x$species == "Picea_abies"])
pi_og  <- unique(x$Ortholog_Group[x$species == "Pinus_sylvestris"])
shared <- intersect(sp_og, pi_og)                       # orthogroups present in BOTH species
sp_g2og <- setNames(x$Ortholog_Group[x$species == "Picea_abies"],     x$gene[x$species == "Picea_abies"])
pi_g2og <- setNames(x$Ortholog_Group[x$species == "Pinus_sylvestris"], x$gene[x$species == "Pinus_sylvestris"])

deg <- function(tag) {
  e <- new.env(); load(file.path(DEGD, sprintf("DE_all_%s_01_2L2FC.RData", tag)), e)
  as.character(get(ls(e)[1], e))
}
sp_deg <- unique(c(deg("SCN"), deg("SCR"), deg("SDN"), deg("SDR")))
pi_deg <- unique(c(deg("PCN"), deg("PCR"), deg("PDN"), deg("PDR")))

single_species <- function(genes, g2og) sum(vapply(genes, function(g) {
  og <- g2og[[g]]; is.null(og) || !(og %in% shared)      # no orthogroup, or a species-exclusive one
}, logical(1)))

out <- data.frame(
  species = c("Picea_abies", "Pinus_sylvestris"),
  n_single_species_genes = c(single_species(sp_deg, sp_g2og), single_species(pi_deg, pi_g2og)),
  n_total_deg_genes = c(length(sp_deg), length(pi_deg)))
out$pct <- round(100 * out$n_single_species_genes / out$n_total_deg_genes, 1)

write.table(out, file.path(INTEG, "single_species_deg_fraction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", file.path(INTEG, "single_species_deg_fraction.tsv"), "\n"); print(out, row.names = FALSE)
