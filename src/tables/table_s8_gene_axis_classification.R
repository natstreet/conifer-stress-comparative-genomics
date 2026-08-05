#!/usr/bin/env Rscript
## table_s8_gene_axis_classification.R
## Comprehensive per-gene conservation-divergence axis table (Supplementary Table S8).
## Base = every gene in the co-expression universe (coex_category_universe.tsv, genome-wide category),
## LEFT-JOINed with the paper's DE-set category and all cleanly per-gene-joinable committed features.
## Every value traces to a committed producer output; blank = not applicable.
## Population-genomic columns are NORWAY-SPRUCE-ONLY (no Scots pine resequencing) — see legend.
suppressPackageStartupMessages(library(data.table))
I       <- Sys.getenv("INTEG_DIR", "results/integration")
POPGEN  <- Sys.getenv("POPGEN_DIR", "data/popgen")
SIGFILE <- Sys.getenv("SD_POPGEN_SIGNALS", "sd_popgen_signals.tsv")
outfile <- Sys.getenv("S8_OUTFILE", "manuscript/supplementary/Supplementary_Table_S8_gene_axis_classification.xlsx")

## ---- BASE: UNION of the co-expression universe and the manuscript DE set -----------------------
## (some stress-DE genes are absent from the co-expression network / universe; they must still appear
##  so the paper's DE-set category counts reconcile in full — their genome-wide columns are blank.)
uni <- fread(file.path(I, "coex_category_universe.tsv"))
setnames(uni, "coex_category", "coex_category_genomewide")
de  <- fread(file.path(I, "de_gene_coex_category.tsv"))          # pa_gene, coex_category, n_de_tissues, stress_type
base <- merge(data.table(pa_gene = sort(unique(c(uni$pa_gene, de$pa_gene)))),
              uni, by = "pa_gene", all.x = TRUE)
base[, conservation_breadth := fifelse(is.na(cold_needle_present), NA_integer_,
       cold_needle_present + cold_root_present + drought_needle_present + drought_root_present)]

## ---- PAPER DE-SET category (the assignment behind the manuscript figures) ----------------------
base <- merge(base, de[, .(pa_gene, coex_category_paper = coex_category,
                           de_n_tissues = n_de_tissues, de_stress_type = stress_type)],
              by = "pa_gene", all.x = TRUE)
base[, stress_responsive_DE := fifelse(!is.na(coex_category_paper), "yes", "no")]

## ---- SD class + copy co-expression status (per gene from the pair table) ------------------------
sdf <- fread(file.path(I, "sd_pair_features.tsv"))
sd_long <- unique(rbindlist(list(
  sdf[, .(pa_gene = pa_gene1, sd_class, sd_pair_coex_class = pair_coex_class)],
  sdf[, .(pa_gene = pa_gene2, sd_class, sd_pair_coex_class = pair_coex_class)])), by = "pa_gene")
base <- merge(base, sd_long, by = "pa_gene", all.x = TRUE)

## ---- dN/dS (1:1 spruce–pine orthologue; paper QC filter dS in (0,5), dN/dS<10; one row per gene) -
yn <- fread(file.path(I, "cross_species_dnds_yn00.tsv"))
yn <- unique(yn[dS > 0 & dS < 5 & dNdS < 10, .(pa_gene, dNdS_spruce_pine = round(dNdS, 4))], by = "pa_gene")
base <- merge(base, yn, by = "pa_gene", all.x = TRUE)

## ---- not_coex DE mechanism (only defined for not_coex stress-DE genes) --------------------------
nc <- fread(file.path(I, "not_coex_de_genes.tsv"))[, .(pa_gene, not_coex_mechanism = mechanism)]
base <- merge(base, nc, by = "pa_gene", all.x = TRUE)

## ---- POPGEN — NORWAY-SPRUCE-ONLY -------------------------------------------------------------
## selection scan (genome-wide iHS / XP-EHH); file is headerless: gene <tab> signal
sel <- fread(file.path(POPGEN, "selection_xpehh_ihs_genes.tsv"), header = FALSE, col.names = c("pa_gene", "sig"))
sel[, selection_scan_spruce := fifelse(grepl("xpehh", sig) & grepl("ihs", sig), "both",
                                fifelse(grepl("xpehh", sig), "xpehh",
                                fifelse(grepl("ihs", sig), "ihs", NA_character_)))]
base <- merge(base, sel[, .(pa_gene, selection_scan_spruce)], by = "pa_gene", all.x = TRUE)
base[is.na(selection_scan_spruce), selection_scan_spruce := "none"]
## PAV / GWAS per-gene membership (the file the manuscript's popgen category enrichment is built on)
sig <- fread(SIGFILE)                                            # gene, chrom, midpos, signal
base <- merge(base, sig[, .(pa_gene = gene,
                            pav_spruce  = fifelse(grepl("pav",  signal), "yes", ""),
                            gwas_spruce = fifelse(grepl("gwas", signal), "yes", ""))],
              by = "pa_gene", all.x = TRUE)
base[is.na(pav_spruce),  pav_spruce  := ""]
base[is.na(gwas_spruce), gwas_spruce := ""]

## ---- tidy: blanks for character NAs, order columns ---------------------------------------------
for (cc in c("coex_category_paper","de_stress_type","sd_class","sd_pair_coex_class","not_coex_mechanism"))
  base[is.na(get(cc)), (cc) := ""]
setcolorder(base, c("pa_gene","OrthoGroup","coex_category_genomewide",
  "cold_needle_present","cold_root_present","drought_needle_present","drought_root_present",
  "conservation_breadth","stress_responsive_DE","coex_category_paper","de_n_tissues","de_stress_type",
  "sd_class","sd_pair_coex_class","dNdS_spruce_pine","not_coex_mechanism",
  "selection_scan_spruce","pav_spruce","gwas_spruce"))
setorder(base, pa_gene)
cat(sprintf("S8 rows: %d (36632 co-expression universe + %d DE-only genes outside it)\n",
            nrow(base), nrow(base) - 36632L))

## ---- reconciliation counts (report both categorizations) ---------------------------------------
cat("\n-- genome-wide category (coex_category_genomewide) --\n")
print(base[, .N, by = coex_category_genomewide][order(coex_category_genomewide)])
cat("\n-- paper DE-set category (coex_category_paper; blank = non-DE) --\n")
print(base[coex_category_paper != "", .N, by = coex_category_paper][order(coex_category_paper)])

## ---- LEGEND -----------------------------------------------------------------------------------
legend <- data.table(Column = character(), Definition = character())
add <- function(c, d) legend <<- rbind(legend, data.table(Column = c, Definition = d))
add("pa_gene", "Norway spruce (Picea abies) gene identifier.")
add("OrthoGroup", "Hierarchical orthogroup (N10) containing the gene.")
add("coex_category_genomewide", "Co-expression conservation category assigned genome-wide (conserved / cold_specific / drought_specific / multi_tissue / not_coex). Source: coex_category_universe.tsv. Blank for genes not in the co-expression universe (stress-DE genes outside the network). NOTE: this genome-wide assignment differs from the DE-restricted category used in the manuscript figures — see coex_category_paper.")
add("cold_needle_present / cold_root_present / drought_needle_present / drought_root_present", "1 if the gene has a conserved cross-species co-expressolog in that stress-tissue comparison, else 0. Source: coex_category_universe.tsv.")
add("conservation_breadth", "Number of the four stress-tissue comparisons with a conserved co-expressolog (0-4); sum of the four *_present columns.")
add("stress_responsive_DE", "yes if the gene is differentially expressed under stress in Norway spruce (in the manuscript DE set), else no. Source: de_gene_coex_category.tsv.")
add("coex_category_paper", "The co-expression category as used in the manuscript figures/Results, defined over stress-responsive DE genes only (blank for non-DE genes). Per-category totals reconcile with the paper. Source: de_gene_coex_category.tsv (from not_coex_de_analysis.R).")
add("de_n_tissues / de_stress_type", "For DE genes: number of stress-tissue comparisons in which the gene is DE, and the stress type(s). Blank for non-DE genes. Source: de_gene_coex_category.tsv.")
add("sd_class", "Structural-duplicate class of the gene's SD pair: shared_SD (duplicated in both species) or spruce_only_SD; blank if the gene is not in an SD pair. Source: sd_pair_features.tsv.")
add("sd_pair_coex_class", "For SD-pair genes: whether both duplicate copies remained co-expressed in the network (both_coex / diverged / both_not_coex). Source: sd_pair_features.tsv.")
add("dNdS_spruce_pine", "Yang-Nielsen dN/dS for the 1:1 Norway spruce-Scots pine orthologue pair (QC filter: dS in (0,5), dN/dS<10). Blank if the gene has no usable 1:1 orthologue. Source: cross_species_dnds_yn00.tsv.")
add("not_coex_mechanism", "For not_coex stress-DE genes only: why the gene lacks a conserved co-expressolog (no_1to1_ortholog / pine_not_DE / diverged_regulation / pine_orth_not_expressed). Blank otherwise. Source: not_coex_de_genes.tsv.")
add("selection_scan_spruce", "NORWAY SPRUCE ONLY (no Scots pine resequencing). Positive-selection scan signal at the gene: ihs, xpehh, both, or none. Source: data/popgen/selection_xpehh_ihs_genes.tsv.")
add("pav_spruce", "NORWAY SPRUCE ONLY. yes if the gene carries a presence/absence-variation signal in the population-genomic set underlying the manuscript's popgen enrichment; blank otherwise. Source: sd_popgen_signals.tsv.")
add("gwas_spruce", "NORWAY SPRUCE ONLY. yes if the gene carries a climate-GWAS signal in that same set; blank otherwise. Source: sd_popgen_signals.tsv.")

caption <- paste0("Supplementary Table S8. Per-gene conservation-divergence axis classification for Norway spruce. ",
  "Each row is one gene placed on the co-expression conservation axis (genome-wide and, for stress-responsive genes, ",
  "the manuscript DE-set category), with per-gene evolutionary, structural-duplicate and population-genomic features. ",
  "Population-genomic columns (selection_scan_spruce, pav_spruce, gwas_spruce) are Norway-spruce-only; there is no ",
  "Scots pine resequencing. Blank cells indicate the feature is not applicable to that gene. See the Legend sheet ",
  "for full column definitions and committed data sources.")

dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
## writexl produces a clean, minimal, universally-openable xlsx (no dangling drawing/printerSettings
## rels, unlike this openxlsx build). Caption + Legend are their own sheets.
library(writexl)
write_xlsx(list(`Supplementary Table S8` = base,
                Legend = legend,
                Caption = data.frame(Caption = caption, check.names = FALSE)),
           outfile)
cat("Wrote (writexl):", outfile, "\n")
## also emit the flat TSV (traceable, diffable)
fwrite(base, file.path(I, "gene_axis_classification_S8.tsv"), sep = "\t")
cat("Wrote:", file.path(I, "gene_axis_classification_S8.tsv"), "\n")
