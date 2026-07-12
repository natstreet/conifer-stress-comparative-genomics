#!/usr/bin/env Rscript
# integration_analysis.R
#
# Integrates co-expressolog tissue-coverage classifications with four datasets:
#
#   Item 1: Spruce SD genes × co-expressolog status
#            Universe = all spruce genes in ComPlEx expression matrix.
#            Tests whether SD-duplicated genes are over/under-represented
#            as cross-species co-expressologs.
#
#   Item 2: SD gene pairs × TE-in-promoter
#            For SD pairs where the two spruce copies differ in co-expressolog
#            category with their pine ortholog, tests whether the diverged copy
#            is more likely to carry a TE-derived promoter element.
#
#   Item 3: Cross-species Ks / KaKs × co-expressolog conservation
#            Backbone = 1:1 PA-PS HOG orthologs (interspecies_aa_identity.tsv).
#            Tests whether co-expression conservation correlates with coding
#            sequence constraint. Includes GO enrichment.
#
#   Item 4: Population-genetics signals (PAV, GWAS, selection) × co-expressolog
#            Tests enrichment of PAV/GWAS/selection genes in co-expressolog
#            categories.
#
# All outputs written to results/integration/ as TSV files.
# Source file provenance is logged in integration_sources.tsv.

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()

PROJ_DIR      <- normalizePath(file.path(SCRIPT_DIR, "../.."))
SD_DIR        <- PROJ_DIR
COMPLEX_RDATA <- file.path(PROJ_DIR, "results/ComPlEx/RData")
EXPR_DIR      <- file.path(PROJ_DIR, "data/expression")
OUTDIR        <- file.path(PROJ_DIR, "results/integration")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

sources <- list()
ROOT_FOR_LOG <- normalizePath(file.path(PROJ_DIR, ".."))   # repo root, for portable provenance paths
log_source <- function(label, path) {
  if (!file.exists(path)) stop(sprintf("Source file missing: %s", path))
  rel <- sub(paste0("^", ROOT_FOR_LOG, "/?"), "", normalizePath(path))
  sources[[label]] <<- list(path = rel, mtime = format(file.info(path)$mtime))
  invisible(path)
}

# ── Load co-expressolog gene-pair data ─────────────────────────────────────────
cat("Loading co-expressolog data...\n")
log_source("weighted_gene_pairs", file.path(COMPLEX_RDATA, "weighted_gene_pairs.tsv"))
wp <- fread(file.path(COMPLEX_RDATA, "weighted_gene_pairs.tsv"))
setnames(wp, c("Species1","Species2"), c("pa_gene","ps_gene"))

# Best pair per (pa_gene, ps_gene)
wp_best <- wp[order(best_pval)][!duplicated(paste(pa_gene, ps_gene))]
cat(sprintf("  %d co-expressolog gene pairs; %d unique (pa,ps) pairs\n",
            nrow(wp), nrow(wp_best)))

# ── Helper: Fisher's exact test ────────────────────────────────────────────────
fisher_enrich <- function(n_grp_pos, n_grp, n_all_pos, n_all, label = "") {
  n_grp_neg <- n_grp - n_grp_pos
  n_bg_pos  <- n_all_pos - n_grp_pos
  n_bg_neg  <- n_all - n_all_pos - n_grp_neg
  mat <- matrix(c(n_grp_pos, n_grp_neg, n_bg_pos, n_bg_neg), nrow = 2)
  ft  <- fisher.test(mat)
  data.table(
    label         = label,
    n_group       = n_grp,
    n_positive    = n_grp_pos,
    rate_positive = round(n_grp_pos / n_grp, 4),
    n_bg          = n_all - n_grp,
    n_bg_positive = n_bg_pos,
    rate_bg       = round(n_bg_pos / (n_all - n_grp), 4),
    odds_ratio    = round(ft$estimate, 3),
    pvalue        = ft$p.value,
    ci_lo         = round(ft$conf.int[1], 3),
    ci_hi         = round(ft$conf.int[2], 3)
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# Item 1: Spruce SD genes × cross-species co-expressolog status
#
# Universe: all spruce genes in the ComPlEx expression matrix (SD + drought
# expression files). These were the genes eligible to be tested.
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Item 1: SD genes × co-expressolog ──────────────────────────────────\n")

# Build universe of all spruce genes tested in ComPlEx
# Use both stress expression files (cold + drought share the same gene set)
log_source("SC_expression", file.path(EXPR_DIR, "SC_expression.txt"))
log_source("SD_expression", file.path(EXPR_DIR, "SD_expression.txt"))
sc_genes <- fread(file.path(EXPR_DIR, "SC_expression.txt"), select = 1)[[1]]
sd_genes_expr <- fread(file.path(EXPR_DIR, "SD_expression.txt"), select = 1)[[1]]
spruce_universe <- unique(c(sc_genes, sd_genes_expr))
cat(sprintf("  Spruce expression universe: %d genes\n", length(spruce_universe)))

# Spruce genes with at least one cross-species co-expressolog (in any tissue)
coex_spruce <- unique(wp_best$pa_gene)
cat(sprintf("  Spruce genes with >=1 co-expressolog: %d\n", length(coex_spruce)))

# SD genes from within-spruce Ka/Ks analysis
log_source("kaks_results", file.path(SD_DIR, "kaks_results.tsv"))
kaks <- fread(file.path(SD_DIR, "kaks_results.tsv"))
sd_gene_tbl <- unique(rbind(
  kaks[, .(pa_gene = gene1, sd_class = category)],
  kaks[, .(pa_gene = gene2, sd_class = category)]
))
# If gene appears in both classes, keep shared_SD
sd_gene_tbl <- sd_gene_tbl[order(sd_class)][!duplicated(pa_gene)]

cat(sprintf("  SD spruce genes: %d  (shared_SD: %d  spruce_only_SD: %d)\n",
            nrow(sd_gene_tbl),
            sum(sd_gene_tbl$sd_class == "shared_SD"),
            sum(sd_gene_tbl$sd_class == "spruce_only_SD")))

# How many SD genes are in the expression universe?
sd_in_universe <- sd_gene_tbl[pa_gene %in% spruce_universe]
cat(sprintf("  SD genes in expression universe: %d\n", nrow(sd_in_universe)))

# Build gene-level table: universe × (SD class, has co-expressolog)
universe_dt <- data.table(pa_gene = spruce_universe)
universe_dt <- merge(universe_dt, sd_gene_tbl, by = "pa_gene", all.x = TRUE)
universe_dt[is.na(sd_class), sd_class := "non_SD"]
universe_dt[, has_coex := pa_gene %in% coex_spruce]

n_total <- nrow(universe_dt)
n_coex  <- sum(universe_dt$has_coex)

cat(sprintf("  Overall co-expressolog rate in universe: %d/%d = %.1f%%\n",
            n_coex, n_total, 100 * n_coex / n_total))

# Rate by SD class
sd_coex_rate <- universe_dt[, .(
  n_total   = .N,
  n_coex    = sum(has_coex),
  rate_coex = round(mean(has_coex), 4)
), by = sd_class][order(sd_class)]
cat("  Co-expressolog rate by SD class:\n")
print(sd_coex_rate)
fwrite(sd_coex_rate, file.path(OUTDIR, "sd_coexpression_rate.tsv"), sep = "\t")

# Fisher enrichment test
sd_enrichment_dt <- rbindlist(lapply(c("shared_SD","spruce_only_SD"), function(cls) {
  n_cls  <- sum(universe_dt$sd_class == cls)
  if (n_cls == 0) return(NULL)
  n_pos  <- sum(universe_dt$sd_class == cls & universe_dt$has_coex)
  fisher_enrich(n_pos, n_cls, n_coex, n_total, label = cls)
}))
if (nrow(sd_enrichment_dt) > 0) {
  cat("  Fisher enrichment (SD co-expressolog rate vs non-SD background):\n")
  print(sd_enrichment_dt[, .(label, n_group, rate_positive, rate_bg, odds_ratio, pvalue)])
  fwrite(sd_enrichment_dt, file.path(OUTDIR, "sd_category_enrichment.tsv"), sep = "\t")
}

# Also look at co-expressolog category for SD genes with co-expressologs
# (best category per spruce gene = based on best-pval pair)
wp_gene_cat <- wp_best[order(best_pval)][!duplicated(pa_gene),
  .(pa_gene, coex_category = ifelse(
    conserved & NegLog10CliqueSum >= 10,  "conserved",  # Methods threshold (cliques_step2 C_SUM)
    ifelse(cold_specific,    "cold_specific",
    ifelse(drought_specific, "drought_specific",
                             "multi_tissue"))))]

universe_dt2 <- merge(universe_dt, wp_gene_cat, by = "pa_gene", all.x = TRUE)
universe_dt2[is.na(coex_category), coex_category := "not_coex"]

sd_category_table <- universe_dt2[, .N, by = .(sd_class, coex_category)][order(sd_class, coex_category)]
cat("  Co-expressolog category distribution by SD class:\n")
print(sd_category_table)
fwrite(sd_category_table, file.path(OUTDIR, "sd_category_counts.tsv"), sep = "\t")


# ═══════════════════════════════════════════════════════════════════════════════
# Item 2: SD gene pairs x co-expressolog status
#
# Classify each segmental-duplicate spruce gene pair by the cross-species
# co-expressolog status of its two copies (both_coex / both_not_coex / diverged).
# This per-pair classification feeds the Fig 7 integrative model (sd_class) and
# te_promoter_family_analysis Analysis B (pair_coex_class), via sd_pair_features.tsv.
# ===============================================================================
cat("\n-- Item 2: SD pairs x co-expressolog status --------------------------\n")

# Assign co-expressolog category per spruce gene
gene_cat <- universe_dt2[, .(pa_gene, coex_category)]

# Build SD pair table: one row per SD pair with co-expressolog class per copy
sd_pairs <- kaks[, .(pa_gene1 = gene1, pa_gene2 = gene2, sd_class = category)]
sd_pairs <- merge(sd_pairs, gene_cat, by.x = "pa_gene1", by.y = "pa_gene", all.x = TRUE)
setnames(sd_pairs, "coex_category", "cat_gene1")
sd_pairs <- merge(sd_pairs, gene_cat, by.x = "pa_gene2", by.y = "pa_gene", all.x = TRUE)
setnames(sd_pairs, "coex_category", "cat_gene2")
sd_pairs[is.na(cat_gene1), cat_gene1 := "not_coex"]
sd_pairs[is.na(cat_gene2), cat_gene2 := "not_coex"]

# Classify each pair: both_coex / both_not_coex / diverged (one coex, one not)
sd_pairs[, pair_coex_class := ifelse(
  cat_gene1 != "not_coex" & cat_gene2 != "not_coex", "both_coex",
  ifelse(cat_gene1 == "not_coex" & cat_gene2 == "not_coex", "both_not_coex",
  "diverged"))]

cat("  SD pair co-expressolog status distribution:\n")
print(sd_pairs[, .N, by = .(sd_class, pair_coex_class)][order(sd_class, pair_coex_class)])
fwrite(sd_pairs[, .N, by = .(sd_class, pair_coex_class)][order(sd_class, pair_coex_class)],
       file.path(OUTDIR, "sd_pair_coexpression_class_counts.tsv"), sep = "\t")

# -- Join stress expression divergence onto SD pairs (for sd_pair_features) -----
log_source("sd_pair_expression_divergence", file.path(SD_DIR, "sd_pair_expression_divergence.tsv"))
te_stress <- fread(file.path(SD_DIR, "sd_pair_expression_divergence.tsv"))
te_stress_pa <- te_stress[species == "PA"]   # one row per (hog_id, species); keep PA

sd_pairs_te <- merge(sd_pairs,
  te_stress_pa[, .(pa_gene1 = gene1, pa_gene2 = gene2,
                   div_cold_needles    = div_cold_needles,
                   div_cold_roots      = div_cold_roots,
                   div_drought_needles = div_drought_needles,
                   div_drought_roots   = div_drought_roots)],
  by = c("pa_gene1","pa_gene2"), all.x = TRUE)

# Per-pair table: sd_class + pair_coex_class feed the Fig 7 integrative model and
# te_promoter_family_analysis Analysis B (cross-species promoter TE Jaccard).
fwrite(sd_pairs_te, file.path(OUTDIR, "sd_pair_features.tsv"), sep = "	")


# ═══════════════════════════════════════════════════════════════════════════════
# Item 3: Cross-species Ks / KaKs × co-expressolog conservation
#
# Backbone: all cross-species pairs in cross_species_ks.tsv (pair_type ==
# "background_1to1"). Despite the label this file includes pairs from
# multi-copy orthogroups — using it directly recovers ~1,400 co-expressolog
# pairs that the old 1:1 HOG backbone excluded. One row per unique
# (pa_gene, ps_gene) pair is retained.
#
# pct_identity is joined from interspecies_aa_identity.tsv where available
# (strictly 1:1 HOG pairs only; NA for multi-copy pairs).
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Item 3: Ks / KaKs × co-expressolog ─────────────────────────────────\n")

log_source("cross_species_ks", file.path(SD_DIR, "cross_species_ks.tsv"))
ks_raw <- fread(file.path(SD_DIR, "cross_species_ks.tsv"))
setnames(ks_raw, c("spruce_gene","pine_gene"), c("pa_gene","ps_gene"))

# Keep cross-species background pairs only; deduplicate to one row per gene pair
backbone <- ks_raw[pair_type == "background_1to1"][
  !duplicated(paste(pa_gene, ps_gene)),
  .(hog_id, pa_gene, ps_gene)]   # cross_species_ks is now a pair list only (no NG86 Ka/Ks)
cat(sprintf("  Backbone: %d cross-species gene pairs (from cross_species_ks.tsv)\n",
            nrow(backbone)))

# Join pct_identity from 1:1 HOG backbone where available
log_source("interspecies_aa_identity",
           file.path(SD_DIR, "interspecies_aa_identity.tsv"))
pct_id <- fread(file.path(SD_DIR, "interspecies_aa_identity.tsv"),
                select = c("pa_gene","ps_gene","pct_identity"))
backbone <- merge(backbone, pct_id, by = c("pa_gene","ps_gene"), all.x = TRUE)
cat(sprintf("  Pairs with pct_identity (strictly 1:1): %d\n",
            sum(!is.na(backbone$pct_identity))))

# Join co-expressolog category
backbone <- merge(backbone, wp_best, by = c("pa_gene","ps_gene"), all.x = TRUE)
backbone[, coex_category := ifelse(
  is.na(n_tissues), "not_coex",
  ifelse(conserved & NegLog10CliqueSum >= 10,  "conserved",  # Methods threshold (cliques_step2 C_SUM)
  ifelse(cold_specific,    "cold_specific",
  ifelse(drought_specific, "drought_specific",
                           "multi_tissue"))))]

cat(sprintf("  Category distribution:\n"))
print(backbone[, .N, by = coex_category][order(-N)])

# Cross-species dN/dS is PAML YN00 (cross_species_dnds_yn00.tsv). The Nei-Gojobori Ka/Ks columns are
# not carried here, to avoid shipping a second, method-inconsistent divergence estimate.
fwrite(backbone[, .(hog_id, pa_gene, ps_gene, pct_identity, coex_category,
                    n_tissues, best_pval, tissue_pattern)],
       file.path(OUTDIR, "integration_backbone_1to1.tsv"), sep = "\t")

# Cross-species Ka/Ks is reported as PAML YN00 dN/dS (cross_species_dnds_yn00.tsv; see table_s6 and
# figure5). The backbone written above (pair list + coex_category) is the load-bearing output here.
coex_cats <- setdiff(unique(backbone$coex_category), "not_coex")   # category list reused by the TE Fisher tests below

# Write gene lists per category for GO enrichment (external tool)
for (cat in unique(backbone$coex_category)) {
  genes <- backbone[coex_category == cat, pa_gene]
  fwrite(data.table(pa_gene = genes),
         file.path(OUTDIR, sprintf("category_genes_%s.tsv", cat)), sep = "\t")
}

# Per-species category counts for manuscript reporting (spruce pa_gene, pine ps_gene).
# Reported per species because orthology is not strictly one-to-one; counts differ
# slightly between the two genomes.
cat_counts_species <- backbone[, .(spruce_genes = uniqueN(pa_gene),
                                    pine_genes   = uniqueN(ps_gene),
                                    pairs        = .N), by = coex_category]
fwrite(cat_counts_species,
       file.path(OUTDIR, "category_counts_by_species.tsv"), sep = "\t")
cat(sprintf("  Gene lists saved for GO enrichment (one file per category)\n"))

# Item 3 extension: co-expression non-conserved genes — TE-promoter and KaKs links
# Compare 'not_coex' to other categories for TE-in-promoter
log_source("pa_ps_promoter_te_results",
           file.path(SD_DIR, "pa_ps_promoter_te_results.tsv"))
prom_te_raw2 <- fread(file.path(SD_DIR, "pa_ps_promoter_te_results.tsv"))
prom_pa2 <- prom_te_raw2[species == "PA",
  .(pa_gene = gene_id, pa_prom_te = as.logical(te_in_promoter))]

backbone_te <- merge(backbone, prom_pa2, by = "pa_gene", all.x = TRUE)
te_by_cat <- backbone_te[!is.na(pa_prom_te), .(
  n               = .N,
  n_te_promoter   = sum(pa_prom_te),
  rate_te_promoter = round(mean(pa_prom_te), 4)
), by = coex_category][order(coex_category)]
cat("  TE-in-promoter rate by co-expressolog category:\n")
print(te_by_cat)
fwrite(te_by_cat, file.path(OUTDIR, "promoter_te_presence_by_category.tsv"), sep = "\t")

te_fisher3 <- rbindlist(lapply(coex_cats, function(cat) {
  n_tot <- backbone_te[!is.na(pa_prom_te), .N]
  n_te  <- backbone_te[!is.na(pa_prom_te), sum(pa_prom_te)]
  n_cat <- backbone_te[coex_category == cat & !is.na(pa_prom_te), .N]
  if (n_cat == 0) return(NULL)
  n_cat_te <- backbone_te[coex_category == cat & !is.na(pa_prom_te), sum(pa_prom_te)]
  fisher_enrich(n_cat_te, n_cat, n_te, n_tot, label = cat)
}))
if (nrow(te_fisher3) > 0) {
  te_fisher3[, padj := p.adjust(pvalue, method = "BH")]
  cat("  Fisher test (TE-in-promoter enrichment per category):\n")
  print(te_fisher3[, .(label, n_group, rate_positive, rate_bg, odds_ratio, pvalue, padj)])
  fwrite(te_fisher3, file.path(OUTDIR, "promoter_te_category_fisher.tsv"), sep = "\t")
}


# ═══════════════════════════════════════════════════════════════════════════════
# Item 4: Population-genetics signals × co-expressolog categories
#
# Universe: all spruce genes in the ComPlEx expression matrix (same as Item 1).
# The popgen signal genes are SD paralogs (within-species duplicates) and are
# therefore absent from the 1:1 ortholog backbone by construction — 270/272
# popgen genes are in the expression universe vs only 41/272 in the backbone.
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Item 4: Popgen signals × co-expressolog ─────────────────────────────\n")

log_source("sd_popgen_signals", file.path(SD_DIR, "sd_popgen_signals.tsv"))
popgen <- fread(file.path(SD_DIR, "sd_popgen_signals.tsv"))
setnames(popgen, "gene", "pa_gene")
cat(sprintf("  sd_popgen_signals: %d genes\n", nrow(popgen)))

# Use expression universe as background (universe_dt2 built in Item 1)
# universe_dt2: one row per gene in expression universe with coex_category
bb_pop <- merge(universe_dt2[, .(pa_gene, coex_category)],
                popgen[, .(pa_gene, signal)],
                by = "pa_gene", all.x = TRUE)
bb_pop[is.na(signal), signal := "none"]
cat(sprintf("  Popgen genes recovered in expression universe: %d / %d\n",
            sum(bb_pop$signal != "none"), nrow(popgen)))

# Co-expressolog rate summary by (compound) signal label
pop_coex_rate <- bb_pop[signal != "none", .(
  n_total        = .N,
  n_coex         = sum(coex_category != "not_coex"),
  n_conserved    = sum(coex_category == "conserved"),
  rate_coex      = round(mean(coex_category != "not_coex"), 4),
  rate_conserved = round(mean(coex_category == "conserved"), 4)
), by = signal][order(signal)]
cat("  Co-expressolog rate by popgen signal (in expression universe):\n")
cat(sprintf("  Background rate (full universe): %.4f\n",
            sum(universe_dt2$coex_category != "not_coex") / nrow(universe_dt2)))
print(pop_coex_rate)
fwrite(pop_coex_rate, file.path(OUTDIR, "popgen_coexpression_rate.tsv"), sep = "\t")

popgen_cat_table <- bb_pop[signal != "none", .N,
  by = .(signal, coex_category)][order(signal, coex_category)]
fwrite(popgen_cat_table, file.path(OUTDIR, "popgen_category_counts.tsv"), sep = "\t")

# Fisher tests: expand compound signal labels (e.g. "gwas,pav") so each gene
# appears once per signal type; test each type against expression universe
popgen_expanded <- rbindlist(lapply(seq_len(nrow(popgen)), function(i) {
  sigs <- trimws(unlist(strsplit(popgen$signal[i], ",")))
  data.table(pa_gene = popgen$pa_gene[i], signal_indiv = sigs)
}))

n_total_pop <- nrow(universe_dt2)
n_coex_pop  <- sum(universe_dt2$coex_category != "not_coex")
signals_indiv <- unique(popgen_expanded$signal_indiv)
pop_fisher <- rbindlist(lapply(signals_indiv, function(sig) {
  genes_sig <- unique(popgen_expanded[signal_indiv == sig, pa_gene])
  n_sig <- length(genes_sig)
  if (n_sig == 0) return(NULL)
  n_pos <- sum(universe_dt2[pa_gene %in% genes_sig, coex_category != "not_coex"])
  fisher_enrich(n_pos, n_sig, n_coex_pop, n_total_pop, label = sig)
}))
if (nrow(pop_fisher) > 0) {
  pop_fisher[, padj := p.adjust(pvalue, method = "BH")]
  cat("  Fisher enrichment (co-expressolog by popgen signal, vs expression universe):\n")
  print(pop_fisher[, .(label, n_group, rate_positive, rate_bg, odds_ratio, pvalue, padj)])
  fwrite(pop_fisher, file.path(OUTDIR, "popgen_category_enrichment.tsv"), sep = "\t")
}


# ── Provenance log ─────────────────────────────────────────────────────────────
prov <- rbindlist(lapply(names(sources), function(nm)
  data.table(source = nm, path = sources[[nm]]$path, mtime = sources[[nm]]$mtime)))
fwrite(prov, file.path(OUTDIR, "integration_sources.tsv"), sep = "\t")

cat(sprintf("\nAll outputs in: %s\n", OUTDIR))
for (f in sort(list.files(OUTDIR, pattern = "\\.tsv$"))) {
  n <- nrow(fread(file.path(OUTDIR, f), showProgress = FALSE))
  cat(sprintf("  %-55s  %d rows\n", f, n))
}
