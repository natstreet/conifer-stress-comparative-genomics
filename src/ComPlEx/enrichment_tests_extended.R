#!/usr/bin/env Rscript
# enrichment_tests_extended.R
#
# Three additional enrichment analyses:
#   A. SD gene class × co-expressolog category  -- Fisher tests (5 categories × 2 SD classes)
#   B. Popgen signal × co-expressolog category  -- Fisher tests (5 categories × 3 signals)

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
PROJ_DIR    <- normalizePath(file.path(SCRIPT_DIR, "../.."))
SD_DIR      <- normalizePath(file.path(PROJ_DIR, ".."))
OUTDIR      <- file.path(PROJ_DIR, "results/integration")

CATEGORIES <- c("conserved", "cold_specific", "drought_specific", "multi_tissue", "not_coex")

# ── A. SD class × co-expressolog category Fisher tests ─────────────────────────
cat("=== A. SD gene class × co-expressolog category ===\n")

sd_cat  <- fread(file.path(OUTDIR, "sd_category_counts.tsv"))
# Expression-universe totals per category (sum over all SD classes)
cat_totals <- sd_cat[, .(n_cat = sum(N)), by = coex_category]
n_universe <- sum(cat_totals$n_cat)
cat(sprintf("  Expression universe: %d genes\n", n_universe))

sd_classes <- c("shared_SD", "spruce_only_SD")
sd_totals  <- sd_cat[, .(n_sd = sum(N)), by = sd_class]
# non_SD total for reference
n_non_sd <- sd_totals[sd_class == "non_SD", n_sd]

results_A <- rbindlist(lapply(sd_classes, function(sdc) {
  n_sdc <- sd_totals[sd_class == sdc, n_sd]
  n_bg  <- n_non_sd   # compare each SD class against non_SD baseline
  rbindlist(lapply(CATEGORIES, function(cat) {
    n_sdc_cat <- sd_cat[sd_class == sdc        & coex_category == cat, N]
    n_bg_cat  <- sd_cat[sd_class == "non_SD"   & coex_category == cat, N]
    if (length(n_sdc_cat) == 0) n_sdc_cat <- 0L
    if (length(n_bg_cat)  == 0) n_bg_cat  <- 0L
    mat <- matrix(c(n_sdc_cat, n_sdc - n_sdc_cat,
                    n_bg_cat,  n_bg  - n_bg_cat), nrow = 2)
    ft <- fisher.test(mat)
    data.table(
      sd_class      = sdc,
      coex_category = cat,
      n_sd          = n_sdc,
      n_in_cat      = n_sdc_cat,
      rate_sd       = round(n_sdc_cat / n_sdc, 4),
      n_nonsd       = n_bg,
      n_nonsd_cat   = n_bg_cat,
      rate_nonsd    = round(n_bg_cat / n_bg, 4),
      OR            = round(ft$estimate, 3),
      pvalue        = ft$p.value
    )
  }))
}))
results_A[, padj := p.adjust(pvalue, method = "BH")]
setorder(results_A, sd_class, pvalue)

cat("\nSD class vs non_SD baseline — Fisher OR by co-expressolog category:\n")
print(results_A[, .(sd_class, coex_category, n_in_cat, rate_sd, rate_nonsd, OR, pvalue, padj)])
fwrite(results_A, file.path(OUTDIR, "sd_category_fisher.tsv"), sep = "\t")
cat(sprintf("  Saved sd_category_fisher.tsv\n"))

# ── B. Popgen signal × co-expressolog category Fisher tests ────────────────────
cat("\n=== B. Popgen signal × co-expressolog category ===\n")

popgen_cat <- fread(file.path(OUTDIR, "popgen_category_counts.tsv"))

# Expand compound signal labels (e.g. "gwas,pav" → two rows)
popgen_cat_exp <- rbindlist(lapply(seq_len(nrow(popgen_cat)), function(i) {
  sigs <- trimws(unlist(strsplit(popgen_cat$signal[i], ",")))
  data.table(signal_indiv = sigs, coex_category = popgen_cat$coex_category[i],
             N = popgen_cat$N[i])
}))[, .(N = sum(N)), by = .(signal_indiv, coex_category)]

# Total per expanded signal type
sig_totals <- popgen_cat_exp[, .(n_sig = sum(N)), by = signal_indiv]
cat("  Expanded signal-type totals:\n")
print(sig_totals)

results_B <- rbindlist(lapply(c("pav", "selection", "gwas"), function(sig) {
  n_sig <- sig_totals[signal_indiv == sig, n_sig]
  if (length(n_sig) == 0 || n_sig == 0) return(NULL)
  rbindlist(lapply(CATEGORIES, function(cat) {
    n_sig_cat <- popgen_cat_exp[signal_indiv == sig & coex_category == cat, N]
    if (length(n_sig_cat) == 0) n_sig_cat <- 0L
    n_cat_all <- cat_totals[coex_category == cat, n_cat]
    if (length(n_cat_all) == 0) n_cat_all <- 0L
    # 2×2: [signal in cat, signal not in cat; non-signal in cat, non-signal not in cat]
    mat <- matrix(c(n_sig_cat,
                    n_sig - n_sig_cat,
                    n_cat_all - n_sig_cat,
                    n_universe - n_cat_all - (n_sig - n_sig_cat)), nrow = 2)
    ft <- fisher.test(mat)
    data.table(
      signal        = sig,
      coex_category = cat,
      n_sig_total   = n_sig,
      n_sig_in_cat  = n_sig_cat,
      rate_sig      = round(n_sig_cat / n_sig, 4),
      rate_bg       = round(n_cat_all / n_universe, 4),
      OR            = round(ft$estimate, 3),
      pvalue        = ft$p.value
    )
  }))
}))

if (!is.null(results_B) && nrow(results_B) > 0) {
  results_B[, padj := p.adjust(pvalue, method = "BH")]
  setorder(results_B, signal, pvalue)
  cat("\nPopgen signal × co-expressolog category enrichment:\n")
  print(results_B[, .(signal, coex_category, n_sig_in_cat, rate_sig, rate_bg, OR, pvalue, padj)])
  fwrite(results_B, file.path(OUTDIR, "popgen_category_fisher.tsv"), sep = "\t")
  cat(sprintf("  Saved popgen_category_fisher.tsv\n"))
} else {
  cat("  No results produced.\n")
}

cat(sprintf("\nOutputs in: %s\n", OUTDIR))
