#!/usr/bin/env Rscript
# integrative_conservation_model.R
#
# What independently predicts the breadth of cross-species co-expression conservation?
# Instead of the many separate univariate enrichment tests, model conservation breadth
# (the per-gene number of stress-tissue comparisons with a significant conserved
# co-expressolog, 0-4) jointly against the genomic and evolutionary features.
#
# Because sequence-constraint (dN/dS) is defined only for 1:1 orthologs while the
# genome-architecture features (segmental duplication, PAV) are defined on the
# duplicate genes that lack 1:1 orthologs, the two feature sets sit on essentially
# disjoint gene sets and are modelled separately (complementary, not competing):
#   Model A (all co-expression-universe genes): degree, expression, SD status, PAV.
#   Model B (1:1 ortholog subset): dN/dS, promoter-TE divergence, degree, expression.
#
# Primary model: proportional-odds ordinal regression (MASS::polr; odds ratios per
# standardised predictor), with the proportional-odds assumption tested (Brant).
# Robustness: standardised linear regression and random-forest permutation importance.
#
# Run from the AbioticStressConifers directory.
# Output: results/integration/integrative_conservation_model.tsv (+ _importance.tsv)

USERLIB <- file.path(Sys.getenv("HOME"), "Library/R/arm64/4.6/library")
if (dir.exists(USERLIB)) .libPaths(c(USERLIB, .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(MASS); library(ranger); library(brant) })

args <- commandArgs(trailingOnly = FALSE)
sp0  <- sub("--file=", "", args[grep("--file=", args)])
PROJ <- normalizePath(file.path(if (length(sp0) > 0) dirname(normalizePath(sp0)) else getwd(), "../.."))
I    <- file.path(PROJ, "results/integration")

# ── Assemble per-gene features ───────────────────────────────────────────────
wp <- fread(file.path(PROJ, "results/ComPlEx/RData/weighted_gene_pairs.tsv"))
wp[, breadth := cold_needle_present + cold_root_present + drought_needle_present + drought_root_present]
d <- wp[, .(breadth = max(breadth)), by = .(pa_gene = Species1)]           # co-expression universe
deg <- fread(file.path(I, "network_degree.tsv")); setnames(deg, c("pa_gene", "degree"))
d <- merge(d, deg, by = "pa_gene", all.x = TRUE)
sc <- fread(file.path(PROJ, "data/expression/SC_expression.txt")); setnames(sc, 1, "pa_gene")
sd_ <- fread(file.path(PROJ, "data/expression/SD_expression.txt")); setnames(sd_, 1, "pa_gene")
ex <- merge(sc[, .(pa_gene, e1 = rowMeans(.SD)), .SDcols = setdiff(names(sc), "pa_gene")],
            sd_[, .(pa_gene, e2 = rowMeans(.SD)), .SDcols = setdiff(names(sd_), "pa_gene")], by = "pa_gene", all = TRUE)
ex[, expr := rowMeans(cbind(e1, e2), na.rm = TRUE)]
d <- merge(d, ex[, .(pa_gene, expr)], by = "pa_gene", all.x = TRUE)
sdf <- fread(file.path(I, "sd_pair_features.tsv"))
sdg <- unique(rbind(sdf[, .(pa_gene = pa_gene1, sd_class)], sdf[, .(pa_gene = pa_gene2, sd_class)]))[, .(sd_class = sd_class[1]), by = pa_gene]
d <- merge(d, sdg, by = "pa_gene", all.x = TRUE); d[is.na(sd_class), sd_class := "none"]
sig <- fread(file.path(PROJ, "sd_popgen_signals.tsv"))
d[, pav := as.integer(pa_gene %in% sig[signal %like% "pav", gene])]
yn <- unique(fread(file.path(I, "cross_species_dnds_yn00.tsv"))[dS > 0 & dS < 5 & dNdS < 10, .(pa_gene, dNdS)], by = "pa_gene")
te <- unique(fread(file.path(I, "cross_species_promoter_te_jaccard.tsv"))[, .(pa_gene, te_jaccard = cross_species_te_jaccard)], by = "pa_gene")
d <- merge(d, yn, by = "pa_gene", all.x = TRUE); d <- merge(d, te, by = "pa_gene", all.x = TRUE)
d <- d[!is.na(degree) & !is.na(expr)]
z <- function(x) as.numeric(scale(x))
d[, `:=`(zdeg = z(degree), zexpr = z(expr),
         shared_SD = as.integer(sd_class == "shared_SD"),
         spruce_SD = as.integer(sd_class == "spruce_only_SD"),
         bf = factor(breadth, levels = 0:4, ordered = TRUE))]

# ── Fit one model set (ordinal + linear + RF + PO test) ──────────────────────
fit_set <- function(dat, preds, label) {
  f <- reformulate(preds, "bf")
  po <- polr(f, data = dat, Hess = TRUE)
  ct <- coef(summary(po))
  co <- ct[seq_along(preds), 1]; se <- ct[seq_along(preds), 2]
  or <- data.table(model = label, term = rownames(ct)[seq_along(preds)],
                   odds_ratio = exp(co), ci_lo = exp(co - 1.96 * se), ci_hi = exp(co + 1.96 * se),
                   p_value = pnorm(abs(ct[seq_along(preds), 3]), lower.tail = FALSE) * 2)
  lm1  <- lm(reformulate(preds, "breadth"), data = dat)
  lsum <- summary(lm1)$coefficients
  lci  <- confint(lm1)
  or[, lm_std_beta := round(coef(lm1)[term], 4)]
  or[, lm_ci_lo   := round(lci[term, 1], 4)]
  or[, lm_ci_hi   := round(lci[term, 2], 4)]
  or[, lm_p       := signif(lsum[term, 4], 3)]
  brant_p <- tryCatch({ b <- brant(po); b[rownames(b) %in% preds, "probability"] }, error = function(e) NA)
  or[, brant_PO_p := round(as.numeric(brant_p[term]), 4)]
  or[, `:=`(odds_ratio = round(odds_ratio, 3), ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
            p_value = signif(p_value, 3))]
  # Random-forest permutation importance (regression on numeric breadth). A random forest is stochastic
  # (bootstrap resampling + random feature selection + the importance permutation), so a single fit gives
  # a noisy, seed-dependent importance. We therefore report the MEAN over 100 forests and quantify the
  # run-to-run (Monte Carlo) spread as the SD, plotted as error bars in Figure S5. Two deliberate choices:
  #   * num.trees = 1000 (raised from a single 500-tree fit) so each forest is closer to converged and the
  #     per-forest Monte Carlo noise is smaller.
  #   * a FIXED, enumerated seed set (1:100) rather than an unseeded/random draw. This is what makes the
  #     result DETERMINISTICALLY REPRODUCIBLE — anyone re-running gets byte-identical importances — while
  #     still exposing the across-run variation. It also means no single "favourable" seed is chosen: the
  #     figure shows the whole ensemble of 100 seeds, not one hand-picked run.
  rf_f <- reformulate(preds, "breadth")
  imp_mat <- vapply(1:100, function(s)
    ranger(rf_f, data = dat, importance = "permutation", num.trees = 1000, seed = s)$variable.importance[preds],
    numeric(length(preds)))
  imp <- data.table(model = label, term = preds,
                    rf_importance_mean = round(rowMeans(imp_mat), 5),
                    rf_importance_sd   = round(apply(imp_mat, 1, sd), 5))
  list(coef = or, imp = imp, n = nrow(dat), lm_r2 = round(summary(lm1)$adj.r.squared, 4))
}

A <- fit_set(d, c("zdeg", "zexpr", "shared_SD", "spruce_SD", "pav"),
             "A: genome-wide (all genes)")
dB <- d[!is.na(dNdS) & !is.na(te_jaccard)]
dB[, `:=`(zdnds = z(dNdS), zte = z(te_jaccard))]
B <- fit_set(dB, c("zdnds", "zte", "zdeg", "zexpr"), "B: 1:1 sequence-constraint")

coefs <- rbind(A$coef, B$coef); imps <- rbind(A$imp, B$imp)
fwrite(coefs, file.path(I, "integrative_conservation_model.tsv"), sep = "\t")
fwrite(imps,  file.path(I, "integrative_conservation_model_importance.tsv"), sep = "\t")
cat(sprintf("Model A: n=%d, linear adj R2=%.3f\nModel B: n=%d, linear adj R2=%.3f\n",
            A$n, A$lm_r2, B$n, B$lm_r2))
cat("\n=== Ordinal odds ratios (primary), linear std-beta + Brant PO-test p ===\n")
print(coefs, row.names = FALSE)
cat("\n=== Random-forest permutation importance (mean +/- SD over 100 forests) ===\n")
print(imps[order(model, -rf_importance_mean)], row.names = FALSE)
