#!/usr/bin/env Rscript
## pine_sd_category_enrichment.R — reciprocal PINE-axis SD enrichment (co-author review; makes the SD
## conservation axis symmetric). Mirrors the spruce analysis (integration_analysis.R + enrichment_tests_
## extended.R section A) exactly, but on the PINE side: the co-expressolog category is taken from the PINE
## gene (Species2) of each cross-species co-expressolog, the universe is the pine expression universe, and
## the SD classes are the pine gene classes from pine_sd_classify.py (shared_SD / pine_only_SD).
##
## Category rule (identical to spruce): conserved if (conserved AND NegLog10CliqueSum >= 10), else
## cold_specific, else drought_specific, else multi_tissue; genes with no co-expressolog are not_coex.
## Fisher: each pine SD class vs the pine non_SD baseline, per co-expression category (BH-adjusted).
##
## Inputs (all committed / deposited): results/integration/pine_sd_gene_class.tsv (pine_sd_classify.py),
##   results/ComPlEx/RData/weighted_gene_pairs.tsv, data/expression/PC_expression.txt + PD_expression.txt.
## Outputs: results/integration/pine_sd_category_counts.tsv, pine_sd_coexpression_rate.tsv,
##          pine_sd_category_fisher.tsv  (the pine [48] ORs/Padj + Fig 5b pine panel).
suppressPackageStartupMessages({ library(data.table) })

INTEG <- "results/integration"
CATEGORIES <- c("conserved", "cold_specific", "drought_specific", "multi_tissue", "not_coex")

## ── pine expression universe (mirror of spruce SC+SD) ──
pu1 <- fread("data/expression/PC_expression.txt", select = 1)[[1]]
pu2 <- fread("data/expression/PD_expression.txt", select = 1)[[1]]
pine_universe <- unique(c(pu1, pu2))
cat(sprintf("  Pine expression universe: %d genes\n", length(pine_universe)))

## ── co-expressolog category per PINE gene (Species2), best-pval pair, same thresholds as spruce ──
wp <- fread("results/ComPlEx/RData/weighted_gene_pairs.tsv")
setnames(wp, c("Species2"), c("ps_gene"))
wp_best <- wp[order(best_pval)][!duplicated(ps_gene)]        # best pair per pine gene
wp_gene_cat <- wp_best[, .(ps_gene, coex_category = fifelse(
    as.logical(conserved) & NegLog10CliqueSum >= 10, "conserved",
    fifelse(as.logical(cold_specific),    "cold_specific",
    fifelse(as.logical(drought_specific), "drought_specific",
                                          "multi_tissue"))))]

## ── pine SD gene classes (shared_SD / pine_only_SD) from the classifier ──
sd_gene_tbl <- fread(file.path(INTEG, "pine_sd_gene_class.tsv"))    # ps_gene, sd_class

## ── universe × SD class × category (mirror integration_analysis universe_dt2) ──
u <- data.table(ps_gene = pine_universe)
u <- merge(u, sd_gene_tbl,  by = "ps_gene", all.x = TRUE)
u[is.na(sd_class), sd_class := "non_SD"]
u <- merge(u, wp_gene_cat, by = "ps_gene", all.x = TRUE)
u[is.na(coex_category), coex_category := "not_coex"]

cat(sprintf("  Pine SD genes in universe: shared_SD=%d  pine_only_SD=%d  non_SD=%d\n",
            sum(u$sd_class == "shared_SD"), sum(u$sd_class == "pine_only_SD"), sum(u$sd_class == "non_SD")))

## co-expressolog rate by class (mirror sd_coexpression_rate.tsv)
u[, has_coex := coex_category != "not_coex"]
rate <- u[, .(n_total = .N, n_coex = sum(has_coex), rate_coex = round(mean(has_coex), 4)), by = sd_class][order(sd_class)]
fwrite(rate, file.path(INTEG, "pine_sd_coexpression_rate.tsv"), sep = "\t")

## counts table (mirror sd_category_counts.tsv)
counts <- u[, .N, by = .(sd_class, coex_category)][order(sd_class, coex_category)]
fwrite(counts, file.path(INTEG, "pine_sd_category_counts.tsv"), sep = "\t")

## ── Fisher: each pine SD class vs non_SD baseline, per category (identical to enrichment_tests_extended A) ──
sd_totals <- counts[, .(n_sd = sum(N)), by = sd_class]
n_non_sd  <- sd_totals[sd_class == "non_SD", n_sd]
res <- rbindlist(lapply(c("shared_SD", "pine_only_SD"), function(sdc) {
  n_sdc <- sd_totals[sd_class == sdc, n_sd]; if (length(n_sdc) == 0) return(NULL)
  rbindlist(lapply(CATEGORIES, function(cat) {
    n_sdc_cat <- counts[sd_class == sdc      & coex_category == cat, N]
    n_bg_cat  <- counts[sd_class == "non_SD" & coex_category == cat, N]
    if (length(n_sdc_cat) == 0) n_sdc_cat <- 0L
    if (length(n_bg_cat)  == 0) n_bg_cat  <- 0L
    mat <- matrix(c(n_sdc_cat, n_sdc - n_sdc_cat, n_bg_cat, n_non_sd - n_bg_cat), nrow = 2)
    ft  <- fisher.test(mat)
    data.table(sd_class = sdc, coex_category = cat, n_sd = n_sdc, n_in_cat = n_sdc_cat,
               rate_sd = round(n_sdc_cat / n_sdc, 4), n_nonsd = n_non_sd, n_nonsd_cat = n_bg_cat,
               rate_nonsd = round(n_bg_cat / n_non_sd, 4), OR = round(ft$estimate, 3), pvalue = ft$p.value)
  }))
}))
res[, padj := p.adjust(pvalue, method = "BH")]
setorder(res, sd_class, pvalue)
fwrite(res, file.path(INTEG, "pine_sd_category_fisher.tsv"), sep = "\t")

cat("\n  Pine SD class vs non_SD — Fisher OR by co-expression category:\n")
print(res[, .(sd_class, coex_category, n_in_cat, rate_sd, rate_nonsd, OR, pvalue = signif(pvalue, 3), padj = signif(padj, 3))])
cat(sprintf("\n  Wrote pine_sd_category_fisher.tsv / _counts.tsv / _coexpression_rate.tsv\n"))
