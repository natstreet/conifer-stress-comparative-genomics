#!/usr/bin/env Rscript
## pav_category_bh.R — PAV enrichment across the five co-expression categories, with the
## multiple-testing correction applied WITHIN the PAV family (Benjamini-Hochberg over the five
## PAV x category Fisher tests), correcting within the family of tests that address the PAV
## hypothesis rather than the broader 15-test popgen sweep.
##
## The per-category Fisher tests (odds ratio, raw p, PAV panel size) are the committed values in
## results/integration/popgen_category_fisher.tsv (producer: enrichment_tests_extended.R). This
## script re-applies BH across ONLY the five PAV rows, so the reported not_coex adjusted P reflects
## correction within the PAV hypothesis family. The underlying Fisher test and background are
## therefore identical to the committed not_coex OR = 1.77, raw p = 0.006332; only the correction
## family differs from the popgen-wide table.
##
## Output: results/integration/pav_category_bh.tsv  (coex_category, n_pav, OR, p_raw, p_bh)

INTEG <- "results/integration"
src   <- file.path(INTEG, "popgen_category_fisher.tsv")
f <- read.delim(src, stringsAsFactors = FALSE)

pav <- f[f$signal == "pav", c("coex_category", "n_sig_total", "OR", "pvalue")]
names(pav) <- c("coex_category", "n_pav", "OR", "p_raw")
pav$p_bh <- p.adjust(pav$p_raw, method = "BH")          # BH within the five PAV tests
pav <- pav[order(pav$p_raw), ]

## Sanity gate: the not_coex raw p must reproduce the committed 0.006332 exactly; if the upstream
## Fisher test or background has changed, STOP rather than emit a silently different value.
nc_raw <- pav$p_raw[pav$coex_category == "not_coex"]
if (!isTRUE(abs(nc_raw - 0.0063322572044659) < 1e-9))
  stop(sprintf("not_coex raw p = %.10g does not reproduce the committed 0.0063322572044659; aborting.", nc_raw))

write.table(pav, file.path(INTEG, "pav_category_bh.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote", file.path(INTEG, "pav_category_bh.tsv"), "\n")
print(pav, row.names = FALSE)
nc_bh <- pav$p_bh[pav$coex_category == "not_coex"]
cat(sprintf("\nnot_coex: OR = %.3f  p_raw = %.6f  p_bh = %.5f  survives@0.05: %s\n",
            pav$OR[pav$coex_category == "not_coex"], nc_raw, nc_bh, nc_bh < 0.05))
