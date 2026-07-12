#!/usr/bin/env Rscript
# table_s7_threshold_sensitivity.R — Supplementary Table S7: threshold-sensitivity of the
# co-expression conservation categories and the dN/dS gradient.
#
# Reads the threshold-sensitivity summary (produced by threshold_sensitivity.py)
# and writes it as the formatted supplementary workbook with the caption row.
#
# Input : results/integration/threshold_sensitivity.tsv
# Output: manuscript/supplementary/Supplementary_Table_S7_threshold_sensitivity.xlsx
suppressPackageStartupMessages({ library(openxlsx) })

IN  <- "results/integration/threshold_sensitivity.tsv"
OUT <- "manuscript/supplementary/Supplementary_Table_S7_threshold_sensitivity.xlsx"

d <- read.table(IN, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Display tweaks used in the cited table: compact the config label and the boolean header.
d$config <- sub("C_SUM_STRESS=", "stress=", d$config, fixed = TRUE)
names(d)[names(d) == "conserved_lowest_and_not_coex_highest"] <- "conserved_lowest&not_coex_highest"
# Format the p-value column compactly (avoid float-expansion artefacts in the cells).
if ("p_cons_vs_not_coex" %in% names(d))
  d$p_cons_vs_not_coex <- formatC(as.numeric(d$p_cons_vs_not_coex), format = "e", digits = 1)

title <- paste0(
  "Table S7. Sensitivity of the co-expression conservation categories and the dN/dS gradient to the ",
  "co-expressolog evidence thresholds. Categories were re-derived from the weighted co-expressolog gene ",
  "pairs across a range of summed -log10 MaxpVal cut-offs (C_SUM for the conserved category; ",
  "C_SUM_STRESS for the stress-specific and multi-tissue categories); the baseline (structural) ",
  "configuration reproduces the published category sizes. dN/dS estimated under Yang & Nielsen (2000). ",
  "Across all configurations the conserved category retains the lowest median dN/dS and the not_coex ",
  "category the highest, the conserved-vs-not_coex difference remains highly significant (two-sided ",
  "Wilcoxon), and shared segmental duplicates remain enriched in the drought-specific category (OR>1) ",
  "and depleted from not_coex (OR<1).")

wb <- createWorkbook(); addWorksheet(wb, "Table S7")
writeData(wb, "Table S7", title, startRow = 1, colNames = FALSE)
writeData(wb, "Table S7", d, startRow = 3, colNames = TRUE)   # blank row 2, matching the cited layout
saveWorkbook(wb, OUT, overwrite = TRUE)
invisible(system2("python3", c("manuscript/supplementary/strip_xlsx_drawings.py", OUT)))  # submission-clean xlsx
cat("Wrote", OUT, "(", nrow(d), "configurations )\n")
