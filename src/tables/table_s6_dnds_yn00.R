#!/usr/bin/env Rscript
# table_s6_dnds_yn00.R — Supplementary Table S6: dN/dS across cross-species
# co-expression conservation categories.
#
# Reads the per-pair YN00 dN/dS values and joins the AUTHORITATIVE co-expression
# category from the current integration backbone (integration_backbone_1to1.tsv),
# exactly as pnps_confound_analysis.py and the Figure 5 script do. (cross_species_dnds_yn00.tsv carries
# no category column; categories are always taken from the backbone here, so they cannot fall out of
# step with the thresholded set.) dN/dS values are unaffected by category assignment. Summarises dN/dS per category
# with a two-sided Wilcoxon rank-sum test of each category against the not_coex
# reference (Benjamini-Hochberg corrected across the four contrasts).
#
# Inputs: results/integration/cross_species_dnds_yn00.tsv     (per-pair dN/dS)
#         results/integration/integration_backbone_1to1.tsv   (authoritative categories)
# Output: manuscript/supplementary/Supplementary_Table_S6_dNdS_YN00.xlsx
suppressPackageStartupMessages({ library(openxlsx) })

IN  <- "results/integration/cross_species_dnds_yn00.tsv"
BB  <- "results/integration/integration_backbone_1to1.tsv"
OUT <- "manuscript/supplementary/Supplementary_Table_S6_dNdS_YN00.xlsx"

kd <- read.table(IN, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                 quote = "")[, c("pa_gene", "ps_gene", "dN", "dS", "dNdS")]
bb <- read.table(BB, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                 quote = "")[, c("pa_gene", "ps_gene", "coex_category")]
d <- merge(bb, kd, by = c("pa_gene", "ps_gene"))
n_total <- nrow(kd)                      # pairs with a YN00 dN/dS estimate
# Analysis filter: dS in (0,5) and dN/dS < 10 (excludes saturated / unreliable estimates)
d <- d[!is.na(d$dNdS) & d$dS > 0 & d$dS < 5 & d$dNdS < 10, ]
n_kept <- nrow(d)

cats <- c(conserved = "Conserved", cold_specific = "Cold-specific",
          drought_specific = "Drought-specific", multi_tissue = "Multi-tissue",
          not_coex = "Not co-expressed")
ref <- d$dNdS[d$coex_category == "not_coex"]

# per-category summary + Wilcoxon vs not_coex
tab <- data.frame(); pvals <- c(); pkeys <- c()
for (k in names(cats)) {
  x <- d$dNdS[d$coex_category == k]
  W <- P <- NA_real_
  if (k != "not_coex") {
    w <- suppressWarnings(wilcox.test(x, ref, alternative = "two.sided", correct = TRUE))
    W <- unname(w$statistic); P <- w$p.value
    pvals <- c(pvals, P); pkeys <- c(pkeys, k)
  }
  tab <- rbind(tab, data.frame(
    category = cats[[k]], n = length(x),
    med = round(median(x), 3), q1 = round(quantile(x, .25), 3), q3 = round(quantile(x, .75), 3),
    W = W, P = P, stringsAsFactors = FALSE))
}
padj <- setNames(p.adjust(pvals, method = "BH"), pkeys)
tab$Padj <- padj[names(cats)]           # NA for not_coex

fmtp <- function(p) ifelse(is.na(p), "—", formatC(p, format = "e", digits = 1))
out <- data.frame(
  `Co-expression category`   = tab$category,
  `n (1:1 pairs)`            = tab$n,
  `Median dN/dS`             = tab$med,
  `IQR lower (Q1)`           = tab$q1,
  `IQR upper (Q3)`           = tab$q3,
  `Wilcoxon W (vs not_coex)` = ifelse(is.na(tab$W), "—", format(round(tab$W), scientific = FALSE)),
  `P (vs not_coex)`          = fmtp(tab$P),
  `P_adj (BH)`               = fmtp(tab$Padj),
  check.names = FALSE, stringsAsFactors = FALSE)

title <- paste0("Table S6. Non-synonymous to synonymous substitution rate (dN/dS, ",
  "Yang & Nielsen 2000) across cross-species co-expression conservation categories.")
footnote <- sprintf(paste0(
  "Footnote: dN/dS estimated for %s Picea abies–Pinus sylvestris 1:1 ortholog pairs under the ",
  "Yang & Nielsen (2000) maximum-likelihood codon model as implemented in PAML yn00 (version 4.10.10); ",
  "pairs with dS outside (0,5) or dN/dS>10 excluded (%s pairs retained). Categories assigned from ",
  "cross-species co-expressolog conservation across the four stress–tissue comparisons. Two-sided ",
  "Wilcoxon rank-sum tests compare each category against the not co-expressed category; P-values ",
  "Benjamini–Hochberg corrected across the four contrasts."),
  format(n_total, big.mark = ","), format(n_kept, big.mark = ","))

wb <- createWorkbook(); addWorksheet(wb, "Table S6")
writeData(wb, "Table S6", title, startRow = 1, colNames = FALSE)      # row 1 title
writeData(wb, "Table S6", out, startRow = 3, colNames = TRUE)         # blank row 2, header row 3, data 4..8
writeData(wb, "Table S6", footnote, startRow = nrow(out) + 4, colNames = FALSE)  # footnote row 9
saveWorkbook(wb, OUT, overwrite = TRUE)
invisible(system2("python3", c("manuscript/supplementary/strip_xlsx_drawings.py", OUT)))  # submission-clean xlsx
cat("Wrote", OUT, "\n"); print(out, row.names = FALSE)
