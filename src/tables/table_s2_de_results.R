#!/usr/bin/env Rscript
# table_s2_de_results.R — Supplementary Table S2: differential-expression results for
# significantly DE genes, one worksheet per stress comparison.
#
# For each of the eight comparisons (species x stress x organ) DESeq2 is run on the
# deposited DESeqDataSet and results are extracted for every timepoint contrast
# (results(dds, name=<Condition_..._vs_20C_0h>, filter=rowMedians(counts(dds,normalized=TRUE)))),
# matching the Figure 1 / DE_all_* definition. A gene is included on a sheet if it is
# significant (padj < 0.01 AND |log2FoldChange| >= 2) in at least one contrast of that
# comparison; the sheet lists its baseMean and, per contrast, log2FoldChange and padj.
#
# Input : data/dds/dds_{SCN,SCR,SDN,SDR,PCN,PCR,PDN,PDR}.rda  (see SOURCES.tsv)
# Output: manuscript/supplementary/Supplementary_Table_S2_DE_results.xlsx
suppressPackageStartupMessages({ library(DESeq2); library(openxlsx) })

DDS_DIR <- "data/dds"
OUT <- "manuscript/supplementary/Supplementary_Table_S2_DE_results.xlsx"

conds <- c(dds_SCN="Spruce_Cold_Needle",   dds_SCR="Spruce_Cold_Root",
           dds_SDN="Spruce_Drought_Needle", dds_SDR="Spruce_Drought_Root",
           dds_PCN="Pine_Cold_Needle",     dds_PCR="Pine_Cold_Root",
           dds_PDN="Pine_Drought_Needle",  dds_PDR="Pine_Drought_Root")

wb <- createWorkbook()
for (f in names(conds)) {
  e <- new.env(); load(file.path(DDS_DIR, paste0(f, ".rda")), envir = e)
  dds <- get(ls(e)[1], e)
  dds <- DESeq(dds, quiet = TRUE)
  filtered <- rowMedians(counts(dds, normalized = TRUE))
  rn <- setdiff(resultsNames(dds), "Intercept")

  tab <- NULL; sig <- character(0); basemean <- NULL
  for (cn in rn) {
    r <- as.data.frame(results(dds, name = cn, filter = filtered))
    if (is.null(basemean)) basemean <- data.frame(gene = rownames(r),
                                                   baseMean = round(r$baseMean, 2))
    short <- sub("^Condition_", "", sub("_vs_.*$", "", cn))
    col <- data.frame(gene = rownames(r),
                      a = round(r$log2FoldChange, 3),
                      b = signif(r$padj, 3))
    names(col) <- c("gene", paste0("log2FC_", short), paste0("padj_", short))
    tab <- if (is.null(tab)) col else merge(tab, col, by = "gene", all = TRUE)
    sig <- union(sig, rownames(r)[which(r$padj < 0.01 & abs(r$log2FoldChange) >= 2)])
  }
  tab <- merge(basemean, tab, by = "gene")
  tab <- tab[tab$gene %in% sig, ]
  tab <- tab[order(tab$gene), ]
  addWorksheet(wb, conds[f])
  writeData(wb, conds[f], tab)
  cat(sprintf("  %-22s %5d significant genes across %d contrasts\n",
              conds[f], nrow(tab), length(rn)))
}
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, OUT, overwrite = TRUE)
invisible(system2("python3", c("manuscript/supplementary/strip_xlsx_drawings.py", OUT)))  # submission-clean xlsx
cat("Wrote", OUT, "\n")
