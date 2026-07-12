#!/usr/bin/env Rscript
# degree_dnds_correlation.R — is network hubness associated with stronger purifying selection?
# Correlates each Norway spruce gene's mean ComPlEx co-expression network degree with the
# cross-species YN00 dN/dS of its 1:1 ortholog pair. Pure re-analysis of two existing outputs
# (network_degree.tsv, cross_species_dnds_yn00.tsv) — no new data.
# Output: results/integration/degree_dnds_correlation.tsv
suppressPackageStartupMessages({ library(data.table); library(stats) })

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
INT <- normalizePath(file.path(SCRIPT_DIR, "../../results/integration"))

degree <- fread(file.path(INT, "network_degree.tsv"))               # Genes, mean_degree
setnames(degree, "Genes", "pa_gene")
yn00   <- fread(file.path(INT, "cross_species_dnds_yn00.tsv"))[, .(pa_gene, dNdS)]

m <- merge(degree, yn00, by = "pa_gene")
# Require finite, non-zero dN/dS (drop pairs with no measurable non-synonymous divergence, dNdS == 0).
m <- m[is.finite(mean_degree) & is.finite(dNdS) & dNdS > 0]
cat(sprintf("Usable genes (mean_degree + finite dN/dS): %d\n", nrow(m)))

sp  <- suppressWarnings(cor.test(m$mean_degree, m$dNdS, method = "spearman"))
rho <- unname(sp$estimate); P <- sp$p.value
cat(sprintf("Spearman rho = %.3f, P = %.2e\n", rho, P))

m[, deg_q := cut(mean_degree, breaks = quantile(mean_degree, 0:4/4),
                 include.lowest = TRUE, labels = paste0("Q", 1:4))]
qmed <- m[, .(n = .N, median_dNdS = round(median(dNdS), 4)), by = deg_q][order(deg_q)]
cat("Median dN/dS by network-degree quartile (Q1 low degree -> Q4 hubs):\n"); print(qmed)

out <- rbind(
  data.table(statistic = c("n_usable", "spearman_rho", "spearman_P"),
             value = c(nrow(m), round(rho, 4), signif(P, 3))),
  data.table(statistic = paste0(qmed$deg_q, "_median_dNdS"), value = qmed$median_dNdS)
)
fwrite(out, file.path(INT, "degree_dnds_correlation.tsv"), sep = "\t")
cat(sprintf("\nWrote %s\n", file.path(INT, "degree_dnds_correlation.tsv")))
