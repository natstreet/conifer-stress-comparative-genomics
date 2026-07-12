#!/usr/bin/env Rscript
## orthogroup_direction_concordance.R — directional concordance of cross-species DEG orthogroups
## (Results [27]). Producer for the "% directionally concordant" statistic and its binomial tests.
##
## Universe: the SAME shared DEG-orthogroup set as the 'overall' row of Figure 2c — orthogroups
## (>=1 Picea AND >=1 Pinus gene in doc/genes_ortholog_categories.tsv) that are DEG-hit in BOTH
## species pooled across all four stress-tissue comparisons.
##
## Direction rule (stated explicitly): for each species, a gene's per-contrast log2FoldChange is
## taken from every non-Intercept DESeq2 coefficient of its stress dds; a gene counts toward an
## orthogroup only in contrasts where it passes the DEG threshold (padj < 0.01 & |log2FC| >= 2), i.e.
## the exact set that defines DE_all_*. Each orthogroup's species direction = sign of the SUM of those
## passing log2FoldChange values over all member genes; a sum of exactly 0 (perfectly mixed) -> NA and
## the orthogroup is EXCLUDED from the concordance denominator.
##
## Tests: (i) binomial vs p = 0.5 (the manuscript's null); (ii) binomial vs the observed marginal null
## — expected concordance if the two species' directions were independent given each species' own
## up/down proportions: E = f_sp_up*f_pi_up + f_sp_down*f_pi_down.
##
## Run from AbioticStressConifers/.  Output: results/integration/orthogroup_direction_concordance.tsv
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats); library(data.table) })

OGX  <- "doc/genes_ortholog_categories.tsv"
DEGD <- "data/DEG_lists"
DDSD <- "data/dds"
INTEG <- "results/integration"; dir.create(INTEG, showWarnings = FALSE, recursive = TRUE)

## ---- universe + gene->orthogroup maps (identical to figure2c_orthogroup_overlap.R) ----
x <- read.table(OGX, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
sp_map <- x[x$species == "Picea_abies",     c("gene", "Ortholog_Group")]
pi_map <- x[x$species == "Pinus_sylvestris", c("gene", "Ortholog_Group")]
universe <- intersect(unique(sp_map$Ortholog_Group), unique(pi_map$Ortholog_Group))
g2og <- function(map) split(map$Ortholog_Group, map$gene)
g2og_sp <- g2og(sp_map); g2og_pi <- g2og(pi_map)

deg_ids <- function(tag) { e <- new.env(); load(file.path(DEGD, sprintf("DE_all_%s_01_2L2FC.RData", tag)), e); as.character(get(ls(e)[1], e)) }
ogs_of  <- function(genes, map) intersect(unique(unlist(map[genes[genes %in% names(map)]], use.names = FALSE)), universe)

## per-gene signed log2FC summed over the contrasts where the gene passes the DEG threshold, across a
## species' stress datasets. Returns a named numeric (gene -> summed passing log2FC).
gene_signed_lfc <- function(tags) {
  acc <- list()
  for (tag in tags) {
    e <- new.env(); load(file.path(DDSD, sprintf("dds_%s.rda", tag)), e)
    dds <- DESeq(e[[ls(e)[1]]], quiet = TRUE)
    filt <- rowMedians(counts(dds, normalized = TRUE))
    for (nm in setdiff(resultsNames(dds), "Intercept")) {
      r <- as.data.frame(results(dds, name = nm, filter = filt))
      hit <- r[!is.na(r$padj) & r$padj < 0.01 & abs(r$log2FoldChange) >= 2, "log2FoldChange", drop = FALSE]
      if (nrow(hit)) acc[[length(acc) + 1]] <- data.frame(gene = rownames(hit), lfc = hit$log2FoldChange)
    }
  }
  d <- rbindlist(acc)
  s <- d[, .(lfc = sum(lfc)), by = gene]        # sum of passing log2FC over all contrasts
  setNames(s$lfc, s$gene)
}

## per-orthogroup species direction. Each member DEG gene has a sign = sign of its summed passing
## log2FC over timepoints. An orthogroup's direction is 'up' only if ALL its member DEG genes are up,
## 'down' only if ALL are down; if the member genes disagree (mixed) it is ambiguous -> NA and the
## orthogroup is excluded from the concordance denominator.
og_direction <- function(ogset, glfc, map) {
  gene_og <- data.table(gene = names(glfc), sgn = sign(as.numeric(glfc)))
  gene_og <- gene_og[gene %in% names(map) & sgn != 0]
  gene_og[, og := vapply(map[gene], function(v) v[1], character(1))]   # 1:1 gene->OG in this table
  agg <- gene_og[og %in% ogset, .(up = sum(sgn > 0), dn = sum(sgn < 0)), by = og]
  d <- setNames(rep(NA_integer_, length(ogset)), ogset)
  d[agg[up > 0 & dn == 0, og]] <-  1L     # all member DEG genes up-regulated
  d[agg[dn > 0 & up == 0, og]] <- -1L     # all member DEG genes down-regulated
  d                                        # mixed (up>0 & dn>0) or no signed member -> NA (excluded)
}

cat("Computing per-gene signed log2FC (spruce SCN/SCR/SDN/SDR)...\n")
sp_lfc <- gene_signed_lfc(c("SCN","SCR","SDN","SDR"))
cat("Computing per-gene signed log2FC (pine PCN/PCR/PDN/PDR)...\n")
pi_lfc <- gene_signed_lfc(c("PCN","PCR","PDN","PDR"))

## shared DEG-orthogroups (overall), same as Fig2c
sp_ogs <- ogs_of(unique(c(deg_ids("SCN"), deg_ids("SCR"), deg_ids("SDN"), deg_ids("SDR"))), g2og_sp)
pi_ogs <- ogs_of(unique(c(deg_ids("PCN"), deg_ids("PCR"), deg_ids("PDN"), deg_ids("PDR"))), g2og_pi)
shared <- intersect(sp_ogs, pi_ogs)
cat("shared DEG-orthogroups (overall) =", length(shared), "\n")

sp_dir <- og_direction(shared, sp_lfc, g2og_sp)
pi_dir <- og_direction(shared, pi_lfc, g2og_pi)

ok <- !is.na(sp_dir) & !is.na(pi_dir)
denom     <- sum(ok)
up_up     <- sum(sp_dir[ok] ==  1 & pi_dir[ok] ==  1)
down_down <- sum(sp_dir[ok] == -1 & pi_dir[ok] == -1)
concord   <- up_up + down_down
discord   <- denom - concord
pct       <- 100 * concord / denom

## marginal (independence) null over the both-unambiguous set
f_sp_up <- mean(sp_dir[ok] == 1); f_pi_up <- mean(pi_dir[ok] == 1)
E_marg  <- f_sp_up * f_pi_up + (1 - f_sp_up) * (1 - f_pi_up)

p_half <- binom.test(concord, denom, p = 0.5)$p.value
p_marg <- binom.test(concord, denom, p = E_marg)$p.value

out <- data.table(
  metric = c("n_both_unambiguous","n_concordant","pct_concordant","n_up_up","n_down_down",
             "n_discordant","down_up_fold","expected_concordance_marginal_null","binom_P_vs_0.5","binom_P_vs_marginal"),
  value  = c(denom, concord, round(pct,1), up_up, down_down, discord, round(down_down/up_up,2),
             round(E_marg,4), p_half, p_marg))
fwrite(out, file.path(INTEG, "orthogroup_direction_concordance.tsv"), sep = "\t")
cat("\n=== directional concordance (shared DEG-orthogroups) ===\n"); print(out)
cat(sprintf("\nWrote %s\n", file.path(INTEG, "orthogroup_direction_concordance.tsv")))
