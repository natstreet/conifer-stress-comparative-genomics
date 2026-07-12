#!/usr/bin/env Rscript
# conserved_coexpression_dynamics.R  — Figure 4 analysis
#
# Biological question (per stress): (a) are there sets of conserved cross-species
# co-expressologs that respond to the stress, and (b) do those genes show conserved
# expression PROFILES between Norway spruce and Scots pine through the stress
# time-course? And are the same or different orthogroups involved for cold vs drought?
#
# Approach (avoids an arbitrary cluster number and tests profile conservation directly):
#   * Set per stress = co-expressologs with strong conserved cross-species co-expression
#     under the stress (DroughtSum >= 60 / ColdSum >= 40) that are differentially
#     expressed in BOTH species under that stress. Specificity is NOT forced, so an
#     orthogroup responsive under both stresses appears in both sets (real overlap).
#   * Profile conservation = Pearson correlation between the spruce and the pine
#     median expression trajectory across the stress stages, per tissue; compared to
#     a shuffled-partner null.
#   * Response direction = sign of (mean late-stress stages - baseline) in spruce,
#     giving up- vs down-regulated groups for functional (GO) interpretation.
#
# Outputs (results/integration/fig4/):
#   conserved_coexpressolog_sets.tsv        one row per orthogroup x stress
#   conserved_coexpression_profile_summary.tsv
#   stress_set_overlap.tsv
# Run from the AbioticStressConifers directory.

suppressPackageStartupMessages({ library(data.table) })

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
SCRIPT_DIR  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
PROJ        <- normalizePath(file.path(SCRIPT_DIR, "../.."))
DEGL <- file.path(PROJ, "data/DEG_lists")
EXPR <- file.path(PROJ, "data/expression/expr_median")
RDAT <- file.path(PROJ, "results/ComPlEx/RData")
OUT  <- file.path(PROJ, "results/integration/fig4")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Stress definitions: co-expression evidence column, threshold, ordered stages,
# baseline stage, late-stress stages, expression-matrix prefixes, DE conditions.
STRESS <- list(
  drought = list(sumcol = "DroughtSum", thr = 60,
                 stages = c("FC80","FC60","FC40","FC30","FC30d7","Collapsed","Collapsed2d","Rehydrated"),
                 base = "FC80", late = c("Collapsed","Collapsed2d"),
                 sp = "SD", pp = "PD", de_s = c("SDN","SDR"), de_p = c("PDN","PDR")),
  cold    = list(sumcol = "ColdSum", thr = 40,
                 stages = c("20C_0h","5C_6h","5C_24h","5C_3d","5C_10d","neg5C_6h","neg5C_24h","neg5C_3d","neg5C_10d"),
                 base = "20C_0h", late = c("neg5C_3d","neg5C_10d"),
                 sp = "SC", pp = "PC", de_s = c("SCN","SCR"), de_p = c("PCN","PCR")))

ld <- function(cond) {
  e <- new.env(); load(file.path(DEGL, paste0("DE_all_", cond, "_01_2L2FC.RData")), envir = e)
  get(ls(e)[grep("DE_all", ls(e))], envir = e)
}
read_med <- function(f) {
  m <- read.table(f, header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)
  rownames(m) <- m[["Genes"]]; m[, -1]
}
wp <- fread(file.path(RDAT, "weighted_gene_pairs.tsv"))

# trajectory of one gene across the stage axis of one matrix/tissue
traj <- function(mat, pref, gene, stages, tissue) {
  cols <- paste0(pref, "_", stages, "_", tissue)
  if (!gene %in% rownames(mat)) return(rep(NA_real_, length(stages)))
  suppressWarnings(as.numeric(mat[gene, cols]))
}

analyse_stress <- function(name) {
  cfg <- STRESS[[name]]
  pa_de <- unique(unlist(lapply(cfg$de_s, ld)))
  ps_de <- unique(unlist(lapply(cfg$de_p, ld)))
  sel <- wp[get(cfg$sumcol) >= cfg$thr & Species1 %in% pa_de & Species2 %in% ps_de]
  sel <- sel[order(-get(cfg$sumcol))][!duplicated(OrthoGroup)]      # one best pair per orthogroup
  cat(sprintf("%s: %d conserved co-expressolog orthogroups (DE both species)\n", name, nrow(sel)))

  S <- read_med(file.path(EXPR, paste0(cfg$sp, "med_expression.txt")))
  P <- read_med(file.path(EXPR, paste0(cfg$pp, "med_expression.txt")))

  prof <- function(gene_pa, gene_ps, tissue) {
    s <- traj(S, cfg$sp, gene_pa, cfg$stages, tissue)
    p <- traj(P, cfg$pp, gene_ps, cfg$stages, tissue)
    if (sum(!is.na(s)) < 4 || sum(!is.na(p)) < 4 ||
        sd(s, na.rm = TRUE) == 0 || sd(p, na.rm = TRUE) == 0) return(NA_real_)
    suppressWarnings(cor(s, p, method = "pearson", use = "complete.obs"))
  }
  corN <- mapply(prof, sel$Species1, sel$Species2, MoreArgs = list(tissue = "N"))
  corR <- mapply(prof, sel$Species1, sel$Species2, MoreArgs = list(tissue = "R"))

  # shuffled-partner null (seeded)
  set.seed(42); shuf <- sample(sel$Species2)
  nullN <- mapply(prof, sel$Species1, shuf, MoreArgs = list(tissue = "N"))

  # response direction from spruce trajectory (mean late-stress - baseline, averaged over tissues)
  dir_val <- sapply(sel$Species1, function(g) {
    sN <- traj(S, cfg$sp, g, cfg$stages, "N"); names(sN) <- cfg$stages
    sR <- traj(S, cfg$sp, g, cfg$stages, "R"); names(sR) <- cfg$stages
    dN <- mean(sN[cfg$late], na.rm = TRUE) - sN[cfg$base]
    dR <- mean(sR[cfg$late], na.rm = TRUE) - sR[cfg$base]
    mean(c(dN, dR), na.rm = TRUE)
  })
  direction <- ifelse(is.na(dir_val), NA, ifelse(dir_val > 0, "up", "down"))

  data.table(stress = name, OrthoGroup = sel$OrthoGroup,
             pa_gene = sel$Species1, ps_gene = sel$Species2,
             coexp_sum = sel[[cfg$sumcol]],
             profile_cor_needle = round(corN, 4), profile_cor_root = round(corR, 4),
             profile_cor_null_needle = round(nullN, 4),
             response_dir_value = round(dir_val, 3), direction = direction)
}

res <- rbindlist(lapply(names(STRESS), analyse_stress))
fwrite(res, file.path(OUT, "conserved_coexpressolog_sets.tsv"), sep = "\t")

# ── Summary: profile conservation vs null, per stress ────────────────────────
summ <- res[, .(
  n = .N,
  median_cor_needle = round(median(profile_cor_needle, na.rm = TRUE), 3),
  pct_positive_needle = round(100 * mean(profile_cor_needle > 0, na.rm = TRUE), 1),
  median_cor_root = round(median(profile_cor_root, na.rm = TRUE), 3),
  pct_positive_root = round(100 * mean(profile_cor_root > 0, na.rm = TRUE), 1),
  median_cor_null = round(median(profile_cor_null_needle, na.rm = TRUE), 3),
  n_up = sum(direction == "up", na.rm = TRUE),
  n_down = sum(direction == "down", na.rm = TRUE)
), by = stress]
# Wilcoxon observed vs null (needle) per stress
summ[, wilcox_p := sapply(stress, function(s) {
  x <- res[stress == s, profile_cor_needle]; y <- res[stress == s, profile_cor_null_needle]
  format(wilcox.test(x, y)$p.value, digits = 3)
})]
fwrite(summ, file.path(OUT, "conserved_coexpression_profile_summary.tsv"), sep = "\t")
cat("\nProfile-conservation summary:\n"); print(summ)

# ── Stress-set overlap (are the same orthogroups involved?) ──────────────────
dOG <- res[stress == "drought", OrthoGroup]; cOG <- res[stress == "cold", OrthoGroup]
ov <- data.table(
  drought_total = length(dOG), cold_total = length(cOG),
  drought_only = length(setdiff(dOG, cOG)),
  shared_both  = length(intersect(dOG, cOG)),
  cold_only    = length(setdiff(cOG, dOG)))
fwrite(ov, file.path(OUT, "stress_set_overlap.tsv"), sep = "\t")
cat("\nStress-set overlap (orthogroups):\n"); print(ov)
cat(sprintf("\nWritten to %s\n", OUT))
