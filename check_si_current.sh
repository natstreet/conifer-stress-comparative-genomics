#!/usr/bin/env bash
# check_si_current.sh — verify that a deposited supplementary-information set matches
# what the committed producers currently generate. Regenerates the SI tables from
# their producers and content-compares them (numeric-tolerant, formatting-agnostic)
# against the copies in TARGET_DIR. Exits non-zero if any table is stale.
#
# Run this before every submission or deposit. It catches the failure mode where a
# producer or an upstream result table was updated but the deposited SI file was not
# rebuilt from it.
#
# Usage:
#   ./check_si_current.sh TARGET_DIR
#
# Note: this checks the tables (pure formatters of committed result tables, where the
# staleness risk lives). Rebuild and visually check Figures S1-S6 with
# ./build_submission_si.sh — figures need an image comparison, not a cell comparison.
#
# Environment: R (default Rscript).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
R="${R:-Rscript}"
TARGET="${1:?usage: ./check_si_current.sh TARGET_DIR}"
TABLES_DIR="manuscript/supplementary"

echo "Regenerating SI tables from producers..."
"$R" src/tables/table_s1_go_enrichment.R              >/dev/null 2>&1
"$R" src/tables/table_s2_de_results.R                 >/dev/null 2>&1
"$R" src/figures/figureS2_timing_deg_fraction.R       >/dev/null 2>&1   # writes Table S3
"$R" src/tables/table_s4_conserved_coexpressolog_go.R >/dev/null 2>&1
"$R" src/tables/table_s5_coexpression_category_go.R   >/dev/null 2>&1
"$R" src/tables/table_s6_dnds_yn00.R                  >/dev/null 2>&1
"$R" src/tables/table_s7_threshold_sensitivity.R      >/dev/null 2>&1
"$R" src/tables/table_s8_gene_axis_classification.R   >/dev/null 2>&1

echo "Comparing regenerated tables against: $TARGET"
# Copy the target tables to a clean temp dir first: some deposit paths (e.g. OneDrive mounts
# with non-ASCII characters) cannot be opened directly by R's xlsx reader, but a plain shell
# copy handles them. The temp dir is removed on exit.
TMP_TGT="$(mktemp -d)"; trap 'rm -rf "$TMP_TGT"' EXIT
cp -f "$TARGET"/Supplementary_Table_S*.xlsx "$TMP_TGT"/ 2>/dev/null || true
"$R" - "$TABLES_DIR" "$TMP_TGT" <<'RS'
suppressMessages(library(readxl)); options(warn = -1)
a <- commandArgs(trailingOnly = TRUE); regen_dir <- a[1]; target_dir <- a[2]
rd <- function(p) tryCatch(as.matrix(as.data.frame(
  read_excel(p, sheet = 1, col_names = FALSE, .name_repair = "minimal"))), error = function(e) NULL)
# numeric-tolerant, formatting-agnostic count of differing cells (Inf = dimensions differ)
tol_diff <- function(A, B) {
  if (is.null(A) || is.null(B)) return(NA)
  if (!identical(dim(A), dim(B))) return(Inf)
  n <- 0
  for (i in seq_len(nrow(A))) for (j in seq_len(ncol(A))) {
    na <- suppressWarnings(as.numeric(A[i, j])); nb <- suppressWarnings(as.numeric(B[i, j]))
    d <- if (!is.na(na) && !is.na(nb)) abs(na - nb) > 1e-6 * max(1, abs(na), abs(nb))
         else !identical(ifelse(is.na(A[i, j]), "", A[i, j]), ifelse(is.na(B[i, j]), "", B[i, j]))
    if (isTRUE(d)) n <- n + 1
  }
  n
}
files <- sort(list.files(regen_dir, pattern = "^Supplementary_Table_S.*\\.xlsx$"))
stale <- 0
for (f in files) {
  r <- rd(file.path(regen_dir, f)); s <- rd(file.path(target_dir, f))
  if (is.null(s)) { cat(sprintf("  %-48s MISSING in target\n", f)); stale <- stale + 1; next }
  d <- tol_diff(r, s)
  if (is.na(d))              cat(sprintf("  %-48s read-fail\n", f))
  else if (d == 0)           cat(sprintf("  %-48s current\n", f))
  else { cat(sprintf("  %-48s STALE (%s)\n", f, if (is.infinite(d)) "dimensions differ" else paste(d, "differing cells")))
         stale <- stale + 1 }
}
cat(sprintf("\n%d of %d SI tables are STALE.\n", stale, length(files)))
quit(status = if (stale > 0) 1 else 0)
RS
