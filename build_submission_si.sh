#!/usr/bin/env bash
# build_submission_si.sh — regenerate every supplementary-information deliverable
# (Supplementary Tables S1-S8 and Supplementary Figures S1-S6) from the committed
# producers, and optionally copy the complete set to a target directory.
#
# The supplementary files are build artifacts, not hand-maintained documents:
# regenerate them here so the deposited SI always matches the current analysis,
# rather than editing or copying individual files by hand. The table producers
# format the result tables under results/; run ./reproduce_paper.sh first if those
# upstream results need regenerating from the raw data.
#
# Usage:
#   ./build_submission_si.sh [TARGET_DIR]      # regenerate, then copy the SI set to TARGET_DIR
#   ./build_submission_si.sh                   # regenerate in place only
#   SKIP_SLOW=1 ./build_submission_si.sh ...    # skip the dN/dS robustness re-run (Figure S6)
#
# Environment: R (default Rscript), PY (default python3).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
R="${R:-Rscript}"; PY="${PY:-python3}"
TARGET="${1:-}"
TABLES_DIR="manuscript/supplementary"
FIG_DIR="results/final_figures"
fail=0
run(){ printf '  %-46s' "$1"; shift; if "$@" >/dev/null 2>&1; then echo "ok"; else echo "FAILED"; fail=1; fi; }

echo "== Supplementary tables (S1-S8) =="
run "Table S1 (GO enrichment)"              "$R" src/tables/table_s1_go_enrichment.R
run "Table S2 (DE results)"                 "$R" src/tables/table_s2_de_results.R
run "Figure S2 + Table S3 (timing)"         "$R" src/figures/figureS2_timing_deg_fraction.R
run "Table S4 (conserved coexpressolog GO)" "$R" src/tables/table_s4_conserved_coexpressolog_go.R
run "Table S5 (coexpression category GO)"   "$R" src/tables/table_s5_coexpression_category_go.R
run "Table S6 (dN/dS YN00)"                 "$R" src/tables/table_s6_dnds_yn00.R
run "Table S7 (threshold sensitivity)"      "$R" src/tables/table_s7_threshold_sensitivity.R
run "Table S8 (gene axis classification)"   "$R" src/tables/table_s8_gene_axis_classification.R

echo "== Supplementary figures (S1-S6) =="
run "Figure S1 (physiology)"                "$R" src/figures/figureS1_physiology.R
run "Figure S5 (integrative model)"         "$R" src/figures/figureS5_integrative_model.R
if [ "${SKIP_SLOW:-0}" = 1 ]; then
  echo "  Figure S6 (dN/dS robustness)                   skipped (SKIP_SLOW=1)"
else
  run "Figure S6 (dN/dS robustness; slow)"    "$PY" src/ComPlEx/dnds_robustness_subsample.py
fi
run "Figures S3, S4 (assemble composites)"  "$R" src/figures/assemble_figures.R

if [ -n "$TARGET" ]; then
  echo "== Copy SI to: $TARGET =="
  mkdir -p "$TARGET"
  cp -f "$TABLES_DIR"/Supplementary_Table_S*.xlsx "$TARGET"/
  for f in FigureS1 FigureS2 FigureS3 FigureS4 FigureS5 FigureS6; do
    cp -f "$FIG_DIR/$f.pdf" "$FIG_DIR/$f.png" "$TARGET"/ 2>/dev/null || true
  done
  echo "  copied the 8 tables + 6 figures (pdf+png)"
fi

[ "$fail" -eq 0 ] && echo "SI build complete." || { echo "SI build finished WITH FAILURES (see above)."; exit 1; }
