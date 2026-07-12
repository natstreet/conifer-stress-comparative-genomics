#!/bin/bash
# run_complex_networks.sh — rebuild the four ComPlEx co-expressolog networks from the VST
# expression matrices using the validated Python port (complex_py.py). The networks are a
# CORE ANALYSIS STAGE, not a deposited input: this script regenerates the per-network
# comparison tables + centrality that cliques_step1/1b consume. Run on every reproduction.
#
# ComPlEx_python location: $COMPLEX_PY_REPO (default: sibling ../ComPlEx_python of the repo root).
set -euo pipefail
HERE="$(cd "$(dirname "$0")/../.." && pwd)"
CPY="${COMPLEX_PY_REPO:-$HERE/../ComPlEx_python}"
if [ ! -f "$CPY/complex_py.py" ]; then
  echo "ERROR: complex_py.py not found at '$CPY'. Set COMPLEX_PY_REPO to the ComPlEx_python checkout" >&2
  echo "       (https://github.com/natstreet/ComPlEx_python)." >&2
  exit 1
fi
DENSITY="${COMPLEX_DENSITY:-0.03}"
# net  spruce-condition  pine-condition
NETS=("cold_needle SCN PCN" "cold_root SCR PCR" "drought_needle SDN PDN" "drought_root SDR PDR")
mkdir -p "$HERE/results/integration/centrality"
for spec in "${NETS[@]}"; do
  set -- $spec; net=$1; s1=$2; s2=$3
  echo "  [ComPlEx] $net  ($s1 vs $s2)  density=$DENSITY"
  python3 "$CPY/complex_py.py" \
    --s1-expr "$HERE/data/expression/${s1}_expression.txt" \
    --s2-expr "$HERE/data/expression/${s2}_expression.txt" \
    --orthologs "$HERE/doc/genes_ortholog_categories.tsv" \
    --s1-name spruce --s2-name pine \
    --out-dir "$HERE/results/ComPlEx/$net" --density "$DENSITY"
  # mirror per-network centrality to the canonical integration path used by
  # build_network_degree.py (spruce) and pine_axis_replication.R (pine)
  cp "$HERE/results/ComPlEx/$net/RData/centrality/centrality_spruce.tsv" \
     "$HERE/results/integration/centrality/centrality_spruce_${net}.tsv"
  cp "$HERE/results/ComPlEx/$net/RData/centrality/centrality_pine.tsv" \
     "$HERE/results/integration/centrality/centrality_pine_${net}.tsv"
done
echo "  [ComPlEx] all four networks rebuilt from VST matrices"
