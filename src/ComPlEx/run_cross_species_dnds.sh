#!/usr/bin/env bash
# Committed driver for cross_species_dnds_yn00.tsv (the cross-species dN/dS backbone
# feeding Fig 5, pN/pS, Fig 7, threshold_sensitivity, single_copy_dnds_breadth,
# not_coex_de, wood_clique). Records the exact kaks_yn00.py invocation so the file is
# turnkey-reproducible from the committed pairs + raw genome CDS.
#
# Raw CDS inputs (FigShare doi:10.17044/scilifelab.28737623) — override via env:
#   DATA_ROOT  dir holding sprucev2/ + pinev1/ CDS FASTAs (default: ./genome_data)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PROJ="$(cd "$HERE/../.." && pwd)"
DATA="${DATA_ROOT:-$PROJ/genome_data}"
PY="${PYTHON:-python3}"

"$PY" "$HERE/kaks_yn00.py" \
  --pairs      "$PROJ/results/integration/integration_backbone_1to1.tsv" \
  --spruce-cds "$DATA/sprucev2/Picab02_230926_at01_longest_no_TE_cds.fa" \
  --pine-cds   "$DATA/pinev1/Pinsy01_240308_at01_longest_no_TE_cds.fa" \
  --out        "${OUT:-$PROJ/results/integration/cross_species_dnds_yn00.tsv}" \
  --workers    "${WORKERS:-4}"
