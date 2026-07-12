#!/usr/bin/env python3
"""
Expression divergence between the two copies of each Norway spruce
segmental-duplicate (SD) gene pair, in each stress-tissue context.

For every SD pair and context, divergence is 1 - Pearson correlation of the two
copies' expression profiles across that context's samples; pairs whose mean
expression is below MIN_EXPR in either copy are left as NA. These per-pair,
per-context divergence values feed the promoter TE-Jaccard vs expression-
divergence analysis (Supplementary Figure S3).

Output:
  sd_pair_expression_divergence.tsv  -- hog_id, gene1, gene2, category, species,
                                         div_{cold_needles, cold_roots,
                                              drought_needles, drought_roots}
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJ   = Path(__file__).resolve().parents[2]   # AbioticStressConifers/
OUTDIR = PROJ                                  # committed intermediates are written here
# Spruce per-context stress expression matrices (variance-stabilised counts, one file per
# cold/drought x needle/root condition), deposited under data/expression/.
EXPR_DIR = PROJ / "data/expression"

MIN_EXPR = 0.5
CONTEXTS = ["cold_needles", "cold_roots", "drought_needles", "drought_roots"]


def load_expr(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df


def expr_divergence(gene1, gene2, ematrix):
    try:
        v1 = ematrix.loc[gene1].values.astype(float)
        v2 = ematrix.loc[gene2].values.astype(float)
        if np.nanmean(v1) < MIN_EXPR or np.nanmean(v2) < MIN_EXPR:
            return np.nan
        r, _ = stats.pearsonr(v1, v2)
        return float(1 - r)
    except Exception:
        return np.nan


# ── Spruce SD pairs (from kaks_results) ───────────────────────────────────────
kaks = pd.read_csv(OUTDIR / "kaks_results.tsv", sep="\t")
all_pairs = kaks[["hog_id", "gene1", "gene2", "category"]].copy()
all_pairs["species"] = "PA"
print(f"SD pairs: {len(all_pairs)} spruce", flush=True)

# ── Expression matrices per stress context ────────────────────────────────────
print("Loading expression data…", flush=True)
expr = {
    "PA": {
        "cold_needles":    load_expr(EXPR_DIR / "SCN_expression.txt"),
        "cold_roots":      load_expr(EXPR_DIR / "SCR_expression.txt"),
        "drought_needles": load_expr(EXPR_DIR / "SDN_expression.txt"),
        "drought_roots":   load_expr(EXPR_DIR / "SDR_expression.txt"),
    },
}

# ── Per-pair expression divergence ────────────────────────────────────────────
rows = []
for _, row in all_pairs.iterrows():
    sp_expr = expr.get(row["species"], {})
    rec = {"hog_id": row["hog_id"], "gene1": row["gene1"], "gene2": row["gene2"],
           "category": row["category"], "species": row["species"]}
    for ctx in CONTEXTS:
        em = sp_expr.get(ctx)
        rec[f"div_{ctx}"] = expr_divergence(row["gene1"], row["gene2"], em) if em is not None else np.nan
    rows.append(rec)
    if len(rows) % 500 == 0:
        print(f"  {len(rows)}/{len(all_pairs)} pairs…", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUTDIR / "sd_pair_expression_divergence.tsv", sep="\t", index=False)
print(f"Saved sd_pair_expression_divergence.tsv ({len(df)} pairs)", flush=True)
