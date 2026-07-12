#!/usr/bin/env python3
"""Threshold sensitivity analysis for the cross-species co-expression categories.

Re-derives the five co-expression conservation categories from the weighted
co-expressolog gene pairs across a range of evidence thresholds (C_SUM for the
conserved category; C_SUM_STRESS for the stress-specific and multi-tissue
categories), and reports — for each configuration — the category sizes, the
median dN/dS (YN00) per category, the conserved-vs-not_coex Wilcoxon test, and
the shared segmental-duplicate enrichment in the drought-specific and not_coex
categories.

Baseline (structural; C_SUM/C_SUM_STRESS not binding) reproduces the published
category sizes 1890/5096/7687/10316/18062.

Output: results/integration/threshold_sensitivity.tsv
Run from the AbioticStressConifers directory.
"""
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, fisher_exact

WP   = "results/ComPlEx/RData/weighted_gene_pairs.tsv"
SC   = "data/expression/SC_expression.txt"
SD   = "data/expression/SD_expression.txt"
YN   = "results/integration/cross_species_dnds_yn00.tsv"
KAKS = "kaks_results.tsv"                       # in the project root (../ from here)
OUT  = "results/integration/threshold_sensitivity.tsv"

import os
KAKS_PATH = KAKS if os.path.exists(KAKS) else os.path.join("..", KAKS)

wp = pd.read_csv(WP, sep="\t").rename(columns={"Species1": "pa_gene", "Species2": "ps_gene"})
universe = set(pd.read_csv(SC, sep="\t", usecols=[0]).iloc[:, 0]) | \
           set(pd.read_csv(SD, sep="\t", usecols=[0]).iloc[:, 0])
yn = pd.read_csv(YN, sep="\t")
yf = yn[(yn.dS > 0) & (yn.dS < 5) & (yn.dNdS < 10)].dropna(subset=["dNdS"])  # manuscript filter (Fig 5a: dS in (0,5), dNdS<10)
dnds = yf.drop_duplicates("pa_gene").set_index("pa_gene")["dNdS"]
kk = pd.read_csv(KAKS_PATH, sep="\t")
shared_sd = set(kk[kk.category == "shared_SD"].gene1) | set(kk[kk.category == "shared_SD"].gene2)


def classify(c_sum=0.0, c_sum_stress=0.0):
    w = wp.copy()
    cons  = w.conserved        & (w.NegLog10CliqueSum >= c_sum)
    cold  = w.cold_specific    & (w.ColdSum            >= c_sum_stress)
    dro   = w.drought_specific & (w.DroughtSum         >= c_sum_stress)
    multi = (w.NegLog10CliqueSum >= c_sum_stress)
    cat = np.where(cons, "conserved",
          np.where(cold, "cold_specific",
          np.where(dro, "drought_specific",
          np.where(multi, "multi_tissue", "none"))))
    w = w.assign(cat=cat)
    w = w[w.cat != "none"].sort_values("best_pval").drop_duplicates("pa_gene")
    return dict(zip(w.pa_gene, w.cat))


def summarise(g2c, label):
    catof = lambda g: g2c.get(g, "not_coex")
    cats = ["conserved", "multi_tissue", "cold_specific", "drought_specific", "not_coex"]
    med = {c: float(np.median([dnds[g] for g in dnds.index if catof(g) == c])) for c in cats}
    vC = [dnds[g] for g in dnds.index if catof(g) == "conserved"]
    vN = [dnds[g] for g in dnds.index if catof(g) == "not_coex"]
    W, p = mannwhitneyu(vC, vN, alternative="two-sided")
    from collections import Counter
    sz = Counter(catof(g) for g in universe)

    def OR(tc):
        a = sum(catof(g) == tc and g in shared_sd for g in universe)
        b = sum(catof(g) == tc and g not in shared_sd for g in universe)
        c = sum(catof(g) != tc and g in shared_sd for g in universe)
        d = sum(catof(g) != tc and g not in shared_sd for g in universe)
        return fisher_exact([[a, b], [c, d]])[0]

    ends = med["conserved"] == min(med.values()) and med["not_coex"] == max(med.values())
    return [label, sz["conserved"], sz["cold_specific"], sz["drought_specific"], sz["multi_tissue"],
            sz["not_coex"], round(med["conserved"], 3), round(med["multi_tissue"], 3),
            round(med["cold_specific"], 3), round(med["drought_specific"], 3), round(med["not_coex"], 3),
            ends, int(W), f"{p:.1e}", round(OR("drought_specific"), 2), round(OR("not_coex"), 2)]


rows = [summarise(classify(0, 0), "Baseline (structural; as published)")]
for cs, css in [(10, 5), (15, 8), (20, 10), (30, 15)]:
    rows.append(summarise(classify(cs, css), f"C_SUM={cs}, C_SUM_STRESS={css}"))

cols = ["config", "n_conserved", "n_cold", "n_drought", "n_multi", "n_not_coex",
        "med_conserved", "med_multi", "med_cold", "med_drought", "med_not_coex",
        "conserved_lowest_and_not_coex_highest", "W_cons_vs_not_coex", "p_cons_vs_not_coex",
        "sharedSD_OR_drought_specific", "sharedSD_OR_not_coex"]
df = pd.DataFrame(rows, columns=cols)
df.to_csv(OUT, sep="\t", index=False)
print(df.to_string(index=False))
