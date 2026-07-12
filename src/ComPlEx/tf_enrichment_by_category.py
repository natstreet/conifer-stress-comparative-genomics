#!/usr/bin/env python3
"""Are transcription factors over-represented in the conserved co-expression core?

Reproducible producer. Tests whether transcription factors (PlantTFDB annotation of
the Picea abies gene models) are enriched across the five co-expression conservation
categories, with Fisher's exact test (each category vs the rest, and conserved vs
not_coex).

Universe (built from the current pipeline inputs):
  * testable universe = spruce genes that (a) are expressed (SC or SD expression matrix)
    AND (b) belong to an ortholog group containing >= 1 Pinus sylvestris gene, per
    ComPlEx's own ortholog table doc/genes_ortholog_categories.tsv. This is the set of
    genes that COULD be scored as a cross-species co-expressolog (100% of the 24,989
    co-expressolog spruce genes are covered by this table).
  * coex_category per gene = best-pval co-expressolog pair in weighted_gene_pairs.tsv,
    with the Methods `conserved` threshold NegLog10CliqueSum >= 10 (matches cliques_step2
    C_SUM and integration_analysis.R).
  * not_coex = testable-universe genes with no co-expressolog in any tissue.

Inputs (run from the AbioticStressConifers directory):
  data/expression/SC_expression.txt, SD_expression.txt      (expression universe)
  doc/genes_ortholog_categories.tsv                          (ComPlEx ortholog groups)
  results/ComPlEx/RData/weighted_gene_pairs.tsv              (co-expressologs + tissue flags)
  data/annotation/Picab02_230926_at01_longest_no_TE_aa_TF_predictions.tsv  (TF annotation)
Output: results/integration/tf_enrichment_by_category.tsv
"""
import csv
import re
from collections import defaultdict
from scipy.stats import fisher_exact

CATS = ["conserved", "multi_tissue", "cold_specific", "drought_specific", "not_coex"]
C_SUM = 10  # Methods conserved threshold on NegLog10CliqueSum (cliques_step2 C_SUM)

EXPR = ["data/expression/SC_expression.txt", "data/expression/SD_expression.txt"]
ORTHO = "doc/genes_ortholog_categories.tsv"
WGP = "results/ComPlEx/RData/weighted_gene_pairs.tsv"
TF = "data/annotation/Picab02_230926_at01_longest_no_TE_aa_TF_predictions.tsv"
OUT = "results/integration/tf_enrichment_by_category.tsv"

_strip = lambda s: s.strip().strip('"')
base = lambda g: re.sub(r'\.mRNA\.\d+$', '', _strip(g))


def expression_universe():
    u = set()
    for fn in EXPR:
        with open(fn) as f:
            next(f)
            for line in f:
                u.add(base(line.split('\t')[0]))
    return u


def ortholog_to_pine():
    """spruce genes in an ortholog group that also contains >=1 pine gene."""
    og_has_pine, spruce_by_og = set(), defaultdict(list)
    with open(ORTHO) as f:
        for d in csv.DictReader(f, delimiter='\t'):
            og = d["Ortholog_Group"]
            if d["species"] == "Pinus_sylvestris":
                og_has_pine.add(og)
            elif d["species"] == "Picea_abies":
                spruce_by_og[og].append(base(d["gene"]))
    return {g for og in og_has_pine for g in spruce_by_og[og]}


def coex_category_by_gene():
    """best-pval co-expressolog pair per spruce gene -> category (threshold applied)."""
    best = {}
    with open(WGP) as f:
        for d in csv.DictReader(f, delimiter='\t'):
            g = base(d["Species1"])
            bp = float(d["best_pval"])
            cons = d["conserved"] in ("TRUE", "True", "1")
            s = float(d["NegLog10CliqueSum"])
            cs = d["cold_specific"] in ("TRUE", "True", "1")
            ds = d["drought_specific"] in ("TRUE", "True", "1")
            cat = ("conserved" if (cons and s >= C_SUM) else
                   "cold_specific" if cs else
                   "drought_specific" if ds else "multi_tissue")
            if g not in best or bp < best[g][0]:
                best[g] = (bp, cat)
    return {g: v[1] for g, v in best.items()}


def tf_set():
    tf = set()
    with open(TF) as f:
        next(f)
        for line in f:
            g = line.split('\t')[0].strip()
            if g:
                tf.add(base(g))
    return tf


def main():
    expr = expression_universe()
    orth = ortholog_to_pine()
    universe = expr & orth
    cat_of = coex_category_by_gene()
    tf = tf_set()
    cat_gene = {g: cat_of.get(g, "not_coex") for g in universe}

    tot = len(universe)
    tot_tf = sum(g in tf for g in universe)
    rows = []
    for c in CATS:
        genes = [g for g in universe if cat_gene[g] == c]
        a, n = sum(g in tf for g in genes), len(genes)
        b = tot_tf - a
        orr, p = fisher_exact([[a, n - a], [b, (tot - n) - b]])
        rows.append({"category": c, "n": n, "n_tf": a, "pct_tf": round(100 * a / n, 1),
                     "odds_ratio_vs_rest": round(orr, 2), "p_value": p})

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    cc = [g for g in universe if cat_gene[g] == "conserved"]
    nc = [g for g in universe if cat_gene[g] == "not_coex"]
    a, b = sum(g in tf for g in cc), sum(g in tf for g in nc)
    orr, p = fisher_exact([[a, len(cc) - a], [b, len(nc) - b]])

    print(f"testable universe: {tot} spruce genes (expressed & pine-ortholog); TF total {tot_tf}")
    for r in rows:
        print(f"  {r['category']:16} n={r['n']:5} TF={r['n_tf']:4} ({r['pct_tf']}%)  "
              f"OR_vs_rest={r['odds_ratio_vs_rest']}  P={r['p_value']:.2e}")
    print(f"\nconserved {100*a/len(cc):.1f}% vs not_coex {100*b/len(nc):.1f}% TF : "
          f"OR={orr:.2f}, P={p:.2e}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
