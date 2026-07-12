#!/usr/bin/env python3
"""build_category_universe.py — per-gene co-expression conservation universe.

Emits one gene-level table that the downstream category consumers read, so the 5-category
scheme is defined in a single place, from the current pipeline inputs.

Universe = spruce genes that (a) are expressed (SC or SD matrix) AND (b) belong to an
ortholog group containing >= 1 Pinus sylvestris gene (ComPlEx's own ortholog table). This
is the set that could be scored as a cross-species co-expressolog.

Per gene:
  coex_category        best-pval co-expressolog pair in weighted_gene_pairs.tsv, with the
                       Methods `conserved` threshold NegLog10CliqueSum >= 10 (cliques_step2
                       C_SUM); not_coex if the gene has no co-expressolog in any tissue.
  *_present            best pair's per-tissue presence flags (0 for not_coex genes).
  OrthoGroup           from weighted_gene_pairs (co-expressologs) or the ortholog table
                       (not_coex); same OrthoFinder namespace as Orthogroups_130323 / PREDEF.

Inputs (run from AbioticStressConifers/):
  data/expression/SC_expression.txt, SD_expression.txt
  doc/genes_ortholog_categories.tsv
  results/ComPlEx/RData/weighted_gene_pairs.tsv
Output: results/integration/coex_category_universe.tsv
"""
import csv
import re
from collections import defaultdict

C_SUM = 10
EXPR = ["data/expression/SC_expression.txt", "data/expression/SD_expression.txt"]
ORTHO = "doc/genes_ortholog_categories.tsv"
WGP = "results/ComPlEx/RData/weighted_gene_pairs.tsv"
OUT = "results/integration/coex_category_universe.tsv"
TISSUES = ["cold_needle", "cold_root", "drought_needle", "drought_root"]

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


def ortholog_tables():
    """spruce genes in a pine-containing OG, and each spruce gene's OrthoGroup."""
    og_has_pine, spruce_by_og, gene_og = set(), defaultdict(list), {}
    with open(ORTHO) as f:
        for d in csv.DictReader(f, delimiter='\t'):
            og = d["Ortholog_Group"]
            if d["species"] == "Pinus_sylvestris":
                og_has_pine.add(og)
            elif d["species"] == "Picea_abies":
                g = base(d["gene"])
                spruce_by_og[og].append(g)
                gene_og[g] = og
    ortho = {g for og in og_has_pine for g in spruce_by_og[og]}
    return ortho, gene_og


def best_coex():
    """best-pval co-expressolog pair per spruce gene -> (category, presence flags, OrthoGroup)."""
    best = {}
    with open(WGP) as f:
        for d in csv.DictReader(f, delimiter='\t'):
            g = base(d["Species1"])
            bp = float(d["best_pval"])
            if g in best and bp >= best[g][0]:
                continue
            cons = d["conserved"] in ("TRUE", "True", "1")
            s = float(d["NegLog10CliqueSum"])
            cs = d["cold_specific"] in ("TRUE", "True", "1")
            ds = d["drought_specific"] in ("TRUE", "True", "1")
            cat = ("conserved" if (cons and s >= C_SUM) else
                   "cold_specific" if cs else
                   "drought_specific" if ds else "multi_tissue")
            pres = {t: d[f"{t}_present"] for t in TISSUES}
            best[g] = (bp, cat, pres, d.get("OrthoGroup", ""))
    return best


def main():
    expr = expression_universe()
    ortho, gene_og = ortholog_tables()
    universe = expr & ortho
    best = best_coex()

    rows = []
    for g in sorted(universe):
        if g in best:
            _, cat, pres, og = best[g]
            row = {"pa_gene": g, "OrthoGroup": og or gene_og.get(g, ""),
                   "coex_category": cat}
            row.update({f"{t}_present": pres[t] for t in TISSUES})
        else:
            row = {"pa_gene": g, "OrthoGroup": gene_og.get(g, ""),
                   "coex_category": "not_coex"}
            row.update({f"{t}_present": "0" for t in TISSUES})
        rows.append(row)

    fields = ["pa_gene", "OrthoGroup", "coex_category"] + [f"{t}_present" for t in TISSUES]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    c = Counter(r["coex_category"] for r in rows)
    print(f"universe: {len(rows)} spruce genes (expressed & pine-ortholog)")
    for k in ["conserved", "multi_tissue", "cold_specific", "drought_specific", "not_coex"]:
        print(f"  {k:16} {c[k]}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
