#!/usr/bin/env python3
"""single_copy_dnds_breadth.py — single-copy dN/dS-vs-conservation-breadth robustness check.

Tests whether the negative dN/dS vs conservation-breadth relationship holds when the analysis is
restricted to single-copy orthologues, under several explicit single-copy definitions, so the
choice of definition is transparent and each reported value is computed here rather than quoted.
The manuscript uses the N10 hierarchical-orthogroup single-copy set (hog_1to1_sp_pi): 11,166 pairs
(82.6% of the 13,519 with a usable dN/dS), Spearman rho = -0.166, P = 2.1e-69.

Method (mirrors figure5_dnds_evolution.R): conservation breadth per Picea gene =
max over its co-expressolog pairs of the number of the four stress-tissue contexts in
which the pair is present (weighted_gene_pairs.tsv); pairs with no co-expressolog get
breadth 0. dN/dS from cross_species_dnds_yn00.tsv (usable = numeric). Spearman via
scipy.stats.spearmanr (matches R cor.test spearman, exact=FALSE, for large n).

Single-copy definitions (reported for all; PRIMARY is 1:1 spruce:pine at OG level):
  og_1to1_sp_pi : orthogroup has exactly 1 Picea abies AND 1 Pinus sylvestris gene
  hog_1to1_sp_pi: N10 HOG has exactly 1 spruce AND 1 pine gene
  allspecies_sc : OrthoFinder single-copy orthogroup (1 gene in every species)

Inputs (deposit): OrthoFinder outputs under Orthogroups/ at the repo root (SOURCES.tsv).
Run from AbioticStressConifers/:  python3 src/ComPlEx/single_copy_dnds_breadth.py
Output: results/integration/single_copy_dnds_breadth.tsv
"""
import csv, gzip, os, re
from pathlib import Path
from collections import defaultdict
from scipy.stats import spearmanr

PROJ = Path(os.getcwd())
OF = PROJ / "Orthogroups"                                    # deposited OrthoFinder tables (SOURCES.tsv)
HOG10 = PROJ / "Orthogroups/Phylogenetic_Hierarchical_Orthogroups_N10.tsv.gz"
DNDS = PROJ / "results/integration/cross_species_dnds_yn00.tsv"
WGP = PROJ / "results/ComPlEx/RData/weighted_gene_pairs.tsv"
OUT = PROJ / "results/integration/single_copy_dnds_breadth.tsv"

base = lambda g: re.sub(r'\.mRNA\.\d+$', '', g.strip())


def breadth_by_gene():
    b = defaultdict(int)
    with open(WGP) as f:
        r = csv.DictReader(f, delimiter='\t')
        for d in r:
            n = sum(int(d[c]) for c in ("cold_needle_present", "cold_root_present",
                                        "drought_needle_present", "drought_root_present"))
            g = base(d["Species1"])
            if n > b[g]:
                b[g] = n
    return b


def usable_pairs():
    # 'Usable' uses the same saturation filter as figure5_dnds_evolution.R (dS in (0,5), dN/dS < 10),
    # so the usable set is defined identically across the paper; this drops dS-saturated pairs (dS=99).
    b = breadth_by_gene()
    rows = []
    with open(DNDS) as f:
        r = csv.DictReader(f, delimiter='\t')   # read by header name, robust to column order/count
        for d in r:
            try:
                dnds = float(d["dNdS"]); ds = float(d["dS"])
            except (ValueError, KeyError, TypeError):
                continue
            if not (0 < ds < 5 and dnds < 10):
                continue
            pa, ps = base(d["pa_gene"]), base(d["ps_gene"])
            rows.append((pa, ps, dnds, b.get(pa, 0)))
    return rows   # (pa_gene, ps_gene, dNdS, breadth)


def pipeline_strict_pairs():
    """Pairs the pipeline itself flags as strictly 1:1 HOG: pct_identity is populated
    only for those (NA for multi-copy pairs) in integration_backbone_1to1.tsv."""
    s = set()
    bb = PROJ / "results/integration/integration_backbone_1to1.tsv"
    with open(bb) as f:
        for d in csv.DictReader(f, delimiter='\t'):
            if d.get("pct_identity", "") not in ("", "NA"):
                s.add((base(d["pa_gene"]), base(d["ps_gene"])))
    return s


def gene_to_og():
    g2og = {}
    with open(OF / "Orthogroups.tsv") as f:
        h = f.readline().rstrip('\n').split('\t')
        si, pi = h.index("Picea_abies"), h.index("Pinus_sylvestris")
        for line in f:
            p = line.rstrip('\n').split('\t')
            for idx in (si, pi):
                if idx < len(p) and p[idx].strip():
                    for g in p[idx].split(','):
                        g2og[base(g)] = p[0]
    return g2og


def og_single_copy():
    sc = set()
    with open(OF / "Orthogroups.GeneCount.tsv") as f:
        r = csv.DictReader(f, delimiter='\t')
        for d in r:
            if int(d["Picea_abies"]) == 1 and int(d["Pinus_sylvestris"]) == 1:
                sc.add(d["Orthogroup"])
    return sc


def hog_single_copy_genes():
    """return set of Picea genes whose N10 HOG is 1 spruce + 1 pine."""
    sc_genes = set()
    with gzip.open(HOG10, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        si, pi = h.index("Picea_abies"), h.index("Pinus_sylvestris")
        for line in f:
            p = line.rstrip('\n').split('\t')
            sp = [base(x) for x in p[si].split(',') if x.strip()] if si < len(p) and p[si].strip() else []
            pn = [base(x) for x in p[pi].split(',') if x.strip()] if pi < len(p) and p[pi].strip() else []
            if len(sp) == 1 and len(pn) == 1:
                sc_genes.add(sp[0])
    return sc_genes


def allspecies_sc_genes(g2og):
    scogs = set(l.strip() for l in open(OF / "Orthogroups_SingleCopyOrthologues.txt") if l.strip())
    return {g for g, og in g2og.items() if og in scogs}


def spear(rows):
    x = [r[3] for r in rows]; y = [r[2] for r in rows]   # breadth, dNdS
    rho, p = spearmanr(x, y)
    return len(rows), rho, p


def main():
    rows = usable_pairs()
    g2og = gene_to_og()
    og_sc = og_single_copy()
    hog_sc = hog_single_copy_genes()
    allsc = allspecies_sc_genes(g2og)
    strict = pipeline_strict_pairs()

    # each predicate takes a row (pa, ps, dNdS, breadth)
    defs = {
        "full_usable": lambda r: True,
        "pipeline_strict_1to1": lambda r: (r[0], r[1]) in strict,   # pipeline's own flag
        "og_1to1_sp_pi": lambda r: g2og.get(r[0]) in og_sc,
        "hog_1to1_sp_pi": lambda r: r[0] in hog_sc,
        "allspecies_sc": lambda r: r[0] in allsc,
    }
    results = []
    for name, keep in defs.items():
        sub = [r for r in rows if keep(r)]
        n, rho, p = spear(sub)
        pct = 100 * n / len(rows)
        results.append((name, n, round(pct, 1), round(rho, 4), p))

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(["definition", "n", "pct_of_usable", "spearman_rho", "spearman_p"])
        for r in results:
            w.writerow(r)

    print(f"usable dN/dS pairs: {len(rows)}\n")
    print(f"{'definition':<16} {'n':>6} {'%':>6} {'rho':>8} {'p':>12}")
    for name, n, pct, rho, p in results:
        print(f"{name:<16} {n:>6} {pct:>6} {rho:>8} {p:>12.2e}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
