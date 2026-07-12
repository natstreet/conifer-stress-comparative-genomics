#!/usr/bin/env python3
"""pN/pS corroboration and confound control for the co-expression-conservation
to coding-sequence-constraint gradient (Figure 5a).

The dN/dS gradient (conserved co-expressologs least diverged, not_coex most)
could arise trivially: broadly conserved genes are also more central in the
co-expression network and more highly expressed, and both predict stronger
purifying selection (Mahler et al. 2017). Two checks address this.

  1. pN/pS by category. Within-species pN/pS (Norway spruce population
     resequencing) is summarised per category as an independent, polymorphism-
     based measure of selection. It reproduces the gradient (median pN/pS 0.228
     in conserved rising to 0.282 in not_coex; Kruskal-Wallis P = 2.3e-13;
     conserved vs not_coex P = 7.0e-13), so the trend is not an artefact of
     between-species dN/dS estimation.

  2. Confound control. A partial Spearman correlation of dN/dS with
     conservation breadth (the per-gene count of stress-tissue comparisons with a
     conserved co-expressolog, 0-4), controlling for mean expression and network
     degree, remains significant (raw rho -0.14; expression-controlled -0.10;
     expression+degree -0.094; all P < 1e-24), so the coupling is robust to these
     correlated properties. Conserved-category genes have ~2x the degree of
     not_coex (median 385 vs 189).

Outputs (results/integration/): pnps_by_category.tsv, confound_control.tsv,
category_degree_expr_medians.tsv

Run from the AbioticStressConifers directory. pnps_sd_genes.tsv is vendored to
data/popgen/ (see SOURCES.tsv); per-gene network degree is read from the deposited
network_degree.tsv (mean across the four spruce ComPlEx networks).
"""
import os, glob, re
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr, rankdata, pearsonr
import statsmodels.api as sm

INTEG = "results/integration"
EXPR  = "data/expression"
# Norway spruce population-resequencing pN/pS (Kalman et al.), vendored to the deposit.
PNPS  = os.path.join("data/popgen", "pnps_sd_genes.tsv")
NETS  = ["cold_needle", "cold_root", "drought_needle", "drought_root"]
CATS  = ["conserved", "multi_tissue", "cold_specific", "drought_specific", "not_coex"]
PRESENT = ["cold_needle_present", "cold_root_present",
           "drought_needle_present", "drought_root_present"]

# Per-gene median pN/pS, and the gene set per co-expression category (the 1:1
# orthologue set for which dN/dS could be estimated).
pnps = pd.read_csv(PNPS, sep="\t")[["Gene", "median_pNpS"]].dropna()
pnps = dict(zip(pnps.Gene, pnps.median_pNpS))
cat_genes = {}
for f in glob.glob(os.path.join(INTEG, "category_genes_*.tsv")):
    name = re.sub(r"category_genes_|\.tsv", "", os.path.basename(f))
    cat_genes[name] = set(pd.read_csv(f, sep="\t").iloc[:, 0])

# 1. pN/pS by category --------------------------------------------------------
vals = {c: [pnps[g] for g in cat_genes.get(c, ()) if g in pnps] for c in CATS}
pnps_tab = pd.DataFrame([{"category": c, "n": len(vals[c]),
                          "median_pNpS": round(float(np.median(vals[c])), 4)} for c in CATS])
kw = kruskal(*[vals[c] for c in CATS])
cn = mannwhitneyu(vals["conserved"], vals["not_coex"], alternative="two-sided")
pnps_tab.to_csv(os.path.join(INTEG, "pnps_by_category.tsv"), sep="\t", index=False)

# Emit the two in-text pN/pS-by-category test statistics to a committed cell (previously computed
# inline but never written): the 5-category Kruskal-Wallis and the conserved-vs-not_coex pairwise test.
pd.DataFrame([{"test": "kruskal_wallis_5cat",              "statistic": kw.statistic, "p_value": kw.pvalue},
              {"test": "mannwhitney_conserved_vs_not_coex", "statistic": cn.statistic, "p_value": cn.pvalue}]
             ).to_csv(os.path.join(INTEG, "pnps_by_category_stats.tsv"), sep="\t", index=False)

# 2. Confound control: dN/dS ~ conservation breadth | degree + expression -----
# Use the categories from the backbone as the reference (the dN/dS file's own coex_category
# column can fall out of step); join on the full (pa_gene, ps_gene) pair, as the
# Figure 5a script does, so degree/confound use the same gene set as the figure.
k = pd.read_csv(os.path.join(INTEG, "cross_species_dnds_yn00.tsv"), sep="\t")[["pa_gene", "ps_gene", "dN", "dS", "dNdS"]]
bb_cat = pd.read_csv(os.path.join(INTEG, "integration_backbone_1to1.tsv"),
                     sep="\t")[["pa_gene", "ps_gene", "coex_category"]]
k = bb_cat.merge(k, on=["pa_gene", "ps_gene"])
# usable dN/dS (same filter as Figure 5a: reliable dS and dN/dS estimates)
k = k[(k.dS > 0) & (k.dS < 5) & (k.dNdS < 10) & k.dNdS.notna() & k.coex_category.isin(CATS)].copy()
# mean expression across the two spruce stress matrices (proxy for expression level)
expr = pd.concat([pd.read_csv(os.path.join(EXPR, f), sep="\t", index_col=0)
                  for f in ("SC_expression.txt", "SD_expression.txt")], axis=1).mean(axis=1)
# mean co-expression network degree across the four spruce networks; precomputed
# from the ComPlEx centrality output and deposited so this reproduces from figshare.
deg = pd.read_csv(os.path.join(INTEG, "network_degree.tsv"), sep="\t").set_index("Genes")["mean_degree"]
k["expr"] = k.pa_gene.map(expr)
k["deg"]  = k.pa_gene.map(deg)
# Conservation breadth = number of the four stress-tissue comparisons in which the
# gene's orthogroup has a significant conserved co-expressolog (0-4), taken as the
# maximum over the gene's co-expressolog pairs; genes with none have breadth 0. This
# is the actual per-gene count, not a category-derived ordinal.
wp = pd.read_csv("results/ComPlEx/RData/weighted_gene_pairs.tsv", sep="\t")
wp["breadth"] = wp[PRESENT].sum(axis=1)
gene_breadth = wp.groupby("Species1")["breadth"].max()
k["breadth"] = k.pa_gene.map(gene_breadth).fillna(0).astype(int)
k = k.dropna(subset=["expr", "deg"])
k["logexpr"], k["logdeg"] = np.log1p(k.expr), np.log1p(k.deg)

def partial_spearman(y, x, covars):
    """Spearman partial correlation of x and y, controlling for covars,
    computed as the Pearson correlation of the rank residuals."""
    Z = sm.add_constant(np.column_stack([rankdata(c) for c in covars]))
    ry = sm.OLS(rankdata(y), Z).fit().resid
    rx = sm.OLS(rankdata(x), Z).fit().resid
    return pearsonr(rx, ry)

cc = pd.DataFrame(
    [("none",) + spearmanr(k.breadth, k.dNdS),
     ("expression",) + partial_spearman(k.dNdS, k.breadth, [k.logexpr]),
     ("degree",) + partial_spearman(k.dNdS, k.breadth, [k.logdeg]),
     ("expression+degree",) + partial_spearman(k.dNdS, k.breadth, [k.logexpr, k.logdeg])],
    columns=["control", "rho", "p"])
cc["rho"] = cc["rho"].round(3)
cc.to_csv(os.path.join(INTEG, "confound_control.tsv"), sep="\t", index=False)
med = k.groupby("coex_category").agg(
    dNdS=("dNdS", "median"), expr=("expr", "median"), deg=("deg", "median"),
    deg_q25=("deg", lambda s: s.quantile(0.25)),
    deg_q75=("deg", lambda s: s.quantile(0.75))).reindex(CATS).round(3)
med.to_csv(os.path.join(INTEG, "category_degree_expr_medians.tsv"), sep="\t")

print(pnps_tab.to_string(index=False))
print("Kruskal-Wallis P=%.2e ; conserved vs not_coex P=%.2e" % (kw.pvalue, cn.pvalue))
print(med.to_string())
print(cc.to_string(index=False))
