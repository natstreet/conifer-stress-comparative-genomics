#!/usr/bin/env python3
"""Build results/integration/network_degree.tsv — the per-gene mean co-expression
network degree across the four Norway spruce ComPlEx stress networks.

Run from the AbioticStressConifers/ directory.

Inputs:
  results/integration/centrality/centrality_spruce_<net>.tsv  for
  <net> in {cold_needle, cold_root, drought_needle, drought_root}
  (columns: Genes, Degree; the ComPlEx Python-port per-network centrality output).

Output:
  results/integration/network_degree.tsv  (columns: Genes, mean_degree)

mean_degree is the mean of Degree across the networks in which the gene is
present (genes absent from a network do not contribute to their own mean). The
file is consumed by pnps_confound_analysis.py (the Figure 5a degree overlay and
the dN/dS confound control). This builder reproduces the deposited
network_degree.tsv value-for-value across all 28,857 genes.
"""
import pandas as pd

NETS = ["cold_needle", "cold_root", "drought_needle", "drought_root"]
CDIR = "results/integration/centrality"
OUT  = "results/integration/network_degree.tsv"

frames = []
for n in NETS:
    d = pd.read_csv(f"{CDIR}/centrality_spruce_{n}.tsv", sep="\t")
    d = d.rename(columns={d.columns[0]: "Genes", d.columns[1]: "Degree"})
    frames.append(d.set_index("Genes")["Degree"])

deg = pd.concat(frames, axis=1).mean(axis=1).rename("mean_degree")
deg = deg.sort_index().reset_index()
deg.columns = ["Genes", "mean_degree"]
deg.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT} : {len(deg)} genes")
