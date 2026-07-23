#!/usr/bin/env python3
"""pine_sd_classify.py — classify N10 orthogroups by the LINEAGE of segmental duplication and emit the
PINE segmental-duplicate gene classes, so the SD conservation axis can be tested symmetrically for pine as
well as spruce (co-author review). Reuses the committed HOG / gene-set loaders and the EXACT per-HOG gene
selection from kaks_within_species_sd.py (the Ka/Ks representative pair: the first two SD genes of the HOG,
falling back to the first two genes when the HOG has <2 SD genes), so the shared_SD / spruce_only_SD gene
sets it reproduces are identical to the existing spruce analysis (built-in anchor: 1,874 / 3,186 in the
spruce expression universe). It then applies the same selection to the pine side and adds pine_only_SD.

  shared_SD      : spruce SD gene present AND pine SD gene present AND >=2 spruce genes in the HOG  (rule unchanged)
  spruce_only_SD : >=2 spruce genes in the spruce SD set AND no pine SD gene                        (rule unchanged)
  pine_only_SD   : >=2 pine   genes in the pine   SD set AND no spruce SD gene   (NEW — exact mirror of spruce_only)

Genes are tagged on each species' OWN axis via the same representative-pair rule the Ka/Ks analysis used, so
pine_only_SD is directly comparable to spruce_only_SD. Output:
  results/integration/pine_sd_gene_class.tsv   (ps_gene, sd_class in {shared_SD, pine_only_SD})
"""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kaks_within_species_sd import load_gene_set, load_hog_table   # committed, proven loaders

DEP = os.environ.get("SPRUCE_PINE_DEPOSIT", os.getcwd())
HOG   = f"{DEP}/Orthogroups/Phylogenetic_Hierarchical_Orthogroups_N10.tsv.gz"
SP_SD = f"{DEP}/Orthogroups/spruce_all_in_duplicated_blocks_PGset.tsv.gz"
PI_SD = f"{DEP}/Orthogroups/pine_all_in_duplicated_blocks_PGset.tsv.gz"
OUT   = "results/integration/pine_sd_gene_class.tsv"

def pair_exact(genes_ordered, sd_set):
    """The Ka/Ks representative pair (kaks_within_species_sd.py): first two SD genes, else first two genes."""
    ins = [g for g in genes_ordered if g in sd_set]
    if len(ins) >= 2:            return ins[:2]
    if len(genes_ordered) >= 2:  return genes_ordered[:2]
    return []

def expr_universe(*files):
    u = set()
    for f in files:
        r = csv.reader(open(f), delimiter="\t"); next(r)
        u |= {row[0] for row in r if row}
    return u

spruce_sd = load_gene_set(SP_SD)
pine_sd   = load_gene_set(PI_SD)
hogs      = load_hog_table(HOG)
SU = expr_universe("data/expression/SC_expression.txt", "data/expression/SD_expression.txt")
PU = expr_universe("data/expression/PC_expression.txt", "data/expression/PD_expression.txt")

sp_shared = set(); sp_only = set()                 # spruce anchor (must reproduce 1,874 / 3,186 in universe)
pine_class = {}                                    # ps_gene -> sd_class
for d in hogs.values():
    sg = d['spruce']; pg = d['pine']
    sg_in = set(g for g in sg if g in spruce_sd)
    pg_in = set(g for g in pg if g in pine_sd)
    if sg_in and pg_in and len(sg) >= 2:                       # shared_SD (rule unchanged)
        sp_shared |= set(pair_exact(sg, spruce_sd))
        for g in pair_exact(pg, pine_sd): pine_class[g] = "shared_SD"
    elif len(sg_in) >= 2 and not pg_in:                        # spruce_only_SD (unchanged; anchor only)
        sp_only |= set(pair_exact(sg, spruce_sd))
    elif len(pg_in) >= 2 and not sg_in:                        # pine_only_SD (NEW — mirror of spruce_only)
        for g in pair_exact(pg, pine_sd): pine_class.setdefault(g, "pine_only_SD")

os.makedirs("results/integration", exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("ps_gene\tsd_class\n")
    for g in sorted(pine_class):
        fh.write(f"{g}\t{pine_class[g]}\n")

# ── built-in anchor gate ──
sp_sh_u, sp_on_u = len(sp_shared & SU), len(sp_only & SU)
ok = (sp_sh_u == 1874 and sp_on_u == 3186)
print(f"  SPRUCE ANCHOR (in expr universe): shared_SD={sp_sh_u}  spruce_only_SD={sp_on_u}   "
      f"{'OK (matches committed 1,874 / 3,186)' if ok else 'MISMATCH — classification logic drifted!'}")
n_sh = sum(v == 'shared_SD' for v in pine_class.values())
n_po = sum(v == 'pine_only_SD' for v in pine_class.values())
print(f"  PINE classes (all):               shared_SD={n_sh}  pine_only_SD={n_po}")
print(f"  PINE classes (in pine universe):  shared_SD={len(set(g for g,v in pine_class.items() if v=='shared_SD') & PU)}"
      f"  pine_only_SD={len(set(g for g,v in pine_class.items() if v=='pine_only_SD') & PU)}")
print(f"  wrote {OUT}")
if not ok:
    sys.exit("ABORT: spruce anchor drifted; do not trust the pine classes.")
