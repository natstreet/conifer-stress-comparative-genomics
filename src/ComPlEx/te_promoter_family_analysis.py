#!/usr/bin/env python3
"""
Two analyses:

A. Within-species: which TE families/classes are enriched in promoters of
   SD gene pairs with diverged co-expression (both_not_coex) vs co-expressed
   (both_coex)?

B. Between-species (1:1 orthologs): does Jaccard similarity of PROMOTER TE
   family composition between orthologous PA/PS genes predict co-expression
   conservation category?

Inputs (all already computed):
  - Picab02_repeats.gff3.gz    (RepeatMasker GFF with family-level Motif)
  - Pinsy01_chromosomes_and_unplaced.all_repeats.gff.gz
  - PA/PS TSS from the representative (longest-CDS, no-TE) gene models
  - integration_backbone_1to1.tsv (cross-species pairs + coex categories)
"""

import bisect, gzip, os, re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr, kruskal
from statsmodels.stats.multitest import multipletests

# Project paths resolved relative to this script (run from anywhere)
PROJ     = Path(__file__).resolve().parents[2]   # AbioticStressConifers/
INTEG    = PROJ / "results/integration"
OUTDIR   = PROJ.parent                           # repository root
# The large genome and repeat annotations are too large to redistribute with the code.
# Set DATA_ROOT to a directory holding them (with sprucev2/ and pinev1/ subdirectories);
# it defaults to ./genome_data under AbioticStressConifers/. See SOURCES.tsv for origins.
DATA_ROOT = Path(os.environ.get("DATA_ROOT", PROJ / "genome_data"))

# Representative (longest-CDS, no-TE) gene models: the promoter TSS is taken from the
# single representative transcript per gene. This avoids spurious upstream exons in
# minor isoforms, which would otherwise displace the gene-feature start (and hence the
# 2 kb promoter window) far upstream of the true transcription start site.
PA_REPR     = DATA_ROOT / "sprucev2/Picab02_230926_at01_longest_no_TE_sorted.gff3"
PA_RM_GFF   = DATA_ROOT / "sprucev2/Picab02_repeats.gff3.gz"
PS_REPR_BED = PROJ / "data/annotation/Pinsy01_240308_at01_longest_no_TE.only_PASA.only_gene.cleaned_IDs.bed"
PS_RM_GFF   = DATA_ROOT / "pinev1/Pinsy01_chromosomes_and_unplaced.all_repeats.gff.gz"

PROMOTER_WIN = 2000

MOTIF_RE  = re.compile(r'Target "Motif:([^"]+)"')
CLASS_MAP = {  # broad class from family name hints
    "LTR": "LTR", "Gypsy": "LTR/Gypsy", "Copia": "LTR/Copia",
    "LINE": "LINE", "SINE": "SINE", "DNA": "DNA",
    "TIR": "DNA/TIR", "Helitron": "DNA/Helitron",
    "RC": "Rolling_circle", "Unknown": "Unknown"
}

def parse_tss(gff_file):
    genes = {}
    op = gzip.open if str(gff_file).endswith(".gz") else open
    with op(gff_file, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene": continue
            chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
            attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
            gid = attrs.get("Name") or re.sub(r"\.(v\d+\.\d+|\d+)$", "", attrs.get("ID",""))
            if gid: genes[gid] = (chrom, (start if strand=="+" else end), strand)
    print(f"  TSS: {len(genes):,} genes from {Path(gff_file).name}")
    return genes

def parse_tss_bed(bed_file):
    """TSS from a BED of representative (longest-CDS) gene models."""
    genes = {}
    with open(bed_file) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6: continue
            chrom, start, end, gid, strand = p[0], int(p[1]), int(p[2]), p[3], p[5]
            genes[gid] = (chrom, (start + 1) if strand == "+" else end, strand)
    print(f"  TSS: {len(genes):,} genes from {Path(bed_file).name}")
    return genes

def load_te_by_family(te_gff, species_prefix):
    """Return dict: chrom -> sorted list of (start, end, family, broad_class)"""
    by_chrom = defaultdict(list)
    op = gzip.open if str(te_gff).endswith(".gz") else open
    n = 0
    with op(te_gff, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9: continue
            try: s, e = int(p[3]), int(p[4])
            except: continue
            attrs_str = p[8]
            # Get family from Motif target
            m = MOTIF_RE.search(attrs_str)
            if m:
                family = m.group(1)
            else:
                attrs = dict(x.split("=",1) for x in attrs_str.split(";") if "=" in x)
                family = attrs.get("Name", attrs.get("ID", p[2]))
            # Broad class: from gff type or family name
            gff_type = p[2]
            if "LTR" in gff_type or "LTR" in family.upper(): broad = "LTR"
            elif "LINE" in gff_type or "LINE" in family.upper(): broad = "LINE"
            elif "SINE" in gff_type or "SINE" in family.upper(): broad = "SINE"
            elif "DNA" in gff_type or "TIR" in family.upper() or "Helitron" in family: broad = "DNA"
            else: broad = "Other/Unknown"
            by_chrom[p[0]].append((s, e, family, broad))
            n += 1
    # Sort by start
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    print(f"  TE intervals: {n:,} from {Path(te_gff).name}")
    return dict(by_chrom)

def promoter_te_families(gene_id, tss_dict, te_dict, win=PROMOTER_WIN):
    """Return (set of family names, set of broad classes) in promoter window."""
    if gene_id not in tss_dict: return set(), set()
    chrom, tss, strand = tss_dict[gene_id]
    if strand == "+": prom_s, prom_e = tss - win, tss
    else:             prom_s, prom_e = tss, tss + win
    if chrom not in te_dict: return set(), set()
    intervals = te_dict[chrom]
    starts = [x[0] for x in intervals]
    lo = bisect.bisect_left(starts, prom_s - win)  # generous search window
    families, broads = set(), set()
    for s, e, fam, broad in intervals[lo:]:
        if s > prom_e: break
        if e >= prom_s and s <= prom_e:
            families.add(fam)
            broads.add(broad)
    return families, broads

# ─────────────────────────────────────────────────────────────────────────────
print("Loading TSS data...")
pa_tss = parse_tss(PA_REPR)
ps_tss = parse_tss_bed(PS_REPR_BED)

print("Loading TE family data...")
pa_te  = load_te_by_family(PA_RM_GFF, "PA")
ps_te  = load_te_by_family(PS_RM_GFF, "PS")


# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS B: Between-species — promoter TE family Jaccard for 1:1 orthologs
#             vs co-expression conservation category
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Analysis B: Between-species promoter TE Jaccard vs co-expression ===")

backbone = pd.read_csv(INTEG / "integration_backbone_1to1.tsv", sep="\t",
                       usecols=["pa_gene","ps_gene","coex_category"])
backbone = backbone.dropna(subset=["coex_category","pa_gene","ps_gene"])
backbone = backbone.drop_duplicates(subset=["pa_gene","ps_gene"])
print(f"  Backbone: {len(backbone)} cross-species pairs")

B_rows = []
n_done = 0
for _, row in backbone.iterrows():
    pa_fam, pa_brd = promoter_te_families(row["pa_gene"], pa_tss, pa_te)
    ps_fam, ps_brd = promoter_te_families(row["ps_gene"], ps_tss, ps_te)
    union_fam = pa_fam | ps_fam
    shared_fam = pa_fam & ps_fam
    jac = len(shared_fam)/len(union_fam) if union_fam else np.nan
    pa_brd_str = "|".join(sorted(pa_brd)) if pa_brd else ""
    ps_brd_str = "|".join(sorted(ps_brd)) if ps_brd else ""
    B_rows.append({
        "pa_gene": row["pa_gene"], "ps_gene": row["ps_gene"],
        "coex_category": row["coex_category"],
        "pa_n_fam": len(pa_fam), "ps_n_fam": len(ps_fam),
        "n_shared_fam": len(shared_fam),
        "cross_species_te_jaccard": round(jac, 4) if not np.isnan(jac) else np.nan,
        "pa_broad_classes": pa_brd_str,
        "ps_broad_classes": ps_brd_str
    })
    n_done += 1
    if n_done % 2000 == 0: print(f"  {n_done}/{len(backbone)} done...", flush=True)

B_df = pd.DataFrame(B_rows)
# Cross-species promoter-TE Jaccard covariate used in the integrative model (Figure 7).
B_df.to_csv(INTEG/"cross_species_promoter_te_jaccard.tsv", sep="\t", index=False)
print(f"  Saved results/integration/cross_species_promoter_te_jaccard.tsv ({len(B_df)} rows)")

# Kruskal-Wallis across categories
sub_B = B_df.dropna(subset=["cross_species_te_jaccard"])
print(f"\n  Pairs with Jaccard: {len(sub_B)} / {len(B_df)}")
print("  Median cross-species TE Jaccard by co-expression category:")
print(sub_B.groupby("coex_category")["cross_species_te_jaccard"].agg(["median","mean","count"]))
kw = kruskal(*[sub_B.loc[sub_B.coex_category==c, "cross_species_te_jaccard"].values
               for c in sub_B.coex_category.unique() if len(sub_B[sub_B.coex_category==c]) > 5])
print(f"\n  Kruskal-Wallis: H={kw.statistic:.2f}, p={kw.pvalue:.4g}")

# Spearman: jaccard vs -log10(best_pval)
if "best_pval" in backbone.columns:
    merged_B = B_df.merge(backbone[["pa_gene","ps_gene","best_pval"]], on=["pa_gene","ps_gene"])
    sub_sp = merged_B.dropna(subset=["cross_species_te_jaccard","best_pval"])
    rho, pval = spearmanr(sub_sp["cross_species_te_jaccard"],
                          -np.log10(sub_sp["best_pval"].clip(1e-300)))
    print(f"  Spearman rho(cross-species TE Jaccard, -log10 co-expression pval): "
          f"{rho:.3f}, p={pval:.4g}")

print("\n=== Done ===")
