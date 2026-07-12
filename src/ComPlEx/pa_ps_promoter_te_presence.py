#!/usr/bin/env python3
"""
Per-gene promoter transposable-element occupancy for Norway spruce (Picab02) and
Scots pine (Pinsy01).

For every gene, flag whether at least one RepeatMasker TE annotation falls within
the 2 kb upstream of the transcription start site. The TSS is taken from the
representative (longest-CDS, no-TE) gene model of each gene, so that spurious
upstream exons in minor isoforms cannot displace the promoter window away from
the true start site.

Output:
  pa_ps_promoter_te_results.tsv  -- gene_id, species, te_in_promoter
"""

import bisect, gzip, os, re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

PROJ   = Path(__file__).resolve().parents[2]   # AbioticStressConifers/
OUTDIR = PROJ                                  # committed intermediates are written here
# External genome and repeat annotations are too large to redistribute with the code.
# Set DATA_ROOT to a directory holding them (with sprucev2/ and pinev1/ subdirectories);
# it defaults to ./genome_data under AbioticStressConifers/. See SOURCES.tsv for origins.
DATA_ROOT = Path(os.environ.get("DATA_ROOT", PROJ / "genome_data"))
SPRUCEV2  = DATA_ROOT / "sprucev2"
PINEV1    = DATA_ROOT / "pinev1"

# Representative (longest-CDS, no-TE) gene models; spruce is a sorted GFF (gene feature),
# pine is a BED of the representative gene models.
PA_REPR     = SPRUCEV2 / "Picab02_230926_at01_longest_no_TE_sorted.gff3"
PA_TE_GFF   = SPRUCEV2 / "Picab02_repeats.gff3.gz"
PS_REPR_BED = PROJ     / "data/annotation/Pinsy01_240308_at01_longest_no_TE.only_PASA.only_gene.cleaned_IDs.bed"
PS_TE_GFF   = PINEV1   / "Pinsy01_chromosomes_and_unplaced.all_repeats.gff.gz"

PROMOTER_WIN = 2000
# Longest full-length LTR retroelements are < 50 kb, so a hit starting more than this far
# upstream cannot overlap a 2 kb promoter window; used to bound the interval scan.
MAX_TE_SPAN = 50000


def parse_gff_tss(gff_file, gene_feature="gene"):
    genes = {}
    _open = gzip.open if str(gff_file).endswith(".gz") else open
    print(f"  Parsing TSS from {Path(gff_file).name} …", flush=True)
    with _open(gff_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != gene_feature:
                continue
            chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
            attrs = dict(item.split("=", 1) for item in p[8].split(";") if "=" in item)
            gid = attrs.get("Name") or re.sub(r"\.(v\d+\.\d+|\d+)$", "", attrs.get("ID", ""))
            if not gid:
                continue
            genes[gid] = (chrom, start if strand == "+" else end, strand)
    print(f"    {len(genes):,} genes", flush=True)
    return genes


def parse_tss_bed(bed_file):
    """TSS from a BED of representative (longest-CDS) gene models."""
    genes = {}
    print(f"  Parsing TSS from {Path(bed_file).name} …", flush=True)
    with open(bed_file) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            chrom, start, end, gid, strand = p[0], int(p[1]), int(p[2]), p[3], p[5]
            # BED start is 0-based; convert the + strand TSS to 1-based coordinates.
            genes[gid] = (chrom, (start + 1) if strand == "+" else end, strand)
    print(f"    {len(genes):,} genes", flush=True)
    return genes


def load_te_intervals(te_gff):
    """{chrom: (starts, ends)} with both arrays co-sorted by start."""
    raw = defaultdict(list)
    _open = gzip.open if str(te_gff).endswith(".gz") else open
    print(f"  Loading TE GFF: {Path(te_gff).name} …", flush=True)
    total = 0
    with _open(te_gff, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            raw[p[0]].append((s, e))
            total += 1
    te_by_chrom = {}
    for ch, pairs in raw.items():
        pairs.sort()
        starts = np.fromiter((x[0] for x in pairs), dtype=np.int64, count=len(pairs))
        ends   = np.fromiter((x[1] for x in pairs), dtype=np.int64, count=len(pairs))
        te_by_chrom[ch] = (starts, ends)
    print(f"    {total:,} TE intervals", flush=True)
    return te_by_chrom


def has_promoter_te(chrom, tss, strand, te_by_chrom, window=PROMOTER_WIN):
    """Does any TE overlap the 2 kb window upstream of the TSS?"""
    if strand == "+":
        prom_s, prom_e = max(1, tss - window), tss
    else:
        prom_s, prom_e = tss, tss + window
    entry = te_by_chrom.get(chrom)
    if entry is None:
        return False
    starts, ends = entry
    lo = bisect.bisect_left(starts, prom_s - MAX_TE_SPAN)
    hi = bisect.bisect_right(starts, prom_e)   # start <= prom_e guaranteed for lo..hi-1
    for i in range(lo, hi):
        if ends[i] >= prom_s:                  # overlap: start <= prom_e and end >= prom_s
            return True
    return False


def compute_presence(genes, te_gff, species_label):
    te_by_chrom = load_te_intervals(te_gff)
    print(f"  Scanning {len(genes):,} {species_label} promoters …", flush=True)
    rows = [{"gene_id": gid, "species": species_label,
             "te_in_promoter": has_promoter_te(chrom, tss, strand, te_by_chrom)}
            for gid, (chrom, tss, strand) in genes.items()]
    df = pd.DataFrame(rows)
    print(f"  {df.te_in_promoter.sum():,} / {len(df):,} "
          f"({100*df.te_in_promoter.mean():.1f}%) carry a TE in the 2 kb promoter", flush=True)
    return df


print("=== PA (P. abies / Picab02) ===", flush=True)
pa = compute_presence(parse_gff_tss(PA_REPR), PA_TE_GFF, "PA")
print("\n=== PS (P. sylvestris / Pinsy01) ===", flush=True)
ps = compute_presence(parse_tss_bed(PS_REPR_BED), PS_TE_GFF, "PS")

out = pd.concat([pa, ps], ignore_index=True)
out.to_csv(OUTDIR / "pa_ps_promoter_te_results.tsv", sep="\t", index=False)
print(f"\nSaved pa_ps_promoter_te_results.tsv ({len(out):,} genes)", flush=True)
