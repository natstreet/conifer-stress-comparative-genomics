#!/usr/bin/env python3
"""
Cross-species 1:1 ortholog PAIR LIST for P. abies vs P. sylvestris (the integration backbone).
Enumerates two pair sets:
  (a) SD gene pairs vs their pine ortholog(s)
  (b) all qualifying 1:1 single-copy background pairs

This emits the pair list only; all dN/dS in the paper is PAML yn00 (cross_species_dnds_yn00.tsv).
The CDS load and pairwise alignment serve as the pair-qualification filter: a pair is emitted iff
both CDS are present and alignable.

Output columns: hog_id, spruce_gene, pine_gene, pair_type, category

Usage:
  python3 cross_species_ks.py

Inputs are resolved under DATA_ROOT (default ./genome_data, as in SOURCES.tsv and the
rest of the pipeline): the genome CDS FASTAs live in DATA_ROOT/sprucev2 and DATA_ROOT/pinev1.
Override any single path with the SPRUCE_CDS / PINE_CDS / DATA_ROOT environment variables.
"""
import gzip
import re
import csv
import os
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJ = Path(__file__).resolve().parents[2]                  # repo root
DATA = os.environ.get("DATA_ROOT", str(PROJ / "genome_data"))
DEP  = os.environ.get("SPRUCE_PINE_DEPOSIT", os.getcwd())   # repo root: holds the deposited intermediates
SPRUCE_CDS_PATH = os.environ.get("SPRUCE_CDS", f"{DATA}/sprucev2/Picab02_230926_at01_longest_no_TE_cds.fa")
PINE_CDS_PATH   = os.environ.get("PINE_CDS", f"{DATA}/pinev1/Pinsy01_240308_at01_longest_no_TE_cds.fa")
ORTHO_FILE      = f"{DEP}/Picea_abies__v__Pinus_sylvestris.tsv"
SD_KAKS_FILE    = os.environ.get("SD_KAKS_FILE", f"{DEP}/kaks_results.tsv")
OUT_FILE        = os.environ.get("CROSS_SPECIES_KS_OUT", f"{DEP}/cross_species_ks.tsv")


def load_fasta(path):
    """Load a FASTA file (optionally gzipped) → dict: gene_id → seq."""
    seqs = {}
    opener = gzip.open if path.endswith('.gz') else open
    current_id = None
    buf = []
    with opener(path, 'rt') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(buf).upper()
                    buf = []
                header = line[1:].split()[0]
                # Strip .mRNA.X suffix to get gene-level ID
                current_id = re.sub(r'\.mRNA\.\d+$', '', header)
            else:
                buf.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(buf).upper()
    return seqs

def align_cds_pairwise(seq1, seq2):
    """Global codon alignment using BioPython PairwiseAligner.

    Used only as the pair-qualification filter: a pair qualifies iff the two CDS can be
    globally aligned. The returned alignment itself is not scored (no Ka/Ks is computed).
    """
    from Bio.Align import PairwiseAligner   # hard dependency: no silent fallback

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5
    try:
        best = next(iter(aligner.align(seq1, seq2)))
    except StopIteration:
        return None, None
    a1 = str(best[0])
    a2 = str(best[1])
    cols1, cols2 = [], []
    for b1, b2 in zip(a1, a2):
        if b1 == '-' and b2 == '-':
            continue
        cols1.append(b1)
        cols2.append(b2)
    aln1 = ''.join(cols1)
    aln2 = ''.join(cols2)
    trim = len(aln1) - (len(aln1) % 3)
    return aln1[:trim], aln2[:trim]

# ── Load CDS ──────────────────────────────────────────────────────────────────
print("Loading spruce CDS...", flush=True)
spruce_cds = load_fasta(SPRUCE_CDS_PATH)
print(f"  {len(spruce_cds):,} spruce CDS sequences", flush=True)

print("Loading pine CDS...", flush=True)
pine_cds = load_fasta(PINE_CDS_PATH)
print(f"  {len(pine_cds):,} pine CDS sequences", flush=True)

# ── Load ortholog table ───────────────────────────────────────────────────────
print("Loading ortholog table...", flush=True)
spruce_to_pine_1to1 = {}   # spruce_gene → pine_gene  (1:1 chromosomal only)
one_to_one_pairs = []

with open(ORTHO_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        og = row["Orthogroup"]
        def clean_genes(raw):
            return [re.sub(r'\.mRNA\.\d+$', '', x.strip())
                    for x in raw.split(",") if x.strip()]
        picea = clean_genes(row["Picea_abies"])
        pinus = clean_genes(row["Pinus_sylvestris"])
        pinus_chr = [g for g in pinus if not g.startswith("PS_sUP")]
        if len(picea) == 1 and len(pinus_chr) == 1:
            spruce_to_pine_1to1[picea[0]] = pinus_chr[0]
            one_to_one_pairs.append((picea[0], pinus_chr[0], og))

print(f"  {len(one_to_one_pairs):,} 1:1 chromosomal pairs", flush=True)

# ── Load SD gene pairs ────────────────────────────────────────────────────────
print("Loading SD pairs...", flush=True)
sd_pairs = []
with open(SD_KAKS_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sd_pairs.append((row["hog_id"], row["gene1"], row["gene2"], row["category"]))
print(f"  {len(sd_pairs):,} SD pairs", flush=True)

# ── Enumerate qualifying cross-species pairs ──────────────────────────────────
results = []

def qualify_pair(sp_gene, pine_gene, pair_type, hog_id="", category=""):
    """Emit a pair record iff both CDS are present AND globally alignable."""
    if sp_gene not in spruce_cds:
        return None
    if pine_gene not in pine_cds:
        return None
    aln1, aln2 = align_cds_pairwise(spruce_cds[sp_gene], pine_cds[pine_gene])
    if aln1 is None:
        return None
    return {
        "hog_id": hog_id,
        "spruce_gene": sp_gene,
        "pine_gene": pine_gene,
        "pair_type": pair_type,
        "category": category,
    }

# SD pairs vs pine orthologs
print(f"\nEnumerating cross-species pairs for {len(sd_pairs)} SD pairs...", flush=True)
done = 0
for hog_id, g1, g2, cat in sd_pairs:
    pine1 = spruce_to_pine_1to1.get(g1)
    pine2 = spruce_to_pine_1to1.get(g2)
    if pine1:
        r = qualify_pair(g1, pine1, "SD_gene1", hog_id, cat)
        if r:
            results.append(r)
    if pine2:
        r = qualify_pair(g2, pine2, "SD_gene2", hog_id, cat)
        if r:
            results.append(r)
    done += 1
    if done % 200 == 0:
        print(f"  {done}/{len(sd_pairs)} SD pairs done...", flush=True)

print(f"  SD pairs: {done} processed, {len(results)} results", flush=True)

# Background 1:1 orthologs: use all qualifying pairs (no sampling cap)
print(f"\nEnumerating background 1:1 pairs...", flush=True)
sd_genes = set()
for _, g1, g2, _ in sd_pairs:
    sd_genes.add(g1)
    sd_genes.add(g2)

bg_candidates = [(sg, pg, og) for sg, pg, og in one_to_one_pairs
                 if sg not in sd_genes
                 and sg in spruce_cds and pg in pine_cds]

print(f"  {len(bg_candidates):,} qualifying background pairs", flush=True)

bg_done = 0
for sg, pg, og in bg_candidates:
    r = qualify_pair(sg, pg, "background_1to1", og, "background")
    if r:
        results.append(r)
    bg_done += 1
    if bg_done % 100 == 0:
        print(f"  {bg_done}/{len(bg_candidates)} background pairs done...", flush=True)

print(f"  Background: {bg_done} processed", flush=True)

# ── Write output ───────────────────────────────────────────────────────────────
fieldnames = ["hog_id", "spruce_gene", "pine_gene", "pair_type", "category"]
with open(OUT_FILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

n_sd = sum(1 for r in results if r["pair_type"].startswith("SD"))
n_bg = sum(1 for r in results if r["pair_type"] == "background_1to1")
print(f"\nWrote {len(results)} rows to {OUT_FILE} (SD: {n_sd}, background_1to1: {n_bg})", flush=True)
