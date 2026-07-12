#!/usr/bin/env python3
"""
Extract upstream (promoter) sequences for all chromosome-anchored protein-coding
genes in P. abies and P. sylvestris for FIMO motif scanning.

For each gene:
  + strand: region from max(previous_gene_end + 1, TSS - MAX_UPSTREAM) to TSS - 1
  - strand: region from TES + 1 to min(next_gene_start - 1, TES + MAX_UPSTREAM)

TSS/TES is taken from the REPRESENTATIVE ISOFORM (longest CDS; longest concatenated
exons when no CDS), not from the gene-level feature. Gene-level coordinates are
retained only for neighbour-boundary clipping.

Uses strand-unaware nearest-neighbour boundaries (conservative: stop at any
adjacent gene regardless of strand) to avoid reading into coding sequence.

Only chromosome-anchored genes (PA_chrXX / PS_chrXX) are included.
Softmasked genomes are used (lowercase = repeat-masked) so FIMO skips repeats.

Outputs (gzip compressed FASTA):
  spruce_upstream_20kb.fa.gz
  pine_upstream_20kb.fa.gz
"""

import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import os
OUTDIR      = Path(os.environ.get("SPRUCE_PINE_DEPOSIT", os.getcwd()))   # repo root (driver sets this)
MAX_UPSTREAM = 20_000  # bp

# Softmasked genome assemblies + all-isoform GFFs (see SOURCES.tsv). Under $DATA_ROOT
# (default ./genome_data, holding sprucev2/ + pinev1/, as elsewhere in the pipeline).
DATA_ROOT = Path(os.environ.get("DATA_ROOT",
                 str(Path(__file__).resolve().parents[2] / "genome_data")))

CONFIGS = [
    dict(
        species   = "PA",
        gff       = DATA_ROOT / "sprucev2/Picab02_230926_at01_all_sorted.gff3.gz",
        genome    = DATA_ROOT / "sprucev2/Picab02_chromosomes_and_unplaced.softmasked.fa.gz",
        chr_prefix= "PA_chr",
        out       = OUTDIR / "spruce_upstream_20kb.fa.gz",
    ),
    dict(
        species   = "PS",
        gff       = DATA_ROOT / "pinev1/Pinsy01_240308_at01_all.gff3.gz",
        genome    = DATA_ROOT / "pinev1/Pinsy01_chromosomes_and_unplaced_softmasked.fasta.gz",
        chr_prefix= "PS_chr",
        out       = OUTDIR / "pine_upstream_20kb.fa.gz",
    ),
]


# ── Parse GFF ────────────────────────────────────────────────────────────────
def parse_genes(gff_path, chr_prefix):
    """Return dict: chrom → sorted list of (gene_start, gene_end, strand, gene_id,
    rep_tss, rep_tes) where rep_tss/rep_tes are from the representative isoform.

    Representative isoform rule (mirrors annotation filtering):
      - If any isoform has CDS: keep the one with the longest total CDS.
      - If no isoforms have CDS: keep the one with the longest concatenated exons.
    TSS = isoform start (+) or isoform end (-).
    Gene-level start/end retained for neighbour-boundary clipping only.
    """
    # First pass: collect gene features
    gene_info = {}   # gene_id -> (chrom, gstart, gend, strand)
    mrna_gene = {}   # mrna_id -> gene_id
    mrna_span = defaultdict(lambda: [10**15, 0])  # mrna_id -> [min_start, max_end]
    cds_len   = defaultdict(int)   # mrna_id -> total CDS bp
    exon_len  = defaultdict(int)   # mrna_id -> total exon bp
    gene_mrnas = defaultdict(list) # gene_id -> [mrna_id, ...]

    with gzip.open(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0]
            if not chrom.startswith(chr_prefix):
                continue
            ftype  = cols[2]
            s, e   = int(cols[3]), int(cols[4])
            strand = cols[6]
            attrs  = cols[8]

            if ftype == "gene":
                m = re.search(r'ID=([^;]+)', attrs)
                if m:
                    gid = m.group(1)
                    gene_info[gid] = (chrom, s, e, strand)

            elif ftype == "mRNA":
                pm = re.search(r'Parent=([^;]+)', attrs)
                im = re.search(r'ID=([^;]+)', attrs)
                if pm and im:
                    mid = im.group(1)
                    gid = pm.group(1)
                    mrna_gene[mid] = gid
                    gene_mrnas[gid].append(mid)
                    mrna_span[mid][0] = min(mrna_span[mid][0], s)
                    mrna_span[mid][1] = max(mrna_span[mid][1], e)

            elif ftype == "CDS":
                pm = re.search(r'Parent=([^;]+)', attrs)
                if pm:
                    cds_len[pm.group(1)] += e - s + 1

            elif ftype == "exon":
                pm = re.search(r'Parent=([^;]+)', attrs)
                if pm:
                    exon_len[pm.group(1)] += e - s + 1

    # Select representative isoform per gene
    genes_by_chrom = defaultdict(list)
    for gid, (chrom, gstart, gend, strand) in gene_info.items():
        mrnas = gene_mrnas.get(gid, [])
        if not mrnas:
            # No mRNA features — fall back to gene-level coordinates
            genes_by_chrom[chrom].append((gstart, gend, strand, gid, gstart, gend))
            continue

        # Pick representative: prefer longest CDS, else longest exon
        has_cds = [m for m in mrnas if cds_len[m] > 0]
        if has_cds:
            rep = max(has_cds, key=lambda m: cds_len[m])
        else:
            rep = max(mrnas, key=lambda m: exon_len[m])

        rep_s = mrna_span[rep][0]
        rep_e = mrna_span[rep][1]
        # TSS/TES from representative isoform; gene bounds for boundary clipping
        genes_by_chrom[chrom].append((gstart, gend, strand, gid, rep_s, rep_e))

    for ch in genes_by_chrom:
        genes_by_chrom[ch].sort()   # sort by gene start
    return genes_by_chrom


# ── Compute upstream intervals ────────────────────────────────────────────────
def upstream_intervals(genes_by_chrom, max_up=MAX_UPSTREAM):
    """Return list of (chrom, up_start, up_end, strand, gene_id) [1-based, inclusive].

    TSS/TES from representative isoform; neighbour boundaries from gene-level coords.
    """
    intervals = []
    for chrom, genes in genes_by_chrom.items():
        n = len(genes)
        for i, (gstart, gend, strand, gid, rep_s, rep_e) in enumerate(genes):
            # Neighbour boundaries from representative isoform coords (not gene-level),
            # so that long artefactual isoforms (e.g. mRNA.1 with a 147 kb intron) do
            # not block the upstream region of adjacent genes.
            prev_end   = genes[i - 1][5] if i > 0 else 0        # previous rep_e
            next_start = genes[i + 1][4] if i < n - 1 else 10**15 # next rep_s

            if strand == "+":
                tss  = rep_s            # representative isoform start
                up_s = max(prev_end + 1, tss - max_up)
                up_e = tss - 1
            else:
                tes  = rep_e            # representative isoform end (highest coord = TSS)
                up_s = tes + 1
                up_e = min(next_start - 1, tes + max_up)

            if up_s > up_e:
                continue
            intervals.append((chrom, up_s, up_e, strand, gid))
    return intervals


# ── Stream genome and extract sequences ─────────────────────────────────────
def extract_sequences(genome_path, intervals, out_path):
    """Stream genome FASTA; for each chromosome, extract all requested intervals."""
    # Index intervals by chromosome
    chrom_ivs = defaultdict(list)
    for iv in intervals:
        chrom_ivs[iv[0]].append(iv)

    n_written = 0
    with gzip.open(genome_path, "rt") as genome_fh, \
         gzip.open(out_path, "wt") as out_fh:

        current_chrom = None
        current_seq   = []
        chroms_done   = set()

        def flush(chrom, seq_parts):
            if chrom is None or chrom not in chrom_ivs:
                return 0
            seq = "".join(seq_parts)  # full chromosome sequence
            written = 0
            for (_, up_s, up_e, strand, gid) in chrom_ivs[chrom]:
                s0 = up_s - 1          # convert to 0-based
                s1 = up_e              # exclusive
                subseq = seq[s0:s1]
                if len(subseq) < 10:   # skip trivially short sequences
                    continue
                if strand == "-":
                    subseq = reverse_complement(subseq)
                out_fh.write(f">{gid}\n")
                # Write in 60-char lines
                for j in range(0, len(subseq), 60):
                    out_fh.write(subseq[j:j+60] + "\n")
                written += 1
            return written

        for line in genome_fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_chrom:
                    n_written += flush(current_chrom, current_seq)
                    chroms_done.add(current_chrom)
                current_chrom = line[1:].split()[0]
                current_seq   = []
            else:
                current_seq.append(line)

        # Flush last chromosome
        if current_chrom:
            n_written += flush(current_chrom, current_seq)

    return n_written


def reverse_complement(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


# ── Main ──────────────────────────────────────────────────────────────────────
for cfg in CONFIGS:
    sp = cfg["species"]
    print(f"\n{'='*60}", flush=True)
    print(f"Processing {sp}...", flush=True)

    print(f"  Parsing GFF (representative isoforms)...", flush=True)
    genes = parse_genes(cfg["gff"], cfg["chr_prefix"])
    n_genes = sum(len(v) for v in genes.values())
    n_chroms = len(genes)
    print(f"  {n_genes:,} genes on {n_chroms} chromosomes", flush=True)

    print(f"  Computing upstream intervals...", flush=True)
    ivs = upstream_intervals(genes)
    print(f"  {len(ivs):,} upstream intervals (max {MAX_UPSTREAM} bp)", flush=True)

    print(f"  Streaming genome and extracting sequences...", flush=True)
    print(f"  (reading {cfg['genome'].name}, this takes ~15 min per species)", flush=True)
    sys.stdout.flush()

    n = extract_sequences(cfg["genome"], ivs, cfg["out"])
    print(f"  Wrote {n:,} sequences → {cfg['out'].name}", flush=True)

print("\nDone.", flush=True)
