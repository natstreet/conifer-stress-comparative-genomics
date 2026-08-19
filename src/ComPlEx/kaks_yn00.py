#!/usr/bin/env python3
"""
kaks_yn00.py — codon-aware Ka/Ks (dN/dS) for cross-species 1:1 ortholog pairs, using a
protein-guided codon alignment and the Yang & Nielsen (2000) maximum-likelihood estimator
as implemented in PAML yn00 (the canonical reference implementation).

Pipeline per pair: translate CDS -> global protein alignment (BLOSUM62) -> back-translate
to a codon alignment (pal2nal-style) -> PAML yn00 -> dN, dS.

METHODS / REPRODUCIBILITY: dN/dS is estimated with PAML yn00 (version 4.10.10) under the
Yang & Nielsen (2000) maximum-likelihood model. yn00 control parameters set (non-default
noted): icode = 0 (standard/universal genetic code); weighting = 0 (no weighting of
substitution pathways); commonf3x4 = 0 (per-pair F3x4 codon frequencies, not a common set);
verbose = 0. Protein alignment: Biopython PairwiseAligner, global mode, BLOSUM62,
gap-open -11, gap-extend -1; the protein alignment is back-translated to the codon
alignment fed to yn00. (PAML yn00 is the canonical reference implementation; the Biopython
Bio.codonalign.cal_dn_ds YN00 is an experimental, non-canonical reimplementation and is not used.)

Usage:
  python3 kaks_yn00.py --pairs integration_backbone_1to1.tsv \
      --spruce-cds Picab02_..._cds.fa --pine-cds Pinsy01_..._cds.fa \
      --out kaks_yn00.tsv [--start 0 --limit N]
"""
import argparse, warnings, os, re, subprocess, tempfile
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
warnings.filterwarnings("ignore")

# PAML yn00 control parameters (Yang & Nielsen 2000). Non-default: none beyond the standard
# pairwise YN00 setup — icode=0 standard code, weighting=0, commonf3x4=0 (per-pair F3x4).
_YN00_CTL = ("seqfile = in.phy\noutfile = out\nverbose = 0\nicode = 0\n"
             "weighting = 0\ncommonf3x4 = 0\n")

def yn00_dn_ds(codon1, codon2):
    """Run PAML yn00 on a 2-sequence codon alignment (equal length, gaps as '-') and
    return (dN, dS) from the Yang & Nielsen (2000) table, or None on failure."""
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "in.phy"), "w") as f:
            f.write(f" 2 {len(codon1)}\nspruce    {codon1}\npine      {codon2}\n")
        with open(os.path.join(d, "yn00.ctl"), "w") as f:
            f.write(_YN00_CTL)
        try:
            subprocess.run(["yn00", "yn00.ctl"], cwd=d, capture_output=True, timeout=120)
            lines = open(os.path.join(d, "out")).read().splitlines()
        except Exception:
            return None
        for i, l in enumerate(lines):
            if re.search(r"\bS\s+N\s+t\s+kappa\s+omega\s+dN", l):
                for j in range(i + 1, min(i + 6, len(lines))):
                    p = lines[j].split()
                    if len(p) >= 11 and p[0].isdigit():
                        try:
                            return float(p[7]), float(p[10])   # dN, dS
                        except ValueError:
                            return None
    return None

def load_cds(path):
    d = {}
    for rec in SeqIO.parse(path, "fasta"):
        gene = None
        for f in rec.description.split():
            if f.startswith("gene="):
                gene = f[5:]; break
        if gene is None:
            gene = rec.id.split(".")[0]
        s = str(rec.seq).upper()
        d[gene] = s[: len(s) - (len(s) % 3)]          # trim to whole codons
    return d

_al = PairwiseAligner()
_al.substitution_matrix = substitution_matrices.load("BLOSUM62")
_al.open_gap_score = -11; _al.extend_gap_score = -1; _al.mode = "global"

def backtranslate(prot_aln, nuc):
    out = []; i = 0
    for aa in prot_aln:
        if aa == "-":
            out.append("---")
        else:
            out.append(nuc[i:i+3]); i += 3
    return "".join(out)

def kaks(cds1, cds2):
    p1 = str(Seq(cds1).translate(to_stop=False)).rstrip("*")
    p2 = str(Seq(cds2).translate(to_stop=False)).rstrip("*")
    if "*" in p1 or "*" in p2 or len(p1) < 20 or len(p2) < 20:
        return None
    aln = _al.align(p1, p2)[0]
    a1, a2 = format(aln).strip().split("\n")[0::2][0], None
    # robust: use aligned sequences directly
    a1 = aln[0]; a2 = aln[1]
    c1 = backtranslate(a1, cds1); c2 = backtranslate(a2, cds2)
    return yn00_dn_ds(c1, c2)          # (dN, dS) or None

def compute_row(task):
    g1, g2, cds1, cds2 = task
    try:
        res = kaks(cds1, cds2)
        if res is None:
            return (g1, g2, "", "", "")          # attempted, failed
        dN, dS = res
        kk = (dN / dS) if (dS and dS > 0) else float("nan")
        return (g1, g2, round(dN, 5), round(dS, 5),
                round(kk, 5) if kk == kk else "")
    except Exception:
        return (g1, g2, "", "", "")

def main():
    import os
    from multiprocessing import Pool
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--spruce-cds", required=True)
    ap.add_argument("--pine-cds", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()

    pairs = pd.read_csv(a.pairs, sep="\t")
    sp = load_cds(a.spruce_cds); pi = load_cds(a.pine_cds)

    # resume: skip pairs already written (by pa_gene+ps_gene key)
    done = set()
    if os.path.exists(a.out):
        prev = pd.read_csv(a.out, sep="\t")
        done = set(zip(prev.pa_gene, prev.ps_gene))

    tasks = []
    for _, r in pairs.iterrows():
        g1, g2 = r["pa_gene"], r["ps_gene"]
        if (g1, g2) in done:
            continue
        if g1 not in sp or g2 not in pi:
            continue
        tasks.append((g1, g2, sp[g1], pi[g2]))

    mode = "a" if done else "w"
    n = 0
    with open(a.out, mode) as fh:
        if not done:
            fh.write("pa_gene\tps_gene\tdN\tdS\tdNdS\n")
        with Pool(a.workers) as pool:
            for row in pool.imap_unordered(compute_row, tasks, chunksize=20):
                fh.write("\t".join(str(x) for x in row) + "\n")
                n += 1
                if n % 500 == 0:
                    fh.flush()
    print(f"wrote {n} new rows (had {len(done)} done) -> {a.out}")

    # Normalise the output to EXACTLY the current --pairs set, deterministically ordered. Resume above
    # appends and never removes, so a stale pre-existing output would otherwise keep orphan pairs from an
    # older --pairs list (e.g. a superseded backbone). Sorting also makes the row order reproducible
    # regardless of the parallel completion order.
    valid = set(zip(pairs["pa_gene"], pairs["ps_gene"]))
    out = pd.read_csv(a.out, sep="\t")
    out = out[[(g1, g2) in valid for g1, g2 in zip(out.pa_gene, out.ps_gene)]]
    out = out.drop_duplicates(subset=["pa_gene", "ps_gene"]).sort_values(["pa_gene", "ps_gene"])
    out.to_csv(a.out, sep="\t", index=False)
    print(f"normalised output to {len(out)} pairs (== current --pairs set, sorted)")

if __name__ == "__main__":
    main()
