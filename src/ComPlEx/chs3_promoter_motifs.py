#!/usr/bin/env python3
"""chs3_promoter_motifs.py — CHS3 panel-D motif comparison.

Computes the comparison at the 2 kb proximal promoter, for both CHS3 SD paralogs
(PA_chr09_G004115, PA_chr09_G004116). The 2 kb window matches the promoter window used
elsewhere in the paper; a wider 20 kb scan is dominated by background — for example, an
"AtMYB4 p = 7.5e-19" hit there is motif MP00480 (PlantTFDB family BBR-BPC, a GAGA-binder)
landing on a soft-masked AT/GA-rich simple repeat ~14 kb upstream of the G004115 TSS, a
low-complexity-composition artefact rather than a proximal cis-element.

Steps, for both paralogs:

  1. Extract the proximal 2 kb of each gene's upstream sequence. The deposited
     *_upstream_20kb.fa.gz is 5'->3' gene-relative with the TSS-proximal base LAST
     (verified against the genome: stored last-40 bp maps to genomic ~1,096,278,065,
     immediately above the G004115 representative-isoform TSS 1,096,278,064; stored
     first-40 bp maps ~20 kb away). So proximal 2 kb = the last 2000 bp.
  2. Run FIMO (PlantTFDB motifs, p < 1e-4) exactly as the genome-wide run did.
  3. Compute per-gene hit/motif counts, shared/union/Jaccard, and the top shared
     motifs (by best p across the pair), annotated with PlantTFDB family.

Outputs (results/integration/chs3/):
  chs3_promoter_2kb.fa            proximal 2 kb promoters (the FIMO input)
  chs3_fimo_2kb.tsv               raw FIMO hits (p < 1e-4)
  chs3_promoter_motif_2kb.tsv     per-motif comparison table (read by the figure)
  chs3_promoter_motif_summary.tsv one-row summary (counts, Jaccard) for the subtitle

Inputs (deposit / external; override the paths below if they live elsewhere):
  UPSTREAM_FA  spruce_upstream_20kb.fa.gz     (extract_upstream_seqs.py output)
  MOTIFS_MEME  PlantTFDB_motifs.meme          (PlantTFDB binding motifs, MEME format)
  ATH_INFO     Ath_TF_binding_motifs_information.txt  (Matrix_id -> Gene_id/Family)
  FIMO         path to the FIMO binary (MEME suite)

Run from the AbioticStressConifers/ project root:  python3 src/ComPlEx/chs3_promoter_motifs.py
"""

import csv
import gzip
import os
import subprocess
import sys
import textwrap
from collections import defaultdict
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
PROJ = Path(os.getcwd())
DEPOSIT = Path(os.environ.get("SPRUCE_PINE_DEPOSIT", os.getcwd()))   # repo root (driver sets this)

UPSTREAM_FA = Path(os.environ.get("UPSTREAM_FA", DEPOSIT / "spruce_upstream_20kb.fa.gz"))
MOTIFS_MEME = Path(os.environ.get("MOTIFS_MEME", DEPOSIT / "PlantTFDB_motifs.meme"))
ATH_INFO = Path(os.environ.get("ATH_INFO",
                DEPOSIT / "Ath_TF_binding_motifs_information.txt"))   # PlantTFDB, see SOURCES.tsv
KAKS_FILE = Path(os.environ.get("KAKS_FILE", DEPOSIT / "kaks_results.tsv"))
FIMO = os.environ.get("FIMO", os.path.expanduser(
    "~/miniforge3/envs/meme_env/bin/fimo"))

GENES = ["PA_chr09_G004115", "PA_chr09_G004116"]
PROX_BP = 2000        # proximal promoter window (matches the rest of the paper)
PVAL_THR = 1e-4       # FIMO significance threshold (as in run_fimo.sh)
N_TOP_SHARED = 6      # top shared motifs to tabulate for the figure

OUTDIR = PROJ / "results/integration/chs3"
OUTDIR.mkdir(parents=True, exist_ok=True)


# ── 1. Extract proximal 2 kb promoters ────────────────────────────────────────
def load_upstream(fa_gz, targets):
    seqs, cur, buf = {}, None, []
    with gzip.open(fa_gz, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in targets:
                    seqs[cur] = "".join(buf)
                cur = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur in targets:
            seqs[cur] = "".join(buf)
    return seqs


def write_proximal_fasta(seqs, out_fa):
    with open(out_fa, "w") as out:
        for g in GENES:
            prox = seqs[g][-PROX_BP:]          # proximal end is LAST (verified — see header)
            if len(prox) < PROX_BP:
                print(f"  WARNING: {g} upstream only {len(prox)} bp (< {PROX_BP})", flush=True)
            out.write(f">{g}\n" + "\n".join(textwrap.wrap(prox, 60)) + "\n")


# ── 2. Run FIMO ───────────────────────────────────────────────────────────────
def run_fimo(fa, out_tsv):
    with open(out_tsv, "w") as out:
        subprocess.run(
            [FIMO, "--thresh", str(PVAL_THR), "--max-stored-scores", "5000000",
             "--skip-matched-sequence", "--text", str(MOTIFS_MEME), str(fa)],
            check=True, stdout=out, stderr=subprocess.DEVNULL)


# ── 3. Parse + compare ────────────────────────────────────────────────────────
def parse_fimo(tsv):
    best = defaultdict(lambda: defaultdict(lambda: 1.0))  # gene -> motif -> best p
    hits = defaultdict(int)
    with open(tsv) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r, None)
        for row in r:
            if len(row) < 8:
                continue
            motif, gene = row[0], row[2]
            try:
                p = float(row[7])
            except ValueError:
                continue
            hits[gene] += 1
            if p < best[gene][motif]:
                best[gene][motif] = p
    return best, hits


def load_family_map(info):
    fam = {}
    with open(info) as fh:
        for d in csv.DictReader(fh, delimiter="\t"):
            fam[d["Matrix_id"]] = (d.get("Family", "?"), d.get("Gene_id", "?"))
    return fam


def main():
    for pth, name in [(UPSTREAM_FA, "UPSTREAM_FA"), (MOTIFS_MEME, "MOTIFS_MEME")]:
        if not Path(pth).exists():
            sys.exit(f"ERROR: {name} not found at {pth}")

    print(f"Extracting proximal {PROX_BP} bp promoters for {', '.join(GENES)} ...", flush=True)
    seqs = load_upstream(UPSTREAM_FA, set(GENES))
    missing = [g for g in GENES if g not in seqs]
    if missing:
        sys.exit(f"ERROR: genes not in {UPSTREAM_FA}: {missing}")
    fa = OUTDIR / "chs3_promoter_2kb.fa"
    write_proximal_fasta(seqs, fa)

    fimo_tsv = OUTDIR / "chs3_fimo_2kb.tsv"
    print(f"Running FIMO (p < {PVAL_THR}) ...", flush=True)
    run_fimo(fa, fimo_tsv)

    best, hits = parse_fimo(fimo_tsv)
    fam = load_family_map(ATH_INFO) if Path(ATH_INFO).exists() else {}
    g1, g2 = GENES
    s1, s2 = set(best[g1]), set(best[g2])
    inter, union = s1 & s2, s1 | s2
    jacc = len(inter) / len(union) if union else float("nan")

    # per-motif comparison table (shared motifs, best p per gene, family)
    def maxp(m):
        return max(best[g1][m], best[g2][m])
    rows = []
    for m in sorted(inter, key=maxp)[:N_TOP_SHARED]:
        family, gid = fam.get(m, ("?", "?"))
        rows.append(dict(motif_id=m, family=family, ath_locus=gid,
                         G004115_best_p=f"{best[g1][m]:.3e}",
                         G004116_best_p=f"{best[g2][m]:.3e}"))
    with open(OUTDIR / "chs3_promoter_motif_2kb.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # Ka/Ks for the pair, read from the YN00 producer output (not hard-coded)
    kaks = "NA"
    if KAKS_FILE.exists():
        with open(KAKS_FILE) as fh:
            for d in csv.DictReader(fh, delimiter="\t"):
                pair = {d.get("gene1"), d.get("gene2")}
                if pair == set(GENES):
                    kaks = f"{float(d['KaKs']):.3f}"
                    break

    # one-row summary for the figure subtitle / counts
    summ = dict(gene1=g1, gene2=g2, window_bp=PROX_BP, pval_thr=PVAL_THR,
                g1_hits=hits[g1], g2_hits=hits[g2],
                g1_motifs=len(s1), g2_motifs=len(s2),
                shared_motifs=len(inter), union_motifs=len(union),
                jaccard_pct=round(100 * jacc, 1),
                g1_unique=len(s1 - s2), g2_unique=len(s2 - s1),
                kaks=kaks)
    with open(OUTDIR / "chs3_promoter_motif_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(summ)

    print("\n2 kb proximal-promoter TF-motif comparison:")
    print(f"  {g1}: {hits[g1]} hits, {len(s1)} motifs;  {g2}: {hits[g2]} hits, {len(s2)} motifs")
    print(f"  shared {len(inter)} / union {len(union)}  ->  Jaccard {summ['jaccard_pct']}% shared")
    print(f"  unique: {summ['g1_unique']} to {g1}, {summ['g2_unique']} to {g2}")
    print(f"  top shared motif: {rows[0]['motif_id']} ({rows[0]['family']}, {rows[0]['ath_locus']}), "
          f"best p {rows[0]['G004115_best_p']}/{rows[0]['G004116_best_p']}")
    print(f"\nWrote: {OUTDIR}/chs3_promoter_motif_2kb.tsv, chs3_promoter_motif_summary.tsv")


if __name__ == "__main__":
    main()
