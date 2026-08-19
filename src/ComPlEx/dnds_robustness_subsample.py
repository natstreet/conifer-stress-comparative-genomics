#!/usr/bin/env python3
"""dnds_robustness_subsample.py

Robustness check for the dN/dS-vs-co-expression-conservation gradient. The primary estimate
(cross_species_dnds_yn00.tsv) uses PAML yn00 on a protein-guided back-translated codon alignment.
This re-estimates a STRATIFIED SUBSAMPLE of 1:1 spruce-pine pairs with two independent improvements:
  (a) codon-aware alignment with PRANK (-codon) instead of protein-guided back-translation, and
  (b) codeml pairwise ML (runmode=-2, model M0, seqtype=1) instead of the yn00 approximation,
then compares per-category medians, the dN/dS-vs-conservation-breadth Spearman rho, and the per-pair
correlation against the yn00 values on the same pairs. Subsample only; fixed seed.

Env overrides: INTEG_DIR, PRANK_BIN, CODEML_BIN, SPRUCE_CDS, PINE_CDS, PER_CAT, SEED, NPROC.
"""
import os, sys, gzip, random, subprocess, tempfile, shutil, re, csv, math
from statistics import median
from multiprocessing import Pool
from Bio import SeqIO

I          = os.environ.get("INTEG_DIR", "results/integration")
PRANK      = os.environ.get("PRANK_BIN", "/Users/rona006/miniforge3/envs/dnds_robust/bin/prank")
CODEML     = os.environ.get("CODEML_BIN", "codeml")
GENOME_DIR = "/Users/rona006/Library/CloudStorage/OneDrive-Umeåuniversitet/work/manuscripts/spruce2/Nature Genetics/RESUBMISSION November 2025/FigShare_NewCompleteSet/Genome"
SPRUCE_CDS = os.environ.get("SPRUCE_CDS", os.path.join(GENOME_DIR, "Picab02_230926_at01_longest_no_TE_cds.fa.gz"))
PINE_CDS   = os.environ.get("PINE_CDS",   os.path.join(GENOME_DIR, "Pinsy01_240308_at01_longest_no_TE_cds.fa.gz"))
PER_CAT    = int(os.environ.get("PER_CAT", "350"))
SEED       = int(os.environ.get("SEED", "42"))
NPROC      = int(os.environ.get("NPROC", "6"))
CATS       = ["conserved", "cold_specific", "drought_specific", "multi_tissue", "not_coex"]

def load_cds(path):
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    d = {}
    with op as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            gene = None
            for f in rec.description.split():
                if f.startswith("gene="):
                    gene = f[5:]; break
            if gene is None:
                gene = rec.id.split(".")[0]
            s = str(rec.seq).upper()
            d[gene] = s[: len(s) - (len(s) % 3)]
    return d

def prank_codon_align(cds1, cds2, wd):
    inp = os.path.join(wd, "in.fa")
    with open(inp, "w") as fh:
        fh.write(">s1\n%s\n>s2\n%s\n" % (cds1, cds2))
    subprocess.run([PRANK, "-d=" + inp, "-o=" + os.path.join(wd, "out"), "-codon", "-F", "-once"],
                   cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
    best = os.path.join(wd, "out.best.fas")
    if not os.path.exists(best):
        return None
    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(best, "fasta")}
    return (seqs.get("s1"), seqs.get("s2")) if len(seqs) == 2 else None

def aln_identity(a1, a2):
    """% identical over aligned columns where BOTH are non-gap (gap-excluded identity)."""
    n = ident = 0
    for x, y in zip(a1, a2):
        if x == "-" or y == "-":
            continue
        n += 1
        if x == y:
            ident += 1
    return (100.0 * ident / n) if n else 0.0, n  # identity%, aligned-codon-columns count*3(nt)

CTL = """seqfile = aln.phy
treefile = tree.nwk
outfile = codeml.out
noisy = 0
verbose = 0
runmode = -2
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 0
"""

def codeml_pairwise(a1, a2, wd):
    L = len(a1)
    with open(os.path.join(wd, "aln.phy"), "w") as fh:
        fh.write(" 2 %d\n" % L)
        fh.write("s1  %s\n" % a1)
        fh.write("s2  %s\n" % a2)
    with open(os.path.join(wd, "tree.nwk"), "w") as fh:
        fh.write("(s1,s2);\n")
    with open(os.path.join(wd, "codeml.ctl"), "w") as fh:
        fh.write(CTL)
    try:
        subprocess.run([CODEML, "codeml.ctl"], cwd=wd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=300)
    except Exception:
        return None
    out = os.path.join(wd, "codeml.out")
    if not os.path.exists(out):
        return None
    txt = open(out).read()
    m = re.search(r"dN/dS\s*=\s*([\d.]+)\s+dN\s*=\s*([\d.]+)\s+dS\s*=\s*([\d.]+)", txt)
    if not m:
        return None
    return float(m.group(1)), float(m.group(2)), float(m.group(3))  # omega, dN, dS

def worker(task):
    pa, ps, cat, ntis, yn, cds1, cds2 = task     # CDS passed in task (spawn-safe, no shared globals)
    if not cds1 or not cds2:
        return (pa, ps, cat, ntis, yn, "", "", "", "no_cds")
    wd = tempfile.mkdtemp(prefix="dndsrb_")
    try:
        al = prank_codon_align(cds1, cds2, wd)
        if not al or not al[0]:
            return (pa, ps, cat, ntis, yn, "", "", "", "aln_fail")
        ident, ncol = aln_identity(al[0], al[1])
        cm = codeml_pairwise(al[0], al[1], wd)
        if cm is None:
            return (pa, ps, cat, ntis, yn, "", round(ident, 1), "", "codeml_fail")
        omega, dN, dS = cm
        return (pa, ps, cat, ntis, yn, omega, round(ident, 1), dS, "ok")
    finally:
        shutil.rmtree(wd, ignore_errors=True)

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs); mx = sum(rx) / n; my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)) * sum((ry[i] - my) ** 2 for i in range(n)))
    return num / den if den else float("nan")

def main():
    global CDS_SP, CDS_PI
    # ---- join yn00 pairs with category + conservation breadth, stratified subsample ----
    yn = {}
    with open(os.path.join(I, "cross_species_dnds_yn00.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                dS = float(row["dS"]); dnds = float(row["dNdS"])
            except (ValueError, KeyError):
                continue
            if 0 < dS < 5 and dnds < 10:
                yn[(row["pa_gene"], row["ps_gene"])] = dnds
    cat_of, ntis_of = {}, {}
    with open(os.path.join(I, "integration_backbone_1to1.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            k = (row["pa_gene"], row["ps_gene"])
            cat_of[k] = row["coex_category"]
            ntis_of[k] = int(row["n_tissues"]) if row.get("n_tissues", "").strip() else 0
    pool = {c: [] for c in CATS}
    for k, dnds in yn.items():
        c = cat_of.get(k)
        if c in pool:
            pool[c].append((k[0], k[1], c, ntis_of[k], dnds))
    rng = random.Random(SEED)
    sample = []
    for c in CATS:
        # sort deterministically first so the seeded draw is independent of the input file's row order
        # (cross_species_dnds_yn00.tsv is written in parallel, so its row order is not stable run-to-run)
        rows = sorted(pool[c], key=lambda r: (r[0], r[1]))
        rng.shuffle(rows)
        sample.extend(rows[:PER_CAT])
    print("stratified subsample (seed %d, target %d/cat):" % (SEED, PER_CAT), flush=True)
    for c in CATS:
        print("  %-18s available %5d  drawn %4d" % (c, len(pool[c]), min(PER_CAT, len(pool[c]))), flush=True)
    print("  TOTAL drawn:", len(sample), flush=True)

    print("loading CDS...", flush=True)
    CDS_SP = load_cds(SPRUCE_CDS); CDS_PI = load_cds(PINE_CDS)
    print("  spruce CDS:", len(CDS_SP), " pine CDS:", len(CDS_PI), flush=True)

    tasks = [(pa, ps, cat, ntis, yn, CDS_SP.get(pa), CDS_PI.get(ps))
             for (pa, ps, cat, ntis, yn) in sample]
    with Pool(NPROC) as p:
        res = p.map(worker, tasks)

    outp = os.path.join(I, "dnds_robustness_subsample.tsv")
    with open(outp, "w") as fh:
        fh.write("pa_gene\tps_gene\tcoex_category\tn_tissues\tyn00_dNdS\tcodeml_dNdS\taln_identity_pct\tcodeml_dS\tstatus\n")
        for row in res:
            fh.write("\t".join(str(x) for x in row) + "\n")
    ok = [r for r in res if r[8] == "ok" and r[5] != "" and float(r[7]) < 5]
    print("\nusable (ok + codeml dS<5): %d / %d" % (len(ok), len(sample)), flush=True)

    # ---- summary ----
    summ = os.path.join(I, "dnds_robustness_summary.tsv")
    with open(summ, "w") as fh:
        fh.write("metric\tcategory\tyn00\tcodeml\n")
        for c in CATS:
            yv = [r[4] for r in ok if r[2] == c]
            cv = [float(r[5]) for r in ok if r[2] == c]
            if yv and cv:
                fh.write("median_dNdS\t%s\t%.4f\t%.4f\n" % (c, median(yv), median(cv)))
        # gradient: dN/dS vs conservation breadth (n_tissues), both methods, on the same usable pairs
        nt = [r[3] for r in ok]
        rho_yn = spearman(nt, [r[4] for r in ok])
        rho_cm = spearman(nt, [float(r[5]) for r in ok])
        rho_pair = spearman([r[4] for r in ok], [float(r[5]) for r in ok])
        fh.write("spearman_dNdS_vs_breadth\tALL\t%.4f\t%.4f\n" % (rho_yn, rho_cm))
        fh.write("per_pair_yn00_vs_codeml_spearman\tALL\t\t%.4f\n" % rho_pair)
    print("wrote", outp, "and", summ, flush=True)
    print("\n=== per-category median dN/dS (yn00 vs codeml) ===", flush=True)
    for c in CATS:
        yv = [r[4] for r in ok if r[2] == c]; cv = [float(r[5]) for r in ok if r[2] == c]
        if yv and cv:
            print("  %-18s yn00 %.3f  codeml %.3f  (n=%d)" % (c, median(yv), median(cv), len(yv)), flush=True)
    print("gradient rho (dN/dS vs breadth): yn00 %.3f  codeml %.3f | per-pair rho %.3f"
          % (rho_yn, rho_cm, rho_pair), flush=True)

    # ---- per-pair scatter -> Supplementary Figure S6 (plain-language legend to match the main text) ----
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        cols = {"conserved": "#1b9e77", "cold_specific": "#7570b3", "drought_specific": "#d95f02",
                "multi_tissue": "#666666", "not_coex": "#e7298a"}
        labels = {"conserved": "conserved", "cold_specific": "cold-specific",
                  "drought_specific": "drought-specific", "multi_tissue": "multi-tissue",
                  "not_coex": "not co-expressed"}
        fig, ax = plt.subplots(figsize=(5, 5))
        for c in CATS:
            xs = [r[4] for r in ok if r[2] == c]; ys = [float(r[5]) for r in ok if r[2] == c]
            ax.scatter(xs, ys, s=8, alpha=0.5, label=labels[c], color=cols[c], edgecolors="none")
        lim = max(max(r[4] for r in ok), max(float(r[5]) for r in ok)) * 1.02
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, label="y = x")
        ax.set(xlim=(0, lim), ylim=(0, lim),
               xlabel="dN/dS  (yn00, protein-guided back-translation)",
               ylabel="dN/dS  (codeml M0, PRANK codon-aware alignment)")
        ax.set_title("Per-pair dN/dS: yn00 vs codeml (n=%d, Spearman rho=%.2f)" % (len(ok), rho_pair), fontsize=9)
        ax.legend(fontsize=6, markerscale=1.5, frameon=False)
        fig.tight_layout()
        figdir = "results/final_figures"; os.makedirs(figdir, exist_ok=True)
        fig.savefig(os.path.join(figdir, "FigureS6.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(figdir, "FigureS6.png"), dpi=150, bbox_inches="tight")
        print("wrote results/final_figures/FigureS6.pdf/.png", flush=True)
    except Exception as e:
        print("scatter skipped:", e, flush=True)

if __name__ == "__main__":
    main()
