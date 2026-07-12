#!/usr/bin/env python3
"""
Ka/Ks (dN/dS) analysis for spruce paralog pairs from segmental duplications.
Uses PAML yn00 (Yang & Nielsen 2000, ML) — the same estimator as the cross-species
dN/dS (kaks_yn00.py -> cross_species_dnds_yn00.tsv), for methodological consistency:
protein-guided codon alignment, then yn00.
Requires the `yn00` binary (PAML) on PATH.
"""

import gzip
import re
import math
import sys
from Bio.Align import PairwiseAligner
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from kaks_yn00 import kaks as paml_yn00_dnds   # protein-guided codon alignment + PAML yn00 -> (dN, dS)

# ─────────────────────────────────────────────────────────────────────────────
# CODON TABLE
# ─────────────────────────────────────────────────────────────────────────────

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

BASES = ['A', 'T', 'G', 'C']


def synonymous_sites_per_codon(codon):
    """
    Count synonymous and non-synonymous sites in a codon (NG86 method).
    Returns (S, N) where S + N = 3.
    """
    codon = codon.upper()
    if codon not in CODON_TABLE or CODON_TABLE[codon] == '*':
        return 0.0, 3.0
    aa = CODON_TABLE[codon]
    syn_fractions = []
    for pos in range(3):
        n_changes = 0
        n_syn = 0
        orig_base = codon[pos]
        for base in BASES:
            if base == orig_base:
                continue
            mutant = codon[:pos] + base + codon[pos+1:]
            if mutant not in CODON_TABLE:
                continue
            n_changes += 1
            mut_aa = CODON_TABLE[mutant]
            if mut_aa == aa:
                n_syn += 1
        if n_changes > 0:
            syn_fractions.append(n_syn / n_changes)
        else:
            syn_fractions.append(0.0)
    S = sum(syn_fractions)
    N = 3.0 - S
    return S, N


def count_differences(c1, c2):
    """
    Count synonymous (sd) and non-synonymous (nd) differences between two codons.
    Uses NG86: for multi-hit codons, average over possible pathways.
    Returns (sd, nd).
    """
    c1, c2 = c1.upper(), c2.upper()
    if c1 == c2:
        return 0.0, 0.0
    if c1 not in CODON_TABLE or c2 not in CODON_TABLE:
        return 0.0, 0.0
    if CODON_TABLE[c1] == '*' or CODON_TABLE[c2] == '*':
        return 0.0, 0.0

    diff_positions = [i for i in range(3) if c1[i] != c2[i]]
    n_diffs = len(diff_positions)

    if n_diffs == 1:
        pos = diff_positions[0]
        aa1 = CODON_TABLE[c1]
        aa2 = CODON_TABLE[c2]
        if aa1 == aa2:
            return 1.0, 0.0
        else:
            return 0.0, 1.0

    elif n_diffs == 2:
        # Two possible pathways; average
        total_sd = 0.0
        total_nd = 0.0
        n_valid_paths = 0
        for first_pos in diff_positions:
            # Intermediate: change first_pos in c1 to match c2
            second_pos = [p for p in diff_positions if p != first_pos][0]
            intermediate = list(c1)
            intermediate[first_pos] = c2[first_pos]
            intermediate = ''.join(intermediate)
            if intermediate not in CODON_TABLE:
                continue
            inter_aa = CODON_TABLE[intermediate]
            # Step 1: c1 -> intermediate
            sd1 = 1.0 if CODON_TABLE[c1] == inter_aa else 0.0
            nd1 = 1.0 - sd1
            # Step 2: intermediate -> c2
            sd2 = 1.0 if inter_aa == CODON_TABLE[c2] else 0.0
            nd2 = 1.0 - sd2
            total_sd += sd1 + sd2
            total_nd += nd1 + nd2
            n_valid_paths += 1
        if n_valid_paths == 0:
            return 0.0, float(n_diffs)
        return total_sd / n_valid_paths, total_nd / n_valid_paths

    else:  # n_diffs == 3
        # Six possible pathways (3! = 6)
        from itertools import permutations
        total_sd = 0.0
        total_nd = 0.0
        n_valid_paths = 0
        for perm in permutations(diff_positions):
            path_sd = 0.0
            path_nd = 0.0
            current = list(c1)
            valid = True
            for pos in perm:
                current[pos] = c2[pos]
                codon_str = ''.join(current)
                if codon_str not in CODON_TABLE:
                    valid = False
                    break
                prev = list(current)
                prev[pos] = c1[pos]
                prev_codon = ''.join(prev)
                if prev_codon not in CODON_TABLE:
                    valid = False
                    break
                if CODON_TABLE[prev_codon] == CODON_TABLE[codon_str]:
                    path_sd += 1.0
                else:
                    path_nd += 1.0
            if valid:
                total_sd += path_sd
                total_nd += path_nd
                n_valid_paths += 1
        if n_valid_paths == 0:
            return 0.0, float(n_diffs)
        return total_sd / n_valid_paths, total_nd / n_valid_paths


def ng86_kaks(seq1, seq2):
    """
    Compute Ka, Ks, Ka/Ks using NG86 method with JC correction.
    seq1 and seq2 must be aligned CDS of the same length (multiple of 3).
    Returns (Ka, Ks, KaKs) or (None, None, None) on failure.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    if len(seq1) != len(seq2) or len(seq1) % 3 != 0:
        return None, None, None

    S_total = 0.0
    N_total = 0.0
    Sd_total = 0.0
    Nd_total = 0.0

    n_codons = len(seq1) // 3
    for i in range(n_codons):
        c1 = seq1[3*i:3*i+3]
        c2 = seq2[3*i:3*i+3]
        # Skip if contains gap or ambiguous
        if '-' in c1 or '-' in c2:
            continue
        if any(b not in 'ATGC' for b in c1 + c2):
            continue
        if CODON_TABLE.get(c1, '*') == '*' or CODON_TABLE.get(c2, '*') == '*':
            continue

        S1, N1 = synonymous_sites_per_codon(c1)
        S2, N2 = synonymous_sites_per_codon(c2)
        S_total += (S1 + S2) / 2.0
        N_total += (N1 + N2) / 2.0

        sd, nd = count_differences(c1, c2)
        Sd_total += sd
        Nd_total += nd

    if S_total < 1.0 or N_total < 1.0:
        return None, None, None

    # JC correction
    p_s = Sd_total / S_total
    p_n = Nd_total / N_total

    # Guard against saturation
    if p_s >= 0.75 or p_n >= 0.75:
        return None, None, None

    try:
        if p_s == 0.0:
            Ks = 0.0
        elif p_s < 0:
            return None, None, None
        else:
            Ks = -0.75 * math.log(1.0 - (4.0 * p_s) / 3.0)

        if p_n == 0.0:
            Ka = 0.0
        elif p_n < 0:
            return None, None, None
        else:
            Ka = -0.75 * math.log(1.0 - (4.0 * p_n) / 3.0)
    except ValueError:
        return None, None, None

    if Ks == 0:
        KaKs = None
    else:
        KaKs = Ka / Ks

    return Ka, Ks, KaKs


# ─────────────────────────────────────────────────────────────────────────────
# ALIGNMENT
# ─────────────────────────────────────────────────────────────────────────────

def align_cds_pairwise(seq1, seq2):
    """
    Align two CDS sequences using BioPython PairwiseAligner.
    Returns (aligned_seq1, aligned_seq2) trimmed to multiples of 3,
    or (None, None) if alignment fails quality checks.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1.upper(), seq2.upper())
    try:
        best = next(iter(alignments))
    except StopIteration:
        return None, None

    # Extract aligned sequences
    a1 = str(best[0])
    a2 = str(best[1])

    # Remove columns where BOTH have gaps
    cols1 = []
    cols2 = []
    for b1, b2 in zip(a1, a2):
        if b1 == '-' and b2 == '-':
            continue
        cols1.append(b1)
        cols2.append(b2)

    aln1 = ''.join(cols1)
    aln2 = ''.join(cols2)

    # Trim to multiple of 3
    trim = len(aln1) - (len(aln1) % 3)
    aln1 = aln1[:trim]
    aln2 = aln2[:trim]

    if len(aln1) < 120:
        return None, None

    # Check identity
    matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != '-')
    total = sum(1 for a, b in zip(aln1, aln2) if a != '-' or b != '-')
    if total == 0 or matches / total < 0.40:
        return None, None

    return aln1, aln2


def align_aa_pairwise(seq1, seq2):
    """
    Align two AA sequences and return % identity.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1.upper(), seq2.upper())
    try:
        best = next(iter(alignments))
    except StopIteration:
        return None

    a1 = str(best[0])
    a2 = str(best[1])

    matches = sum(1 for b1, b2 in zip(a1, a2) if b1 == b2 and b1 != '-')
    total = sum(1 for b1, b2 in zip(a1, a2) if not (b1 == '-' and b2 == '-'))
    if total == 0:
        return None
    return 100.0 * matches / total


# ─────────────────────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────

def load_fasta(path, strip_suffix=True):
    """
    Load FASTA file into dict {base_id: sequence}.
    If strip_suffix=True, strips .mRNA.N suffix from IDs.
    """
    seqs = {}
    current_id = None
    current_seq = []

    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                tx_id = line[1:].split()[0]
                if strip_suffix:
                    base_id = re.sub(r'\.mRNA\.\d+$', '', tx_id)
                else:
                    base_id = tx_id
                current_id = base_id
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(current_seq)
    return seqs


def load_gene_set(path):
    """Load a gzipped TSV of gene IDs (first column, skip header)."""
    genes = set()
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        next(f)  # skip header
        for line in f:
            gene_id = line.strip().split('\t')[0]
            if gene_id:
                genes.add(re.sub(r'\.mRNA\.\d+$', '', gene_id))
    return genes


def load_hog_table(path, spruce_col_name='Picea_abies', pine_col_name='Pinus_sylvestris'):
    """
    Load N10 HOG table. Returns dict {hog_id: {'spruce': [base_ids], 'pine': [base_ids]}}.
    """
    hogs = {}
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        print(f"  HOG table columns: {header}", file=sys.stderr)

        # Find column indices
        try:
            spruce_idx = header.index(spruce_col_name)
        except ValueError:
            raise ValueError(f"Column '{spruce_col_name}' not found. Available: {header}")
        try:
            pine_idx = header.index(pine_col_name)
        except ValueError:
            raise ValueError(f"Column '{pine_col_name}' not found. Available: {header}")

        print(f"  Spruce col index: {spruce_idx}, Pine col index: {pine_idx}", file=sys.stderr)

        for line in f:
            parts = line.rstrip('\n').split('\t')
            hog_id = parts[0]

            def parse_genes(idx):
                if idx < len(parts) and parts[idx].strip():
                    return [re.sub(r'\.mRNA\.\d+$', '', g.strip())
                            for g in parts[idx].split(',') if g.strip()]
                return []

            hogs[hog_id] = {
                'spruce': parse_genes(spruce_idx),
                'pine': parse_genes(pine_idx),
            }
    return hogs


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=== Ka/Ks Analysis for Spruce Paralog Pairs ===\n", flush=True)

    # ── File paths (env-configurable; see SOURCES.tsv for raw-input provenance) ──
    # DATA_ROOT holds the genome CDS/AA FASTAs and orthogroup tables (FigShare
    # doi:10.17044/scilifelab.28737623). DEP is the deposit root (this repo's parent).
    import os
    _PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    DATA = os.environ.get("DATA_ROOT", os.path.join(_PROJ, "genome_data"))
    DEP  = os.environ.get("SPRUCE_PINE_DEPOSIT", os.getcwd())   # repo root: Orthogroups/ + deposited intermediates
    SPRUCE_CDS_PATH = os.environ.get("SPRUCE_CDS",
                          f"{DATA}/sprucev2/Picab02_230926_at01_longest_no_TE_cds.fa")
    SPRUCE_AA_PATH  = os.environ.get("SPRUCE_AA",
                          f"{DATA}/sprucev2/Picab02_230926_at01_longest_no_TE_aa.fa")
    PINE_AA_PATH    = os.environ.get("PINE_AA",
                          f"{DATA}/pinev1/Pinsy01_240308_at01_longest_no_TE_aa.fa")
    HOG_PATH        = f"{DEP}/Orthogroups/Phylogenetic_Hierarchical_Orthogroups_N10.tsv.gz"
    SPRUCE_SD_PATH  = f"{DEP}/Orthogroups/spruce_all_in_duplicated_blocks_PGset.tsv.gz"
    PINE_SD_PATH    = f"{DEP}/Orthogroups/pine_all_in_duplicated_blocks_PGset.tsv.gz"
    KAKS_OUT        = os.environ.get("KAKS_OUT", f"{DEP}/kaks_results.tsv")
    AA_IDENT_OUT    = os.environ.get("AA_IDENT_OUT", f"{DEP}/interspecies_aa_identity.tsv")

    # ── Step 1: Load data ────────────────────────────────────────────────────
    print("Step 1: Loading data...", flush=True)

    print("  Loading spruce CDS...", flush=True)
    spruce_cds = load_fasta(SPRUCE_CDS_PATH)
    print(f"  Loaded {len(spruce_cds):,} spruce CDS sequences", flush=True)

    print("  Loading spruce AA...", flush=True)
    spruce_aa = load_fasta(SPRUCE_AA_PATH)
    print(f"  Loaded {len(spruce_aa):,} spruce AA sequences", flush=True)

    print("  Loading pine AA...", flush=True)
    pine_aa = load_fasta(PINE_AA_PATH)
    print(f"  Loaded {len(pine_aa):,} pine AA sequences", flush=True)

    print("  Loading spruce SD gene set...", flush=True)
    spruce_sd = load_gene_set(SPRUCE_SD_PATH)
    print(f"  Loaded {len(spruce_sd):,} spruce SD genes", flush=True)

    print("  Loading pine SD gene set...", flush=True)
    pine_sd = load_gene_set(PINE_SD_PATH)
    print(f"  Loaded {len(pine_sd):,} pine SD genes", flush=True)

    print("  Loading HOG table...", flush=True)
    hogs = load_hog_table(HOG_PATH)
    print(f"  Loaded {len(hogs):,} HOGs", flush=True)

    # ── Step 2: Classify HOGs ────────────────────────────────────────────────
    print("\nStep 2: Classifying HOGs...", flush=True)

    shared_SD = []
    spruce_only_SD = []
    single_copy_orthologs = []

    for hog_id, data in hogs.items():
        sg = set(data['spruce'])
        pg = set(data['pine'])
        sg_in_sd = sg & spruce_sd
        pg_in_sd = pg & pine_sd

        if sg_in_sd and pg_in_sd and len(sg) >= 2:
            shared_SD.append(hog_id)
        elif len(sg_in_sd) >= 2 and not pg_in_sd:
            spruce_only_SD.append(hog_id)
        elif len(sg) == 1 and len(pg) == 1:
            single_copy_orthologs.append(hog_id)

    print(f"  Shared SD HOGs:          {len(shared_SD):,}", flush=True)
    print(f"  Spruce-only SD HOGs:     {len(spruce_only_SD):,}", flush=True)
    print(f"  Single-copy orthologs:   {len(single_copy_orthologs):,}", flush=True)

    # ── Steps 3-5: Compute Ka/Ks ─────────────────────────────────────────────
    print("\nSteps 3-5: Computing Ka/Ks and AA identity...", flush=True)

    results = []  # (hog_id, gene1, gene2, category, Ka, Ks, KaKs)

    def process_kaks_hogs(hog_list, category):
        done = 0
        skipped = 0
        for hog_id in hog_list:
            data = hogs[hog_id]
            sg = data['spruce']
            sg_in_sd = [g for g in sg if g in spruce_sd]

            if len(sg_in_sd) < 2:
                # Fall back: use any two spruce genes from the HOG
                if len(sg) >= 2:
                    pair = sg[:2]
                else:
                    continue
            else:
                pair = sg_in_sd[:2]

            g1, g2 = pair[0], pair[1]

            if g1 not in spruce_cds or g2 not in spruce_cds:
                skipped += 1
                continue

            res = paml_yn00_dnds(spruce_cds[g1], spruce_cds[g2])   # PAML yn00 (dN, dS); same as cross-species
            if res is None:
                skipped += 1
                continue
            Ka, Ks = res
            KaKs = (Ka / Ks) if (Ks and Ks > 0) else None
            results.append((hog_id, g1, g2, category, Ka, Ks, KaKs))
            done += 1

            if done % 100 == 0:
                print(f"    {category}: {done} pairs processed...", flush=True)

        print(f"  {category}: {done} pairs computed, {skipped} skipped", flush=True)

    process_kaks_hogs(shared_SD, 'shared_SD')
    process_kaks_hogs(spruce_only_SD, 'spruce_only_SD')

    # ── Calibration: interspecies AA identity ────────────────────────────────
    print("\n  Computing interspecies AA identity (calibration)...", flush=True)
    aa_results = []  # (hog_id, pa_gene, ps_gene, pct_identity)

    done_calib = 0
    skipped_calib = 0
    for hog_id in single_copy_orthologs:
        data = hogs[hog_id]
        pa_gene = data['spruce'][0] if data['spruce'] else None
        ps_gene = data['pine'][0] if data['pine'] else None

        if pa_gene is None or ps_gene is None:
            continue
        if pa_gene not in spruce_aa or ps_gene not in pine_aa:
            skipped_calib += 1
            continue

        pct_id = align_aa_pairwise(spruce_aa[pa_gene], pine_aa[ps_gene])
        if pct_id is not None:
            aa_results.append((hog_id, pa_gene, ps_gene, pct_id))
            done_calib += 1
        else:
            skipped_calib += 1

        if done_calib % 50 == 0 and done_calib > 0:
            print(f"    Calibration: {done_calib} pairs...", flush=True)

    print(f"  Calibration: {done_calib} pairs computed, {skipped_calib} skipped", flush=True)

    # ── Step 6: Output ───────────────────────────────────────────────────────
    print("\nStep 6: Writing output files...", flush=True)

    with open(KAKS_OUT, 'w') as f:
        f.write("hog_id\tgene1\tgene2\tcategory\tKa\tKs\tKaKs\n")
        for row in results:
            hog_id, g1, g2, cat, Ka, Ks, KaKs = row
            Ka_str  = f"{Ka:.6f}"  if Ka  is not None else "NA"
            Ks_str  = f"{Ks:.6f}"  if Ks  is not None else "NA"
            KaKs_str = f"{KaKs:.6f}" if KaKs is not None else "NA"
            f.write(f"{hog_id}\t{g1}\t{g2}\t{cat}\t{Ka_str}\t{Ks_str}\t{KaKs_str}\n")

    print(f"  Written: {KAKS_OUT}", flush=True)

    with open(AA_IDENT_OUT, 'w') as f:
        f.write("hog_id\tpa_gene\tps_gene\tpct_identity\n")
        for row in aa_results:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]:.2f}\n")

    print(f"  Written: {AA_IDENT_OUT}", flush=True)

    # ── Summary statistics ───────────────────────────────────────────────────
    def summarise(rows, category):
        subset = [r for r in rows if r[3] == category]
        valid_Ka   = [r[4] for r in subset if r[4] is not None]
        valid_Ks   = [r[5] for r in subset if r[5] is not None]
        valid_KaKs = [r[6] for r in subset if r[6] is not None]

        def mean(lst):
            return sum(lst) / len(lst) if lst else float('nan')
        def median(lst):
            if not lst: return float('nan')
            s = sorted(lst)
            n = len(s)
            return (s[n//2] + s[(n-1)//2]) / 2

        print(f"\n  [{category}]")
        print(f"    Pairs with valid Ka/Ks: {len(valid_KaKs)} / {len(subset)}")
        if valid_Ka:
            print(f"    Ka  — mean={mean(valid_Ka):.4f}, median={median(valid_Ka):.4f}")
        if valid_Ks:
            print(f"    Ks  — mean={mean(valid_Ks):.4f}, median={median(valid_Ks):.4f}")
        if valid_KaKs:
            print(f"    Ka/Ks — mean={mean(valid_KaKs):.4f}, median={median(valid_KaKs):.4f}")
            n_pos = sum(1 for v in valid_KaKs if v > 1.0)
            n_neg = sum(1 for v in valid_KaKs if v < 1.0)
            print(f"    Ka/Ks > 1 (positive selection): {n_pos}")
            print(f"    Ka/Ks < 1 (purifying selection): {n_neg}")

    print("\n=== SUMMARY STATISTICS ===")
    summarise(results, 'shared_SD')
    summarise(results, 'spruce_only_SD')

    if aa_results:
        pct_ids = [r[3] for r in aa_results]
        mean_id = sum(pct_ids) / len(pct_ids)
        s = sorted(pct_ids)
        n = len(s)
        med_id = (s[n//2] + s[(n-1)//2]) / 2
        print(f"\n  [Single-copy orthologs — AA identity calibration]")
        print(f"    Pairs:  {len(aa_results)}")
        print(f"    Mean AA identity:   {mean_id:.1f}%")
        print(f"    Median AA identity: {med_id:.1f}%")

    print("\n=== DONE ===", flush=True)


if __name__ == '__main__':
    main()
