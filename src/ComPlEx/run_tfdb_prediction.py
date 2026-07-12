#!/usr/bin/env python3
"""
PlantTFDB transcription-factor prediction pipeline (provenance script)
=====================================================================
Submits gymnosperm proteome FASTA files to the PlantTFDB prediction service
(planttfdb.gao-lab.org) for transcription-factor annotation, then concatenates
the per-chunk results into one TSV per species.

PROVENANCE
----------
This is the script that produced the transcription-factor annotation used in the
manuscript. Run on the Picea abies and Pinus sylvestris longest-isoform,
TE-filtered protein sets, it generated:

    data/annotation/Picab02_230926_at01_longest_no_TE_aa_TF_predictions.tsv
    data/annotation/Pinsy01_240308_at01_longest_no_TE_aa_TF_predictions.tsv

which are the inputs to tf_enrichment_by_category.py, plot_tf_figures.R and
figure3_coexpressolog_tf.R. The two committed files are byte-for-byte the output of
this script.

NOTE ON REPRODUCIBILITY
-----------------------
This step queries a live external web service (PlantTFDB), so it is provenance
rather than a deterministic, offline-reproducible pipeline stage: re-running it
depends on the remote server and may change if PlantTFDB updates its models. For
that reason the prediction TSVs are committed as static inputs; this script
documents exactly how they were produced. The only TF-annotation tool used was
PlantTFDB (Jin et al. 2017); no other predictor was involved.

The input proteomes are the gene-model protein sets (longest isoform per gene,
transposable-element-derived models removed) and are not redistributed here.

Usage:
    python3 run_tfdb_prediction.py        # after setting INPUT_DIR below

Requirements:
    pip install requests beautifulsoup4 lxml
"""

import requests
import urllib3
import re
import io
import os
import time
import sys
from bs4 import BeautifulSoup
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
# INPUT_DIR must contain the proteome .fa files (not shipped in the repo).
# OUTPUT_DIR defaults to the repo annotation directory where the committed
# prediction TSVs live; existing outputs are skipped, so it will not overwrite.
INPUT_DIR  = Path(os.environ.get("TFDB_INPUT_DIR", "data/annotation/proteomes"))
OUTPUT_DIR = Path(os.environ.get("TFDB_OUTPUT_DIR", "data/annotation"))

CHUNK_SIZE  = 4000   # sequences per submission (keeps well within server timeout)
MAX_RETRIES = 4      # retries per chunk on failure
RETRY_DELAY = 20     # seconds between retries
BASE_URL    = 'https://planttfdb.gao-lab.org'

# Suppress SSL warnings if certificate validation is turned off
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# ── Core functions ─────────────────────────────────────────────────────────────

def count_sequences(fasta_path: Path) -> int:
    with open(fasta_path, encoding='utf-8', errors='replace') as f:
        return sum(1 for line in f if line.startswith('>'))


def get_chunks(fasta_path: Path, chunk_size: int):
    """Yield FASTA chunks as UTF-8 bytes, each containing at most chunk_size sequences."""
    current = []
    chunk   = []
    with open(fasta_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if current:
                    chunk.append(''.join(current))
                    if len(chunk) == chunk_size:
                        yield ''.join(chunk).encode('utf-8')
                        chunk = []
                current = [line]
            else:
                current.append(line)
        if current:
            chunk.append(''.join(current))
        if chunk:
            yield ''.join(chunk).encode('utf-8')


def get_fresh_userid(session: requests.Session) -> str:
    r = session.get(f'{BASE_URL}/prediction.php', timeout=30, verify=False)
    r.raise_for_status()
    m = re.search(r'prediction_result\.php\?userid=(\w+)', r.text)
    if not m:
        raise ValueError("Could not extract userid from prediction page.")
    return m.group(1)


def submit_chunk(session: requests.Session, chunk_bytes: bytes,
                 species: str, chunk_idx: int) -> str:
    """Submit one FASTA chunk; return the TSV result text (empty string if no TFs found)."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            userid = get_fresh_userid(session)
            r = session.post(
                f'{BASE_URL}/prediction_result.php?userid={userid}',
                files={'input_file': ('chunk.fa', io.BytesIO(chunk_bytes), 'text/plain')},
                data={'best1_ath': 'on', 'input_seq': ''},
                timeout=90,
                verify=False
            )
            r.raise_for_status()

            soup = BeautifulSoup(r.text, 'lxml')

            # Locate the TF-results download link
            tf_link = None
            for a in soup.find_all('a', href=True):
                if 'TF_and_best1' in a['href']:
                    tf_link = BASE_URL + '/' + a['href'].lstrip('/')
                    break

            if tf_link is None:
                return ''   # valid: no TFs in this chunk

            r2 = session.get(tf_link, timeout=30, verify=False)
            r2.raise_for_status()
            return r2.text

        except Exception as e:
            print(f"\n    ⚠  Attempt {attempt}/{MAX_RETRIES} failed: {e}")
            if attempt < MAX_RETRIES:
                print(f"    Retrying in {RETRY_DELAY}s...", flush=True)
                time.sleep(RETRY_DELAY)

    raise RuntimeError(
        f"All {MAX_RETRIES} attempts failed for {species} chunk {chunk_idx}."
    )


def process_species(fasta_path: Path) -> None:
    species   = fasta_path.stem.replace('.aa', '')
    out_path  = OUTPUT_DIR / f'{species}_TF_predictions.tsv'

    if out_path.exists():
        print(f'\n  ✔  Skipping {species} — output already exists.')
        return

    total_seqs = count_sequences(fasta_path)
    n_chunks   = (total_seqs + CHUNK_SIZE - 1) // CHUNK_SIZE

    print(f'\n{"="*64}')
    print(f'  Species : {species}')
    print(f'  Seqs    : {total_seqs:,}   Chunks: {n_chunks} × {CHUNK_SIZE}')
    print(f'{"="*64}')

    session     = requests.Session()
    all_results = []
    header      = 'Gene_ID\tFamily\tBest_hit_Ath\tE_value\tDescription\n'
    t0_species  = time.time()

    for idx, chunk_bytes in enumerate(get_chunks(fasta_path, CHUNK_SIZE), 1):
        print(f'  Chunk {idx}/{n_chunks}  ({len(chunk_bytes):,} bytes) ... ',
              end='', flush=True)
        t0 = time.time()
        result = submit_chunk(session, chunk_bytes, species, idx)
        elapsed = time.time() - t0
        n_tfs = len([l for l in result.strip().splitlines() if l.strip()])
        print(f'✓  {n_tfs} TFs  [{elapsed:.1f}s]')
        if result.strip():
            all_results.append(result.strip())

    with open(out_path, 'w') as fh:
        fh.write(header)
        if all_results:
            fh.write('\n'.join(all_results) + '\n')

    total_tfs     = sum(len([l for l in r.strip().splitlines() if l]) for r in all_results)
    species_mins  = (time.time() - t0_species) / 60
    print(f'\n  → {total_tfs:,} TFs saved to {out_path.name}  [{species_mins:.1f} min]')


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(
        f for f in INPUT_DIR.glob('*.fa')
        if not f.name.startswith('._')
    )

    if not fasta_files:
        print(f'\nNo .fa files found in {INPUT_DIR}')
        print('Set TFDB_INPUT_DIR (or edit INPUT_DIR) to the proteome directory.')
        sys.exit(1)

    print(f'\nPlantTFDB Prediction Pipeline')
    print(f'Input  : {INPUT_DIR}')
    print(f'Output : {OUTPUT_DIR}')
    print(f'Species to process ({len(fasta_files)}):')
    for f in fasta_files:
        print(f'  - {f.stem}')

    t0_total = time.time()
    for fasta in fasta_files:
        process_species(fasta)

    total_mins = (time.time() - t0_total) / 60
    print(f'\n{"="*64}')
    print(f'All species complete in {total_mins:.1f} minutes.')
    print(f'\nSummary:')
    for tsv in sorted(OUTPUT_DIR.glob('*_TF_predictions.tsv')):
        with open(tsv) as fh:
            n = sum(1 for line in fh) - 1  # subtract header row
        print(f'  {tsv.stem}: {n:,} TFs')
    print(f'\nResults saved to: {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
