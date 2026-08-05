#!/usr/bin/env bash
###############################################################################
# reproduce_paper.sh — single driver that regenerates every figure, table and
# reported number in van Zalen et al. from the committed scripts + the deposited
# data (FigShare + ENA genomes). Replaces run_pipeline.sh.
#
# HOW TO REPRODUCE (see README.md + SOURCES.tsv for the file list and DOIs):
#   1. Clone this repository (the scripts).
#   2. Download the FigShare deposits listed in SOURCES.tsv and extract them in
#      place, so the repo contains data/expression, data/dds, data/annotation,
#      data/popgen … , data/annotation/gene_annotation.tsv.gz and results/ComPlEx/ (the deposited
#      ComPlEx co-expression network output).
#   3. Place the genome CDS/protein/repeat files under $DATA_ROOT (default
#      ./genome_data) as sprucev2/ and pinev1/ (ENA assemblies; see SOURCES.tsv).
#   4. Run:  DATA_ROOT=/path/to/genome_data ./reproduce_paper.sh
#
# ENVIRONMENT:
#   DATA_ROOT  genome CDS/AA/repeat/TSS dir (default ./genome_data)
#   R          Rscript                       (default: Rscript)
#   PY_SCI     python with scipy+pandas       (default: python3)
#   PY_BIO     python with Biopython+scipy    (default: python3; cross-species Ka/Ks)
#
# Stages run in dependency order; each is logged to reproduce_paper.log.
# STOP_ON_FAIL=1 aborts on the first failure (default: continue, report at end).
###############################################################################
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
export SPRUCE_PINE_DEPOSIT="$HERE"
export DATA_ROOT="${DATA_ROOT:-$HERE/genome_data}"
R="${R:-Rscript}"; PY_SCI="${PY_SCI:-python3}"; PY_BIO="${PY_BIO:-python3}"
export PYTHONUNBUFFERED=1
STOP_ON_FAIL="${STOP_ON_FAIL:-0}"
LOG="$HERE/reproduce_paper.log"; : > "$LOG"
ok=0; fail=0; failed_stages=""

# START_AT=<stage-name>: skip every stage before the named one (they still print as SKIP). For
# re-running only a downstream segment after an upstream change, when the upstream outputs are known
# current. With no START_AT the whole pipeline runs. Does NOT relax the no-cache rule for stages that
# do run — each still executes its producer.
START_AT="${START_AT:-}"; _reached=1; [ -n "$START_AT" ] && _reached=0
stage(){ local name="$1"; shift; local t0=$SECONDS
  if [ "$_reached" = 0 ]; then
    if [ "$name" = "$START_AT" ]; then _reached=1; else
      printf '%s%-34s SKIP (before START_AT=%s)\n' '-- ' "$name" "$START_AT"; return 0; fi
  fi
  printf '[%s] %s%-34s ' "$(date '+%Y-%m-%d %H:%M:%S')" '-- ' "$name"
  if "$@" >>"$LOG" 2>&1; then printf 'OK   (%ds)\n' $((SECONDS-t0)); ok=$((ok+1))
  else printf 'FAIL (%ds)  see reproduce_paper.log\n' $((SECONDS-t0))
    fail=$((fail+1)); failed_stages="$failed_stages $name"
    [ "$STOP_ON_FAIL" = 1 ] && { echo "STOP_ON_FAIL set — aborting."; exit 1; }
  fi; }
# NOTE: every stage runs its producer on every invocation — there is no output-cache skip. The
# genome-input producers (Ka/Ks, dN/dS, promoter FASTA/TE, FIMO) therefore REQUIRE the genome files
# under DATA_ROOT and a FIMO binary; without them those stages fail loudly (as they should) rather than
# being silently skipped. The deposited copies of their outputs are a cross-check target, not a trusted
# input. (A previous stage_cached skip-when-present wrapper masked the cold_root orphan and a
# hardcoded-path bug from every "clean-room" run; it has been removed.)
rmd(){ "$R" -e "rmarkdown::render('$1', quiet=TRUE)"; }

# Output directory skeleton — created up front so a from-empty run does not fail on a missing
# target dir (e.g. the DE renders save() into data/DEG_lists). Producers that dir.create their
# own subdir still work; this just guarantees the common output roots exist.
mkdir -p data/DEG_lists \
         results/ComPlEx/RData \
         results/integration/centrality results/integration/chs3 \
         results/integration/fig4 results/integration/figures results/integration/go \
         results/final_figures

echo "=== reproduce_paper.sh  (DATA_ROOT=$DATA_ROOT) ==="

# ── DIFFERENTIAL EXPRESSION — per condition DE tables (feed Fig 1c, Fig 2, S3) ──
for f in SpruceColdNeedle SpruceColdRoot SpruceDroughtNeedle SpruceDroughtRoot \
         PineColdNeedle PineColdRoot PineDroughtNeedle PineDroughtRoot; do
  stage "DE_$f"                        rmd "src/BioQA_DE/DE_$f.Rmd"
done
stage "build_deg_overlap_counts"       "$R"      src/BioQA_DE/build_deg_overlap_counts.R
stage "ComPlExDataPrep"                "$R"      src/ComPlEx/ComPlExDataPrep.R   # per-condition VST (Fig 1 PCA + ComPlEx inputs)

# ── CO-EXPRESSION NETWORKS — rebuilt from the VST matrices every run (NOT a deposited input). ─
# The four ComPlEx co-expressolog networks are a core analysis stage and are regenerated here from the
# count-derived VST matrices; cliques_step1/1b then build weighted_gene_pairs (the conservation-breadth
# table). Requires the ComPlEx_python checkout (COMPLEX_PY_REPO; https://github.com/natstreet/ComPlEx_python).
# Producers intentionally NOT staged here (every other src/ producer IS a stage above/below):
#   src/ComPlEx/ComPlEx.R           — R REFERENCE network implementation; the validated Python port
#                                     (run_complex_networks.sh) builds the networks the pipeline uses,
#                                     and no staged script consumes ComPlEx.R's outputs.
#   src/ComPlEx/run_tfdb_prediction.py — provenance one-off: submits proteomes to the LIVE PlantTFDB
#                                     web service (non-deterministic/external); its TF annotation is
#                                     deposited and read by figure3. Not part of automated reproduction.
#   src/ComPlEx/kaks_yn00.py        — runs via run_cross_species_dnds.sh (the cross_species_dnds stage).
#   src/lib/*.R                     — sourced helpers (go_enrichment.R, fig_palette.R), not standalone.
stage "build_complex_networks"         bash      src/ComPlEx/run_complex_networks.sh   # 4 networks via the Python port
stage "cliques_step1"                  "$R"      src/ComPlEx/cliques_step1.R           # networks -> co_expressologs
stage "cliques_step1b"                 "$R"      src/ComPlEx/cliques_step1b.R          # co_expressologs -> weighted_gene_pairs

# ── EVOLUTIONARY BACKBONES — within/cross-species Ka/Ks, expression divergence ─
stage "kaks_within_species_sd"         "$PY_BIO" src/ComPlEx/kaks_within_species_sd.py
stage "cross_species_ks"               "$PY_BIO" src/ComPlEx/cross_species_ks.py
stage "sd_pair_expression_divergence"  "$PY_SCI" src/ComPlEx/sd_pair_expression_divergence.py

# ── CO-EXPRESSION CONSERVATION — clique categories from the freshly rebuilt weighted_gene_pairs ─
stage "cliques_step2"                  "$R"      src/ComPlEx/cliques_step2.R

# ── PROMOTER TF / TE — FASTA extraction, TE presence, CHS3 motifs ──────────────
stage "extract_upstream_seqs"          "$PY_SCI" src/ComPlEx/extract_upstream_seqs.py   # needs softmasked genome assemblies
stage "pa_ps_promoter_te_presence"     "$PY_SCI" src/ComPlEx/pa_ps_promoter_te_presence.py
stage "chs3_promoter_motifs"           "$PY_SCI" src/ComPlEx/chs3_promoter_motifs.py   # Fig 6 motif stat (needs FIMO + upstream FASTA)


# ── INTEGRATION HUB — backbone_1to1 + category gene lists + TE covariate ───────
stage "integration_analysis"           "$R"      src/ComPlEx/integration_analysis.R
stage "te_promoter_family_analysis"    "$PY_SCI" src/ComPlEx/te_promoter_family_analysis.py  # Fig 7 zte (needs repeat GFFs)

# ── CROSS-SPECIES dN/dS (YN00) — needs backbone_1to1 ──────────────────────────
PYTHON="$PY_BIO" stage "cross_species_dnds" bash src/ComPlEx/run_cross_species_dnds.sh

# ── CATEGORY UNIVERSE, ENRICHMENT, dN/dS-vs-breadth, popgen, PAV ───────────────
stage "build_category_universe"        "$PY_SCI" src/ComPlEx/build_category_universe.py
stage "build_network_degree"           "$PY_SCI" src/ComPlEx/build_network_degree.py
stage "tf_enrichment_by_category"      "$PY_SCI" src/ComPlEx/tf_enrichment_by_category.py     # Fig 3b context
stage "pnps_confound_analysis"         "$PY_SCI" src/ComPlEx/pnps_confound_analysis.py
stage "single_copy_dnds_breadth"       "$PY_SCI" src/ComPlEx/single_copy_dnds_breadth.py      # Fig 5a
stage "degree_dnds_correlation"        "$R"      src/ComPlEx/degree_dnds_correlation.R        # Fig 7: hubness vs purifying selection
stage "pine_axis_replication"          "$R"      src/ComPlEx/pine_axis_replication.R          # reciprocal Scots pine: degree-dN/dS + Model A/B
stage "threshold_sensitivity"          "$PY_SCI" src/ComPlEx/threshold_sensitivity.py         # Supp Table S1
stage "not_coex_de_analysis"           "$R"      src/ComPlEx/not_coex_de_analysis.R           # Fig S mechanism inputs (spruce)
stage "not_coex_de_pine"               "$R"      src/ComPlEx/not_coex_de_pine.R               # reciprocal Scots pine not_coex confirmation
stage "wood_clique_and_conditioned_pav" "$PY_SCI" src/ComPlEx/wood_clique_and_conditioned_pav.py
stage "stress_selection_enrichment"    "$R"      src/ComPlEx/stress_selection_enrichment.R    # selection enrichment (supp fig)
stage "enrichment_tests_extended"      "$R"      src/ComPlEx/enrichment_tests_extended.R
stage "pine_sd_classify"               "$PY_BIO" src/ComPlEx/pine_sd_classify.py               # pine SD gene classes (reciprocal axis)
stage "pine_sd_category_enrichment"    "$R"      src/ComPlEx/pine_sd_category_enrichment.R     # pine SD x category enrichment (Fig 5b pine)
stage "pav_category_bh"                "$R"      src/ComPlEx/pav_category_bh.R                # PAV BH within the 5-category family [not_coex p_bh]
stage "gwas_deg_overlap"               "$R"      src/ComPlEx/gwas_deg_overlap.R
stage "orthogroup_direction_concordance" "$R"    src/ComPlEx/orthogroup_direction_concordance.R  # Fig 2 directional concordance [27]
stage "single_species_deg_fraction"    "$R"      src/ComPlEx/single_species_deg_fraction.R       # single-species DEG fraction 8.5/10.9

# ── CONSERVED DYNAMICS (Fig 4) + GO ENRICHMENT (single shared method) ─────────
stage "conserved_coexpression_dynamics" "$R"     src/ComPlEx/conserved_coexpression_dynamics.R
stage "coex_category_go"               "$R"      src/ComPlEx/coex_category_go.R               # Results-text GO
stage "overlap_set_go"                 "$R"      src/BioQA_DE/overlap_set_go.R                # Results-text GO
stage "subset_go"                      "$R"      src/ComPlEx/subset_go.R                      # Results-text GO
stage "deg_directional_go"             "$R"      src/ComPlEx/deg_directional_go.R             # Supp Table S2
stage "table1_2_conserved_dynamics_go" "$R"      src/tables/table1_2_conserved_dynamics_go.R # Table 1 (drought) + Table 2 (cold)

# ── INTEGRATIVE MODEL (feeds Supp Fig S5 regression) ──────────────────────────────────────────────────────────────────────────
stage "integrative_conservation_model" "$R"      src/ComPlEx/integrative_conservation_model.R

# ── FIGURES — component producers + final assembly ────────────────────────────
stage "figure1ab_pca"                  "$R"      src/figures/figure1ab_pca.R                 # Fig 1a/b
stage "figure1c_deg_counts"            "$R"      src/figures/figure1c_deg_counts.R           # Fig 1c
stage "figure2ab_deg_overlap"          "$R"      src/figures/figure2ab_deg_overlap.R         # Fig 2a/b
stage "figure2c_orthogroup_overlap"    "$R"      src/figures/figure2c_orthogroup_overlap.R   # Fig 2c
stage "figure3_coexpressolog_tf"       "$R"      src/figures/figure3_coexpressolog_tf.R      # Fig 3 (co-expressolog UpSet + TF)
stage "figure4_conserved_dynamics"     "$R"      src/figures/figure4_conserved_dynamics.R    # Fig 4
stage "figure5_dnds_evolution"         "$R"      src/figures/figure5_dnds_evolution.R        # Fig 5 (dN/dS vs conservation)
stage "figure6_chs3_case_study"        "$R"      src/figures/figure6_chs3_case_study.R       # Fig 6 (a/b/c)
stage "figureS5_integrative_model"      "$R"      src/figures/figureS5_integrative_model.R     # Supp Fig S5 (regression, demoted from main Fig 7)
stage "figureS1_physiology"            "$R"      src/figures/figureS1_physiology.R           # Supp Fig S1
stage "figureS2_timing"                "$R"      src/figures/figureS2_timing_deg_fraction.R  # Supp Fig S2 + Supp Table S3 (timing)
stage "figureS3a_notcoex_mechanism"    "$R"      src/figures/figureS3a_notcoex_mechanism.R   # Supp Fig S3a
stage "figureS3b_pav_enrichment"       "$R"      src/figures/figureS3b_pav_enrichment.R      # Supp Fig S3b
stage "figureS4_chs_phylogeny"         "$R"      src/figures/figureS4_chs_phylogeny.R        # Supp Fig S4
stage "assemble_figures"               "$R"      src/figures/assemble_figures.R              # composes final Figure*/FigureS* PDFs

# ── TABLES — main-text (Table 1/2 built above by table1_2_...) + supplementary ─
stage "table_s1_go"                    "$R"      src/tables/table_s1_go_enrichment.R                # Supp Table S1 (GO of DEGs)
stage "table_s2_de"                    "$R"      src/tables/table_s2_de_results.R                   # Supp Table S2 (DE lists; S3 timing = figureS2 stage)
stage "table_s4_conserved_go"          "$R"      src/tables/table_s4_conserved_coexpressolog_go.R   # Supp Table S4 (conserved co-expressolog GO)
stage "table_s5_category_go"           "$R"      src/tables/table_s5_coexpression_category_go.R     # Supp Table S5 (co-expression category GO)
stage "table_s6_dnds"                  "$R"      src/tables/table_s6_dnds_yn00.R                    # Supp Table S6 (dN/dS YN00)
stage "table_s7_threshold"             "$R"      src/tables/table_s7_threshold_sensitivity.R        # Supp Table S7 (threshold sensitivity)
stage "table_s8_gene_axis"             "$R"      src/tables/table_s8_gene_axis_classification.R     # Supp Table S8 (per-gene conservation-divergence axis)

echo "=========================================================================="
echo "reproduce_paper.sh finished: $ok OK, $fail FAILED.${failed_stages:+  failed:$failed_stages}"
echo "Full log: reproduce_paper.log"

# NOTE: Figure 8 is a hand-drawn schematic (no producer); its underlying numbers come
# from selection_category_enrichment.tsv (stress_selection_enrichment). Table 1/Table 2
# are both produced by conserved_dynamics_go (drought / cold rows). The ComPlEx network
# layer (the four networks + co_expressologs + weighted_gene_pairs) is now REGENERATED here
# from the VST matrices (build_complex_networks -> cliques_step1 -> cliques_step1b), not loaded
# from the deposit — the deposited copies are a cross-check target, not a trusted input.
[ "$fail" -eq 0 ]
