## Lens status
- **probe:** applied — two passes (blast-radius Mash on 5 cohorts; reshaped hard-dedup → disclosure layer).
- **brainstorm:** applied — Claude-only (codex was hook-blocked); surfaced C1-C3 + M1-M4, all folded into the Steps below.
- **review:** not applied.

## Problem Statement
The provenance-disjoint report card reports raw-isolate sens/spec and a tier string implying biological independence the data lacks: every SCORED R class is clonally dominated (<=17 effective lineages at Mash 0.001; Klebsiella cipro R = 30 BioProjects but 3 lineages). Hard-dedup to >=20/class would empty the card. Goal: a lineage-disclosure layer — per-cell raw_N + effective_lineage_N@{0.001,0.005} + cluster-weighted sens/spec with Wilson CI + effective-N + a graded lineage annotation + a reframed honest tier — keeping every cell visible with honest dual N, never hard-demoting. Whole-cohort greedy-representative clustering (chaining-resistant); discordant clones surfaced, not voted away.

## Codebase Context
- scripts/provenance_disjoint_validate.py — select_disjoint = pools[lab][:per_class] raw order, no dedup; genomes download after; artifact metrics from independent_cohort_validate._conf; carries registry_organism/drug. Per-strain pred = call_resistance(main.tsv, drug, organism=).
- scripts/build_validation_report_card.py — _key(org,drug)=lower; load_scored/census/registry + 7-state classify() + main(); PROV_TIER string to reframe; emits decoder_validation_report_card.{md,json}.
- dna_decode/eval/phylogeny.compute_mash_distances(strain_genomes, use_docker=True) -> DistanceMatrix.
- scripts/mash_cluster_n147.py::{score_threshold,pick_best_threshold,per_clade_label_balance} — clustering logic to lift (verify its linkage method in Step 1; the new module uses greedy-representative, not single-linkage).
- dna_decode/eval/amr_rules.call_resistance — per-strain recompute.
- dna_decode/models/cache.verify_complete — the project's completeness-gate pattern to mirror for genome integrity.
- data/raw/<slug>_provdisjoint_<drug>/{selected.tsv, amrfinder_runs/<acc>/main.tsv, refseq/<acc>/genome.fna} — 5 cohorts genome-complete; kleb gent (3/60) + tet (33/60) need re-fetch.

### Reusable-Code Survey
- dna_decode/eval/phylogeny.compute_mash_distances — distance engine (Docker-proven).
- scripts/mash_cluster_n147.py clustering — lift, but swap single-linkage -> greedy-representative (C2).
- scripts/independent_cohort_validate._conf — raw sens/spec/confusion; the weighted metric mirrors its shape + adds Wilson CI.
- dna_decode/models/cache.verify_complete — completeness-gate pattern for genome integrity (M1).
- dna_decode/eval/amr_rules.call_resistance — per-strain pred recompute (M4).
- dna_decode/data/refseq.py — genome fetch-by-accession (confirm entrypoint in Step 1).
- None — searched dna_decode/eval/, scripts/, graphify-out/ (absent) for an existing cluster-weighted-metric / Wilson-CI / greedy-dedup helper; none exists.

## Pre-Change Baseline
- Report card (committed): 7 SCORED cells, raw-isolate sens/spec only; no lineage columns; PROV_TIER implies independence.
- Measured inflation (5/7 cohorts): every R class <20 effective lineages @0.001.
- Tests: tests/test_build_validation_report_card.py (12) + test_amr_rules* green.
- Toolkit validation of THIS plan: validated at save time (tools restored after the sentinel-guard hook ACL fix).

## Verification Signal
- wiki/provdisjoint_lineage_metrics.json emitted for all 7 SCORED cells: raw_N, effective_lineage_N@{0.001,0.005}, cluster_weighted_{sens,spec} each with Wilson CI + effective-N, n_discordant_lineages, partial/n_genomes_missing.
- Report card renders new columns; no cluster-weighted point estimate appears without its CI+N (emitter asserts this).
- Recomputed raw sens/spec reconciles with each cell's existing artifact metrics (M4 assert) before the weighted number is trusted.
- Tier string reframed; no SCORED cell removed or demoted; cells with missing genomes show "lineage: incomplete (k/N)".
- New tests green; existing report-card + amr_rules tests unchanged.

## Implementation Steps

### Step 1: Clonality clustering, discordance, and cluster-weighted metric module
Files: dna_decode/eval/clonality.py
Depends on: none

**What changes:**
- New module. greedy_representative_clusters(genomes: dict[str,Path], threshold: float, use_docker=True) -> dict[str,int] — whole-cohort Mash via compute_mash_distances, then greedy-representative dedup (pick a strain, drop all within threshold, repeat — chaining-resistant, C2). Returns strain_id->cluster_id. Sort strain ids before picking for determinism.
- cluster_class(member_labels) -> "R"|"S"|"DISCORDANT" — a cluster mixing R+S labels is DISCORDANT (C1: not majority-voted).
- cluster_weighted_confusion(preds, labels, clusters) -> dict — one vote per same-label cluster (pred=majority, label=cluster_class); DISCORDANT clusters excluded from sens/spec, counted as n_discordant_lineages; returns {tp,fp,tn,fn,sens,spec,n_clusters_R,n_clusters_S,n_discordant} (C1/C3).
- wilson_ci(k:int, n:int) -> (lo,hi) — Wilson interval (C3).
- effective_lineage_n(clusters, labels, cls) -> int — count same-label clusters in class.
- Pure-logic (clustering math, aggregation, CI) split from the Docker Mash call for offline unit tests.

**Test strategy:**
- Unit: synthetic distance matrices -> greedy clusters resist a chain that single-linkage would merge; cluster_class on a mixed cluster -> DISCORDANT; cluster_weighted_confusion excludes discordant + collapses a clone to one vote; wilson_ci(8,8) ~ [0.63,1.0].

### Step 2: Per-cohort lineage-metrics recompute script
Files: scripts/compute_lineage_metrics.py, dna_decode/data/cell_key.py
Depends on: Step 1

**What changes:**
- New dna_decode/data/cell_key.py::canonical_cell_key(organism, drug) — the SINGLE canonical join key, imported by this script AND the report card (M2).
- New script. For each data/raw/*_provdisjoint_* cohort: ensure genomes with a completeness gate (exists + non-empty + valid FASTA header + sane contig count, mirroring cache.verify_complete); fetch missing via dna_decode/data/refseq (restartable, skip-complete); on fetch failure mark cohort partial + record n_genomes_missing and do not emit a lineage tier for it (M1).
- Read the cohort's original artifact registry_organism/drug; recompute per-strain call_resistance with exactly those args and assert recomputed raw sens/spec == artifact metrics before proceeding (M4).
- greedy_representative_clusters whole-cohort @ {0.001,0.005}; compute raw (_conf) + cluster_weighted_confusion (with wilson_ci); derive a graded lineage_N annotation (buckets, not a binary >=20 tier — M3).
- Persist wiki/provdisjoint_lineage_metrics.json (provdisjoint-lineage-metrics-v1), keyed by canonical_cell_key, checkpointed per-cohort (Docker-wedge-safe), idempotent upsert.

**Test strategy:**
- Unit: tiny on-disk fixture cohort (no network) -> sidecar schema; raw-vs-weighted divergence on a clustered fixture; M4 reconcile-assert fires on a mismatched-args mock; idempotent upsert; partial-cohort path sets partial+n_genomes_missing and emits no tier.

### Step 3: Report card consumes lineage metrics + CI rendering + reframed tier
Files: scripts/build_validation_report_card.py
Depends on: Step 2

**What changes:**
- Replace _key with imported canonical_cell_key (M2); load_lineage_metrics() reads the sidecar; main() adds per SCORED cell: raw_N, effective_lineage_N@{0.001,0.005}, cluster_weighted sens/spec rendered inline with Wilson CI + effective-N (C3 — an emitter assertion refuses a weighted value lacking CI+N), n_discordant_lineages, graded lineage_N annotation.
- A SCORED cell with no lineage row (or partial) renders raw + "lineage: not computed / incomplete (k/N)" — never silently blank (M1/M2).
- Reframe PROV_TIER -> "isolate-level provenance-disjoint stress test; R classes clonally dominated; lineage-effective N + cluster-weighted metrics (with CI) disclosed; NOT lineage-independent external validation." Keep no-aggregate-headline. lineage data augments, never demotes, the 7-state machine.

**Test strategy:**
- Unit: lineage-metrics fixture -> report-card JSON/MD carry new fields + CI + reframed tier; emitter assertion raises on a weighted value without CI; a partial cell renders "incomplete"; SCORED cells not removed; existing 12 tests unchanged.

### Step 4: Tests + docs
Files: tests/test_clonality.py, tests/test_compute_lineage_metrics.py, CLAUDE.md
Depends on: Step 1, Step 2, Step 3

**What changes:**
- New test files for the clonality module + recompute script; extend tests/test_build_validation_report_card.py for lineage columns/CI.
- Document the lineage-disclosure layer + the honest reframe + the clonality-inflation finding + greedy-vs-single-linkage rationale in the CLAUDE.md decoder-suite gotchas.

**Test strategy:**
- uv run pytest tests/ -q -> 0 regressions + all new tests green.

## Execution Preview
- Wave 0: Step 1.
- Wave 1: Step 2 (depends Step 1; also creates cell_key.py).
- Wave 2: Step 3 (depends Step 2 — sidecar schema + cell_key).
- Wave 3: Step 4.
- Total waves: 4. Max parallelism: 1 (linear data flow). Critical path: 1->2->3->4 (length 4).

## Risk Flags
- **Greedy-representative is order-sensitive** — sort strains by id before picking so the result is reproducible.
- **Genome re-fetch (kleb gent/tet, ~84)** — network-bound; M1 completeness-gate + partial-flag handle it; cohorts that can't complete emit no tier rather than a wrong one.
- **Single-linkage in mash_cluster_n147** — [unverified] confirm its linkage in Step 1; the new module deliberately uses greedy-representative regardless.
- **Cluster-weighted N is tiny (8-17) -> wide CIs** — that's the honest point; C3 mandates CI+N rendering.
- **Class-d cross-cutting + shared wiki/ namespace** — brainstorm done (Claude-only); a codex pass is advisable once convenient but not mandatory to proceed.

## Open Questions
- Graded lineage_N bucket boundaries (e.g. >=15 / 8-14 / <8) — pick at implementation; cosmetic.
- Whether DISCORDANT clusters get a tiny dedicated report-card sub-line or just the n_discordant_lineages count — lean count-only for v1.

## Verification
- uv run pytest tests/ -q -> 0 regressions + new tests green.
- Manual: scripts/compute_lineage_metrics.py -> sidecar has all 7 cells (or partial where genomes missing) with CI+effective-N.
- Manual: rebuild report card -> every SCORED cell shows raw + lineage cols + CI; tier reframed; no cell removed.

## Save-time amendments

Captured at: 2026-06-11
Source: /save-plan arguments

Audit-notes-only: these brainstorm findings are already folded into the Implementation Steps above (this plan was regenerated by /technical-plan after the brainstorm). Listed here for provenance.

- C1 whole-cohort clustering + DISCORDANT_LINEAGE category (not majority-vote mixed clones)
- C2 greedy-representative dedup not single-linkage (chaining)
- C3 cluster-weighted sens/spec MUST render Wilson CI + effective-N inline
- M1 genome completeness-gate + partial-cohort flag
- M2 canonical join key across 3 sidecars
- M3 graded lineage_N not constant >=20 tier
- M4 recompute with original rule-path args + assert raw reconciles

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
