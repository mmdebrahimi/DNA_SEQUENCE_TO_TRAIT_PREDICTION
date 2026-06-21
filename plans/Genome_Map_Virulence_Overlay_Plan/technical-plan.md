# Genome-map virulence-determinant overlay tier — Technical Plan (v2, brainstorm-hardened)

## Lens status
Inputs: `wiki/genome_map_virulence_overlay_feasibility_2026-06-19.md` (probe-grade GREEN) + the pre-exec `/brainstorm` grounding in `features/genome-map/feature.md` (3 criticals C1/C2/C3 + M1/M2 + the empirically-inverted coord-join finding) + the shipped genome-map (`dna_decode/genome_map/`, 102 tests) + the pathotype resolver (`dna_decode/pathotype/`, RULES_VERSION pathotype-rules-v0.2.0).
Degradations: repo-index unconfigured (direct reading); no `/feature-design` record (the 5 design Q's + Q5 carried in Open Questions); no `project-rules.md`/`DESIGN.md`. This is the v2 post-/brainstorm regeneration — the criticals are folded into the Steps, NOT stacked as amendments.

## Problem Statement
Add a 5th `virulence-determinant` overlay tier to the genome-map honesty tool: where a curated VirulenceFinder (VF) allele is present in a single E. coli genome, surface it behind the SAME coordinate-join integrity gate + a presence-only wall as the AMR `determinant-phenotype` tier (presence of a curated determinant, NEVER a learned pathogenicity claim), with the deterministic pathotype-resolver call shown SEPARATELY as the genome-level overlay (the virulence analog of the AMR R/S call). E. coli/Shigella-scoped in v1 (the committed DB is `data/virulencefinder_db/virulence_ecoli.fsa`). Offline-degrades without blastn.

The pre-exec `/brainstorm` folded in (all [grounded]): **C1** `genome_pathotype_call` must mirror the deployed `resolve_call` contract (qc_pass + the ExPEC cross-axis support inputs) or it silently over-claims `COMMENSAL_LOW_MARKER_BURDEN`/CONFIDENT on a low-QC genome; **C2** the tier must surface ALL called VF alleles (the DB has 4942, cluster mapping keeps ~23) — only clustered hits feed the resolver, + capture the DB SHA; **C3** the canonical VF caller does best-hit-per-allele + `-max_target_seqs 5`, losing tandem copies — add an all-hits coord-retaining mode + interval dedup; **inverted Issue 3** (empirical): `makeblastdb` WITHOUT `-parse_seqids` yields `sseqid` = the exact FASTA header = the SAME token AMRFinder reports → the shared `build_contig_name_map` is valid for both overlays (the plan must NOT add `-parse_seqids`); **M1** detect non-unique contig first-tokens as degraded join quality; **M2** isolate VF join-quality metrics from the AMR `all_joins_symbol_fallback` + the AMR GO/NO-GO spike gate.

Non-goals: a learned virulence/pathogenicity predictor; non-E.coli VF DBs; folding virulence into the AMR tier; a virulence GO/NO-GO spike gate (the spike stays AMR-only in v1); ANY edit to the frozen AMR surface (`dna_decode/eval/amr_rules.py` + `calibrated_amr_rules.json`); ANY change to `build_vf_diff`'s cluster-audit semantics.

## Codebase Context
- `dna_decode/pathotype/vf_runner.py::run_canonical_vf` (`:79`) — runs `blastn` (VF DB vs assembly) with outfmt `6 qseqid pident length qlen` (`:109`, NO coords), best-hit-per-allele (`:130-131`), `-max_target_seqs 5` (`:111`); `_cluster_for_allele` (`:69`) maps only `CLUSTER_MARKERS` prefixes; `makeblastdb` called WITHOUT `-parse_seqids` (`:103`). `find_blastn` (`:45`) resolves native `C:/Users/Farshad/ncbi-blast/bin`. `NON_INDEPENDENCE_CAVEAT` (`:172`). NOT frozen.
- `dna_decode/pathotype/resolve.py::resolve_call(profile, *, partial_clusters, qc_pass, support_gene_count, cross_axis_support)` (`:21`) — `:158` returns `COMMENSAL_LOW_MARKER_BURDEN`/CONFIDENT; the ExPEC rescue is gated on `cross_axis_support` (`:138`), not `support_gene_count` alone.
- `dna_decode/pathotype/detect.py::assembly_qc(contigs)` (`:111`) — returns `{n_contigs, total_bp, n50, qc_verdict}` (key is `qc_verdict`, NOT `pass`).
- `dna_decode/pathotype/expec_score.py` — `support_gene_count(pergene_cov)` (`:29`) + `meets_cross_axis_support(pergene_cov)` (`:45`); `pergene_cov` = `{lowercase-gene-prefix: coverage∈[0,1]}`.
- `dna_decode/pathotype/cli.py` (`:32-37`) — the reference contract: `assembly_qc(contigs)` → `resolve_call(profile, qc_pass=(qc_verdict!="FAIL"))`.
- `dna_decode/genome_map/phenotype_overlay.py` — `DeterminantHit` (`:61`) + `join_hits` (`:protein-id → coord-overlap → symbol-fallback`, region-exclusion + smallest-span tie-break) + `build_contig_name_map` (length reconciliation) + `all_joins_symbol_fallback`. Reused verbatim.
- `dna_decode/genome_map/build_map.py::build_genome_map` (`:106`) + `__init__.py` tier constants + `TIER_PRECEDENCE` (`:37`). The assembler + vocabulary.
- `scripts/genome_map_spike.py::run_genome_map_for` + `scripts/genome_map.py` (CLI) — the wiring points.

### Reusable-Code Survey
- `phenotype_overlay.{DeterminantHit, join_hits, build_contig_name_map, all_joins_symbol_fallback}` — coord-join machinery, reused verbatim.
- `vf_runner.run_canonical_vf` — extended with an all-hits coord-retaining mode (C3); `_cluster_for_allele` reused; `build_vf_diff` UNTOUCHED.
- `pathotype.{resolve_call, detect.assembly_qc, expec_score}` — the genome-level call + its full honesty contract.
- `build_map.build_genome_map` + the `__init__` constants — assembler + tier vocabulary.
- None from: graphify-out/GRAPH_REPORT.md (absent), src/lib/utils dirs (absent). Grounding consulted: the brainstorm @ 2026-06-19 entry + `wiki/genome_map_virulence_overlay_feasibility_2026-06-19.md`.

## Pre-Change Baseline
- The genome-map has 4 tiers + an AMR-determinant overlay only; 102 genome-map tests; full suite 1472 passed.
- `run_canonical_vf` emits identity+coverage, NO coords, best-hit-per-allele; the pathotype resolver is standalone (never wired into the genome-map).
- Frozen AMR surface byte-frozen — this feature touches NONE of it (vf_runner + genome_map + pathotype non-frozen).
- EMPIRICAL (this host): `makeblastdb` without `-parse_seqids` → `sseqid` = the exact FASTA header (`CP021689.1`); WITH it → mangled (`gb|CP021689.1|`). Minus-strand HSPs emit `sstart>send`.

## Verification Signal
- On an E. coli genome with VF alleles: a `virulence-determinant` tier appears ONLY behind a HIGH-confidence coordinate join (symbol-fallback VF joins are visible secondary evidence, never the tier); each carries gene + cluster (if any) + the DB SHA + the non-independence caveat (presence, never "pathogenic"). ALL called VF alleles surface (not just the 23 resolver clusters).
- The genome-level pathotype call appears SEPARATELY, computed via the FULL `resolve_call` honesty contract (qc_pass + cross-axis support); a low-QC genome yields `AMBIGUOUS_LOW_QC` / `insufficient_context`, NEVER a confident commensal call.
- VF join-quality metrics are ISOLATED — they never touch the AMR `all_joins_symbol_fallback` nor the AMR GO/NO-GO spike gate. Non-unique contig first-tokens flagged as degraded join quality.
- Non-E.coli/Shigella → no virulence overlay (`SKIPPED_NON_ECOLI`). No blastn → `virulence_status=UNAVAILABLE_NO_BLASTN`, map still emits. Frozen AMR surface byte-unchanged. `uv run pytest tests/test_genome_map_*.py tests/test_pathotype_*.py -q` green; full suite 0 regressions.

## Implementation Steps

### Step 1: Coordinate-retaining VF caller (all-hits mode + per_hit + the no-parse_seqids pin)
Files: dna_decode/pathotype/vf_runner.py, tests/test_pathotype_vf_diff.py, tests/test_vf_runner_coords.py
Depends on: none

**What changes:**
- Add an `all_hits=False` param to `run_canonical_vf`. Extend the blastn outfmt to `6 qseqid sseqid sstart send pident length qlen`; parse `sseqid/sstart/send` per HSP; normalize `start=min(sstart,send)`, `stop=max(sstart,send)`, `strand` from order. Raise `-max_target_seqs` (e.g. 10000) when `all_hits=True` so tandem/multi-copy alleles are retained.
- Add a NEW `per_hit` list to the return: `{allele_id, vf_gene, cluster (None if unclustered — C2), sseqid, start, stop, strand, percent_identity, percent_coverage, called}` — built from ALL called HSPs (NOT cluster-filtered). Interval-dedup per `allele_id + sseqid` (collapse overlapping HSPs of one copy; keep distinct copies at distinct coords). Keep `per_gene`/`per_cluster` byte-identical (cluster-scoped, best-hit — `build_vf_diff` UNCHANGED). Offline-degrade emits `per_hit: []`.
- **Do NOT add `-parse_seqids`** (a comment pins WHY: it mangles `sseqid` to `gb|…|` and breaks the contig-name map; verified empirically 2026-06-19). Capture the VF DB SHA256 into the result (`db_sha`).

**Test strategy:**
- Synthetic blastn-outfmt stdout (incl. a minus-strand reversed pair + two tandem copies of one allele at distinct coords + an unclustered allele) → `per_hit` parsed with normalized coords, both copies retained, the unclustered allele present with `cluster=None`; `per_gene`/`per_cluster` byte-identical vs pre-change (regression pin in `test_pathotype_vf_diff.py`). LIVE integration fixture (`test_vf_runner_coords.py`, skipped when blastn absent): a synthetic reverse-strand allele → `makeblastdb`+`blastn` → assert exact `sseqid` = the FASTA first-token, normalized coords, AND `-parse_seqids` absent from the command.

### Step 2: Virulence-determinant tier constant
Files: dna_decode/genome_map/__init__.py, tests/test_genome_map_virulence_overlay.py
Depends on: none

**What changes:**
- Add `TIER_VIRULENCE_DETERMINANT = "virulence-determinant"`; insert into `TIER_PRECEDENCE` immediately AFTER `TIER_DETERMINANT_PHENOTYPE` (AMR wins if a feature is both; near-zero overlap), before `TIER_CURATED_FUNCTION`; export in `__all__`. `PHENOTYPE_TIER` (the AMR wall constant) unchanged.

**Test strategy:**
- `virulence-determinant` is at `TIER_PRECEDENCE` index 1; `PHENOTYPE_TIER` unchanged.

### Step 3: VF overlay module — adapter + join + the full-contract pathotype call + organism scope
Files: dna_decode/genome_map/virulence_overlay.py, tests/test_genome_map_virulence_overlay.py
Depends on: Step 1, Step 2

**What changes:**
- `parse_virulence_hits(canonical_vf_result) -> list[DeterminantHit]`: map each `per_hit` to a `DeterminantHit` (symbol=`vf_gene`, cls="VIRULENCE", subclass=`cluster or ""`, method="blastn", protein_id=None, contig=`sseqid`, start/stop). Only `called` hits.
- `join_virulence(features, hits, contig_name_map) -> (joined, counts)`: thin reuse of `phenotype_overlay.join_hits` (protein-id never fires → coord-overlap; symbol-fallback retained but excluded from the tier). Emit `virulence_join_quality` (`n_vf_rows/n_high_confidence_join/n_symbol_fallback/n_unjoined`) + `all_virulence_joins_symbol_fallback`. **M1:** detect non-unique contig first-tokens (a helper over the contig map) → a `n_ambiguous_contig` count surfaced in the counts (degraded join quality, never silent).
- **C1** `genome_pathotype_call(canonical_vf_result, contigs) -> dict`: build the cluster→bool profile from `per_cluster[*].called` (clustered only); build `pergene_cov` `{lowercase-prefix: max(percent_coverage)/100}` from `per_gene` (mapping allele→prefix); compute `qc = detect.assembly_qc(contigs)` + pass `qc_pass=(qc["qc_verdict"]!="FAIL")`, `support_gene_count=expec_score.support_gene_count(pergene_cov)`, `cross_axis_support=expec_score.meets_cross_axis_support(pergene_cov)` into `resolve_call`; carry `NON_INDEPENDENCE_CAVEAT` + the DB SHA. Stamp `status="insufficient_context"` ONLY when `contigs`/the VF result is unavailable/error (never when support is merely absent).
- `virulence_organism_in_scope(organism) -> bool`: a pure helper — True for `Escherichia`/`Escherichia_coli`/`Escherichia coli`/`Escherichia_coli_Shigella`/`Shigella*` (case/underscore-insensitive, genus-token based), False otherwise.

**Test strategy:**
- `per_hit` with coords → high-confidence coord join → DeterminantHits incl. tandem copies; a coords-absent/ambiguous-length hit → symbol-fallback (excluded), `all_virulence_joins_symbol_fallback` flagged; `n_ambiguous_contig` counted on a duplicate-first-token contig set. `genome_pathotype_call`: a low-QC contig set → `AMBIGUOUS_LOW_QC` (never COMMENSAL/CONFIDENT); a STX+LEE profile → EHEC-class with the caveat + DB SHA; an iron+capsule pergene_cov → cross-axis ExPEC rescue fires; unavailable contigs → `insufficient_context`. `virulence_organism_in_scope`: Escherichia/Shigella variants True; Klebsiella/Salmonella False.

### Step 4: Assemble the virulence tier + the genome-level pathotype overlay (metric isolation)
Files: dna_decode/genome_map/build_map.py, tests/test_genome_map_virulence_overlay.py
Depends on: Step 2, Step 3

**What changes:**
- `build_genome_map` gains optional `virulence_joined_hits` + `virulence_join_counts` + `pathotype_call`. Per-feature tier precedence: high-confidence AMR join → `determinant-phenotype`; ELSE high-confidence virulence join → `virulence-determinant`; else `classify_feature_tier`. The virulence WALL: a `virulence` field populated ONLY on a `virulence-determinant` feature (`vf_gene` + `cluster` + the pathotype the cluster contributes to (clustered only) + the non-independence caveat + DB SHA); symbol-fallback VF hits → `secondary_evidence` only.
- Metrics gain `virulence_determinant_feature_count` + `virulence_join_quality` (incl. `n_ambiguous_contig`) + `all_virulence_joins_symbol_fallback` + a top-level `genome_pathotype_call` block. **M2:** these are SEPARATE keys — `build_genome_map` must NOT let any virulence metric feed `all_joins_symbol_fallback` (the AMR key) nor the spike gate (`gate.py` reads only the AMR `join_quality`/`all_joins_symbol_fallback` — unchanged here; a test pins that an all-symbol-fallback VF set leaves the AMR gate + verdict untouched).

**Test strategy:**
- Synthetic features + a high-confidence virulence join + no AMR → `virulence-determinant` with the `virulence` field; symbol-fallback VF → NOT the tier, visible as secondary; the AMR phenotype wall + the virulence wall both hold; `genome_pathotype_call` surfaced separately; raw fields retained. **M2 isolation test:** an all-symbol-fallback VF set + a clean AMR overlay → `all_joins_symbol_fallback` (AMR) is False + `evaluate_gate`/`aggregate_spike_verdict` unchanged; `all_virulence_joins_symbol_fallback` True.

### Step 5: Wire the live VF run into the single-genome path + CLI + render
Files: scripts/genome_map_spike.py, scripts/genome_map.py, tests/test_genome_map_cli.py
Depends on: Step 4

**What changes:**
- `run_genome_map_for` (+ the CLI): when a FASTA is present AND `virulence_organism_in_scope(organism)` AND `--no-virulence` not set, run `vf_runner.run_canonical_vf(fasta, VF_DB, all_hits=True)` → `parse_virulence_hits` → `join_virulence` (REUSING the SAME `build_contig_name_map` already built for the AMR join) → `genome_pathotype_call(result, contigs)`; pass into `build_genome_map`. Stamp top-level `virulence_status` ∈ {FULL, UNAVAILABLE_NO_BLASTN, SKIPPED_NON_ECOLI, SKIPPED_USER}.
- Offline/degrade: no blastn → `UNAVAILABLE_NO_BLASTN`, map still emits (no virulence tier, exit 0). Non-E.coli/Shigella → `SKIPPED_NON_ECOLI`. `--no-virulence` → `SKIPPED_USER`. Render surfaces `virulence_status` + the virulence-determinant count + the genome pathotype call + the DB SHA.

**Test strategy:**
- CLI offline (monkeypatch `run_canonical_vf` → unavailable) → `UNAVAILABLE_NO_BLASTN`, map emits, exit 0; non-E.coli organism → `SKIPPED_NON_ECOLI`; `--no-virulence` → `SKIPPED_USER`; a monkeypatched blastn-present path → a virulence-determinant feature + the pathotype call + DB SHA in the JSON + md. Frozen AMR surface byte-unchanged.

## Execution Preview

Wave 0 (2):  Step 1 — coord-retaining VF caller, Step 2 — virulence tier constant
Wave 1 (1):  Step 3 — VF overlay module (adapter + join + full-contract pathotype call + organism scope)
Wave 2 (1):  Step 4 — build_map virulence tier + metric isolation
Wave 3 (1):  Step 5 — live VF wiring + CLI + render

Critical path: Step 1 → 3 → 4 → 5 (4 waves)
Max parallelism: 2 agents (Wave 0)

Note: Steps 1 + 2 touch disjoint files; Steps 3/4 share `tests/test_genome_map_virulence_overlay.py` — sequenced by `Depends on`.

## Risk Flags
- Severity: medium — **VF coord-join is the AMR integrity crux** (the gene-symbol trap). Mitigated by reusing `join_hits` + `build_contig_name_map` + symbol-fallback exclusion + the all-symbol-fallback guard verbatim + the LIVE integration fixture pinning `sseqid`/coords. [grounded — brainstorm]
- Severity: medium — **Honesty-rail (C1):** `genome_pathotype_call` over-claims `COMMENSAL`/CONFIDENT on a low-QC genome unless it passes the full `resolve_call` contract (qc_pass + cross_axis_support). Folded into Step 3 + a low-QC test. [grounded — brainstorm]
- Severity: medium — **Coverage over-claim (C2):** a cluster-scoped tier silently drops ~4900 of 4942 VF alleles. Folded: `per_hit` = all called hits; only clustered feed the resolver; DB SHA captured. [grounded — brainstorm]
- Severity: medium — **External-tool surface:** the blastn `-outfmt 6` column order + minus-strand `sstart>send` + the NO-`-parse_seqids` contract are load-bearing. Empirically verified on this host; pinned by the integration fixture. Native blastn present on THIS host only → other machines hit the offline-degrade path. [grounded]
- Severity: low — **Metric poisoning (M2):** VF metrics must not feed the AMR spike gate. Folded: separate keys + an isolation test. [grounded — brainstorm]
- Severity: low — **Claim-surface widening (Q5):** the genome-level pathotype call moves the tool from "AMR resistance" to "pathotype" — an AUTHORITY decision (Open Question A). Drafted YES with the C1 honesty contract; scope-out path documented. [grounded]
- Severity: low — **File overlap:** Steps 2/3/4 share `tests/test_genome_map_virulence_overlay.py`; sequenced by `Depends on`.

## Open Questions
- A. **(AUTHORITY — ratify before /execute-plan)** Q5: the genome-level overlay = the pathotype resolver call (EPEC/EHEC/ETEC/UPEC/EAEC/ExPEC/commensal). DRAFTED YES (with the C1 QC-gated honesty contract — a low-QC genome can't yield a confident commensal). To SCOPE OUT: drop `genome_pathotype_call` from Steps 3-5 (the per-feature `virulence-determinant` tier stands alone as "located VF determinant present"); everything else is unchanged. Confirm.
- B. Tier precedence when a feature is BOTH an AMR and a virulence determinant — DRAFTED: AMR `determinant-phenotype` wins (near-zero overlap). Confirm.
- C. Resolved (feasibility memo + brainstorm): own tier (not folded into AMR); presence-only wall; E.coli/Shigella scope v1; VF non-independence caveat + DB SHA carried; all-called-hits coverage; no `-parse_seqids`.

## Verification
1. `uv run pytest tests/test_pathotype_vf_diff.py tests/test_vf_runner_coords.py tests/test_genome_map_virulence_overlay.py tests/test_genome_map_cli.py -q` — all green.
2. `uv run pytest tests/ -q` (excluding `tests/test_models_foundation.py`) — 0 regressions; frozen AMR surface byte-unchanged; the AMR GO/NO-GO spike gate behavior unchanged (M2).
3. A live single-genome run on a cached E. coli FASTA (D:-free: native blastn + committed VF DB) → a `virulence-determinant` tier behind a high-confidence coord join + the separate genome pathotype call (full QC contract) + `virulence_status=FULL` + DB SHA; symbol-fallback VF joins excluded from the tier; tandem copies each get a feature.
4. Offline (no blastn) → `virulence_status=UNAVAILABLE_NO_BLASTN`, map still emits; non-E.coli/Shigella → `SKIPPED_NON_ECOLI`; a low-QC genome → pathotype call is `AMBIGUOUS_LOW_QC`/`insufficient_context`, never confident commensal.

## Save-time amendments

Captured at: 2026-06-19
Source: `/save-plan` arguments

> Audit-notes-only: `/execute-plan` reads ONLY `## Implementation Steps`. These amendments are provenance — the fixes below are ALREADY folded into the Steps above (this is the v2 post-/brainstorm regeneration, not a v1+amendments stack).

- genome-map virulence-determinant overlay v2 (brainstorm-hardened): coord-retaining VF caller + all-called-hits + no-parse_seqids pin, full-contract genome_pathotype_call (QC + cross-axis ExPEC), metric isolation from the AMR spike gate, E.coli/Shigella scope, DB SHA, offline-degrades
- Saved to plans/ (folder-per-plan) NOT features/genome-map/technical-plan.md — that slot holds the EXECUTED AMR plan; this avoids the silent-overwrite trap. Feature grounding stays in features/genome-map/feature.md.
- Q5 (genome-level pathotype call = AUTHORITY decision) remains OPEN — drafted YES with the C1 QC-gated honesty contract; scope-out path documented in Open Questions A.


<!-- toolkit: check=clean waves=clean gate=fired:open-questions mode=feature:genome-map-virulence -->
