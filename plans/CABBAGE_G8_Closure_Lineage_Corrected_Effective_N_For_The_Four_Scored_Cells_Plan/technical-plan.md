# CABBAGE G8 Closure — Lineage-Corrected Effective-N for the Four Scored Cells — Technical Plan

> Close the one open label gate on the project's strongest label substrate by computing minority-class effective lineage counts for the four scored CABBAGE cells and feeding them into the existing code-owned G8 gate — in a namespace-separate sidecar, because two of the four cells collide exactly with existing provdisjoint keys.

## Lens status
Inputs: conversation
Degradations: repo-index absent (unindexed — no Layer 2 to refuse); `graphify-out/GRAPH_REPORT.md` absent; `project-rules.md` absent so `## Project Rules Check` is omitted; `DESIGN.md` absent but plan touches no UI; sentrux not installed; `python3` on this host is a Microsoft Store stub so the toolkit was run with `C:/ProgramData/Anaconda3/python.exe`; the EMBL-EBI portal `genomes/` tree URL schema is `[unverified]` (not fetched from this host); this laptop cannot execute the plan (no `D:` data tree, no `uv`, no Docker) — execution target is server2.

## Problem Statement

CABBAGE / EMBL-EBI AMR Portal is the strongest label substrate this project has identified: 1.7M lab-measured genome-phenotype pairs, CC-BY 4.0, 8.3 MB parquet, MD5-verified. Four cells were scored from it on 2026-10-10/11. Every one of them carries, in its own `honest_limits`, the sentence *"NOT CLONALITY-CORRECTED … The lineage layer (G8) is the open gate on this substrate and is NOT closed by this run."*

G8 is the last open gate on this substrate, and it is one of exactly three items in the mission's `amr-label-substrate-validated` gap. The gate itself is already code-owned: `rejection_gates.g8_dedup_collapses_balance` returns `TRIP if n < 3 else PASS` against `MIN_EFFECTIVE_LINEAGES = 3`, and returns `INSUFFICIENT_DATA` when no `min_effective_lineages` is supplied. Nothing supplies it. The arithmetic to produce it also already exists: `clonality.effective_lineage_n(clusters, labels, cls)` counts distinct same-label clusters of a class.

So this is not a new capability. It is wiring an existing measurement into an existing gate for four existing cells, plus the honesty layers that measurement requires.

Two facts make it non-trivial rather than mechanical:

- **The key collision is real, not hypothetical.** The CABBAGE cells key on `amrfinder_organism` — `Campylobacter`, `Escherichia_coli_Shigella`, `Salmonella`, `Salmonella`. The existing `wiki/provdisjoint_lineage_metrics.json` already holds `Campylobacter x ciprofloxacin` and `Escherichia_coli_Shigella x gentamicin`. Writing CABBAGE lineage metrics into that file under `canonical_cell_key` would silently overwrite two provdisjoint cells — the documented shared-key trap, and the reason the external-cohort arm has Fix C.
- **Clustering by determinant profile would be circular.** The portal ships a pre-computed AMRFinderPlus tree, which is the cheap-looking lineage proxy. It is not a lineage: the decoder's prediction is a function of those same determinants, so correcting a determinant-based decoder with a determinant-derived partition would manufacture agreement. Lineages must come from genome distance (Mash), not from the determinant tree.

## Codebase Context

- `dna_decode/eval/rejection_gates.py` (19,598 bytes) — `MIN_EFFECTIVE_LINEAGES = 3` at line 51; `g8_dedup_collapses_balance` at line 246 consumes `ev["min_effective_lineages"]`; `LABEL_GATES` includes G8. The gate needs a number, not a rewrite.
- `dna_decode/eval/clonality.py` (9,394 bytes) — `greedy_representative_clusters_from_matrix(distance_matrix, threshold)` (chaining-resistant greedy-representative, explicitly NOT `phylogeny.cluster_by_ani`'s single-linkage union-find), `effective_lineage_n(clusters, labels, cls)`, `cluster_weighted_confusion`, `cluster_members`, `cluster_class`, `wilson_ci`.
- `scripts/compute_lineage_metrics.py` (15,529 bytes) — the precedent that already solved this shape for the provdisjoint arm: `THRESHOLDS = (0.001, 0.005)` fixed at line 50, `GRADE_THRESHOLD = 0.005`, `reconcile_raw_metrics` (line 135) as an M4 gate, a genome-completeness gate emitting `partial` + `n_genomes_missing` and **no lineage tier** when incomplete, `canonical_cell_key` keying, per-cohort checkpointing, idempotent upsert.
- `scripts/cabbage_score_cell.py` (18,995 bytes) — the four cells' producer. Hardcodes `PHENO_PARQUET = "D:/dna_decode_cache/amr_portal/phenotype.parquet"` (line 42) and `CACHE = Path("D:/dna_decode_cache/amr_portal/amrfinder")` (line 43). Carries the fail-closed leakage gate via `cohort_manifest.build_manifest` / `prior_accessions`, and already names the registered constraint `cabbage-genomes-tree-coverage-uncensused`.
- Parquet columns used by that scorer: `species`, `antibiotic_name`, `assembly_ID`, `measurement`, `resistance_phenotype`, `database`, `collection_year`, `country`. **No MLST/ST column and no BioProject column** — so there is no cheap in-table lineage proxy, which is what forces Mash.
- `dna_decode/eval/phylogeny.py` — `compute_mash_distances(..., use_docker=True)`, batched to 2 Docker calls (1 sketch + 1 dist) rather than N(N−1)/2.
- `scripts/screen_candidate_gates.py` (13,207 bytes) — `--verify` re-derives committed hand verdicts. It currently reproduces four: `pear`, `hbv`, `oxford`, `phenosense`. **`grep -i cabbage` returns nothing** — CABBAGE is absent from the mechanical screen entirely, so its "9 of 10 gates pass" status is memo-only.
- Tests: `tests/` with the established offline-skip pattern (real-data tests skip when gitignored data is absent).

### Reusable-Code Survey
- `scripts/compute_lineage_metrics.py` — the whole shape: fixed thresholds, reconcile gate, completeness gate, checkpointing, idempotent upsert. Mirror it; do not re-derive it.
- `dna_decode/eval/clonality.py` — `effective_lineage_n` is literally the quantity G8 wants; `greedy_representative_clusters_from_matrix` is the chaining-resistant clusterer and `wilson_ci` the interval.
- `dna_decode/eval/rejection_gates.py` — `g8_dedup_collapses_balance` already owns the verdict; supply evidence, never re-implement the threshold.
- `dna_decode/eval/phylogeny.py` — `compute_mash_distances(use_docker=True)` for the distance matrix.
- `scripts/cabbage_score_cell.py` — the portal assembly fetcher + the fail-closed leakage gate to reuse verbatim.

## Pre-Change Baseline

Every number below was **measured this session by running code**, not recalled — per the 2026-10-12 lesson that a plan's baseline claims are measurements and three such claims in one plan were false.

Isolate-level confusion for the four scored cells, read from the committed JSONs:

| cell | tp | fn | tn | fp | n_scored | sens | spec | minority class (R) |
|---|---|---|---|---|---|---|---|---|
| Campylobacter jejuni × ciprofloxacin | 997 | 21 | 3236 | 11 | 4265 | 0.9794 | 0.9966 | 1018 |
| Escherichia coli × gentamicin | 55 | 6 | 7848 | 39 | 7948 | 0.9016 | 0.9951 | 61 |
| Salmonella enterica × ceftriaxone | 591 | 36 | 7190 | 15 | 7832 | 0.9426 | 0.9979 | 627 |
| Salmonella enterica × ciprofloxacin | 55 | 12 | 7421 | 346 | 7834 | 0.8209 | 0.9555 | 67 |

- `min_effective_lineages` is **absent for all four cells**; `g8_dedup_collapses_balance` therefore returns `INSUFFICIENT_DATA` today.
- All four cells' `honest_limits` contain the explicit "G8 … is NOT closed by this run" sentence.
- `wiki/cabbage_lineage_metrics.json` does not exist.
- `scripts/screen_candidate_gates.py --verify` reproduces **4 of 4** committed hand verdicts (pear, hbv, oxford, phenosense) and contains **no CABBAGE case**.
- `wiki/provdisjoint_lineage_metrics.json` exists and holds 10+ cells including the two colliding keys.
- Minority-class assemblies to resolve: 1018 + 61 + 627 + 67 = **1,773**.

## Verification Signal

- `g8_dedup_collapses_balance` returns `PASS` or `TRIP` (never `INSUFFICIENT_DATA`) for each of the four cells, or the cell is marked `partial` and emits no lineage number at all.
- `wiki/cabbage_lineage_metrics.json` exists with schema `cabbage-lineage-metrics-v1` and one entry per scored cell, each carrying `effective_lineage_N` per class, `largest_cluster_fraction`, `discordant_isolate_fraction`, and a Wilson CI beside every weighted point.
- `wiki/provdisjoint_lineage_metrics.json` is **byte-identical** before and after (sha256 compared, not `git status`).
- The four cell JSONs are **augment-only**: identical top-level key set, and the only differing field is the new `lineage` block plus the one rewritten `honest_limits` sentence. Verified by diff, not asserted.
- `scripts/screen_candidate_gates.py --verify` reproduces **5 of 5** hand verdicts with CABBAGE added.
- Frozen AMR surface (`amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`, `cohort_manifest.py`) byte-unchanged and the 2026-08-31 prospective lock re-verifies `ok=True`.
- `tests/test_cabbage_lineage.py` green offline with no network and no Docker.

## Implementation Steps

### Step 1: Namespace-separate sidecar, merge-or-refuse semantics, frozen verdict rule
Files: dna_decode/eval/cabbage_lineage.py, wiki/cabbage_g8_acceptance_bar.json
Depends on: none

**What changes:**
- New module `dna_decode/eval/cabbage_lineage.py` owning schema `cabbage-lineage-metrics-v1` and the sidecar path `wiki/cabbage_lineage_metrics.json`. It must NOT write to `wiki/provdisjoint_lineage_metrics.json` and must not be matched by that file's loader glob.
- **Merge-or-refuse semantics defined before the first increment**, per the 2026-10-14 lesson that staged writes over a shared artifact are silently lossy: `upsert_cell` merges by cell key and REFUSES when the incoming run would drop cells already present. A run that computed one cell may not emit a one-cell artifact replacing a four-cell one.
- `wiki/cabbage_g8_acceptance_bar.json` frozen before any lineage number exists: the bar is the **existing** `MIN_EFFECTIVE_LINEAGES = 3`, imported from `rejection_gates`, never restated as a literal. Records that **TRIP is a valid pre-committed outcome, not a failure to fix** — a substrate whose minority class collapses to <3 lineages is a finding about the substrate.
- Records the definitional decision: `min_effective_lineages` is the **minority class's** effective lineage count, because G8's own text is "smallest class has N effective lineages". States that the majority class is not separately clustered, that this is ~100× larger, and that lineage count is not monotone in isolate count so this is a scoped reading rather than a proof.

**Test strategy:**
- Unit: `upsert_cell` refuses a cell-dropping write; accepts an additive one; is idempotent on replay.
- Unit: the module's sidecar path is not the provdisjoint path, and the provdisjoint loader glob does not match it.
- Unit: the bar is read from `rejection_gates.MIN_EFFECTIVE_LINEAGES`, so changing that constant moves the bar (pins against a restated literal).

### Step 2: Minority-class assembly resolution with a completeness gate
Files: scripts/cabbage_lineage_metrics.py
Depends on: Step 1

**What changes:**
- New script. Reuses `cabbage_score_cell.py`'s portal fetcher verbatim and resolves only the **minority-class** assemblies per cell (1,773 total), not all ~27,879 scored isolates.
- M1 genome-completeness gate mirroring `compute_lineage_metrics.py`: a cell whose assemblies cannot all be resolved is marked `partial` with `n_genomes_missing` and emits **no lineage tier and no `min_effective_lineages`**. A partial cell leaves G8 `INSUFFICIENT_DATA` rather than reporting a number from a subset.
- Discharges the registered constraint `cabbage-genomes-tree-coverage-uncensused` by measuring coverage instead of assuming it.
- Reuses the fail-closed leakage gate: an incomplete `cohort_manifest` aborts rather than scoring.

**Test strategy:**
- Unit: a stubbed resolver missing 1 of N assemblies yields `partial=True`, `n_genomes_missing=1`, and no lineage number.
- Unit: a complete resolver yields `partial=False` and proceeds.
- Unit: incomplete manifest aborts unless the explicit degraded flag is passed, and then stamps `manifest_degraded=True`.

### Step 3: Raw reconciliation, then Mash distances, then greedy clustering at fixed thresholds
Files: scripts/cabbage_lineage_metrics.py
Depends on: Step 2

**What changes:**
- **M4 reconcile gate first:** recompute each cell's raw confusion through the deployed `call_resistance` with the artifact's own `amrfinder_organism` and drug, and assert it equals the committed `confusion` block before any weighted number is computed. Mismatch raises and reports rather than proceeding — the 2026-10-06 lesson that a reconcile gate against an already-published number is what reveals silent disagreement.
- `phylogeny.compute_mash_distances(..., use_docker=True)` over each cell's minority-class assemblies (batched 2 Docker calls per cell).
- `clonality.greedy_representative_clusters_from_matrix` at the **fixed** `(0.001, 0.005)` thresholds inherited from `compute_lineage_metrics.py`. **No threshold search.** The documented failure is that `variance_ratio = intra/inter` → 0 as clusters → singletons, so minimizing it is monotonically biased toward over-splitting; the coarse grid in `mash_cluster_n147.py` was load-bearing by accident.
- Determinant-profile clustering is explicitly forbidden here (circular — see Problem Statement).

**Test strategy:**
- Unit: a planted confusion mismatch raises and names the differing cell.
- Unit: clustering is called with exactly `(0.001, 0.005)` and no search path exists (grep-style guard on the module for a threshold-sweep call).
- Unit (pure, offline): a synthetic distance matrix with three obvious clades yields 3 clusters under greedy-rep, and a chaining case that single-linkage would merge stays separate.

### Step 4: Cluster-structure metrics, not just cluster counts
Files: dna_decode/eval/cabbage_lineage.py
Depends on: Step 1

**What changes:**
- `effective_lineage_n` per class via the existing `clonality` helper.
- **`largest_cluster_fraction` and `discordant_isolate_fraction`** alongside the count, per the TB lesson: a lineage-collapse metric must report cluster STRUCTURE, because 43 discordant clusters there hid 97% of the isolates and a count alone looked fine.
- Wilson CI beside every weighted point, and a `_assert_weighted_renderable`-style refusal to emit a weighted point without its interval — weighted N is small, so the interval is the point.
- Discordant (mixed-label) clusters are excluded from both class counts and counted separately; never majority-voted.

**Test strategy:**
- Unit: a partition where one cluster holds 77% of isolates reports `largest_cluster_fraction ≈ 0.77` even though the cluster count looks healthy.
- Unit: a mixed-label cluster lands in `discordant`, is excluded from both class counts, and is reported.
- Unit: emitting a weighted point without a CI raises.

### Step 5: Supply G8 its evidence and make the CABBAGE screen mechanical
Files: scripts/cabbage_lineage_metrics.py, scripts/screen_candidate_gates.py
Depends on: Step 3, Step 4

**What changes:**
- Pass `min_effective_lineages` into `g8_dedup_collapses_balance` per cell and record the returned `GateResult` verbatim in the sidecar. The verdict stays code-owned; the script supplies evidence only.
- Add CABBAGE as a case in `screen_candidate_gates.py` with its committed hand verdict, so `--verify` re-derives it mechanically. Today CABBAGE is absent from that harness, so the "9 of 10 gates pass" claim is prose-only — the 2026-10-13 lesson is that code-gated numbers survived and prose-only numbers did not, with no exceptions either way.
- G1 and G3 remain human-judgment gates requiring an explicit evidence string; the screen must still refuse a verdict without them rather than defaulting to pass.

**Test strategy:**
- Unit: supplying `min_effective_lineages = 2` yields `TRIP`; `3` yields `PASS`; omitting it yields `INSUFFICIENT_DATA`.
- Integration: `screen_candidate_gates.py --verify` reproduces 5 of 5 hand verdicts.
- Unit: the CABBAGE screen case refuses a whole verdict when the G1/G3 evidence strings are absent.

### Step 6: Augment-only re-render of the four cells and the synthesis
Files: wiki/cabbage_scored_cell_campylobacter_ciprofloxacin_2026-10-10.json, wiki/cabbage_scored_cell_escherichia_gentamicin_2026-10-11.json, wiki/cabbage_scored_cell_salmonella_ceftriaxone_2026-10-11.json, wiki/cabbage_scored_cell_salmonella_ciprofloxacin_2026-10-11.json, wiki/cabbage_four_cell_synthesis_2026-10-11.md
Depends on: Step 5

**What changes:**
- Attach a `lineage` block to each cell under its own key. Existing `metrics`, `confusion`, `leakage` and `source_concentration_of_scored_set` blocks are untouched — the layer augments and never replaces, exactly as the lineage / source-concentration / prospective / species-composition layers do on the provdisjoint card.
- Rewrite the one `honest_limits` sentence that asserts G8 is open, replacing it with the measured verdict. Leaving it would make the artifact contradict its own new block.
- Update `cabbage_four_cell_synthesis_2026-10-11.md` with the four G8 verdicts and the structure figures.
- Verify augment-only **by diff**: identical top-level key set per cell, and the only differing fields are `lineage` plus that one sentence.

**Test strategy:**
- Integration: load each cell before and after; assert equal key sets and that the symmetric difference of changed fields is exactly `{lineage, honest_limits}`.
- Integration: sha256 of `wiki/provdisjoint_lineage_metrics.json` unchanged.
- Manual: read one rendered cell and confirm the G8 line and the structure figures agree with the sidecar.

### Step 7: Offline regression suite
Files: tests/test_cabbage_lineage.py
Depends on: Step 5

**What changes:**
- Collision guard proven **non-vacuous** with deliberately different values, the way the prospective anti-overwrite test does: write a CABBAGE cell and a provdisjoint cell that share a key but carry different numbers, and assert neither clobbers the other. Simulate the merge bug and confirm the test fails.
- Merge-or-refuse: a one-cell run cannot replace a four-cell artifact.
- `partial` emits no lineage number and leaves G8 `INSUFFICIENT_DATA`.
- Both `PASS` and `TRIP` render correctly, so the TRIP branch is exercised before it is ever needed.
- No network, no Docker, no `D:` dependency; real-data assertions skip when the gitignored parquet is absent, following the established pattern.

**Test strategy:**
- The step is the tests. Run `pytest tests/test_cabbage_lineage.py -q` plus the full suite to confirm zero regressions against the 6,409-passing baseline.

## Execution Preview

Wave 0 (1 parallel):  Step 1 — Namespace-separate sidecar, merge-or-refuse semantics, frozen verdict rule
Wave 1 (2 parallel):  Step 2 — Minority-class assembly resolution with a completeness gate, Step 4 — Cluster-structure metrics, not just cluster counts
Wave 2 (1 parallel):  Step 3 — Raw reconciliation, then Mash distances, then greedy clustering at fixed thresholds
Wave 3 (1 parallel):  Step 5 — Supply G8 its evidence and make the CABBAGE screen mechanical
Wave 4 (2 parallel):  Step 6 — Augment-only re-render of the four cells and the synthesis, Step 7 — Offline regression suite

Critical path: Step 1 → Step 2 → Step 3 → Step 5 → Step 6 (5 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: critical — **Shared-key silent overwrite.** Two of the four CABBAGE cells (`Campylobacter × ciprofloxacin`, `Escherichia_coli_Shigella × gentamicin`) collide exactly with keys already in `wiki/provdisjoint_lineage_metrics.json`. Measured, not hypothetical. Resolved by Step 1's separate file + schema and Step 7's non-vacuous collision guard, with a sha256 check in Verification Signal.
- Severity: high — **This plan cannot run on the host it was written on.** The laptop has no `D:` data tree, no `uv`, and no Docker; `scripts/cabbage_score_cell.py` hardcodes two `D:/dna_decode_cache/amr_portal/` paths. Execution target is server2. Writing it here is fine; executing it here is not.
- Severity: high — **The repo copy on this host has no git** (`.git` is an empty directory, so git commands resolve to the home repo). Save and execute this plan where git exists, or the artifacts become another unversioned single-disk divergence.
- Severity: high — **G8 may TRIP, and that is a valid outcome.** The two small cells have minority classes of 61 and 67 isolates; the provdisjoint Klebsiella cipro R class was 60 isolates and collapsed to roughly 3 lineages. Step 1 pre-commits TRIP as a finding about the substrate, so a trip cannot be re-read as a bug to tune away.
- Severity: medium — **Portal `genomes/` tree coverage is unknown** (registered constraint `cabbage-genomes-tree-coverage-uncensused`), and `cabbage_score_cell.py` already notes fetch failures for absent assemblies. A large coverage gap forces cells to `partial` and yields no G8 number — a real possibility, handled rather than assumed.
- Severity: medium — **External tool integrations.** Mash via `phylogeny.compute_mash_distances(use_docker=True)`, verified against `dna_decode/eval/phylogeny.py` source (batched 1 sketch + 1 dist, pinned image `quay.io/biocontainers/mash:2.3--hb105d93_10` per CLAUDE.md). `g8_dedup_collapses_balance` signature verified by reading `rejection_gates.py:246-255`. `clonality` API verified by reading `clonality.py:59,234`. Fixed thresholds verified at `compute_lineage_metrics.py:50`. The EMBL-EBI portal `genomes/` URL schema is `[unverified]` — not fetched from this host; Step 2's completeness gate is the verification step.
- Severity: medium — **File overlaps.** Steps 2, 3 and 5 all write `scripts/cabbage_lineage_metrics.py`; resolved by making 3 depend on 2 and 5 depend on 3. Steps 1 and 4 both write `dna_decode/eval/cabbage_lineage.py`; resolved by making 4 depend on 1. No intra-wave overlap remains.
- Severity: medium — **Transitive imports not in Files: lists.** Steps 3 and 5 import `rejection_gates`, `clonality`, `phylogeny`, `cohort_manifest` and `call_resistance` without modifying them. `cohort_manifest.py` is one of the five sha256-pinned frozen-surface files, so a change there would retire the prospective lock — read-only use only, asserted in Verification Signal.
- Severity: low — **Minority-class-only clustering is a scoped reading of G8.** The majority class is ~100× larger and is not separately clustered; lineage count is not monotone in isolate count, so this is defensible but not a proof. Recorded in Step 1 and surfaced in Open Questions.
- Severity: low — **Restructuring applied.** None beyond the dependency additions noted above for file overlap.

## Open Questions

- Should the majority (S) class also be clustered, making `min_effective_lineages` a true minimum over both classes rather than the minority class's count? I have resolved this as minority-class-only for now, on G8's own "smallest class" wording and because the majority class would add ~26,000 assemblies. Reversing it is a scope decision, not a technical one.
- Is a `TRIP` on any cell grounds to withdraw that cell from the trust surface, or to disclose it the way the lineage/source-concentration layers disclose rather than demote? The project's standing pattern is disclose-never-demote, and Step 6 follows it — but if a trip should actually demote a SCORED cell, that is an authority call about the trust surface.

## Verification

1. On server2, with the repo at a committed state: `uv run pytest tests/test_cabbage_lineage.py -q`, then the full suite, confirming zero regressions against the 6,409-passing baseline.
2. `uv run python scripts/cabbage_lineage_metrics.py` for each of the four cells; confirm each emits either a `min_effective_lineages` with a G8 `PASS`/`TRIP` or `partial=True` with no lineage number.
3. `uv run python scripts/screen_candidate_gates.py --verify` prints 5 of 5 hand verdicts reproduced.
4. `sha256sum wiki/provdisjoint_lineage_metrics.json` identical to its pre-change value.
5. Diff each of the four cell JSONs pre/post; confirm identical key sets and that only `lineage` and the one `honest_limits` sentence changed.
6. Re-verify the frozen five are byte-unchanged and `prospective_lock.verify_lock` returns `ok=True`.
7. Read `wiki/cabbage_four_cell_synthesis_2026-10-11.md` and confirm its four G8 verdicts and structure figures match the sidecar.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high -->
