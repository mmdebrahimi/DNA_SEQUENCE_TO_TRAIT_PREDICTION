# Cipro Decision Bundle Plan

> Tier-0 cheap decision bundle on N=38 + BV-BRC MIC census, fired BEFORE any Databricks burst or per-gene NT diagnostic. Replaces the rejected binary Path A vs Path B framing after SUSPEND_CONDITION_4 verdict landed 2026-05-17.

> **Executable scope:** see `plans/Cipro_Decision_Bundle_Technical_Plan.md`. This document is now the reconciled rationale (design-decision narrative + open tradeoffs context). All steps + file metadata + waves live in the technical plan. Last reconciled 2026-05-17 to absorb /review reductions: Bakta smoke (D7) dropped; mechanism-completeness differential-test (D8) dropped; 3 of 5 manifest variants dropped (only `label_original` + `label_exclude_no_mic_borderline` remain); `--ignore-gate` CLI flag dropped; mean+max preflight v3 deferred to a conditional follow-up plan.

---

## Problem Statement

The N=38 cipro cohort has been declared structurally unusable for testing PIVOT TRIGGER condition 4. The pre-curated-baseline gate fired `SUSPEND_CONDITION_4` because:

- Cohort signal quality = 0.17 (clean_count = 7 of 40; only 7 HIGH_R / 0 HIGH_S strains)
- Opacity_count = 0 → AMRFinder is NOT the bottleneck; the BV-BRC labels themselves are unreliable.
- 6 SUSPECT_S strains carry a primary cipro mechanism + borderline MIC = almost certainly mislabeled.

State of evidence (2026-05-17):
1. Stage 1 N=38 mean-pool — FAIL (NT-XGB 0.568 vs k-mer 0.648).
2. Stage 1b N=38 mean+max + scaled NT-LR — FAIL (NT-LR 0.673 vs k-mer 0.648, gap +2.5 pp).
3. Cipro attribution preflight v2 — INCONCLUSIVE_MISS (zero QRDR + expanded-set loci in top-K=20).
4. Experiment 1 (AMRFinder all 38) — QRDR_DOMINANT (18/20 R have QRDR; 7/20 plasmid; 2/20 no mechanism).
5. Experiment 2 (raw BV-BRC AST MIC rejoin) — NOISY (7 HIGH_R / 0 HIGH_S; 9 R have no MIC; 12 S borderline).
6. Merge (Experiment 1 × 2) — NOISE_DOMINATES; gate fires SUSPEND_CONDITION_4.

The original binary decision (Path A = N=150 strict-MIC cohort expansion via Databricks burst; Path B = per-gene NT windows on 7 CLEAN_R strains) is rejected: both fail core feasibility checks. Path A spends burst on an unverified strict-MIC universe; Path B tests a different architecture (per-locus) than the one that failed (whole-genome pooling).

The right next move is a free/cheap Tier-0 decision bundle that produces an actionable verdict on whether burst spend or per-gene diagnostic is justified, BEFORE committing to either.

## Design Decisions

### D1: Replace binary A/B with a 4-tier decision bundle

**Decision:** Build a tiered experiment bundle with Tier 0 (cheap, free, ≤1 hr each), Tier 1 (small refactor + run), Tier 2 (medium investment gated on Tier 0), Tier 3 (large spend gated on Tier 0+1+2).

**Rationale:** Both Path A and Path B as originally framed answer questions other than the load-bearing one. The load-bearing question after SUSPEND_CONDITION_4 is "can the cipro phenotype be cleanly labeled at any reasonable cohort size in BV-BRC?" — this is a feasibility question, not a model question.

**Trade-off:** A linear "fire Path A now" approach is faster if burst is justified; the tiered approach takes 1-2 days longer but avoids burst-spend on an uncensused universe.

### D2: BV-BRC census must cover multiple phenotype policies

**Decision:** The Tier-0 census reports strain counts under at LEAST three phenotype policies: (a) HIGH_R / HIGH_S extreme-only (MIC ≥ 8.0 / ≤ 0.125), (b) CLSI R / CLSI S excluding I and borderline, (c) EUCAST R / EUCAST S excluding CLSI/EUCAST disagreement. Each policy is cross-tabulated with assembly_accession presence, MLST presence, assembly_quality, and lineage-confounding metrics (MLST counts + largest-clade fraction).

**Rationale:** Today's MIC audit shows HIGH_S = 0 because HIGH_S requires MIC ≤ 0.125 (1/4× CLSI-S). Many real S strains sit at 0.25 (EUCAST boundary, CLSI S) — the current threshold under-counts feasibility. Without multi-policy reporting, the census may falsely kill expansion.

**Trade-off:** Multi-policy reporting is more script complexity but the alternative — running census once per policy — is harder to keep consistent.

### D3: Frozen label-policy manifest before any relabeling experiment

**Decision:** Emit a persisted `cipro_label_manifest_<date>.parquet` artifact with columns: strain_id × accession × original_label × MIC tier × mechanism class × 2 variant labels × inclusion flag. All downstream label-sensitivity experiments consume this manifest; no hand-coded label overrides. Variants RETAINED: `label_original`, `label_exclude_no_mic_borderline`. Variants DROPPED 2026-05-17 per /review reductions: `label_mechanism_derived` (circular by construction — uses mechanism to label, then mechanism as feature); `label_suspect_s_as_r` (post-hoc-risky, too tempting to relabel-until-model-works); `label_exclude_no_mic` (standalone — subsumed by `label_exclude_no_mic_borderline`).

**Rationale:** `scripts/stage1_n40_cipro.py:158-167` builds `labels_by_strain` directly from `s.ast_labels`. Without a frozen manifest, label-sensitivity reruns are post-hoc relabel search — confirmation-bias-risky. The manifest forces pre-declaration. Minimum viable variant set is 2: original (regression baseline) + exclude_no_mic_borderline (the noise-stress-test).

### D4: Primary estimand for relabeling = per-strain error concentration + rank-order stability, NOT max AUROC

**Decision:** Pre-declare that the relabeling variants are noise-stress-tests. The decision-bearing output is (a) do errors cluster on the noise-class subset, (b) does the rank ordering of held-out predictions stay stable across relabelings. Max AUROC across variants is informational only, NOT evidence of model ability.

**Rationale:** Relabeling SUSPECT_S as R is defensible only as a phenotype-noise stress test. Reporting "best AUROC under best relabeling" = optimizing the label policy to fit the model, a known confirmation-bias pattern.

**Trade-off:** Strict discipline costs a real shipping outcome — if relabeling actually helps, we cannot report it as Phase 1 success without adopting the new phenotype target formally.

### D5: Mean+max attribution preflight v3 is a closeout falsifier, not a fork-decider

**Decision:** Refactor `gene_level_mutagenesis` at `dna_decode/interp/mutagenesis.py:124,142` to accept `aggregation` kwarg (default "mean" for backwards compatibility). Update `scripts/cipro_attribution_preflight.py` to train a mean+max NT-XGBoost classifier + pass `aggregation="mean+max"` to the mutagenesis call. Run as v3 preflight. Schedule as Tier 1, AFTER Tier 0 lands.

**Rationale:** Preflight v2's own verdict says DAMNING_MISS requires mean+max attribution. The audit trail leaves an objection until v3 runs. But: even if v3 misses, strongest conclusion is "current head's gene-level attribution misses" — not "NT can't represent QRDR." So it's a closeout falsifier, not the next decision-driver.

### D6: Curated baseline informational run needs a 3rd verdict field

**Decision:** Add `given_suspended_gate: INFORMATIONAL_ONLY` field to `cipro_curated_baseline.py` JSON payload. The script's existing 2-layer verdict (original_condition_4 + amended_condition_4) is kept; the 3rd field signals readers that a PASS is informational unless the cohort-label audit is rescued. The originally-proposed `--ignore-gate` CLI flag DROPPED 2026-05-17 per /review reductions: the curated baseline never reads the merge gate, so the flag would be a no-op. The "informational only" annotation is a JSON-field constant.

**Rationale:** Without the 3rd field, readers may confuse "amended condition 4 PASS" with "actionable result." With SUSPEND fired, any curated-baseline verdict on N=38 is descriptive only.

### D7 (DROPPED 2026-05-17 per /review reductions)

**Decision:** Bakta 4-strain smoke was originally specified as a Tier 2 step with a 4-strain selection (2 CLEAN_R + 1 SUSPECT_S_QRDR + 1 borderline_S no-mech) to estimate false-positive gene-symbol leakage. /review surfaced that this is Phase 2 prep folded into a Tier 0 decision plan; its output doesn't change the Path A / Path B decision. Lifted out of this plan; will return in a standalone Phase 2 Bakta annotation prep plan if Path A is alive.

### D8 (DROPPED 2026-05-17 per /review reductions)

**Decision:** Mechanism completeness differential-test (AMRFinder-first; manual blastn only on discordant rows) was originally specified as a Tier 2 step. /review surfaced that the underlying question (is AMRFinder missing QRDR?) is already answered by today's `opacity_count = 0` finding — none of the HIGH_R strains lack a primary mechanism. The step is redundant given the audit outputs already in hand.

### D9 (Open tradeoff): Phase 1 deliverable framing affects bundle weights

**Tradeoff:** The bundle weights split on the unresolved Phase 1 EP1 deliverable definition:
- "Publish at end of Phase 1" (defensible negative/indeterminate): audit hygiene wins (census + manifest + curated w/ given_suspended_gate + mean+max v3 + mechanism discordance review).
- "Ship a working classifier": data acquisition wins (census first, then strict-MIC cohort build if feasible, then phenotype-policy selection, then classical baselines before more NT attribution).

**Status:** Decision NOT made — user to lock before bundle Tier 2/3 fires. Tier 0 + Tier 1 are equally load-bearing under both framings.

## Implementation Plan

### Tier 0 — fire first (parallel-eligible after Step 1)

**Step 1. BV-BRC MIC feasibility census with multiple phenotype policies**
- Create `scripts/bvbrc_cipro_feasibility_census.py`.
- Inputs: `C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv` (AST) + `BVBRC_genome (1).csv` (assembly metadata). Reuses parsers from `dna_decode/data/ast_data.py` + `dna_decode/data/bvbrc_genome.py` where possible.
- For each cipro AST row: tier (HIGH_R / HIGH_S / DECISIVE / BORDERLINE / AMBIGUOUS / CONFLICT / NO_MIC) under 3 phenotype policies.
- Cross with assembly_accession non-null, MLST non-null, contig_count + N50 thresholds.
- Report strata: method (broth_microdilution / disk / etest), source/testing-standard, duplicate-row policy sensitivity.
- Lineage-confounding sub-census: per policy, MLST count + largest-clade fraction + R/S separability by MLST alone.
- Output: `wiki/cipro_bvbrc_feasibility_census_<date>.{md,json}`.
- Output gate: does any phenotype policy reach ≥75R/≥75S decisive strains with MLST + downloadable assembly + assembly_quality pass? If yes, Path A is alive UNDER THAT POLICY. If no, condition 4 is structurally untestable at N=150 from BV-BRC alone.
- Runtime: ~5 min.

**Step 2. Frozen N=38 label-policy manifest + per-strain error audit**
- See technical plan Step 7 (manifest subcommand) for the executable spec. Folded into the single `scripts/cipro_feasibility_and_label_audit.py` script per /review's accepted reduction (shared `_confidence_tier` + algorithm overlap with Step 1 census).
- Manifest variants RETAINED (2 of 5): `label_original` (unchanged BV-BRC binary) + `label_exclude_no_mic_borderline` (drop NO_MIC + BORDERLINE + AMBIGUOUS in one go).
- Manifest variants DROPPED 2026-05-17 per /review reductions: `label_exclude_no_mic` (standalone — subsumed by combined), `label_exclude_borderline` (standalone — subsumed by combined), `label_mechanism_derived` (circular by construction), `label_suspect_s_as_r` (post-hoc-risky).
- Per-strain error audit: see technical plan Step 10.5 (`scripts/cipro_error_audit.py`) for the executable spec. Consumes Step 8's new JSON sidecar (`wiki/stage1_n40_cipro_*.predictions.json`).
- Decision rule (unchanged): strong clustering of errors on noise_class subset = label noise is the bottleneck; uniform errors = architecture bottleneck.

**Step 3. Pre-declared single label-sensitivity LOSO**
- Run Stage 1b runner with the single retained noise-stress-test variant: `label_exclude_no_mic_borderline`.
- Output: `wiki/stage1_n40_cipro_mean-plus-max_<date>_label_exclude_no_mic_borderline.md` + matching JSON sidecar.
- Primary estimand: per-strain error concentration + rank-order stability vs original. NOT max AUROC.

### Tier 1 (DEFERRED 2026-05-17) — Mean+max attribution preflight v3

**Step 4. Mean+max attribution preflight v3** — DEFERRED to a separate plan, conditional on the decision-cell `True_low_threshold` outcome from Step 11. Even if v3 misses, strongest conclusion is "current head's gene-level attribution misses" — not "NT can't represent QRDR." So it's a closeout falsifier, not the next decision-driver. The `gene_level_mutagenesis` aggregation kwarg refactor still ships (technical plan Step 2) so the refactor is ready when needed.

### Tier 2 (DROPPED 2026-05-17 per /review reductions) — Bakta + mechanism completeness

Both steps dropped from this plan. Bakta 4-strain smoke was Phase 2 prep folded into a Tier 0 decision plan; lift to standalone Phase 2 Bakta annotation prep plan if Path A is alive. Mechanism completeness differential-test was redundant given today's `opacity_count = 0` finding.

### Tier 3 — gated on Tier 0 census passing + D9 framing locked

**Step 7. Conditional Path A (N=150 strict-MIC cohort build)**
- ONLY fire if Step 1 census shows ≥75R/≥75S decisive strains under some phenotype policy.
- Extend `scripts/build_stage2_n150_cohort.py` to accept the chosen phenotype policy as input.
- Databricks burst spend justified IFF the strict-MIC cohort is achievable AND D9 framing is "ship a working classifier."

**Step 8. Conditional Path B (per-gene NT windows on 7 CLEAN_R strains)** — DEFERRED as a representation sanity check; not pivot-deciding.

## Verification

See `plans/Cipro_Decision_Bundle_Technical_Plan.md` §Verification for the executable verification gates. Summary:

- Census artifact lists 3 phenotype policies × strain counts × MLST coverage; HIGH_R under `high_extreme` policy matches today's audit (7 in N=38).
- Manifest parquet has 38 rows × 2 variant labels (`label_original` + `label_exclude_no_mic_borderline`) + matching inclusion flags. Per-strain error audit shows whether errors are noise-class-clustered.
- Stage 1b runner accepts the manifest, runs LOSO under the chosen variant, emits JSON sidecar + markdown packet with the pre-declared estimand (NOT max AUROC).
- `cipro_decision_cell.py` refuses to write the runtime artifact without the PC1/PC2 pre-conditions markdown.

**End-to-end gate (after Tier 0 lands):** the runtime decision artifact `wiki/cipro_decision_bundle_runtime_<date>.json` contains `decision_cell` + `pc1_framing`; Step 12's runtime check is mechanical. The user can answer "should Tier 3 (Path A burst) fire?" with grounded yes/no.
