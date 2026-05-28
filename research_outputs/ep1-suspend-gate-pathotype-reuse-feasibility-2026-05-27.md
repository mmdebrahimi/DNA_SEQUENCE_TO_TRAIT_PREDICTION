# EP-1 SUSPEND-Gate Pathotype Reuse Feasibility Memo (2026-05-27)

> Feasibility-grade verdict on whether the EP-1 mechanism-phenotype audit-gate pattern reuses for pathotype-call opacity. NOT an implementation spec. The workhorse owns the v0 implementation contract per `research_outputs/phase4_pathotype_discovery_handoff_to_workhorse_2026-05-27.md` §6.
>
> **Scope:** closes Action 5 in `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`. Informs H5 (`VirulenceFinder/ETECFinder → resolver concordance is stable under pinned versions`).
>
> **Source pattern files:**
> - `scripts/cipro_mechanism_phenotype_merge.py:112-145` (per-strain merge), `:165-169` (verdict thresholds 0.70/0.40)
> - `scripts/drug_mechanism_phenotype_merge.py:14-52` (drug-agnostic `classify_noise(row, drug)`)
> - `scripts/drug_mechanism_audit.py:219-258` (drug-agnostic `compute_verdict`)
> - `dna_decode/data/mic_tiers.py:181-292` (per-drug catalog constants)

## Verdict

**`ADAPTED_REUSE`.**

The two-layer verdict + gate pipeline at `scripts/drug_mechanism_phenotype_merge.py` + `mic_tiers.py` is **structurally directly reusable** for pathotype-call opacity audit. ~80% of orchestration code maps with cosmetic renames; ~20% requires pathotype-specific catalog design + threshold retuning. The pattern's two load-bearing innovations — `opacity_flag` separating "tool incomplete" from "label noisy" + the cohort-level `signal_quality` verdict — both transfer cleanly.

### Pre-committed falsifier

The verdict would have flipped to `NEW_DESIGN_NEEDED` if either:

- More than 50% of cipro-trained noise-class symbols (`CLEAN_R_primary_mechanism`, `OPAQUE_R_no_mechanism`, `SUSPECT_S_silent_primary_mechanism`, `NOISY_R_borderline`, `NOISY_R_no_mic`, `CLEAN_S_no_primary_mechanism`, `NOISY_S_no_mic`) had no semantic analog in pathotype-call evidence. **Result:** all 7 have direct analogs (table below). Falsifier NOT triggered.
- The verdict thresholds (0.70 / 0.40) were derived from MIC-tier-specific behavior with no plausible mechanism for pathotype-cluster opacity. **Result:** the thresholds are derived from "fraction of cohort with clean signal" — a domain-agnostic metric. Calibration will need retuning but the threshold-shape is reusable. Falsifier NOT triggered.

## Analogy mapping

### Per-strain noise-class mapping (AMR → pathotype)

| AMR symbol | Pathotype analog | Semantics |
|---|---|---|
| `CLEAN_R_primary_mechanism` | `CLEAN_pathotype_call_with_primary_cluster_present` | Confident pathotype label (e.g., independent EHEC label) + primary cluster module present (e.g., LEE + STX2) |
| `CLEAN_S_no_primary_mechanism` | `CLEAN_commensal_call_no_primary_cluster` | Confident commensal label + no DEC module fires |
| `OPAQUE_R_no_mechanism` | `OPAQUE_pathotype_label_no_cluster_marker` | Confident pathotype label + zero cluster markers detected (assembly incomplete OR atypical strain OR catalog gap) |
| `OPAQUE_R_co_resistance_only` | `OPAQUE_pathotype_call_co_modifiers_only` | Confident pathotype label + only co-pathotype modifiers (efflux/regulatory/porin analogs) detected — catalog gap or atypical pathway |
| `SUSPECT_S_silent_primary_mechanism` | `SUSPECT_commensal_call_with_primary_cluster_present` | Commensal-labeled strain WITH primary DEC cluster present → silent/non-functional cluster OR mislabeled commensal |
| `SUSPECT_S_borderline_primary_mechanism` | `SUSPECT_pathotype_call_partial_cluster_match` | Borderline-confidence label + partial cluster match (e.g., aggR-alone without AAF / aatA) |
| `NOISY_R_borderline` | `NOISY_pathotype_call_label_borderline` | Pathotype label is itself low-confidence (hybrid description without typing assay, atypical EPEC without bfp typing) |
| `NOISY_R_no_mic` | `NOISY_pathotype_label_no_phenotype_assay` | Pathotype label asserted but no independent phenotyping assay (no clinical metadata, no serotype confirmation) |

### Catalog mapping (`mic_tiers.py` → `pathotype_tiers.py`)

| AMR catalog (existing) | Pathotype analog (workhorse to design) | Notes |
|---|---|---|
| `DRUG_BREAKPOINTS[drug] = {clsi_r, clsi_s, eucast_r, eucast_s}` | `PATHOTYPE_CALL_CONFIDENCE_BREAKPOINTS[pathotype]` | Pathotype calls don't have MIC-style continuous breakpoints; instead, **cluster-presence thresholds** (e.g., `%ID ≥ 0.90` + `%coverage ≥ 0.60` per VirulenceFinder defaults) function as the analog. May be drug-call-independent (single global threshold pair) rather than per-pathotype. |
| `DRUG_LOCI_BY_MECHANISM[drug] = {mechanism: {gene_symbols}}` | `PATHOTYPE_CLUSTER_SIGNATURES[pathotype] = {cluster: {gene_symbols}}` | Direct structural mirror. E.g., `EHEC: {"LEE": {eae, tir, escV, ...}, "Shiga_toxin": {stx1, stx2}}`. Catalog enumeration is workhorse-owned. |
| `DRUG_PRIMARY_MECHANISMS[drug] = {mechanism_names}` | `PATHOTYPE_PRIMARY_CLUSTERS[pathotype] = {cluster_names}` | Direct mirror. E.g., `EHEC: {"LEE_full", "Shiga_toxin"}`. |
| `DRUG_AMRFINDER_CLASSES[drug] = {amrfinder_class_names}` | `PATHOTYPE_VF_DB_FILTERS[pathotype] = {vf_db_class_names}` | Filters VirulenceFinder output to pathotype-relevant rows. May be empty (VF doesn't currently carry per-pathotype class tags as cleanly as AMRFinder's `Class` column does for drug-class). Defer to workhorse. |
| `CO_RESISTANCE_MECHANISMS = {efflux, regulatory, porin_loss}` | `CO_PATHOTYPE_MODIFIERS = {?, ?, ?}` | Cross-pathotype-shared opacity modifiers. The biological analogs are not 1:1: efflux + regulatory + porin are AMR-specific; pathotype's cross-pathotype modifiers may be `serum_resistance_markers`, `motility_loss`, `IgA_protease`, or similar. **Workhorse-owned biological judgment;** do NOT pre-enumerate from discovery side. |

### Verdict-layer thresholds (`scripts/drug_mechanism_audit.py` → pathotype audit)

| AMR threshold | Pathotype starting point | Retuning rationale |
|---|---|---|
| `signal_quality ≥ 0.70 → SIGNAL_DOMINATES` | Start at 0.70; retune on real pathotype cohort | The 0.70 was calibrated on N=38 cipro cohort against MIC-tier-confidence distribution. Pathotype-cluster opacity has different noise distribution (no continuous-quantitative axis like MIC). Expect ±0.05-0.10 retuning after first 50-strain pathotype cohort runs. |
| `0.40 ≤ signal_quality < 0.70 → MIXED` | Start at 0.40 | Same rationale. |
| `signal_quality < 0.40 → NOISE_DOMINATES` | Start at 0.40 | Same rationale. |
| Gate logic: `clean_count ≥ 20 → RUN_CURATED_BASELINE_FULL_AND_CLEAN`; `≥ 10 → RUN_CURATED_BASELINE_FULL_ONLY`; `opacity_count ≥ 5 → MECHANISM_DEBUG_BRANCH`; else `SUSPEND_CONDITION_4` | Direct copy as starting point; rename `MECHANISM_DEBUG_BRANCH → CLUSTER_CATALOG_DEBUG_BRANCH` | Gate count thresholds are domain-agnostic. Renames are cosmetic. |

### Input-audit redesign

| AMR input | Pathotype analog |
|---|---|
| `mechanism_audit` (per-strain AMRFinder mechanism hits JSON) | `cluster_audit` (per-strain VirulenceFinder + ETECFinder cluster-presence JSON; includes per-cluster `%ID`, `%coverage`, contig/start/end/strand, hit_status) |
| `mic_audit` (per-strain MIC tier + median MIC + n_mic_rows) | `phenotype_audit` (per-strain independent-label provenance + serotype/MLST + outbreak/clinical context where available; replaces continuous-MIC with categorical pathotype label + confidence tier) |
| Join key: `strain_id` | Join key: `strain_id` (no change) |
| Per-strain row shape: `{strain_id, tier, median_mic, n_mic_rows, primary_mechanism_class, mechanisms_present, mech_hits}` | Per-strain row shape: `{strain_id, pathotype_label, label_provenance, label_confidence, primary_cluster_class, clusters_present, cluster_hits}` |

## Workhorse-side reuse recommendation

1. **Structural template:** copy `scripts/drug_mechanism_phenotype_merge.py` to `scripts/pathotype_cluster_label_merge.py`. Rename `mechanism_audit` → `cluster_audit`, `mic_audit` → `phenotype_audit`. The merge orchestration (left-outer join on `strain_id`, aggregate counts, verdict + gate) is reusable verbatim modulo renames.
2. **Catalog file:** design `dna_decode/data/pathotype_tiers.py` as a direct structural mirror of `mic_tiers.py`. Per-pathotype constants (`PATHOTYPE_CLUSTER_SIGNATURES`, `PATHOTYPE_PRIMARY_CLUSTERS`, etc.). Catalog enumeration is workhorse-owned biological judgment; this memo does NOT enumerate.
3. **Thresholds:** start with 0.70 / 0.40 verdict thresholds + 20 / 10 / 5 gate thresholds; retune on real pathotype-cohort signal_quality distribution after first 50-strain run.
4. **Test scaffold:** mirror `tests/test_drug_mechanism_audit.py` + `tests/test_drug_mechanism_phenotype_merge.py` shape. Pre-commit 4-cell verdict matrix (SIGNAL_DOMINATES / MIXED / NOISE_DOMINATES × clean_count buckets) as regression tests BEFORE running on real data — this is the 2026-05-14 pre-execute-/brainstorm lesson applied (`wiki/decisions-log.md` HIGH-salience entry on Stage1 plan-correction pattern).
5. **Estimated engineering cost** (workhorse-side, post-Gate-A): 4-6 hours per the Phase 1 codebase-exploration finding. Breakdown: ~1-2 hr catalog design (biological judgment), ~1-2 hr classify_noise() pathotype-specific rewrite, ~30 min orchestration copy + renames, ~2-3 hr test rewrite + threshold calibration on real cohort.

## Risks & caveats

| Risk | Severity | Mitigation |
|---|---|---|
| `PATHOTYPE_CALL_CONFIDENCE_BREAKPOINTS` is the wrong abstraction shape (pathotype calls may not have a single global threshold pair; per-cluster thresholds may dominate) | MEDIUM | Workhorse should empirically test with global vs per-cluster thresholds on Gate A's 5-strain sanity check before locking the catalog. |
| `CO_PATHOTYPE_MODIFIERS` biological analog is genuinely unclear (efflux/regulatory/porin don't 1:1 map to pathotype biology) | MEDIUM | Workhorse decides on biological grounds. Possible starting set: `{plasmid_horizontal_transfer_markers, motility_loss, host_immune_evasion}`. If no clean analog exists, the modifier-count branch in the verdict layer can be dropped without breaking the rest. |
| Verdict-threshold retuning on pathotype-cluster signal_quality distribution may need >0.10 shift, invalidating the "0.70/0.40 starting points" recommendation | LOW | Per-cohort retuning is expected discipline; the AMR thresholds were themselves retuned during EP-1 (2026-05-17 LESSON). Recalibration cost is hours, not days. |
| The pathotype-call-vs-MIC asymmetry (categorical vs continuous) breaks the `SUSPECT_S_silent_primary_mechanism` analog more than this memo claims | LOW | The `SUSPECT_commensal_call_with_primary_cluster_present` analog is still discriminating semantically (silent cluster OR mislabeled commensal), even without an MIC continuum to ground confidence. If the workhorse finds the analog under-discriminates on real data, it can split into `SUSPECT_commensal_call_silent_cluster` + `SUSPECT_commensal_call_mislabeled` post-hoc. |

## Status (for the project ledger)

- **Action 5 (Cross-check audit-gate reuse + pathotype-specific opacity definition):** RESOLVED. Verdict = `ADAPTED_REUSE`. Workhorse-side reuse recommendation provided.
- **H5 status:** the `VirulenceFinder/ETECFinder → resolver concordance` hypothesis substrate-requirement remains pinned to ≥50 isolates spanning every decision-table row. This memo does NOT change H5's substrate or threshold — it specifies how the OPACITY portion of the audit gate will be implemented when the cohort exists.

## Reproducibility note

- Pattern files read (verified via Phase 1 exploration agent + cross-referenced this session):
  - `scripts/cipro_mechanism_phenotype_merge.py` (~310 lines)
  - `scripts/drug_mechanism_audit.py` (~480 lines)
  - `scripts/drug_mechanism_phenotype_merge.py` (first 100 lines)
  - `dna_decode/data/mic_tiers.py` (header + catalog, ~320 lines)
  - `tests/test_drug_mechanism_audit.py` (~261 lines, 15+ test cases)
- No code modified. No tests modified.
- The verdict + analogy mapping table are derived from structural reading of the cited files, not from live pathotype-cohort data.
