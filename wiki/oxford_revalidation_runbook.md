# Oxford external re-validation — runbook

End-to-end procedure for re-validating the FROZEN v0.5.0 AMR decoder on the Oxford
E. coli cohort (`PRJNA604975`), using the ingestion layer shipped 2026-06-15.

> **NOTE.** `cohort_manifest_external_<run_id>.json` is a per-run ARTIFACT (the exact
> scored-cohort definition). It is UNRELATED to the frozen `dna_decode/eval/cohort_manifest.py`
> leakage-registry MODULE. Do not route the artifact through the frozen module.

## 0. Prerequisites (the empirical gates — do these FIRST)
1. **Fetch the per-isolate MIC table.** Confirm it is openly downloadable (not MTA-gated).
   Its schema is UNVERIFIED until inspected — do not guess the columns.
2. **W0 probe** — pin the schema + crosswalk feasibility before anything codes against it:
   ```
   python scripts/oxford_w0_probe.py --project PRJNA604975 --mic-table <table> \
     --key-col <isolate-key-col> --drug-col CIP=ciprofloxacin --drug-col CRO=ceftriaxone --drug-col CN=gentamicin
   ```
   Read `wiki/oxford_w0_probe_<date>.json`: which candidate field the MIC key matches
   (run/sample_alias/secondary/BioSample), the MIC-key→BioSample resolution rate, the
   operator/censoring distribution, and duplicate-row counts. **Pin** `--key-col` +
   `--drug-col` map from this. A low resolution rate or many conflicts ⇒ fix the key
   choice before proceeding.

## 1. One-command run (after the gates)
```
python scripts/run_oxford_revalidation.py --project PRJNA604975 --mic-table <table> \
  --key-col <key> --drug-col CIP=ciprofloxacin --drugs ciprofloxacin
```
Chains: W0 probe (advisory) → `build_oxford_labels` (→ `cohort_manifest_external_<run_id>.json`,
ABORTS on crosswalk hard conflict) → exact-set `external_cohort_preflight --cohort-manifest`
(abort unless PASS) → `external_cohort_revalidate` per drug → run-scoped roll-up.

## 2. Gate + exit-code semantics
- **Preflight (exact-set):** `scored_gate=true`; FAILs on tuning-BioSample overlap, >5%
  unresolved, incomplete manifest (unless `--allow-degraded`), or zero FREE BioSamples.
  The whole-project path (`--project` without `--cohort-manifest`) is a non-scored
  `external_project_probe_*` diagnostic only.
- **Scorer exit codes:** `3` powering hard-fail (`n_scored==0` or strict-scored R/S below
  `--min-per-class` default 10) · `1` degraded (indeterminate fraction >20% of attempted-FREE)
  unless `--allow-degraded` · `2` gate refusal · `0` clean. Drift guard: `selected.tsv`
  BioSamples must be ⊆ the manifest.
- **Driver exit** = worst child by severity **3 > 1 > 2 > 0**. Roll-up runs ONLY IF every
  drug run is acceptable (0, or 1 with `--allow-degraded`).
- **Roll-up:** run-scoped (`--run-id`); refuses glob-all without `--allow-unscoped-glob`;
  skips `powering.hard_fail` cells; publishes `run_degraded` cells only with `--allow-degraded`
  (stamped, never headline).

## 3. Censoring
MIC operators are modeled (`external_mic_labels.MicValue`): a `<=`/`<` upper bound can only
support S (CENSORED_HIGH_S) and is NEVER called R; a `>`/`>=` lower bound can only support R
(CENSORED_HIGH_R) and is NEVER called S; mid-range censored → CENSORED_EXCLUDED.

## 4. Frozen invariant
The decoder + the 5 FROZEN files (`amr_rules.py`, `mic_tiers.py`, `cohort_manifest.py`,
`build_validation_report_card.py`, `compute_lineage_metrics.py`) stay byte-unchanged — this
arm validates the tool, it does not modify it.
