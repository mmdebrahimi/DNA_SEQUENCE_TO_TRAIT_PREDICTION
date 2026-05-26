# Post-Codex-Push Operational Checklist

> Exact command sequence to run on this side (GTX 860M laptop) the moment Codex pushes the outstanding v0.1 cipro slice 1 + EP-0 close artifacts from Precision 7780 to origin. Pre-staged 2026-05-26 to remove "what was the sequence?" friction at the critical handoff moment.

**Anchor:** `Downloads/response_to_codex_v0_1_handoff_2026-05-26.md` was sent to Codex with the explicit "PUSH ORIGIN" request. When that lands on `origin/main`, run this checklist top-to-bottom.

---

## ⚡ FAST PATH (post-2026-05-26 — use the handoff gate runner)

The `plans/Two_Machine_Operating_Contract.md` §2.1 spec is now implemented as
`scripts/handoff_gate.py`. For the cef audit-aware closeout + every future
Codex push, the new fast path is:

```bash
cd C:/Users/Farshad/PythonProjects/dna_decode
git pull --ff-only origin main
uv run python -m scripts.handoff_gate
```

Verdict interpretation:
- **ACCEPTED** (exit 0): all 5 contract checks PASS. Handoff is integrated. Proceed to bundle-content-specific verification (which is what the rest of THIS file covers — Steps 3-9 are still relevant for the v0.1 cipro context they were drafted for).
- **PROVISIONAL** (exit 1): ≥1 check FAIL. Per contract §2.1, work depending on this artifact is provisional. Read the rendered `reports/handoff_gate_<DATE>.md` for per-check suggested actions.

The fast path subsumes steps 1-2 + 6 of the longer checklist below. Steps 3-5 + 7-9 (bundle-content verification + tagging + plans-index update + cef kickoff) remain manual per-handoff work.

---

## 0. Trigger

Codex on Precision 7780 reports `git push origin main` succeeded. Or user receives some signal (Slack, email, calendar nudge, etc.) that the push happened.

Optionally: poll with `git fetch origin && git log --oneline HEAD..origin/main | head -5` every N hours; if non-empty, the push landed.

## 1. Pull + sync verify (target: 1 minute)

```bash
cd C:/Users/Farshad/PythonProjects/dna_decode
git fetch origin
git log --oneline HEAD..origin/main | head -20   # confirm Codex's commits arrived
git pull --ff-only origin main
```

Expected NEW files (per Codex's 2026-05-25 handoff):
- `reports/dna_decoder_v0_release_candidate_2026-05-24.md`
- `reports/dna_decoder_v0_release_candidate_example_2026-05-24.{md,json}`
- `reports/cipro_v0_scope_limit_decision_2026-05-23.md`
- `reports/dna_decoder_v0_1_genome_input_release_candidate_2026-05-25.md`
- `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.{md,json}`
- `reports/dna_decoder_v0_1_genome_input_validation_2026-05-25.{md,json}`
- `tests/test_pipeline_predict_genome_input_e2e.py` (new test file)
- Updated `tests/test_data_annotations.py` (26 tests)
- Updated `scripts/pipeline.py` (`--genome-fasta`, `--annotations`, `--sample-id`, `--allow-missing-audit` flags + `cv_strategy` train-time persistence)
- Updated `dna_decode/data/annotations.py` (CDS extraction for genome-input)
- Retrained model: `data/processed/models/ciprofloxacin_nucleotide_transformer.pkl` (leakage-safe `leave_one_accession_out` CV; AUROC 0.8697)
- Bounded falsifier results: `wiki/cipro_bounded_falsifier_results_2026-05-23.{md,json}` (or similar naming)

If `git pull --ff-only` fails with non-fast-forward: STOP. Codex made commits on a branch that diverged from origin/main. Investigate before merging — DO NOT force-push.

## 2. Cross-machine sync diagnostic (target: 30 seconds)

```bash
uv run python scripts/cross_machine_sync_check.py --skip-pytest
```

Expected output:
- `commit-gap ahead=0 behind=0` (we're in lockstep with origin)
- `spec-divergence missing markers: 0/<N>` (all KNOWN_DIVERGENCE_TARGETS satisfied)
- `Downloads/ recent artifacts: <small N>` (only files newer than the latest commit)
- Working tree clean OR a small number of new untracked items only

If the sync check still reports DRIFT DETECTED after the push: investigate which markers are still missing; could indicate Codex pushed PARTIAL work.

## 3. Verify retrained model loads (target: 2 minutes)

```bash
uv run python -c "
import pickle
from pathlib import Path
p = Path('data/processed/models/ciprofloxacin_nucleotide_transformer.pkl')
with open(p, 'rb') as f:
    bundle = pickle.load(f)
print('keys:', sorted(bundle.keys() if isinstance(bundle, dict) else ['(not a dict)']))
prov = bundle.get('provenance', {}) if isinstance(bundle, dict) else {}
print('cv_strategy:', prov.get('cv_strategy'))
print('cv_auroc:', prov.get('cv_auroc'))
print('trained_on:', prov.get('trained_on'))
print('n_strains:', prov.get('n_strains'))
"
```

Expected:
- `cv_strategy: leave_one_accession_out`
- `cv_auroc: 0.8697...` (or similar)
- `trained_on: 2026-05-22` (or thereabouts)
- `n_strains: 146` (147 cohort minus the 1 dedup'd duplicate)

If the pickle has `cv_strategy: None` or `cv_auroc: None`: Codex's train-time persistence didn't activate; flag back to Codex.

## 4. Run the leakage-safe v0 predict smoke (target: 5 minutes)

Quick sanity smoke against a held-out cipro strain via the cached-strain path (the original v0 surface):

```bash
# Pick a strain_id from the N=147 cohort that wasn't in any test set
# (or use one Codex named in the v0 release example).
uv run python -m scripts.pipeline predict \
  --drug ciprofloxacin \
  --strain-id <held-out-strain-id> \
  --model-path data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --cache D:/dna_decode_cache/embeddings/nt_n147_cipro.h5 \
  --annotations <path/to/that-strain.gff3> \
  --audit-merge-json wiki/cipro_mechanism_phenotype_merge_2026-05-17.json \
  --output reports/post_push_cipro_smoke_2026-05-26.json
```

Verify the output JSON has:
- `prediction` (R or S)
- `calibrated_probability` (float)
- `confidence_tier` (HIGH / MEDIUM / LOW)
- `attribution_scope_confidence` (HIGH / PARTIAL / INDETERMINATE)
- `top_k_attribution` (list)
- `audit_verdict` (object with `suspend_gate_fired` boolean)
- `provenance.cv_strategy = leave_one_accession_out`
- `provenance.cv_auroc = 0.8697...`

If the output JSON is missing the v0.1 provenance fields: Codex's runtime `pipeline.py` changes didn't land cleanly; flag back.

## 5. Run the v0.1 genome-input smoke (target: 10-15 minutes — GPU bound)

This is the actual v0.1 ship-validation. **Requires Precision 7780 GPU** — must run there, NOT this laptop (GTX 860M can't fit NT v2 100M).

If running on this side somehow (CPU fallback, etc.): expect very slow embedding. Otherwise, ask Codex to re-run their `reports/dna_decoder_v0_1_genome_input_example_2026-05-25.{md,json}` smoke from the freshly-pulled origin state to confirm reproducibility post-push.

## 6. Run the full test suite (target: 1-2 minutes)

```bash
uv run pytest tests/ -q --tb=line
```

Expected count: **695 + Codex's new tests**. Per Codex's handoff:
- `tests/test_data_annotations.py` → 26 (was lower; ~22?; +N from Codex)
- `tests/test_pipeline_predict_genome_input_e2e.py` → new file; ~?? tests
- `tests/test_pipeline_predict_v0.py` → existing 31; Codex may have added more
- `tests/test_pipeline_cli.py` → existing 9; may have additions

Total expected: somewhere between 720-740 tests passing. Investigate any failures.

## 7. Tag v0.0-cipro + v0.1-cipro (target: 1 minute)

Once Steps 1-6 all pass cleanly:

```bash
# v0.0-cipro: tag at the EP-0 close commit (where retrained pickle + release packet landed)
git log --oneline --grep="v0 release candidate" | head -3
git tag v0.0-cipro <that-commit-hash>

# v0.1-cipro: tag at the genome-input slice 1 commit
git log --oneline --grep="genome-input" | head -3
git tag v0.1-cipro <that-commit-hash>

# Local-only initially. Don't push tags until user confirms.
git tag --list
```

Per `plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md` R5 risk flag, tags are LOCAL-ONLY until user explicitly approves push. Reverting a pushed tag requires force-push to the tag, which the project's commit pattern hasn't required before.

## 8. Update plans-index + ledger (target: 5 minutes)

```bash
# Add row to wiki/plans-index.md for any new Codex-shipped plans
# Update project_state/dna-decode-2026-05-11.md Bellman Current State
# Append Action Log rows recording the v0.1 cipro slice 1 landing + sync close
```

Specifically:
- Bellman Current State: "v0.0-cipro + v0.1-cipro tagged on origin/main; cross-machine sync IN SYNC; cef slice next per plans/Cef_Cached_Strain_V0_1_Slice_Plan.md."
- Action Log: 1 row for "Codex pushed v0.1 slice 1 + EP-0 close artifacts" + 1 row for "Claude pulled + tagged + verified".

Commit + push.

## 9. Kick off cef slice (target: starts execution; ~2-3 hr compute on Precision 7780)

Per `plans/Cef_Cached_Strain_V0_1_Slice_Plan.md` Step 1-6. Codex executes on Precision 7780.

Claude side: monitor; pull cef artifacts when Codex pushes; tag `v0.1-cef`.

---

## What to do if something breaks at any step

| Symptom | Action |
|---|---|
| `git pull --ff-only` fails | DO NOT force-pull. Investigate divergence; fetch + inspect `git log HEAD..origin/main` + `git log origin/main..HEAD`; possibly merge or rebase manually. Flag to user before acting. |
| Sync check reports drift even after pull | Codex pushed partial work; identify missing markers + ask Codex to push remainder. |
| Retrained pickle missing `cv_strategy` | Codex's train-time persistence didn't activate; runtime pipeline.py changes may have shipped without the cmd_train code change. |
| Predict v0 smoke output is missing fields | Codex's runtime pipeline.py changes broke the schema. Compare against `wiki/decoder_v0_ux_and_success_criterion.md` RELOCKED spec. |
| Test suite drops below 695 | Regression. Investigate which tests broke; binary-search the commit that introduced it. |
| Tag commands find no matching commits | Codex's commit messages may not match the grep patterns; use `git log --oneline -n 30` and pick manually. |

## What this checklist deliberately does NOT cover

- Re-running the EP-1A same-strain parity test (Codex's 4/4 + max delta 0.011599 already validated it)
- Re-running the EP-1B external benchmark (separate EP)
- Cef slice execution (`plans/Cef_Cached_Strain_V0_1_Slice_Plan.md` covers it)
- EP-1.5 architecture POC (gated on EP-2 / EP-3 multi-drug expansion needs)
- Mass cleanup of `Downloads/` (one-off artifacts; orthogonal to the push)

## Bottom line

8 numbered steps + ~25 min wallclock if everything is clean. Worst case (~30-60 min) if Codex's push has subtleties to investigate. Either way, this checklist replaces "what was the sequence again?" with a single read-through.
