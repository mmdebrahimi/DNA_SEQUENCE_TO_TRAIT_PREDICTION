# Codex Active Handoff Packet — Bounded Falsifier Execution

**Generated:** 2026-05-23 (overnight session wrap by Claude on GTX 860M laptop)
**For:** Codex CLI session on Precision 7780 (RTX 3500 Ada)
**Repo:** `mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`
**Target commit:** `43bf94d` on `main` (origin pushed 2026-05-23)

---

## 0. Pre-execution: pull latest

Your repo on Precision 7780 needs the 11 new commits from the overnight session:

```bash
cd <your-dna_decode-checkout>
git fetch origin
git log --oneline HEAD..origin/main | head -15   # confirm 11 commits ahead
git pull --ff-only origin main
git log --oneline -1   # should show 43bf94d
```

If `git pull` fails with conflicts: STOP. There are local commits on Precision that conflict — escalate to user before proceeding.

---

## 1. Run the leakage check FIRST (< 5 s, blocking gate)

```bash
uv run python scripts/leakage_check_dup_accession.py
```

Expected:
- Reads `data/processed/stage2_n150_cipro_cohort.parquet` + `data/processed/models/ciprofloxacin_nucleotide_transformer.pkl`.
- Writes `reports/cipro_leakage_check_dup_accession_2026-05-22.{json,md}` (or the new run date if you want — the runner accepts whatever path you pass via `--leakage-check-json`).
- **If `loso_leakage_present == True`** in the output: ABORT the falsifier. Halt + send the JSON back to me. We dedup + retrain (Step L of `plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md`).
- **If `loso_leakage_present == False`**: proceed to step 2.

The duplicate accession in question: `GCA_025200635.1` shared by strain_ids `562.109860` AND `562.111036`.

---

## 2. Diff the runner draft against your spec

```bash
git show 1e2ad95 -- scripts/cipro_bounded_falsifier.py | head -50
# OR
cat scripts/cipro_bounded_falsifier.py
```

Authority split (per `wiki/cipro_bounded_falsifier_coordination_plan_2026-05-22.md` Section 1):
- **Codex owns runner mechanics:** model invocation, ISM call ordering, classifier-load idiom, cache `bulk_get` arg shape, exception handling.
- **Claude owns:** `StrainResult` schema, per-bucket pass criteria, verdict matrix → exit codes, diagnostic exports schema, `--leakage-check-json` gate.

If your runtime version diverges from Claude's draft in a way that breaks the schema, surface the divergence in the results MD before I dispatch. The 15 contract tests at `tests/test_cipro_bounded_falsifier.py` pin the verdict-matrix semantics — they should pass against your runner output too.

---

## 3. Run the bounded falsifier

```bash
uv run python scripts/cipro_bounded_falsifier.py \
  --cohort data/processed/stage2_n150_cipro_cohort.parquet \
  --model data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --cache "C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_repo/data/processed/cache/nt_n150_cipro.h5" \
  --refseq-cache "C:/Users/b0652085/OneDrive - Bombardier/Apps/Stress-DNA Project/dna_decode_repo/data/cache/refseq" \
  --subset wiki/cipro_bounded_falsifier_subset_2026-05-22.json \
  --leakage-check-json reports/cipro_leakage_check_dup_accession_2026-05-22.json \
  --drug ciprofloxacin \
  --output-prefix wiki/cipro_bounded_falsifier_results_2026-05-23
```

**Adjust `--cache` and `--refseq-cache`** to your actual local paths (the example above is from the 2026-05-21 audit JSON).

**Expected wallclock:** ~20 min (12 strains × ~95 s per ISM). If any single strain takes > 30 min: ESCALATE per the coordination plan Section 1 escalation triggers.

**Per-strain progress is logged to stdout:** `[falsifier] A_ERS 562.7693 (GCA_001284665.1) ...`

---

## 4. Expected outputs (transfer back to Claude)

Two artifact pairs land at:
- `wiki/cipro_bounded_falsifier_results_2026-05-23.{md,json}`
- `reports/cipro_leakage_check_dup_accession_2026-05-22.{md,json}` (from step 1)

Either:
- **Push to origin:** `git add wiki/ reports/ && git commit -m "artifact(cipro-falsifier): results from Precision 7780 run" && git push origin main`. Then Claude pulls them on the GTX 860M.
- **OR Gmail-transfer them** to `Downloads/` on the other laptop (rename `.json` → `.json.tab` only if Gmail blocks; `.md` should pass).

**Claude is waiting on these 4 files.** Anything else (e.g., diagnostic plots, runtime logs) is bonus.

---

## 5. Verdict matrix — what Claude does next on each outcome

| Verdict | Claude's next action |
|---|---|
| `PASS` | Step P of `plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md` — Mash-cluster N=147 (PASS-path artifact: `scripts/mash_cluster_n147.py` already in repo) + spec-tighten + tag `v0.0-cipro`. |
| `FAIL` | Step F — fill scope-limit doc, populate `attribution_scope_confidence` field with locus-tag-prefix proxy, ship v0 with documented limitation, tag `v0.0-cipro`. |
| `RUNNER_REGRESSION` | Step R — halt; diff Codex runner vs Claude draft at `reports/cipro_runner_diff_<DATE>.md`; re-run. |
| `INDETERMINATE_BUCKET_C` | Step V — ship as FAIL with "Bucket C status indeterminate at falsifier-tested fidelity" + queue follow-up plan. |

Step 0 verdict-reconciliation gate fires on every path: Claude re-computes the verdict from raw `per_known_locus` data locally; if local verdict ≠ your emitted verdict → escalate before dispatch.

---

## 6. Heartbeat protocol (if run exceeds 2 hr wallclock)

If the falsifier hasn't completed by 2 hr:
- Emit a 1-line heartbeat at `Downloads/codex_heartbeat_<HHMM>.txt` with: current bucket, current strain, per-strain elapsed time, any error.
- Claude on the GTX 860M will read it on next session start.

If the runner is failing or stuck on per-strain ISM:
- Try `--max-strains 2` (NOT yet implemented in runner — would need to add). Stand-in: edit `wiki/cipro_bounded_falsifier_subset_2026-05-22.json` to keep only 1 strain per bucket for fast iteration. Document the edit in the results MD.

---

## 7. What's already done on Claude's side (no action needed from you)

These 11 commits ship infrastructure + plans + tests that ride alongside the falsifier:

```
43bf94d docs(overnight): status report + CLAUDE.md gotchas + README Current state
d2003e9 draft(mash-cluster): N=147 orchestration script + 8 pure-logic tests
2755385 tests(cipro-falsifier): 15 contract tests + tightened _ranked_by docstring
773e8b0 feat(predict-v0): attribution_scope_confidence field + locus-tag-prefix proxy
5b0eae0 fix(cohort): assert assembly_accession uniqueness in build_cohort + regression tests
827387c tests(predict-v0 + mic-tiers): commit existing 98 tests previously untracked
bd6fb4f plan(cipro-post-falsifier): apply /review + /brainstorm edits in-place
1e72c2c plan(cipro-post-falsifier): verdict-conditional ship-path technical plan
1e2ad95 draft(cipro-falsifier): runner skeleton for Codex to diff against
a7424b0 ledger(cipro-falsifier): 3 LESSONS + Bellman refresh + decisions log + plans index
ae97ca9 plan(cipro-falsifier): coordination plan + 12-strain subset + leakage check snippet
```

664 tests green; +295 from session start; zero regressions.

---

## 8. Quick reference (1-screen run summary)

```
git pull --ff-only origin main
uv run python scripts/leakage_check_dup_accession.py
# (if loso_leakage_present=False)
uv run python scripts/cipro_bounded_falsifier.py \
  --cohort <cohort.parquet> --model <model.pkl> \
  --cache <cache.h5> --refseq-cache <refseq-root> \
  --subset wiki/cipro_bounded_falsifier_subset_2026-05-22.json \
  --leakage-check-json reports/cipro_leakage_check_dup_accession_2026-05-22.json \
  --output-prefix wiki/cipro_bounded_falsifier_results_2026-05-23
# (transfer the 4 result files back to Claude)
```

**Estimated runtime: ~25 min total (5 s leakage check + ~20 min falsifier).**
