# Soraya --until-mvp — J2 Phase 1 free-GPU runner (2026-07-04-1834)

## Directive
"use your best judgement and choose one, or perhaps even distribute it between those for faster response.
whatever you think best, decide and move forward" (execute-mode). Context: prior two turns established J2
Phase 1 (ESM2-650M inference → the real ~0.48 field number) can run FREE on Kaggle/Colab. The user asked me
to pick a free-GPU path or distribute across them, and move forward.

## Decision (execute-mode, draft-then-ratify)
Build the **J2 Phase 1 free-GPU runner** — Kaggle primary + Colab fallback + **built-in sharding** so the
user can "distribute between those" (shard 0/2 Kaggle + shard 1/2 Colab in parallel → merge → ~half
wall-clock). Reuse the EXISTING validated scorer (`scripts/esm_zeroshot_dms.py`), do not reinvent it.

## MVP bar (frozen at run start) — ALL MET
| # | Criterion | Kind | Status |
|---|---|---|---|
| C1 | `notebooks/j2_phase1_esm2_proteingym.py` — self-contained sharded runner | file-exists | ✅ MET |
| C2 | `notebooks/J2_PHASE1_RUNBOOK.md` — Kaggle + Colab + distribute recipe | file-exists | ✅ MET |
| C3 | `tests/test_j2_phase1_notebook.py` passes offline (no GPU/net) | test-exit-0 | ✅ MET (10/10) |

## What shipped
- **Self-contained runner** (upload ONE file to a fresh Kaggle/Colab kernel): ESM2-650M masked-marginals
  over joinable ProteinGym DMS assays → Spearman + shuffled control + median; `--shard i/n` (deterministic
  strided, disjoint, covers-all) + `--out JSON` + `--merge` (no-GPU pooling) + `--fetch` (reference file).
- **Drift guard:** the runner's pure scoring core (spearman / parse_variant / constants) is asserted
  byte-equal to the canonical `scripts/esm_zeroshot_dms.py` by `test_drift_guard_scoring_core_matches_canonical`
  — the self-contained copy cannot silently diverge from the validated pipeline.
- **Runbook:** data-into-kernel (upload D:/dna_decode_cache/proteingym as a Kaggle Dataset = zero URL
  fragility, uses the exact validated data) + single-kernel path + distribute-across-two path + expected
  output + honest scope.

## Verification (verify-in-batch)
- Offline unit tests: **10/10 green** (spearman monotonic/shuffle, parse, load_dms multi-mutant filter,
  shard disjoint+covers-all+deterministic, merge dedup/pool, drift guard, constants match).
- **Real-surface plumbing check** (against the actual cached D:/dna_decode_cache/proteingym): `load_reference`
  found all **217** assays, all with CSV + target_seq present; shard 0/2=109 + 1/2=108, disjoint + covers-all;
  sample paths resolved. The GPU/model half is NOT runnable here (860M laptop; won't pull 650M to tight C:) —
  that is the user's Kaggle step by design (login-gated, like OpenGWAS/UKB).

## Honest scope / wall classification
- **The NUMBER is not produced yet** — the deliverable is the ready-to-run, offline-tested, real-plumbing-
  verified runner. The ~0.48 lands when the user runs it on a free GPU. External wall = **free-GPU kernel +
  the user's Google/Kaggle login** (not code-closable here — cloud login is user identity).
- Phase 1 mostly RE-CONFIRMS a known field number (the falsifier already proved the substrate at 0.417 via
  AlphaMissense). The real bet is Phase 2 (beat 0.48 via training) — deliberately NOT built here; gated on
  "does beating it power a decoder feature."
- Frozen AMR surface (amr_rules.py + calibrated_amr_rules.json) byte-unchanged — purely additive.

## Terminal
`mvp-reached` (all 3 criteria live-MET). No spillover into Phase 2 (it's the GPU/training bet, an authority
+ resource fork — not adjacent reversible work). Next real move is the user's: run the runbook on Kaggle,
paste the median back.

## No-resume
Run complete. Re-invoke after the user runs the runbook (paste the number) → then the Phase 2 go/no-go fork.
