# EP-4 v0.1 ExPEC-recall — Corrected Round-2 Commitments (user-authored 2026-06-03)

> BINDING overrides for the minister maiden run `ep4-v01-expec-recall`. The as-committed
> round-2 answers were the riskier ones; these corrected answers REPLACE them. Apply these
> BEFORE `--until-mvp`. Recorded from this (wrong-cwd) session so the real run launched from
> a `dna_decode/`-rooted session uses them verbatim. NOTE: no journal existed yet — this is a
> fresh run, so these are the round-2 commitments of record, not an amendment to recorded receipts.

## Q1 — rescue scope
Rescue **JSMY + JSLG ONLY** → 11/12 = **0.917** (clears ≥0.85).
**Do NOT rescue capsule-only JSPG.**

## Q2 — support-branch threshold (cross-axis, no hand-tuned K)
Support rule = **"≥1 iron gene AND ≥1 capsule gene"** (cross-axis).
- Structurally excludes capsule-only JSPG (no hand-tuned pooled-K).
- **JSMY** rescued via this support rule.
- **JSLG** rescued via the existing **strong-adhesin (P-fimbriae) lever** (not the support rule).

## Q3 — DEC-module precedence (keep as committed)
DEC-module gate stays **above** the ExPEC support branch.
- ADD a regression test: an **LEE + full-support synthetic profile returns EPEC (not ExPEC)**.

## Pinning + test contract
- Pin the support rule as a **tested constant** (no magic numbers inline).
- The recall test must assert: **recall ≥ 0.85 AND confident-supported-call precision == 1.0**.
- Report the **lowest non-ExPEC support-gene-count margin** in the cohort (how close any
  non-ExPEC came to tripping the support rule — the overfit canary).

## Pre-authorized STOP condition
If **JSLG cannot be cleanly rescued** without weakening the strong-marker rule:
**STOP at `blocked:strategy-cap` (≈0.833)** rather than reaching back for JSPG.
A clean 0.833 beats an overfit 0.85.

## Launch (from a dna_decode-rooted session ONLY)
```
cd C:\Users\Farshad\PythonProjects\dna_decode
claude                      # fresh; do NOT resume a wrong-cwd session
/soraya minister ep4-v01-expec-recall
# then, per the maiden-run wiki:
SORAYA_RUN_ID=2026-06-03-1309-ep4-v01-expec-recall python minister_drive_ep4.py
```
