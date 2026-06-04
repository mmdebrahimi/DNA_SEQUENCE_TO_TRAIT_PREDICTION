# Soraya minister run — ep4-v01-vf-diff — RESULT

- **Run id:** 2026-06-04-0008-ep4-v01-vf-diff
- **Mission shape:** in_place (improve the incumbent v0 pathotype resolver CLI; no family generation)
- **Verdict:** `mission-mvp-reached`
- **Φ (mission potential):** 2 → 2 (untouched; in_place generates nothing → finiteness trivially preserved)
- **Frontier:** [ep4-v01-vf-diff] (incumbent bound)

## Gate (attended audit-evidence, T1=(c))
- Round 1 PARKED on the interrogation gate (`blocked:user-only`) — gate precedes execution (ID3). Incumbent NOT run.
- 2 interrogation receipts recorded for `ep4-v01-vf-diff` after user ratification via AskUserQuestion (draft-then-ratify):
  - **Q1 canonical-VF path → native BLAST+ install** (ratified). NCBI BLAST+ 2.17.0 installed to `C:/Users/Farshad/ncbi-blast/bin` (blastn + makeblastdb).
  - **Q2 honesty bar → per-gene concordance + `caller_is_independent_baseline: false` + same-DB non-independence caveat** (ratified). No bare headline agreement %.
- Round 2 (receipts present) executed the incumbent's --until-mvp → both endpoints live-MET → terminal.

## What shipped (the --until-mvp work)
- `dna_decode/pathotype/vf_runner.py` — canonical VirulenceFinder caller via real `blastn` of the VF allele DB vs the assembly (VF's own method: identity + coverage thresholds, 90% / 60%), + `build_vf_diff` side-by-side builder. Offline-safe: degrades to `status: unavailable` (section retained) when blastn absent.
- `dna_decode/pathotype/cli.py` — `vf_diff` section added to the provenance JSON + a summary line; `--no-vf-diff` flag.
- `tests/test_pathotype_vf_diff.py` — 4 tests (shape, non-independence flag, offline-safe degradation, real-blastn STX2 both_present agreement on a synthetic stx2 assembly built from the DB).

## Endpoints (both MET, live-rechecked)
1. `test-exit-0 python -m pytest tests/test_pathotype_vf_diff.py` → 4 passed (system-python gated runner, exit 0).
2. `project-state-row …ecoli-pathotype-prediction-cli-2026-05-26.md:VF side-by-side diff shipped` → AC9 Action Log row 21 present.

## Verification
- vf_diff suite 4/4; full pathotype suite 31/31; full repo suite **789 passed, 0 failures**.
- Real blastn run on synthetic stx2 assembly: canonical STX2 97.3% id / 100% cov; resolver STX2 100% cov; agreement `both_present`; per-cluster concordance 1.0 (1 both+, 22 both-, 0 disagree).

## Honesty residual (carried into every diff)
Both callers match the SAME VF DB → high concordance is EXPECTED. The diff is an AUDIT of the fast
k-mer-seed caller against canonical blastn over one DB, NOT independent validation. (ledger Decision T1.)

## Next (locked sequence, §3 of the handoff)
v0 is now COMPLETE. Next: (2) freeze v0 as shipped + scope-capped; (3) de-confounded labels frontier
(within-study / lineage-matched ExPEC+EPEC) — a research/substrate epoch, NOT modeling.
