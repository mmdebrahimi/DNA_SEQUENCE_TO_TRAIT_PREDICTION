# Decisions pending, 2026-11-28

Moved out of `NEXT.md` the same day it was written. That file is declared transient and sits
under a 9,000-word ceiling (`tests/test_claude_md_growth_ceiling.py`); it stood at 8,606 words
with 394 of headroom, and this block cost 704. The ceiling's own test file records the
precedent for exactly this situation -- *"a ceiling that moves whenever it binds is not a
ceiling"*, and both a previous raise and its addition were reverted -- while compressing other
sessions' notes to make room is a user authority call. So the content lives here and `NEXT.md`
keeps a pointer, which is the pattern that file already uses elsewhere.

## DECISIONS NEEDED 2026-11-28 — three authority calls, none of them mine to make

**1. A frozen bar and its grader have drifted apart.**
`dna_decode/glm/frame_verdict.CHECK_ORDER` holds **nine** outcomes;
`wiki/promoteratlas_frame_acceptance_bar.json` names **eight**. The grader gained
`INDETERMINATE_CONTROL_PREMISE_UNMET` at index 4, so **it can now return a verdict the
pre-registered bar gives no reading for**. Both files are committed; this predates 2026-11-28.
Reconciling means either adding an outcome to a pre-registered bar after the fact or removing a
refusal state from the grader — both authority calls, and the standing rule is to prove
outcome-neutrality before editing a frozen bar. Caught only by
`test_freeze_frame_bar.py::test_the_check_order_MATCHES_the_graders_check_order_exactly`;
`test_the_check_order_lists_every_outcome_exactly_once` PASSES because it checks the artifact
against its own outcomes, so a self-consistency check cannot see drift between two files.

**2. Ratify (or reject) the 2-of-5 baseline-panel definability split.**
`wiki/glm_baseline_panel_ce_2026-11-28.json` measures the two mandated baselines that ARE defined
for a per-position cross-entropy readout (site-independent, shuffled null), records one-hot ridge as
defined on the ACTIVITY criterion and already measured there, and records mutation-count and
BLOSUM62 as structurally UNDEFINED with reasons. **Nothing was written to the budget-4 gap's
`file-exists` artifact**, because that path is also the self-firing expiry argument of
`north-star-step4-animals-deferred` and one of five non-vacuity fixtures in
`tests/test_mission_gap_predicates.py` — writing it would close a mission gap, re-open the north
star's step 4 and flip a test fixture at once. If the split is ratified, the gate can be met by that
artifact plus a random-proposal null. Until then the gap stands open, correctly.

**3. The working-tree `wiki/data_inventory.{json,md}` is a DEGRADED census — do not commit it.**
73 datasets / 5 orphans against the committed 257 / 49. Cause: `D:/dna_decode_cache` exists and
holds `amr_portal`, `mpra`, `mpra_human`, `promcalc`, `promoteratlas_positional` but is **missing
the required `refseq`**, so the probe correctly reports `fit: false` and refuses to substitute — yet
the artifact still writes `roots_missing: 0`, so **it does not announce its own degradation**. The
orphan-signal test is the only thing that catches it, and the orphan signal is the inventory's whole
purpose. Restore from HEAD rather than commit. Left untouched: another session's in-flight change.

### Measured 2026-11-28 — the leak DiD at 500 windows (the arms ledger is UNCHANGED at 0.6535)

The clean half had 519 windows available and only 120 were scored. Re-scored all three checkpoints
on both halves at 500 (six runs, ~25 min, always with explicit `--out`; the committed 11-22 and
11-24 artifacts verified byte-untouched afterwards).

- **`run2_reseed_seed1` is the CONTROL, not a second re-pretrain** (`ckpt_dir
  data/glm_ckpt_reseed_control`, a reseed on the corpus that still CONTAINS the draft). Nothing in
  its label says so. A first pass misread it and reported a cross-seed sign disagreement that does
  not exist; retracted in `800de831`.
- Treatment DiD **−0.000844** (left_only) / **−0.001970** (bidirectional). Floor, measured at the
  same n, **0.000824** / **0.000922** — down from 0.00281 at 120 windows.
- The treatment now **EXCEEDS its floor in both arms**, which it did NOT at 120 windows. The
  grader's gate still fails on its CI half, so the verdict stays **LEAK_NOT_IDENTIFIABLE**.
- The pre-registered expected **−0.0054 is now excluded** in both arms — the 120-window run could
  not establish that. Its own provenance was NOT re-derived, so this is a statement about the
  recorded expectation, not a re-measured leak. One control draw is still one draw of trajectory
  noise.
- Artifacts: `wiki/glm_leak_did_n500_synthesis_2026-11-28.json` + four per-cell + six scored.
