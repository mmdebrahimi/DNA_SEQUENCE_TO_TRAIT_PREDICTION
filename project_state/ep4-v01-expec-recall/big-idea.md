# Big Idea — EP-4 v0.1 ExPEC-recall hardening

> Bounded minister mission (maiden real run, 2026-06-03). Scope: push the shipped v0 pathotype resolver's
> ExPEC recall past its current 0.75 via per-gene ExPEC scoring, WITHOUT regressing confident-call
> precision (currently 1.0). This is a bounded sub-mission of the open-world pathotype project — NOT the
> apex. One gap; the minister generates a candidate family to close it, gates it through interrogation,
> seeds its ledger, and runs --until-mvp to retire the gap.

## Mission Gaps

### gap: improve-expec-recall
<!-- strategy-budget: 2 -->
- test-exit-0 python -m pytest tests/test_pathotype_expec_recall.py
- project-state-row project_state/ecoli-pathotype-prediction-cli-2026-05-26.md:ExPEC recall hardened

## Endpoint Criteria
- test-exit-0 python -m pytest tests/test_pathotype_expec_recall.py
- project-state-row project_state/ecoli-pathotype-prediction-cli-2026-05-26.md:ExPEC recall hardened
