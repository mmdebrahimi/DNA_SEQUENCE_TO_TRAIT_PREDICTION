<!-- recommendation.md for 2026-06-03-0000-ep4-v0-cohort-calibrate -->
# Recommendation — next moves (for the user on waking)

## Done overnight (committed + pushed)
v0 resolver cohort-evaluated + ExPEC-calibrated to H4 ship-gate (precision 1.0, abstention 0.08, EPEC recall 1.0, ExPEC recall 0.75). Coverage cache makes future resolver tuning instant. VF diff documented as gated (needs aligner binaries).

## Ranked next actions
1. **VF side-by-side diff (gated — needs you)** — install BLAST+ or KMA (see `research_outputs/pathotype_vf_sidebyside_feasibility_2026-06-03.md`), then I can wire `scripts/` to run real VirulenceFinder + diff its calls vs our `marker_hits`. This is the last open ledger v0 requirement.
2. **Per-gene ExPEC scoring (v0.1)** — split SIDEROPHORES/CAPSULE clusters into per-gene presence to recover the screen's 0.882-style granularity and push ExPEC recall past 0.75 without touching precision. Cheap now (coverage cache + ~5 min re-detect with a finer catalog).
3. **Expand the cohort to ETEC** — add von Mentzer ETEC genomes (8 LR refs in `research_outputs/etec_vonmentzer_collection_gca_2026-05-30.csv`) to exercise the ETEC arm of the 11-class table + report 3-class (ExPEC/EPEC/ETEC) supported-surface metrics.
4. **Break study==class** — the standing confound: get within-study or ST/lineage-matched genomes so a learned comparison is interpretable. Not needed for the resolver (it's confound-immune) but needed before any learned-vs-classical claim.

## EMIT (user-only `/project-state`)
```
/project-state ecoli-pathotype-prediction-cli-2026-05-26 --append-action --class build \
  "v0 resolver cohort-evaluated (H4): confident precision 1.0, EPEC recall 1.0; ExPEC calibrated 0.25->0.75 recall / 0.58->0.08 abstention via LOW_CONFIDENCE ExPEC_COMPATIBLE (1 strong + 2 support). 20/20 tests. Both H4 ship-gate targets met. VF side-by-side diff gated (needs KMA/BLAST+). Advances Goal 5."
```

## Held
- VF dep-install — declined unattended (fail-safe). All pure-Python work done + pushed.
