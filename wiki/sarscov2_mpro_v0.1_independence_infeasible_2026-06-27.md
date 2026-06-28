# SARS-CoV-2 Mpro v0.1 independent number — FREE PATH EXHAUSTED (verified 2026-06-27)

**Conclusion: the v0.1 "independent number" for the SARS-CoV-2 Mpro cell cannot be produced from the
free CoV-RDB surface. The catalog and ALL Mpro-inhibitor fold-change records come from the SAME
studies, so there is NO held-out partition to score. An independent number is an EXTERNAL wall
(non-CoV-RDB clinical-isolate fold), NOT a code or search-effort problem.** Same structural pattern as
the TB gold-set exhaustion (`wiki/tb_goldset_public_source_exhaustion_2026-06-22.md`): the binding
constraint is the independent LABEL, not the pipeline.

## What v0.1 needed
The v0 number is `COV_RDB_IN_DISTRIBUTION_KNOWLEDGE_BASELINE` — the Mpro catalog
(`sarscov2_amr.MPRO_MAJOR_DRMS`) is built from CoV-RDB `invitro_selection_results`, and the validation
fold-change is ALSO CoV-RDB `rx_fold`. A truly-independent number needs fold from a source NOT used to
build the catalog (the honesty rail in `scripts/sarscov2_mpro_validate.py`). The cheapest free attempt:
a **held-out-by-study** partition — score the frozen catalog on fold records whose `ref_name` is disjoint
from the selection studies that built the catalog.

## The census (verified against the committed CoV-RDB payload)
`covid-drdb-payload` (D: cache, 11M). Mpro selection studies (catalog provenance): 10 refs incl.
`Bouzidi23, FDA23-NTV, Heilmann22, Iketani22c, Krismer23, Takashita23, Zhou22c` (nirmatrelvir).
`rx_fold` Mpro-inhibitor records by study:

| drug | fold records | studies | held-out studies (ref ∉ selection) |
|---|---|---|---|
| nirmatrelvir | 49 | FDA23-NTV (42), Krismer23 (7) | **0** |
| ensitrelvir | 7 | Krismer23 (7) | **0** |
| lufotrelvir | 0 | — | **0** |

**Every Mpro-inhibitor fold record comes from a study that is already in the catalog-building selection
set.** The held-out partition is EMPTY across all three drugs. (The `rx_fold` table is dominated by
monoclonal-antibody drugs — CB6/Vir-7831/REGN… — which are irrelevant to the Mpro target-site cell.)

## Why this is not "try the next source"
CoV-RDB (`hivdb/covid-drdb-payload`) IS the field's aggregator of published SARS-CoV-2 resistance
fold-change; it already pulls in the public Mpro-inhibitor in-vitro selection + fold studies. There is no
*free, public, per-isolate, MEASURED Mpro-inhibitor fold from a study CoV-RDB has not already absorbed*.
An independent number requires clinical-isolate nirmatrelvir fold OUTSIDE CoV-RDB (author contact / a
clinical cohort / a DUA) — the acquisition wall.

## Disposition
- **v0 stands** as the honest in-distribution baseline (sens 0.68 / spec 0.0, UNDERPOWERED 37R/5S,
  `wiki/sarscov2_mpro_validation_result_2026-06-23.md`). Do NOT relabel it independent.
- **v0.1 = BLOCKED:external** (acquisition), with proof — joins TB on the label wall. The
  `sarscov2_mpro_validate.py` scorer is ready the instant an external fold table lands (it already
  partitions by `ref_name` candidates and censors operators) — the blocker is exclusively the label.
- The other deferred SARS items (RdRp/remdesivir nsp12 with its −1 ribosomal frameshift; per-drug
  differential) are separate, heavier, and also need data not on the free surface.

## Artifacts
Census reproducible from `scripts/sarscov2_mpro_validate.py::load_records` + the committed payload.
Companion walls: `wiki/tb_goldset_public_source_exhaustion_2026-06-22.md`,
`wiki/negative_results_map_2026-06-13.md`, `wiki/reproducibility_freeze_2026-06-13.md`.
