# AMR Portal new-cell exploratory curation — cipro/QRDR transfer (2026-06-28)

Tier 3/4 of the unscored-cell triage (`wiki/amr_portal_unscored_triage_2026-06-28.md`): per-cell curation,
not a sweep. Scored **ciprofloxacin via the conserved `qrdr_point` rule** (≥2 gyrA/parC/parE target-site
point mutations — the one drug whose rule legitimately transfers across Gram-negatives) on 7 candidate new
organisms from the AMR Portal (free, provenance-disjoint, measured AST). Scorer:
`scripts/amr_portal_explore_newcells.py`. Artifact: `wiki/amr_portal_newcell_exploratory_2026-06-28.json`.
EXPLORATORY — these organisms have NO deployed claim, so they are NAMESPACE-SEPARATE from the deployed
`amr_portal_independent` card (scope contract + shared-key trap). FROZEN AMR surface untouched (read-only).

## Result + the key finding

| cell | nR / nS | sens | spec | acc | verdict |
|---|---|---|---|---|---|
| Neisseria gonorrhoeae | 5618 / 6406 | 0.846 | 0.993 | 0.924 | **TRANSFERS_CLEAN** |
| Staphylococcus aureus | 1563 / 2095 | 0.944 | 0.892 | 0.914 | **TRANSFERS_CLEAN** |
| Acinetobacter baumannii | 1276 / 221 | 0.904 | 0.995 | 0.918 | **TRANSFERS_CLEAN** |
| Pseudomonas aeruginosa | 797 / 434 | 0.577 | 0.991 | 0.723 | TRANSFERS_PARTIAL |
| Enterobacter cloacae | 87 / 48 | 0.264 | 1.000 | 0.526 | NEEDS_CURATION |
| Enterobacter hormaechei | 29 / 52 | 0.000 | 1.000 | 0.642 | NEEDS_CURATION |
| Serratia marcescens | 14 / 36 | 0.000 | 1.000 | 0.720 | NEEDS_CURATION |

**THE FINDING (inverts the naive prior — an R2-discipline catch): cipro-QRDR transfer is MECHANISM-specific,
NOT taxonomy-specific.** The expectation going in was "Enterobacterales transfer (E. coli-like QRDR);
Gram-positive / efflux organisms don't." The data says the OPPOSITE on both counts:
- **It transfers wherever cipro resistance is actually driven by ≥2 QRDR target-site point mutations** —
  Neisseria gonorrhoeae, Staphylococcus aureus (Gram-positive!), Acinetobacter baumannii all score cleanly.
- **It fails wherever cipro R runs through a DIFFERENT mechanism** — Pseudomonas (MexAB efflux → sens 0.58,
  the efflux isolates lack 2 QRDR mutations), and the Enterobacterales (Enterobacter/Serratia: qnr plasmid
  genes + AmpC-derepression + single-mutation R → qrdr_point catches almost none; sens 0.0–0.26).
- Spec is ~1.0 everywhere → when 2 QRDR mutations ARE present the call is right; the failures are all
  UNDER-calls (missed R), never over-calls. The rule's precision holds; its recall is mechanism-bounded.

## Honest disposition
- **3 new TRANSFERS_CLEAN powered cells** (Neisseria gonorrhoeae / Staphylococcus aureus / Acinetobacter
  baumannii cipro) — independently validated on free measured AST. These are PROMOTION CANDIDATES: adding
  each as a `<organism>|ciprofloxacin` calibrated-registry entry is a USER-ratified frozen-surface amendment
  (the `calibrated_amr_rules.json` is sha-pinned) — NOT auto-promoted here. Staph is Gram-positive, so confirm
  the AMRFinder QRDR-determinant identity (gyrA/grlA) before promotion (the number is clean — sens 0.944 /
  spec 0.892 discriminates, not a degenerate all-R call — but identity-confirm is good hygiene).
- **NEEDS_CURATION (Enterobacter ×2, Serratia)** — the conserved rule does NOT transfer; their cipro R is
  qnr/AmpC/efflux/single-mutation. An organism-specific rule would need qnr-gene + single-QRDR logic,
  developed + validated on a SEPARATE cohort. **NOT done by lowering the threshold on THIS cohort** (that
  would overfit the independent test set — the threshold-vs-null trap). Honest ABSTAIN-until-curated.
- **Pseudomonas cipro PARTIAL** — qrdr_point is a high-precision floor (spec 0.99) but recall-bounded by
  efflux; an organism-specific efflux-aware rule is the curation path, or accept the partial floor.

## Scope + integrity
- EXPLORATORY, not deployed; invisible to the frozen report card by design; namespace-separate artifact.
- FROZEN AMR surface (`amr_rules.py` + `calibrated_amr_rules.json`) byte-unchanged; leak guard 9/9.
- Acquired-gene drugs (tet/cef/gent/mero) on these organisms were deliberately NOT scored — the
  intrinsic-gene guardrail (their determinants vary with each organism's intrinsic flora; cipro/QRDR is the
  sole conserved-mechanism transfer).
- Tests: `tests/test_amr_portal_explore_newcells.py` (pure verdict-thresholds; synthetic).

## Companion
`wiki/amr_portal_unscored_triage_2026-06-28.md` (the board), `wiki/frontier_reassessment_2026-06-28.md`
(the corrected frontier), `wiki/amr_portal_independent_report_card.md` (the deployed-cell card, untouched).
