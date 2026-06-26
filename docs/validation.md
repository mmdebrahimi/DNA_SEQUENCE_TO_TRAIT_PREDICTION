# Validation — what the tiers mean and where the numbers come from

`dna-decode` attaches a **validation tier + provenance** to every call. Tiers are **not averaged** into a
single score — a surface validated against a free independent wet-lab label is a different claim from one
that is faithful to a reference tool, and the tool keeps them distinct on purpose. The authoritative,
always-current per-cell view is `dna-decode list` and the report cards under `wiki/`.

## Tier glossary

| Tier | Meaning |
|---|---|
| `INDEPENDENT_MEASURED` | Scored on a free, **independent measured-AST** cohort (submitter/lab/country provenance-disjoint; non-circular). The strongest bacterial-AMR tier. |
| `INDEPENDENT_WETLAB` | Scored on a free, **independent isolate-level wet-lab** label (e.g. Stanford HIVDB PhenoSense fold-change — *not* the tool's own interpretation, so non-circular). |
| `near_independent` | Provenance-disjoint stress test or consensus panel (e.g. GeT-RM for PGx) — strong but not a fresh measured cohort. |
| `faithful_to_tool` | Deterministic re-implementation of a curated reference database/method (PlasmidFinder, SerotypeFinder, …); validated that it *matches the tool*, not that it beats an independent label. |
| `knowledge_baseline` | A literature/catalogue rule scored in-distribution (e.g. SARS-CoV-2 CoV-RDB) — informative but not externally independent. |
| `NOT_CENSUSED` / `no_free_source` | CLI-routable but not yet scored / no free isolate-level phenotype label exists. Surfaced explicitly, never hidden. |

## Headline numbers (see the cited cards for the full, caveated tables)

- **Bacterial AMR** (E. coli / Klebsiella / Pseudomonas / S. aureus × cipro/cef/tet/gent/meropenem):
  in-cohort + held-out + cross-source + cross-organism; e.g. cipro acc 0.925 (held-out 0.862), cef 0.933,
  gent 0.945, tet 0.833, mero 0.867. Every per-drug rule beats naive AMRFinder.
  → `wiki/amr_multiorganism_capstone_2026-06-07.md`, `wiki/amr_portal_independent_report_card.md`.
- **M. tuberculosis** (WHO-2023 catalogue rule on the CRyPTIC compendium): rifampicin / isoniazid on
  measured-AST isolates. → `wiki/tb_report_card.md`.
- **HIV-1** (Stanford HIVDB PhenoSense, free independent wet-lab fold-change): NNRTI AUC ≈ 0.96, plus NRTI /
  PI / INSTI / CAI surfaces at their own honest tiers. → `wiki/hiv_decoder_report_card.md`.
- **C. auris fungal AMR** (de-confounded WGS+MIC, Gate G1): fluconazole sens 1.0 across clades,
  label-limited specificity. → `wiki/fungal_ep7_g1_closeout_2026-06-08.md`.
- **One legible cross-kingdom view** (tiers preserved, never averaged): `wiki/cross_kingdom_validation_summary.md`.

## Honest limitations (the tool reports these per call — preserve them)

- An **S call cannot rule out** resistance via efflux / porin-loss / regulatory mechanisms absent from the
  curated determinant databases.
- Raw per-isolate sens/spec on clonally-dominated R classes is **clonality-inflated**; the report cards
  disclose lineage-effective N + cluster-weighted metrics with confidence intervals.
- Typing/finder calls are **faithful-to-tool**, not independent baselines.
- The learned foundation-model (embedding) track is a **closed negative** — the shipped product is the
  deterministic decoder. See `wiki/negative_results_map_2026-06-13.md`.

This file is a human summary. The machine-readable evidence is the committed `wiki/*report_card*.json`
and the Evidence-Contract Registry (`dna_decode/data/cell_registry.py`).
