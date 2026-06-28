# Cross-kingdom validation summary — the whole validated decoder surface

One legible view of every validation surface. **No aggregate headline** — each surface has a DIFFERENT independence construction + honesty tier, preserved verbatim (never averaged). See each card for per-cell detail.

| Kingdom | Surface | Headline | Independence tier (honest) | Card |
|---|---|---|---|---|
| Bacteria | NCBI-PD provenance-disjoint (frozen) | 10 SCORED / 27 cells | isolate-level provenance-disjoint (different submitter/lab/country); NOT methodology-independent | `wiki/decoder_validation_report_card.md` |
| Bacteria | EBI AMR Portal independent (measured AST) | 23 SCORED_INDEPENDENT / 27 cells | BioSample/GCA disjoint vs CRyPTIC + tuning cohorts; measured wet-lab AST (non-circular); E.coli/Salmonella/Klebsiella/Shigella 0.83–0.995 acc | `wiki/amr_portal_independent_report_card.md` |
| Bacteria (M. tuberculosis) | EBI AMR Portal INDEPENDENT (WHO rule) | RIF acc 0.937; ISO acc 0.914 (N~2,845, full cohort) | BioSample-resolution-checked (ENA-side biosample-grade; NCBI-side 0/30 cross-archive); RAW headline (homoplasy -> lineage is disclosure); measured AST | `wiki/tb_report_card.md` |
| Virus (HIV) | Stanford HIVDB PhenoSense | 25 cells (NNRTI/NRTI/PI/INSTI/CAI) | PhenoSense fold-change is NOT HIVDB's own Sierra interpretation (non-circular) | `wiki/hiv_decoder_report_card.md` |
| Bacteria (external cohort) | external-cohort revalidation arm | 2 cells | external clinical re-validation of the frozen decoder; strict-tier is the primary metric, relaxed secondary; raw sens/spec is clonality-inflated — see the cluster-weighted block. Separate from the frozen decoder report card. | `wiki/external_validation_report_card.md` |

## The arc (what this surface represents)
- **The binding constraint of the project — a FREE, independent, measured-phenotype label — is broken across bacteria AND M. tuberculosis.** The deterministic decoder is independently validated on E. coli / Salmonella / Klebsiella / Shigella (EBI AMR Portal, measured AST) **and** M. tuberculosis (WHO rule, N~2,845), all FREE (no DUA, no author-contact).
- **HIV** is validated against a free independent wet-lab fold-change (Stanford HIVDB PhenoSense) — the project's first independent-label win.
- **SARS-CoV-2 / influenza / fungal** cells are in-distribution or no-free-phenotype-source (their honesty tiers say so on their own cards).
- **Learned-embedding expansion is a CLOSED 0-for-4 negative** (cipro within-lineage, pathotype, Arabidopsis ×2) — `wiki/negative_results_map_2026-06-13.md`. The validated shippable artifact is the DETERMINISTIC decoder suite, not the foundation-model embedding bet.

## Honesty discipline (preserved, not flattened)
Each surface's independence is a DIFFERENT construction: NCBI-PD provenance-disjoint (submitter/lab/country) ≠ EBI AMR Portal accession-disjoint measured-AST ≠ TB BioSample-resolution-checked ≠ HIV free wet-lab fold-change ≠ external-cohort measured-MIC. This summary lists them side by side; it does NOT compute a single cross-kingdom accuracy (that would be a category error).

Rebuild: `uv run python scripts/build_cross_kingdom_summary.py`.
