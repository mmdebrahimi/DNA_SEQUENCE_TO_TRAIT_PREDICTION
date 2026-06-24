# M. tuberculosis — validation report card

Standing TB validation surface, **namespace-separate** from the frozen NCBI-PD / AMR-Portal / HIV / external cards (TB is the non-frozen `organism_rules` cell; its phenotype is BMD-MIC / genomic-VCF and its independence is constructed differently — it must never be keyed into the frozen `canonical_cell_key` card). **Headline = RAW per-isolate sens/spec** (TB resistance is HOMOPLASIC, so a lineage-majority collapse measures the wrong question; the lineage figure is a clonality DISCLOSURE, never the headline — mirrors the bacterial disclose-not-demote discipline).

## INDEPENDENT — provenance-disjoint, measured AST (the gold-set saga's resolved holy grail)
Free (no DUA, no author-contact). WHO-2023 catalogue rule applied UNCHANGED; isolates the rule was never tuned on.
| Drug | n (R/S) | RAW sens (95% CI) | RAW spec (95% CI) | RAW acc | lineage disclosure (NOT headline) |
|---|---|---|---|---|---|
| rifampicin | 1452/1388 | 0.920 [0.905, 0.933] | 0.955 [0.943, 0.965] | 0.937 | sens 0.444 / spec 0.979 |
| isoniazid | 1643/1199 | 0.879 [0.862, 0.894] | 0.962 [0.949, 0.971] | 0.914 | sens 0.321 / spec 0.972 |

Source: EBI AMR Portal (CABBAGE) provenance-disjoint cohort + assembly→H37Rv-VCF→WHO-rule pipeline (`scripts/run_tb_independent_amr_portal.py`). Memo `wiki/tb_independent_number_2026-06-23.md`; homoplasy/disclosure rationale `wiki/tb_independent_lineage_finding_2026-06-23.md`.

## IN-DISTRIBUTION — knowledge baseline (NOT independent)
The WHO catalogue was built partly FROM CRyPTIC, so a CRyPTIC-scored number is in-distribution. Shown for comparison; the independent number above is the real external test.
| Drug | n | RAW sens / spec | lineage (disclosure) |
|---|---|---|---|
| rifampicin | 8955 | 0.918 / 0.974 | 0.410 / 0.991 |
| isoniazid | 9518 | 0.889 / 0.989 | 0.349 / 0.995 |

Source: `wiki/tb_{rif,inh}_cryptic_parquet_baseline_*.json` (`wiki/tb_cryptic_parquet_baseline_2026-06-22.md`).

## Honesty rails
- **Independence is BioSample-resolution-CHECKED (upgraded 2026-06-23, `wiki/tb_independence_biosample_check.json`).** The ENA-side disjoint isolates (1,364 with an `ERS` accession) are already BioSample-grade — their `ERS` is string-matched DIRECTLY against CRyPTIC's `ENA_SAMPLE` (ERS), the same namespace. The NCBI-side (1,480 `SAMN`) are the only cross-archive risk vs the European CRyPTIC set; a bounded ENA-portal probe found **0/30** of their ENA-mirror accessions in CRyPTIC. The one irreducible residual is genomic RE-SUBMISSION (an isolate sequenced twice as distinct BioSample records) — which needs Mash genomic dedup, NOT accession resolution.
- **Measured phenotype = non-circular** (BMD-MIC / measured DST); **WHO rule applied UNCHANGED.**
- **RAW is the headline; lineage is disclosure** (homoplasy). The independent lineage figures (~0.44 / 0.32) match the in-distribution lineage figures (0.41 / 0.349), confirming the clonal structure of the R classes.
- FROZEN bacterial AMR surface byte-unchanged; this card is READ-only. Rebuild: `uv run python scripts/build_tb_report_card.py`.
