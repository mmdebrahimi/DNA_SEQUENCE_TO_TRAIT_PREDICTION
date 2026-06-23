# The independent M. tuberculosis number — the gold-set saga's holy grail (2026-06-23)

The first genuinely-INDEPENDENT TB RIF/INH number in the project's history: the WHO-2023-catalogue rule
scored on a PROVENANCE-DISJOINT, MEASURED-AST cohort — **free, no DUA, no author contact.** Ran on this host
the moment Docker came up; built end-to-end (assembly → H37Rv variant call → WHO rule), not a literature
re-score.

## Result — N=60 smoke (provenance-disjoint, measured AST, WHO rule UNCHANGED)
| Drug | n (R / S) | sens (95% CI) | spec (95% CI) | accuracy |
|---|---|---|---|---|
| **Rifampicin** | 40 / 20 | 0.900 [0.769, 0.960] | 0.950 [0.764, 0.991] | **0.917** (TP36 FP1 TN19 FN4) |
| **Isoniazid** | 48 / 12 | 0.854 [0.728, 0.928] | 1.000 [0.757, 1.000] | **0.883** (TP41 FP0 TN12 FN7) |

Sane, strong, and INDEPENDENT — consistent with the WHO catalogue's expected determinant-rule performance,
on isolates the rule was never tuned on. **The label wall that walled this number for the entire gold-set saga
(5 sources, all author-request / DUA / circular) is gone** — the EBI AMR Portal supplied free measured
phenotypes + accessions, and the assembly→VCF→score pipeline supplied the genotype.

## What ran (pipeline, all on this host)
1. **Cohort:** `data/raw/tb_goldset/amr_portal_tb_disjoint_cohort.tsv` (2,845 disjoint isolates with a
   fetchable GCA assembly + measured RIF/INH; the N=60 smoke = the first 60).
2. **Fetch:** GCA assembly FASTA via the NCBI Datasets API (HTTPS — no Docker).
3. **Variant call:** minimap2 `-cx asm5 --cs` to H37Rv NC_000962.3 + `paftools.js call` (Docker biocontainer
   `quay.io/biocontainers/minimap2:2.28`), FILTER `.`→PASS masking (the asm5 confident-difference set).
4. **Score:** `tb_vcf.parse_masked_calls` → `tb_amr.score_drug` (RIF + INH) vs the measured label, via the
   FROZEN WHO-catalogue rule (`organism_rules/tb_amr`, catalogue pinned `0bb39143`).
   - VERIFY-IN-BATCH that proved it: the first R/R isolate (SAMEA1015921) carried **761155 C>T = rpoB S450L**
     (the catalogue's flagship RIF determinant) → RIF=R ✓, INH=R ✓.

Runner: `scripts/run_tb_independent_amr_portal.py` (CHECKPOINTED + restartable; skips fetched/aligned/scored
isolates on resume; per-isolate result → `data/raw/tb_indep/results.jsonl`).

## Honest rails
- **N=60 smoke, not the full cohort.** Underpowered relative to the 2,845 available; the CIs reflect it. The
  full run is staged to D: (C: is at 95% — see below) and is the scale step.
- **Raw sens/spec (not yet lineage-collapsed).** TB R classes are clonally dominated; the honest headline for
  the full number is the lineage-collapsed metric (Mash clustering + `clonality.cluster_weighted_confusion`,
  as the in-distribution baseline already does). The N=60 raw number is the proof, not the final figure.
- **Independent at the ACCESSION level (upper bound).** BioSample/GCA disjoint vs CRyPTIC + our cohorts;
  BioSample cross-archive resolution would only tighten it.
- **Callability unassessed (no regeno VCF):** a determinant non-match = S (documented `tb_amr` behavior).
- **Non-circular:** measured wet-lab AST, and the WHO catalogue applied unchanged (this is a genuine
  genotype→independent-phenotype test, unlike scoring AMRFinder's own TB calls).

## Scale step (full 2,845) — staged to D:
C: is at 95% (13 GB free); 2,845 assemblies ≈ 3.7 GB → the full fetch must stage on D: (4.4 TB free), NOT the
repo disk. The runner is checkpointed, so the full run is `--max 0` with the work dir on D:; it resumes if
interrupted. The lineage-collapsed full number is then the publication-grade independent TB figure.

## Provenance
Runner `scripts/run_tb_independent_amr_portal.py`; cohort `scripts/amr_portal_tb_cohort.py`; WHO rule
`organism_rules/tb_amr` + `data/tb_who_catalogue`. Result `wiki/tb_independent_amr_portal_scores.json`.
Runbook `wiki/amr_portal_tb_independent_runbook_2026-06-23.md`; feasibility
`wiki/amr_portal_feasibility_result_2026-06-23.md`. The in-distribution comparator: RIF raw 0.916/0.974,
INH 0.889/0.989 (`wiki/tb_cryptic_parquet_baseline_2026-06-22.md`) — the independent N=60 lands in the same
range, as it should.
