# The independent M. tuberculosis number — the gold-set saga's holy grail (2026-06-23)

The first genuinely-INDEPENDENT TB RIF/INH number in the project's history: the WHO-2023-catalogue rule
scored on a PROVENANCE-DISJOINT, MEASURED-AST cohort — **free, no DUA, no author contact.** Ran on this host
the moment Docker came up; built end-to-end (assembly → H37Rv variant call → WHO rule), not a literature
re-score.

## Result — FULL N=2,845 (provenance-disjoint, measured AST, WHO rule UNCHANGED) ✅ COMPLETE
| Drug | n (R / S) | sens (95% CI) | spec (95% CI) | accuracy |
|---|---|---|---|---|
| **Rifampicin** | 1,452 / 1,388 | **0.920** [0.905, 0.933] | **0.955** [0.943, 0.965] | **0.937** (TP1336 FP62 TN1326 FN116) |
| **Isoniazid** | 1,643 / 1,199 | **0.879** [0.862, 0.894] | **0.962** [0.949, 0.971] | **0.914** (TP1444 FP46 TN1153 FN199) |

Fully powered, tight CIs, INDEPENDENT — the WHO determinant rule on 2,845 isolates it was never tuned on.
**The label wall that walled this number for the entire gold-set saga (5 sources, all author-request / DUA /
circular) is GONE** — the EBI AMR Portal supplied free measured phenotypes + accessions, and the
assembly→VCF→score pipeline supplied the genotype. The number lands right where the in-distribution CRyPTIC
baseline predicted (RIF raw 0.916/0.974, INH 0.889/0.989), as a faithful independent test should — and the
N=60 smoke (RIF 0.917 / INH 0.883) held up at full power (RIF 0.937 / INH 0.914).

### Clonality disclosure (full N — NOT a demoted headline; see the homoplasy finding)
Per `wiki/tb_independent_lineage_finding_2026-06-23.md`, TB resistance is HOMOPLASIC (acquired independently
within sublineages), so a lineage-majority collapse measures the wrong question and is reported as a
DISCLOSURE, not the headline. At full N the lineage-collapsed figures are RIF sens 0.444 (20 R-lineages / 47
S-lineages / 43 discordant) and INH 0.321 (30 / 36 / 44) — matching the in-distribution baseline's lineage
figures (0.41 / 0.349), confirming the clonal structure. **The honest headline is the RAW per-isolate number
above; the lineage figures disclose that the R classes are clonally structured (resistance is not a clonal
property).**

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
- **Full cohort N=2,845 — COMPLETE + fully powered.** (The N=60 smoke was the de-risking proof; it held up.)
- **RAW per-isolate is the honest headline for TB AMR** (not lineage-collapsed) — because TB resistance is
  homoplasic; the lineage figures are a clonality DISCLOSURE (above), mirroring the project's bacterial
  disclose-not-demote discipline.
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
