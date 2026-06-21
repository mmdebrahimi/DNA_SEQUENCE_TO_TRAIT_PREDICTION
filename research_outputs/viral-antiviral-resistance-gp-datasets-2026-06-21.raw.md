# Free viral genotype–phenotype susceptibility datasets for antiviral-resistance interpretation (raw research, 2026-06-21)

> Captured 2026-06-21. Source: Claude Code (`/research` orchestrator, v0.5). Topic: "free, redistributable, isolate-level laboratory genotype–phenotype susceptibility datasets for antiviral resistance interpretation (HIV reverse-transcriptase/protease/integrase, HCV, SARS-CoV-2 protease), screened against the project's 8 rejection gates + license terms". Slug: viral-antiviral-resistance-gp-datasets-2026-06-21.
> Web via Claude Code WebSearch + WebFetch, single-pass. **Thin yield + degraded web layer:** broad multi-pathogen queries (HCV / SARS-CoV-2) were repeatedly blocked by a WebSearch usage-policy filter this session; narrow domain-restricted HIV queries succeeded. HIVDB pages are JS-rendered (WebFetch read only the title). The HIV findings are solid; HCV/SARS-CoV-2 are honest gaps, NOT absence-of-data.

## Audit table (verbatim, all candidate rows)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section / table / figure | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Stanford HIVDB genotype-phenotype dataset is publicly available, isolate-level, lab-measured | qualitative | — | Genotype-phenotype datasets — HIV Drug Resistance Database | Stanford HIVDB | 2026 | Genotype-phenotype datasets page | https://hivdb.stanford.edu/pages/genopheno.dataset.html | 2026-06-21 | "genotype-phenotype datasets ... isolates on which in vitro susceptibility tests were performed using the PhenoSense assay" | Page exists + is a public dataset download page; the phenotype is a wet-lab in-vitro susceptibility assay | primary (DB) | medium |
| HIVDB genotype-phenotype dataset size | 2,167 | isolates | (study citing the HIVDB genotype-phenotype correlation dataset) | via HIVDB | 2026 | dataset description | https://hivdb.stanford.edu/pages/genopheno.dataset.html | 2026-06-21 | "12,442 phenotype results from 2,167 isolates" | Isolate-level genotype↔phenotype pairs at usable scale | secondary (search synthesis) | medium |
| HIVDB genotype-phenotype dataset — phenotype result count | 12,442 | phenotype results | (study citing the HIVDB genotype-phenotype correlation dataset) | via HIVDB | 2026 | dataset description | https://hivdb.stanford.edu/pages/genopheno.dataset.html | 2026-06-21 | "12,442 phenotype results from 2,167 isolates" | Multiple drugs × isolates → many genotype-phenotype pairs | secondary (search synthesis) | medium |
| Phenotype assay = PhenoSense (Monogram) — a fold-change in-vitro susceptibility measurement, independent of any genotype-interpretation algorithm | qualitative | — | HIVDB genotype-phenotype | Stanford HIVDB / Monogram | 2026 | assay description | https://hivdb.stanford.edu/pages/genotype-phenotype.html | 2026-06-21 | "in vitro susceptibility tests were performed using the PhenoSense assay (Monogram ...)" | Label is a WET-LAB measurement, NOT derived from a genotype rule → clears the circular-label gate | primary (DB) | medium |
| HIVDB underlying publications are openly licensed (CC BY 4.0); the dataset's own terms page was not directly confirmed | qualitative | — | HIV-1 Protease, Reverse Transcriptase, and Integrase Variation (Rhee et al.) | Rhee et al. / J Virol | 2016 | copyright/licence statement | https://hivdb.stanford.edu/_wrapper/pages/pdf/Rhee.2016.JVI.pdf | 2026-06-21 | "open-access article distributed under the terms of the Creative Commons Attribution 4.0 International license" | License of the associated PUBLICATIONS is CC BY 4.0; data often under CC0 — but the dataset-download terms page itself was not fetched | primary (paper) | low |
| HIVDB sequence corpus scale (context for cohort construction) | 119,000 PR / 128,000 RT / 13,000 IN | sequences | HIV-1 Protease, RT, and Integrase Variation (Rhee et al.) | Rhee et al. / J Virol | 2016 | results | https://hivdb.stanford.edu/_wrapper/pages/pdf/Rhee.2016.JVI.pdf | 2026-06-21 | "119,000 PR, 128,000 RT, and 13,000 IN sequences from 132,000 individuals in 143 countries" | Large public sequence base; 85% in GenBank | primary (paper) | medium |

## Highest-confidence rows (top 5)

1. Row 1 — HIVDB genotype-phenotype dataset is a real, public, isolate-level, lab-measured dataset (the core GO finding).
2. Row 4 — the phenotype is a PhenoSense lab fold-change INDEPENDENT of any genotype-interpretation rule → clears the circular-label gate that bounded the bacterial learned arm.
3. Rows 2-3 — dataset scale (2,167 isolates / 12,442 phenotype results) is sufficient for cohort construction with a de-confound precondition.
4. Row 6 — large public sequence corpus (context).
5. Row 5 — permissive license signal (CC BY 4.0 / CC0), pending a direct terms-page confirmation.

## Low-confidence rows

- Row 5 (license): inferred from associated journal-article licenses, NOT the authoritative HIVDB dataset terms-of-use page (JS-rendered, not fetchable this session). VERIFY directly before any redistribution claim.

## 8-rejection-gate screen (the specific ask — applied to the HIVDB genotype-phenotype dataset)

| Gate | Verdict for HIVDB gp-dataset | Note |
|---|---|---|
| circular-label (label derived from the tool you'd compete against) | **CLEARS** | phenotype = PhenoSense wet-lab fold-change, independent of genotype rules — the decisive win vs AMRFinder-derived bacterial labels |
| study==class | **CLEARS** | continuous lab measurement, not a study/category axis |
| sampling-defined (label IS the sampling context) | **CLEARS** | lab assay, not a clinical-site/isolation-source category (the gate that killed pathotype) |
| surveillance-domination | likely clears | clinical+lab isolates, not a single surveillance program; confirm submitter spread |
| assembly-attrition | **CLEARS** | short single-gene sequences (PR/RT/IN); no whole-genome assembly attrition |
| MIC-censoring | **PARTIAL** | PhenoSense fold-change has assay dynamic-range cutoffs → operator-aware/censored handling needed (reuse `MicValue`) |
| provenance-not-separable | likely clears | isolate-level with assay provenance |
| dedup-collapses-balance | **NEEDS DE-CONFOUND CHECK** | HIV subtype/clonal structure is the lineage-confound analog → apply within-subtype concordance (the standing discipline) before any AUROC claim |

**Net:** HIVDB clears the LABEL gates that bound the bacterial learned arm; the residual work is the project's standard de-confound precondition (within-subtype) + censoring-aware fold-change handling — NOT a label-availability blocker.

## Honest gaps

- **HCV (NS5A/NS3/NS5B) and SARS-CoV-2 (Mpro/nirmatrelvir) datasets:** NOT assessed — every broad multi-pathogen WebSearch this session was blocked by a usage-policy filter. This is a TOOLING gap (blocked search), not evidence of absence. Re-run as separate narrow domain-restricted queries (e.g. restricted to `hcv.geno2pheno.org`, `covdb.stanford.edu`).
- **Authoritative HIVDB dataset terms-of-use / redistribution license:** not directly confirmed (JS-rendered page). The CC BY 4.0 / CC0 signal is from associated publications. A human/manual check of the official terms + citation page is the one verification-needed item before treating the data as redistributable.
- **Exact file format + whether RT/protease/integrase sequences ship alongside the fold-change values in the download:** not confirmed (page not fetchable); standard HIVDB practice is downloadable TSV/Excel with sequences, but verify.
