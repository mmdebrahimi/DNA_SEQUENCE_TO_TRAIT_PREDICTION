# Free viral genotype–phenotype susceptibility datasets for antiviral-resistance interpretation — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-06-21. Source: Claude (`/research` WebSearch synthesis). Slug: viral-antiviral-resistance-gp-datasets-2026-06-21.
> Audit floor: 5 of 5 locators per row. Mapping floor: rationale → quantity. Banned-phrase + cite-token + source-identity scans applied.
> **Degraded web layer:** broad multi-pathogen WebSearch queries were blocked by a usage-policy filter this session; HIVDB pages are JS-rendered (WebFetch read only the title). ALL rows below derive from WebSearch summary synthesis (not a direct page fetch) → every row carries a websearch-summary provenance caveat; HCV/SARS-CoV-2 are honest gaps.

## Research Context

- **Problem:** free, redistributable, isolate-level laboratory genotype–phenotype susceptibility datasets for antiviral resistance interpretation (HIV RT/protease/integrase, HCV, SARS-CoV-2 protease), screened against the project's 8 rejection gates + license terms
- **Why:** Wave B of the bacteria/virus phenotype→trait tool — the project's viral determinant decoder (influenza NA) exists in code but is unvalidated for lack of a free isolate-level label. This run is the GO/NO-GO on a free, de-confounded viral genotype↔phenotype substrate.
- **Captured:** 2026-06-21
- **Schema:** memo-schema 0.4

## Audit table (verbatim, supported rows only)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Stanford HIVDB has a public, isolate-level genotype-phenotype dataset measured by a wet-lab assay | qualitative | — | Genotype-phenotype datasets — HIVDB | Stanford HIVDB | 2026 | dataset page | https://hivdb.stanford.edu/pages/genopheno.dataset.html | 2026-06-21 | "in vitro susceptibility tests were performed using the PhenoSense assay" | A public download page whose phenotype is a wet-lab in-vitro susceptibility measurement | primary (DB) | medium |
| HIVDB genotype-phenotype dataset — isolates | 2,167 | isolates | study citing HIVDB gp-dataset | via HIVDB | 2026 | dataset description | https://hivdb.stanford.edu/pages/genopheno.dataset.html | 2026-06-21 | "12,442 phenotype results from 2,167 isolates" | Isolate-level genotype↔phenotype pairs at usable scale | secondary (search synthesis) | medium |
| HIVDB genotype-phenotype dataset — phenotype results | 12,442 | phenotype results | study citing HIVDB gp-dataset | via HIVDB | 2026 | dataset description | https://hivdb.stanford.edu/pages/genopheno.dataset.html | 2026-06-21 | "12,442 phenotype results from 2,167 isolates" | Many drug × isolate fold-change measurements | secondary (search synthesis) | medium |
| Phenotype = PhenoSense fold-change, independent of any genotype-interpretation algorithm | qualitative | — | HIVDB genotype-phenotype | Stanford HIVDB / Monogram | 2026 | assay description | https://hivdb.stanford.edu/pages/genotype-phenotype.html | 2026-06-21 | "in vitro susceptibility tests were performed using the PhenoSense assay (Monogram ...)" | Label is a wet-lab measurement, NOT a genotype rule → clears the circular-label gate | primary (DB) | medium |
| HIVDB public sequence corpus scale (cohort-construction context) | 119,000 PR / 128,000 RT / 13,000 IN | sequences | HIV-1 PR/RT/Integrase Variation (Rhee et al.) | Rhee et al. / J Virol | 2016 | results | https://hivdb.stanford.edu/_wrapper/pages/pdf/Rhee.2016.JVI.pdf | 2026-06-21 | "119,000 PR, 128,000 RT, and 13,000 IN sequences from 132,000 individuals in 143 countries" | Large public sequence base; 85% in GenBank | primary (paper) | medium |

## Source-Locator Coverage

- Total rows submitted: 6
- Survived audit floor: 6
- Survived mapping floor: 6
- Survived banned-phrase scan: 6
- Final supported: 5 (1 row dropped on low-confidence: the license/terms-page row)
- Survival rate: 5 / 6 (83%)

## Caveats per row

- **All rows:** `provenance: websearch-summary (re-verify against the direct HIVDB page before high-confidence use)` — the HIVDB pages are JS-rendered and were not directly WebFetched; rows derive from WebSearch synthesis. None were rated `high`, so no confidence downgrade was applied, but a direct-page confirmation is the discipline before any uplift.
- **Rows 2-3 (2,167 / 12,442):** the exact counts come from a study citing the dataset, not the dataset page's own header — confirm against the live download.

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| HIVDB has a public isolate-level genotype-phenotype dataset (PhenoSense fold-change) | 2,167 | isolates | https://hivdb.stanford.edu/pages/genopheno.dataset.html | **Candidate use:** the substrate for the FIRST validated viral cell — validate the HIV RT/protease/integrase determinant decoder against an INDEPENDENT wet-lab fold-change label (clears the circular-label + sampling-defined gates that bound the bacterial learned arm). **Verification needed:** download the actual dataset + confirm it pairs sequences with fold-change per drug; confirm isolate count + drug coverage on the live file. | medium |
| Phenotype is a PhenoSense lab fold-change, independent of genotype rules | qualitative | — | https://hivdb.stanford.edu/pages/genotype-phenotype.html | **Candidate use:** this independence is the decisive property — it is why HIV clears the circular-label gate (unlike AMRFinder-derived bacterial labels). **Verification needed:** confirm the fold-change is the deposited label (not a derived/predicted susceptibility), and capture the assay's lower/upper dynamic-range cutoffs (censoring → reuse `MicValue`). | medium |
| HIVDB data is permissively licensed (publications CC BY 4.0; data often CC0) | qualitative | — | https://hivdb.stanford.edu/_wrapper/pages/pdf/Rhee.2016.JVI.pdf | **Candidate use:** confirms redistribution is likely permitted → the data can ship as a committed validation cohort. **Verification needed (LOAD-BEARING):** the authoritative HIVDB dataset-download terms-of-use + citation page (NOT the associated papers) — this row is `low` confidence and is the single must-verify item before treating the data as redistributable. | low (see unsupported) |
| HCV (NS5A/NS3/NS5B) + SARS-CoV-2 (Mpro) free genotype-phenotype datasets | qualitative | — | (blocked) | **Candidate use:** secondary viral substrates if HIV proves the pattern. **Verification needed:** re-run as separate narrow domain-restricted `/research` queries — every broad query was usage-policy-blocked this session (tooling gap, NOT absence of data). Candidate sources: geno2pheno[hcv], Stanford CoVDB. | — (honest gap) |

additional candidates exist; review the full memo + the raw memo's 8-gate screen.

## Verification trace (Mission Control L1)

This intake was invoked as part of Mission Control run `2026-06-21-1314-research-viral-antiviral-resistance-gp-datasets`. The parent Intent Contract is at `mission-control-runs/2026-06-21-1314-research-viral-antiviral-resistance-gp-datasets/intent-contract.md`.

**Validation steps applied:**
- Audit floor (Step 2): 6 pass / 0 fail
- Mapping floor (Step 3): 6 pass / 0 fail
- Banned-phrase scan (Step 4): 0 hard-reject / 0 soft-warn
- Cite-token noise scan (Step 5): 0 flagged
- Source-text identity advisory (Step 5.5): 6 provenance-flag (all websearch-summary) / 0 author-identity-uncertain / 0 quote-shape-table-cell

**Verification result for parent run's sub-task "Intake validation":**
- Status: PASS
- Criterion: rows pass audit + mapping floors + banned-phrase + cite-token + source-identity advisory
- Evidence: this memo (5 supported rows) + `research_outputs/viral-antiviral-resistance-gp-datasets-2026-06-21_unsupported.md` (1 rejected row)

## Promotion Gate reminder

INPUT to the 4-step Promotion Gate, NOT an approval. Before treating HIVDB as a committed validation cohort: (1) the dataset resolves + downloads at the cited URL; (2) the genotype↔phenotype pairing + fold-change label is confirmed on the live file; (3) the authoritative terms-of-use/license is read directly; (4) the project's de-confound discipline (within-subtype concordance) is applied — overall AUROC on HIV will conflate subtype structure + mechanism exactly as it did for bacterial lineage.
