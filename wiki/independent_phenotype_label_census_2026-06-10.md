# Independent phenotype-label census — can a shipped decoder be validated against a different-LAB label? — 2026-06-10

> Anchor-3 first deliverable (class-(e) decision memo). The shipped deterministic AMR decoder has THREE
> validation layers — BV-BRC held-out (`wiki/dna_amr_external_validation_2026-06-05.md`), NCBI-PD
> cross-SOURCE (`wiki/dna_amr_xsource_validation_2026-06-07.md`), wider-organism transfer
> (`wiki/wider_amr_transferability_synthesis_2026-06-08.md`) — and ALL THREE explicitly name the same unmet
> gap: a **different-lab / different-methodology** check (NARMS / NI-assay / EUCAST). This census asks the
> binding feasibility question before any harness is built: **does a free, machine-readable, isolate-level,
> genome-linkable, methodology-independent + provenance-disjoint phenotype source exist for any shipped
> decoder?**

## VERDICT: NO — not on free, systematic data. The free isolate-level AST ecosystem is integrated.

A truly methodology-independent, provenance-disjoint, free, genome-linkable label **does not exist as a
systematic source** for any shipped decoder. The deepest *free* independence achievable is cross-SOURCE
(already done via NCBI PD), **not** cross-methodology/different-lab.

## Evidence (census table)

| Candidate source | free + machine-readable? | isolate-level phenotype? | genome-linkable? | provenance-DISJOINT from NCBI-PD/BV-BRC? | verdict |
|---|---|---|---|---|---|
| **NARMS Now Integrated Data** (bacterial: Salmonella/Campylobacter/E. coli/Enterococcus) | YES (downloadable spreadsheets, **numeric MIC** + genes) | YES (per-isolate MIC) | YES (genomes uploaded to NCBI) | **NO** | **disqualified on independence** |
| **CDC/GISAID influenza NI-assay IC50** (flu NA decoder) | NO (aggregate/surveillance; GISAID-gated, not free-open) | mostly aggregate | weak | (moot) | **disqualified on availability** |
| Per-study supplementary MIC + accession tables | per-paper, NOT systematic/machine-readable | yes (case-by-case) | yes (case-by-case) | YES (if MIC not submitted as an NCBI antibiogram) | **out of scope** (manual curation, not a census-able source) |

### Why NARMS is NOT provenance-independent (the decisive finding)
- **NARMS WGS isolates are submitted to NCBI Pathogen Detection weekly** ([NCBI PD submit-data](https://www.ncbi.nlm.nih.gov/pathogens/submit-data/)), and the NCBI-PD `AST_phenotypes` column **IS the submitted antibiogram** (MICs / disk diffusion) per [NCBI PD FAQ](https://www.ncbi.nlm.nih.gov/pathogens/faq/). The wider-AMR validation already scores against that exact `AST_phenotypes` field (`scripts/organism_drug_validate.py:58`). So a "NARMS validation" reuses labels from the **same pipeline** as the NCBI-PD validation already done — it is source-distinct in *portal*, not in *measurement provenance*.
- More structurally: **genome-linkable + free ⇒ the genome is in NCBI ⇒ its antibiogram is very likely already in NCBI-PD's `AST_phenotypes`.** The free isolate-level AST ecosystem (BV-BRC ↔ NCBI-PD ↔ NARMS) is one integrated system; you cannot get a free downloadable genome whose phenotype is *also* from a provenance-disjoint lab, systematically.

### Flu NA decoder
Highest scientific value — it has **zero** phenotype validation (it was validated by genotype self-consistency: H275Y presence on real field isolates, `wiki/antiviral_na_4th_kingdom_2026-06-10.md`), and true NAI phenotype is an NI-assay IC50 fold-change (`dna_decode/data/antiviral_amr.py`). But public NI-assay data is aggregate/surveillance-first and the isolate-level fold-changes live behind GISAID (not free-open) — **no free, isolate-level, accession-linked IC50 dataset was found.** (Two targeted web queries were blocked by a usage-policy filter on the assay phrasing; the availability assessment rests on the prior grounding + the aggregate nature of public flu antiviral surveillance.)

## The one salvageable FREE move (clearly NOT gap-closing)
NARMS Now carries **numeric MIC values**, whereas NCBI-PD `AST_phenotypes` carries R/S interpretations. So a NARMS run would enable a **strict-MIC-tier** label (HIGH_R/HIGH_S at a 4× safety margin via `dna_decode/data/mic_tiers.py`) that the NCBI-PD R/S validation could not. That is a **stricter-LABEL, portal-distinct** check — it adds label-quality rigor — but it is **NOT a different-lab / methodology-independent** check (NARMS AST feeds NCBI-PD). It would be incremental over the existing cross-source validation, and must NOT be reported as "independent-lab validation" (that would repeat the just-corrected denominator/independence overclaim).

## Disposition (recommended)
1. **Bank the bound.** The honest ceiling on free validation is: *the deterministic AMR decoders are validated in-cohort + out-of-cohort + cross-SOURCE (NCBI PD); a different-lab / methodology-independent validation is NOT achievable on free, systematic data because the free isolate-level AST ecosystem (BV-BRC/NCBI-PD/NARMS) is integrated.* Propagate this one-line bound into the decoder caveats / README so no future reader over-reads "validated."
2. **Flu phenotype validation = blocked on a free isolate-level NI-assay source** (GISAID-gated). Record as "highest value, not free-source feasible"; revisit only if a free isolate-level IC50+accession dataset appears.
3. **Optional, user's call — NARMS strict-MIC-tier re-validation** (a separate `/idea-anchor`'d build): a stricter-LABEL, portal-distinct check (NOT different-lab). Value is label-quality rigor over the existing R/S validation; cost is a NARMS-spreadsheet→NCBI-assembly linkage build (the `assembly_accession` wall is the risk). Only worth it if "strict 4×-margin MIC labels" is a rigor bar you want; it does NOT close the different-lab gap.

## Bottom line
Anchor 3's premise was right that the gap is real and thrice-named — but the census's answer is that the gap is **not crossable on free systematic data**. That is itself the deliverable: it bounds every "validated" claim the project makes. The deterministic decoders are as independently validated as free data allows (cross-source); going further needs either paid/controlled different-lab assays or per-study manual MIC+accession curation, neither of which is a census-able free source. No harness should be built to chase an independent label that doesn't exist; the strongest available *free* increment is a stricter-label NARMS re-run, explicitly not a different-lab check.
