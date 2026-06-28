> **⚠️ CORRECTED (2026-06-28, same day).** This scan's headline "no NEW free, public, isolate-level,
> MEASURED-phenotype DB exists beyond the Stanford family" is WRONG — it MISSED the **EBI AMR Portal**
> (`wiki/amr_portal_feasibility_result_2026-06-23.md`), which is exactly that (1.71 M measured-AST rows,
> 74 powered provenance-disjoint cells) and was already in the repo + exploited 5 days earlier. A
> wrong-denominator / under-search miss. The scan's per-candidate gate analysis (FungAMR / MalariaGEN
> verdicts) stands; only the "no free source exists" conclusion is retracted. See
> `wiki/frontier_reassessment_2026-06-28.md`.

# Free independent-label scan — hunting the next HIVDB-equivalent (2026-06-28)

**Goal (Track A):** find a NEW free, public, **isolate-level, MEASURED-phenotype** genotype↔phenotype
database — the pattern that made HIV the project's first genuinely-independent validated cell (Stanford
HIVDB PhenoSense fold-change). One such label is transformative: it converts a built-but-unvalidated caller
into an independently-validated cell. The project found HIVDB + CoV-RDB opportunistically; this is the first
*deliberate* broad sweep across the currently-`NO_FREE_PHENOTYPE_SOURCE` cells (fungal / antiviral /
antimalarial — 11 cells on the report card) plus adjacent pathogens.

**The bar (must clear ALL):** free + public + downloadable; per-ISOLATE (one isolate's genome/sequence ↔ its
OWN measured phenotype); phenotype is a wet-lab/clinical MEASUREMENT (not a genome-derived prediction);
clears the 8 rejection gates (`wiki/negative_results_map_2026-06-13.md`). The decisive distinction is
**isolate-level measured label** (HIVDB pattern) vs **mutation-curation catalog** (WHO-TB/CARD pattern — G1
circular for a determinant caller).

## Candidate scan

| candidate (cell) | source | structure | verdict | gate |
|---|---|---|---|---|
| **Fungal** (C. auris ERG11/FKS1) | **FungAMR** (Nat Microbiol 2025; CARD-integrated) | **mutation→resistance CATALOG** (35,792 entries = gene+mutation+drug+annotation curated from 501 studies; confidence scores). NO per-isolate measured-MIC↔genome table. | **NOT an independent label** | G1 (catalog-vs-catalog if used to validate a determinant caller) |
| **Malaria** (Pf kelch13/crt) | MalariaGEN **Pf6/Pf7** (~7k/20k samples, open) | labels are **resistance INFERRED FROM MARKERS** (pfcrt/kelch13 genotype rules), not measured | **BLOCKED** | G1 (genome-derived label) |
| **Malaria** (measured IC50) | genetic-cross studies (7G8×GB4, HB3×Dd2) + scattered ex-vivo field supplements; WWARN aggregator | real per-sample IC50 exists but is **tiny + lab-cross** (n≈32 progeny) or scattered; WWARN is access-request-gated | **BLOCKED** | G8 (within-class N) + access-gated |
| **Influenza NA** (oseltamivir/zanamivir/peramivir) | WHO GISRS antiviral surveillance / NA-inhibition IC50 | per-isolate NA-inhibition fold exists but is largely **GISAID-gated** (not free-public per-isolate); surveillance reports are aggregate | **NEEDS-VERIFICATION** (search-blocked this session) | likely G7/access; unconfirmed |
| **HCV / HBV** (DAA / NRTI) | geno2pheno[HCV]/[HBV], HCV-GLUE | these are interpretation **TOOLS**; underlying per-isolate measured fold-change typically clinical-not-public | **NEEDS-VERIFICATION** (search-blocked this session) | likely G1/G7; unconfirmed |
| Bacterial beyond E. coli/Kleb | BV-BRC / NCBI-PD / Oxford | already the project's main track (provenance-disjoint + Oxford independent MIC) | ALREADY-EXPLOITED | — |
| Other Stanford-family antiviral DBs | Stanford HIVDB program | only **HIVDB + CoV-RDB** exist; both already exploited | ALREADY-EXPLOITED | — |

## Verdict — decisive negative on the main question (+ 1 adjacent positive)
**No NEW free, public, isolate-level, MEASURED-phenotype genotype↔phenotype database surfaced beyond the
already-exploited Stanford-DB family (HIVDB→HIV, CoV-RDB→SARS).** Each candidate trips a gate by the same
structural pattern the negative-results map predicted:
- the big, open, downloadable resources publish **genome-derived / marker-inferred labels** (MalariaGEN) or
  **curated mutation catalogs** (FungAMR) — both G1-circular for a determinant caller;
- the genuinely-measured phenotype data is either **tiny/lab-only** (Pf crosses, G8) or **access-gated**
  (WWARN, GISAID/influenza) — G7/access, not free-public per-isolate.

This **confirms the banked thesis** (`wiki/reproducibility_freeze_2026-06-13.md`): the binding constraint is
LABELS, and the Stanford antiviral-DB family is a rare exception, not a repeatable well. **For a new
INDEPENDENT cell, ACQUISITION (the user-authority track) is the unambiguous path** — there is no free
shortcut hiding on the open surface.

### The one actionable POSITIVE (adjacent, not an independent cell)
**FungAMR is a high-value CATALOG-ENRICHMENT for the fungal caller** — a comprehensive, curated,
CARD-integrated fungal determinant list (95 species, 246 genes) far larger than the project's hand-curated
`fungal_amr.py` catalog. Integrating it would **improve the fungal ERG11/FKS1 determinant coverage** — but
that is a *caller upgrade*, NOT an independent validation (it cannot score the caller non-circularly). It's
a real, separate, code-closable opportunity (a better catalog), distinct from Track A's goal.

## Honest scope limits
- **2 cells unverified (influenza NA, HCV/HBV):** a content filter blocked every WebSearch phrasing for
  viral-resistance terms this session, and the primary article hosts (Nature/bioRxiv) were paywalled/403.
  These are recorded **NEEDS-VERIFICATION**, not claimed-negative. To close them: check the WHO GISRS
  antiviral-susceptibility data-release terms (influenza) + HCV-GLUE / geno2pheno data-availability
  statements (HCV/HBV) in a session without the filter, against the same isolate-level-measured bar.
- The verified candidates (fungal, malaria) are decisive negatives for an *independent label*.

## Recommended disposition
- **Track A main question: CLOSED (decisive negative)** — no free shortcut to a new independent cell;
  do not keep hunting public DBs as the path to independence. → the lever is **acquisition (Track C)**.
- **New adjacent option surfaced:** FungAMR catalog-enrichment of the fungal caller (a caller upgrade).
- **2 NEEDS-VERIFICATION cells** (influenza, HCV/HBV) parked for a filter-free session.

## Sources (verified 2026-06-28)
- FungAMR: [Nat Microbiol 2025](https://www.nature.com/articles/s41564-025-02084-7) (paywalled) · [bioRxiv 2024.10.07.617009](https://www.biorxiv.org/content/10.1101/2024.10.07.617009v1) · [Springer Communities summary](https://communities.springernature.com/posts/teaming-up-to-investigate-the-landscape-of-antimicrobial-resistance-mutations-in-fungi) · [CARD FungAMR](https://card.mcmaster.ca/fungamrhome)
- MalariaGEN: [Pf7 open dataset](https://www.researchgate.net/publication/367186406_Pf7_an_open_dataset_of_Plasmodium_falciparum_genome_variation_in_20000_worldwide_samples) · [AAT1 CQ-resistance (Nat Microbiol 2023)](https://www.nature.com/articles/s41564-023-01377-z) · [7G8×GB4/HB3×Dd2 IC50 (PNAS)](https://www.pnas.org/doi/10.1073/pnas.0911317106)
- Gates + thesis: `wiki/negative_results_map_2026-06-13.md`, `wiki/reproducibility_freeze_2026-06-13.md`; the HIVDB precedent: `wiki/hiv_decoder_report_card.md`.
