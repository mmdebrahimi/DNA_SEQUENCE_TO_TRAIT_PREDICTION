# Data-acquisition sweep — the 4 free-independent-measured-label reservoirs (2026-06-25)

Ran the 4 reservoirs in sequential order to get a **new free independent measured-label validation number**
(the "labels not models" wall). Outcome: **reservoir 1 is a clear WIN** (a real independent number, the first
for a non-AMR typing trait). Reservoirs 2–4 are gated/closed (documented walls, not unexplored).

## Reservoir 1 — Phenotypic serology typing → **WIN**
The free independent measured label was hiding in plain sight: the **GPS pipeline paper (Nat Commun 2025)
Supplementary Data 1** carries per-isolate **`Phenotypic_serotype`** (method QUELLUNG/phenotypic) + `ERR`/`ERS`
accessions for **260 Poland isolates** (32 serotypes), and **ENA serves the GPS-deposited assemblies**
(`ERS → analysis ERZ → contig.fa.gz`). Both the label (wet-lab serology) and the assembly (GPS-deposited) are
**independent of our caller** → clears the circularity rail (this is NOT the in-silico Monocle field).

- Runner: `scripts/pneumo_gps_quellung_validate.py` (native blastn, no Docker, checkpointed; contigs deleted per-isolate).
- **Interim n=32: serogroup concordance 0.969, exact-serotype 0.688** (full n≤260 accruing). exact misses are
  WITHIN-serogroup (9A↔9V, 6B↔6E, 15B↔15C) = the documented v0 single-best-reference ceiling, not a bug.
- **Honest headline: serogroup ~0.97 is the v0's real independent resolution; exact ~0.69 is the lower bound
  motivating a v0.1** (allele-level within-serogroup typing). Artifact: `wiki/pneumo_serotype_cohort_validation.json`.
- **This is the answer to "the data we need":** a genuine independent-measured-label number for a freshly-built cell.

**Bonus discovery:** the SAME Supplementary Data 1 also carries **measured AST** (broth-microdilution MIC +
S/R) for ~16 antibiotics (penicillin, ceftriaxone, meropenem, erythromycin, levofloxacin, tetracycline, …)
on these pneumococcal isolates — a NEW free measured **AMR** phenotype source for a NEW organism (S. pneumoniae
is not currently an AMR cell). A candidate future cell, screened-clean of the bacterial-AMR gates by being a
new organism on measured MIC.

**Salmonella half** (the other reservoir-1 candidate): still gated on building the real SeqSero2 antigen DB
(the data-engineering documented in `wiki/salm_serovar_report_card.md`) + a wet-lab-serotyped (not SISTR-
predicted) cohort. Deferred behind the pneumo win.

## Reservoir 2 — TB CRyPTIC measured MIC → gated (in-distribution baseline already exists)
CRyPTIC BMD-MIC is a real measurement, and the TB cell already has its **in-distribution baseline**
(`wiki/tb_cryptic_parquet_baseline_2026-06-22`: RIF raw 0.916/0.974, lineage-collapsed 0.41; INH similar).
The remaining gap is **independence**: the WHO catalogue (the rule) was built partly from CRyPTIC, so scoring
it on CRyPTIC is in-distribution. A near-independent number needs a **hand-curated post-2023 gold set**
(`wiki/tb_independent_goldset_acquisition_2026-06-17.md`) — NOT autonomously available tonight. No new free
data from this reservoir without that curation. **Status: blocked on gold-set acquisition.**

## Reservoir 3 — Viral measured-phenotype DBs (HCV/HBV) → no free consolidated measured dataset
The HIV win came from Stanford HIVDB publishing the actual **PhenoSense measurements** as a downloadable
genotype↔measured-phenotype dataset. There is **no equivalent free consolidated measured-phenotype dataset
for HCV or HBV**: the available resources (**geno2pheno[HCV]/[HBV]**, HCV-GLUE) are **rule/annotation tools**
(literature-curated) = faithful-to-tool, not independent measured labels. HCV/HBV replicon-assay fold-change
data exists but is scattered across primary DAA-development papers, not a single free downloadable repository.
**Status: a HCV/HBV cell would be faithful-to-tool only (not the independent win); skip unless a measured
dataset surfaces.** (Primary-literature mining is a large manual effort with uncertain payoff.)

## Reservoir 4 — Bacterial AMR MIC (NCBI-PD/BV-BRC) → closed for public expansion
The 10 SCORED provenance-disjoint cells + the **8 rejection gates** (`wiki/negative_results_map_2026-06-13.md`)
establish public bacterial-AMR R/S is saturated. Only **non-public acquisition** (collaborator/biobank/clinical
wet-lab AST) reopens it — a USER decision, no autonomous data. **Status: closed (gates), pending acquisition.**

## Verdict
| Reservoir | Outcome | New free independent number? |
|---|---|---|
| 1 Phenotypic serology (pneumo Quellung) | **WIN** | **YES** — serogroup 0.97 / exact 0.69 (interim n=32; full n≤260 accruing) |
| 2 TB CRyPTIC | gated | No (independent gap needs a hand-curated post-2023 gold set) |
| 3 HCV/HBV | no measured dataset | No (geno2pheno = rule tool; faithful-to-tool only) |
| 4 Bacterial AMR | closed | No (8 gates; non-public acquisition only) |

**We got what we needed from reservoir 1.** The next-best *new* data leads, both surfaced here: (a) the
**pneumococcal measured-AST** bonus in the same GPS supplement (a new-organism AMR cell), and (b) the standing
**non-public acquisition** path (reservoir 4) that the user controls.
