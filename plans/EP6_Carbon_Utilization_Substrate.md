# EP-6 — Carbon-Source Utilization Substrate (the embedding-frontier candidate)

> **Status:** SCOPING 2026-06-06. Infrastructure scaffolded + green; awaiting a real BacDive export to fire the feasibility gate.
> **Anchors on:** `plans/Trait_Decoding_Roadmap.md` Phase 4 (non-AMR bacterial phenotypes); `research_outputs/ecoli-bacterial-phenotype-decoder-substrate-feasibility-2026-06-05.md` (the survey that picked this substrate); `dna_decode/eval/cohort_deconfound.py` (the de-confound precondition).
> **Why now:** AMR + pathotype both FAILED the embedding thesis — AMR lost to mechanism features (cipro NT 0.914 < QRDR-POINT 0.943; within-lineage NT=chance); pathotype labels were sampling-confounded + unvalidatable. Carbon-utilization is the first candidate that clears BOTH failure modes.

---

## Why carbon-utilization (the two-half test the prior substrates failed)

The embedding niche requires BOTH halves; AMR/pathotype each failed one:

| Requirement | AMR (cipro) | Pathotype | Carbon-utilization |
|---|---|---|---|
| **Sampling-INDEPENDENT label** (lab measurement, not a sampling-context category) | ✓ (MIC) | ✗ (isolation source == label) | ✓ (Biolog/API growth assay) |
| **NO curated mechanism catalog to lose to** | ✗ (AMRFinder QRDR-POINT won) | ✓ (no catalog) | ✓ (gene-content RF exists but is no oracle) |

Carbon-utilization is the first substrate that is YES on both. Source: `research_outputs/ecoli-bacterial-phenotype-decoder-substrate-feasibility-2026-06-05.md` (Li et al. 2023, PMC10729968: BacDive 4397 strains × 58 carbon sources; embeddings never tried).

## The pre-named trap (why this could still fail)

Li et al. 2023 found carbon-utilization prediction is **largely phylogenetic** — gene-content RF nails it in-clade but **fails out-of-clade until scale rescues it** (tryptophan out-of-clade 92.2% only at the 4397-strain scale). That is the *exact* lineage-vs-mechanism crux that killed cipro NT (within-lineage = chance → learned lineage, not mechanism). So:

- The de-confound gate (`cohort_deconfound.py`) is a **hard precondition**, not an afterthought.
- The within-lineage diagnostic (`scripts/within_lineage_diagnostic.py`) is the **promotion gate** — same as AMR.
- Ranking favors **balanced/harder** carbon sources (high minority fraction): the easy ones are already near-ceiling for gene-content RF out-of-clade, so embeddings can only add value on the hard ones.

---

## Terminal claim

`pipeline.py predict --phenotype carbon_utilization:<source> --strain-id <E. coli strain>` runs end-to-end AND, on a DE_CONFOUNDED cohort under `leave_one_accession_out` CV, **NT-frozen-mean-pool embeddings beat gene-content (KO presence/absence) by ≥3 pp AUROC on ≥1 carbon source**, with **within-lineage concordance > chance** (the mechanism-not-lineage gate).

A clean FAIL (NT does not beat gene-content, OR within-lineage = chance) is an acceptable terminal too — it would make the embedding thesis **0-for-3** and justify pivoting the decoder to deterministic/gene-content heads permanently. Per north star: honest, failure-tolerant.

## Substrate

- **Organism:** E. coli (first); BacDive is multi-species so an E. coli-only slice is required.
- **Phenotype:** binary carbon-source utilization (grows / doesn't) per source.
- **Labels:** BacDive curated Biolog/API assays (sampling-independent lab measurements).
- **Cohort floor:** ≥100 E. coli strains per carbon source with downloadable assembly_accession AND a DE_CONFOUNDED within-MLST contrast.

## Architecture

- **Baseline (the bar to beat):** gene-content / KO presence-absence → RF or XGBoost (Li et al.'s method; reuse `dna_decode/eval/loso_kmer.py` shape + a KO/gene-presence feature builder).
- **Candidate:** NT v2 100M frozen mean-pool → XGBoost (same as AMR Phase 0; reuse `dna_decode/models/` + `eval/cv.py`).
- **v0 schema extension:** `phenotype` field generalizes the existing `drug` field in `pipeline.py predict`.

## Dataset prerequisite

1. A BacDive carbon-utilization export (long-format: strain × carbon-source × utilization). Acquisition = BacDive API/registration OR the Li et al. supplementary. **This is the only blocker** — all parsing/feasibility infra is built + green.
2. MLST per strain (for the de-confound gate) — from BacDive metadata or `mlst` on the assemblies.
3. assembly_accession per strain (NT-cacheable genomes) — BacDive↔NCBI crosswalk.

## Falsifier

- **Substrate falsifier:** the feasibility census (`scripts/bacdive_carbon_util_feasibility.py`) returns **0 FEASIBLE carbon sources** (no source has ≥100 de-confoundable E. coli strains with accessions) → EP-6 is **substrate-infeasible**; document + drop (same shape as the 2026-05-18 BV-BRC strict-MIC census that killed tet/gent).
- **Architecture falsifier:** on a FEASIBLE source, NT does NOT beat gene-content by ≥3 pp OR within-lineage concordance = chance → the embedding thesis is 0-for-3; pivot the decoder permanently to gene-content/deterministic heads.

## Go/No-Go gates (ordered — each STOPS)

1. **Acquire** a BacDive E. coli carbon export → `data/raw/`.
2. **Feasibility census** (`bacdive_carbon_util_feasibility.py --export X --mlst Y --min-strains 100`). Exit 0 = ≥1 FEASIBLE source. **No-Go if 0.**
3. **Build cohort** on the top-ranked FEASIBLE source (per-class MLST-balanced, assembly-accession-unique — reuse `build_stage2_n150_cohort.py` shape).
4. **Populate NT cache** (Databricks burst for N≥100; GTX 860M won't scale — see CLAUDE.md).
5. **Falsifier** (`amr_falsifier.py`-shaped runner adapted for carbon-util: NT-XGB vs gene-content-XGB, paired bootstrap CI) + **within-lineage diagnostic**. This is the terminal.

## What is built + green (2026-06-06)

- `dna_decode/data/bacdive.py` — carbon-util loader (tolerant column map, multi-vocab binarize, ambiguous-drop, organism filter, majority-vote, per-source census). 8 tests `tests/test_bacdive.py`.
- `scripts/bacdive_carbon_util_feasibility.py` — layered census (count → balance → accession → de-confound gate) reusing `cohort_deconfound.py`. Smoke-validated: confounded source → BLOCKED_CONFOUNDED; de-confounded → FEASIBLE.
- Reused as-is: `cohort_deconfound.py` (de-confound precondition), `within_lineage_diagnostic.py` (promotion gate), `amr_falsifier.py` (CI-aware runner template).

## Honest scope

Single substrate, single organism. A FEASIBLE verdict only proves a de-confoundable cohort *exists* — necessary, not sufficient. The architecture falsifier is the real test. This EP exists to give the embedding thesis its fairest possible shot (sampling-independent label, no curated catalog) and to take the answer — pass or fail — as the decoder's architectural input.
