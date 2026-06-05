# HANDOFF — AMR Phase 2 embedding-validation frontier (workhorse-gated) — 2026-06-04

> Written from the GTX 860M laptop. REVISED 2026-06-04 to the **GPU-only** workhorse model (user
> directive: the workhorse is for its GPU, nothing else — no personal repo, no git, no dev there).
> The ONLY workhorse job is producing the NT embedding cache from PUBLIC data; **the falsifier + all
> repo/ledger/analysis work runs on the laptop.**
>
> ## ⚠ Division of machines (GPU-only)
> ```
> PUBLIC DATA (NCBI genomes + HuggingFace NT model)
>     │  workhorse fetches these itself — public, NOT the personal repo
>     ▼
> WORKHORSE (GPU only): run populate_cache.py -> ONE output file: the NT embedding cache (.h5)
>     │  ← the only artifact that crosses back (numeric embeddings of PUBLIC genomes; low DLP profile)
>     ▼
> LAPTOP (project home): git/ledgers + the falsifier (NT-XGBoost on cached embeddings vs k-mer vs
>                        AMRFinder — all CPU) + classical baselines (BLAST+/Docker already here)
> ```
> This keeps the personal repo OFF the corporate machine (Zscaler/Sentinel DLP), and prevents the
> bundle-drift that already hit `cv.py` (no dev inside the workhorse bundle).
>
> **For cef specifically the GPU work is likely already DONE** — the cef cached-strain result (AUROC
> 0.895) means the cef NT embeddings already exist on the workhorse (`nt_gate_b_cohort_67.h5`). So the
> only ask to the workhorse is: **return that one .h5 cache file.** Then the entire cef falsifier runs
> on the laptop (script staged at `scripts/cef_falsifier.py`).

## Why this is the next epoch (the strategic pivot, evidence-backed)

This session closed the pathotype embedding path: pathotype labels are sampling-defined categories →
study==class confound is intrinsic → no de-confounded discrimination cohort is buildable
(`research_outputs/horesh-f1-label-provenance-audit-2026-06-04.md`). The lesson generalizes:
**embeddings can only be honestly validated where the label is a sampling-INDEPENDENT lab measurement.**
AMR MIC is exactly that.

**De-risk evidence (computed on the local N=147 cipro cohort, 2026-06-04):**
- R and S **co-occur within 6 shared MLST lineages (35 R + 10 S strains)** — because cipro resistance is
  acquired by point mutations WITHIN a lineage (susceptible ST131 → resistant ST131). The MIC label is a
  broth-microdilution assay result, orthogonal to lineage/sampling. Contrast pathotype, where ExPEC vs
  EPEC essentially never share a both-class lineage (the label IS the sampling site).
- The cipro embedding signal already **survives leakage-safe CV**: `leave_one_accession_out` AUROC
  **0.8697** (vs pathotype k-mer 0.514 = chance). The substrate is validatable.

⇒ AMR is where the v0.1 embedding-vs-resolver/classical falsifier has an honest answer.

## Current AMR state (per CLAUDE.md + repo, verify on workhorse)
- cipro: v0 cached + v0.1 genome-input LANDED (Codex, Precision 7780). leakage-safe retrain 0.8697.
- cef: cached-strain + genome-input LANDED (AUROC 0.895, N=49, dup-accession PASS). Currently DEBUG-MODE
  (`--allow-missing-audit --no-attribution`). **Audit-aware packet IN FLIGHT** —
  `plans/Cef_Audit_Aware_Packet_Design.md` (5 artifacts + `scripts/drug_mechanism_phenotype_merge.py`).
- tet/gent: deferred (distributed-mechanism / substrate-infeasible per the 2026-05-18 census).
- AMR ledger: `project_state/dna-decode-2026-05-11.md`.

## The epoch — terminal claim + falsifier (roadmap Phase 2)

**Terminal claim:** `pipeline.py predict --drug X` for X in {cef} ships an audit-aware packet AND the
embedding model beats the best classical baseline (k-mer / AMRFinder-gene-presence) by **≥3 pp AUROC**
under leakage-safe CV — OR a documented mechanism-class scope-limit (EP-1.5 finding) names why not.

**Falsifier:** embeddings do NOT beat classical on cef (concentrated β-lactamase mechanism) → the
NT-mean-pool advantage is absent even on a concentrated mechanism at clean labels → re-evaluate the
embedding architecture before any further drug.

## ⚠ UPDATE 2026-06-04 — cef cohort is CONFOUNDED; clean substrate is cipro N=147

A `/brainstorm` proved `gate_b_cohort` (cef) is geography/lineage confounded (R≈USA, S≈Africa/India,
1 shared MLST) → an NT verdict there measures batch, not biology. The cef run is BLOCKED
(`scripts/cef_falsifier.py` now exits 4 via the new de-confound gate). The prior cef AUROC 0.895 is
non-citable. See `plans/cef_falsifier_brainstorm.md`.

**New tooling:** `dna_decode/eval/cohort_deconfound.py` (de-confound gate, a PRECONDITION) +
`scripts/amr_falsifier.py` (drug-agnostic; gate baked in; CI-aware verdict).

**The decisive clean test = cipro N=147** (6 shared R/S lineages → DE_CONFOUNDED). One-file GPU-only ask
to the workhorse: return `nt_n147_cipro.h5`. Then on the laptop:
```
uv run python scripts/amr_falsifier.py --drug ciprofloxacin \
  --cohort data/processed/stage2_n150_cipro_cohort.parquet \
  --nt-cache data/processed/embeddings/nt_n147_cipro.h5
```
(cipro genomes for the k-mer baseline are public-NCBI fetchable on the laptop, same as cef.)

## Sequenced next steps (GPU-only split)

### Workhorse — ONE job (GPU only)
- **Return the cef NT embedding cache** `nt_gate_b_cohort_67.h5` to the laptop (it already exists from the
  cef cached-strain run). If for some reason it must be regenerated: `scripts/populate_cache.py` on the
  cef cohort (`data/processed/gate_b_cohort.parquet`, 60 cef-labelled strains) — pulling genomes from
  public NCBI + the public NT model. No git, no repo, no ledger work on the workhorse.
- (Workhorse housekeeping, your call: commit its stranded `reports/dna_decoder_dual_machine_handoff_2026-06-04.md` to main, or just relay it — but per the GPU-only model, ideally nothing on that machine touches git.)

### Laptop — everything else (CPU; no GPU)
1. Drop the returned cef cache at `data/processed/embeddings/nt_gate_b_cohort_67.h5` (or pass `--nt-cache`).
   Fetch the cef genome FASTAs from public NCBI for the k-mer baseline (the script prints the command if missing).
2. **Run the falsifier:** `uv run python scripts/cef_falsifier.py` — NT-XGBoost + NT-logreg + k-mer-XGB
   under `leave_one_accession_out` CV (cef has no duplicate accessions, so this equals strain-out; canonical
   + future-proof). Optional second pass under Mash-clade-out via `scripts/mash_cluster_n147.py` (Docker).
   Gate: max(NT) − k-mer-XGB ≥ 3 pp → embeddings beat classical on a concentrated β-lactamase mechanism.
3. **Cef audit-aware packet** (`plans/Cef_Audit_Aware_Packet_Design.md`; `scripts/drug_mechanism_phenotype_merge.py`):
   CPU-only — runs here. NOTE: the workhorse reports the packet design is NO LONGER the blocker; audit it
   for stale DEBUG-MODE drift rather than rebuilding.
4. **Update both ledgers** (AC9) here: AMR `dna-decode-2026-05-11` + roadmap Phase 2 row.

## What stays on the laptop (no GPU)
- The falsifier (`scripts/cef_falsifier.py`), the k-mer + AMRFinder baselines, the audit/merge scripts,
  all ledger ops, all git, all doc reconcile. The workhorse does NOT touch any of it.

## Provenance
De-risk numbers: `data/processed/stage2_n150_cipro_cohort.parquet` (R/S × MLST overlap, inline pandas).
Strategic basis: `research_outputs/horesh-f1-label-provenance-audit-2026-06-04.md` +
[[feedback_sampling_defined_phenotype_intrinsic_confound]]. Roadmap: `plans/Trait_Decoding_Roadmap.md`
Phase 2.
