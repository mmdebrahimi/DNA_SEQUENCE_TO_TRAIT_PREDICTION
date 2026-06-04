# HANDOFF — AMR Phase 2 embedding-validation frontier (workhorse-gated) — 2026-06-04

> Written from the GTX 860M laptop. The embedding core (NT cache populate + train) needs the **Precision
> 7780 GPU** — the 860M orphans on multi-minute GPU jobs. Per the two-machine rule, this is a genuine
> hand-off, not laptop work. This doc is the workhorse-executable next epoch.

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

## Sequenced next steps (workhorse)
1. **Finish the cef audit-aware packet** (`plans/Cef_Audit_Aware_Packet_Design.md`): run
   `scripts/drug_mechanism_phenotype_merge.py` for cef; flip cef out of DEBUG-MODE (audit verdict +
   attribution in the v0 packet). Laptop-doable PARTS: the merge + audit (no GPU) — but cef cohort/cache
   live on the workhorse, so co-locate.
2. **Embedding-vs-classical falsifier on cef** under `leave_one_accession_out` + Mash-clade-out CV:
   NT-XGBoost vs k-mer-XGB vs AMRFinder-gene-presence. Reuse `dna_decode/eval/cv.py` +
   `dna_decode/eval/loso_kmer.py`. Record the ≥3 pp gap (or scope-limit).
3. **Update both ledgers** (AC9): AMR ledger `dna-decode-2026-05-11` + the roadmap Phase 2 row.

## What stays on the laptop (no GPU)
- Cohort/label de-risk + audit analysis (like the cipro de-risk above), the merge/audit scripts, ledger
  ops, doc reconcile. NOT the NT populate/train.

## Provenance
De-risk numbers: `data/processed/stage2_n150_cipro_cohort.parquet` (R/S × MLST overlap, inline pandas).
Strategic basis: `research_outputs/horesh-f1-label-provenance-audit-2026-06-04.md` +
[[feedback_sampling_defined_phenotype_intrinsic_confound]]. Roadmap: `plans/Trait_Decoding_Roadmap.md`
Phase 2.
