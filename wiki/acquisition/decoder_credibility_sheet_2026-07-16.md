# dna-decode — what it is, what it's validated on, and what we'd do with your data

*One page. Every number below is machine-generated from a committed artifact in the repo — no aggregate
"X% accurate" claim, because the honest answer is per-cell. Research use only; not a clinical tool.*

> **Audience: technical** (bioinformatician / clinical microbiologist / PI). For a non-specialist reader — a
> lab manager, an admin who gates data access, or anyone deciding whether this is worth a meeting — send
> **`decoder_credibility_sheet_PLAIN_2026-07-16.md`** instead (same facts + numbers, no jargon). When in
> doubt, send **both**: the plain sheet first, this one as the appendix.

## What it is

An **open-source, deterministic genotype→phenotype decoder**. Point it at a genome (or a VCF) and it calls
curated, mechanism-level traits — antibiotic resistance, pathotype, plasmid/serotype/MLST typing,
pharmacogenomics, and (new) visible traits. It is **rule-based over curated determinant catalogues**, not a
black-box model: every call names the determinant it fired on, and **abstains** when the mechanism isn't
decodable rather than guessing.

`pip install dna-decode` (v0.7.0, PyPI) · CLI: `dna-decode amr | pathotype | plasmid | serotype | mlst | pgx | …`

## What it's validated on (the honest grid)

**27 decoder cells** on the standing trust surface: **10 SCORED · 3 UNDERPOWERED · 2 ABSTAINS_BY_DESIGN ·
1 LABEL_CONFOUNDED · 11 NO_FREE_PHENOTYPE_SOURCE**. Highlights, each with its own caveat attached:

| cell | validation | number |
|---|---|---|
| **Bacterial AMR** (E. coli / Klebsiella / Campylobacter × cipro/cef/gent/tet/meropenem) | **provenance-disjoint** NCBI-PD cohorts (submitters *outside* NARMS/CDC/FDA/GenomeTrakr/PulseNet) | 10 SCORED cells; raw sens/spec disclosed **alongside clonality-corrected**, lineage-effective N + Wilson CI |
| **M. tuberculosis** RIF / INH | **independent** EBI AMR-Portal cohort, provenance-disjoint from the WHO-catalogue build set (N=2,845) | RIF raw 0.920/0.955 · **lineage-collapsed 0.444 [0.246–0.663] / 0.979** |
| **HIV-1** (5 drug classes, 4 genes) | **free, independent, isolate-level wet-lab labels** (Stanford PhenoSense fold-change — *not* HIVDB's own interpretation, so non-circular) | efavirenz **AUC 0.962** · nevirapine **0.985** |
| SARS-CoV-2 Mpro, *C. auris* azoles | in-distribution / label-limited — **stated as such** | — |

## Why the honesty posture is the point

- **We publish the number that makes us look worse.** TB raw sens is 0.920; the clonality-corrected number is
  0.444, and *that* is our headline, because 2,845 isolates collapse to ~67 lineages and raw counts one vote
  per isolate.
- **We record negative results — and the scope of each is exact.** **ZERO-SHOT / frozen** genomic-FM
  embeddings failed under de-confounding on every substrate tried (0-for-5, across the bacterial→eukaryote
  boundary); zero-shot ESM2 scored *below chance* where the curated catalogue hits 0.926. **That negative is
  about zero-shot embeddings — it is NOT a verdict on learned models.**
- **The shipped architecture is a HYBRID, not a rules purist.** A **supervised** sequence model *rescues the
  catalogue's structural blind spot* (catalogue-negative resistant isolates): **0.81 leave-STUDY-out**
  (deployable) vs **0.449** for zero-shot ESM on the same gap; **GENERAL_RESCUE** across 8/11 HIV RT drugs
  (all 5 NNRTIs). Catalogue fold-in was tested and **rejected** (hard rules trade sens for spec, −0.006
  bal-acc) — the value is the continuous weighting, so it ships as an **offline complement**:
  `dna_decode/data/hiv_supervised_complement.py` + 3 committed model JSONs. Deterministic rule leads +
  abstains; learned layer covers its blind spot.
- **The regime boundary is the real finding** (`wiki/supervised_learned_layer_synthesis_2026-07-12.md`):
  supervised rescue works where resistance is **CONVERGENT** (HIV — same mutations recur across the
  phylogeny) and fails where it is **CLONALLY CONFOUNDED** (TB RIF — plain 5-fold said 0.66, but
  leave-one-lineage-out revealed **0.51 = chance**; it had learned lineage, so we did not ship it).
  **The de-confounded split is the whole ballgame** — and the HIV win is trustworthy *because* it survived
  the equivalent de-confound that killed the TB one.
- **Cells abstain by design** when a mechanism isn't gene-decodable (efflux/regulatory/intrinsic), rather than
  over-calling.
- **Frozen + reproducible:** the deployed rule surface is sha256-pinned and byte-verified on every run; a
  prospective-lock manifest makes any future score provably leakage-free by date.

## What we're missing — and what we'd do with your data

**The binding constraint is labels, not models.** Public genotype↔phenotype pairs are saturated for the
organisms above; what we cannot get freely is **measured, isolate-level phenotype paired with sequence** for
new organisms/drugs.

**Given a paired WGS + measured-AST (MIC or disk) set, we would:**
1. Score the **existing frozen decoder** on it as a genuinely external cohort — no retraining, no tuning, so
   the result is a real out-of-distribution test (we've built this arm already: leakage preflight at the
   BioSample level, operator-aware MIC censoring, powering gate).
2. Report **per-cell sens/spec with clonality correction + CI**, including the failures.
3. Where a mechanism is missing, **extend the curated catalogue** and say exactly what it can't see.
4. Return the analysis to you — **your data stays yours**; we publish only what you approve, and can work
   under an MTA / data-sharing agreement. Nothing is redistributed without permission.

**What we do NOT need:** patient identifiers, clinical metadata, or anything beyond sequence + the AST result
and its method.

## Provenance

Every claim here traces to a committed artifact: `wiki/decoder_validation_report_card.{md,json}` ·
`wiki/tb_independent_lineage_collapsed_result_2026-07-02.md` · `wiki/hiv_decoder_report_card.md` ·
`wiki/reproducibility_freeze_2026-06-13.md` · `wiki/negative_results_map_2026-06-13.md`.
Repo: `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`
