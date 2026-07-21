# MaveDB scoping for the `forward` variant-effect cell (2026-07-21)

**Question:** what would a MaveDB pull add to the `forward` (R2 molecular / protein variant-effect) cell,
currently validated on ProteinGym v1.1? **Verdict:** MaveDB is a genuine, free, CC0, human-inclusive
expansion — but NOT as a ProteinGym replacement (heavy overlap + loses curation). Its two high-value uses are
**targeted human-gene pulls** and a **leakage-free prospective validation tier**. This is a scoping memo, not
a build.

## Grounded facts (verified in-session via api.mavedb.org, read-only)
- **License = CC0 (public domain)** on 98/100 sampled score sets (2 are CC BY-NC-SA). Fully usable, no DUA,
  no access gate. (This is stronger than the ChatGPT report's vague "open-access".)
- **Human-heavy:** 69/100 sampled score sets target *Homo sapiens* (rest: artificial sequences + assorted
  bacteria). MaveDB directly serves the "move to human" north star in the ONE regime where learned models
  work.
- **Recent:** sampled published dates span 2021-04 → 2026-06; **40% are ≥2024** — newer than ProteinGym
  v1.1's collection cutoff.
- **All protein_coding** in the sample; 224,055 variants across 100 score sets. The full catalog is larger
  (the search endpoint page-caps at 100; MaveDB's own site reports thousands of score sets — not verified
  in-session).
- Record model: `POST https://api.mavedb.org/api/v1/score-sets/search` → `{scoreSets:[{urn, numVariants,
  license, targetGenes:[{name, category, targetSequence.taxonomy.organismName}], publishedDate, experiment,
  ...}]}`. Per-score-set variant CSV downloadable by URN.

## Current baseline (what the forward cell already has)
- **ProteinGym v1.1 = 217 substitution DMS assays** (Zenodo 13936340), forward cell validated across all 217
  (`wiki/proteingym_esm2_650m_full_2026-07-09.md`; hybrid ESM2+ProSST(+GEMME) is the deployed scorer).
- Committed catalog at `wiki/proteingym_v1.1_substitutions_catalog.tsv` — the dedup key against MaveDB.

## Relationship (the load-bearing fact)
**ProteinGym is CURATED FROM MaveDB + other sources.** ProteinGym adds benchmark discipline MaveDB lacks:
standardized reference sequences, normalized DMS score direction (higher = more fit), quality filtering, and
a train/test split. MaveDB is the raw upstream superset. So:
- Pulling MaveDB wholesale mostly re-imports what ProteinGym already curated (**overlap**), minus the
  curation — you'd redo reference alignment + score-direction normalization + dedup.
- The marginal value is the assays ProteinGym did NOT include: newer submissions, non-substitution MAVEs
  (indels/regulatory — NOT yet consumable by the substitution-only forward cell), and specific target genes.

## Gate screen (R2 molecular regime) — clean
- **G1 circular-label:** CLEAR. The label is a wet-lab functional readout, not a predictor's output.
- **G3 sampling-defined:** CLEAR. A designed saturation-mutagenesis library, not field isolates.
- **Population-structure / clonality confound:** **NONE BY CONSTRUCTION** — a within-protein mutant library
  has no lineage/ancestry axis. This is the cleanest regime the project has; the confound that killed
  Arabidopsis + cipro-within-lineage cannot arise here.
- **G8 / leakage (the ONE real trap):** MaveDB ⊇ ProteinGym, so any MaveDB pull MUST be deduped against
  `proteingym_v1.1_substitutions_catalog.tsv` (by target gene + assay/DOI) before scoring, or the hybrid's
  ProteinGym numbers leak into a "new" MaveDB eval.

## Recommendation — two high-value moves, one free near-term
1. **(Highest VOI, free, novel) Prospective-lock the molecular cell.** MaveDB score sets published AFTER
   ProteinGym v1.1's cutoff, deduped against the committed catalog, are a **leakage-free-by-time** validation
   set for the FROZEN hybrid — the exact prospective-lock discipline the AMR track already uses
   (`prospective_lock.py`), now applied to R2. ~40% of the sample is ≥2024 → a real, growing, free holdout.
   No new training, no money, no confound. This ADDS a validation tier rather than de-risking a build.
2. **(On-demand) Targeted human-gene cells.** When a specific human clinical gene is wanted (BRCA1, PTEN,
   TP53, MSH2, …), pull its MaveDB DMS directly (CC0, already human) rather than waiting for a benchmark.
   This is how the forward cell "does human" concretely.
- **NOT recommended:** wholesale MaveDB re-ingest to "widen training" — mostly ProteinGym overlap, loses
  curation, and the substitution-only cell can't consume the non-substitution assays that are MaveDB's true
  additive content.

## Concrete first step (if the user picks lane R2)
`scripts/mavedb_prospective_holdout.py`: search score-sets → filter `publishedDate >= <ProteinGym cutoff>` +
`organismName == Homo sapiens` (or any) + `category == protein_coding` → dedup vs
`proteingym_v1.1_substitutions_catalog.tsv` → download variant CSVs → score with the frozen hybrid → report
Spearman per assay with a Wilson/bootstrap CI. Mirrors the AMR prospective-lock harness. Deferred pending the
lane decision (below).
