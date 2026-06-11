# Klebsiella meropenem provdisjoint FN dissection — the "blind spot" is a clonal single-submitter batch — 2026-06-11

The carbapenemase rule scored **sens 0.467** (FN=16/30 R) on the provenance-disjoint Klebsiella-meropenem
cohort vs **sens 1.0 in-cohort (N=30)**. The "mechanism-complete decoder" idea anchored on this gap as a
diverse missing-mechanism class (porin-loss + efflux + AmpC hyperproduction) worth building scanners for.
The idea-anchor + probe demanded the FN be dissected FIRST. It was. The gap is a phantom.

## Three grounded findings (all from on-disk data — no re-fetch)

1. **16/16 FN carry ZERO acquired determinants.** Re-applying `call_resistance(Klebsiella, meropenem)` to each
   FN's on-disk AMRFinder `main.tsv`: none carry a carbapenemase, ESBL, AmpC, or porin marker AMRFinder reports.
   So the meropenem-R label is unexplained by any acquired/point determinant.

2. **15/16 FN are ONE BioProject.** Parsing each assembly's `assembly_data_report.jsonl`:
   - **15 of 16 = PRJNA504784, University of Rochester, Rochester NY, mostly 2017.**
   - 1 of 16 = PRJNA234117 (Broad Institute, Boston, 2013).
   Half of the cohort's entire R class (15 of 30) came from a single clinical submitter's collection.

3. **13/16 FN are ONE near-clonal lineage.** Mash on the 16 FN genomes: the Rochester strains have median
   pairwise distance **0.0001** (>99.99% ANI), 78/105 pairs near-clonal (<0.001). Single-linkage at Mash<0.005
   yields **2 clusters** — one of 13 near-clonal strains + one of 3. The "16 false negatives" are effectively
   **~2 independent biological observations**, one outbreak clone sequenced ~13×.

## Verdict — the premise is invalid; do NOT build the scanners

The meropenem **sens 0.467 is clonal / single-submitter batch-inflated.** The rule did not miss 16 diverse
carbapenem-R strains — it missed ~2 effective lineages, one of which one hospital sequenced 13 times. Building
a porin-loss / efflux / AmpC-hyperproduction scanner against this FN set would **overfit to one clone from one
institution** (the study==class / validate-on-independent-data traps this project has hit repeatedly). "Recovering"
the 13-clone would lift the metric by recovering a single lineage sequenced 13× — not generalizable mechanism
detection. The honest recoverable-sensitivity ceiling on a *clonally-deduplicated* cohort is unknown but the FN
denominator collapses from 16 to ~2.

Not determined (and not decision-relevant): whether the clone is genuine non-carbapenemase CRE (porin-LoF) or a
mislabeled/borderline-AST batch. An ompK35/36 BLAST would refine the *why* but cannot change the verdict —
one clone either way.

## The REAL finding — a trust-surface hardening gap (high value, generalizable)

The provenance-disjoint cohort selector (`scripts/provenance_disjoint_validate.py::select_disjoint`) excludes the
surveillance ECOSYSTEM but does **NOT deduplicate by clonality** — so one outbreak clone can dominate a class and
inflate/deflate sens/spec on ANY cell, silently. This is the genuinely valuable next move that the dissection
surfaced:

- **Add clonality-aware dedup** (Mash-cluster the selected cohort, keep ≤1 per near-clonal cluster per class)
  BEFORE scoring, so the report-card metrics measure independent lineages, not sequencing-depth of one outbreak.
- **Re-score meropenem** after dedup → the honest sens estimate (the 13-clone collapses to 1 effective R).
- This improves EVERY cell's metric honesty — it is a report-card capability, not a decoder mechanism gap.

## Recommendation

KILL the "mechanism-complete decoder" framing as scoped. Replace the next move with **clonality-aware cohort
dedup in the provenance-disjoint validator** (+ re-score the affected cells). The "decoder misses non-carbapenemase
CRE" claim is unproven on independent data and cannot be earned from this single-clone batch.
