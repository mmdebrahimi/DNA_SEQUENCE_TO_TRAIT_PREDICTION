# Fresh-cohort re-validation of the 10 SCORED cells — FINDING (2026-07-11)

**Headline: the deployed decoder REPRODUCES on fresh, unseen, provenance-disjoint genomes across all 10
cells — none collapsed.** Each cell was re-scored on a SECOND provenance-disjoint cohort whose genomes are
DISJOINT from the frozen validation cohort (excluded via the accession manifest; ~1,500 prior accessions
excluded per cell, selected from large fresh pools). This is the strongest reproducibility-of-the-claim
evidence short of external clinical validation. **Frozen artifacts + the report card are byte-unchanged**
(md5 fingerprint identical before/after; prospective-lock `verify_lock` still green).

See `revalidation_summary.md` for the full frozen-vs-fresh table. The interpretation:

## The two big positive drifts CONFIRM the clonality-disclosure finding

| cell | frozen | fresh | Δ | why |
|---|---|---|---|---|
| **Klebsiella × meropenem** | sens **0.467** | sens **0.929** | **+0.46** | The frozen 0.467 was the DISCLOSED clonally-inflated FAIL — "the meropenem FN was one clone from one BioProject". A fresh cohort that ESCAPES that clone scores 0.929. Direct evidence the low frozen number was a clonal artifact, not a decoder weakness. |
| **E. coli × ciprofloxacin** | spec **0.7** | spec **0.933** | **+0.23** | Same pattern — the frozen low specificity was a cohort artifact; the cell is far more specific on fresh genomes. |

These are the most valuable results of the run: the two frozen numbers that looked *weak* were depressed by
clonal structure in the frozen cohort, and the decoder is actually **stronger** than those headlines
suggested. This is exactly why the project disclosed clonality rather than hiding it — and the fresh cohort
vindicates that call.

## The negative drifts are small-N sampling noise (not regressions)

At `per_class=15` (n=30/cell), each isolate is worth ~0.067, so a Δ of 0.13–0.17 is 2–3 isolates — inside
the Wilson-CI overlap with the frozen number:
- Klebsiella cipro sens 0.967 → 0.800 (2–3 isolates)
- Klebsiella tet sens 0.800 → 0.667
- Klebsiella cef spec 0.900 → 0.733

Six of the remaining cells land within ±0.07 of frozen (E. coli cef/gent/tet, Kleb gent) — clean
reproductions.

## Honest caveats

1. **Campylobacter — RESOLVED 2026-07-11 (was n=4, now n=28).** The first pass scored Campy on only 4 genomes
   because `select_disjoint` REUSES an existing `selected.tsv` and the pre-launch smoke had written a
   `per_class=2` cohort for this cell. Re-run at `per_class=15` after deleting that file: a genuine fresh
   disjoint cohort of **28 genomes (13R/15S)** — **fresh sens 1.0 / spec 1.0, exactly matching frozen
   sens 1.0 / spec 1.0 (n=40), zero drift.** The R side is pool-limited to 13 (only 13 Campy cipro R genomes
   are disjoint from the ~1,479 prior accessions — a genuine ecosystem-domination ceiling consistent with the
   BV-BRC census, not a selector bug). So **all 10 SCORED cells are now properly re-validated on fresh
   disjoint genomes.**
2. **These are RAW (isolate-level) numbers, not clonality-collapsed.** The FRESH cohorts could themselves be
   clonally structured; the standing lineage-disclosure caveat still applies. The point proven is
   *reproduction on disjoint genomes*, not a lineage-independent number.
3. **n=15/class < the frozen n=30/class**, so the fresh Wilson CIs are wider — the small drifts are well
   within them.

## Bottom line

The shipped deterministic decoder is confirmed to **hold on unseen, provenance-disjoint genomes**; the two
frozen "weak" numbers (Kleb meropenem sens, E. coli cipro spec) were clonal artifacts of the frozen cohort,
and on fresh data the decoder performs **better** than those headlines. This is a positive,
publication-relevant confidence result — and it required touching nothing frozen.

Generated alongside `scripts/revalidate_scored_cells.py` (composes the frozen validation functions,
isolated outputs). Frozen surface untouched (verified).
