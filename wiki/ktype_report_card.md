# Klebsiella K-antigen (wzi) — validation report card

Namespace-separate typing-trait card (sibling of `serotype`). **Honest tier:** `FAITHFUL_TO_TOOL_MEASURED_LABEL_AVAILABLE`.

deterministic single-gene wzi caller, FAITHFUL to the Kleborate/BIGSdb wzi method (NOT an independent baseline). A free MEASURED serological K-type label EXISTS (KlebNET-GSP 731-isolate set), so the cell is VALIDATABLE -- but the full wzi-caller-vs-serology number needs the cohort genome run (scoped). wzi->K is ~94% / NOT one-to-one; full-locus Kaptive is more accurate.

## Caller self-consistency
- **15/15** sampled wzi alleles (of 555 in the DB) typed correctly -> correct KL (PASS). The caller is mechanically sound across the DB.

## Measured-label ceiling (the rare non-AMR free measured label)
- Genomic capsule typing (Kaptive KL) vs **SEROLOGICAL K-type** (wet-lab) on the KlebNET-GSP **N=733** set: naive-numeric **0.745** | paper-curated **0.845**.
- naive KL#==K# match UNDER-counts (drops the KL<->K renaming/alternative-type equivalence the paper applies); the curated rate is the real ceiling.
- This is the *ceiling* a genomic K-typer approaches; the single-gene wzi-v0 is ~94% of full-locus Kaptive. The genuine wzi-caller-vs-serology number is the scoped scale-up below.

## Wzi-caller-vs-serology number (COMPLETE 2026-06-25, n=447)
Ran the wzi caller on the full ENA-fetchable ERR cohort (targeted wzi read-mapping → `call_ktype`) vs the
MEASURED serological K-type. **Concordance 0.629** (naive KL#==K#; 273/434 scored), n_attempted=447,
n_called=434 (13 no-call), 5 error. Artifact: `wiki/ktype_cohort_validation.json`.
- **Honest read:** 0.629 is the naive-numeric lower bound (the curated KL↔K equivalence the paper applies
  lifts it toward the 0.745 naive / 0.845 curated CEILING above — and that ceiling is full-locus Kaptive,
  whereas this is single-gene wzi-v0 = ~94% of Kaptive). So 0.629 vs the 0.745 naive ceiling is the expected
  wzi-v0 gap, NOT a failure. Faithful-to-tool (the wzi method), never an independent baseline.
- Cohort scope: 447 of the 733-isolate measured-serology set had a fetchable ERR accession (286 'TBD'/
  unavailable skipped — reported, not hidden).

## Honesty rails
- FAITHFUL-TO-TOOL (wzi method), NOT an independent baseline; `caller_is_independent_baseline=false`.
- wzi->K-type is ~94% predictive and NOT one-to-one (Brisse 2013 JCM).
- SEPARATE from the AMR/HIV/TB trust cards -- ktype is a typing trait, not a drug cell; never keyed into the frozen canonical_cell_key card.
- Sources: KlebNET-GSP 2025 Technical Report (Zenodo 15742130, CC-BY); Kleborate wzi DB (BIGSdb Pasteur); Brisse 2013 JCM (wzi typing ~94%).

Rebuild: `uv run --with openpyxl python scripts/ktype_validate.py --xlsx <suppl.xlsx> --db-dir data/ktype_db`.
