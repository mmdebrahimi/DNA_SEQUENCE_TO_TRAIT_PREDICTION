# Klebsiella K-antigen (wzi) — validation report card

Namespace-separate typing-trait card (sibling of `serotype`). **Honest tier:** `FAITHFUL_TO_TOOL_MEASURED_LABEL_AVAILABLE`.

deterministic single-gene wzi caller, FAITHFUL to the Kleborate/BIGSdb wzi method (NOT an independent baseline). A free MEASURED serological K-type label EXISTS (KlebNET-GSP 731-isolate set), so the cell is VALIDATABLE -- but the full wzi-caller-vs-serology number needs the cohort genome run (scoped). wzi->K is ~94% / NOT one-to-one; full-locus Kaptive is more accurate.

## Caller self-consistency
- **15/15** sampled wzi alleles (of 555 in the DB) typed correctly -> correct KL (PASS). The caller is mechanically sound across the DB.

## Measured-label ceiling (the rare non-AMR free measured label)
- Genomic capsule typing (Kaptive KL) vs **SEROLOGICAL K-type** (wet-lab) on the KlebNET-GSP **N=733** set: naive-numeric **0.745** | paper-curated **0.845**.
- naive KL#==K# match UNDER-counts (drops the KL<->K renaming/alternative-type equivalence the paper applies); the curated rate is the real ceiling.
- This is the *ceiling* a genomic K-typer approaches; the single-gene wzi-v0 is ~94% of full-locus Kaptive. The genuine wzi-caller-vs-serology number is the scoped scale-up below.

## Pending scale-up
- run the wzi caller on the 731 genomes (ENA run reads -> targeted wzi mapping) -> the genuine wzi-caller-vs-measured-serology number (namespace-separate).

## Honesty rails
- FAITHFUL-TO-TOOL (wzi method), NOT an independent baseline; `caller_is_independent_baseline=false`.
- wzi->K-type is ~94% predictive and NOT one-to-one (Brisse 2013 JCM).
- SEPARATE from the AMR/HIV/TB trust cards -- ktype is a typing trait, not a drug cell; never keyed into the frozen canonical_cell_key card.
- Sources: KlebNET-GSP 2025 Technical Report (Zenodo 15742130, CC-BY); Kleborate wzi DB (BIGSdb Pasteur); Brisse 2013 JCM (wzi typing ~94%).

Rebuild: `uv run --with openpyxl python scripts/ktype_validate.py --xlsx <suppl.xlsx> --db-dir data/ktype_db`.
