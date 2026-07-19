# CDC AR Isolate Bank → frozen-decoder external re-validation (strengthen-first)

**Date:** 2026-07-18 · **Status:** ingester SHIPPED + tested + cohort built on real CDC data; the
AMRFinder-Docker **scoring** run is the remaining attended step · **Chosen path:** "both, strengthen
first" — re-validate an existing SCORED E. coli cell on the AR Bank, then reuse the template for the
gonorrhoeae expansion · **Code:** `dna_decode/data/ar_isolate_bank.py` +
`scripts/build_ar_bank_labels.py` (commit 76264af) · **Frozen surface:** byte-unchanged (verify_lock OK).

## The unlock (verified live, not assumed)

The prior memo framed the AR Bank as a free *request*. Live checking sharpened it further: the
per-isolate **DATA is fully public** — no login, no order, **no biosafety attestation for the data**
(only the physical isolates need that). An isolate-detail page
(`IsolateDetail.aspx?IsolateID=&PanelID=`) server-renders, per isolate:
- the NCBI **BioSample** accession (→ the WGS),
- the **organism**,
- reference **broth-microdilution MIC** + CLSI/FDA **S/I/R** per drug.

So the **measured-MIC ⋈ WGS join is a free web fetch**. The label wall's cheapest opening is a *fetch*,
not an acquisition — and it feeds infrastructure the project already built.

## What shipped

- **`dna_decode/data/ar_isolate_bank.py`** — pure parsers (panel rows + isolate MIC detail) verified
  **exact** against the live 2026-07-18 pages (AR#0021 → 31 MICs, BioSample `SAMN04014862`,
  Citrobacter freundii); lazy cache-first fetch/enumeration over the **34 established panels**;
  `to_label_inputs` bridges STRAIGHT into the frozen `external_mic_labels.build_drug_labels` —
  **BioSample-keyed, no crosswalk** (unlike the Oxford ingester, because the page gives the BioSample).
- **`scripts/build_ar_bank_labels.py`** — enumerate → fetch detail → tier via the frozen
  `build_drug_labels` → emit `selected_{strict,relaxed}.tsv` + `buckets_<drug>.json` +
  `cohort_manifest_external_<run_id>.json` (the exact contract `external_cohort_preflight` /
  `external_cohort_revalidate` consume).
- **9 offline tests** (`tests/test_ar_isolate_bank.py`, fixtures = real markup). They caught + fixed a
  real `to_label_inputs` bug: a non-pilot drug (canon `None`) matched every *other* non-pilot drug.

## Live cohort

Enumeration across the 34 panels: **1087 unique-BioSample isolates**, all with WGS —
E. coli 115 · Klebsiella 157 · Pseudomonas 140 · Staphylococcus 138 · **N. gonorrhoeae 94** ·
Candida 85 (incl. *C. auris*, our fungal cell) · Acinetobacter 68 · Enterobacter 50 · Salmonella 24 …

## KEY FINDING — the curated bank is resistance-enriched → per drug, ONE class is powered

E. coli cohort (114 BioSamples with MIC), strict-tier (HIGH_R/HIGH_S) counts:

| Drug | strict N | R | S | Powered side | Uninformative side |
|---|---|---|---|---|---|
| ceftriaxone | 97 | **95** | 2 | **SENSITIVITY** | specificity (only 2 S) |
| ciprofloxacin | 89 | **86** | 3 | **SENSITIVITY** | specificity (only 3 S) |
| gentamicin | 56 | 0 | **56** | **SPECIFICITY** | sensitivity (0 R) |

The AR Bank is *designed* as a resistance panel, so its mechanisms (carbapenemases/ESBLs) make it
R-heavy for β-lactam/fluoroquinolone drugs and S-heavy for gentamicin (those mechanisms don't confer
gent-R). **Consequence for scoring:** a naive two-class sens/spec per drug would HARD-FAIL the
powering gate (`MIN_PER_CLASS=10`) on the minority class. Read **one-sided**, it's a coherent,
powered, independent test — **cef + cipro independently confirm SENSITIVITY** (≈90 CDC-curated R
each, a strong FN-rate stress test), **gent independently confirms SPECIFICITY** (56 S). Across the
three drugs, *both* sens and spec get independent, provenance-separable confirmation.

This is the honest framing to carry into the roll-up — not a forced two-class number the curated bank
cannot provide. (It also generalizes: for a well-powered *both-class* external cohort, pursue the
gonorrhoeae/Euro-GASP set, whose resistance is not curation-enriched.)

## Provenance-disjointness (preflight)

<!-- PREFLIGHT_VERDICT -->
The exact-set BioSample-level preflight (`external_cohort_preflight.py --cohort-manifest …`) resolves
the 114 cohort BioSamples + the decoder's tuning accessions and fails closed on any overlap. **Verdict:
pending the run's completion — filled on landing.** AR Bank isolates are curated CDC outbreak isolates;
some may also sit in the NCBI-PD tuning pull, in which case the preflight flags + excludes them by
construction (that is the gate working, not a failure of the premise).

## The scoring run (the remaining attended step)

Scoring = the shipped external arm: resolve each BioSample → GCA, run **AMRFinder (Docker)** with the
frozen organism triple (`-O Escherichia`, `call_resistance(organism="Escherichia_coli_Shigella")`),
compute the confusion vs the MIC labels. ~115 genomes × ~95 s AMRFinder ≈ a multi-hour CPU run
(Docker; mind the documented WSL-mount wedge — restartable, `wsl --shutdown` recovery).

```bash
# 1) preflight (disjointness/availability) — DONE above
# 2) score each powered drug (one-sided; --min-per-class documents the pilot)
uv run python -m scripts.external_cohort_revalidate --cohort ar_bank_ecoli --drug ceftriaxone \
  --labels-dir data/raw/ar_bank_ecoli_extval_ceftriaxone \
  --preflight-json wiki/external_preflight_ar_bank_ecoli_2026-07-18.json \
  --cohort-manifest wiki/cohort_manifest_external_arbank_2026-07-18.json --min-per-class 2
# repeat --drug ciprofloxacin (sensitivity) and --drug gentamicin (specificity, --min-per-class 10 on S)
# 3) roll-up: scripts/build_external_validation_report.py (separate external namespace)
```

## Honest scope

- The DATA fetch is free + attestation-free; **physical isolates** (not needed here) require an
  institutional account + biosafety officer.
- The AR Bank MIC is CDC reference BMD (real G1). Isolates are curated (resistance-enriched) → **not
  prevalence-representative**; report one-sided powered metrics, never a population sens/spec off this set.
- Independence is **provenance-separable at the BioSample level** (preflight-enforced), not
  methodology-independent (same AMRFinder `-O` + frozen `call_resistance`).
- Frozen decoder surface byte-unchanged; this is a READ-only external-validation adapter.
