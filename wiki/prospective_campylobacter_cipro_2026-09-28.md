# Prospective-lock accrual — Campylobacter × ciprofloxacin (2026-09-28)

**Artifact:** `wiki/prospective_lock_validation_Campylobacter_ciprofloxacin_2026-09-28.json`
**Lock:** `wiki/prospective_lock_manifest_2026-08-31.json` (the LIVE v2 lock; `verify_lock` OK, 5 files pinned, `frozen_commit=edb4c5e`)

## Result

| | |
|---|---|
| N scored | **426** (149 R / 277 S), 0 excluded, 0 abstain |
| acc | **0.998** |
| sens | **0.993** (148/149) |
| spec | **1.000** (277/277, 0 FP) |
| powering | `POWERED` |

Leakage-free BY CONSTRUCTION: every isolate became public `2026-09-03`, strictly after the
2026-08-31 lock, so the decoder cannot have been tuned to it. This is the first prospective cell
scored against the **live v2 surface** — both E. coli cells still read
`superseded_by_surface_change` (their numbers describe the retired v1 rule and stay WITHHELD).

## What this discharges — the load-bearing point

The rule that fired is **not** the frozen default. It is the opt-in registry entry, and its own
text states the precondition:

> `calibrated_organism_v1 (Campylobacter|ciprofloxacin; IN-SAMPLE N=30, LOO bal-acc 1.0).`
> `Registry is opt-in; needs an independent cohort before becoming a default.`

This cohort **is** that independent cohort: N=426 (14× the calibration set), a different
BioProject, a different country, a different surveillance network. That is the stated condition
met, measured rather than argued.

## The confusion reproduces independently

Recomputing from the cached AMRFinder runs by calling the deployed `call_resistance` directly:
**426 scored / 1 FN / 0 FP** — identical to the artifact. The number is not an artifact-writing
accident.

## The single FN is a DECLARED blind spot, not drift

`GCA_060803015.1` / `SAMN62387178` carries **zero determinants** — no `gyrA_T86I`, nothing. The
rule's own output names exactly this case:

- `undetectable_mechanisms: ['efflux', 'porin_loss', 'regulatory']`
- `caveat: "An S call cannot rule out resistance via efflux, porin_loss, regulatory."`

So the one miss falls inside a limit the decoder states at call time. Per the standing rule
(diagnose the FN before reading it as decay) this is **abstention-by-design leaking into a hard S
call**, not catalog drift and not a new gap.

## Limits — read these with the 0.998

1. **SINGLE-SOURCE. All 426 isolates are one BioProject** — `PRJNA563067`, *PulseNet Canada
   Campylobacter whole genome sequencing*, Public Health Agency of Canada (umbrella
   `PRJNA563656`). Measured, not inferred: 40 GCAs sampled stratified across the entire accession
   range resolve **40/40** to that one project. Largest-source share **1.000**, which **FAILS this
   repo's own diversity bar** (`source_diverse_validate`: ≥5 BioProjects AND ≤0.60 share) by the
   maximum possible margin. A single-source cohort cannot reveal a determinant-family blind spot
   it happens not to contain — the mechanism that made the E. coli gentamicin `rmt` gap
   structurally invisible.
   Corroborating layout signal: the 426 BioSamples occupy 426 of 428 *consecutive* SAMN slots
   (span 427, density 0.995, one cluster at gap>50) on a single release date — one submission batch.

2. **It tests essentially ONE mechanism.** Determinants across the cohort are overwhelmingly
   `gyrA_T86I` (the QRDR point mutation). Near-perfection on a single well-catalogued point
   mutation is *mechanistically expected*; this is not evidence of broad catalog coverage.

3. **NOT clonality-corrected.** Prospective rows are not lineage-collapsed (the lineage table
   covers the provdisjoint cohorts only), and a national surveillance network is clonally
   structured. Treat 0.998 as isolate-level.

4. In-sample calibration origin still applies to the *rule*; this cohort tests it out-of-sample but
   does not re-derive it.

## What it is NOT, and what it IS

Not: source-diverse validation, lineage-independent validation, or clinical validation.

Is: an independent-source, temporally leakage-free, powered out-of-sample test of an opt-in rule
whose own text asked for exactly that — and it passes, **on a different institution and country
from the cohort it was calibrated against**. The provenance-disjoint Campylobacter cipro cell is
`PRJNA560409` (University of Melbourne, N=40, sens 1.000 / spec 1.000). Two internally-concentrated
but mutually independent sources now agree on near-perfect Campylobacter cipro performance. Each
alone fails the diversity bar; that they are *different* sources is the part worth keeping.

## Reproduce

```bash
uv run python -m scripts.prospective_lock_validate \
  --cohort-tsv "D:/dna_decode_cache/data files donwload/prospective_cohort.tsv" \
  --organism Campylobacter --amrfinder-organism Campylobacter --drug ciprofloxacin
uv run python scripts/build_validation_report_card.py
```
