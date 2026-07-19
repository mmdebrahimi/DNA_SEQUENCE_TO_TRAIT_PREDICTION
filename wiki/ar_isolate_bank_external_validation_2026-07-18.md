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

## Provenance-disjointness — the gate worked, and a methodological catch

The exact-set BioSample-level preflight (`external_cohort_preflight.py --cohort-manifest …`) FAILED
closed on the full 114-isolate cohort — correctly — for two reasons:

1. **14 confirmed same-isolate leaks.** 14 AR Bank E. coli BioSamples are ALSO in the decoder's
   tuning set (expected: CDC deposits AR Bank isolates into NCBI-PD, which the frozen cells tuned on).
   This is the GCA-vs-BioSample same-isolate overlap that accession-string matching misses — the
   BioSample-level gate caught it.
2. **40.6% of the 1509 tuning accessions unresolved** to BioSamples — a transient NCBI-throttling
   artifact during the run (1508/1509 are proper GCA/GCF, so it is NOT structural; the resolver
   caches successes and a warm re-run would retry the failures).

**The catch (and a preflight-hardening insight):** because 40% of the tuning GCAs were throttle-
unresolved, the BioSample-level check *could not see* whether those tuning isolates matched the
cohort. A complementary **resolution-free** check — cohort assembly accession-base ∩ tuning GCA
accession-base (needs no Entrez resolution, cannot be throttled) — found **9 ADDITIONAL leaks** the
BioSample preflight structurally missed. So the fail-closed behavior was doubly right: forcing a
"pass" would have let those 9 same-assembly isolates leak in. **Recommended hardening:**
`external_cohort_preflight` should add this direct cohort-assembly-base ∩ tuning-GCA-base test
alongside the resolution-dependent one (surfaced here, not yet applied — it touches the shared arm
used by Oxford too).

**Full leaked set = 14 (BioSample) + 9 (assembly) = 23**, excluded via the new
`build_ar_bank_labels.py --exclude-biosamples`. The rebuilt cohort is disjoint at BOTH the BioSample
and assembly level (`wiki/cohort_manifest_external_arbank_disjoint_2026-07-18.json`, 91 BioSamples;
`leakage_excluded` recorded in the manifest).

## The scorable disjoint cohort

Only 25 of the 91 fully-disjoint isolates have a downloadable assembly (66 are reads-only /
ASSEMBLY-REQUIRED — AR Bank deposits SRA reads but not always an assembled GCA). The FREE ∩
fully-disjoint scorable set, strict-tier:

| Drug | R | S | Powered side (informative n) |
|---|---|---|---|
| ceftriaxone | 19 | 2 | **SENSITIVITY** (n=19 ✓) |
| ciprofloxacin | 18 | 1 | **SENSITIVITY** (n=18 ✓) |
| gentamicin | 0 | 15 | **SPECIFICITY** (n=15 ✓) |

Each drug's informative class clears the ≥10 powering floor — a legitimate one-sided independent
external test survives even after removing all 23 leaks and the reads-only isolates.

## The scoring run — BLOCKED on Docker (external infra, not code)

Scoring = the shipped external arm: resolve each of the 25 BioSamples → GCA, run **AMRFinder
(Docker)** with the frozen organism triple (`-O Escherichia`,
`call_resistance(organism="Escherichia_coli_Shigella")`), compute the confusion vs the MIC labels.
25 genomes × ~95 s ≈ ~40 min (far less than the earlier 115-isolate estimate).

**Status 2026-07-18: `docker ps` → DOCKER_DOWN.** Docker Desktop is not running on this host, so
AMRFinder cannot execute. This is the classic WSL-mount wedge (restartable; `wsl --shutdown` +
relaunch Docker Desktop recovers it). The gap is **external** (start Docker), not code-closable —
the ingester, disjoint cohort, labels, and manifest are all built and committed.

When Docker is up, score each powered drug (one-sided). The fully-disjoint cohort was already
excluded of all 23 leaks, so `--allow-degraded` is justified (the preflight's residual FAIL is the
throttle-induced unresolved-fraction, and disjointness was independently verified at the assembly
level above); prefer a clean warm-preflight PASS if the re-run resolves the unresolved fraction.

```bash
# ensure Docker Desktop is running first (wsl --shutdown to recover a wedge)
M=wiki/cohort_manifest_external_arbank_disjoint_2026-07-18.json
uv run python -m scripts.external_cohort_revalidate --cohort ar_bank_ecoli --drug ceftriaxone \
  --labels-dir data/raw/ar_bank_ecoli_extval_ceftriaxone --cohort-manifest $M \
  --allow-degraded --min-per-class 2            # sensitivity arm (19 R)
# repeat --drug ciprofloxacin (18 R, sensitivity) and --drug gentamicin (15 S, specificity)
# then roll-up: scripts/build_external_validation_report.py (separate external namespace)
```

## Honest scope

- The DATA fetch is free + attestation-free; **physical isolates** (not needed here) require an
  institutional account + biosafety officer.
- The AR Bank MIC is CDC reference BMD (real G1). Isolates are curated (resistance-enriched) → **not
  prevalence-representative**; report one-sided powered metrics, never a population sens/spec off this set.
- Independence is **provenance-separable at the BioSample level** (preflight-enforced), not
  methodology-independent (same AMRFinder `-O` + frozen `call_resistance`).
- Frozen decoder surface byte-unchanged; this is a READ-only external-validation adapter.
