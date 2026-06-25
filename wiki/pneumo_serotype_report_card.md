# S. pneumoniae capsular serotype — decoder report card

**Decoder:** `dna-pneumo-serotype` (also `dna-decode pneumoserotype`) — deterministic cps-reference blastn caller.
**Trait class:** capsular serotype ("their look" / antigenic identity). Sibling of `dna-serotype` (E. coli O:H) + `dna-ktype` (Klebsiella capsule).
**Date:** 2026-06-24. **Status:** caller SHIPPED + offline-safe; validated on a synthetic cps control; full-cohort number = a runnable step (see below).

## GREEN-cell gate (from `plans/Non_AMR_GREEN_Cell_Triage_Round2_2026-06-24.md`)
| Gate | Result |
|---|---|
| **1. Determinant catalog exists?** | ✓ curated cps-locus reference sets — PneumoCaT (`phe-bioinformatics/PneumoCaT`) / SeroBA (`sanger-pathogens/seroba`) / the GPS pipeline; ~90–107 serotypes. |
| **2. FREE, independent, MEASURED isolate-level label?** | ✓✓ **Global Pneumococcal Sequencing (GPS): 11,810 genomes with phenotypic Quellung serotype** records, ENA-public. The richest free measured label in this project (vs the Klebsiella K-antigen cell's 733). |
| **Verdict** | **GREEN-VALIDATED candidate** — passes both gates; the gold-standard shape (ship + score vs the measured label + trust badge). |

## Honesty tier
- **`caller_is_independent_baseline = False`** — the v0 is FAITHFUL to the cps-reference typing method (blastn best-match over the curated reference DB). Validating it against an in-silico tool (SeroBA/PneumoCaT) would be in-distribution, NOT the GREEN-VALIDATED tier.
- **The GREEN-VALIDATED number must be scored vs the wet-lab MEASURED Quellung label**, not vs another genomic tool. Circularity rail (load-bearing; from the VF-diff lesson).
- **v0 resolution ceiling:** single-best-reference resolves **serogroup** reliably; within-serogroup pairs that differ by a single locus/SNP (6A/6B at wciP, 19A/19F) need the allele-level logic the full tools add. Published **in-silico-vs-Quellung concordance ≈ 89.3%** (GPS pipeline, n=10,549/11,810 concordant; Nat Commun 2025) — the ceiling for any cps-based caller to match.

## Validation status
- **Synthetic control (committed, offline-safe):** `tests/test_pneumoserotype.py` — a synthetic 19F cps fixture → caller returns `19F` (real blastn) + the offline-safe degrade + pure-logic parsers. Always-green in CI without the real DB.
- **Real cps DB:** NOT committed (gitignored-class external DB). Build path: derive a per-serotype `cps_references.fasta` (header `serotype__<ST>__<id>`) from PneumoCaT's `pneumo_capsular_locus_references` or SeroBA's database.
- **Full-cohort GREEN-VALIDATED number (PENDING — runnable):** `scripts/serotype_cohort_validate.py --cell pneumo` — per GPS isolate with a Quellung label: fetch assembly → `call_pneumo_serotype` → concordance vs the phenotypic Quellung serotype. Multi-hour cohort op (best on D: / a long window; native blastn, no Docker). Reports serogroup-level + exact-serotype concordance separately (honest, given the v0 within-serogroup ceiling).

## Provenance / reproducibility
- Caller: `dna_decode/pneumoserotype/{runner,cli}.py`; thresholds identity 90 / coverage 70.
- Record schema: `pneumo-serotype-call-v0` (carries `caller_is_independent_baseline=False` + the caveat).
- FROZEN AMR surface (`amr_rules.py` + `calibrated_amr_rules.json`) byte-unchanged — this is an additive typing cell.
- **NOT a clinical tool.**

Sources: [GPS Pipeline, Nat Commun 2025](https://www.nature.com/articles/s41467-025-64018-5) · [GPS Pipeline (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12460886/).
