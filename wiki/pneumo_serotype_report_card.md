# S. pneumoniae capsular serotype — decoder report card

**Decoder:** `dna-pneumo-serotype` (also `dna-decode pneumoserotype`) — deterministic cps-reference blastn caller.
**Trait class:** capsular serotype ("their look" / antigenic identity). Sibling of `dna-serotype` (E. coli O:H) + `dna-ktype` (Klebsiella capsule).
**Date:** 2026-06-24 (Quellung validation 2026-06-25). **Status:** caller SHIPPED; **real-DB reference-control 4/4 exact**; **INDEPENDENT phenotypic-Quellung validation IN PROGRESS** (GPS Poland cohort, n≤260; pilot below).

## INDEPENDENT measured-label validation vs phenotypic Quellung (2026-06-25) — the real number
The independent-measured-label win (the HIV pattern, for a typing trait): scored the deterministic caller vs
the **wet-lab phenotypic serotype** from the GPS pipeline paper (Nat Commun 2025, **Supplementary Data 1**,
column `Phenotypic_serotype`, method QUELLUNG/phenotypic) on **GPS-deposited ENA assemblies** — label AND
assembly both independent of our caller (clears the circularity rail; NOT the in-silico Monocle field).
Cohort: 260 Poland isolates (50 explicit-QUELLUNG + 210 phenotypic), 32 serotypes. Runner:
`scripts/pneumo_gps_quellung_validate.py` (ENA ERS→ERZ→contig.fa.gz, native blastn, checkpointed).

**Pilot (n=10):** serogroup concordance **0.9**, exact-serotype **0.4**. The exact misses are systematically
WITHIN-serogroup (9A↔9V, 6B↔6E, 15B↔15C) — the documented v0 ceiling (single-best cps reference resolves
serogroup; within-serogroup pairs need the allele-level logic the full tools add), NOT a bug. **Honest
headline: serogroup-level concordance is the v0's real resolution; exact-serotype is the lower bound that
motivates a v0.1** (allele-level within-serogroup typing). Full n≤260 number: `wiki/pneumo_serotype_cohort_validation.json` (run in progress).

## Real-DB reference-control (2026-06-24) — a REAL number
Built the real cps DB from PneumoCaT's Stage-1 reference (95 serotypes; `scripts/build_pneumo_cps_db.py`)
and ran the caller (native blastn) on **4 textbook-known reference genomes** (strain names header-verified
from ENA — this caught + excluded 4 wrong-accession fetches that returned *Beijerinckia*/*Korarchaeum*/*Leptothrix*):

| strain | measured serotype | predicted | exact | %id / %cov |
|---|---|---|---|---|
| TIGR4 | 4 | 04→4 | ✓ | 100 / 100 |
| D39 | 2 | 02→2 | ✓ | 100 / 100 |
| Hungary19A-6 | 19A | 19A | ✓ | 98.6 / 100 |
| ATCC 700669 (Spain23F-1) | 23F | 23F | ✓ | 99.7 / 100 |

**Reference-control: exact 4/4, serogroup 4/4.** This validates the real DB + caller integration end-to-end
on real genomes with INDEPENDENT textbook labels. **Honest scope:** n=4 reference-control, NOT the full GPS
Quellung cohort (the headline GREEN-VALIDATED number still needs the 11,810-isolate cohort run). Artifact:
`wiki/pneumo_serotype_reference_control_2026-06-24.json`.

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
- **Real cps DB:** NOT committed (gitignored-class external DB). Build path: derive a per-serotype `cps_references.fasta` (header `serotype__<ST>__<id>`) from PneumoCaT's `pneumo_capsular_locus_references` or SeroBA's database. **Acquisition note (verified 2026-06-24):** `git clone https://github.com/phe-bioinformatics/PneumoCaT` → its `streptococcus-pneumoniae-ctvdb` / capsular reference FASTAs → concatenate per-serotype references under the `serotype__<ST>__<id>` header convention. This is the MORE tractable of the two new cells' DBs (a flat per-serotype reference set, vs SeqSero2's algorithmic serovar logic).
- **Full-cohort GREEN-VALIDATED number (PENDING — runnable):** `scripts/serotype_cohort_validate.py --cell pneumo` — per GPS isolate with a Quellung label: fetch assembly → `call_pneumo_serotype` → concordance vs the phenotypic Quellung serotype. Multi-hour cohort op (best on D: / a long window; native blastn, no Docker). Reports serogroup-level + exact-serotype concordance separately (honest, given the v0 within-serogroup ceiling).

## Provenance / reproducibility
- Caller: `dna_decode/pneumoserotype/{runner,cli}.py`; thresholds identity 90 / coverage 70.
- Record schema: `pneumo-serotype-call-v0` (carries `caller_is_independent_baseline=False` + the caveat).
- FROZEN AMR surface (`amr_rules.py` + `calibrated_amr_rules.json`) byte-unchanged — this is an additive typing cell.
- **NOT a clinical tool.**

Sources: [GPS Pipeline, Nat Commun 2025](https://www.nature.com/articles/s41467-025-64018-5) · [GPS Pipeline (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12460886/).
