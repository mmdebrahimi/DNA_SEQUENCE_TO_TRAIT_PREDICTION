# Salmonella enterica serovar — decoder report card

**Decoder:** `dna-salmserovar` (also `dna-decode salmserovar`) — deterministic antigen-blastn + Kauffmann-White formula caller.
**Trait class:** serovar (the canonical Salmonella identity / "their look"). Sibling of `dna-serotype` (E. coli O:H) + `dna-ktype`.
**Date:** 2026-06-24. **Status:** caller SHIPPED + offline-safe; validated on a synthetic Typhimurium control; full-cohort number = a runnable step (see below).

## GREEN-cell gate (from `plans/Non_AMR_GREEN_Cell_Triage_Round2_2026-06-24.md`)
| Gate | Result |
|---|---|
| **1. Determinant catalog exists?** | ✓ antigen-formula DBs — SeqSero2 (`denglab/SeqSero2`) / SISTR; the White-Kauffmann-Le Minor scheme (O-group + H1=fliC + H2=fljB → serovar). |
| **2. FREE, independent, MEASURED isolate-level label?** | ✓ traditional Kauffmann-White slide-agglutination serotyping, reported for large public sets (NARMS / PulseNet / the SeqSero2 + Frontiers evaluation cohorts) — **must filter to wet-lab-serotyped, NOT tool-predicted, isolates**. |
| **Verdict** | **GREEN-VALIDATED candidate (with circularity filter)** — passes both gates; the canonical Salmonella identity trait. |

## Honesty tier
- **`caller_is_independent_baseline = False`** — FAITHFUL to the SeqSero2 / Kauffmann-White method (blastn over the antigen allele DB + formula lookup).
- **The GREEN-VALIDATED number must be scored vs the wet-lab MEASURED serovar** (traditional serotyping), not vs SeqSero2/SISTR predictions — the dominant circularity trap here, since many public "serovars" ARE tool-predicted. The validation cohort must be filtered to lab-serotyped isolates.
- **v0 resolution:** the serovar is reported **only when the O:H1:H2 formula resolves uniquely** in the White-Kauffmann table (else the antigenic FORMULA is reported with serovar=None — like O?/H?). Shared formulas / phase-incomplete genomes → formula-only, honestly.

## Validation status
- **Synthetic control (committed, offline-safe):** `tests/test_salmserovar.py` — a synthetic Typhimurium fixture (O=4 / H1=i / H2=1,2) → formula `4:i:1,2` → `Typhimurium` (real blastn) + offline-safe degrade + pure-logic parsers. Always-green in CI without the real DB.
- **Real antigen DB:** NOT committed (gitignored-class external DB). Build path: derive `salmonella_antigens.fasta` (headers `O__<g>__id` / `H1__<a>__id` / `H2__<a>__id`) + `serovar_table.tsv` (the White-Kauffmann-Le Minor formula table) from the SeqSero2 database. **Acquisition note (verified 2026-06-24):** SeqSero2 is **bioconda-only (NOT on PyPI)** — get the DB via `git clone https://github.com/denglab/SeqSero2` (or `conda install -c bioconda seqsero2`); the DB→`O__/H1__/H2__` + `serovar_table.tsv` adaptation is the remaining data-engineering step (SeqSero2's serovar logic is an algorithm, not a flat formula table, so the table needs deriving).
- **Full-cohort GREEN-VALIDATED number (PENDING — runnable):** `scripts/serotype_cohort_validate.py --cell salm` — per lab-serotyped isolate: fetch assembly → `call_serovar` → concordance vs the wet-lab serovar. Native blastn, no Docker. Reports formula-resolved-rate + serovar concordance separately.

## Provenance / reproducibility
- Caller: `dna_decode/salmserovar/{runner,cli}.py`; thresholds identity 90 / coverage 80.
- Record schema: `serovar-call-v0` (carries `caller_is_independent_baseline=False` + the caveat).
- FROZEN AMR surface byte-unchanged — additive typing cell.
- **NOT a clinical tool.**

Sources: [SeqSero2, ASM AEM](https://journals.asm.org/doi/10.1128/aem.01746-19) · [WGS Salmonella serotyping vs gold-standard, Frontiers 2025](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1685741/full).
