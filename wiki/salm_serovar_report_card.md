# Salmonella enterica serovar — decoder report card

**Decoder:** `dna-salmserovar` (also `dna-decode salmserovar`) — deterministic antigen-blastn + Kauffmann-White formula caller.
**Trait class:** serovar (the canonical Salmonella identity / "their look"). Sibling of `dna-serotype` (E. coli O:H) + `dna-ktype`.
**Date:** 2026-06-24 (real DB 2026-06-25). **Status:** caller SHIPPED + **REAL DB BUILT** + verified on a real reference genome (S. Typhimurium LT2 → **"Typhimurium"** 4:i:1,2, all 100%); full-cohort number = a runnable step.

## Real DB BUILT (2026-06-25) — the deferred data-engineering, done
The previously-deferred "real SeqSero2 DB" is built: `scripts/build_salmserovar_db.py` derives BOTH artifacts
from a SeqSero2 clone — **`serovar_table.tsv`** (2365 White-Kauffmann-Le Minor formulas, from `Initial_Conditions.py`'s
`phaseO`/`phase1`/`phase2`/`sero` parallel lists) + **`salmonella_antigens.fasta`** (360 alleles: 201 H1=fliC,
97 H2=fljB, 62 O-group wzx/wzy, reformatted to the `<axis>__<antigen>__<id>` convention). DB at
`data/salmserovar_db/` (gitignored; rebuild from a clone).
- **Real-genome verification:** `dna-salmserovar` on S. Typhimurium LT2 (ENA GCA_000006945.2) → **Typhimurium**
  (formula 4:i:1,2; O/H1/H2 all 100% id/cov). End-to-end correct.
- **Caller bug FIXED (was live in shipped 0.5.2):** `_best_per_axis` selected by coverage-only → flagellin
  alleles cross-hybridize at full coverage, so it picked the WRONG H antigen (LT2 gave 4:r:1,5,7). Fixed to
  **identity-primary** selection (the true antigen is the ~100%-identity hit). Regression test added.

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
- **Real antigen DB:** NOT committed (gitignored-class external DB). Build path: derive `salmonella_antigens.fasta` (headers `O__<g>__id` / `H1__<a>__id` / `H2__<a>__id`) + `serovar_table.tsv` (the White-Kauffmann-Le Minor formula table) from the SeqSero2 database. **Acquisition note (verified 2026-06-24 by cloning + inspecting the DB):** SeqSero2 is **bioconda-only (NOT on PyPI)**; `git clone https://github.com/denglab/SeqSero2` → `seqsero2_db/`. Its DB is **`H_and_O_and_specific_genes.fasta`** (368 seqs; headers encode the H1=fliC / H2=fljB antigen, e.g. `fliC_g,m_...`, `fljB_1,2_...` — cleanly extractable) **+ `antigens.pickle` (a per-allele k-mer DETECTION index, NOT a serovar formula table)**. So unlike PneumoCaT (which ships a flat per-serotype reference FASTA → the pneumo cell's real DB built in minutes), SeqSero2 has **no flat allele-DB + formula-TSV to adapt**: the O-antigen is detected via specific genes and the serovar is resolved by an ALGORITHM (an internal White-Kauffmann-Le Minor table in their Python). Building this cell's real DB = genuine data-engineering: (1) extract fliC/fljB H-antigen alleles (easy, from the FASTA), (2) assemble an O-antigen allele set + (3) source the antigenic-formula→serovar TSV (e.g. from the published White-Kauffmann-Le Minor scheme / SISTR). Deferred — best done attended.
- **Full-cohort GREEN-VALIDATED number (PENDING — runnable):** `scripts/serotype_cohort_validate.py --cell salm` — per lab-serotyped isolate: fetch assembly → `call_serovar` → concordance vs the wet-lab serovar. Native blastn, no Docker. Reports formula-resolved-rate + serovar concordance separately.

## Provenance / reproducibility
- Caller: `dna_decode/salmserovar/{runner,cli}.py`; thresholds identity 90 / coverage 80.
- Record schema: `serovar-call-v0` (carries `caller_is_independent_baseline=False` + the caveat).
- FROZEN AMR surface byte-unchanged — additive typing cell.
- **NOT a clinical tool.**

Sources: [SeqSero2, ASM AEM](https://journals.asm.org/doi/10.1128/aem.01746-19) · [WGS Salmonella serotyping vs gold-standard, Frontiers 2025](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1685741/full).
