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
- **v0 resolution (CORRECTED 2026-09-09 — this line asserted uniqueness and was FALSE as shipped):** the serovar is reported when the O:H1:H2 formula **RESOLVES**, which is **not** the same as resolving uniquely. The table is built with `ambiguous_policy="first"`, so for the **195 formulas carried by >1 Kauffmann-White serovar** (after canonical naming) it records an **arbitrary winner** instead of abstaining — kept deliberately, because the principled `omit` alternative was measured at **0 hits gained and 3 lost**. On top of that a **phase-incomplete fallback** matches `O+H1` ignoring H2 and reports a serovar iff that collapses to one name, which may itself be such a winner. So a returned serovar is one of: uniquely resolved · an arbitrary winner on a contested formula · fallback-resolved and possibly winner-contaminated — and **v0 does not yet disclose which**. serovar=None still means the formula did not resolve at all.

## Validation status
- **Synthetic control (committed, offline-safe):** `tests/test_salmserovar.py` — a synthetic Typhimurium fixture (O=4 / H1=i / H2=1,2) → formula `4:i:1,2` → `Typhimurium` (real blastn) + offline-safe degrade + pure-logic parsers. Always-green in CI without the real DB.
- **Real antigen DB:** NOT committed (gitignored-class external DB). Build path: derive `salmonella_antigens.fasta` (headers `O__<g>__id` / `H1__<a>__id` / `H2__<a>__id`) + `serovar_table.tsv` (the White-Kauffmann-Le Minor formula table) from the SeqSero2 database. **Acquisition note (verified 2026-06-24 by cloning + inspecting the DB):** SeqSero2 is **bioconda-only (NOT on PyPI)**; `git clone https://github.com/denglab/SeqSero2` → `seqsero2_db/`. Its DB is **`H_and_O_and_specific_genes.fasta`** (368 seqs; headers encode the H1=fliC / H2=fljB antigen, e.g. `fliC_g,m_...`, `fljB_1,2_...` — cleanly extractable) **+ `antigens.pickle` (a per-allele k-mer DETECTION index, NOT a serovar formula table)**. So unlike PneumoCaT (which ships a flat per-serotype reference FASTA → the pneumo cell's real DB built in minutes), SeqSero2 has **no flat allele-DB + formula-TSV to adapt**: the O-antigen is detected via specific genes and the serovar is resolved by an ALGORITHM (an internal White-Kauffmann-Le Minor table in their Python). Building this cell's real DB = genuine data-engineering: (1) extract fliC/fljB H-antigen alleles (easy, from the FASTA), (2) assemble an O-antigen allele set + (3) source the antigenic-formula→serovar TSV (e.g. from the published White-Kauffmann-Le Minor scheme / SISTR). Deferred — best done attended.
- **Full-cohort GREEN-VALIDATED number — LANDED 2026-09-04, refined 2026-09-09 (this line read "PENDING — runnable" until 2026-09-09; under-claiming is as much a trust-surface falsehood as over-claiming).** Scored on **200 wet-lab slide-agglutination-labelled isolates** (74 serovars / 29 BioProjects / largest share 0.125): **accuracy 0.8400** (168 hit / 11 miss / 21 no-call) after the O-procedure port + name-table rebuild, up from 0.7050. Pinned **SeqSero2 1.3.2** scored on the SAME isolates through the SAME equivalence function reaches **0.8800** → **delta −0.0400**. **Quote the delta with more confidence than the levels** (the per-serovar cap of 12 over-weights rare serovars and depresses every caller — SeqSero2's own published figure is ~0.95). Artifacts: `wiki/salmserovar_validation_2026-09-04.{md,json}` · `salmserovar_seqsero2_2026-09-08.{md,json}` · `salmserovar_o_fix_result_2026-09-09.{md,json}` · `salmserovar_name_table_2026-09-09.md`. **STRATIFIED 2026-09-09** (`wiki/salmserovar_resolution_strata_2026-09-09.{md,json}`) — the headline SURVIVES: **165 of 168 hits (98.2%) are uniquely resolved**; only **3** rest on an arbitrary tie-break, and the `O+H1` fallback produced **zero** calls. Tie-break-independent accuracy = **0.8250**. Quote 0.8400 with the stratum disclosed, or 0.8250 for the tie-break-independent figure. **Recorded against it:** the arbitrary-winner policy gains those 3 hits by converting **5 abstentions into confidently-wrong calls** (hit rate 3/8 within that stratum) — which is the losing side of the asymmetry this same cell pre-registered on 2026-09-04 (*"an abstention is recoverable by a human, a confident wrong serovar is not"*). Whether that reading should govern is an **open acceptance-bar question, not settled**.

## Provenance / reproducibility
- Caller: `dna_decode/salmserovar/{runner,cli}.py`; thresholds identity 90 / **coverage 40** (lowered from 80 on 2026-09-04 against a pre-registered bar; this line still said 80 until 2026-09-09 — the same prose-drift that let the CLI ship a stale default, surviving here because the drift guard resolves argparse defaults and cannot see a threshold restated in prose).
- Record schema: `serovar-call-v0` (carries `caller_is_independent_baseline=False` + the caveat).
- FROZEN AMR surface byte-unchanged — additive typing cell.
- **NOT a clinical tool.**

Sources: [SeqSero2, ASM AEM](https://journals.asm.org/doi/10.1128/aem.01746-19) · [WGS Salmonella serotyping vs gold-standard, Frontiers 2025](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1685741/full).
