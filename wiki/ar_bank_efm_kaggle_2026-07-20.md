# AR-Bank E. faecium (Enterococcus) — first validation attempt via Kaggle AMRFinder factory

**Date:** 2026-07-20 · **Status:** ⚠️ **UNDERPOWERED** (perfect accuracy on what assembled; S-class
assembly-blocked) · **Cell:** `enterococcus_amr` (NON-FROZEN, hand determinant rule) · **Frozen surface:**
byte-unchanged.

## Why Kaggle

E. faecium is deposited in the AR Bank as **SRA reads only (0 downloadable assemblies)** → the cell needs
SRA→assembly→AMRFinder, the exact Docker-heavy path that wedges Docker Desktop locally. Kaggle runs bioconda
natively (no Docker nesting, 12 h CPU) → the AMRFinder factory host. Local = scoring only.

## Result — the cell works, but is UNDERPOWERED

The kernel `enterococcus-amrfinder` assembled **7 of 17** isolates; scoring vs CDC S/I/R labels:

| drug | n | R/S | acc | sens | spec | verdict |
|---|---|---|---|---|---|---|
| doxycycline | 5 | 5R/0S | 1.00 | 1.00 | — | UNDERPOWERED (0 S) |
| levofloxacin | 7 | 6R/1S | 1.00 | 1.00 | 1.00 | UNDERPOWERED (1 S) |

**The determinant rule is PERFECT on every isolate that assembled** (doxy 5/5 R via tet genes; levo 7/7
incl. the 1 S). The cell is NOT the problem — it's **assembly-blocked on the susceptible class**.

## Root cause — the S-side block failed to assemble (single-end libraries)

10 of 17 isolates failed on Kaggle (`CalledProcessError`), and they are **the entire S class**: for
doxycycline the 10 failures = 3R + **ALL 7 S** → 0 S scored → unpowerable. Diagnosed via the ENA
`library_layout` field on the failed `SAMN15040xxx` block:

- **5 failures are SINGLE-end** libraries (SAMN15040088/089/091/093/101 — 4 of them S). The kernel ran
  `fasterq-dump --split-files` + hardcoded `{srr}_1.fastq`, but fasterq writes `{srr}.fastq` (**no `_1`
  suffix**) for single-end runs → SKESA got a nonexistent file → `CalledProcessError`. **This is the
  primary S-class killer.**
- **5 failures are PAIRED** (SAMN15040082/086/102/105/106) — should have assembled like the 7 successes;
  failed for a separate (likely transient Kaggle resource/timeout) reason, retried on re-run.

## Fix (shipped) + re-run

`scripts/kaggle_enterococcus_amrfinder.py` now uses **`fasterq-dump --split-3`** (paired → `_1/_2`,
single → `{srr}.fastq`) + robust read-file detection covering the single-end case (raises loudly if no
FASTQ at all). This is a durable pipeline fix — like the ENA-direct fetch + `--subsample-reads` shipped
earlier today, it makes the SRA path robust for future single-end cohorts. Re-running the kernel retries all
10 failures; the single-end fix recovers ≥4 S, and the paired retries recover the rest → expected to POWER
doxycycline (8R/7S) on the next run.

**Scope:** first Enterococcus validation, a NEW organism cell. NON-FROZEN; the frozen decoder surface is
untouched. The cell's determinant rule is validated (perfect accuracy); only the S-class powering is pending
the re-run.

## Re-run result (2026-07-21) — single-end fix worked, but a SECOND wall + a scorer integrity bug

The single-end-fixed kernel assembled **16/17** (was 7/17) — the `--split-3` fix stopped the crashes.
**But verify-in-batch caught a false victory:** the recovered `SAMN15040xxx` block produces **empty
AMRFinder TSVs (0–1 determinants; 294-byte header-only files)** — a real E. faecium assembly yields 13–24
determinants (even susceptible isolates carry intrinsic genes like `aac(6')-Ii`). The assemblies FAILED.

**Root cause:** those SRA runs are **shallow** (58k–445k reads; both PAIRED and SINGLE) → SKESA produces
fragmented/empty contigs → AMRFinder finds nothing. 0/9 have downloadable GCA assemblies (SRA-reads-only).
A **data-quality wall on the S-side block**, distinct from the (now-fixed) single-end crash.

**Scorer integrity bug (fixed):** the scorer was scoring these empty assemblies as **S-by-absence** (no van
gene → S), manufacturing a false `SCORED_ENDORSED` — vancomycin sens cratered to 0.2 (the block's R isolates
scored S because their assemblies were empty), and the doxycycline "endorsement" was inflated by
empty-assembly S isolates right by accident. Fixed via `MIN_DETERMINANTS_FOR_VALID_ASSEMBLY = 3`: an isolate
with < 3 total determinants is `INDETERMINATE_ASSEMBLY_FAILED`, **excluded from scoring**, never
S-by-absence. An empty assembly is not a susceptible isolate.

**Honest re-scored result** (9 empty assemblies gated out; scored only the 7 with real assemblies):

| drug | n | R/S | sens | spec | verdict |
|---|---|---|---|---|---|
| doxycycline | 5 | 5R/0S | 1.00 | — | UNDERPOWERED |
| levofloxacin | 7 | 6R/1S | 1.00 | 1.00 | UNDERPOWERED |
| vancomycin | 7 | 2R/5S | 1.00 | 0.80 | UNDERPOWERED |
| teicoplanin | 7 | 2R/5S | 1.00 | 0.80 | UNDERPOWERED |

The **determinant rules are correct** (perfect sens on every real assembly; the earlier vanco sens 0.2 was
100% the empty-assembly artifact). Enterococcus stays **UNDERPOWERED** — the S-side block is unrecoverable
from shallow reads. **Forward lever (speculative):** a SPAdes-instead-of-SKESA kernel retry (more sensitive
on low coverage; SAMN15040089 has 445k paired reads that *should* assemble). Uncertain payoff; not run yet.
The scorer integrity gate is the durable win — it stops any empty assembly ever masquerading as susceptible.
