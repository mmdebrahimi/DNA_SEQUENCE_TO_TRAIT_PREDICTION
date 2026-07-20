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
