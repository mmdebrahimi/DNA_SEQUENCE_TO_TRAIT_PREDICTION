# Why the gentamicin gap survived validation: the cell that should have caught it was 95% one BioProject

The report card discloses that the SCORED cells are **clonally** dominated. This measures the sibling
property one level up — **source** concentration — and it turns out to explain a finding this project
spent three runs chasing.

## The measurement

| cell | n | BioProjects | largest share | dominant source |
|---|---:|---:|---:|---|
| campylobacter × ciprofloxacin | 40 | **1** | **100%** | PRJNA560409 |
| e. coli × ciprofloxacin | 60 | 2 | **97%** | PRJNA278886 |
| e. coli × gentamicin | 60 | 4 | **95%** | PRJNA278886 |
| e. coli × ceftriaxone | 60 | 4 | 55% | PRJNA278886 |
| e. coli × tetracycline | 60 | 6 | 53% | PRJNA662792 |
| klebsiella × ceftriaxone | 60 | 2 | 50% | PRJNA278886 |
| klebsiella × gentamicin | 60 | 32 | 48% | PRJNA278886 |
| klebsiella × tetracycline | 60 | 2 | 50% | PRJNA278886 |
| klebsiella × meropenem | 60 | 15 | 37% | PRJNA504784 |
| klebsiella × ciprofloxacin | 60 | **36** | **18%** | PRJNA529587 |

**3 of 10 cells are effectively single-source.** Others are genuinely diverse — Klebsiella cipro draws on
36 BioProjects across 60 isolates. The trust surface currently does not distinguish these.

## The consequence, measured

The E. coli × gentamicin cell reports **sens 0.893**. Two independent source-diverse measurements of the
same organism and drug with the same frozen rule report **0.429** (prospective accrual) and **0.523**
(131 accession-disjoint isolates across 8 BioProjects).

Here is why:

| | n | `rmt`-family carriers |
|---|---:|---:|
| the SCORED validation cohort (95% one BioProject) | 60 | **0** |
| a source-diverse disjoint set | 131 | **37 (28%)** |

**The catalog gap was not missed — it was structurally invisible.** A validation cohort containing zero
carriers of an entire determinant family cannot detect a rule that is blind to that family. The 0.893 was
never wrong about its cohort; it was a statement about one hospital's isolates.

## The direction of the error is not predictable

Single-source cells are unreliable **both ways**, which is the reason to disclose rather than to correct:

- **e. coli × gentamicin** (95% one source) reads **optimistic** — 0.893 vs 0.523 on diverse data.
- **e. coli × ciprofloxacin** (97% one source) reads **pessimistic** — spec 0.700 vs 0.988 on an
  8-BioProject set.

A single-site estimate is not biased in a knowable direction; it is simply an estimate of that site.

## What this is NOT

- **Not a demotion.** The cells remain provenance-**disjoint** from the tuning data, which is what was
  claimed and what was verified. This measures whether they are also provenance-**diverse** — a different
  property that was never claimed and is now measured.
- **Not a re-scoring.** No published number changes. `sens 0.893` remains the correct number for that
  cohort.
- **Not a bug in the cohort builder.** Provenance-disjointness was the design goal and it was met.

## A hypothesis of mine died on the way here

I expected my 131-isolate set to score better because it was *easier* — accession-disjoint but
in-distribution. Provenance said the opposite: mine spans **8 BioProjects**, the report-card cohort spans
**2**, and their BioProject overlap is **1**. The comparison set was the more diverse one, which is what
makes its numbers informative rather than suspect.

## Recommendation

Source concentration belongs on the trust surface next to clonality — the same disclose-don't-demote
shape, in its own namespace. A reader should be able to see that Klebsiella cipro rests on 36 sources and
Campylobacter cipro on one, because those two numbers do not mean the same thing.

Reproduce: `uv run python scripts/provdisjoint_source_concentration.py` (network; writes
`wiki/provdisjoint_source_concentration.json`).
