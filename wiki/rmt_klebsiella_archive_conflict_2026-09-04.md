# The two archives disagree about Klebsiella — and the diverse one shows no over-call

Yesterday's Klebsiella over-call (BV-BRC, PPV 0.475) failed this project's own source-diversity bar. The
obvious next question — **what does the other archive say about Klebsiella specifically?** — had never
been asked, because the NCBI-PD sweep was read as a pooled number and the per-organism split sat
unexamined in the committed artifact.

It says the opposite.

| archive | R | S | PPV | sources | largest share | clears our 0.60 bar? |
|---|---|---|---|---|---|---|
| **NCBI-PD** *(artifact project excluded)* | 53 | **0** | **1.000** | **12** | **0.264** | **YES** |
| BV-BRC *(dominant study kept)* | 58 | 64 | 0.475 | 6 | 0.664 | no |

**The archive that clears the diversity bar contradicts the over-call.** If the true Klebsiella PPV were
0.475, observing zero susceptible carriers among 12 independent BioProjects has probability
**p = 1.3 × 10⁻⁴**.

> **Quote the conservative number.** The per-isolate figure is 7.7 × 10⁻¹⁸, but it assumes carriers are
> independent, which clonality makes false. Collapsing to **one vote per BioProject** throws away all
> within-project replication and asks only whether 12 independent sources all came back clean. That is
> the defensible bound, and it is still decisive.

## Why the exclusions are not special pleading

This is the part a reader should be most suspicious of, so it goes first. Both archives had a dominant
source. Both got the **same pre-registered `aac(3)` control**, written before either result was seen:

- **PD's `PRJNA1322038` FAILED it** — calls `aac(3)` carriers R **2%** of the time vs **97%** elsewhere,
  and calls isolates with no determinant R 86% of the time. Excluded as a label artifact.
- **BV-BRC's pmid 36801013 PASSED it** — 99% vs 83%, and 4% for no-determinant isolates. Kept.

The asymmetry is the **control's output**, not a choice made after seeing which way it cut. Had the
verdicts been reversed, the exclusions would have been reversed with them.

## What this changes

**The over-call is not retracted.** The warning stays and the safe direction stands. What changes is that
it is now **contradicted by the more source-diverse evidence** and rests on a single study. The shipped
`organism_scope` warning carries both sides — `one_line` refuses to state the PPV without naming the
contradicting archive, pinned by test.

**The E. coli scope is untouched and further reinforced:** 70R/0S in PD, 12R/0S in BV-BRC, 12/12 and
146/146 across the earlier sweeps, plus <1 in 4,979 prevalence in Oxford.

## What this hinges on, and how it could flip

**PD's 42 susceptible Klebsiella carriers all sit inside the excluded project.** So PD's zero is a
*zero-after-exclusion*. **If that exclusion were wrong, PD reads 53R/42S = PPV 0.558 and would
corroborate the over-call instead of contradicting it.** The entire result turns on one control call.

Other limits:

- **Absence is weaker than presence.** 53R/0S bounds the over-call; it does not prove the rule safe in
  Klebsiella.
- **Clonality is not corrected within projects.** The per-source collapse is a blunt proxy for a
  lineage-collapsed analysis, not a substitute for one.
- **Different label provenance** — PD's `AST_phenotypes` versus BV-BRC's publication-curated MICs. A
  systematic difference between those two label sources is not excluded, and would be the natural
  explanation if both archives are internally sound.

## The methodological point

The source-diversity bar was applied to the evidence *against* the rule first, because that is the
evidence that had just been published. Applying it only in that direction would be motivated reasoning,
so it was then applied to the evidence *for* the rule: **17 sources, largest share 0.273 — it passes.**
The bar is not being used selectively; it simply cuts differently on the two sides.

## Reproduce

```bash
uv run python scripts/rmt_klebsiella_archive_conflict.py
```

Offline — reads two committed artifacts. No frozen file is touched.
