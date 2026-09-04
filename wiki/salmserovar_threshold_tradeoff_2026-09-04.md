# The threshold change pays for itself — measured against a bar set before the run

Three causal claims about this cell were asserted and then measured wrong in a row: phase-2 flagellin
as the dominant defect, then DB coverage as the real priority, then *nearly* "so lower the threshold".
Each was reached by reasoning from the previous result rather than by measuring. **So this one fixed its
verdict rule before reading any number.**

## Pre-registered rule (frozen before the sweep)

```
ADOPT     if rescued_correct >= 7  AND  newly_wrong <= 2  AND  net >= +5
REJECT    if newly_wrong > 2, regardless of how much is rescued
NO_CHANGE otherwise
```

The asymmetry is deliberate: **an abstention is recoverable by a human, a confident wrong serovar is
not.** The rule demands many rescues and tolerates almost no new errors.

## Result — 200 isolates, coverage cut 80 → 40

| baseline → relaxed | n | |
|---|---|---|
| correct → correct | 98 | |
| **abstain → correct** | **36** | **gain** |
| wrong → wrong | 35 | |
| abstain → abstain | 20 | |
| wrong → correct | 7 | bonus the rule didn't count |
| abstain → wrong | 3 | |
| **correct → wrong** | **1** | **loss** |

**36 rescued against 1 new error, net +35.** Every bar cleared. **Verdict: ADOPT.**

## Why it is safe, not merely favourable

The usual selective-classification trade is coverage *against* accuracy — abstain less, be wrong more.
That did not happen here:

| | coverage | accuracy on covered | forced-call accuracy |
|---|---|---|---|
| deployed (80) | 0.705 | 0.702 | 0.495 |
| **relaxed (40)** | **0.900** | **0.783** | **0.705** |

**Both improved.** The abstentions being converted were overwhelmingly correct calls being discarded,
not coin-flips being admitted — consistent with the earlier measurement that 14 of 21 sub-threshold O
hits named the *correct* O group with zero wrong.

## What moved, in cell terms

| | before | after |
|---|---|---|
| accuracy vs wet-lab label | 0.702 | **0.783** |
| abstention rate | 29.5% | **10.0%** |
| delta vs NCBI-PD's in-silico field (0.925) | −0.223 | **−0.142** |

Still behind the incumbent field. **The cell is better, not fixed** — and "incumbent" here means NCBI's
published production call, not a pinned reference-tool run (see the comparator correction in the
[validation memo](salmserovar_validation_2026-09-04.md)).

## Why coverage and not identity

Identity stays at 90. The discarded O-antigen hits were at **near-perfect identity (median 99.8) and
partial coverage (median 58.4, max 78.9)** — a partial-alignment / allele-length mismatch concentrated
on the O7 `wzx/wzy` reference. Identity was never the failing axis, and relaxing it would admit
genuinely different alleles.

## Honest limits

- **Only the coverage cut moved.** A joint identity×coverage sweep was not run.
- One cohort (N=200, reference-lab-filtered), one antigen-DB build.
- Equivalence uses the same notation + White-Kauffmann rule as the original validation, applied
  identically to both settings, so no leniency can favour either.
- **40 is not shown to be optimal** — it is shown to clear a pre-registered bar. A sweep over several
  cuts would find the best value; this answers "does relaxing pay?", not "what is the best cut?".
- **ADOPT does not mean the O-antigen defect is closed.** The allele-length hypothesis — 11 of 14
  sub-threshold hits were a single O group — is still untested, and 20 isolates still abstain.

## Reproduce

```bash
uv run python scripts/salmserovar_threshold_tradeoff.py
```

Needs blastn + the cached assemblies. Frozen AMR surface byte-unchanged — typing cell.
