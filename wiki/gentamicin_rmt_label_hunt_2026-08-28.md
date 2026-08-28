# The missing measurement, found: 62/63 of `rmt`-family carriers are gentamicin-R on public labels

Last run left one thing open — the candidate's over-calling risk was **UNTESTED, not zero**, because every
methyltransferase carrier in the local labelled data was R. The binding constraint was S-labelled carriers.
This run went and got the labels: NCBI Pathogen Detection has `AST_phenotypes` for many of the 109 carriers
whose genomes and AMRFinder output are already cached here. Free, and the same source the frozen SCORED
cells came from.

## What the sweep found (COMPLETE — all four organism groups streamed cleanly)

| | |
|---|---:|
| carriers with cached AMRFinder | 109 |
| of those, with a public gentamicin call | **63** |
| **R** | **62** |
| **S** | **1** |

**PPV of a 16S methyltransferase for gentamicin-R = 62/63 = 0.984**, on independent public labels. That is
what the mechanism predicts, and it is now measured rather than asserted.

## The diagnosis is sharper than the memo's, and narrower

I expected AMRFinder's subclass filing to be *inconsistent per record*. **It is not** — it is perfectly
consistent per **gene**:

| gene | Subclass | rows | visible to the frozen rule? |
|---|---|---:|---|
| `armA` | `GENTAMICIN` | 24 / 24 | **yes — already counted** |
| `rmtB1`,`rmtC`,`rmtE`,`rmtE1`,`rmtF`,`rmtF1`,`rmtG` | `AMINOGLYCOSIDE` | 134 / 134 | no |

So the frozen rule already sees **15%** of methyltransferase rows, and **the gap is exactly the `rmt*`
family**. The accrual memo's "`rmtE1`/`rmtE`/`armA`-family" lumps in a gene that was never missing.

The candidate has been **narrowed to the measured gap** (`rmt*`/`npmA`, not `armA`). Scores are
byte-identical afterwards — which is the proof that the `armA` clause was a no-op, exactly as the per-gene
table predicts.

## The one S-labelled carrier does not implicate the candidate

`GCA_020406995.1` carries **`armA` only**, filed as `GENTAMICIN` — so **the frozen rule already calls it R**.
It is an existing false positive, and the candidate adds nothing to it.

Which leaves the honest statement: **0 S-labelled `rmt` carriers in a complete public sweep.** That is
stronger than last run (the sweep finished, and it is public data rather than this repo's local slice) but
it is still an **absence, not a bound**. 62 R and 0 S means the over-calling rate is consistent with zero
and cannot be distinguished from small.

## Two hypotheses of mine died this run

1. **"Compound subclasses naming gentamicin are a second blind spot"** — refuted last run; the rule matches
   by token, so they were already counted.
2. **"AMRFinder files `rmt` inconsistently"** — refuted here; it is perfectly consistent per gene, which is
   *why* `armA` was already covered and `rmt` never was.

Both were checked before they reached a memo.

## A defect in my own tooling, fixed

The label hunt printed a confident *"still ZERO S-labelled carriers … now a measured property"* immediately
after printing **INCOMPLETE: 4 group(s) failed**. A failed sweep's zero means *the sweep failed*, not *none
exist*. It now prints **NO CONCLUSION** on an incomplete sweep, and a test pins that. Same
reassuring-verdict pattern as the vacuous-specificity line fixed last run.

## Status

- The candidate **rescues 1 of 1** rescuable false negative, adds **0** false positives across 150 local
  isolates and the 1 public S-labelled carrier, and is now scoped to the gap that actually exists.
- **Still not deployed.** Changing the frozen surface invalidates the prospective lock and the
  reproducibility freeze; a v2 lock is a user authority call. A test asserts the frozen files carry no
  trace of the candidate.
- What would still sharpen it: an S-labelled **`rmt`** carrier. None is reachable in public data from here,
  and 62/63 suggests they are genuinely rare rather than merely unsampled.

Reproduce: `uv run python scripts/gentamicin_rmt_label_hunt.py` (network) then
`uv run python scripts/gentamicin_rmt_candidate.py` (offline).
