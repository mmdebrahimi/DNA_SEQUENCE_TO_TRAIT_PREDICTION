# A third candidate through the rejection screen — it fires a gate nothing had reached, and exposes one real limit

Run: `uv run python scripts/screen_candidate_gates.py --verify` (offline, seconds — all four hand
verdicts re-derive) or `--candidate phenosense`.

## Why a third

`dna_decode/eval/rejection_gates.py` was derived from **n=2 worked examples**, and its own memo names
that as the schema's stated limit. The evidence-surface ledger carried "screen a THIRD candidate" as an
open candidate for three weeks. The useful third is not any candidate — it is one that **reaches a path
the first three never did**.

| candidate | layer | verdict | gate that decided it |
|---|---|---|---|
| PEAR | L4 | CLEARS | — |
| HBV | L1 | REJECTED | G1 (circular label) |
| Oxford | L1 | INCOMPLETE | G7 trips (provenance not separable) |
| **PhenoSense (new)** | **L4** | **REJECTED** | **G6 — assay degeneracy** |

**No candidate had ever tripped G6**, the layer-dispatched gate that asks whether a *continuous* label is
too censored to be usable. It was built, tested on synthetic input, and never exercised by a real
candidate. Now it has been.

## The result

Screening the Stanford PhenoSense fold-change set (lamivudine, the worst of its 24 fold columns) returns
**`REJECTED`, G6 `trip`: 45.2% of values at the mode over 138 distinct levels** — against bars of 25% and
20, imported from the shipped `assay_degeneracy` rather than restated.

**That reproduces, mechanically, what was measured by hand the same morning**
(`wiki/hiv_fold_censoring_audit_2026-09-22.md`), by a completely separate code path. The screen was not
tuned to agree — the packet carries only measured fields, and the verdict fell out.

Six gates return `insufficient_data` because only measurable fields were supplied. That is the screen
behaving correctly: filling them to make the table look complete would be fabricating evidence, which is
the failure the screen exists to prevent.

## The limit it exposes — measured, not asserted

**The screen treats a candidate as monolithic, and a multi-drug substrate is not.** Same dataset, same
packet, only the drug column varying:

| drug | mode-share | distinct | G6 | verdict |
|---|---|---|---|---|
| lamivudine | 0.452 | 138 | **trip** | **REJECTED** |
| nevirapine | 0.330 | 175 | **trip** | **REJECTED** |
| zidovudine | 0.137 | 190 | pass | INCOMPLETE |
| etravirine | 0.056 | 146 | pass | INCOMPLETE |
| rilpivirine | 0.051 | 102 | pass | INCOMPLETE |

So "is the PhenoSense dataset a usable L4 substrate?" **has no single answer** — it is usable for
etravirine and not for lamivudine. The schema has no notion of a candidate with per-arm answers, and
screening such a substrate once, on whichever arm you happen to pick, yields a verdict that is true of
that arm and silently presented as true of the dataset.

**Deliberately NOT fixed here.** Adding a per-arm axis changes the schema for all four committed
candidates and their pinned verdicts; doing that off one example would repeat the n=2 generalisation the
third screen exists to test. The honest move is to record it with the evidence and screen the *worst*
arm meanwhile — which is what the committed packet does, and says so.

## Honest limits

- One new candidate. n=3 is not "validated"; it is one more than two, and it bought one gate-path
  exercise plus one named limit.
- The PhenoSense packet asserts only G1/G3/G6/G9/G10. Its G2/G4/G5/G7/G8 are genuinely unmeasured, so
  this screen does **not** establish that the substrate clears those.
- A `REJECTED` here bounds the **label**, nothing else — and specifically it does not retract anything:
  the HIV cells are validated at a fold *cutoff* of 3, where censoring at 100 is immaterial. G6 rejects
  the dataset as a **continuous** target space, which is a different question from its use as a binary
  label.

## Reusable

**The valuable next worked example is the one that reaches an unexercised path, not the one that is
easiest to write up.** Three candidates had between them exercised G1, G7 and a clean pass; the fourth
was chosen because it would land on G6, and it did.
