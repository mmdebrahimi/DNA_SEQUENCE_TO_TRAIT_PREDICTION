# The ten rejection gates are runnable now — and the first thing they did was correct me twice

**`dna_decode/eval/rejection_gates.py` + `scripts/screen_candidate_gates.py`.** The gate set that has been
applied by hand since 2026-06-13 is a function. It screens a candidate dataset and **refuses a verdict**
when the parts that cannot be computed have not been supplied by a human.

`uv run python scripts/screen_candidate_gates.py --verify` re-derives both committed hand verdicts in
under a second, offline.

---

## Why now, and not earlier

This was deliberately deferred twice. Coding a judgment checklist before there is anything to screen
produces trust-layer theatre: a polished artifact that encodes opinions as booleans and teaches nobody
anything. The stated precondition was *scope it to one concrete candidate first*.

That is now over-satisfied — **two** candidates have been screened by hand and committed:
[PEAR](pear_substrate_screen_2026-08-31.md) (2026-08-31) and [HBV](hbv_cell_gate_screen_2026-09-01.md)
(2026-09-01). Two is enough to derive the schema and few enough that it stays grounded in real cases.

## The split that makes it honest: 8 mechanical, 2 judgment

Reading the ten definitions rather than remembering them, **eight already carry a countable rule** —
G4 `<20/class`, G6 majority-censored, G8 `<~3 effective lineages`, G9/G10 a majority, plus G2's
contingency, G5's fetchable-accession count and G7's field population.

**Two do not**, and no amount of code changes that:

- **G1** — is the label wet-lab/clinical, or produced by a genomic tool the decoder would compete against?
- **G3** — is the label an assay reading, or a description of where and why the isolate was collected?

Those are readings of a methods section. So they take a **human evidence string plus an explicit
assertion**, and a screen missing either returns `NEEDS_HUMAN_EVIDENCE` — at which point the overall
verdict is **REFUSED**, not "pass". Prose alone does not satisfy them either; narrative without the
assertion still refuses. *A screen that defaulted its uncomputable gates to pass would be exactly the
theatre this was deferred to avoid.*

## The layer comes first, and G6 is where that bites

The same dataset screens differently depending on what it is **for**, so `intended_layer` is required and
never inferred (`L1_AMR_RS` or `L4_forward_continuous`).

G6 is the subtle one. Its **letter** is MIC interval-censoring at a clinical breakpoint — meaningless for a
continuous fitness readout. Its **spirit** is wider: *the quantitative readout cannot separate where it
matters*. Under L4 that failure arrives as an **assay floor** — every dead variant scoring the same value.

So G6 is layer-dispatched: L1 evaluates censoring, L4 evaluates degeneracy on the shipped
`assay_degeneracy` thresholds (mode-share > 25%, fewer than 20 distinct levels), imported rather than
restated and pinned by a drift test. Returning `not_applicable` for L4 would have cleared a candidate on
**the exact failure that already bit this repo once** — CcdB posted the forward/inverse sweep's *best*
number because 79.3% of it was tied at the ceiling.

This is also why the PEAR memo was right to score G6 OPEN on a continuous assay, and my first cut of the
code was wrong to call it n/a.

## Two corrections it produced immediately

**1. A real bug in my own gate code.** `--verify` reported `G2: memo says not_applicable, screen says
insufficient_data`. Cause: G2 demanded its measurement *before* checking applicability, unlike its four
sibling gates. `insufficient_data` reads as *go measure this* — misleading when the honest answer is
*this cannot apply*. **Applicability must precede input requirement**, now pinned by test.

**2. An overstated headline in a committed memo.** The mechanical screen returns **INCOMPLETE** for PEAR,
not "clears". PEAR's own gate table and honest-limits section both record G6 as unscreened; only the
headline said it cleared everything. The memo headline is corrected. I did not tune the screen to agree
with it — the disagreement was the finding.

Neither would have surfaced without running the screen against verdicts that already existed. **A checker
with nothing to check is unfalsifiable.**

## What a PASS does and does not mean

Shipped in the result object, so it travels with every verdict:

> A PASS bounds only whether a usable LABEL exists and whether the rule is scoreable. It is NOT a build
> recommendation: artifact reachability, regime fit (`eval/regime.py`) and worth-doing are separate
> questions.

PEAR is the standing proof. It clears every applicable gate and is still not buildable here — its
processed data ships as serialized ggplot2 objects and R is not installed. **The gates bound the label
question only.**

## Honest limits

- **n=2 worked examples.** The schema is derived from two candidates, one per layer. A third could well
  need a field that does not exist yet.
- **G2's threshold is not in the memo.** The 60% dominant-source bar is imported from the
  source-concentration arm (`scripts/source_diverse_validate.py`), not from the gate definition. It is
  stated in the code so it is arguable rather than hidden.
- **Decision support, not an oracle.** The mechanical gates are only as good as the numbers fed in, and
  G1/G3 remain human judgments that the tool records rather than makes.
- **It cannot screen what it is not given.** A candidate nobody supplies evidence for returns REFUSED —
  correctly, and unhelpfully.

## Run it

```bash
uv run python scripts/screen_candidate_gates.py --verify              # reproduce both hand verdicts
uv run python scripts/screen_candidate_gates.py --candidate pear      # full per-gate JSON
uv run python scripts/screen_candidate_gates.py --packet mine.json    # screen a new candidate
```

18 tests at `tests/test_rejection_gates.py`, including a non-vacuity check that corrupting an expected
verdict makes the reproduction check fail. Frozen AMR surface untouched — this reads nothing and scores
nothing.
