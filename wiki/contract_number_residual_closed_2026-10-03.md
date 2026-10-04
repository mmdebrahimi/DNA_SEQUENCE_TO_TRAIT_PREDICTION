# The contract-number residual is CLOSED — and the candidate that asked for a new exemption kind rested on a retired matcher (2026-10-03)

**Verdict: `n_numbers_unverifiable` 4 → 0.** Both remaining cells now cite the artifact that holds
their numbers. **No new provenance kind was added** — the work the ledger queued turned out to be
unnecessary, because its founding measurement had been invalidated three days earlier and nobody
re-derived it.

## What the queued candidate asked for

`project_state/evidence-surface-2026-08-31.md` carried this as its next action:

> cyp2d6 (`0.62`,`1.0`) + pigment (`0.468`,`1.0`) are the ONLY remaining residual (4 numbers, 2
> cells) and are unverifiable **BY MEASUREMENT** — **730 and 737 of ~1,200 wiki files contain them**,
> so a citation would be TRUE and VACUOUS. Still open: an explicit `citation-would-be-vacuous` kind
> would stop the audit's own legend … from reading as "go cite it", which is FALSE for these.

The reasoning is sound *given* 730/1200 and 737/1200. A number present in 61% of the corpus cannot
be evidence of anything, so a citation really would be theatre, and a kind that says so is the right
fix.

## The premise is false, and the cause is already in the repo's record

Re-measured under the **whole-number matcher that replaced the loose substring one on 2026-09-28**:

| cell | numbers | artifacts containing ALL of them | claimed |
|---|---|---|---|
| `pgx:human:cyp2d6` | `0.62`, `1.0` | **13 of 1,293 (1.0%)** | 730 of ~1,200 |
| `typing:human:pigment` | `0.468`, `1.0` | **4 of 1,293 (0.3%)** | 737 of ~1,200 |

The 730/737 figures were produced by the matcher where `v in haystack` let `1.0` match inside `1.04`,
`21.0` and `0.100` — so nearly every artifact "contained" it. That defect is already documented, and
its own correction notice says the matcher *"inflated every rate 4-10x and MANUFACTURED the entire LOW
bucket"*. **This candidate was one more thing it manufactured, and it survived the fix because it had
been copied out of the measurement into prose.**

So these numbers were never unverifiable-by-measurement. They were **uncited**.

## Both citations, and what carries each

**`pgx:human:cyp2d6` → `wiki/cyp2d6_hybrid_identity_2026-07-06.md`**, which states verbatim:

> SNP diplotype 46/47 · copy-number \*5/\*xN 26/26 · **hybrid-presence sens 0.62/spec 1.0** · …

**`typing:human:pigment` → `wiki/pigment_1000g_population_2026-07-29.md`**, whose table carries
`EUR P(blue) = 0.4676` and `EAS P(brown) = 1.0`. The contract's `0.468` is **0.4676 ROUNDED**, so it
matches only through the audit's `_matches_by_rounding` path — and the contract now says so in words,
because telling a reader a value is "in the artifact" when what is there is `0.4676` would be a small
lie of exactly the kind this audit exists to catch.

## I mis-cited cyp2d6 first, and the semantic test caught it

The first citation pointed at `wiki/cyp2d6_structural_2026-07-06.md`. That file **does** contain
`0.62` — in these rows:

```
| NA18855 | `*1/*5` | 0.62 | 1 | deletion | deletion | OK |
| NA18992 | `*1/*5` | 0.62 | 1 | deletion | deletion | OK |
```

It is a **per-sample read-depth ratio**, not a sensitivity. The right VALUE in the WRONG QUANTITY —
the audit's own documented open limitation (*"whole-number matching is not SEMANTIC"*), hit live. The
citation would have passed the audit while pointing at the wrong measurement. It is now pinned by
`test_the_cyp2d6_near_miss_is_recorded_a_value_match_in_the_WRONG_quantity`, which asserts both halves:
the near-miss artifact really carries `0.62`, and really does not carry the sens/spec claim.

## The decoy control, reported rather than rounded up

The single-file containment rates above are **not** the audit's control. Its decoy trials
**concatenate several artifacts**, which makes a full match much easier:

| cell | `decoy_match_rate` | grade |
|---|---|---|
| `typing:human:pigment` | **0.0049** | HIGH |
| `pgx:human:cyp2d6` | **0.2642** | **MODERATE** |

Both clear the standing citation rule (resolve + numbers present + discrimination beats LOW), but
**cyp2d6's citation rests on the hand-read verbatim string, not on its decoy rate** — at 26% a decoy
group containing both `0.62` and `1.0` is unremarkable. Said plainly rather than quoting the 1.0%
single-file figure, which would flatter it.

## Before / after

| | 2026-09-28 (committed) | 2026-10-03 |
|---|---|---|
| `n_numbers_unverifiable` | 4 | **0** |
| `cells_unverifiable` | `{cyp2d6, pigment}` | **`{}`** |
| `n_cells_audited` | 20 | 22 |
| `n_numbers_checked` | 185 | 190 |
| `n_candidate_drift` | 0 | **0** |
| `verdict` | NOTHING_TO_ADJUDICATE | NOTHING_TO_ADJUDICATE |

The +5 checked numbers are exactly these two cells' numbers entering the checked pool now that they
are cited; `cells removed` is empty, so nothing lost coverage.

**A stale figure in `CLAUDE.md` is corrected in passing:** it records "9 still unverifiable in 3
cells", while the committed 09-28 artifact says 4 in 2 cells. Same lesson one layer up — re-derive
from the artifact, not from the prose about it.

## Tests

6 new, and **10 pre-existing tests rewritten rather than number-bumped**. They had pinned exact
live-corpus counts (`n_numbers_checked == 185`, `n_numbers_extracted == 196`, `n_cells_audited == 20`)
and fired on a legitimate citation — the documented exact-anchor-against-a-live-corpus trap. Each now
asserts the invariant it was protecting (totals reconcile with the per-cell rows; per-kind counts sum
to the extracted total) with the historical numbers kept in the docstrings as provenance.

One was renamed: `test_the_cells_the_rule_refused_still_have_no_citation` →
`test_no_cell_is_refused_any_more_and_the_refusal_reason_was_RETRACTED`. Bumping it would have left a
green test asserting a claim the repo has retracted.

`test_the_disclosure_block_is_a_TRIPWIRE_now_that_the_residual_is_empty` replaces an
`n_numbers_unverifiable > 0` assertion that now fails **by achievement**. Flipping it to `== 0` would
pass forever without exercising the block again, so it is a tripwire instead: an empty residual must
be reported consistently, and any cell re-entering it must appear in both places.

**Non-vacuity proven by injection:** removing the pigment citation turns 3 distinct tests red.

## Honest limits

- This closes CITATION coverage. It does **not** show either measurement is sound — the audit verifies
  **presence**, and the semantic step was my own reading of two files.
- `0.468` is reached only by rounding. Disclosed, not hidden.
- cyp2d6 is MODERATE-discrimination; its citation is carried by the verbatim string.
- The residual reaching 0 is a property of today's corpus and contracts, not a permanent invariant —
  hence the tripwire rather than a frozen zero.
- Frozen AMR surface untouched; this is registry prose + an audit script's tests only.
