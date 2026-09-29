# The audit's variant matcher TRUNCATED, and one artifact token was verifying two different claims

**2026-09-28.** `scripts/contract_number_audit.py` · `_number_variants` · third and largest instalment of
one defect class.

## What was wrong

`_number_variants(tok)` expands a cited number into forms an artifact might legitimately store it in
(zero-padding, percent↔fraction). It generated those with `f"{f:.{nd}f}"` for `nd` in 1..4, which **pads**
above the token's own precision and **truncates** below it:

```
_number_variants("0.939")  ->  ['0.9', '0.939', '0.9390', '0.94', '93.9', '94']
_number_variants("2.13")   ->  ['0.0213', '2.1', '2.13', '2.130', '2.1300']
```

It feeds the **fast path**, so a truncated hit was recorded as an exact match and never reached
`_matches_by_rounding` — the branch that exists to disclose exactly this. `matched_only_after_rounding`,
the disclosure field, stayed **empty** for every one of them.

## The case that settled it

`pgx:human:ugt1a1` cites **both** `0.939` and `0.941`. The single artifact token **`0.9`** was matching
**both**. One artifact value cannot be evidence for two different claims — so the check was not
discriminating between the numbers at all.

Others found the same way: `finder:any:forward`'s `0.454` passing on **`0.5`**, its `0.926` on **`93`**,
`typing:Escherichia_coli:serotype`'s `0.167` on **`17`**.

## Measured before fixing

**7 of 166 verified numbers passed ONLY on a truncation.** The blast radius on the *decoy control* is far
larger, because the decoy haystack is hundreds of unrelated artifacts:

| pinned-pool decoy matches | before | after |
|---|---|---|
| `finder:Escherichia_coli:pointfinder` | **175 / 599** | **0 / 599** |
| `typing:Klebsiella:ktype` | **29 / 599** | **0 / 599** |

pointfinder cites `0.9923`, whose truncating variants included **`1.0`** — so it matched any artifact
carrying that token. **The 2026-09-27 boundary-strictness fix took pointfinder 59% → 29% and reported 29%
as the corrected state; 29% was still almost entirely artifact.**

## The direction is the whole point

`_matches_by_rounding` already had it right: the **artifact** must round to the **cited** value at the
cited precision, so `0.9` correctly fails against `0.939`. Truncating the cited value down to the
artifact's coarser one is that same test run **backwards**. Dropping the truncations does not lose
legitimate rounding — it **routes it to the disclosing path**, where a reader can see it happened.

`_number_variants` now filters to value-preserving forms only (same number, optionally rescaled by 100).

## Result

| | before | after |
|---|---|---|
| candidate drift | 0 | **0** |
| `matched_only_after_rounding` (disclosed) | **0** for every affected cell | **19** |
| `unverifiable` numbers | 9 | **4** |
| cells with no citation at all | 3 | **2** |

Two genuine drift candidates surfaced and were resolved, not suppressed:

- **`finder:any:forward` `0.926`** — a real measured number (the HIV curated-catalogue AUROC) living in an
  artifact the cell did not cite. Citation added: `wiki/hiv_esm_vs_catalog_2026-07-09.md` carries both it
  and the paired `0.454`.
- **`typing:Streptococcus_pneumoniae:pneumoserotype` `2.13`** — mine, introduced in the same session while
  documenting a derivation. No artifact carries it; removed from the prose.

## The sibling checker was already correct

`scripts/contract_number_semantics.py::value_matches` does `round(value, decimals) == float(token)` —
artifact rounded to cited precision, the right direction. The two checkers disagreed and the **newer** one
was right. Nothing changed there.

## What the fix cost, reported not buried

**The discrimination GRADE has saturated at `HIGH` for every cell**, because the bands (`<5%` / `<50%`) were
calibrated against the inflated rates. That is the bands going stale, not the control breaking: the
underlying **rate still spans an order of magnitude** (0 → 0.0365) and still orders cells the way the
control intends — `ugt1a1` (3 round figures) highest, `salmserovar` (30 measured figures) at zero.

Inventing new bands off 20 cells would be one more asserted threshold, so **the rate is what ships** and
the grade's saturation is recorded rather than engineered around. Consequently:

- `test_discrimination_actually_discriminates` now asserts the **rate** separates, not the grade.
- The per-cell `0.0 < rate` guard was wrong once matching became strict — a 30-number fingerprint appearing
  in an unrelated artifact would be extraordinary, so **0.0 is the correct expected value for most cells**.
  Non-vacuity moved to the **pool** level (some cell must be nonzero, some must be zero, rates must vary).
- All four pinned-pool anchors became 0, which would make that test **vacuous** — four zeros cannot tell a
  correct matcher from an inert one. Two nonzero anchors (`ugt1a1` 11, `slco1b1` 10) were added for exactly
  that reason.

## Honest limits

- This bounds **presence**, never whether a measurement was sound.
- The remaining 4 `unverifiable` numbers (`cyp2d6`, `pigment`) are unverifiable **by measurement** — 730+
  wiki files contain their values, so no citation could discriminate. That is a property of the numbers,
  not a missing to-do.
- Matching is still **not semantic** on this side; a cited number can match the right VALUE at the wrong
  QUANTITY. `contract_number_semantics.py` is the check for that, and it is separate.
- Six pre-registered anchors were re-baselined. Each carries its previous value and the reason inline.

## Knock-on the same day: the vocabulary lever was not exhausted after all

Earlier that day a probe (`scripts/vocab_expansion_probe.py`) measured **zero addable `LABEL_VOCAB`
entries** across the 141 remaining unlabeled numbers, and that was recorded as the lever being exhausted.

Citing pneumoserotype's two artifacts (above) moved 5 numbers into the cited set, and **`concordance`
immediately surfaced as a clean candidate** — 2 near / 2 would-confirm / **0 would-mismatch**. Added:
`measurable_by_heuristic` 17 → 18, `label_unlabeled` 145 → 144, **0 new false leads**.

**So the vocabulary lever's headroom is a function of WHICH CELLS ARE CITED and re-opens whenever coverage
grows.** A measurement over "all remaining numbers" is a measurement over the *cited* ones — and citing was
the other thing being done that day. The probe's own test was the tripwire and it fired on the next full
suite run.

Two further results from adding it:

- **The probe is an UPPER BOUND on conversions.** It counts a word as converting when it sits in the window
  AND the artifact path matches, but `label_for` applies leading-label-wins plus inter-number clipping, so at
  most **one** word wins per number. `concordance` was estimated at 2 and delivered 1. The `would_mismatch`
  cost column is a ceiling for the same reason, which makes a marginal candidate even less attractive.
- **`concordance` is deliberately NOT an alias of `acc`, and `agreement` was refused.** Concordance is
  agreement between two *callers* where neither is truth; accuracy is measured against a *label*. That
  distinction is the basis of this repo's `FAITHFUL_TO_TOOL` vs `INDEPENDENT_MEASURED` split, so collapsing
  them would let a tool-agreement figure verify a prose claim of accuracy — the right-value/wrong-quantity
  error the semantic check exists for. Pinned in both directions. `agreement` scored 2 conversions against
  **1 false lead** and was refused, since a false lead is strictly worse than an honestly-unlabeled number.

The test that encoded the retracted ranking was **rewritten, not re-baselined**:
`test_a_vocabulary_entry_has_MORE_leverage_than_a_binding` →
`test_vocabulary_and_bindings_are_COMPLEMENTARY_not_ranked`. Adjusting its number would have left a green
test asserting a retracted claim, which is worse than a red one.
