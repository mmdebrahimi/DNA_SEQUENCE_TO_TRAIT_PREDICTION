# The salmserovar 0.705 collision is resolved — and one artifact holds the same fraction as two different quantities (2026-10-03)

**Verdict: `MISMATCH 1 → 0`, `field_confirmed_by_binding 5 → 6`, zero binding defects.** The last open
lead in `scripts/contract_number_semantics.py` is closed.

**Credit where it belongs: the DIAGNOSIS was already correct and recorded.** Evidence-surface ledger
row 35 (2026-09-28) names the exact field and the exact cause:

> salmserovar coverage 0.705 is a CITATION-COVERAGE lead — the coverage genuinely exists at
> `salmserovar_threshold_tradeoff:selective_classification.deployed.coverage` but that artifact is not
> cited, while the cited ones carry 0.705 as `ours.accuracy`, incidentally confirming that cell's
> accuracy/coverage 0.705 collision is a real coincidence.

What was missing was the **fix**. This run independently reproduced that diagnosis and then closed it.

## The collision is real, and one block carries both readings

`wiki/salmserovar_threshold_tradeoff_2026-09-04.json`, `selective_classification`:

```
deployed: { correct:  99, wrong: 42, abstain: 59, n: 200, coverage: 0.705, accuracy_forced_call: 0.495 }
relaxed:  { correct: 141, wrong: 39, abstain: 20, n: 200, coverage: 0.900, accuracy_forced_call: 0.705 }
```

**0.705 appears twice, as two different quantities, and both are 141/200:**

- `deployed.coverage` — 99 + 42 = **141 of 200 GOT A CALL** (before the threshold fix)
- `relaxed.accuracy_forced_call` — **141 of 200 were CORRECT** (after it)

Same fraction; the numerator means something different each time. The contract quotes both: `Coverage
0.705->0.900` and, separately, `ours 0.7050 (141/39/20)`.

## What was actually wrong

The parent contract-number audit reported `0.705` as **`artifact`-verified** — against
`ours.accuracy` in three artifacts. The prose's claim at that token is **coverage**. Right value,
wrong quantity, reported as verified. Its own sibling in the same sentence, `0.900`, was correctly
typed `derived` with stated arithmetic, so the pair was internally inconsistent.

The coverage field existed all along. It sat in an artifact the prose **never cited** — so the
designed resolution (a `metric_bindings` pointer) could not resolve against it either.

## The fix needs BOTH halves, and that is pinned

1. **Cite** `wiki/salmserovar_threshold_tradeoff_2026-09-04.json` in the prose (it holds that whole
   paragraph — 36 rescued, net +35, coverage 0.705→0.900 — and was uncited).
2. **Bind** `0.705 → selective_classification.deployed.coverage`, quantity `coverage`.

Citing without binding leaves the audit matching an accuracy; binding without citing cannot resolve
(measured: the binding-only attempt reported `binding_declared_but_path_unresolved`). Both halves are
pinned by test.

`0.705` and `0.7050` stay **separate tokens with different quantities** — merging them, e.g. by
normalising a trailing zero, would let an accuracy verify a coverage claim, which is the exact error
this check exists to catch.

## A second defect, found because my own run caused the harm

`contract_number_semantics.py` wrote its output to the **literal** path
`wiki/contract_number_semantics_2026-09-28.json`. So every later run overwrote a file whose filename
asserts a date it no longer holds — and this is not hypothetical: **my run replaced the artifact that
ledger row 35 cites as its evidence.**

| | committed 09-28 | after my run |
|---|---|---|
| `n_numbers` | 185 | 190 |
| `n_mismatch` | **1** | 0 |
| `n_bindings_declared` | 5 | 6 |

The 09-28 file was **restored from git**, the writer now date-stamps (as its sibling
`contract_number_audit.py` already did), and `vocab_expansion_probe.py` — which read the same pinned
name — now resolves the **newest** artifact instead. Because these artifacts now *accrue* rather than
overwrite, a test asserts every one of them stays inside the parent audit's own-family circularity
exclusion, or the audit would start grading citations against its own output.

## Tests

7 new in `tests/test_contract_number_semantics.py`, plus **3 pre-existing exact-count anchors
rewritten** — `n_bindings_declared == 5`, `measurable_by_binding == 5`, `measurable_by_heuristic == 18`.
Adding one legitimate binding broke all three. Same class as the 10 rewritten yesterday.

The heuristic count moving 18 → 17 while the total held at 23 is worth keeping: `field_mismatch`
already counted as measurable, so binding the number **reclassified** it rather than converting a new
one. Asserting the literal would have read a clean reclassification as a regression, so the tests now
pin the conservation law (`n_measurable == by_binding + by_heuristic`) and zero binding defects.

**Non-vacuity proven by injection:** removing the binding turns 2 tests red. One of them originally
died with `StopIteration` on a bare `next()` — a crash dressed as a finding — and now asserts first.

## Honest limits

- This resolves **one** lead. The checker remains **weakly powered by construction**: 147 of 190
  numbers carry no quantity word near them and are `label_unlabeled`.
- A binding is a **provenance pointer**, not a metric. It makes the prose checkable; it does not
  establish that the sentence faithfully summarises the artifact.
- `field_confirmed_by_binding` is still only 6 of 190 numbers. Coverage is linear in authoring effort.
- The diagnosis was not mine — only the fix and the date-stamp defect were.
