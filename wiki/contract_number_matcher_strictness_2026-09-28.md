# The contract-number audit was strict on one side and loose on the other (2026-09-28)

WP1 of the high-VOI pass over the evidence surface. Fixes the two coverage defects the 2026-09-27 run
found and reported but deliberately left unfixed — and, in the course of fixing them, **corrects that
run's own headline**, which mis-apportioned the blame between two real defects.

## What was wrong

`scripts/contract_number_audit.py` checks that a number in a cell contract appears in the artifact the
contract cites. The **prose** side was always boundary-guarded (`DECIMAL_RE` carries `(?<![\w.])` and a
trailing guard). The **artifact** side was a raw substring test — `any(v in haystack for v in variants)`.

So a short token matched *inside* a longer number and counted as provenance:

| cited | "matched" | actually |
|---|---|---|
| `2.1` | inside `72.1` | a coverage metric in a pneumo artifact |
| `0.51` (variant of `0.5115`) | inside `0.5189` | a different protein's Spearman |
| `0.9` (variant of `0.900`) | inside `0.9217877…` | a different accuracy |

Separately, the trailing guard `(?![\d.])` did two jobs: reject a dotted continuation (`1.2.3`, a version
string) and — accidentally — **hide any number sitting before a sentence-final period**. Five numbers were
invisible to the audit entirely.

## The correction to the previous run's headline

The 2026-09-27 run replaced a 10-artifact sampled decoy control with a deterministic full-pool scan and
reported: *"the sample was underpowered and graded the worst cell best — `essentiality:any:essentiality`
matched 0 of 10 while 44 of 600 unrelated artifacts contain all its numbers."*

**Both of those figures came from the loose matcher.** Under whole-number matching the true count is
**5 of 600 (0.83%)**, and the cell correctly grades **HIGH**. The 10-sample's original `HIGH` verdict was
*right*. "The control graded the worst cell best" was an artifact of the matcher, not of the sample size.

Holding the matcher constant at STRICT and re-running the k=10 seeded sample against the full pool
isolates the two errors:

- **k=10 vs full pool agree on 14 of 18 cell grades.** All four disagreements sit near a band boundary
  (`cyp4f2` 0/10 HIGH vs 7% MODERATE · `hla:b5701` 0/10 HIGH vs 6% MODERATE · `essentiality` 1/10 MODERATE
  vs 3% HIGH · `forward` 1/10 MODERATE vs 4% HIGH). Sampling error is real but **modest** — it is boundary
  noise, not a wrong verdict.
- **The matcher was the dominant defect**, inflating every rate 4–10× and *manufacturing the entire LOW
  bucket*.

Two real defects pointed in opposite directions, and the previous run attributed both to sampling.

## Measured effect of whole-number matching

| cell | loose | strict | grade |
|---|---|---|---|
| `pgx:human:ugt1a1` | 249/299 (83%) | 64/299 (21%) | MODERATE |
| `pgx:human:cyp4f2` | 204/299 (68%) | 20/299 (7%) | MODERATE |
| `pgx:human:slco1b1` | 392/599 (65%) | 114/599 (19%) | MODERATE |
| `typing:Escherichia_coli:serotype` | 131/199 (66%) | 16/199 (8%) | MODERATE |
| `finder:Escherichia_coli:pointfinder` | 354/599 (59%) | 175/599 (29%) | MODERATE |
| `typing:Klebsiella:ktype` | 307/599 (51%) | **29/599 (5%)** | **HIGH** |
| `typing:bacteriophage:phage` | 245/599 (41%) | 22/599 (4%) | HIGH |
| `finder:any:forward` | 72/299 (24%) | 8/199 (4%) | HIGH |
| `typing:bacteria:mlst` | 6/599 (1%) | **0/599** | HIGH |
| `finder:bacteria:resfinder` | 1/599 (0%) | **0/599** | HIGH |

**LOW-discrimination cells: 6 → 0.** Grades are now HIGH 12 / MODERATE 6 / n-a 1. A citation cannot be
judged weak on evidence that one of its numbers appeared inside a different, longer number.

**Six cells reach exactly zero** (`mlst`, `resfinder`, `flowering`, `inverse`, `klebsiella:kleb`,
`salmserovar`) — including `salmserovar`, which carries thirty measured figures.

### A band justification that has to be withdrawn

The 2026-09-27 run moved the HIGH band off `rate == 0` to `< 5%`, on the stated grounds that *no cell
reaches zero, not even salmserovar*. **That measurement was taken with the loose matcher and the stated
reason is wrong** — six cells reach zero once matching is strict.

The band **stays** at `< 5%`, on a different and now-measured argument: it separates the 0–5% group from
the 5–50% group, whereas `rate == 0` would file the dog cells at 0.3% beside `pointfinder` at 29%. The
exact `decoy_match_rate` ships per cell either way, so no reader is limited to the band. Both bars remain
**asserted, not derived**, and the code says so.

## The 4 substring-only numbers, each given the provenance it actually has

Not a blanket exemption — the four split three ways, and the split is the point:

| cell | number | true status | action |
|---|---|---|---|
| `finder:any:forward` | `0.5115` | real, in an **uncited sibling** artifact | added the citation `wiki/forward_inverse_sweep_2026-07-17.{md,json}` (bounded-verified) |
| `typing:bacteriophage:phage` | `0.862` | **derived** — prose states `25/29 called = 0.862` | new `derived` kind |
| `typing:bacteriophage:phage` | `0.291` | **derived** — prose states `25/86=0.291` | new `derived` kind |
| `typing:Salmonella:salmserovar` | `0.900` | **derived** — prose states 200 isolates and `abstention 10.0%`, so 180/200 | new `derived` kind |

**A `derived` reason must carry the arithmetic.** "It's derived" with no numerator and denominator is an
unaudited exemption; with them a reader re-does the division and the claim is checkable *without any
artifact at all*. A test re-computes each stated fraction and asserts it rounds to the cited value, so a
wrong arithmetic claim fails rather than sits there.

**The verdict stays `NOTHING_TO_ADJUDICATE` — legitimately.** The fix was expected to flip it to
`ADJUDICATION_REQUIRED` with 4 candidates; once each number is given the provenance it genuinely has,
zero drift remains. A clean verdict earned this way is the goal; a clean verdict from loose matching was
not.

## Standing figures

`NOTHING_TO_ADJUDICATE 0/180 across 19 cells` · **196** numbers extracted (was 192; +4 formerly hidden
by the lookahead, all of which match) · 9 still `unverifiable` in 3 cells · **0** LOW-discrimination.

Per-kind: `artifact` 171 · `unverifiable` 9 · `structural-non-measurement` 5 · `threshold` 4 · `derived` 3
· `enforced-by-test` 2 · `external-reference` 1 · `superseded-value` 1. `audit()` raises unless these sum
to 196.

## Performance, because correctness that costs 20× is a real cost

Whole-number matching sends far more numbers down the `_matches_by_rounding` fallback, which linearly
scanned every numeric token in the artifact. The audit went **14s → 7m46s**. Two guessed optimizations
(a substring pre-filter; removing a per-trial list copy) each changed **nothing** — so it was profiled
instead of guessed at a third time: **7 regex searches consumed 0.814 of 1.031 seconds**, 0.116s each,
because a lookbehind forces a full scan of a multi-megabyte artifact.

`_bounded_in` is now a `str.find` loop checking boundary characters directly — C-level, O(1) per
occurrence. `_artifact_numbers` returns a sorted list and `_matches_by_rounding` bisects three ranges
instead of scanning. **Back to 19.8s.** The find loop was verified equivalent to the regex across 560
(artifact, token) pairs on the real corpus: **0 disagreements**.

## Limits

- **Whole-number matching is still not semantic.** A number can now only match as a whole number, but
  nothing checks it means the same *quantity* — a cited sensitivity may match a specificity of equal value.
- **The `derived` and declaration tables are human-written.** Each entry carries a reason and the
  arithmetic is machine-checked, but the *claim that a number is derived* is judgment, capped at 5 entries
  per table so the tables stay auditable.
- **9 numbers in 3 cells remain unverifiable** (`cyp2d6`, `pneumoserotype`, `pigment`) — unchanged by this
  work; two of them are vacuous-citation cases by measurement, one is a derived no-call rate.
- **`0/180` is not "every contract number is verified."** It covers cells whose prose names a `wiki/`
  artifact, and it verifies presence, not correctness of the measurement.
