# CABBAGE scored cell — *E. coli* × gentamicin, N=7,948, and the `rmt` specificity question is answered

**2026-10-11.** The v2 gentamicin lock's own honest limit read: *"zero S-labelled `rmt` carriers exist
in any dataset checked, so the unchanged specificity is an ABSENCE of counter-examples, not a bound on
over-calling."* This cohort contains counter-examples, so that limit is now superseded by measurement.
Artifacts: `wiki/cabbage_scored_cell_escherichia_gentamicin_2026-10-11.json`,
`wiki/cabbage_fn_anatomy_escherichia_gentamicin_2026-10-11.json`.
Reproduce: `uv run python scripts/cabbage_score_cell.py --species "Escherichia coli" --drug gentamicin --organism Escherichia_coli_Shigella --max-isolates 0`
then `uv run python scripts/cabbage_fn_anatomy.py --species "Escherichia coli" --drug gentamicin --organism Escherichia_coli_Shigella`.

## Result

Frozen `eval/amr_rules.call_resistance` (byte-unchanged), labels from portal MIC through this repo's
`external_mic_labels.tier_for_isolate`, STRICT tiers, operator-aware; fail-closed leakage gate excluded
260 prior accessions.

| | |
|---|---|
| with assembly + MIC | 13,032 |
| strict-labelled | 8,259 (67 R / 8,192 S) |
| scoreable after leakage gate | 7,999 |
| **scored** | **7,948** |
| confusion | tp 55 / fn 6 / tn 7,848 / fp 39 |
| **accuracy** | **0.9943** |
| **sensitivity** | **0.9016** (55/61) |
| **specificity** | **0.9951** (7,848/7,887) |
| fetch failures | 51 (0.6%) |

**This cohort PASSES this repo's own source-diversity bar — verdict `source_diverse`.** On the SCORED
set: **21 countries, largest 0.5844** (USA); **10 database labels, largest 0.3534** (NCBI_antibiogram);
collection years 2000–2024. The bar is ≥5 sources AND largest ≤0.60, so country clears it by 0.016 —
narrowly, and worth saying so rather than rounding it into comfort. (The label census before scoring
reads 23 countries / 0.580 and 14 labels / 0.353; the scored-set figures above are the ones that
describe the measurement.) This is rare here and is the main reason the cell is worth more than its
size: the Campylobacter CABBAGE cell fails the same bar at 96.4% one country.

**The v2 rescue replicates.** The v2 rule was justified on N=131 leakage-gated isolates at sens 0.523 →
0.892. On this independent, ~60× larger, source-diverse cohort sens is **0.9016** — the same number,
measured again on different isolates.

## The `rmt` question: counter-examples exist, and the rescue is still the best-performing rule component

Per-determinant PPV over the whole scored set (marginal, see limits):

| determinant | AMRFinder Subclass | counted by | carriers | R | S | PPV | Wilson CI |
|---|---|---|---|---|---|---|---|
| `aac(3)-IId` | GENTAMICIN | v1 | 40 | 20 | 20 | **0.500** | [0.352, 0.648] |
| `rmtB1` | AMINOGLYCOSIDE | **v2 rescue only** | 22 | 20 | 2 | **0.909** | [0.722, 0.975] |
| `aac(3)-IVa` | GENTAMICIN | v1 | 18 | 15 | 3 | 0.833 | [0.608, 0.942] |
| `aac(3)-IIe` | GENTAMICIN | v1 | 9 | 6 | 3 | 0.667 | [0.354, 0.879] |
| `armA` | GENTAMICIN | v1 | 8 | 8 | 0 | 1.000 | [0.676, 1.0] |
| `aac(3)-VIa` | GENTAMICIN | v1 | 8 | 0 | 8 | **0.000** | [0.0, 0.324] |
| `aac(6')-Ib` | GENTAMICIN | v1 | 2 | 0 | 2 | 0.000 | [0.0, 0.658] |
| `aac(3)-IV` | GENTAMICIN | v1 | 1 | 1 | 0 | 1.000 | [0.207, 1.0] |
| `ant(2'')-Ia` | GENTAMICIN | v1 | 1 | 1 | 0 | 1.000 | [0.207, 1.0] |
| `aac(6')-Ib4` | GENTAMICIN | v1 | 1 | 0 | 1 | 0.000 | [0.0, 0.793] |

**Two S-labelled `rmtB1` carriers now exist, so the "absence of counter-examples" caveat is discharged —
and the measured PPV is 0.909, the highest of any multi-carrier determinant in the set.** `armA`, which
v1 already counted because AMRFinder files it under Subclass GENTAMICIN, is 8/8.

**The over-call is NOT the v2 widening, and that is the finding a reader would most likely get backwards.**
Attribution verified per isolate through the deployed rule rather than inferred: `aac(3)-IId` and
`aac(3)-VIa` both carry Subclass **GENTAMICIN** and are matched by v1's `subclass_any={GENTAMICIN}`,
while `rmtB1` carries Subclass **AMINOGLYCOSIDE** and is reachable *only* through v2's
`symbol_rescue=^(rmt[A-H]\d*|npmA\d*)$`. Splitting the 39 false positives by which rule version admits
them:

- **37 of 39 come from v1-counted determinants** — `aac(3)-IId` 20, `aac(3)-VIa` 8, `aac(3)-IVa` 3,
  `aac(3)-IIe` 3, `aac(6')-Ib` 2, `aac(6')-Ib4` 1.
- **2 of 39 come from the v2 rescue** (`rmtB1`).

So the widening that reset the lock contributes 5% of the over-call while fixing the blind spot it was
built for. The precision problem that remains is `aac(3)-IId` at PPV 0.500 on 40 carriers (CI excludes
0.9) and `aac(3)-VIa` at 0/8 — determinants v1 counted all along.

**All 6 false negatives carry ZERO counted determinants**, the same declared blind spot as every other
cell measured this way (`undetectable_mechanisms` is stated at call time).

**Acting on any of this is a USER AUTHORITY call.** Changing what the rule counts edits the sha256-pinned
frozen surface and retires the active 2026-08-31 lock. Nothing was changed; the frozen five are
byte-unchanged and `verify_lock` re-verified ok.

## Honest limits

- **PPV is MARGINAL per determinant, not causal.** A carrier may hold several determinants, so credit or
  blame is not attributable to one symbol from these counts alone — `aac(3)-IId`'s 20 S carriers may
  co-carry something else, and conversely.
- **The R class is small and that is what bounds sensitivity.** 61 scored R against 7,887 S: sens 0.9016
  rests on 61 isolates while spec 0.9951 rests on 7,887. Read the two with very different confidence.
- **Label censoring does real work here.** Of 13,032 isolates with MIC, 3,706 are `CENSORED_EXCLUDED`
  and 3,132 `CENSORED_HIGH_S`; the R class is 15 `HIGH_R` + 52 `CENSORED_HIGH_R`. Operator-aware
  censoring is what makes those 52 usable rather than discarded or miscalled.
- **NOT clonality-corrected** (G8), as with every CABBAGE cell —
  `cabbage-clonality-g8-blocked-no-mash`. Isolate-level counts; *E. coli* collections are clonally
  structured, and 23 countries bounds source concentration, not lineage redundancy.
- **Source-diverse on the metadata axes measured** (country, database label), and the country axis
  clears the 0.60 bar at 0.5844 — a margin of 0.016, so a different cohort draw could fail it.
  `database` remains an overlapping aggregator label rather than a BioProject, so this is not a
  per-BioProject diversity claim.
- **The portal's AMRFinderPlus run is not under our control.** Its organism setting is EBI's; the
  Campylobacter cell shows the portal does run organism-aware (point mutations appear), but for
  gentamicin the rule keys on acquired genes rather than point mutations, so this matters less here
  than it would for a target-site drug.
- 51 fetch failures (0.6%) are excluded and counted, never scored wrong — a far lower rate than the
  Campylobacter cell's 8.3%, which is a datapoint toward the registered
  `cabbage-genomes-tree-coverage-uncensused` constraint, not a census.
- Not temporally locked: these isolates are **not** prospective relative to the active 2026-08-31 lock.
